function patchmeister()
%
% Created by Marcel Goldschen-Ohm <goldschen-ohm@utexas.edu, marcel.goldschen@gmail.com>

%% data structures
% Data is organized into a collection of time series data traces in one or
% more channels. The traces are arranged as a 2D struct array where each
% row is a sweep and each column is a channel. See the templates below for
% details on this data structure.
%
% Note: Struct arrays of individual traces were chosen over lumping data
% into higher dimensional arrays or cell arrays primarily because I feel
% that this provides a more intuitive interface to the data.
%
% !!! All traces in a sweep have the same x data and repeated sweeps have
% the same x/y labels and units. Despite this, each trace is assigned its
% own x and y data fields as well as their label and unit specifications.
% The disadvantage of this approach is that it is possible to group traces
% in a way that violates the idea of sequences of recorded sweeps. However,
% the advantages of this approach are that it provides an intuitive
% interface for trace manipulation (e.g. every trace knows about it's own
% data) and provides for very flexible trace grouping (e.g. can handle
% inclusion of sweeps where recording was stopped prematurely, or similar
% stimulus sequences with slightly different recording lengths).
% Furthermore, it is fairly simple to enforce that all sweeps in a group
% have the same x data, channels and units (e.g. the deal() function can be
% used to easily assign labels/units to all traces in a group). Also, this
% approach does not necessarily mean redundant copying of shared x data.
% For example, assigning the x values of one trace to other traces will
% avoid copying the x data via MATLAB's copy on write rules (this will
% remain true so long as the original raw x values are never edited, see
% below). Thus, we can have the best of both worlds, intuitive access to x
% data from within each trace and shared memory to x data for all traces in
% a group.
%
% !!! The original raw (x,y) data in each trace should almost NEVER be
% altered unless you are absolutely sure of what you are doing. Instead,
% offsets and scale factors, as well as masked, zeroed and interpolated
% regions are stored separately, and can be used to reconstruct the
% adjusted (e.g. baselined, scaled, masked, etc.) trace from the raw data.
% For such reconstruction, see the interface functions described below.
% This also means that the original raw data can always be viewed at any
% time if desired.
%
% !!! The function getXY_(...) defined below provides a consistent
% interface to the trace data, and should be used in all cases rather than
% accessing the x, y data fields directly.

% A time series trace of (x,y) data points is the basic fundamental unit
% that we will deal with in this program.
template.trace.x = []; % 1x1 sample interval OR Tx1 sample points, e.g. time
template.trace.y = []; % Tx1 sample points, e.g. current, voltage, etc.
% basic info about the (x,y) data
template.trace.xlabel = "Time";
template.trace.xunit = "s";
template.trace.ylabel = ""; % e.g. "Current"
template.trace.yunit = ""; % e.g. "A"
template.trace.timestamp = NaT; % timestamp for recorded data
% basic offsets and scaling
template.trace.x0 = 0; % 1x1 time zero offset
template.trace.y0 = 0; % 1x1 (uniform) OR Tx1 (nonuniform) baseline offset
template.trace.yscale = 1; % 1x1 (uniform) OR Tx1 (nonuniform) scale factor
% optional masking, zeroing and interpolation
template.trace.ismasked = false; % 1x1 logical flag for masking entire trace
template.trace.masked = []; % Tx1 logical for masked data points
template.trace.zeroed = []; % Tx1 logical for zeroed data points
template.trace.interpolated = []; % Tx1 logical for interpolated segments

% A series is a 2D struct array of traces where each column is a channel
% and each row is a sweep. The UI is setup to view a single series.
template.series.traces = repmat(template.trace, [0,0]); % SxC data traces, e.g. traces(sweep,channel).y
% Sweeps can optionally be grouped as desired.
template.series.groupids = []; % Sx1 group indices for each sweep
template.series.grouplabels = repmat("", [0,0]); % Gx1 group labels
% metadata
template.series.meta.date = datestr(now, 'yyyy-mm-dd');
template.series.meta.patchid = ''; % cell/patch ID
template.series.meta.construct = ''; % construct, e.g. denote subunit composition and mutations
template.series.meta.experiment = ''; % very short one-line summary of experiment (more detail can go in notes)
template.series.meta.notes = '';

%% init
% for debugging (comment out if you don't want command line access to the
% data)
global data;
global ui;

% data for a single series of sweeps
data = template.series;

% UI
ui = struct();
initUI_();

%% user init
% Put any user specific initialization code here...

%% test data
% This is just some dummy data for testing purposes.

% data.traces = repmat(template.trace, [5,2]);
% for i_ = 1:size(data.traces,1)
%     for j_ = 1:size(data.traces,2)
%         data.traces(i_,j_).x = [0:99]';
%         data.traces(i_,j_).y = rand(100,1) .* 10000;
%     end
% end
% [data.traces(:,1).ylabel] = deal("Current");
% [data.traces(:,2).ylabel] = deal("Voltage");
% [data.traces(:,1).yunit] = deal("pA");
% [data.traces(:,2).yunit] = deal("mV");
% data.groupids = [1; 1; 1; 2; 2];
% data.grouplabels = ["Group 1"; "Group 2"];
% updateUI_();

%% trace data interface
% You should almost always use these functions to access trace (x,y) data
% rather than directly referencing the trace.x and trace.y data fields.
% This provides a consistent interface to the data (see description above).

    function [x,y] = getXY_(trace, israw)
        x = [];
        y = [];
        % the trace's original raw (x,y) data
        if numel(trace.x) == 1
            x = reshape((0:length(trace.y)-1) .* trace.x, size(trace.y));
        else
            x = trace.x;
        end
        y = trace.y;
        if exist('israw', 'var') && israw
            return;
        end
        % offset and scale (x,y) data
        x = x - trace.x0;
        y = (y - trace.y0) .* trace.yscale;
        % interpolate segments
        % first and last pts of each interpolated segment are unchanged,
        % inbetween pts are linearly interpolated based on x values
        if any(trace.interpolated)
            pts = find(trace.interpolated);
            segments = breakIntoContiguousSegments_(pts);
            for i = 1:numel(segments)
                a = segments{i}(1);
                b = segments{i}(end);
                if b-a > 1
                    fracb = (x(a+1:b-1) - x(a)) ./ (x(b) - x(a));
                    y(a+1:b-1) = (1-fracb) .* y(a) + fracb .* y(b);
                end
            end
        end
        % zero segments
        if any(trace.zeroed)
            y(trace.zeroed) = 0;
        end
        % denoise
        if denoise_()
            y = wdenoise(y);
        end
        % decimate
        n = str2num(ui.decimateEdit.String);
        if ~isempty(n) && n > 1
            y = decimate(y, n);
            x = downsample(x, n);
            % set masked segments to NaN?
            masked = downsample(trace.masked, n);
            if any(masked) && ~showMasked_()
                y(masked) = nan;
            end
        else
            % set masked segments to NaN?
            if any(trace.masked) && ~showMasked_()
                y(trace.masked) = nan;
            end
        end
    end

    function [x,y] = getOverlappingXY_(traces, israw)
        % return [x, y1, y2,...] for region where all traces overlap
        % !!! traces should all have the same sample interval
        x = [];
        y = [];
        if isempty(traces); return; end
        if ~exist('israw', 'var'); israw = false; end
        [x,y] = getXY_(traces(1), israw);
        dx = min(diff(x));
        epsilon = 0.01 * dx;
        for i = 2:numel(traces)
            [xi,yi] = getXY_(traces(i), israw);
            dxi = min(diff(xi));
            epsiloni = 0.01 * dxi;
            ind = intersect(find(x >= xi(1) - epsilon), find(x <= xi(end) + epsilon));
            indi = intersect(find(xi >= x(1) - epsiloni), find(xi <= x(end) + epsiloni));
            if numel(ind) ~= numel(indi)
                error('getOverlappingXY_: Invalid overlapping indices.');
                x = [];
                y = [];
                return;
            end
            x = x(ind);
            y = y(ind,:);
            xi = xi(indi);
            yi = yi(indi);
            if any(abs(x - xi) > epsilon)
                error('getOverlappingXY_: Traces have different sample intervals.');
                x = [];
                y = [];
                return;
            end
            y = [y yi];
        end
    end

    function [x,y] = getMeanXY_(traces, israw)
        if ~exist('israw', 'var'); israw = false; end
        [x,y] = getOverlappingXY_(traces, israw);
        y = nanmean(y,2);
    end

    function segments = breakIntoContiguousSegments_(ind)
        % return cell array of numeric arrays of contiguous indices
        segments = {};
        segment = [ind(1)];
        for i = 2:numel(ind)
            if ind(i) == segment(end)+1
                segment(end+1) = ind(i);
            else
                segments{end+1} = segment;
                segment = [ind(i)];
            end
        end
        if ~isempty(segment)
            segments{end+1} = segment;
        end
    end

%% per axes visible trace manipulations
% These functions apply manipulations (e.g. offsets, scaling, masking,
% ect.) to all visible traces in a given axes or list of axes. Most of
% these manipulations require that data points within the traces be
% selected using MATLAB's brush tool prior to applying the manipulation
% (e.g. baseline requires a selected region to which to baseline).

    function autoscale_(ax, xy, same)
        if ~exist('ax', 'var') || isempty(ax)
            ax = vertcat(ui.groups.ax);
        end
        if ~exist('xy', 'var'); xy = "xy"; end
        if ~exist('same', 'var'); same = false; end
        if xy == "x"
            lims = [];
            for i = 1:numel(ax)
                lims = [lims; axesXYLims_(ax(i))];
            end
            lims = [min(lims(:,1)), max(lims(:,2))];
            for i = 1:numel(ax)
                ax(i).XLim = lims;
            end
        elseif xy == "y"
            lims = [];
            for i = 1:numel(ax)
                lims = [lims; axesXYLims_(ax(i))];
            end
            if same
                i = ~isinf(lims(:,3));
                if any(i)
                    lims(:,3) = min(lims(i,3));
                end
                i = ~isinf(lims(:,4));
                if any(i)
                    lims(:,4) = max(lims(i,4));
                end
            end
            for i = 1:numel(ax)
                ax(i).YLim = lims(i,3:4);
            end
        elseif xy == "xy"
            lims = [];
            for i = 1:numel(ax)
                lims = [lims; axesXYLims_(ax(i))];
            end
            i = ~isinf(lims(:,1));
            if any(i)
                lims(:,1) = min(lims(i,1));
            end
            i = ~isinf(lims(:,2));
            if any(i)
                lims(:,2) = max(lims(i,2));
            end
            if same
                i = ~isinf(lims(:,3));
                if any(i)
                    lims(:,3) = min(lims(i,3));
                end
                i = ~isinf(lims(:,4));
                if any(i)
                    lims(:,4) = max(lims(i,4));
                end
            end
            for i = 1:numel(ax)
                axis(ax(i), lims(i,:));
            end
        end
    end

    function xylims = axesXYLims_(ax)
        xylims = [-inf, inf, -inf, inf];
        for i = 1:numel(ax.Children)
            if ax.Children(i).Type == "line"
                xmin = min(ax.Children(i).XData);
                xmax = max(ax.Children(i).XData);
                ymin = min(ax.Children(i).YData);
                ymax = max(ax.Children(i).YData);
                if isinf(xylims)
                    xylims = [xmin, xmax, ymin, ymax];
                else
                    xylims = [xylims; xmin, xmax, ymin, ymax];
                end
            end
        end
        if size(xylims, 1) > 1
            xylims(1,1) = min(xylims(:,1));
            xylims(1,2) = max(xylims(:,2));
            xylims(1,3) = min(xylims(:,3));
            xylims(1,4) = max(xylims(:,4));
            xylims(2:end,:) = [];
        end
        if ~isinf(xylims)
            % add margins around y data
            xylims(3) = xylims(3) - 0.05 * diff(xylims(3:4));
            xylims(4) = xylims(4) + 0.05 * diff(xylims(3:4));
        end
    end

    function setBaselineVisibleTraces_(ax, y0)
        if ~exist('ax', 'var') || isempty(ax)
            ax = vertcat(ui.groups.ax);
        end
        if ~exist('y0', 'var')
            answer = inputdlg('Baseline:', 'Set Baseline', 1, {'0'});
            if isempty(answer); return; end
            y0 = str2num(answer{1});
        end
        israw = showRaw_();
        for i = 1:numel(ax)
            rows = ax(i).UserData.rows;
            channel = ax(i).UserData.channel;
            if israw
                % set baseline
                [data.traces(rows,channel).y0] = deal(y0);
            else
                % shift baseline (i.e. set based on current view)
                for j = 1:numel(rows)
                    data.traces(rows(j),channel).y0 = ...
                        data.traces(rows(j),channel).y0 - y0;
                end
            end
        end
        redraw_();
    end

    function baselineFlatVisibleTraces_(ax)
        if ~exist('ax', 'var') || isempty(ax)
            ax = vertcat(ui.groups.ax);
        end
        israw = showRaw_();
        didit = false;
        for i = 1:numel(ax)
            for j = 1:numel(ax(i).UserData.traces)
                brushdata = logical(ax(i).UserData.traces(j).BrushData);
                if any(brushdata)
                    row = ax(i).UserData.rows(j);
                    channel = ax(i).UserData.channel;
                    trace = data.traces(row,channel);
                    [x,y] = getXY_(trace, israw);
                    if israw
                        data.traces(row,channel).y0 = mean(y(brushdata));
                    else
                        data.traces(row,channel).y0 = ...
                            trace.y0 + mean(y(brushdata));
                    end
                    didit = true;
                end
            end
        end
        if didit
            redraw_();
        else
            msgbox('Select data to baseline to with the brush tool.', ...
                'Baseline Flat');
        end
    end

    function baselineSlopingVisibleTraces_(ax)
        if ~exist('ax', 'var') || isempty(ax)
            ax = vertcat(ui.groups.ax);
        end
        israw = showRaw_();
        didit = false;
        for i = 1:numel(ax)
            for j = 1:numel(ax(i).UserData.traces)
                brushdata = logical(ax(i).UserData.traces(j).BrushData);
                if any(brushdata)
                    selpts = find(brushdata);
                    dsel = find(diff(selpts) > 1);
                    if numel(dsel) == 1
                        selpts1 = selpts(1:dsel);
                        selpts2 = selpts(dsel+1:end);
                        row = ax(i).UserData.rows(j);
                        channel = ax(i).UserData.channel;
                        trace = data.traces(row,channel);
                        [x,y] = getXY_(trace, israw);
                        x1 = mean(x(selpts1));
                        y1 = mean(y(selpts1));
                        x2 = mean(x(selpts2));
                        y2 = mean(y(selpts2));
                        m = (y2 - y1) / (x2 - x1);
                        if israw
                            data.traces(row,channel).y0 = m .* (trace.x - x1) + y1;
                        else
                            data.traces(row,channel).y0 = ...
                                data.traces(row,channel).y0 + m .* (x - x1) + y1;
                        end
                        didit = true;
                    end
                end
            end
        end
        if didit
            redraw_();
        else
            msgbox('Select two contiguous regions of data to baseline to with the brush tool.', ...
                'Baseline Sloping Two Region');
        end
    end

    function baselineSplineVisibleTraces_(ax, nsegments)
        if ~exist('ax', 'var') || isempty(ax)
            ax = vertcat(ui.groups.ax);
        end
        if ~exist('nsegments', 'var') || isempty(nsegments)
            if isfield(ui, 'baselineSplineNumSegmentsEdit') && isvalid(ui.baselineSplineNumSegmentsEdit)
                nsegments = max(1, str2num(ui.baselineSplineNumSegmentsEdit.String));
            else
                answer = inputdlg('# Spline Segments:', 'Baseline Spline', 1, {'3'});
                if isempty(answer); return; end
                nsegments = str2num(answer{1});
            end
        end
        didit = false;
        for i = 1:numel(ax)
            for j = 1:numel(ax(i).UserData.traces)
                brushdata = logical(ax(i).UserData.traces(j).BrushData);
                if any(brushdata)
                    selpts = find(brushdata);
                    row = ax(i).UserData.rows(j);
                    channel = ax(i).UserData.channel;
                    trace = data.traces(row,channel);
                    % !!! ALWAYS apply spline fit to the raw data
                    % This allows repeated spline trial fits without
                    % accruing distortions from bogus fits
                    [x,y] = getXY_(trace, true);
                    try
                        pp = splinefit(x(selpts), y(selpts), nsegments);
                        data.traces(row,channel).y0 = ppval(pp, x);
                    catch
                        msgbox("!!! Requires package 'splinefit'. See Add-On Explorer.", ...
                            'Baseline Spline');
                        return;
                    end
                    didit = true;
                end
            end
        end
        if didit
            redraw_();
        else
            msgbox('Select data to baseline to with the brush tool.', ...
                'Baseline Spline');
        end
    end

    function baselineSplineVisibleTracesInteractive_(ax)
        if ~exist('ax', 'var')
            ax = gobjects(0);
        end
        if isfield(ui, 'baselineSplineWindow') && isvalid(ui.baselineSplineWindow)
            figure(ui.baselineSplineWindow);
            return;
        end
        ui.baselineSplineWindow = figure('Name', 'Baseline Spline', ...
            'Menubar', 'none', ...
            'Toolbar', 'none', ...
            'numbertitle', 'off');
        ui.baselineSplineWindow.Position(3) = 250;
        ui.baselineSplineWindow.Position(4) = 75;
        uicontrol(ui.baselineSplineWindow, ...
            'Style', 'text', ...
            'String', '# Baseline Spline Segments', ...
            'Units', 'normalized', ...
            'Position', [0, 0.66, 1, 0.33]);
        ui.baselineSplineNumSegmentsEdit = uicontrol(ui.baselineSplineWindow, ...
            'Style', 'edit', ...
            'String', '3', ...
            'Units', 'normalized', ...
            'Position', [0, 0.33, 1, 0.33]);
        uicontrol(ui.baselineSplineWindow, ...
            'Style', 'pushbutton', ...
            'String', 'Fit Spline to Baseline', ...
            'Units', 'normalized', ...
            'Position', [0, 0, 1, 0.33], ...
            'Callback', @(s,e) baselineSplineVisibleTraces_(ax));
%         if ui.showRawBtn.Checked == "off"
%             ui.showRawBtn.Checked = 'on';
%             ui.showBaselineBtn.Checked = 'on';
%             refresh_();
%             autoscale_();
%         elseif ui.showBaselineBtn.Checked == "off"
%             ui.showBaselineBtn.Checked = 'on';
%             redraw_();
%         end
    end

    function scaleVisibleTraces_(ax, yscale)
        if ~exist('ax', 'var') || isempty(ax)
            ax = vertcat(ui.groups.ax);
        end
        if ~exist('yscale', 'var')
            answer = inputdlg('Y Scale:', 'Set Y Scale', 1, {'1'});
            if isempty(answer); return; end
            yscale = str2num(answer{1});
        end
        israw = showRaw_();
        for i = 1:numel(ax)
            rows = ax(i).UserData.rows;
            channel = ax(i).UserData.channel;
            if israw
                % set scale
                [data.traces(rows,channel).yscale] = deal(yscale);
            else
                % multiply current scale
                for j = 1:numel(rows)
                    data.traces(rows(j),channel).yscale = ...
                        data.traces(rows(j),channel).yscale .* yscale;
                end
            end
        end
        redraw_();
    end

    function normalizePositiveVisibleTraces_(ax)
        if ~exist('ax', 'var') || isempty(ax)
            ax = vertcat(ui.groups.ax);
        end
        israw = showRaw_();
        didit = false;
        for i = 1:numel(ax)
            for j = 1:numel(ax(i).UserData.traces)
                brushdata = logical(ax(i).UserData.traces(j).BrushData);
                if any(brushdata)
                    row = ax(i).UserData.rows(j);
                    channel = ax(i).UserData.channel;
                    trace = data.traces(row,channel);
                    [x,y] = getXY_(trace, israw);
                    if israw
                        data.traces(row,channel).yscale = 1.0 / max(y(brushdata));
                    else
                        data.traces(row,channel).yscale = trace.yscale ./ max(y(brushdata));
                    end
                    didit = true;
                end
            end
        end
        if didit
            redraw_();
        else
            msgbox('Select data to normalize to with the brush tool.', ...
                'Normalize Positive');
        end
    end

    function normalizeNegativeVisibleTraces_(ax)
        if ~exist('ax', 'var') || isempty(ax)
            ax = vertcat(ui.groups.ax);
        end
        israw = showRaw_();
        didit = false;
        for i = 1:numel(ax)
            for j = 1:numel(ax(i).UserData.traces)
                brushdata = logical(ax(i).UserData.traces(j).BrushData);
                if any(brushdata)
                    row = ax(i).UserData.rows(j);
                    channel = ax(i).UserData.channel;
                    trace = data.traces(row,channel);
                    [x,y] = getXY_(trace, israw);
                    if israw
                        data.traces(row,channel).yscale = 1.0 / abs(min(y(brushdata)));
                    else
                        data.traces(row,channel).yscale = trace.yscale ./ abs(min(y(brushdata)));
                    end
                    didit = true;
                end
            end
        end
        if didit
            redraw_();
        else
            msgbox('Select data to normalize to with the brush tool.', ...
                'Normalize Negative');
        end
    end

    function normalizeAbsPeakVisibleTraces_(ax)
        if ~exist('ax', 'var') || isempty(ax)
            ax = vertcat(ui.groups.ax);
        end
        israw = showRaw_();
        didit = false;
        for i = 1:numel(ax)
            for j = 1:numel(ax(i).UserData.traces)
                brushdata = logical(ax(i).UserData.traces(j).BrushData);
                if any(brushdata)
                    row = ax(i).UserData.rows(j);
                    channel = ax(i).UserData.channel;
                    trace = data.traces(row,channel);
                    [x,y] = getXY_(trace, israw);
                    if israw
                        data.traces(row,channel).yscale = 1.0 / max(abs(y(brushdata)));
                    else
                        data.traces(row,channel).yscale = trace.yscale ./ max(abs(y(brushdata)));
                    end
                    didit = true;
                end
            end
        end
        if didit
            redraw_();
        else
            msgbox('Select data to normalize to with the brush tool.', ...
                'Normalize Abs Peak');
        end
    end

    function setTimeZeroVisibleTraces_(ax, x0)
        if ~exist('ax', 'var') || isempty(ax)
            ax = vertcat(ui.groups.ax);
        end
        if ~exist('x0', 'var')
            answer = inputdlg('Time Zero:', 'Set Time Zero', 1, {'0'});
            if isempty(answer); return; end
            x0 = str2num(answer{1});
        end
        israw = showRaw_();
        for i = 1:numel(ax)
            rows = ax(i).UserData.rows;
            channel = ax(i).UserData.channel;
            if israw
                % set time zero
                [data.traces(rows,channel).x0] = deal(x0);
            else
                % shift time zero
                for j = 1:numel(rows)
                    data.traces(rows(j),channel).x0 = ...
                        data.traces(rows(j),channel).x0 + x0;
                end
            end
        end
        redraw_();
    end

    function alignToOnsetVisibleTraces_(ax, xSD, direction)
        if ~exist('ax', 'var') || isempty(ax)
            ax = vertcat(ui.groups.ax);
        end
        if ~exist('xSD', 'var') || isempty(xSD)
            answer = inputdlg('x Standard Deviation:', 'Align to Onset', 1, {'8'});
            if isempty(answer); return; end
            xSD = str2num(answer{1});
        end
        if ~exist('direction', 'var') || isempty(direction)
            choices = {'Absolute Value', 'Positive', 'Negative'};
            ind = listdlg('ListString', choices, 'SelectionMode', 'single');
            if isempty(ind); return; end
            direction = choices{ind};
        end
        israw = showRaw_();
        didit = false;
        for i = 1:numel(ax)
            for j = 1:numel(ax(i).UserData.traces)
                brushdata = logical(ax(i).UserData.traces(j).BrushData);
                if any(brushdata)
                    selpts = find(brushdata);
                    row = ax(i).UserData.rows(j);
                    channel = ax(i).UserData.channel;
                    trace = data.traces(row,channel);
                    [x,y] = getXY_(trace, israw);
                    y = y - mean(y(brushdata));
                    threshold = xSD * std(y(brushdata));
                    if direction == "Absolute Value"
                        ind = find(abs(y(selpts(end)+1:end)) > threshold, 1);
                    elseif direction == "Positive"
                        ind = find(y(selpts(end)+1:end) > threshold, 1);
                    elseif direction == "Negative"
                        ind = find(y(selpts(end)+1:end) < -threshold, 1);
                    end
                    if israw
                        data.traces(row,channel).x0 = x(selpts(end) + ind);
                    else
                        data.traces(row,channel).x0 = ...
                            data.traces(row,channel).x0 + x(selpts(end) + ind);
                    end
                    didit = true;
                end
            end
        end
        if didit
            redraw_();
        else
            msgbox('Select baseline data just prior to peak onset with the brush tool.', ...
                'Align To Onset');
        end
    end

    function maskSelectedRegions_(ax)
        if ~exist('ax', 'var') || isempty(ax)
            ax = vertcat(ui.groups.ax);
        end
        didit = false;
        for i = 1:numel(ax)
            for j = 1:numel(ax(i).UserData.traces)
                brushdata = logical(ax(i).UserData.traces(j).BrushData);
                if any(brushdata)
                    row = ax(i).UserData.rows(j);
                    channel = ax(i).UserData.channel;
                    if isempty(data.traces(row,channel).masked)
                        data.traces(row,channel).masked = brushdata';
                    else
                        data.traces(row,channel).masked = ...
                            data.traces(row,channel).masked | brushdata';
                    end
                    didit = true;
                end
            end
        end
        if didit
            redraw_();
        else
            msgbox('Select data to mask with the brush tool.', ...
                'Mask Selected Region');
        end
    end

    function zeroSelectedRegions_(ax)
        if ~exist('ax', 'var') || isempty(ax)
            ax = vertcat(ui.groups.ax);
        end
        didit = false;
        for i = 1:numel(ax)
            for j = 1:numel(ax(i).UserData.traces)
                brushdata = logical(ax(i).UserData.traces(j).BrushData);
                if any(brushdata)
                    row = ax(i).UserData.rows(j);
                    channel = ax(i).UserData.channel;
                    if isempty(data.traces(row,channel).zeroed)
                        data.traces(row,channel).zeroed = brushdata';
                    else
                        data.traces(row,channel).zeroed = ...
                            data.traces(row,channel).zeroed | brushdata';
                    end
                    didit = true;
                end
            end
        end
        if didit
            redraw_();
        else
            msgbox('Select data to zero with the brush tool.', ...
                'Zero Selected Region');
        end
    end

    function interpolateSelectedRegions_(ax)
        if ~exist('ax', 'var') || isempty(ax)
            ax = vertcat(ui.groups.ax);
        end
        didit = false;
        for i = 1:numel(ax)
            for j = 1:numel(ax(i).UserData.traces)
                brushdata = logical(ax(i).UserData.traces(j).BrushData);
                if any(brushdata)
                    row = ax(i).UserData.rows(j);
                    channel = ax(i).UserData.channel;
                    if isempty(data.traces(row,channel).interpolated)
                        data.traces(row,channel).interpolated = brushdata';
                    else
                        data.traces(row,channel).interpolated = ...
                            data.traces(row,channel).interpolated | brushdata';
                    end
                    didit = true;
                end
            end
        end
        if didit
            redraw_();
        else
            msgbox('Select data to interpolate with the brush tool.', ...
                'Interpolate Selected Region');
        end
    end

    function measure_(ax)
        if ~exist('ax', 'var') || isempty(ax)
            ax = vertcat(ui.groups.ax);
        end
        israw = showRaw_();
        values = [];
        for i = 1:numel(ax)
            for j = 1:numel(ax(i).UserData.traces)
                brushdata = logical(ax(i).UserData.traces(j).BrushData);
                if any(brushdata)
                    row = ax(i).UserData.rows(j);
                    channel = ax(i).UserData.channel;
                    trace = data.traces(row,channel);
                    [x,y] = getXY_(trace, israw);
                    values = [values; y(brushdata)];
                end
            end
        end
        if ~isempty(values)
            msgbox({['Mean: ' num2str(mean(values))]; ...
                ['Std: ' num2str(std(values))]; ...
                ['Var: ' num2str(var(values))]}, ...
                'Measure');
        else
            msgbox('Select data to measure with the brush tool.', ...
                'Measure');
        end
    end

    function setIsMaskedVisibleTraces_(ax, ismasked)
        if ~exist('ax', 'var') || isempty(ax)
            ax = vertcat(ui.groups.ax);
        end
        didit = false;
        for i = 1:numel(ax)
            rows = ax(i).UserData.rows;
            channel = ax(i).UserData.channel;
            if ~isempty(rows)
                [data.traces(rows,channel).ismasked] = deal(ismasked);
                didit = true;
            end
        end
        if didit
            refresh_();
        else
            msgbox('No visible traces. To unmask, first show masked traces.', ...
                'Mask/UnMask');
        end
    end

    function revertToRawVisibleTraces_(ax)
        if questdlg('Revert visible traces to raw data?', ...
                'Revert to Raw', ...
                'OK', 'Cancel', 'Cancel') == "Cancel"
            return;
        end
        if ~exist('ax', 'var') || isempty(ax)
            ax = vertcat(ui.groups.ax);
        end
        for i = 1:numel(ax)
            rows = ax(i).UserData.rows;
            channel = ax(i).UserData.channel;
            [data.traces(rows,channel).x0] = deal(0);
            [data.traces(rows,channel).y0] = deal(0);
            [data.traces(rows,channel).yscale] = deal(1);
            [data.traces(rows,channel).masked] = deal([]);
            [data.traces(rows,channel).zeroed] = deal([]);
            [data.traces(rows,channel).interpolated] = deal([]);
        end
        refresh_();
    end

    function deleteVisibleTraces_(ax)
        if questdlg('Delete visible sweeps? Traces for all channels will be deleted.', ...
                'Delete', ...
                'OK', 'Cancel', 'Cancel') == "Cancel"
            return;
        end
        if ~exist('ax', 'var') || isempty(ax)
            ax = vertcat(ui.groups.ax);
        end
        rows = unique(vertcat(ax.UserData.rows));
        data.traces(rows,:) = [];
        data.groupids(rows) = [];
        updateUI_();
    end

%% UI

    function initUI_()
        ui.mainWindow = figure( ...
            'Name', 'Patch Meister', ...
            'Units', 'normalized', ...
            'Position', [0 0 1 1], ...
            'numbertitle', 'off');
        ui.mainWindow.Units = 'pixels';

        loadIcons_();
        
        ui.menu = createMainMenu_();
        
        % visible channels
        ui.channelspanel = uipanel(ui.mainWindow, ...
            'Units', 'pixels', ...
            'Position', [0, 0, 100, 60]);
        uicontrol(ui.channelspanel, ...
            'Style', 'text', ...
            'String', 'Channels', ...
            'HorizontalAlignment', 'left', ...
            'ForegroundColor', [0 0 1], ...
            'Units', 'pixels', ...
            'Position', [0, 45, 100, 15]);
        ui.visibleChannels = uicontrol(ui.channelspanel, ...
            'Style', 'listbox', ...
            'Value', [], ...
            'Units', 'pixels', ...
            'Position', [0, 0, 100, 45], ...
            'Callback', @(s,e) updateUI_());
        
        % sweep traversal
        ui.sweepspanel = uipanel(ui.mainWindow, ...
            'Units', 'pixels', ...
            'Position', [0, 0, 100, 70]);
        uicontrol(ui.sweepspanel, ...
            'Style', 'text', ...
            'String', 'Sweeps', ...
            'HorizontalAlignment', 'left', ...
            'ForegroundColor', [0 0 1], ...
            'Units', 'pixels', ...
            'Position', [0, 55, 100, 15]);
        ui.prevSweepBtn = uicontrol(ui.sweepspanel, ...
            'Style', 'pushbutton', ...
            'Tooltip', 'Previous Sweep', ...
            'String', '<', ...
            'Units', 'pixels', ...
            'Position', [0, 20, 50, 35], ...
            'Callback', @(s,e) selectPrevSweep_());
        ui.nextSweepBtn = uicontrol(ui.sweepspanel, ...
            'Style', 'pushbutton', ...
            'Tooltip', 'Next Sweep', ...
            'String', '>', ...
            'Units', 'pixels', ...
            'Position', [50, 20, 50, 35], ...
            'Callback', @(s,e) selectNextSweep_());
        ui.selectedSweepsEdit = uicontrol(ui.sweepspanel, ...
            'Style', 'edit', ...
            'Tooltip', 'Selected Sweeps, e.g. 1 3,5-12', ...
            'String', '1', ...
            'Units', 'pixels', ...
            'Position', [0, 0, 60, 20], ...
            'Callback', @(s,e) refresh_());
        ui.numSweepsBtn = uicontrol(ui.sweepspanel, ...
            'Style', 'pushbutton', ...
            'Tooltip', 'Select All Sweeps', ...
            'String', '/ N', ...
            'Units', 'pixels', ...
            'Position', [60, 0, 40, 20], ...
            'Callback', @(s,e) selectAllSweeps_());
        
        % buttons
        ui.btnspanel = uipanel(ui.mainWindow, ...
            'Units', 'pixels', ...
            'Position', [0, 0, 100, 40]);
        uicontrol(ui.btnspanel, ...
            'Style', 'text', ...
            'String', 'Autoscale', ...
            'HorizontalAlignment', 'left', ...
            'ForegroundColor', [0 0 1], ...
            'Units', 'pixels', ...
            'Position', [0, 25, 100, 15]);
        uicontrol(ui.btnspanel, ...
            'Style', 'pushbutton', ...
            'Tooltip', 'Autoscale X', ...
            'Units', 'pixels', ...
            'Position', [0, 0, 25, 25], ...
            'CData', imresize(ui.icons('autoscaleX'), [23,23]), ...
            'Callback', @(s,e) autoscale_([], "x"));
        uicontrol(ui.btnspanel, ...
            'Style', 'pushbutton', ...
            'Tooltip', 'Autoscale Y', ...
            'Units', 'pixels', ...
            'Position', [25, 0, 25, 25], ...
            'CData', imresize(ui.icons('autoscaleY'), [23,23]), ...
            'Callback', @(s,e) autoscale_([], "y"));
        uicontrol(ui.btnspanel, ...
            'Style', 'pushbutton', ...
            'Tooltip', 'Autoscale Y Same', ...
            'Units', 'pixels', ...
            'Position', [50, 0, 25, 25], ...
            'CData', imresize(ui.icons('autoscaleYSame'), [23,23]), ...
            'Callback', @(s,e) autoscale_([], "y", true));
        uicontrol(ui.btnspanel, ...
            'Style', 'pushbutton', ...
            'Tooltip', 'Autoscale XY', ...
            'Units', 'pixels', ...
            'Position', [75, 0, 25, 25], ...
            'CData', imresize(ui.icons('autoscaleXY'), [23,23]), ...
            'Callback', @(s,e) autoscale_([], "xy"));
        
        % display options
        h = 18;
        ui.optionspanel = uipanel(ui.mainWindow, ...
            'Units', 'pixels', ...
            'Position', [0, 0, 100, 3*h+45]);
        uicontrol(ui.optionspanel, ...
            'Style', 'text', ...
            'String', 'Display Options', ...
            'HorizontalAlignment', 'left', ...
            'ForegroundColor', [0 0 1], ...
            'Units', 'pixels', ...
            'Position', [0, 30+3*h, 100, 15]);
        ui.showRawBtn = uicontrol(ui.optionspanel, ...
            'Style', 'togglebutton', ...
            'String', 'Raw', ...
            'Tooltip', 'Show Raw', ...
            'Units', 'pixels', ...
            'Position', [0, 30+2*h, 50, h], ...
            'Callback', @(s,e) refresh_(true));
        ui.showBaselineBtn = uicontrol(ui.optionspanel, ...
            'Style', 'togglebutton', ...
            'String', 'Baseline', ...
            'Tooltip', 'Show Baseline', ...
            'Units', 'pixels', ...
            'Position', [50, 30+2*h, 50, h], ...
            'Callback', @(s,e) refresh_());
        ui.showMaskedBtn = uicontrol(ui.optionspanel, ...
            'Style', 'togglebutton', ...
            'String', 'Masked', ...
            'Tooltip', 'Show Masked', ...
            'Units', 'pixels', ...
            'Position', [0, 30+h, 50, h], ...
            'Callback', @(s,e) refresh_());
        ui.showAverageBtn = uicontrol(ui.optionspanel, ...
            'Style', 'togglebutton', ...
            'String', 'Average', ...
            'Tooltip', 'Show Average', ...
            'Units', 'pixels', ...
            'Position', [50, 30+h, 50, h], ...
            'Callback', @(s,e) refresh_());
        ui.showAverageOnlyBtn = uicontrol(ui.optionspanel, ...
            'Style', 'togglebutton', ...
            'String', 'Average ONLY', ...
            'Tooltip', 'Show Average ONLY', ...
            'Units', 'pixels', ...
            'Position', [0, 30, 100, h], ...
            'Callback', @(s,e) refresh_());
        ui.denoiseBtn = uicontrol(ui.optionspanel, ...
            'Style', 'togglebutton', ...
            'String', 'Auto Denoise', ...
            'Tooltip', 'Denoise', ...
            'Units', 'pixels', ...
            'Position', [0, 15, 100, h], ...
            'Callback', @(s,e) redraw_());
        uicontrol(ui.optionspanel, ...
            'Style', 'text', ...
            'String', 'decimate:', ...
            'HorizontalAlignment', 'right', ...
            'Units', 'pixels', ...
            'Position', [0, 0, 50, 15]);
        ui.decimateEdit = uicontrol(ui.optionspanel, ...
            'Style', 'edit', ...
            'String', '', ...
            'Units', 'pixels', ...
            'Position', [50, 0, 50, 15], ...
            'Callback', @(s,e) refresh_());
%         uicontrol(ui.optionspanel, ...
%             'Style', 'text', ...
%             'String', 'filter type:', ...
%             'HorizontalAlignment', 'right', ...
%             'Units', 'pixels', ...
%             'Position', [0, 15, 50, 15]);
%         ui.lowpassFilterType = uicontrol(ui.optionspanel, ...
%             'Style', 'popupmenu', ...
%             'String', {'Bessel', 'Butterworth'}, ...
%             'Units', 'pixels', ...
%             'Position', [50, 15, 50, 15], ...
%             'Callback', @(s,e) refresh_());
%         uicontrol(ui.optionspanel, ...
%             'Style', 'text', ...
%             'String', 'filter (Hz):', ...
%             'HorizontalAlignment', 'right', ...
%             'Units', 'pixels', ...
%             'Position', [0, 0, 50, 15]);
%         ui.lowpassFilterCutoffFreqHzEdit = uicontrol(ui.optionspanel, ...
%             'Style', 'edit', ...
%             'String', '', ...
%             'Units', 'pixels', ...
%             'Position', [50, 0, 50, 15], ...
%             'Callback', @(s,e) refresh_());
        
        ui.groups = repmat(struct(), [0,0]);
        
        set(ui.mainWindow, 'KeyPressFcn', @(s,e) keyPress_(e));
        set(ui.mainWindow, 'SizeChangedFcn', @(s,e) resizeUI_());
        ui.path = pwd();

        resizeUI_();
    end

    function updateUI_()
        ngroups = getNumGroups_();
        nchannels = getNumChannels_();
        % update visible channels from ui list
        ui.visibleChannels.String = cellstr(horzcat(data.traces(1,:).ylabel));
        ui.visibleChannels.Min = 0;
        ui.visibleChannels.Max = nchannels;
        if isempty(ui.visibleChannels.Value)
            ui.visibleChannels.Value = 1:nchannels;
        end
        ui.visibleChannels.Value(ui.visibleChannels.Value > nchannels) = [];
        nvischannels = numel(ui.visibleChannels.Value);
        % delete all current group UIs
        if isfield(ui.groups, 'panel')
            delete(vertcat(ui.groups.panel));
        end
        ui.groups = repmat(struct(), [0,0]);
        % create new group UIs
        for i = 1:ngroups % group index
            ui.groups(i).panel = uipanel(ui.mainWindow, ...
                'Units', 'pixels');
            rows = getTracesInGroup_(i);
            for j = 1:nvischannels
                channel = ui.visibleChannels.Value(j); % channel index
                ui.groups(i).ax(j,1) = axes(ui.groups(i).panel, ...
                    'Units', 'pixels', ...
                    'TickLength', [0.005, 0.01]);
                ax = ui.groups(i).ax(j);
                ax.UserData.groupid = i;
                ax.UserData.channel = channel;
                if j ~= nvischannels
                    ax.XTickLabels = [];
                end
                if i ~= ngroups || j ~= nvischannels
                    xlabel(ax, '');
                else
                    xlabel(ax, [char(data.traces(rows(1),1).xlabel) ' (' char(data.traces(rows(1),1).xunit) ')']);
                end
                ylabel(ax, {getGroupLabel_(i); ...
                    [char(data.traces(rows(1),channel).ylabel) ' (' char(data.traces(rows(1),channel).yunit) ')']});
                % axes context menu
                ax.UserData.menu = createAxesMenu_(ax);
                ax.UIContextMenu = ax.UserData.menu;
                ax.UserData.menuBtn = uicontrol(ui.groups(i).panel, ...
                    'Style', 'pushbutton', ...
                    'Units', 'pixels', ...
                    'Position', [2, 0, 20, 20], ...
                    'CData', imresize(ui.icons('menu'), [19, 19]), ...
                    'UserData', ax, ...
                    'Callback', @popupAxesMenu_);
            end
        end
        linkaxes(vertcat(ui.groups.ax), 'x');
        ui.numSweepsBtn.String = ['/ ' num2str(getMaxNumSweepsPerGroup_())];
        resizeUI_();
        refresh_();
        autoscale_();
        updateMainWindowTitle_();
    end

    function resizeUI_()
        w = ui.mainWindow.Position(3);
        h = ui.mainWindow.Position(4);
        margin = 2;
        sep = 2;
        
        x = margin;
        y = h-margin;
        ui.channelspanel.Position(1) = x;
        ui.channelspanel.Position(2) = y - ui.channelspanel.Position(4);
        y = y-ui.channelspanel.Position(4)-sep;
        
        ui.sweepspanel.Position(1) = x;
        ui.sweepspanel.Position(2) = y - ui.sweepspanel.Position(4);
        y = y-ui.sweepspanel.Position(4)-sep;
        
        ui.btnspanel.Position(1) = x;
        ui.btnspanel.Position(2) = y - ui.btnspanel.Position(4);
        y = y-ui.btnspanel.Position(4)-sep;
        
        ui.optionspanel.Position(1) = x;
        ui.optionspanel.Position(2) = y - ui.optionspanel.Position(4);
        y = y-ui.optionspanel.Position(4)-sep;
        
%         btnsz = [100,100]./3.3;
%         ui.btnControls.autoscaleXBtn.CData = imresize(ui.icons('autoscaleX'), btnsz);
%         ui.btnControls.autoscaleYBtn.CData = imresize(ui.icons('autoscaleY'), btnsz);
%         ui.btnControls.autoscaleXYBtn.CData = imresize(ui.icons('autoscaleXY'), btnsz);
%         ui.btnControls.baselineFlatBtn.CData = imresize(ui.icons('baselineFlat'), btnsz);
%         ui.btnControls.baselineSlopingBtn.CData = imresize(ui.icons('baselineSloping'), btnsz);
%         ui.btnControls.baselineSplineBtn.CData = imresize(ui.icons('baselineSpline'), btnsz);
%         ui.btnControls.normalizePositiveBtn.CData = imresize(ui.icons('normalizePositive'), btnsz);
%         ui.btnControls.normalizeNegativeBtn.CData = imresize(ui.icons('normalizeNegative'), btnsz);
%         ui.btnControls.normalizeAbsPeakBtn.CData = imresize(ui.icons('normalizeAbsPeak'), btnsz);
%         ui.btnControls.timeZeroBtn.CData = imresize(ui.icons('timeZero'), btnsz);
%         ui.btnControls.alignToOnsetBtn.CData = imresize(ui.icons('alignToOnset'), btnsz);
%         ui.btnControls.zeroRegionBtn.CData = imresize(ui.icons('zeroRegion'), btnsz);
%         ui.btnControls.interpolateRegionBtn.CData = imresize(ui.icons('interpolateRegion'), btnsz);
%         ui.btnControls.maskRegionBtn.CData = imresize(ui.icons('maskRegion'), btnsz);
        
        x = x+100+margin;
        ngroups = numel(ui.groups);
        gh = floor(double(h-2*margin-15) / ngroups);
        for i = 1:ngroups
            y = h-margin-i*gh;
            if i < ngroups
                ui.groups(i).panel.Position = [x, y, w-margin-x, gh];
            else
                ui.groups(i).panel.Position = [x, y-15, w-margin-x, gh+15];
            end
            pw = ui.groups(i).panel.Position(3);
            ph = ui.groups(i).panel.Position(4);
            nvischannels = numel(ui.groups(i).ax);
            for j = 1:nvischannels
                ax = ui.groups(i).ax(j);
                if i < ngroups
                    axh = floor((ph-27) / nvischannels);
                    ax.Position = [80, ph-7-j*axh, pw-100, axh-10];
                else
                    axh = floor((ph-27-15) / nvischannels);
                    ax.Position = [80, ph-7-j*axh, pw-100, axh-10];
                end
                ax.UserData.menuBtn.Position(2) = ...
                    ax.Position(2) + ax.Position(4) - ax.UserData.menuBtn.Position(4);
            end
        end
    end

    function redraw_(refresh)
        if ~exist('refresh', 'var'); refresh = false; end
        vischannels = getVisibleChannels_();
        selsweeps = getSelectedSweeps_();
        groupuids = unique(data.groupids);
        for i = 1:numel(ui.groups)
            grouprows = find(data.groupids == groupuids(i));
            for j = 1:numel(ui.groups(i).ax)
                ax = ui.groups(i).ax(j);
                channel = vischannels(j);
                if refresh
                    cla(ax);
                    ax.UserData.groupid = groupuids(i);
                    ax.UserData.channel = channel;
                    ax.UserData.rows = [];
                    ax.UserData.traces = gobjects(0);
                    ax.UserData.baselines = gobjects(0);
                    ax.UserData.avgtrace = gobjects(0);
                    hold(ax, 'on');
                    cmap = colormap(lines());
                    ncolors = size(cmap,1);
                    groupselsweeps = intersect(1:numel(grouprows), selsweeps);
                    if ~isempty(groupselsweeps)
                        rows = grouprows(groupselsweeps);
                        nrows = numel(rows);
                        if showRaw_() % show raw data
                            if ~showAverageOnly_() % show traces
                                for k = 1:nrows
                                    trace = data.traces(rows(k),channel);
                                    [x,y] = getXY_(trace, showRaw_());
                                    if ~trace.ismasked % not masked
                                        ax.UserData.traces = [ax.UserData.traces; ...
                                            plot(ax, x, y, 'color', cmap(1+mod(k-1,ncolors),:))];
                                        ax.UserData.rows = [ax.UserData.rows; rows(k)];
                                    elseif showMasked_() % show masked traces
                                        ax.UserData.traces = [ax.UserData.traces; ...
                                            plot(ax, x, y, 'color', [0.5, 0.5, 0.5])];
                                        ax.UserData.rows = [ax.UserData.rows; rows(k)];
                                    end
                                    if showBaseline_() ... % show baseline
                                            && (~trace.ismasked || showMasked_())
                                        y0 = trace.y0;
                                        if numel(y0) == 1
                                            y0 = repmat(y0, size(y));
                                        end
                                        ax.UserData.baselines = [ax.UserData.baselines; ...
                                            plot(ax, x, y0, 'k--')];
                                    end
                                end
                            end
                            if showAverage_() || showAverageOnly_() % show average
                                [x,y] = getMeanXY_(data.traces(rows,channel), showRaw_());
                                ax.UserData.avgtrace = plot(ax, x, y, 'k-');
                            end
                        else % show offset and scaled data
                            if ~showAverageOnly_() % show traces
                                for k = 1:nrows
                                    trace = data.traces(rows(k),channel);
                                    [x,y] = getXY_(trace, showRaw_());
                                    if ~trace.ismasked
                                        ax.UserData.traces = [ax.UserData.traces; ...
                                            plot(ax, x, y, 'color', cmap(1+mod(k-1,ncolors),:))];
                                        ax.UserData.rows = [ax.UserData.rows; rows(k)];
                                    elseif showMasked_() % show masked traces
                                        ax.UserData.traces = [ax.UserData.traces; ...
                                            plot(ax, x, y, 'color', [0.5, 0.5, 0.5])];
                                        ax.UserData.rows = [ax.UserData.rows; rows(k)];
                                    end
                                end
                            end
                            if showAverage_() || showAverageOnly_() % show average
                                [x,y] = getMeanXY_(data.traces(rows,channel), showRaw_());
                                ax.UserData.avgtrace = plot(ax, x, y, 'k-');
                            end
                            if showBaseline_() % show baseline
                                plot(ax, ax.XLim', zeros([2,1]), 'k--');
                            end
                        end
                    end
                else % refresh
                    for k = 1:numel(ax.UserData.traces)
                        trace = data.traces(ax.UserData.rows(k),channel);
                        [x,y] = getXY_(trace, showRaw_());
                        ax.UserData.traces(k).XData(:) = x;
                        ax.UserData.traces(k).YData(:) = y;
                        if numel(ax.UserData.baselines) >= k && showRaw_()
                            y0 = trace.y0;
                            if numel(y0) == 1
                                y0 = repmat(y0, size(y));
                            end
                            ax.UserData.baselines(k).XData(:) = x;
                            ax.UserData.baselines(k).YData(:) = y0;
                        end
                    end
                    if ~isempty(ax.UserData.avgtrace)
                        [x,y] = getMeanXY_(data.traces(ax.UserData.rows,channel), showRaw_());
                        ax.UserData.avgtrace.XData(:) = x;
                        ax.UserData.avgtrace.YData(:) = y;
                    end
                end
            end
        end
    end

    function refresh_(autoscale)
        if ~exist('autoscale', 'var'); autoscale = false; end
        redraw_(true);
        if autoscale
            autoscale_();
        end
    end

    function loadIcons_()
        ui.icons = containers.Map();
        ui.icons('menu') = imread('icons/menu.png');
        ui.icons('autoscaleX') = imread('icons/autoscaleX.png');
        ui.icons('autoscaleY') = imread('icons/autoscaleY.png');
        ui.icons('autoscaleYSame') = imread('icons/autoscaleYSame.png');
        ui.icons('autoscaleXY') = imread('icons/autoscaleXY.png');
%         ui.icons('baselineFlat') = imread('icons/baselineFlat.png');
%         ui.icons('baselineSloping') = imread('icons/baselineSloping.png');
%         ui.icons('baselineSpline') = imread('icons/baselineSpline.png');
%         ui.icons('normalizePositive') = imread('icons/normalizePositive.png');
%         ui.icons('normalizeNegative') = imread('icons/normalizeNegative.png');
%         ui.icons('normalizeAbsPeak') = imread('icons/normalizeAbsPeak.png');
%         ui.icons('timeZero') = imread('icons/timeZero.png');
%         ui.icons('alignToOnset') = imread('icons/alignToOnset.png');
%         ui.icons('zeroRegion') = imread('icons/zeroRegion.png');
%         ui.icons('interpolateRegion') = imread('icons/interpolateRegion.png');
%         ui.icons('maskRegion') = imread('icons/maskRegion.png');
%         ui.icons('chainlinked') = imread('icons/chainlinked.png');
%         ui.icons('chainunlinked') = imread('icons/chainunlinked.png');
    end

    function menu = createMainMenu_()
        menu = uimenu(ui.mainWindow, ...
            'Text', 'PatchMeister');
        
        uimenu(menu, ...
            'Text', 'About', ...
            'MenuSelectedFcn', @(s,e) about_());

        uimenu(menu, ...
            'Separator', 'on', ...
            'Text', 'Load Data', ...
            'MenuSelectedFcn', @(s,e) loadData_());
        uimenu(menu, ...
            'Text', 'Save Data', ...
            'MenuSelectedFcn', @(s,e) saveData_());

        uimenu(menu, ...
            'Separator', 'on', ...
            'Text', 'Import HEKA Data File', ...
            'MenuSelectedFcn', @(s,e) loadHEKA_());
        uimenu(menu, ...
            'Text', 'Import Axon ABF Data File', ...
            'MenuSelectedFcn', @(s,e) loadABF_());
%         uimenu(menu, ...
%             'Text', 'Import Axograph Data File', ...
%             'MenuSelectedFcn', @(s,e) loadAxograph_());

        uimenu(menu, ...
            'Separator', 'on', ...
            'Text', 'Edit Patch Info', ...
            'MenuSelectedFcn', @(s,e) editInfo_());
        uimenu(menu, ...
            'Text', 'Show/Edit Notes', ...
            'MenuSelectedFcn', @(s,e) editNotes_());
        uimenu(menu, ...
            'Text', 'Edit Channels', ...
            'MenuSelectedFcn', @(s,e) editChannels_());
        
        uimenu(menu, ...
            'Separator', 'on', ...
            'Text', 'Merge Groups', ...
            'MenuSelectedFcn', @(s,e) mergeGroups_());
        uimenu(menu, ...
            'Text', 'Group Interleaved', ...
            'MenuSelectedFcn', @(s,e) groupInterleaved_());
        uimenu(menu, ...
            'Text', 'Group Blocks', ...
            'MenuSelectedFcn', @(s,e) groupBlocks_());
        uimenu(menu, ...
            'Text', 'Edit Groups', ...
            'MenuSelectedFcn', @(s,e) editGroups_());
        uimenu(menu, ...
            'Text', 'Edit Group Labels', ...
            'MenuSelectedFcn', @(s,e) editGroupLabels_());

        % buttons are now available for these
%         uimenu(menu, ...
%             'Separator', 'on', ...
%             'Label', 'Autoscale All XY', ...
%             'Callback', @(s,e) autoscale_([], "xy"));
%         uimenu(menu, ...
%             'Label', 'Autoscale All X', ...
%             'Callback', @(s,e) autoscale_([], "x"));
%         uimenu(menu, ...
%             'Label', 'Autoscale All Y', ...
%             'Callback', @(s,e) autoscale_([], "y"));
%         uimenu(menu, ...
%             'Label', 'Autoscale All Same Y', ...
%             'Callback', @(s,e) autoscale_([], "y", true));

        % buttons are now available for these
%         ui.showRawBtn = uimenu(menu, ...
%             'Separator', 'on', ...
%             'Label', 'Show Raw', ...
%             'Callback', @toggleDisplayOptionBtn_);
%         ui.showBaselineBtn = uimenu(menu, ...
%             'Label', 'Show Baseline', ...
%             'Callback', @toggleDisplayOptionBtn_);
%         ui.showMaskedBtn = uimenu(menu, ...
%             'Label', 'Show Masked', ...
%             'Callback', @toggleDisplayOptionBtn_);
%         ui.showAverageBtn = uimenu(menu, ...
%             'Label', 'Show Average', ...
%             'Callback', @toggleDisplayOptionBtn_);
%         ui.showAverageOnlyBtn = uimenu(menu, ...
%             'Label', 'Show Average ONLY', ...
%             'Callback', @toggleDisplayOptionBtn_);
    end

    function menu = createAxesMenu_(ax)
        menu = uicontextmenu(ui.mainWindow);
        uimenu(menu, ...
            'Label', 'Autoscale XY', ...
            'Callback', @(s,e) autoscale_(ax));
        uimenu(menu, ...
            'Label', 'Autoscale X', ...
            'Callback', @(s,e) autoscale_(ax, 'x'));
        uimenu(menu, ...
            'Label', 'Autoscale Y', ...
            'Callback', @(s,e) autoscale_(ax, 'y'));
        
        uimenu(menu, ...
            'Separator', 'on', ...
            'Label', 'Set Baseline', ...
            'Callback', @(s,e) setBaselineVisibleTraces_(ax));
        uimenu(menu, ...
            'Label', 'Baseline Flat', ...
            'Callback', @(s,e) baselineFlatVisibleTraces_(ax));
        uimenu(menu, ...
            'Label', 'Baseline Sloping Two Region', ...
            'Callback', @(s,e) baselineSlopingVisibleTraces_(ax));
        uimenu(menu, ...
            'Label', 'Baseline Spline', ...
            'Callback', @(s,e) baselineSplineVisibleTracesInteractive_(ax));
        
        uimenu(menu, ...
            'Separator', 'on', ...
            'Label', 'Set Y Scale', ...
            'Callback', @(s,e) scaleVisibleTraces_(ax));
        uimenu(menu, ...
            'Label', 'Normalize Positive', ...
            'Callback', @(s,e) normalizePositiveVisibleTraces_(ax));
        uimenu(menu, ...
            'Label', 'Normalize Negative', ...
            'Callback', @(s,e) normalizeNegativeVisibleTraces_(ax));
        uimenu(menu, ...
            'Label', 'Normalize Abs Peak', ...
            'Callback', @(s,e) normalizeAbsPeakVisibleTraces_(ax));
        
        uimenu(menu, ...
            'Separator', 'on', ...
            'Label', 'Set Time Zero', ...
            'Callback', @(s,e) setTimeZeroVisibleTraces_(ax));
        uimenu(menu, ...
            'Label', 'Align to Onset', ...
            'Callback', @(s,e) alignToOnsetVisibleTraces_(ax));
        
        uimenu(menu, ...
            'Separator', 'on', ...
            'Label', 'Zero Selected Regions', ...
            'Callback', @(s,e) zeroSelectedRegions_(ax));
        uimenu(menu, ...
            'Label', 'Interpolate Selected Regions', ...
            'Callback', @(s,e) interpolateSelectedRegions_(ax));
        uimenu(menu, ...
            'Label', 'Mask Selected Regions', ...
            'Callback', @(s,e) maskSelectedRegions_(ax));
        
        uimenu(menu, ...
            'Separator', 'on', ...
            'Label', 'Measure', ...
            'Callback', @(s,e) measure_(ax));
        
        uimenu(menu, ...
            'Separator', 'on', ...
            'Label', 'Mask', ...
            'Callback', @(s,e) setIsMaskedVisibleTraces_(ax, true));
        uimenu(menu, ...
            'Label', 'UnMask', ...
            'Callback', @(s,e) setIsMaskedVisibleTraces_(ax, false));
        
        uimenu(menu, ...
            'Separator', 'on', ...
            'Label', 'Revert to Raw', ...
            'Callback', @(s,e) revertToRawVisibleTraces_(ax));
        
        uimenu(menu, ...
            'Separator', 'on', ...
            'Label', 'Delete', ...
            'Callback', @(s,e) deleteVisibleTraces_(ax));
    end

    function popupAxesMenu_(btn, varargin)
        ax = btn.UserData;
        btnPos = getpixelposition(btn, true);
        ax.UserData.menu.Position(1) = btnPos(1);
        ax.UserData.menu.Position(2) = btnPos(2);
        ax.UserData.menu.Visible = 'on';
    end

function UNUSED_toggleMenuDisplayOptionBtn_(btn, varargin)
        if btn.Checked == "on"
            btn.Checked = 'off';
        else
            btn.Checked = 'on';
        end
        refresh_();
        if btn.Text == "Show Raw"
            autoscale_();
        end
    end

function UNUSED_createViewControls_()
    ui.viewControls.panel = uipanel(ui.mainWindow, ...
        'Units', 'pixels');
    ui.viewControls.showOptionsText = uicontrol(ui.viewControls.panel, ...
        'Style', 'text', ...
        'String', 'Display Options:', ...
        'Units', 'normalized', ...
        'Position', [0.00, 0.75, 1.00, 0.25]);
    ui.viewControls.showRawBtn = uicontrol(ui.viewControls.panel, ...
        'Style', 'togglebutton', ...
        'Tooltip', 'Show Raw Data', ...
        'String', 'Raw', ...
        'Units', 'normalized', ...
        'Position', [0.00, 0.50, 0.50, 0.25], ...
        'Callback', @(s,e) refresh_());
    ui.viewControls.showBaselineBtn = uicontrol(ui.viewControls.panel, ...
        'Style', 'togglebutton', ...
        'Tooltip', 'Show Baseline', ...
        'String', 'Baseline', ...
        'Units', 'normalized', ...
        'Position', [0.50, 0.50, 0.50, 0.25], ...
        'Callback', @(s,e) refresh_());
    ui.viewControls.showMaskedBtn = uicontrol(ui.viewControls.panel, ...
        'Style', 'togglebutton', ...
        'Tooltip', 'Show Masked Sweeps', ...
        'String', 'Masked', ...
        'Units', 'normalized', ...
        'Position', [0.00, 0.25, 0.50, 0.25], ...
        'Callback', @(s,e) refresh_());
    ui.viewControls.showAverageBtn = uicontrol(ui.viewControls.panel, ...
        'Style', 'togglebutton', ...
        'Tooltip', 'Show Group Average', ...
        'String', 'Average', ...
        'Units', 'normalized', ...
        'Position', [0.50, 0.25, 0.50, 0.25], ...
        'Callback', @(s,e) refresh_());
    ui.viewControls.showAverageOnlyBtn = uicontrol(ui.viewControls.panel, ...
        'Style', 'togglebutton', ...
        'Tooltip', 'Show Group Average Only', ...
        'String', 'Average Only', ...
        'Units', 'normalized', ...
        'Position', [0.00, 0.00, 1.00, 0.25], ...
        'Callback', @(s,e) refresh_());
end

function UNUSED_createBtnControls_()
    nrows = 5;
    ncols = 3;
    w = 1.0 / ncols;
    h = 1.0 / nrows;
    ui.btnControls.panel = uipanel(ui.mainWindow, ...
        'Units', 'pixels');
    y = 1-h;
    ui.btnControls.autoscaleXBtn = uicontrol(ui.btnControls.panel, ...
        'Style', 'pushbutton', ...
        'Tooltip', 'Autoscale X', ...
        'Units', 'normalized', ...
        'Position', [0, y, w, h], ...
        'Callback', @(s,e) autoscale_([], 'x'));
    ui.btnControls.autoscaleYBtn = uicontrol(ui.btnControls.panel, ...
        'Style', 'pushbutton', ...
        'Tooltip', 'Autoscale Y', ...
        'Units', 'normalized', ...
        'Position', [w, y, w, h], ...
        'Callback', @(s,e) autoscale_([], 'y'));
    ui.btnControls.autoscaleXYBtn = uicontrol(ui.btnControls.panel, ...
        'Style', 'pushbutton', ...
        'Tooltip', 'Autoscale XY', ...
        'Units', 'normalized', ...
        'Position', [1-w, y, w, h], ...
        'Callback', @(s,e) autoscale_([]));
    y = y-h;
    ui.btnControls.baselineFlatBtn = uicontrol(ui.btnControls.panel, ...
        'Style', 'pushbutton', ...
        'Tooltip', 'Baseline Flat', ...
        'Units', 'normalized', ...
        'Position', [0, y, w, h], ...
        'Callback', @(s,e) baselineFlat_());
    ui.btnControls.baselineSlopingBtn = uicontrol(ui.btnControls.panel, ...
        'Style', 'pushbutton', ...
        'Tooltip', 'Baseline Sloping Two Region', ...
        'cdata', [], ...
        'Units', 'normalized', ...
        'Position', [w, y, w, h], ...
        'Callback', @(s,e) baselineSloping_());
    ui.btnControls.baselineSplineBtn = uicontrol(ui.btnControls.panel, ...
        'Style', 'pushbutton', ...
        'Tooltip', 'Baseline Spline', ...
        'Units', 'normalized', ...
        'Position', [1-w, y, w, h], ...
        'Callback', @(s,e) baselineSplineInteractive_());
    y = y-h;
    ui.btnControls.normalizePositiveBtn = uicontrol(ui.btnControls.panel, ...
        'Style', 'pushbutton', ...
        'Tooltip', 'Normalize Positive', ...
        'Units', 'normalized', ...
        'Position', [0, y, w, h], ...
        'Callback', @(s,e) normalizePositive_());
    ui.btnControls.normalizeNegativeBtn = uicontrol(ui.btnControls.panel, ...
        'Style', 'pushbutton', ...
        'Tooltip', 'Normalize Negative', ...
        'Units', 'normalized', ...
        'Position', [w, y, w, h], ...
        'Callback', @(s,e) normalizeNegative_());
    ui.btnControls.normalizeAbsPeakBtn = uicontrol(ui.btnControls.panel, ...
        'Style', 'pushbutton', ...
        'Tooltip', 'Normalize Abs Peak', ...
        'Units', 'normalized', ...
        'Position', [1-w, y, w, h], ...
        'Callback', @(s,e) normalizeAbsPeak_());
    y = y-h;
    ui.btnControls.timeZeroBtn = uicontrol(ui.btnControls.panel, ...
        'Style', 'pushbutton', ...
        'Tooltip', 'Set Time Zero', ...
        'Units', 'normalized', ...
        'Position', [0, y, w, h]);
    ui.btnControls.alignToOnsetBtn = uicontrol(ui.btnControls.panel, ...
        'Style', 'pushbutton', ...
        'Tooltip', 'Align to Onset', ...
        'Units', 'normalized', ...
        'Position', [w, y, w, h]);
    y = y-h;
    ui.btnControls.zeroRegionBtn = uicontrol(ui.btnControls.panel, ...
        'Style', 'pushbutton', ...
        'Tooltip', 'Zero Selected Regions', ...
        'Units', 'normalized', ...
        'Position', [0, y, w, h]);
    ui.btnControls.interpolateRegionBtn = uicontrol(ui.btnControls.panel, ...
        'Style', 'pushbutton', ...
        'Tooltip', 'Interpolate Selected Regions', ...
        'Units', 'normalized', ...
        'Position', [w, y, w, h]);
    ui.btnControls.maskRegionBtn = uicontrol(ui.btnControls.panel, ...
        'Style', 'pushbutton', ...
        'Tooltip', 'Mask Selected Regions', ...
        'Units', 'normalized', ...
        'Position', [1-w, y, w, h]);
end

    function keyPress_(event)
        key = event.Key; %disp(key);
        if strcmp(key, 'leftarrow') disp('L');
        elseif strcmp(key, 'rightarrow') disp('R');
        elseif strcmp(key, 'uparrow') disp('U');
        elseif strcmp(key, 'downarrow') disp('D');
        end
    end

    function updateMainWindowTitle_()
        if ~isempty(data.traces)
            ui.mainWindow.Name = [data.meta.date ...
                '  [ ' data.meta.patchid ...
                ' ]  [ ' data.meta.construct ...
                ' ]  [ ' data.meta.experiment ' ]'];
        else
            ui.mainWindow.Name = 'Patch Meister';
        end
    end

    function nchannels = getNumChannels_()
        nchannels = size(data.traces,2);
    end

    function vischannels = getVisibleChannels_()
        vischannels = ui.visibleChannels.Value;
        vischannels(vischannels < 1) = [];
        vischannels(vischannels > getNumChannels_()) = [];
    end

    function ngroups = getNumGroups_()
        ngroups = max(1, numel(unique(data.groupids)));
    end

    function grouplabel = getGroupLabel_(groupind)
        if numel(data.grouplabels) >= groupind && ~isempty(data.grouplabels{groupind})
        	grouplabel = data.grouplabels{groupind};
        else
            grouplabel = ['Group ' num2str(groupind)];
        end
    end

    function traceind = getTracesInGroup_(groupind)
        groupuids = unique(data.groupids);
        traceind = find(data.groupids == groupuids(groupind));
    end

    function nsweepsmax = getMaxNumSweepsPerGroup_()
        groupuids = unique(data.groupids);
        ngroups = numel(groupuids);
        nsweepsmax = 0;
        for i = 1:ngroups
            nsweepsmax = max(nsweepsmax, numel(find(data.groupids == groupuids(i))));
        end
    end

    function tf = showRaw_(); tf = ui.showRawBtn.Value; end
    function tf = showBaseline_(); tf = ui.showBaselineBtn.Value; end
    function tf = showMasked_(); tf = ui.showMaskedBtn.Value; end
    function tf = showAverage_(); tf = ui.showAverageBtn.Value; end
    function tf = showAverageOnly_(); tf = ui.showAverageOnlyBtn.Value; end
    function tf = denoise_(); tf = ui.denoiseBtn.Value; end

    function about_()
        msgbox({ ...
            'Created by Marcel Goldschen-Ohm'; ...
            ''; ...
            'Email: goldschen-ohm@utexas.edu'; ...
            'Alt. Email: marcel.goldschen@gmail.com'; ...
            ''; ...
            'GitHub: https://github.com/marcel-goldschen-ohm/PatchMeister'; ...
            '' ...
            }, 'About Patch Meister');
    end

%% sweep selection
% Interface for traversing through sweeps.

    function ind = getSelectedSweeps_()
        ind = str2ind_(ui.selectedSweepsEdit.String);
        if isempty(ind)
            ind = 1:getMaxNumSweepsPerGroup_();
        end
        ind(ind < 1) = [];
        ind(ind > getMaxNumSweepsPerGroup_()) = [];
    end

    function setSelectedSweeps_(ind)
        ui.selectedSweepsEdit.String = ind2str_(ind);
        redraw_(true);
    end

    function selectPrevSweep_()
        ind = getSelectedSweeps_();
        if numel(ind) == 1 && ind > 1
            setSelectedSweeps_(ind - 1);
        elseif numel(ind) > 1
            setSelectedSweeps_(ind(end));
        end
    end

    function selectNextSweep_()
        ind = getSelectedSweeps_();
        if numel(ind) == 1 && ind < getMaxNumSweepsPerGroup_()
            setSelectedSweeps_(ind + 1);
        elseif numel(ind) > 1
            setSelectedSweeps_(ind(1));
        end
    end

    function selectAllSweeps_()
        setSelectedSweeps_(1:getMaxNumSweepsPerGroup_());
    end

    function ind = str2ind_(str)
        ind = [];
        fields = split(string(str), [",", " "]);
        for i = 1:numel(fields)
            field = strtrim(fields{i});
            if ~isempty(field)
                if ~isempty(strfind(field, "-"))
                    subfields = split(field, "-");
                    if numel(subfields) == 2
                        first = str2num(subfields{1});
                        last = str2num(subfields{2});
                        for j = first:last
                            ind = [ind, j];
                        end
                    end
                else
                    ind = [ind, str2num(field)];
                end
            end
        end
        ind = unique(ind);
    end

    function str = ind2str_(ind)
        str = '';
        ind = unique(ind);
        range = [];
        for i = 1:numel(ind)
            if isempty(range) || ind(i) == range(end) + 1
                range = [range, ind(i)];
            else
                if numel(range) == 1
                    if isempty(str)
                        str = num2str(range);
                    else
                        str = [str ',' num2str(range)];
                    end
                elseif numel(range) > 1
                    if isempty(str)
                        str = [num2str(range(1)) '-' num2str(range(end))];
                    else
                        str = [str ',' num2str(range(1)) '-' num2str(range(end))];
                    end
                end
                range = [ind(i)];
            end
        end
        if numel(range) == 1
            if isempty(str)
                str = num2str(range);
            else
                str = [str ',' num2str(range)];
            end
        elseif numel(range) > 1
            if isempty(str)
                str = [num2str(range(1)) '-' num2str(range(end))];
            else
                str = [str ',' num2str(range(1)) '-' num2str(range(end))];
            end
        end
    end

%% group sweeps

    function mergeGroups_()
        data.groupids(:) = 1;
        updateUI_();
        autoscale_();
        if isfield(ui, 'groupSweepsTable') && isvalid(ui.groupSweepsTable)
            ui.groupSweepsTable.Data(:,1) = data.groupids;
        end
    end

    function groupInterleaved_(ngroups)
        if ~exist('ngroups', 'var')
            answer = inputdlg('Group every:', 'Group Interleaved', 1, {'2'});
            if isempty(answer); return; end
            ngroups = str2num(answer{1});
        end
        for i = 1:ngroups
            data.groupids(i:ngroups:end) = i;
        end
        updateUI_();
        autoscale_();
        if isfield(ui, 'groupSweepsTable') && isvalid(ui.groupSweepsTable)
            ui.groupSweepsTable.Data(:,1) = data.groupids;
        end
    end

    function groupBlocks_(ntracesPerBlock)
        if ~exist('ntracesPerBlock', 'var')
            answer = inputdlg('Group blocks of:', 'Group Blocks', 1, {'2'});
            if isempty(answer); return; end
            ntracesPerBlock = str2num(answer{1});
        end
        ntraces = size(data.traces,1);
        gidx = 1;
        for i = 1:ntracesPerBlock:ntraces
            if i + ntracesPerBlock - 1 <= ntraces
                data.groupids(i:i+ntracesPerBlock-1) = gidx;
            else
                data.groupids(i:ntraces) = deal(gidx);
            end
            gidx = gidx + 1;
        end
        updateUI_();
        autoscale_();
        if isfield(ui, 'groupSweepsTable') && isvalid(ui.groupSweepsTable)
            ui.groupSweepsTable.Data(:,1) = data.groupids;
        end
    end

    function editGroups_()
        if isfield(ui, 'groupSweepsWindow') && isvalid(ui.groupSweepsWindow)
            figure(ui.groupSweepsWindow);
            return;
        end
        ui.groupSweepsWindow = figure('Name', 'Group Sweeps', ...
            'Menu', 'none', 'Toolbar', 'none', 'numbertitle', 'off');
        ui.groupSweepsTable = uitable(ui.groupSweepsWindow, ...
            'Data', data.groupids, ...
            'ColumnName', {'Group'}, ...
            'ColumnEditable', true, ...
            'Units', 'normalized', ...
            'Position', [0, 0, 1, 1], ...
            'CellEditCallback', @groupSweepsTableCellEdited_);
    end

    function groupSweepsTableCellEdited_(~, celldata)
        row = celldata.Indices(1);
        col = celldata.Indices(2);
        data.groupids(row) = celldata.NewData;
        updateUI_();
    end

    function editGroupLabels_()
        ngroups = getNumGroups_();
        for i = 1:ngroups
            fieldlabels{i} = ['Group ' num2str(i)];
            grouplabels{i} = getGroupLabel_(i);
        end
        answer = inputdlg(fieldlabels, 'Group Labels', 1, grouplabels);
        if isempty(answer); return; end
        data.grouplabels = answer;
        updateUI_();
    end

%% edit info

    function editInfo_()
        answer = inputdlg({'Date:', 'Patch/Cell ID:', 'Construct:', 'Experiment:'}, ...
            'Info', 1, ...
            {data.meta.date, data.meta.patchid, data.meta.construct, data.meta.experiment});
        if isempty(answer); return; end
        data.meta.date = answer{1};
        data.meta.patchid = answer{2};
        data.meta.construct = answer{3};
        data.meta.experiment = answer{4};
        updateMainWindowTitle_();
    end

    function editNotes_()
        if isfield(ui, 'notesWindow') && isvalid(ui.notesWindow)
            figure(ui.notesWindow);
            return;
        end
        ui.notesWindow = figure('Name', 'Notes', ...
            'Menu', 'none', 'Toolbar', 'none', 'numbertitle', 'off');
        ui.notesEdit = uicontrol(ui.notesWindow, ...
            'Style', 'edit', ...
            'Min', 0, ...
            'Max', 2, ...
            'String', data.meta.notes, ...
            'Units', 'normalized', ...
            'Position', [0, 0, 1, 1], ...
            'HorizontalAlignment', 'left', ...
            'Callback', @(s,e) updateNotesFromUI_());
    end

    function updateNotesFromUI_()
        if isfield(ui, 'notesEdit') && isvalid(ui.notesEdit)
            data.meta.notes = ui.notesEdit.String;
        end
    end

    function editChannels_()
        labels = {};
        channels = {};
        units = {};
        nchannels = getNumChannels_();
        for i = 1:nchannels
            labels{i} = ['Channel ' num2str(i)];
            channels{i} = char(data.traces(1,i).ylabel);
            units{i} = char(data.traces(1,i).yunit);
        end
        newchannels = inputdlg(labels, 'Channel Labels', 1, channels);
        if isempty(newchannels); return; end
        for i = 1:nchannels
            [data.traces(:,i).ylabel] = deal(string(newchannels{i}));
            labels{i} = [labels{i} ' -  ' newchannels{i}];
        end
        newunits = inputdlg(labels, 'Channel Units', 1, units);
        if isempty(newunits); return; end
        for i = 1:nchannels
            if string(units{i}) ~= string(newunits{i})
                [data.traces(:,i).yunit] = deal(string(newunits{i}));
                scaleFactor = inputdlg( ...
                    {['Scale factor for conversion from ' units{i} ' -> ' newunits{i}]}, ...
                    'Unit Conversion Factor', 1, {'1'});
                if ~isempty(scaleFactor)
                    scaleFactor = str2num(scaleFactor{1});
                    if scaleFactor ~= 1
                        for row = 1:size(data.traces,1)
                            data.traces(row,i).y = data.traces(row,i).y .* scaleFactor;
                        end
                    end
                end
            end
        end
        updateUI_();
    end

function UNUSED_editTraces_()
    if isfield(ui, 'editTracesWindow') && isvalid(ui.editTracesWindow)
        figure(ui.editTracesWindow);
        return;
    end
    ui.editTracesWindow = figure('Name', 'Traces', ...
        'Menu', 'none', 'Toolbar', 'none', 'numbertitle', 'off');

    tabledata = [num2cell([vertcat(data.groupids), ...
        vertcat(data.traces.x0), ...
        vertcat(data.traces.y0), ...
        vertcat(data.traces.yscale), ...
        ]), num2cell(vertcat(data.traces.ismasked))];

    uicontrol(ui.editTracesWindow, ...
        'Style', 'pushbutton', ...
        'String', 'Merge Groups', ...
        'Units', 'normalized', ...
        'Position', [0, 0.9, 0.2, 0.1], ...
        'Callback', @(s,e) mergeGroups_());

    uicontrol(ui.editTracesWindow, ...
        'Style', 'pushbutton', ...
        'String', 'Group Interleaved', ...
        'Units', 'normalized', ...
        'Position', [0.2, 0.9, 0.2, 0.1], ...
        'Callback', @(s,e) groupInterleaved_());

    uicontrol(ui.editTracesWindow, ...
        'Style', 'pushbutton', ...
        'String', 'Group Blocks', ...
        'Units', 'normalized', ...
        'Position', [0.4, 0.9, 0.2, 0.1], ...
        'Callback', @(s,e) groupBlocks_());

    ui.tracesInfoTable = uitable(ui.editTracesWindow, ...
        'Data', tabledata, ...
        'ColumnName', {'Group', 'X0', 'Y0', 'YScale', 'Masked'}, ...
        'ColumnEditable', true, ...
        'Units', 'normalized', ...
        'Position', [0, 0, 0.6, 0.9], ...
        'CellEditCallback', @tracesInfoTableCellEdited_, ...
        'CellSelectionCallback', @tracesInfoTableCellSelected_);

    ui.traceXYTable = uitable(ui.editTracesWindow, ...
        'Data', [], ...
        'ColumnName', {'X', 'Y'}, ...
        'ColumnEditable', true, ...
        'Units', 'normalized', ...
        'Position', [0.6, 0, 0.4, 1], ...
        'CellEditCallback', @traceXYTableCellEdited_);
end

function UNUSED_tracesInfoTableCellEdited_(~, celldata)
    i = celldata.Indices(1);
    j = celldata.Indices(2);
    switch j
        case 1
            data.traces(i).groupid = celldata.NewData;
            updateUI_();
        case 2
            data.traces(i).x0 = celldata.NewData;
        case 3
            data.traces(i).y0 = celldata.NewData;
        case 4
            data.traces(i).yscale = celldata.NewData;
        case 5
            data.traces(i).ismasked = celldata.NewData;
    end
    redraw_();
end

function UNUSED_tracesInfoTableCellSelected_(~, selection)
    if isfield(ui, 'traceXYTable') && isvalid(ui.traceXYTable)
        i = selection.Indices(1);
        ui.traceXYTable.UserData = i;
        ui.traceXYTable.Data = [data.traces(i).x, data.traces(i).y];
    end
end

function UNUSED_traceXYTableCellEdited_(~, celldata)
    trace = ui.editTraces.traceXYTable.UserData;
    i = celldata.Indices(1);
    j = celldata.Indices(2);
    switch j
        case 1
            data.traces(trace).x(i) = celldata.NewData;
        case 2
            data.traces(trace).y(i) = celldata.NewData;
    end
    redraw_();
end

%% i/o

    function clearData_()
        data = template.series;
    end

    function loadData_(filepath)
        if ~exist('filepath', 'var') || isempty(filepath)
            [file, path] = uigetfile(fullfile(ui.path, '*.mat'), 'Open data file.');
            if isequal(file, 0); return; end
            filepath = fullfile(path, file);
        end
        [ui.path, file, ext] = fileparts(filepath);
        if ext ~= ".mat"
            warndlg('Requires a *.mat file.', 'ERROR');
            return
        end
        wb = waitbar(0, 'Loading data file...');
        tmp = load(filepath);
        close(wb);
        file = strrep(file, '_', ' ');
        if ~isfield(tmp, 'data')
            warndlg('Unknown data format.', 'ERROR');
            return
        end
        clearData_();
        % a lot of this has checks for legacy data structures
        data.traces = repmat(template.trace, size(tmp.data.traces));
        for i = 1:size(tmp.data.traces,1)
            for j = 1:size(tmp.data.traces,2)
                data.traces(i,j).x = tmp.data.traces(i,j).x;
                data.traces(i,j).y = tmp.data.traces(i,j).y;
                data.traces(i,j).x0 = tmp.data.traces(i,j).x0;
                data.traces(i,j).y0 = tmp.data.traces(i,j).y0;
                data.traces(i,j).yscale = tmp.data.traces(i,j).yscale;
                data.traces(i,j).ismasked = tmp.data.traces(i,j).ismasked;
                if isfield(tmp.data.traces, 'masked')
                    data.traces(i,j).masked = tmp.data.traces(i,j).masked;
                end
                if isfield(tmp.data.traces, 'zeroed')
                    data.traces(i,j).zeroed = tmp.data.traces(i,j).zeroed;
                end
                if isfield(tmp.data.traces, 'interpolated')
                    data.traces(i,j).interpolated = tmp.data.traces(i,j).interpolated;
                end
                if isfield(tmp.data.traces, 'xlabel')
                    data.traces(i,j).xlabel = tmp.data.traces(i,j).xlabel;
                elseif isfield(tmp.data, 'xlabels')
                    data.traces(i,j).xlabel = string(tmp.data.xlabels{1});
                end
                if isfield(tmp.data.traces, 'xunit')
                    data.traces(i,j).xunit = tmp.data.traces(i,j).xunit;
                elseif isfield(tmp.data, 'xlabels')
                    data.traces(i,j).xunit = string(tmp.data.xlabels{2});
                end
                if isfield(tmp.data.traces, 'ylabel')
                    data.traces(i,j).ylabel = tmp.data.traces(i,j).ylabel;
                elseif isfield(tmp.data, 'ylabels')
                    data.traces(i,j).ylabel = string(tmp.data.ylabels{j,1});
                end
                if isfield(tmp.data.traces, 'yunit')
                    data.traces(i,j).yunit = tmp.data.traces(i,j).yunit;
                elseif isfield(tmp.data, 'ylabels')
                    data.traces(i,j).yunit = string(tmp.data.ylabels{j,2});
                end
            end
        end
        if isfield(tmp.data, 'groupids')
            data.groupids = tmp.data.groupids;
        elseif isfield(tmp.data.traces, 'groupid')
            data.groupids = vertcat(tmp.data.traces.groupid);
        else
            data.groupids = [1:size(data.traces,1)]';
        end
        if isfield(tmp.data, 'grouplabels')
            data.grouplabels = string(tmp.data.grouplabels);
        elseif isfield(tmp.data, 'groupnames')
            data.grouplabels = string(tmp.data.groupnames);
        else
            data.grouplabels = repmat("", [0,0]);
        end
        if isfield(tmp.data, 'meta') && isfield(tmp.data.meta, 'date')
            data.meta.date = tmp.data.meta.date;
        elseif isfield(tmp.data, 'info') && isKey(tmp.data.info, 'date')
            data.meta.date = tmp.data.info('date');
        end
        if isfield(tmp.data, 'meta') && isfield(tmp.data.meta, 'patchid')
            data.meta.patchid = tmp.data.meta.patchid;
        elseif isfield(tmp.data, 'info') && isKey(tmp.data.info, 'patchid')
            data.meta.patchid = tmp.data.info('patchid');
        end
        if isfield(tmp.data, 'meta') && isfield(tmp.data.meta, 'construct')
            data.meta.construct = tmp.data.meta.construct;
        elseif isfield(tmp.data, 'info') && isKey(tmp.data.info, 'construct')
            data.meta.construct = tmp.data.info('construct');
        end
        if isfield(tmp.data, 'meta') && isfield(tmp.data.meta, 'experiment')
            data.meta.experiment = tmp.data.meta.experiment;
        elseif isfield(tmp.data, 'info') && isKey(tmp.data.info, 'experiment')
            data.meta.experiment = tmp.data.info('experiment');
        end
        if isfield(tmp.data, 'meta') && isfield(tmp.data.meta, 'notes')
            data.meta.notes = tmp.data.meta.notes;
        elseif isfield(tmp.data, 'info') && isKey(tmp.data.info, 'notes')
            data.meta.notes = tmp.data.info('notes');
        end
%         else % assume PatchMaster .mat file export
%             traceLabels = fieldnames(tmp);
%             ntraces = numel(traceLabels);
%             data.traces = repmat(template.trace, [ntraces,1]);
%             for i = 1:ntraces
%                 xy = tmp.(traceLabels{i});
%                 data.traces(i,1).x = xy(:,1);
%                 data.traces(i,1).y = xy(:,2) .* 1e12; % A -> pA
%             end
%             data.xlabels = {'Time', 's'};
%             data.ylabels = {'Current', 'pA'};
%             data.groupids = [1:ntraces]';
%             data.grouplabels = {};
%             ind = strfind(file, ' ');
%             ind = [ind, repmat(length(file) + 1, [1, 3])];
%             data.meta.date = file(1:ind(1) - 1);
%             data.meta.patchid = file(ind(1) + 1:ind(2) - 1);
%             data.meta.construct = file(ind(2) + 1:ind(3) - 1);
%             data.meta.experiment = file(ind(3) + 1:end);
%         end
        updateUI_();
        setSelectedSweeps_(1);
        autoscale_();
    end

    function loadHEKA_(filepath)
        if ~exist('filepath', 'var') || isempty(filepath)
            [file, path] = uigetfile(fullfile(ui.path, '*.*'), 'Open HEKA data file.');
            if isequal(file, 0); return; end
            filepath = fullfile(path, file);
        end
        [ui.path, file, ext] = fileparts(filepath);
        % load HEKA data
        try
            heka = HEKA_Importer(filepath);
        catch
            warndlg("!!! Requires package 'HEKA Patchmaster Importer' by Christian Keine. Find in MATLAB's Add-On Explorer.", ...
                'HEKA file loader');
            return
        end
        nrecordings = size(heka.RecTable,1);
        % info for each recording are in the nonempty leaves of dataTree
        recdata = heka.trees.dataTree(:,end);
        clip = [];
        for i = 1:numel(recdata)
            if isempty(recdata{i}); clip = [clip i]; end
        end
        recdata(clip) = [];
        if numel(recdata) ~= nrecordings
            warndlg('Unexpected data structure. Please report this error.');
            return
        end
        if nrecordings > 1
            % Ask which recordings to load. Loading multiple recordings is
            % allowed provided they have the same channels.
            stimuli = {};
            for rec = 1:nrecordings
                stimuli{rec} = heka.RecTable.Stimulus{rec};
            end
            selrec = listdlg('ListString', stimuli, ...
                'PromptString', 'Select recordings to load:');
        else
            selrec = 1;
        end
        nchannels = numel(heka.RecTable.dataRaw{selrec(1)});
        for i = 2:numel(selrec)
            if nchannels ~= numel(heka.RecTable.dataRaw{selrec(i)}) ...
                    || ~isequal(heka.RecTable.ChName{selrec(1)}, heka.RecTable.ChName{selrec(i)})
                warndlg('Selected recordings do NOT have the same channels, and cannot be loaded together.');
                return
            end
        end
        clearData_();
        groupid = 1;
        for i = 1:numel(selrec)
            rec = selrec(i);
            if ~isempty(data.traces) && size(data.traces,2) ~= nchannels
                clearData_();
                warndlg('Selected recordings do NOT have the same channels, and cannot be loaded together.');
                return
            end
            nsweeps = size(heka.RecTable.dataRaw{rec}{1},2);
            npts = size(heka.RecTable.dataRaw{rec}{1},1);
            traces = repmat(template.trace, [nsweeps,nchannels]);
            for sweep = 1:nsweeps
                for channel = 1:nchannels
                    if sweep == 1 && channel == 1
                        traces(sweep,channel).x = [0:npts-1]' .* recdata{rec}.TrXInterval;
                    else
                        traces(sweep,channel).x = traces(1,1).x;
                    end
                    traces(sweep,channel).y = heka.RecTable.dataRaw{rec}{channel}(:,sweep);
                    traces(sweep,channel).xlabel = "Time";
                    traces(sweep,channel).xunit = string(heka.RecTable.TimeUnit{rec}{channel});
                    traces(sweep,channel).ylabel = string(heka.RecTable.ChName{rec}{channel}); 
                    traces(sweep,channel).yunit = string(heka.RecTable.ChUnit{rec}{channel});
                    traces(sweep,channel).timestamp = heka.RecTable.TimeStamp{rec}(sweep);
                end
            end
            data.traces = [data.traces; traces];
            data.groupids = [data.groupids; repmat(groupid, [nsweeps,1])];
            data.grouplabels(groupid) = string(heka.RecTable.Stimulus{rec});
            groupid = groupid + 1;
        end
        updateUI_();
        setSelectedSweeps_(1);
        autoscale_();
    end

    function loadABF_(filepath)
        if ~exist('filepath', 'var') || isempty(filepath)
            [file, path] = uigetfile(fullfile(ui.path, '*.*'), 'Open ABF or ABF2 data file.');
            if isequal(file, 0); return; end
            filepath = fullfile(path, file);
        end
        [ui.path, file, ext] = fileparts(filepath);
        try
            [d,si,h] = abfload(filepath); % data, sample interval (us), header
        catch
            warndlg("!!! Requires package 'fcollman/abfload' by Forrest Collman. Find in MATLAB's Add-On Explorer.", ...
                'pCLAMP ABF file loader');
            return
        end
        if iscell(d)
            nsweeps = numel(d);
            nchannels = size(d{1},2); % should be the same for all rows
            data.traces = repmat(template.trace, [nsweeps,nchannels]);
            for sweep = 1:nsweeps
                npts = size(d{sweep},1);
                for channel = 1:nchannels
                    if channel == 1
                        data.traces(sweep,channel).x = [0:npts-1]' .* (si * 1e-6); % time (s)
                    else
                        data.traces(sweep,channel).x = data.traces(sweep,1).x;
                    end
                    data.traces(sweep,channel).y = d{sweep}(:,channel);
                    data.traces(sweep,channel).xlabel = "Time";
                    data.traces(sweep,channel).xunit = "s";
                    data.traces(sweep,channel).ylabel = string(h.recChNames(channel));
                    data.traces(sweep,channel).yunit = string(h.recChUnits(channel));
                end
            end
        else
            nsweeps = size(d,3);
            nchannels = size(d,2);
            npts = size(d,1);
            data.traces = repmat(template.trace, [nsweeps,nchannels]);
            for sweep = 1:nsweeps
                for channel = 1:nchannels
                    if sweep == 1 && channel == 1
                        data.traces(sweep,channel).x = [0:npts-1]' .* (si * 1e-6); % time (s)
                    else
                        data.traces(sweep,channel).x = data.traces(1,1).x;
                    end
                    data.traces(sweep,channel).y = d(:,channel,sweep);
                    data.traces(sweep,channel).xlabel = "Time";
                    data.traces(sweep,channel).xunit = "s";
                    data.traces(sweep,channel).ylabel = string(h.recChNames(channel));
                    data.traces(sweep,channel).yunit = string(h.recChUnits(channel));
                end
            end
        end
        data.groupids = ones([nsweeps,1], 'uint32');
        data.grouplabels = repmat("", [0,0]);
        data.meta.date = '';
        data.meta.patchid = '';
        data.meta.construct = '';
        data.meta.experiment = '';
        data.meta.notes = '';
        updateUI_();
        setSelectedSweeps_(1);
        autoscale_();
    end

    function loadAxograph_(filepath)
        if ~exist('filepath', 'var') || isempty(filepath)
            [file, path] = uigetfile(fullfile(ui.path, '*.*'), 'Open Axograph or AxographX data file.');
            if isequal(file, 0); return; end
            filepath = fullfile(path, file);
        end
        [ui.path, file, ext] = fileparts(filepath);
        try
            axgr = readAxograph(filepath);
        catch
            warndlg("!!! Failed to load Axograph or AxographX file. Make sure readAxograph.m is in your MATLAB path.", ...
                'Axograph file loader');
            return
        end
        ncolumns = numel(axgr.columnData);
        % TODO: load data from axgr
        
        
        updateUI_();
        setSelectedSweeps_(1);
        autoscale_();
    end

    function saveData_(filepath)
        if ~exist('filepath', 'var') || isempty(filepath)
            default = fullfile(ui.path, [data.meta.date ...
                ' ' data.meta.patchid ...
                ' ' data.meta.construct ...
                ' ' data.meta.experiment ...
                '.mat']);
            [file, path] = uiputfile(default, 'Save data to file.');
            if isequal(file, 0); return; end
            filepath = fullfile(path, file);
        end
        [ui.path, file, ext] = fileparts(filepath);
        wb = waitbar(0, 'Saving data to file...');
        save(filepath, 'data');
        close(wb);
    end

end
