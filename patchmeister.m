function patchmeister()
%
% Created by Marcel Goldschen-Ohm <goldschen-ohm@utexas.edu, marcel.goldschen@gmail.com>

%% init
% for debugging
global data;
global ui;

% A trace is a list of (x,y) data points for some measurement vs. time,
% e.g. as obtained during episodic acquisition in an electrophysiology
% experiment.
%
% !!! The original raw (x,y) data should NEVER be altered unless you are
% absolutely sure of what you are doing. Instead, offsets and scale
% factors, as well as masked, zeroed and interpolated regions are stored
% separately, and can be used to reconstruct the adjusted (e.g. baselined,
% scaled, masked, etc.) trace from the raw data.
%
% Each trace is alloted its own x data to ease grouping traces with
% differing x values. Note, however, that assigning the x values of one
% trace to other traces will avoid copying the x data via Matlab's copy on
% write rules.
%
template.trace.x = []; % Tx1, e.g. time
template.trace.y = []; % Tx1, e.g. current, voltage, etc.
template.trace.x0 = 0; % 1x1 time zero offset
template.trace.y0 = 0; % 1x1 (uniform) OR Tx1 (nonuniform) baseline offset
template.trace.yscale = 1; % 1x1 (uniform) OR Tx1 (nonuniform) scale factor
template.trace.ismasked = false; % 1x1 logical flag for masking entire trace
template.trace.masked = []; % Tx1 logical for masked data points
template.trace.zeroed = []; % Tx1 logical for zeroed data points
template.trace.interpolated = []; % Tx1 logical for interpolated segments

% Each column is a channel (e.g. current, voltage, etc.)
% Each row is a sweep - a single data trace for each channel.
data.traces = repmat(template.trace, [0,0]); % sweeps x channels
data.xlabels = {}; % e.g. {'Time', 'ms'}
data.ylabels = {}; % (channels x 2), e.g. {'Current', 'pA'; 'Voltage', 'mV'}

% group sweeps
data.groupids = []; % (sweeps x 1) group index for each sweep
data.grouplabels = {}; % labels for each group

% experiment info
data.info = containers.Map;
data.info('date') = datestr(now, 'yyyy-mm-dd');
data.info('patchid') = ''; % cell/patch ID
data.info('construct') = ''; % construct, e.g. denote subunit composition and mutations
data.info('experiment') = ''; % very short one-line summary of experiment (more detail can go in notes)
data.info('notes') = '';

% UI
ui = struct();
initUI_();

%% custom user init
% loadData_('/Users/marcel/Box Sync/Goldschen-Ohm Lab/People/Wagner Nors/data/rGABAAR a1L263T-b2-g2/HEK293T/2019-02-01 T9-30 a1L263T-b2-g2 500ms 1mM PTX alt pipes A and B.mat');

%% test data
% data.traces = repmat(template.trace, [5,2]);
% data.groupids = [1; 1; 1; 2; 2];
% data.grouplabels = {'Group 1', 'Group 2'};
% for i_ = 1:size(data.traces,1)
%     for j_ = 1:size(data.traces,2)
%         data.traces(i_,j_).x = [1:100]';
%         data.traces(i_,j_).y = rand(100,1) .* 10000;
%     end
% end
% data.xlabels = {'Time', 'ms'};
% data.ylabels = {'Current', 'pA'; 'Voltage', 'mV'};
% updateUI_();

%% trace data
    function [x,y] = getTrace_(trace)
        % return the trace after applying offsets and scaling
        x = trace.x - trace.x0;
        y = (trace.y - trace.y0) .* trace.yscale;
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
        if any(trace.zeroed)
            y(trace.zeroed) = 0;
        end
        if any(trace.masked) && ui.showMaskedBtn.Checked == "off"
            y(trace.masked) = nan;
        end
    end

    function [x,y] = getAvgTrace_(traces, raw)
        % return average trace only for x range in which all traces overlap
        x = [];
        y = [];
        if isempty(traces); return; end
        if exist('raw', 'var') && raw
            x = traces(1).x;
            y = traces(1).y;
        else
            [x,y] = getTrace_(traces(1));
        end
        n = ones(size(y));
        n(isnan(y)) = 0;
        dx = min(diff(x));
        epsilon = 0.01 * dx;
        for i = 2:numel(traces)
            if exist('raw', 'var') && raw
                xi = traces(i).x;
                yi = traces(i).y;
            else
                [xi,yi] = getTrace_(traces(i));
            end
            dxi = min(diff(xi));
            if abs(dx - dxi) < epsilon
                % assume sweeps have the same sample interval
                ind = intersect(find(x >= xi(1) - epsilon), find(x <= xi(end) + epsilon));
                indi = intersect(find(xi >= x(1) - epsilon), find(xi <= x(end) + epsilon));
                if numel(ind) == numel(indi)
                    x = x(ind);
                    n = n(ind);
                    y = y(ind);
                    yi = yi(indi);
                    ynan = isnan(y);
                    yinan = isnan(yi);
                    notnan = ~ynan & ~yinan;
                    y(ynan) = yi(ynan);
                    n(ynan) = 1;
                    y(notnan) = y(notnan) + yi(notnan);
                    n(notnan) = n(notnan) + 1;
                    n(isnan(y)) = 0;
                else
                    disp('ERROR: Invalid overlapping indices for sweep average.');
                end
            else
                disp('ERROR: Sweeps have different sample intervals.');
            end
        end
        y = y ./ n;
    end

    function segments = breakIntoContiguousSegments_(pts)
        segments = {};
        segment = [pts(1)];
        for i = 2:numel(pts)
            if pts(i) == segment(end)+1
                segment(end+1) = pts(i);
            else
                segments{end+1} = segment;
                segment = [pts(i)];
            end
        end
        if ~isempty(segment)
            segments{end+1} = segment;
        end
    end

%% per axes visible trace manipulations
    function autoscale_(ax, xy)
        if ~exist('ax', 'var') || isempty(ax)
            ax = vertcat(ui.groups.ax);
        end
        if exist('xy', 'var')
            if xy == 'x'
                lims = [];
                for i = 1:numel(ax)
                    lims = [lims; axesXYLims_(ax(i))];
                end
                lims = [min(lims(:,1)), max(lims(:,2))];
                for i = 1:numel(ax)
                    ax(i).XLim = lims;
                end
            elseif xy == 'y'
                for i = 1:numel(ax)
                    lims = axesXYLims_(ax(i));
                    ax(i).YLim = lims(3:4);
                end
            end
            return
        end
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
        for i = 1:numel(ax)
            axis(ax(i), lims(i,:));
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
        for i = 1:numel(ax)
            rows = ax(i).UserData.rows;
            channel = ax(i).UserData.channel;
            if ui.showRawBtn.Checked == "on"
                [data.traces(rows,channel).y0] = deal(y0);
            else
                for j = 1:numel(rows)
                    data.traces(rows(j),channel).y0 = ...
                        data.traces(rows(j),channel).y0 - y0;
                end
            end
        end
        redraw_();
    end

function UNUSED_shiftBaselineVisibleTraces_(ax, dy0)
        if ~exist('ax', 'var') || isempty(ax)
            ax = vertcat(ui.groups.ax);
        end
        if ~exist('dy0', 'var')
            answer = inputdlg('Baseline:', 'Shift Baseline', 1, {'0'});
            if isempty(answer); return; end
            dy0 = str2num(answer{1});
        end
        for i = 1:numel(ax)
            rows = ax(i).UserData.rows;
            channel = ax(i).UserData.channel;
            for j = 1:numel(rows)
                data.traces(rows(j),channel).y0 = data.traces(rows(j),channel).y0 - dy0;
            end
        end
        redraw_();
    end

    function baselineFlatVisibleTraces_(ax)
        if ~exist('ax', 'var') || isempty(ax)
            ax = vertcat(ui.groups.ax);
        end
        didit = false;
        for i = 1:numel(ax)
            for j = 1:numel(ax(i).UserData.traces)
                brushdata = logical(ax(i).UserData.traces(j).BrushData);
                selpts = find(brushdata);
                if ~isempty(selpts)
                    row = ax(i).UserData.rows(j);
                    channel = ax(i).UserData.channel;
                    trace = data.traces(row,channel);
                    if ui.showRawBtn.Checked == "on"
                        data.traces(row,channel).y0 = mean(trace.y(brushdata));
                    else
                        [x,y] = getTrace_(trace);
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
        didit = false;
        for i = 1:numel(ax)
            for j = 1:numel(ax(i).UserData.traces)
                brushdata = logical(ax(i).UserData.traces(j).BrushData);
                selpts = find(brushdata);
                if ~isempty(selpts)
                    dsel = find(diff(selpts) > 1);
                    if numel(dsel) == 1
                        selpts1 = selpts(1:dsel);
                        selpts2 = selpts(dsel+1:end);
                        row = ax(i).UserData.rows(j);
                        channel = ax(i).UserData.channel;
                        trace = data.traces(row,channel);
                        if ui.showRawBtn.Checked == "on"
                            x1 = mean(trace.x(selpts1));
                            y1 = mean(trace.y(selpts1));
                            x2 = mean(trace.x(selpts2));
                            y2 = mean(trace.y(selpts2));
                            m = (y2 - y1) / (x2 - x1);
                            npts = numel(trace.y);
                            data.traces(row,channel).y0 = m .* (trace.x - x1) + y1;
                        else
                            [x,y] = getTrace_(trace);
                            x1 = mean(x(selpts1));
                            y1 = mean(y(selpts1));
                            x2 = mean(x(selpts2));
                            y2 = mean(y(selpts2));
                            m = (y2 - y1) / (x2 - x1);
                            npts = numel(y);
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
                selpts = find(brushdata);
                if ~isempty(selpts)
                    row = ax(i).UserData.rows(j);
                    channel = ax(i).UserData.channel;
                    trace = data.traces(row,channel);
                    %[x,y] = getTrace_(trace);
                    npts = numel(trace.y);
                    try
                        pp = splinefit(trace.x(selpts), trace.y(selpts), nsegments);
                        data.traces(row,channel).y0 = ppval(pp, trace.x);
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
        if ui.showRawBtn.Checked == "off"
            ui.showRawBtn.Checked = 'on';
            ui.showBaselineBtn.Checked = 'on';
            refresh_();
            autoscale_();
        elseif ui.showBaselineBtn.Checked == "off"
            ui.showBaselineBtn.Checked = 'on';
            redraw_();
        end
    end

    function setScaleVisibleTraces_(ax, yscale)
        if ~exist('ax', 'var') || isempty(ax)
            ax = vertcat(ui.groups.ax);
        end
        if ~exist('yscale', 'var')
            answer = inputdlg('Y Scale:', 'Set Y Scale', 1, {'1'});
            if isempty(answer); return; end
            yscale = str2num(answer{1});
        end
        for i = 1:numel(ax)
            rows = ax(i).UserData.rows;
            channel = ax(i).UserData.channel;
            if ui.showRawBtn.Checked == "on"
                [data.traces(rows,channel).yscale] = deal(yscale);
            else
                for j = 1:numel(rows)
                    data.traces(rows(j),channel).yscale = ...
                        data.traces(rows(j),channel).yscale .* yscale;
                end
            end
        end
        redraw_();
    end

function UNUSED_scaleVisibleTraces_(ax, yscale)
        if ~exist('ax', 'var') || isempty(ax)
            ax = vertcat(ui.groups.ax);
        end
        if ~exist('yscale', 'var')
            answer = inputdlg('Y Scale:', 'Set Y Scale', 1, {'1'});
            if isempty(answer); return; end
            yscale = str2num(answer{1});
        end
        for i = 1:numel(ax)
            rows = ax(i).UserData.rows;
            for j = 1:numel(rows)
                data.traces(rows(j),channel).yscale = data.traces(rows(j),channel).yscale .* yscale;
            end
        end
        redraw_();
    end

    function normalizePositiveVisibleTraces_(ax)
        if ~exist('ax', 'var') || isempty(ax)
            ax = vertcat(ui.groups.ax);
        end
        didit = false;
        for i = 1:numel(ax)
            for j = 1:numel(ax(i).UserData.traces)
                brushdata = logical(ax(i).UserData.traces(j).BrushData);
                selpts = find(brushdata);
                if ~isempty(selpts)
                    row = ax(i).UserData.rows(j);
                    channel = ax(i).UserData.channel;
                    trace = data.traces(row,channel);
                    if ui.showRawBtn.Checked == "on"
                        data.traces(row,channel).yscale = 1.0 / max(trace.y(selpts));
                    else
                        [x,y] = getTrace_(trace);
                        data.traces(row,channel).yscale = trace.yscale ./ max(y(selpts));
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
        didit = false;
        for i = 1:numel(ax)
            for j = 1:numel(ax(i).UserData.traces)
                brushdata = logical(ax(i).UserData.traces(j).BrushData);
                selpts = find(brushdata);
                if ~isempty(selpts)
                    row = ax(i).UserData.rows(j);
                    channel = ax(i).UserData.channel;
                    trace = data.traces(row,channel);
                    if ui.showRawBtn.Checked == "on"
                        data.traces(row,channel).yscale = 1.0 / abs(min(trace.y(selpts)));
                    else
                        [x,y] = getTrace_(trace);
                        data.traces(row,channel).yscale = trace.yscale ./ abs(min(y(selpts)));
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
        didit = false;
        for i = 1:numel(ax)
            for j = 1:numel(ax(i).UserData.traces)
                brushdata = logical(ax(i).UserData.traces(j).BrushData);
                selpts = find(brushdata);
                if ~isempty(selpts)
                    row = ax(i).UserData.rows(j);
                    channel = ax(i).UserData.channel;
                    trace = data.traces(row,channel);
                    if ui.showRawBtn.Checked == "on"
                        data.traces(row,channel).yscale = 1.0 / max(abs(trace.y(selpts)));
                    else
                        [x,y] = getTrace_(trace);
                        data.traces(row,channel).yscale = trace.yscale ./ max(abs(y(selpts)));
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
        for i = 1:numel(ax)
            rows = ax(i).UserData.rows;
            channel = ax(i).UserData.channel;
            if ui.showRawBtn.Checked == "on"
                [data.traces(rows,channel).x0] = deal(x0);
            else
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
        didit = false;
        for i = 1:numel(ax)
            for j = 1:numel(ax(i).UserData.traces)
                brushdata = logical(ax(i).UserData.traces(j).BrushData);
                selpts = find(brushdata);
                if ~isempty(selpts)
                    row = ax(i).UserData.rows(j);
                    channel = ax(i).UserData.channel;
                    trace = data.traces(row,channel);
                    if ui.showRawBtn.Checked == "on"
                        x = trace.x;
                        y = trace.y;
                    else
                        [x,y] = getTrace_(trace);
                    end
                    y = y - mean(y(selpts));
                    threshold = xSD * std(y(selpts));
                    if direction == "Absolute Value"
                        ind = find(abs(y(selpts(end)+1:end)) > threshold, 1);
                    elseif direction == "Positive"
                        ind = find(y(selpts(end)+1:end) > threshold, 1);
                    elseif direction == "Negative"
                        ind = find(y(selpts(end)+1:end) < -threshold, 1);
                    end
                    data.traces(row,channel).x0 = trace.x(selpts(end) + ind);
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

    function maskSelectedRegions_(ax)
        if ~exist('ax', 'var') || isempty(ax)
            ax = vertcat(ui.groups.ax);
        end
        didit = false;
        for i = 1:numel(ax)
            for j = 1:numel(ax(i).UserData.traces)
                brushdata = logical(ax(i).UserData.traces(j).BrushData);
                selpts = find(brushdata);
                if ~isempty(selpts)
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
                selpts = find(brushdata);
                if ~isempty(selpts)
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
                selpts = find(brushdata);
                if ~isempty(selpts)
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
            msgbox('Select data to mask with the brush tool.', ...
                'Interpolate Selected Region');
        end
    end

    function setIsMaskedVisibleTraces_(ax, ismasked)
        if ~exist('ax', 'var') || isempty(ax)
            ax = vertcat(ui.groups.ax);
        end
        for i = 1:numel(ax)
            rows = ax(i).UserData.rows;
            channel = ax(i).UserData.channel;
            [data.traces(rows,channel).ismasked] = deal(ismasked);
        end
        refresh_();
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
        
        ui.visibleChannels = uicontrol(ui.mainWindow, ...
            'Style', 'listbox', ...
            'Value', [], ...
            'Callback', @(s,e) updateUI_());

        loadIcons_();
        
        ui.menu = createMainMenu_();
        
        % sweep traversal
        ui.sweepspanel = uipanel(ui.mainWindow, ...
            'Units', 'pixels', ...
            'Position', [0, 0, 100, 70]);
        ui.prevSweepBtn = uicontrol(ui.sweepspanel, ...
            'Style', 'pushbutton', ...
            'Tooltip', 'Previous Sweep', ...
            'String', '<', ...
            'Units', 'pixels', ...
            'Position', [0, 20, 50, 50], ...
            'Callback', @(s,e) selectPrevSweep_());
        ui.nextSweepBtn = uicontrol(ui.sweepspanel, ...
            'Style', 'pushbutton', ...
            'Tooltip', 'Next Sweep', ...
            'String', '>', ...
            'Units', 'pixels', ...
            'Position', [50, 20, 50, 50], ...
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
        
        ui.groups = repmat(struct(), [0,0]);
        
        set(ui.mainWindow, 'KeyPressFcn', @(s,e) keyPress_(e));
        set(ui.mainWindow, 'SizeChangedFcn', @(s,e) resizeUI_());
        ui.path = pwd();

        resizeUI_();
    end

    function updateUI_()
        ngroups = getNumGroups_();
        nchannels = getNumChannels_();
        % update visible channels
        ui.visibleChannels.String = data.ylabels(:,1);
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
                    xlabel(ax, [data.xlabels{1} ' (' data.xlabels{2} ')']);
                end
                ylabel(ax, {getGroupLabel_(i); ...
                    [data.ylabels{channel,1} ' (' data.ylabels{channel,2} ')']});
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
        
        x = margin;
        y = h-margin;
        ui.visibleChannels.Position = [x, y-40, 100, 40];
        y = y-40-margin;
%         ui.viewControls.panel.Position = [5, y-70, 100, 70];
%         y = y-70;
        ui.sweepspanel.Position(1) = x;
        ui.sweepspanel.Position(2) = y - ui.sweepspanel.Position(4);
        y = y-ui.sweepspanel.Position(4)-margin;
%         ui.btnControls.panel.Position = [5, y-100*5.0/3, 100, 100*5.0/3];
        
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
                        if ui.showRawBtn.Checked == "on" % show raw data
                            if ui.showAverageOnlyBtn.Checked == "off" % show traces
                                for k = 1:nrows
                                    trace = data.traces(rows(k),channel);
                                    if ~trace.ismasked % not masked
                                        ax.UserData.traces = [ax.UserData.traces; ...
                                            plot(ax, trace.x, trace.y, 'color', cmap(1+mod(k-1,ncolors),:))];
                                        ax.UserData.rows = [ax.UserData.rows; rows(k)];
                                    elseif ui.showMaskedBtn.Checked == "on" % show masked traces
                                        ax.UserData.traces = [ax.UserData.traces; ...
                                            plot(ax, trace.x, trace.y, 'color', [0.5, 0.5, 0.5])];
                                        ax.UserData.rows = [ax.UserData.rows; rows(k)];
                                    end
                                    if ui.showBaselineBtn.Checked == "on" ... % show baseline
                                            && (~trace.ismasked || ui.showMaskedBtn.Checked == "on")
                                        y0 = trace.y0;
                                        if numel(y0) == 1
                                            y0 = repmat(y0, size(trace.y));
                                        end
                                        ax.UserData.baselines = [ax.UserData.baselines; ...
                                            plot(ax, trace.x, y0, 'k--')];
                                    end
                                end
                            end
                            if ui.showAverageBtn.Checked == "on" ...
                                    || ui.showAverageOnlyBtn.Checked == "on" % show average
                                [x, y] = getAvgTrace_(data.traces(rows, channel), true);
                                ax.UserData.avgtrace = plot(ax, x, y, 'k-');
                            end
                        else % show offset and scaled data
                            if ui.showAverageOnlyBtn.Checked == "off" % show traces
                                for k = 1:nrows
                                    trace = data.traces(rows(k),channel);
                                    if ~trace.ismasked
                                        [x, y] = getTrace_(trace);
                                        ax.UserData.traces = [ax.UserData.traces; ...
                                            plot(ax, x, y, 'color', cmap(1+mod(k-1,ncolors),:))];
                                        ax.UserData.rows = [ax.UserData.rows; rows(k)];
                                    elseif ui.showMaskedBtn.Checked == "on" % show masked traces
                                        [x, y] = getTrace_(trace);
                                        ax.UserData.traces = [ax.UserData.traces; ...
                                            plot(ax, x, y, 'color', [0.5, 0.5, 0.5])];
                                        ax.UserData.rows = [ax.UserData.rows; rows(k)];
                                    end
                                end
                            end
                            if ui.showAverageBtn.Checked == "on" ...
                                    || ui.showAverageOnlyBtn.Checked == "on" % show average
                                [x, y] = getAvgTrace_(data.traces(rows,channel));
                                ax.UserData.avgtrace = plot(ax, x, y, 'k-');
                            end
                            if ui.showBaselineBtn.Checked == "on" % show baseline
                                plot(ax, ax.XLim', zeros([2,1]), 'k--');
                            end
                        end
                    end
                else % refresh
                    for k = 1:numel(ax.UserData.traces)
                        trace = data.traces(ax.UserData.rows(k),channel);
                        if ui.showRawBtn.Checked == "on" % show raw data
                            ax.UserData.traces(k).XData(:) = trace.x;
                            ax.UserData.traces(k).YData(:) = trace.y;
                        else % show offset and scaled data
                            [x, y] = getTrace_(trace);
                            ax.UserData.traces(k).XData(:) = x;
                            ax.UserData.traces(k).YData(:) = y;
                        end
                    end
                    for k = 1:numel(ax.UserData.baselines)
                        trace = data.traces(ax.UserData.rows(k),channel);
                        if ui.showRawBtn.Checked == "on" % show raw data
                            y0 = trace.y0;
                            if numel(y0) == 1
                                y0 = repmat(y0, size(trace.y));
                            end
                            ax.UserData.baselines(k).XData(:) = trace.x;
                            ax.UserData.baselines(k).YData(:) = y0;
                        end
                    end
                    if ~isempty(ax.UserData.avgtrace)
                        if ui.showRawBtn.Checked == "on" % show raw data
                            [x, y] = getAvgTrace_(data.traces(ax.UserData.rows,channel), true);
                            ax.UserData.avgtrace.XData(:) = x;
                            ax.UserData.avgtrace.YData(:) = y;
                        else % show offset and scaled data
                            [x, y] = getAvgTrace_(data.sweeps(ax.UserData.rows,channel));
                            ax.UserData.avgtrace.XData(:) = x;
                            ax.UserData.avgtrace.YData(:) = y;
                        end
                    end
                end
            end
        end
    end

    function refresh_()
        redraw_(true);
    end

    function loadIcons_()
        ui.icons = containers.Map();
        ui.icons('menu') = imread('icons/menu.png');
%         ui.icons('autoscaleX') = imread('icons/autoscaleX.png');
%         ui.icons('autoscaleY') = imread('icons/autoscaleY.png');
%         ui.icons('autoscaleXY') = imread('icons/autoscaleXY.png');
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
            'Text', 'Load Data', ...
            'MenuSelectedFcn', @(s,e) loadData_());
        uimenu(menu, ...
            'Text', 'Save Data', ...
            'MenuSelectedFcn', @(s,e) saveData_());

        uimenu(menu, ...
            'Separator', 'on', ...
            'Text', 'Edit Patch Info', ...
            'MenuSelectedFcn', @(s,e) editInfo_());
        uimenu(menu, ...
            'Text', 'Show/Edit Notes', ...
            'MenuSelectedFcn', @(s,e) editNotes_());
        
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

        uimenu(menu, ...
            'Separator', 'on', ...
            'Label', 'Autoscale All XY', ...
            'Callback', @(s,e) autoscale_([]));
        uimenu(menu, ...
            'Label', 'Autoscale All X', ...
            'Callback', @(s,e) autoscale_([], 'x'));
        uimenu(menu, ...
            'Label', 'Autoscale All Y', ...
            'Callback', @(s,e) autoscale_([], 'y'));

        ui.showRawBtn = uimenu(menu, ...
            'Separator', 'on', ...
            'Label', 'Show Raw', ...
            'Callback', @toggleDisplayOptionBtn_);
        ui.showBaselineBtn = uimenu(menu, ...
            'Label', 'Show Baseline', ...
            'Callback', @toggleDisplayOptionBtn_);
        ui.showMaskedBtn = uimenu(menu, ...
            'Label', 'Show Masked', ...
            'Callback', @toggleDisplayOptionBtn_);
        ui.showAverageBtn = uimenu(menu, ...
            'Label', 'Show Average', ...
            'Callback', @toggleDisplayOptionBtn_);
        ui.showAverageOnlyBtn = uimenu(menu, ...
            'Label', 'Show Average ONLY', ...
            'Callback', @toggleDisplayOptionBtn_);
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
%         uimenu(menu, ...
%             'Label', 'Shift Baseline', ...
%             'Callback', @(s,e) shiftBaselineVisibleTraces_(ax));
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
            'Callback', @(s,e) setScaleVisibleTraces_(ax));
%         uimenu(menu, ...
%             'Label', 'Multiply Y By Scale', ...
%             'Callback', @(s,e) scaleVisibleTraces_(ax));
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

    function toggleDisplayOptionBtn_(btn, varargin)
        if btn.Checked == "on"
            btn.Checked = 'off';
        else
            btn.Checked = 'on';
        end
        refresh_();
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
            ui.mainWindow.Name = [data.info('date') ...
                '  [ ' data.info('patchid') ...
                ' ]  [ ' data.info('construct') ...
                ' ]  [ ' data.info('experiment') ' ]'];
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

    function traceind = getTracesInGroup(groupind)
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

%% sweep selection
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

function tf = UNUSED_isSweepMaskedForAllGroups_(sweepind)
        tf = true;
        groupuids = unique(data.groupids);
        for i = 1:numel(ui.groups)
            groupsweepind = find(data.groupids == groupuids(i));
            if numel(groupsweepind) >= sweepind ...
                    && ~isempty(find(horzcat(~data.traces(groupsweepind(sweepind),:).ismasked)))
                tf = false;
                return;
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
            {data.info('date'), data.info('patchid'), data.info('construct'), data.info('experiment')});
        if isempty(answer); return; end
        data.info('date') = answer{1};
        data.info('patchid') = answer{2};
        data.info('construct') = answer{3};
        data.info('experiment') = answer{4};
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
            'String', data.info('notes'), ...
            'Units', 'normalized', ...
            'Position', [0, 0, 1, 1], ...
            'HorizontalAlignment', 'left', ...
            'Callback', @(s,e) updateNotesFromUI_());
    end

    function updateNotesFromUI_()
        if isfield(ui, 'notesEdit') && isvalid(ui.notesEdit)
            data.info('notes') = ui.notesEdit.String;
        end
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
    function loadData_(filepath)
        if ~exist('filepath', 'var') || isempty(filepath)
            [file, path] = uigetfile(fullfile(ui.path, '*.*'), 'Open data file.');
            if isequal(file, 0); return; end
            filepath = fullfile(path, file);
        end
        [ui.path, file, ext] = fileparts(filepath);
        if ext == ".mat"
            wb = waitbar(0, 'Loading data file...');
            tmp = load(filepath);
            close(wb);
            file = strrep(file, '_', ' ');
            if isfield(tmp, 'data') % assume PatchMeister .mat file
                data.info = tmp.data.info;
                if ~isKey(data.info, 'notes')
                    data.info('notes') = '';
                end
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
                    end
                end
                if isfield(tmp.data, 'xlabels')
                    data.xlabels = tmp.data.xlabels;
                elseif isfield(tmp.data, 'units')
                    data.xlabels = {'Time', tmp.data.units{1}};
                else
                    data.xlabels = {'Time', 's'};
                end
                if isfield(tmp.data, 'ylabels')
                    data.ylabels = tmp.data.ylabels;
                elseif isfield(tmp.data, 'units')
                    data.ylabels = {'Current', tmp.data.units{2}};
                else
                    data.ylabels = {'Current', 'pA'};
                end
                if isfield(tmp.data, 'groupids')
                    data.groupids = tmp.data.groupids;
                elseif isfield(tmp.data.traces, 'groupid')
                    data.groupids = vertcat(tmp.data.traces.groupid);
                else
                    data.groupids = [1:size(data.traces,1)]';
                end
                if isfield(tmp.data, 'grouplabels')
                    data.grouplabels = tmp.data.grouplabels;
                elseif isfield(tmp.data, 'groupnames')
                    data.grouplabels = tmp.data.groupnames;
                else
                    data.grouplabels = {};
                end
            else % assume PatchMaster .mat file export
                traceLabels = fieldnames(tmp);
                ntraces = numel(traceLabels);
                data.traces = repmat(template.trace, [ntraces, 1]);
                for i = 1:ntraces
                    xy = tmp.(traceLabels{i});
                    data.traces(i,1).x = xy(:,1);
                    data.traces(i,1).y = xy(:,2) .* 1e12; % A -> pA
                end
                data.xlabels = {'Time', 's'};
                data.ylabels = {'Current', 'pA'};
                data.groupids = [1:ntraces]';
                data.grouplabels = {};
                ind = strfind(file, ' ');
                ind = [ind, repmat(length(file) + 1, [1, 3])];
                data.info('date') = file(1:ind(1) - 1);
                data.info('patchid') = file(ind(1) + 1:ind(2) - 1);
                data.info('construct') = file(ind(2) + 1:ind(3) - 1);
                data.info('experiment') = file(ind(3) + 1:end);
            end
        end
        updateUI_();
        setSelectedSweeps_(1);
        autoscale_();
    end

    function saveData_(filepath)
        if ~exist('filepath', 'var') || isempty(filepath)
            default = fullfile(ui.path, [data.info('date') ...
                ' ' data.info('patchid') ...
                ' ' data.info('construct') ...
                ' ' data.info('experiment') ...
                '.mat']);
            [file, path] = uiputfile(default, 'Save data to file.');
            if isequal(file, 0); return; end
            filepath = fullfile(path, file);
            ui.path = path;
        end
        wb = waitbar(0, 'Saving data to file...');
        save(filepath, 'data');
        close(wb);
    end

end
