function patch_meister()
% Created by Marcel Goldschen-Ohm <goldschen-ohm@utexas.edu>

%% init
template.trace = struct();
template.trace.x = [];
template.trace.y = [];
template.trace.groupid = 1;
template.trace.x0 = 0;
template.trace.y0 = 0;
template.trace.yscale = 1;
template.trace.ismasked = false;
template.trace.rois = {};

function [x, y] = getTrace(trace)
    x = trace.x - trace.x0;
	y = (trace.y - trace.y0) .* trace.yscale;
end

function [x, y] = getAverageTrace(traces)
    x = [];
    y = [];
    if isempty(traces); return; end
    [x, y] = getTrace(traces(1));
    dx = min(diff(x));
    dx = dx(1);
    epsilon = 0.01 * dx;
    for i = 2:numel(traces)
        [tx, ty] = getTrace(traces(i));
        idx = intersect(find(x >= tx(1) - epsilon), find(x <= tx(end) + epsilon));
        tidx = intersect(find(tx >= x(1) - epsilon), find(tx <= x(end) + epsilon));
        if numel(idx) == numel(tidx)
            x = x(idx);
            y = y(idx) + ty(tidx);
        else
            disp('ERROR: Invalid overlapping indexes for average.');
        end
    end
    y = y ./ numel(traces);
end

template.data = struct();
template.data.info = containers.Map;
template.data.info('date') = '';
template.data.info('patchid') = '';
template.data.info('construct') = '';
template.data.info('experiment') = '';
template.data.traces = [];
template.data.units = {'sec', 'pA'};
template.data.groupnames = {};

data = template.data;
olddata = data; % copy of data for undoing last operation

ui = struct();
initUi();
set(ui.mainWindow, 'KeyPressFcn', @keyPress);
set(ui.mainWindow, 'SizeChangedFcn', @resizeUi);
ui.path = pwd();

%% user init

%% UI
function initUi()
    ui.mainWindow = figure( ...
        'Name', 'Patch Meister', ...
        'Units', 'normalized', ...
        'Position', [0 0 1 1]);
    ui.mainWindow.Units = 'pixels';
    
    ui.menu = uimenu(ui.mainWindow, ...
        'Text', 'PatchMeister');
    
    uimenu(ui.menu, ...
        'Text', 'Load Data', ...
        'MenuSelectedFcn', @loadData);
    uimenu(ui.menu, ...
        'Text', 'Save Data', ...
        'MenuSelectedFcn', @saveData);
    
    uimenu(ui.menu, ...
        'Separator', 'on', ...
        'Text', 'Edit Patch Info', ...
        'MenuSelectedFcn', @editInfo);
    
    uimenu(ui.menu, ...
        'Separator', 'on', ...
        'Text', 'Undo Last Operation', ...
        'MenuSelectedFcn', @undo);
    uimenu(ui.menu, ...
        'Text', 'Baseline', ...
        'MenuSelectedFcn', @baseline);
    uimenu(ui.menu, ...
        'Text', 'Normalize', ...
        'MenuSelectedFcn', @normalize);
    uimenu(ui.menu, ...
        'Text', 'Align To Onset', ...
        'MenuSelectedFcn', @alignToOnset);
    
    ui.showTracesBtn = uimenu(ui.menu, ...
        'Separator', 'on', ...
        'Text', 'Show Traces', ...
        'Checked', 'on', ...
        'MenuSelectedFcn', @toggleBtnCheckedOption);
    ui.showMaskedTracesBtn = uimenu(ui.menu, ...
        'Text', 'Show Masked Traces', ...
        'Checked', 'off', ...
        'MenuSelectedFcn', @toggleBtnCheckedOption);
    ui.showSelectedTracesOnlyBtn = uimenu(ui.menu, ...
        'Text', 'Show Selected Traces Only', ...
        'Checked', 'on', ...
        'MenuSelectedFcn', @toggleBtnCheckedOption);
    ui.linkGroupSelectedTracesBtn = uimenu(ui.menu, ...
        'Text', 'Link Group Selected Traces', ...
        'Checked', 'on', ...
        'MenuSelectedFcn', @toggleBtnCheckedOption);
    ui.showAverageTraceBtn = uimenu(ui.menu, ...
        'Text', 'Show Average Trace', ...
        'Checked', 'off', ...
        'MenuSelectedFcn', @toggleBtnCheckedOption);
    
    uimenu(ui.menu, ...
        'Separator', 'on', ...
        'Text', 'Merge Groups', ...
        'MenuSelectedFcn', @mergeGroups);
    uimenu(ui.menu, ...
        'Text', 'Group Interleaved', ...
        'MenuSelectedFcn', @groupInterleavedTraces);
    uimenu(ui.menu, ...
        'Text', 'Group Adjacent', ...
        'MenuSelectedFcn', @groupAdjacentTraces);
    
    ax = axes('Parent', ui.mainWindow, ...
        'Units', 'pixels', ...
        'XLimMode', 'manual', ...
        'YLimMode', 'manual', ...
        'UserData', 1);
    plot(ax, 0, 0);
    ui.plots = [ax];
    ui.traces = gobjects(0);
    ui.visibleTraces = [];
    ui.selectedTraces = [];
    ui.plotcontrols = [getPlotControls()];
    ui.plotcontrols(end).panel.UserData = 1;
    
    ui.autoscaleXBtn = uicontrol(ui.mainWindow, ...
        'Style', 'pushbutton', ...
        'String', 'autoX', ...
        'Units', 'pixels', ...
        'Position', [10, 10, 50, 16], ...
        'CallBack', @(src, event) autoscalePlots(src, event, 'x'));
    ui.autoscaleYBtn = uicontrol(ui.mainWindow, ...
        'Style', 'pushbutton', ...
        'String', 'autoY', ...
        'Units', 'pixels', ...
        'Position', [60, 10, 50, 16], ...
        'CallBack', @(src, event) autoscalePlots(src, event, 'y'));
    ui.autoscaleXYBtn = uicontrol(ui.mainWindow, ...
        'Style', 'pushbutton', ...
        'String', 'autoXY', ...
        'Units', 'pixels', ...
        'Position', [110, 10, 50, 16], ...
        'CallBack', @autoscalePlots);
    
    resizeUi();
end

function toggleBtnCheckedOption(src, event)
    if strcmp(src.Checked, 'on')
        src.Checked = 'off';
    else
        src.Checked = 'on';
    end
    redraw();
end

function updateUi(src, event)
    if isempty(data.traces) || ~isfield(data.traces, 'groupid'); return; end
    ngroups = numel(unique(vertcat(data.traces.groupid)));
    while numel(ui.plots) < ngroups
        groupid = numel(ui.plots) + 1;
        ax = axes('Parent', ui.mainWindow, ...
            'Units', 'pixels', ...
            'XLimMode', 'manual', ...
            'YLimMode', 'manual', ...
            'UserData', groupid);
        plot(ax, 0, 0);
        ui.plots = [ui.plots; ax];
        ui.plotcontrols = [ui.plotcontrols; getPlotControls()];
        ui.plotcontrols(end).panel.UserData = groupid;
    end
    while numel(ui.plots) > ngroups
        delete(ui.plots(end));
        ui.plots(end) = [];
        delete(ui.plotcontrols(end).panel);
        ui.plotcontrols(end) = [];
    end
    for i = 1:ngroups
        ax = ui.plots(i);
        xlabel(ax, ['Time (' data.units{1} ')']);
        ylabel(ax, ['Current (' data.units{2} ')']);
        if numel(data.groupnames) >= i && ~isempty(data.groupnames{i})
            ui.plotcontrols(i).groupNameEdit.String = data.groupnames{i};
        else 
            ui.plotcontrols(i).groupNameEdit.String = ['Group ' num2str(i)];
        end
        ntraces = numel(find(vertcat(data.traces.groupid) == i));
        ui.plotcontrols(i).numTracesText.String = ['/' num2str(ntraces)];
    end
    linkaxes(ui.plots, 'x');
    resizeUi();
    redraw();
end

function resizeUi(src, event)
    winbox = ui.mainWindow.Position;
    winleft = 10;
    winbottom = 10;
    winwidth = winbox(3) - 20;
    winheight = winbox(4) - 20;
    wintop = winbottom + winheight;
    nplots = numel(ui.plots);
    plotsleft = winleft + 200;
    plotswidth = winwidth - 200;
    plotstop = wintop - 20;
    plotsbottom = winbottom + 50;
    plotsheight = plotstop - plotsbottom;
    plotheight = floor(double(plotsheight) / nplots);
    for row = 1:nplots
        ui.plotcontrols(row).panel.Position(1) = winleft;
        ui.plotcontrols(row).panel.Position(2) = plotstop - (row - 0.5) * plotheight ...
            - ui.plotcontrols(row).panel.Position(4) / 2;
        ui.plots(row).Position = [plotsleft, plotstop - row * plotheight, plotswidth, plotheight];
        ui.plots(row).TickLength = [0.005, 0.01];
    end
end

function redraw(src, event)
    ui.traces = gobjects(numel(data.traces), 1);
    ui.visibleTraces = [];
    for i = 1:numel(ui.plots)
        ax = ui.plots(i);
        cla(ax);
        hold(ax, 'on');
        cmap = colormap(lines());
        traceIdxs = find(vertcat(data.traces.groupid) == i);
        selectedTraceIdxs = traceIdxs(find(ui.selectedTraces(traceIdxs)));
        nonselectedTraceIdxs = setxor(traceIdxs, selectedTraceIdxs);
        if isempty(selectedTraceIdxs)
            ui.plotcontrols(i).selectedTracesEdit.String = ['1-' num2str(numel(traceIdxs))];
            ui.plotcontrols(i).maskTraceCheckBox.Visible = false;
            %ui.plotcontrols(i).editTraceInfoBtn.Visible = false;
        elseif numel(selectedTraceIdxs) == 1
            ui.plotcontrols(i).selectedTracesEdit.String = num2str(find(traceIdxs == selectedTraceIdxs));
            ui.plotcontrols(i).maskTraceCheckBox.Value = data.traces(selectedTraceIdxs).ismasked;
            ui.plotcontrols(i).maskTraceCheckBox.Visible = true;
            %ui.plotcontrols(i).editTraceInfoBtn.Visible = true;
        else
            ui.plotcontrols(i).maskTraceCheckBox.Visible = false;
            %ui.plotcontrols(i).editTraceInfoBtn.Visible = false;
        end
        if strcmp(ui.showTracesBtn.Checked, 'on') ...
                && (strcmp(ui.showSelectedTracesOnlyBtn.Checked, 'off') || isempty(selectedTraceIdxs))
            % non-selected traces
            for j = 1:numel(nonselectedTraceIdxs)
                traceIdx = nonselectedTraceIdxs(j);
                tidx = find(traceIdxs == traceIdx, 1);
                [x, y] = getTrace(data.traces(traceIdx));
                if ~data.traces(traceIdx).ismasked
                    ui.traces(traceIdx) = plot(ax, x, y, 'color', cmap(tidx, :));
                    ui.visibleTraces = [ui.visibleTraces; traceIdx];
                elseif strcmp(ui.showMaskedTracesBtn.Checked, 'on')
                    ui.traces(traceIdx) = plot(ax, x, y, 'color', [0.5, 0.5, 0.5]);
                    ui.visibleTraces = [ui.visibleTraces; traceIdx];
                end
            end
        end
        if strcmp(ui.showAverageTraceBtn.Checked, 'on')
            % trace average
            avgTraceIdxs = traceIdxs;
            masked = find(vertcat(data.traces(avgTraceIdxs).ismasked));
            avgTraceIdxs(masked) = [];
            if ~isempty(avgTraceIdxs)
                [x, y] = getAverageTrace(data.traces(avgTraceIdxs));
                plot(ax, x, y, 'color', [0, 0, 0]);
            end
        end
        if strcmp(ui.showTracesBtn.Checked, 'on')
            % selected traces
            for j = 1:numel(selectedTraceIdxs)
                traceIdx = selectedTraceIdxs(j);
                tidx = find(traceIdxs == traceIdx, 1);
                [x, y] = getTrace(data.traces(traceIdx));
                if ~data.traces(traceIdx).ismasked
                    ui.traces(traceIdx) = plot(ax, x, y, 'color', cmap(tidx, :), 'linewidth', 1.5);
                    ui.visibleTraces = [ui.visibleTraces; traceIdx];
                elseif strcmp(ui.showMaskedTracesBtn.Checked, 'on')
                    ui.traces(traceIdx) = plot(ax, x, y, 'color', [0.5, 0.5, 0.5], 'linewidth', 1.5);
                    ui.visibleTraces = [ui.visibleTraces; traceIdx];
                end
            end
        end
    end
end

function autoscalePlots(src, event, xy)
    if exist('xy', 'var')
        if xy == 'x'
            lims = [];
            for i = 1:numel(ui.plots)
                ax = ui.plots(i);
                lims = [lims; axesXYLims(ax)];
            end
            lims = [min(lims(:,1)), max(lims(:,2))];
            for i = 1:numel(ui.plots)
                ax = ui.plots(i);
                ax.XLim = lims;
            end
        elseif xy == 'y'
            for i = 1:numel(ui.plots)
                ax = ui.plots(i);
                lims = axesXYLims(ax);
                ax.YLim = lims(3:4);
            end
        end
        return
    end
    lims = [];
    for i = 1:numel(ui.plots)
        ax = ui.plots(i);
        lims = [lims; axesXYLims(ax)];
    end
    lims(:,1) = min(lims(:,1));
    lims(:,2) = max(lims(:,2));
    for i = 1:numel(ui.plots)
        ax = ui.plots(i);
        axis(ax, lims(i,:));
    end
end

%% Plot control panels
function controls = getPlotControls()
    lh = 16;
    bh = 40;
    w = 150;
    h = 10 + 3 * lh + bh;
    controls.panel = uipanel(ui.mainWindow, ...
        'Units', 'pixels', ...
        'Position', [0, 0, w, h]);
    t = h - 5;
    w = w - 10;
    controls.groupNameEdit = uicontrol(controls.panel, ...
        'Style', 'edit', ...
        'Units', 'pixels', ...
        'Position', [5 t-lh, w, lh], ...
        'CallBack', @groupNameEditChanged);
    t = t - lh;
    controls.selectedTracesEdit = uicontrol(controls.panel, ...
        'Style', 'edit', ...
        'Units', 'pixels', ...
        'Position', [5, t-lh, 70, lh], ...
        'CallBack', @selectedTracesEditChanged);
    controls.numTracesText = uicontrol(controls.panel, ...
        'Style', 'text', ...
        'String', '/0', ...
        'HorizontalAlignment', 'left', ...
        'Units', 'pixels', ...
        'Position', [75, t-lh, 35, lh]);
    controls.showAllTracesBtn = uicontrol(controls.panel, ...
        'Style', 'pushbutton', ...
        'String', 'all', ...
        'Units', 'pixels', ...
        'Position', [110, t-lh, 35, lh], ...
        'CallBack', @allTraces);
    t = t - lh;
    controls.firstTraceBtn = uicontrol(controls.panel, ...
        'Style', 'pushbutton', ...
        'String', '<<', ...
        'Units', 'pixels', ...
        'Position', [5, t-bh, 20, bh], ...
        'CallBack', @firstTrace);
    controls.prevTraceBtn = uicontrol(controls.panel, ...
        'Style', 'pushbutton', ...
        'String', '<', ...
        'Units', 'pixels', ...
        'Position', [25, t-bh, 50, bh], ...
        'CallBack', @prevTrace);
    controls.nextTraceBtn = uicontrol(controls.panel, ...
        'Style', 'pushbutton', ...
        'String', '>', ...
        'Units', 'pixels', ...
        'Position', [75, t-bh, 50, bh], ...
        'CallBack', @nextTrace);
    controls.lastTraceBtn = uicontrol(controls.panel, ...
        'Style', 'pushbutton', ...
        'String', '>>', ...
        'Units', 'pixels', ...
        'Position', [125, t-bh, 20, bh], ...
        'CallBack', @lastTrace);
    t = t - bh;
    controls.maskTraceCheckBox = uicontrol(controls.panel, ...
        'Style', 'checkbox', ...
        'String', 'Masked', ...
        'Units', 'pixels', ...
        'Position', [5, t-lh, w, lh], ...
        'CallBack', @toggleMaskTrace);
    controls.editTraceInfoBtn = uicontrol(controls.panel, ...
        'Style', 'pushbutton', ...
        'String', 'info', ...
        'Units', 'pixels', ...
        'Position', [110, t-lh, 35, lh], ...
        'CallBack', @editTraceInfo);
    t = t - lh;
end

function groupNameEditChanged(src, event)
    groupid = src.Parent.UserData;
    data.groupnames{groupid} = src.String;
    if isempty(src.String)
        src.String = ['Group ' num2str(groupid)];
    end
end

function selectedTracesEditChanged(src, event)
    groupid = src.Parent.UserData;
    gidx = find(vertcat(data.traces.groupid) == groupid);
    idx = str2idx(src.String);
    ui.selectedTraces(gidx) = false;
    idx(find(idx < 1)) = [];
    idx(find(idx > numel(gidx))) = [];
    if numel(idx) < numel(gidx) && numel(idx) > 0
        ui.selectedTraces(gidx(idx)) = true;
    end
    if strcmp(ui.linkGroupSelectedTracesBtn.Checked, 'on')
        propagateSelectedTracesToAllOtherGroups(groupid);
    end
    redraw();
end

function prevTrace(src, event)
    groupid = src.Parent.UserData;
    gidx = find(vertcat(data.traces.groupid) == groupid);
    if strcmp(ui.showMaskedTracesBtn.Checked, 'off')
        masked = vertcat(data.traces(gidx).ismasked);
        idx = gidx(find(~masked));
    else
        idx = gidx;
    end
    if isempty(idx); return; end
    hdx = find(ui.selectedTraces(idx));
    ui.selectedTraces(gidx) = false;
    if isempty(hdx)
        idx = idx(end);
    elseif hdx(1) == 1
        idx = idx(1);
    else
        idx = idx(hdx(1) - 1);
    end
    ui.selectedTraces(idx) = true;
    if strcmp(ui.linkGroupSelectedTracesBtn.Checked, 'on')
        propagateSelectedTracesToAllOtherGroups(groupid);
    end
    redraw();
end

function nextTrace(src, event)
    groupid = src.Parent.UserData;
    gidx = find(vertcat(data.traces.groupid) == groupid);
    if strcmp(ui.showMaskedTracesBtn.Checked, 'off')
        masked = vertcat(data.traces(gidx).ismasked);
        idx = gidx(find(~masked));
    else
        idx = gidx;
    end
    if isempty(idx); return; end
    hdx = find(ui.selectedTraces(idx));
    ui.selectedTraces(gidx) = false;
    if isempty(hdx)
        idx = idx(1);
    elseif hdx(end) == numel(idx)
        idx = idx(end);
    else
        idx = idx(hdx(end) + 1);
    end
    ui.selectedTraces(idx) = true;
    if strcmp(ui.linkGroupSelectedTracesBtn.Checked, 'on')
        propagateSelectedTracesToAllOtherGroups(groupid);
    end
    redraw();
end

function firstTrace(src, event)
    groupid = src.Parent.UserData;
    gidx = find(vertcat(data.traces.groupid) == groupid);
    if strcmp(ui.showMaskedTracesBtn.Checked, 'off')
        masked = vertcat(data.traces(gidx).ismasked);
        idx = gidx(find(~masked));
    else
        idx = gidx;
    end
    if isempty(idx); return; end
    ui.selectedTraces(gidx) = false;
    ui.selectedTraces(idx(1)) = true;
    if strcmp(ui.linkGroupSelectedTracesBtn.Checked, 'on')
        propagateSelectedTracesToAllOtherGroups(groupid);
    end
    redraw();
end

function lastTrace(src, event)
    groupid = src.Parent.UserData;
    gidx = find(vertcat(data.traces.groupid) == groupid);
    if strcmp(ui.showMaskedTracesBtn.Checked, 'off')
        masked = vertcat(data.traces(gidx).ismasked);
        idx = gidx(find(~masked));
    else
        idx = gidx;
    end
    if isempty(idx); return; end
    ui.selectedTraces(gidx) = false;
    ui.selectedTraces(idx(end)) = true;
    if strcmp(ui.linkGroupSelectedTracesBtn.Checked, 'on')
        propagateSelectedTracesToAllOtherGroups(groupid);
    end
    redraw();
end

function allTraces(src, event)
    groupid = src.Parent.UserData;
    gidx = find(vertcat(data.traces.groupid) == groupid);
    ui.selectedTraces(gidx) = false;
    if strcmp(ui.linkGroupSelectedTracesBtn.Checked, 'on')
        ui.selectedTraces(:) = false;
    end
    redraw();
end

function toggleMaskTrace(src, event)
    if strcmp(ui.showTracesBtn.Checked, 'off'); return; end
    groupid = src.Parent.UserData;
    gidx = find(vertcat(data.traces.groupid) == groupid);
    if strcmp(ui.showMaskedTracesBtn.Checked, 'off')
        masked = vertcat(data.traces(gidx).ismasked);
        idx = gidx(find(~masked));
    else
        idx = gidx;
    end
    if isempty(idx); return; end
    hdx = find(ui.selectedTraces(idx));
    if numel(hdx) == 1
        idx = idx(hdx);
        data.traces(idx).ismasked = ~data.traces(idx).ismasked;
        redraw();
    end
end

function editTraceInfo(src, event)
    if strcmp(ui.showTracesBtn.Checked, 'off'); return; end
    groupid = src.Parent.UserData;
    gidx = find(vertcat(data.traces.groupid) == groupid);
    if strcmp(ui.showMaskedTracesBtn.Checked, 'off')
        masked = vertcat(data.traces(gidx).ismasked);
        idx = gidx(find(~masked));
    else
        idx = gidx;
    end
    if isempty(idx); return; end
    hdx = find(ui.selectedTraces(idx));
    if numel(hdx) == 1
        idx = idx(hdx);
        x0 = num2str(data.traces(idx).x0);
        y0 = data.traces(idx).y0;
        if numel(y0) > 1; y0 = ''; else y0 = num2str(y0); end
        yscale = data.traces(idx).yscale;
        if numel(yscale) > 1; yscale = ''; else yscale = num2str(yscale); end
        answer = inputdlg({'X zero:', 'Y zero:', 'Y scale:'}, ...
            ['Trace ' num2str(hdx)], 1, ...
            {x0, y0, yscale});
        x0 = answer{1};
        y0 = answer{2};
        yscale = answer{3};
        if ~isempty(x0); data.traces(idx).x0 = str2num(x0); end
        if ~isempty(y0); data.traces(idx).y0 = str2num(y0); end
        if ~isempty(yscale); data.traces(idx).yscale = str2num(yscale); end
        redraw();
    elseif isempty(hdx)
        % apply to all traces in group
        x0 = unique(vertcat(data.traces(idx).x0));
        if numel(x0) ~= 1; x0 = ''; else x0 = num2str(x0); end
        y0 = [];
        for i = 1:numel(idx)
            if numel(data.traces(idx(i)).y0) == 1
                y0 = [y0, data.traces(idx(i)).y0];
            else
                y0 = [];
                break
            end
        end
        y0 = unique(y0);
        if numel(y0) ~= 1; y0 = ''; else y0 = num2str(y0); end
        yscale = [];
        for i = 1:numel(idx)
            if numel(data.traces(idx(i)).yscale) == 1
                yscale = [yscale, data.traces(idx(i)).yscale];
            else
                yscale = [];
                break
            end
        end
        yscale = unique(yscale);
        if numel(yscale) ~= 1; yscale = ''; else yscale = num2str(yscale); end
        answer = inputdlg({'X zero:', 'Y zero:', 'Y scale:'}, ...
            ['Traces 1-' num2str(numel(idx))], 1, ...
            {x0, y0, yscale});
        x0 = str2num(answer{1});
        y0 = str2num(answer{2});
        yscale = str2num(answer{3});
        for i = 1:numel(idx)
            if ~isempty(x0); data.traces(idx(i)).x0 = data.traces(idx(i)).x0 + x0; end
            if ~isempty(y0); data.traces(idx(i)).y0 = data.traces(idx(i)).y0 + y0; end
            if ~isempty(yscale); data.traces(idx(i)).yscale = data.traces(idx(i)).yscale .* yscale; end
        end
        redraw();
    end
end

%% I/O
function loadData(src, event, filepath)
    if ~exist('filepath', 'var') || isempty(filepath)
        [file, path] = uigetfile(fullfile(ui.path, '*.mat'), 'Open data file.');
        if isequal(file, 0); return; end
        filepath = fullfile(path, file);
    end
    wb = waitbar(0, 'Loading data file...');
    tmp = load(filepath);
    close(wb);
    [ui.path, file, ext] = fileparts(filepath);
    file = strrep(file, '_', ' ');
    if isfield(tmp, 'data') % assume PatchMeister .mat file
        data = tmp.data;
    else % assume PatchMaster .mat file export
        data = template.data;
        traceLabels = fieldnames(tmp);
        ntraces = length(traceLabels);
        data.traces = repmat(template.trace, [ntraces, 1]);
        for i = 1:ntraces
            xy = tmp.(traceLabels{i});
            data.traces(i).x = xy(:,1);
            data.traces(i).y = xy(:,2) .* 1e12; % A -> pA
        end
        idx = strfind(file, ' ');
        idx = [idx, repmat(length(file) + 1, [1, 3])];
        data.info('date') = file(1:idx(1) - 1);
        data.info('patchid') = file(idx(1) + 1:idx(2) - 1);
        data.info('construct') = file(idx(2) + 1:idx(3) - 1);
        data.info('experiment') = file(idx(3) + 1:end);
    end
    olddata = data;
    figure(ui.mainWindow);
    title(ui.plots(1), [data.info('date') ...
        ' | PatchID: ' data.info('patchid') ...
        ' | Construct: ' data.info('construct') ...
        ' | Exp: ' data.info('experiment')]);
    ui.selectedTraces = false(size(data.traces));
	updateUi();
    autoscalePlots();
end

function saveData(src, event, filepath)
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

function editInfo(src, event)
    answer = inputdlg({'Date:', 'Patch ID:', 'Construct:', 'Experiment:'}, ...
        'Info', 1, ...
        {data.info('date'), data.info('patchid'), data.info('construct'), data.info('experiment')});
    data.info('date') = answer{1};
    data.info('patchid') = answer{2};
    data.info('construct') = answer{3};
    data.info('experiment') = answer{4};
    figure(ui.mainWindow);
    title(ui.plots(1), [data.info('date') ...
        ' | PatchID: ' data.info('patchid') ...
        ' | Construct: ' data.info('construct') ...
        ' | Exp: ' data.info('experiment')]);
end

%% Groups
function mergeGroups(src, event)
    olddata = data;
    [data.traces(:).groupid] = deal(1);
    ui.selectedTraces(:) = false;
	updateUi();
    autoscalePlots();
end

function groupInterleavedTraces(src, event, ngroups)
    if ~exist('ngroups', 'var')
        answer = inputdlg('Group every:', 'Group Interleaved', 1, {'2'});
        ngroups = str2num(answer{1});
    end
    olddata = data;
    for i = 1:ngroups
        [data.traces(i:ngroups:end).groupid] = deal(i);
    end
    ui.selectedTraces(:) = false;
	updateUi();
    autoscalePlots();
end

function groupAdjacentTraces(src, event, blocksize)
    if ~exist('blocksize', 'var')
        answer = inputdlg('Group blocks of:', 'Group Adjacent', 1, {'2'});
        blocksize = str2num(answer{1});
    end
    olddata = data;
    ntraces = numel(data.traces);
    gidx = 1;
    for i = 1:blocksize:ntraces
        if i + blocksize - 1 <= ntraces
            [data.traces(i:i + blocksize - 1).groupid] = deal(gidx);
        else
            [data.traces(i:end).groupid] = deal(gidx);
        end
        gidx = gidx + 1;
    end
    ui.selectedTraces(:) = false;
	updateUi();
    autoscalePlots();
end

%% Data manipulations
function undo(src, event)
    data = olddata;
    updateUi();
end

function baseline(src, event, type)
    if ~exist('type', 'var')
        list = {'Flat', 'Sloping Two Region'};
        idx = listdlg('ListString', list, 'SelectionMode', 'single');
        if isempty(idx); return; end
        type = list{idx};
    end
    olddata = data;
    for vi = 1:numel(ui.visibleTraces)
        i = ui.visibleTraces(vi);
        sel = find(logical(ui.traces(i).BrushData));
        if ~isempty(sel)
            if strcmp(type, 'Flat')
                data.traces(i).y0 = mean(data.traces(i).y(sel));
            elseif strcmp(type, 'Sloping Two Region')
                dsel = find(diff(sel) > 1);
                if numel(dsel) == 1
                    sel1 = sel(1:dsel);
                    sel2 = sel(dsel + 1:end);
                    x1 = mean(data.traces(i).x(sel1));
                    y1 = mean(data.traces(i).y(sel1));
                    x2 = mean(data.traces(i).x(sel2));
                    y2 = mean(data.traces(i).y(sel2));
                    m = (y2 - y1) / (x2 - x1);
                    data.traces(i).y0 = m .* (data.traces(i).x - x1) + y1;
                end
            end
        end
    end
    redraw();
end

function normalize(src, event, direction)
    if ~exist('direction', 'var')
        list = {'Absolute Value', 'Positive', 'Negative'};
        idx = listdlg('ListString', list, 'SelectionMode', 'single');
        if isempty(idx); return; end
        direction = list{idx};
    end
    olddata = data;
    for vi = 1:numel(ui.visibleTraces)
        i = ui.visibleTraces(vi);
        sel = find(logical(ui.traces(i).BrushData));
        if ~isempty(sel)
            y = data.traces(i).y - data.traces(i).y0;
            if strcmp(direction, 'Absolute Value')
                peak = max(abs(y(sel)));
            elseif strcmp(direction, 'Positive')
                peak = max(y(sel));
            elseif strcmp(direction, 'Negative')
                peak = abs(min(y(sel)));
            end
            data.traces(i).yscale = 1.0 / peak(1);
        end
    end
    redraw();
end

function alignToOnset(src, event, xSD, direction)
    if ~exist('xSD', 'var')
        answer = inputdlg('x Standard Deviation:', 'Align to Onset', 1, {'8'});
        xSD = str2num(answer{1});
    end
    if ~exist('direction', 'var')
        list = {'Absolute Value', 'Positive', 'Negative'};
        idx = listdlg('ListString', list, 'SelectionMode', 'single');
        if isempty(idx); return; end
        direction = list{idx};
    end
    olddata = data;
    for vi = 1:numel(ui.visibleTraces)
        i = ui.visibleTraces(vi);
        sel = find(logical(ui.traces(i).BrushData));
        if ~isempty(sel)
            [x, y] = getTrace(data.traces(i));
            y = y - mean(y(sel));
            bs = std(y(sel));
            if strcmp(direction, 'Absolute Value')
                idx = find(abs(y(sel(end) + 1:end)) > xSD * bs, 1);
            elseif strcmp(direction, 'Positive')
                idx = find(y(sel(end) + 1:end) > xSD * bs, 1);
            elseif strcmp(direction, 'Negative')
                idx = find(y(sel(end) + 1:end) < -xSD * bs, 1);
            end
            if idx
                data.traces(i).x0 = data.traces(i).x(sel(end) + idx);
            end
        end
    end
    redraw();
end

%% helper functions
function idx = str2idx(str)
    idx = [];
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
                        idx = [idx, j];
                    end
                end
            else
                idx = [idx, str2num(field)];
            end
        end
    end
    idx = unique(idx);
end

function str = idx2str(idx)
    idx = unique(idx);
    strs = [];
    a = [];
    for i = 1:numel(idx)
        if isempty(a) || idx(i) == a(end) + 1
            a = [a, idx(i)];
        else
            if numel(a) == 1
                strs = [strs, num2str(a(1))];
            else
                strs = [strs, num2str(a(1)) + "-" + num2str(a(end))];
            end
        end
    end
    str = join(strs, ",");
end

function xylims = axesXYLims(ax)
    xylims = [-inf, inf, -inf, inf];
    for i = 1:numel(ax.Children)
        if strcmp(ax.Children(i).Type, "line")
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

function propagateSelectedTracesToAllOtherGroups(groupid)
    groupids = vertcat(data.traces.groupid);
    gidx = find(groupids == groupid);
    sdx = find(ui.selectedTraces(gidx));
    ngroups = numel(unique(groupids));
    for i = 1:ngroups
        if i ~= groupid
            gidx = find(groupids == i);
            ui.selectedTraces(gidx) = false;
            idx = sdx;
            idx(find(idx > numel(gidx))) = [];
            ui.selectedTraces(gidx(idx)) = true;
        end
    end
end

%% keyboard input
function keyPress(src, event)
    key = event.Key; %disp(key);
	if strcmp(key, 'leftarrow') disp('L');
	elseif strcmp(key, 'rightarrow') disp('R');
	elseif strcmp(key, 'uparrow') disp('U');
	elseif strcmp(key, 'downarrow') disp('D');
    end
end

end
