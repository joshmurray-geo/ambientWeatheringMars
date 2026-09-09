function [figHandle, patchHandles, textHandles, fieldLabels] =...
                                    plumx(minModesTable, varargin)
%perpetrate Takes a mineral modes table and produces a phase diagram. The
% Name-Value arguments can optionally specify: the column to serve as the x
% axis; the column to serve as the y axis; an array of colours for each
% field; cell array of preferred name for each column; value for the
% smoothing of the boundary; value for the factor by which each field is
% expanded (as a fraction of the increment in the x and y direction; the
% width of the outline of the phase fields; the colour of the outline of
% the phase fields; the order in which to draw the fields (specify 'ascend'
% or 'descend' so they draw from smallest to largest and vice versa.
% Additional optional parameters for the labelling are: the font; the
% maximum font size; the minimum font size; the maximum number of line
% breaks in the labels; the colour of the text.

    % ------ parse all the inputs
    p = inputParser;
    
    % replace underscores for the sake of the MATLAB text parser. use the
    % 'MineralLabels' name-value pair if you want to specify the labels
    colLabels = minModesTable.Properties.VariableNames;
    colLabels = strrep(colLabels, '_', '.');
    
    % patch parameters
    addParameter(p, 'XAxisColumn', 1, @isnumeric);
    addParameter(p, 'YAxisColumn', 2, @isnumeric);
    addParameter(p, 'PhaseColors', nan, @isnumeric);
    addParameter(p, 'MineralLabels', colLabels);
    addParameter(p, 'SmoothValue', 0.2, @isnumeric);
    addParameter(p, 'ExpansionValue', 0.6, @isnumeric);
    addParameter(p, 'LineWidth', 3, @isnumeric);
    addParameter(p, 'LineColor', [1 1 1], @isnumeric);
    addParameter(p, 'FieldOrder', 'descend', @ischar);
    
    % label parameters
    addParameter(p, 'Numerate', false, @islogical);
    addParameter(p, 'FontName', 'Monospaced', @ischar);
    addParameter(p, 'MaxFontSize', 30, @isnumeric);
    addParameter(p, 'MinFontSize', 4, @isnumeric);
    addParameter(p, 'LineBreaks', Inf, @isnumeric);
    addParameter(p, 'FontColor', [0 0 0], @isnumeric);
    addParameter(p, 'TextOverflow', true, @islogical);
    
    parse(p, varargin{:});

    xInd = p.Results.XAxisColumn;
    yInd = p.Results.YAxisColumn;
    cs = p.Results.PhaseColors;
    minNames = p.Results.MineralLabels;
    smoothVal = p.Results.SmoothValue;
    expVal = p.Results.ExpansionValue;
    lw = p.Results.LineWidth;
    lc = p.Results.LineColor;
    plotOrder = p.Results.FieldOrder;
    
    numField = p.Results.Numerate;
    maxFont = p.Results.MaxFontSize;
    minFont = p.Results.MinFontSize;
    maxLineBreaks = p.Results.LineBreaks;
    fontName = p.Results.FontName;
    txtColor = p.Results.FontColor;
    overBool = p.Results.TextOverflow;
    
    %% -------- Pull out the labels and extract the x and y variable
    n = width(minModesTable);
    colNums = 1:n;
    colInds = setdiff(colNums, [xInd, yInd]);
    
    minModes = minModesTable(:,colInds);
    minNames = minNames(colInds);
    xs = minModesTable{:,xInd};
    ys = minModesTable{:,yInd};
    
    m = sqrt(length(xs));
    if m ~= floor(m)
        error(['perpetrate assumes a table of height such that it can be',...
                'drawn as an n x n square']);
    end
    
    
    % ----- convert mineral modes into assemblage index
    assIn = minModes{:,:};
    assIn = ~isnan(assIn);
    [assCols, ~, assIn] = unique(assIn, 'rows', 'stable');
    
    % ----- extract text label for each field
    fieldLabels = cell(height(assCols), 1);
    for i = 1:length(fieldLabels)
        assMins = minNames(assCols(i,:));
        for j = 1:length(assMins)
            if j == 1
                fieldLabels{i} = assMins{j};
            else
                fieldLabels{i} = [fieldLabels{i}, ' + ', assMins{j}];
            end
        end
    end
    
    % ----- locate boundaries of the perplex points
    zs = assIn;
    bds = cell(max(zs, [], 'omitnan'), 1);
    xRng = [min(xs, [], 'omitnan'), max(xs, [], 'omitnan')];
    yRng = [min(ys, [], 'omitnan'), max(ys, [], 'omitnan')];
    xInc = (xRng(2) - xRng(1)) / m;
    yInc = (yRng(2) - yRng(1)) / m;
    
    for i = 1:length(bds)
        % find all the points of that index
        ind = zs == i;
        tmpX = xs(ind);
        tmpY = ys(ind);
        
        % draw a boundary around those points
        tmpB = boundary(tmpX, tmpY, smoothVal);
        
        % edge case where field is just a line
        if isempty(tmpB)
            tmpX = repmat(tmpX, 4, 1);
            tmpY = repmat(tmpY, 4, 1);
            tmpX = tmpX + (rand(size(tmpX)) - 0.5) * xInc / 3;
            tmpY = tmpY + (rand(size(tmpX)) - 0.5) * xInc / 3;
            tmpB = boundary(tmpX, tmpY, smoothVal);
        end
        
        res = [tmpX(tmpB), tmpY(tmpB)];
        resCtr = mean(res);
        resDiff = res - resCtr;
        resDiff = resDiff ./ max(resDiff);
        res = res + [resDiff(:,1) * xInc, resDiff(:,2) * yInc] *...
                                            expVal;

        bds{i} = res;        
    end
    
    % if no colour was put in, then default to grayscale
    if isnan(cs)
        cs = gray(length(bds) + 6);
        cs = cs(6:end-1,:);
    end
    
    % sort the boundaries by area
    areas = nan(size(bds));
    for i = 1:length(bds)
        tmp = bds{i};
        areas(i) = polyarea(tmp(:,1), tmp(:,2));
    end
    [~, ind] = sort(areas, plotOrder);
    bds = bds(ind);
    fieldLabels = fieldLabels(ind);
    cs = cs(ind,:);
    
    
    % new figure. Big and square.
    f = gcf;
    set(f, 'Position', [100, 100, 1000, 1000]);
    
    % plot the fields for each mineral assemblage
    ps = gobjects(length(bds), 1);
    for i = 1:length(bds)
        tmp = bds{i};
        if ~isempty(tmp)
            ps(i) = patch(tmp(:,1), tmp(:,2), cs(i,:),...
                                            'LineStyle', 'none');
            hold on;
            plot(tmp(:,1), tmp(:,2), 'LineWidth', lw, 'Color', lc);
        end
    end
    
    % make axes square
    pbaspect([1 1 1]);
    xlim(xRng);
    ylim(yRng);
        
    figureAesthetics();
    
    ts = cell(size(bds));
    % go through each field and label it
    if numField
        numCount = 1;
    end
    % labelPolygon is a custom function for labels within fields.
    for i = 1:length(bds)
        if numField
            [t, ot] = labelPolygon(fieldLabels{i}, bds{i}, 'Numerate', numCount,...
                    'FontName', fontName,...
                    'LineBreaks', maxLineBreaks, 'MaxFontSize', maxFont,...
                    'MinFontSize', minFont, 'FontColor', txtColor,...
                    'TextOverflow', overBool);
            if ot
                numCount = numCount + 1;
            end
        else
            t = labelPolygon(fieldLabels{i}, bds{i}, 'FontName', fontName,...
                    'LineBreaks', maxLineBreaks, 'MaxFontSize', maxFont,...
                    'MinFontSize', minFont, 'FontColor', txtColor,...
                    'TextOverflow', overBool);
        end
        if ~isempty(t)
            ts{i} = t;
        end
    end
    
    % asign array to return
    figHandle = f;
    textHandles = ts;
    patchHandles = ps;
    
    function [] = figureAesthetics()
    % custom plotting spec. Alters the current figure. 
        axCol = [0 0 0];
        backCol = [1 1 1];
        ax = gca;
        fMod = gcf;
        fMod.Position = [70 70 850 850];
        set(ax, 'color', backCol);
        set(gcf, 'color', backCol);
        grid on;
        ax.GridColor = [0.95, 0.95, 0.95];
        ax.GridAlpha = 0;
        x = get(ax, 'xlabel');
        y = get(ax, 'ylabel');
        ax.FontSize = 20;
        tit = ax.Title;
        tit.FontWeight = 'normal';
        tit.FontSize = 35;
        x.FontSize = 35;
        y.FontSize = 35;
        ax.FontName = 'Helvetica';
        ax.LineWidth = 3;
        if (isprop(ax, 'Legend'))
            if (~isempty(ax.Legend))
                ax.Legend.LineWidth = 3;
            end
        end
        ax.XColor = axCol;
        ax.YColor = axCol;
        left_color = axCol;
        right_color = axCol;
        set(fMod,'defaultAxesColorOrder',[left_color; right_color]);
        box off;
        yt = yticks;
        if (length(yt) > 4)
            yticks(yt(1:2:end));
        end
        xt = xticks;
        if (length(xt) > 5)
            xticks(xt(1:2:end));
        end
        set(ax, 'TickLength', [0.03, 0.2])
        set(ax, 'YMinorTick', 'off');
        set(ax, 'XMinorTick', 'off');
        % make legend edge black
        if ~isempty(ax.Legend)
            ax.Legend.EdgeColor = [0 0 0];
        end

        box on;
        fMod.Renderer='Painters';
        %mlabel off; plabel off; gridm off;
        pbaspect([1 1 1]);
        set(ax, 'layer', 'top')
        shg;
    end
end

