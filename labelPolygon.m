function [textHandle, overflowFlag] = labelPolygon(inputText,polygonNodes,varargin)
%labelPolygon Takes a piece of text and a polygon and adds the text label
%WITHIN the polygon. It allows the user to dictate the maximum font size,
%the minimum font size, and the maximum number of line breaks. If you
%change the window size or axis limits after adding the text, the font size
%may be wrong. Best to plot all your shapes first and then label them. 
    p = inputParser;
    defMaxFont = 30;
    defMinFont = 4;
    defMaxLineBreaks = Inf;
    defFontName = 'Monospaced';
    
    addParameter(p, 'Numerate', nan, @isnumeric);
    addParameter(p, 'FontName', defFontName, @ischar);
    addParameter(p, 'LineBreaks', defMaxLineBreaks, @isnumeric);
    addParameter(p, 'MaxFontSize', defMaxFont, @isnumeric);
    addParameter(p, 'MinFontSize', defMinFont, @isnumeric);
    addParameter(p, 'FontColor', [0 0 0], @isnumeric);
    addParameter(p, 'TextOverflow', true, @islogical);
    
    parse(p, varargin{:});
    
    numField = p.Results.Numerate;
    maxFont = p.Results.MaxFontSize;
    minFont = p.Results.MinFontSize;
    maxLineBreaks = p.Results.LineBreaks;
    fontName = p.Results.FontName;
    txtColor = p.Results.FontColor;
    overBool = p.Results.TextOverflow;
    
    
    % size of polygon
    [~, m] = size(polygonNodes);
    if m ~= 2
        error('polygonNodes should be an n x 2 array of points (x, y)');
    end

    % can't have more line breaks than there are spaces
    maxLineBreaks = min([maxLineBreaks, sum(isspace(inputText))]);
    
    shp = polygonNodes;
    [cX, cY] = centroid(polyshape(unique(shp, 'rows', 'stable'),...
                        'KeepCollinearPoints', true));
                    
    overflowFlag = false;
    
    % draw the text at max size
    t = text(cX, cY, inputText, 'FontSize', maxFont, 'FontName',...
                    fontName, 'HorizontalAlignment', 'center',...
                    'Color', txtColor);
    
    % extent of the text
    tmp = t.Extent;
    bdX = [tmp(1), tmp(1), tmp(1) + tmp(3), tmp(1) + tmp(3)];
    bdY = [tmp(2), tmp(2) + tmp(4), tmp(2) + tmp(4), tmp(2)];
    
    % text variants across multiple lines
    spTexts = cell(maxLineBreaks+1, 1);
    spTexts{1} = inputText;
    for i = 1:maxLineBreaks
        spTexts{i+1} = multiline(inputText, i);
    end
    
    fSize = maxFont;
    % slowly decrease font size trying monoline and multiline variants of
    % the label until it fits.
    while fSize >= minFont
        t.FontSize = fSize;
        for i = 1:length(spTexts)
            t.String = spTexts{i};

            tmp = t.Extent;

            bdX = [tmp(1), tmp(1), tmp(1) + tmp(3), tmp(1) + tmp(3)];
            bdY = [tmp(2), tmp(2) + tmp(4), tmp(2) + tmp(4), tmp(2)];
        
            if inpolygon(bdX, bdY, shp(:,1), shp(:,2))
                break;
            else
                fSize = fSize - 1;
            end
        end
        if inpolygon(bdX, bdY, shp(:,1), shp(:,2))
                break;
        end
    end
    
    % if it never fits then flag it and optionally add a number to it. 
    if fSize < minFont
        delete(t);
        warning(['text ''', inputText, ''' unable to fit inside polygon.',....
            'Consider decreasing MinFontSize or shortening the label']);
        overflowFlag = true;
        if overBool
            t = text(cX, cY, multiline(inputText, 1), 'FontSize', minFont, 'FontName',...
                    fontName, 'HorizontalAlignment', 'center',...
                    'Color', txtColor);
        end
        if ~isnan(numField)
             numSize = (maxFont + minFont) / 2;
             t = text(cX, cY, num2str(numField), 'FontSize', numSize, 'FontName',...
                    fontName, 'HorizontalAlignment', 'center',...
                    'Color', txtColor);
        end
    end
    textHandle = t;
    
    % function for generating multi-line versions of the label 
    function cString = multiline(txt, lines)
        nls = linspace(1, length(txt), lines + 2);
        nls = nls(2:end-1);
        spInds = find(isspace(txt));
        diff = abs(nls - spInds');
        [a, ~] = find(diff == min(diff));
        repInd = spInds(a);
        txt(repInd) = newline;
        cString = txt;
    end
    
end

