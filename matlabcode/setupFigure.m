function h = setupFigure(ifig, figsize, varargin)
%%
% h = setupFigure(ifig,figsize) creates a figure and set some standard properties
% figsize: vector [ posx posy width height]
% ifig: number of vector
% h: handle of figure
%

if nargin > 2
    Visible = varargin{1};
else
    Visible = 'on';
end

if nargin > 3
    units = varargin{2};
else
    units = 'pixel';
end

if nargin >4
    fontSize = varargin{3};
else
    fontSize = 12;
end
if nargin > 5
    box = varargin{4};
else
    box = 'off';
end

try
    if ifig> 0
        h = figure(ifig);
        set(gcf, 'Visible', Visible);
    else
        h = figure('Visible', Visible);
        set(gcf, 'Visible', Visible);
    end
catch
    if ifig> 0
        h = figure(ifig);
        set(gcf, 'HandleVisibility', Visible);
    else
        h = figure('HandleVisibility', Visible);
        set(gcf, 'Visible', Visible);
    end
end

clf
hold('on');
%set(h,'OuterPosition',figsize);
set(gcf,'Position',figsize);
try
set(gcf,'ActivePositionProperty','Position')
end
set(gcf, 'units', units, 'pos', figsize)
%set(gcf, 'Color', [1 1 1]); % Sets figure background
%set(gca, 'Color', [1 1 1]); % Sets axes background
%set(gcf, 'Renderer', 'painters');
set(gcf, 'PaperPositionMode', 'auto');
set(gca,'Box',box,'FontSize',fontSize);
end
