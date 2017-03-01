function [G, Ax] = fwsvgdetails(G, C, T)
%FOODWEBLAYOUTDETAILS Convert details from foodwebforce.js SVG to plot
%
% [G2, Ax] = fwsvgdetails(G, C, T)
% [G2, Ax] = fwsvgdetails(G, C)
%
% This function applies the circle/text info from a d3.foodweblayout-
% generated SVG to the Node properties of a digraph object (assuming that
% the nodes correspond to those used to call d3.foodweblayout).  It also
% returns information regarding the size and position of the SVG canvas
% that can be used to recreate the axis in a Matlab figure.
%
% Input variables:
%
%   G:      food web digraph object to which details will be added. Nodes
%           table should include a Name property, and the names should
%           correspond to the strings in the label field of the C struture
%           input (though order does not need to match).
%
%   C:      structure with details extracted from SVG circle objects (see
%           extractsvg.m).
%
%   T:      structure with details extracted from SVG text objects (see
%           extractsvg.m).  This is optional, and if not included (or if
%           all fields are empty), text labels will be assigned centered on
%           each node, using the font and fontsize specified by the
%           DefaultAxesFontName and DefaultAxesFontSize, respectively, for
%           the Matlab session. 
%
% Output variables:
%
%   G2:     digraph object, a copy of input G with the following properties
%           added to the Nodes table: 
%           y:          y-coordinate of node center, in pixel units
%           x:          x-coordinate of node center, in pixel units
%           r:          radius of node, in pixel units
%           tx:         x-coordinate of node label, in pixel units
%           ty:         y-coordinate of node label, in pixel units
%           th:         horizontal alignment of node label
%           tv:         vertical alignment of node label
%
%   Ax:     structure holding information about the axis that are required
%           to replicate the svg canvas
%           ylim:       y limits (scale: 1 unit = 1 pixel).  Note that this
%                       scale assumes a reverse y-orientation like those
%                       used by images. 
%           xlim:       x limits (scale: 1 unit = 1 pixel)
%           figpos:     position vector for figure (in pixels)
%           axpos:      position vector for axis (normalized)
%           fontsize:   fontsize for text labels
%           fontname:   fontname for text labels

% Copyright 2016 Kelly Kearney

% Parse inputs


% Assign circle details to appropriate nodes in the graph object

G.Nodes.x  = vlookup(C, G.Nodes.Name, 'x', 'label');
G.Nodes.y  = vlookup(C, G.Nodes.Name, 'y', 'label');
G.Nodes.r  = vlookup(C, G.Nodes.Name, 'r', 'label');

% Add details about text, if relevant

tempty = nargin < 3 || isempty(T.x);

if tempty
    
    % Defaults: center of nodes
    
    G.Nodes.tx = G.Nodes.x;
    G.Nodes.ty = G.Nodes.y;
    G.Nodes.th = repmat({'center'}, numnodes(G), 1);
    G.Nodes.tv = repmat({'middle'}, numnodes(G), 1);
else
    
    % Locations based on SVG text objects
    
    G.Nodes.tx = vlookup(T, G.Nodes.Name, 'x', 'label');
    G.Nodes.ty = vlookup(T, G.Nodes.Name, 'y', 'label');

    G.Nodes.th = vlookup(T, G.Nodes.Name, 'anchor', 'label');
    G.Nodes.th = regexprep(G.Nodes.th, {'start', 'middle'}, {'left', 'center'});

    [~, loc] = ismember(G.Nodes.th, {'left', 'center'});
    vert = {'base'; 'middle'};

    G.Nodes.tv = vert(loc);
    
end

% Extract details of figure and axis

Ax.ylim = [C.axpos(2) C.axpos(2)+C.axpos(4)]; % pixel units
Ax.xlim = [C.axpos(1) C.axpos(1)+C.axpos(3)]; % pixel units

Ax.figpos = [0 0 C.svgsz]; % pixels

w = C.axpos(3);
h = C.axpos(4);
Ax.axpos = [C.axpos(1) C.svgsz(2)-C.axpos(2)-h w h]./[C.svgsz C.svgsz]; % normalized

if tempty
    Ax.fontsize = get(0, 'DefaultAxesFontSize');
    Ax.fontname = get(0, 'DefaultAxesFontName');
else
    Ax.fontsize = T.fontsize(1);
    Ax.fontname = T.font{1};
end






