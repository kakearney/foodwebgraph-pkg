function h = plotfoodweb(G, Ax, varargin)
%PLOTFOODWEB Plots a trophic-group/bundled-edge style food web
%
% h = plotfoodweb(G, Ax)
%
% Plots a food web, using plotdeb.m for edges and circular patches for
% nodes.  This is intended to plot the final product of the
% trophic-group-trophic-level layout algorithm (i.e. foodweblayout.m)
% combined with edge bundling.
%
% Input variables:
%
%   G:          digraph object, with node and text label positions included
%               in the Node property table properties (see foodweblayout)  
%
%   Ax:         structure holding axis info (see foodweblayout)
%
% Optional input variables (passed as parameter/value pairs):
%
%   nodecolor:  1 x 3 RGB color vector, color for node faces [0.8 0.8 0.8]
%
%   edgecolor:  method used to color edges:
%               gradient:   edge color varies from source to target, same
%                           for each edge (expects axis color limits to
%                           range from 0-1)
%               byEdge:     color reflects edge number (best with colormap
%                           with nedge discreet colors and color limits 
%                           [0 nedge]+0.5)
%               bySource:   color reflect source node number (best with
%                           colormap with nnode discreet colors and color
%                           limits [0 nnode]+0.5
%               byTarget:   color reflect target node number (best with
%                           colormap with nnode discreet colors and color
%                           limits [0 nnode]+0.5
%               [gradient]
%
%   Also accepts all parameters applicable to plotdeb.m
%
% Output variables:
%
%   h:          scalar structure holding graphics handles:
%
%               fig:    1 x 1, handle to figure   
%
%               ax:     1 x 1, handle to axis
%
%               edg:    1 x 1, handle to edges patch.  Note that all edges
%                       are one multifaceted patch; each row of the Faces
%                       property and each column of the XData/YData/CData
%                       properties correspond to a single edge.
%
%               nd:     nnode x 1, handles to node patches
%
%               txt:    nnode x 1, handles to node text labels

% Copyright 2016 Kelly Kearney

% Parse input

p = inputParser;
p.KeepUnmatched = true;
p.addParameter('nodecolor', ones(1,3)*0.8,    @(x) validateattributes(x, {'numeric'}, {'size', [1 3], '>=', 0, '<=', 1}));
p.addParameter('edgecolor', 'gradient', @(x) validateattributes(x, {'char'},{}));

p.parse(varargin{:});

Opt = p.Results;

Opt.edgecolor = validatestring(Opt.edgecolor, {'gradient', 'byEdge', 'bySource', 'byTarget'});

% Assume unmatched options apply to the plotting of edges, and let
% plotdeb.m throw an error if not

debOpt = p.Unmatched;

% Set up axis

h.fig = figure('position', Ax.figpos, 'color', 'w');
h.ax = axes('position', Ax.axpos, ...
            'xlim', Ax.xlim, ...
            'ylim', Ax.ylim, ...
            'xcolor', 'none');
hold on;

% Plot edges, and alter color if necessary

h.edg = plotdeb(G, debOpt);
set(h.edg, 'clipping', 'off');


nedge = numedges(G);
nvert = size(h.edg.CData,1);
switch Opt.edgecolor
    case 'byEdge'
        h.edg.CData = ones(nvert,1)*(1:nedge);
    case 'bySource'
        h.edg.CData = ones(nvert,1)*findnode(G, G.Edges.EndNodes(:,1))';
    case 'byTarget'
        h.edg.CData = ones(nvert,1)*findnode(G, G.Edges.EndNodes(:,2))';
end


% Plot nodes

th = linspace(0,2*pi,50);
h.nd = arrayfun(@(x,y,r) patch(r.*cos(th)+x, r.*sin(th)+y, 'w'), ...
    G.Nodes.x, G.Nodes.y, G.Nodes.r, 'uni', 0);
h.nd = cat(1, h.nd{:});
set(h.nd, 'facecolor', ones(1,3)*0.8, 'clipping', 'off', 'edgecolor' ,'none');

% Plot label

h.txt = text(G.Nodes.tx, G.Nodes.ty, G.Nodes.Name, 'fontsize', Ax.fontsize, 'fontname', Ax.fontname, 'interpreter', 'none');
set(h.txt, {'horizontalalignment', 'VerticalAlignment'}, [G.Nodes.th G.Nodes.tv]);

% Label axis

ylabel('Trophic Level');

