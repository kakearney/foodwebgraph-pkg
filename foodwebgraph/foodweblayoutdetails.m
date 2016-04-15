function [G, Ax] = foodweblayoutdetails(G, C, T, tllim)
%FOODWEBLAYOUTDETAILS Convert details from foodwebforce.js SVG to plot
%
% [G, Ax] = foodweblayoutdetails(G, C, T)
% [G, Ax] = foodweblayoutdetails(G, C, T, tllim)
%
% This function converts the pixel-based details returned by
% renderfoodwebforce.m in 'extract' mode to a more plotting-friendly
% coordinate system.
%
% Input variables:
%
%   G:      graph object to which details will be added.  In the usual
%           foodweblayout.js Matlab workflow, this will be the object
%           passed *into* trophicgroupgraph.m, as opposed to the one passed
%           to renderfoodwebforce.m; however, any graph object holding
%           Nodes with names corresponding to those in the C and T
%           structures can be used.
%
%   C:      structure with circle element data.  See renderfoodwebforce.m
%           for details.
%
%   T:      structure with text element data.  See renderfoodwebforce.m
%           for details.
%
%   tllim:  trophic level limits used to generate the C and T data
%           structures.  If not included, will assume the default (extent
%           of TL values in the graph object G).
%
% Output variables:
%
%   G:      graph object, same as input with the following parameters added
%           to the Nodes table:
%           y:  y-coordinate of node center, in trophic level units
%           x:  x-coordinate of node center, in trophic level units.  The
%               units in the x-direction are arbitrary, and are kept the
%               same as the y-units to make plotting of circles simpler.
%           r:  radius of node, in trophic level units
%           tx: x-coordinate of node label, in trophic level units
%           ty: y-coordinate of node label, in trophic level units
%           th: horizontal alignment of node label
%           tv: vertical alignment of node label
%
%   Ax:     structure holding information about the axis
%           ylim:   y limits
%           xlim:   x limits
%           dy:     height of axis, in pixels
%           dx:     width of axis, in pixels

% Copyright 2016 Kelly Kearney


% Trophic level limits assigned to reference points

if nargin < 4
    tllim = [min([G.Nodes.TL]) max([G.Nodes.TL])];
end

% Convert pixel coordinates to trophic level coordinates

m = (tllim(1) - tllim(2))./(C.yref(2)-C.yref(1));
fun = @(x) m*(x - C.yref(2)) + tllim(1);

fac = diff(C.yref)./diff(tllim);

if isempty(T.x)
    error('No text object found; did you forget to label?');
end

C.x = C.x./fac;
C.y = fun(C.y);
C.r = C.r./fac;
T.x = T.x./fac;
T.y = fun(T.y);

% Assign details to appropriate nodes in the graph object

G.Nodes.x  = vlookup(C, G.Nodes.Name, 'x', 'label');
G.Nodes.y  = vlookup(C, G.Nodes.Name, 'y', 'label');
G.Nodes.r  = vlookup(C, G.Nodes.Name, 'r', 'label');
G.Nodes.tx = vlookup(T, G.Nodes.Name, 'x', 'label');
G.Nodes.ty = vlookup(T, G.Nodes.Name, 'y', 'label');

G.Nodes.th = vlookup(T, G.Nodes.Name, 'anchor', 'label');
G.Nodes.th = regexprep(G.Nodes.th, {'start', 'middle'}, {'left', 'center'});

[tf, loc] = ismember(G.Nodes.th, {'left', 'center'});
vert = {'base'; 'middle'};

G.Nodes.tv = vert(loc);

% Calculate figure/axis size necessary to get 

Ax.ylim = fun(C.yref([2 1]));
Ax.xlim = C.xref./fac;
Ax.figpos = [0 0 C.svgsz];

w = diff(C.xref);
h = diff(C.yref);
Ax.axpos = [C.trans(1) C.svgsz(2)-C.trans(1)-h w h]./[C.svgsz C.svgsz];


% Ax.dx = diff(C.xref);
% Ax.dy = diff(C.yref);




