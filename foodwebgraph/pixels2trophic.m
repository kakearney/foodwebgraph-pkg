function [G, Ax] = pixels2trophic(G, Ax, tllim)
%PIXELS2TROPHIC Convert food web graph to trophic level units
%
% [G2, Ax2] = pixels2trophic(G, Ax, tllim)
% [G2, Ax2] = pixels2trophic(G, Ax)
%
% The d3-foodweb functions position nodes based on pixels, with an internal
% calculation for the trophic-level-to-pixel positioning needed for the y
% axis.  When plotting in Matlab it's easier to recast everything (x- and
% y-coordinates as well as node radius) in terms of these trophic level
% units.  This allows us to maintain a 1:1 axis ratio for easy circle
% plotting while keeping the y-axis in intuitive trophic level units.
%
% Input variables:
%
%   G:      food web digraph object.  The Nodes table should include node
%           position, size, and text labels details, as returned by the
%           foodweblayout.m function.  The Edges table can include
%           coordinates for the edge pathways (x and y variables), but
%           these are optional.
%
%   Ax:     structure including figure and axis size and limit details, as
%           returned by foodweblayout or fwsvgdetails.
%
%   tllim:  trophic level limits used to calculate the positions of the
%           nodes in the graph. If not included, will assume the default
%           (extent of TL values in the graph object G).   
%
% Output variables:
%
%   G2:     food web digraph object, identical to input G except the x, y,
%           tx, ty, and r Nodes properties, and x and y Edge properties
%           (if present) are converted from image-style pixel coordinates
%           to equivalent trophic level coordinates.  This includes a flip
%           in direction for the y-coordinates to match Matlab's 'normal'
%           axis direction convention.
%
%   Ax2:    structure identical to the input Ax but with the xlim and ylim
%           fields modified to match the new coordinate system.

% Copyright 2017 Kelly Kearney

if nargin < 3
    tllim = [min([G.Nodes.TL]) max([G.Nodes.TL])];
end

% Convert pixel coordinates to trophic level coordinates

axpospix = Ax.axpos .* [Ax.figpos(3:4) Ax.figpos(3:4)];

m = (tllim(1) - tllim(2))./axpospix(4);
mart = Ax.figpos(4) - axpospix(4) - axpospix(2);
fun = @(x) m*(x - (axpospix(4)+mart)) + tllim(1);

fac = axpospix(4)./diff(tllim);

% Change node and text coordinates to trophic level units

G.Nodes.x = G.Nodes.x./fac;
G.Nodes.tx = G.Nodes.tx./fac;
G.Nodes.r = G.Nodes.r./fac;

G.Nodes.y = fun(G.Nodes.y);
G.Nodes.ty = fun(G.Nodes.ty);

% If edge path coordinates exist, convert these too

haspath = all(ismember({'x','y'}, G.Edges.Properties.VariableNames));
if haspath
    G.Edges.x = cellfun(@(x) x./fac, G.Edges.x, 'uni', 0);
    G.Edges.y = cellfun(fun, G.Edges.y, 'uni', 0);
end

% Convert x- and y-axis limits

Ax.ylim = fliplr(fun(Ax.ylim));
Ax.xlim = Ax.xlim./fac;
