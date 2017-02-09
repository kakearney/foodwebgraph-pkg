function gen37_foodweb
%GEN37_FOODWEB An example script for the food web layout/edge bundling

%% Step 1: Start with an Ecopath model

pth = '~/Documents/Research/Working/EcopathModels/';
EM = mdb2ecopathmodel(fullfile(pth, 'Generic_37.EwEmdb'));

%% Step 2: Calculate trophic groups

[grp, S, G] = trophicgroup(EM);
% [grp, S, G] = trophicgroup(EM, 'method', 'simulatedannealing', 'verbose', true);

%% Step 3: Position nodes with interactive viewer

% First, a small cosmetic adjustment

G.Nodes.B = min(G.Nodes.B, 25); % detrital node overshadows, so shrink it

% Now run the layout tool

[G, Ax] = foodweblayout(G, grp);

%% Step 4: Edge bundling

wlim = minmax(log10(G.Edges.Weight), 'expand', 0.01);
wtfun = @(x) log10(x) - wlim(1);
G = debundle(G, 'edgefun', wtfun, 'l', 100);

%% Step 5: Plot

% The default: solid-colored nodes and gradient edges

h = plotfoodweb(G, Ax, 'p', 1/3, 'w', 50, 'rloop', 60, 'cthresh', 1);

h.cb = colorbar('south');
set(h.cb, 'position', [h.ax.Position(1:2)+[0.03 0] h.cb.Position(3:4).*[0.25 1]]);
set(h.cb, 'ticks', [0 1], 'ticklabels', {'Source','Target'});

% Color by source or target node

nnode = numnodes(G);
cmap = jet(nnode);

h = plotfoodweb(G, Ax, 'p', 1/3, 'w', 50, 'rloop', 60, 'cthresh', 1, ...
    'edgecolor', 'bySource');
colormap(cmap);
set(gca, 'clim', [0 nnode] + 0.5);
set(h.nd, {'facecolor'}, num2cell(cmap,2));

h = plotfoodweb(G, Ax, 'p', 1/3, 'w', 50, 'rloop', 60, 'cthresh', 1, ...
    'edgecolor', 'byTarget');
colormap(cmap);
set(gca, 'clim', [0 nnode] + 0.5);
set(h.nd, {'facecolor'}, num2cell(cmap,2));

% Color by detritus fraction

detfrac = sourcefraction(G, 'Detritus');

h = plotfoodweb(G, Ax, 'p', 1/3, 'w', 50, 'rloop', 60, 'cthresh', 1, ...
    'edgecolor', 'bySource');
h.edg.CData = detfrac(h.edg.CData);
set(h.nd, {'cdata'}, num2cell(detfrac));
set(h.nd, 'facecolor', 'flat');

h.cb = colorbar('south');
set(h.cb, 'position', [h.ax.Position(1:2)+[0.03 0] h.cb.Position(3:4).*[0.25 1]]);
title(h.cb, 'Fraction from detritus');
set(h.txt(findnode(G, 'Phytoplankton')), 'color', 'w');








