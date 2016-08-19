function h = trophicgroupdendro(links, cutoff, energy, names)
%TROPHICGROUPDENDRO Interactive dendrogram of hierarchical trophic groups
%
% trophicgroupdendro(links, cutoff, energy, names)
% 
% Plot interactive dendrogram of hierachical trophic grouping results (see
% trophicgroup.m).
%
% Input variables:
%
%   links:  cluster tree links, output of linkage function
%
%   cutoff: cutoff values for color threshhold
%
%   energy: trophic group index value at each cutoff value
%
%   names:  food web node names
%
% Output variables:
%
%   h:      structure holding graphics handles.

% Copyright 2016 Kelly Kearney

elim = [0 max(energy)*1.1];

plotdendro = @(c) dendrogram(links, 0, 'orientation', 'right', ...
    'colorthreshold', c, 'labels', names);

% h = plotgrid('setup', cell(1), [],[],'mar', 0.05, 'ml', 0.2);
% h.ax = subgrid(h.ax, [0.1 0.9], 1);
% h.ax = cat(1, h.ax{:});

h.fig = figure;
h.ax = gobjects(2,1);
h.ax(1) = axes('position', [0.2 0.14 0.75 0.81]);
h.ax(2) = axes('position', [0.2 0.05 0.75 0.09]);

axes(h.ax(2));
h.stm = stem(cutoff, energy, '.');
h.ref = line(nan(1,2), elim, 'color', 'r');

h.stm.ButtonDownFcn = {@stemclick, h.ref, plotdendro, h.ax(1)}

axes(h.ax(1));
plotdendro(Inf);

xlim = minmax(h.ax(1).XLim, 'expand', 0.02);

set(h.ax, 'xlim', xlim);

set(h.ax(1), 'xaxisloc', 'top', 'box', 'off', 'TickLabelInterpreter', 'none');
set(h.ax(2), 'box', 'off');

% Subfunction:

function stemclick(h,ed, href, dfun, dax)

xyz = get(ancestor(h, 'axes'), 'CurrentPoint');
[dx, imin] = min(abs(xyz(1) - h.XData));

thresh = h.XData(imin);

set(href, 'xdata', ones(1,2)*thresh);

props = {'xlim' 'xaxisloc', 'box', 'TickLabelInterpreter'};
vals = get(dax, props);
pv = [props; vals];

cla(dax);
axes(dax);
dfun(thresh);

set(dax, pv{:});
title(dax, sprintf('Cutoff = %.2f (index %d)', thresh, imin));