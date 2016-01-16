function gen37_foodweb

%% Step 1: Start with an Ecopath model

Ewein = mdb2ewein('Generic_37.EwEmdb');

%% Step 2: Convert to a graph object

G = ecopath2graph(Ewein);

%% Step 3: Prepare graph object for trophic grouping

G2 = G;

detid = G2.Nodes.Name(G2.Nodes.type == 2);
oosid = G2.Nodes.Name(G2.Nodes.type == 4);

G2 = rmnode(G2, oosid);
isdet = ismember(G2.Edges.EndNodes(:,2), detid);
G2 = rmedge(G2, find(isdet));

%% Step 4: Calculate trophic groups

grp = gauzensgroupep(G2, 1, 'type', 'trophicgroup');

%% Step 5: Position nodes

[H3, G3] = trophicgroupgraph(G2, grp);

opt = {'drawmode', 'force', ...
    'init', 'cpack', ...
    'cpfac', 0.6, ...
    'totheight', 600, ...
    'totwidth', 500, ...
    'fontsz', 8, ...
    'charge', -100, ...
    'padding', 10, ...
    'tlnudge', 2.0, ...
    'tllim', [0.5 5]};

% [hw, stat] = renderfoodwebforce(G3, 'mode', 'test', opt{:});

[C,T,F] = renderfoodwebforce(G3, 'mode', 'extract', opt{:});

[G2, Ax] = foodweblayoutdetails(G2, C, T, [0.5 5]);


%% Step 8: Run edge bundling

wlim = minmax(log10(G2.Edges.Weight), 'expand', 0.01);
wtfun = @(x) log10(x) - wlim(1);
G2 = debundle(G2, 'edgefun', wtfun, 'l', 100);

%% Step 9: Plot

mar = 50;
h.fig = figure('position', [10 10 Ax.dx+mar*2 Ax.dy+mar*2], 'color', 'w');
h.ax = axes('units', 'pixels', ...
            'position', [mar mar Ax.dx Ax.dy], ...
            'xlim', Ax.xlim, ...
            'ylim', Ax.ylim, ...
            'xcolor', 'none');
hold on;



h.edg = plotdeb(G2, 'p', 1/3, 'w', 50, 'rloop', 50);
set(h.edg, 'clipping', 'off');

th = linspace(0,2*pi,50);
h.nd = arrayfun(@(x,y,r) patch(r.*cos(th)+x, r.*sin(th)+y, 'w'), ...
    G2.Nodes.x, G2.Nodes.y, G2.Nodes.r, 'uni', 0);
h.nd = cat(1, h.nd{:});
set(h.nd, 'facecolor', ones(1,3)*0.8, 'clipping', 'off', 'edgecolor' ,'none');

h.txt = text(G2.Nodes.tx, G2.Nodes.ty, G2.Nodes.Name, 'fontsize', 8, 'fontname', 'Arial');
set(h.txt, {'horizontalalignment', 'VerticalAlignment'}, [G2.Nodes.th G2.Nodes.tv]);

ylabel('Trophic Level');
colormap(usercolormap([1 0 0], [0 0 1]));

% export_fig('foodwebcode', gcf, '-png', '-r150', '-nocrop');
