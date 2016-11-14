function [H, G] = trophicgroupgraph(G, tg)
%TROPHICGROUPGRAPH Create graph object resulting from trophic grouping
%
% [H, G] = trophicgroupgraph(G, tg)
%
% Creates a graph object with the same nodes as G plus pseudo-nodes
% representing trophic groups.  Edges connect nodes with the appropriate
% pseudo-nodes.  See gauzensgroup.m or trophicgroup.m for details of
% grouping algorithms.
%
% Input variables:
%
%   G:  Ecopath graph object
%
%   tg: nnode x nlev of indices.  Nodes that share an index are grouped
%       together.  Multiple columns indicate multiple levels of grouping.
%       For example, for a 5-node graph
%       tg = [...
%             1   1
%             1   1
%             2   2 
%             2   2 
%             3   1];
%       indicates 3 node clusters, and then two clusters-of-clusters (i.e.
%       clusters 1 and 3 from the first level are clustered at the second
%       level).  An example of this would be multiple columns from the
%       S.tgmetricDetails.partitions output of the trophicgroup.m
%       hierarchical algorithm.
%       
% Output variables:
%
%   H:  graph object. Nodes are the same as those in the input graph (but
%       stripped of Nodes fields except Name), with the addition of several
%       grouping nodes.  Edges indicate the trophic-group connections,
%       resulting in a dendrogram-like graph.
%
%   G:  graph object. Nodes are the same as H, but this time with all Node
%       properties intact and two additional ones added: 'parent' lists the
%       name of the parent node in the trophic grouping hierarchy, and 'TG'
%       lists the trophic group index for each node (at the lowest grouping
%       level).  Edges are the same as in the input graph object. 

% Copyright 2015-2016 Kelly Kearney

%--------------------------
% Trophic group nodes
%--------------------------

% Check grouping indices

nlev = size(tg,2);
for ii = 1:nlev-1
    [~,chk] = aggregate(tg(:,ii), tg(:,ii+1), @(x) length(unique(x)));
    chk = cat(1, chk{:});
    if any(chk > 1)
        error('Nodes grouped at one level must stay together at the next');
    end
end

% Create dendrogram-like graph

nnode = numnodes(G);
tmp = [(1:nnode)' tg];

cidx = ones(nnode,1)*(1:nlev);
nodename = [G.Nodes.Name ...
    arrayfun(@(idx,lev) sprintf('lev%02d_%02d', lev, idx), tg, cidx, 'uni', 0), ...
    repmat({'web'}, nnode, 1)]; 

src = nodename(:,1:end-1);
tar = nodename(:,2:end);
src = src(:);
tar = tar(:);

unqnode = unique([src tar]);
[~,snum] = ismember(src, unqnode);
[~,tnum] = ismember(tar, unqnode);
[~,ia] = unique([snum tnum], 'rows');

Htmp = digraph(src(ia), tar(ia), ones(size(ia)));

% Get rid of intermediate nodes, and connect top-level nodes to an overall
% parent node (web)

nin = indegree(Htmp);
newsrc = Htmp.Nodes.Name(nin ~= 1);
newtar = cell(size(newsrc));

isweb = strcmp(newsrc, 'web');
for ii = 1:length(newsrc)
    if ~isweb(ii)
        pth = shortestpath(Htmp, newsrc{ii}, 'web');
        idx = find(ismember(pth, newsrc), 2);
        newtar{ii} = pth{idx(2)};
    end
end
newsrc = newsrc(~isweb);
newtar = newtar(~isweb);

H = digraph(newsrc, newtar, ones(size(newsrc)));

% Add new nodes to the original graph

newnode = setdiff(H.Nodes.Name, G.Nodes.Name);

wstat = warning('off', 'MATLAB:table:RowsAddedExistingVars');
G = addnode(G, newnode);
warning(wstat);
isn = ismember(G.Nodes.Name, newnode);
G.Nodes.B(isn) = 0;
G.Nodes.type(isn) = 5;
G.Nodes.TL(isn) = NaN;

widx = findnode(G, 'web');
G.Nodes.type(widx) = 6;

% Names for nodes and their parents

[tf, loc] = ismember(G.Nodes.Name, newsrc);
parent = cell(numnodes(G),1);
[parent{~tf}] = deal('null');
parent(tf) = newtar(loc(tf));

G.Nodes.parent = parent;

% Add TG field

G.Nodes.TG = zeros(numnodes(G),1);
G.Nodes.TG(1:nnode) = tg(:,1);





