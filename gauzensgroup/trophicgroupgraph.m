function [H, G] = trophicgroupgraph(G, tg)
%TROPHICGROUPGRAPH Create graph object resulting from trophic grouping
%
% [H, G] = trophicgroupgraph(G, tg)
%
% Creates a graph object with the same nodes as G plus pseudo-nodes
% representing trophic groups.  Edges connect nodes with the appropriate
% pseudo-nodes.  See gauzensgroup.m for details of trophic grouping.
%
% Input variables:
%
%   G:  Ecopath graph object
%
%   tg: nnode x 1 vector of trophic group indices corresponding to each
%       node, or n x 1 cell array of holding multiple levels of grouping
%       indices (for the latter, see gauzensgroupep.m).
%
% Output variables:
%
%   H:  graph object. Nodes are the same as those in the input graph (but
%       stripped of Nodes fields except Name), with the addition of several
%       grouping nodes.  Edges indicate the trophic-group connections 
%
%   G:  graph object. Nodes are the same as H, but this time with all Node
%       properties intact and two additional ones added: 'parent' lists the
%       name of the parent node in the trophic grouping hierarchy, and 'TG'
%       lists the trophic group index for each node.  Edges are the same as
%       in the input graph object.

% Copyright 2015-2016 Kelly Kearney

%--------------------------
% Trophic group nodes
%--------------------------

% Names for nodes and their parents

if ~iscell(tg)
    tg = {tg};
end

[nm1, nm2] = deal(cell(size(tg)));
for it = 1:length(tg)
    if it == 1
        nm1{it} = G.Nodes.Name;
    else
        nm1{it} = arrayfun(@(x) sprintf('lev%d_%02d',it-1,x), (1:length(tg{it}))', 'uni', 0);
    end
    nm2{it} = arrayfun(@(x) sprintf('lev%d_%02d',it,x), tg{it}, 'uni', 0);
end
nm1 = cat(1, nm1{:});
nm2 = cat(1, nm2{:});

extranm = setdiff([nm2; G.Nodes.Name], nm1);
nm1 = [nm1; extranm];
nm2 = [nm2; repmat({'web'}, size(extranm))];

% Add new nodes to the graph

newnode = setdiff([nm1; nm2], G.Nodes.Name);
nnew = length(newnode);
new = table(newnode, zeros(nnew,1), ones(nnew,1)*5, nan(nnew,1), ...
    'variableNames', {'Name', 'B', 'type', 'TL'});

G = addnode(G, new);
widx = findnode(G, 'web');
G.Nodes.type(widx) = 6;

% Trophic group hierarchy adjacency matrix

ntot = numnodes(G);
hadj = zeros(ntot);

cidx = findnode(G, nm1);
pidx = findnode(G, nm2);
idx = sub2ind([ntot ntot], cidx, pidx);
hadj(idx) = 1;

H = digraph(hadj, G.Nodes.Name);

% Add parent fields to original graph nodes

parent = cell(ntot,1);
parent(cidx) = G.Nodes.Name(pidx);
isemp = cellfun('isempty', parent);
[parent{isemp}] = deal('null');

G.Nodes.parent = parent;

% Add TG field

G.Nodes.TG = zeros(numnodes(G),1);
G.Nodes.TG(1:length(tg{1})) = tg{1};

% % Add mass flux links as "imports" field, as needed for d3 layout-parsers
% 
% madj = adjacency(G);
% imports = cell(ntot,1);
% for ii = 1:ntot
%     imports{ii} = G.Nodes.Name(madj(:,ii) > 0);
% end
% G.Nodes.imports = imports;
