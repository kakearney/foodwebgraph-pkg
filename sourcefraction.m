function [frac, A] = sourcefraction(G, tagnodes)
%SOURCEFRACTION Calculates fraction of fluxes at each node from source
%
% [frac, A] = sourcefraction(G, tagnodes)
%
% This function calculates the fraction of the fluxes reaching a particular
% node in a graph that originate from certain specified source nodes.  A
% source node is defined as a node in the web that has no predecessors.
% Tagged source nodes are assigned a value of 1, and untagged ones a value
% of 0; all other nodes are assigned the average value of their predecessor
% nodes, weighted by the value of the edges connecting them to those
% predecessor nodes (values between 0 and 1).
%
% This calculation is based on the inverse method of Levine (1980) for
% calculating trophic position:
%
% Levine S (1980) Several measures of trophic structure applicable to
% complex food webs. J Theor Biol 83:195?207 
% DOI: 10.1016/0022-5193(80)90288-X
%
% Input variables:
%
%   G:          a digraph object
%
%   tagnodes:   string, or cell array of strings.  Names of node(s) to tag;
%               must be source nodes, i.e. nodes with no predecessors. 
%
% Output variables:
%
%   frac:       nnode x 1 array holding fraction of fluxes each node
%               receives from the tagged node(s)
%
%   A:          1 x 1 structure with the following fields:
%
%               srcnodes:   nsrc x 1 cell array of strings, names of all
%                           source nodes
%               srcval:     nsrc x 1 vector, value assigned to each source
%                           node 

% Copyright 2016 Kelly Kearney

% Set value for each node, 1 for tagged, 0 for others

nnode = numnodes(G);

bidx = findnode(G, tagnodes);
if ~all(bidx)
    str = sprintf('%s,', tagnodes{bidx==0});
    error('Specified tagged nodes not found (%s)', str(1:end-1));
end

val = zeros(nnode,1);
val(bidx) = 1;

% Build transition matrix

edgeidx = reshape(findnode(G, G.Edges.EndNodes), [], 2);

adj = full(sparse(edgeidx(:,1), edgeidx(:,2), G.Edges.Weight, nnode, nnode));
adj = bsxfun(@rdivide, adj, sum(adj,1));
adj(isnan(adj)) = 0;

T = adj';

% Sort so source nodes come first

issrc = sum(T,2) == 0;
srcidx = find(issrc);

tf = ismember(bidx, srcidx);
if ~all(tf)
    str = sprintf('%s,', tagnodes{~tf});
    error('Specified tagged node (%s) is not a source node', str(1:end-1));
end

nsrc = length(srcidx);
T(sub2ind([nnode nnode], srcidx, srcidx)) = 1;

[~,isrt] = sort(issrc, 'descend');
T = T(isrt,isrt);

sval = val(isrt);
sval = sval(1:nsrc);

% Solve for y (i.e. bfrac)

Itmp = T(1:nsrc,1:nsrc); % should be eye(nsrc)
Q = T(nsrc+1:end,nsrc+1:end);
R = T(nsrc+1:end,1:nsrc);
ztmp = T(1:nsrc,nsrc+1:end); % should be 0s

bfrac = (Q - eye(nnode-nsrc))\(-R*sval);

frac = nan(nnode,1);
frac(isrt) = [sval; bfrac];

A.srcnodes = G.Nodes.Name(isrt(1:nsrc));
A.srcval = sval;


