function fwgraph2json(G, tg, filename, copyjs)
%FWGRAPH2JSON Saves food web graph object to a JSON file
%
% fwgraph2json(G, tg, filename)
% fwgraph2json(G, tg, filename, copyjs)
%
% This function exports an ecopathmodel graph object to a JSON file,
% allowing it to be used as input for the d3-foodweb plotting tools.
%
% Input variables:
%
%   G:          food web graph object, with the same properties as an
%               ecopathmodel graph object. Properties expected are:
%               Node table:
%                   Name:   string, value used for node IDs
%                   B:      value used for node size scaling, typically
%                           biomass 
%                   TL:     value used for y-positioning, typically trophic
%                           level
%                   type:   optional, used for node color scaling if no
%                           cval property is found 
%                   cval:   optional, value used for node color scaling
%               Edge table:
%                   EndNodes:   source and target node names for each edge
%                   xpath:      optional, x-coordinates of edge path for
%                               d3.foodwebstatic
%                   ypath:      optional, y-coordinates of edge path for
%                               d3.foodwebstatic
%                   Weight:     optional, value used to scale edge width in
%                               d3.foodwebstatic
%
%   tg:         nnode x 1 array, trophic group indices, used to add TG
%               property to Node property table.  If empty, will assume
%               each node is in its own group.  
%
%   filename:   path to new file, with or without .json extension.
%
%   copyjs:     logical scalar, if true, copy to the same folder as the new
%               .json file all files necessary to run the foodwebgraph
%               layout utility in the same location.  Default if not
%               included is false.
%   

% Copyright 2017 Kelly Kearney

% Check input

validateattributes(G, {'digraph'}, {'scalar'});
ng = numnodes(G);

validateattributes(tg, {'numeric'}, {'vector', 'nonnegative','integer', 'numel', ng});
tg = tg(:);

validateattributes(filename, {'char'}, {});

if nargin < 4
    copyjs = false;
end
validateattributes(copyjs, {'logical'}, {'scalar'});

% Parse file name to find folder, and add .json if necessary

[pth,~,ex] = fileparts(filename);
if isempty(pth)
    pth = '.';
end

if ~strcmp(ex, '.json')
    filename = [filename '.json'];
end

if ~exist(pth, 'dir')
    error('Specified folder for new file (%s) does not exist', pth);
end

% Add trophic group and id properties to node table

G.Nodes.id = G.Nodes.Name;
G.Nodes.TG = tg;

% Construct trophic group linkages graph

adj = bsxfun(@eq, tg, tg');
adj = triu(adj .* ~eye(ng));
G2 = digraph(adj, G.Nodes.id);

idx = findnode(G, G2.Edges.EndNodes(:,1));

% Construct node table, edge table, and trophic links table

J.nodes = table2struct(G.Nodes);
J.nodes = rmfield(J.nodes, 'Name');

edgedata = {'source', G.Edges.EndNodes(:,1), 'target', G.Edges.EndNodes(:,2)};
if all(ismember({'x','y'}, G.Edges.Properties.VariableNames))
    xpath = cellfun(@(x) x(:)', G.Edges.x, 'uni', 0);
    ypath = cellfun(@(x) x(:)', G.Edges.y, 'uni', 0);
    edgedata = [edgedata, {'xpath', xpath, 'ypath', ypath}];
end
if ismember('Weight', G.Edges.Properties.VariableNames)
    edgedata = [edgedata, {'Weight', num2cell(G.Edges.Weight)}];
end

J.links1 = struct(edgedata{:});

% if all(ismember({'x','y'}, G.Edges.Properties.VariableNames))
%     xpath = cellfun(@(x) x(:)', G.Edges.x, 'uni', 0);
%     ypath = cellfun(@(x) x(:)', G.Edges.y, 'uni', 0);
%     J.links1 = struct(...
%         'source', G.Edges.EndNodes(:,1),  ...
%         'target', G.Edges.EndNodes(:,2), ...
%         'xpath', xpath, ...
%         'ypath', ypath, ...
%         'Weight', num2cell(G.Edges.Weight));
% else
%     J.links1 = struct(...
%         'source', G.Edges.EndNodes(:,1),  ...
%         'target', G.Edges.EndNodes(:,2));
% end
J.links2 = struct(...
    'source', G2.Edges.EndNodes(:,1), ...
    'target', G2.Edges.EndNodes(:,2), ...
    'TG', num2cell(G.Nodes.TG(idx)));

% Write to file

Jopt = struct('Compact', 0, ...
              'FileName', filename, ...
              'NoRowBracket', 1); 
savejson('', J, Jopt);

% Set up sandbox environment in the targeted folder

if copyjs

    thisfile = mfilename('fullpath');
    jsdir = fileparts(thisfile);
    
    copyfile(fullfile(jsdir, 'foodweblayouttool.html'), fullfile(pth, 'index.html'));
    copyfile(fullfile(jsdir, 'd3-foodweb.min.js'), pth);
    
end


