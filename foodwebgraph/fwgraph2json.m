function fwgraph2json(G, tg, filename, copyjs)
%FWGRAPH2JSON Saves food web graph object to a JSON file
%
% graph2json(G, tg, filename, folder, copyjs)
%
% Input variables:
%
%   G:          food web graph object, with the same properties as an
%               ecopathmodel graph object (i.e Node table includes Name, B,
%               type, and TL propeties)  
%
%   tg:         nnode x 1 array, trophic group indices.  If empty, will
%               assume each node is in its own group.
%
%   filename:   name for file, with or without .json extension
%
%   copyjs:     logical scalar, if true, copy to the same folder as the new
%               .json file all files necessary to run the foodwebgraph
%               layout utility in the same location
%   

% Create JSON file

[pth,fl,ex] = fileparts(filename);
if isempty(pth)
    pth = '.';
end

if ~strcmp(ex, '.json')
    filename = [filename '.json'];
end

G.Nodes.id = G.Nodes.Name;
G.Nodes.TG = tg;

ng = numnodes(G);
adj = bsxfun(@eq, tg, tg');
adj = triu(adj .* ~eye(ng));
G2 = digraph(adj, G.Nodes.id);

idx = findnode(G, G2.Edges.EndNodes(:,1));

J.nodes = table2struct(G.Nodes);
J.nodes = rmfield(J.nodes, 'Name');
if all(ismember({'x','y'}, G.Edges.Properties.VariableNames))
    xpath = cellfun(@(x) x(:)', G.Edges.x, 'uni', 0);
    ypath = cellfun(@(x) x(:)', G.Edges.y, 'uni', 0);
    J.links1 = struct('source', G.Edges.EndNodes(:,1),  'target', G.Edges.EndNodes(:,2), ...
        'xpath', xpath, 'ypath', ypath, 'Weight', num2cell(G.Edges.Weight));
else
    J.links1 = struct('source', G.Edges.EndNodes(:,1),  'target', G.Edges.EndNodes(:,2));
end
J.links2 = struct('source', G2.Edges.EndNodes(:,1), 'target', G2.Edges.EndNodes(:,2), 'TG', num2cell(G.Nodes.TG(idx)));

Jopt = struct('Compact', 0, ...
              'FileName', filename, ...
              'NoRowBracket', 1); 
savejson('', J, Jopt);

% Set up sandbox environment in the targeted folder

if copyjs
    jsfiles = {'index.html', 'foodweblayout.js', 'labeler.js'};

    thisfile = mfilename('fullpath');
    jsdir = fileparts(thisfile);
    
    for ii = 1:length(jsfiles)
        copyfile(fullfile(jsdir, jsfiles{ii}), pth);
    end
    
end


