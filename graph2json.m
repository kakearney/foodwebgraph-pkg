function graph2json(G, filename, copyjs)
%GRAPH2JSON Saves food web graph object to a JSON file
%
% graph2json(G, filename, folder, copyjs)
%
% Input variables:
%
%   G:          food web graph object, see second output of
%               trophocgroupgraph.m 
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

E = struct('source', G.Edges.EndNodes(:,1), 'target', G.Edges.EndNodes(:,2));

Jopt = struct('Compact', 0, ...
              'FileName', filename, ...
              'NoRowBracket', 1); 
Json.nodes = table2struct(G.Nodes)';
Json.edges = E;
savejson('', Json, Jopt);

% Set up sandbox environment in the targeted folder

if copyjs
    jsfiles = {'index.html', 'foodweblayout.js', 'labeler.js', ...
        'packages.js', 'seedrandom.min.js'};

    thisfile = mfilename('fullpath');
    jsdir = fileparts(thisfile);
    
    for ii = 1:length(jsfiles)
        copyfile(fullfile(jsdir, jsfiles{ii}), pth);
    end
    
end


