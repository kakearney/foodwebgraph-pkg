function [G, Ax, Param, Svg] = foodweblayout(G, tg, varargin)
%FOODWEBLAYOUT Wrapper to call d3.foodweblayout.js
%
% This function provides a wrapper around the d3-foodweb d3.foodweblayout
% tool.  It prepares a properly-formatted .json file and then launches a
% webpage that allows one to interactively test parameters for the node
% positioning algorithm; users can also make manual adjustments by clicking
% on nodes and dragging them to a desired location. Once the user has
% decided on a set of parameters and made any desired manual adjustments to
% the layout, details of the layout (node position, node size, text label
% position and alignment) are extracted from the SVG and added to the Node
% properties table of the input graph.
%
% [Gout, Ax, Param, Svg] = foodweblayout(Gin, tg)
% [Gout, Ax, Param, Svg] = foodweblayout(Gin, tg, 'folder', folder)
% 
% Input variables:
%
%   Gin:    digraph object, where nodes represent functional groups and
%           edges represent the fluxes between them.  Should be consistent
%           with that returned by the ecopathmodel graph method.
%
%   tg:     nnode x n array, trophic group indices.  Nodes sharing an index
%           will be linked to each other in the force layout tool.
%
% Optional input variables (passed as parameter/value pairs):
%
%   folder: path name to folder where html file (and supporting javascript
%           files) should be saved.  Providing this detail is necessary if
%           you wish to view the file in an external browser.  If not
%           included, files will be added to a folder in the temporary
%           directory.
%
% Output variables:
%
%   Gout:   digraph object, same as input digraph with the following
%           parameters added to the Nodes table:  
%           y:          y-coordinate of node center, in trophic level units
%           x:          x-coordinate of node center, in trophic level
%                       units.  The units in the x-direction are arbitrary,
%                       and use the same pixels-to-trophic level units
%                       conversion as the y-units to maintain a 1:1 aspect
%                       ratio and make plotting of circles simpler.
%           r:          radius of node, in trophic level units
%           tx:         x-coordinate of node label, in trophic level units
%           ty:         y-coordinate of node label, in trophic level units
%           th:         horizontal alignment of node label
%           tv:         vertical alignment of node label
%
%   Ax:     structure holding information about the axis that is required
%           to replicate the svg canvas:
%           ylim:       y limits
%           xlim:       x limits
%           figpos:     position vector for figure (in pixels)
%           axpos:      position vector for axis (normalized)
%           fontsize:   fontsize for text labels
%           fontname:   fontname for text labels
%
%   Param:  structure with user-set layout parameters (these correspond to
%           the input arguments for d3.foodweblayout, and can be used to
%           replicate the SVG canvas setup and node scaling details with
%           d3.foodwebstatic). 
%
%   Svg:    structure with details of elements extracted from the svg.  See
%           extractsvg.m for details.

% Copyright 2016-2017 Kelly Kearney

p = inputParser;
p.addParameter('folder', tempname, @(x) validateattributes(x, {'char'}, {}));

p.parse(varargin{:});
Opt = p.Results;

if nargin < 2
    tg = (1:numnodes(G))';
end

% Create folder with copies of necessary files

if ~exist(Opt.folder, 'dir')
    mkdir(Opt.folder);
end

fwgraph2json(G, tg, fullfile(Opt.folder, 'foodweb.json'), true);

[stat, h] = web(fullfile(Opt.folder, 'index.html'));

% Wait for the user to check the DONE box, then dump the html from the
% browser into a file

tobj = timer('TimerFcn', {@checkstatus, h}, ...
             'ExecutionMode', 'fixedRate', ...
             'StartDelay', 1, ...
             'TasksToExecute', 3600);
start(tobj);
wait(tobj);

if strcmp(get(tobj, 'UserData'), 'browserclosed')
    error('Web browser closed prematurely');
end

txt = get(h, 'HtmlText');
close(h);

file = [tempname '.html'];
fid = fopen(file, 'wt');
fprintf(fid, '%s', txt);
fclose(fid);

% Extract the parameter settings set by the user

linemarker = '<div id="textparams"';
txt = regexp(txt, '\n', 'split');
isparam = strncmp(strtrim(txt), linemarker, length(linemarker));
paramtxt = txt{isparam};
idx1 = strfind(paramtxt, 'marl');
idx2 = strfind(paramtxt, '<h1>DONE');
paramtxt = paramtxt(idx1:(idx2-1));
pv = regexp(paramtxt, '<br>', 'split');
isemp = cellfun('isempty', pv);
pv = pv(~isemp);
pv = regexp(pv, '=', 'split');
pv = strtrim(cat(1, pv{:}));
Param = cell2struct(pv(:,2), pv(:,1));
Param = structfun(@str2double, Param, 'uni', 0);

% Extract relevant details from the html

if isnan(Param.tlmin)
    Param.tlmin = min(G.Nodes.TL);
end
if isnan(Param.tlmax)
    Param.tlmax = max(G.Nodes.TL);
end 

[C,T] = extractsvg(file);
[G, Ax] = fwsvgdetails(G, C, T);

Svg.C = C;
Svg.T = T;

% Subfunction: Parse html and check for DONE marker

function checkstatus(obj, ~, hweb)
txt = get(hweb, 'HtmlText');
if isempty(txt)
    set(obj, 'UserData', 'browserclosed');
    isdone = true;
else
    idx1 = strfind(txt, '<div id="textparams"');
    if isempty(idx1)
        isdone = false;
    else
        idx2 = strfind(txt, '</div>');
        idx2 = idx2(find(idx2 > idx1,1)) + 6;
        divtext = txt(idx1:idx2);
        isdone = ~isempty(strfind(divtext, 'DONE'));
    end
end
if isdone
    stop(obj);
end





