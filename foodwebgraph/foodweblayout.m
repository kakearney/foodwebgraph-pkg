function [G, Ax, Param, Svg] = foodweblayout(G, tg, varargin)
%FOODWEBLAYOUT Wrapper to call foodweblayout.js
%
% This function provides a wrapper around the d3-foodweb foodweblayouttool,
% preparing the input .json file and then extracting details from the
% user-configured food web.
%
% [G, Ax] = foodweblayout(G, tg)
% [G, Ax] = foodweblayout(G, tg, 'folder', folder)
% 
% Input variables:
%
%   G:      digraph object, where nodes represent functional groups and
%           edges represent the fluxes between them.  Should be consistent
%           with that returned by the ecopathmodel graph method.
%
%   tg:     nnode x n array, trophic group indices.  Nodes sharing an index
%           will be linked to each other in the force layout tool.  See
%           trophicgroupgraph for further details.
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
%   G:      digraph object, same as input with the following parameters
%           added to the Nodes table: 
%           y:          y-coordinate of node center, in trophic level units
%           x:          x-coordinate of node center, in trophic level
%                       units.  The units in the x-direction are arbitrary,
%                       and are kept the same as the y-units to make
%                       plotting of circles simpler.
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
%   Param:  structure with user-set layout parameters
%
%   Svg:    structure with details of elements extracted from the svg.  See
%           extractsvg.m for details.

% Copyright 2016 Kelly Kearney

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
[G, Ax] = foodweblayoutdetails(G, C, T, [Param.tlmin Param.tlmax]);

Svg.C = C;
Svg.T = T;

% Subfunction: Parse html and check for DONE marker

function checkstatus(obj, ev, hweb)
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





