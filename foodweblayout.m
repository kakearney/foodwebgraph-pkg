function [G, Ax] = foodweblayout(G, tg, varargin)
%FOODWEBLAYOUT Wrapper to call foodweblayout.js
%
%

p = inputParser;
p.addParameter('folder', tempname, @(x) validateattributes(x, {'char'}, {}));

p.parse(varargin{:});
Opt = p.Results;


if nargin < 2
    tg = (1:numnodes(G))';
end

[Htg, Gtg] = trophicgroupgraph(G, tg);

% Create folder with copies of necessary files

if ~exist(Opt.folder, 'dir')
    mkdir(Opt.folder);
end

graph2json(Gtg, fullfile(Opt.folder, 'foodweb.json'), true);

[stat, h] = web(fullfile(Opt.folder, 'index.html'));

% Wait for the user to check the DONE box, then dump the html from the
% browser into a file

while ~checkstatus(h)
    pause(1);
end

txt = get(h, 'HtmlText');
close(h);

file = [tempname '.html'];
fid = fopen(file, 'wt');
fprintf(fid, '%s', txt);
fclose(fid);

% Extract relevant details from the html

[C,T] = extractsvg(file);
[G, Ax] = foodweblayoutdetails(G, C, T);

% Subfunction: Parse html and check for DONE marker

function isdone = checkstatus(h)

txt = get(h, 'HtmlText');
if isempty(txt)
    isdone = false;
else
    idx1 = strfind(txt, '<div id="textparams">');
    if isempty(idx1)
        isdone = false;
    else
        idx2 = strfind(txt, '</div>');
        idx2 = idx2(find(idx2 > idx1,1)) + 6;
        divtext = txt(idx1:idx2);
        isdone = ~isempty(strfind(divtext, 'DONE'));
    end
end







