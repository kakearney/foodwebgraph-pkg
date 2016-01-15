function varargout = renderfoodwebforce(G, varargin)
%RENDERFOODWEBFORCE Render food web diagram via foodweblayout.js
%
% [stat, h, File] = renderfoodwebforce(G, 'mode', 'test', p1, v1, ...)
% [C, T, File]    = renderfoodwebforce(G, 'mode', 'extract', p1, v1, ...)
%
% This function provides a wrapper environment for running the
% foodweblayout.js script.  It creates a simple web page to call the script
% with the appropriate input and parameters, and either shows that page
% using the Matlab web browser, or runs the page through phantomJS to
% render and extract details from the page.  For the latter option,
% phantomJS (http://phantomjs.org/) must be installed locally.
%
% Input variables:
%
%   G:          A graph object.  The graph's Nodes table must include the
%               following variables:
%               Name:   name of each node
%               B:      biomass value of each node, used for size scaling.
%               type:   type of node
%                       0 = consumer
%                       1 = producer
%                       2 = detritus
%                       3 = fleet
%                       4 = out-of-system
%                       5 = trophic group node
%                       6 = top-level web grouping node
%               TL:     trophic level, used for trophic nudging.  If value
%                       is NaN, no nudging will be applied to that node,
%               parent: name of node above each one in the trophic grouping
%                       hierarchy
%               TG:     index of trophic group to which each node belongs.
%                       Value is 0 for non-leaf nodes (i.e. trophic group
%                       nodes).
%
% Optional input variables (passed as parameter/value pairs):
%
%   mode:       mode to run ['test']
%               'test':     render the resulting webpage via Matlab's web
%                           browser 
%               'extract':  render via phantomJS and extract the svg
%                           elements 
%               
%   jsdir:      folder where javascript libraries reside.  If empty,
%               assumed to be in the same location as this function. []
% 
%   marl:       left margin (pixels) for main svg element grouping [50]
%
%   marr:       right margin (pixels) for main svg element grouping [50]
%
%   mart:       top margin (pixels) for main svg element grouping [50]
%
%   marb:       bottom margin (pixels) for main svg element grouping [50]
%
%   totwidth:   width (pixels) of svg canvas [800]
%
%   totheight:  height (pixels) of svg canvas [800]
%
%   padding:    minimum separation between nodes [10]
%
%   tlnudge:    strength of trophic level nudging [2.0]
%
%   linkdist:   link distance for force-directed layout [5]
%
%   charge:     charge for force-directed layout [-100]
%
%   hlineop:    opacity for temporary trophic group lines [0.5]
%
%   fontsz:     font size for node labels [8]
%
%   seed:       string used to seed the random number generator ['hello']
%
%   cpfac:      scaling factor for the box used for initial circle packing,
%               relative to svg canvas size [0.8]
%
%   drawmode:   drawing mode ['cpack']:
%               'cpack':    plot initial circle-packing layout only
%               'force':    plot force-directed layout
%
%   init:       initialization mode ['cpack']:
%               'cpack':    initialize using circle-packing layout
%               'coords':   initialize using coordinates in the JSON file
%
% Output variables:
%
%   stat:       status of web command (0 = success, 1 = browser not found,
%               2 = browser found but couldn't be launched).  'test' mode
%               only.
%
%   h:          handle to the active web browser. 'test' mode only
%
%   C:          structure holding details of circle elements found in the
%               resulting svg image, with the following fields
%               x:      x-coordinate of nodes, in pixels (measured from
%                       left to right)
%               y:      y-coordinate of nodes, in pixels (measured from top
%                       to bottom)
%               xref:   x-coordinate of corner reference dots, in pixels
%               yref:   y-coordinate of corner reference dots, in pixels
%               r:      radius of nodes, in pixels
%               label:  node name associated with each x, y, and r value
%
%   T:          structure holding details of the text elements found in the
%               resulting svg image, with the following fields:
%               x:      x-coordinate of the text, in pixels
%               y:      y-coordinate of text, in pixels
%               label:  text string
%               anchor: anchor point described by the x/y coordinates.  If
%                       'middle', the coordinates are centered both
%                       vertically and horizontally.  If 'start', the
%                       coordinates reference the lower left corner of the
%                       text extent.
%
%   F:          structure holding path names to the temporary files created
%               by this function.
%               json:   the JSON file used as input, holding the G.Nodes
%                       table
%               html:   html file, includes a single <div> object and a
%                       script calling foodweblayout.js with the
%                       appropriate input data and parameters.  Note that
%                       most browsers will require a local server to be
%                       running in order to render this file due to the
%                       cross-origin call to the D3 library (Matlab's
%                       browser does not have this restriction).
%               html2:  The html file created after passing the previous
%                       file through a phantomJS waitFor() script.
%                       Includes the final rendered SVG element within the
%                       <div> element.
%               svg:    The SVG elements from the previous file, extracted
%                       into a file as a standalone image.

% Copyright 2015 Kelly Kearney

% Parse input

Opt.mode      = 'test';
Opt.jsdir     = []; 

Opt.marl      = 50;
Opt.marr      = 50;
Opt.mart      = 50;
Opt.marb      = 50;
Opt.totwidth  = 800;
Opt.totheight = 600;
Opt.padding   = 10;
Opt.tlnudge   = 2.0;
Opt.linkdist  = 5;
Opt.charge    = -100;
Opt.hlineop   = 0.5;
Opt.fontsz    = 8;
Opt.seed      = 'hello';
Opt.cpfac     = 0.8;
Opt.tllim     = [NaN NaN];
Opt.drawmode  = 'force';
Opt.init      = 'cpack';

Opt = parsepv(Opt, varargin);

% Where are javascript scripts located?

if isempty(Opt.jsdir)
    thisfile = mfilename('fullpath');
    Opt.jsdir = fileparts(thisfile);
end

% Save Nodes structure to JSON file

File.json = [tempname '.json'];

Jopt = struct('Compact', 0, ...
              'FileName', File.json, ...
              'NoRowBracket', 1); 
Json.nodes = table2struct(G.Nodes)';
savejson('', Json, Jopt);


% Create html file

htmltxt = {...
        '<!DOCTYPE html>'                                                       
        '<html>'                                                                 
        '<head>'                                                                 
        '    <meta http-equiv="content-Type" content="text/html; charset=utf-8"/>'
        '</head>'                                                                
        ''                                                                       
        '<!--'                                                                   
        '***************'                                                        
        'Main page'                                                              
        '***************'                                                        
        '-->'                                                                    
        ''                                                                       
        '<!-- Location for SVG canvas -->'                                       
        ''                                                                      
        '<div id="fwf"></div>'                                                   
        ''                                                                      
        '<!--'                                                               
        '***************'                                                        
        'Scripts'                                                                
        '***************'                                                        
        '-->'                                                                    
        ''                                                                       
        '<script src="http://d3js.org/d3.v3.min.js" charset="utf-8"></script>'   
sprintf('<script src="%s"></script>', fullfile(Opt.jsdir, 'packages.js'))                                   
sprintf('<script src="%s"></script>', fullfile(Opt.jsdir, 'foodweblayout.js'))                               
sprintf('<script src="%s"></script>', fullfile(Opt.jsdir, 'seedrandom.min.js'))                             
sprintf('<script src="%s"></script>', fullfile(Opt.jsdir, 'labeler.js'))                           
        ''                                                                       
        '<script>'                                                               
        ''                                                                       
sprintf('d3.json("%s", function(data) {', File.json)                                  
        '    d3.select("#fwf")'                                                  
        '        .datum(data)'                                                   
        '    .call(foodweblayout()'                                              
sprintf('	    .marl(%d)',      Opt.marl)
sprintf('	    .marr(%d)',      Opt.marr) 
sprintf('	    .mart(%d)',      Opt.mart) 
sprintf('	    .marb(%d)',      Opt.marb)
sprintf('	    .totwidth(%d)',  Opt.totwidth) 
sprintf('	    .totheight(%d)', Opt.totheight) 
sprintf('	    .padding(%d)',   Opt.padding) 
sprintf('	    .tlnudge(%.2f)', Opt.tlnudge) 
sprintf('	    .linkdist(%.1f)',Opt.linkdist)
sprintf('	    .charge(%d)',    Opt.charge) 
sprintf('	    .hlineop(%.2f)', Opt.hlineop)
sprintf('	    .fontsz(%d)',    Opt.fontsz) 
sprintf('	    .seed(''%s'')',      Opt.seed)
sprintf('	    .cpfac(%.2f)',   Opt.cpfac) 
sprintf('	    .tllim([%.2f,%.2f])',   Opt.tllim) 
sprintf('       .init(''%s'')', Opt.init)
sprintf('	    .drawmode(''%s''));',  Opt.drawmode)                                           
        '});'                                                                    
        ''                                                                       
        '</script>'
};

% Render webpage

File.html = [tempname '.html'];
printtextarray(htmltxt, File.html);

switch Opt.mode
    case 'test'
        
        % Preview in Matlab browser
        
        [stat, h] = web(File.html);
        tmp = {stat, h, File};
        varargout = tmp(1:nargout);
        
    case 'extract'
        
        % Render svg file using phantomJS
        
        File.html2 = [tempname '.html'];
        File.svg = [tempname '.svg'];
        
        
        cmd = sprintf('phantomjs %s %s > %s', ...
            fullfile(Opt.jsdir, 'waitforfoodweb.js'), ...
            File.html, File.html2);
        [s,r] = system(cmd);
        if s
            error('Error generating SVG via phantomJS: %s', r);
        end
        
        extractsvg(File.html2, File.svg);
        
        % Extract node positions, radii, and labels, plus reference point
        % locations

        [C,T] = parsesvg(File.svg);
        
        tmp = {C, T, File};
        varargout = tmp(1:nargout);
end

%-----------------------------------
% Parse circle coordinates from
% svg file
%-----------------------------------

function [C,T] = parsesvg(file)

xdoc = xmlread(file);
circ = xdoc.getElementsByTagName('circle');

n = circ.getLength;

[C.x ,C.y, C.xref, C.yref, C.r] = deal(nan(n,1));
C.label = cell(n,1);

for k = 1:n
    c = circ.item(k-1);
    
    if c.hasAttribute('cx')
        C.xref(k) = str2double(char(c.getAttribute('cx')));
        C.yref(k) = str2double(char(c.getAttribute('cy')));
    else
        C.r(k) = str2double(char(c.getAttribute('r')));
        t = c.getElementsByTagName('title');
        p = c.getParentNode;
        
        C.label{k} = char(t.item(0).getFirstChild.getData);
        tform = char(p.getAttribute('transform'));
        
        xy = regexp(tform, 'translate\(([\d\.-]+),([\d\.-]+)\)', 'tokens', 'once');
        xy = cellfun(@str2double, xy);
        C.x(k) = xy(1);
        C.y(k) = xy(2);
        
    end
        
end
    
isref = isnan(C.x);
C.xref = C.xref(isref);
C.yref = C.yref(isref);
C.x = C.x(~isref);
C.y = C.y(~isref);
C.r = C.r(~isref);
C.label = C.label(~isref);

txt = xdoc.getElementsByTagName('text');
nt = txt.getLength;

[T.x, T.y] = deal(nan(nt,1));
[T.label, T.anchor] = deal(cell(nt,1));
for it = 1:nt
    t = txt.item(it-1);
    T.x(it) = str2double(char(t.getAttribute('x')));
    T.y(it) = str2double(char(t.getAttribute('y'))); 
    T.anchor{it} = char(t.getAttribute('text-anchor')); 
    if t.getLength > 0
        T.label{it} = char(t.item(0).getData);
    else
        T.label{it} = '';
    end
end

%-----------------------------------
% Pull out the svg portion of file
%-----------------------------------

function extractsvg(htmlfile, svgfile)

txt = fileread(htmlfile);
txt = regexp(txt, '<svg.*</svg>', 'match', 'once');
    
newtxt = '<svg xmlns="http://www.w3.org/2000/svg" xmlns:xlink="http://www.w3.org/1999/xlink"';
txt = regexprep(txt, '<svg', newtxt);

fid = fopen(svgfile, 'wt');
fprintf(fid, '%s', txt);
fclose(fid);
