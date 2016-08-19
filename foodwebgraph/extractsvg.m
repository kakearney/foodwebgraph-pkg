function [C,T] = extractsvg(file)
%EXTRACTSVG Parse food web graphics information out of an svg element
%
% [C,T] = extractsvg(file)
%
% Input variables:
%
%   file:       name of either an svg file, or an html file holding a
%               single svg element.  The svg should contain a
%               foodweblayout.js-generated food web diagram. 
%
% Output variables:
%
%   C:          structure holding details of circle elements found in the
%               resulting svg image, with the following fields
%               x:      x-coordinate of nodes, in pixels (measured from
%                       left to right, with 0 at the left of the main
%                       grouping element)
%               y:      y-coordinate of nodes, in pixels (measured from top
%                       to bottom, with 0 at the top of the main grouping
%                       element) 
%               xref:   x-coordinate of corner reference dots, in pixels
%                       (these define the main grouping element)
%               yref:   y-coordinate of corner reference dots, in pixels
%               r:      radius of nodes, in pixels
%               label:  node name associated with each x, y, and r value
%               trans:  transform (i.e. shift left, shift down) from upper
%                       right of svg canvas to upper right of main grouping
%                       element 
%               svgsz:  width and height of svg canvas
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

% Copyright 2016 Kelly Kearney

% First, check to see if input is a pure SVG file or HTML-with-svg.  If the
% latter, pull out the svg element and put it in its own file temporarily
% so it can be easily parsed via xmlread.

txt = fileread(file);
ishtml = ~isempty(strfind(txt, '<!DOCTYPE html'));
issvg  = ~isempty(strfind(txt, '<!DOCTYPE svg'));

if ishtml
    idx1 = strfind(txt, '<svg');
    idx2 = strfind(txt, '</svg>');
    if length(idx1) > 1
        error('Found multiple SVG elements in this file; only one allowed');
    end
    svgcontent = txt(idx1:(idx2+5));
    file = [tempname '.svg'];
    fid = fopen(file, 'wt');
    fprintf(fid, '%s', svgcontent);
    fclose(fid);
end
    
% Extract size and position of the circle nodes (as well as the corner dots
% that mark the edge of the "axis" g element.

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

% Extract the position, alignment, and font details from the text labels

txt = xdoc.getElementsByTagName('text');
nt = txt.getLength;

[T.x, T.y] = deal(nan(nt,1));
[T.label, T.anchor] = deal(cell(nt,1));
for it = 1:nt
    t = txt.item(it-1);
    T.x(it) = str2double(char(t.getAttribute('x')));
    T.y(it) = str2double(char(t.getAttribute('y'))); 
    T.anchor{it} = char(t.getAttribute('text-anchor'));
    T.fontsize{it} =  str2double(char(t.getAttribute('font-size')));
    T.font{it} = char(t.getAttribute('font-family'));
    if t.getLength > 0
        T.label{it} = char(t.item(0).getData);
    else
        T.label{it} = '';
    end
end

% Extract the location and size of the main svg element and "axis" grouping
% element

g = xdoc.getElementsByTagName('g');
gmain = g.item(0);
trans = char(gmain.getAttribute('transform'));
C.trans = cellfun(@str2double, regexp(trans, 'translate\(([0-9]*),([0-9]*)\)', 'tokens', 'once'));

svg = xdoc.getElementsByTagName('svg');
svg = svg.item(0);
C.svgsz = [str2double(char(svg.getAttribute('width'))), ...
           str2double(char(svg.getAttribute('height')))];
