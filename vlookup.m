function tval = vlookup(tbl, val, col, lookcol, fillval)
%VLOOKUP Lookup table
%
% Return variables in one matrix based on matching values, lookup table
% style.
%
% Input variables:
%
%   tbl:        table array
%
%   val:        elements of tbl.lookcol to look for
%
%   col:        name of column to return values from
%
%   lookcol:    name of column to look for e values
%
%   fillval:    value to return if an element of e isn't found
%
% Output variables:
%
%   tval:       array same size as e, holding corresponding values from
%               column col

% Copyright 2015 Kelly Kearney

nv = numel(val);

if isnumeric(tbl.(col))
    if nargin < 5
        fillval = NaN;
    end
    tval = ones(size(val))*fillval;
else
    tval = cell(size(val));
    if nargin >= 5
        [tval{:}] = deal(fillval);
    end
end


if isnumeric(val) && isnumeric(tbl.(lookcol))
    [tf,loc] = ismember(val, tbl.(lookcol));
elseif iscell(val) && iscell(tbl.(lookcol)) && all(cellfun(@ischar,val)) && all(cellfun(@ischar,tbl.(lookcol)))
    [tf,loc] = ismember(val, tbl.(lookcol));
else
    % TODO
    error('Haven''t written for other non-numeric, cell array of strings data');
    
end

tval(tf) = tbl.(col)(loc(tf));














