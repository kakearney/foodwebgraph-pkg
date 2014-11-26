function printtextarray(array, filename, removeLastEndline)
%PRINTTEXTARRAY Prints text in a cell array to a file
%
% printtextarray(array, filename)
%
% This program prints the data in a cell array of strings to a file.
%
% Input variables:
%
%   array:      n x 1 or 1 x n cell array of strings
%
%   filename:   name of text file to print to

% Copyright 2004-2005 Kelly Kearney

if ~isvector(array)
    error('Cell array must be a vector cell array of strings');
end

fid = fopen(filename, 'wt');
for i = 1:length(array)-1
    fprintf(fid, '%s\n', array{i});
end
if nargin == 3  && removeLastEndline
    fprintf(fid, '%s', array{end});
else
    fprintf(fid, '%s\n', array{end});
end
    
fclose(fid);