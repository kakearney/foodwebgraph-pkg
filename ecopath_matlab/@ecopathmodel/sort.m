function A = sort(A, varargin)
%SORT Sort groups and fleets in an ecopathmodel object
%
% B = sort(A, grporder)
% B = sort(A, grporder, fltorder)
%
% This function rearranges the tables in an ecopathmodel object based on a
% prescribed order. This function allows any order, but note that most
% ecopath-related calculations require detrital groups to be listed after
% live groups; a warning will be issued if the prescribed group order
% breaks this rule.
%
% Input variables:
%
%   A:          ecopathmodel object
%
%   grporder:   any permutation of the integers 1:A.ngroup.  Values
%               correspond to the indices of groups in the current model.
%
%   fltorder:   any permutation of the integers 1:A.ngear. Values
%               correspond to the indices of fleets in the current model.
%
% Ouptut variables:
%
%   B:          ecopathmodel object, same as A but with groups and fleets
%               rearranged as requested

% Copyright 2016 Kelly Kearney


% Parse and check input

p = inputParser;
p.addOptional('grporder', 1:A.ngroup, @(x) validateattributes(x, {'numeric'}, {'vector', 'numel', A.ngroup}));
p.addOptional('fltorder', 1:A.ngear,  @(x) validateattributes(x, {'numeric'}, {'vector', 'numel', A.ngear}));
p.parse(varargin{:});

grporder = p.Results.grporder;
fltorder = p.Results.fltorder;

% Check detritus group locations

detidx = find(A.groupdata.pp == 2);

[~, detloc] = ismember(detidx, grporder);
[~, detorder] = sort(detloc);

if any(detloc < A.nlive)
    warning('Detrital groups must follow live groups for ecopath-related calculations');
end

% Resort tables

A.groupdata = A.groupdata(grporder,:);
A.dc = A.dc(grporder,grporder);
A.landing = A.landing(grporder, fltorder);
A.discard = A.discard(grporder, fltorder);
A.df = A.df(grporder, detorder);
A.discardFate(fltorder, detorder);

A.name = A.name(grporder);
A.fleet = A.fleet(fltorder);

% Change row/column indices of pedigree table to match

if ~isempty(A.pedigree)

    [~,rowgrp] = ismember(A.pedigree.row, grporder);
    [~,rowflt] = ismember(A.pedigree.row, fltorder);
    [~,colgrp] = ismember(A.pedigree.column, grporder);
    [~,coldet] = ismember(A.pedigree.column, detorder);
    [~,colflt] = ismember(A.pedigree.column, fltorder);

    Ped = A.pedigree;

    isgrp = strcmp(Ped.property, 'groupdata');
    Ped.row(isgrp) = rowgrp(isgrp);

    isdc = strcmp(Ped.property, 'dc');
    Ped.row(isdc)    = rowgrp(isdc);
    Ped.column(isdc) = colgrp(isdc);

    iscatch = ismember(Ped.property, {'landing', 'discard'});
    Ped.row(iscatch)    = rowgrp(iscatch);
    Ped.column(iscatch) = colflt(iscatch);

    isdf = strcmp(Ped.property, 'df');
    Ped.row(isdf)    = rowgrp(isdf);
    Ped.column(isdf) = coldet(isdf);

    isdis = strcmp(Ped.property, 'discardFate');
    Ped.row(isdis)    = rowflt(isdis);
    Ped.column(isdis) = coldet(isdis);

    A.pedigree = Ped;
end



