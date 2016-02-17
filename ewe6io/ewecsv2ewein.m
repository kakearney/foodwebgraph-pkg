function S = ewecsv2ewein(file)
%EWECSV2EWEIN Converts a file with EwE model data to an EwE input structure
%
% S = ewecsv2ewein(file)
%
% This function converts cut-and-pasted data from the EwE 6 GUI into an EwE
% input structure.  The model data should be saved in a comma-delimited
% file.  Please see eweModelData.xls for the format of this file.  The EwE
% input structure holds all of this data in a Matlab structure, and is the
% required input format for my EwE Matlab functions.
%
% Note: This function was designed for EwE version 5, and updated a bit for
% the early releases of Ewe6. But it hasn't been tested against recent
% versions. 
%
% Input variables:
%
%   file:       filename of comma-delimited file
%
% Output variables:
%
%   S:          EwE input structure, with the following fields:
%
%               ngroup:         1 x 1 array, number of functional groups in
%                               the model  
%
%               nlive:          1 x 1 array, number of live (non-detrital)
%                               groups in the model
%
%               ngear:          1 x 1, array number of fishing gear types
%                               in the model 
%                          
%               name:           ngroup x 1 cell array, names of groups 
%
%               pp:             ngroup x 1 array, fraction of diet
%                               consisting of primary production, pp = 2
%                               indicates detritus 
%
%               areafrac:       ngroup x 1 array, fraction of habitat area
%                               occupied by each group  
%
%               b:              ngroup x 1 array, biomass, -9999 indicates
%                               unknown value 
%
%               pb:             ngroup x 1 array, produdction over biomass
%                               ratios, -9999 indicates unknown value 
%   
%               qb:             ngroup x 1 array, consumption over biomass
%                               ratios, -9999 indicates unknown value
%
%               ee:             ngroup x 1 array, ecotrophic efficiencies,
%                               -9999 indicates unknown value
% 
%               ge:             ngroup x 1 array, gross efficiency, i.e.
%                               production over consumption ratio, -9999
%                               indicates unknown value
%
%               gs:             ngroup x 1 array, fraction of consumed food
%                               that is not assimilated
%
%               dtImp:          ngroup x 1 array, detritus import (must be
%                               zero for all non-detrital groups)  
%
%               bh:             ngroup x 1 array,  habitat biomass, i.e.
%                               biomass per unit area (b/areafrac) 
%
% 
%               dc:             ngroup x ngroup array, diet composition,
%                               dc(i,j) tells fraction predator j's diet
%                               consisting of prey i   
% 
%               df:             ngroup x (ngroup - nlive) array, fraction
%                               of each group that goes to each detrital
%                               group due to other mortality and egestion  
%
%               immig:          ngroup x 1 array, immigration into area
%
%               emig:           ngroup x 1 array, emigration out of area
%
%               emigRate:       ngroup x 1 array, emigration per unit
%                               biomass 
% 
%               ba:             ngroup x 1 array, biomass accumulation
%
%               baRate:         ngroup x 1 array, biomass accumulation per
%                               unit biomass 
%
%               landing:        ngroup x ngear array, landings of each
%                               group by each gear type 
%
%               discard:        ngroup x ngear array, discards of each
%                               group by each gear type 
%
%               discardFate:    ngear x (ngroup - nlive) array, fraction of
%                               discards from each gear type that go to
%                               each detritus group  
%
%               maxrelpb:       nlive x 1 array, maximum p/b ratio,
%                               relative to initial p/b, possible if prey
%                               is abundant (prevents unrealistically high
%                               consumption)
%               
%               maxrelfeed:     nlive x 1 array, maximum amount of feeding
%                               time, relative to initial amount, a
%                               predator will spend searching for prey
%                               (prevents unrealistically high time when
%                               prey is scarce)
%
%               feedadj:        nlive x 1 array, factor determining how
%                               fast organisms adjust their feeding time, 0
%                               = remains constant, 1 = adjusts very
%                               quickly to food changes   
%
%               fracsens:       nlive x 1 array, fraction of other
%                               mortality that is sensitive to feeding time 
%
%               predeffect:     nlive x 1 array, allows for direct response
%                               of feeding time to predator abundance  
%
%               densecatch:     nlive x 1 array, allows for changes in
%                               fishing rate based on biomass  
%
%               qbmaxqb0:       nlive x 1 array, something to do with
%                               handling time
%
%               switchpower:    nlive x 1 array, ?

% Copyright 2007 Kelly Kearney

[data, Result] = readtext(file, ',', '', '"');

% Groups and gears

S.ngroup = data{1,2};
S.nlive  = data{2,2};
S.ngear  = data{3,2};

% Group names

irow = find(strcmp('"Group numbers, names, and primary producer indicators"', data(:,1)));
S.name = data(irow+1:irow+S.ngroup,2);
for iname = 1:length(S.name)
    if regexp(S.name{iname}, '^"') & regexp(S.name{iname}, '"$')
        S.name{iname} = regexprep(S.name{iname}, '^"' ,'');
        S.name{iname} = regexprep(S.name{iname}, '"$', '');
    end
end

S.pp = cell2matfilled(data(irow+1:irow+S.ngroup,3), 0);

% Basic input

irow = find(strcmp('Basic input', data(:,1)));
basic = cell2matfilled(data(irow+1:irow+S.ngroup,3:10), NaN);

S.areafrac = basic(:,1);
S.bh       = basic(:,2);
S.pb       = basic(:,3);
S.qb       = basic(:,4);
S.ee       = basic(:,5);
S.ge       = basic(:,6);
S.gs       = basic(:,7);
S.dtImp    = basic(:,8);

S.b = S.bh .* S.areafrac;

% temp = S.areafrac > 0 & S.b > 0;
% S.bh = zeros(S.ngroup,1);
% S.bh(temp) = S.b(temp) ./ S.areafrac(temp);

% Diet composition

nconsumer = length(find(S.pp < 1));
irow = find(strcmp('Diet composition', data(:,1)));
dc = cell2matfilled(data(irow+1:irow+S.ngroup,3:nconsumer+2), 0);
S.dc = zeros(S.ngroup);
S.dc(:, ~S.pp) = dc(1:S.ngroup,:);

% Detritus fate

ndet = S.ngroup - S.nlive;
irow = find(strcmp('Detritus fate', data(:,1)));
S.df = cell2matfilled(data(irow+1:irow+S.ngroup,3:ndet+2), 0);

% Other production

irow = find(strcmp('Other production', data(:,1)));
other = cell2matfilled(data(irow+1:irow+S.ngroup,3:7), 0);
S.immig    = other(:,1);
S.emig     = other(:,2);
S.emigRate = other(:,3);
S.ba       = other(:,4);
S.baRate   = other(:,5);

% Definition of fleets

if S.ngear > 0
    irow = find(strcmp('Definition of fleets', data(:,1)));
    S.fleetname = data(irow+1:irow+S.ngear,2);
else
    S.fleetname = cell(0,1);
end

% Landings

irow = find(strcmp('Landings', data(:,1)));
if S.ngear > 0
    S.landing = cell2mat(data(irow+1:irow+S.ngroup,3:S.ngear+2));
else
    S.landing = zeros(S.ngroup,S.ngear);
end

% Discards

irow = find(strcmp('Discards', data(:,1)));
if S.ngear > 0
    S.discard = cell2mat(data(irow+1:irow+S.ngroup,3:S.ngear+2));
else
    S.discard = zeros(S.ngroup,S.ngear);
end

% Discard Fate

irow = find(strcmp('Discard Fate', data(:,1)));
if S.ngear > 0
    S.discardFate = cell2mat(data(irow+1:irow+S.ngear,3:ndet+2));
else
    S.discardFate = zeros(S.ngear,ndet);
end

% Group Info

irow = find(strcmp('Group Info', data(:,1)));
groupinfo = cell2mat(data(irow+1:S.nlive+irow,3:10));
S.maxrelpb    = groupinfo(:,1);
S.maxrelfeed  = groupinfo(:,2);
S.feedadj     = groupinfo(:,3);
S.fracsens    = groupinfo(:,4);
S.predeffect  = groupinfo(:,5);
S.densecatch  = groupinfo(:,6);
S.qbmaxqb0    = groupinfo(:,7);
S.switchpower = groupinfo(:,8);

% Vulnerabilities

irow = find(strcmp('Vulnerabilities', data(:,1)));
kv = cell2matfilled(data(irow+1:irow+S.ngroup,3:nconsumer+2), 2);
S.kv = zeros(S.ngroup);
S.kv(:, ~S.pp) = kv(1:S.ngroup,:);
S.kv(~S.dc) = 0;

% Subfunction: Convert cell array to matrix, filling in empty cells

function mat = cell2matfilled(cell, empty)

mat = zeros(size(cell));
for icell = 1:numel(cell)
    if isempty(cell{icell})
        mat(icell) = empty;
    else
        mat(icell) = cell{icell};
    end
end



