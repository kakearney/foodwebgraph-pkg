function [Ewein, A] = mdb2ewein(file)
%MDB2EWEIN Convert EwE6 files to Ecopath input structures
%
% Ewein = mdb2ewein(file)
% [Ewein, A] = mdb2ewein(file)
%
% This function reads data from EwE6 database files and converts it into
% the structure format required by ecopathlite.  It is only designed to
% work with EwE version 6 files; older files must be converted via the
% conversion utility that comes with EwE6.
%
% This function relies on the mdbtools utilities to read the MS Access
% database files (http://mdbtools.sourceforge.net/); these work on
% Linux/Unix and Mac systems, but I have no idea if they can be installed
% on a Windows operating system.
%
% Input variables:
%
%   file:   name of Ecopath file.  Ecopath 6 uses Microsoft Access 2003
%           database files to store its data; by default they have .ewemdb
%           extensions
%
% Output variables:
%
%   Ewein:  ecopathlite input structure (see ecopathlite.m).  Also includes
%           a few extra fields just for reference:
%
%           name:   ngroup x 1 cell array of strings, names of each
%                   functional group
%
%           Info:   Structure with Ecopath model information (e.g. model
%                   name, description, units, etc.).
%
%   A:      structure of database arrays holding all data from the file in
%           its original format

% Copyright 2012-2014 Kelly Kearney

%----------------------------
% Extract data from .mdb file
%----------------------------

if ~exist(file, 'file')
    error('Could not find file: %s', file);
end

% Translate file names such that they can be passed via command line

file = regexprep(file, '([\s,<>|:\(\)&;\?\*])', '\\$1');

% Read table names

cmd = sprintf('mdb-tables -d, %s', file);
[s,r] = system(cmd);
if s
    error('Error reading table names: %s', r);
end
    
tables = regexp(r, ',', 'split');
tables = regexprep(tables, '\n', '');
isemp = cellfun('isempty', tables);
tables = tables(~isemp);    

% isecopath = strncmp(tables, 'Ecopath', 7);
% tables = tables(isecopath);

% Crazy regular expression for Excel-style comma-delimited stuff... find
% commas that aren't imbedded within quotes

pattern = ',(?=(?:[^\"]*\"[^\"]*\")*(?![^\"]*\"))';

% Read all Ecopath table data

err = cell(0,2);
for it = 1:length(tables)
    cmd = sprintf('mdb-export %s %s > mdbtmp.txt', file, tables{it});
    [s,r] = system(cmd);
    if s
        err = [err; {tables{it} r}];
        continue
    end
    
    try
        A.(tables{it}) = dataset('File', 'mdbtmp.txt', 'Delimiter', ',');  
        
        vname = get(A.(tables{it}), 'VarNames');
        for iv = 1:length(vname)
            if isnumeric(A.(tables{it}).(vname{iv}))
                A.(tables{it}).(vname{iv})(A.(tables{it}).(vname{iv}) == -9999) = NaN;
            else
                A.(tables{it}).(vname{iv}) = regexprep(A.(tables{it}).(vname{iv}), {'(^")|("$)', '""'}, {'', '"'}); % Strip out beginning/end quotes and escaped quotes
            end
        end
        
    catch
        
        try
            % Dataset doesn't like Access-style commas-in-quotes fields, so
            % parse those manually

            [data, s] = readtext('mdbtmp.txt', ',', '', '"');
            data(s.stringMask) = regexprep(data(s.stringMask), {'(^")|("$)', '""'}, {'', '"'});
            
            tmp = cell2mat(data(s.numberMask));
            tmp(tmp == -9999) = NaN;
            data(s.numberMask) = num2cell(tmp);
            
            A.(tables{it}) = cell2dataset(data);
            
        catch
            
            % if still no good, just return the text (and warn if the user
            % is getting this info)
            
            if nargout > 1
                warning('MDBEWEIN:parse', 'Could not parse table: %s', tables{it});
            end
            A.(tables{it}) = r;
            
        end
    end
end

%----------------------------
% Build Ewein structure
%----------------------------

% Reorder data if necessary (had to make things difficult, didn't you?).

[A.EcopathGroup, isrtg] = sortrows(A.EcopathGroup, 'Sequence');
[A.EcopathFleet, isrtf] = sortrows(A.EcopathFleet, 'Sequence');

% Group, living group, and gear/fleet number

Ewein.ngroup = size(A.EcopathGroup,1);
Ewein.nlive = sum(A.EcopathGroup.Type ~= 2);
Ewein.ngear  = size(A.EcopathFleet,1);

% Basic info

Ewein.areafrac = A.EcopathGroup.Area;
Ewein.b        = A.EcopathGroup.Biomass;
Ewein.pb       = A.EcopathGroup.ProdBiom;
Ewein.qb       = A.EcopathGroup.ConsBiom;
Ewein.ee       = A.EcopathGroup.EcoEfficiency;
Ewein.ge       = A.EcopathGroup.ProdCons;
Ewein.gs       = A.EcopathGroup.Unassim;
Ewein.dtImp    = A.EcopathGroup.DtImports;
Ewein.bh       = nan(Ewein.ngroup,1); % Not in file?
Ewein.pp       = A.EcopathGroup.Type;

% Diet


[tf, seqi] = ismember(A.EcopathDietComp.PreyID, A.EcopathGroup.GroupID);
[tf, seqj] = ismember(A.EcopathDietComp.PredID, A.EcopathGroup.GroupID);

Ewein.dc = full(sparse(seqi, seqj, A.EcopathDietComp.Diet, Ewein.ngroup, Ewein.ngroup));

% Detritus fate

detidx = A.EcopathGroup.Sequence(A.EcopathGroup.Type == 2);

Ewein.df = full(sparse(seqi, seqj, A.EcopathDietComp.DetritusFate, Ewein.ngroup, Ewein.ngroup));
Ewein.df = Ewein.df(detidx,:)'; % Why is this backwards?... so confusing

% Immigration/emigration

Ewein.immig    = A.EcopathGroup.Immigration;
Ewein.emig     = A.EcopathGroup.Emigration;
Ewein.emigRate = A.EcopathGroup.EmigRate;

% Biomass accumulation

Ewein.ba       = A.EcopathGroup.BiomAcc;
Ewein.baRate   = A.EcopathGroup.BiomAccRate;

% Landings and discards

[tf, seqci] = ismember(A.EcopathCatch.GroupID, A.EcopathGroup.GroupID);
[tf, seqcf] = ismember(A.EcopathCatch.FleetID, A.EcopathFleet.FleetID);

Ewein.landing = full(sparse(seqci, seqcf, A.EcopathCatch.Landing, Ewein.ngroup, Ewein.ngear));
Ewein.discard = full(sparse(seqci, seqcf, A.EcopathCatch.Discards, Ewein.ngroup, Ewein.ngear));

% Discard fate

[tf, seqdi] = ismember(A.EcopathDiscardFate.GroupID, A.EcopathGroup.GroupID);
[tf, seqdf] = ismember(A.EcopathDiscardFate.FleetID, A.EcopathFleet.FleetID);

Ewein.discardFate = full(sparse(seqdi, seqdf, A.EcopathDiscardFate.DiscardFate, Ewein.ngroup, Ewein.ngear));
Ewein.discardFate = Ewein.discardFate(detidx,:)';

% Names and Model information

Ewein.name = A.EcopathGroup.GroupName;
Ewein.fleet = A.EcopathFleet.FleetName;

Ewein.Info = A.EcopathModel;

% Multi-stanza details

[tf, seqs] = ismember(A.StanzaLifeStage.GroupID, A.EcopathGroup.GroupID);

Ewein.stanza = zeros(Ewein.ngroup,1);
Ewein.stanza(seqs) = A.StanzaLifeStage.StanzaID;

Ewein.ageStart = zeros(Ewein.ngroup,1);
Ewein.ageStart(seqs) = A.StanzaLifeStage.AgeStart;

Ewein.stanzadata = A.Stanza;
Ewein.vbK = A.EcopathGroup.vbK;





        
