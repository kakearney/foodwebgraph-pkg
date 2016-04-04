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
% database files (https://github.com/brianb/mdbtools).  Mac users can get
% it via either MacPorts or Homebrew if they don't want to compile from
% source. Please make sure this utility is properly compiled prior to
% calling mdb2ewein.m.  Windows users: you will need a C compiler (I
% believe you can get one for free via Visual Studio Express).
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
%           name:       ngroup x 1 cell array of strings, names of each
%                       functional group
%
%           Info:       dataset array with Ecopath model information (e.g.
%                       model name, description, units, etc.). 
%
%           pedigree:   ngroup x 7 pedigree array.  Columns correspond to
%                       B, PB, QB, DC, EE, GE, and catch, with NaNs used as
%                       placeholders where no pedigree data is defined.
%
%   A:      structure of database arrays holding all data from the file in
%           its original format

% Copyright 2012-2015 Kelly Kearney

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

% Crazy regular expression for Excel-style comma-delimited stuff... find
% commas that aren't imbedded within quotes

pattern = ',(?=(?:[^\"]*\"[^\"]*\")*(?![^\"]*\"))';

% Read all Ecopath table data

err = cell(0,2);
for it = 1:length(tables)
    if regexpfound(tables{it}, '\s')
        tbl = regexprep(tables{it}, '([\s,<>|:\(\)&;\?\*])', '\\$1');
        cmd = sprintf('mdb-export %s %s > mdbtmp.txt', file, tbl);
        tables{it} = regexprep(tables{it}, '\s', '_');
    else
        cmd = sprintf('mdb-export %s %s > mdbtmp.txt', file, tables{it});
    end
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

% Imports

Ewein.import = A.EcopathGroup.ImpVar;

% Diet

[tf, seqi] = ismember(A.EcopathDietComp.PreyID, A.EcopathGroup.GroupID);
[tf, seqj] = ismember(A.EcopathDietComp.PredID, A.EcopathGroup.GroupID);

Ewein.dc = full(sparse(seqi, seqj, A.EcopathDietComp.Diet, Ewein.ngroup, Ewein.ngroup));
Ewein.dc(:,Ewein.pp>1) = 0; % Some files seem to store flow to det fractions here

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

% Pedigree (At the moment, EwE6 allows pedigree values for B, PB, QB, DC,
% and catch, while I focus on B, PB, QB, DC, EE, and GE)

Ewein.pedigree = nan(Ewein.ngroup, 7);

if ~isempty(A.EcopathGroupPedigree)

    cols = {'BiomassAreaInput', 'PBInput', 'QBInput', 'DietComp', 'ee','ge','TCatchInput', 'Biomass'};

    [tf, ridx] = ismember(A.EcopathGroupPedigree.GroupID, A.EcopathGroup.GroupID);
    [tf, cidx] = ismember(A.EcopathGroupPedigree.VarName, cols);
    [tf, lidx] = ismember(A.EcopathGroupPedigree.LevelID, A.Pedigree.LevelID);

    cidx(cidx == 8) = 1; % TODO: Is this the same?  Found in Albatross Bay model, which was upgraded from older version of EwE
    
    idx = sub2ind(size(Ewein.pedigree), ridx, cidx);
    Ewein.pedigree(idx) = A.Pedigree.Confidence(lidx)./100;
end





        
