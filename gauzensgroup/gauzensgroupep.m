function tgroup = gauzensgroupep(G, nrep, varargin)
%GAUZENSGROUPMULTI Run the Gauzens grouping algorithm on an Ecopath model
%
% tgroup = gauzensgroupep(G, nrep, p1, v1)
%
% This function applies the Gauzens et al. modularity or trophic grouping
% algorithm to an Ecopath model.  For the most part, the grouping algorithm
% is the same as in gauzensgroup.m, with a few modifications:
%
% - For trophic grouping, certain types of nodes are separated from each
%   other, even if the grouping metric places them together.  The metric
%   tends to put fishing fleets and detrital groups together; this is
%   prevented here.
%
% - The trophic group indices (for the first pass) are sorted by trophic
%   level.  In the original algorithm, the actual numbers assigned to
%   trophic groups are irrelevant; the key feature of the output is simply
%   seeing which input nodes are grouped together.  Renumbering so the
%   indices are sequential and continuous (numbers can be skipped in the
%   original) makes for slightly easier bookkeeping.
%
% - Multiple iterations of the algorithm can be performed.  This works by
%   running the grouping once, then creating a new food web consisting of
%   the consolidated functional groups and running the algorithm again.
%   This process is a bit experimental, but might be a useful way to build
%   a dendrogram-like hierarchy of functional groups when starting with a
%   very species-specific food web.
%
% Input variables:
%
%   G:      ecopath graph object (see ecopath2graph)
%
%   nrep:   number of grouping repetitions
%
% Optional input varialbes:
%
%   See gauzensgroup.m, all optional parameters are the same.
%
% Output variables:
%
%   tgroup: 1 x nrep cell array, holding group indices calculated at each
%           iteration

tgroup = cell(1, nrep);

Opt.tini = 0.9995;      % Initial temperature
Opt.tfin = 1e-13;       % Final temperature
Opt.ts = 0.95;          % Cooling rate
Opt.type = 'modularity';
Opt.plot = false;
Opt.limit = 10;
Opt.etol = 1e-5;
Opt.split = 'random';

Opt = parsepv(Opt, varargin);
p = fieldnames(Opt);
v = struct2cell(Opt);
pv = [p v]';

for it = 1:nrep
    if it == 1
        [tgroup{it}, Tg] = gauzensgroup(adjacency(G), pv{:});
        
        if strcmp(Opt.type, 'trophicgroup') && ismember('type', G.Nodes.Properties.VariableNames)
            
            % Break up groups that include certain type combos
            
            [tgtmp, typ] = aggregate(tgroup{it}, G.Nodes.type);
            for ii = 1:length(tgtmp)
                isdetfish = all(ismember([0 3], typ{ii})); % Detritus with fleet
                if isdetfish
                    tgnew = max(tgroup{it}) + 1;
                    change = G.Nodes.type == 3 & tgroup{1} == tgtmp(ii);
                    tgroup{1}(change) = tgnew;
                end                    
            end
        end
        
        if ismember('TL', G.Nodes.Properties.VariableNames)
        
            % Grouping indices are arbitrary; in many cases, it's useful to add
            % some order.  For now, renumbering so groups run from lowest
            % average trophic level to highest

            [gidx, tl] = aggregate(tgroup{it}, G.Nodes.TL, @nanmean);
            [~, isrt] = sort(cat(1, tl{:}));

            [~, tgroup{it}] = ismember(tgroup{it}, gidx(isrt));
        end
        
    else
        [tgroup{it}, Tg] = gauzensgroup(adj, pv{:});
    end
    
    % Set up adjacency for next repetition
    
    [ii,jj] = find(Tg.adj);
    ij = unique([tgroup{it}(ii) tgroup{it}(jj)], 'rows');
    ng = max(tgroup{it});
    adj = full(sparse(ij(:,1), ij(:,2), ones(size(ij,1),1), ng, ng));
end

