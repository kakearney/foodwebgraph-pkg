function [tgpart, Out, G] = trophicgroup(G, varargin)
%TROPHICGROUP Partition the nodes in a food web into trophic groups
%
% tgpart = trophicgroup(G)
% tgpart = trophicgroup(EM)
% [tgpart, Out, G] = trophicgroup(..., p1, v1, ...)
%
% Based on Gauzens et al. 2014, this function attempts to partition a food
% web network graph such that the resulting node clusters maximize the
% trophic similarity index (notated as G(E) in the Gauzens et al. paper,
% defined in Eq 2.2) of the graph. 
%
%   G:          a digraph object, where nodes represent functional groups
%               in the ecosystem and edges reresent predator-prey links.
%   
%   EM:         an ecopathmodel object.  Nodes in the graph will consist of
%               all groups and fleets in the model, and edges will be all
%               predator-prey and fisheries catch links.
%
% Optional input variables (passed as parameter/value pairs):
%
%   method:     algorithm to use when calculating optimal partition
%               'hierarchical': Hierarchical clustering, using 1 - trophic
%                   similarity as the distance metric.  Optimal cutoff
%                   value is chosen based on the one that maximizes the
%                   partition trophic similarity index.  This method is
%                   fast, and usually gets to a decent partition, with cost
%                   ~90% of the optimal.
%               'greedy': A greedy algorithm that begins with all nodes
%                   in their own clusters, then iteratively combines the
%                   two groups that result in the highest-scoring
%                   partition. This method is very fast but usually results
%                   in a partition with cost value much lower than the
%                   optimal one, ~30% of optimal.   
%               'simulatedannealing': The simulated annealing option
%                   used by Gauzen et al.  This method is the slowest, but
%                   gets closest to the optimal value, usually within 99%.
%               ['hierarchical']
%
%   splittype:  logical scalar, true to split apart any clusters that
%               combine nodes of different types (as defined by the graph's
%               type node property; for an ecopathmodel object, this refers
%               to the producer/consumer/detrital/fleet designation).  Note
%               that this splitting is done after the optimization is
%               complete. [true for ecopathmodel input, false for digraph
%               input]
%
%   sort:       logical scalar, true to sort resulting clusters from lowest
%               to highest average trophic level.  The cluster indices
%               returned by all algorithms are arbitrary, so it can be
%               useful to put them into some sort of order [true for
%               ecopathmodel input, false for digraph input]
%
%   (For hierachical method only)
%
%   linkmethod: Algorithm used to create the hierachical cluster tree.
%               'single':   nearest distance (default)
%               'complete': furthest distance
%               'average':  unweighted average distance (UPGMA) (also known as group average)
%               'weighted': weighted average distance (WPGMA)
%               'centroid': unweighted center of mass distance (UPGMC)
%               'median':   weighted center of mass distance (WPGMC)
%               'ward':     inner squared distance (min variance algorithm)
%               ['weighted'] 
%
%   (For simulated annealing method only)
%
%   tini:       Initial temperature for simulated annealing [0.9995]
%
%   tfin:       Final temperature for simulated annealing [1e-13]
%
%   ts:         Cooling rate [0.95]
%
%   limit:      The simulated annealing algorithm is terminated when the
%               change in the cost function is less than a specified error
%               tolerance  over this many consecutive temperature time
%               steps. [10]     
%
%   etol:       error tolerance used for simulated annealing stopping
%               criteria [1e-5] 
%
%   split:      method used to split a cluster into two potential new
%               clusters 
%               'random':   split at random (highly recommended)
%               'optimal':  split to maximize modularity or trophic
%                           similarity, via additional simulated annealing
%                           (slow)  
%               'optimal2': split to maximize modularity or trophic
%                           similarity, via brute-force checking all
%                           potential split-partitions (slow).  
%               ['random']
%
%   verbose:    Print progress to screen as simulation proceeds.  Includes
%               the temperature (t, from 1 to 0), the change in energy
%               metric over the last set of moves (dE), and the current
%               energy value (E). [true]   
%
% Output variables:
%
%   tgpart:     nnode x 1 array, cluster indices corresponding to each node
%               in the final partitioning
%
%   Out:        structure with some additional info:
%
%               tgmetricMax:        trophic group index value achieved
%                                   by the returned partitioning
%               tgmetricDetails:    structure with additional data, fields
%                                   vary by algorithm
%
%                   (Hierarchical)
%
%                   costAtCutoff:   metric at each potential cutoff value
%                   cutoff:         cutoff values
%                   partitions:     partitions resulting at each cutoff
%                                   point
%
%                   (Greedy)
%
%                   costAtIteration:metric after each iteration
%
%                   (Simulated Annealing)
%
%                   costAtMove:     metric from each proposed move
%                   costAtAcceptedMove: metric after move has been accepted
%                                   or rejected
%
%   G:          graph object used for calculations.

% Copyright 2016 Kelly Kearney

%--------------------
% Parse input
%--------------------

validateattributes(G, {'digraph', 'ecopathmodel'}, {});
if isa(G, 'ecopathmodel')
    G = G.graph('det', false, 'oos', false);
    splitdefault = true;
    sortdefault = true;
else
    splitdefault = false;
    sortdefault = false;
end 

% Overall parameters

p = inputParser;
p.KeepUnmatched = true;
p.addParameter('method',    'hierarchical', @(x) validateattributes(x, {'char'}, {}));
p.addParameter('splittype', splitdefault,   @(x) validateattributes(x, {'logical'}, {'scalar'}));
p.addParameter('sort',      sortdefault,    @(x) validateattributes(x, {'logical'}, {'scalar'}));

p.parse(varargin{:});
Opt = p.Results;

validatestring(Opt.method, {'hierarchical', 'greedy', 'simulatedannealing'});

% Parameters for hierarchical clustering option

p1 = inputParser;
p1.KeepUnmatched = true;
p1.addParameter('linkmethod', 'weighted', @(x) validateattributes(x, {'char'}, {}));

p1.parse(varargin{:});
Hopt = p1.Results;

% Parameters for simulated annealing option:

p2 = inputParser;
p2.KeepUnmatched = true;

p2.addParameter('tini',    0.9995,       @(x) validateattributes(x, {'numeric'}, {'scalar'}));
p2.addParameter('tfin',    1e-13,        @(x) validateattributes(x, {'numeric'}, {'scalar'}));
p2.addParameter('ts',      0.95,         @(x) validateattributes(x, {'numeric'}, {'scalar'}));
p2.addParameter('limit',   10,           @(x) validateattributes(x, {'numeric'}, {'scalar'}));
p2.addParameter('etol',    1e-5,         @(x) validateattributes(x, {'numeric'}, {'scalar'}));
p2.addParameter('plot',    false,        @(x) validateattributes(x, {'logical'}, {'scalar'}));
p2.addParameter('verbose', false,        @(x) validateattributes(x, {'logical'}, {'scalar'}));
p2.addParameter('split',   'random',     @(x) validateattributes(x, {'char'}, {}));


p2.parse(varargin{:});
SAopt = p2.Results;

% Parameters for greedy option:

%--------------------
% Setup calculations
% for all methods
%--------------------

% Basic graph properties

nnode = numnodes(G);
adj = adjacency(G);

nedge = nnz(adj);                      % # links in the food web
degree = full(sum(adj>0,1)' + sum(adj>0,2));  % degree of each node

% Connected components

cidx = conncomp(graph(adj | adj'));

% Trophic similarity and expected trophic similarity

[tsim, etsim] = trophicsimilarity(adj);
[ir,ic] = meshgrid(1:nnode);

% Difference between tsim and etsim.  Zero out elements where ir >= ic, to
% reduce calculation times later.

tminuse = triu(tsim - etsim, 1);

% Cost function: return cost value for a given partition.  This is the
% value to be maximized.

costfun = @(x) costmetric(x, tminuse, ir);

%---------------------
% Calculate partitions
%---------------------

switch Opt.method
    
    case 'hierarchical'
        
        % Build linkage matrix based on inverse of trophic similarity

        dists = squareform(1 - tsim, 'tovector');
        links = linkage(dists, Hopt.linkmethod);

        nlink = size(links,1);
        
        % Identify clusters for each possible partitioning (i.e. at cutoff
        % values just above each link in the dendrogram).
        
        idx = cell(nlink,1);

        for ii = 1:nlink
            if links(ii,1) <= nnode
                tmp1 = links(ii,1);
            else
                tmp1 = idx{links(ii,1) - nnode};
            end
            if links(ii,2) <= nnode
                tmp2 = links(ii,2);
            else
                tmp2 = idx{links(ii,2) - nnode};
            end
            idx{ii} = [tmp1 tmp2];
        end
        
        % Create lists of cluster indices for each partition

        linkdist = unique(links(:,3));
        cutoff = mean([linkdist(1:end-1) linkdist(2:end)], 2);
        cutoff = [-0.1; cutoff];
        npart = length(cutoff);

        pidx = zeros(nnode,npart);

        for iw = 1:npart
            if iw == 1
                pidx(:,iw) = 1:nnode;
            else
                pidx(:,iw) = cluster(links, 'cutoff', cutoff(iw), 'criterion', 'distance'); 
            end
        end

        % Renumber so identify which cluster the indices refer to, keeping
        % things consistent across partitions (not entirely necessary, but
        % it makes bookkeeping easier)  

        tg = zeros(size(pidx)); 
        for iw = 1:npart
            [xcon, yagg] = aggregate(pidx(:,iw), (1:nnode)');
            for ii = 1:length(xcon)
                if isscalar(yagg{ii})
                    tg(pidx(:,iw)==xcon(ii),iw) = yagg{ii};
                else
                    isc = cellfun(@(x) isequal(sort(yagg{ii}'), sort(x)), idx);
                    tg(pidx(:,iw)==xcon(ii),iw) = find(isc) + 47;
                end
            end
        end

        % Calculate energy for each potential partition

        energy = zeros(npart,1);
        for ii = 1:npart
            energy(ii) = costfun(tg(:,ii));
        end
        
        % Choose the partition that minimizes the trophic group index
        
        [~,imax] = max(energy);
        
        tgpart = pidx(:,imax);
        
        Out.tgmetricMax = max(energy);
        Out.tgmetricDetails = struct('costAtCutoff', energy, 'cutoff', ...
            cutoff, 'partitions', tg, 'links', links);
        
    case 'greedy'
        
        % Start with all nodes in their own clusters
        
        pidx = (1:nnode)';
        
        energy = costfun(pidx);
        
        % Iteratively combine clusters
        
        etrack = nan(nnode-1,1);
        
        for irep = 1:nnode-1
            
            % Calculate energy for all potential combinations of existing
            % clusters
            
            ncluster = max(pidx);
            
            energynew = zeros(ncluster);
            for ii = 1:ncluster
                for jj = (ii+1):ncluster
                    tmp = pidx;
                    tmp(pidx == jj) = pidx(ii);
                    energynew(ii,jj) = costfun(tmp);
                end
            end
            
            % Choose the combo that maximizes energy.  If none of them
            % raise it, we're done.
            
            delta = triu(energynew - energy,1);
            if all(delta <= 0)
                break
            end
            
            [~,imax] = max(delta(:));
            [c1,c2] = ind2sub([ncluster ncluster], imax);
            pidx(pidx == c2) = c1;
            
            pidx(pidx > c2) = pidx(pidx > c2) - 1; % Fill the gap
            
            etrack(irep) = energynew(imax);
            energy = energynew(imax);
            
        end
        
        tgpart = pidx;
        
        Out.tgmetricMax = max(etrack);
        Out.tgmetricDetails = struct('costAtIteration', etrack);
        
        
    case 'simulatedannealing'
        
        [tgpart, etrack] = simanpartition(costfun, nnode, cidx, SAopt);
        
        Out.tgmetricMax = etrack(end,2); %max(etrack(:,2));
        Out.tgmetricDetails = struct('costAtMove', etrack(:,1), 'costAtAcceptedMove', etrack(:,2));
        
end

%---------------------
% Ecopath-specific 
% modifications
%---------------------   

% Disallow grouping across types

if Opt.splittype
    if ~ismember('type', G.Nodes.Properties.VariableNames)
        warning('splittype flag ignored because graph does not include type as a node property');
    else
        [~,~,tgpart] = unique([tgpart G.Nodes.type], 'rows');
    end
end

% Grouping indices are arbitrary, but it can be useful to have some sort of
% order to the groups.  Sort by trophic level.

if Opt.sort
    if ~ismember('TL', G.Nodes.Properties.VariableNames)
        warning('sort flag ignored because graph does not include TL as a node property');
    else
        [gidx, tl] = aggregate(tgpart, G.Nodes.TL, @nanmean);
        [~, isrt] = sort(cat(1, tl{:}));

        [~, tgpart] = ismember(tgpart, gidx(isrt));
    end
end



    
%---------------------
% Trophic group index
%---------------------    

% This function calculates the trophic group index for a given partitioning
% of the graph.  The tminuse matrix should be the upper triangular portion
% of the tsim - etsim matrix.

function energy = costmetric(x, tminuse, ir)

insamegrp = bsxfun(@eq, x, x');

if any(insamegrp(:))
    n = histc(x, 1:max(x));
    grpnum = x(ir(insamegrp));
    
    dt = accumarray(grpnum, tminuse(insamegrp), [max(x) 1]);
    energy = nansum(dt./n);
    
else
    energy = 0;
end






%*********






% %--------------------
% % Viewer
% %--------------------
% 
% if Opt.plot
% 
%     elim = [0 max(energy)*1.1];
% 
%     plotdendro = @(c) dendrogram(links, 0, 'orientation', 'right', ...
%         'colorthreshold', c, 'labels', G.Nodes.Name);
% 
%     h = plotgrid('setup', cell(1), [],[],'mar', 0.05, 'ml', 0.2);
%     h.ax = subgrid(h.ax, [0.1 0.9], 1);
%     h.ax = cat(1, h.ax{:});
% 
%     axes(h.ax(2));
%     h.stm = stem(cutoff, energy, '.');
%     h.ref = line(nan(1,2), elim, 'color', 'r');
% 
%     h.stm.ButtonDownFcn = {@stemclick, h.ref, plotdendro, h.ax(1)}
% 
%     axes(h.ax(1));
%     plotdendro(Inf);
% 
%     xlim = minmax(h.ax(1).XLim, 'expand', 0.02);
% 
%     set(h.ax, 'xlim', xlim);
% 
%     set(h.ax(1), 'xaxisloc', 'top', 'box', 'off', 'TickLabelInterpreter', 'none');
%     set(h.ax(2), 'box', 'off');
% 
% end
% 
% 
% %--------------------
% % Subfunctions
% %--------------------
% 
% % Calculate trophic similarity index (G(E), i.e. Eq 2.2 from Gauzen et al.)
% 
% function energy = computetg(x, tsim, etsim, ir, ic)
% 
% insamegrp = x(ir) == x(ic) & ic > ir;
% 
% if any(insamegrp(:))
%     n = histc(x, 1:max(x));
%     grpnum = x(ir(insamegrp));
%     tsimgrp = tsim(insamegrp);
%     etsimgrp = etsim(insamegrp);
%     tt = accumarray(grpnum, tsimgrp, [max(x) 1]);
%     et = accumarray(grpnum, etsimgrp, [max(x) 1]);
%     energy = nansum((tt - et)./n);
% else
%     energy = 0;
% end
% 
% function stemclick(h,ed, href, dfun, dax)
% 
% xyz = get(ancestor(h, 'axes'), 'CurrentPoint');
% [dx, imin] = min(abs(xyz(1) - h.XData));
% 
% thresh = h.XData(imin);
% 
% set(href, 'xdata', ones(1,2)*thresh);
% 
% props = {'xlim' 'xaxisloc', 'box', 'TickLabelInterpreter'};
% vals = get(dax, props);
% pv = [props; vals];
% 
% cla(dax);
% axes(dax);
% dfun(thresh);
% 
% set(dax, pv{:});
% title(dax, sprintf('Cutoff = %.2f (index %d)', thresh, imin));





