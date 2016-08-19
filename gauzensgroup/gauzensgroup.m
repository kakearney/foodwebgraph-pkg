function [grp, A] = gauzensgroup(adj, varargin)
%GAUZENSGROUP Trophic groups/modularity calculation
% 
% [grp, A] = gauzensgroup(adj, p1, v1, ...)
%
% This function applies the trophic group and maximum modularity algorithms
% described in Gauzens et al. 2014 (J. R. Soc. Interface, volume 12).  It
% is based on the code developed by Gauzens, with a few minor tweaks for my
% own applications.
%
% Input variables:
%
%   adj:    n x n predator-prey adjacency matrix, where rows
%           indicate the prey group and columns the predator group, and a
%           non-zero value indicates predation of the prey by the predator.
%
% Optional input variables, passed as parameter/value pairs:
%
%   tini:   Initial temperature for simulated annealing [0.9995]
%
%   tfin:   Final temperature for simulated annealing [1e-13]
%
%   ts:     Cooling rate [0.95]
%
%   type:   Grouping algorithm to apply ['modularity']
%           'modularity':   group based on maximum modularity
%           'trophicgroup': group based on trophic simlarity
%
%   limit:  The simulated annealing algorithm is terminated when the change
%           in the cost function is less than a specified error tolerance
%           over this many consecutive temperature time steps. [10]  
%
%   etol:   error tolerance used for simulated annealing stopping criteria
%           [1e-5]
%
%   split:  method used to split a cluster into two potential new clusters
%           ['random']
%           'random':   split at random (recommended)
%           'optimal':  split to maximize modularity or trophic similarity,
%                       via additional simulated annealing
%                       (used in original code, but I found this added to
%                       the simulation time without making a significant
%                       difference in the end results).  
%           'optimal2': split to maximize modularity or trophic similarity,
%                       via brute-force checking all potential
%                       split-partitions (experiental, to see if this might
%                       be more efficient than the previous optimal choice,
%                       but still very slow).
%
%   verbose: Print progress to screen as simulation proceeds.  Includes the
%            temperature (t, from 1 to 0), the change in energy metric over
%            the last set of moves (dE), and the current energy value (E).
%            [true]
%
% Ouptut variables:
%
%   grp:    n x 1 array, indices of trophic groups to which each
%           predator/prey group belongs
%
%   A:      1 x 1 structure
%
%           adj:    Adjacency matrix used, either from user input or
%                   calculated from Ecopath model 
%
%           h:      handle to plotted figure, if plot = true
%
%           etrack: metric used for simulated annealing, over time

% Copyright 2015 Kelly Kearney

% Options

p = inputParser;
p.addParameter('tini',    0.9995,       @(x) validateattributes(x, {'numeric'}, {'scalar'}));
p.addParameter('tfin',    1e-13,        @(x) validateattributes(x, {'numeric'}, {'scalar'}));
p.addParameter('ts',      0.95,         @(x) validateattributes(x, {'numeric'}, {'scalar'}));
p.addParameter('type',    'modularity', @(x) validateattributes(x, {'char'}, {}));
p.addParameter('limit',   10,           @(x) validateattributes(x, {'numeric'}, {'scalar'}));
p.addParameter('etol',    1e-5,         @(x) validateattributes(x, {'numeric'}, {'scalar'}));
p.addParameter('split',   'random',     @(x) validateattributes(x, {'char'}, {}));
p.addParameter('plot',    false,        @(x) validateattributes(x, {'logical'}, {'scalar'}));
p.addParameter('verbose', false,        @(x) validateattributes(x, {'logical'}, {'scalar'}));

p.parse(varargin{:});
Opt = p.Results;

Opt.type = validatestring(Opt.type, {'modularity', 'trophicgroup'}, '', 'type');
Opt.split = validatestring(Opt.split, {'random', 'optimal'}, '', 'type');

% Setup

if Opt.verbose
    fprintf('Initializing...\n');
end

nnode = size(adj,1);                    % Number of groups in original food web
nlinks = nnz(adj);                      % # links in the food web
degree = sum(adj>0,1)' + sum(adj>0,2);  % degree of each node

% Connected components

cidx = conncomp(graph(adj | adj'));

% Initial plot

% if Opt.plot
%     h.fig = figure;
%     h.ax(1) = subplot(2,1,1);
%     h.ax(2) = subplot(2,1,2);
%     if verLessThan('matlab', 'R2015b')
%         axes(h.ax(1));
%         text(0.5, 0.5, 'Graph plot requires R2015b or later', 'horiz', 'center');
%     else
%         cmap = [...
%           0.12157      0.46667      0.70588
%           0.68235      0.78039       0.9098
%                 1      0.49804     0.054902
%                 1      0.73333      0.47059
%           0.17255      0.62745      0.17255
%           0.59608      0.87451      0.54118
%           0.48627      0.15294      0.15686
%                 1      0.59608      0.58824
%           0.58039      0.40392      0.74118
%           0.77255       0.6902      0.83529
%           0.54902      0.33725      0.29412
%           0.76863      0.61176      0.58039
%            0.8902      0.46667      0.76078
%           0.96863      0.71373      0.82353
%           0.49804      0.49804      0.49804
%           0.78039      0.78039      0.78039
%           0.73725      0.74118      0.13333
%           0.85882      0.85882      0.55294
%          0.090196       0.7451          0.8
%           0.61961       0.8549      0.89804]; % d3cat20
%         
%         axes(h.ax(1));
%         
%         G = digraph(adj);
%         h.gp = plot(G, 'layout', 'layered', 'edgecolor', ones(1,3)*0.5, ...
%             'direction', 'up', 'nodecdata', cidx, 'nodecolor', 'flat');
%         colormap(cmap);
%         set(h.ax(1), 'clim', [0.5 20.5]);
%         drawnow;
%     end
% end

% Cost functions for trophic grouping or modularity

if strcmp(Opt.type, 'trophicgroup')
    
    [tsim, etsim] = trophicsimilarity(adj);
    [ir,ic] = meshgrid(1:nnode);
    
    tsimdiff = triu(tsim - etsim,1);
    
    efun = @(grp) computetg(grp, tsimdiff, ir);
    
else
    
    efun = @(grp) computemod(grp, nlinks, degree, adj);
    
end

% Run simulated annealing

SAopt = rmfield(Opt, {'type'});

[grp, etrack] = simanpartition(efun, nnode, cidx, SAopt);


% 
% 
% 
% % Initialization
% 
% grp = (1:nnode)';   % start with nodes in separate groups
% 
% switch Opt.type
%     case 'modularity'
%         efun = @(grp) computemod(grp, nlinks, degree, adj);
%     case 'trophicgroup'
%         efun = @(grp) computetg(grp, tsimdiff, ir);
% end
%         
% energyant = efun(grp);
%         
% % energyant = computeenergy(Opt.type, grp, nlinks, degree, adj);
% 
% ntmax = ceil(log(Opt.tfin/Opt.tini)./log(Opt.ts));
% nmoves = ntmax * (nnode^2 + nnode);
% movecount = 0;
% etrack = nan(nmoves,2);
% 
% % Tracking plot
% 
% if Opt.plot
% 
%     npt = 10000;
%     h.en = plot(h.ax(2), 1:npt, nan(npt,2));
%     set(h.ax(2), 'xlim', [0 npt]);
%     if strcmp(Opt.type, 'modularity')
%         set(h.ax(2), 'ylim', [-0.5 1]);
%     else
%         set(h.ax(2), 'ylim', [0 25]);
%     end
%     drawnow;
% end
% 
% t = Opt.tini;
% count = 0;
% eprev = -Inf;
% 
% if Opt.verbose
%     tmp = repmat(' ', 1, 44);
%     fprintf('Beginning simulated annealing:\n%s\n', tmp);
% end
% 
% while t > Opt.tfin && count < Opt.limit
%     
%     % Individual movements
%     
%     try
%         npergrp = accumarray(grp, grp, [nnode 1], @length);
%     catch ME
%         sstr = sprintf('%d ', grp);
%         fprintf('accumarray issue:\n  subs = %s\n  sz = %d', sstr, nnode);
%         rethrow(ME);
%     end
%     openidx = find(npergrp > 0);
%     ng = length(openidx);
%     
%     for ii = 1:nnode^2
%         srcidx = ceil(rand(1)*nnode);
%         taridx = openidx(ceil(rand(1)*ng));
%         grpn = grp;
%         grpn(srcidx) = taridx;
%         
%         energy = efun(grpn);
% %         test = efun2(grpn);
% %         if abs(test-energy)>1e-15
% %             blah = 1;
% %         end
%         
%         if energy > energyant || exp((energy - energyant)/t) > rand(1)
%             grp = grpn;
%             energyant = energy;
%         end
%         
%         movecount = movecount + 1;
%         etrack(movecount,1) = energy;
%         etrack(movecount,2) = energyant;
%     end
%     
%     % Cleanup: get rid of any now-empty groups
%     
%     [~,~,grp] = unique(grp);
%     
%     % Merges and splits
%     
%     split = rand(nnode,1) > 0.5;
%     for ii = 1:nnode
%         ng = max(grp);
%         
%         grpn = grp;
%         
%         if split(ii) % Split one group into two
%             
%             srcidx = ceil(rand(1)*ng);
%             isin = grp == srcidx;
%             nsub = sum(isin);
%             
%             if nsub > 1
%                 subidx = find(isin);
%                 if nsub == 2 % Only two nodes, so just split evenly
%                     grpn(subidx(1)) = ng + 1;
%                 else  
%                     unqc = unique(cidx(subidx));
%                     if length(unqc) == 2 % Only two connected components, split there
%                         move = cidx(subidx) == unqc(1);
%                         grpn(subidx(move)) = ng + 1;
%                     else % Otherwise, optimize the split
%                          
%                         switch Opt.split
%                             case 'random'
%                         
%                                 move = zeros(size(subidx)); 
%                                 while ~all(move) && ~any(move)
%                                     move = rand(size(subidx)) > 0.5; % Split randomly
%                                 end
% 
%                                 grpn(subidx(move)) = ng + 1;
%                        
%                             case 'optimal' % TODO update for tg/mod diffs
%                         
%                                 nsub = sum(isin);
%                                 subgrp = rand(nsub,1) > 0.5;
%                                 subadj = adj(isin,isin);
%                                 subdeg = sum(subadj,1)'+sum(subadj,2);
%                                 subL   = nnz(subadj);
%                                 tsub = tini;
% 
%                                 esplit = computeenergy(Opt.type, subgrp, subL, subdeg, subadj);
%                                 while tsub >= t
%                                     for is = 1:nsub^2
%                                         idx = ceil(rand(1)*nsub);
%                                         subgrpn = subgrp;
%                                         subgrpn(idx) = ~subgrpn(idx);
% 
%                                         esplitn = computeenergy(Opt.type, subgrpn, subL, subdeg, subadj);
% 
%                                         if esplitn > esplit || exp((esplitn - esplit)/tsub) > rand(1)
%                                             subgrp = subgrpn;
%                                             esplit = esplitn;
%                                         end 
%                                     end
%                                     tsub = tsub * ts;
% 
%                                 end
% 
%                                 subidx = find(isin);
%                                 grpn(subidx(subgrp)) = ng + 1;
%                         end
%                     
%                         
%                     end 
%                 end
%             end
%            
%         else % Merge two groups into one
%             
%             srcidx = ceil(rand(1,2)*ng);
%             tf = ismember(grpn, srcidx);
%             grpn(tf) = srcidx(1);
%             grpn(grp > srcidx(2)) = grp(grp > srcidx(2)) - 1; % Fill in the gap
%         end
%         
%         
%         energy = efun(grpn);
%         
%         if energy > energyant || exp((energy - energyant)/t) > rand(1)
%             grp = grpn;
%             energyant = energy;
%         end
%         
%         movecount = movecount + 1;
%         etrack(movecount,1) = energy;
%         etrack(movecount,2) = energyant;
%         
%         if Opt.plot
%             if movecount <= npt
%                 set(h.en, {'ydata'}, num2cell(etrack(1:npt,:),1)');
%             else
%                 set(h.en, {'ydata'}, num2cell(etrack(movecount-(npt-1):movecount,:),1)');
%             end
%             drawnow;
%         end
%     end
%     
%     % Termination criteria
%     
%     eprev(eprev == 0) = eps;
%     de = abs(energyant - eprev)/abs(eprev);
%     if  de < Opt.etol
%         count = count + 1;
%     else
%         count = 0;
%     end
%     eprev = energyant;
%     
%     if Opt.verbose
%         fprintf(repmat('\b', 1, 44));
%         fprintf('t = %7.4f, dE = %10.6f, E = %10.6f', t, de, energyant); 
%     end
%            
%     t = t * Opt.ts;
% 
% end
% if Opt.verbose
%     fprintf('\nTermination criteria reached\n');
% end
% 
% % Plot 
% 
% if Opt.plot
%     if ~verLessThan('matlab', 'R2015b')
%         set(h.gp, 'NodeCData', grp);
%     end
% end

if nargout == 2
    isn = any(isnan(etrack),2);
    A.etrack = etrack(~isn,:);
    A.adj = adj;
end

% Energy calculation

function energy = computemod(x, L, d, adj)

gidx = uniquefast(x);
ng = length(gidx);
        
m = zeros(ng,1);
for ig = 1:ng
    madj = adj(x == gidx(ig), x == gidx(ig));
    ls = nnz(madj);
    ds = sum(d(x == gidx(ig)));
    m(ig) = ls./L - (ds./(2*L)).^2;
end
energy = sum(m);

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

function energy = computetg(x, tminuse, ir)

insamegrp = bsxfun(@eq, x, x');

if any(insamegrp(:))
    n = histc(x, 1:max(x));
    grpnum = x(ir(insamegrp));
    
    dt = accumarray(grpnum, tminuse(insamegrp), [max(x) 1]);
    energy = nansum(dt./n);
    
else
    energy = 0;
end


function x = uniquefast(x)

xsrt = sort(x(:));
dx = diff(xsrt);
isin = [true; dx~=0];
x = xsrt(isin);
