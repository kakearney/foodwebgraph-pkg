function [pidx, etrack] = simanpartition(costfun, nnode, cidx, SAopt)
%SIMANPARTITION Partition a graph using simulated annealing
%
% This function implements the version of simulated annealing used by
% Gauzens et al, 2014 to partition a graph in node clusters such that a
% specified metric is maximized.
%
% Input variables:
% 
%   costfun:    function handle of function of the form cost = f(x), where
%               x is a nnode x 1 vector of cluster indices corresponding to
%               each node in a given partition, and cost is a scalar metric
%               to be maximized by the simulated annealing algorithm
%
%   nnode:      scalar, number of nodes in the graph
%
%   cidx:       nnode x 1 array, indices of connected components to which
%               each node belongs
%
%   SAopt:      structure of simulated annealing options

% Copyright 2016 Kelly Kearney

% Initialize with partition where each node is in its own cluster

pidx = (1:nnode)';
energyant = costfun(pidx);


% Initialize annealing 

ntmax = ceil(log(SAopt.tfin/SAopt.tini)./log(SAopt.ts)); % Maximum # temperature steps
nmoves = ntmax * (nnode^2 + nnode);                % Maximum # of moves
movecount = 0;                                     % Number of moves performed
etrack = nan(nmoves,2);                            % Track energy over time

t = SAopt.tini;
count = 0;
eprev = -Inf;


% Begin simulated annealing

if SAopt.verbose
    fmt = '  t = %7.4f, dE = %10.6f, E = %10.6f';
    fmt2 = '  t = %7.4f, dE = %10.6g, E = %10.6f';
    progressstr = sprintf(fmt, t, NaN, NaN);
    len = length(progressstr);
    
    fprintf('%s\n%s', 'Beginning simulated annealing', progressstr);
end

while t > SAopt.tfin && count < SAopt.limit

    % Individual movements: Choose a node at random and move it
    % to a different cluster.  Repeat nnode^2 times.

    try
        npercluster = accumarray(pidx, pidx, [nnode 1], @length);
    catch ME
        sstr = sprintf('%d ', pidx);
        fprintf('accumarray issue:\n  subs = %s\n  sz = %d', sstr, nnode);
        rethrow(ME);
    end
    openidx = find(npercluster > 0);
    ng = length(openidx);

    for ii = 1:nnode^2

        % Make move, calculate new energy metric

        srcidx = ceil(rand(1)*nnode);       % which node to switch
        taridx = openidx(ceil(rand(1)*ng)); % cluster to move it to
        grpn = pidx;
        grpn(srcidx) = taridx;

        energy = costfun(grpn);

        % Choose whether to accept the move

        if energy > energyant || exp((energy - energyant)/t) > rand(1)
            pidx = grpn;
            energyant = energy;
        end

        % For reference: energy of this move, energy of last
        % accepted move

        movecount = movecount + 1;
        etrack(movecount,1) = energy;
        etrack(movecount,2) = energyant;
    end

    % Cleanup: get rid of any now-empty clusters in the
    % partition (such that cluster indices remain consecutive)

    [~,~,pidx] = unique(pidx);

    % Split/merge movements: Either choose two clusters at
    % random, and try combining them, or choose one cluster and
    % split it into two.  Repeat nnode times.

    splitflag = rand(nnode,1) > 0.5; % choose between splitting and merging
    for ii = 1:nnode

        ncluster = max(pidx); % # clusters in current partition
        grpn = pidx;          % New partition (equal to current to start)

        if splitflag(ii) % Split one cluster into two

            srcidx = ceil(rand(1)*ncluster); % Which cluster to split
            isin = pidx == srcidx;           % Mask for nodes in that cluster
            nsub = sum(isin);                % # nodes in cluster

            if nsub > 1

                % Figure out how to split the cluster:
                % If there are only two nodes in the cluster,
                % move one to the new cluster.  If there are
                % more than two nodes, but they fall into
                % exactly two connected-component groups, then
                % split along these lines.  Otherwise, split
                % either at random.

                subidx = find(isin); 

                if nsub == 2 % Only two nodes, so just split evenly
                    grpn(subidx(1)) = ncluster + 1;
                else  
                    unqc = unique(cidx(subidx)); % Connected components to which each node belong
                    if length(unqc) == 2 % Only two connected components, split there
                        move = cidx(subidx) == unqc(1);
                        grpn(subidx(move)) = ng + 1;
                    else 
                        switch SAopt.split
                            case 'random'
                        
                                move = zeros(size(subidx)); 
                                while ~all(move) && ~any(move)
                                    move = rand(size(subidx)) > 0.5; % Split randomly
                                end

                                grpn(subidx(move)) = ng + 1;
                            case 'optimal'
 
                                movetest = rand(nsub,1) > 0.5;
                                subgrp = grpn;
                                subgrp(subidx(movetest)) = ng + 1;
                                esplit = costfun(subgrp);
                                
                                tsub = SAopt.tini; 
                                
                                c1 = ng+1;
                                c2 = srcidx;
                                
                                while tsub >= t
                                    for is = 1:nsub^2
                                        idx = ceil(rand(1)*nsub);
                                        
                                        subgrpn = subgrp;
                                        if subgrpn(subidx(idx)) == c1
                                            subgrpn(subidx(idx)) = c2;
                                        else
                                            subgrpn(subidx(idx)) = c1;
                                        end
                                        
                                        esplitn = costfun(subgrp);

                                        if esplitn > esplit || exp((esplitn - esplit)/tsub) > rand(1)
                                            subgrp = subgrpn;
                                            esplit = esplitn;
                                        end 
                                    end
                                    tsub = tsub * SAopt.ts;

                                end

                                grpn = subgrp;
                                
                            case 'optimal2'
                                
                                possiblesplits = partitions(subidx,2);
                                nposs = length(possiblesplits);
                                grptmp = grpn * ones(1,nposs);
                                for ii = 1:nposs
                                    grptmp(possiblesplits{ii}{1}) = ng + 1;
                                end
                                tmp = cellfun(costfun, num2cell(grptmp,1));
                                [~,imax] = max(tmp);
                                
                                grpn = grptmp(:,imax);
                                
                        end
                    end 
                end
            end

        else % Merge two groups into one

            srcidx = ceil(rand(1,2)*ng);
            tf = ismember(grpn, srcidx);
            grpn(tf) = srcidx(1);
            grpn(pidx > srcidx(2)) = pidx(pidx > srcidx(2)) - 1; % Fill in the gap
        end
        
        % Calculate new energy metric after split/merge

        energy = costfun(grpn);

        % Choose whether to accept the move

        if energy > energyant || exp((energy - energyant)/t) > rand(1)
            pidx = grpn;
            [~,~,pidx] = unique(pidx); % Clean up, fill in gaps
            energyant = energy;
        end
        
        % For reference: energy of this move, energy of last
        % accepted move

        movecount = movecount + 1;
        etrack(movecount,1) = energy;
        etrack(movecount,2) = energyant;

    end

    % Termination criteria: If change between current energy value
    % and the previous one is small enough, increment the counter.

    eprev(eprev == 0) = eps;
    de = abs(energyant - eprev)/abs(eprev);
    if  de < SAopt.etol
        count = count + 1;
    else
        count = 0;
    end

    % Print and plot progress
    
    if SAopt.verbose
        progressstr = sprintf(fmt, t, de, energyant); 
        if length(progressstr) ~= len
            progressstr = sprintf(fmt2, t, de, energyant); 
        end

        fprintf(repmat('\b', 1, len));
        fprintf('%s', progressstr);
        
    end
    
    eprev = energyant;
    
    % Decrease temperature

    t = t * SAopt.ts;

end

if SAopt.verbose
    fprintf('\nTermination criteria reached\n');
end

isn = all(isnan(etrack),2);
etrack = etrack(~isn,:);
