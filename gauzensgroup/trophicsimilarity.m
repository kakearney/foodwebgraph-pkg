function [tsim, etsim] = trophicsimilarity(adj)
%TROPHICSIMILARILTY Calculate trophic similarity
%
% tsim = trophicsimilarity(adj)
%
% Calculates trophic similarity between groups in a food web, based on the
% equation from Gauzens et al. 2014 (J. R. Soc. Interface, volume 12).
% Trophic similarity is defined as the ratio between the number of predator
% and prey species interacting with both of two node groups in the web
% versus the number of predator and prey groups interacting with either of
% those two groups (i.e intersection/union).
%
% Input variables:
%
%   adj:    ngroup x ngroup adjacency matrix.  adj(i,j) = 1 if group j
%           preys upon group i, and 0 otherwise.
%
% Output variables:
%
%   tsim:   ngroup x ngroup array, trophic similarity matrix.
%
%   etsim:  ngroup x ngroup array, expected trophic similarity, i.e.
%           trophic similarity in a random graph.

% Copyright 2016 Kelly Kearney

nnode = size(adj,1);

npred = sum(adj,2);
nprey = sum(adj,1);

etsim = zeros(nnode);
tsim  = zeros(nnode);
    
% The hypergeometric pdf calculation is pretty time-consuming.  Here, we
% calculate it only the minimum number of times necessary, for each combo
% of parameters

[minpred, maxpred, minprey, maxprey] = deal(zeros(nnode));
for ii = 1:nnode
    for jj = 1:nnode
        minpred(ii,jj) = min(npred([ii,jj]));
        maxpred(ii,jj) = max(npred([ii,jj]));
        minprey(ii,jj) = min(nprey([ii,jj]));
        maxprey(ii,jj) = max(nprey([ii,jj]));
    end
end

[combos, ~, cidx] = unique([minpred(:) maxpred(:) minprey(:) maxprey(:)], 'rows');
nc = size(combos,1);
[NN,nn,hg1,hg2] = deal(cell(nc,1));
for ic = 1:nc
    N = 0:combos(ic,1); % 0:minpred
    n = 0:combos(ic,3); % 0:minprey

    [NN{ic},nn{ic}] = meshgrid(N,n);

    hg1tmp = hygepdf(N, nnode, combos(ic,1), combos(ic,2));
    hg2tmp = hygepdf(n, nnode, combos(ic,3), combos(ic,4));
    [hg1{ic}, hg2{ic}] = meshgrid(hg1tmp, hg2tmp);
    
end
cidx = reshape(cidx, nnode, nnode);

% Calculate trophic similarity for each pairwise node combo

for ii = 1:nnode
    for jj = 1:nnode

        % Observed

        tsim(ii,jj) = (sum(adj(:,ii) & adj(:,jj)) + sum(adj(ii,:) & adj(jj,:))) ./ ...
                      (sum(adj(:,ii) | adj(:,jj)) + sum(adj(ii,:) | adj(jj,:)));

        % Expected in random graph
        
%         N = 0:min(npred([ii jj]));
%         n = 0:min(nprey([ii jj]));
%         
%         [NNo,nno] = meshgrid(N,n);
% 
%         hg1o = hygepdf(N, nnode, min(npred([ii,jj])), max(npred([ii,jj])));
%         hg2o = hygepdf(n, nnode, min(nprey([ii,jj])), max(nprey([ii,jj])));
%         [hg1o, hg2o] = meshgrid(hg1o, hg2o);
% 
%         tmp = ((NNo + nno)./(npred(ii) + npred(jj) - NNo + nprey(ii) + nprey(jj) - nno)) .* ...
%                 hg1o .* hg2o; 
%         etsim(ii,jj) = nansum(tmp(:));
        
        % Expected in random graph
        
        c = cidx(ii,jj);
        
        tmp = ((NN{c} + nn{c})./(npred(ii) + npred(jj) - NN{c} + nprey(ii) + nprey(jj) - nn{c})) .* ...
                hg1{c} .* hg2{c}; 
        etsim(ii,jj) = nansum(tmp(:));

    end
end
tsim(isnan(tsim)) = 0; % Occurs only for disconnected groups  