function [tsim, etsim] = trophicsimilarity(adj)
%TROPHICSIMILARILTY Calculate trophic similarity
%
% tsim = trophicsimilarity(adj)
%
% Calculates trophic similarity between groups in a food web, based on the
% equation from Gauzens et al. 2014 (J. R. Soc. Interface, volume 12).
%
% Input variables:
%
%   adj:    ngroup x ngroup adjacency matrix, 1 if column j eats row i
%
% Output variables:
%
%   tsim:   ngroup x ngroup trophic similarity matrix

% Copyright 2016 Kelly Kearney

nnode = size(adj,1);

npred = sum(adj,2);
nprey = sum(adj,1);
etsim = zeros(nnode);
    
tsim = zeros(nnode);
for ii = 1:nnode
    for jj = 1:nnode

        % Observed

        tsim(ii,jj) = (sum(adj(:,ii) & adj(:,jj)) + sum(adj(ii,:) & adj(jj,:))) ./ ...
                      (sum(adj(:,ii) | adj(:,jj)) + sum(adj(ii,:) | adj(jj,:)));

        % Expected in random graph

        N = 0:min(npred([ii jj]));
        n = 0:min(nprey([ii jj]));
        [NN,nn] = meshgrid(N,n);

        hg1 = hygepdf(N, nnode, min(npred([ii,jj])), max(npred([ii,jj])));
        hg2 = hygepdf(n, nnode, min(nprey([ii,jj])), max(nprey([ii,jj])));
        [hg1, hg2] = meshgrid(hg1, hg2);

        tmp = ((NN + nn)./(npred(ii) + npred(jj) - NN + nprey(ii) + nprey(jj) - nn)) .* ...
                hg1 .* hg2; 
        etsim(ii,jj) = nansum(tmp(:));

    end
end
tsim(isnan(tsim)) = 0; % Occurs only for disconnected groups  