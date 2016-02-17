function [tsim, etsim] = trophicsimilarity(adj)
%TROPHICSIMILARILTY Calculate trophic similarity
%
% tsim = trophicsimilarity(adj)

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