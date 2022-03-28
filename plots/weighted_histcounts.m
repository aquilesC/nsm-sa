function [N, edges] = weighted_histcounts (X,n,nbins)

if sum(n==1)==length(n)
    X2=X;
else
    X2=[];
    for i = 1:length(X)
        X2 = [X2, repmat(X(i),1,n(i))];
    end
end

if length(nbins)==1
    
    [N, edges0] = histcounts (X2,nbins);
    
else
    
    edges0=(nbins(1:end-1) + nbins(2:end))/2;
    edges0=[2*edges0(1) - edges0(2), edges0,2*edges0(end)-edges0(end-1)];
    [N, edges0] = histcounts (X2,edges0);
    
end

edges=(edges0(1:end-1) + edges0(2:end))/2;
