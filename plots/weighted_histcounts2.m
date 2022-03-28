function [N, edgesX, edgesY] = weighted_histcounts2 (X,Y,n,nbinsX,nbinsY)

if sum(n==1)==length(n)
    X2=X;
    Y2=Y;
else
    
    X2=[]; Y2=[];
    for i = 1:length(X)
        X2 = [X2, repmat(X(i),1,n(i))];
        Y2 = [Y2, repmat(Y(i),1,n(i))];
    end
end

if nargin==4
    
    [N, edgesX0, edgesY0] = histcounts2 (X2,Y2,nbinsX);
    
else
    
    edgesX0=(nbinsX(1:end-1) + nbinsX(2:end))/2;
    edgesX0=[2*edgesX0(1) - edgesX0(2), edgesX0,2*edgesX0(end)-edgesX0(end-1)];
    edgesY0=(nbinsY(1:end-1) + nbinsY(2:end))/2;
    edgesY0=[2*edgesY0(1) - edgesY0(2), edgesY0,2*edgesY0(end)-edgesY0(end-1)];
    [N, edgesX0, edgesY0] = histcounts2 (X2,Y2,edgesX0,edgesY0);
    
end

edgesX=(edgesX0(1:end-1) + edgesX0(2:end))/2;
edgesY=(edgesY0(1:end-1) + edgesY0(2:end))/2;
