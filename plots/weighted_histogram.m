function [iOC_edges, Deff_edges, iOC_edges2, Deff_edges2, iOC, Deff, iOCxDeff] = weighted_histogram (collection,iOC_edges_step,Deff_edges_step, fun, iOC_limit, Deff_limit)


%fun = 1:   each trajectory is included as 2 dimensional gaussian iOCxDeff
%           with widths iOC_error and Deff_error. The height is set to fulfill
%           sum(sum(iOCxDeff)=N where N is number of frames of the trajectory

%fun =2:    each trajectory is included as N points with defined iOC (integrated optical contrast) and Deff (diffusivity). N is number of the frames of the trajectory 

%collection.N - number of frames of each trajectory
%collection.Deff - diffusivity of the trajectory
%collection.iOC - integrated optical contrast of the trajectory

%1D histogram of iOC can be plotted as bar(iOC_edges, iOC)
%1D histogram of diffusivity can be plotter as bar(Deff_edges, Deff)
%2D histogram of iOC x diffusivity can be plotter as surf(iOC_edges2, Deff_edges2, iOCxDeff)

if fun==1

    if nargin==6
        iOC_edges=iOC_limit(1):iOC_edges_step:iOC_limit(2);
        Deff_edges=Deff_limit(1):Deff_edges_step:Deff_limit(2);
    
    else
    
        iOC_edges=min(collection.iOC_mean)-iOC_edges_step:iOC_edges_step:max(collection.iOC_mean)+iOC_edges_step;
        Deff_edges=min(collection.Deff_mean)-Deff_edges_step:Deff_edges_step:max(collection.Deff_mean)+Deff_edges_step;
    end

    [iOC_edges2,Deff_edges2]=meshgrid(iOC_edges,Deff_edges);

    iOCxDeff=zeros(size(iOC_edges2));
    iOC=zeros(size(iOC_edges));
    Deff=zeros(size(Deff_edges));
    for i=1:length(collection.iOC_mean)
        iOCxDeff = iOCxDeff + exp(-0.5*(((Deff_edges2-collection.Deff_mean(i))/collection.Deff_error(i)).^2 + ((iOC_edges2-collection.iOC_mean(i))/collection.iOC_error(i)).^2))/collection.Deff_error(i)/collection.iOC_error(i)/pi/2*...
            collection.N(i)*iOC_edges_step*Deff_edges_step;
        iOC = iOC + exp(-0.5*(((iOC_edges-collection.iOC_mean(i))/collection.iOC_error(i)).^2))/collection.iOC_error(i)/sqrt(pi*2)*collection.N(i)*iOC_edges_step;
        Deff = Deff + exp(-0.5*(((Deff_edges-collection.Deff_mean(i))/collection.Deff_error(i)).^2))/collection.Deff_error(i)/sqrt(pi*2)*collection.N(i)*Deff_edges_step;
    end
    
elseif fun==2
    
    if nargin==6
        iOC_edges=iOC_limit(1):iOC_edges_step:iOC_limit(2);
        Deff_edges=Deff_limit(1):Deff_edges_step:Deff_limit(2);
    
    else
    
        iOC_edges=min(collection.iOC_mean)-2*iOC_edges_step:iOC_edges_step:max(collection.iOC_mean)+2*iOC_edges_step;
        Deff_edges=min(collection.Deff_mean)-2*Deff_edges_step:Deff_edges_step:max(collection.Deff_mean)+2*Deff_edges_step;
    end
    
    iOC_all=[]; Deff_all=[];
    for i=1:length(collection.iOC_mean)
        iOC_all=[iOC_all,repmat(collection.iOC_mean(i),1,collection.N(i))];
        Deff_all=[Deff_all,repmat(collection.Deff_mean(i),1,collection.N(i))];
    end
    
    [Deff,Deff_edges]=histcounts(Deff_all,Deff_edges);
    [iOC,iOC_edges]=histcounts(iOC_all,iOC_edges);
    [iOCxDeff,iOC_edges,Deff_edges]=histcounts2(iOC_all,Deff_all,iOC_edges,Deff_edges);
    
    iOC_edges=(iOC_edges(1:end-1)+iOC_edges(2:end))/2;
    Deff_edges=(Deff_edges(1:end-1)+Deff_edges(2:end))/2; 
    [iOC_edges2,Deff_edges2]=meshgrid(iOC_edges,Deff_edges);
    
end
    
