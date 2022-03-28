function iOC = integrateOpticalContrast(data,timeFrame, position, DLS, DLSW)

X11=-(DLSW-1)/2:(DLSW-1)/2;
N=length(timeFrame);
    
if N>=3
   
    
    X0=repmat(1:size(data.Im,2),N,1);
    Y0=repmat(timeFrame,1,size(X0,2));
    Z0=data.Im(timeFrame,:);
    X1=repmat(X11,N,1);
    X1=X1+repmat(position,1,size(X1,2));
    Y1=repmat(timeFrame,1,size(X1,2));
    
    Z1=interp2(X0,Y0,Z0,X1,Y1,'spline');
    Z1m=median(Z1,1)-1;
    
    fitresult=fit(X11',Z1m','gauss1',...
            'Lower',[-Inf,0,DLS],...
            'Upper',[0,0,3*DLS],...
            'TolFun',1e-9);
    iOC=sqrt(pi)*fitresult.a1*fitresult.c1*(data.Yum(2)-data.Yum(1));
    
else
        
    iOC=NaN;
    
end




