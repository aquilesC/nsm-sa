function [opticalContrast]=intensityMoments(Yum,Im,Yav)

% finds positions of minimas and calculates their intensity moments
% based on Crocker, J.C., Grier, D.G., 1996. Methods of digital video microscopy for colloidal studies. J. Coll. Interface Sci. 179, 298?310.
% or I.F. Sbalzarini, P. Koumoutsakos / Journal of Structural Biology 151 (2005) 182?195

% Yum - space axis [um], recommended: Yum=data.Yum; 
% Im - matrix of intensities, recommended: Im=data.Im;
% Yav - width of the expected dips [pixel]

% opticalContrast.timeFrame - serie of number of the time frames correspondng to found minimas [frame #]
% opticalContrast.position - corresponding minima positions [pixel]
% opticalContrast.positionUm - corresponding minima positions [um]
% opticalContrast.m0 - correspondng 0-th intensity moment [pixel]

 
X=1:size(Im,2);
t=1:size(Im,1);
X2=repmat(X,length(t),1);
t2=repmat(t',1,length(X));
Y=Im-1;

m0=conv2(Y,ones(1,Yav),'same');
m1=conv2(Y.*X,ones(1,Yav),'same')./m0;
%m2=conv2(Y,[-(Yav-1)/2:(Yav-1)/2].^2,'same')./m0;
m1_difference=abs(X-m1);
 
    
A =  -Y;
Mask=ones(1,Yav); Mask((Yav+1)/2)=0;
B = imdilate(A,Mask);
dips = A > B;
dips_position=X2(dips);
opticalContrast.timeFrame=t2(dips);
 
%refinement
for i=1:length(dips_position)
        
     ff=max([dips_position(i)-(Yav-1)/2,1]):min([dips_position(i)+(Yav-1)/2,size(Y,2)]);
     [a,b]=min(m1_difference(opticalContrast.timeFrame(i),ff));
     opticalContrast.position(i,:)=m1(opticalContrast.timeFrame(i),ff(b));
     %opticalContrast.m2(i,:)=m2(opticalContrast.timeFrame(i),ff(b));
     opticalContrast.m0(i,:)=m0(opticalContrast.timeFrame(i),ff(b));

end
    
opticalContrast.positionUm=(opticalContrast.position-1)*(Yum(2)-Yum(1))+Yum(1); 
              
    

 

