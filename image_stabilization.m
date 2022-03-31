function [dIm, x3, kk]=image_stabilization (Im,Im0,Yav)

% image stabilization  - correction for vibration - spatial shift of an image  
% more details in "B. Špačková et al.: Label-Free Nanofluidic Scattering Microscopy of Size and Mass of Single Diffusing Molecules and Nanoparticles", SI section 5  

%Im - one-dimensional matrix of intensities (intensities along the space) at (i+1)th frame
%Im0 - one-dimensional matrix of intensities (intensities along the space) at ith frame
%Yav - moving average parameter

%Im_temp - Im stabilized in terms of vibrations
%x3 - found space shift in (i=1)th frame [pixels]
%kk - number of interations
%dIm - ratio between Im and Im0 which is stabilized in Y position
%and normalized by moving mean
 
dx = 1e-2; %initial guess of space shift in the (i+1)th frame
threshold = 1e-4; %threshold value for the Newton-Rapson interation

Y = 1:length(Im0);

%% two initial guesses of spatial shift for Newton-Rapson interation
x1 = 0;
Im_temp=interp1(Y,Im,Y+x1,'spline');
dIm1 = Im_temp./Im0;
dIm1=dIm1./smooth(dIm1,Yav)';
        
x2 = dx;
Im_temp=interp1(Y,Im,Y+x2,'spline');
dIm2 = Im_temp./Im0;
dIm2=dIm2./smooth(dIm2,Yav)';

%% Newton-Rapson iteration
kk=1;
while abs(x2-x1)>threshold

            x3 = -(x2-x1)./(dIm2(2:end-1)-dIm1(2:end-1)+eps).*(dIm1(2:end-1)-1)+x1; 
            W = abs(dIm2(2:end-1)-dIm1(2:end-1)).*(Im(2:end-1)/max(Im));
            
            ff=isnan(x3)==0 & abs(x3)<Inf & isnan(W)==0 & abs(W)<Inf;
            x3 = sum(x3(ff).*W(ff))/sum(W(ff));
            Im_temp=interp1(Y,Im,Y+x3,'spline');

            dIm = Im_temp./Im0;
            
            I_mean=smooth(dIm,Yav)';
            dIm=dIm./I_mean;
            


            x1=x2; x2=x3;
            dIm1=dIm2; dIm2=dIm; 

            kk=kk+1;
            if kk>100
                disp('image_stabilizationY8: does not converge, kk>100'); return
            end
end
   
        

