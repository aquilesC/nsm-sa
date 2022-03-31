function Iout = intensity_stabilization(I, denoise_setting, fun)

% image stabilization and background subtration
% more details in "B. Špačková et al.: Label-Free Nanofluidic Scattering Microscopy of Size and Mass of Single Diffusing Molecules and Nanoparticles", SI section 5  
%
% I - raw data = 1D image in time = matrix of intensities [time frame, pixel]
% denoise_setting.Yav_drift - span of the moving average in spatial coordinate for backround estimation [px]
% denoise_setting.tav_drift - span of the moving average in time coordinate for backround estimation [time frames]
% fun - type of background estimation; 
%   fun = 'median' is based on moving median; 
%   fun = 'mean' is based on moving mean

% Iout - stabilized image with background subtracted

Yav_drift=denoise_setting.Yav_drift;
tav_drift=denoise_setting.tav_drift;

%% image stabilization
disp('image stabilization')
for i=1:size(I,1)-1 
%parfor i=1:size(I,1)-1 
        [Yc(i,:),x3(i),kk(i)]=image_stabilization (I(i,:),I(i+1,:),Yav_drift); %Y position stabilized intesity ratio between i+1 and i frame
end

Ic=ones(1,size(I,2));
for i=1:size(Yc,1)
    Ic(i+1,:)=Ic(i,:)./Yc(i,:);
end

%% subtration of the background
disp('drifts stabilization')
if strcmp(fun,'median')==1
    I_mean = movmedian(Ic, tav_drift,1);
elseif strcmp(fun,'mean')
    I_mean = conv2(Ic-1,ones(tav_drift,1)/tav_drift,'same');
    I_mean=I_mean+1;
end

% I_mean = moving_median(Ic, Yav_drift,2);
% I=Ic./I_mean;

Iout=Ic./I_mean;



