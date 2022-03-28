function data=intensity_stabilization(data,fun)

%corrects for vibration along Y axis (long axis of a nanochannel), short-term intensity instabilities and slow drifts

X=1:length(data.Yum);
I=data.Im;
Yav_drift=data.Yav_drift;
tav_drift=data.tav_drift;


disp('image stabilization')
parfor i=1:size(I,1)-1          
%                      Yc=image_stabilizationY8 (X,I(i,:),I(round(size(I,1)/2),:),Yav_drift); %Y position stabilized intesity ratio between i+1 and i frame
%                      Ic(i,:)=I(round(size(I,1)/2),:).*Yc;
%     if tav>1
%         Yc=image_stabilization (X,mean(I(i:i+tav-1,:),1),mean(I(i+tav:i+2*tav-1,:),1),Yav_drift); %Y position stabilized intesity ratio between i+1 and i frame
%     else
        [Yc(i,:),x3(i),kk(i)]=image_stabilization (X,I(i,:),I(i+1,:),Yav_drift); %Y position stabilized intesity ratio between i+1 and i frame
   

end
Ic=ones(1,size(I,2));
for i=1:size(Yc,1)
    Ic(i+1,:)=Ic(i,:)./Yc(i,:);
end
I=Ic;


disp('drifts stabilization')
if strcmp(fun,'median')==1
    I_mean = movmedian(I, tav_drift,1);
    %I_mean = moving_median(I, tav_drift,1);
elseif strcmp(fun,'mean')
    I_mean = conv2(I-1,ones(tav_drift,1)/tav_drift,'same');
    I_mean=I_mean+1;
end
% Iplus=[I(1,:)-flipud(I(2:tav_drift+1,:)-I(1,:)); I; I(end,:)-(flipud(I(end-tav_drift:end-1,:))-I(end,:))];
% I_mean = conv2(Iplus,ones(tav_drift,1)/tav_drift,'same');
% I_mean = I_mean(tav_drift+1:end-tav_drift,:);
I=I./I_mean;
%I1=I;
% I_mean = moving_median(I, Yav_drift,2);
% I=I./I_mean;

data.Im=I;