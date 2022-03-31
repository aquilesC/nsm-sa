function [weight,calibration] = iOCToWeight (iOC, channel_parameters, material)

%iOC  - integrated optical contrast [um]
%weight - molecular weight [Da]
%calibration = weight/iOC [Da/um]

%FOR BARE NANOCHANNEL
%channel_parameters = A;

% FOR COATED NANOCHANNEL
%channel_parameters -  [A, A_i, I_rel]

%A - area of a nanochannel [um2]
%A_i - are of the nanochannel without the coating [um2]
%I_rel - difference between the intensity before and after the deposition - I_without/I_with

%material - default: protein; other option: lipids

n_i=1.33; %RI of the inside of a nanochannel
n_o=1.46; %RI of the outside of a nanochannel


n_TE = 2*n_i^2/(n_i^2 - n_o^2);
n_TM = (n_i^2 + n_o^2)/(n_i^2 - n_o^2);
n_mean = 0.5*(n_TE + n_TM);

if length(channel_parameters)==1
    
    A=channel_parameters;
    
else 
    
    A=channel_parameters(1);
    A_i=channel_parameters(2);
    I_rel=channel_parameters(3);
    
    f=A_i/A;
        
      n_s = sqrt((n_i^2-n_o^2)*(sqrt(1/I_rel)-f)/(1-f)+n_o^2); %for TE
      %n_s = %for TM
      
      n_mean=n_mean*sqrt(I_rel);
    
end

if nargin == 2 %protein
    alpha_MW = 0.461e-12; %[mum2/Da]
else
    if strcmp(material,'lipid')==1
        alpha_MW = 0.461e-12/1.85*1.6; %[mum2/Da]
    end
end
        
calibration = A/alpha_MW/n_mean;
weight = iOC*calibration;
