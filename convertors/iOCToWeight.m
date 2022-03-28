function [weight,calibration] = iOCToWeight (iOC, A)

%clear;
%addpath('/Users/spackova/Box Sync/_project/SPR')

%protein, BSA
% lambda=0.6;
% r_bsa=3.7e-3; %marek rika 3.7e-3
% n_bsa=1.44; %marek rika 1.44
% n_m=1.33;
% weight_bsa=66e3; %Da %marek rika 65kDa

% [alpha_bsa,L] = alpha_comp2('es_withoutLW', [r_bsa,r_bsa,r_bsa], 0, n_m^2, n_m^2, n_bsa^2, lambda);
% alpha_bsa=alpha_bsa(1)*4*pi;
%alpha_bsa = 2.8986e-08;

%nanochannel
n_m=1.33;
n_c=1.46; %??
% r_c=50e-3;
% A=pi*r_c^2;
%A=0.1*0.17;

%weight=logspace(4,8); %Da

% contrast_bsa_TM=alpha_bsa/A*(2*n_m^2)/(n_m^2 - n_c^2);
% contrast_bsa_TE=alpha_bsa/A*(n_m^2 + n_c^2)/(n_m^2 - n_c^2);
% contrast_bsa=alpha_bsa/A*(3*n_m^2 + n_c^2)/(n_m^2 - n_c^2)/2;
% 
% weight=iOC./contrast_bsa*weight_bsa;

% figure;
% loglog(weight,- contrast_int);
% xlabel('Molecular weight (Da)');
% ylabel('\int Contrast')

n_mean = (3*n_m^2 + n_c^2)/(n_m^2 - n_c^2)/2;
alpha_MW = 0.461e-12; %[mum2/Da]
calibration = A/alpha_MW/n_mean

weight = iOC*calibration;