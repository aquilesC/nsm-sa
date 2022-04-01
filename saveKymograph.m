% main file for the standard analysis (SA) for Nanofluidic Scattering Microscopy (NSM), 
% described in "B. Špačková et al.: Label-Free Nanofluidic Scattering Microscopy of Size and Mass of Single Diffusing Molecules and Nanoparticles". 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INPUT file = data collected by the camera (saved as strcat(ExperimentTimeStamp,'_M.mat'))
% data.Im - raw data = 1D image in time = matrix of intensities [time frame, pixel] 
% data.time - temporal coordinate [s]
% data.Yum - spatial coordinate [um]


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% OUTPUT files 
% = data - struct containing the processed images (saved as strcat(ExperimentTimeStamp,'_C.mat'))
% data.Im - kymograph = 1D image in time = matrix of intensities [time frame, pixel] 
% data.time - temporal coordinate [s]
% data.Yum - spatial coordinate [um]

% = trajectory containing characteristics of a trajectory (saved as strcat(ExperimentTimeStamp,'_D.mat'))
% trajectory.timeFrame - temporal coordinates of found trajectories [number of frame]
% trajectroy.position - spatial coordinates of found trajectories [pixel]
% trajectroy.positionUm - spatial coordinates of found trajectories [um]
% trajectory.ID_bound - indexes of found minimas corresponding to bound particles (not diffusing)
% trajectory.N - temporal length of foudn trajectory
% trajectory.Deff_mean - diffusivity of a particle [um2/s]
% trajectory.iOC_mean - integrated optical contrast of a particle [um]


clear;
close all;

ExperimentTimeStamp={'./data/24-08-20_16-06-15_D10_iOC0.0005'}; %names of the INPUT files (without the letter at the end) to be anlyzed


% SETINGS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% integrated optical contrast setting
iOC_setting.DLS=20; %width of the diffraction limited spot  - initial guess for gauss fitting [px]
iOC_setting.DLSW=2*iOC_setting.DLS+1; %number of selected of points around minima [px]
iOC_setting.Yav2=iOC_setting.DLS+1; %aperture for calculating the intensity moments [px]
iOC_setting.DLSW2=2.5*iOC_setting.DLS; %aperture for masking the found responces [px]

% denoising setting    
denoise_setting.Yav=iOC_setting.DLS/2 - 1; %span of the gaussian moving average of the signal in spatial coordinate (imgaussfilt3) [px]
denoise_setting.tav = 1; %span of the moving average of the signal in time coordinate [frame]
denoise_setting.Yav_drift = 200; %span of the moving average in spatial coordinate for backround estimation [px]
denoise_setting.tav_drift = 200; %span of the moving average in time coordinate for backround estimation [time frames]
denoise_setting.PODM = 1; %iterative condition for masking procedure; put inf for quick and unprecise analysis, put 1 for slow precise analysis 

% particle tracking settings
PT_setting.FACTOR_PARTICLE_LIMIT = 4; %factor setting maximal values of STD(m) of a trajectory that are considred to be a particle
PT_setting.FACTOR_NOISE_LIMIT8 = 3; %factor setting maximal values of -MEAN(m) of a 8-frame-trajectory that are considered to be noise 
PT_setting.FACTOR_NOISE_LIMIT15 = 4; %factor setting maximal values of -MEAN(m) of a 15-frame-trajectory that are considered to be noise
PT_setting.DIFF_COEFFUm=50; %D_limit - the highest diffusivity of a particle that are found by the algorithm [um^2/s]
PT_setting.MinimalTemporalLength=40; %all found trajectroeis whose temporal length is lower than PT_setting.MinimalTemporalLength are considered to be noise and are discarded

% find bound molecule setting
filterBoundMolecule_settings.N=30;
filterBoundMolecule_settings.A=0.9;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

addpath('convertors')
addpath('plots')
        
for itest=1:length(ExperimentTimeStamp)
    
    clear data trajectory
    tic
    
    %% load raw data
    load(strcat(ExperimentTimeStamp{itest},'_M.mat'));
    if size(data.Im,2)==length(data.time)
        disp('changing dimensions of data.Im!')
        data.Im=data.Im';
    end
    
    %% exclude zero intensity frames
    ff=find(mean(data.Im,2)==0);
    if length(ff>0)
        display(strcat('darkness at frame no',num2str(ff)));
        data.time(ff)=[];
        data.Im(ff,:)=[];
    end
    
    %% add SA settings to data file 
    data.denoise_setting=denoise_setting;
    data.iOC_setting=iOC_setting; 
    data.PT_setting=PT_setting;
    data.filterBoundMolecule_settings=filterBoundMolecule_settings;
    
%     %% select time frames
%     data.Im=data.Im(1:1000,:);
%     data.time=data.time(1:1000);

    %% convolution of the optical responce with gaussian profile 
    data.Im=imgaussfilt3(data.Im,[1e-5,data.denoise_setting.Yav,1e-5]);
           
    %% initial values for iteration
    Ycut=1+data.denoise_setting.Yav_drift/2:size(data.Im,2)-data.denoise_setting.Yav_drift/2;
    data0=data;
    Icut=data.Im(:,Ycut);
    W=ones(size(data.Im));
    podm=Inf; kk=0; 
     
    %% estimation of STD of intenity profile
    STD_profile0=1./sqrt(Icut(1,:))./mean(1./sqrt(Icut(1,:)));
    STD_profile0_fit=polyfit(Ycut,STD_profile0,3); 
    a=polyval(STD_profile0_fit,1:length(data.Yum));
    STD_profile0_fit=STD_profile0_fit/mean(a);
     
     while podm>=data.denoise_setting.PODM && kk<20
         
         kk=kk+1;
         data=data0;
         I0cut=Icut;

         %% image stabilization and background subtraction
         data.Im=data.Im./W;
         data.Im=intensity_stabilization(data.Im,denoise_setting,'median');
         data.Im=data.Im.*W;
         if data.denoise_setting.tav>1
            data.Im=conv2(data.Im-1,ones(data.denoise_setting.tav,1)/data.denoise_setting.tav,'same');
            data.Im=data.Im+1;
         end
         data.Im=imgaussfilt3(data.Im,[1e-4,data.denoise_setting.Yav,1e-4]);
         Icut=data.Im(:,Ycut);
         
         %% estimation of noise limits 
         if kk==1 
             
             % estimation of noise limits for m0 (zero intensity moments): M0_std_polyfit, M0_mean_polyfit = std and mean values of m0 
             data_inv=data;
             data_inv.Im=2-data_inv.Im;
             OC_inv=intensityMoments(data_inv.Yum,data_inv.Im,data.iOC_setting.Yav2);
             bin=ceil(OC_inv.position);
             M0=OC_inv.m0;
             for i=1:length(Ycut)
                    M0_temp=M0(bin==Ycut(i));
                    M0_mean_temp=median(M0_temp);
                    M0_std_temp=sqrt(sum((M0_temp-M0_mean_temp).^2)/(length(M0_temp)-1));
                    ff=find(abs(M0_temp-M0_mean_temp)>3*M0_std_temp);
                    while length(ff)>0
                        M0_temp(ff)=[];
                        M0_mean_temp=median(M0_temp);
                        M0_std_temp=sqrt(sum((M0_temp-M0_mean_temp).^2)/(length(M0_temp)-1));
                        ff=find(abs(M0_temp-M0_mean_temp)>3*M0_std_temp);
                    end
                    M0_std(i)=M0_std_temp;
                    M0_mean(i)=M0_mean_temp;
             end
             M0_std_polyfit=polyfit(Ycut,M0_std,3);
             M0_mean_polyfit=polyfit(Ycut,M0_mean,3);
                
             % estimation of noise limits for intensity: STD_profile = std of intensity
             a=Icut(Icut>=1);
             STD=sqrt(sum((a-1).^2)/(length(a)-1));
             STD_profile_fit=STD.*STD_profile0_fit;
             STD_profile=repmat(polyval(STD_profile_fit,1:size(W,2)),size(W,1),1);

%                 % plot noise limits 
%                 figure
%                 subplot(2,2,1);
%                 plot(Ycut,M0_std); hold on
%                 plot(Ycut,polyval(M0_std_polyfit,Ycut))
%                 ylabel('M0 std')
%                 subplot(2,2,2);
%                 plot(Ycut,M0_mean); hold on
%                 plot(Ycut,polyval(M0_mean_polyfit,Ycut))
%                 ylabel('M0 mean')
%                 subplot(2,2,3)
%                 plot(Ycut,STD_profile(1,Ycut))
%                 ylabel('STD')
              
         end
         
         %% calculation of intensity moments
         OC=intensityMoments(data.Yum,data.Im,data.iOC_setting.Yav2);

         %% definition of thresholds for particle tracking
         M_std = polyval(M0_std_polyfit,OC.position);
         M_mean = polyval(M0_mean_polyfit,OC.position);
         number_of_frames = 8;
         OC.VAR_m8_threshold = (M_std*(1+0.72*PT_setting.FACTOR_PARTICLE_LIMIT/sqrt(number_of_frames))).^2; %(sigma_limit)^2 for 8-frame-trajectories
         OC.MEAN_m8_threshold = M_mean - PT_setting.FACTOR_NOISE_LIMIT8*M_std/sqrt(number_of_frames); %m_limit for 8-frame-trajectories
         number_of_frames = 15;
         OC.VAR_m15_threshold = (M_std*(1+0.72*PT_setting.FACTOR_PARTICLE_LIMIT/sqrt(number_of_frames))).^2; %(sigma_limit)^2 for 15-frame-trajectories
         OC.MEAN_m15_threshold = M_mean - PT_setting.FACTOR_NOISE_LIMIT15*M_std/sqrt(number_of_frames); %%m_limit for 15-frame-trajectories
         PT_setting.DIFF_COEFF=PT_setting.DIFF_COEFFUm *(data.time(2)-data.time(1))/(data.Yum(2)-data.Yum(1))^2; %D_limit [pixels^2/frame]
    
         
         %% particle tracking 
          display('finding trajectories')
          [MI_decided, MI_trajectory] = findTrajectory (OC, PT_setting.DIFF_COEFF, PT_setting.FACTOR_PARTICLE_LIMIT);

          %% mask 
          % mask from found particles
          trajectory_positionUm = OC.positionUm(MI_decided);
          trajectory_timeFrame = OC.timeFrame(MI_decided);
          W0=zeros(size(data.Im));
          for i=1:length(trajectory_positionUm)
              W0(trajectory_timeFrame(i),data.Yum>=trajectory_positionUm(i)-data.iOC_setting.DLSW2/2*(data.Yum(2)-data.Yum(1)) & data.Yum<=trajectory_positionUm(i)+data.iOC_setting.DLSW2/2*(data.Yum(2)-data.Yum(1)))=1;
          end

          % exclude parts with not equidistant timing (when collection of data stopped)
          W1=ones(size(data.Im));
          W1(1:data.denoise_setting.tav_drift/2,:)=0;
          W1(end-data.denoise_setting.tav_drift/2:end,:)=0;
          W1(:,1:data.denoise_setting.Yav_drift/2)=0;
          W1(:,end-data.denoise_setting.Yav_drift/2:end)=0;
          dtime=diff(data.time);
          a=median(dtime);
          timeFrame_jump=find(dtime>1.5*a);
          for i=1:length(timeFrame_jump)
                W1(max([1,timeFrame_jump(i)-data.denoise_setting.tav_drift/2]):min([timeFrame_jump(i)+data.denoise_setting.tav_drift/2],size(W1,1)),:)=0;
          end
          
          % combine mask with estimation based on intensities relative to STD of noise
          W=data.Im-1;
          ff=abs(W)<3*STD_profile;
          W(ff)=W(ff).*(W(ff)./(3*STD_profile(ff))).^2;
          W(W1==0)=0;
          ff=W0==1;
          W(ff)=data.Im(ff)-1;
          W=W+1;
            
          podm=max(max(abs(Icut-I0cut)./STD_profile(:,Ycut))); %the indicator of difference between the iteration steps
          disp([num2str(100*kk./20) ' % done... '])
     end
     disp('100% Done')
     data.STD_profile_fit=STD_profile_fit;
     data.M0_std_polyfit=M0_std_polyfit;
     data.M0_mean_polyfit=M0_mean_polyfit;

    %% trajectory characteristics
    itra=0;
    trajectory=[];
    for i=1:length(MI_trajectory)
        if length(MI_trajectory{i})>=data.PT_setting.MinimalTemporalLength
              itra=itra+1;
              trajectory(itra).timeFrame=OC.timeFrame(MI_trajectory{i});
              trajectory(itra).position=OC.position(MI_trajectory{i});
              trajectory(itra).positionUm=(trajectory(itra).position-1)*(data.Yum(2)-data.Yum(1))+data.Yum(1);
              trajectory(itra).ID_bound = filterBoundMolecules(trajectory(itra).position,filterBoundMolecule_settings.N,filterBoundMolecule_settings.A);
              trajectory(itra).N=length(trajectory(itra).position)-length(trajectory(itra).ID_bound);
              ff=setdiff(1:length(trajectory(itra).timeFrame),trajectory(itra).ID_bound);
              [trajectory(itra).Deff_mean, diffusion_coefficient_correction,trajectory(itra).Deff_std] = trajectoriesToDiffusivity(trajectory(itra).positionUm(ff)', data.time(trajectory(itra).timeFrame(ff)'), trajectory(itra).timeFrame(ff)',data.denoise_setting.tav);
              trajectory(itra).iOC_mean = integrateOpticalContrast(data,trajectory(itra).timeFrame(ff),trajectory(itra).position(ff), iOC_setting.DLS, iOC_setting.DLSW);
        end
    end

    %% save results
    save(strcat(ExperimentTimeStamp{itest},'_C.mat'),'data');
    save(strcat(ExperimentTimeStamp{itest},'_D.mat'),'trajectory');
    
    toc
end


                
            