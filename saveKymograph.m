clear;
addpath('convertors')
addpath('plots')
analysis='C11';
analysis_trajectory='D11';

% ffolder={'Z:\AstraZeneca-protein_ladder\2021-06-02-protein_ladder_HEPESplusNaCl\protein_ladder_3000nm\channel_30x85_2_flow_in_microchannels\',...
%     'Z:\AstraZeneca-protein_ladder\2021-06-02-protein_ladder_HEPESplusNaCl\protein_ladder_3000nm\channel_30x85_2\'}; 
%ffolder={'Z:\simulated_tests\velocity0_distance20_timesteps10000/'};
ffolder={'C:\Users\lab\Desktop\results\experimental_noise\'};
% integrated optical contrast setting
iOC_setting.DLS=20; %width of the diffraction limited spot  - initial guess for gauss fitting [px]
iOC_setting.DLSW=2*iOC_setting.DLS+1; %number of selected of points around minima [px]

% denoising setting    
Yav=iOC_setting.DLS/2 - 1; %span of the gaussian moving average of the signal in spatial coordinate (imgaussfilt3) [px]
tav = 1; %span of the moving average of the signal in time coordinate [frame]
Yav_drift = 200; %span of the moving average in spatial coordinate for backround estimation [px]
tav_drift = 200; %span of the moving average in time coordinate for backround estimation [px]
Yav2=iOC_setting.DLS+1; %aperture for calculating the intensity moments [px]
PODM = 1; %iterative condition for masking procedure
W_width=2.5*iOC_setting.DLS; %aperture for masking the found responces [px]

% particle tracking settings
PT_setting.FACTOR_PARTICLE_LIMIT = 4; %for particcletracking
PT_setting.FACTOR_NOISE_LIMIT8 = 3;%for particcletracking
PT_setting.FACTOR_NOISE_LIMIT15 = 4;%for particcletracking
PT_setting.DIFF_COEFFUm=50; %[um^2/s]

% find boud molecule setting
filterBoundMolecule_settings.N=30;
filterBoundMolecule_settings.A=0.9;

fun_plot=1;

for ifolder=1:length(ffolder)
    folder=ffolder{ifolder};

    %% list of files
    d=dir(folder);
    [a1,a2]=cellfun(@size,strfind({d.name},'_S.mat')); 
    
    ff=find(a1==1);
    ttest=[]; kk=0;
    for i=1:length(ff)
        test0=d(ff(i)).name(1:end-6);
          [b1,b2]=cellfun(@size,strfind({d.name},strcat(test0,'_C11.mat')));
          if sum(b1)==0
            kk=kk+1;
            ttest{kk}=test0;
          end
    end
    
    %ttest={'24-08-20_16-06-15_D50_iOC0.002'};
        
for itest=1:length(ttest)
    test=ttest{itest}
    data=[];
    
    load(strcat(folder,'/',test,'_S.mat'));
    %load(strcat(folder,'/',test,'_M.mat'));
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
    
    
    data.Im=imgaussfilt3(data.Im,[1e-5,Yav,1e-5]);
    
     data.Yav_drift=Yav_drift;
     data.tav_drift=tav_drift;
     data.Yav2 = Yav2;
     data.Yav = Yav;
     data.tav = tav;
     data.PODM = PODM;
     data.W_width = W_width;
     data.iOC_setting=iOC_setting; 
     data.PT_setting=PT_setting;
     data.filterBoundMolecule_settings=filterBoundMolecule_settings;
     
    tic
    %cut time
%     data.Im=data.Im(3000:4000,:);
%     data.time=data.time(3000:4000);
           
     %% image & intensity stabilization
     DIFF_COEFF=PT_setting.DIFF_COEFFUm *(data.time(2)-data.time(1))/(data.Yum(2)-data.Yum(1))^2; %limiting Diff coeff [pixels^2/frame]
     Ycut=1+Yav_drift/2:size(data.Im,2)-Yav_drift/2;
     
     data0=data;
     Icut=data.Im(:,Ycut);
     W=ones(size(data.Im));
     podm=Inf; kk=0; 
     
     %STD_profile
             STD_profile0=1./sqrt(Icut(1,:))./mean(1./sqrt(Icut(1,:)));
             STD_profile0_fit=polyfit(Ycut,STD_profile0,3); 
             a=polyval(STD_profile0_fit,1:length(data.Yum));
             STD_profile0_fit=STD_profile0_fit/mean(a);
     
     while podm>=PODM && kk<20
         
         kk=kk+1;
         data=data0;
         I0cut=Icut;
         data.Im=data.Im./W;
%          if kk==1
            data=intensity_stabilization(data,'median');
%          else 
%              data=intensity_stabilization_fast3(data,'mean');
%          end
         data.Im=data.Im.*W;
         if data.tav>1
            data.Im=conv2(data.Im-1,ones(data.tav,1)/data.tav,'same');
            data.Im=data.Im+1;
         end
         data.Im=imgaussfilt3(data.Im,[1e-4,data.Yav,1e-4]);
         Icut=data.Im(:,Ycut);
         
         %setting the noise levels - STD_profile, M0_std_polyfit, M0_mean_polyfit
         if kk==1 
             
             % M0_std_polyfit, M0_mean_polyfit
             data_inv=data;
             data_inv.Im=2-data_inv.Im;
             OC_inv=intensityMoments(data_inv,data_inv.Yav2);
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
                
             
             a=Icut(Icut>=1);
             STD=sqrt(sum((a-1).^2)/(length(a)-1));
             STD_profile_fit=STD.*STD_profile0_fit;
             STD_profile=repmat(polyval(STD_profile_fit,1:size(W,2)),size(W,1),1);
%              
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
         
         %% intensity moments
         OC=intensityMoments(data,data.Yav2);
         OC.Yav2=data.Yav2;
         M_std = polyval(M0_std_polyfit,OC.position);
         M_mean = polyval(M0_mean_polyfit,OC.position);
         number_of_frames = 8;
         OC.VAR_m8_threshold = (M_std.^2*(1+0.72*PT_setting.FACTOR_PARTICLE_LIMIT/sqrt(number_of_frames))).^2;
         OC.MEAN_m8_threshold = M_mean - PT_setting.FACTOR_NOISE_LIMIT8*M_std/sqrt(number_of_frames);
         number_of_frames = 15;
         OC.VAR_m15_threshold = (M_std.^2*(1+0.72*PT_setting.FACTOR_PARTICLE_LIMIT/sqrt(number_of_frames))).^2;
         OC.MEAN_m15_threshold = M_mean - PT_setting.FACTOR_NOISE_LIMIT15*M_std/sqrt(number_of_frames);
         
         %% mask
          
          % mask from trajectories 
          display('finding trajectories')
          [MI_decided, MI_trajectory] = findTrajectory (OC, data,DIFF_COEFF,PT_setting.FACTOR_PARTICLE_LIMIT);
          trajectory_positionUm = OC.positionUm(MI_decided);
          trajectory_timeFrame = OC.timeFrame(MI_decided);
          W0=zeros(size(data.Im));
          for i=1:length(trajectory_positionUm)
              W0(trajectory_timeFrame(i),data.Yum>=trajectory_positionUm(i)-W_width/2*(data.Yum(2)-data.Yum(1)) & data.Yum<=trajectory_positionUm(i)+W_width/2*(data.Yum(2)-data.Yum(1)))=1;
          end

          % excluding mask
          W1=ones(size(data.Im));
          W1(1:data.tav_drift/2,:)=0;
          W1(end-data.tav_drift/2:end,:)=0;
          W1(:,1:data.Yav_drift/2)=0;
          W1(:,end-data.Yav_drift/2:end)=0;
          dtime=diff(data.time);
          a=median(dtime);
          timeFrame_jump=find(dtime>1.5*a);
          for i=1:length(timeFrame_jump)
                W1(max([1,timeFrame_jump(i)-data.tav_drift/2]):min([timeFrame_jump(i)+data.tav_drift/2],size(W1,1)),:)=0;
          end
          
          % mask from relative to STD combined with mask from trajectory
          W=data.Im-1;
          ff=abs(W)<3*STD_profile;
          W(ff)=W(ff).*(W(ff)./(3*STD_profile(ff))).^2;
          W(W1==0)=0;
          ff=W0==1;
          W(ff)=data.Im(ff)-1;
          W=W+1;
            
          podm=max(max(abs(Icut-I0cut)./STD_profile(:,Ycut)))
     end
      
     data.STD_profile_fit=STD_profile_fit;
     data.M0_std_polyfit=M0_std_polyfit;
     data.M0_mean_polyfit=M0_mean_polyfit;

    %% trajectory characteristics
    itra=0;
    trajectory=[];
    for i=1:length(MI_trajectory)
        if length(MI_trajectory{i})>=3
              itra=itra+1;
              trajectory(itra).timeFrame=OC.timeFrame(MI_trajectory{i});
              trajectory(itra).position=OC.position(MI_trajectory{i});
              trajectory(itra).positionUm=(trajectory(itra).position-1)*(data.Yum(2)-data.Yum(1))+data.Yum(1);
              trajectory(itra).ID_bound = filterBoundMolecules(trajectory(itra).position,filterBoundMolecule_settings.N,filterBoundMolecule_settings.A);
              trajectory(itra).N=length(trajectory(itra).position)-length(trajectory(itra).ID_bound);
              ff=setdiff(1:length(trajectory(itra).timeFrame),trajectory(itra).ID_bound);
              [trajectory(itra).Deff_mean, diffusion_coefficient_correction,trajectory(itra).Deff_std] = trajectoriesToDiffusivity(trajectory(itra).positionUm(ff)', data.time(trajectory(itra).timeFrame(ff)'), trajectory(itra).timeFrame(ff)',data.tav,'');
              trajectory(itra).iOC_mean = integrateOpticalContrast(data,trajectory(itra).timeFrame(ff),trajectory(itra).position(ff), iOC_setting.DLS, iOC_setting.DLSW);
        end
    end


if fun_plot==1
    
    %% plot kymograph
    %close all;
    figure('Position',[1 200 2000 500]);
    surf(1:size(data.Im,1),data.Yum,data.Im'); shading flat; view(2); colorbar; hold on
    colormap(bone)
    ylim([data.Yum(1) data.Yum(end)])
%     xlim([data.time(1) data.time(end)])
    xlim([1 size(data.Im,1)])
    ylabel('Position (\mum)')
%     xlabel('Time (s)')
    xlabel('Time frame')
    caxis([1-3e-4 1+3e-4])
    drawnow;
    shg;
    
    %% plot trajectory
    for itra=1:length(trajectory)
         plot3(trajectory(itra).timeFrame,trajectory(itra).positionUm,2*ones(size(trajectory(itra).timeFrame)),'Marker','.')
    end

    %% plot trajectories characteristics
    if length(trajectory)>0
        collection.Deff_mean=[trajectory.Deff_mean];
        collection.iOC_mean=[trajectory.iOC_mean];
        collection.N=[trajectory.N];

        %filter short trajectories
        a=collection.N>20;
        fname=fieldnames(collection);
        for i=1:length(fname)
            collection.(fname{i})= collection.(fname{i})(a);
        end
        
        %iOC x Deff histogram
        [N, edgesX, edgesY] = weighted_histcounts2 (collection.Deff_mean,collection.iOC_mean,collection.N,100);
        figure;
        surf(edgesY*1000,edgesX,N); view(2); shading flat
        xlabel('iOC (nm)')
        ylabel('Diffusivity (\mum^2/s)')
        
        % transform from iOC to MW and transform from Deff to HR
        nanochannelArea = list_of_nanochannels ('channel_30x85');
        HR = EffectiveDiffusivityToSize (collection.Deff_mean, nanochannelArea);
        MW = iOCToWeight (-collection.iOC_mean, nanochannelArea);
        
        % MW x HR histogram
        [N, edgesX, edgesY] = weighted_histcounts2 (HR,MW,collection.N,300);
        figure;
        surf(edgesY/1000,edgesX*1000,N); view(2); shading flat
        xlabel('Molecular weight (kDa)')
        ylabel('Hydrodynamic radius (nm)')
    end
    
    
end
    
    %%
    toc
    tic
    save(strcat(folder,'/',test,'_',analysis,'.mat'),'data');
    save(strcat(folder,'/',test,'_',analysis_trajectory,'.mat'),'trajectory');
    clear data trajectory
    toc
end
end

                
            