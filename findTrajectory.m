function [MI_decided, MI_trajectory] = findTrajectory (OC, DIFF_COEFF, FACTOR_PARTICLE_LIMIT)

%particle tracking algorithm, described in "B. Špačková et al.: Label-Free Nanofluidic Scattering Microscopy of Size and Mass of Single Diffusing Molecules and Nanoparticles". 

% OC.timeFrame - number of time frames correspondng to found minimas [frame #]
% OC.position - corresponding minima positions [pixel]
% OC.m0 - correspondng 0-th intensity moment [pixel]

% FACTOR_PARTICLE_LIMIT

% MI_decided - indexes of minima that were attributed to particles
% MI_trajectory - groups of indexes of minima that were attributed to a trajectory {index of trajectory}(index of minima)


%% D_limits
DISTANCE_LIMIT=sqrt(2*DIFF_COEFF) * 3;
number_of_frames = 8;
VAR_Dx8_threshold =  DIFF_COEFF*(1+0.72*FACTOR_PARTICLE_LIMIT/sqrt(number_of_frames-1)).^2;
number_of_frames=15;
VAR_Dx15_threshold =  DIFF_COEFF*(1+0.72*FACTOR_PARTICLE_LIMIT/sqrt(number_of_frames-1)).^2;

                
%% indexes of minima
MI_0=1:length(OC.timeFrame);
MI=[];
MI{1}=MI_0(OC.timeFrame==1);
for it=2:max(OC.timeFrame)
     MI{it}=MI_0(OC.timeFrame==it);
end
                     

T=841; %for parallel computing: number of frames analyzed by a single core
%parfor iT=1:ceil((length(MI)-8)/T)
for iT=1:ceil((length(MI)-8)/T) %iT is the number of core
               
               MI_variation=[];
               SUM_Dx2 = [];
               SUM_m = [];
               SUM_m2 = [];
               SUM_m_withoutend = [];
               SUM_m2_withoutend = [];
               ID=[];
               cost=[];
               
                %% 8-frame-trajectories
                number_of_frames = 8;
                for it=1:7:min([T,length(MI)-((iT-1)*T+7)]) %it=(iT-1)*T+7:min([iT*T+6,length(MI)-1])
                    it0=(it-1)/7+1;
                
                    % MI allowed for 2-frame-trajectories 
                    MI_variation12=combvec(MI{(iT-1)*T+it},MI{(iT-1)*T+it+1});
                    MI_variation12=MI_variation12(:,abs(diff(OC.position(MI_variation12),1,1))<DISTANCE_LIMIT);
                    MI_variation34=combvec(MI{(iT-1)*T+it+2},MI{(iT-1)*T+it+3});
                    MI_variation34=MI_variation34(:,abs(diff(OC.position(MI_variation34),1,1))<DISTANCE_LIMIT);
                    MI_variation56=combvec(MI{(iT-1)*T+it+4},MI{(iT-1)*T+it+5});
                    MI_variation56=MI_variation56(:,abs(diff(OC.position(MI_variation56),1,1))<DISTANCE_LIMIT);
                    MI_variation78=combvec(MI{(iT-1)*T+it+6},MI{(iT-1)*T+it+7});
                    MI_variation78=MI_variation78(:,abs(diff(OC.position(MI_variation78),1,1))<DISTANCE_LIMIT);

                    % MI allowed for 4-frame-trajectories 
                    MI_variation14=combvec(MI_variation12,MI_variation34);
                    MI_variation14 = MI_variation14(:,abs(diff(OC.position(MI_variation14(end-2:end-1,:)),1,1))<DISTANCE_LIMIT);
                    MI_variation58=combvec(MI_variation56,MI_variation78);
                    MI_variation58 = MI_variation58(:,abs(diff(OC.position(MI_variation58(end-2:end-1,:)),1,1))<DISTANCE_LIMIT);

                    % MI allowed for 8-frame-trajectories 
                    MI_variation0=combvec(MI_variation14,MI_variation58);
                    MI_variation0 = MI_variation0(:,abs(diff(OC.position(MI_variation0(end-2:end-1,:)),1,1))<DISTANCE_LIMIT);
                    
                    % diffusivity (VAR_Dx) for 8-frame-trajectories 
                    Dx = diff(OC.position(MI_variation0),1,1);
                    SUM_Dx20 = sum(Dx.^2,1);
                    VAR_Dx = SUM_Dx20/(number_of_frames-1);
                    
                    % filter 8-frame-trajectories with too high diffusivities 
                    ff = VAR_Dx < VAR_Dx8_threshold;
                    if sum(ff)>0
                        MI_variation0 = MI_variation0(:,ff);
                        SUM_Dx20 = SUM_Dx20(ff);
                        
                        m = OC.m0(MI_variation0);
                        SUM_m0 = sum(m,1);
                        SUM_m20 = sum(m.^2,1);
                        MEAN_m = SUM_m0/8;
                        MEAN_m8_threshold = mean(OC.MEAN_m8_threshold(MI_variation0),1);
                        
                        % filter 8-frame-trajectories with too low mean of zero order intesity moments
                        ff = MEAN_m < MEAN_m8_threshold;
                        if sum(ff)>0
                            MI_variation0 = MI_variation0(:,ff);
                            SUM_Dx20 = SUM_Dx20(ff);
                            SUM_m0 = SUM_m0(ff);
                            SUM_m20 = SUM_m20(ff);
                        
                            VAR_m = (SUM_m20 - SUM_m0.^2/number_of_frames)/number_of_frames;
                            VAR_m8_threshold = mean(OC.VAR_m8_threshold(MI_variation0),1);

                            % filter 8-frame-trajectories with too high STD of zero order intesity moments
                            ff = VAR_m < VAR_m8_threshold;
                            if sum(ff)>0

                                % properties of filtered 8-frame-trajectories
                                MI_variation{it0} = MI_variation0(:,ff);
                                SUM_Dx2{it0} = SUM_Dx20(ff);
                                SUM_m{it0} = SUM_m0(ff);
                                SUM_m2{it0} = SUM_m20(ff);
                                SUM_m_withoutend{it0} = SUM_m{it0} - OC.m0(MI_variation{it0}(end,:))';
                                SUM_m2_withoutend{it0} = SUM_m2{it0} - (OC.m0(MI_variation{it0}(end,:))').^2;
                            end
                        end
                    end
                end
                
                %% 15-frame-trajectories
                len=1;
                number_of_frames = 15;
                for it=1:length(MI_variation)-1
                    
                    if size(MI_variation{it},2)>0 & size(MI_variation{it+1},2)>0

                        % MI of 15-frame-trajectories created out of two 8-frame-trajectories
                        ID0=combvec(1:size(MI_variation{it},2),1:size(MI_variation{it+1},2));
                        MI_variation0=combvec(MI_variation{it},MI_variation{it+1});
                        ff=MI_variation0(8,:)==MI_variation0(9,:);
                        MI_variation0=MI_variation0(:,ff);
                        MI_variation0(8,:)=[];
                        ID0=ID0(:,ff);
                        
                        % properties of 15-frame-trajectories
                        VAR_Dx = (SUM_Dx2{it}(ID0(1,:)) + SUM_Dx2{it+1}(ID0(2,:)))/(number_of_frames-1);
                        MEAN_m = (SUM_m_withoutend{it}(ID0(1,:)) + SUM_m{it+1}(ID0(2,:)))/number_of_frames;
                        VAR_m = (SUM_m2_withoutend{it}(ID0(1,:)) + SUM_m2{it+1}(ID0(2,:)) - ...
                            (SUM_m_withoutend{it}(ID0(1,:)) + SUM_m{it+1}(ID0(2,:))).^2/number_of_frames)/number_of_frames;
                        MEAN_m15_threshold = mean(OC.MEAN_m15_threshold(MI_variation0),1);
                        VAR_m15_threshold = mean(OC.VAR_m15_threshold(MI_variation0),1);
                        
                        % filter 15-frame-trajectories that have too low mean or too high diffusivity 
                        ff = MEAN_m < MEAN_m15_threshold  & VAR_Dx < VAR_Dx15_threshold;% & VAR_m < VAR_m15_threshold

                         ID{len}{it} = [];
                         cost{len}{it} = [];
                         cost0=VAR_Dx(ff).*VAR_m(ff); %cost function
                         MI_variation0=MI_variation0(:,ff);
                         ID0=ID0(:,ff);

                         % from trajectories that have the same minima at the beggining and at the end, select the ones with minimal cost
                         [a1,b1,c1]=unique(MI_variation0(1,:));
                         [a2,b2,c2]=unique(MI_variation0(end,:));
                          for i=1:length(a1)
                                    for j=1:length(a2)
                                        f=find(c1==i & c2==j);
                                        [d,e]=min(cost0(f));
                                        cost{len}{it}=[cost{len}{it},d];
                                        ID{len}{it}=[ID{len}{it},ID0(:,f(e))];
                                    end
                           end
                    end
                end
                
              %% [(len+1)*7+1]-frame-trajectories  
              while len<T/7 && length(ID)==len

                	len = len+1;
                    number_of_frames=(len+1)*7+1;

                for it=1:min([size(ID{len-1},2),size(ID{1},2)-(len-1)])
                    
                    if size(ID{len-1}{it},2)>0 && size(ID{1}{it+len-1},2)>0

                        % MI of [(len+1)*7+1]-frame-trajectories created out of one [(len)*7+1]-frame trajetory and 8-frame-trajectories
                        ID0=combvec(ID{len-1}{it},ID{1}{it+len-1});
                        ff=ID0(end-2,:)==ID0(end-1,:);
                        ID0=ID0(:,ff);
                        ID0(end-1,:)=[];
                        
                        if len>2
                            %delete shorter trajectories that are contained in [(len+1)*7+1]-frame-trajectories
                            IDD0=combvec(1:size(ID{len-1}{it},2),1:size(ID{1}{it+len-1},2));
                            IDD0=unique(IDD0(1,ff));
                            ID{len-1}{it}(:,IDD0)=[];
                            cost{len-1}{it}(:,IDD0)=[];
                        end
                        
                      if sum(ff)>0  
                          % properties/cost function of [(len+1)*7+1]-frame-trajectories
                          SUM_Dx20=0;
                          SUM_m20=0;
                          SUM_m0=0;
                          for i=1:len
                              SUM_Dx20 = SUM_Dx20 + SUM_Dx2{it+i-1}(ID0(i,:));
                              SUM_m20 = SUM_m20 + SUM_m2_withoutend{it+i-1}(ID0(i,:));
                              SUM_m0 = SUM_m0 + SUM_m_withoutend{it+i-1}(ID0(i,:));
                          end
                            i=len+1;
                            SUM_Dx20 = SUM_Dx20 + SUM_Dx2{it+i-1}(ID0(i,:));
                              SUM_m20 = SUM_m20 + SUM_m2{it+i-1}(ID0(i,:));
                              SUM_m0 = SUM_m0 + SUM_m{it+i-1}(ID0(i,:));
                          
                        VAR_Dx = SUM_Dx20/(number_of_frames-1);
                        VAR_m = (SUM_m20 - ...
                            (SUM_m0).^2/number_of_frames)/number_of_frames;

                        cost0  = VAR_Dx.*VAR_m;
                        ID{len}{it}=[];
                        cost{len}{it}=[];

                        % from trajectories that have the same minima at the end, select the ones with minimal cost
                        [a2,b2,c2]=unique(ID0(end,:));
                        for j=1:length(a2)
                              f=find(c2==j);
                              [d,e]=min(cost0(f));
                              cost{len}{it}=[cost{len}{it},d];
                              ID{len}{it}=[ID{len}{it},ID0(:,f(e))];
                        end
                                
                        
                      end
                    end
                end
              end
              
              %% select the longest trajectories with minimal cost
                MI_decided_par{iT}=[];
                MI_trajectory_par{iT}=[];
                kk=0;
                for len=1:length(ID)
                    ID0=ID{length(ID)-len+1};
                    cost0=cost{length(ID)-len+1};
                    cost0=cell2mat(cost0);
                    MI_trajectory0=[];
                    for it=1:length(ID0)
                        if size(ID0{it},2)>0
                          for i=1:size(ID0{it},2)
                            MI_trajectory00=[];
                            for j=1:size(ID0{it},1)-1
                                MI_trajectory00=[MI_trajectory00;MI_variation{it+j-1}(1:end-1,ID0{it}(j,i))];
                            end
                            MI_trajectory00=[MI_trajectory00;MI_variation{it+j}(:,ID0{it}(j+1,i))];
                            MI_trajectory0=[MI_trajectory0,MI_trajectory00];
                          end
                        end
                    end
                    
                    [cost0,a]=sort(cost0);
                    MI_trajectory0=MI_trajectory0(:,a);
                    for i=1:size(MI_trajectory0,2)
                        if length(intersect(MI_trajectory0(:,i),MI_decided_par{iT}))==0
                            MI_decided_par{iT}=[MI_decided_par{iT};MI_trajectory0(:,i)];
                            kk=kk+1;
                            MI_trajectory_par{iT}{kk}=MI_trajectory0(:,i);
                        end
                    end
                end
        
end
           
%% combine trajectories found from all the cores
for iT=1:length(MI_decided_par)-1
    a=intersect(MI_decided_par{iT},MI_decided_par{iT+1});
    MI_trajectory_pair=zeros(2,length(a));
    for i=1:length(a)
        for kk=1:length(MI_trajectory_par{iT})
            if sum(MI_trajectory_par{iT}{kk}==a(i))>0
                MI_trajectory_pair(1,i)=kk;
            end
        end
        for kk=1:length(MI_trajectory_par{iT+1})
            if sum(MI_trajectory_par{iT+1}{kk}==a(i))>0
                MI_trajectory_pair(2,i)=kk;
            end
        end
    end
        MI_trajectory_pair=MI_trajectory_pair(:,MI_trajectory_pair(1,:)>0 & MI_trajectory_pair(2,:)>0);
        [c,ia,ib]=unique(MI_trajectory_pair','rows');
        for i=1:length(ia)
            if length(ib==ia(i))>=4
                ff1=OC.timeFrame(MI_trajectory_par{iT+1}{MI_trajectory_pair(2,ia(i))}(1)); %first frame of second
                ff=OC.timeFrame(MI_trajectory_par{iT}{MI_trajectory_pair(1,ia(i))})<ff1; %frame of the first that are lower than ff1
                MI_trajectory_par{iT+1}{MI_trajectory_pair(2,ia(i))}=[MI_trajectory_par{iT}{MI_trajectory_pair(1,ia(i))}(ff);MI_trajectory_par{iT+1}{MI_trajectory_pair(2,ia(i))}];
                MI_trajectory_par{iT}{MI_trajectory_pair(1,ia(i))}=[];
            end
        end

end
    
MI_decided=MI_decided_par{1};
MI_trajectory=MI_trajectory_par{1};
for iT=2:length(MI_decided_par)
    MI_decided=[MI_decided; MI_decided_par{iT}];
    MI_trajectory=[MI_trajectory MI_trajectory_par{iT}];
end
MI_decided=unique(MI_decided); 
 