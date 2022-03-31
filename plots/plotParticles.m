% the script plots the result of the particle trajectory analysis, saved in the file ExperimentTimeStamp_C.m

%% defined values 
n = [trajectory.N]; %temporal length of trajectories
iOC = [trajectory.iOC_mean]; %integrated optical contrast of particles
D = [trajectory.Deff_mean]; %diffusivity of particles

%% filter (relevant only for SA)
a=n>=40;
n=n(a);
iOC=iOC(a);
D=D(a);

%% transform from iOC to MW and transform from D to HR
channel_name='channel_30x85';
nanochannelArea = list_of_nanochannels (channel_name);
HR = EffectiveDiffusivityToSize (D, nanochannelArea); %hydrodynamic radius
[MW, calibration] = iOCToWeight (iOC, nanochannelArea); %molecular weight
MW = abs(MW);
 
%% change units to iOC [nm], D[um2/s], MW [kDa], HR [nm]
iOC = iOC*1000;
MW = MW/1000;
HR = HR*1000;

%% iOC histogram
[N, edges] = weighted_histcounts (iOC,n,100);
figure;
bar(edges,N);
xlabel('iOC (nm)')
ylabel('Counts')

%% Deff histogram
[N, edges] = weighted_histcounts (D,n,100);
figure;
bar(edges,N);
xlabel('Diffusivity (\mum^2/s)')
ylabel('Counts')


%% MW histogram
[N, edges] = weighted_histcounts (MW,n,100);
figure;
bar(edges,N,...
    'FaceColor',[0.8 0.8 0.8],...
    'EdgeColor','none'); hold on
xlabel('Molecular weight (kDa)')
ylabel('Counts')
xlim([edges(1) edges(end)])


%% HR histogram
[N, edges] = weighted_histcounts (HR,n,100);
figure;
bar(edges,N,...
    'FaceColor',[0.8 0.8 0.8],...
    'EdgeColor','none'); hold on
xlabel('Hydrodynamic radius (nm)')
ylabel('Counts')
xlim([0, edges(end)])

if exist('number_of_peaks')==1
    for i=1:number_of_peaks
        a=MW>peak_MWxHR(i).X_mean-3*peak_MWxHR(i).X_std & MW<peak_MWxHR(i).X_mean+3*peak_MWxHR(i).X_std...
            & HR>peak_MWxHR(i).Y_mean-3*peak_MWxHR(i).Y_std & HR<peak_MWxHR(i).Y_mean+3*peak_MWxHR(i).Y_std;
        [N, edges] = weighted_histcounts (HR(a),n(a),HR_edges);
        bar(edges,N,...
            'FaceColor',BasicColor(i),...
            'EdgeColor','none');
        plot(edges,peak_MWxHR(i).Y_N*exp(-0.5*((edges-peak_MWxHR(i).Y_mean)/(peak_MWxHR(i).Y_std)).^2),...
          'Color',BasicColor(i)) 
       text(peak_MWxHR(i).Y_mean,peak_MWxHR(i).Y_N,...
           strcat(num2str(peak_MWxHR(i).Y_mean)),...
            'Color',BasicColor(i));
    end
end

%% scatter plot MW x HR
figure;
for i=1:length(MW)
                scatter(MW(i),HR(i),'MarkerFaceAlpha',n(i)/max(n),...
                    'MarkerFaceColor','black',...
                    'MarkerEdgeColor','none'); hold on
end
xlabel('Molecular weight (kDa)')
ylabel('Hydrodynamic radius (nm)')

%% scatter plot MW x HR
figure;
for i=1:length(MW)
                scatter(iOC(i),D(i),'MarkerFaceAlpha',n(i)/max(n),...
                    'MarkerFaceColor','black',...
                    'MarkerEdgeColor','none'); hold on
end
xlabel('iOC (nm)')
ylabel('Diffusivity (\mum^2/s)')





