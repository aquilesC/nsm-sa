%% defined values
n = collection.N;
iOC = collection.iOC_mean;
D = collection.Deff_mean;

channel_name='channel_30x85_AZ';

%% filter
a=n>40;
n=n(a);
iOC=iOC(a);
D=D(a);

%% transform from iOC to MW and transform from Deff to HR
nanochannelArea = list_of_nanochannels (channel_name);
HR = EffectiveDiffusivityToSize (D, nanochannelArea);
[MW, calibration] = iOCToWeight (-iOC, nanochannelArea);
 
%% change units to iOC [nm], D[um2/s], MW [kDa], HR [nm]
iOC = iOC*1000;
MW = MW/1000;
HR = HR*1000;

%% define edges
MW_step=5; %[kDa]
HR_step=0.1; %[nm]
MW_edges=0:MW_step:max(MW);
HR_edges=0:HR_step:max(HR);

%% iOC histogram
[N, edges] = weighted_histcounts (iOC,n,1000);
figure;
bar(edges,N);
xlabel('iOC (nm)')
ylabel('Counts')

%% Deff histogram
[N, edges] = weighted_histcounts (D,n,1000);
figure;
bar(edges,N);
xlabel('Diffusivity (\mum^2/s)')
ylabel('Counts')

%% iOC x Deff histogram
[N, edgesX, edgesY] = weighted_histcounts2 (D,iOC,n,500);
figure;
surf(edgesY,edgesX,N); view(2); shading flat
xlabel('iOC (nm)')
ylabel('Diffusivity (\mum^2/s)')

%% MW x HR histogram
%[N, edgesX, edgesY] = weighted_histcounts2 (HR,MW,n,300);

[N, edgesY, edgesX] = weighted_histcounts2 (HR,MW,n,HR_edges,MW_edges);
HR_globular = massToHR (edgesX*1e3,'globular')*1000;
figure;
surf(edgesX,edgesY,N); view(2); shading flat; hold on
plot3(edgesX,HR_globular,repmat(max(max(N))*2,size(edgesX)),'Color','white')
xlabel('Molecular weight (kDa)')
ylabel('Hydrodynamic radius (nm)')

%% definition of peaks
number_of_peaks = input('Number_of_peaks?');
for i=1:number_of_peaks
    peak_MWxHR(i) = find_peak_from_2D_histogram(N, edgesY, edgesX);
end

%% add MW x HR peaks
if exist('peak_MWxHR')==1
    for i=1:length(peak_MWxHR)
        theta = linspace(0,2*pi);
        x=peak_MWxHR(i).X_mean + cos(theta)*3*peak_MWxHR(i).X_std;
        y=peak_MWxHR(i).Y_mean + sin(theta)*3*peak_MWxHR(i).Y_std;
        plot3(x,y,repmat(max(max(N))*2,size(x)),'Color',BasicColor(i));
        text(peak_MWxHR(i).X_mean,peak_MWxHR(i).Y_mean,...
            strcat('[',num2str(peak_MWxHR(i).X_mean),',',num2str(peak_MWxHR(i).Y_mean),']'),...
            'Color',BasicColor(i));
    end
end

%% add iOC and D axis to MW and HR histogram
ax1=gca;
box off
hold on

YTickLabel0=fliplr([0.2,0.5,1,2,5,10,20,50]);
YTickLabel=[];
for i=1:length(YTickLabel0)
    YTickLabel{i}=YTickLabel0(i);
end
YTick = EffectiveDiffusivityToSize (YTickLabel0, nanochannelArea);

ax2=axes('Position', ax1.Position,'XAxisLocation', 'top','YAxisLocation','right','color','none',...
    'YTick',YTick*1000,'YTickLabel',YTickLabel);
line(-[ax1.XLim(1), ax1.XLim(2)]/calibration*1e6,[ax1.YLim(1), ax1.YLim(2)], 'Color','none')
xlim(-[ax1.XLim(1), ax1.XLim(2)]/calibration*1e6)
ylim([ax1.YLim(1), ax1.YLim(2)])
xlabel('iOC (nm)')
ylabel('Diffusivity (\mum^2/s)')


%% MW histogram
%[N, edges] = weighted_histcounts (MW,n,1000);
[N, edges] = weighted_histcounts (MW,n,MW_edges);
figure;
bar(edges,N,...
    'FaceColor',[0.8 0.8 0.8],...
    'EdgeColor','none'); hold on
xlabel('Molecular weight (kDa)')
ylabel('Counts')

if exist('peak_MWxHR')==1
    for i=1:number_of_peaks
        a=MW>peak_MWxHR(i).X_mean-3*peak_MWxHR(i).X_std & MW<peak_MWxHR(i).X_mean+3*peak_MWxHR(i).X_std...
            & HR>peak_MWxHR(i).Y_mean-3*peak_MWxHR(i).Y_std & HR<peak_MWxHR(i).Y_mean+3*peak_MWxHR(i).Y_std;
        [N, edges] = weighted_histcounts (MW(a),n(a),MW_edges);
        bar(edges,N,...
            'FaceColor',BasicColor(i),...
            'EdgeColor','none');
        plot(edges,peak_MWxHR(i).X_N*exp(-((edges-peak_MWxHR(i).X_mean)/(2*peak_MWxHR(i).X_std)).^2),...
          'Color',BasicColor(i)) 
       text(peak_MWxHR(i).X_mean,peak_MWxHR(i).X_N,...
           strcat(num2str(peak_MWxHR(i).X_mean)),...
            'Color',BasicColor(i));
    end
end

%% add iOC axis to MW histogram
ax1=gca;
box off
hold on

ax2=axes('Position', ax1.Position,'XAxisLocation', 'top','YAxisLocation','right','color','none',...
    'YTick',[]);
line(-[ax1.XLim(1), ax1.XLim(2)]/calibration*1e6,[ax1.YLim(1), ax1.YLim(2)], 'Color','none')
xlim(-[ax1.XLim(1), ax1.XLim(2)]/calibration*1e6)
ylim([ax1.YLim(1), ax1.YLim(2)])
xlabel('iOC (nm)')



%% HR histogram
%[N, edges] = weighted_histcounts (HR,n,1000);
[N, edges] = weighted_histcounts (HR,n,HR_edges);
figure;
bar(edges,N,...
    'FaceColor',[0.8 0.8 0.8],...
    'EdgeColor','none'); hold on
xlabel('Hydrodynamic radius (nm)')
ylabel('Counts')

if exist('number_of_peaks')==1
    for i=1:number_of_peaks
        a=MW>peak_MWxHR(i).X_mean-3*peak_MWxHR(i).X_std & MW<peak_MWxHR(i).X_mean+3*peak_MWxHR(i).X_std...
            & HR>peak_MWxHR(i).Y_mean-3*peak_MWxHR(i).Y_std & HR<peak_MWxHR(i).Y_mean+3*peak_MWxHR(i).Y_std;
        [N, edges] = weighted_histcounts (HR(a),n(a),HR_edges);
        bar(edges,N,...
            'FaceColor',BasicColor(i),...
            'EdgeColor','none');
        plot(edges,peak_MWxHR(i).Y_N*exp(-((edges-peak_MWxHR(i).Y_mean)/(2*peak_MWxHR(i).Y_std)).^2),...
          'Color',BasicColor(i)) 
       text(peak_MWxHR(i).Y_mean,peak_MWxHR(i).Y_N,...
           strcat(num2str(peak_MWxHR(i).Y_mean)),...
            'Color',BasicColor(i));
    end
end


%% add D axis on HR histogram
box off
hold on

XTickLabel0=fliplr([0.2,0.5,1,2,5,10,20,50]);
XTickLabel=[];
for i=1:length(XTickLabel0)
    XTickLabel{i}=XTickLabel0(i);
end
XTick = EffectiveDiffusivityToSize (XTickLabel0, nanochannelArea);

ax1=gca;
ax2=axes('Position', ax1.Position,'XAxisLocation', 'top','YAxisLocation','right','color','none',...
    'XTick',XTick*1000,'XTickLabel',XTickLabel,'YTick',[]);
line([0, edges(end)*1000],[0 max(N)*1.2], 'Color','none')
xlim([0, edges(end)*1000])
ylim([0 max(N)*1.2])
xlabel('Diffusivity (\mum^2/s)')



