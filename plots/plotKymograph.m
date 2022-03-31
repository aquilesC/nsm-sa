%the script plots the kymograph, saved in the file ExperimentTimeStamp_C.m
%and adds the particle trajectory, saved in the file ExperimentTimeStamp_D.m

%% plot kymograph
figure('Position',[1 200 2000 500]);
%surf(1:size(data.Im,1),1:size(data.Im,2),data.Im'); shading flat; view(2);  hold on; %colorbar;
imagesc(1:size(data.Im,1),1:size(data.Im,2),data.Im');
colormap(bone)
ylim([1 size(data.Im,2)])
xlim([1 size(data.Im,1)])
ylabel('Pixel')
xlabel('Time frame')
caxis([1-3e-4 1+3e-4])
box off
hold on


%% plot trajectory
fun_label = 1; %1 - iOC and D, 2 - MW and HR
if exist('trajectory') == 1
    
    for itra=1:length(trajectory)
         plot3(trajectory(itra).timeFrame,trajectory(itra).position,2*ones(size(trajectory(itra).timeFrame)),...
             'Color',BasicColor(mod(itra,8)+1),'Marker','.')
         
         if fun_label==1
             text(trajectory(itra).timeFrame(1),trajectory(itra).position(1),2,...
                     {strcat('iOC=',num2str(-trajectory(itra).iOC_mean*1e3)),strcat('D=',num2str(trajectory(itra).Deff_mean))},...
                     'Color',BasicColor(mod(itra,8)+1), 'FontSize',14);%,...
                     %'BackgroundColor',[1 1 1],'EdgeColor',[1 1 1])
         end
    end

end

%% add secondary axis
% ax1=gca;
% ax2=axes('Position', ax1.Position,'XAxisLocation', 'top','YAxisLocation','right','color','none');
% line([data.time(1), data.time(end)],[data.Yum(1), data.Yum(end)], 'Color','none')
% xlim([data.time(1), data.time(end)])
% ylim([data.Yum(1), data.Yum(end)])
% xlabel('Time (s)')
% ylabel('Position (/mum)')



