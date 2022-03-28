function opticalContrast=integrateOpticalContrast1(data)

%calculates the integrated optical contrast for the all found maxima and
%minima in the kymograph and returns their position in time and space

X=data.Yum; %space axis
Y=data.Im; % matrix of intensities
XPX=1:length(data.Yum);
%Yav %span of the moving average
 
% opticalContrast.timeFrame=[]; % serie of number of the time frames correspondng to found maximas or minimas
% opticalContrast.positionUm=[]; % corresponding serie of positions 
% opticalContrast.iOC=[]; %corresponding serie of integrated optical contrast
% opticalContrast.edge_left=[];
% opticalContrast.edge_right=[];
% opticalContrast.length=[];

dY=(Y(:,1:end-1)-1).*(Y(:,2:end)-1);
%dY(:,1)=0; dY(:,end)=0;
[edge,timeFrame]=find(dY'<=0);
dtimeFrame=diff(timeFrame);
a=find(dtimeFrame==0);
opticalContrast.edgePX_left=edge(a);
opticalContrast.edgePX_right=edge(a+1);
opticalContrast.lengthPX=opticalContrast.edgePX_right-opticalContrast.edgePX_left;
opticalContrast.timeFrame=timeFrame(a);
opticalContrast.iOCpx=zeros(size(opticalContrast.timeFrame));

for i=1:length(opticalContrast.timeFrame)
    x=XPX(opticalContrast.edgePX_left(i):opticalContrast.edgePX_right(i));
    y=Y(opticalContrast.timeFrame(i),opticalContrast.edgePX_left(i):opticalContrast.edgePX_right(i))-1;
    opticalContrast.iOCpx(i,1)=trapz(x,y);
    opticalContrast.position(i,1)=trapz(x,x.*y)./opticalContrast.iOCpx(i,1);
    
end

opticalContrast.edge_left=(opticalContrast.edgePX_left+1)*(data.Yum(2)-data.Yum(1)) - data.Yum(1);
opticalContrast.edge_right=(opticalContrast.edgePX_right+1)*(data.Yum(2)-data.Yum(1)) - data.Yum(1);
opticalContrast.positionUm=opticalContrast.position*(data.Yum(2)-data.Yum(1));
%opticalContrast.length=opticalContrast.lengthPX*(data.Yum(2)-data.Yum(1));
opticalContrast.iOC=opticalContrast.iOCpx*(data.Yum(2)-data.Yum(1));

% for i=1:size(opticalContrast.timeFrame)
%     hold off
%     plot(data.Yum,data.Im(opticalContrast.timeFrame(i),:)); hold on
%     plot(opticalContrast.edge_right(i),1,'.'); hold on
%     plot(opticalContrast.edge_left(i),1,'o');
% end
        

% for it=1:size(data.Im,1)
%    
%        
%         %Ys=(smooth(Y(it,:),Yav,'lowess'))';
%         Ys=Y(it,:);
%         %dYs=Ys(2:end)-Ys(1:end-1);
%         ff=find((Ys(1:end-1)-1).*(Ys(2:end)-1)<=0);
%         
% %         hold off
% %         plot(X,Ys-1); hold on
% %         plot(X,zeros(size(X)))
% %         plot(X(ff), Ys(ff)-1,'.','MarkerSize',20)
%         
%         iOC=zeros(1,length(ff)-1);
%         positionUm=zeros(1,length(ff)-1);
%         edge_left=zeros(1,length(ff)-1);
%         edge_right=zeros(1,length(ff)-1);
% 
%         for i=1:length(ff)-1
%             iOC(i)=trapz(X(ff(i):ff(i+1)),Y(it,ff(i):ff(i+1))-1);
%             positionUm(i)=trapz(X(ff(i):ff(i+1)),X(ff(i):ff(i+1)).*(Y(it,ff(i):ff(i+1))-1))./iOC(i);
%             edge_left(i)=X(ff(i));
%             edge_right(i)=X(ff(i+1));
%         end
%         opticalContrast.timeFrame=[opticalContrast.timeFrame, repmat(it,1,length(ff)-1)];
%         opticalContrast.iOC=[opticalContrast.iOC, iOC];
%         opticalContrast.positionUm=[opticalContrast.positionUm,positionUm];
%         opticalContrast.edge_left=[opticalContrast.edge_left, edge_left];
%         opticalContrast.edge_right=[opticalContrast.edge_right, edge_right];
%         opticalContrast.length=[opticalContrast.length, edge_right-edge_left];
%             
%          %plot(positionUm,iOC,'o','MarkerFaceColor','red')
% end

% for it=1:size(data.Im,1)
%    
%        
%         %Ys=(smooth(Y(it,:),Yav,'lowess'))';
%         Ys=Y(it,:);
%         dYs=Ys(2:end)-Ys(1:end-1);
%         ff_zero=find((Ys(1:end-1)-1).*(Ys(2:end)-1)<=0);
%         
% %         hold off
% %         plot(X,Ys-1); hold on
% %         plot(X,zeros(size(X)))
% %         plot(X(ff_zero), Ys(ff_zero)-1,'.','MarkerSize',20)
%         
%         ff_extreme=[1,1,find(dYs(2:end).*dYs(1:end-1)<0)+1,length(Ys),length(Ys)];
%         
% %         plot(X(ff_extreme), Ys(ff_extreme)-1,'x')
%         
%         A=(Ys(ff_extreme(1:end-1))-1).*(Ys(ff_extreme(2:end))-1); %extrmes which does not cross zero
%         B=abs(Ys(ff_extreme(1:end-1))-1)-abs(Ys(ff_extreme(2:end))-1); %from the couple those who are closer to the zero
%         ff_extreme1=find(A>0 & B<0);
%         
% %         plot(X(ff_extreme(ff_extreme1)), Ys(ff_extreme(ff_extreme1))-1,'o')
%         
%         A=(Ys(ff_extreme(ff_extreme1))-1)./(Ys(ff_extreme(ff_extreme1-1))-1); %enough contrast
%         B=(Ys(ff_extreme(ff_extreme1))-1)./(Ys(ff_extreme(ff_extreme1+1))-1);
%         ff_extreme2 = find(A < 0.74 & B<0.74);
%         
% %         plot(X(ff_extreme(ff_extreme1(ff_extreme2))), Ys(ff_extreme(ff_extreme1(ff_extreme2)))-1,'o','MarkerFaceColor',[0 0 0])
%         
%         ff=[ff_zero,(ff_extreme(ff_extreme1(ff_extreme2)))];
%         ff=sort(unique(ff)); %bounderies for trapz
%         
% %         plot(X(ff),Ys(ff)-1,'o','MarkerFaceColor',[0 0 0])
%         
%         iOC=zeros(1,length(ff)-1);
%         positionUm=zeros(1,length(ff)-1);
%         edge_left=zeros(1,length(ff)-1);
%         edge_right=zeros(1,length(ff)-1);
% 
%         for i=1:length(ff)-1
%             iOC(i)=trapz(X(ff(i):ff(i+1)),Y(it,ff(i):ff(i+1))-1);
%             positionUm(i)=trapz(X(ff(i):ff(i+1)),X(ff(i):ff(i+1)).*(Y(it,ff(i):ff(i+1))-1))./iOC(i);
%             edge_left(i)=X(ff(i));
%             edge_right(i)=X(ff(i+1));
%         end
%         opticalContrast.timeFrame=[opticalContrast.timeFrame, repmat(it,1,length(ff)-1)];
%         opticalContrast.iOC=[opticalContrast.iOC, iOC];
%         opticalContrast.positionUm=[opticalContrast.positionUm,positionUm];
%         opticalContrast.edge_left=[opticalContrast.edge_left, edge_left];
%         opticalContrast.edge_right=[opticalContrast.edge_right, edge_right];
%             
% %         plot(positionUm,iOC,'o','MarkerFaceColor','red')
% end
 


