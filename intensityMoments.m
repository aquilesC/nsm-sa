function [opticalContrast]=intensityMoments(data,Yav)
 
X=1:size(data.Im,2);
t=1:size(data.Im,1);
X2=repmat(X,length(t),1);
t2=repmat(t',1,length(X));
Y=data.Im-1;
 
   tic 
    m0=conv2(Y,ones(1,Yav),'same');
    m1=conv2(Y.*X,ones(1,Yav),'same')./m0;
    %m2=conv2(Y,[-(Yav-1)/2:(Yav-1)/2].^2,'same')./m0;
%     for i=1:(Yav2-1)/2
%             OC(:,:,i)=conv2(Y,ones(1,2*i-1),'same');
%     end
%     OC(:,:,i+1)=m0;
    
    m1_difference=abs(X-m1);
 
    
    A =  -Y;
    Mask=ones(1,Yav); Mask((Yav+1)/2)=0;
    B = imdilate(A,Mask);
    dips = A > B;
    dips_position=X2(dips);
    opticalContrast.timeFrame=t2(dips);
    
%     A =  Y;
%     Mask=ones(1,Yav); Mask((Yav+1)/2)=0;
%     B = imdilate(A,Mask);
%     peaks = A > B;
%     peaks_position = X2(peaks);
%     peaks_timeFrame = t2(peaks);
% 
%     cross = Y(:,1:end-1).*Y(:,2:end)<=0;
%     cross_position = X2(cross);
%     cross_timeFrame = t2(cross);
%     
%     cross_position=[cross_position; peaks_position];
%     cross_timeFrame=[cross_timeFrame; peaks_timeFrame];
 
    %refinement
    for i=1:length(dips_position)
        
        ff=max([dips_position(i)-(Yav-1)/2,1]):min([dips_position(i)+(Yav-1)/2,size(Y,2)]);
        [a,b]=min(m1_difference(opticalContrast.timeFrame(i),ff));
        opticalContrast.position(i,:)=m1(opticalContrast.timeFrame(i),ff(b));
%         opticalContrast.OC(:,i)=OC(opticalContrast.timeFrame(i),ff(b),:);
        %opticalContrast.m2(i,:)=m2(opticalContrast.timeFrame(i),ff(b));
        opticalContrast.m0(i,:)=m0(opticalContrast.timeFrame(i),ff(b));
        
%         cross_position0=cross_position(cross_timeFrame==opticalContrast.timeFrame(i));
%         cross_position0_right=cross_position0(cross_position0>=opticalContrast.position(i));
%         cross_position0_left=cross_position0(cross_position0<opticalContrast.position(i));
%         if length(cross_position0_right)>0
%             [opticalContrast.position_right(i),b]=min(cross_position0_right);
%         else
%             opticalContrast.position_right(i)=X(end);
%         end
%         if length(cross_position0_left)>0
%             [opticalContrast.position_left(i),b]=max(cross_position0_left);
%         else
%             opticalContrast.position_left(i)=X(1);
%         end
        
%         if opticalContrast.position_left(i)==opticalContrast.position_right(i)
%             opticalContrast.iOC(i)=0;
%         else
%             opticalContrast.iOC(i)=trapz(data.Yum(opticalContrast.position_left(i):opticalContrast.position_right(i)),Y(opticalContrast.timeFrame(i),opticalContrast.position_left(i):opticalContrast.position_right(i)));
%         end    
        
%         if fun_plot==1 %& opticalContrast.position(i)>50
%             hold off
%             plot(Y(opticalContrast.timeFrame(i),:)); hold on
%             plot(opticalContrast.position(i),opticalContrast.OC(1,i),'.')
%             plot(opticalContrast.position_right(i),0,'o');
%             plot(opticalContrast.position_left(i),0,'*');
%         end
    end
    
    opticalContrast.positionUm=(opticalContrast.position-1)*(data.Yum(2)-data.Yum(1))+data.Yum(1); 
              
    
%     for i=1:size(Y,1)
%         hold off
%         plot(Y(i,:)); hold on
%         ff=find(opticalContrast_all.timeFrame==i);
%         plot(opticalContrast_all.position0(ff),Y(i,opticalContrast_all.position0(ff)),'o')
%         plot(opticalContrast_all.position(ff),opticalContrast_all.OC(1,ff),'x')
%     end
 

