function [dIm, x3, kk]=image_stabilization (Y,Im,Im0,Yav)

%dIm is a ratio between Im and Im0 which is stabilized in Y position
%and normalized by fitted polynom 

funplot = 0;
dx = 1e-2;

threshold = 1e-4;


    
%         x0 = 0;
%         dIm0 = Im./Im0;
%         pf=polyfit(Y,dIm0,porder);
%         dIm0=dIm0./polyval(pf,Y);
%         dImsum_start=sum(abs(dIm0-1));
        
        x1 = 0;%-dx;
        Im_temp=interp1(Y,Im,Y+x1,'spline');
        dIm1 = Im_temp./Im0;
        dIm1=dIm1./smooth(dIm1,Yav)';
        
        x2 = dx;
        Im_temp=interp1(Y,Im,Y+x2,'spline');
        dIm2 = Im_temp./Im0;
        dIm2=dIm2./smooth(dIm2,Yav)';
        
        if funplot == 1
            hold off
            plot(dIm1); hold on
            plot(dIm2);
        end
        
        kk=1;
        while abs(x2-x1)>threshold

            %x4 = (x1^2*(y2 - y3) + x2^2 * (y3 - y1) + x3^2 * (y1 - y2))/(2* (x1 *(y2 - y3) + x2 *(y3 - y1) + x3* (y1 - y2)));
            %x4 = (x1^2*(dIm2 - dIm3) + x2^2 * (dIm3 - dIm1) + x3^2 * (dIm1 - dIm2))./(2* (x1 *(dIm2 - dIm3) + x2 *(dIm3 - dIm1) + x3* (dIm1 - dIm2)));
            
%             x3 = -(x2-x1)./(dIm2(1+Yav/2:end-Yav/2)-dIm1(1+Yav/2:end-Yav/2)).*(dIm1(1+Yav/2:end-Yav/2)-1)+x1;
%             W = abs(dIm2(1+Yav/2:end-Yav/2)-dIm1(1+Yav/2:end-Yav/2)).^2;
            
            x3 = -(x2-x1)./(dIm2(2:end-1)-dIm1(2:end-1)+eps).*(dIm1(2:end-1)-1)+x1;
            %W = abs(dIm2-dIm1).^2.*(Im/max(Im)).^2;
            W = abs(dIm2(2:end-1)-dIm1(2:end-1)).*(Im(2:end-1)/max(Im));
            
            ff=isnan(x3)==0 & abs(x3)<Inf & isnan(W)==0 & abs(W)<Inf;
            x3 = sum(x3(ff).*W(ff))/sum(W(ff));
            Im_temp=interp1(Y,Im,Y+x3,'spline');
%             if length(find(isnan(Im_temp)==1))>0
%                 disp('');
%             end
            dIm = Im_temp./Im0;
            
            I_mean=smooth(dIm,Yav)';
            dIm=dIm./I_mean;
            if funplot == 1
                plot(dIm,'.'); hold on
            end
%             if length(find(isnan(dIm)==1))>0
%                 disp('');
%             end
           
            %y3=std(dIm);
%             plot(x4,y4,'.')

            x1=x2; x2=x3; 
            %y1=y2; y2=y3; 
            dIm1=dIm2; dIm2=dIm; 

            kk=kk+1;
            if kk>100
                disp('image_stabilizationY8: does not converge, kk>100'); return
            end
        end
   
   
%         dImsum_stop=sum(abs(dIm-1));
%         if dImsum_stop>dImsum_start
%             disp('image_stabilizationY2: higher noise at the end!'); return
%         end
        

