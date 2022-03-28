function [diffusion_coefficient, diffusion_coefficient_correction, diffusion_coefficient_error] = trajectoriesToDiffusivity(position, time, frame,tav, fun)

if length(position)>2
% Estimation of the diffusion coefficient based on Vestergaard, C. L.;  Blainey, P. C.; Flyvbjerg, H., Optimal estimation of diffusion coefficients from single-particle trajectories. Phys Rev E 2014, 89 (2).

% position - a time serie of measured positions of a particle
% time - corresponding time of the measurement
% frame - no. of the frame
% tav - time smoothen factor

% diffusion_coefficient - estimated diffusion coefficient corrected for the motion blur and localization error
% diffusion_coefficient_correction - correction term 
% diffusion_coefficient_error - estimated error of the value

    if length(position)/tav ~= round(length(position)/tav)
        %disp('trajectoriesToDiffusivity: length(position)/tav ~= round(length(position)/tav)');
        a=floor(length(position)/tav)*tav;
        position=position(1:a);
        time=time(1:a);
        frame=frame(1:a);
    end

    if tav>1
        position=reshape(position,tav,[])'; position=position(:);
        time=reshape(time,tav,[])'; time=time(:);
        frame=repmat((1:length(position)/tav)',1,tav); frame=frame(:);
    end
    
    Dtime=diff(time);
    Dframe=diff(frame);
    Dposition=diff(position);
    ff=find(Dframe==1);
    Dtime=Dtime(ff);
    Dposition=Dposition(ff);
    
    diffusion_coefficient0 = mean(Dposition.^2./(2*Dtime));
    diffusion_coefficient_correction = mean(Dposition(1:end-1).*Dposition(2:end)./Dtime(1:end-1));
    
    
    if strcmp(fun,'fit')==1
        
        a = Dposition./sqrt(2*Dtime);
        edges_step = diffusion_coefficient0/10;
        edges0 = -max(abs(a)):edges_step:max(abs(a));
        if length(edges0<50)
            edges0 = linspace(-max(abs(a)),max(abs(a)),50);
        end
        edges=(edges0(2:end)+edges0(1:end-1))/2;
        [N, edges0]=histcounts(a,edges0);
        
        [fitresult1, gof1] = fit(edges', N', 'a * exp(-0.5*(x)^2/b) + c * exp(-0.5*(x)^2/d)',...
             'Upper',[Inf,Inf, Inf, 3*abs(diffusion_coefficient_correction)],...
             'Lower',[0,abs(diffusion_coefficient_correction), 0, 0],...                                                                                                                                                                                                                                                                                                                                                                                                                              ,0,diffusion_coefficient_correction],...
             'StartPoint',[max(a),diffusion_coefficient0, 0, diffusion_coefficient_correction],...
             'TolFun',1e-9);
%              'Robust','LAR',...
        [fitresult2, gof2] = fit(edges', N', 'a * exp(-0.5*(x)^2/b)',...
             'Upper',[Inf,Inf],...
             'Lower',[0,0],...                                                                                                                                                                                                                                                                                                                                                                                                                              ,0,diffusion_coefficient_correction],...
             'StartPoint',[max(a),diffusion_coefficient0],...
             'TolFun',1e-9);
         
%          fitresult1.b
%          fitresult1.d
%          fitresult2.b
         
% %          if gof1.rsquare < gof2.rsquare
%            
%              if fitresult1.b>fitresult1.d
%                 diffusion_coefficient1=fitresult1.b;
%              else diffusion_coefficient1=fitresult1.d;
%              end
%              %disp('interaction with surface')
%              
%          %else 
%          diffusion_coefficient2=fitresult2.b;
%              disp('only diffusion')
%          %end
         
%          figure
%          plot(edges,N,'.'); hold on
%          plot(edges, fitresult1.a * exp(-0.5*(edges).^2/fitresult1.b));
%          plot(edges, fitresult1.c * exp(-0.5*(edges).^2/fitresult1.d));
%          plot(edges, fitresult2.a * exp(-0.5*(edges).^2/fitresult2.b));

            diffusion_coefficient0=[max([fitresult1.b, fitresult1.d]), fitresult2.b];
         
    end
    
    diffusion_coefficient = diffusion_coefficient0  + diffusion_coefficient_correction;
    
    epsilon = diffusion_coefficient_correction./diffusion_coefficient;
    diffusion_coefficient_variance = diffusion_coefficient.^2.*( (6 + 4*epsilon + 2*epsilon.^2)/length(Dposition) + 4*(1+epsilon).^2/(length(Dposition)).^2);
    diffusion_coefficient_error = sqrt(diffusion_coefficient_variance); 
    
else
    
    diffusion_coefficient=NaN;
    diffusion_coefficient_correction=NaN;
    diffusion_coefficient_error=NaN;
    
end
    


