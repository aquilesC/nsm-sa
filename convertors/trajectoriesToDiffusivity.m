function [diffusion_coefficient, diffusion_coefficient_correction, diffusion_coefficient_error] = trajectoriesToDiffusivity(position, time, frame,tav)


% Estimation of the diffusion coefficient based on Vestergaard, C. L.;  Blainey, P. C.; Flyvbjerg, H., Optimal estimation of diffusion coefficients from single-particle trajectories. Phys Rev E 2014, 89 (2).

% position - a time serie of measured positions of a particle
% time - corresponding time of the measurement
% frame - no. of the frame
% tav - time smoothen factor

% diffusion_coefficient - estimated diffusion coefficient corrected for the motion blur and localization error [um^2/s]
% diffusion_coefficient_correction - correction term [um^2/s]
% diffusion_coefficient_error - estimated error of the value [um^2/s]

if length(position)>2

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
    diffusion_coefficient = diffusion_coefficient0  + diffusion_coefficient_correction;
    
    epsilon = diffusion_coefficient_correction./diffusion_coefficient;
    diffusion_coefficient_variance = diffusion_coefficient.^2.*( (6 + 4*epsilon + 2*epsilon.^2)/length(Dposition) + 4*(1+epsilon).^2/(length(Dposition)).^2);
    diffusion_coefficient_error = sqrt(diffusion_coefficient_variance); 
    
else
    
    diffusion_coefficient=NaN;
    diffusion_coefficient_correction=NaN;
    diffusion_coefficient_error=NaN;
    
end
    


