function nanochannelArea = list_of_nanochannels (channel)

%in micrometers

switch channel
    case 'channel_30x85' %Channel I - 200505A
        nanochannelArea=2.7e-3;
    case 'channel_50x70' %Channel II - 190314A
        nanochannelArea=7.9e-3;
    case 'channel_20x85' %Channel III - 200505A
        nanochannelArea=1.48e-3;
    case 'channel_30x135' %Channel IV - 200505A
        nanochannelArea=3.92e-3; 
    case 'channel_200x200' %Channel V - 2120121
        nanochannelArea=0.225*0.2;
    case 'channel_200x200_coated' %Channel V - coated by LUV - 2120121 
        nanochannelArea=[0.225*0.2,(0.225-4.9e-3*2)*(0.2-4.9e-3*2),1.12]; %[the whole area, area without the shell, difference between the intensity before and after the deposition]
        
    case 'channel_B' %Channel VI - 210618
        nanochannelArea=0.08*0.04; 
        
    case 'channel_B_coated' %Channel V - coated by LUV - 210618
        nanochannelArea=[0.08*0.04,0.072*0.032,1.7]; %[the whole area, area without the shell, difference between the intensity before and after the deposition]
            
    
end
