function nanochannelArea = list_of_nanochannels (channel)

%in micrometers

switch channel
    case 'channel_30x85' %Channel I - 200505A
        nanochannelArea=2.7e-3;
    case 'channel_50x70' %Channel II - 190314A
        nanochannelArea=7.9e-3;
    case 'channel_20x85' %Channel IIi - 200505A
        nanochannelArea=1.48e-3;
    case 'channel_30x135' %Channel IV - 200505A
        nanochannelArea=3.92e-3; 
    case 'channel_200x200' %2120121
        nanochannelArea=0.2024*0.2213;
        
    case 'channel_30x85_AZ' %Channel I - 200505A
        nanochannelArea=2.7e-3/0.8542;    
    
end