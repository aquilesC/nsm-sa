function rp = DiffusivityToRadius(D)


%rp in micrometers

%from /Users/spackova/Box Sync/_project/transport_simulation/diffusivity
kB=1.38064852*1e-23;	%[J?K?1]= [ kg?m2?s?2?K?1]
T=273.15 + 25 ;%absolute temperature [K]
ny=  8.9*1e-4 ;%viscosity [Pa?s] [(N?s)/m2 = kg/(s?m)]


%D=kB*T./(6*pi*ny*rp*1e-6) *1e12; %diffusivity [m2/s * 1e-12 = um2/s]
rp = kB*T./(6*pi*ny*D*1e-6) *1e12;