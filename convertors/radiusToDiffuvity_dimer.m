
clear;

L = 12;%length
d = 6;%diamater
r = L/d;

kB=1.38064852*1e-23;	%[J?K?1]= [ kg?m2?s?2?K?1]
T=273.15 + 25 ;%absolute temperature [K]
ny=  8.9*1e-4 ;%viscosity [Pa?s] [(N?s)/m2 = kg/(s?m)]

beta = acosh(r)/r/sqrt(r^2-1);
alpha_parallel=2*(r^2*beta-1)/(r^2-1);
alpha_perpedicular=r^2*(1-beta)/(r^2-1);
f_parallel=8*pi*ny*L/r^2/(2*beta+alpha_parallel);
f_perpedicular=8*pi*ny*L/(2*r^2*beta+alpha_perpedicular);

D_parallel=kB*T/f_parallel;
D_perpendicular=kB*T/f_perpedicular;

D_ellipsoid=1/3*(2*D_perpendicular+D_parallel)*1e18

rp=d/2;
D_sphere=kB*T./(6*pi*ny*rp*1e-6) *1e12

1/(D_ellipsoid/D_sphere)

