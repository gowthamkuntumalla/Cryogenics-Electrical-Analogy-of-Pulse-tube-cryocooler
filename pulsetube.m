%CRYOGENICS SEMINAR 
%PULSE TUBE ELECTRICAL ANALOGY
% Figure 15
%%%%%%constants%%%%%%
clc;
clear;
f=44;
rho=9.2769;%density kg/m3
phi=.621; %porosity
A=(pi/4)*((15.5)^2)*1e-6;%Cross-section area
rg=10e5; %resistance in regenerator
gam=1.6615;%gamma
pm=33e5;% mean pressure
g=12.6;% Check this value
a=2*pi*f*A/gam/pm;%some constant
%%%%%%%solve linear equation%%%%%%%%%
%syms u(x) 
flowo=1.16e-4*(cosd(-34.5)+1i*sind(-34.5));% This can be obtained from regenerator flow rate at x=.079
x = (0:0.005:0.1);
preso=.65e5*(cosd(-64)+1i*sind(-64));
u=flowo-a*1i*preso*x;
%flowvec=subs(u);
flow = abs(u);
phase = angle(u);
xo = (0:0.005:0.1);
plot(xo,flow)
%plot(xo,phase*57.1)
