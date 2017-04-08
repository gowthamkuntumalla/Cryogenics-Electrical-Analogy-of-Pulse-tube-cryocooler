%Cryogenics all objects combined
clc;
clear;
%parameters%
f=44;
rho=9.2769;%density kg/m3
phi=.621; %porosity
A=(pi/4)*((15.5)^2)*1e-6;%Cross-section area
S=0.0104;%perimeter of intertance tube cross-section*
mu=13.637e-6;%dynamic viscosity*
nu=1.47e-6;%kinematic viscosity*
delv=sqrt(nu/(pi*f));%viscous penetration depth* 
d=mu*S/(A^2*delv);%some constant*
rg=6.87e9; %resistance in regenerator
gam=1.6615;%gamma
pm=33e5;% mean pressure
g=16.73;% Check this value
c=2*pi*f*A/gam/pm;%some constant
a=2*pi*f*rho/(phi*A); %a,b = some combined constants for simplicity
b=2*pi*f*phi*A/(gam*pm);
ro=2.06*10^9;% orifice resistance
%%%%%% Regenerator %%%%%%%
%syms Ur(x) Pr(x)
%Pr(x)= pm-int((a*1i+rg)*Ur(x),x,0,0.079);
%Ur(x)= Ur(0)-int(b*1i*Pr(x)+g*Ur(x),x,0,0.079);
syms u1(x1) p1(x1)
ode1 = diff(p1)+(a*1i+rg)*u1 == 0; % simultaneous differential equations
ode2 = diff(u1)+b*1i*p1+g*u1 == 0; % simultaneous differential equations
odes = [ode1; ode2]; % make a vector
cond1=p1(0)==3.15e5*(cosd(12.5)+1i*sind(12.5)); 
cond2=u1(0)==0.00042*(cosd(24)+1i*sind(24));
conds=[cond1; cond2];%initial conditions
[psol(x1), usol(x1)] = dsolve(odes,conds); % main solver
x1 = (0:0.005:0.079);
flowvec=subs(usol);
flow = abs(flowvec);
phase = angle(flowvec);
xo = (0:0.005:0.079);
%pressure
p_subs=subs(psol);
p_abs=abs(p_subs);
p_angle=angle(p_subs)*57.1;


%%%%%% Pulse Tube %%%%%%%
%preso=0.65e5*(cosd(-64)+1i*sind(-64));
syms u2(x2)
p2=psol(.079);
u2(x2)=usol(0.079)-c*1i*p2*x2;% p2= pressure in pulse tube(constant along length), u2 = volume flow rate 
x2 = (0:0.005:0.1);
usol2=subs(u2);
psol2=subs(p2);
phase2 = angle(usol2)*57.1;


%%%%%% Orifice %%%%%%%
p3 =p2-ro*u2(0.1);

%%%%%% Inertance Tube %%%%%
syms p1i(xi) ui(xi)
ode1i = diff(p1i)+(a*phi*1i+d)*ui == 0; % simultaneous differential equations
ode2i = diff(ui)+c*1i*p1i == 0; % simultaneous differential equations
odesi = [ode1i; ode2i]; % make a vector
cond1i=p1i(0)==p2; 
cond2i=ui(0)==u2(0.1);
condsi=[cond1i; cond2i];%initial conditions
[psoli(xi), usoli(xi)] = dsolve(odesi,condsi); % main solver
xi = (0:0.4:4.4);
flowveci=subs(usoli);
flowi = abs(flowveci);
phasei = angle(flowveci)*57.1;
presveci=subs(psoli);
presi = abs(presveci);
presphasei = angle(presveci)*57.1;

%%%%%% Reservoir %%%%%%%
syms u4(x4)
u4=u2(.1)-c*1i*p3*x4;