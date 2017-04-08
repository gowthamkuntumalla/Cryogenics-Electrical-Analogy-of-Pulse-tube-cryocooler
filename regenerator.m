%CRYOGENICS SEMINAR 
%RGENERATOR ELECTRICAL ANALOGY
%http://in.mathworks.com/help/symbolic/solve-a-system-of-differential-equations.html
% Asish ' Hostel 13 C wing 104
%%%%%%constants%%%%%%
f=44;
rho=9.2769;%density kg/m3
phi=.621; %porosity
A=(pi/4)*(15.5)^2*1e-6;%Cross-section area
rg=10e5; %resistance in regenerator
gam=1.6615;%gamma
pm=33e5;% mean pressure
g=12.6;% Check this value
a=2*pi*f*rho/(phi*A); %a,b = some combined constants for simplicity
b=2*pi*f*phi/(gam*pm);

%%%%%%%solve ode%%%%%%%%%
syms u(x) p(x)
ode1 = diff(p)+(a*1i+rg)*u == 0; % simultaneous differential equations
ode2 = diff(u)+b*1i*A*p+g*u == 0; % simultaneous differential equations
odes = [ode1; ode2]; % make a vector
cond1=p(0)==3.3e5*(cosd(12.5)+1i*sind(12.5)); 
cond2=u(0)==0.00042*(cosd(24)+1i*sind(24));
conds=[cond1; cond2];%initial conditions
[psol(x), usol(x)] = dsolve(odes,conds); % main solver
x = (0:0.005:0.079);

%volume flowrate
flowvec=subs(usol);
flow = abs(flowvec);
phase = angle(flowvec);
%pressure
p_subs=subs(psol);
p_abs=abs(p_subs);
p_angle=angle(p_subs);
%%%%%%%plots%%%%%%
xo = (0:0.005:0.079);
%plot(xo,flow);
%plot(xo,phase*57.1);

%plot(xo,p_abs);
%plot(xo,p_angle*57.1);