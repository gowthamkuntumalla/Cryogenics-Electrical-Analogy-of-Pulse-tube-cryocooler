
%function valueofk = k;
%valueofk=1;
f=44;%frequency
rg=193.216;
A=(pi/4)*(15.5)^2*1e-4;
k=.0167;%constant= omega^2 * rho / (gamma*pm)
l=(3.134e-5)*rg*A;%constant= omega*porosity*rg/(gamma * pm)*A
g=21;%In regenerator % This value was fixed upon 21 after iterating by verifying with paper

%********Volume**********%
syms qin(x); % q = flowrate, x = distance variable (metre)
cond1 = qin(0)==0.00042*(cosd(24)+1i*sind(24));
cond2 = qin(.079)==.00012*(cosd(-6)+1i*sind(-6));
conds = [cond1; cond2]; % initial conditions
q=dsolve(diff(diff(qin,x),x)+(k-l*1i)*qin+g*diff(qin,x)==0,conds);
x = (0:0.005:0.079);
flowvec=subs(q);
flow = abs(flowvec);
phase = angle(flowvec);
%***plots***%
xo = (0:0.005:0.079);
%plot(xo,flow);
%plot(xo,phase*57.1);
grid on

%********Pressure**********%
syms p(x);% p= pressure
a=2.83e6;%omega*rhp/(porosity*area)
p_sol = dsolve(diff(p,x)+(rg+a*1i)*q==0,p(0)==33e5);
x = (0:0.005:0.079);
p_subs=subs(p_sol);
p_abs=abs(p_subs);
p_angle=angle(p_subs);
%plot(xo,p_abs);
%plot(xo,p_angle*57.1);
