%% dsolve attempt

syms ona(t) offa(t) a(t) bba(t) unb(t) onb(t) offb(t) b(t)

ode1 = diff(ona) == r_onbasal_A1*offa - r_off_A1*ona;
ode2 = diff(offa) == -r_onbasal_A1*offa + r_off_A1*ona;
ode3 = diff(a) == r_prodon_A1*ona - r_deg_A1*a;
ode4 = diff(bba) == (r_bind_byA1_B1*a^n_A1/(k_A1^n_A1 + a^n_A1))*unb - r_unbind_byA1_B1*bba;
ode5 = diff(unb) == -(r_bind_byA1_B1*a^n_A1/(k_A1^n_A1 + a^n_A1))*unb + r_unbind_byA1_B1*bba;
ode6 = diff(onb) == r_bound_byA1_B1*offb - r_off_B1*onb;
ode7 = diff(offb) == -r_bound_byA1_B1*offb + r_off_B1*onb;
ode8 = diff(b) == r_prodon_B1*onb - r_deg_B1*b;

%%
odes = [ode1;ode2;ode3;ode4;ode5;ode6;ode7;ode8];

%%
S = dsolve(odes);
% explicit solution could not be found

%% ODE45 version

% vector: y: 
% y(1) = ona
% y(2) = offa
% y(3) = a 
% y(4) = bba
% y(5) = unb
% y(6) = onb
% y(7) = offb 
% y(8) = b
% function dydt = odefun1(t,y)
% dydt = zeros(8,1);
% dydt(1) = r_onbasal_A1*y(2) - r_off_A1*y(1);
% dydt(2) = -r_onbasal_A1*y(2) + r_off_A1*y(1);
% dydt(3) = r_prodon_A1*y(1) - r_deg_A1*y(3);
% dydt(4) = (r_bind_byA1_B1*y(3)^n_A1/(k_A1^n_A1 + y(3)^n_A1))*y(5) - r_unbind_byA1_B1*y(4);
% dydt(5) = -(r_bind_byA1_B1*y(3)^n_A1/(k_A1^n_A1 + y(3)^n_A1))*y(5) + r_unbind_byA1_B1*y(4);
% dydt(6) = r_bound_byA1_B1*y(7) - r_off_B1*y(6);
% dydt(7) = -r_bound_byA1_B1*y(7) + r_off_B1*y(6);
% dydt(8) = r_prodon_B1*y(6) - r_deg_B1*y(8);
% end
%%
ts = [0 100];
ic = [0;1;0;0;1;0;1;0];
%%
[t,y] = ode45(@(t,y) odefun1(t,y,ra_1), ts, ic);

t1=t; y1=y;

%% diploid wt/wt
ts_wt = [0 100];
ic_wt = [0;1;0;1;0;0;1;0;1;0;1;0;1;0];

[t,y] = ode45(@(t,y) odefun_wtwt(t,y,ra_1), ts_wt, ic_wt);

%% diploid wt/mut
ts_wt = [0 100];
ic_wt = [0;1;0;1;0;0;1;0;1;0;1;0;1;0];

[t,y] = ode45(@(t,y) odefun_wtmut(t,y,ra_1), ts_wt, ic_wt);