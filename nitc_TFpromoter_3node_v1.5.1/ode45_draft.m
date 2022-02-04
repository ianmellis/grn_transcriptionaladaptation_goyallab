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

%
t_wtwt = t;
y_wtwt = y;
%% diploid wt/mut
% vector: y: 
% y(1) = ona_1 (wt)
% y(2) = offa_1 (wt)
% y(3) = ona_2 (mut)
% y(4) = offa_2 (mut)
% y(5) = a 
% y(6) = anons
% y(7) = onap_1 
% y(8) = offap_1 
% y(9) = onap_2 
% y(10) = offap_2 
% y(11) = aprim 
% y(12) = bba_1
% y(13) = bbap_1
% y(14) = unb_1
% y(15) = bba_2
% y(16) = bbap_2
% y(17) = unb_2
% y(18) = onb_1
% y(19) = offb_1
% y(20) = onb_2
% y(21) = offb_2
% y(22) = b
ts_wtmut = [0 1000];
ic_wtmut = [0;1;0;1;0;0;...
    0;1;0;1;0;...
    0;0;1;0;0;1;...
    0;1;0;1;0];

[t,y] = ode45(@(t,y) odefun_wtmut(t,y,ra_1), ts_wt, ic_wt);

t_wtmut = t;
y_wtmut = y;

%% het for all 100 paramsets

for i = 1:100
    
    rng(8574);
    
    % rate ratios
    % r_onbasal_A1/r_nitc = relative on-rates of wt A1 and NITC-regulated
    % alleles
    basal_nitc_on_ratio = latinhyp(i, 1);
    
    % r_onbasal_A1/r_off = on-off ratio, for speed of burst-off. Assume same
    % off rate for all other alleles
    onbasalA1_off_ratio = latinhyp(i, 2);
    
    % r_bound_byA1_B1/r_bound_byAprime1_B1 = relative B1 on-rates caused by bound
    % A1 vs Aprime1
    A1_Aprime1_boundOn_ratio = latinhyp(i, 3);
    
    % r_bind_byA1_B1/r_bind_byAprime1_B1 = relative promoter-binding rates of
    % A1 vs Aprime1
    A1_Aprime_binding_ratio = latinhyp(i, 4);
    
    % r_bind_byA1_B1/r_unbind_byA1_B1 = relative bind-unbind rates. Assume same
    % unbind rate for Aprime regardless of Aprime bind rate
    bindbyA1_unbindbyA1_ratio = latinhyp(i, 5);
    
    paramsetnum = i;
    
    % rates
    r_prod_basal = 0;
    r_prod_on = latinhyp(i, 6);
    r_deg = 0.1;
    r_onbasal_A1 = 1;
    r_onbasal_other = 0;
    r_nitc_byAnonsense1_A1 = r_onbasal_A1/basal_nitc_on_ratio;
    r_nitc_byAnonsense1_Anonsense1 = r_nitc_byAnonsense1_A1;
    r_nitc_byAnonsense1_Aprime1 = r_nitc_byAnonsense1_A1;
    r_bound_byA1_B1 = latinhyp(i, 7);
    r_bound_byAprime1_B1 = r_bound_byA1_B1/A1_Aprime1_boundOn_ratio;
    r_off = 1;
    n_all = latinhyp(i, 9);
    x = 0.5;
    r_bind_byA1_B1 = latinhyp(i, 8);
    r_bind_byAprime1_B1 = r_bind_byA1_B1/A1_Aprime_binding_ratio;
    r_unbind_byA1_B1 = r_bind_byA1_B1/bindbyA1_unbindbyA1_ratio;
    r_unbind_byAprime1_B1 = r_unbind_byA1_B1;
    
    k_all = x*(r_prod_basal + r_prod_on)/r_deg;
    
    r_prodbasal_A1 = r_prod_basal;
    r_prodbasal_Anonsense1 = r_prod_basal;
    r_prodbasal_Aprime1 = r_prod_basal;
    r_prodbasal_B1 = r_prod_basal;
    r_prodon_A1 = r_prod_on;
    r_prodon_Anonsense1 = r_prod_on;
    r_prodon_Aprime1 = r_prod_on;
    r_prodon_B1 = r_prod_on;
    r_deg_A1 = r_deg;
    r_deg_Anonsense1 = r_deg;
    r_deg_Aprime1 = r_deg;
    r_deg_B1 = r_deg;
    r_onbasal_A1 = r_onbasal_A1;
    r_onbasal_Anonsense1 = r_onbasal_A1;
    r_onbasal_Aprime1 = r_onbasal_other;
    r_onbasal_B1 = r_onbasal_other;
    r_nitc_byAnonsense1_A1 = r_nitc_byAnonsense1_A1;
    r_nitc_byAnonsense1_Anonsense1 = r_nitc_byAnonsense1_Anonsense1;
    r_nitc_byAnonsense1_Aprime1 = r_nitc_byAnonsense1_Aprime1;
    r_bound_byA1_B1 = r_bound_byA1_B1;
    r_bound_byAprime1_B1 = r_bound_byAprime1_B1;
    r_off_A1 = r_off;
    r_off_Anonsense1 = r_off;
    r_off_Aprime1 = r_off;
    r_off_B1 = r_off;
    r_bind_byA1_B1 = r_bind_byA1_B1;
    r_bind_byAprime1_B1 = r_bind_byAprime1_B1;
    k_A1 = k_all;
    k_Anonsense1 = k_all;
    k_Aprime1 = k_all;
    k_B1 = k_all;
    n_A1 = n_all;
    n_Anonsense1 = n_all;
    n_Aprime1 = n_all;
    n_B1 = n_all;
    r_unbind_byA1_B1 = r_unbind_byA1_B1;
    r_unbind_byAprime1_B1 = r_unbind_byAprime1_B1;
    
    ra_1 = [r_prodbasal_A1,...
        r_prodbasal_Anonsense1,...
        r_prodbasal_Aprime1,...
        r_prodbasal_B1,...
        r_prodon_A1,...
        r_prodon_Anonsense1,...
        r_prodon_Aprime1,...
        r_prodon_B1,...
        r_deg_A1,...
        r_deg_Anonsense1,...
        r_deg_Aprime1,...
        r_deg_B1,...
        r_onbasal_A1,...
        r_onbasal_Anonsense1,...
        r_onbasal_Aprime1,...
        r_onbasal_B1,...
        r_nitc_byAnonsense1_A1,...
        r_nitc_byAnonsense1_Anonsense1,...
        r_nitc_byAnonsense1_Aprime1,...
        r_bound_byA1_B1,...
        r_bound_byAprime1_B1,...
        r_off_A1,...
        r_off_Anonsense1,...
        r_off_Aprime1,...
        r_off_B1,...
        r_bind_byA1_B1,...
        r_bind_byAprime1_B1,...
        k_A1,...
        k_Anonsense1,...
        k_Aprime1,...
        k_B1,...
        n_A1,...
        n_Anonsense1,...
        n_Aprime1,...
        n_B1,...
        r_unbind_byA1_B1,...
        r_unbind_byAprime1_B1];
    
    ts_wtmut = [0 500];
    ic_wtmut = [0;1;0;1;0;0;...
        0;1;0;1;0;...
        0;0;1;0;0;1;...
        0;1;0;1;0];
    
    [t,y] = ode45(@(t,y) odefun_wtmut(t,y,ra_1), ts_wt, ic_wt);
    
    t_wtmut = t;
    y_wtmut = y;
    
end

