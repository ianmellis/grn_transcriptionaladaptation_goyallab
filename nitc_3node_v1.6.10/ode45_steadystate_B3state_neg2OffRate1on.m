% ODE45-based steady state analysis of NITC_3node_B3state_v1.6.8 (inhib)
% simulations with regulation of target off-rate

% runs ode45 for finding steady-state expression values of all species in
% all 3 ancestral genotypes (wt/wt, wt/mut, mut/mut). Uses the same 100 LHS
% parameter sets from v1.6.8

% edit as needed
outdir = '~/code/grn_nitc/nitc_3node_v1.6.10/';


%% LHS
% Larsson et al., 2019 estimated burst kinetics transcriptome-wide using
% single-cell RNA-seq data. They inferred r_on, r_off, and r_prod relative
% to r_deg for individual alleles. 
% 
% Fixed parameter values: r_deg = 1, x= 0.5. Can do sweep over x later.
% 
% search over: basal_nitc_on_ratio, onbasalA1_off_ratio,
% A1_Aprime1_addon_ratio, A1_Aprime_prodon_ratio,
% r_prod_on, r_addon_byA1_B1, n (Hill coefficient n), r_onbasal_A1
%
% Ratio definitions:
% - basal_nitc_on_ratio = r_onbasal_A1/r_nitc = relative on-rates of wt A1 and NITC-regulated
% - onbasalA1_off_ratio = r_onbasal_A1/r_off = on-off ratio, for speed of burst-off
% - A1_Aprime1_addon_ratio = r_addon_byA1_B1/r_addon_byAprime1_B1 = scalar
% multiple strength decrease of Aprime1 relative to A1 in causing B1-on
% - A1_Aprime_prodon_ratio = d_Aprime_B1 = scalar multiple descrease in
% production rate of B1 burst when turned on by para instead of orig
%  - r_onbasal_Aprime_ratio = r_onbasal_A1/r_onbasal_Aprime1 = scalar
%
rng(8563);

nruns = 10000;

fseeds = 358239:4:3409333;
simseeds = fseeds(1:nruns);

min_range = [repmat(0.1,1,1),...        %basal_nitc_on_ratio (0.1,10) - is 10 high enough?
    repmat(0.01,1,1),...                %onbasalA1_off_ratio (0.01,2) - per Larsson this is the bulk of the distribution (very quick off-burst)
    repmat(0.1,1,1),...                 %A1_Aprime1_addon_ratio (0.1,10)
    repmat(0.1,1,1),...                 %A1_Aprime_prodon_ratio (0.1,10)
    repmat(1,1,1),...                   %r_prod_on (1,1000)
    repmat(0.1,1,1),...                 %r_addon_byA1_B1 (0.1,10)
    repmat(0.1,1,1),...                 %n (Hill coefficient n) (0.1,5) too large a range? I think so, but having some non-linearity is good. Perhaps ok with 10000 paramsets.
    repmat(0.1,1,1),...                 %r_onbasal_A1 (0.1, 10)
    repmat(1,1,1)];                     %r_onbasal_Aprime_ratio (1, 100)
%                                       

max_range = [repmat(10,1,1),...         %basal_nitc_on_ratio (0.1,10) - is 10 high enough?
    repmat(2,1,1),...                   %onbasalA1_off_ratio (0.01,2) - per Larsson this is the bulk of the distribution (very quick off-burst)
    repmat(10,1,1),...                  %A1_Aprime1_addon_ratio (0.1,10)
    repmat(10,1,1),...                  %A1_Aprime_prodon_ratio (0.1,10)
    repmat(1000,1,1),...                %r_prod_on (1,1000)
    repmat(10,1,1),...                  %r_addon_byA1_B1 (0.1,10)
    repmat(5,1,1),...                   %n (Hill coefficient n) (0.1,5) too large a range?
    repmat(10,1,1),...                  %r_onbasal_A1 (0.1, 10)
    repmat(100,1,1)];                   %r_onbasal_Aprime_ratio (1, 100)

latinhyp = lhsdesign_modified(nruns, min_range, max_range);

%%
% % rate ratios
% % r_onbasal_A1/r_nitc = relative on-rates of wt A1 and NITC-regulated
% % alleles
% basal_nitc_on_ratios = [10, 10, 10, 2, 10];
% 
% % r_onbasal_A1/r_off = on-off ratio, for speed of burst-off. Assume same
% % off rate for all other alleles
% onbasalA1_off_ratios = [1, 1, 1, 1, 1];
% 
% % r_bound_byA1_B1/r_bound_byAprime1_B1 = relative B1 on-rates caused by bound
% % A1 vs Aprime1
% A1_Aprime1_boundOn_ratios = [5, 10, 1, 5, 5];
% 
% % r_bind_byA1_B1/r_bind_byAprime1_B1 = relative promoter-binding rates of
% % A1 vs Aprime1
% A1_Aprime_binding_ratios = [5, 1, 10, 5, 5];
% 
% % r_bind_byA1_B1/r_unbind_byA1_B1 = relative bind-unbind rates. Assume same
% % unbind rate for Aprime regardless of Aprime bind rate
% bindbyA1_unbindbyA1_ratios = [100, 100, 100, 100, 10];
lhs_1_f = [psc_outdir, 'latinhyp_sampledSets.csv'];
lhs_1_f2 = [trace_outdir, 'latinhyp_sampledSets.csv'];

lhs_1_s = array2table(latinhyp);
lhs_1_s.Properties.VariableNames = {'basal_nitc_on_ratio',... 
    'onbasalA1_off_ratio',...
    'A1_Aprime1_addon_ratio',... 
    'A1_Aprime_prodon_ratio',...
    'r_prod_on',...
    'r_addon_byA1_B1',...
    'Hill_coefficient_n',...
    'r_onbasal_A1',...
    'r_onbasal_Aprime_ratio'};
%% het for all 100 paramsets
ssSP_ww = zeros(100,5);
ssSP_wm = zeros(100,5);
ssSP_mm = zeros(100,5);
for i = 1:100
    
    % Set seed
    %
    rng(8363);
    
    % rate ratios
    % r_onbasal_A1/r_nitc = relative on-rates of wt A1 and NITC-regulated
    % alleles
    basal_nitc_on_ratio = latinhyp(i, 1);
    
    % r_onbasal_A1/r_off = on-off ratio, for speed of burst-off. Assume same
    % off rate for all other alleles
    onbasalA1_off_ratio = latinhyp(i, 2);
    
    % r_bound_byA1_B1/r_bound_byAprime1_B1 = relative B1 on-rates caused by bound
    % A1 vs Aprime1
    A1_Aprime1_addon_ratio = latinhyp(i, 3);
    
    % r_bind_byA1_B1/r_bind_byAprime1_B1 = relative promoter-binding rates of
    % A1 vs Aprime1
    A1_Aprime_onrate_ratio = latinhyp(i, 4);
    
    paramsetnum = i;
    
    % rates
    r_prod_basal = 0;
    r_prod_on = latinhyp(i, 5);
    r_deg = 1;
    r_onbasal_A1 = latinhyp(i, 8);
    r_onbasal_Aprime = r_onbasal_A1/latinhyp(i, 9);
    r_onbasal_other = r_onbasal_A1; % basal B1 expression same as A1
    r_nitc_byAnonsense1_A1 = r_onbasal_A1/basal_nitc_on_ratio;
    r_nitc_byAnonsense1_Anonsense1 = r_nitc_byAnonsense1_A1;
    r_nitc_byAnonsense1_Aprime1 = r_nitc_byAnonsense1_A1;
    r_addon_byA1_B1 = latinhyp(i, 6);
    r_addon_byAprime1_B1 = r_addon_byA1_B1/A1_Aprime1_addon_ratio;
    r_off = r_onbasal_A1/onbasalA1_off_ratio;
    n_all = latinhyp(i, 7);
    x = 0.5;
    
    k_all = x*(r_prod_basal + r_prod_on)/r_deg;
    
    r_prodbasal_A1 = r_prod_basal;
    r_prodbasal_Anonsense1 = r_prod_basal;
    r_prodbasal_Aprime1 = r_prod_basal;
    r_prodbasal_B1 = r_prod_basal;
    r_prodon_A1 = r_prod_on;
    r_prodon_Anonsense1 = r_prod_on;
    r_prodon_Aprime1 = r_prod_on;
    r_prodon_B1 = r_prod_on;
    
    d_Aprime1_B1 = A1_Aprime_onrate_ratio; % multiplied in sim
    
    r_deg_A1 = r_deg;
    r_deg_Anonsense1 = r_deg;
    r_deg_Aprime1 = r_deg;
    r_deg_B1 = r_deg;
    r_onbasal_A1 = r_onbasal_A1;
    r_onbasal_Anonsense1 = r_onbasal_A1;
    r_onbasal_Aprime1 = r_onbasal_Aprime;
    r_onbasal_B1 = r_onbasal_other;
    
    r_nitc_byAnonsense1_A1 = r_nitc_byAnonsense1_A1;
    r_nitc_byAnonsense1_Anonsense1 = r_nitc_byAnonsense1_Anonsense1;
    r_nitc_byAnonsense1_Aprime1 = r_nitc_byAnonsense1_Aprime1;
    r_addon_byA1_B1 = r_addon_byA1_B1;
    r_addon_byAprime1_B1 = r_addon_byAprime1_B1;
    
    r_off_A1 = r_off;
    r_off_Anonsense1 = r_off;
    r_off_Aprime1 = r_off;
    r_offorig_B1 = r_off;
    r_offpara_B1 = 0; % basal off-para rate
    
    k_A1 = k_all;
    k_Anonsense1 = k_all;
    k_Aprime1 = k_all;
    k_B1 = k_all;
    n_A1 = n_all;
    n_Anonsense1 = n_all;
    n_Aprime1 = n_all;
    n_B1 = n_all;
    
    
    ra_1 = [r_prodbasal_A1,...
        r_prodbasal_Anonsense1,...
        r_prodbasal_Aprime1,...
        r_prodbasal_B1,...
        r_prodon_A1,...
        r_prodon_Anonsense1,...
        r_prodon_Aprime1,...
        r_prodon_B1,...
        d_Aprime1_B1,...
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
        r_addon_byA1_B1,...
        r_addon_byAprime1_B1,...
        r_off_A1,...
        r_off_Anonsense1,...
        r_off_Aprime1,...
        r_offorig_B1,...
        r_offpara_B1,...
        k_A1 ,...
        k_Anonsense1 ,...
        k_Aprime1 ,...
        k_B1 ,...
        n_A1 ,...
        n_Anonsense1 ,...
        n_Aprime1 ,...
        n_B1];
    
    
    ra_1_s = array2table(ra_1);
    ra_1_s.Properties.VariableNames = {'r_prodbasal_A1',...
        'r_prodbasal_Anonsense1',...
        'r_prodbasal_Aprime1',...
        'r_prodbasal_B1',...
        'r_prodon_A1',...
        'r_prodon_Anonsense1',...
        'r_prodon_Aprime1',...
        'r_prodon_B1',...
        'd_Aprime1_B1',...
        'r_deg_A1',...
        'r_deg_Anonsense1',...
        'r_deg_Aprime1',...
        'r_deg_B1',...
        'r_onbasal_A1',...
        'r_onbasal_Anonsense1',...
        'r_onbasal_Aprime1',...
        'r_onbasal_B1',...
        'r_nitc_byAnonsense1_A1',...
        'r_nitc_byAnonsense1_Anonsense1',...
        'r_nitc_byAnonsense1_Aprime1',...
        'r_addon_byA1_B1',...
        'r_addon_byAprime1_B1',...
        'r_off_A1',...
        'r_off_Anonsense1',...
        'r_off_Aprime1',...
        'r_offorig_B1',...
        'r_offpara_B1',...
        'k_A1',...
        'k_Anonsense1',...
        'k_Aprime1',...
        'k_B1',...
        'n_A1',...
        'n_Anonsense1',...
        'n_Aprime1',...
        'n_B1'};
    
    ts = [0 500];
    ic_wtwt = [0;1;0;1;0;0;...
        0;1;0;1;0;...
        1;0;0;1;0;0;0];
    
    [t,y] = ode45(@(t,y) odefun_wtwt_B3state_neg2OffRate1on(t,y,ra_1), ts, ic_wtwt); 
    
    ss_A1 = real(y(size(y,1),5));
    ss_Anons1 = real(y(size(y,1),6));
    ss_Aprim1 = real(y(size(y,1),11));
    ss_B1 = real(y(size(y,1),18));
    
    ssSP_ww(i,:) = [i,ss_A1,ss_Anons1,ss_Aprim1,ss_B1];
    
    ic_wtmut = ic_wtwt;
    
    [t,y] = ode45(@(t,y) odefun_wtmut_B3state_neg2OffRate1on(t,y,ra_1), ts, ic_wtmut);
    
    ss_A1 = real(y(size(y,1),5));
    ss_Anons1 = real(y(size(y,1),6));
    ss_Aprim1 = real(y(size(y,1),11));
    ss_B1 = real(y(size(y,1),18));
    
    ssSP_wm(i,:) = [i,ss_A1,ss_Anons1,ss_Aprim1,ss_B1];
    
    ic_mutmut = ic_wtwt;
    
    [t,y] = ode45(@(t,y) odefun_mutmut_B3state_neg2OffRate1on(t,y,ra_1), ts, ic_mutmut);
    
    ss_A1 = real(y(size(y,1),5));
    ss_Anons1 = real(y(size(y,1),6));
    ss_Aprim1 = real(y(size(y,1),11));
    ss_B1 = real(y(size(y,1),18));
    
    ssSP_mm(i,:) = [i,ss_A1,ss_Anons1,ss_Aprim1,ss_B1];
    

    
end

%%
ss_file = [outdir, 'steady_state_ODE45_wtwt_B3state_neg2OffRate1onODE.csv'];
ss_table = array2table(ssSP_ww, 'VariableNames', {'paramset', 'ss_A1','ss_Anons1','ss_Aprim1','ss_B1'});
writetable(ss_table, ss_file, 'Delimiter', ',')

ss_file = [outdir, 'steady_state_ODE45_wtmut_B3state_neg2OffRate1onODE.csv'];
ss_table = array2table(ssSP_wm, 'VariableNames', {'paramset', 'ss_A1','ss_Anons1','ss_Aprim1','ss_B1'});
writetable(ss_table, ss_file, 'Delimiter', ',')

ss_file = [outdir, 'steady_state_ODE45_mutmut_B3state_neg2OffRate1onODE.csv'];
ss_table = array2table(ssSP_mm, 'VariableNames', {'paramset', 'ss_A1','ss_Anons1','ss_Aprim1','ss_B1'});
writetable(ss_table, ss_file, 'Delimiter', ',')

