% ODE45-based steady state analysis of NITC-TFpromoter_3node_v1.5.1
% simulations

% runs ode45 for finding steady-state expression values of all species in
% all 3 ancestral genotypes (wt/wt, wt/mut, mut/mut). Uses the same 100 LHS
% parameter sets from v1.5.1

% edit as needed
outdir = '~/code/grn_nitc/nitc_TFpromoter_3node_v1.5.1/';

%% LHS
% Fixed parameter values: r_deg = 0.1, r_onbasal_A1 = 1, x= 0.5
% 
% search over: basal_nitc_on_ratio, onbasalA1_off_ratio,
% A1_Aprime1_boundOn_ratio, A1_Aprime_binding_ratio,
% bindbyA1_unbindbyA1_ratio, r_prod_on, r_bound_byA1_B1, r_bind_byA1_B1, 
% r_unbind_byA1_B1
rng(8734);

nruns = 100;

min_range = [repmat(0.5,1,1),...       %basal_nitc_on_ratio
    repmat(0.1,1,1),...                 %onbasalA1_off_ratio
    repmat(1,1,1),...                   %A1_Aprime1_boundOn_ratio
    repmat(0.1,1,1),...                 %A1_Aprime_binding_ratio
    repmat(1,1,1),...                   %bindbyA1_unbindbyA1_ratio
    repmat(1,1,1),...                   %r_prod_on
    repmat(0.1,1,1),...                 %r_bound_byA1_B1
    repmat(0.5,1,1),...                 %r_bind_byA1_B1
    repmat(0.1,1,1)];                   %n (Hill coefficient n)

max_range = [repmat(25,1,1),...          %basal_nitc_on_ratio
    repmat(10,1,1),...                  %onbasalA1_off_ratio
    repmat(25,1,1),...                  %A1_Aprime1_boundOn_ratio
    repmat(10,1,1),...                  %A1_Aprime_binding_ratio
    repmat(100,1,1),...                 %bindbyA1_unbindbyA1_ratio
    repmat(100,1,1),...                 %r_prod_on
    repmat(5,1,1),...                   %r_bound_byA1_B1
    repmat(50,1,1),...                  %r_bind_byA1_B1
    repmat(10,1,1)];                    %n (Hill coefficient n)

latinhyp = lhsdesign_modified(nruns, min_range, max_range);

% %%
% latinhyp_prod = repmat(latinhyp(:,1),1,n_species);
% latinhyp_deg = repmat(latinhyp(:,2),1,n_species);
% latinhyp_onbasal = repmat(latinhyp(:,3),1,n_species_upstr);
% latinhyp_n = repmat(latinhyp(:,4),1,n_species);
% latinhyp_ondep = repmat(latinhyp(:,5),1,n_species_downstr);
% latinhyp_off = repmat(latinhyp(:,6),1,n_species);
% latinhyp_proddiff = repmat(latinhyp(:,7),1,n_species);
% 
% %     new parameters
% latinhyp_onbasal_multiple_aprime = repmat(latinhyp(:,8),1,n_species_paralog); % change this to only apply to paralogs
% latinhyp_onbasal_multiple_b = repmat(latinhyp(:,9),1,n_species_downstr);
% latinhyp_ondep_multiple_prime = repmat(latinhyp(:,10),1,n_species_downstr);
% latinhyp_nitc = repmat(latinhyp(:,11),1,n_species_paralog);
% 
% x = 0.5;
% 
% %% wrap for multiple paramsets
% 
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
% lhs_1_f = [outdir, 'latinhyp_sampledSets.csv'];

lhs_1_s = array2table(latinhyp);
lhs_1_s.Properties.VariableNames = {'basal_nitc_on_ratio',...
    'onbasalA1_off_ratio',...
    'A1_Aprime1_boundOn_ratio',...
    'A1_Aprime_binding_ratio',...
    'bindbyA1_unbindbyA1_ratio',...
    'r_prod_on',...
    'r_bound_byA1_B1',...
    'r_bind_byA1_B1',...
    'Hill_coefficient_n'};
%% LHS
% Larsson et al., 2019 estimated burst kinetics transcriptome-wide using
% single-cell RNA-seq data. They inferred r_on, r_off, and r_prod relative
% to r_deg for individual alleles. 
% 
% Fixed parameter values: r_deg = 1, x= 0.5. 
% 
% search over: basal_nitc_on_ratio, onbasalA1_off_ratio,
% A1_Aprime1_addon_ratio, A1_Aprime_prodon_ratio,
% r_prod_on, r_addon_byA1_B1, n (Hill coefficient n), r_onbasal_A1
rng(8734);

nruns = 100;

min_range = [repmat(0.1,1,1),...        %basal_nitc_on_ratio (0.1,10) - is 10 high enough?
    repmat(0.01,1,1),...                %onbasalA1_off_ratio (0.01,1) - per Larsson this is the bulk of the distribution (very quick off-burst)
    repmat(0.1,1,1),...                 %A1_Aprime1_addon_ratio (0.1,10)
    repmat(0.1,1,1),...                 %A1_Aprime_prodon_ratio (0.1,10)
    repmat(1,1,1),...                   %r_prod_on (1,1000)
    repmat(0.1,1,1),...                 %r_addon_byA1_B1 (0.1,10)
    repmat(0.1,1,1),...                 %n (Hill coefficient n) (0.1,10) too large a range?
    repmat(0.1,1,1)];                   %r_onbasal_A1 (0.1, 10)
%                                       

max_range = [repmat(10,1,1),...         %basal_nitc_on_ratio (0.1,10) - is 10 high enough?
    repmat(1,1,1),...                   %onbasalA1_off_ratio (0.01,1) - per Larsson this is the bulk of the distribution (very quick off-burst)
    repmat(10,1,1),...                  %A1_Aprime1_addon_ratio (0.1,10)
    repmat(10,1,1),...                  %A1_Aprime_prodon_ratio (0.1,10)
    repmat(1000,1,1),...                %r_prod_on (1,1000)
    repmat(10,1,1),...                  %r_addon_byA1_B1 (0.1,10)
    repmat(10,1,1),...                  %n (Hill coefficient n) (0.1,10) too large a range?
    repmat(10,1,1)];                    %r_onbasal_A1 (0.1, 10)

latinhyp = lhsdesign_modified(nruns, min_range, max_range);

% %%
% latinhyp_prod = repmat(latinhyp(:,1),1,n_species);
% latinhyp_deg = repmat(latinhyp(:,2),1,n_species);
% latinhyp_onbasal = repmat(latinhyp(:,3),1,n_species_upstr);
% latinhyp_n = repmat(latinhyp(:,4),1,n_species);
% latinhyp_ondep = repmat(latinhyp(:,5),1,n_species_downstr);
% latinhyp_off = repmat(latinhyp(:,6),1,n_species);
% latinhyp_proddiff = repmat(latinhyp(:,7),1,n_species);
% 
% %     new parameters
% latinhyp_onbasal_multiple_aprime = repmat(latinhyp(:,8),1,n_species_paralog); % change this to only apply to paralogs
% latinhyp_onbasal_multiple_b = repmat(latinhyp(:,9),1,n_species_downstr);
% latinhyp_ondep_multiple_prime = repmat(latinhyp(:,10),1,n_species_downstr);
% latinhyp_nitc = repmat(latinhyp(:,11),1,n_species_paralog);
% 
% x = 0.5;
% 
% %% wrap for multiple paramsets
% 
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
lhs_1_f = [outdir, 'latinhyp_sampledSets.csv'];

lhs_1_s = array2table(latinhyp);
lhs_1_s.Properties.VariableNames = {'basal_nitc_on_ratio',... 
    'onbasalA1_off_ratio',...
    'A1_Aprime1_addon_ratio',... 
    'A1_Aprime_prodon_ratio',...
    'r_prod_on',...
    'r_addon_byA1_B1',...
    'Hill_coefficient_n',...
    'r_onbasal_A1'};


%% het for all 100 paramsets
ssSP_ww = zeros(100,5);
ssSP_wm = zeros(100,5);
ssSP_mm = zeros(100,5);
for i = 1:100
    
    % Set seed
    %
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
    A1_Aprime1_addon_ratio = latinhyp(i, 3);
    
    % r_bind_byA1_B1/r_bind_byAprime1_B1 = relative promoter-binding rates of
    % A1 vs Aprime1
    A1_Aprime1_prodon_ratio = latinhyp(i, 4);
    
    paramsetnum = i;
    
    % rates
    r_prod_basal = 0;
    r_prod_on = latinhyp(i, 5);
    r_deg = 1;
    r_onbasal_A1 = latinhyp(i, 8);
    r_onbasal_other = 0;
    r_nitc_byAnonsense1_A1 = r_onbasal_A1/basal_nitc_on_ratio;
    r_nitc_byAnonsense1_Anonsense1 = r_nitc_byAnonsense1_A1;
    r_nitc_byAnonsense1_Aprime1 = r_nitc_byAnonsense1_A1;
    r_addon_byA1_B1 = latinhyp(i, 6);
    r_addon_byAprime1_B1 = r_bound_byA1_B1/A1_Aprime1_addon_ratio;
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
    
    d_Aprime1_B1 = 1/A1_Aprime1_prodon_ratio; % multiplied in sim, so reciprocal taken here
    
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
    r_addon_byA1_B1 = r_addon_byA1_B1;
    r_addon_byAprime1_B1 = r_addon_byAprime1_B1/A1_Aprime1_addon_ratio;
    
    r_off_A1 = r_off;
    r_off_Anonsense1 = r_off;
    r_off_Aprime1 = r_off;
    r_offorig_B1 = r_off;
    r_offpara_B1 = r_off;
    
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
        r_prodbasal_Aprime1,...2
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
        'r_prodbasal_Aprime1',...2
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
    
    % Initialize species
    
    A1 = 0;
    Anonsense1 = 0;
    Aprime1 = 0;
    B1 = 0;
    
    Burst1_onorig_targ_allele1 = 0;
    Burst1_onpara_targ_allele1 = 0;
    Burst1_off_targ_allele1 = 1;
    Burst1_onorig_targ_allele2 = 0;
    Burst1_onpara_targ_allele2 = 0;
    Burst1_off_targ_allele2 = 1;
    Burst1_on_para_allele1 = 0;
    Burst1_on_para_allele2 = 0;
    Burst1_on_orig_allele1 = 0;
    Burst1_on_orig_allele2 = 0;
    Burst1_is_mutated_allele1 = 0;
    Burst1_is_mutated_allele2 = 0;
    
    sp_1 = [A1,...
    Anonsense1,...
    Aprime1,...
    B1,...
    Burst1_onorig_targ_allele1,...
    Burst1_onpara_targ_allele1,...
    Burst1_off_targ_allele1,...
    Burst1_onorig_targ_allele2,...
    Burst1_onpara_targ_allele2,...
    Burst1_off_targ_allele2,...
    Burst1_on_para_allele1,...
    Burst1_on_para_allele2,...
    Burst1_on_orig_allele1,...
    Burst1_on_orig_allele2,...
    Burst1_is_mutated_allele1,...
    Burst1_is_mutated_allele2];
    
    ts_wtmut = [0 500];
    ic_wtmut = [0;1;0;1;0;0;...
        0;1;0;1;0;...
        0;0;1;0;0;1;...
        0;1;0;1;0];
    
    [t,y] = ode45(@(t,y) odefun_wtwt_B3state(t,y,ra_1), ts_wtmut, ic_wtmut); %odefun_wtwt was slightly different format but also worked
    
    ss_A1 = real(y(size(y,1),5));
    ss_Anons1 = real(y(size(y,1),6));
    ss_Aprim1 = real(y(size(y,1),11));
    ss_B1 = real(y(size(y,1),22));
    
    ssSP_ww(i,:) = [i,ss_A1,ss_Anons1,ss_Aprim1,ss_B1];
    
    [t,y] = ode45(@(t,y) odefun_wtmut_B3state(t,y,ra_1), ts_wtmut, ic_wtmut);
    
    ss_A1 = real(y(size(y,1),5));
    ss_Anons1 = real(y(size(y,1),6));
    ss_Aprim1 = real(y(size(y,1),11));
    ss_B1 = real(y(size(y,1),22));
    
    ssSP_wm(i,:) = [i,ss_A1,ss_Anons1,ss_Aprim1,ss_B1];
    
    [t,y] = ode45(@(t,y) odefun_mutmut_B3state(t,y,ra_1), ts_wtmut, ic_wtmut);
    
    ss_A1 = real(y(size(y,1),5));
    ss_Anons1 = real(y(size(y,1),6));
    ss_Aprim1 = real(y(size(y,1),11));
    ss_B1 = real(y(size(y,1),22));
    
    ssSP_mm(i,:) = [i,ss_A1,ss_Anons1,ss_Aprim1,ss_B1];
    

    
end

%%
ss_file = [outdir, 'steady_state_ODE45_wtwt_B3state.csv'];
ss_table = array2table(ssSP_ww, 'VariableNames', {'paramset', 'ss_A1','ss_Anons1','ss_Aprim1','ss_B1'});
writetable(ss_table, ss_file, 'Delimiter', ',')

ss_file = [outdir, 'steady_state_ODE45_wtmut_B3state.csv'];
ss_table = array2table(ssSP_wm, 'VariableNames', {'paramset', 'ss_A1','ss_Anons1','ss_Aprim1','ss_B1'});
writetable(ss_table, ss_file, 'Delimiter', ',')

ss_file = [outdir, 'steady_state_ODE45_mutmut_B3state.csv'];
ss_table = array2table(ssSP_mm, 'VariableNames', {'paramset', 'ss_A1','ss_Anons1','ss_Aprim1','ss_B1'});
writetable(ss_table, ss_file, 'Delimiter', ',')

%% het for all 100 paramsets - consecutive starting conditions
ssSP_ww = zeros(100,5);
ssSP_wm = zeros(100,5);
ssSP_mm = zeros(100,5);
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
    
    [t,y] = ode45(@(t,y) odefun_wtwt2(t,y,ra_1), ts_wtmut, ic_wtmut); %odefun_wtwt was slightly different format but also worked
    
    ss_A1 = real(y(size(y,1),5));
    ss_Anons1 = real(y(size(y,1),6));
    ss_Aprim1 = real(y(size(y,1),11));
    ss_B1 = real(y(size(y,1),22));
    
    ssSP_ww(i,:) = [i,ss_A1,ss_Anons1,ss_Aprim1,ss_B1];
    
    ic_wtmut = real(y(size(y,1),:));
    
    [t,y] = ode45(@(t,y) odefun_wtmut(t,y,ra_1), ts_wtmut, ic_wtmut);
    
    ss_A1 = real(y(size(y,1),5));
    ss_Anons1 = real(y(size(y,1),6));
    ss_Aprim1 = real(y(size(y,1),11));
    ss_B1 = real(y(size(y,1),22));
    
    ssSP_wm(i,:) = [i,ss_A1,ss_Anons1,ss_Aprim1,ss_B1];
    
    ic_wtmut = real(y(size(y,1),:));
    
    [t,y] = ode45(@(t,y) odefun_mutmut(t,y,ra_1), ts_wtmut, ic_wtmut);
    
    ss_A1 = real(y(size(y,1),5));
    ss_Anons1 = real(y(size(y,1),6));
    ss_Aprim1 = real(y(size(y,1),11));
    ss_B1 = real(y(size(y,1),22));
    
    ssSP_mm(i,:) = [i,ss_A1,ss_Anons1,ss_Aprim1,ss_B1];
    

    
end

%%
ss_file = [outdir, 'steady_state_ODE45_wtwt3.csv'];
ss_table = array2table(ssSP_ww, 'VariableNames', {'paramset', 'ss_A1','ss_Anons1','ss_Aprim1','ss_B1'});
writetable(ss_table, ss_file, 'Delimiter', ',')

ss_file = [outdir, 'steady_state_ODE45_wtmut3.csv'];
ss_table = array2table(ssSP_wm, 'VariableNames', {'paramset', 'ss_A1','ss_Anons1','ss_Aprim1','ss_B1'});
writetable(ss_table, ss_file, 'Delimiter', ',')

ss_file = [outdir, 'steady_state_ODE45_mutmut3.csv'];
ss_table = array2table(ssSP_mm, 'VariableNames', {'paramset', 'ss_A1','ss_Anons1','ss_Aprim1','ss_B1'});
writetable(ss_table, ss_file, 'Delimiter', ',')