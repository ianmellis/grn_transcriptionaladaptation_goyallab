%% Generate data - TF-promoter-Bonly model
% diploid, alleles independent, NITC affects ancestral and paralog
% substrate competition between ancestral and paralog at target promoter

outdir = '/Volumes/IAMYG1/grn_nitc_data/v1.5.2/';

if ~exist(outdir, 'dir')
    mkdir(outdir)
end

%% LHS
% Larsson et al., 2019 estimated burst kinetics transcriptome-wide using
% single-cell RNA-seq data. They inferred r_on, r_off, and r_prod relative
% to r_deg for individual alleles. 
% 
% Fixed parameter values: r_deg = 1, x= 0.5. What to do about r_onbasal?
% 
% search over: basal_nitc_on_ratio, onbasalA1_off_ratio,
% A1_Aprime1_boundOn_ratio, A1_Aprime_binding_ratio,
% bindbyA1_unbindbyA1_ratio, r_prod_on, r_bound_byA1_B1, r_bind_byA1_B1, 
% r_unbind_byA1_B1
rng(8734);

nruns = 100;

min_range = [repmat(0.1,1,1),...        %basal_nitc_on_ratio (0.1,10) - is 10 high enough?
    repmat(0.01,1,1),...                %onbasalA1_off_ratio (0.01,1) - per Larsson this is the bulk of the distribution (very quick off-burst)
    repmat(0.1,1,1),...                 %A1_Aprime1_boundOn_ratio (0.1,10)
    repmat(0.1,1,1),...                 %A1_Aprime_binding_ratio (0.1,10)
    repmat(1,1,1),...                   %bindbyA1_unbindbyA1_ratio (1,100)
    repmat(1,1,1),...                   %r_prod_on (1,1000)
    repmat(0.1,1,1),...                 %r_bound_byA1_B1 (0.1,10)
    repmat(0.1,1,1),...                 %r_bind_byA1_B1 (0.1,10)
    repmat(0.1,1,1),...                 %n (Hill coefficient n) (0.1,10) too large a range?
    repmat(0.1,1,1)];                   %r_onbasal_A1 (0.1, 10)
%                                       

max_range = [repmat(10,1,1),...         %basal_nitc_on_ratio (0.1,10) - is 10 high enough?
    repmat(1,1,1),...                   %onbasalA1_off_ratio (0.01,1) - per Larsson this is the bulk of the distribution (very quick off-burst)
    repmat(10,1,1),...                  %A1_Aprime1_boundOn_ratio (0.1,10)
    repmat(10,1,1),...                  %A1_Aprime_binding_ratio (0.1,10)
    repmat(100,1,1),...                 %bindbyA1_unbindbyA1_ratio (1,100)
    repmat(1000,1,1),...                %r_prod_on (1,1000)
    repmat(10,1,1),...                  %r_bound_byA1_B1 (0.1,10)
    repmat(10,1,1),...                  %r_bind_byA1_B1 (0.1,10)
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
    'A1_Aprime1_boundOn_ratio',...
    'A1_Aprime_binding_ratio',...
    'bindbyA1_unbindbyA1_ratio',...
    'r_prod_on',...
    'r_bound_byA1_B1',...
    'r_bind_byA1_B1',...
    'Hill_coefficient_n'};

%%
writetable(lhs_1_s,lhs_1_f,'Delimiter',',')
%%
tic
parfor i = 1:nruns
    
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
    r_onbasal_A1 = latinhyp(i, 10);
    r_onbasal_other = 0;
    r_nitc_byAnonsense1_A1 = r_onbasal_A1/basal_nitc_on_ratio;
    r_nitc_byAnonsense1_Anonsense1 = r_nitc_byAnonsense1_A1;
    r_nitc_byAnonsense1_Aprime1 = r_nitc_byAnonsense1_A1;
    r_bound_byA1_B1 = latinhyp(i, 7);
    r_bound_byAprime1_B1 = r_bound_byA1_B1/A1_Aprime1_boundOn_ratio;
    r_off = r_onbasal_A1/onbasalA1_off_ratio;
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
    
    d_Aprime_B1 = A1_Aprime1_on_ratio;
    
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
    
    
    ra_1 = [r_prodbasal_A1,...0];
        r_prodbasal_Anonsense1,...1];
        r_prodbasal_Aprime1,...2];
        r_prodbasal_B1,...3];
        
        r_prodon_A1,...4];
        r_prodon_Anonsense1,...5];
        r_prodon_Aprime1,...6];
        r_prodon_B1,...7];
        
        d_Aprime1_B1,...8];
        
        r_deg_A1,...9];
        r_deg_Anonsense1,...10];
        r_deg_Aprime1,...11];
        r_deg_B1,...12];
        
        r_onbasal_A1,...13];
        r_onbasal_Anonsense1,...14];
        r_onbasal_Aprime1,...15];
        r_onbasal_B1,...16];
        
        r_nitc_byAnonsense1_A1,...17];
        r_nitc_byAnonsense1_Anonsense1,...18];
        r_nitc_byAnonsense1_Aprime1,...19];
        r_addon_byA1_B1,...20];
        r_addon_byAprime1_B1,...21];
        
        r_off_A1,...22];
        r_off_Anonsense1,...23];
        r_off_Aprime1,...24];
        r_offorig_B1,...25];
        r_offpara_B1,...26];
        
        k_A1 ,...27];
        k_Anonsense1 ,...28];
        k_Aprime1 ,...29];
        k_B1 ,...30];
        
        n_A1 ,...31];
        n_Anonsense1 ,...32];
        n_Aprime1 ,...33];
        n_B1];%  = rates[34]];
    
    
    ra_1_s = array2table(ra_1);
    ra_1_s.Properties.VariableNames = {'r_prodbasal_A1',...
        'r_prodbasal_Anonsense1',...
        'r_prodbasal_Aprime1',...
        'r_prodbasal_B1',...
        'r_prodon_A1',...
        'r_prodon_Anonsense1',...
        'r_prodon_Aprime1',...
        'r_prodon_B1',...
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
        'r_bound_byA1_B1',...
        'r_bound_byAprime1_B1',...
        'r_off_A1',...
        'r_off_Anonsense1',...
        'r_off_Aprime1',...
        'r_off_B1',...
        'r_bind_byA1_B1',...
        'r_bind_byAprime1_B1',...
        'k_A1',...
        'k_Anonsense1',...
        'k_Aprime1',...
        'k_B1',...
        'n_A1',...
        'n_Anonsense1',...
        'n_Aprime1',...
        'n_B1',...
        'r_unbind_byA1_B1',...
        'r_unbind_byAprime1_B1'};
    
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
    
    % Initialize propensities
    % force production of ancestral gene as first reaction, arbitrarily
    pr_1 = [1,zeros(1,27)];
        
    % Run sim
    % runs infinitely when both alleles mutated
    
    maxgillespie = 300000;
    
    %
    
    [times,savespecies] = gillespie_burstshistomex_TFpromoterB(0,sp_1,ra_1,pr_1,sum(clock*100),maxgillespie,maxgillespie);
    
    sp_q300 = savespecies(:,[599:300:100000, 100599:300:200000, 200599:300:300000]);
    sp_q300(19,:) = [599:300:100000, 100599:300:200000, 200599:300:300000];
    
    % Save results
%     ddir = './nitc_TFpromoter_3node_v1.5/';
    
    sp_f = [outdir, 'initialsim_species', num2str(paramsetnum), '.csv'];
    sp_q300_f = [outdir, 'initialsim_species', num2str(paramsetnum), '_q300.csv'];
    ra_f = [outdir, 'initialsim_rates', num2str(paramsetnum), '.csv'];
    
    writetable(array2table(transpose(savespecies),'VariableNames',{'A1',...
        'Anonsense1',...
        'Aprime1',...
        'B1',...
        'Promoter1_unbound_targ_allele1',...
        'Promoter1_boundbyorig_targ_allele1',...
        'Promoter1_boundbypara_targ_allele1',...
        'Promoter1_unbound_targ_allele2',...
        'Promoter1_boundbyorig_targ_allele2',...
        'Promoter1_boundbypara_targ_allele2',...
        'Burst1_on_targ_allele1',...
        'Burst1_on_targ_allele2',...
        'Burst1_on_para_allele1',...
        'Burst1_on_para_allele2',...
        'Burst1_on_orig_allele1',...
        'Burst1_on_orig_allele2',...
        'Burst1_is_mutated_allele1',...
        'Burst1_is_mutated_allele2'}), sp_f, 'Delimiter', ',') 
    writetable(ra_1_s,ra_f,'Delimiter',',')
    
    writetable(array2table(transpose(sp_q300),'VariableNames',{'A1',...
        'Anonsense1',...
        'Aprime1',...
        'B1',...
        'Promoter1_unbound_targ_allele1',...
        'Promoter1_boundbyorig_targ_allele1',...
        'Promoter1_boundbypara_targ_allele1',...
        'Promoter1_unbound_targ_allele2',...
        'Promoter1_boundbyorig_targ_allele2',...
        'Promoter1_boundbypara_targ_allele2',...
        'Burst1_on_targ_allele1',...
        'Burst1_on_targ_allele2',...
        'Burst1_on_para_allele1',...
        'Burst1_on_para_allele2',...
        'Burst1_on_orig_allele1',...
        'Burst1_on_orig_allele2',...
        'Burst1_is_mutated_allele1',...
        'Burst1_is_mutated_allele2',...
        'time'}), sp_q300_f, 'Delimiter', ',') 
    %
    
end
toc;

%% LHS from Lea
% % search over: r_prod, r_deg
% min_range = [repmat(0.01,1,1),...       %production rate - same for all genes
%     repmat(0.001,1,1),...               %degradation rate - same for all genes
%     repmat(0.001,1,1),...               %basal on-rate of burst - original regulator
%     repmat(0.1,1,1),...                 %Hill coefficient n (Hill function)
%     repmat(0.1,1,1),...                 %additional on-rate due to dependency to other node
%     repmat(0.01,1,1),...                %off rate of burst
%     repmat(2,1,1),...                   %proddiff
%     repmat(0.001,1,1),...               %basal on-rate of burst - paralogs; a multiple of the basal on-rate of original gene
%     repmat(0.001,1,1),...               %basal on-rate of burst - targets; a multiple of the basal on-rate of original regulator gene
%     repmat(0.001,1,1),...               %additional on-rate due to dependency to other node - paralog; a multiple of the r_add of original gene
%     repmat(0.01,1,1)];                  %additional on-rate of paralog due to NITC
% 
% max_range = [repmat(1,1,1),...          %production rate - same for all genes
%     repmat(0.1,1,1),...                 %degradation rate - same for all genes
%     repmat(0.1,1,1),...                 %basal on-rate of burst - original regulator
%     repmat(10,1,1),...                  %Hill coefficient n (Hill function)
%     repmat(1,1,1),...                   %additional on-rate due to dependency to other node
%     repmat(0.1,1,1),...                 %off rate of burst
%     repmat(100,1,1),...                 %proddiff
%     repmat(1,1,1),...                   %basal on-rate of burst - paralogs; a multiple of the basal on-rate of original gene
%     repmat(1,1,1),...                   %basal on-rate of burst - targets; a multiple of the basal on-rate of original regulator gene
%     repmat(2,1,1),...                   %additional on-rate due to dependency to other node - paralog; a multiple of the r_add of original gene
%     repmat(1,1,1)];                     %additional on-rate of paralog due to NITC
% 
% 
% latinhyp = lhsdesign_modified(nruns, min_range, max_range);
% 
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
% k = repmat(x*latinhyp(:,1)./latinhyp(:,2).*latinhyp(:,7),1,n_species); %0.95 of stationary on-state of system
% % do a search over x