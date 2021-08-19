% This is a hard-coded version of generateData_nitc_hap for running
% simulations of GRN output from a haploid 3-gene model {A, A', B} where A
% is the reference regulator, A' is a paralog of A, and B is the downstream
% target gene whose output is of interest. This model simulates nonsense-
% induced transcriptional compensation. It is based on the symmetric 
% network wrapper function generateData.m by Lea Schuh and Yogesh Goyal, 
% 2020.
% Ian Mellis, 2021

%This script produces Gillespie simulations according to random
%(latin hyper cube sampled) parameter sets
%
%INPUT:
%
%nruns:                 number of different runs (= number of parameter sets)
%n_species_upstr:       number of reference upstream regulators (A); hard coded to 1
%n_species_paralog:     number of paralog regulators (A'); hard coded to 1
%n_species_downstr:     number of downstream targets (B); hard coded to 1
%maxgillespie:          number of Gillespie simulated time units
%gen:               'yes' = generate new parameters (or 'no' = load parameters)
%init:              initial values for all species and initial values bursts
%                   (where 0 = 'off' and 1 = 'on')
%data:              if gen = 'no', which parameter matrix to use
%type:              normal - normal parameter space; constrained - constrained
%                   parameter space; asym - if asymmetric architecture;
%                   hard coded to normal; nitc - NITC architecture, an
%                   asymmetric subtype
%
%OUTPUT:
%
%S_outpar:           matrix of size 3*n_species times maxgillespie, where
%                   the columns give the gene produt count and DNA state
%                   (on/off) per time

%%
function generateData_3node_nitc_hap(nruns,n_species_upstr, n_species_paralog, n_species_downstr, maxgillespie,gen,init,data,type)

n_species_upstr = 1;
n_species_paralog = 1;
n_species_downstr = 1;
type = 'rigid_3';
gen = 'yes';

%overall specifications (for up to 10 species)
set_spec_orig = cell(1,1);
set_spec_para = cell(1,1);
set_Burston = cell(1,1);
set_Burstoff = cell(1,1);
set_prod = cell(1,1);
set_proddiff = cell(1,1);
set_deg = cell(1,1);
set_onbasal = cell(1,1);
set_ondep = cell(1,1);
set_off = cell(1,1);
set_onbasal_aprime = cell(1,1);
set_onbasal_b = cell(1,1);
set_ondep_prime = cell(1,1);
set_nitc = cell(1,1);
set_spec_targ = cell(1,1);

for iname = 1:1
    
    set_spec_orig{iname} = sprintf('A%d', iname);              %species
    set_Burston_orig{iname} = sprintf('Burst%d_on_orig', iname);                      %burst 'on' original regulator
    set_Burstoff_orig{iname} = sprintf('Burst%d_off_orig', iname);                    %burst 'off' original regulator
    set_prod{iname} = sprintf('prod%d', iname);                     %production rate
    set_proddiff{iname} = sprintf('proddiff%d', iname);             %difference in production rate between burst 'on' and burst 'off' (> 1)
    set_deg{iname} = sprintf('deg%d', iname);                       %degradation rate
    set_onbasal{iname} = sprintf('onbasal_a%d', iname);               %basal on-rate of burst
    set_ondep{iname} = sprintf('ondep%d', iname);                   %additional on-rate of target due to dependency to original regulator
    set_off{iname} = sprintf('off%d', iname);                       %off rate of burst
    set_onbasal_aprime{iname} = sprintf('onbasal_aprime%d', iname); %basal on-rate of burst of paralog
    set_onbasal_b{iname} = sprintf('onbasal_b%d', iname);           %basal on-rate of burst of downstream target
    set_ondep_prime{iname} = sprintf('ondep_prime', iname);                %additional on-rate of target due to dependency on paralog
    set_nitc{iname} = sprintf('nitc%d', iname);                            %additional on-rate of paralog due to dependency on NITC after mutation
    set_spec_para{iname} = sprintf('Aprime%d', iname);
    set_spec_targ{iname} = sprintf('B%d', iname); 
    set_Burston_para{iname} = sprintf('Burst%d_on_para', iname);                      %burst 'on' paralog
    set_Burstoff_para{iname} = sprintf('Burst%d_off_para', iname);                    %burst 'off' paralog
    set_Burston_targ{iname} = sprintf('Burst%d_on_targ', iname);                      %burst 'on' target
    set_Burstoff_targ{iname} = sprintf('Burst%d_off_targ', iname);                    %burst 'off' target

end

%set the initial values for Burst_off
for isetoff = 1: n_species_upstr
    if init.Burston_orig(isetoff) == 1
        init.Burstoff_orig(isetoff) = 0;
    else
        init.Burstoff_orig(isetoff) = 1;
    end
end

for isetoff = 1: n_species_paralog
    if init.Burston_para(isetoff) == 1
        init.Burstoff_para(isetoff) = 0;
    else
        init.Burstoff_para(isetoff) = 1;
    end
end

for isetoff = 1: n_species_downstr
    if init.Burston_targ(isetoff) == 1
        init.Burstoff_targ(isetoff) = 0;
    else
        init.Burstoff_targ(isetoff) = 1;
    end
end

%load all weakly-connected non-isomorphic symmetric networks of size
%n_species

if isequal(type,'normal') == 1
    %     load_M_iso = sprintf('/Volumes/MELANOMAII/Data/M_iso%d',n_species);
    load_M_iso = sprintf('./Example/M_iso%d',n_species);
    load(load_M_iso);
    nstruc = length(M_iso);                     %number of all possible networks considered
elseif isequal(type,'asym') == 1
    %     load('/Volumes/MELANOMAII/Example/AsymArchMat')
    load('./Example/AsymArchMat')
    M_iso{1} = Mat;
    nstruc = length(M_iso);
elseif isequal(type,'constrained') == 1
    %     load_M_iso = sprintf('/Volumes/MELANOMAII/Data/M_iso%d',n_species);
    load_M_iso = sprintf('./Example/M_iso%d',n_species);
    load(load_M_iso);
    nstruc = length(M_iso);
elseif isequal(type, 'nitc') == 1
    load('./M_iso_nitc_A1_Aprime1_B1'); % path check
    M_iso = M_iso_nitc_A1_Aprime1_B1;
    nstruc = length(M_iso);
end

M_orig_targ = [1]; % A->B edges; rows are A
M_para_targ = [1]; % A'->B edges; rows are A'
M_nitc_para = [1]; % Post-NITC A->A' edges; rows are NITC-A

if isequal(gen,'yes') == 1                  %generate new parameters by latin hypercube sampling method
    rng('shuffle');                         %receive different values upon new start
    
    min_range = [repmat(0.01,1,1),...       %production rate - same for all genes
        repmat(0.001,1,1),...               %degradation rate - same for all genes
        repmat(0.001,1,1),...               %basal on-rate of burst - original regulator
        repmat(0.1,1,1),...                 %Hill coefficient n (Hill function)
        repmat(0.1,1,1),...                 %additional on-rate due to dependency to other node
        repmat(0.01,1,1),...                %off rate of burst
        repmat(2,1,1),...                   %proddiff
        repmat(0.001,1,1),...               %basal on-rate of burst - paralogs
        repmat(0.001,1,1),...               %basal on-rate of burst - targets
        repmat(0.1,1,1),...                 %additional on-rate due to dependency to other node - paralog
        repmat(0.1,1,1)];                   %additional on-rate of paralog due to NITC
    
    max_range = [repmat(1,1,1),...          %production rate - same for all genes
        repmat(0.1,1,1),...                 %degradation rate - same for all genes
        repmat(0.1,1,1),...                 %basal on-rate of burst - original regulator
        repmat(10,1,1),...                  %Hill coefficient n (Hill function)
        repmat(1,1,1),...                   %additional on-rate due to dependency to other node
        repmat(0.1,1,1),...                 %off rate of burst
        repmat(100,1,1),...                 %proddiff
        repmat(0.1,1,1),...                 %basal on-rate of burst - paralogs
        repmat(0.1,1,1),...                 %basal on-rate of burst - targets
        repmat(1,1,1),...                   %additional on-rate due to dependency to other node - paralog
        repmat(1,1,1)];                     %additional on-rate of paralog due to NITC
    
    
    latinhyp = lhsdesign_modified(nruns, min_range, max_range);
    
    latinhyp_prod = repmat(latinhyp(:,1),1,n_species);
    latinhyp_deg = repmat(latinhyp(:,2),1,n_species);
    latinhyp_onbasal = repmat(latinhyp(:,3),1,n_species);
    latinhyp_n = repmat(latinhyp(:,4),1,n_species);
    latinhyp_ondep = repmat(latinhyp(:,5),1,n_species);
    latinhyp_off = repmat(latinhyp(:,6),1,n_species);
    latinhyp_proddiff = repmat(latinhyp(:,7),1,n_species);
    
%     new parameters
    latinhyp_onbasal_aprime = repmat(latinhyp(:,8),1,n_species);
    latinhyp_onbasal_b = repmat(latinhyp(:,9),1,n_species);
    latinhyp_ondep_prime = repmat(latinhyp(:,10),1,n_species);
    latinhyp_nitc = repmat(latinhyp(:,11),1,n_species);
    
    k = repmat(0.95*latinhyp(:,1)./latinhyp(:,2).*latinhyp(:,7),1,n_species); %0.95 of stationary on-state of system
    % do a search over k
    
else
    load(data)  %load parameter set
    
    if isequal(type,'normal') == 1
        R_outpar_par = Data50;
    elseif isequal(type,'asym') == 1
        R_outpar_par = Data50;
    elseif isequal(type,'constrained') == 1
        R_outpar_par = Data50C;
    end
    latinhyp_prod = repmat(R_outpar_par(:,1),1,n_species);
    latinhyp_deg = repmat(R_outpar_par(:,2),1,n_species);
    latinhyp_onbasal = repmat(R_outpar_par(:,6),1,n_species);
    latinhyp_n = repmat(R_outpar_par(:,8),1,n_species);
    latinhyp_ondep = repmat(R_outpar_par(:,3),1,n_species);
    latinhyp_off = repmat(R_outpar_par(:,4),1,n_species);
    latinhyp_proddiff = repmat(R_outpar_par(:,5),1,n_species);
    
    k = repmat(R_outpar_par(:,7),1,n_species);
    
    %     new parameters
    latinhyp_onbasal_aprime = repmat(R_outpar_par(:,9),1,n_species);
    latinhyp_onbasal_b = repmat(R_outpar_par(:,10),1,n_species);
    latinhyp_ondep_prime = repmat(R_outpar_par(:,11),1,n_species);
    latinhyp_nitc = repmat(R_outpar_par(:,11),1,n_species);
    
end


for istruc = 1:nstruc
    
    fclose all;
    clear functions
    
    istruc
    
    clear T_outpar S_outpar
    
    n_species = n_species_upstr + n_species_paralog + n_species_downstr;
    
    % needs to change? ***
    T_outpar = cell(1,nruns);
    S_outpar = cell(1,nruns);
%     S = zeros(nruns,3*n_species);
%     R = zeros(nruns,8*n_species);
%     P = zeros(nruns,4*n_species);
    S = zeros(nruns,3*n_species);
    R = zeros(nruns,11*n_species);
    P = zeros(nruns,4*n_species);
    
    for iruns = 1:nruns
        
        fclose all;
        clear functions
%         clearvars -except nruns maxgillespie n_species init gen...
%             set_spec set_Bon set_Boff set_prod set_proddiff set_deg...
%             set_onbasal set_ondep set_off M_iso nstruc...
%             latinhyp_prod latinhyp_deg latinhyp_onbasal latinhyp_n...
%             latinhyp_ondep latinhyp_off latinhyp_proddiff k...
%             iruns istruc S P R type

        clearvars -except nruns maxgillespie n_species init gen...
            n_species_downstr n_species_upstr n_species_paralog...
            set_prod set_proddiff set_deg...
            set_spec_orig set_Burston_orig set_Burstoff_orig...
            set_onbasal set_ondep set_off set_onbasal_aprime...
            set_onbasal_b set_ondep_prime set_nitc set_spec_para...
            set_spec_targ set_Burston_para set_Burstoff_para...
            set_Burston_targ set_Burstoff_targ M_iso nstruc...
            M_orig_targ M_para_targ M_nitc_para...
            latinhyp_prod latinhyp_deg latinhyp_onbasal latinhyp_n...
            latinhyp_ondep latinhyp_off latinhyp_proddiff k...
            iruns istruc S P R type
        
        parmat = M_iso{istruc};
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %create txt file for specific network/sampled parameter set
        %         fileID = fopen('/Volumes/MELANOMAII/Example/gillespie_bursts.txt','w');
        fileID = fopen('./Example/gillespie_bursts.txt','w');
        
        fprintf(fileID,'%d\n', maxgillespie);
        
%         set_prod for each species
        for prod_txt = 1:n_species_upstr
            fprintf(fileID,'A_%s = %f : = %s \n', set_prod{prod_txt}, latinhyp_prod(iruns,prod_txt), set_spec_orig{prod_txt}) ;
        end
        
        for prod_txt = 1:n_species_paralog
            fprintf(fileID,'Aprime_%s = %f : = %s \n', set_prod{prod_txt}, latinhyp_prod(iruns,prod_txt), set_spec_para{prod_txt}) ;
        end
        
        for prod_txt = 1:n_species_downstr
            fprintf(fileID,'B_%s = %f : = %s \n', set_prod{prod_txt}, latinhyp_prod(iruns,prod_txt), set_spec_targ{prod_txt}) ;
        end
        
%         set_deg for each species
        for deg_txt = 1:n_species_upstr
            fprintf(fileID,'A_%s = %f : %s = \n', set_deg{deg_txt}, latinhyp_deg(iruns,deg_txt), set_spec_orig{deg_txt}) ;
        end
        
        for deg_txt = 1:n_species_paralog
            fprintf(fileID,'Aprime_%s = %f : %s = \n', set_deg{deg_txt}, latinhyp_deg(iruns,deg_txt), set_spec_para{deg_txt}) ;
        end
        
        for deg_txt = 1:n_species_downstr
            fprintf(fileID,'B_%s = %f : %s = \n', set_deg{deg_txt}, latinhyp_deg(iruns,deg_txt), set_spec_targ{deg_txt}) ;
        end
        
%         ondep for each downstream target B regulated by A 
        for ondep_txt = 1:n_species_downstr
            fprintf(fileID,'B_%s = %f : %s = %s\n', set_ondep{ondep_txt}, latinhyp_ondep(iruns,ondep_txt), set_Burstoff_targ{ondep_txt}, set_Burston_targ{ondep_txt});
        end
        
%         ondep for each downstream target B regulated by Aprime
        for ondep_txt = 1:n_species_downstr
            fprintf(fileID,'B_%s = %f : %s = %s\n', set_ondep_prime{ondep_txt}, latinhyp_ondep_prime(iruns,ondep_txt), set_Burstoff_targ{ondep_txt}, set_Burston_targ{ondep_txt});
        end
        
%         NITC-based ondep for Aprime regulated by NITC-type A
        for nitc_txt = 1:nspecies
            fprintf(fileID,'Aprime%s = %f : %s = %s\n', set_nitc{nitc_txt}, latinhyp_nitc(iruns,nitc_txt), set_Burstoff_para{ondep_txt}, set_Burston_para{ondep_txt});
        end
        
%         off for each species
        for on_txt = 1:n_species_upstr
            fprintf(fileID,'A_%s = %f : %s = %s\n', set_off{on_txt}, latinhyp_off(iruns,on_txt), set_Bon{on_txt}, set_Boff{on_txt});
        end
        
        for on_txt = 1:n_species_paralog
            fprintf(fileID,'Aprime_%s = %f : %s = %s\n', set_off{on_txt}, latinhyp_off(iruns,on_txt), set_Bon{on_txt}, set_Boff{on_txt});
        end
        
        for on_txt = 1:n_species_downstr
            fprintf(fileID,'B_%s = %f : %s = %s\n', set_off{on_txt}, latinhyp_off(iruns,on_txt), set_Bon{on_txt}, set_Boff{on_txt});
        end
        
        
        fprintf(fileID,'%s\n', 'Production difference');
        
        for proddiff_txt = 1:n_species_upstr
            fprintf(fileID,'A_%s = %d\n', set_proddiff {proddiff_txt}, latinhyp_proddiff(iruns,proddiff_txt));
        end
        
        for proddiff_txt = 1:n_species_paralog
            fprintf(fileID,'Aprime_%s = %d\n', set_proddiff {proddiff_txt}, latinhyp_proddiff(iruns,proddiff_txt));
        end
        
        for proddiff_txt = 1:n_species_downstr
            fprintf(fileID,'B_%s = %d\n', set_proddiff {proddiff_txt}, latinhyp_proddiff(iruns,proddiff_txt));
        end
        
        
        fprintf(fileID,'%s\n', 'Basal values');
        
        for onbasak_txt = 1:n_species
            fprintf(fileID,'%s = %d\n', set_onbasal{onbasak_txt}, latinhyp_onbasal(iruns,onbasak_txt));
        end
        
        for onbasal_aprime_txt = 1:n_species
            fprintf(fileID,'%s = %d\n', set_onbasal_aprime{onbasal_aprime_txt}, latinhyp_onbasal_aprime(iruns,onbasal_aprime_txt));
        end
        
        for onbasal_b_txt = 1:n_species
            fprintf(fileID,'%s = %d\n', set_onbasal_b{onbasal_b_txt}, latinhyp_onbasal_b(iruns,onbasal_b_txt));
        end
        
        
        fprintf(fileID,'%s\n', 'Hill function k');
%         how to handle?
        for spec_txt = 1:n_species
            fprintf(fileID,'k%s = %d\n', set_spec{spec_txt}, k(iruns,spec_txt));
        end
        
        fprintf(fileID,'%s\n', 'Hill function n');
        
        for n_txt = 1:n_species
            fprintf(fileID,'n%s = %d\n', set_spec{n_txt}, latinhyp_n(iruns,n_txt));
        end
        
        
        fprintf(fileID,'%s\n', 'Initial values');
        
        for initspec_txt = 1:n_species_upstr
            fprintf(fileID,'%s = %d\n', set_spec_orig{initspec_txt}, init.spec_A(initspec_txt));
        end
        
        for initspec_txt = 1:n_species_paralog
            fprintf(fileID,'%s = %d\n', set_spec_para{initspec_txt}, init.spec_Aprime(initspec_txt));
        end
        
        for initspec_txt = 1:n_species_downstr
            fprintf(fileID,'%s = %d\n', set_spec_targ{initspec_txt}, init.spec_B(initspec_txt));
        end
        
        for initB_txt = 1:n_species_upstr
            fprintf(fileID,'%s = %d\n', set_Burstoff_orig{initB_txt}, init.Burstoff_A(initB_txt));
            fprintf(fileID,'%s = %d\n', set_Burston_orig{initB_txt}, init.Burston_A(initB_txt));
        end
        
        for initB_txt = 1:n_species_paralog
            fprintf(fileID,'%s = %d\n', set_Burstoff_para{initB_txt}, init.Burstoff_Aprime(initB_txt));
            fprintf(fileID,'%s = %d\n', set_Burston_para{initB_txt}, init.Burston_Aprime(initB_txt));
        end
        
        for initB_txt = 1:n_species_downstr
            fprintf(fileID,'%s = %d\n', set_Burstoff_targ{initB_txt}, init.Burstoff_B(initB_txt));
            fprintf(fileID,'%s = %d\n', set_Burston_targ{initB_txt}, init.Burston_B(initB_txt));
        end
        
        fprintf(fileID,'%s\n', 'Network');
        for net_txt = 1:n_species_upstr
            fprintf(fileID,'%d\n',M_orig_targ(net_txt,:));
        end
        for net_txt = 1:n_species_paralog
            fprintf(fileID,'%d\n',M_para_targ(net_txt,:));
        end
        for net_txt = 1:n_species_upstr
            fprintf(fileID,'%d\n',M_nitc_para(net_txt,:));
        end
        fclose(fileID);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        fclose all;
        clear functions
        
        %         tok = make_param_mex_bursts('/Volumes/MELANOMAII/Example/gillespie_bursts.txt',...
        %             '/Volumes/MELANOMAII/Example/gillespie_bursts',n_species,iruns);
        tok = make_param_mex_bursts('./Example/gillespie_bursts.txt',...
            './Example/gillespie_bursts',n_species,iruns);
        
        gillespie_burstsparams  %load parameters set by txt
        
        S(iruns,:) = transpose(species);
        R(iruns,:) = transpose(rates);
        P(iruns,:) = transpose(propensity);
        
    end
    
    if isequal(type,'normal') == 1
        %         R_save = sprintf('/Volumes/MELANOMAII/Example/R_outpar%d_%d',n_species,istruc);
        %         S_save = sprintf('/Volumes/MELANOMAII/Example/S_outpar%d_%d',n_species,istruc);
        %         P_save = sprintf('/Volumes/MELANOMAII/Example/P_outpar%d_%d',n_species,istruc);
        R_save = sprintf('./Example/R_outpar%d_%d',n_species,istruc);
        S_save = sprintf('./Example/S_outpar%d_%d',n_species,istruc);
        P_save = sprintf('./Example/P_outpar%d_%d',n_species,istruc);
    elseif isequal(type,'asym') == 1
        %         R_save = sprintf('/Volumes/MELANOMAII/Example/R_outparAsymArch%d_%d',n_species,istruc);
        %         S_save = sprintf('/Volumes/MELANOMAII/Example/S_outparAsymArch%d_%d',n_species,istruc);
        %         P_save = sprintf('/Volumes/MELANOMAII/Example/P_outparAsymArch%d_%d',n_species,istruc);
        R_save = sprintf('./Example/R_outparAsymArch%d_%d',n_species,istruc);
        S_save = sprintf('./Example/S_outparAsymArch%d_%d',n_species,istruc);
        P_save = sprintf('./Example/P_outparAsymArch%d_%d',n_species,istruc);
    elseif isequal(type,'constrained') == 1
        %         R_save = sprintf('/Volumes/MELANOMAII/Example/R_outparC%d_%d',n_species,istruc);
        %         S_save = sprintf('/Volumes/MELANOMAII/Example/S_outparC%d_%d',n_species,istruc);
        %         P_save = sprintf('/Volumes/MELANOMAII/Example/P_outparC%d_%d',n_species,istruc);
        R_save = sprintf('./Example/R_outparC%d_%d',n_species,istruc);
        S_save = sprintf('./Example/S_outparC%d_%d',n_species,istruc);
        P_save = sprintf('./Example/P_outparC%d_%d',n_species,istruc);
    end
    
    save(R_save,'R','-v7.3');
    save(S_save,'S','-v7.3');
    save(P_save,'P','-v7.3');
    
    for kmem = 1:1
        
        parfor jruns = 1:nruns
            %         for jruns = 1:nruns
            
            par_spec = S(100*kmem-100+jruns,:);
            par_rates = R(100*kmem-100+jruns,:);
            par_prop = P(100*kmem-100+jruns,:);
            
            [times,savespecies] = gillespie_burstshistomex(0,par_spec,par_rates,par_prop,sum(clock*100),maxgillespie,maxgillespie);
            
            S_outpar{jruns} = savespecies;
            
            times = [];
            savespecies = [];
            
        end
        
        if isequal(type,'normal') == 1
            %             S_save = sprintf('/Volumes/MELANOMAII/Example/S_outpar%d_%d_%d',n_species,istruc,kmem);
            S_save = sprintf('./Example/S_outpar%d_%d_%d',n_species,istruc,kmem);
        elseif isequal(type,'asym') == 1
            %             S_save = sprintf('/Volumes/MELANOMAII/Example/S_outparAsymArch%d_%d_%d',n_species,istruc,kmem);
            S_save = sprintf('./Example/S_outparAsymArch%d_%d_%d',n_species,istruc,kmem);
        elseif isequal(type,'constrained') == 1
            %             S_save = sprintf('/Volumes/MELANOMAII/Example/S_outparC%d_%d_%d',n_species,istruc,kmem);
            S_save = sprintf('./Example/S_outparC%d_%d_%d',n_species,istruc,kmem);
        end
        
        save(S_save,'S_outpar','-v7.3');
        
        clear S_outpar
    end
end

end
