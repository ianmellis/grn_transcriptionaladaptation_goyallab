% Parameter file


% Simulation parameter values

% Maximum number of Gillespie steps
maxgillespiesteps = 10000;

% Initial time
currT = 0;


% Initial values for species

A1 = 1;
Aprime1 = 1;
B1 = 1;
Burst1_off_orig = 0;
Burst1_on_orig = 1;
Burst1_off_para = 1;
Burst1_on_para = 0;
Burst1_off_targ = 1;
Burst1_on_targ = 0;

% EVERYTHING BELOW IS BOOKKEEPING; DO NOT ALTER!


% Reaction rates

% Rxn: = A1 
A_prod1  = 0.831786;
% Rxn: = Aprime1 
Aprime_prod1  = 0.831786;
% Rxn: = B1 
B_prod1  = 0.831786;
% Rxn: A1 = 
A_deg1  = 0.041163;
% Rxn: Aprime1 = 
Aprime_deg1  = 0.041163;
% Rxn: B1 = 
B_deg1  = 0.041163;
% Rxn: Burst1_off_targ = Burst1_on_targ
B_ondep1  = 0.891703;
% Rxn: Burst1_off_targ = Burst1_on_targ
B_ondep_prime  = 0.0685844;
% Rxn: Burst1_off_para = Burst1_on_para
Aprimenitc1  = 0.247986;
% Rxn: Burst1_on_orig = Burst1_off_orig
A_off1  = 0.035732;
% Rxn: Burst1_on_para = Burst1_off_para
Aprime_off1  = 0.035732;
% Rxn: Burst1_on_targ = Burst1_off_targ
B_off1  = 0.035732;
A_proddiff1  = 45.9486;
Aprime_proddiff1  = 45.9486;
B_proddiff1  = 45.9486;
onbasal_a1  = 0.0677509;
onbasal_aprime1  = 0.0801551;
onbasal_b1  = 0.0939215;
kA1  = 464.244;
kAprime1  = 464.244;
kB1  = 464.244;
nA1  = 9.20147;
nAprime1  = 9.20147;
nB1  = 9.20147;
numrxns = 12;
nspecies = 9;
species = zeros(nspecies,1);
species(1) = A1;
species(2) = Aprime1;
species(3) = B1;
species(4) = Burst1_off_targ;
species(5) = Burst1_on_targ;
species(6) = Burst1_off_para;
species(7) = Burst1_on_para;
species(8) = Burst1_on_orig;
species(9) = Burst1_off_orig;
rates = zeros(numrxns,1);
rates(1) = A_prod1 ;
rates(2) = Aprime_prod1 ;
rates(3) = B_prod1 ;
rates(4) = A_deg1 ;
rates(5) = Aprime_deg1 ;
rates(6) = B_deg1 ;
rates(7) = B_ondep1 ;
rates(8) = B_ondep_prime ;
rates(9) = Aprimenitc1 ;
rates(10) = A_off1 ;
rates(11) = Aprime_off1 ;
rates(12) = B_off1 ;
rates(13) = A_proddiff1 ;
rates(14) = Aprime_proddiff1 ;
rates(15) = B_proddiff1 ;
rates(16) = onbasal_a1 ;
rates(17) = onbasal_aprime1 ;
rates(18) = onbasal_b1 ;
rates(19) = kA1 ;
rates(20) = kAprime1 ;
rates(21) = kB1 ;
rates(22) = nA1 ;
rates(23) = nAprime1 ;
rates(24) = nB1 ;
y0 = zeros(nspecies,1);
y0(1) = A1;
y0(2) = Aprime1;
y0(3) = B1;
y0(4) = Burst1_off_targ;
y0(5) = Burst1_on_targ;
y0(6) = Burst1_off_para;
y0(7) = Burst1_on_para;
y0(8) = Burst1_on_orig;
y0(9) = Burst1_off_orig;
% Intialize the propensities...
propensity(1) = A_prod1 *A_proddiff1 *Burst1_on_targ+A_prod1 *Burst1_off_targ;
propensity(2) = Aprime_prod1 *Aprime_proddiff1 *Burst1_on_para+Aprime_prod1 *Burst1_off_para;
propensity(3) = B_prod1 *B_proddiff1 *Burst1_off_orig+B_prod1 *Burst1_on_orig;
propensity(4) = A_deg1 *A1;
propensity(5) = Aprime_deg1 *Aprime1;
propensity(6) = B_deg1 *B1;
% needs on propensity for A and a single on propensity for B (as well as
% updated Aprime on propensity
propensity(7) = B_ondep1 *((A1^nA1 )/(kA1 ^nA1 +A1^nA1 ))*Burst1_off_targ + B_ondep_prime *((Aprime1^nAprime1 )/(kAprime1 ^nAprime1 +Aprime1^nAprime1 ))*Burst1_off_targ + onbasal_a1 *Burst1_off_targ;
propensity(8) = B_ondep_prime *((A1^nA1 )/(kA1 ^nA1 +A1^nA1 ))*Burst1_off_para+onbasal_aprime1 *Burst1_off_para;
propensity(9) = Aprimenitc1 *((A1^nA1 )/(kA1 ^nA1 +A1^nA1 ))*Burst1_on_orig+onbasal_b1 *Burst1_on_orig;
propensity(10) = A_off1 *Burst1_on_targ;
propensity(11) = Aprime_off1 *Burst1_on_para;
propensity(12) = B_off1 *Burst1_off_orig;

