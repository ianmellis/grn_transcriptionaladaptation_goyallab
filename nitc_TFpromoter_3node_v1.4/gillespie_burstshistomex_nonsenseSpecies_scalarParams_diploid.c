#include <math.h>
#include <matrix.h>

//Here's the stuff for the random number generator
static unsigned long jz,jsr=123456789;

#define SHR3 (jz=jsr, jsr^=(jsr<<13), jsr^=(jsr>>17), jsr^=(jsr<<5),jz+jsr)
#define UNI (.5 + (signed) SHR3*.2328306e-9)
#define IUNI SHR3

static long hz;
static unsigned long iz, kn[128], ke[256];
static double wn[128],fn[128], we[256],fe[256];

#define RNOR (hz=SHR3, iz=hz&127, (fabs(hz)<kn[iz])? hz*wn[iz] : nfix())
#define REXP (jz=SHR3, iz=jz&255, (    jz <ke[iz])? jz*we[iz] : efix())


double nfix(void);
double efix(void);
double ran_expo(double lambda);
void zigset(unsigned long jsrseed);



//END random number generator

#ifndef __STDIO_H__
#include <stdio.h>
#endif

#ifndef mex_h
#include "mex.h"
#endif

#define TIMES_OUT     plhs[0]
#define SPECIES_OUT   plhs[1]

#define CURRT_IN          prhs[0]
#define SPECIES_IN        prhs[1]
#define RATES_IN          prhs[2]
#define PROPENSITIES_IN   prhs[3]
#define SEED_IN           prhs[4]
#define M_IN              prhs[5]
#define TOTALT_IN         prhs[6]

//INCLUDE NSPECIES HERE
#define NSPECIES 20

//INCLUDE NUMRXNS HERE
#define NUMRXNS 24


//void gillespie(long m, double *times_out,double *species_out,double currT,double *species,double *rates,double *propensities,double *rand1,double *rand2)

void gillespie(long m, double *times_out,double *species_out,double currT,double *species,double *rates,double *propensities,long seed, double totalt)
{
  double alpha, deltaT, p;
  double cumpropensities[NUMRXNS];

  double deltaTsave;
  long savecount;
  long savecountcheck;
  long totaliterations;

//INSERT ALL VARIABLE DECLARATIONS HERE
  long A1,A1_nonsense,Aprime1,B1,Burst1_off_targ_allele1,Burst1_on_targ_allele1,Burst1_off_para_allele1,Burst1_on_para_allele1,Burst1_on_orig_allele1,Burst1_off_orig_allele1,Burst1_is_mutated_allele1,Burst1_not_mutated_allele1,Burst1_off_targ_allele2,Burst1_on_targ_allele2,Burst1_off_para_allele2,Burst1_on_para_allele2,Burst1_on_orig_allele2,Burst1_off_orig_allele2,Burst1_is_mutated_allele2,Burst1_not_mutated_allele2;
  double A_prod1 ,Anonsense_prod1 , Aprime_prod1 ,B_prod1 ,A_deg1 ,Anonsense_deg1 ,Aprime_deg1 ,B_deg1 ,B_ondep1 ,B_ondep_prime ,Aprimenitc1 ,A_off1 ,Aprime_off1 ,B_off1 ,A_proddiff1 ,Aprime_proddiff1 ,B_proddiff1 ,onbasal_a1 ,onbasal_aprime1 ,onbasal_b1 ,kA1 ,kA1_nonsense ,kAprime1 ,kB1 ,nA1 ,nA1_nonsense ,nAprime1 ,nB1 ;

  long i, j, k;
//UNPACK ALL SPECIES HERE
  A1 = (long)species[0];
  A1_nonsense = (long)species[1];
  Aprime1 = (long)species[2];
  B1 = (long)species[3];
  
  Promoter1_unbound_targ_allele1 = (long)species[4];
  Promoter1_boundbyorig_targ_allele1 = (long)species[4];
  Promoter1_boundbypara_targ_allele1 = (long)species[4];
  Promoter1_unbound_targ_allele2 = (long)species[4];
  Promoter1_boundbyorig_targ_allele2 = (long)species[4];
  Promoter1_boundbypara_targ_allele2 = (long)species[4];
  Promoter1_unbound_para_allele1 = (long)species[4];
  Promoter1_boundbynonsense_para_allele1 = (long)species[4];
  Promoter1_unbound_para_allele2 = (long)species[4];
  Promoter1_boundbynonsense_para_allele2 = (long)species[4];
  Promoter1_unbound_orig_allele1 = (long)species[4];
  Promoter1_boundbynonsense_orig_allele1 = (long)species[4];
  Promoter1_unbound_orig_allele2 = (long)species[4];
  Promoter1_boundbynonsense_orig_allele2 = (long)species[4];

  Burst1_off_targ_allele1 = (long)species[4];
  Burst1_on_targ_allele1 = (long)species[5];
  Burst1_off_targ_allele2 = (long)species[6];
  Burst1_on_targ_allele2 = (long)species[7];
  Burst1_off_para_allele1 = (long)species[8];
  Burst1_on_para_allele1 = (long)species[9];
  Burst1_off_para_allele2 = (long)species[10];
  Burst1_on_para_allele2 = (long)species[11];
  Burst1_on_orig_allele1 = (long)species[12];
  Burst1_off_orig_allele1 = (long)species[13];
  Burst1_on_orig_allele2 = (long)species[14];
  Burst1_off_orig_allele2 = (long)species[15];
  Burst1_is_mutated_allele1 = (long)species[16];
  Burst1_not_mutated_allele1 = (long)species[17];
  Burst1_is_mutated_allele2 = (long)species[18];
  Burst1_not_mutated_allele2 = (long)species[19];

//UNPACK ALL RATES HERE
  r_prod_A1 = rates[0]
  r_prod_Anonsense1 = rates[1]
  r_prod_Aprime1 = rates[2]
  r_prod_B1 = rates[3]
  
  d_on_A1 = rates[4]
  d_on_Anonsense1 = rates[5]
  d_on_Aprime1 = rates[6]
  d_on_B1 = rates[7]

  r_deg_A1 = rates[8]
  r_deg_Anonsense1 = rates[9]
  r_deg_Aprime1 = rates[10]
  r_deg_B1 = rates[11]

  r_onbasal_A1 = rates[12]
  r_onbasal_Anonsense1 = rates[13]
  r_onbasal_Aprime1 = rates[14]
  r_onbasal_B1 = rates[15]

  d_bound_byAnonsense1_A1 = rates[16]
  d_bound_byAnonsense1_Anonsense1 = rates[17]
  d_bound_byAnonsense1_Aprime1 = rates[18]
  d_bound_byA1_B1 = rates[19] 
  d_bound_byAprime1_B1 = rates[20] 

  r_off_A1 = rates[21]
  r_off_Anonsense1 = rates[22]
  r_off_Aprime1 = rates[23]
  r_off_B1 = rates[24]

  r_bind_byAnonsense1_A1 = rates[25]
  r_bind_byAnonsense1_Anonsense1 = rates[26]
  r_bind_byAnonsense1_Aprime1 = rates[27]
  r_bind_byA1_B1 = rates[28] 
  r_bind_byAprime1_B1 = rates[29] 

  k_A1  = rates[30];
  k_Anonsense1  = rates[31];
  k_Aprime1  = rates[32];
  k_B1  = rates[33];

  n_A1  = rates[34];
  n_Anonsense1  = rates[35];
  n_Aprime1  = rates[36];
  n_B1  = rates[37];

  r_unbind_byAnonsense1_A1 = rates[38]
  r_unbind_byAnonsense1_Anonsense1 = rates[39]
  r_unbind_byAnonsense1_Aprime1 = rates[40]
  r_unbind_byA1_B1 = rates[41] 
  r_unbind_byAprime1_B1 = rates[42] 

 

  zigset(seed);
  totaliterations = 0;

  savecount = 0;
  savecountcheck = 0;
  deltaTsave = totalt/(m-1);


  //  for (i=0; i<m; i++) {
  while (currT<totalt) {
    if (currT > 50000) {
      Burst1_is_mutated_allele1 = 1;
      Burst1_not_mutated_allele1 = 0;
    }
    if (currT > 100000) {
      Burst1_is_mutated_allele1 = 1;
      Burst1_not_mutated_allele1 = 0;
      Burst1_is_mutated_allele2 = 1;
      Burst1_not_mutated_allele2 = 0;
    }
    totaliterations++;
    cumpropensities[0] = propensities[0];
    for (j=1; j<NUMRXNS; j++)
      cumpropensities[j] = cumpropensities[j-1]+propensities[j];

    alpha = cumpropensities[NUMRXNS-1];
    //deltaT = -1.0/alpha*log(rand1[i]);
    deltaT = -1.0/alpha*log(UNI);
    currT += deltaT;
    //p = rand2[i]*alpha;
    p = UNI*alpha;

//if(savecount < 1000){
    while (currT > savecount*deltaTsave) {
        if (savecount < m){
//INSERT SAVE HERE
species_out[savecount*NSPECIES+0] = A1;
species_out[savecount*NSPECIES+1] = A1_nonsense;
species_out[savecount*NSPECIES+2] = Aprime1;
species_out[savecount*NSPECIES+3] = B1;
species_out[savecount*NSPECIES+4] = Burst1_off_targ_allele1;
species_out[savecount*NSPECIES+5] = Burst1_on_targ_allele1;
species_out[savecount*NSPECIES+6] = Burst1_off_targ_allele2;
species_out[savecount*NSPECIES+7] = Burst1_on_targ_allele2;
species_out[savecount*NSPECIES+8] = Burst1_off_para_allele1;
species_out[savecount*NSPECIES+9] = Burst1_on_para_allele1;
species_out[savecount*NSPECIES+10] = Burst1_off_para_allele2;
species_out[savecount*NSPECIES+11] = Burst1_on_para_allele2;
species_out[savecount*NSPECIES+12] = Burst1_on_orig_allele1;
species_out[savecount*NSPECIES+13] = Burst1_off_orig_allele1;
species_out[savecount*NSPECIES+14] = Burst1_on_orig_allele2;
species_out[savecount*NSPECIES+15] = Burst1_off_orig_allele2;
species_out[savecount*NSPECIES+16] = Burst1_is_mutated_allele1;
species_out[savecount*NSPECIES+17] = Burst1_not_mutated_allele1;
species_out[savecount*NSPECIES+18] = Burst1_is_mutated_allele2;
species_out[savecount*NSPECIES+19] = Burst1_not_mutated_allele2;
            times_out[savecount] = savecount*deltaTsave;
            savecount++;
            savecountcheck++;
            }
		else{savecount++;
		}
    }

//INSERT IF STATEMENT HERE
if (p<cumpropensities[0]) {
  // rxn: = A1 from allele1
  if (Burst1_is_mutated_allele1 == 0){ // extra check for mutation ?necessary?
    A1=A1 + 1;
  }
  //update propensity for = A1 from allele1
  propensities[0] = Burst1_not_mutated_allele1*(A_prod1 *A_proddiff1 *Burst1_on_orig_allele1+A_prod1 *Burst1_off_orig_allele1);
  //update propensity for = A1_nonsense from allele1
  propensities[1] = Burst1_is_mutated_allele1*(Anonsense_prod1 *A_proddiff1 *Burst1_on_orig_allele1+Anonsense_prod1 *Burst1_off_orig_allele1);
  //update propensity for = A1 from allele2
  propensities[2] = Burst1_not_mutated_allele2*(A_prod1 *A_proddiff1 *Burst1_on_orig_allele2+A_prod1 *Burst1_off_orig_allele2);
  //update propensity for = A1_nonsense from allele2
  propensities[3] = Burst1_is_mutated_allele2*(Anonsense_prod1 *A_proddiff1 *Burst1_on_orig_allele2+Anonsense_prod1 *Burst1_off_orig_allele2);
  //update propensity for = Aprime1 from allele1
  propensities[4] = Aprime_prod1 *Aprime_proddiff1 *Burst1_on_para_allele1+Aprime_prod1 *Burst1_off_para_allele1;
  //update propensity for = Aprime1 from allele2
  propensities[5] = Aprime_prod1 *Aprime_proddiff1 *Burst1_on_para_allele2+Aprime_prod1 *Burst1_off_para_allele2;
  //update propensity for = B1 from allele1
  propensities[6] = B_prod1 *B_proddiff1 *Burst1_on_targ_allele1+B_prod1 *Burst1_off_targ_allele1;
  //update propensity for = B1 from allele2
  propensities[7] = B_prod1 *B_proddiff1 *Burst1_on_targ_allele1+B_prod1 *Burst1_off_targ_allele1;
  //update propensity for A1 = 
  propensities[8] = A_deg1 *A1;
  //update propensity for A1_nonsense = 
  propensities[9] = Anonsense_deg1 *A1_nonsense;
  //update propensity for Aprime1 = 
  propensities[10] = Aprime_deg1 *Aprime1;
  //update propensity for B1 = 
  propensities[11] = B_deg1 *B1;
  //update propensity for Burst1_off_orig = Burst1_on_orig for allele1
  propensities[12] = onbasal_a1 *Burst1_off_orig_allele1;
  //update propensity for Burst1_off_orig = Burst1_on_orig for allele2
  propensities[13] = onbasal_a1 *Burst1_off_orig_allele2;
  //update propensity for Burst1_off_para = Burst1_on_para for allele1
  propensities[14] = Aprimenitc1 *(pow(A1_nonsense,nA1_nonsense )/(pow(kA1_nonsense ,nA1_nonsense )+pow(A1_nonsense,nA1_nonsense )))*Burst1_off_para_allele1+onbasal_aprime1 *Burst1_off_para_allele1;
  //update propensity for Burst1_off_para = Burst1_on_para for allele2
  propensities[15] = Aprimenitc1 *(pow(A1_nonsense,nA1_nonsense )/(pow(kA1_nonsense ,nA1_nonsense )+pow(A1_nonsense,nA1_nonsense )))*Burst1_off_para_allele2+onbasal_aprime1 *Burst1_off_para_allele2;
  //update propensity for Burst1_off_targ = Burst1_on_targ for allele1
  propensities[16] = B_ondep1 *(pow(A1,nA1 )/(pow(kA1 ,nA1 )+pow(A1,nA1 )))*Burst1_off_targ_allele1+ B_ondep_prime *(pow(Aprime1,nAprime1 )/(pow(kAprime1 ,nAprime1 )+pow(Aprime1,nAprime1 )))*Burst1_off_targ_allele1 +onbasal_b1 *Burst1_off_targ_allele1;
  //update propensity for Burst1_off_targ = Burst1_on_targ for allele2
  propensities[17] = B_ondep1 *(pow(A1,nA1 )/(pow(kA1 ,nA1 )+pow(A1,nA1 )))*Burst1_off_targ_allele2+ B_ondep_prime *(pow(Aprime1,nAprime1 )/(pow(kAprime1 ,nAprime1 )+pow(Aprime1,nAprime1 )))*Burst1_off_targ_allele2 +onbasal_b1 *Burst1_off_targ_allele2;
  //update propensity for Burst1_on_orig = Burst1_off_orig for allele1
  propensities[18] = A_off1 *Burst1_on_orig_allele1;
  //update propensity for Burst1_on_orig = Burst1_off_orig for allele2
  propensities[19] = A_off1 *Burst1_on_orig_allele2;
  //update propensity for Burst1_on_para = Burst1_off_para for allele1
  propensities[20] = Aprime_off1 *Burst1_on_para_allele1;
  //update propensity for Burst1_on_para = Burst1_off_para for allele2
  propensities[21] = Aprime_off1 *Burst1_on_para_allele2;
  //update propensity for Burst1_on_targ = Burst1_off_targ for allele1
  propensities[22] = B_off1 *Burst1_on_targ_allele1;
  //update propensity for Burst1_on_targ = Burst1_off_targ for allele2
  propensities[23] = B_off1 *Burst1_on_targ_allele2;
} else if (p<cumpropensities[1]) {
  // rxn: = A1_nonsense 
  if (Burst1_is_mutated_allele1 == 1){
    A1_nonsense=A1_nonsense + 1;
  }
  //update propensity for = A1 from allele1
  propensities[0] = Burst1_not_mutated_allele1*(A_prod1 *A_proddiff1 *Burst1_on_orig_allele1+A_prod1 *Burst1_off_orig_allele1);
  //update propensity for = A1_nonsense from allele1
  propensities[1] = Burst1_is_mutated_allele1*(Anonsense_prod1 *A_proddiff1 *Burst1_on_orig_allele1+Anonsense_prod1 *Burst1_off_orig_allele1);
  //update propensity for = A1 from allele2
  propensities[2] = Burst1_not_mutated_allele2*(A_prod1 *A_proddiff1 *Burst1_on_orig_allele2+A_prod1 *Burst1_off_orig_allele2);
  //update propensity for = A1_nonsense from allele2
  propensities[3] = Burst1_is_mutated_allele2*(Anonsense_prod1 *A_proddiff1 *Burst1_on_orig_allele2+Anonsense_prod1 *Burst1_off_orig_allele2);
  //update propensity for = Aprime1 from allele1
  propensities[4] = Aprime_prod1 *Aprime_proddiff1 *Burst1_on_para_allele1+Aprime_prod1 *Burst1_off_para_allele1;
  //update propensity for = Aprime1 from allele2
  propensities[5] = Aprime_prod1 *Aprime_proddiff1 *Burst1_on_para_allele2+Aprime_prod1 *Burst1_off_para_allele2;
  //update propensity for = B1 from allele1
  propensities[6] = B_prod1 *B_proddiff1 *Burst1_on_targ_allele1+B_prod1 *Burst1_off_targ_allele1;
  //update propensity for = B1 from allele2
  propensities[7] = B_prod1 *B_proddiff1 *Burst1_on_targ_allele1+B_prod1 *Burst1_off_targ_allele1;
  //update propensity for A1 = 
  propensities[8] = A_deg1 *A1;
  //update propensity for A1_nonsense = 
  propensities[9] = Anonsense_deg1 *A1_nonsense;
  //update propensity for Aprime1 = 
  propensities[10] = Aprime_deg1 *Aprime1;
  //update propensity for B1 = 
  propensities[11] = B_deg1 *B1;
  //update propensity for Burst1_off_orig = Burst1_on_orig for allele1
  propensities[12] = onbasal_a1 *Burst1_off_orig_allele1;
  //update propensity for Burst1_off_orig = Burst1_on_orig for allele2
  propensities[13] = onbasal_a1 *Burst1_off_orig_allele2;
  //update propensity for Burst1_off_para = Burst1_on_para for allele1
  propensities[14] = Aprimenitc1 *(pow(A1_nonsense,nA1_nonsense )/(pow(kA1_nonsense ,nA1_nonsense )+pow(A1_nonsense,nA1_nonsense )))*Burst1_off_para_allele1+onbasal_aprime1 *Burst1_off_para_allele1;
  //update propensity for Burst1_off_para = Burst1_on_para for allele2
  propensities[15] = Aprimenitc1 *(pow(A1_nonsense,nA1_nonsense )/(pow(kA1_nonsense ,nA1_nonsense )+pow(A1_nonsense,nA1_nonsense )))*Burst1_off_para_allele2+onbasal_aprime1 *Burst1_off_para_allele2;
  //update propensity for Burst1_off_targ = Burst1_on_targ for allele1
  propensities[16] = B_ondep1 *(pow(A1,nA1 )/(pow(kA1 ,nA1 )+pow(A1,nA1 )))*Burst1_off_targ_allele1+ B_ondep_prime *(pow(Aprime1,nAprime1 )/(pow(kAprime1 ,nAprime1 )+pow(Aprime1,nAprime1 )))*Burst1_off_targ_allele1 +onbasal_b1 *Burst1_off_targ_allele1;
  //update propensity for Burst1_off_targ = Burst1_on_targ for allele2
  propensities[17] = B_ondep1 *(pow(A1,nA1 )/(pow(kA1 ,nA1 )+pow(A1,nA1 )))*Burst1_off_targ_allele2+ B_ondep_prime *(pow(Aprime1,nAprime1 )/(pow(kAprime1 ,nAprime1 )+pow(Aprime1,nAprime1 )))*Burst1_off_targ_allele2 +onbasal_b1 *Burst1_off_targ_allele2;
  //update propensity for Burst1_on_orig = Burst1_off_orig for allele1
  propensities[18] = A_off1 *Burst1_on_orig_allele1;
  //update propensity for Burst1_on_orig = Burst1_off_orig for allele2
  propensities[19] = A_off1 *Burst1_on_orig_allele2;
  //update propensity for Burst1_on_para = Burst1_off_para for allele1
  propensities[20] = Aprime_off1 *Burst1_on_para_allele1;
  //update propensity for Burst1_on_para = Burst1_off_para for allele2
  propensities[21] = Aprime_off1 *Burst1_on_para_allele2;
  //update propensity for Burst1_on_targ = Burst1_off_targ for allele1
  propensities[22] = B_off1 *Burst1_on_targ_allele1;
  //update propensity for Burst1_on_targ = Burst1_off_targ for allele2
  propensities[23] = B_off1 *Burst1_on_targ_allele2;
} else if (p<cumpropensities[2]) {
    // rxn: = A1 from allele1
  if (Burst1_is_mutated_allele2 == 0){ // extra check for mutation ?necessary?
    A1=A1 + 1;
  }
  //update propensity for = A1 from allele1
  propensities[0] = Burst1_not_mutated_allele1*(A_prod1 *A_proddiff1 *Burst1_on_orig_allele1+A_prod1 *Burst1_off_orig_allele1);
  //update propensity for = A1_nonsense from allele1
  propensities[1] = Burst1_is_mutated_allele1*(Anonsense_prod1 *A_proddiff1 *Burst1_on_orig_allele1+Anonsense_prod1 *Burst1_off_orig_allele1);
  //update propensity for = A1 from allele2
  propensities[2] = Burst1_not_mutated_allele2*(A_prod1 *A_proddiff1 *Burst1_on_orig_allele2+A_prod1 *Burst1_off_orig_allele2);
  //update propensity for = A1_nonsense from allele2
  propensities[3] = Burst1_is_mutated_allele2*(Anonsense_prod1 *A_proddiff1 *Burst1_on_orig_allele2+Anonsense_prod1 *Burst1_off_orig_allele2);
  //update propensity for = Aprime1 from allele1
  propensities[4] = Aprime_prod1 *Aprime_proddiff1 *Burst1_on_para_allele1+Aprime_prod1 *Burst1_off_para_allele1;
  //update propensity for = Aprime1 from allele2
  propensities[5] = Aprime_prod1 *Aprime_proddiff1 *Burst1_on_para_allele2+Aprime_prod1 *Burst1_off_para_allele2;
  //update propensity for = B1 from allele1
  propensities[6] = B_prod1 *B_proddiff1 *Burst1_on_targ_allele1+B_prod1 *Burst1_off_targ_allele1;
  //update propensity for = B1 from allele2
  propensities[7] = B_prod1 *B_proddiff1 *Burst1_on_targ_allele1+B_prod1 *Burst1_off_targ_allele1;
  //update propensity for A1 = 
  propensities[8] = A_deg1 *A1;
  //update propensity for A1_nonsense = 
  propensities[9] = Anonsense_deg1 *A1_nonsense;
  //update propensity for Aprime1 = 
  propensities[10] = Aprime_deg1 *Aprime1;
  //update propensity for B1 = 
  propensities[11] = B_deg1 *B1;
  //update propensity for Burst1_off_orig = Burst1_on_orig for allele1
  propensities[12] = onbasal_a1 *Burst1_off_orig_allele1;
  //update propensity for Burst1_off_orig = Burst1_on_orig for allele2
  propensities[13] = onbasal_a1 *Burst1_off_orig_allele2;
  //update propensity for Burst1_off_para = Burst1_on_para for allele1
  propensities[14] = Aprimenitc1 *(pow(A1_nonsense,nA1_nonsense )/(pow(kA1_nonsense ,nA1_nonsense )+pow(A1_nonsense,nA1_nonsense )))*Burst1_off_para_allele1+onbasal_aprime1 *Burst1_off_para_allele1;
  //update propensity for Burst1_off_para = Burst1_on_para for allele2
  propensities[15] = Aprimenitc1 *(pow(A1_nonsense,nA1_nonsense )/(pow(kA1_nonsense ,nA1_nonsense )+pow(A1_nonsense,nA1_nonsense )))*Burst1_off_para_allele2+onbasal_aprime1 *Burst1_off_para_allele2;
  //update propensity for Burst1_off_targ = Burst1_on_targ for allele1
  propensities[16] = B_ondep1 *(pow(A1,nA1 )/(pow(kA1 ,nA1 )+pow(A1,nA1 )))*Burst1_off_targ_allele1+ B_ondep_prime *(pow(Aprime1,nAprime1 )/(pow(kAprime1 ,nAprime1 )+pow(Aprime1,nAprime1 )))*Burst1_off_targ_allele1 +onbasal_b1 *Burst1_off_targ_allele1;
  //update propensity for Burst1_off_targ = Burst1_on_targ for allele2
  propensities[17] = B_ondep1 *(pow(A1,nA1 )/(pow(kA1 ,nA1 )+pow(A1,nA1 )))*Burst1_off_targ_allele2+ B_ondep_prime *(pow(Aprime1,nAprime1 )/(pow(kAprime1 ,nAprime1 )+pow(Aprime1,nAprime1 )))*Burst1_off_targ_allele2 +onbasal_b1 *Burst1_off_targ_allele2;
  //update propensity for Burst1_on_orig = Burst1_off_orig for allele1
  propensities[18] = A_off1 *Burst1_on_orig_allele1;
  //update propensity for Burst1_on_orig = Burst1_off_orig for allele2
  propensities[19] = A_off1 *Burst1_on_orig_allele2;
  //update propensity for Burst1_on_para = Burst1_off_para for allele1
  propensities[20] = Aprime_off1 *Burst1_on_para_allele1;
  //update propensity for Burst1_on_para = Burst1_off_para for allele2
  propensities[21] = Aprime_off1 *Burst1_on_para_allele2;
  //update propensity for Burst1_on_targ = Burst1_off_targ for allele1
  propensities[22] = B_off1 *Burst1_on_targ_allele1;
  //update propensity for Burst1_on_targ = Burst1_off_targ for allele2
  propensities[23] = B_off1 *Burst1_on_targ_allele2;
} else if (p<cumpropensities[3]) {
  // rxn: = A1_nonsense 
  if (Burst1_is_mutated_allele2 == 1){
    A1_nonsense=A1_nonsense + 1;
  }
  //update propensity for = A1 from allele1
  propensities[0] = Burst1_not_mutated_allele1*(A_prod1 *A_proddiff1 *Burst1_on_orig_allele1+A_prod1 *Burst1_off_orig_allele1);
  //update propensity for = A1_nonsense from allele1
  propensities[1] = Burst1_is_mutated_allele1*(Anonsense_prod1 *A_proddiff1 *Burst1_on_orig_allele1+Anonsense_prod1 *Burst1_off_orig_allele1);
  //update propensity for = A1 from allele2
  propensities[2] = Burst1_not_mutated_allele2*(A_prod1 *A_proddiff1 *Burst1_on_orig_allele2+A_prod1 *Burst1_off_orig_allele2);
  //update propensity for = A1_nonsense from allele2
  propensities[3] = Burst1_is_mutated_allele2*(Anonsense_prod1 *A_proddiff1 *Burst1_on_orig_allele2+Anonsense_prod1 *Burst1_off_orig_allele2);
  //update propensity for = Aprime1 from allele1
  propensities[4] = Aprime_prod1 *Aprime_proddiff1 *Burst1_on_para_allele1+Aprime_prod1 *Burst1_off_para_allele1;
  //update propensity for = Aprime1 from allele2
  propensities[5] = Aprime_prod1 *Aprime_proddiff1 *Burst1_on_para_allele2+Aprime_prod1 *Burst1_off_para_allele2;
  //update propensity for = B1 from allele1
  propensities[6] = B_prod1 *B_proddiff1 *Burst1_on_targ_allele1+B_prod1 *Burst1_off_targ_allele1;
  //update propensity for = B1 from allele2
  propensities[7] = B_prod1 *B_proddiff1 *Burst1_on_targ_allele1+B_prod1 *Burst1_off_targ_allele1;
  //update propensity for A1 = 
  propensities[8] = A_deg1 *A1;
  //update propensity for A1_nonsense = 
  propensities[9] = Anonsense_deg1 *A1_nonsense;
  //update propensity for Aprime1 = 
  propensities[10] = Aprime_deg1 *Aprime1;
  //update propensity for B1 = 
  propensities[11] = B_deg1 *B1;
  //update propensity for Burst1_off_orig = Burst1_on_orig for allele1
  propensities[12] = onbasal_a1 *Burst1_off_orig_allele1;
  //update propensity for Burst1_off_orig = Burst1_on_orig for allele2
  propensities[13] = onbasal_a1 *Burst1_off_orig_allele2;
  //update propensity for Burst1_off_para = Burst1_on_para for allele1
  propensities[14] = Aprimenitc1 *(pow(A1_nonsense,nA1_nonsense )/(pow(kA1_nonsense ,nA1_nonsense )+pow(A1_nonsense,nA1_nonsense )))*Burst1_off_para_allele1+onbasal_aprime1 *Burst1_off_para_allele1;
  //update propensity for Burst1_off_para = Burst1_on_para for allele2
  propensities[15] = Aprimenitc1 *(pow(A1_nonsense,nA1_nonsense )/(pow(kA1_nonsense ,nA1_nonsense )+pow(A1_nonsense,nA1_nonsense )))*Burst1_off_para_allele2+onbasal_aprime1 *Burst1_off_para_allele2;
  //update propensity for Burst1_off_targ = Burst1_on_targ for allele1
  propensities[16] = B_ondep1 *(pow(A1,nA1 )/(pow(kA1 ,nA1 )+pow(A1,nA1 )))*Burst1_off_targ_allele1+ B_ondep_prime *(pow(Aprime1,nAprime1 )/(pow(kAprime1 ,nAprime1 )+pow(Aprime1,nAprime1 )))*Burst1_off_targ_allele1 +onbasal_b1 *Burst1_off_targ_allele1;
  //update propensity for Burst1_off_targ = Burst1_on_targ for allele2
  propensities[17] = B_ondep1 *(pow(A1,nA1 )/(pow(kA1 ,nA1 )+pow(A1,nA1 )))*Burst1_off_targ_allele2+ B_ondep_prime *(pow(Aprime1,nAprime1 )/(pow(kAprime1 ,nAprime1 )+pow(Aprime1,nAprime1 )))*Burst1_off_targ_allele2 +onbasal_b1 *Burst1_off_targ_allele2;
  //update propensity for Burst1_on_orig = Burst1_off_orig for allele1
  propensities[18] = A_off1 *Burst1_on_orig_allele1;
  //update propensity for Burst1_on_orig = Burst1_off_orig for allele2
  propensities[19] = A_off1 *Burst1_on_orig_allele2;
  //update propensity for Burst1_on_para = Burst1_off_para for allele1
  propensities[20] = Aprime_off1 *Burst1_on_para_allele1;
  //update propensity for Burst1_on_para = Burst1_off_para for allele2
  propensities[21] = Aprime_off1 *Burst1_on_para_allele2;
  //update propensity for Burst1_on_targ = Burst1_off_targ for allele1
  propensities[22] = B_off1 *Burst1_on_targ_allele1;
  //update propensity for Burst1_on_targ = Burst1_off_targ for allele2
  propensities[23] = B_off1 *Burst1_on_targ_allele2;
} else if (p<cumpropensities[4]) {
  // rxn: = Aprime1 from allele1
  Aprime1=Aprime1 + 1;

  //update propensity for = A1 from allele1
  propensities[0] = Burst1_not_mutated_allele1*(A_prod1 *A_proddiff1 *Burst1_on_orig_allele1+A_prod1 *Burst1_off_orig_allele1);
  //update propensity for = A1_nonsense from allele1
  propensities[1] = Burst1_is_mutated_allele1*(Anonsense_prod1 *A_proddiff1 *Burst1_on_orig_allele1+Anonsense_prod1 *Burst1_off_orig_allele1);
  //update propensity for = A1 from allele2
  propensities[2] = Burst1_not_mutated_allele2*(A_prod1 *A_proddiff1 *Burst1_on_orig_allele2+A_prod1 *Burst1_off_orig_allele2);
  //update propensity for = A1_nonsense from allele2
  propensities[3] = Burst1_is_mutated_allele2*(Anonsense_prod1 *A_proddiff1 *Burst1_on_orig_allele2+Anonsense_prod1 *Burst1_off_orig_allele2);
  //update propensity for = Aprime1 from allele1
  propensities[4] = Aprime_prod1 *Aprime_proddiff1 *Burst1_on_para_allele1+Aprime_prod1 *Burst1_off_para_allele1;
  //update propensity for = Aprime1 from allele2
  propensities[5] = Aprime_prod1 *Aprime_proddiff1 *Burst1_on_para_allele2+Aprime_prod1 *Burst1_off_para_allele2;
  //update propensity for = B1 from allele1
  propensities[6] = B_prod1 *B_proddiff1 *Burst1_on_targ_allele1+B_prod1 *Burst1_off_targ_allele1;
  //update propensity for = B1 from allele2
  propensities[7] = B_prod1 *B_proddiff1 *Burst1_on_targ_allele1+B_prod1 *Burst1_off_targ_allele1;
  //update propensity for A1 = 
  propensities[8] = A_deg1 *A1;
  //update propensity for A1_nonsense = 
  propensities[9] = Anonsense_deg1 *A1_nonsense;
  //update propensity for Aprime1 = 
  propensities[10] = Aprime_deg1 *Aprime1;
  //update propensity for B1 = 
  propensities[11] = B_deg1 *B1;
  //update propensity for Burst1_off_orig = Burst1_on_orig for allele1
  propensities[12] = onbasal_a1 *Burst1_off_orig_allele1;
  //update propensity for Burst1_off_orig = Burst1_on_orig for allele2
  propensities[13] = onbasal_a1 *Burst1_off_orig_allele2;
  //update propensity for Burst1_off_para = Burst1_on_para for allele1
  propensities[14] = Aprimenitc1 *(pow(A1_nonsense,nA1_nonsense )/(pow(kA1_nonsense ,nA1_nonsense )+pow(A1_nonsense,nA1_nonsense )))*Burst1_off_para_allele1+onbasal_aprime1 *Burst1_off_para_allele1;
  //update propensity for Burst1_off_para = Burst1_on_para for allele2
  propensities[15] = Aprimenitc1 *(pow(A1_nonsense,nA1_nonsense )/(pow(kA1_nonsense ,nA1_nonsense )+pow(A1_nonsense,nA1_nonsense )))*Burst1_off_para_allele2+onbasal_aprime1 *Burst1_off_para_allele2;
  //update propensity for Burst1_off_targ = Burst1_on_targ for allele1
  propensities[16] = B_ondep1 *(pow(A1,nA1 )/(pow(kA1 ,nA1 )+pow(A1,nA1 )))*Burst1_off_targ_allele1+ B_ondep_prime *(pow(Aprime1,nAprime1 )/(pow(kAprime1 ,nAprime1 )+pow(Aprime1,nAprime1 )))*Burst1_off_targ_allele1 +onbasal_b1 *Burst1_off_targ_allele1;
  //update propensity for Burst1_off_targ = Burst1_on_targ for allele2
  propensities[17] = B_ondep1 *(pow(A1,nA1 )/(pow(kA1 ,nA1 )+pow(A1,nA1 )))*Burst1_off_targ_allele2+ B_ondep_prime *(pow(Aprime1,nAprime1 )/(pow(kAprime1 ,nAprime1 )+pow(Aprime1,nAprime1 )))*Burst1_off_targ_allele2 +onbasal_b1 *Burst1_off_targ_allele2;
  //update propensity for Burst1_on_orig = Burst1_off_orig for allele1
  propensities[18] = A_off1 *Burst1_on_orig_allele1;
  //update propensity for Burst1_on_orig = Burst1_off_orig for allele2
  propensities[19] = A_off1 *Burst1_on_orig_allele2;
  //update propensity for Burst1_on_para = Burst1_off_para for allele1
  propensities[20] = Aprime_off1 *Burst1_on_para_allele1;
  //update propensity for Burst1_on_para = Burst1_off_para for allele2
  propensities[21] = Aprime_off1 *Burst1_on_para_allele2;
  //update propensity for Burst1_on_targ = Burst1_off_targ for allele1
  propensities[22] = B_off1 *Burst1_on_targ_allele1;
  //update propensity for Burst1_on_targ = Burst1_off_targ for allele2
  propensities[23] = B_off1 *Burst1_on_targ_allele2;
} else if (p<cumpropensities[5]) {
  // rxn: = Aprime1 from allele2
  Aprime1=Aprime1 + 1;

  //update propensity for = A1 from allele1
  propensities[0] = Burst1_not_mutated_allele1*(A_prod1 *A_proddiff1 *Burst1_on_orig_allele1+A_prod1 *Burst1_off_orig_allele1);
  //update propensity for = A1_nonsense from allele1
  propensities[1] = Burst1_is_mutated_allele1*(Anonsense_prod1 *A_proddiff1 *Burst1_on_orig_allele1+Anonsense_prod1 *Burst1_off_orig_allele1);
  //update propensity for = A1 from allele2
  propensities[2] = Burst1_not_mutated_allele2*(A_prod1 *A_proddiff1 *Burst1_on_orig_allele2+A_prod1 *Burst1_off_orig_allele2);
  //update propensity for = A1_nonsense from allele2
  propensities[3] = Burst1_is_mutated_allele2*(Anonsense_prod1 *A_proddiff1 *Burst1_on_orig_allele2+Anonsense_prod1 *Burst1_off_orig_allele2);
  //update propensity for = Aprime1 from allele1
  propensities[4] = Aprime_prod1 *Aprime_proddiff1 *Burst1_on_para_allele1+Aprime_prod1 *Burst1_off_para_allele1;
  //update propensity for = Aprime1 from allele2
  propensities[5] = Aprime_prod1 *Aprime_proddiff1 *Burst1_on_para_allele2+Aprime_prod1 *Burst1_off_para_allele2;
  //update propensity for = B1 from allele1
  propensities[6] = B_prod1 *B_proddiff1 *Burst1_on_targ_allele1+B_prod1 *Burst1_off_targ_allele1;
  //update propensity for = B1 from allele2
  propensities[7] = B_prod1 *B_proddiff1 *Burst1_on_targ_allele1+B_prod1 *Burst1_off_targ_allele1;
  //update propensity for A1 = 
  propensities[8] = A_deg1 *A1;
  //update propensity for A1_nonsense = 
  propensities[9] = Anonsense_deg1 *A1_nonsense;
  //update propensity for Aprime1 = 
  propensities[10] = Aprime_deg1 *Aprime1;
  //update propensity for B1 = 
  propensities[11] = B_deg1 *B1;
  //update propensity for Burst1_off_orig = Burst1_on_orig for allele1
  propensities[12] = onbasal_a1 *Burst1_off_orig_allele1;
  //update propensity for Burst1_off_orig = Burst1_on_orig for allele2
  propensities[13] = onbasal_a1 *Burst1_off_orig_allele2;
  //update propensity for Burst1_off_para = Burst1_on_para for allele1
  propensities[14] = Aprimenitc1 *(pow(A1_nonsense,nA1_nonsense )/(pow(kA1_nonsense ,nA1_nonsense )+pow(A1_nonsense,nA1_nonsense )))*Burst1_off_para_allele1+onbasal_aprime1 *Burst1_off_para_allele1;
  //update propensity for Burst1_off_para = Burst1_on_para for allele2
  propensities[15] = Aprimenitc1 *(pow(A1_nonsense,nA1_nonsense )/(pow(kA1_nonsense ,nA1_nonsense )+pow(A1_nonsense,nA1_nonsense )))*Burst1_off_para_allele2+onbasal_aprime1 *Burst1_off_para_allele2;
  //update propensity for Burst1_off_targ = Burst1_on_targ for allele1
  propensities[16] = B_ondep1 *(pow(A1,nA1 )/(pow(kA1 ,nA1 )+pow(A1,nA1 )))*Burst1_off_targ_allele1+ B_ondep_prime *(pow(Aprime1,nAprime1 )/(pow(kAprime1 ,nAprime1 )+pow(Aprime1,nAprime1 )))*Burst1_off_targ_allele1 +onbasal_b1 *Burst1_off_targ_allele1;
  //update propensity for Burst1_off_targ = Burst1_on_targ for allele2
  propensities[17] = B_ondep1 *(pow(A1,nA1 )/(pow(kA1 ,nA1 )+pow(A1,nA1 )))*Burst1_off_targ_allele2+ B_ondep_prime *(pow(Aprime1,nAprime1 )/(pow(kAprime1 ,nAprime1 )+pow(Aprime1,nAprime1 )))*Burst1_off_targ_allele2 +onbasal_b1 *Burst1_off_targ_allele2;
  //update propensity for Burst1_on_orig = Burst1_off_orig for allele1
  propensities[18] = A_off1 *Burst1_on_orig_allele1;
  //update propensity for Burst1_on_orig = Burst1_off_orig for allele2
  propensities[19] = A_off1 *Burst1_on_orig_allele2;
  //update propensity for Burst1_on_para = Burst1_off_para for allele1
  propensities[20] = Aprime_off1 *Burst1_on_para_allele1;
  //update propensity for Burst1_on_para = Burst1_off_para for allele2
  propensities[21] = Aprime_off1 *Burst1_on_para_allele2;
  //update propensity for Burst1_on_targ = Burst1_off_targ for allele1
  propensities[22] = B_off1 *Burst1_on_targ_allele1;
  //update propensity for Burst1_on_targ = Burst1_off_targ for allele2
  propensities[23] = B_off1 *Burst1_on_targ_allele2;
} else if (p<cumpropensities[6]) {
  // rxn: = B1 from allele1
  B1=B1 + 1;

   //update propensity for = A1 from allele1
  propensities[0] = Burst1_not_mutated_allele1*(A_prod1 *A_proddiff1 *Burst1_on_orig_allele1+A_prod1 *Burst1_off_orig_allele1);
  //update propensity for = A1_nonsense from allele1
  propensities[1] = Burst1_is_mutated_allele1*(Anonsense_prod1 *A_proddiff1 *Burst1_on_orig_allele1+Anonsense_prod1 *Burst1_off_orig_allele1);
  //update propensity for = A1 from allele2
  propensities[2] = Burst1_not_mutated_allele2*(A_prod1 *A_proddiff1 *Burst1_on_orig_allele2+A_prod1 *Burst1_off_orig_allele2);
  //update propensity for = A1_nonsense from allele2
  propensities[3] = Burst1_is_mutated_allele2*(Anonsense_prod1 *A_proddiff1 *Burst1_on_orig_allele2+Anonsense_prod1 *Burst1_off_orig_allele2);
  //update propensity for = Aprime1 from allele1
  propensities[4] = Aprime_prod1 *Aprime_proddiff1 *Burst1_on_para_allele1+Aprime_prod1 *Burst1_off_para_allele1;
  //update propensity for = Aprime1 from allele2
  propensities[5] = Aprime_prod1 *Aprime_proddiff1 *Burst1_on_para_allele2+Aprime_prod1 *Burst1_off_para_allele2;
  //update propensity for = B1 from allele1
  propensities[6] = B_prod1 *B_proddiff1 *Burst1_on_targ_allele1+B_prod1 *Burst1_off_targ_allele1;
  //update propensity for = B1 from allele2
  propensities[7] = B_prod1 *B_proddiff1 *Burst1_on_targ_allele1+B_prod1 *Burst1_off_targ_allele1;
  //update propensity for A1 = 
  propensities[8] = A_deg1 *A1;
  //update propensity for A1_nonsense = 
  propensities[9] = Anonsense_deg1 *A1_nonsense;
  //update propensity for Aprime1 = 
  propensities[10] = Aprime_deg1 *Aprime1;
  //update propensity for B1 = 
  propensities[11] = B_deg1 *B1;
  //update propensity for Burst1_off_orig = Burst1_on_orig for allele1
  propensities[12] = onbasal_a1 *Burst1_off_orig_allele1;
  //update propensity for Burst1_off_orig = Burst1_on_orig for allele2
  propensities[13] = onbasal_a1 *Burst1_off_orig_allele2;
  //update propensity for Burst1_off_para = Burst1_on_para for allele1
  propensities[14] = Aprimenitc1 *(pow(A1_nonsense,nA1_nonsense )/(pow(kA1_nonsense ,nA1_nonsense )+pow(A1_nonsense,nA1_nonsense )))*Burst1_off_para_allele1+onbasal_aprime1 *Burst1_off_para_allele1;
  //update propensity for Burst1_off_para = Burst1_on_para for allele2
  propensities[15] = Aprimenitc1 *(pow(A1_nonsense,nA1_nonsense )/(pow(kA1_nonsense ,nA1_nonsense )+pow(A1_nonsense,nA1_nonsense )))*Burst1_off_para_allele2+onbasal_aprime1 *Burst1_off_para_allele2;
  //update propensity for Burst1_off_targ = Burst1_on_targ for allele1
  propensities[16] = B_ondep1 *(pow(A1,nA1 )/(pow(kA1 ,nA1 )+pow(A1,nA1 )))*Burst1_off_targ_allele1+ B_ondep_prime *(pow(Aprime1,nAprime1 )/(pow(kAprime1 ,nAprime1 )+pow(Aprime1,nAprime1 )))*Burst1_off_targ_allele1 +onbasal_b1 *Burst1_off_targ_allele1;
  //update propensity for Burst1_off_targ = Burst1_on_targ for allele2
  propensities[17] = B_ondep1 *(pow(A1,nA1 )/(pow(kA1 ,nA1 )+pow(A1,nA1 )))*Burst1_off_targ_allele2+ B_ondep_prime *(pow(Aprime1,nAprime1 )/(pow(kAprime1 ,nAprime1 )+pow(Aprime1,nAprime1 )))*Burst1_off_targ_allele2 +onbasal_b1 *Burst1_off_targ_allele2;
  //update propensity for Burst1_on_orig = Burst1_off_orig for allele1
  propensities[18] = A_off1 *Burst1_on_orig_allele1;
  //update propensity for Burst1_on_orig = Burst1_off_orig for allele2
  propensities[19] = A_off1 *Burst1_on_orig_allele2;
  //update propensity for Burst1_on_para = Burst1_off_para for allele1
  propensities[20] = Aprime_off1 *Burst1_on_para_allele1;
  //update propensity for Burst1_on_para = Burst1_off_para for allele2
  propensities[21] = Aprime_off1 *Burst1_on_para_allele2;
  //update propensity for Burst1_on_targ = Burst1_off_targ for allele1
  propensities[22] = B_off1 *Burst1_on_targ_allele1;
  //update propensity for Burst1_on_targ = Burst1_off_targ for allele2
  propensities[23] = B_off1 *Burst1_on_targ_allele2;
} else if (p<cumpropensities[7]) {
  // rxn: = B1 from allele2
  B1=B1 + 1;

   //update propensity for = A1 from allele1
  propensities[0] = Burst1_not_mutated_allele1*(A_prod1 *A_proddiff1 *Burst1_on_orig_allele1+A_prod1 *Burst1_off_orig_allele1);
  //update propensity for = A1_nonsense from allele1
  propensities[1] = Burst1_is_mutated_allele1*(Anonsense_prod1 *A_proddiff1 *Burst1_on_orig_allele1+Anonsense_prod1 *Burst1_off_orig_allele1);
  //update propensity for = A1 from allele2
  propensities[2] = Burst1_not_mutated_allele2*(A_prod1 *A_proddiff1 *Burst1_on_orig_allele2+A_prod1 *Burst1_off_orig_allele2);
  //update propensity for = A1_nonsense from allele2
  propensities[3] = Burst1_is_mutated_allele2*(Anonsense_prod1 *A_proddiff1 *Burst1_on_orig_allele2+Anonsense_prod1 *Burst1_off_orig_allele2);
  //update propensity for = Aprime1 from allele1
  propensities[4] = Aprime_prod1 *Aprime_proddiff1 *Burst1_on_para_allele1+Aprime_prod1 *Burst1_off_para_allele1;
  //update propensity for = Aprime1 from allele2
  propensities[5] = Aprime_prod1 *Aprime_proddiff1 *Burst1_on_para_allele2+Aprime_prod1 *Burst1_off_para_allele2;
  //update propensity for = B1 from allele1
  propensities[6] = B_prod1 *B_proddiff1 *Burst1_on_targ_allele1+B_prod1 *Burst1_off_targ_allele1;
  //update propensity for = B1 from allele2
  propensities[7] = B_prod1 *B_proddiff1 *Burst1_on_targ_allele1+B_prod1 *Burst1_off_targ_allele1;
  //update propensity for A1 = 
  propensities[8] = A_deg1 *A1;
  //update propensity for A1_nonsense = 
  propensities[9] = Anonsense_deg1 *A1_nonsense;
  //update propensity for Aprime1 = 
  propensities[10] = Aprime_deg1 *Aprime1;
  //update propensity for B1 = 
  propensities[11] = B_deg1 *B1;
  //update propensity for Burst1_off_orig = Burst1_on_orig for allele1
  propensities[12] = onbasal_a1 *Burst1_off_orig_allele1;
  //update propensity for Burst1_off_orig = Burst1_on_orig for allele2
  propensities[13] = onbasal_a1 *Burst1_off_orig_allele2;
  //update propensity for Burst1_off_para = Burst1_on_para for allele1
  propensities[14] = Aprimenitc1 *(pow(A1_nonsense,nA1_nonsense )/(pow(kA1_nonsense ,nA1_nonsense )+pow(A1_nonsense,nA1_nonsense )))*Burst1_off_para_allele1+onbasal_aprime1 *Burst1_off_para_allele1;
  //update propensity for Burst1_off_para = Burst1_on_para for allele2
  propensities[15] = Aprimenitc1 *(pow(A1_nonsense,nA1_nonsense )/(pow(kA1_nonsense ,nA1_nonsense )+pow(A1_nonsense,nA1_nonsense )))*Burst1_off_para_allele2+onbasal_aprime1 *Burst1_off_para_allele2;
  //update propensity for Burst1_off_targ = Burst1_on_targ for allele1
  propensities[16] = B_ondep1 *(pow(A1,nA1 )/(pow(kA1 ,nA1 )+pow(A1,nA1 )))*Burst1_off_targ_allele1+ B_ondep_prime *(pow(Aprime1,nAprime1 )/(pow(kAprime1 ,nAprime1 )+pow(Aprime1,nAprime1 )))*Burst1_off_targ_allele1 +onbasal_b1 *Burst1_off_targ_allele1;
  //update propensity for Burst1_off_targ = Burst1_on_targ for allele2
  propensities[17] = B_ondep1 *(pow(A1,nA1 )/(pow(kA1 ,nA1 )+pow(A1,nA1 )))*Burst1_off_targ_allele2+ B_ondep_prime *(pow(Aprime1,nAprime1 )/(pow(kAprime1 ,nAprime1 )+pow(Aprime1,nAprime1 )))*Burst1_off_targ_allele2 +onbasal_b1 *Burst1_off_targ_allele2;
  //update propensity for Burst1_on_orig = Burst1_off_orig for allele1
  propensities[18] = A_off1 *Burst1_on_orig_allele1;
  //update propensity for Burst1_on_orig = Burst1_off_orig for allele2
  propensities[19] = A_off1 *Burst1_on_orig_allele2;
  //update propensity for Burst1_on_para = Burst1_off_para for allele1
  propensities[20] = Aprime_off1 *Burst1_on_para_allele1;
  //update propensity for Burst1_on_para = Burst1_off_para for allele2
  propensities[21] = Aprime_off1 *Burst1_on_para_allele2;
  //update propensity for Burst1_on_targ = Burst1_off_targ for allele1
  propensities[22] = B_off1 *Burst1_on_targ_allele1;
  //update propensity for Burst1_on_targ = Burst1_off_targ for allele2
  propensities[23] = B_off1 *Burst1_on_targ_allele2;
} else if (p<cumpropensities[8]) {
  // rxn: A1 = 
  A1=A1 + -1;

   //update propensity for = A1 from allele1
  propensities[0] = Burst1_not_mutated_allele1*(A_prod1 *A_proddiff1 *Burst1_on_orig_allele1+A_prod1 *Burst1_off_orig_allele1);
  //update propensity for = A1_nonsense from allele1
  propensities[1] = Burst1_is_mutated_allele1*(Anonsense_prod1 *A_proddiff1 *Burst1_on_orig_allele1+Anonsense_prod1 *Burst1_off_orig_allele1);
  //update propensity for = A1 from allele2
  propensities[2] = Burst1_not_mutated_allele2*(A_prod1 *A_proddiff1 *Burst1_on_orig_allele2+A_prod1 *Burst1_off_orig_allele2);
  //update propensity for = A1_nonsense from allele2
  propensities[3] = Burst1_is_mutated_allele2*(Anonsense_prod1 *A_proddiff1 *Burst1_on_orig_allele2+Anonsense_prod1 *Burst1_off_orig_allele2);
  //update propensity for = Aprime1 from allele1
  propensities[4] = Aprime_prod1 *Aprime_proddiff1 *Burst1_on_para_allele1+Aprime_prod1 *Burst1_off_para_allele1;
  //update propensity for = Aprime1 from allele2
  propensities[5] = Aprime_prod1 *Aprime_proddiff1 *Burst1_on_para_allele2+Aprime_prod1 *Burst1_off_para_allele2;
  //update propensity for = B1 from allele1
  propensities[6] = B_prod1 *B_proddiff1 *Burst1_on_targ_allele1+B_prod1 *Burst1_off_targ_allele1;
  //update propensity for = B1 from allele2
  propensities[7] = B_prod1 *B_proddiff1 *Burst1_on_targ_allele1+B_prod1 *Burst1_off_targ_allele1;
  //update propensity for A1 = 
  propensities[8] = A_deg1 *A1;
  //update propensity for A1_nonsense = 
  propensities[9] = Anonsense_deg1 *A1_nonsense;
  //update propensity for Aprime1 = 
  propensities[10] = Aprime_deg1 *Aprime1;
  //update propensity for B1 = 
  propensities[11] = B_deg1 *B1;
  //update propensity for Burst1_off_orig = Burst1_on_orig for allele1
  propensities[12] = onbasal_a1 *Burst1_off_orig_allele1;
  //update propensity for Burst1_off_orig = Burst1_on_orig for allele2
  propensities[13] = onbasal_a1 *Burst1_off_orig_allele2;
  //update propensity for Burst1_off_para = Burst1_on_para for allele1
  propensities[14] = Aprimenitc1 *(pow(A1_nonsense,nA1_nonsense )/(pow(kA1_nonsense ,nA1_nonsense )+pow(A1_nonsense,nA1_nonsense )))*Burst1_off_para_allele1+onbasal_aprime1 *Burst1_off_para_allele1;
  //update propensity for Burst1_off_para = Burst1_on_para for allele2
  propensities[15] = Aprimenitc1 *(pow(A1_nonsense,nA1_nonsense )/(pow(kA1_nonsense ,nA1_nonsense )+pow(A1_nonsense,nA1_nonsense )))*Burst1_off_para_allele2+onbasal_aprime1 *Burst1_off_para_allele2;
  //update propensity for Burst1_off_targ = Burst1_on_targ for allele1
  propensities[16] = B_ondep1 *(pow(A1,nA1 )/(pow(kA1 ,nA1 )+pow(A1,nA1 )))*Burst1_off_targ_allele1+ B_ondep_prime *(pow(Aprime1,nAprime1 )/(pow(kAprime1 ,nAprime1 )+pow(Aprime1,nAprime1 )))*Burst1_off_targ_allele1 +onbasal_b1 *Burst1_off_targ_allele1;
  //update propensity for Burst1_off_targ = Burst1_on_targ for allele2
  propensities[17] = B_ondep1 *(pow(A1,nA1 )/(pow(kA1 ,nA1 )+pow(A1,nA1 )))*Burst1_off_targ_allele2+ B_ondep_prime *(pow(Aprime1,nAprime1 )/(pow(kAprime1 ,nAprime1 )+pow(Aprime1,nAprime1 )))*Burst1_off_targ_allele2 +onbasal_b1 *Burst1_off_targ_allele2;
  //update propensity for Burst1_on_orig = Burst1_off_orig for allele1
  propensities[18] = A_off1 *Burst1_on_orig_allele1;
  //update propensity for Burst1_on_orig = Burst1_off_orig for allele2
  propensities[19] = A_off1 *Burst1_on_orig_allele2;
  //update propensity for Burst1_on_para = Burst1_off_para for allele1
  propensities[20] = Aprime_off1 *Burst1_on_para_allele1;
  //update propensity for Burst1_on_para = Burst1_off_para for allele2
  propensities[21] = Aprime_off1 *Burst1_on_para_allele2;
  //update propensity for Burst1_on_targ = Burst1_off_targ for allele1
  propensities[22] = B_off1 *Burst1_on_targ_allele1;
  //update propensity for Burst1_on_targ = Burst1_off_targ for allele2
  propensities[23] = B_off1 *Burst1_on_targ_allele2;
} else if (p<cumpropensities[9]) {
  // rxn: Aprime1 = 
  A1_nonsense=A1_nonsense + -1;

   //update propensity for = A1 from allele1
  propensities[0] = Burst1_not_mutated_allele1*(A_prod1 *A_proddiff1 *Burst1_on_orig_allele1+A_prod1 *Burst1_off_orig_allele1);
  //update propensity for = A1_nonsense from allele1
  propensities[1] = Burst1_is_mutated_allele1*(Anonsense_prod1 *A_proddiff1 *Burst1_on_orig_allele1+Anonsense_prod1 *Burst1_off_orig_allele1);
  //update propensity for = A1 from allele2
  propensities[2] = Burst1_not_mutated_allele2*(A_prod1 *A_proddiff1 *Burst1_on_orig_allele2+A_prod1 *Burst1_off_orig_allele2);
  //update propensity for = A1_nonsense from allele2
  propensities[3] = Burst1_is_mutated_allele2*(Anonsense_prod1 *A_proddiff1 *Burst1_on_orig_allele2+Anonsense_prod1 *Burst1_off_orig_allele2);
  //update propensity for = Aprime1 from allele1
  propensities[4] = Aprime_prod1 *Aprime_proddiff1 *Burst1_on_para_allele1+Aprime_prod1 *Burst1_off_para_allele1;
  //update propensity for = Aprime1 from allele2
  propensities[5] = Aprime_prod1 *Aprime_proddiff1 *Burst1_on_para_allele2+Aprime_prod1 *Burst1_off_para_allele2;
  //update propensity for = B1 from allele1
  propensities[6] = B_prod1 *B_proddiff1 *Burst1_on_targ_allele1+B_prod1 *Burst1_off_targ_allele1;
  //update propensity for = B1 from allele2
  propensities[7] = B_prod1 *B_proddiff1 *Burst1_on_targ_allele1+B_prod1 *Burst1_off_targ_allele1;
  //update propensity for A1 = 
  propensities[8] = A_deg1 *A1;
  //update propensity for A1_nonsense = 
  propensities[9] = Anonsense_deg1 *A1_nonsense;
  //update propensity for Aprime1 = 
  propensities[10] = Aprime_deg1 *Aprime1;
  //update propensity for B1 = 
  propensities[11] = B_deg1 *B1;
  //update propensity for Burst1_off_orig = Burst1_on_orig for allele1
  propensities[12] = onbasal_a1 *Burst1_off_orig_allele1;
  //update propensity for Burst1_off_orig = Burst1_on_orig for allele2
  propensities[13] = onbasal_a1 *Burst1_off_orig_allele2;
  //update propensity for Burst1_off_para = Burst1_on_para for allele1
  propensities[14] = Aprimenitc1 *(pow(A1_nonsense,nA1_nonsense )/(pow(kA1_nonsense ,nA1_nonsense )+pow(A1_nonsense,nA1_nonsense )))*Burst1_off_para_allele1+onbasal_aprime1 *Burst1_off_para_allele1;
  //update propensity for Burst1_off_para = Burst1_on_para for allele2
  propensities[15] = Aprimenitc1 *(pow(A1_nonsense,nA1_nonsense )/(pow(kA1_nonsense ,nA1_nonsense )+pow(A1_nonsense,nA1_nonsense )))*Burst1_off_para_allele2+onbasal_aprime1 *Burst1_off_para_allele2;
  //update propensity for Burst1_off_targ = Burst1_on_targ for allele1
  propensities[16] = B_ondep1 *(pow(A1,nA1 )/(pow(kA1 ,nA1 )+pow(A1,nA1 )))*Burst1_off_targ_allele1+ B_ondep_prime *(pow(Aprime1,nAprime1 )/(pow(kAprime1 ,nAprime1 )+pow(Aprime1,nAprime1 )))*Burst1_off_targ_allele1 +onbasal_b1 *Burst1_off_targ_allele1;
  //update propensity for Burst1_off_targ = Burst1_on_targ for allele2
  propensities[17] = B_ondep1 *(pow(A1,nA1 )/(pow(kA1 ,nA1 )+pow(A1,nA1 )))*Burst1_off_targ_allele2+ B_ondep_prime *(pow(Aprime1,nAprime1 )/(pow(kAprime1 ,nAprime1 )+pow(Aprime1,nAprime1 )))*Burst1_off_targ_allele2 +onbasal_b1 *Burst1_off_targ_allele2;
  //update propensity for Burst1_on_orig = Burst1_off_orig for allele1
  propensities[18] = A_off1 *Burst1_on_orig_allele1;
  //update propensity for Burst1_on_orig = Burst1_off_orig for allele2
  propensities[19] = A_off1 *Burst1_on_orig_allele2;
  //update propensity for Burst1_on_para = Burst1_off_para for allele1
  propensities[20] = Aprime_off1 *Burst1_on_para_allele1;
  //update propensity for Burst1_on_para = Burst1_off_para for allele2
  propensities[21] = Aprime_off1 *Burst1_on_para_allele2;
  //update propensity for Burst1_on_targ = Burst1_off_targ for allele1
  propensities[22] = B_off1 *Burst1_on_targ_allele1;
  //update propensity for Burst1_on_targ = Burst1_off_targ for allele2
  propensities[23] = B_off1 *Burst1_on_targ_allele2;
} else if (p<cumpropensities[10]) {
  // rxn: Aprime1 = 
  Aprime1=Aprime1 + -1;

   //update propensity for = A1 from allele1
  propensities[0] = Burst1_not_mutated_allele1*(A_prod1 *A_proddiff1 *Burst1_on_orig_allele1+A_prod1 *Burst1_off_orig_allele1);
  //update propensity for = A1_nonsense from allele1
  propensities[1] = Burst1_is_mutated_allele1*(Anonsense_prod1 *A_proddiff1 *Burst1_on_orig_allele1+Anonsense_prod1 *Burst1_off_orig_allele1);
  //update propensity for = A1 from allele2
  propensities[2] = Burst1_not_mutated_allele2*(A_prod1 *A_proddiff1 *Burst1_on_orig_allele2+A_prod1 *Burst1_off_orig_allele2);
  //update propensity for = A1_nonsense from allele2
  propensities[3] = Burst1_is_mutated_allele2*(Anonsense_prod1 *A_proddiff1 *Burst1_on_orig_allele2+Anonsense_prod1 *Burst1_off_orig_allele2);
  //update propensity for = Aprime1 from allele1
  propensities[4] = Aprime_prod1 *Aprime_proddiff1 *Burst1_on_para_allele1+Aprime_prod1 *Burst1_off_para_allele1;
  //update propensity for = Aprime1 from allele2
  propensities[5] = Aprime_prod1 *Aprime_proddiff1 *Burst1_on_para_allele2+Aprime_prod1 *Burst1_off_para_allele2;
  //update propensity for = B1 from allele1
  propensities[6] = B_prod1 *B_proddiff1 *Burst1_on_targ_allele1+B_prod1 *Burst1_off_targ_allele1;
  //update propensity for = B1 from allele2
  propensities[7] = B_prod1 *B_proddiff1 *Burst1_on_targ_allele1+B_prod1 *Burst1_off_targ_allele1;
  //update propensity for A1 = 
  propensities[8] = A_deg1 *A1;
  //update propensity for A1_nonsense = 
  propensities[9] = Anonsense_deg1 *A1_nonsense;
  //update propensity for Aprime1 = 
  propensities[10] = Aprime_deg1 *Aprime1;
  //update propensity for B1 = 
  propensities[11] = B_deg1 *B1;
  //update propensity for Burst1_off_orig = Burst1_on_orig for allele1
  propensities[12] = onbasal_a1 *Burst1_off_orig_allele1;
  //update propensity for Burst1_off_orig = Burst1_on_orig for allele2
  propensities[13] = onbasal_a1 *Burst1_off_orig_allele2;
  //update propensity for Burst1_off_para = Burst1_on_para for allele1
  propensities[14] = Aprimenitc1 *(pow(A1_nonsense,nA1_nonsense )/(pow(kA1_nonsense ,nA1_nonsense )+pow(A1_nonsense,nA1_nonsense )))*Burst1_off_para_allele1+onbasal_aprime1 *Burst1_off_para_allele1;
  //update propensity for Burst1_off_para = Burst1_on_para for allele2
  propensities[15] = Aprimenitc1 *(pow(A1_nonsense,nA1_nonsense )/(pow(kA1_nonsense ,nA1_nonsense )+pow(A1_nonsense,nA1_nonsense )))*Burst1_off_para_allele2+onbasal_aprime1 *Burst1_off_para_allele2;
  //update propensity for Burst1_off_targ = Burst1_on_targ for allele1
  propensities[16] = B_ondep1 *(pow(A1,nA1 )/(pow(kA1 ,nA1 )+pow(A1,nA1 )))*Burst1_off_targ_allele1+ B_ondep_prime *(pow(Aprime1,nAprime1 )/(pow(kAprime1 ,nAprime1 )+pow(Aprime1,nAprime1 )))*Burst1_off_targ_allele1 +onbasal_b1 *Burst1_off_targ_allele1;
  //update propensity for Burst1_off_targ = Burst1_on_targ for allele2
  propensities[17] = B_ondep1 *(pow(A1,nA1 )/(pow(kA1 ,nA1 )+pow(A1,nA1 )))*Burst1_off_targ_allele2+ B_ondep_prime *(pow(Aprime1,nAprime1 )/(pow(kAprime1 ,nAprime1 )+pow(Aprime1,nAprime1 )))*Burst1_off_targ_allele2 +onbasal_b1 *Burst1_off_targ_allele2;
  //update propensity for Burst1_on_orig = Burst1_off_orig for allele1
  propensities[18] = A_off1 *Burst1_on_orig_allele1;
  //update propensity for Burst1_on_orig = Burst1_off_orig for allele2
  propensities[19] = A_off1 *Burst1_on_orig_allele2;
  //update propensity for Burst1_on_para = Burst1_off_para for allele1
  propensities[20] = Aprime_off1 *Burst1_on_para_allele1;
  //update propensity for Burst1_on_para = Burst1_off_para for allele2
  propensities[21] = Aprime_off1 *Burst1_on_para_allele2;
  //update propensity for Burst1_on_targ = Burst1_off_targ for allele1
  propensities[22] = B_off1 *Burst1_on_targ_allele1;
  //update propensity for Burst1_on_targ = Burst1_off_targ for allele2
  propensities[23] = B_off1 *Burst1_on_targ_allele2;
} else if (p<cumpropensities[11]) {
  // rxn: B1 = 
  B1=B1 + -1;

   //update propensity for = A1 from allele1
  propensities[0] = Burst1_not_mutated_allele1*(A_prod1 *A_proddiff1 *Burst1_on_orig_allele1+A_prod1 *Burst1_off_orig_allele1);
  //update propensity for = A1_nonsense from allele1
  propensities[1] = Burst1_is_mutated_allele1*(Anonsense_prod1 *A_proddiff1 *Burst1_on_orig_allele1+Anonsense_prod1 *Burst1_off_orig_allele1);
  //update propensity for = A1 from allele2
  propensities[2] = Burst1_not_mutated_allele2*(A_prod1 *A_proddiff1 *Burst1_on_orig_allele2+A_prod1 *Burst1_off_orig_allele2);
  //update propensity for = A1_nonsense from allele2
  propensities[3] = Burst1_is_mutated_allele2*(Anonsense_prod1 *A_proddiff1 *Burst1_on_orig_allele2+Anonsense_prod1 *Burst1_off_orig_allele2);
  //update propensity for = Aprime1 from allele1
  propensities[4] = Aprime_prod1 *Aprime_proddiff1 *Burst1_on_para_allele1+Aprime_prod1 *Burst1_off_para_allele1;
  //update propensity for = Aprime1 from allele2
  propensities[5] = Aprime_prod1 *Aprime_proddiff1 *Burst1_on_para_allele2+Aprime_prod1 *Burst1_off_para_allele2;
  //update propensity for = B1 from allele1
  propensities[6] = B_prod1 *B_proddiff1 *Burst1_on_targ_allele1+B_prod1 *Burst1_off_targ_allele1;
  //update propensity for = B1 from allele2
  propensities[7] = B_prod1 *B_proddiff1 *Burst1_on_targ_allele1+B_prod1 *Burst1_off_targ_allele1;
  //update propensity for A1 = 
  propensities[8] = A_deg1 *A1;
  //update propensity for A1_nonsense = 
  propensities[9] = Anonsense_deg1 *A1_nonsense;
  //update propensity for Aprime1 = 
  propensities[10] = Aprime_deg1 *Aprime1;
  //update propensity for B1 = 
  propensities[11] = B_deg1 *B1;
  //update propensity for Burst1_off_orig = Burst1_on_orig for allele1
  propensities[12] = onbasal_a1 *Burst1_off_orig_allele1;
  //update propensity for Burst1_off_orig = Burst1_on_orig for allele2
  propensities[13] = onbasal_a1 *Burst1_off_orig_allele2;
  //update propensity for Burst1_off_para = Burst1_on_para for allele1
  propensities[14] = Aprimenitc1 *(pow(A1_nonsense,nA1_nonsense )/(pow(kA1_nonsense ,nA1_nonsense )+pow(A1_nonsense,nA1_nonsense )))*Burst1_off_para_allele1+onbasal_aprime1 *Burst1_off_para_allele1;
  //update propensity for Burst1_off_para = Burst1_on_para for allele2
  propensities[15] = Aprimenitc1 *(pow(A1_nonsense,nA1_nonsense )/(pow(kA1_nonsense ,nA1_nonsense )+pow(A1_nonsense,nA1_nonsense )))*Burst1_off_para_allele2+onbasal_aprime1 *Burst1_off_para_allele2;
  //update propensity for Burst1_off_targ = Burst1_on_targ for allele1
  propensities[16] = B_ondep1 *(pow(A1,nA1 )/(pow(kA1 ,nA1 )+pow(A1,nA1 )))*Burst1_off_targ_allele1+ B_ondep_prime *(pow(Aprime1,nAprime1 )/(pow(kAprime1 ,nAprime1 )+pow(Aprime1,nAprime1 )))*Burst1_off_targ_allele1 +onbasal_b1 *Burst1_off_targ_allele1;
  //update propensity for Burst1_off_targ = Burst1_on_targ for allele2
  propensities[17] = B_ondep1 *(pow(A1,nA1 )/(pow(kA1 ,nA1 )+pow(A1,nA1 )))*Burst1_off_targ_allele2+ B_ondep_prime *(pow(Aprime1,nAprime1 )/(pow(kAprime1 ,nAprime1 )+pow(Aprime1,nAprime1 )))*Burst1_off_targ_allele2 +onbasal_b1 *Burst1_off_targ_allele2;
  //update propensity for Burst1_on_orig = Burst1_off_orig for allele1
  propensities[18] = A_off1 *Burst1_on_orig_allele1;
  //update propensity for Burst1_on_orig = Burst1_off_orig for allele2
  propensities[19] = A_off1 *Burst1_on_orig_allele2;
  //update propensity for Burst1_on_para = Burst1_off_para for allele1
  propensities[20] = Aprime_off1 *Burst1_on_para_allele1;
  //update propensity for Burst1_on_para = Burst1_off_para for allele2
  propensities[21] = Aprime_off1 *Burst1_on_para_allele2;
  //update propensity for Burst1_on_targ = Burst1_off_targ for allele1
  propensities[22] = B_off1 *Burst1_on_targ_allele1;
  //update propensity for Burst1_on_targ = Burst1_off_targ for allele2
  propensities[23] = B_off1 *Burst1_on_targ_allele2;
} else if (p<cumpropensities[12]) {
  // rxn: Burst1_off_orig_allele1 = Burst1_on_orig_allele1
  
  if(Burst1_off_orig_allele1 == 1) {
    Burst1_off_orig_allele1 = 0;
    Burst1_on_orig_allele1 = 1;
  }

  //Burst1_off_orig=Burst1_off_orig + -1;
  //Burst1_on_orig=Burst1_on_orig + 1;
  
   //update propensity for = A1 from allele1
  propensities[0] = Burst1_not_mutated_allele1*(A_prod1 *A_proddiff1 *Burst1_on_orig_allele1+A_prod1 *Burst1_off_orig_allele1);
  //update propensity for = A1_nonsense from allele1
  propensities[1] = Burst1_is_mutated_allele1*(Anonsense_prod1 *A_proddiff1 *Burst1_on_orig_allele1+Anonsense_prod1 *Burst1_off_orig_allele1);
  //update propensity for = A1 from allele2
  propensities[2] = Burst1_not_mutated_allele2*(A_prod1 *A_proddiff1 *Burst1_on_orig_allele2+A_prod1 *Burst1_off_orig_allele2);
  //update propensity for = A1_nonsense from allele2
  propensities[3] = Burst1_is_mutated_allele2*(Anonsense_prod1 *A_proddiff1 *Burst1_on_orig_allele2+Anonsense_prod1 *Burst1_off_orig_allele2);
  //update propensity for = Aprime1 from allele1
  propensities[4] = Aprime_prod1 *Aprime_proddiff1 *Burst1_on_para_allele1+Aprime_prod1 *Burst1_off_para_allele1;
  //update propensity for = Aprime1 from allele2
  propensities[5] = Aprime_prod1 *Aprime_proddiff1 *Burst1_on_para_allele2+Aprime_prod1 *Burst1_off_para_allele2;
  //update propensity for = B1 from allele1
  propensities[6] = B_prod1 *B_proddiff1 *Burst1_on_targ_allele1+B_prod1 *Burst1_off_targ_allele1;
  //update propensity for = B1 from allele2
  propensities[7] = B_prod1 *B_proddiff1 *Burst1_on_targ_allele1+B_prod1 *Burst1_off_targ_allele1;
  //update propensity for A1 = 
  propensities[8] = A_deg1 *A1;
  //update propensity for A1_nonsense = 
  propensities[9] = Anonsense_deg1 *A1_nonsense;
  //update propensity for Aprime1 = 
  propensities[10] = Aprime_deg1 *Aprime1;
  //update propensity for B1 = 
  propensities[11] = B_deg1 *B1;
  //update propensity for Burst1_off_orig = Burst1_on_orig for allele1
  propensities[12] = onbasal_a1 *Burst1_off_orig_allele1;
  //update propensity for Burst1_off_orig = Burst1_on_orig for allele2
  propensities[13] = onbasal_a1 *Burst1_off_orig_allele2;
  //update propensity for Burst1_off_para = Burst1_on_para for allele1
  propensities[14] = Aprimenitc1 *(pow(A1_nonsense,nA1_nonsense )/(pow(kA1_nonsense ,nA1_nonsense )+pow(A1_nonsense,nA1_nonsense )))*Burst1_off_para_allele1+onbasal_aprime1 *Burst1_off_para_allele1;
  //update propensity for Burst1_off_para = Burst1_on_para for allele2
  propensities[15] = Aprimenitc1 *(pow(A1_nonsense,nA1_nonsense )/(pow(kA1_nonsense ,nA1_nonsense )+pow(A1_nonsense,nA1_nonsense )))*Burst1_off_para_allele2+onbasal_aprime1 *Burst1_off_para_allele2;
  //update propensity for Burst1_off_targ = Burst1_on_targ for allele1
  propensities[16] = B_ondep1 *(pow(A1,nA1 )/(pow(kA1 ,nA1 )+pow(A1,nA1 )))*Burst1_off_targ_allele1+ B_ondep_prime *(pow(Aprime1,nAprime1 )/(pow(kAprime1 ,nAprime1 )+pow(Aprime1,nAprime1 )))*Burst1_off_targ_allele1 +onbasal_b1 *Burst1_off_targ_allele1;
  //update propensity for Burst1_off_targ = Burst1_on_targ for allele2
  propensities[17] = B_ondep1 *(pow(A1,nA1 )/(pow(kA1 ,nA1 )+pow(A1,nA1 )))*Burst1_off_targ_allele2+ B_ondep_prime *(pow(Aprime1,nAprime1 )/(pow(kAprime1 ,nAprime1 )+pow(Aprime1,nAprime1 )))*Burst1_off_targ_allele2 +onbasal_b1 *Burst1_off_targ_allele2;
  //update propensity for Burst1_on_orig = Burst1_off_orig for allele1
  propensities[18] = A_off1 *Burst1_on_orig_allele1;
  //update propensity for Burst1_on_orig = Burst1_off_orig for allele2
  propensities[19] = A_off1 *Burst1_on_orig_allele2;
  //update propensity for Burst1_on_para = Burst1_off_para for allele1
  propensities[20] = Aprime_off1 *Burst1_on_para_allele1;
  //update propensity for Burst1_on_para = Burst1_off_para for allele2
  propensities[21] = Aprime_off1 *Burst1_on_para_allele2;
  //update propensity for Burst1_on_targ = Burst1_off_targ for allele1
  propensities[22] = B_off1 *Burst1_on_targ_allele1;
  //update propensity for Burst1_on_targ = Burst1_off_targ for allele2
  propensities[23] = B_off1 *Burst1_on_targ_allele2;
} else if (p<cumpropensities[13]) {
  // rxn: Burst1_off_orig_allele2 = Burst1_on_orig_allele2
  
  if(Burst1_off_orig_allele2 == 1) {
    Burst1_off_orig_allele2 = 0;
    Burst1_on_orig_allele2 = 1;
  }

  //Burst1_off_orig=Burst1_off_orig + -1;
  //Burst1_on_orig=Burst1_on_orig + 1;
  
   //update propensity for = A1 from allele1
  propensities[0] = Burst1_not_mutated_allele1*(A_prod1 *A_proddiff1 *Burst1_on_orig_allele1+A_prod1 *Burst1_off_orig_allele1);
  //update propensity for = A1_nonsense from allele1
  propensities[1] = Burst1_is_mutated_allele1*(Anonsense_prod1 *A_proddiff1 *Burst1_on_orig_allele1+Anonsense_prod1 *Burst1_off_orig_allele1);
  //update propensity for = A1 from allele2
  propensities[2] = Burst1_not_mutated_allele2*(A_prod1 *A_proddiff1 *Burst1_on_orig_allele2+A_prod1 *Burst1_off_orig_allele2);
  //update propensity for = A1_nonsense from allele2
  propensities[3] = Burst1_is_mutated_allele2*(Anonsense_prod1 *A_proddiff1 *Burst1_on_orig_allele2+Anonsense_prod1 *Burst1_off_orig_allele2);
  //update propensity for = Aprime1 from allele1
  propensities[4] = Aprime_prod1 *Aprime_proddiff1 *Burst1_on_para_allele1+Aprime_prod1 *Burst1_off_para_allele1;
  //update propensity for = Aprime1 from allele2
  propensities[5] = Aprime_prod1 *Aprime_proddiff1 *Burst1_on_para_allele2+Aprime_prod1 *Burst1_off_para_allele2;
  //update propensity for = B1 from allele1
  propensities[6] = B_prod1 *B_proddiff1 *Burst1_on_targ_allele1+B_prod1 *Burst1_off_targ_allele1;
  //update propensity for = B1 from allele2
  propensities[7] = B_prod1 *B_proddiff1 *Burst1_on_targ_allele1+B_prod1 *Burst1_off_targ_allele1;
  //update propensity for A1 = 
  propensities[8] = A_deg1 *A1;
  //update propensity for A1_nonsense = 
  propensities[9] = Anonsense_deg1 *A1_nonsense;
  //update propensity for Aprime1 = 
  propensities[10] = Aprime_deg1 *Aprime1;
  //update propensity for B1 = 
  propensities[11] = B_deg1 *B1;
  //update propensity for Burst1_off_orig = Burst1_on_orig for allele1
  propensities[12] = onbasal_a1 *Burst1_off_orig_allele1;
  //update propensity for Burst1_off_orig = Burst1_on_orig for allele2
  propensities[13] = onbasal_a1 *Burst1_off_orig_allele2;
  //update propensity for Burst1_off_para = Burst1_on_para for allele1
  propensities[14] = Aprimenitc1 *(pow(A1_nonsense,nA1_nonsense )/(pow(kA1_nonsense ,nA1_nonsense )+pow(A1_nonsense,nA1_nonsense )))*Burst1_off_para_allele1+onbasal_aprime1 *Burst1_off_para_allele1;
  //update propensity for Burst1_off_para = Burst1_on_para for allele2
  propensities[15] = Aprimenitc1 *(pow(A1_nonsense,nA1_nonsense )/(pow(kA1_nonsense ,nA1_nonsense )+pow(A1_nonsense,nA1_nonsense )))*Burst1_off_para_allele2+onbasal_aprime1 *Burst1_off_para_allele2;
  //update propensity for Burst1_off_targ = Burst1_on_targ for allele1
  propensities[16] = B_ondep1 *(pow(A1,nA1 )/(pow(kA1 ,nA1 )+pow(A1,nA1 )))*Burst1_off_targ_allele1+ B_ondep_prime *(pow(Aprime1,nAprime1 )/(pow(kAprime1 ,nAprime1 )+pow(Aprime1,nAprime1 )))*Burst1_off_targ_allele1 +onbasal_b1 *Burst1_off_targ_allele1;
  //update propensity for Burst1_off_targ = Burst1_on_targ for allele2
  propensities[17] = B_ondep1 *(pow(A1,nA1 )/(pow(kA1 ,nA1 )+pow(A1,nA1 )))*Burst1_off_targ_allele2+ B_ondep_prime *(pow(Aprime1,nAprime1 )/(pow(kAprime1 ,nAprime1 )+pow(Aprime1,nAprime1 )))*Burst1_off_targ_allele2 +onbasal_b1 *Burst1_off_targ_allele2;
  //update propensity for Burst1_on_orig = Burst1_off_orig for allele1
  propensities[18] = A_off1 *Burst1_on_orig_allele1;
  //update propensity for Burst1_on_orig = Burst1_off_orig for allele2
  propensities[19] = A_off1 *Burst1_on_orig_allele2;
  //update propensity for Burst1_on_para = Burst1_off_para for allele1
  propensities[20] = Aprime_off1 *Burst1_on_para_allele1;
  //update propensity for Burst1_on_para = Burst1_off_para for allele2
  propensities[21] = Aprime_off1 *Burst1_on_para_allele2;
  //update propensity for Burst1_on_targ = Burst1_off_targ for allele1
  propensities[22] = B_off1 *Burst1_on_targ_allele1;
  //update propensity for Burst1_on_targ = Burst1_off_targ for allele2
  propensities[23] = B_off1 *Burst1_on_targ_allele2;
} else if (p<cumpropensities[14]) {
  // rxn: Burst1_off_para_allele1 = Burst1_on_para_allele1
  
if(Burst1_off_para_allele1 == 1) {
    Burst1_off_para_allele1 = 0;
    Burst1_on_para_allele1 = 1;
  }

  //Burst1_off_para=Burst1_off_para + -1;
  //Burst1_on_para=Burst1_on_para + 1;

  //update propensity for = A1 from allele1
  propensities[0] = Burst1_not_mutated_allele1*(A_prod1 *A_proddiff1 *Burst1_on_orig_allele1+A_prod1 *Burst1_off_orig_allele1);
  //update propensity for = A1_nonsense from allele1
  propensities[1] = Burst1_is_mutated_allele1*(Anonsense_prod1 *A_proddiff1 *Burst1_on_orig_allele1+Anonsense_prod1 *Burst1_off_orig_allele1);
  //update propensity for = A1 from allele2
  propensities[2] = Burst1_not_mutated_allele2*(A_prod1 *A_proddiff1 *Burst1_on_orig_allele2+A_prod1 *Burst1_off_orig_allele2);
  //update propensity for = A1_nonsense from allele2
  propensities[3] = Burst1_is_mutated_allele2*(Anonsense_prod1 *A_proddiff1 *Burst1_on_orig_allele2+Anonsense_prod1 *Burst1_off_orig_allele2);
  //update propensity for = Aprime1 from allele1
  propensities[4] = Aprime_prod1 *Aprime_proddiff1 *Burst1_on_para_allele1+Aprime_prod1 *Burst1_off_para_allele1;
  //update propensity for = Aprime1 from allele2
  propensities[5] = Aprime_prod1 *Aprime_proddiff1 *Burst1_on_para_allele2+Aprime_prod1 *Burst1_off_para_allele2;
  //update propensity for = B1 from allele1
  propensities[6] = B_prod1 *B_proddiff1 *Burst1_on_targ_allele1+B_prod1 *Burst1_off_targ_allele1;
  //update propensity for = B1 from allele2
  propensities[7] = B_prod1 *B_proddiff1 *Burst1_on_targ_allele1+B_prod1 *Burst1_off_targ_allele1;
  //update propensity for A1 = 
  propensities[8] = A_deg1 *A1;
  //update propensity for A1_nonsense = 
  propensities[9] = Anonsense_deg1 *A1_nonsense;
  //update propensity for Aprime1 = 
  propensities[10] = Aprime_deg1 *Aprime1;
  //update propensity for B1 = 
  propensities[11] = B_deg1 *B1;
  //update propensity for Burst1_off_orig = Burst1_on_orig for allele1
  propensities[12] = onbasal_a1 *Burst1_off_orig_allele1;
  //update propensity for Burst1_off_orig = Burst1_on_orig for allele2
  propensities[13] = onbasal_a1 *Burst1_off_orig_allele2;
  //update propensity for Burst1_off_para = Burst1_on_para for allele1
  propensities[14] = Aprimenitc1 *(pow(A1_nonsense,nA1_nonsense )/(pow(kA1_nonsense ,nA1_nonsense )+pow(A1_nonsense,nA1_nonsense )))*Burst1_off_para_allele1+onbasal_aprime1 *Burst1_off_para_allele1;
  //update propensity for Burst1_off_para = Burst1_on_para for allele2
  propensities[15] = Aprimenitc1 *(pow(A1_nonsense,nA1_nonsense )/(pow(kA1_nonsense ,nA1_nonsense )+pow(A1_nonsense,nA1_nonsense )))*Burst1_off_para_allele2+onbasal_aprime1 *Burst1_off_para_allele2;
  //update propensity for Burst1_off_targ = Burst1_on_targ for allele1
  propensities[16] = B_ondep1 *(pow(A1,nA1 )/(pow(kA1 ,nA1 )+pow(A1,nA1 )))*Burst1_off_targ_allele1+ B_ondep_prime *(pow(Aprime1,nAprime1 )/(pow(kAprime1 ,nAprime1 )+pow(Aprime1,nAprime1 )))*Burst1_off_targ_allele1 +onbasal_b1 *Burst1_off_targ_allele1;
  //update propensity for Burst1_off_targ = Burst1_on_targ for allele2
  propensities[17] = B_ondep1 *(pow(A1,nA1 )/(pow(kA1 ,nA1 )+pow(A1,nA1 )))*Burst1_off_targ_allele2+ B_ondep_prime *(pow(Aprime1,nAprime1 )/(pow(kAprime1 ,nAprime1 )+pow(Aprime1,nAprime1 )))*Burst1_off_targ_allele2 +onbasal_b1 *Burst1_off_targ_allele2;
  //update propensity for Burst1_on_orig = Burst1_off_orig for allele1
  propensities[18] = A_off1 *Burst1_on_orig_allele1;
  //update propensity for Burst1_on_orig = Burst1_off_orig for allele2
  propensities[19] = A_off1 *Burst1_on_orig_allele2;
  //update propensity for Burst1_on_para = Burst1_off_para for allele1
  propensities[20] = Aprime_off1 *Burst1_on_para_allele1;
  //update propensity for Burst1_on_para = Burst1_off_para for allele2
  propensities[21] = Aprime_off1 *Burst1_on_para_allele2;
  //update propensity for Burst1_on_targ = Burst1_off_targ for allele1
  propensities[22] = B_off1 *Burst1_on_targ_allele1;
  //update propensity for Burst1_on_targ = Burst1_off_targ for allele2
  propensities[23] = B_off1 *Burst1_on_targ_allele2;
} else if (p<cumpropensities[15]) {
  // rxn: Burst1_off_para_allele2 = Burst1_on_para_allele2
  
if(Burst1_off_para_allele2 == 1) {
    Burst1_off_para_allele2 = 0;
    Burst1_on_para_allele2 = 1;
  }

  //Burst1_off_para=Burst1_off_para + -1;
  //Burst1_on_para=Burst1_on_para + 1;

  //update propensity for = A1 from allele1
  propensities[0] = Burst1_not_mutated_allele1*(A_prod1 *A_proddiff1 *Burst1_on_orig_allele1+A_prod1 *Burst1_off_orig_allele1);
  //update propensity for = A1_nonsense from allele1
  propensities[1] = Burst1_is_mutated_allele1*(Anonsense_prod1 *A_proddiff1 *Burst1_on_orig_allele1+Anonsense_prod1 *Burst1_off_orig_allele1);
  //update propensity for = A1 from allele2
  propensities[2] = Burst1_not_mutated_allele2*(A_prod1 *A_proddiff1 *Burst1_on_orig_allele2+A_prod1 *Burst1_off_orig_allele2);
  //update propensity for = A1_nonsense from allele2
  propensities[3] = Burst1_is_mutated_allele2*(Anonsense_prod1 *A_proddiff1 *Burst1_on_orig_allele2+Anonsense_prod1 *Burst1_off_orig_allele2);
  //update propensity for = Aprime1 from allele1
  propensities[4] = Aprime_prod1 *Aprime_proddiff1 *Burst1_on_para_allele1+Aprime_prod1 *Burst1_off_para_allele1;
  //update propensity for = Aprime1 from allele2
  propensities[5] = Aprime_prod1 *Aprime_proddiff1 *Burst1_on_para_allele2+Aprime_prod1 *Burst1_off_para_allele2;
  //update propensity for = B1 from allele1
  propensities[6] = B_prod1 *B_proddiff1 *Burst1_on_targ_allele1+B_prod1 *Burst1_off_targ_allele1;
  //update propensity for = B1 from allele2
  propensities[7] = B_prod1 *B_proddiff1 *Burst1_on_targ_allele1+B_prod1 *Burst1_off_targ_allele1;
  //update propensity for A1 = 
  propensities[8] = A_deg1 *A1;
  //update propensity for A1_nonsense = 
  propensities[9] = Anonsense_deg1 *A1_nonsense;
  //update propensity for Aprime1 = 
  propensities[10] = Aprime_deg1 *Aprime1;
  //update propensity for B1 = 
  propensities[11] = B_deg1 *B1;
  //update propensity for Burst1_off_orig = Burst1_on_orig for allele1
  propensities[12] = onbasal_a1 *Burst1_off_orig_allele1;
  //update propensity for Burst1_off_orig = Burst1_on_orig for allele2
  propensities[13] = onbasal_a1 *Burst1_off_orig_allele2;
  //update propensity for Burst1_off_para = Burst1_on_para for allele1
  propensities[14] = Aprimenitc1 *(pow(A1_nonsense,nA1_nonsense )/(pow(kA1_nonsense ,nA1_nonsense )+pow(A1_nonsense,nA1_nonsense )))*Burst1_off_para_allele1+onbasal_aprime1 *Burst1_off_para_allele1;
  //update propensity for Burst1_off_para = Burst1_on_para for allele2
  propensities[15] = Aprimenitc1 *(pow(A1_nonsense,nA1_nonsense )/(pow(kA1_nonsense ,nA1_nonsense )+pow(A1_nonsense,nA1_nonsense )))*Burst1_off_para_allele2+onbasal_aprime1 *Burst1_off_para_allele2;
  //update propensity for Burst1_off_targ = Burst1_on_targ for allele1
  propensities[16] = B_ondep1 *(pow(A1,nA1 )/(pow(kA1 ,nA1 )+pow(A1,nA1 )))*Burst1_off_targ_allele1+ B_ondep_prime *(pow(Aprime1,nAprime1 )/(pow(kAprime1 ,nAprime1 )+pow(Aprime1,nAprime1 )))*Burst1_off_targ_allele1 +onbasal_b1 *Burst1_off_targ_allele1;
  //update propensity for Burst1_off_targ = Burst1_on_targ for allele2
  propensities[17] = B_ondep1 *(pow(A1,nA1 )/(pow(kA1 ,nA1 )+pow(A1,nA1 )))*Burst1_off_targ_allele2+ B_ondep_prime *(pow(Aprime1,nAprime1 )/(pow(kAprime1 ,nAprime1 )+pow(Aprime1,nAprime1 )))*Burst1_off_targ_allele2 +onbasal_b1 *Burst1_off_targ_allele2;
  //update propensity for Burst1_on_orig = Burst1_off_orig for allele1
  propensities[18] = A_off1 *Burst1_on_orig_allele1;
  //update propensity for Burst1_on_orig = Burst1_off_orig for allele2
  propensities[19] = A_off1 *Burst1_on_orig_allele2;
  //update propensity for Burst1_on_para = Burst1_off_para for allele1
  propensities[20] = Aprime_off1 *Burst1_on_para_allele1;
  //update propensity for Burst1_on_para = Burst1_off_para for allele2
  propensities[21] = Aprime_off1 *Burst1_on_para_allele2;
  //update propensity for Burst1_on_targ = Burst1_off_targ for allele1
  propensities[22] = B_off1 *Burst1_on_targ_allele1;
  //update propensity for Burst1_on_targ = Burst1_off_targ for allele2
  propensities[23] = B_off1 *Burst1_on_targ_allele2;
} else if (p<cumpropensities[16]) {
  // rxn: Burst1_off_targ_allele1 = Burst1_on_targ_allele1
  
  if(Burst1_off_targ_allele1 == 1) {
    Burst1_off_targ_allele1 = 0;
    Burst1_on_targ_allele1 = 1;
  }

  //Burst1_off_targ=Burst1_off_targ + -1;
  //Burst1_on_targ=Burst1_on_targ + 1;

  //update propensity for = A1 from allele1
  propensities[0] = Burst1_not_mutated_allele1*(A_prod1 *A_proddiff1 *Burst1_on_orig_allele1+A_prod1 *Burst1_off_orig_allele1);
  //update propensity for = A1_nonsense from allele1
  propensities[1] = Burst1_is_mutated_allele1*(Anonsense_prod1 *A_proddiff1 *Burst1_on_orig_allele1+Anonsense_prod1 *Burst1_off_orig_allele1);
  //update propensity for = A1 from allele2
  propensities[2] = Burst1_not_mutated_allele2*(A_prod1 *A_proddiff1 *Burst1_on_orig_allele2+A_prod1 *Burst1_off_orig_allele2);
  //update propensity for = A1_nonsense from allele2
  propensities[3] = Burst1_is_mutated_allele2*(Anonsense_prod1 *A_proddiff1 *Burst1_on_orig_allele2+Anonsense_prod1 *Burst1_off_orig_allele2);
  //update propensity for = Aprime1 from allele1
  propensities[4] = Aprime_prod1 *Aprime_proddiff1 *Burst1_on_para_allele1+Aprime_prod1 *Burst1_off_para_allele1;
  //update propensity for = Aprime1 from allele2
  propensities[5] = Aprime_prod1 *Aprime_proddiff1 *Burst1_on_para_allele2+Aprime_prod1 *Burst1_off_para_allele2;
  //update propensity for = B1 from allele1
  propensities[6] = B_prod1 *B_proddiff1 *Burst1_on_targ_allele1+B_prod1 *Burst1_off_targ_allele1;
  //update propensity for = B1 from allele2
  propensities[7] = B_prod1 *B_proddiff1 *Burst1_on_targ_allele1+B_prod1 *Burst1_off_targ_allele1;
  //update propensity for A1 = 
  propensities[8] = A_deg1 *A1;
  //update propensity for A1_nonsense = 
  propensities[9] = Anonsense_deg1 *A1_nonsense;
  //update propensity for Aprime1 = 
  propensities[10] = Aprime_deg1 *Aprime1;
  //update propensity for B1 = 
  propensities[11] = B_deg1 *B1;
  //update propensity for Burst1_off_orig = Burst1_on_orig for allele1
  propensities[12] = onbasal_a1 *Burst1_off_orig_allele1;
  //update propensity for Burst1_off_orig = Burst1_on_orig for allele2
  propensities[13] = onbasal_a1 *Burst1_off_orig_allele2;
  //update propensity for Burst1_off_para = Burst1_on_para for allele1
  propensities[14] = Aprimenitc1 *(pow(A1_nonsense,nA1_nonsense )/(pow(kA1_nonsense ,nA1_nonsense )+pow(A1_nonsense,nA1_nonsense )))*Burst1_off_para_allele1+onbasal_aprime1 *Burst1_off_para_allele1;
  //update propensity for Burst1_off_para = Burst1_on_para for allele2
  propensities[15] = Aprimenitc1 *(pow(A1_nonsense,nA1_nonsense )/(pow(kA1_nonsense ,nA1_nonsense )+pow(A1_nonsense,nA1_nonsense )))*Burst1_off_para_allele2+onbasal_aprime1 *Burst1_off_para_allele2;
  //update propensity for Burst1_off_targ = Burst1_on_targ for allele1
  propensities[16] = B_ondep1 *(pow(A1,nA1 )/(pow(kA1 ,nA1 )+pow(A1,nA1 )))*Burst1_off_targ_allele1+ B_ondep_prime *(pow(Aprime1,nAprime1 )/(pow(kAprime1 ,nAprime1 )+pow(Aprime1,nAprime1 )))*Burst1_off_targ_allele1 +onbasal_b1 *Burst1_off_targ_allele1;
  //update propensity for Burst1_off_targ = Burst1_on_targ for allele2
  propensities[17] = B_ondep1 *(pow(A1,nA1 )/(pow(kA1 ,nA1 )+pow(A1,nA1 )))*Burst1_off_targ_allele2+ B_ondep_prime *(pow(Aprime1,nAprime1 )/(pow(kAprime1 ,nAprime1 )+pow(Aprime1,nAprime1 )))*Burst1_off_targ_allele2 +onbasal_b1 *Burst1_off_targ_allele2;
  //update propensity for Burst1_on_orig = Burst1_off_orig for allele1
  propensities[18] = A_off1 *Burst1_on_orig_allele1;
  //update propensity for Burst1_on_orig = Burst1_off_orig for allele2
  propensities[19] = A_off1 *Burst1_on_orig_allele2;
  //update propensity for Burst1_on_para = Burst1_off_para for allele1
  propensities[20] = Aprime_off1 *Burst1_on_para_allele1;
  //update propensity for Burst1_on_para = Burst1_off_para for allele2
  propensities[21] = Aprime_off1 *Burst1_on_para_allele2;
  //update propensity for Burst1_on_targ = Burst1_off_targ for allele1
  propensities[22] = B_off1 *Burst1_on_targ_allele1;
  //update propensity for Burst1_on_targ = Burst1_off_targ for allele2
  propensities[23] = B_off1 *Burst1_on_targ_allele2;
} else if (p<cumpropensities[17]) {
  // rxn: Burst1_off_targ_allele2 = Burst1_on_targ_allele2
  
  if(Burst1_off_targ_allele2 == 1) {
    Burst1_off_targ_allele2 = 0;
    Burst1_on_targ_allele2 = 1;
  }

  //Burst1_off_targ=Burst1_off_targ + -1;
  //Burst1_on_targ=Burst1_on_targ + 1;

  //update propensity for = A1 from allele1
  propensities[0] = Burst1_not_mutated_allele1*(A_prod1 *A_proddiff1 *Burst1_on_orig_allele1+A_prod1 *Burst1_off_orig_allele1);
  //update propensity for = A1_nonsense from allele1
  propensities[1] = Burst1_is_mutated_allele1*(Anonsense_prod1 *A_proddiff1 *Burst1_on_orig_allele1+Anonsense_prod1 *Burst1_off_orig_allele1);
  //update propensity for = A1 from allele2
  propensities[2] = Burst1_not_mutated_allele2*(A_prod1 *A_proddiff1 *Burst1_on_orig_allele2+A_prod1 *Burst1_off_orig_allele2);
  //update propensity for = A1_nonsense from allele2
  propensities[3] = Burst1_is_mutated_allele2*(Anonsense_prod1 *A_proddiff1 *Burst1_on_orig_allele2+Anonsense_prod1 *Burst1_off_orig_allele2);
  //update propensity for = Aprime1 from allele1
  propensities[4] = Aprime_prod1 *Aprime_proddiff1 *Burst1_on_para_allele1+Aprime_prod1 *Burst1_off_para_allele1;
  //update propensity for = Aprime1 from allele2
  propensities[5] = Aprime_prod1 *Aprime_proddiff1 *Burst1_on_para_allele2+Aprime_prod1 *Burst1_off_para_allele2;
  //update propensity for = B1 from allele1
  propensities[6] = B_prod1 *B_proddiff1 *Burst1_on_targ_allele1+B_prod1 *Burst1_off_targ_allele1;
  //update propensity for = B1 from allele2
  propensities[7] = B_prod1 *B_proddiff1 *Burst1_on_targ_allele1+B_prod1 *Burst1_off_targ_allele1;
  //update propensity for A1 = 
  propensities[8] = A_deg1 *A1;
  //update propensity for A1_nonsense = 
  propensities[9] = Anonsense_deg1 *A1_nonsense;
  //update propensity for Aprime1 = 
  propensities[10] = Aprime_deg1 *Aprime1;
  //update propensity for B1 = 
  propensities[11] = B_deg1 *B1;
  //update propensity for Burst1_off_orig = Burst1_on_orig for allele1
  propensities[12] = onbasal_a1 *Burst1_off_orig_allele1;
  //update propensity for Burst1_off_orig = Burst1_on_orig for allele2
  propensities[13] = onbasal_a1 *Burst1_off_orig_allele2;
  //update propensity for Burst1_off_para = Burst1_on_para for allele1
  propensities[14] = Aprimenitc1 *(pow(A1_nonsense,nA1_nonsense )/(pow(kA1_nonsense ,nA1_nonsense )+pow(A1_nonsense,nA1_nonsense )))*Burst1_off_para_allele1+onbasal_aprime1 *Burst1_off_para_allele1;
  //update propensity for Burst1_off_para = Burst1_on_para for allele2
  propensities[15] = Aprimenitc1 *(pow(A1_nonsense,nA1_nonsense )/(pow(kA1_nonsense ,nA1_nonsense )+pow(A1_nonsense,nA1_nonsense )))*Burst1_off_para_allele2+onbasal_aprime1 *Burst1_off_para_allele2;
  //update propensity for Burst1_off_targ = Burst1_on_targ for allele1
  propensities[16] = B_ondep1 *(pow(A1,nA1 )/(pow(kA1 ,nA1 )+pow(A1,nA1 )))*Burst1_off_targ_allele1+ B_ondep_prime *(pow(Aprime1,nAprime1 )/(pow(kAprime1 ,nAprime1 )+pow(Aprime1,nAprime1 )))*Burst1_off_targ_allele1 +onbasal_b1 *Burst1_off_targ_allele1;
  //update propensity for Burst1_off_targ = Burst1_on_targ for allele2
  propensities[17] = B_ondep1 *(pow(A1,nA1 )/(pow(kA1 ,nA1 )+pow(A1,nA1 )))*Burst1_off_targ_allele2+ B_ondep_prime *(pow(Aprime1,nAprime1 )/(pow(kAprime1 ,nAprime1 )+pow(Aprime1,nAprime1 )))*Burst1_off_targ_allele2 +onbasal_b1 *Burst1_off_targ_allele2;
  //update propensity for Burst1_on_orig = Burst1_off_orig for allele1
  propensities[18] = A_off1 *Burst1_on_orig_allele1;
  //update propensity for Burst1_on_orig = Burst1_off_orig for allele2
  propensities[19] = A_off1 *Burst1_on_orig_allele2;
  //update propensity for Burst1_on_para = Burst1_off_para for allele1
  propensities[20] = Aprime_off1 *Burst1_on_para_allele1;
  //update propensity for Burst1_on_para = Burst1_off_para for allele2
  propensities[21] = Aprime_off1 *Burst1_on_para_allele2;
  //update propensity for Burst1_on_targ = Burst1_off_targ for allele1
  propensities[22] = B_off1 *Burst1_on_targ_allele1;
  //update propensity for Burst1_on_targ = Burst1_off_targ for allele2
  propensities[23] = B_off1 *Burst1_on_targ_allele2;
} else if (p<cumpropensities[18]) {
  // rxn: Burst1_on_orig_allele1   = Burst1_off_orig_allele1

  if(Burst1_on_orig_allele1 == 1) {
    Burst1_off_orig_allele1 = 1;
    Burst1_on_orig_allele1 = 0;
  }
  //Burst1_on_orig=Burst1_on_orig + -1;
  //Burst1_off_orig=Burst1_off_orig + 1;

  //update propensity for = A1 from allele1
  propensities[0] = Burst1_not_mutated_allele1*(A_prod1 *A_proddiff1 *Burst1_on_orig_allele1+A_prod1 *Burst1_off_orig_allele1);
  //update propensity for = A1_nonsense from allele1
  propensities[1] = Burst1_is_mutated_allele1*(Anonsense_prod1 *A_proddiff1 *Burst1_on_orig_allele1+Anonsense_prod1 *Burst1_off_orig_allele1);
  //update propensity for = A1 from allele2
  propensities[2] = Burst1_not_mutated_allele2*(A_prod1 *A_proddiff1 *Burst1_on_orig_allele2+A_prod1 *Burst1_off_orig_allele2);
  //update propensity for = A1_nonsense from allele2
  propensities[3] = Burst1_is_mutated_allele2*(Anonsense_prod1 *A_proddiff1 *Burst1_on_orig_allele2+Anonsense_prod1 *Burst1_off_orig_allele2);
  //update propensity for = Aprime1 from allele1
  propensities[4] = Aprime_prod1 *Aprime_proddiff1 *Burst1_on_para_allele1+Aprime_prod1 *Burst1_off_para_allele1;
  //update propensity for = Aprime1 from allele2
  propensities[5] = Aprime_prod1 *Aprime_proddiff1 *Burst1_on_para_allele2+Aprime_prod1 *Burst1_off_para_allele2;
  //update propensity for = B1 from allele1
  propensities[6] = B_prod1 *B_proddiff1 *Burst1_on_targ_allele1+B_prod1 *Burst1_off_targ_allele1;
  //update propensity for = B1 from allele2
  propensities[7] = B_prod1 *B_proddiff1 *Burst1_on_targ_allele1+B_prod1 *Burst1_off_targ_allele1;
  //update propensity for A1 = 
  propensities[8] = A_deg1 *A1;
  //update propensity for A1_nonsense = 
  propensities[9] = Anonsense_deg1 *A1_nonsense;
  //update propensity for Aprime1 = 
  propensities[10] = Aprime_deg1 *Aprime1;
  //update propensity for B1 = 
  propensities[11] = B_deg1 *B1;
  //update propensity for Burst1_off_orig = Burst1_on_orig for allele1
  propensities[12] = onbasal_a1 *Burst1_off_orig_allele1;
  //update propensity for Burst1_off_orig = Burst1_on_orig for allele2
  propensities[13] = onbasal_a1 *Burst1_off_orig_allele2;
  //update propensity for Burst1_off_para = Burst1_on_para for allele1
  propensities[14] = Aprimenitc1 *(pow(A1_nonsense,nA1_nonsense )/(pow(kA1_nonsense ,nA1_nonsense )+pow(A1_nonsense,nA1_nonsense )))*Burst1_off_para_allele1+onbasal_aprime1 *Burst1_off_para_allele1;
  //update propensity for Burst1_off_para = Burst1_on_para for allele2
  propensities[15] = Aprimenitc1 *(pow(A1_nonsense,nA1_nonsense )/(pow(kA1_nonsense ,nA1_nonsense )+pow(A1_nonsense,nA1_nonsense )))*Burst1_off_para_allele2+onbasal_aprime1 *Burst1_off_para_allele2;
  //update propensity for Burst1_off_targ = Burst1_on_targ for allele1
  propensities[16] = B_ondep1 *(pow(A1,nA1 )/(pow(kA1 ,nA1 )+pow(A1,nA1 )))*Burst1_off_targ_allele1+ B_ondep_prime *(pow(Aprime1,nAprime1 )/(pow(kAprime1 ,nAprime1 )+pow(Aprime1,nAprime1 )))*Burst1_off_targ_allele1 +onbasal_b1 *Burst1_off_targ_allele1;
  //update propensity for Burst1_off_targ = Burst1_on_targ for allele2
  propensities[17] = B_ondep1 *(pow(A1,nA1 )/(pow(kA1 ,nA1 )+pow(A1,nA1 )))*Burst1_off_targ_allele2+ B_ondep_prime *(pow(Aprime1,nAprime1 )/(pow(kAprime1 ,nAprime1 )+pow(Aprime1,nAprime1 )))*Burst1_off_targ_allele2 +onbasal_b1 *Burst1_off_targ_allele2;
  //update propensity for Burst1_on_orig = Burst1_off_orig for allele1
  propensities[18] = A_off1 *Burst1_on_orig_allele1;
  //update propensity for Burst1_on_orig = Burst1_off_orig for allele2
  propensities[19] = A_off1 *Burst1_on_orig_allele2;
  //update propensity for Burst1_on_para = Burst1_off_para for allele1
  propensities[20] = Aprime_off1 *Burst1_on_para_allele1;
  //update propensity for Burst1_on_para = Burst1_off_para for allele2
  propensities[21] = Aprime_off1 *Burst1_on_para_allele2;
  //update propensity for Burst1_on_targ = Burst1_off_targ for allele1
  propensities[22] = B_off1 *Burst1_on_targ_allele1;
  //update propensity for Burst1_on_targ = Burst1_off_targ for allele2
  propensities[23] = B_off1 *Burst1_on_targ_allele2;
  } else if (p<cumpropensities[19]) {
  // rxn: Burst1_on_orig_allele2   = Burst1_off_orig_allele2

  if(Burst1_on_orig_allele2 == 1) {
    Burst1_off_orig_allele2 = 1;
    Burst1_on_orig_allele2 = 0;
  }
  //Burst1_on_orig=Burst1_on_orig + -1;
  //Burst1_off_orig=Burst1_off_orig + 1;

  //update propensity for = A1 from allele1
  propensities[0] = Burst1_not_mutated_allele1*(A_prod1 *A_proddiff1 *Burst1_on_orig_allele1+A_prod1 *Burst1_off_orig_allele1);
  //update propensity for = A1_nonsense from allele1
  propensities[1] = Burst1_is_mutated_allele1*(Anonsense_prod1 *A_proddiff1 *Burst1_on_orig_allele1+Anonsense_prod1 *Burst1_off_orig_allele1);
  //update propensity for = A1 from allele2
  propensities[2] = Burst1_not_mutated_allele2*(A_prod1 *A_proddiff1 *Burst1_on_orig_allele2+A_prod1 *Burst1_off_orig_allele2);
  //update propensity for = A1_nonsense from allele2
  propensities[3] = Burst1_is_mutated_allele2*(Anonsense_prod1 *A_proddiff1 *Burst1_on_orig_allele2+Anonsense_prod1 *Burst1_off_orig_allele2);
  //update propensity for = Aprime1 from allele1
  propensities[4] = Aprime_prod1 *Aprime_proddiff1 *Burst1_on_para_allele1+Aprime_prod1 *Burst1_off_para_allele1;
  //update propensity for = Aprime1 from allele2
  propensities[5] = Aprime_prod1 *Aprime_proddiff1 *Burst1_on_para_allele2+Aprime_prod1 *Burst1_off_para_allele2;
  //update propensity for = B1 from allele1
  propensities[6] = B_prod1 *B_proddiff1 *Burst1_on_targ_allele1+B_prod1 *Burst1_off_targ_allele1;
  //update propensity for = B1 from allele2
  propensities[7] = B_prod1 *B_proddiff1 *Burst1_on_targ_allele1+B_prod1 *Burst1_off_targ_allele1;
  //update propensity for A1 = 
  propensities[8] = A_deg1 *A1;
  //update propensity for A1_nonsense = 
  propensities[9] = Anonsense_deg1 *A1_nonsense;
  //update propensity for Aprime1 = 
  propensities[10] = Aprime_deg1 *Aprime1;
  //update propensity for B1 = 
  propensities[11] = B_deg1 *B1;
  //update propensity for Burst1_off_orig = Burst1_on_orig for allele1
  propensities[12] = onbasal_a1 *Burst1_off_orig_allele1;
  //update propensity for Burst1_off_orig = Burst1_on_orig for allele2
  propensities[13] = onbasal_a1 *Burst1_off_orig_allele2;
  //update propensity for Burst1_off_para = Burst1_on_para for allele1
  propensities[14] = Aprimenitc1 *(pow(A1_nonsense,nA1_nonsense )/(pow(kA1_nonsense ,nA1_nonsense )+pow(A1_nonsense,nA1_nonsense )))*Burst1_off_para_allele1+onbasal_aprime1 *Burst1_off_para_allele1;
  //update propensity for Burst1_off_para = Burst1_on_para for allele2
  propensities[15] = Aprimenitc1 *(pow(A1_nonsense,nA1_nonsense )/(pow(kA1_nonsense ,nA1_nonsense )+pow(A1_nonsense,nA1_nonsense )))*Burst1_off_para_allele2+onbasal_aprime1 *Burst1_off_para_allele2;
  //update propensity for Burst1_off_targ = Burst1_on_targ for allele1
  propensities[16] = B_ondep1 *(pow(A1,nA1 )/(pow(kA1 ,nA1 )+pow(A1,nA1 )))*Burst1_off_targ_allele1+ B_ondep_prime *(pow(Aprime1,nAprime1 )/(pow(kAprime1 ,nAprime1 )+pow(Aprime1,nAprime1 )))*Burst1_off_targ_allele1 +onbasal_b1 *Burst1_off_targ_allele1;
  //update propensity for Burst1_off_targ = Burst1_on_targ for allele2
  propensities[17] = B_ondep1 *(pow(A1,nA1 )/(pow(kA1 ,nA1 )+pow(A1,nA1 )))*Burst1_off_targ_allele2+ B_ondep_prime *(pow(Aprime1,nAprime1 )/(pow(kAprime1 ,nAprime1 )+pow(Aprime1,nAprime1 )))*Burst1_off_targ_allele2 +onbasal_b1 *Burst1_off_targ_allele2;
  //update propensity for Burst1_on_orig = Burst1_off_orig for allele1
  propensities[18] = A_off1 *Burst1_on_orig_allele1;
  //update propensity for Burst1_on_orig = Burst1_off_orig for allele2
  propensities[19] = A_off1 *Burst1_on_orig_allele2;
  //update propensity for Burst1_on_para = Burst1_off_para for allele1
  propensities[20] = Aprime_off1 *Burst1_on_para_allele1;
  //update propensity for Burst1_on_para = Burst1_off_para for allele2
  propensities[21] = Aprime_off1 *Burst1_on_para_allele2;
  //update propensity for Burst1_on_targ = Burst1_off_targ for allele1
  propensities[22] = B_off1 *Burst1_on_targ_allele1;
  //update propensity for Burst1_on_targ = Burst1_off_targ for allele2
  propensities[23] = B_off1 *Burst1_on_targ_allele2;
} else if (p<cumpropensities[20]) {
  // rxn: Burst1_on_para_allele1  = Burst1_off_para_allele1
  
  if(Burst1_on_para_allele1 == 1) {
    Burst1_off_para_allele1 = 1;
    Burst1_on_para_allele1 = 0;
  }
  //Burst1_off_para=Burst1_off_para + 1;
  //Burst1_on_para=Burst1_on_para + -1;

  //update propensity for = A1 from allele1
  propensities[0] = Burst1_not_mutated_allele1*(A_prod1 *A_proddiff1 *Burst1_on_orig_allele1+A_prod1 *Burst1_off_orig_allele1);
  //update propensity for = A1_nonsense from allele1
  propensities[1] = Burst1_is_mutated_allele1*(Anonsense_prod1 *A_proddiff1 *Burst1_on_orig_allele1+Anonsense_prod1 *Burst1_off_orig_allele1);
  //update propensity for = A1 from allele2
  propensities[2] = Burst1_not_mutated_allele2*(A_prod1 *A_proddiff1 *Burst1_on_orig_allele2+A_prod1 *Burst1_off_orig_allele2);
  //update propensity for = A1_nonsense from allele2
  propensities[3] = Burst1_is_mutated_allele2*(Anonsense_prod1 *A_proddiff1 *Burst1_on_orig_allele2+Anonsense_prod1 *Burst1_off_orig_allele2);
  //update propensity for = Aprime1 from allele1
  propensities[4] = Aprime_prod1 *Aprime_proddiff1 *Burst1_on_para_allele1+Aprime_prod1 *Burst1_off_para_allele1;
  //update propensity for = Aprime1 from allele2
  propensities[5] = Aprime_prod1 *Aprime_proddiff1 *Burst1_on_para_allele2+Aprime_prod1 *Burst1_off_para_allele2;
  //update propensity for = B1 from allele1
  propensities[6] = B_prod1 *B_proddiff1 *Burst1_on_targ_allele1+B_prod1 *Burst1_off_targ_allele1;
  //update propensity for = B1 from allele2
  propensities[7] = B_prod1 *B_proddiff1 *Burst1_on_targ_allele1+B_prod1 *Burst1_off_targ_allele1;
  //update propensity for A1 = 
  propensities[8] = A_deg1 *A1;
  //update propensity for A1_nonsense = 
  propensities[9] = Anonsense_deg1 *A1_nonsense;
  //update propensity for Aprime1 = 
  propensities[10] = Aprime_deg1 *Aprime1;
  //update propensity for B1 = 
  propensities[11] = B_deg1 *B1;
  //update propensity for Burst1_off_orig = Burst1_on_orig for allele1
  propensities[12] = onbasal_a1 *Burst1_off_orig_allele1;
  //update propensity for Burst1_off_orig = Burst1_on_orig for allele2
  propensities[13] = onbasal_a1 *Burst1_off_orig_allele2;
  //update propensity for Burst1_off_para = Burst1_on_para for allele1
  propensities[14] = Aprimenitc1 *(pow(A1_nonsense,nA1_nonsense )/(pow(kA1_nonsense ,nA1_nonsense )+pow(A1_nonsense,nA1_nonsense )))*Burst1_off_para_allele1+onbasal_aprime1 *Burst1_off_para_allele1;
  //update propensity for Burst1_off_para = Burst1_on_para for allele2
  propensities[15] = Aprimenitc1 *(pow(A1_nonsense,nA1_nonsense )/(pow(kA1_nonsense ,nA1_nonsense )+pow(A1_nonsense,nA1_nonsense )))*Burst1_off_para_allele2+onbasal_aprime1 *Burst1_off_para_allele2;
  //update propensity for Burst1_off_targ = Burst1_on_targ for allele1
  propensities[16] = B_ondep1 *(pow(A1,nA1 )/(pow(kA1 ,nA1 )+pow(A1,nA1 )))*Burst1_off_targ_allele1+ B_ondep_prime *(pow(Aprime1,nAprime1 )/(pow(kAprime1 ,nAprime1 )+pow(Aprime1,nAprime1 )))*Burst1_off_targ_allele1 +onbasal_b1 *Burst1_off_targ_allele1;
  //update propensity for Burst1_off_targ = Burst1_on_targ for allele2
  propensities[17] = B_ondep1 *(pow(A1,nA1 )/(pow(kA1 ,nA1 )+pow(A1,nA1 )))*Burst1_off_targ_allele2+ B_ondep_prime *(pow(Aprime1,nAprime1 )/(pow(kAprime1 ,nAprime1 )+pow(Aprime1,nAprime1 )))*Burst1_off_targ_allele2 +onbasal_b1 *Burst1_off_targ_allele2;
  //update propensity for Burst1_on_orig = Burst1_off_orig for allele1
  propensities[18] = A_off1 *Burst1_on_orig_allele1;
  //update propensity for Burst1_on_orig = Burst1_off_orig for allele2
  propensities[19] = A_off1 *Burst1_on_orig_allele2;
  //update propensity for Burst1_on_para = Burst1_off_para for allele1
  propensities[20] = Aprime_off1 *Burst1_on_para_allele1;
  //update propensity for Burst1_on_para = Burst1_off_para for allele2
  propensities[21] = Aprime_off1 *Burst1_on_para_allele2;
  //update propensity for Burst1_on_targ = Burst1_off_targ for allele1
  propensities[22] = B_off1 *Burst1_on_targ_allele1;
  //update propensity for Burst1_on_targ = Burst1_off_targ for allele2
  propensities[23] = B_off1 *Burst1_on_targ_allele2;
} else if (p<cumpropensities[21]) {
  // rxn: Burst1_on_para_allele2  = Burst1_off_para_allele2
  
  if(Burst1_on_para_allele2 == 1) {
    Burst1_off_para_allele2 = 1;
    Burst1_on_para_allele2 = 0;
  }
  //Burst1_off_para=Burst1_off_para + 1;
  //Burst1_on_para=Burst1_on_para + -1;

  //update propensity for = A1 from allele1
  propensities[0] = Burst1_not_mutated_allele1*(A_prod1 *A_proddiff1 *Burst1_on_orig_allele1+A_prod1 *Burst1_off_orig_allele1);
  //update propensity for = A1_nonsense from allele1
  propensities[1] = Burst1_is_mutated_allele1*(Anonsense_prod1 *A_proddiff1 *Burst1_on_orig_allele1+Anonsense_prod1 *Burst1_off_orig_allele1);
  //update propensity for = A1 from allele2
  propensities[2] = Burst1_not_mutated_allele2*(A_prod1 *A_proddiff1 *Burst1_on_orig_allele2+A_prod1 *Burst1_off_orig_allele2);
  //update propensity for = A1_nonsense from allele2
  propensities[3] = Burst1_is_mutated_allele2*(Anonsense_prod1 *A_proddiff1 *Burst1_on_orig_allele2+Anonsense_prod1 *Burst1_off_orig_allele2);
  //update propensity for = Aprime1 from allele1
  propensities[4] = Aprime_prod1 *Aprime_proddiff1 *Burst1_on_para_allele1+Aprime_prod1 *Burst1_off_para_allele1;
  //update propensity for = Aprime1 from allele2
  propensities[5] = Aprime_prod1 *Aprime_proddiff1 *Burst1_on_para_allele2+Aprime_prod1 *Burst1_off_para_allele2;
  //update propensity for = B1 from allele1
  propensities[6] = B_prod1 *B_proddiff1 *Burst1_on_targ_allele1+B_prod1 *Burst1_off_targ_allele1;
  //update propensity for = B1 from allele2
  propensities[7] = B_prod1 *B_proddiff1 *Burst1_on_targ_allele1+B_prod1 *Burst1_off_targ_allele1;
  //update propensity for A1 = 
  propensities[8] = A_deg1 *A1;
  //update propensity for A1_nonsense = 
  propensities[9] = Anonsense_deg1 *A1_nonsense;
  //update propensity for Aprime1 = 
  propensities[10] = Aprime_deg1 *Aprime1;
  //update propensity for B1 = 
  propensities[11] = B_deg1 *B1;
  //update propensity for Burst1_off_orig = Burst1_on_orig for allele1
  propensities[12] = onbasal_a1 *Burst1_off_orig_allele1;
  //update propensity for Burst1_off_orig = Burst1_on_orig for allele2
  propensities[13] = onbasal_a1 *Burst1_off_orig_allele2;
  //update propensity for Burst1_off_para = Burst1_on_para for allele1
  propensities[14] = Aprimenitc1 *(pow(A1_nonsense,nA1_nonsense )/(pow(kA1_nonsense ,nA1_nonsense )+pow(A1_nonsense,nA1_nonsense )))*Burst1_off_para_allele1+onbasal_aprime1 *Burst1_off_para_allele1;
  //update propensity for Burst1_off_para = Burst1_on_para for allele2
  propensities[15] = Aprimenitc1 *(pow(A1_nonsense,nA1_nonsense )/(pow(kA1_nonsense ,nA1_nonsense )+pow(A1_nonsense,nA1_nonsense )))*Burst1_off_para_allele2+onbasal_aprime1 *Burst1_off_para_allele2;
  //update propensity for Burst1_off_targ = Burst1_on_targ for allele1
  propensities[16] = B_ondep1 *(pow(A1,nA1 )/(pow(kA1 ,nA1 )+pow(A1,nA1 )))*Burst1_off_targ_allele1+ B_ondep_prime *(pow(Aprime1,nAprime1 )/(pow(kAprime1 ,nAprime1 )+pow(Aprime1,nAprime1 )))*Burst1_off_targ_allele1 +onbasal_b1 *Burst1_off_targ_allele1;
  //update propensity for Burst1_off_targ = Burst1_on_targ for allele2
  propensities[17] = B_ondep1 *(pow(A1,nA1 )/(pow(kA1 ,nA1 )+pow(A1,nA1 )))*Burst1_off_targ_allele2+ B_ondep_prime *(pow(Aprime1,nAprime1 )/(pow(kAprime1 ,nAprime1 )+pow(Aprime1,nAprime1 )))*Burst1_off_targ_allele2 +onbasal_b1 *Burst1_off_targ_allele2;
  //update propensity for Burst1_on_orig = Burst1_off_orig for allele1
  propensities[18] = A_off1 *Burst1_on_orig_allele1;
  //update propensity for Burst1_on_orig = Burst1_off_orig for allele2
  propensities[19] = A_off1 *Burst1_on_orig_allele2;
  //update propensity for Burst1_on_para = Burst1_off_para for allele1
  propensities[20] = Aprime_off1 *Burst1_on_para_allele1;
  //update propensity for Burst1_on_para = Burst1_off_para for allele2
  propensities[21] = Aprime_off1 *Burst1_on_para_allele2;
  //update propensity for Burst1_on_targ = Burst1_off_targ for allele1
  propensities[22] = B_off1 *Burst1_on_targ_allele1;
  //update propensity for Burst1_on_targ = Burst1_off_targ for allele2
  propensities[23] = B_off1 *Burst1_on_targ_allele2;
} else if (p<cumpropensities[22]) {
  // rxn: Burst1_on_targ_allele1 = Burst1_off_targ_allele1

  if(Burst1_on_targ_allele1 == 1) {
    Burst1_off_targ_allele1 = 1;
    Burst1_on_targ_allele1 = 0;
  }
  //Burst1_off_targ=Burst1_off_targ + 1;
  //Burst1_on_targ=Burst1_on_targ + -1;

  //update propensity for = A1 from allele1
  propensities[0] = Burst1_not_mutated_allele1*(A_prod1 *A_proddiff1 *Burst1_on_orig_allele1+A_prod1 *Burst1_off_orig_allele1);
  //update propensity for = A1_nonsense from allele1
  propensities[1] = Burst1_is_mutated_allele1*(Anonsense_prod1 *A_proddiff1 *Burst1_on_orig_allele1+Anonsense_prod1 *Burst1_off_orig_allele1);
  //update propensity for = A1 from allele2
  propensities[2] = Burst1_not_mutated_allele2*(A_prod1 *A_proddiff1 *Burst1_on_orig_allele2+A_prod1 *Burst1_off_orig_allele2);
  //update propensity for = A1_nonsense from allele2
  propensities[3] = Burst1_is_mutated_allele2*(Anonsense_prod1 *A_proddiff1 *Burst1_on_orig_allele2+Anonsense_prod1 *Burst1_off_orig_allele2);
  //update propensity for = Aprime1 from allele1
  propensities[4] = Aprime_prod1 *Aprime_proddiff1 *Burst1_on_para_allele1+Aprime_prod1 *Burst1_off_para_allele1;
  //update propensity for = Aprime1 from allele2
  propensities[5] = Aprime_prod1 *Aprime_proddiff1 *Burst1_on_para_allele2+Aprime_prod1 *Burst1_off_para_allele2;
  //update propensity for = B1 from allele1
  propensities[6] = B_prod1 *B_proddiff1 *Burst1_on_targ_allele1+B_prod1 *Burst1_off_targ_allele1;
  //update propensity for = B1 from allele2
  propensities[7] = B_prod1 *B_proddiff1 *Burst1_on_targ_allele1+B_prod1 *Burst1_off_targ_allele1;
  //update propensity for A1 = 
  propensities[8] = A_deg1 *A1;
  //update propensity for A1_nonsense = 
  propensities[9] = Anonsense_deg1 *A1_nonsense;
  //update propensity for Aprime1 = 
  propensities[10] = Aprime_deg1 *Aprime1;
  //update propensity for B1 = 
  propensities[11] = B_deg1 *B1;
  //update propensity for Burst1_off_orig = Burst1_on_orig for allele1
  propensities[12] = onbasal_a1 *Burst1_off_orig_allele1;
  //update propensity for Burst1_off_orig = Burst1_on_orig for allele2
  propensities[13] = onbasal_a1 *Burst1_off_orig_allele2;
  //update propensity for Burst1_off_para = Burst1_on_para for allele1
  propensities[14] = Aprimenitc1 *(pow(A1_nonsense,nA1_nonsense )/(pow(kA1_nonsense ,nA1_nonsense )+pow(A1_nonsense,nA1_nonsense )))*Burst1_off_para_allele1+onbasal_aprime1 *Burst1_off_para_allele1;
  //update propensity for Burst1_off_para = Burst1_on_para for allele2
  propensities[15] = Aprimenitc1 *(pow(A1_nonsense,nA1_nonsense )/(pow(kA1_nonsense ,nA1_nonsense )+pow(A1_nonsense,nA1_nonsense )))*Burst1_off_para_allele2+onbasal_aprime1 *Burst1_off_para_allele2;
  //update propensity for Burst1_off_targ = Burst1_on_targ for allele1
  propensities[16] = B_ondep1 *(pow(A1,nA1 )/(pow(kA1 ,nA1 )+pow(A1,nA1 )))*Burst1_off_targ_allele1+ B_ondep_prime *(pow(Aprime1,nAprime1 )/(pow(kAprime1 ,nAprime1 )+pow(Aprime1,nAprime1 )))*Burst1_off_targ_allele1 +onbasal_b1 *Burst1_off_targ_allele1;
  //update propensity for Burst1_off_targ = Burst1_on_targ for allele2
  propensities[17] = B_ondep1 *(pow(A1,nA1 )/(pow(kA1 ,nA1 )+pow(A1,nA1 )))*Burst1_off_targ_allele2+ B_ondep_prime *(pow(Aprime1,nAprime1 )/(pow(kAprime1 ,nAprime1 )+pow(Aprime1,nAprime1 )))*Burst1_off_targ_allele2 +onbasal_b1 *Burst1_off_targ_allele2;
  //update propensity for Burst1_on_orig = Burst1_off_orig for allele1
  propensities[18] = A_off1 *Burst1_on_orig_allele1;
  //update propensity for Burst1_on_orig = Burst1_off_orig for allele2
  propensities[19] = A_off1 *Burst1_on_orig_allele2;
  //update propensity for Burst1_on_para = Burst1_off_para for allele1
  propensities[20] = Aprime_off1 *Burst1_on_para_allele1;
  //update propensity for Burst1_on_para = Burst1_off_para for allele2
  propensities[21] = Aprime_off1 *Burst1_on_para_allele2;
  //update propensity for Burst1_on_targ = Burst1_off_targ for allele1
  propensities[22] = B_off1 *Burst1_on_targ_allele1;
  //update propensity for Burst1_on_targ = Burst1_off_targ for allele2
  propensities[23] = B_off1 *Burst1_on_targ_allele2;
} else if (p<cumpropensities[23]) {
  // rxn: Burst1_on_targ_allele2 = Burst1_off_targ_allele2

  if(Burst1_on_targ_allele2 == 1) {
    Burst1_off_targ_allele2 = 1;
    Burst1_on_targ_allele2 = 0;
  }
  //Burst1_off_targ=Burst1_off_targ + 1;
  //Burst1_on_targ=Burst1_on_targ + -1;

  //update propensity for = A1 from allele1
  propensities[0] = Burst1_not_mutated_allele1*(A_prod1 *A_proddiff1 *Burst1_on_orig_allele1+A_prod1 *Burst1_off_orig_allele1);
  //update propensity for = A1_nonsense from allele1
  propensities[1] = Burst1_is_mutated_allele1*(Anonsense_prod1 *A_proddiff1 *Burst1_on_orig_allele1+Anonsense_prod1 *Burst1_off_orig_allele1);
  //update propensity for = A1 from allele2
  propensities[2] = Burst1_not_mutated_allele2*(A_prod1 *A_proddiff1 *Burst1_on_orig_allele2+A_prod1 *Burst1_off_orig_allele2);
  //update propensity for = A1_nonsense from allele2
  propensities[3] = Burst1_is_mutated_allele2*(Anonsense_prod1 *A_proddiff1 *Burst1_on_orig_allele2+Anonsense_prod1 *Burst1_off_orig_allele2);
  //update propensity for = Aprime1 from allele1
  propensities[4] = Aprime_prod1 *Aprime_proddiff1 *Burst1_on_para_allele1+Aprime_prod1 *Burst1_off_para_allele1;
  //update propensity for = Aprime1 from allele2
  propensities[5] = Aprime_prod1 *Aprime_proddiff1 *Burst1_on_para_allele2+Aprime_prod1 *Burst1_off_para_allele2;
  //update propensity for = B1 from allele1
  propensities[6] = B_prod1 *B_proddiff1 *Burst1_on_targ_allele1+B_prod1 *Burst1_off_targ_allele1;
  //update propensity for = B1 from allele2
  propensities[7] = B_prod1 *B_proddiff1 *Burst1_on_targ_allele1+B_prod1 *Burst1_off_targ_allele1;
  //update propensity for A1 = 
  propensities[8] = A_deg1 *A1;
  //update propensity for A1_nonsense = 
  propensities[9] = Anonsense_deg1 *A1_nonsense;
  //update propensity for Aprime1 = 
  propensities[10] = Aprime_deg1 *Aprime1;
  //update propensity for B1 = 
  propensities[11] = B_deg1 *B1;
  //update propensity for Burst1_off_orig = Burst1_on_orig for allele1
  propensities[12] = onbasal_a1 *Burst1_off_orig_allele1;
  //update propensity for Burst1_off_orig = Burst1_on_orig for allele2
  propensities[13] = onbasal_a1 *Burst1_off_orig_allele2;
  //update propensity for Burst1_off_para = Burst1_on_para for allele1
  propensities[14] = Aprimenitc1 *(pow(A1_nonsense,nA1_nonsense )/(pow(kA1_nonsense ,nA1_nonsense )+pow(A1_nonsense,nA1_nonsense )))*Burst1_off_para_allele1+onbasal_aprime1 *Burst1_off_para_allele1;
  //update propensity for Burst1_off_para = Burst1_on_para for allele2
  propensities[15] = Aprimenitc1 *(pow(A1_nonsense,nA1_nonsense )/(pow(kA1_nonsense ,nA1_nonsense )+pow(A1_nonsense,nA1_nonsense )))*Burst1_off_para_allele2+onbasal_aprime1 *Burst1_off_para_allele2;
  //update propensity for Burst1_off_targ = Burst1_on_targ for allele1
  propensities[16] = B_ondep1 *(pow(A1,nA1 )/(pow(kA1 ,nA1 )+pow(A1,nA1 )))*Burst1_off_targ_allele1+ B_ondep_prime *(pow(Aprime1,nAprime1 )/(pow(kAprime1 ,nAprime1 )+pow(Aprime1,nAprime1 )))*Burst1_off_targ_allele1 +onbasal_b1 *Burst1_off_targ_allele1;
  //update propensity for Burst1_off_targ = Burst1_on_targ for allele2
  propensities[17] = B_ondep1 *(pow(A1,nA1 )/(pow(kA1 ,nA1 )+pow(A1,nA1 )))*Burst1_off_targ_allele2+ B_ondep_prime *(pow(Aprime1,nAprime1 )/(pow(kAprime1 ,nAprime1 )+pow(Aprime1,nAprime1 )))*Burst1_off_targ_allele2 +onbasal_b1 *Burst1_off_targ_allele2;
  //update propensity for Burst1_on_orig = Burst1_off_orig for allele1
  propensities[18] = A_off1 *Burst1_on_orig_allele1;
  //update propensity for Burst1_on_orig = Burst1_off_orig for allele2
  propensities[19] = A_off1 *Burst1_on_orig_allele2;
  //update propensity for Burst1_on_para = Burst1_off_para for allele1
  propensities[20] = Aprime_off1 *Burst1_on_para_allele1;
  //update propensity for Burst1_on_para = Burst1_off_para for allele2
  propensities[21] = Aprime_off1 *Burst1_on_para_allele2;
  //update propensity for Burst1_on_targ = Burst1_off_targ for allele1
  propensities[22] = B_off1 *Burst1_on_targ_allele1;
  //update propensity for Burst1_on_targ = Burst1_off_targ for allele2
  propensities[23] = B_off1 *Burst1_on_targ_allele2;
}
  //}
  }

  printf("Total iterations = %d\n",totaliterations);
  //printf("savecount = %d\n",savecount);
  //printf("savecountcheck = %d\n",savecountcheck);
}
/*-------------------------------------------------------------------*/

void mexFunction(int nlhs, mxArray *plhs[], int nrhs,
                 const mxArray *prhs[])
{
  long i;

  double *rand1, *rand2, *propensities_in, *species, *rates;
  double propensities[NUMRXNS];
  double *currTpr;
  double *seedpr, *mpr, *totaltpr;
  double seed, m, totalt;

  double *times_out, *species_out;

  double currT;

  //Load scalars first...

  currTpr = mxGetPr(CURRT_IN);
  currT = currTpr[0];

  seedpr = mxGetPr(SEED_IN);
  seed = seedpr[0];
  mpr = mxGetPr(M_IN);
  m = mpr[0];

  totaltpr = mxGetPr(TOTALT_IN);
  totalt = totaltpr[0];

  //Now load matrices/vectors...
  species = mxGetPr(SPECIES_IN);
  propensities_in = mxGetPr(PROPENSITIES_IN);
  rates = mxGetPr(RATES_IN);

  //Now let's create the output matrices
  TIMES_OUT = mxCreateDoubleMatrix(1,(long) m,mxREAL);
  times_out = mxGetPr(TIMES_OUT);
  SPECIES_OUT = mxCreateDoubleMatrix(NSPECIES,(long) m,mxREAL);
  species_out = mxGetPr(SPECIES_OUT);

  // Copy propensities_in into propensities
  for (i = 0; i<NUMRXNS; i++){
        propensities[i] = propensities_in[i];
  }

  gillespie((long)m,times_out,species_out,currT,species,rates,propensities,(long)seed, totalt);

} /* end function mexFunction */

/*-------------------------------------------------------------------*/
// Random number stuff below


/* nfix() generates variates from the residue when rejection in RNOR occurs. */

double nfix(void)
{
const double r = 3.442620f;     /* The start of the right tail */
static double x, y;
 for(;;)
  {  x=hz*wn[iz];      /* iz==0, handles the base strip */
     if(iz==0)
       { do{ x=-log(UNI)*0.2904764; y=-log(UNI);}	/* .2904764 is 1/r */
        while(y+y<x*x);
        return (hz>0)? r+x : -r-x;
       }
                         /* iz>0, handle the wedges of other strips */
      if( fn[iz]+UNI*(fn[iz-1]-fn[iz]) < exp(-.5*x*x) ) return x;

     /* initiate, try to exit for(;;) for loop*/
      hz=SHR3;
      iz=hz&127;
      if(fabs(hz)<kn[iz]) return (hz*wn[iz]);
  }

}

/* efix() generates variates from the residue when rejection in REXP occurs. */
double efix(void)
{ double x;
 for(;;)
  {  if(iz==0) return (7.69711-log(UNI));          /* iz==0 */
     x=jz*we[iz]; if( fe[iz]+UNI*(fe[iz-1]-fe[iz]) < exp(-x) ) return (x);

      /* initiate, try to exit for(;;) loop */
   jz=SHR3;
   iz=(jz&255);
   if(jz<ke[iz]) return (jz*we[iz]);
  }
}
/*--------This procedure sets the seed and creates the tables------*/

void zigset(unsigned long jsrseed)
{  const double m1 = 2147483648.0, m2 = 4294967296.;
   double dn=3.442619855899,tn=dn,vn=9.91256303526217e-3, q;
   double de=7.697117470131487, te=de, ve=3.949659822581572e-3;
   int i;
   jsr^=jsrseed;

/* Set up tables for RNOR */
   q=vn/exp(-.5*dn*dn);
   kn[0]=(dn/q)*m1;
   kn[1]=0;

   wn[0]=q/m1;
   wn[127]=dn/m1;

   fn[0]=1.;
   fn[127]=exp(-.5*dn*dn);

    for(i=126;i>=1;i--)
    {dn=sqrt(-2.*log(vn/dn+exp(-.5*dn*dn)));
     kn[i+1]=(dn/tn)*m1;
     tn=dn;
     fn[i]=exp(-.5*dn*dn);
     wn[i]=dn/m1;
    }

/* Set up tables for REXP */
    q = ve/exp(-de);
    ke[0]=(de/q)*m2;
    ke[1]=0;

    we[0]=q/m2;
    we[255]=de/m2;

    fe[0]=1.;
    fe[255]=exp(-de);

   for(i=254;i>=1;i--)
  {de=-log(ve/de+exp(-de));
   ke[i+1]= (de/te)*m2;
   te=de;
   fe[i]=exp(-de);
   we[i]=de/m2;
  }
}
double ran_expo(double lambda){
    double u;
    u = rand() / (RAND_MAX + 1.0);
    return -log(1- u) / lambda;
}
