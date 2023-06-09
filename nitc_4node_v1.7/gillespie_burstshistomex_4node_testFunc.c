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
#define NSPECIES 21

//INCLUDE NUMRXNS HERE
#define NUMRXNS 39

void updateSpecies(double *species, long A1, long Anonsense1, long Aprime1, long Aprime2, long B1, long Burst1_onorig_targ_allele1, long Burst1_onpara1_targ_allele1,  long Burst1_onpara2_targ_allele1,long Burst1_off_targ_allele1, long Burst1_onorig_targ_allele2, long Burst1_onpara1_targ_allele2, long Burst1_onpara2_targ_allele2, long Burst1_off_targ_allele2, long Burst1_on_para1_allele1, long Burst1_on_para1_allele2, long Burst1_on_para2_allele1, long Burst1_on_para2_allele2, long Burst1_on_orig_allele1, long Burst1_on_orig_allele2, long Burst1_is_mutated_allele1, long Burst1_is_mutated_allele2)
{
    species[0] = A1;
    species[1] = Anonsense1;
    species[2] = Aprime1;
    species[3] = Aprime2;
    species[4] = B1;

    species[5] = Burst1_onorig_targ_allele1;
    species[6] = Burst1_onpara1_targ_allele1;
    species[7] = Burst1_onpara2_targ_allele1;
    species[8] = Burst1_off_targ_allele1;
    species[9] = Burst1_onorig_targ_allele2;
    species[10] = Burst1_onpara1_targ_allele2;
    species[11] = Burst1_onpara2_targ_allele2;
    species[12] = Burst1_off_targ_allele2;
    species[13] = Burst1_on_para1_allele1;
    species[14] = Burst1_on_para1_allele2;
    species[15] = Burst1_on_para2_allele1;
    species[16] = Burst1_on_para2_allele2;
    species[17] = Burst1_on_orig_allele1;
    species[18] = Burst1_on_orig_allele2;
    species[19] = Burst1_is_mutated_allele1;
    species[20] = Burst1_is_mutated_allele2;

     // for (int i = 0; i<16; i++)
     //  printf(" species %i %f \n", i, species[i]);
}

void updatePropensities(double *species, double *rates,double *propensities)
{
  long A1, Anonsense1, Aprime1, Aprime2, B1, Burst1_onorig_targ_allele1, Burst1_onpara1_targ_allele1, Burst1_onpara2_targ_allele1,  Burst1_off_targ_allele1, Burst1_onorig_targ_allele2, Burst1_onpara1_targ_allele2, Burst1_onpara2_targ_allele2, Burst1_off_targ_allele2, Burst1_on_para1_allele1, Burst1_on_para1_allele2, Burst1_on_para2_allele1, Burst1_on_para2_allele2, Burst1_on_orig_allele1, Burst1_on_orig_allele2, Burst1_is_mutated_allele1, Burst1_is_mutated_allele2;
  double r_prodbasal_A1, r_prodbasal_Anonsense1, r_prodbasal_Aprime1, r_prodbasal_Aprime2, r_prodbasal_B1, r_prodon_A1, r_prodon_Anonsense1, r_prodon_Aprime1, r_prodon_Aprime2, r_prodon_B1, d_Aprime1_B1, d_Aprime2_B1, r_deg_A1, r_deg_Anonsense1, r_deg_Aprime1, r_deg_Aprime2, r_deg_B1, r_onbasal_A1, r_onbasal_Anonsense1, r_onbasal_Aprime1, r_onbasal_Aprime2, r_onbasal_B1, r_nitc_byAnonsense1_A1, r_nitc_byAnonsense1_Anonsense1, r_nitc_byAnonsense1_Aprime1, r_nitc_byAnonsense1_Aprime2, r_addon_byA1_B1, r_addon_byAprime1_B1, r_addon_byAprime2_B1, r_off_A1, r_off_Anonsense1, r_off_Aprime1, r_off_Aprime2, r_offorig_B1, r_offpara1_B1, r_offpara2_B1, k_A1, k_Anonsense1, k_Aprime1, k_Aprime2, k_B1, n_A1, n_Anonsense1, n_Aprime1, n_Aprime2, n_B1;

  A1 = (long)species[0];
  Anonsense1 = (long)species[1];
  Aprime1 = (long)species[2];
  Aprime2 = (long)species[3];
  B1 = (long)species[4];

  Burst1_onorig_targ_allele1 = (long)species[5];
  Burst1_onpara1_targ_allele1 = (long)species[6];
  Burst1_onpara2_targ_allele1 = (long)species[7];
  Burst1_off_targ_allele1 = (long)species[8];
  Burst1_onorig_targ_allele2 = (long)species[9];
  Burst1_onpara1_targ_allele2 = (long)species[10];
  Burst1_onpara2_targ_allele2 = (long)species[11];
  Burst1_off_targ_allele2 = (long)species[12];
  Burst1_on_para1_allele1 = (long)species[13];
  Burst1_on_para1_allele2 = (long)species[14];
  Burst1_on_para2_allele1 = (long)species[15];
  Burst1_on_para2_allele2 = (long)species[16];
  Burst1_on_orig_allele1 = (long)species[17];
  Burst1_on_orig_allele2 = (long)species[18];
  Burst1_is_mutated_allele1 = (long)species[19];
  Burst1_is_mutated_allele2 = (long)species[20];

//UNPACK ALL RATES HERE
  r_prodbasal_A1 = rates[0];
  r_prodbasal_Anonsense1 = rates[1];
  r_prodbasal_Aprime1 = rates[2];
  r_prodbasal_Aprime2 = rates[3];
  r_prodbasal_B1 = rates[4];
  
  r_prodon_A1 = rates[5];
  r_prodon_Anonsense1 = rates[6];
  r_prodon_Aprime1 = rates[7];
  r_prodon_Aprime2 = rates[8];
  r_prodon_B1 = rates[9];

  d_Aprime1_B1 = rates[10];
  d_Aprime2_B1 = rates[11];

  r_deg_A1 = rates[12];
  r_deg_Anonsense1 = rates[13];
  r_deg_Aprime1 = rates[14];
  r_deg_Aprime2 = rates[15];
  r_deg_B1 = rates[16];

  r_onbasal_A1 = rates[17];
  r_onbasal_Anonsense1 = rates[18];
  r_onbasal_Aprime1 = rates[19];
  r_onbasal_Aprime2 = rates[20];
  r_onbasal_B1 = rates[21];

  r_nitc_byAnonsense1_A1 = rates[22];
  r_nitc_byAnonsense1_Anonsense1 = rates[23];
  r_nitc_byAnonsense1_Aprime1 = rates[24];
  r_nitc_byAnonsense1_Aprime2 = rates[25];
  r_addon_byA1_B1 = rates[26];
  r_addon_byAprime1_B1 = rates[27];
  r_addon_byAprime2_B1 = rates[28];

  r_off_A1 = rates[29];
  r_off_Anonsense1 = rates[30];
  r_off_Aprime1 = rates[31];
  r_off_Aprime2 = rates[32];
  r_offorig_B1 = rates[33];
  r_offpara1_B1 = rates[34];
  r_offpara2_B1 = rates[35];

  k_A1  = rates[36];
  k_Anonsense1  = rates[37];
  k_Aprime1  = rates[38];
  k_Aprime2  = rates[39];
  k_B1  = rates[40];

  n_A1  = rates[41];
  n_Anonsense1  = rates[42];
  n_Aprime1  = rates[43];
  n_Aprime2  = rates[44];
  n_B1  = rates[45];


  //update propensity for = A1 from allele1
  propensities[0] = (1-Burst1_is_mutated_allele1)*(r_prodon_A1 * Burst1_on_orig_allele1 + r_prodbasal_A1 * (1-Burst1_on_orig_allele1));
  //update propensity for = A1_nonsense from allele1
  propensities[1] = Burst1_is_mutated_allele1*(r_prodon_Anonsense1 * Burst1_on_orig_allele1 + r_prodbasal_Anonsense1 * (1-Burst1_on_orig_allele1));
  //update propensity for = A1 from allele2
  propensities[2] = (1-Burst1_is_mutated_allele2)*(r_prodon_A1 * Burst1_on_orig_allele2 + r_prodbasal_A1 * (1-Burst1_on_orig_allele2));
  //update propensity for = A1_nonsense from allele2
  propensities[3] = Burst1_is_mutated_allele2*(r_prodon_Anonsense1 * Burst1_on_orig_allele2 + r_prodbasal_Anonsense1 * (1-Burst1_on_orig_allele2));
  //update propensity for = Aprime1 from allele1
  propensities[4] = r_prodon_Aprime1 *Burst1_on_para1_allele1 + r_prodbasal_Aprime1 * (1-Burst1_on_para1_allele1);
  //update propensity for = Aprime1 from allele2
  propensities[5] = r_prodon_Aprime1 *Burst1_on_para1_allele2 + r_prodbasal_Aprime1 * (1-Burst1_on_para1_allele2);
  //update propensity for = Aprime2 from allele1
  propensities[6] = r_prodon_Aprime2 *Burst1_on_para2_allele1 + r_prodbasal_Aprime2 * (1-Burst1_on_para2_allele1);
  //update propensity for = Aprime2 from allele2
  propensities[7] = r_prodon_Aprime2 *Burst1_on_para2_allele2 + r_prodbasal_Aprime2 * (1-Burst1_on_para2_allele2);
  //update propensity for = B1 from allele1
  propensities[8] = r_prodon_B1 *Burst1_onorig_targ_allele1 + d_Aprime1_B1*r_prodon_B1*Burst1_onpara1_targ_allele1 + d_Aprime2_B1*r_prodon_B1*Burst1_onpara2_targ_allele1 + r_prodbasal_B1 *Burst1_off_targ_allele1;
  //update propensity for = B1 from allele2
  propensities[9] = r_prodon_B1 *Burst1_onorig_targ_allele2 + d_Aprime1_B1*r_prodon_B1*Burst1_onpara1_targ_allele2+ d_Aprime2_B1*r_prodon_B1*Burst1_onpara2_targ_allele2 + r_prodbasal_B1 *Burst1_off_targ_allele2;

  //update propensity for A1 = 
  propensities[10] = r_deg_A1 *A1;
  //update propensity for A1_nonsense = 
  propensities[11] = r_deg_Anonsense1 *Anonsense1;
  //update propensity for Aprime1 = 
  propensities[12] = r_deg_Aprime1 *Aprime1;
  //update propensity for Aprime2 = 
  propensities[13] = r_deg_Aprime2 *Aprime2;
  //update propensity for B1 = 
  propensities[14] = r_deg_B1 *B1;

  //update propensity for Burst1_off_orig = Burst1_on_orig for allele1
  propensities[15] = (1-Burst1_on_orig_allele1)*((1-Burst1_is_mutated_allele1)*(r_onbasal_A1 + r_nitc_byAnonsense1_A1*(pow(Anonsense1,n_Anonsense1 )/(pow(k_Anonsense1 ,n_Anonsense1 )+pow(Anonsense1,n_Anonsense1 )))) + Burst1_is_mutated_allele1*(r_onbasal_Anonsense1 + r_nitc_byAnonsense1_Anonsense1*(pow(Anonsense1,n_Anonsense1 )/(pow(k_Anonsense1 ,n_Anonsense1 )+pow(Anonsense1,n_Anonsense1 )))));
  //update propensity for Burst1_off_orig = Burst1_on_orig for allele2
  propensities[16] = (1-Burst1_on_orig_allele2)*((1-Burst1_is_mutated_allele2)*(r_onbasal_A1 + r_nitc_byAnonsense1_A1*(pow(Anonsense1,n_Anonsense1 )/(pow(k_Anonsense1 ,n_Anonsense1 )+pow(Anonsense1,n_Anonsense1 )))) + Burst1_is_mutated_allele2*(r_onbasal_Anonsense1 + r_nitc_byAnonsense1_Anonsense1*(pow(Anonsense1,n_Anonsense1 )/(pow(k_Anonsense1 ,n_Anonsense1 )+pow(Anonsense1,n_Anonsense1 )))));
  //update propensity for Burst1_off_para1 = Burst1_on_para1 for allele1
  propensities[17] = (1-Burst1_on_para1_allele1)*(r_onbasal_Aprime1 + r_nitc_byAnonsense1_Aprime1*(pow(Anonsense1,n_Anonsense1 )/(pow(k_Anonsense1 ,n_Anonsense1 )+pow(Anonsense1,n_Anonsense1 ))));
  //update propensity for Burst1_off_para1 = Burst1_on_para1 for allele2
  propensities[18] = (1-Burst1_on_para1_allele2)*(r_onbasal_Aprime1 + r_nitc_byAnonsense1_Aprime1*(pow(Anonsense1,n_Anonsense1 )/(pow(k_Anonsense1 ,n_Anonsense1 )+pow(Anonsense1,n_Anonsense1 ))));
  //update propensity for Burst1_off_para2 = Burst1_on_para2 for allele1
  propensities[19] = (1-Burst1_on_para2_allele1)*(r_onbasal_Aprime2 + r_nitc_byAnonsense1_Aprime2*(pow(Anonsense1,n_Anonsense1 )/(pow(k_Anonsense1 ,n_Anonsense1 )+pow(Anonsense1,n_Anonsense1 ))));
  //update propensity for Burst1_off_para2 = Burst1_on_para2 for allele2
  propensities[20] = (1-Burst1_on_para2_allele2)*(r_onbasal_Aprime2 + r_nitc_byAnonsense1_Aprime2*(pow(Anonsense1,n_Anonsense1 )/(pow(k_Anonsense1 ,n_Anonsense1 )+pow(Anonsense1,n_Anonsense1 ))));
  //update propensity for Burst1_off_targ = Burst1_onorig_targ for allele1
  propensities[21] = (Burst1_off_targ_allele1)*(r_onbasal_B1 + r_addon_byA1_B1*(pow(A1,n_A1 )/(pow(k_A1 ,n_A1 )+pow(A1,n_A1 ))));
  //update propensity for Burst1_off_targ = Burst1_onpara1_targ for allele1
  propensities[22] = (Burst1_off_targ_allele1)*(r_onbasal_B1 + r_addon_byAprime1_B1*(pow(Aprime1,n_Aprime1 )/(pow(k_Aprime1 ,n_Aprime1 )+pow(Aprime1,n_Aprime1 ))));
  //update propensity for Burst1_off_targ = Burst1_onpara2_targ for allele1
  propensities[23] = (Burst1_off_targ_allele1)*(r_onbasal_B1 + r_addon_byAprime2_B1*(pow(Aprime2,n_Aprime2 )/(pow(k_Aprime2 ,n_Aprime2 )+pow(Aprime2,n_Aprime2 ))));
  //update propensity for Burst1_off_targ = Burst1_onorig_targ for allele2
  propensities[24] = (Burst1_off_targ_allele2)*(r_onbasal_B1 + r_addon_byA1_B1*(pow(A1,n_A1 )/(pow(k_A1 ,n_A1 )+pow(A1,n_A1 ))));
  //update propensity for Burst1_off_targ = Burst1_onpara1_targ for allele2
  propensities[25] = (Burst1_off_targ_allele2)*(r_onbasal_B1 + r_addon_byAprime1_B1*(pow(Aprime1,n_Aprime1 )/(pow(k_Aprime1 ,n_Aprime1 )+pow(Aprime1,n_Aprime1 ))));
  //update propensity for Burst1_off_targ = Burst1_onpara2_targ for allele2
  propensities[26] = (Burst1_off_targ_allele2)*(r_onbasal_B1 + r_addon_byAprime2_B1*(pow(Aprime2,n_Aprime2 )/(pow(k_Aprime2 ,n_Aprime2 )+pow(Aprime2,n_Aprime2 ))));


  //update propensity for Burst1_on_orig = Burst1_off_orig for allele1
  propensities[27] = Burst1_on_orig_allele1*(Burst1_is_mutated_allele1*r_off_Anonsense1 + (1-Burst1_is_mutated_allele1)*r_off_A1);
  //update propensity for Burst1_on_orig = Burst1_off_orig for allele2
  propensities[28] = Burst1_on_orig_allele2*(Burst1_is_mutated_allele2*r_off_Anonsense1 + (1-Burst1_is_mutated_allele2)*r_off_A1);
  //update propensity for Burst1_on_para1 = Burst1_off_para1 for allele1
  propensities[29] = Burst1_on_para1_allele1*r_off_Aprime1;
  //update propensity for Burst1_on_para1 = Burst1_off_para1 for allele2
  propensities[30] = Burst1_on_para1_allele2*r_off_Aprime1;
  //update propensity for Burst1_on_para2 = Burst1_off_para2 for allele1
  propensities[31] = Burst1_on_para2_allele1*r_off_Aprime2;
  //update propensity for Burst1_on_para2 = Burst1_off_para2 for allele2
  propensities[32] = Burst1_on_para2_allele2*r_off_Aprime2;
  //update propensity for Burst1_onorig_targ = Burst1_off_targ for allele1
  propensities[33] = Burst1_onorig_targ_allele1*r_offorig_B1;
  //update propensity for Burst1_onpara1_targ = Burst1_off_targ for allele1
  propensities[34] = Burst1_onpara1_targ_allele1*r_offpara1_B1;
  //update propensity for Burst1_onpara2_targ = Burst1_off_targ for allele1
  propensities[35] = Burst1_onpara2_targ_allele1*r_offpara2_B1;
  //update propensity for Burst1_onorig_targ = Burst1_off_targ for allele2
  propensities[36] = Burst1_onorig_targ_allele2*r_offorig_B1;
  //update propensity for Burst1_onpara1_targ = Burst1_off_targ for allele2
  propensities[37] = Burst1_onpara1_targ_allele2*r_offpara1_B1;
  //update propensity for Burst1_onpara2_targ = Burst1_off_targ for allele2
  propensities[38] = Burst1_onpara2_targ_allele2*r_offpara2_B1;


  species[0] = A1;
  species[1] = Anonsense1;
  species[2] = Aprime1;
  species[3] = Aprime2;
  species[4] = B1;

  species[5] = Burst1_onorig_targ_allele1;
  species[6] = Burst1_onpara1_targ_allele1;
  species[7] = Burst1_onpara2_targ_allele1;
  species[8] = Burst1_off_targ_allele1;
  species[9] = Burst1_onorig_targ_allele2;
  species[10] = Burst1_onpara1_targ_allele2;
  species[11] = Burst1_onpara2_targ_allele2;
  species[12] = Burst1_off_targ_allele2;
  species[13] = Burst1_on_para1_allele1;
  species[14] = Burst1_on_para1_allele2;
  species[15] = Burst1_on_para2_allele1;
  species[16] = Burst1_on_para2_allele2;
  species[17] = Burst1_on_orig_allele1;
  species[18] = Burst1_on_orig_allele2;
  species[19] = Burst1_is_mutated_allele1;
  species[20] = Burst1_is_mutated_allele2;
  // updateSpecies(species, A1, Anonsense1, Aprime1, B1, Burst1_onorig_targ_allele1, Burst1_onpara_targ_allele1, Burst1_off_targ_allele1, Burst1_onorig_targ_allele2, Burst1_onpara_targ_allele2, Burst1_off_targ_allele2, Burst1_on_para_allele1, Burst1_on_para_allele2, Burst1_on_orig_allele1, Burst1_on_orig_allele2, Burst1_is_mutated_allele1, Burst1_is_mutated_allele2);

  // for (int i = 0; i<28; i++)
  //     printf(" prop %i %f \n", i, propensities[i]);
  // need to return both species and propensities
}

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
  long A1, Anonsense1, Aprime1, Aprime2, B1, Burst1_onorig_targ_allele1, Burst1_onpara1_targ_allele1, Burst1_onpara2_targ_allele1,  Burst1_off_targ_allele1, Burst1_onorig_targ_allele2, Burst1_onpara1_targ_allele2, Burst1_onpara2_targ_allele2, Burst1_off_targ_allele2, Burst1_on_para1_allele1, Burst1_on_para1_allele2, Burst1_on_para2_allele1, Burst1_on_para2_allele2, Burst1_on_orig_allele1, Burst1_on_orig_allele2, Burst1_is_mutated_allele1, Burst1_is_mutated_allele2;
  double r_prodbasal_A1, r_prodbasal_Anonsense1, r_prodbasal_Aprime1, r_prodbasal_Aprime2, r_prodbasal_B1, r_prodon_A1, r_prodon_Anonsense1, r_prodon_Aprime1, r_prodon_Aprime2, r_prodon_B1, d_Aprime1_B1, d_Aprime2_B1, r_deg_A1, r_deg_Anonsense1, r_deg_Aprime1, r_deg_Aprime2, r_deg_B1, r_onbasal_A1, r_onbasal_Anonsense1, r_onbasal_Aprime1, r_onbasal_Aprime2, r_onbasal_B1, r_nitc_byAnonsense1_A1, r_nitc_byAnonsense1_Anonsense1, r_nitc_byAnonsense1_Aprime1, r_nitc_byAnonsense1_Aprime2, r_addon_byA1_B1, r_addon_byAprime1_B1, r_addon_byAprime2_B1, r_off_A1, r_off_Anonsense1, r_off_Aprime1, r_off_Aprime2, r_offorig_B1, r_offpara1_B1, r_offpara2_B1, k_A1, k_Anonsense1, k_Aprime1, k_Aprime2, k_B1, n_A1, n_Anonsense1, n_Aprime1, n_Aprime2, n_B1;
long i, j, k;
//UNPACK ALL SPECIES HERE
  A1 = (long)species[0];
  Anonsense1 = (long)species[1];
  Aprime1 = (long)species[2];
  Aprime2 = (long)species[3];
  B1 = (long)species[4];

  Burst1_onorig_targ_allele1 = (long)species[5];
  Burst1_onpara1_targ_allele1 = (long)species[6];
  Burst1_onpara2_targ_allele1 = (long)species[7];
  Burst1_off_targ_allele1 = (long)species[8];
  Burst1_onorig_targ_allele2 = (long)species[9];
  Burst1_onpara1_targ_allele2 = (long)species[10];
  Burst1_onpara2_targ_allele2 = (long)species[11];
  Burst1_off_targ_allele2 = (long)species[12];
  Burst1_on_para1_allele1 = (long)species[13];
  Burst1_on_para1_allele2 = (long)species[14];
  Burst1_on_para2_allele1 = (long)species[15];
  Burst1_on_para2_allele2 = (long)species[16];
  Burst1_on_orig_allele1 = (long)species[17];
  Burst1_on_orig_allele2 = (long)species[18];
  Burst1_is_mutated_allele1 = (long)species[19];
  Burst1_is_mutated_allele2 = (long)species[20];

//UNPACK ALL RATES HERE
  r_prodbasal_A1 = rates[0];
  r_prodbasal_Anonsense1 = rates[1];
  r_prodbasal_Aprime1 = rates[2];
  r_prodbasal_Aprime2 = rates[3];
  r_prodbasal_B1 = rates[4];
  
  r_prodon_A1 = rates[5];
  r_prodon_Anonsense1 = rates[6];
  r_prodon_Aprime1 = rates[7];
  r_prodon_Aprime2 = rates[8];
  r_prodon_B1 = rates[9];

  d_Aprime1_B1 = rates[10];
  d_Aprime2_B1 = rates[11];

  r_deg_A1 = rates[12];
  r_deg_Anonsense1 = rates[13];
  r_deg_Aprime1 = rates[14];
  r_deg_Aprime2 = rates[15];
  r_deg_B1 = rates[16];

  r_onbasal_A1 = rates[17];
  r_onbasal_Anonsense1 = rates[18];
  r_onbasal_Aprime1 = rates[19];
  r_onbasal_Aprime2 = rates[20];
  r_onbasal_B1 = rates[21];

  r_nitc_byAnonsense1_A1 = rates[22];
  r_nitc_byAnonsense1_Anonsense1 = rates[23];
  r_nitc_byAnonsense1_Aprime1 = rates[24];
  r_nitc_byAnonsense1_Aprime2 = rates[25];
  r_addon_byA1_B1 = rates[26];
  r_addon_byAprime1_B1 = rates[27];
  r_addon_byAprime2_B1 = rates[28];

  r_off_A1 = rates[29];
  r_off_Anonsense1 = rates[30];
  r_off_Aprime1 = rates[31];
  r_off_Aprime2 = rates[32];
  r_offorig_B1 = rates[33];
  r_offpara1_B1 = rates[34];
  r_offpara2_B1 = rates[35];

  k_A1  = rates[36];
  k_Anonsense1  = rates[37];
  k_Aprime1  = rates[38];
  k_Aprime2  = rates[39];
  k_B1  = rates[40];

  n_A1  = rates[41];
  n_Anonsense1  = rates[42];
  n_Aprime1  = rates[43];
  n_Aprime2  = rates[44];
  n_B1  = rates[45];

 

  zigset(seed);
  totaliterations = 0;

  savecount = 0;
  savecountcheck = 0;
  deltaTsave = totalt/(m-1);


  //  for (i=0; i<m; i++) {
  while (currT<totalt) {
    if (currT > 100000) {
      Burst1_is_mutated_allele1 = 1;
      Burst1_is_mutated_allele2 = 0;
    }
    if (currT > 200000) {
      Burst1_is_mutated_allele1 = 1;
      Burst1_is_mutated_allele2 = 1;
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
    // printf("alpha %f \n", alpha);
    // printf("deltaT %f \n", deltaT);
    // printf("p %f \n", p);
    // printf("total cumulativeprops %f \n", cumpropensities[27]);

    // if(totaliterations == 5) {
    //   break;
    // }

//if(savecount < 1000){
    while (currT > savecount*deltaTsave) {
        if (savecount < m){
//INSERT SAVE HERE
species_out[savecount*NSPECIES+0] = A1;
species_out[savecount*NSPECIES+1] = Anonsense1;
species_out[savecount*NSPECIES+2] = Aprime1;
species_out[savecount*NSPECIES+3] = Aprime2;
species_out[savecount*NSPECIES+4] = B1;

species_out[savecount*NSPECIES+5] = Burst1_onorig_targ_allele1;
species_out[savecount*NSPECIES+6] = Burst1_onpara1_targ_allele1;
species_out[savecount*NSPECIES+7] = Burst1_onpara2_targ_allele1;
species_out[savecount*NSPECIES+8] = Burst1_off_targ_allele1;
species_out[savecount*NSPECIES+9] = Burst1_onorig_targ_allele2;
species_out[savecount*NSPECIES+10] = Burst1_onpara1_targ_allele2;
species_out[savecount*NSPECIES+11] = Burst1_onpara2_targ_allele2;
species_out[savecount*NSPECIES+12] = Burst1_off_targ_allele2;
species_out[savecount*NSPECIES+13] = Burst1_on_para1_allele1;
species_out[savecount*NSPECIES+14] = Burst1_on_para1_allele2;
species_out[savecount*NSPECIES+15] = Burst1_on_para2_allele1;
species_out[savecount*NSPECIES+16] = Burst1_on_para2_allele2;
species_out[savecount*NSPECIES+17] = Burst1_on_orig_allele1;
species_out[savecount*NSPECIES+18] = Burst1_on_orig_allele2;
species_out[savecount*NSPECIES+19] = Burst1_is_mutated_allele1;
species_out[savecount*NSPECIES+20] = Burst1_is_mutated_allele2;
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
 
  updateSpecies(species, A1, Anonsense1, Aprime1, Aprime2, B1, Burst1_onorig_targ_allele1, Burst1_onpara1_targ_allele1, Burst1_onpara2_targ_allele1,  Burst1_off_targ_allele1, Burst1_onorig_targ_allele2, Burst1_onpara1_targ_allele2, Burst1_onpara2_targ_allele2, Burst1_off_targ_allele2, Burst1_on_para1_allele1, Burst1_on_para1_allele2, Burst1_on_para2_allele1, Burst1_on_para2_allele2, Burst1_on_orig_allele1, Burst1_on_orig_allele2, Burst1_is_mutated_allele1, Burst1_is_mutated_allele2);
  updatePropensities(species, rates, propensities);

} else if (p<cumpropensities[1]) {
  // rxn: = A1_nonsense 
  if (Burst1_is_mutated_allele1 == 1){
    Anonsense1=Anonsense1 + 1;
  }
  
  updateSpecies(species, A1, Anonsense1, Aprime1, Aprime2, B1, Burst1_onorig_targ_allele1, Burst1_onpara1_targ_allele1, Burst1_onpara2_targ_allele1,  Burst1_off_targ_allele1, Burst1_onorig_targ_allele2, Burst1_onpara1_targ_allele2, Burst1_onpara2_targ_allele2, Burst1_off_targ_allele2, Burst1_on_para1_allele1, Burst1_on_para1_allele2, Burst1_on_para2_allele1, Burst1_on_para2_allele2, Burst1_on_orig_allele1, Burst1_on_orig_allele2, Burst1_is_mutated_allele1, Burst1_is_mutated_allele2);
  updatePropensities(species, rates, propensities);

} else if (p<cumpropensities[2]) {
    // rxn: = A1 from allele2
  if (Burst1_is_mutated_allele2 == 0){ // extra check for mutation ?necessary?
    A1=A1 + 1;
  }
  
  updateSpecies(species, A1, Anonsense1, Aprime1, Aprime2, B1, Burst1_onorig_targ_allele1, Burst1_onpara1_targ_allele1, Burst1_onpara2_targ_allele1,  Burst1_off_targ_allele1, Burst1_onorig_targ_allele2, Burst1_onpara1_targ_allele2, Burst1_onpara2_targ_allele2, Burst1_off_targ_allele2, Burst1_on_para1_allele1, Burst1_on_para1_allele2, Burst1_on_para2_allele1, Burst1_on_para2_allele2, Burst1_on_orig_allele1, Burst1_on_orig_allele2, Burst1_is_mutated_allele1, Burst1_is_mutated_allele2);
  updatePropensities(species, rates, propensities);

} else if (p<cumpropensities[3]) {
  // rxn: = A1_nonsense from allele2
  if (Burst1_is_mutated_allele2 == 1){
    Anonsense1=Anonsense1 + 1;
  }

  updateSpecies(species, A1, Anonsense1, Aprime1, Aprime2, B1, Burst1_onorig_targ_allele1, Burst1_onpara1_targ_allele1, Burst1_onpara2_targ_allele1,  Burst1_off_targ_allele1, Burst1_onorig_targ_allele2, Burst1_onpara1_targ_allele2, Burst1_onpara2_targ_allele2, Burst1_off_targ_allele2, Burst1_on_para1_allele1, Burst1_on_para1_allele2, Burst1_on_para2_allele1, Burst1_on_para2_allele2, Burst1_on_orig_allele1, Burst1_on_orig_allele2, Burst1_is_mutated_allele1, Burst1_is_mutated_allele2);
  updatePropensities(species, rates, propensities);

} else if (p<cumpropensities[4]) {
  // rxn: = Aprime1 from allele1
  Aprime1=Aprime1 + 1;

  updateSpecies(species, A1, Anonsense1, Aprime1, Aprime2, B1, Burst1_onorig_targ_allele1, Burst1_onpara1_targ_allele1, Burst1_onpara2_targ_allele1,  Burst1_off_targ_allele1, Burst1_onorig_targ_allele2, Burst1_onpara1_targ_allele2, Burst1_onpara2_targ_allele2, Burst1_off_targ_allele2, Burst1_on_para1_allele1, Burst1_on_para1_allele2, Burst1_on_para2_allele1, Burst1_on_para2_allele2, Burst1_on_orig_allele1, Burst1_on_orig_allele2, Burst1_is_mutated_allele1, Burst1_is_mutated_allele2);
  updatePropensities(species, rates, propensities);

} else if (p<cumpropensities[5]) {
  // rxn: = Aprime1 from allele2
  Aprime1=Aprime1 + 1;

  updateSpecies(species, A1, Anonsense1, Aprime1, Aprime2, B1, Burst1_onorig_targ_allele1, Burst1_onpara1_targ_allele1, Burst1_onpara2_targ_allele1,  Burst1_off_targ_allele1, Burst1_onorig_targ_allele2, Burst1_onpara1_targ_allele2, Burst1_onpara2_targ_allele2, Burst1_off_targ_allele2, Burst1_on_para1_allele1, Burst1_on_para1_allele2, Burst1_on_para2_allele1, Burst1_on_para2_allele2, Burst1_on_orig_allele1, Burst1_on_orig_allele2, Burst1_is_mutated_allele1, Burst1_is_mutated_allele2);
  updatePropensities(species, rates, propensities);

} else if (p<cumpropensities[6]) {
  // rxn: = Aprime2 from allele1
  Aprime2=Aprime2 + 1;

  updateSpecies(species, A1, Anonsense1, Aprime1, Aprime2, B1, Burst1_onorig_targ_allele1, Burst1_onpara1_targ_allele1, Burst1_onpara2_targ_allele1,  Burst1_off_targ_allele1, Burst1_onorig_targ_allele2, Burst1_onpara1_targ_allele2, Burst1_onpara2_targ_allele2, Burst1_off_targ_allele2, Burst1_on_para1_allele1, Burst1_on_para1_allele2, Burst1_on_para2_allele1, Burst1_on_para2_allele2, Burst1_on_orig_allele1, Burst1_on_orig_allele2, Burst1_is_mutated_allele1, Burst1_is_mutated_allele2);
  updatePropensities(species, rates, propensities);

} else if (p<cumpropensities[7]) {
  // rxn: = Aprime2 from allele2
  Aprime2=Aprime2 + 1;

  updateSpecies(species, A1, Anonsense1, Aprime1, Aprime2, B1, Burst1_onorig_targ_allele1, Burst1_onpara1_targ_allele1, Burst1_onpara2_targ_allele1,  Burst1_off_targ_allele1, Burst1_onorig_targ_allele2, Burst1_onpara1_targ_allele2, Burst1_onpara2_targ_allele2, Burst1_off_targ_allele2, Burst1_on_para1_allele1, Burst1_on_para1_allele2, Burst1_on_para2_allele1, Burst1_on_para2_allele2, Burst1_on_orig_allele1, Burst1_on_orig_allele2, Burst1_is_mutated_allele1, Burst1_is_mutated_allele2);
  updatePropensities(species, rates, propensities);

} else if (p<cumpropensities[8]) {
  // rxn: = B1 from allele1
  B1=B1 + 1;

  updateSpecies(species, A1, Anonsense1, Aprime1, Aprime2, B1, Burst1_onorig_targ_allele1, Burst1_onpara1_targ_allele1, Burst1_onpara2_targ_allele1,  Burst1_off_targ_allele1, Burst1_onorig_targ_allele2, Burst1_onpara1_targ_allele2, Burst1_onpara2_targ_allele2, Burst1_off_targ_allele2, Burst1_on_para1_allele1, Burst1_on_para1_allele2, Burst1_on_para2_allele1, Burst1_on_para2_allele2, Burst1_on_orig_allele1, Burst1_on_orig_allele2, Burst1_is_mutated_allele1, Burst1_is_mutated_allele2);
  updatePropensities(species, rates, propensities);

} else if (p<cumpropensities[9]) {
  // rxn: = B1 from allele2
  B1=B1 + 1;

  updateSpecies(species, A1, Anonsense1, Aprime1, Aprime2, B1, Burst1_onorig_targ_allele1, Burst1_onpara1_targ_allele1, Burst1_onpara2_targ_allele1,  Burst1_off_targ_allele1, Burst1_onorig_targ_allele2, Burst1_onpara1_targ_allele2, Burst1_onpara2_targ_allele2, Burst1_off_targ_allele2, Burst1_on_para1_allele1, Burst1_on_para1_allele2, Burst1_on_para2_allele1, Burst1_on_para2_allele2, Burst1_on_orig_allele1, Burst1_on_orig_allele2, Burst1_is_mutated_allele1, Burst1_is_mutated_allele2);
  updatePropensities(species, rates, propensities);

} else if (p<cumpropensities[10]) {
  // rxn: A1 = 
  A1=A1 + -1;

  updateSpecies(species, A1, Anonsense1, Aprime1, Aprime2, B1, Burst1_onorig_targ_allele1, Burst1_onpara1_targ_allele1, Burst1_onpara2_targ_allele1,  Burst1_off_targ_allele1, Burst1_onorig_targ_allele2, Burst1_onpara1_targ_allele2, Burst1_onpara2_targ_allele2, Burst1_off_targ_allele2, Burst1_on_para1_allele1, Burst1_on_para1_allele2, Burst1_on_para2_allele1, Burst1_on_para2_allele2, Burst1_on_orig_allele1, Burst1_on_orig_allele2, Burst1_is_mutated_allele1, Burst1_is_mutated_allele2);
  updatePropensities(species, rates, propensities);

} else if (p<cumpropensities[11]) {
  // rxn: Aprime1 = 
  Anonsense1=Anonsense1 + -1;

  updateSpecies(species, A1, Anonsense1, Aprime1, Aprime2, B1, Burst1_onorig_targ_allele1, Burst1_onpara1_targ_allele1, Burst1_onpara2_targ_allele1,  Burst1_off_targ_allele1, Burst1_onorig_targ_allele2, Burst1_onpara1_targ_allele2, Burst1_onpara2_targ_allele2, Burst1_off_targ_allele2, Burst1_on_para1_allele1, Burst1_on_para1_allele2, Burst1_on_para2_allele1, Burst1_on_para2_allele2, Burst1_on_orig_allele1, Burst1_on_orig_allele2, Burst1_is_mutated_allele1, Burst1_is_mutated_allele2);
  updatePropensities(species, rates, propensities);

} else if (p<cumpropensities[12]) {
  // rxn: Aprime1 = 
  Aprime1=Aprime1 + -1;

  updateSpecies(species, A1, Anonsense1, Aprime1, Aprime2, B1, Burst1_onorig_targ_allele1, Burst1_onpara1_targ_allele1, Burst1_onpara2_targ_allele1,  Burst1_off_targ_allele1, Burst1_onorig_targ_allele2, Burst1_onpara1_targ_allele2, Burst1_onpara2_targ_allele2, Burst1_off_targ_allele2, Burst1_on_para1_allele1, Burst1_on_para1_allele2, Burst1_on_para2_allele1, Burst1_on_para2_allele2, Burst1_on_orig_allele1, Burst1_on_orig_allele2, Burst1_is_mutated_allele1, Burst1_is_mutated_allele2);
  updatePropensities(species, rates, propensities);

} else if (p<cumpropensities[13]) {
  // rxn: Aprime1 = 
  Aprime2=Aprime2 + -1;

  updateSpecies(species, A1, Anonsense1, Aprime1, Aprime2, B1, Burst1_onorig_targ_allele1, Burst1_onpara1_targ_allele1, Burst1_onpara2_targ_allele1,  Burst1_off_targ_allele1, Burst1_onorig_targ_allele2, Burst1_onpara1_targ_allele2, Burst1_onpara2_targ_allele2, Burst1_off_targ_allele2, Burst1_on_para1_allele1, Burst1_on_para1_allele2, Burst1_on_para2_allele1, Burst1_on_para2_allele2, Burst1_on_orig_allele1, Burst1_on_orig_allele2, Burst1_is_mutated_allele1, Burst1_is_mutated_allele2);
  updatePropensities(species, rates, propensities);

} else if (p<cumpropensities[14]) {
  // rxn: B1 = 
  B1=B1 + -1;

  updateSpecies(species, A1, Anonsense1, Aprime1, Aprime2, B1, Burst1_onorig_targ_allele1, Burst1_onpara1_targ_allele1, Burst1_onpara2_targ_allele1,  Burst1_off_targ_allele1, Burst1_onorig_targ_allele2, Burst1_onpara1_targ_allele2, Burst1_onpara2_targ_allele2, Burst1_off_targ_allele2, Burst1_on_para1_allele1, Burst1_on_para1_allele2, Burst1_on_para2_allele1, Burst1_on_para2_allele2, Burst1_on_orig_allele1, Burst1_on_orig_allele2, Burst1_is_mutated_allele1, Burst1_is_mutated_allele2);
  updatePropensities(species, rates, propensities);

} else if (p<cumpropensities[15]) {
  // rxn: Burst1_off_orig_allele1 = Burst1_on_orig_allele1

  if(Burst1_on_orig_allele1 == 1) { // this should never happen, so a value of 2 is an error flag
    Burst1_on_orig_allele1 = 2; 
  }  
  if(Burst1_on_orig_allele1 == 0) {
    Burst1_on_orig_allele1 = 1;
  }

  updateSpecies(species, A1, Anonsense1, Aprime1, Aprime2, B1, Burst1_onorig_targ_allele1, Burst1_onpara1_targ_allele1, Burst1_onpara2_targ_allele1,  Burst1_off_targ_allele1, Burst1_onorig_targ_allele2, Burst1_onpara1_targ_allele2, Burst1_onpara2_targ_allele2, Burst1_off_targ_allele2, Burst1_on_para1_allele1, Burst1_on_para1_allele2, Burst1_on_para2_allele1, Burst1_on_para2_allele2, Burst1_on_orig_allele1, Burst1_on_orig_allele2, Burst1_is_mutated_allele1, Burst1_is_mutated_allele2);
  updatePropensities(species, rates, propensities);

} else if (p<cumpropensities[16]) {
  // rxn: Burst1_off_orig_allele2 = Burst1_on_orig_allele2
  
  if(Burst1_on_orig_allele2 == 1) { // this should never happen, so a value of 2 is an error flag
    Burst1_on_orig_allele2 = 2; 
  }
  if(Burst1_on_orig_allele2 == 0) {
    Burst1_on_orig_allele2 = 1;
  }

  updateSpecies(species, A1, Anonsense1, Aprime1, Aprime2, B1, Burst1_onorig_targ_allele1, Burst1_onpara1_targ_allele1, Burst1_onpara2_targ_allele1,  Burst1_off_targ_allele1, Burst1_onorig_targ_allele2, Burst1_onpara1_targ_allele2, Burst1_onpara2_targ_allele2, Burst1_off_targ_allele2, Burst1_on_para1_allele1, Burst1_on_para1_allele2, Burst1_on_para2_allele1, Burst1_on_para2_allele2, Burst1_on_orig_allele1, Burst1_on_orig_allele2, Burst1_is_mutated_allele1, Burst1_is_mutated_allele2);
  updatePropensities(species, rates, propensities);

} else if (p<cumpropensities[17]) {
  // rxn: Burst1_off_para1_allele1 = Burst1_on_para1_allele1
  
  if(Burst1_on_para1_allele1 == 1) { // this should never happen, so a value of 2 is an error flag
    Burst1_on_para1_allele1 = 2; 
  }  
  if(Burst1_on_para1_allele1 == 0) {
    Burst1_on_para1_allele1 = 1;
  }

  updateSpecies(species, A1, Anonsense1, Aprime1, Aprime2, B1, Burst1_onorig_targ_allele1, Burst1_onpara1_targ_allele1, Burst1_onpara2_targ_allele1,  Burst1_off_targ_allele1, Burst1_onorig_targ_allele2, Burst1_onpara1_targ_allele2, Burst1_onpara2_targ_allele2, Burst1_off_targ_allele2, Burst1_on_para1_allele1, Burst1_on_para1_allele2, Burst1_on_para2_allele1, Burst1_on_para2_allele2, Burst1_on_orig_allele1, Burst1_on_orig_allele2, Burst1_is_mutated_allele1, Burst1_is_mutated_allele2);
  updatePropensities(species, rates, propensities);

} else if (p<cumpropensities[18]) {
  // rxn: Burst1_off_para1_allele2 = Burst1_on_para1_allele2

  if(Burst1_on_para1_allele2 == 1) { // this should never happen, so a value of 2 is an error flag
    Burst1_on_para1_allele2 = 2; 
  }  
  if(Burst1_on_para1_allele2 == 0) {
    Burst1_on_para1_allele2 = 1;
  }

  updateSpecies(species, A1, Anonsense1, Aprime1, Aprime2, B1, Burst1_onorig_targ_allele1, Burst1_onpara1_targ_allele1, Burst1_onpara2_targ_allele1,  Burst1_off_targ_allele1, Burst1_onorig_targ_allele2, Burst1_onpara1_targ_allele2, Burst1_onpara2_targ_allele2, Burst1_off_targ_allele2, Burst1_on_para1_allele1, Burst1_on_para1_allele2, Burst1_on_para2_allele1, Burst1_on_para2_allele2, Burst1_on_orig_allele1, Burst1_on_orig_allele2, Burst1_is_mutated_allele1, Burst1_is_mutated_allele2);
  updatePropensities(species, rates, propensities);

} else if (p<cumpropensities[19]) {
  // rxn: Burst1_off_para2_allele1 = Burst1_on_para2_allele1
  
  if(Burst1_on_para2_allele1 == 1) { // this should never happen, so a value of 2 is an error flag
    Burst1_on_para2_allele1 = 2; 
  }  
  if(Burst1_on_para2_allele1 == 0) {
    Burst1_on_para2_allele1 = 1;
  }

  updateSpecies(species, A1, Anonsense1, Aprime1, Aprime2, B1, Burst1_onorig_targ_allele1, Burst1_onpara1_targ_allele1, Burst1_onpara2_targ_allele1,  Burst1_off_targ_allele1, Burst1_onorig_targ_allele2, Burst1_onpara1_targ_allele2, Burst1_onpara2_targ_allele2, Burst1_off_targ_allele2, Burst1_on_para1_allele1, Burst1_on_para1_allele2, Burst1_on_para2_allele1, Burst1_on_para2_allele2, Burst1_on_orig_allele1, Burst1_on_orig_allele2, Burst1_is_mutated_allele1, Burst1_is_mutated_allele2);
  updatePropensities(species, rates, propensities);

} else if (p<cumpropensities[20]) {
  // rxn: Burst1_off_para2_allele2 = Burst1_on_para2_allele2

  if(Burst1_on_para2_allele2 == 1) { // this should never happen, so a value of 2 is an error flag
    Burst1_on_para2_allele2 = 2; 
  }  
  if(Burst1_on_para2_allele2 == 0) {
    Burst1_on_para2_allele2 = 1;
  }

  updateSpecies(species, A1, Anonsense1, Aprime1, Aprime2, B1, Burst1_onorig_targ_allele1, Burst1_onpara1_targ_allele1, Burst1_onpara2_targ_allele1,  Burst1_off_targ_allele1, Burst1_onorig_targ_allele2, Burst1_onpara1_targ_allele2, Burst1_onpara2_targ_allele2, Burst1_off_targ_allele2, Burst1_on_para1_allele1, Burst1_on_para1_allele2, Burst1_on_para2_allele1, Burst1_on_para2_allele2, Burst1_on_orig_allele1, Burst1_on_orig_allele2, Burst1_is_mutated_allele1, Burst1_is_mutated_allele2);
  updatePropensities(species, rates, propensities);

} else if (p<cumpropensities[21]) {
  // rxn: Burst1_off_targ_allele1 = Burst1_onorig_targ_allele1
 
  Burst1_off_targ_allele1 = 0;
  Burst1_onorig_targ_allele1 = 1;
  Burst1_onpara1_targ_allele1 = 0;
  Burst1_onpara2_targ_allele1 = 0;

  updateSpecies(species, A1, Anonsense1, Aprime1, Aprime2, B1, Burst1_onorig_targ_allele1, Burst1_onpara1_targ_allele1, Burst1_onpara2_targ_allele1,  Burst1_off_targ_allele1, Burst1_onorig_targ_allele2, Burst1_onpara1_targ_allele2, Burst1_onpara2_targ_allele2, Burst1_off_targ_allele2, Burst1_on_para1_allele1, Burst1_on_para1_allele2, Burst1_on_para2_allele1, Burst1_on_para2_allele2, Burst1_on_orig_allele1, Burst1_on_orig_allele2, Burst1_is_mutated_allele1, Burst1_is_mutated_allele2);
  updatePropensities(species, rates, propensities);

} else if (p<cumpropensities[22]) {
  // rxn: Burst1_off_targ_allele1 = Burst1_onpara1_targ_allele1
  
  Burst1_off_targ_allele1 = 0;
  Burst1_onorig_targ_allele1 = 0;
  Burst1_onpara1_targ_allele1 = 1;
  Burst1_onpara2_targ_allele1 = 0;

  updateSpecies(species, A1, Anonsense1, Aprime1, Aprime2, B1, Burst1_onorig_targ_allele1, Burst1_onpara1_targ_allele1, Burst1_onpara2_targ_allele1,  Burst1_off_targ_allele1, Burst1_onorig_targ_allele2, Burst1_onpara1_targ_allele2, Burst1_onpara2_targ_allele2, Burst1_off_targ_allele2, Burst1_on_para1_allele1, Burst1_on_para1_allele2, Burst1_on_para2_allele1, Burst1_on_para2_allele2, Burst1_on_orig_allele1, Burst1_on_orig_allele2, Burst1_is_mutated_allele1, Burst1_is_mutated_allele2);
  updatePropensities(species, rates, propensities);

} else if (p<cumpropensities[23]) {
  // rxn: Burst1_off_targ_allele1 = Burst1_onpara1_targ_allele1
  
  Burst1_off_targ_allele1 = 0;
  Burst1_onorig_targ_allele1 = 0;
  Burst1_onpara1_targ_allele1 = 0;
  Burst1_onpara2_targ_allele1 = 1;

  updateSpecies(species, A1, Anonsense1, Aprime1, Aprime2, B1, Burst1_onorig_targ_allele1, Burst1_onpara1_targ_allele1, Burst1_onpara2_targ_allele1,  Burst1_off_targ_allele1, Burst1_onorig_targ_allele2, Burst1_onpara1_targ_allele2, Burst1_onpara2_targ_allele2, Burst1_off_targ_allele2, Burst1_on_para1_allele1, Burst1_on_para1_allele2, Burst1_on_para2_allele1, Burst1_on_para2_allele2, Burst1_on_orig_allele1, Burst1_on_orig_allele2, Burst1_is_mutated_allele1, Burst1_is_mutated_allele2);
  updatePropensities(species, rates, propensities);

} else if (p<cumpropensities[24]) {
  // rxn: Burst1_off_targ_allele2 = Burst1_onorig_targ_allele2
 
  Burst1_off_targ_allele2 = 0;
  Burst1_onorig_targ_allele2 = 1;
  Burst1_onpara1_targ_allele2 = 0;
  Burst1_onpara2_targ_allele2 = 0;

  updateSpecies(species, A1, Anonsense1, Aprime1, Aprime2, B1, Burst1_onorig_targ_allele1, Burst1_onpara1_targ_allele1, Burst1_onpara2_targ_allele1,  Burst1_off_targ_allele1, Burst1_onorig_targ_allele2, Burst1_onpara1_targ_allele2, Burst1_onpara2_targ_allele2, Burst1_off_targ_allele2, Burst1_on_para1_allele1, Burst1_on_para1_allele2, Burst1_on_para2_allele1, Burst1_on_para2_allele2, Burst1_on_orig_allele1, Burst1_on_orig_allele2, Burst1_is_mutated_allele1, Burst1_is_mutated_allele2);
  updatePropensities(species, rates, propensities);

} else if (p<cumpropensities[25]) {
  // rxn: Burst1_off_targ_allele2 = Burst1_onpara1_targ_allele2
  
  Burst1_off_targ_allele2 = 0;
  Burst1_onorig_targ_allele2 = 0;
  Burst1_onpara1_targ_allele2 = 1;
  Burst1_onpara2_targ_allele2 = 0;

  updateSpecies(species, A1, Anonsense1, Aprime1, Aprime2, B1, Burst1_onorig_targ_allele1, Burst1_onpara1_targ_allele1, Burst1_onpara2_targ_allele1,  Burst1_off_targ_allele1, Burst1_onorig_targ_allele2, Burst1_onpara1_targ_allele2, Burst1_onpara2_targ_allele2, Burst1_off_targ_allele2, Burst1_on_para1_allele1, Burst1_on_para1_allele2, Burst1_on_para2_allele1, Burst1_on_para2_allele2, Burst1_on_orig_allele1, Burst1_on_orig_allele2, Burst1_is_mutated_allele1, Burst1_is_mutated_allele2);
  updatePropensities(species, rates, propensities);

} else if (p<cumpropensities[26]) {
  // rxn: Burst1_off_targ_allele2 = Burst1_onpara1_targ_allele2
  
  Burst1_off_targ_allele2 = 0;
  Burst1_onorig_targ_allele2 = 0;
  Burst1_onpara1_targ_allele2 = 0;
  Burst1_onpara2_targ_allele2 = 1;

  updateSpecies(species, A1, Anonsense1, Aprime1, Aprime2, B1, Burst1_onorig_targ_allele1, Burst1_onpara1_targ_allele1, Burst1_onpara2_targ_allele1,  Burst1_off_targ_allele1, Burst1_onorig_targ_allele2, Burst1_onpara1_targ_allele2, Burst1_onpara2_targ_allele2, Burst1_off_targ_allele2, Burst1_on_para1_allele1, Burst1_on_para1_allele2, Burst1_on_para2_allele1, Burst1_on_para2_allele2, Burst1_on_orig_allele1, Burst1_on_orig_allele2, Burst1_is_mutated_allele1, Burst1_is_mutated_allele2);
  updatePropensities(species, rates, propensities);

} else if (p<cumpropensities[27]) {
  // rxn: Burst1_on_orig_allele1  = Burst1_off_orig_allele1

  if(Burst1_on_orig_allele1 == 0) {
    Burst1_on_orig_allele1 = 2; // should never happen, so a 2 value is an error flag
  }  
  if(Burst1_on_orig_allele1 == 1) {
    Burst1_on_orig_allele1 = 0;
  }

  updateSpecies(species, A1, Anonsense1, Aprime1, Aprime2, B1, Burst1_onorig_targ_allele1, Burst1_onpara1_targ_allele1, Burst1_onpara2_targ_allele1,  Burst1_off_targ_allele1, Burst1_onorig_targ_allele2, Burst1_onpara1_targ_allele2, Burst1_onpara2_targ_allele2, Burst1_off_targ_allele2, Burst1_on_para1_allele1, Burst1_on_para1_allele2, Burst1_on_para2_allele1, Burst1_on_para2_allele2, Burst1_on_orig_allele1, Burst1_on_orig_allele2, Burst1_is_mutated_allele1, Burst1_is_mutated_allele2);
  updatePropensities(species, rates, propensities);

} else if (p<cumpropensities[28]) {
  // rxn: Burst1_on_orig_allele2  = Burst1_off_orig_allele2
  
  if(Burst1_on_orig_allele2 == 0) {
    Burst1_on_orig_allele2 = 2; // should never happen, so a 2 value is an error flag
  }
  if(Burst1_on_orig_allele2 == 1) {
    Burst1_on_orig_allele2 = 0;
  }

  updateSpecies(species, A1, Anonsense1, Aprime1, Aprime2, B1, Burst1_onorig_targ_allele1, Burst1_onpara1_targ_allele1, Burst1_onpara2_targ_allele1,  Burst1_off_targ_allele1, Burst1_onorig_targ_allele2, Burst1_onpara1_targ_allele2, Burst1_onpara2_targ_allele2, Burst1_off_targ_allele2, Burst1_on_para1_allele1, Burst1_on_para1_allele2, Burst1_on_para2_allele1, Burst1_on_para2_allele2, Burst1_on_orig_allele1, Burst1_on_orig_allele2, Burst1_is_mutated_allele1, Burst1_is_mutated_allele2);
  updatePropensities(species, rates, propensities);

} else if (p<cumpropensities[29]) {
  // rxn: Burst1_on_para1_allele1 = Burst1_off_para1_allele1

  if(Burst1_on_para1_allele1 == 0) {
    Burst1_on_para1_allele1 = 2; // should never happen, so a 2 value is an error flag
  }
  if(Burst1_on_para1_allele1 == 1) {
    Burst1_on_para1_allele1 = 0;
  }

  updateSpecies(species, A1, Anonsense1, Aprime1, Aprime2, B1, Burst1_onorig_targ_allele1, Burst1_onpara1_targ_allele1, Burst1_onpara2_targ_allele1,  Burst1_off_targ_allele1, Burst1_onorig_targ_allele2, Burst1_onpara1_targ_allele2, Burst1_onpara2_targ_allele2, Burst1_off_targ_allele2, Burst1_on_para1_allele1, Burst1_on_para1_allele2, Burst1_on_para2_allele1, Burst1_on_para2_allele2, Burst1_on_orig_allele1, Burst1_on_orig_allele2, Burst1_is_mutated_allele1, Burst1_is_mutated_allele2);
  updatePropensities(species, rates, propensities);

} else if (p<cumpropensities[30]) {
  // rxn: Burst1_on_para1_allele2 = Burst1_off_para1_allele2

  if(Burst1_on_para1_allele2 == 0) {
    Burst1_on_para1_allele2 = 2; // should never happen, so a 2 value is an error flag
  }
  if(Burst1_on_para1_allele2 == 1) {
    Burst1_on_para1_allele2 = 0;
  }
  
  updateSpecies(species, A1, Anonsense1, Aprime1, Aprime2, B1, Burst1_onorig_targ_allele1, Burst1_onpara1_targ_allele1, Burst1_onpara2_targ_allele1,  Burst1_off_targ_allele1, Burst1_onorig_targ_allele2, Burst1_onpara1_targ_allele2, Burst1_onpara2_targ_allele2, Burst1_off_targ_allele2, Burst1_on_para1_allele1, Burst1_on_para1_allele2, Burst1_on_para2_allele1, Burst1_on_para2_allele2, Burst1_on_orig_allele1, Burst1_on_orig_allele2, Burst1_is_mutated_allele1, Burst1_is_mutated_allele2);
  updatePropensities(species, rates, propensities);

} else if (p<cumpropensities[31]) {
  // rxn: Burst1_on_para2_allele1 = Burst1_off_para2_allele1

  if(Burst1_on_para2_allele1 == 0) {
    Burst1_on_para2_allele1 = 2; // should never happen, so a 2 value is an error flag
  }
  if(Burst1_on_para2_allele1 == 1) {
    Burst1_on_para2_allele1 = 0;
  }

  updateSpecies(species, A1, Anonsense1, Aprime1, Aprime2, B1, Burst1_onorig_targ_allele1, Burst1_onpara1_targ_allele1, Burst1_onpara2_targ_allele1,  Burst1_off_targ_allele1, Burst1_onorig_targ_allele2, Burst1_onpara1_targ_allele2, Burst1_onpara2_targ_allele2, Burst1_off_targ_allele2, Burst1_on_para1_allele1, Burst1_on_para1_allele2, Burst1_on_para2_allele1, Burst1_on_para2_allele2, Burst1_on_orig_allele1, Burst1_on_orig_allele2, Burst1_is_mutated_allele1, Burst1_is_mutated_allele2);
  updatePropensities(species, rates, propensities);

} else if (p<cumpropensities[32]) {
  // rxn: Burst1_on_para2_allele2 = Burst1_off_para2_allele2

  if(Burst1_on_para2_allele2 == 0) {
    Burst1_on_para2_allele2 = 2; // should never happen, so a 2 value is an error flag
  }
  if(Burst1_on_para2_allele2 == 1) {
    Burst1_on_para2_allele2 = 0;
  }
  
  updateSpecies(species, A1, Anonsense1, Aprime1, Aprime2, B1, Burst1_onorig_targ_allele1, Burst1_onpara1_targ_allele1, Burst1_onpara2_targ_allele1,  Burst1_off_targ_allele1, Burst1_onorig_targ_allele2, Burst1_onpara1_targ_allele2, Burst1_onpara2_targ_allele2, Burst1_off_targ_allele2, Burst1_on_para1_allele1, Burst1_on_para1_allele2, Burst1_on_para2_allele1, Burst1_on_para2_allele2, Burst1_on_orig_allele1, Burst1_on_orig_allele2, Burst1_is_mutated_allele1, Burst1_is_mutated_allele2);
  updatePropensities(species, rates, propensities);

} else if (p<cumpropensities[33]) {
  // rxn: Burst1_onorig_targ_allele1   = Burst1_off_targ_allele1

  Burst1_off_targ_allele1 = 1;
  Burst1_onorig_targ_allele1 = 0;
  Burst1_onpara1_targ_allele1 = 0;
  Burst1_onpara2_targ_allele1 = 0;

  updateSpecies(species, A1, Anonsense1, Aprime1, Aprime2, B1, Burst1_onorig_targ_allele1, Burst1_onpara1_targ_allele1, Burst1_onpara2_targ_allele1,  Burst1_off_targ_allele1, Burst1_onorig_targ_allele2, Burst1_onpara1_targ_allele2, Burst1_onpara2_targ_allele2, Burst1_off_targ_allele2, Burst1_on_para1_allele1, Burst1_on_para1_allele2, Burst1_on_para2_allele1, Burst1_on_para2_allele2, Burst1_on_orig_allele1, Burst1_on_orig_allele2, Burst1_is_mutated_allele1, Burst1_is_mutated_allele2);
  updatePropensities(species, rates, propensities);

} else if (p<cumpropensities[34]) {
  // rxn: Burst1_onpara1_targ_allele1   = Burst1_off_targ_allele1

  Burst1_off_targ_allele1 = 1;
  Burst1_onorig_targ_allele1 = 0;
  Burst1_onpara1_targ_allele1 = 0;
  Burst1_onpara2_targ_allele1 = 0;

  updateSpecies(species, A1, Anonsense1, Aprime1, Aprime2, B1, Burst1_onorig_targ_allele1, Burst1_onpara1_targ_allele1, Burst1_onpara2_targ_allele1,  Burst1_off_targ_allele1, Burst1_onorig_targ_allele2, Burst1_onpara1_targ_allele2, Burst1_onpara2_targ_allele2, Burst1_off_targ_allele2, Burst1_on_para1_allele1, Burst1_on_para1_allele2, Burst1_on_para2_allele1, Burst1_on_para2_allele2, Burst1_on_orig_allele1, Burst1_on_orig_allele2, Burst1_is_mutated_allele1, Burst1_is_mutated_allele2);
  updatePropensities(species, rates, propensities);

} else if (p<cumpropensities[35]) {
  // rxn: Burst1_onpara2_targ_allele1   = Burst1_off_targ_allele1

  Burst1_off_targ_allele1 = 1;
  Burst1_onorig_targ_allele1 = 0;
  Burst1_onpara1_targ_allele1 = 0;
  Burst1_onpara2_targ_allele1 = 0;

  updateSpecies(species, A1, Anonsense1, Aprime1, Aprime2, B1, Burst1_onorig_targ_allele1, Burst1_onpara1_targ_allele1, Burst1_onpara2_targ_allele1,  Burst1_off_targ_allele1, Burst1_onorig_targ_allele2, Burst1_onpara1_targ_allele2, Burst1_onpara2_targ_allele2, Burst1_off_targ_allele2, Burst1_on_para1_allele1, Burst1_on_para1_allele2, Burst1_on_para2_allele1, Burst1_on_para2_allele2, Burst1_on_orig_allele1, Burst1_on_orig_allele2, Burst1_is_mutated_allele1, Burst1_is_mutated_allele2);
  updatePropensities(species, rates, propensities);

} else if (p<cumpropensities[36]) {
  // rxn: Burst1_onorig_targ_allele2   = Burst1_off_targ_allele2

  Burst1_off_targ_allele2 = 1;
  Burst1_onorig_targ_allele2 = 0;
  Burst1_onpara1_targ_allele2 = 0;
  Burst1_onpara2_targ_allele2 = 0;

  updateSpecies(species, A1, Anonsense1, Aprime1, Aprime2, B1, Burst1_onorig_targ_allele1, Burst1_onpara1_targ_allele1, Burst1_onpara2_targ_allele1,  Burst1_off_targ_allele1, Burst1_onorig_targ_allele2, Burst1_onpara1_targ_allele2, Burst1_onpara2_targ_allele2, Burst1_off_targ_allele2, Burst1_on_para1_allele1, Burst1_on_para1_allele2, Burst1_on_para2_allele1, Burst1_on_para2_allele2, Burst1_on_orig_allele1, Burst1_on_orig_allele2, Burst1_is_mutated_allele1, Burst1_is_mutated_allele2);
  updatePropensities(species, rates, propensities);

  } else if (p<cumpropensities[37]) {
  // rxn: Burst1_onpara1_targ_allele2   = Burst1_off_targ_allele2

  Burst1_off_targ_allele2 = 1;
  Burst1_onorig_targ_allele2 = 0;
  Burst1_onpara1_targ_allele2 = 0;
  Burst1_onpara1_targ_allele2 = 0;

  updateSpecies(species, A1, Anonsense1, Aprime1, Aprime2, B1, Burst1_onorig_targ_allele1, Burst1_onpara1_targ_allele1, Burst1_onpara2_targ_allele1,  Burst1_off_targ_allele1, Burst1_onorig_targ_allele2, Burst1_onpara1_targ_allele2, Burst1_onpara2_targ_allele2, Burst1_off_targ_allele2, Burst1_on_para1_allele1, Burst1_on_para1_allele2, Burst1_on_para2_allele1, Burst1_on_para2_allele2, Burst1_on_orig_allele1, Burst1_on_orig_allele2, Burst1_is_mutated_allele1, Burst1_is_mutated_allele2);
  updatePropensities(species, rates, propensities);
} else if (p<cumpropensities[38]) {
  // rxn: Burst1_onpara2_targ_allele2   = Burst1_off_targ_allele2

  Burst1_off_targ_allele2 = 1;
  Burst1_onorig_targ_allele2 = 0;
  Burst1_onpara1_targ_allele2 = 0;
  Burst1_onpara1_targ_allele2 = 0;

  updateSpecies(species, A1, Anonsense1, Aprime1, Aprime2, B1, Burst1_onorig_targ_allele1, Burst1_onpara1_targ_allele1, Burst1_onpara2_targ_allele1,  Burst1_off_targ_allele1, Burst1_onorig_targ_allele2, Burst1_onpara1_targ_allele2, Burst1_onpara2_targ_allele2, Burst1_off_targ_allele2, Burst1_on_para1_allele1, Burst1_on_para1_allele2, Burst1_on_para2_allele1, Burst1_on_para2_allele2, Burst1_on_orig_allele1, Burst1_on_orig_allele2, Burst1_is_mutated_allele1, Burst1_is_mutated_allele2);
  updatePropensities(species, rates, propensities);
}
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
