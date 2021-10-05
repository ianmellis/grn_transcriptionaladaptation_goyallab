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
#define NSPECIES 9

//INCLUDE NUMRXNS HERE
#define NUMRXNS 12


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
  long A1,Aprime1,B1,Burst1_off_targ,Burst1_on_targ,Burst1_off_para,Burst1_on_para,Burst1_on_orig,Burst1_off_orig;
  double A_prod1 ,Aprime_prod1 ,B_prod1 ,A_deg1 ,Aprime_deg1 ,B_deg1 ,B_ondep1 ,B_ondep_prime ,Aprimenitc1 ,A_off1 ,Aprime_off1 ,B_off1 ,A_proddiff1 ,Aprime_proddiff1 ,B_proddiff1 ,onbasal_a1 ,onbasal_aprime1 ,onbasal_b1 ,kA1 ,kAprime1 ,kB1 ,nA1 ,nAprime1 ,nB1 ;

  long i, j, k;
//UNPACK ALL SPECIES HERE
  A1 = (long)species[0];
  Aprime1 = (long)species[1];
  B1 = (long)species[2];
  Burst1_off_targ = (long)species[3];
  Burst1_on_targ = (long)species[4];
  Burst1_off_para = (long)species[5];
  Burst1_on_para = (long)species[6];
  Burst1_on_orig = (long)species[7];
  Burst1_off_orig = (long)species[8];
//UNPACK ALL RATES HERE
  A_prod1  = rates[0];
  Aprime_prod1  = rates[1];
  B_prod1  = rates[2];
  A_deg1  = rates[3];
  Aprime_deg1  = rates[4];
  B_deg1  = rates[5];
  B_ondep1  = rates[6];
  B_ondep_prime  = rates[7];
  Aprimenitc1  = rates[8];
  A_off1  = rates[9];
  Aprime_off1  = rates[10];
  B_off1  = rates[11];
  A_proddiff1  = rates[12];
  Aprime_proddiff1  = rates[13];
  B_proddiff1  = rates[14];
  onbasal_a1  = rates[15];
  onbasal_aprime1  = rates[16];
  onbasal_b1  = rates[17];
  kA1  = rates[18];
  kAprime1  = rates[19];
  kB1  = rates[20];
  nA1  = rates[21];
  nAprime1  = rates[22];
  nB1  = rates[23];

  zigset(seed);
  totaliterations = 0;

  savecount = 0;
  savecountcheck = 0;
  deltaTsave = totalt/(m-1);


  //  for (i=0; i<m; i++) {
  while (currT<totalt) {
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
species_out[savecount*NSPECIES+1] = Aprime1;
species_out[savecount*NSPECIES+2] = B1;
species_out[savecount*NSPECIES+3] = Burst1_off_targ;
species_out[savecount*NSPECIES+4] = Burst1_on_targ;
species_out[savecount*NSPECIES+5] = Burst1_off_para;
species_out[savecount*NSPECIES+6] = Burst1_on_para;
species_out[savecount*NSPECIES+7] = Burst1_on_orig;
species_out[savecount*NSPECIES+8] = Burst1_off_orig;
            times_out[savecount] = savecount*deltaTsave;
            savecount++;
            savecountcheck++;
            }
		else{savecount++;
		}
    }

//INSERT IF STATEMENT HERE
if (p<cumpropensities[0]) {
  // rxn: = A1 
  A1=A1 + 1;

  //update propensity for = A1 
  propensities[0] = A_prod1 *A_proddiff1 *Burst1_on_targ+A_prod1 *Burst1_off_targ;
  //update propensity for = Aprime1 
  propensities[1] = Aprime_prod1 *Aprime_proddiff1 *Burst1_on_para+Aprime_prod1 *Burst1_off_para;
  //update propensity for = B1 
  propensities[2] = B_prod1 *B_proddiff1 *Burst1_off_orig+B_prod1 *Burst1_on_orig;
  //update propensity for A1 = 
  propensities[3] = A_deg1 *A1;
  //update propensity for Aprime1 = 
  propensities[4] = Aprime_deg1 *Aprime1;
  //update propensity for B1 = 
  propensities[5] = B_deg1 *B1;
  //update propensity for Burst1_off_orig = Burst1_on_orig
  propensities[6] = onbasal_aprime1 *Burst1_off_orig;
  //update propensity for Burst1_off_para = Burst1_on_para
  propensities[7] = Aprimenitc1 *(pow(A1,nA1 )/(pow(kA1 ,nA1 )+pow(A1,nA1 )))*Burst1_off_para+onbasal_aprime1 *Burst1_off_para;
  //update propensity for Burst1_off_targ = Burst1_on_targ
  propensities[8] = B_ondep1 *(pow(A1,nA1 )/(pow(kA1 ,nA1 )+pow(A1,nA1 )))*Burst1_off_targ+ B_ondep_prime *(pow(Aprime1,nAprime1 )/(pow(kAprime1 ,nAprime1 )+pow(Aprime1,nAprime1 )))*Burst1_off_targ +onbasal_a1 *Burst1_off_targ;
  //update propensity for Burst1_on_orig = Burst1_off_orig
  propensities[9] = A_off1 *Burst1_on_targ;
  //update propensity for Burst1_on_para = Burst1_off_para
  propensities[10] = Aprime_off1 *Burst1_on_para;
  //update propensity for Burst1_on_targ = Burst1_off_targ
  propensities[11] = B_off1 *Burst1_off_orig;
} else if (p<cumpropensities[1]) {
  // rxn: = Aprime1 
  Aprime1=Aprime1 + 1;

  //update propensity for = A1 
  propensities[0] = A_prod1 *A_proddiff1 *Burst1_on_targ+A_prod1 *Burst1_off_targ;
  //update propensity for = Aprime1 
  propensities[1] = Aprime_prod1 *Aprime_proddiff1 *Burst1_on_para+Aprime_prod1 *Burst1_off_para;
  //update propensity for = B1 
  propensities[2] = B_prod1 *B_proddiff1 *Burst1_off_orig+B_prod1 *Burst1_on_orig;
  //update propensity for A1 = 
  propensities[3] = A_deg1 *A1;
  //update propensity for Aprime1 = 
  propensities[4] = Aprime_deg1 *Aprime1;
  //update propensity for B1 = 
  propensities[5] = B_deg1 *B1;
  //update propensity for Burst1_off_targ = Burst1_on_targ
  propensities[6] = B_ondep1 *(pow(A1,nA1 )/(pow(kA1 ,nA1 )+pow(A1,nA1 )))*Burst1_off_targ+onbasal_a1 *Burst1_off_targ;
  //update propensity for Burst1_off_targ = Burst1_on_targ
  propensities[7] = B_ondep_prime *(pow(A1,nA1 )/(pow(kA1 ,nA1 )+pow(A1,nA1 )))*Burst1_off_para+onbasal_aprime1 *Burst1_off_para;
  //update propensity for Burst1_off_para = Burst1_on_para
  propensities[8] = Aprimenitc1 *(pow(A1,nA1 )/(pow(kA1 ,nA1 )+pow(A1,nA1 )))*Burst1_on_orig+onbasal_b1 *Burst1_on_orig;
  //update propensity for Burst1_on_orig = Burst1_off_orig
  propensities[9] = A_off1 *Burst1_on_targ;
  //update propensity for Burst1_on_para = Burst1_off_para
  propensities[10] = Aprime_off1 *Burst1_on_para;
  //update propensity for Burst1_on_targ = Burst1_off_targ
  propensities[11] = B_off1 *Burst1_off_orig;
} else if (p<cumpropensities[2]) {
  // rxn: = B1 
  B1=B1 + 1;

  //update propensity for = A1 
  propensities[0] = A_prod1 *A_proddiff1 *Burst1_on_targ+A_prod1 *Burst1_off_targ;
  //update propensity for = Aprime1 
  propensities[1] = Aprime_prod1 *Aprime_proddiff1 *Burst1_on_para+Aprime_prod1 *Burst1_off_para;
  //update propensity for = B1 
  propensities[2] = B_prod1 *B_proddiff1 *Burst1_off_orig+B_prod1 *Burst1_on_orig;
  //update propensity for A1 = 
  propensities[3] = A_deg1 *A1;
  //update propensity for Aprime1 = 
  propensities[4] = Aprime_deg1 *Aprime1;
  //update propensity for B1 = 
  propensities[5] = B_deg1 *B1;
  //update propensity for Burst1_off_targ = Burst1_on_targ
  propensities[6] = B_ondep1 *(pow(A1,nA1 )/(pow(kA1 ,nA1 )+pow(A1,nA1 )))*Burst1_off_targ+onbasal_a1 *Burst1_off_targ;
  //update propensity for Burst1_off_targ = Burst1_on_targ
  propensities[7] = B_ondep_prime *(pow(A1,nA1 )/(pow(kA1 ,nA1 )+pow(A1,nA1 )))*Burst1_off_para+onbasal_aprime1 *Burst1_off_para;
  //update propensity for Burst1_off_para = Burst1_on_para
  propensities[8] = Aprimenitc1 *(pow(A1,nA1 )/(pow(kA1 ,nA1 )+pow(A1,nA1 )))*Burst1_on_orig+onbasal_b1 *Burst1_on_orig;
  //update propensity for Burst1_on_orig = Burst1_off_orig
  propensities[9] = A_off1 *Burst1_on_targ;
  //update propensity for Burst1_on_para = Burst1_off_para
  propensities[10] = Aprime_off1 *Burst1_on_para;
  //update propensity for Burst1_on_targ = Burst1_off_targ
  propensities[11] = B_off1 *Burst1_off_orig;
} else if (p<cumpropensities[3]) {
  // rxn: A1 = 
  A1=A1 + -1;

  //update propensity for = A1 
  propensities[0] = A_prod1 *A_proddiff1 *Burst1_on_targ+A_prod1 *Burst1_off_targ;
  //update propensity for = Aprime1 
  propensities[1] = Aprime_prod1 *Aprime_proddiff1 *Burst1_on_para+Aprime_prod1 *Burst1_off_para;
  //update propensity for = B1 
  propensities[2] = B_prod1 *B_proddiff1 *Burst1_off_orig+B_prod1 *Burst1_on_orig;
  //update propensity for A1 = 
  propensities[3] = A_deg1 *A1;
  //update propensity for Aprime1 = 
  propensities[4] = Aprime_deg1 *Aprime1;
  //update propensity for B1 = 
  propensities[5] = B_deg1 *B1;
  //update propensity for Burst1_off_targ = Burst1_on_targ
  propensities[6] = B_ondep1 *(pow(A1,nA1 )/(pow(kA1 ,nA1 )+pow(A1,nA1 )))*Burst1_off_targ+onbasal_a1 *Burst1_off_targ;
  //update propensity for Burst1_off_targ = Burst1_on_targ
  propensities[7] = B_ondep_prime *(pow(A1,nA1 )/(pow(kA1 ,nA1 )+pow(A1,nA1 )))*Burst1_off_para+onbasal_aprime1 *Burst1_off_para;
  //update propensity for Burst1_off_para = Burst1_on_para
  propensities[8] = Aprimenitc1 *(pow(A1,nA1 )/(pow(kA1 ,nA1 )+pow(A1,nA1 )))*Burst1_on_orig+onbasal_b1 *Burst1_on_orig;
  //update propensity for Burst1_on_orig = Burst1_off_orig
  propensities[9] = A_off1 *Burst1_on_targ;
  //update propensity for Burst1_on_para = Burst1_off_para
  propensities[10] = Aprime_off1 *Burst1_on_para;
  //update propensity for Burst1_on_targ = Burst1_off_targ
  propensities[11] = B_off1 *Burst1_off_orig;
} else if (p<cumpropensities[4]) {
  // rxn: Aprime1 = 
  Aprime1=Aprime1 + -1;

  //update propensity for = A1 
  propensities[0] = A_prod1 *A_proddiff1 *Burst1_on_targ+A_prod1 *Burst1_off_targ;
  //update propensity for = Aprime1 
  propensities[1] = Aprime_prod1 *Aprime_proddiff1 *Burst1_on_para+Aprime_prod1 *Burst1_off_para;
  //update propensity for = B1 
  propensities[2] = B_prod1 *B_proddiff1 *Burst1_off_orig+B_prod1 *Burst1_on_orig;
  //update propensity for A1 = 
  propensities[3] = A_deg1 *A1;
  //update propensity for Aprime1 = 
  propensities[4] = Aprime_deg1 *Aprime1;
  //update propensity for B1 = 
  propensities[5] = B_deg1 *B1;
  //update propensity for Burst1_off_targ = Burst1_on_targ
  propensities[6] = B_ondep1 *(pow(A1,nA1 )/(pow(kA1 ,nA1 )+pow(A1,nA1 )))*Burst1_off_targ+onbasal_a1 *Burst1_off_targ;
  //update propensity for Burst1_off_targ = Burst1_on_targ
  propensities[7] = B_ondep_prime *(pow(A1,nA1 )/(pow(kA1 ,nA1 )+pow(A1,nA1 )))*Burst1_off_para+onbasal_aprime1 *Burst1_off_para;
  //update propensity for Burst1_off_para = Burst1_on_para
  propensities[8] = Aprimenitc1 *(pow(A1,nA1 )/(pow(kA1 ,nA1 )+pow(A1,nA1 )))*Burst1_on_orig+onbasal_b1 *Burst1_on_orig;
  //update propensity for Burst1_on_orig = Burst1_off_orig
  propensities[9] = A_off1 *Burst1_on_targ;
  //update propensity for Burst1_on_para = Burst1_off_para
  propensities[10] = Aprime_off1 *Burst1_on_para;
  //update propensity for Burst1_on_targ = Burst1_off_targ
  propensities[11] = B_off1 *Burst1_off_orig;
} else if (p<cumpropensities[5]) {
  // rxn: B1 = 
  B1=B1 + -1;

  //update propensity for = A1 
  propensities[0] = A_prod1 *A_proddiff1 *Burst1_on_targ+A_prod1 *Burst1_off_targ;
  //update propensity for = Aprime1 
  propensities[1] = Aprime_prod1 *Aprime_proddiff1 *Burst1_on_para+Aprime_prod1 *Burst1_off_para;
  //update propensity for = B1 
  propensities[2] = B_prod1 *B_proddiff1 *Burst1_off_orig+B_prod1 *Burst1_on_orig;
  //update propensity for A1 = 
  propensities[3] = A_deg1 *A1;
  //update propensity for Aprime1 = 
  propensities[4] = Aprime_deg1 *Aprime1;
  //update propensity for B1 = 
  propensities[5] = B_deg1 *B1;
  //update propensity for Burst1_off_targ = Burst1_on_targ
  propensities[6] = B_ondep1 *(pow(A1,nA1 )/(pow(kA1 ,nA1 )+pow(A1,nA1 )))*Burst1_off_targ+onbasal_a1 *Burst1_off_targ;
  //update propensity for Burst1_off_targ = Burst1_on_targ
  propensities[7] = B_ondep_prime *(pow(A1,nA1 )/(pow(kA1 ,nA1 )+pow(A1,nA1 )))*Burst1_off_para+onbasal_aprime1 *Burst1_off_para;
  //update propensity for Burst1_off_para = Burst1_on_para
  propensities[8] = Aprimenitc1 *(pow(A1,nA1 )/(pow(kA1 ,nA1 )+pow(A1,nA1 )))*Burst1_on_orig+onbasal_b1 *Burst1_on_orig;
  //update propensity for Burst1_on_orig = Burst1_off_orig
  propensities[9] = A_off1 *Burst1_on_targ;
  //update propensity for Burst1_on_para = Burst1_off_para
  propensities[10] = Aprime_off1 *Burst1_on_para;
  //update propensity for Burst1_on_targ = Burst1_off_targ
  propensities[11] = B_off1 *Burst1_off_orig;
} else if (p<cumpropensities[6]) {
  // rxn: Burst1_off_targ = Burst1_on_targ
  Burst1_off_targ=Burst1_off_targ + -1;
  Burst1_on_targ=Burst1_on_targ + 1;

  //update propensity for = A1 
  propensities[0] = A_prod1 *A_proddiff1 *Burst1_on_targ+A_prod1 *Burst1_off_targ;
  //update propensity for = Aprime1 
  propensities[1] = Aprime_prod1 *Aprime_proddiff1 *Burst1_on_para+Aprime_prod1 *Burst1_off_para;
  //update propensity for = B1 
  propensities[2] = B_prod1 *B_proddiff1 *Burst1_off_orig+B_prod1 *Burst1_on_orig;
  //update propensity for A1 = 
  propensities[3] = A_deg1 *A1;
  //update propensity for Aprime1 = 
  propensities[4] = Aprime_deg1 *Aprime1;
  //update propensity for B1 = 
  propensities[5] = B_deg1 *B1;
  //update propensity for Burst1_off_targ = Burst1_on_targ
  propensities[6] = B_ondep1 *(pow(A1,nA1 )/(pow(kA1 ,nA1 )+pow(A1,nA1 )))*Burst1_off_targ+onbasal_a1 *Burst1_off_targ;
  //update propensity for Burst1_off_targ = Burst1_on_targ
  propensities[7] = B_ondep_prime *(pow(A1,nA1 )/(pow(kA1 ,nA1 )+pow(A1,nA1 )))*Burst1_off_para+onbasal_aprime1 *Burst1_off_para;
  //update propensity for Burst1_off_para = Burst1_on_para
  propensities[8] = Aprimenitc1 *(pow(A1,nA1 )/(pow(kA1 ,nA1 )+pow(A1,nA1 )))*Burst1_on_orig+onbasal_b1 *Burst1_on_orig;
  //update propensity for Burst1_on_orig = Burst1_off_orig
  propensities[9] = A_off1 *Burst1_on_targ;
  //update propensity for Burst1_on_para = Burst1_off_para
  propensities[10] = Aprime_off1 *Burst1_on_para;
  //update propensity for Burst1_on_targ = Burst1_off_targ
  propensities[11] = B_off1 *Burst1_off_orig;
} else if (p<cumpropensities[7]) {
  // rxn: Burst1_off_targ = Burst1_on_targ
  Burst1_off_targ=Burst1_off_targ + -1;
  Burst1_on_targ=Burst1_on_targ + 1;

  //update propensity for = A1 
  propensities[0] = A_prod1 *A_proddiff1 *Burst1_on_targ+A_prod1 *Burst1_off_targ;
  //update propensity for = Aprime1 
  propensities[1] = Aprime_prod1 *Aprime_proddiff1 *Burst1_on_para+Aprime_prod1 *Burst1_off_para;
  //update propensity for = B1 
  propensities[2] = B_prod1 *B_proddiff1 *Burst1_off_orig+B_prod1 *Burst1_on_orig;
  //update propensity for A1 = 
  propensities[3] = A_deg1 *A1;
  //update propensity for Aprime1 = 
  propensities[4] = Aprime_deg1 *Aprime1;
  //update propensity for B1 = 
  propensities[5] = B_deg1 *B1;
  //update propensity for Burst1_off_targ = Burst1_on_targ
  propensities[6] = B_ondep1 *(pow(A1,nA1 )/(pow(kA1 ,nA1 )+pow(A1,nA1 )))*Burst1_off_targ+onbasal_a1 *Burst1_off_targ;
  //update propensity for Burst1_off_targ = Burst1_on_targ
  propensities[7] = B_ondep_prime *(pow(A1,nA1 )/(pow(kA1 ,nA1 )+pow(A1,nA1 )))*Burst1_off_para+onbasal_aprime1 *Burst1_off_para;
  //update propensity for Burst1_off_para = Burst1_on_para
  propensities[8] = Aprimenitc1 *(pow(A1,nA1 )/(pow(kA1 ,nA1 )+pow(A1,nA1 )))*Burst1_on_orig+onbasal_b1 *Burst1_on_orig;
  //update propensity for Burst1_on_orig = Burst1_off_orig
  propensities[9] = A_off1 *Burst1_on_targ;
  //update propensity for Burst1_on_para = Burst1_off_para
  propensities[10] = Aprime_off1 *Burst1_on_para;
  //update propensity for Burst1_on_targ = Burst1_off_targ
  propensities[11] = B_off1 *Burst1_off_orig;
} else if (p<cumpropensities[8]) {
  // rxn: Burst1_off_para = Burst1_on_para
  Burst1_off_para=Burst1_off_para + -1;
  Burst1_on_para=Burst1_on_para + 1;

  //update propensity for = A1 
  propensities[0] = A_prod1 *A_proddiff1 *Burst1_on_targ+A_prod1 *Burst1_off_targ;
  //update propensity for = Aprime1 
  propensities[1] = Aprime_prod1 *Aprime_proddiff1 *Burst1_on_para+Aprime_prod1 *Burst1_off_para;
  //update propensity for = B1 
  propensities[2] = B_prod1 *B_proddiff1 *Burst1_off_orig+B_prod1 *Burst1_on_orig;
  //update propensity for A1 = 
  propensities[3] = A_deg1 *A1;
  //update propensity for Aprime1 = 
  propensities[4] = Aprime_deg1 *Aprime1;
  //update propensity for B1 = 
  propensities[5] = B_deg1 *B1;
  //update propensity for Burst1_off_targ = Burst1_on_targ
  propensities[6] = B_ondep1 *(pow(A1,nA1 )/(pow(kA1 ,nA1 )+pow(A1,nA1 )))*Burst1_off_targ+onbasal_a1 *Burst1_off_targ;
  //update propensity for Burst1_off_targ = Burst1_on_targ
  propensities[7] = B_ondep_prime *(pow(A1,nA1 )/(pow(kA1 ,nA1 )+pow(A1,nA1 )))*Burst1_off_para+onbasal_aprime1 *Burst1_off_para;
  //update propensity for Burst1_off_para = Burst1_on_para
  propensities[8] = Aprimenitc1 *(pow(A1,nA1 )/(pow(kA1 ,nA1 )+pow(A1,nA1 )))*Burst1_on_orig+onbasal_b1 *Burst1_on_orig;
  //update propensity for Burst1_on_orig = Burst1_off_orig
  propensities[9] = A_off1 *Burst1_on_targ;
  //update propensity for Burst1_on_para = Burst1_off_para
  propensities[10] = Aprime_off1 *Burst1_on_para;
  //update propensity for Burst1_on_targ = Burst1_off_targ
  propensities[11] = B_off1 *Burst1_off_orig;
} else if (p<cumpropensities[9]) {
  // rxn: Burst1_on_orig = Burst1_off_orig
  Burst1_on_orig=Burst1_on_orig + -1;
  Burst1_off_orig=Burst1_off_orig + 1;

  //update propensity for = A1 
  propensities[0] = A_prod1 *A_proddiff1 *Burst1_on_targ+A_prod1 *Burst1_off_targ;
  //update propensity for = Aprime1 
  propensities[1] = Aprime_prod1 *Aprime_proddiff1 *Burst1_on_para+Aprime_prod1 *Burst1_off_para;
  //update propensity for = B1 
  propensities[2] = B_prod1 *B_proddiff1 *Burst1_off_orig+B_prod1 *Burst1_on_orig;
  //update propensity for A1 = 
  propensities[3] = A_deg1 *A1;
  //update propensity for Aprime1 = 
  propensities[4] = Aprime_deg1 *Aprime1;
  //update propensity for B1 = 
  propensities[5] = B_deg1 *B1;
  //update propensity for Burst1_off_targ = Burst1_on_targ
  propensities[6] = B_ondep1 *(pow(A1,nA1 )/(pow(kA1 ,nA1 )+pow(A1,nA1 )))*Burst1_off_targ+onbasal_a1 *Burst1_off_targ;
  //update propensity for Burst1_off_targ = Burst1_on_targ
  propensities[7] = B_ondep_prime *(pow(A1,nA1 )/(pow(kA1 ,nA1 )+pow(A1,nA1 )))*Burst1_off_para+onbasal_aprime1 *Burst1_off_para;
  //update propensity for Burst1_off_para = Burst1_on_para
  propensities[8] = Aprimenitc1 *(pow(A1,nA1 )/(pow(kA1 ,nA1 )+pow(A1,nA1 )))*Burst1_on_orig+onbasal_b1 *Burst1_on_orig;
  //update propensity for Burst1_on_orig = Burst1_off_orig
  propensities[9] = A_off1 *Burst1_on_targ;
  //update propensity for Burst1_on_para = Burst1_off_para
  propensities[10] = Aprime_off1 *Burst1_on_para;
  //update propensity for Burst1_on_targ = Burst1_off_targ
  propensities[11] = B_off1 *Burst1_off_orig;
} else if (p<cumpropensities[10]) {
  // rxn: Burst1_on_para = Burst1_off_para
  Burst1_off_para=Burst1_off_para + 1;
  Burst1_on_para=Burst1_on_para + -1;

  //update propensity for = A1 
  propensities[0] = A_prod1 *A_proddiff1 *Burst1_on_targ+A_prod1 *Burst1_off_targ;
  //update propensity for = Aprime1 
  propensities[1] = Aprime_prod1 *Aprime_proddiff1 *Burst1_on_para+Aprime_prod1 *Burst1_off_para;
  //update propensity for = B1 
  propensities[2] = B_prod1 *B_proddiff1 *Burst1_off_orig+B_prod1 *Burst1_on_orig;
  //update propensity for A1 = 
  propensities[3] = A_deg1 *A1;
  //update propensity for Aprime1 = 
  propensities[4] = Aprime_deg1 *Aprime1;
  //update propensity for B1 = 
  propensities[5] = B_deg1 *B1;
  //update propensity for Burst1_off_targ = Burst1_on_targ
  propensities[6] = B_ondep1 *(pow(A1,nA1 )/(pow(kA1 ,nA1 )+pow(A1,nA1 )))*Burst1_off_targ+onbasal_a1 *Burst1_off_targ;
  //update propensity for Burst1_off_targ = Burst1_on_targ
  propensities[7] = B_ondep_prime *(pow(A1,nA1 )/(pow(kA1 ,nA1 )+pow(A1,nA1 )))*Burst1_off_para+onbasal_aprime1 *Burst1_off_para;
  //update propensity for Burst1_off_para = Burst1_on_para
  propensities[8] = Aprimenitc1 *(pow(A1,nA1 )/(pow(kA1 ,nA1 )+pow(A1,nA1 )))*Burst1_on_orig+onbasal_b1 *Burst1_on_orig;
  //update propensity for Burst1_on_orig = Burst1_off_orig
  propensities[9] = A_off1 *Burst1_on_targ;
  //update propensity for Burst1_on_para = Burst1_off_para
  propensities[10] = Aprime_off1 *Burst1_on_para;
  //update propensity for Burst1_on_targ = Burst1_off_targ
  propensities[11] = B_off1 *Burst1_off_orig;
} else if (p<cumpropensities[11]) {
  // rxn: Burst1_on_targ = Burst1_off_targ
  Burst1_off_targ=Burst1_off_targ + 1;
  Burst1_on_targ=Burst1_on_targ + -1;

  //update propensity for = A1 
  propensities[0] = A_prod1 *A_proddiff1 *Burst1_on_targ+A_prod1 *Burst1_off_targ;
  //update propensity for = Aprime1 
  propensities[1] = Aprime_prod1 *Aprime_proddiff1 *Burst1_on_para+Aprime_prod1 *Burst1_off_para;
  //update propensity for = B1 
  propensities[2] = B_prod1 *B_proddiff1 *Burst1_off_orig+B_prod1 *Burst1_on_orig;
  //update propensity for A1 = 
  propensities[3] = A_deg1 *A1;
  //update propensity for Aprime1 = 
  propensities[4] = Aprime_deg1 *Aprime1;
  //update propensity for B1 = 
  propensities[5] = B_deg1 *B1;
  //update propensity for Burst1_off_targ = Burst1_on_targ
  propensities[6] = B_ondep1 *(pow(A1,nA1 )/(pow(kA1 ,nA1 )+pow(A1,nA1 )))*Burst1_off_targ+onbasal_a1 *Burst1_off_targ;
  //update propensity for Burst1_off_targ = Burst1_on_targ
  propensities[7] = B_ondep_prime *(pow(A1,nA1 )/(pow(kA1 ,nA1 )+pow(A1,nA1 )))*Burst1_off_para+onbasal_aprime1 *Burst1_off_para;
  //update propensity for Burst1_off_para = Burst1_on_para
  propensities[8] = Aprimenitc1 *(pow(A1,nA1 )/(pow(kA1 ,nA1 )+pow(A1,nA1 )))*Burst1_on_orig+onbasal_b1 *Burst1_on_orig;
  //update propensity for Burst1_on_orig = Burst1_off_orig
  propensities[9] = A_off1 *Burst1_on_targ;
  //update propensity for Burst1_on_para = Burst1_off_para
  propensities[10] = Aprime_off1 *Burst1_on_para;
  //update propensity for Burst1_on_targ = Burst1_off_targ
  propensities[11] = B_off1 *Burst1_off_orig;
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
