/* This program is for a 2 to 3 kinematics*/



#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "cuba.h"
#define pi M_PI

/* The values of mass and energy*/

#define k3 2
#define k4 3
#define k5 4
#define m3 2
#define m4 4
#define m5 5
#define s k3+k4+k5
#define sigma s-k5
#define m_plus m3+m4
#define m_min m3-m4

long double k5max = s/2-(((m3*m3+m4*m4)-m5*m5)/2*s);
long double t = sigma*sigma-(k5*k5-m5*m5);
long double k3max = (sigma*(t+m_plus*m_min)+sqrt((k5*k5-m5*m5)*(t-m_plus*m_plus)*(t-m_min*m_min)))/2*t;
long double k3min = (sigma*(t+m_plus*m_min)-sqrt((k5*k5-m5*m5)*(t-m_plus*m_plus)*(t-m_min*m_min)))/2*t;

/*FILE* fp;
fp = fopen("datafile.dat","w");*/
static int Integrand(const int *ndim, const cubareal xx[],
  const int *ncomp, cubareal ff[], void *userdata) {
  
#define x xx[0]  // represents k5_zero
#define y xx[1]  // represents k3_zero
#define z xx[2]  // represents theta
#define q xx[3]  // represents n 
#define f ff[0]

f= (2*pi*pi*(k5max-m5)*(k3max-k3min)*-1*sin(pi*z))/(8*16*pi*pi*pi*pi);

return 0;
  
  }
  
#define NDIM 5
#define NCOMP 1
#define USERDATA NULL
#define NVEC 1
#define EPSREL 1e-3
#define EPSABS 1e-7
#define VERBOSE 1
#define LAST 4
#define SEED 0
#define MINEVAL 0
#define MAXEVAL 50000

/*-----------------------For Vegas------------------------*/

#define NSTART 1000
#define NINCREASE 500
#define NBATCH 1000
#define GRIDNO 0
#define STATEFILE NULL
#define SPIN NULL
#define RETAINSTATEFILE false

/*----------------------For Suave-------------------------*/

#define NNEW 1000
#define NMIN 2
#define FLATNESS 50

/*----------------------For Divonne-----------------------*/

#define KEY1 47
#define KEY2 1
#define KEY3 1
#define MAXPASS 5
#define BORDER 1e-6
#define MAXCHISQ 10
#define MINDEVIATION .25
#define NGIVEN 0
#define LDXGIVEN NDIM
#define NEXTRA 0


/*---------------------For Cuhre--------------------------*/

#define KEY 0

int main() {
  int comp, nregions, neval, fail;
  cubareal integral[NCOMP], error[NCOMP], prob[NCOMP];
  
  /*printf("Enter the values of m1, m2, m3, m4 and P_in respectively:")
  scanf("%d %d %d %d %d"&m1,&m2,&m3,&m4,&pin);*/
  

#if 1
  printf("-------------------- Vegas test --------------------\n");

  Vegas(NDIM, NCOMP, Integrand, USERDATA, NVEC,
    EPSREL, EPSABS, VERBOSE, SEED,
    MINEVAL, MAXEVAL, NSTART, NINCREASE, NBATCH,
    GRIDNO, STATEFILE, SPIN,
    &neval, &fail, integral, error, prob);

  printf("VEGAS RESULT:\tneval %d\tfail %d\n",
    neval, fail);
  for( comp = 0; comp < NCOMP; ++comp )
    printf("VEGAS RESULT:\t%.8f +- %.8f\tp = %.3f\n",
      (double)integral[comp], (double)error[comp], (double)prob[comp]);
#endif

#if 1
  printf("\n-------------------- Suave test --------------------\n");

  Suave(NDIM, NCOMP, Integrand, USERDATA, NVEC,
    EPSREL, EPSABS, VERBOSE | LAST, SEED,
    MINEVAL, MAXEVAL, NNEW, NMIN, FLATNESS,
    STATEFILE, SPIN,
    &nregions, &neval, &fail, integral, error, prob);

  printf("SUAVE RESULT:\tnregions %d\tneval %d\tfail %d\n",
    nregions, neval, fail);
  for( comp = 0; comp < NCOMP; ++comp )
    printf("SUAVE RESULT:\t%.8f +- %.8f\tp = %.3f\n",
      (double)integral[comp], (double)error[comp], (double)prob[comp]);
#endif

#if 1
  printf("\n------------------- Divonne test -------------------\n");

  Divonne(NDIM, NCOMP, Integrand, USERDATA, NVEC,
    EPSREL, EPSABS, VERBOSE, SEED,
    MINEVAL, MAXEVAL, KEY1, KEY2, KEY3, MAXPASS,
    BORDER, MAXCHISQ, MINDEVIATION,
    NGIVEN, LDXGIVEN, NULL, NEXTRA, NULL,
    STATEFILE, SPIN,
    &nregions, &neval, &fail, integral, error, prob);

  printf("DIVONNE RESULT:\tnregions %d\tneval %d\tfail %d\n",
    nregions, neval, fail);
  for( comp = 0; comp < NCOMP; ++comp )
    printf("DIVONNE RESULT:\t%.8f +- %.8f\tp = %.3f\n",
      (double)integral[comp], (double)error[comp], (double)prob[comp]);
#endif

#if 1
  printf("\n-------------------- Cuhre test --------------------\n");

  Cuhre(NDIM, NCOMP, Integrand, USERDATA, NVEC,
    EPSREL, EPSABS, VERBOSE | LAST,
    MINEVAL, MAXEVAL, KEY,
    STATEFILE, SPIN,
    &nregions, &neval, &fail, integral, error, prob);

  printf("CUHRE RESULT:\tnregions %d\tneval %d\tfail %d\n",
    nregions, neval, fail);
  for( comp = 0; comp < NCOMP; ++comp )
    printf("CUHRE RESULT:\t%.8f +- %.8f\tp = %.3f\n",
      (double)integral[comp], (double)error[comp], (double)prob[comp]);
#endif

  return 0;
}




