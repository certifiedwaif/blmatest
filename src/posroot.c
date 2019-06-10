#include <R.h>
#include <math.h>

#include "posroot.h"

void posroot(double a, double b, double c, double *root, double *status)
{ /* this computes the real roots of a cubic polynomial; in the end, if
     status==1, root stores the nonegative root; if status is not one,
     then status is the total number of nonegative roots and root is
     useless
  */

  int i;
  double Q,R,disc,Q3,A,B,aux,x[3];

  *root = 0.;
  *status=0.;

  Q=(R_pow_di(a,2)-3.*b)/9.;
  R=(2*R_pow_di(a,3)-9*a*b+27.*c)/54.;
  Q3=R_pow_di(Q,3);

  disc=R_pow_di(R,2)-Q3;

  if(disc>=0.){
    if(R>=0) A=-cbrt(R+sqrt(disc));
    else A=-cbrt(R-sqrt(disc));

    if(A==0.) B=0.;
    else B=Q/(A);
    *root=(A+B)-a/3.;
    if(*root>=0) *status=1.;
  }
  else{
    A=acos(R/sqrt(Q3));
    aux= 2. * sqrt(Q);
    x[0]=-aux * cos(A/3.);
    x[1]=-aux * cos((A+4.*asin(1.))/3.);
    x[2]=-aux * cos((A-4.*asin(1.))/3.);
    aux=a/3.;
    for(i=0;i<3;i++) x[i]=x[i]-aux;
    for(i=0;i<3;i++){
      if (x[i]>=0.){
	*status=*status+1.;
	*root=x[i];
      }
    }
  }
}

