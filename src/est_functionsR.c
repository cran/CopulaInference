#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include<string.h>
#include <R.h>


#define cte = 0.6366198



/*************************************************************/
/* This file contains functions to compute estimators        */
/* of copula-based models for arbitrary data                 */
/*                                                           */
/*                                                           */
/****************************************************** ******/

/* Gauss-Legendre integration with 16 points */


double  gaussint16( double (*f)(double), double *a, double *b)
{
  double x[8], w[8], sum=0.0;
  double aa,bb;
  int j;

  aa = 0.5*(b[0]-a[0]);
  bb = 0.5*(a[0]+b[0]);

  x[0] = 0.095012509837637;
  x[1] = 0.281603550779259;
  x[2] = 0.458016777657227;
  x[3] = 0.617876244402644;
  x[4] = 0.755404408355003;
  x[5] = 0.865631202387832;
  x[6] = 0.944575023073233;
  x[7] = 0.989400934991650;


  w[0] = 0.189450610455069;
  w[1] = 0.182603415044924;
  w[2] = 0.169156519395003;
  w[3] = 0.149595988816577;
  w[4] = 0.124628971255534;
  w[5] = 0.095158511682493;
  w[6] = 0.062253523938648;
  w[7] = 0.027152459411754;


  for(j=0;j<8;j++)
  {
    sum += w[j]*f(bb-aa*x[j]) + w[j]*f(bb+aa*x[j]);

  }
  sum = aa*sum;
  return sum;
}

/* With parameters in the function */
double  gaussint16v2( double (*f)(double, double *), double *cpar, double *a, double *b)
{
  double x[8], w[8], sum=0.0;
  double aa,bb;
  int j;

  aa = 0.5*(b[0]-a[0]);
  bb = 0.5*(a[0]+b[0]);

  x[0] = 0.095012509837637;
  x[1] = 0.281603550779259;
  x[2] = 0.458016777657227;
  x[3] = 0.617876244402644;
  x[4] = 0.755404408355003;
  x[5] = 0.865631202387832;
  x[6] = 0.944575023073233;
  x[7] = 0.989400934991650;


  w[0] = 0.189450610455069;
  w[1] = 0.182603415044924;
  w[2] = 0.169156519395003;
  w[3] = 0.149595988816577;
  w[4] = 0.124628971255534;
  w[5] = 0.095158511682493;
  w[6] = 0.062253523938648;
  w[7] = 0.027152459411754;


  for(j=0;j<8;j++)
  {
    sum += w[j]*f(bb-aa*x[j], cpar) + w[j]*f(bb+aa*x[j],cpar);

  }
  sum = aa*sum;
  return sum;
}


/* Double integral With parameters in the function */
double  gaussint162d( double (*f)(double, double, double *), double *cpar, double *a, double *b)
{
  double x[8], w[8], sum=0.0;
  double aa,bb,xx,ww;
  int j,k;

  aa = 0.5*(b[0]-a[0]);
  bb = 0.5*(a[0]+b[0]);

  x[0] = 0.095012509837637;
  x[1] = 0.281603550779259;
  x[2] = 0.458016777657227;
  x[3] = 0.617876244402644;
  x[4] = 0.755404408355003;
  x[5] = 0.865631202387832;
  x[6] = 0.944575023073233;
  x[7] = 0.989400934991650;


  w[0] = 0.189450610455069;
  w[1] = 0.182603415044924;
  w[2] = 0.169156519395003;
  w[3] = 0.149595988816577;
  w[4] = 0.124628971255534;
  w[5] = 0.095158511682493;
  w[6] = 0.062253523938648;
  w[7] = 0.027152459411754;


  for(j=0;j<8;j++)
  {
    xx = bb-aa*x[j];
    ww = w[j];
    for(k=0;k<8;k++)
      {
        sum += w[k]*f(xx, bb-aa*x[k], cpar) + w[k]*f(xx, bb+aa*x[k],cpar);
      }
    sum *= ww;

  }
  sum *= aa*aa;
  return sum;
}

/* Kendall's tau for Frank copula */

double Debye(double x)
{
  return x/(exp(x)-1.0);
}

void taufrank(double *cpar, double *tau)
{
  double sum =0.0;

  double *a = calloc(1,sizeof(double));
  double *b = calloc(1,sizeof(double));

  a[0] = 0.0;
  b[0] = fabs(cpar[0]);

  if(cpar[0]==0.0){tau[0]=0.0;}
  else{
       sum = gaussint16(Debye, a,b);
       tau[0] = 1.0-4.0/b[0]+4.0/(b[0]*b[0])*sum;

      if(cpar[0]<0.0) tau[0]=-tau[0];
  }
 free(a); free(b);
}

/* Kendall's tau for Joe copula */

double KJ(double x, double *cpar)
{
  return x*log(x)*pow(1.0-x, 2.0*(1.0/cpar[0]-1.0) );
}

void taujoe(double *cpar, double *tau)
{
  double sum =0.0;

  double *a = calloc(1,sizeof(double));
  double *b = calloc(1,sizeof(double));

  a[0] = 0.0;
  b[0] = 1.0;

  if(cpar[0]==1.0){tau[0]= 0.0;}
  else{
    sum = gaussint16v2(KJ, cpar, a,b);
    tau[0] = 1.0+4.0/(cpar[0]*cpar[0])*sum;

    }
  free(a); free(b);
}


/* Kendall's tau for bb6 copula */

double Kbb6(double x, double *cpar)
{
  double theta,  tt1, tt2;

  theta = cpar[0];
  tt1 = -log(1.0-pow(1.0-x,theta));
  tt2 = (1.0-x)*(1.0-pow(1.0-x,-theta));
  return tt1*tt2;
}

void taubb6(double *cpar, double *tau)
{
  double sum =0.0;

  double *a = calloc(1,sizeof(double));
  double *b = calloc(1,sizeof(double));

  a[0] = 0.0;
  b[0] = 1.0;


    sum = gaussint16v2(Kbb6, cpar, a,b);
    tau[0] = 1.0+4.0*sum/(cpar[0]*cpar[1]);

  free(a); free(b);
}
/* Kendall's tau for bb7 copula */

double Kbb7(double x, double *cpar)
{
  double theta, delta, tt0, tt1, tt2;

  theta = cpar[0]; delta = cpar[1];

  tt0 = 1.0-pow(1.0-x,theta);
  tt1 = -1.0+ pow(tt0,-delta);
  tt2 = pow(1.0-x,theta-1.0)*pow(tt0,-1.0-delta);
  return tt1/tt2;
}


void taubb7(double *cpar, double *tau)
{
  double sum =0.0;

  double *a = calloc(1,sizeof(double));
  double *b = calloc(1,sizeof(double));

  a[0] = 0.0;
  b[0] = 1.0;

   sum = gaussint16v2(Kbb7, cpar, a,b);
    tau[0] = 1.0-4.0*sum/(cpar[0]*cpar[1]);


  free(a); free(b);
}
/* Kendall's tau for bb8 copula */

double Kbb8(double x, double *cpar)
{
  double theta, delta, tt0, tt1, tt2;

  theta = cpar[0]; delta = cpar[1];

  tt0 = 1.0-x*delta;
  tt1 = -log( (-1.0 + pow(tt0,theta)/(pow(1.0-delta,theta)-1.0) ));
  tt2 =  tt0*(1.0-pow(tt0,-theta));
  return  tt1*tt2;
}

void taubb8(double *cpar, double *tau)
{
  double sum =0.0;

  double *a = calloc(1,sizeof(double));
  double *b = calloc(1,sizeof(double));

  a[0] = 0.0;
  b[0] = 1.0;


    sum = gaussint16v2(Kbb8, cpar, a,b);
    tau[0] = 1.0+4.0*sum/(cpar[0]*cpar[1]);


  free(a); free(b);
}

void taucopula(int *family_number, int *rotation, double *cpar, double *tau)
{
  int fnumber;
  fnumber = (int)family_number[0];
  /*printf("\n cpar0 = %f\n",cpar[0]);*/
  switch(fnumber)
  {
  case 1:  tau[0] = 0.6366198*asin(cpar[0]); break;
  case 2:  tau[0] = 0.6366198*asin(cpar[0]); break;
  case 3:  tau[0] = cpar[0]/(2.0 + cpar[0]);
                    if(rotation[0]==90 || rotation[0]==270)
                    {tau[0] = -tau[0];} break;
  case 4 : tau[0] = 1.0 - 1.0/cpar[0];
                    if(rotation[0]==90 || rotation[0]==270)
                     {tau[0] = -tau[0];} break;
  case 5 : taufrank(cpar,tau); break;
  case 6 : taujoe(cpar,tau);
            if(rotation[0]==90 || rotation[0]==270)
             {tau[0] = -tau[0];} break;

  case 7 : tau[0]=1.0-2.0/(2.0+cpar[0])/cpar[1];
    if(rotation[0]==90 || rotation[0]==270)
    {tau[0] = -tau[0];} break;

  case 8 : taubb6(cpar,tau);
    if(rotation[0]==90 || rotation[0]==270)
    {tau[0] = -tau[0];} break;

  case 9 : taubb7(cpar,tau);
    if(rotation[0]==90 || rotation[0]==270)
    {tau[0] = -tau[0];} break;

  case 10 : taubb8(cpar,tau);
    if(rotation[0]==90 || rotation[0]==270)
    {tau[0] = -tau[0];} break;
  }
}



/* Kendall's tau and Spearman's rho for arbitrary data*/
void estdep(double *x, double *y, int *n, double *tau, double *rho, double *std, double *Fx, double *Fxm, double *Fy, double *Fym, int *Ix, int *Iy)
{
  int i,j;
  double x0,y0, c1, c2, s1, s2, sum3;
  double n1;
  int sum , sum11, sum12, sum21, sum22,  a11, a21, a12, a22;


  sum  = 0;
  sum3 = 0.0;
  s1 = 0.0;
  s2 = 0.0;
  /*n1 = 1.0*n[0];*/
   n1 = 1.0+n[0];

  for(i=0;i<n[0];i++)
  {
    sum11 = 0;
    sum21 = 0;
    sum12 = 0;
    sum22 = 0;

    x0 = x[i];
    y0 = y[i];

    for(j=0;j<n[0];j++)
    {
      a11 = ( x[j] <= x0) ;
      a12 = ( x[j] < x0)  ;
      a21 = ( y[j] <= y0) ;
      a22 = ( y[j] < y0)  ;
      sum  +=  (a11+a12)*(a21+a22);
      sum11 += a11;
      sum22 += a22;
      sum12 += a12;
      sum21 += a21;

    }
    c1 = ((double)(sum11+sum12))/((double)n[0]) -1.0;
    c2 = ((double)(sum21+sum22))/((double)n[0]) -1.0;
    Fx[i] = ((double)sum11)/((double) n1);
    Fxm[i] = ((double)sum12)/((double) n1);
    Fy[i] = ((double)sum21)/((double) n1);
    Fym[i] = ((double)sum22)/((double) n1);
    Ix[i] = (sum11-sum12-1>0);
    Iy[i] = (sum21-sum22-1>0);
    sum3 += c1*c2;
    s1 +=  c1*c1;
    s2 +=  c2*c2;
  }
  s1 =  s1/((double)n[0]);
  s2 =  s2/((double)n[0]);
  std[0] = sqrt(s1*s2);
  tau[0] = -1.0+ ((double)sum)/((double)n[0]*n[0]);
  rho[0] = ((double)sum3)/((double)n[0]) /std[0];
}






   double maxi(double u, double v)

   {
      if( u> v)
         return u;
      else
         return v;

   }
   int maxint(int u, int v)

   {
      if( u> v)
         return u;
      else
         return v;

   }

   double mini(double u, double v)

   {
      if( u> v)
         return v;
      else
         return u;

   }




void quick_sort(double *t, int lo, int hi)
{
    int i,j;
    double mid;
    double tmp;


    i=lo; j=hi;
    mid = t[  (lo + hi) >> 1 ];


    do
    {
        while (t[i] < mid) i++;
        while (mid < t[j]) j--;

        if (i <= j)
        {

            tmp = t[i];
            t[i++] = t[j];
            t[j--] = tmp;
        }
    } while (i <=j);

    if (lo < j) quick_sort(t, lo, j);
    if (i < hi) quick_sort(t, i, hi);
}





void unique(double *x, int *n, double *values, int *m)
{

int i,n1,k;
double *y = calloc(n[0],sizeof(double));
for(i=0;i<n[0];i++)
y[i]=x[i];
n1 = n[0]-1;

quick_sort(y,0,n1);
k=0;
values[0]=y[0];
for(i=1;i<n[0];i++)
  {
     if(y[i]>values[k])
       {
         k++;
         values[k]=y[i];
       }
  }
m[0]=k+1;
free(y);
}
void prepare_data(double *x, int *n, double *values, int *m, double *Fn, double *fn)
   {
        int i,k, somme;
        double v;
        for(k=0;k<m[0];k++){
            somme=0;
            v = values[k];
            for(i=0;i<n[0];i++){
                    somme += (x[i]<= v);
            }
            Fn[k] = ((double) somme)/((double)n[0]);
        }
        fn[0]=Fn[0];
        for(k=1;k<m[0];k++){
            fn[k] = Fn[k]-Fn[k-1];
        }

   }




   void rank(double *x, double *r, int n)

   {



      int i, j;
      int count;

      for(i=0;i<n;i++)
      {
         count=0;
         for(j=0;j<n;j++)
         {
            if(x[j] <= x[i])
               ++count;
         }
         r[i] = (double)count;

      }
   }



   double mean(double *x, int n)

   {
      int i;
      double sum = 0.0;

      for(i=0;i<n;i++)
         sum += x[i];

      return sum/((double) n);
   }



   double sumBR(double *x, int n)

   {
      int i;
      double s = 0.0;

      for(i=0;i<n;i++)
         s += x[i];


      return s;
   }

double stdev(double *x, int n)

{
  int i;
  double m, sum = 0.0;

  m = mean(x,n);

  for(i=0;i<n;i++)
    sum += (x[i]-m)*(x[i]-m);

  return sqrt(sum/((double) n));
}




double maxvec(double *x, int n)

   {
      int i;
      double y, s;

      s = 0.0;
      for(i=0;i<n;i++)
      {
         y = fabs(x[i]);

         if(s <y)
            s = y;
        /* printf("s = %f\n",s); */


      }

      return s;
   }


   void multvec(double *x, double *y, double *xy, int n)

   {
      int i;


      for(i=0;i<n;i++)
         xy[i] = x[i]*y[i];


   }


void hpla(double *uu, double *vv, int *n, int *cond_var, double *cpar, double *hval)
{

  int j;
  double eta, B, tem0, u, v;

  eta  =cpar[0]-1.0;
  for(j=0;j<n[0];j++)
  {
     u = uu[j];
     v = vv[j];
     B =  1. +eta*(v+u);
     tem0 = sqrt(B*B -4.*cpar[0]*eta*v*u);
     if(cond_var[0]==2)/* derivative with respect to v */
     {
       if(u==0.0){hval[j]=0.0;}
       else if(u==1.0){hval[j]=1.0;}
       else{ hval[j] = 0.5*(1.0-(B-2.0*cpar[0]*u)/tem0);}
     }
     else
       { /* derivative with respect to u */
         if(v==0.0){hval[j]=0.0;}
         else if(v==1.0){hval[j]=1.0;}
         else{ hval[j] = 0.5*(1.0-(B-2.0*cpar[0]*v)/tem0);}
       }
  }
}


void dpla(double *uu, double *vv, int *n, double *cpar, double *pdf)
{
  int j;
  double eta, B, tem0,  u, v;

  eta  =cpar[0]-1.0;
  for(j=0;j<n[0];j++)
  {
    u = uu[j];
    v = vv[j];
    B =  1. +eta*(v+u);
    tem0 = sqrt(B*B -4.*cpar[0]*eta*v*u);
    pdf[j] = cpar[0]*( B - 2*eta*u*v)/(tem0*tem0*tem0);

  }
}

void ppla(double *uu, double *vv, int *n, double *cpar, double *cdf)
{
  int j;
  double eta, B, tem0, tem1, u, v;

  eta  =cpar[0]-1.0;

  for(j=0;j<n[0];j++)
  {
    u = uu[j];
    v = vv[j];
    if(u==0.0 || v==0.0){cdf[j]=0.0;}
    else if(u==1.0){cdf[j]=v;}
    else if(v==1.0){cdf[j]=u;}
    else{
      B =  1. +eta*(v+u);
      tem0 = sqrt(B*B -4.*cpar[0]*eta*v*u);
      tem1 = B + tem0;
      cdf[j] = 2.0*cpar[0]*u*v/tem1;
    }
    }


  }


void taupla(double *cpar, double*tau)
{
  double x[16], w[16];
  double aa, bb, ww, sum=0.0;
  int j,k,nn;
  double *pdf = calloc(16,sizeof(double));
  double *cdf = calloc(16,sizeof(double));
  int *n = calloc(1,sizeof(int));
  double *u = calloc(16,sizeof(double));
  double *v = calloc(16,sizeof(double));
  double *prod = calloc(16,sizeof(double));
  double *prod1 = calloc(16,sizeof(double));
  n[0] = 16;
 nn=16;
  aa = 0.5;
  bb = 0.5;

  x[0] = -0.095012509837637;
  x[1] = -0.281603550779259;
  x[2] = -0.458016777657227;
  x[3] = -0.617876244402644;
  x[4] = -0.755404408355003;
  x[5] = -0.865631202387832;
  x[6] = -0.944575023073233;
  x[7] = -0.989400934991650;
  x[8] = 0.095012509837637;
  x[9] = 0.281603550779259;
  x[10] = 0.458016777657227;
  x[11] = 0.617876244402644;
  x[12] = 0.755404408355003;
  x[13] = 0.865631202387832;
  x[14] = 0.944575023073233;
  x[15] = 0.989400934991650;


  w[0] = 0.189450610455069;
  w[1] = 0.182603415044924;
  w[2] = 0.169156519395003;
  w[3] = 0.149595988816577;
  w[4] = 0.124628971255534;
  w[5] = 0.095158511682493;
  w[6] = 0.062253523938648;
  w[7] = 0.027152459411754;
  w[8] = 0.189450610455069;
  w[9] = 0.182603415044924;
  w[10] = 0.169156519395003;
  w[11] = 0.149595988816577;
  w[12] = 0.124628971255534;
  w[13] = 0.095158511682493;
  w[14] = 0.062253523938648;
  w[15] = 0.027152459411754;

  if(cpar[0]==1.0) tau[0]=0.0;
  else{
     for(j=0;j<16;j++)
      {
        ww=x[j];
        for(k=0;k<16;k++)
        {
          u[k] = bb+aa*ww;
          v[k] = bb+aa*x[k];
        }
        for(k=0;k<16;k++)
        {
          dpla(u,v,n, cpar, pdf);
          ppla(u,v,n, cpar, cdf);
          multvec(pdf,cdf,prod,nn);
          multvec(prod,w,prod1,nn);
        }
          sum += w[j]*sumBR(prod1,nn);
        }

        tau[0] = sum-1.0;
  }
  free(u); free(v); free(prod); free(prod1); free(n); free(pdf); free(cdf);
}


void bi_emp_cdf(double *X1, double *X2, int *n, int *n1, int *n2, double *y1, double *y2, double*cdf)
{
 int i,j,k,l;

 double sum, y1new,y2new;

 l=0;

  for(k=0;k<n2[0];k++)
  {
    y2new = y2[k];

    for(j=0;j<n1[0];j++)
    {
      y1new = y1[j];

      sum=0.0;
      for(i=0;i<n[0];i++)
      {
       if(X1[i] <= y1new &&    X2[i] <= y2new)
        { sum = sum +1.0;}
             }
      cdf[l] = sum/(1.0*n[0]);
      l++;
    }
  }

}

void emp_cdf(double *x, int *n, double *Fx, double *Fxm, int *Ix)
{
  int i,j;
  double x0, n1;
  int sum1, sum2, a1, a2;



  n1 = 1.0+ n[0];

  for(i=0;i<n[0];i++)
  {
    sum1 = 0;
    sum2 = 0;
    x0 = x[i];


    for(j=0;j<n[0];j++)
    {
      a1 = ( x[j] <= x0) ;
      a2 = ( x[j] < x0)  ;
      sum2 += a2;
      sum1 += a1;

    }

    Fx[i] = ((double)sum1)/((double) n1);
    Fxm[i] = ((double)sum2)/((double) n1);
    Ix[i] = (sum1-sum2>1);

  }

}


/* Kendall's tau and Spearman's rho for arbitrary data using margins*/
void est_dep(double *x, double *y,int *n, double *tau, double *rho)
{
  int i,j;
  double x0, y0, s1, s2, c1, c2, sum3;
  int sum , sum11, sum12, sum21, sum22, a11, a21, a12, a22;


  sum  = 0;
  sum3 = 0.0;
  s1 = 0.0;
  s2 = 0.0;

  for(i=0;i<n[0];i++)
  {
    sum11 = 0;
    sum21 = 0;
    sum12 = 0;
    sum22 = 0;

    x0  = x[i];
    y0  = y[i];


    for(j=0;j<n[0];j++)
    {
      a11 = ( x[j]  <= x0) ;
      a12 = ( x[j]  < x0)  ;
      a21 = ( y[j]  <= y0) ;
      a22 = ( y[j]  < y0)  ;
      sum  +=  (a11+a12)*(a21+a22);
      sum11 += a11;
      sum22 += a22;
      sum12 += a12;
      sum21 += a21;

    }
    c1 = ((double)(sum11+sum12))/((double)n[0]) -1.0;
    c2 = ((double)(sum21+sum22))/((double)n[0]) -1.0;
    s1 +=  c1*c1;
    s2 +=  c2*c2;
    sum3 += c1*c2;
    }
  s1 =  s1/((double)n[0]);
  s2 =  s2/((double)n[0]);

  tau[0] = -1.0+ ((double)sum)/((double)n[0]*n[0]);
  rho[0] = ((double)sum3)/((double)n[0]) /sqrt(s1*s2);
}


