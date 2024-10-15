#include "pluto.h"

/* *************************************************************** */
void ComputeUserVar (const Data *d, Grid *grid)
/*
 *
 *  PURPOSE
 *
 *    Define user-defined output variables
 *
 *
 *
 ***************************************************************** */
{  
  int  i, j, k;
  static int first_call = 1;
  double *dx1, *dx2, *dx3;
  double *r, *rp, *th, *thp;
  double s, sp;
  double ***Jx1, ***Jx2, ***Jx3;
  double ***J1, ***J2, ***J3;
  double ***DivB, ***DB;
  double dx1_Bx2 = 0.0, dx2_Bx1 = 0.0;
  double dx2_Bx3 = 0.0, dx3_Bx2 = 0.0;
  double dx1_Bx3 = 0.0, dx3_Bx1 = 0.0;
  double d12, d21, d13, d31, d23, d32;
  double h3;
  static double ***Bx1, ***Bx2, ***Bx3;
  static double ***a23Bx3, ***a13Bx3, ***a12Bx2;
  RBox box;

  J1 = GetUserVar("J1");
  J2 = GetUserVar("J2");
  J3 = GetUserVar("J3");
  // DivB = GetUserVar("DivB");


/* --------------------------------------------------------
   1. Set pointer shortcuts.
      The primary magnetic field will be the staggered one
      or the cell-centered field.
      
      Note: in 2.5 dimensions, we also need Jx1 and Jx2 which
            contribute to the total energy equation but not
            to the induction equation.
            For this reason, we set  Bx3 to be equal to the \b
            cell centered value.
   ---------------------------------------------------------------- */

#if (defined STAGGERED_MHD)
  Bx1 = d->Vc[BX1];
  Bx2 = d->Vc[BX2];
  Bx3 = d->Vc[BX3];

  DIM_EXPAND(Bx1 = d->Vs[BX1s];  , 
           Bx2 = d->Vs[BX2s];  ,
           Bx3 = d->Vs[BX3s];)
#else
  if (first_call){
    Bx1 = ARRAY_3D(NX3_MAX,NX2_MAX,NX1_MAX,double);
    Bx2 = ARRAY_3D(NX3_MAX,NX2_MAX,NX1_MAX,double);
    Bx3 = ARRAY_3D(NX3_MAX,NX2_MAX,NX1_MAX,double);
  }

/* -- Construct interface values -- */

  for (k = 0; k < NX3_TOT; k++){
  for (j = 0; j < NX2_TOT; j++){
  for (i = 0; i < NX1_TOT-INCLUDE_IDIR; i++){
    Bx1[k][j][i] = AVERAGE_X(d->Vc[BX1],k,j,i);
  }}}

  for (k = 0; k < NX3_TOT; k++){
  for (j = 0; j < NX2_TOT-INCLUDE_JDIR; j++){
  for (i = 0; i < NX1_TOT; i++){
    Bx2[k][j][i] = AVERAGE_Y(d->Vc[BX2],k,j,i);
  }}}

  for (k = 0; k < NX3_TOT-INCLUDE_KDIR; k++){
  for (j = 0; j < NX2_TOT; j++){
  for (i = 0; i < NX1_TOT; i++){ 
    Bx3[k][j][i] = AVERAGE_Z(d->Vc[BX3],k,j,i);
  }}}

#endif

  Jx1 = d->J[IDIR];
  Jx2 = d->J[JDIR];
  Jx3 = d->J[KDIR];

  dx1 = grid->dx[IDIR];
  dx2 = grid->dx[JDIR];
  dx3 = grid->dx[KDIR];

/* ----------------------------------------------------------------
   2. Allocate static memory areas for the three arrays
      a12Bx2, a23Bx3, a13Bx3 if necessary. Otherwise, they will
      just be shortcuts to the default magnetic field arrays.
   ---------------------------------------------------------------- */

  if (first_call){
    a12Bx2 = Bx2;  
    a23Bx3 = Bx3;
    a13Bx3 = Bx3;
    #if GEOMETRY == CYLINDRICAL
    a13Bx3 = ARRAY_3D(NX3_MAX, NX2_MAX, NX1_MAX, double);
    #elif GEOMETRY == POLAR
    a12Bx2 = ARRAY_3D(NX3_MAX, NX2_MAX, NX1_MAX, double);
    #elif GEOMETRY == SPHERICAL
    a12Bx2 = ARRAY_3D(NX3_MAX, NX2_MAX, NX1_MAX, double); 
    a13Bx3 = ARRAY_3D(NX3_MAX, NX2_MAX, NX1_MAX, double);
    a23Bx3 = ARRAY_3D(NX3_MAX, NX2_MAX, NX1_MAX, double);
    #endif
  }

/* --------------------------------------------------------
   3. Construct goemetrical coefficients and arrays
   -------------------------------------------------------- */

#if GEOMETRY == CYLINDRICAL
  r   = grid->x[IDIR]; rp  = grid->xr[IDIR];
  TOT_LOOP(k,j,i) a13Bx3[k][j][i] = Bx3[k][j][i]*fabs(r[i]);
#elif GEOMETRY == POLAR
  r   = grid->x[IDIR]; rp  = grid->xr[IDIR];
  TOT_LOOP(k,j,i) a12Bx2[k][j][i] = Bx2[k][j][i]*fabs(r[i]);
#elif GEOMETRY == SPHERICAL
  r   = grid->x[IDIR]; rp  = grid->xr[IDIR];
  th  = grid->x[JDIR]; thp = grid->xr[JDIR];
  TOT_LOOP(k,j,i) {
    a12Bx2[k][j][i] = Bx2[k][j][i]*r[i];
    a13Bx3[k][j][i] = Bx3[k][j][i]*r[i];
    a23Bx3[k][j][i] = Bx3[k][j][i]*fabs(sin(th[j]));
  }
#endif

/* --------------------------------------------------------
   4a. Compute the three components of currents at cell
       edges, so that the different components of J have
       different spatial locations:

       Jx  at  (i    , j+1/2, k+1/2)
       Jy  at  (i+1/2, j    , k+1/2)
       Jz  at  (i+1/2, j+1/2, k    )
   -------------------------------------------------------- */

  d21 = d23 = d32 = d31 = 0.0;
  for (k = 0; k < NX3_TOT-INCLUDE_KDIR; k++){
  for (j = 0; j < NX2_TOT-INCLUDE_JDIR; j++){
  for (i = 0; i < NX1_TOT-INCLUDE_IDIR; i++){

    DIM_EXPAND(d12 = d13 = 1.0/dx1[i];  ,
               d21 = d23 = 1.0/dx2[j];  ,
               d31 = d32 = 1.0/dx3[k];)

    #if GEOMETRY == CYLINDRICAL
    d32 = d31 = 0.0;
    d13 = 2.0/(fabs(r[i+1])*r[i+1] - fabs(r[i])*r[i]);  /* = 1.0/(rp*dr) */
    #elif GEOMETRY == POLAR
    DIM_EXPAND(d12  = 2.0/(fabs(r[i+1])*r[i+1] - fabs(r[i])*r[i]); ,  /* = 1.0/(rp*dr) */
             d23 /= r[i];  d21 /= rp[i];                         ,
                                       ; )
    #elif GEOMETRY == SPHERICAL
    s  = fabs(sin(th[j]));
    sp = 0.5*(fabs(sin(th[j])) + fabs(sin(th[j+1]))); 
    DIM_EXPAND(d12 /= rp[i];   d13 /= rp[i];     ,
             d21 /= rp[i];   d23 /= r[i]*sp;   ,
             d32 /= r[i]*sp; d31 /= rp[i]*s;)
    #endif

    Jx1[k][j][i] = FDIFF_X2(a23Bx3,k,j,i)*d23 - FDIFF_X3(   Bx2,k,j,i)*d32;
    Jx2[k][j][i] = FDIFF_X3(   Bx1,k,j,i)*d31 - FDIFF_X1(a13Bx3,k,j,i)*d13;
    Jx3[k][j][i] = FDIFF_X1(a12Bx2,k,j,i)*d12 - FDIFF_X2(Bx1,k,j,i)*d21;

    // DB[k][j][i] = FDIFF_X1(Bx1,k,j,i)*d12;
    // DB[k][j][i] = FDIFF_X1(Bx1,k,j,i)*d12 + FDIFF_X2(Bx2,k,j,i)*d21 + FDIFF_X3(Bx2,k,j,i)*d31;


    
  }}}
  


  // int i, j, k;  
  
  DOM_LOOP(k,j,i){
    // printf("HI");

    // J1[k][j][i] = 4;
    
    J1[k][j][i] = Jx1[k][j][i];
    J2[k][j][i] = Jx2[k][j][i];
    J3[k][j][i] = Jx3[k][j][i];
    // DivB[k][j][i] = DB[k][j][i];
    
 

  }
}
/* ************************************************************* */
void ChangeOutputVar ()
/* 
 *
 * 
 *************************************************************** */
{ 
  Image *image;

#if PARTICLES
  //SetOutputVar ("energy",PARTICLES_FLT_OUTPUT, NO);
//  SetOutputVar ("x1",    PARTICLES_FLT_OUTPUT, NO);
  //SetOutputVar ("vx1",   PARTICLES_FLT_OUTPUT, NO);
#endif

}





