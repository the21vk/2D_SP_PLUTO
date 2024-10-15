/* ///////////////////////////////////////////////////////////////////// */
/*!
  \file
  \brief Reconnection test (Harris sheet) in 2D.
 
  In this setup - see also Section 5.3 of Mignone et al., ApJS (2012) 198:7 -
  we reproduce a 2D Harris current sheet with magnetic field profile given by 
  \f[
    B_x(y) = B_0 \tanh(y/l)
  \f] 
  where \c l is the half thickness of the layer.
  The density profile is given by
  \f[ 
    \rho(y) = \rho_0 \cosh^{-2}(y/l) + \rho_{\infty}
  \f]
  We use \f$ \rho_0=1 \f$ and \f$  \rho_{\infty}  = 0.2 \f$,
  following the guidelines of Birn et al., 2001, while
  \c l is user supplied. \n
  In order to achieve equilibrium with the magnetic pressure, 
  the thermal pressure is chosen to be \f$ p = c_s^2 \rho \f$, where
  \f$ c_s^2 = \frac{B_0^2}{2\rho_0} \f$.
  The initial equilibrium is pertubed by an additional magnetic field 
  defined as
  \f[
     \begin{array}{lcl}
      B_x(x,y) &=& \DS -\Psi_0\frac{\pi}{L_y}\cos\left(\frac{2\pi x}{L_x}\right)
                    \sin\left(\frac{\pi y}{L_y}\right),  \\ \noalign{\medskip}
     B_y(x,y) &=& \DS +\Psi_0 \frac{2\pi}{L_x}\sin\left(\frac{2\pi x}{L_x}\right)
                    \cos\left(\frac{\pi y}{L_y}\right).      
     \end{array}
  \f]  
  
  The Lundquist number \f$ S \f$ of a plasma is defined as 
  \f[
    S = \frac{v_A L}{\eta}
  \f]
  where \f$ v_A \f$ is the Alfvén velocity, 
  \f$ v_A = \DS \frac{B}{\sqrt{\rho}}\f$, \f$ L \f$ is a typical lenght scale, 
  and \f$ \eta \f$ the plasma resistivity.
  The reconnection rate \f$\mathcal{E} = \DS \frac{v_{in}}{v_{out}}\f$, with
  \f$ v_{in} \f$ and \f$ v_{out}\f$ the plasma inflow and outflow velocities,
  follows the Sweet-Parker scaling law 
  \f$\mathcal{E} \sim \frac{\delta}{L} \sim \frac{1}{\sqrt{S}}\f$.
  In this example several values of the resitivity \f$ \eta \f$, 
  that correspond to different values of the Lundquist number \f$ S \f$, 
  are provided. 
  The reconnection rate, calculated as the ratio \f$ \frac{\delta}{L} \f$
  (see Mignone et al., 2012) verifies the Sweet-Parker scaling in the range 
  \f$ \eta = 10^{-2} - 10^{-4} \f$ (see the first figure below).

  The input parameters (read from \c pluto.ini) for this test problem are:

  - <tt> g_inputParam[ETA]</tt>:   sets the value of resistivity \f$ \eta \f$;
  - <tt> g_inputParam[WIDTH]</tt>: sets the layer width \c l;
  - <tt> g_inputParam[PSI0]</tt>:  sets the amplitude of perturbation \f$\Psi_0\f$.

  \note
  - Configuration #02 employs a small width (l -\> 0, current-sheet)
    large resistivity (test passes only with the new implementation of 
    the resistive-CT module in PLUTO 4.1. Crash with PLUTO 4.0).
  - Configuratation #09 employs adaptive mesh refinement as in the original
    PLUTO-Chombo paper (ApJS 2012).
   

  \image html fig1.png "Computed Sweet-Parker scaling for different values of eta with a resolution of 512x256."
  \image html vx_rho_plot__137.png "Density map and magnetic field lines for eta = 2.e-3 at t = 137."

  \authors E. Striani (edoardo.striani@iaps.inaf.it)\n
           A. Mignone (mignone@ph.unito.it)
  \date    March 02, 2017
  
  \b Reference
     - "The PLUTO Code for Adaptive Mesh Computations in
        Astrophysical FLuid Dynamics" Mignone et al., ApJS (2012) 198,7
     - "Geospace Environmental Modeling (GEM) magnetic
        reconnection challenge" Birn et al., JGR (2001) 106, 3715
*/
/* ///////////////////////////////////////////////////////////////////// */
#include "pluto.h"

/* ********************************************************************* */
void Init (double *v, double x, double y, double x3)
/*
 *********************************************************************** */
{
  double cs2 = 1, b0 = 1.0, l, Psi0;
  double Lx, Ly, kx, ky;
  double rnd, beta;
  g_isoSoundSpeed = 3.0;
  beta    = g_inputParam[BETA];
  //beta = 5.0;
  l = g_inputParam[WIDTH];

#if PHYSICS == MHD || PHYSICS == RMHD
  v[BX2] = tanh(x/l) * ( 1 + cos(pow(CONST_PI * y / 2 + 0.01,2))) /2;
  
  v[BX1] = ((50 * l * CONST_PI * CONST_PI * y + CONST_PI * l) * sin(pow((CONST_PI * y / 2 + 0.01), 2)) * log(cosh(x/l))) / 200;

//  v[RHO] = 2.0 - v[BX2]*v[BX2]/(2*cs2);
  v[RHO] = 2 ;  /* v[PRS] = 0.5 in the original PLUTO (2012) paper. */
  // v[RHO] = v[PRS]/cs2;
  // if (x < 0){
  //   v[VX1] = 0.05;    
  //   }

  // if (x > 0){
  //   v[VX1] = -0.05;    
  //   }
  
  v[VX1] = 0.0;
  v[VX2] = 0.0;

  #if USE_RANDOM_PERTURBATION == YES
    rnd    = RandomNumber(-1,1);
    v[BX1] = 1.e-2*rnd*exp(-x*x*200.0);
  #endif
  Lx = g_domEnd[IDIR] - g_domBeg[IDIR]; kx = 2*CONST_PI/Lx;
  Ly = g_domEnd[JDIR] - g_domBeg[JDIR]; ky = 2*CONST_PI/Ly;

  Psi0    = g_inputParam[PSI0];
  v[BX2] += 0;
  v[BX1] += Psi0*sin(ky*y+0.1);
  v[BX3]  = 0.0;

  v[AX1] = 0.0;
  v[AX2] = 0.0;
  // v[AX3] = Psi0*cos(2.0*ky*y);

#endif
}

/* ********************************************************************* */
void InitDomain (Data *d, Grid *grid)
/*! 
 * Assign initial condition by looping over the computational domain.
 * Called after the usual Init() function to assign initial conditions
 * on primitive variables.
 * Value assigned here will overwrite those prescribed during Init().
 *
 *
 *********************************************************************** */
{
}

/* ********************************************************************* */
void Analysis (const Data *d, Grid *grid)
/* *********************************************************************** */
{

}
#if PHYSICS == MHD
/* ********************************************************************* */
void BackgroundField (double x1, double x2, double x3, double *B0)
/* *********************************************************************** */
{
   B0[0] = 0.0;
   B0[1] = 0.0;
   B0[2] = 0.0;
}
#endif

/* ********************************************************************* */
void UserDefBoundary (const Data *d, RBox *box, int side, Grid *grid) 
/* *********************************************************************** */
{
  int   i, j, k, nv;
  double  x, y, l;
  double *x1, *x2, *x3;

  x1 = grid->x[IDIR];
  x2 = grid->x[JDIR];
  x3 = grid->x[KDIR];

  l = g_inputParam[WIDTH];


  if (side == X1_BEG){
    BOX_LOOP(box,k,j,i){
      x  = x1[i];
      y  = x2[j];
      // d->Vc[VX1][k][j][i] = 0.01;     
      // d->Vc[VX2][k][j][i] = 0.0;    
      d->Vc[BX2][k][j][i] = tanh(x1[i]/l) * ( 1 + cos(pow(CONST_PI * x2[j] / 2 + 0.01,2))) /2; 
      d->Vc[BX1][k][j][i] = ((50 * l * CONST_PI * CONST_PI * x2[j] + CONST_PI * l) * sin(pow((CONST_PI * x2[j] / 2 + 0.01), 2)) * log(cosh(x1[i]/l))) / 200;
      d->Vc[RHO][k][j][i] = 2;    
    }
  }

  if (side == X1_END){
    BOX_LOOP(box,k,j,i){
      x  = x1[i];
      y  = x2[j];
      // d->Vc[VX1][k][j][i] = -0.01;    
      // d->Vc[VX2][k][j][i] = 0.0;    
      d->Vc[BX2][k][j][i] = tanh(x1[i]/l) * ( 1 + cos(pow(CONST_PI * x2[j] / 2 + 0.01,2))) /2; 
      d->Vc[BX1][k][j][i] = ((50 * l * CONST_PI * CONST_PI * x2[j] + CONST_PI * l) * sin(pow((CONST_PI * x2[j] / 2 + 0.01), 2)) * log(cosh(x1[i]/l))) / 200;
      
      
      // d->Vc[BX1][k][j][i] = d->Vc[BX2][k][j][i] * y / x;
      d->Vc[RHO][k][j][i] = 2;    
    }
  }
}

#if BODY_FORCE != NO
/* ********************************************************************* */
void BodyForceVector(double *v, double *g, double x1, double x2, double x3)
/* *********************************************************************** */
{
  g[IDIR] = 0.0;
  g[JDIR] = 0.0;
  g[KDIR] = 0.0;
}
/* ********************************************************************* */
double BodyForcePotential(double x1, double x2, double x3)
/* *********************************************************************** */
{
  return 0.0;
}
#endif


