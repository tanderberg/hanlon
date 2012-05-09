/* evolve_frw.cpp : Euler, RK4 and symplectic (default) integrators (FRW).
 * 
 * Copyright (C) 2007-2012, Tommy Anderberg.
 * 
 * Part of the HaNLON program package. For more information, visit 
 * http://simplicial.net/hanlon/
 * 
 * This software is distributed under the GNU General Public License (GPL) 3.0
 * with the following attribution requirements (GPL Section 7):
 * 
 * - you agree to retain in this software and in any derivative thereof the 
 *   copyright, author attribution, URL and licensing information provided 
 *   in this file;
 * 
 * - if you use this software or any derivative thereof in work to be 
 *   published or otherwise publicly distributed, you agree to acknowledge 
 *   Tommy Anderberg with the relevant citations.
 * 
 * This software is distributed in the hope that it will be useful, but 
 * WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
 * or FITNESS FOR A PARTICULAR PURPOSE. See the GPL for more details.
 * 
 * You should have received a copy of the GPL with this software. If not, see
 * http://www.gnu.org/licenses/
 * 
 * To contact the author, send email to
 * Tommy dot Anderberg at simplicial dot net
 * (with obvious substitutions of anti-spam address elements).
 */

#include "stdafx.h"
#include "core_frw.cpp"

const double DEFAULT_DAOA_0                 = -0.0001;     /* Default initial da/a */
const double DEFAULT_EPSILON                = 1e-6;        /* Default "negligible" error */
const double DEFAULT_H_0                    = 71;          /* Default initial Hubble constant, in km/s/Mpc */
const double MPC_TO_KM                      = 3.086e19;    /* Conversion factor from Mpc to km */
const double RHO_SIM_TO_MKS                 = 1.407e25;    /* Conversion factor from simulated (Minkowski) rho to kg/m^3, using nu = 246.3 GeV */
const double SEC_TO_SIM                     = 3.741954e26; /* Conversion factor from seconds to simulated (unitless, lattice) time using nu = 246.3 GeV. */
const double G_FACTOR                       = 5.59143e-10; /* 8 pi G_N / 3 in m^3/(kg s^2) */
const int DEFAULT_MAX_ITERATIONS            = 100;         /* Default max number of iterations in leapfrog solver */
const int DEFAULT_REPORT_STEPS              = 0;           /* Report after this number of steps */
const int DEFAULT_SAVE_STEPS                = 50;          /* Save after this number of steps */

#ifdef RK4
TFields dotFields[3][MAX_NSIDE][MAX_NSIDE][MAX_NSIDE];     /* Max allowed number of cells always allocated. Ugly but fast. */
TFields dotMomenta[3][MAX_NSIDE][MAX_NSIDE][MAX_NSIDE];    /* Ditto */
#else
TFields dotMomenta[MAX_NSIDE][MAX_NSIDE][MAX_NSIDE];       /* Ditto */
TFields dotFields[MAX_NSIDE][MAX_NSIDE][MAX_NSIDE];        /* Ditto */
#endif

class TRunLattice : public TLattice
{
   double compute_rho_static_core(TFields *fields, double twoaR2, double twoaR4,
                                  double h_00, double h_01, double h_02, double h_11, double h_12, double h_22, 
                                  double dx_ph_0, double dx_ph_1, double dx_ph_2, double dx_ph_3, 
                                  double dx_A_2, double dx_A_3, 
                                  double dy_ph_0, double dy_ph_1, double dy_ph_2, double dy_ph_3, 
                                  double dy_A_1, double dy_A_3, 
                                  double dz_ph_0, double dz_ph_1, double dz_ph_2, double dz_ph_3, 
                                  double dz_A_1, double dz_A_2,
                                  bool includeLambdaTerm);

   double compute_rho(int i, int j, int k, double aR2, double aR4,
                      TFields fields[MAX_NSIDE][MAX_NSIDE][MAX_NSIDE], 
                      TFields momenta[MAX_NSIDE][MAX_NSIDE][MAX_NSIDE],
                      bool includeLambdaTerm);

public:

#ifdef RK4
   void rk4(double dtn);
#else
   int sv(double dtn, double epsilon, int nMaxIterations);
   void euler(double dtn);
#endif

   double compute_rho(TFields fields[MAX_NSIDE][MAX_NSIDE][MAX_NSIDE], 
                      TFields momenta[MAX_NSIDE][MAX_NSIDE][MAX_NSIDE],
   bool includeLambdaTerm);
   void init_a(double h0, double o_sim, double o_d, double o_v, double rho_tol);
   void evolve_a(double dtn);
};

#ifdef RK4

void TRunLattice::rk4(double dtn)
/* IMPORTANT: Assumes field derivatives have already been computed!  */
{
   int substep_1 = (active + 1) % 3;
   int substep_2 = (substep_1 + 1) % 3;

   double dtHalf = dtn / 2;
   double dtThird = dtn / 3;
   double dtSixth = dtThird / 2;

   // First Euler substep
   
//   compute_derivatives(fields[active]);

#ifdef _OPENMP
#pragma omp parallel
{
#pragma omp for
#endif
   for (int i = 0; i <= nSideM1; i++)
       for (int j = 0; j <= nSideM1; j++)
           for (int k = 0; k <= nSideM1; k++)
           {
              compute_dots(i, j, k, fields[active], &(momenta[active][i][j][k]), dotFields[active][i][j][k], dotMomenta[active][i][j][k]);

              fields[substep_1][i][j][k].ph_0 = fields[active][i][j][k].ph_0 + dtHalf * dotFields[active][i][j][k].ph_0;
              fields[substep_1][i][j][k].ph_1 = fields[active][i][j][k].ph_1 + dtHalf * dotFields[active][i][j][k].ph_1;
              fields[substep_1][i][j][k].ph_2 = fields[active][i][j][k].ph_2 + dtHalf * dotFields[active][i][j][k].ph_2;
              fields[substep_1][i][j][k].ph_3 = fields[active][i][j][k].ph_3 + dtHalf * dotFields[active][i][j][k].ph_3;

              fields[substep_1][i][j][k].A_1 = fields[active][i][j][k].A_1 + dtHalf * dotFields[active][i][j][k].A_1;
              fields[substep_1][i][j][k].A_2 = fields[active][i][j][k].A_2 + dtHalf * dotFields[active][i][j][k].A_2;
              fields[substep_1][i][j][k].A_3 = fields[active][i][j][k].A_3 + dtHalf * dotFields[active][i][j][k].A_3;

              momenta[substep_1][i][j][k].ph_0 = momenta[active][i][j][k].ph_0 + dtHalf * dotMomenta[active][i][j][k].ph_0;
              momenta[substep_1][i][j][k].ph_1 = momenta[active][i][j][k].ph_1 + dtHalf * dotMomenta[active][i][j][k].ph_1;
              momenta[substep_1][i][j][k].ph_2 = momenta[active][i][j][k].ph_2 + dtHalf * dotMomenta[active][i][j][k].ph_2;
              momenta[substep_1][i][j][k].ph_3 = momenta[active][i][j][k].ph_3 + dtHalf * dotMomenta[active][i][j][k].ph_3;

              momenta[substep_1][i][j][k].A_1 = momenta[active][i][j][k].A_1 + dtHalf * dotMomenta[active][i][j][k].A_1;
              momenta[substep_1][i][j][k].A_2 = momenta[active][i][j][k].A_2 + dtHalf * dotMomenta[active][i][j][k].A_2;
              momenta[substep_1][i][j][k].A_3 = momenta[active][i][j][k].A_3 + dtHalf * dotMomenta[active][i][j][k].A_3;
           }
#ifdef _OPENMP
#pragma omp barrier
}
#endif

   // Update time derivatives
   
   compute_derivatives(fields[substep_1]);

#ifdef _OPENMP
#pragma omp parallel
{
#pragma omp for
#endif
   for (int i = 0; i <= nSideM1; i++)
       for (int j = 0; j <= nSideM1; j++)
           for (int k = 0; k <= nSideM1; k++)
           {
               compute_dots(i, j, k, fields[substep_1], &(momenta[substep_1][i][j][k]), dotFields[substep_1][i][j][k], dotMomenta[substep_1][i][j][k]);
           }

#ifdef _OPENMP
#pragma omp barrier
}
#endif

   // Second Euler substep

#ifdef _OPENMP
#pragma omp parallel
{
#pragma omp for
#endif
   for (int i = 0; i <= nSideM1; i++)
       for (int j = 0; j <= nSideM1; j++)
           for (int k = 0; k <= nSideM1; k++)
           {
              fields[substep_1][i][j][k].ph_0 = fields[active][i][j][k].ph_0 + dtHalf * dotFields[substep_1][i][j][k].ph_0;
              fields[substep_1][i][j][k].ph_1 = fields[active][i][j][k].ph_1 + dtHalf * dotFields[substep_1][i][j][k].ph_1;
              fields[substep_1][i][j][k].ph_2 = fields[active][i][j][k].ph_2 + dtHalf * dotFields[substep_1][i][j][k].ph_2;
              fields[substep_1][i][j][k].ph_3 = fields[active][i][j][k].ph_3 + dtHalf * dotFields[substep_1][i][j][k].ph_3;

              fields[substep_1][i][j][k].A_1 = fields[active][i][j][k].A_1 + dtHalf * dotFields[substep_1][i][j][k].A_1;
              fields[substep_1][i][j][k].A_2 = fields[active][i][j][k].A_2 + dtHalf * dotFields[substep_1][i][j][k].A_2;
              fields[substep_1][i][j][k].A_3 = fields[active][i][j][k].A_3 + dtHalf * dotFields[substep_1][i][j][k].A_3;

              momenta[substep_1][i][j][k].ph_0 = momenta[active][i][j][k].ph_0 + dtHalf * dotMomenta[substep_1][i][j][k].ph_0;
              momenta[substep_1][i][j][k].ph_1 = momenta[active][i][j][k].ph_1 + dtHalf * dotMomenta[substep_1][i][j][k].ph_1;
              momenta[substep_1][i][j][k].ph_2 = momenta[active][i][j][k].ph_2 + dtHalf * dotMomenta[substep_1][i][j][k].ph_2;
              momenta[substep_1][i][j][k].ph_3 = momenta[active][i][j][k].ph_3 + dtHalf * dotMomenta[substep_1][i][j][k].ph_3;

              momenta[substep_1][i][j][k].A_1 = momenta[active][i][j][k].A_1 + dtHalf * dotMomenta[substep_1][i][j][k].A_1;
              momenta[substep_1][i][j][k].A_2 = momenta[active][i][j][k].A_2 + dtHalf * dotMomenta[substep_1][i][j][k].A_2;
              momenta[substep_1][i][j][k].A_3 = momenta[active][i][j][k].A_3 + dtHalf * dotMomenta[substep_1][i][j][k].A_3;
           }
#ifdef _OPENMP
#pragma omp barrier
}
#endif

   // Update time derivatives
   
   compute_derivatives(fields[substep_1]);

#ifdef _OPENMP
#pragma omp parallel
{
#pragma omp for
#endif
   for (int i = 0; i <= nSideM1; i++)
       for (int j = 0; j <= nSideM1; j++)
           for (int k = 0; k <= nSideM1; k++)
           {
               compute_dots(i, j, k, fields[substep_1], &(momenta[substep_1][i][j][k]), dotFields[substep_2][i][j][k], dotMomenta[substep_2][i][j][k]);
           }

#ifdef _OPENMP
#pragma omp barrier
}
#endif

   // Third Euler substep

#ifdef _OPENMP
#pragma omp parallel
{
#pragma omp for
#endif
   for (int i = 0; i <= nSideM1; i++)
       for (int j = 0; j <= nSideM1; j++)
           for (int k = 0; k <= nSideM1; k++)
           {
              fields[substep_1][i][j][k].ph_0 = fields[active][i][j][k].ph_0 + dtHalf * dotFields[substep_2][i][j][k].ph_0;
              fields[substep_1][i][j][k].ph_1 = fields[active][i][j][k].ph_1 + dtHalf * dotFields[substep_2][i][j][k].ph_1;
              fields[substep_1][i][j][k].ph_2 = fields[active][i][j][k].ph_2 + dtHalf * dotFields[substep_2][i][j][k].ph_2;
              fields[substep_1][i][j][k].ph_3 = fields[active][i][j][k].ph_3 + dtHalf * dotFields[substep_2][i][j][k].ph_3;

              fields[substep_1][i][j][k].A_1 = fields[active][i][j][k].A_1 + dtHalf * dotFields[substep_2][i][j][k].A_1;
              fields[substep_1][i][j][k].A_2 = fields[active][i][j][k].A_2 + dtHalf * dotFields[substep_2][i][j][k].A_2;
              fields[substep_1][i][j][k].A_3 = fields[active][i][j][k].A_3 + dtHalf * dotFields[substep_2][i][j][k].A_3;

              momenta[substep_1][i][j][k].ph_0 = momenta[active][i][j][k].ph_0 + dtHalf * dotMomenta[substep_2][i][j][k].ph_0;
              momenta[substep_1][i][j][k].ph_1 = momenta[active][i][j][k].ph_1 + dtHalf * dotMomenta[substep_2][i][j][k].ph_1;
              momenta[substep_1][i][j][k].ph_2 = momenta[active][i][j][k].ph_2 + dtHalf * dotMomenta[substep_2][i][j][k].ph_2;
              momenta[substep_1][i][j][k].ph_3 = momenta[active][i][j][k].ph_3 + dtHalf * dotMomenta[substep_2][i][j][k].ph_3;

              momenta[substep_1][i][j][k].A_1 = momenta[active][i][j][k].A_1 + dtHalf * dotMomenta[substep_2][i][j][k].A_1;
              momenta[substep_1][i][j][k].A_2 = momenta[active][i][j][k].A_2 + dtHalf * dotMomenta[substep_2][i][j][k].A_2;
              momenta[substep_1][i][j][k].A_3 = momenta[active][i][j][k].A_3 + dtHalf * dotMomenta[substep_2][i][j][k].A_3;
           }
#ifdef _OPENMP
#pragma omp barrier
}
#endif

   // Final Euler step
   
   compute_derivatives(fields[substep_1]);

#ifdef _OPENMP
#pragma omp parallel
{
#pragma omp for
#endif
   for (int i = 0; i <= nSideM1; i++)
       for (int j = 0; j <= nSideM1; j++)
           for (int k = 0; k <= nSideM1; k++)
           {
              TFields df, dm;
              compute_dots(i, j, k, fields[substep_1], &(momenta[substep_1][i][j][k]), df, dm);
           
              fields[active][i][j][k].ph_0 += dtSixth * (dotFields[active][i][j][k].ph_0 + df.ph_0) +
                                              dtThird * (dotFields[substep_1][i][j][k].ph_0 + dotFields[substep_2][i][j][k].ph_0);
              fields[active][i][j][k].ph_1 += dtSixth * (dotFields[active][i][j][k].ph_1 + df.ph_1) +
                                              dtThird * (dotFields[substep_1][i][j][k].ph_1 + dotFields[substep_2][i][j][k].ph_1);
              fields[active][i][j][k].ph_2 += dtSixth * (dotFields[active][i][j][k].ph_2 + df.ph_2) +
                                              dtThird * (dotFields[substep_1][i][j][k].ph_2 + dotFields[substep_2][i][j][k].ph_2);
              fields[active][i][j][k].ph_3 += dtSixth * (dotFields[active][i][j][k].ph_3 + df.ph_3) +
                                              dtThird * (dotFields[substep_1][i][j][k].ph_3 + dotFields[substep_2][i][j][k].ph_3);

              fields[active][i][j][k].A_1 += dtSixth * (dotFields[active][i][j][k].A_1 + df.A_1) +
                                             dtThird * (dotFields[substep_1][i][j][k].A_1 + dotFields[substep_2][i][j][k].A_1);
              fields[active][i][j][k].A_2 += dtSixth * (dotFields[active][i][j][k].A_2 + df.A_2) +
                                             dtThird * (dotFields[substep_1][i][j][k].A_2 + dotFields[substep_2][i][j][k].A_2);
              fields[active][i][j][k].A_3 += dtSixth * (dotFields[active][i][j][k].A_3 + df.A_3) +
                                             dtThird * (dotFields[substep_1][i][j][k].A_3 + dotFields[substep_2][i][j][k].A_3);

              momenta[active][i][j][k].ph_0 += dtSixth * (dotMomenta[active][i][j][k].ph_0 + dm.ph_0) +
                                               dtThird * (dotMomenta[substep_1][i][j][k].ph_0 + dotMomenta[substep_2][i][j][k].ph_0);
              momenta[active][i][j][k].ph_1 += dtSixth * (dotMomenta[active][i][j][k].ph_1 + dm.ph_1) +
                                               dtThird * (dotMomenta[substep_1][i][j][k].ph_1 + dotMomenta[substep_2][i][j][k].ph_1);
              momenta[active][i][j][k].ph_2 += dtSixth * (dotMomenta[active][i][j][k].ph_2 + dm.ph_2) +
                                               dtThird * (dotMomenta[substep_1][i][j][k].ph_2 + dotMomenta[substep_2][i][j][k].ph_2);
              momenta[active][i][j][k].ph_3 += dtSixth * (dotMomenta[active][i][j][k].ph_3 + dm.ph_3) +
                                               dtThird * (dotMomenta[substep_1][i][j][k].ph_3 + dotMomenta[substep_2][i][j][k].ph_3);

              momenta[active][i][j][k].A_1 += dtSixth * (dotMomenta[active][i][j][k].A_1 + dm.A_1) +
                                              dtThird * (dotMomenta[substep_1][i][j][k].A_1 + dotMomenta[substep_2][i][j][k].A_1);
              momenta[active][i][j][k].A_2 += dtSixth * (dotMomenta[active][i][j][k].A_2 + dm.A_2) +
                                              dtThird * (dotMomenta[substep_1][i][j][k].A_2 + dotMomenta[substep_2][i][j][k].A_2);
              momenta[active][i][j][k].A_3 += dtSixth * (dotMomenta[active][i][j][k].A_3 + dm.A_3) +
                                              dtThird * (dotMomenta[substep_1][i][j][k].A_3 + dotMomenta[substep_2][i][j][k].A_3);
           }
#ifdef _OPENMP
#pragma omp barrier
}
#endif

   // Done!
}

#else

int TRunLattice::sv(double dtn, double epsilon, int nMaxIterations)
/* IMPORTANT: Assumes field derivatives have already been computed! */
{
   int scratch = (active + 1) % 3;
   int inActive = (scratch + 1) % 3;

   int err = 0;

   // Precompute auxiliary data
   
   double dtnHalf = dtn / 2;
   double dtnHalfOaR3 = dtnHalf / (a*a*a);
   
//   compute_derivatives(fields[active]);

   // Initial scratch fields guess: forward Euler

#ifdef _OPENMP
#pragma omp parallel
{
#pragma omp for
#endif
   for (int i = 0; i <= nSideM1; i++)
       for (int j = 0; j <= nSideM1; j++)
           for (int k = 0; k <= nSideM1; k++)
           {
              TFields df;
              compute_dotFields(&(fields[active][i][j][k]), &(momenta[active][i][j][k]), df);

              fields[scratch][i][j][k].ph_0 = fields[active][i][j][k].ph_0 + dtnHalf * df.ph_0;
              fields[scratch][i][j][k].ph_1 = fields[active][i][j][k].ph_1 + dtnHalf * df.ph_1;
              fields[scratch][i][j][k].ph_2 = fields[active][i][j][k].ph_2 + dtnHalf * df.ph_2;
              fields[scratch][i][j][k].ph_3 = fields[active][i][j][k].ph_3 + dtnHalf * df.ph_3;
              
              fields[scratch][i][j][k].A_1 = fields[active][i][j][k].A_1 + dtnHalf * df.A_1;
              fields[scratch][i][j][k].A_2 = fields[active][i][j][k].A_2 + dtnHalf * df.A_2;
              fields[scratch][i][j][k].A_3 = fields[active][i][j][k].A_3 + dtnHalf * df.A_3;
           }
#ifdef _OPENMP
#pragma omp barrier
}
#endif

#ifndef _OPENMP
   double local_max_rel_dy;
#else
   double global_max_rel_dy;
#endif

   // Implicit scratch fields solution
   
   int iteration = 0;

   do
   {
      // Precompute auxiliary data
       
      compute_derivatives(fields[scratch]);

#ifdef _OPENMP
      global_max_rel_dy = 0.0;

#pragma omp parallel
{
      double local_max_rel_dy = 0.0;
#else
      local_max_rel_dy = 0.0;
#endif

      TFields df, f;
      double A[NCOMPS][NCOMPS];
      double b[NCOMPS];
      int indx[NCOMPS];
      double deviation;

#ifdef _OPENMP
#pragma omp for
#endif
      for (int i = 0; i <= nSideM1; i++)
          for (int j = 0; j <= nSideM1; j++)
              for (int k = 0; k <= nSideM1; k++)
              {
                 compute_dotFields_and_dDotFields_dFields_aR3(&(fields[scratch][i][j][k]), &(momenta[active][i][j][k]), df, A);

                 A[0][0] *= -dtnHalfOaR3; A[0][0] += 1;
                 A[0][1] *= -dtnHalfOaR3;
                 A[0][2] *= -dtnHalfOaR3;
                 A[0][3] *= -dtnHalfOaR3;
                 A[0][4] *= -dtnHalfOaR3;
                 A[0][5] *= -dtnHalfOaR3;
                 A[0][6] *= -dtnHalfOaR3;

                 A[1][0] *= -dtnHalfOaR3;
                 A[1][1] *= -dtnHalfOaR3; A[1][1] += 1;
                 A[1][2] *= -dtnHalfOaR3;
                 A[1][3] *= -dtnHalfOaR3;
                 A[1][4] *= -dtnHalfOaR3;
                 A[1][5] *= -dtnHalfOaR3;
                 A[1][6] *= -dtnHalfOaR3;

                 A[2][0] *= -dtnHalfOaR3;
                 A[2][1] *= -dtnHalfOaR3;
                 A[2][2] *= -dtnHalfOaR3; A[2][2] += 1;
                 A[2][3] *= -dtnHalfOaR3;
                 A[2][4] *= -dtnHalfOaR3;
                 A[2][5] *= -dtnHalfOaR3;
                 A[2][6] *= -dtnHalfOaR3;

                 A[3][0] *= -dtnHalfOaR3;
                 A[3][1] *= -dtnHalfOaR3;
                 A[3][2] *= -dtnHalfOaR3;
                 A[3][3] *= -dtnHalfOaR3; A[3][3] += 1;
                 A[3][4] *= -dtnHalfOaR3;
                 A[3][5] *= -dtnHalfOaR3;
                 A[3][6] *= -dtnHalfOaR3;

                 A[4][4] = 1;
                 A[5][5] = 1;
                 A[6][6] = 1;
                   
                 b[0] = (fields[scratch][i][j][k].ph_0 - fields[active][i][j][k].ph_0) - dtnHalf * df.ph_0;
                 b[1] = (fields[scratch][i][j][k].ph_1 - fields[active][i][j][k].ph_1) - dtnHalf * df.ph_1;
                 b[2] = (fields[scratch][i][j][k].ph_2 - fields[active][i][j][k].ph_2) - dtnHalf * df.ph_2;
                 b[3] = (fields[scratch][i][j][k].ph_3 - fields[active][i][j][k].ph_3) - dtnHalf * df.ph_3;
                   
                 b[4] = (fields[scratch][i][j][k].A_1 - fields[active][i][j][k].A_1) - dtnHalf * df.A_1;
                 b[5] = (fields[scratch][i][j][k].A_2 - fields[active][i][j][k].A_2) - dtnHalf * df.A_2;
                 b[6] = (fields[scratch][i][j][k].A_3 - fields[active][i][j][k].A_3) - dtnHalf * df.A_3;

                 solve7special(&(A[0][0]), &(indx[0]), &(b[0]));

                 f.ph_0 = fields[scratch][i][j][k].ph_0 - b[0];
                 f.ph_1 = fields[scratch][i][j][k].ph_1 - b[1];
                 f.ph_2 = fields[scratch][i][j][k].ph_2 - b[2];
                 f.ph_3 = fields[scratch][i][j][k].ph_3 - b[3];

                 f.A_1 = fields[scratch][i][j][k].A_1 - b[4];
                 f.A_2 = fields[scratch][i][j][k].A_2 - b[5];
                 f.A_3 = fields[scratch][i][j][k].A_3 - b[6];

                 deviation = fabs(b[0]) / (fabs(f.ph_0) + 1); if (local_max_rel_dy < deviation) local_max_rel_dy = deviation;
                 deviation = fabs(b[1]) / (fabs(f.ph_1) + 1); if (local_max_rel_dy < deviation) local_max_rel_dy = deviation;
                 deviation = fabs(b[2]) / (fabs(f.ph_2) + 1); if (local_max_rel_dy < deviation) local_max_rel_dy = deviation;
                 deviation = fabs(b[3]) / (fabs(f.ph_3) + 1); if (local_max_rel_dy < deviation) local_max_rel_dy = deviation;

                 deviation = fabs(b[4]) / (fabs(f.A_1) + 1); if (local_max_rel_dy < deviation) local_max_rel_dy = deviation;
                 deviation = fabs(b[5]) / (fabs(f.A_2) + 1); if (local_max_rel_dy < deviation) local_max_rel_dy = deviation;
                 deviation = fabs(b[6]) / (fabs(f.A_3) + 1); if (local_max_rel_dy < deviation) local_max_rel_dy = deviation;

                 fields[inActive][i][j][k] = f;
              }

#ifdef _OPENMP
#pragma omp barrier
#pragma omp critical
{
      if (local_max_rel_dy > global_max_rel_dy) global_max_rel_dy = local_max_rel_dy;
}
#pragma omp barrier
}
      // Leave updated fields in scratch array

      swap(scratch, inActive);

   } while ((global_max_rel_dy > epsilon) && (++iteration < nMaxIterations));
   
   if (global_max_rel_dy > epsilon) err = 1;

#else

      // Leave updated fields in scratch array

      swap(scratch, inActive);
     
   } while ((local_max_rel_dy > epsilon) && (++iteration < nMaxIterations));
   
   if (local_max_rel_dy > epsilon) err = 1;
   
#endif

   // Precompute auxiliary data and intermediate dot momenta; initialize momenta by forward Euler
   
   compute_derivatives(fields[scratch]);

#ifdef _OPENMP
#pragma omp parallel
{
#pragma omp for
#endif
   for (int i = 0; i <= nSideM1; i++)
       for (int j = 0; j <= nSideM1; j++)
           for (int k = 0; k <= nSideM1; k++)
           {
              compute_dotMomenta(i, j, k, fields[scratch], &(momenta[active][i][j][k]), dotMomenta[i][j][k]);

              momenta[scratch][i][j][k].ph_0 = momenta[active][i][j][k].ph_0 + dtn * dotMomenta[i][j][k].ph_0;
              momenta[scratch][i][j][k].ph_1 = momenta[active][i][j][k].ph_1 + dtn * dotMomenta[i][j][k].ph_1;
              momenta[scratch][i][j][k].ph_2 = momenta[active][i][j][k].ph_2 + dtn * dotMomenta[i][j][k].ph_2;
              momenta[scratch][i][j][k].ph_3 = momenta[active][i][j][k].ph_3 + dtn * dotMomenta[i][j][k].ph_3;
              
              momenta[scratch][i][j][k].A_1 = momenta[active][i][j][k].A_1 + dtn * dotMomenta[i][j][k].A_1;
              momenta[scratch][i][j][k].A_2 = momenta[active][i][j][k].A_2 + dtn * dotMomenta[i][j][k].A_2;
              momenta[scratch][i][j][k].A_3 = momenta[active][i][j][k].A_3 + dtn * dotMomenta[i][j][k].A_3;
           }

#ifdef _OPENMP
#pragma omp barrier
}
#endif

   // Implicit momenta solution

   int from_p = scratch;
   int to_p = inActive;
   
   iteration = 0;
   
   do
   {
#ifdef _OPENMP
      global_max_rel_dy = 0.0;

#pragma omp parallel
{
      double local_max_rel_dy = 0.0;
#else
      local_max_rel_dy = 0.0;
#endif

      TFields dm, m;
      double A[NCOMPS][NCOMPS];
      double b[NCOMPS];
      int indx[NCOMPS];
      double deviation;

#ifdef _OPENMP
#pragma omp for
#endif
      for (int i = 0; i <= nSideM1; i++)
          for (int j = 0; j <= nSideM1; j++)
              for (int k = 0; k <= nSideM1; k++)
              {
                 compute_dotMomenta_and_dDotFields_dFields_aR3(i, j, k, fields[scratch], &(momenta[from_p][i][j][k]), dm, A);

                 A[0][0] *= dtnHalfOaR3; A[0][0] += 1;
                 A[0][1] *= dtnHalfOaR3;
                 A[0][2] *= dtnHalfOaR3;
                 A[0][3] *= dtnHalfOaR3;
                 A[0][4] *= dtnHalfOaR3;
                 A[0][5] *= dtnHalfOaR3;
                 A[0][6] *= dtnHalfOaR3;

                 A[1][0] *= dtnHalfOaR3;
                 A[1][1] *= dtnHalfOaR3; A[1][1] += 1;
                 A[1][2] *= dtnHalfOaR3;
                 A[1][3] *= dtnHalfOaR3;
                 A[1][4] *= dtnHalfOaR3;
                 A[1][5] *= dtnHalfOaR3;
                 A[1][6] *= dtnHalfOaR3;

                 A[2][0] *= dtnHalfOaR3;
                 A[2][1] *= dtnHalfOaR3;
                 A[2][2] *= dtnHalfOaR3; A[2][2] += 1;
                 A[2][3] *= dtnHalfOaR3;
                 A[2][4] *= dtnHalfOaR3;
                 A[2][5] *= dtnHalfOaR3;
                 A[2][6] *= dtnHalfOaR3;

                 A[3][0] *= dtnHalfOaR3;
                 A[3][1] *= dtnHalfOaR3;
                 A[3][2] *= dtnHalfOaR3;
                 A[3][3] *= dtnHalfOaR3; A[3][3] += 1;
                 A[3][4] *= dtnHalfOaR3;
                 A[3][5] *= dtnHalfOaR3;
                 A[3][6] *= dtnHalfOaR3;

                 A[4][4] = 1;
                 A[5][5] = 1;
                 A[6][6] = 1;
 
                 b[0] = (momenta[from_p][i][j][k].ph_0 - momenta[active][i][j][k].ph_0) - dtnHalf * (dotMomenta[i][j][k].ph_0 + dm.ph_0);
                 b[1] = (momenta[from_p][i][j][k].ph_1 - momenta[active][i][j][k].ph_1) - dtnHalf * (dotMomenta[i][j][k].ph_1 + dm.ph_1);
                 b[2] = (momenta[from_p][i][j][k].ph_2 - momenta[active][i][j][k].ph_2) - dtnHalf * (dotMomenta[i][j][k].ph_2 + dm.ph_2);
                 b[3] = (momenta[from_p][i][j][k].ph_3 - momenta[active][i][j][k].ph_3) - dtnHalf * (dotMomenta[i][j][k].ph_3 + dm.ph_3);
                   
                 b[4] = (momenta[from_p][i][j][k].A_1 - momenta[active][i][j][k].A_1) - dtnHalf * (dotMomenta[i][j][k].A_1 + dm.A_1);
                 b[5] = (momenta[from_p][i][j][k].A_2 - momenta[active][i][j][k].A_2) - dtnHalf * (dotMomenta[i][j][k].A_2 + dm.A_2);
                 b[6] = (momenta[from_p][i][j][k].A_3 - momenta[active][i][j][k].A_3) - dtnHalf * (dotMomenta[i][j][k].A_3 + dm.A_3);

                 solve7special(&(A[0][0]), &(indx[0]), &(b[0]));

                 m.ph_0 = momenta[from_p][i][j][k].ph_0 - b[0];
                 m.ph_1 = momenta[from_p][i][j][k].ph_1 - b[1];
                 m.ph_2 = momenta[from_p][i][j][k].ph_2 - b[2];
                 m.ph_3 = momenta[from_p][i][j][k].ph_3 - b[3];

                 m.A_1 = momenta[from_p][i][j][k].A_1 - b[4];
                 m.A_2 = momenta[from_p][i][j][k].A_2 - b[5];
                 m.A_3 = momenta[from_p][i][j][k].A_3 - b[6];

                 deviation = fabs(b[0]) / (fabs(m.ph_0) + 1); if (local_max_rel_dy < deviation) local_max_rel_dy = deviation;
                 deviation = fabs(b[1]) / (fabs(m.ph_1) + 1); if (local_max_rel_dy < deviation) local_max_rel_dy = deviation;
                 deviation = fabs(b[2]) / (fabs(m.ph_2) + 1); if (local_max_rel_dy < deviation) local_max_rel_dy = deviation;
                 deviation = fabs(b[3]) / (fabs(m.ph_3) + 1); if (local_max_rel_dy < deviation) local_max_rel_dy = deviation;

                 deviation = fabs(b[4]) / (fabs(m.A_1) + 1); if (local_max_rel_dy < deviation) local_max_rel_dy = deviation;
                 deviation = fabs(b[5]) / (fabs(m.A_2) + 1); if (local_max_rel_dy < deviation) local_max_rel_dy = deviation;
                 deviation = fabs(b[6]) / (fabs(m.A_3) + 1); if (local_max_rel_dy < deviation) local_max_rel_dy = deviation;
                   
                 momenta[to_p][i][j][k] = m;
              }

#ifdef _OPENMP
#pragma omp barrier
#pragma omp critical
{
      if (local_max_rel_dy > global_max_rel_dy) global_max_rel_dy = local_max_rel_dy;
}
#pragma omp barrier
}

      swap(from_p, to_p);

   } while ((global_max_rel_dy > epsilon) && (++iteration < nMaxIterations));
   
   if (global_max_rel_dy > epsilon) err |= 2;

#else

      swap(from_p, to_p);

   } while ((local_max_rel_dy > epsilon) && (++iteration < nMaxIterations));
   
   if (local_max_rel_dy > epsilon) err |= 2;
   
#endif

   // Make sure final momenta are in scratch array

   if (scratch != from_p)
   {
      memcpy(&(momenta[scratch][0][0][0]), &(momenta[from_p][0][0][0]), (nSideM1+1) * (nSideM1+1) * (nSideM1+1) * sizeof(double) * NCOMPS);
   }

   // Compute final field values in scratch array

#ifdef _OPENMP
#pragma omp parallel
{
#pragma omp for
#endif
   for (int i = 0; i <= nSideM1; i++)
       for (int j = 0; j <= nSideM1; j++)
           for (int k = 0; k <= nSideM1; k++)
           {
              TFields df;
              compute_dotFields(&(fields[scratch][i][j][k]), &(momenta[scratch][i][j][k]), df);

              fields[scratch][i][j][k].ph_0 += dtnHalf * df.ph_0;
              fields[scratch][i][j][k].ph_1 += dtnHalf * df.ph_1;
              fields[scratch][i][j][k].ph_2 += dtnHalf * df.ph_2;
              fields[scratch][i][j][k].ph_3 += dtnHalf * df.ph_3;

              fields[scratch][i][j][k].A_1 += dtnHalf * df.A_1;
              fields[scratch][i][j][k].A_2 += dtnHalf * df.A_2;
              fields[scratch][i][j][k].A_3 += dtnHalf * df.A_3;
           }
#ifdef _OPENMP
#pragma omp barrier
}
#endif

   // Done, set scratch arrays to active

   active = scratch;

   return(err);
}

void TRunLattice::euler(double dtn)
/* IMPORTANT: Assumes field derivatives have already been computed! */
{
   int inActive;

   if (0 == active) inActive = 1; else inActive = 0;

   // Forward Euler
   
//   compute_derivatives(fields[active]);

#ifdef _OPENMP
#pragma omp parallel
{
#pragma omp for
#endif
   for (int i = 0; i <= nSideM1; i++)
       for (int j = 0; j <= nSideM1; j++)
           for (int k = 0; k <= nSideM1; k++)
           {
              TFields df, dm;
              compute_dots(i, j, k, fields[active], &(momenta[active][i][j][k]), df, dm);
           
              fields[inActive][i][j][k].ph_0 = fields[active][i][j][k].ph_0 + dtn * df.ph_0;
              fields[inActive][i][j][k].ph_1 = fields[active][i][j][k].ph_1 + dtn * df.ph_1;
              fields[inActive][i][j][k].ph_2 = fields[active][i][j][k].ph_2 + dtn * df.ph_2;
              fields[inActive][i][j][k].ph_3 = fields[active][i][j][k].ph_3 + dtn * df.ph_3;

              fields[inActive][i][j][k].A_1 = fields[active][i][j][k].A_1 + dtn * df.A_1;
              fields[inActive][i][j][k].A_2 = fields[active][i][j][k].A_2 + dtn * df.A_2;
              fields[inActive][i][j][k].A_3 = fields[active][i][j][k].A_3 + dtn * df.A_3;

              momenta[inActive][i][j][k].ph_0 = momenta[active][i][j][k].ph_0 + dtn * dm.ph_0;
              momenta[inActive][i][j][k].ph_1 = momenta[active][i][j][k].ph_1 + dtn * dm.ph_1;
              momenta[inActive][i][j][k].ph_2 = momenta[active][i][j][k].ph_2 + dtn * dm.ph_2;
              momenta[inActive][i][j][k].ph_3 = momenta[active][i][j][k].ph_3 + dtn * dm.ph_3;

              momenta[inActive][i][j][k].A_1 = momenta[active][i][j][k].A_1 + dtn * dm.A_1;
              momenta[inActive][i][j][k].A_2 = momenta[active][i][j][k].A_2 + dtn * dm.A_2;
              momenta[inActive][i][j][k].A_3 = momenta[active][i][j][k].A_3 + dtn * dm.A_3;
           }
#ifdef _OPENMP
#pragma omp barrier
}
#endif

   active = inActive;
}

#endif

double TRunLattice::compute_rho_static_core(TFields *fields, double twoaR2, double twoaR4,
                                            double h_00, double h_01, double h_02, double h_11, double h_12, double h_22, 
                                            double dx_ph_0, double dx_ph_1, double dx_ph_2, double dx_ph_3, 
                                            double dx_A_2, double dx_A_3, 
                                            double dy_ph_0, double dy_ph_1, double dy_ph_2, double dy_ph_3, 
                                            double dy_A_1, double dy_A_3, 
                                            double dz_ph_0, double dz_ph_1, double dz_ph_2, double dz_ph_3, 
                                            double dz_A_1, double dz_A_2,
                                            bool includeLambdaTerm)
{
   double AcrossGradPh_am[4][3];

   AcrossGradPh_am[0][0] = fields->A_2*dz_ph_0 - fields->A_3*dy_ph_0;
   AcrossGradPh_am[0][1] = fields->A_3*dx_ph_0 - fields->A_1*dz_ph_0;
   AcrossGradPh_am[0][2] = fields->A_1*dy_ph_0 - fields->A_2*dx_ph_0;
   
   AcrossGradPh_am[1][0] = fields->A_2*dz_ph_1 - fields->A_3*dy_ph_1;
   AcrossGradPh_am[1][1] = fields->A_3*dx_ph_1 - fields->A_1*dz_ph_1;
   AcrossGradPh_am[1][2] = fields->A_1*dy_ph_1 - fields->A_2*dx_ph_1;
   
   AcrossGradPh_am[2][0] = fields->A_2*dz_ph_2 - fields->A_3*dy_ph_2;
   AcrossGradPh_am[2][1] = fields->A_3*dx_ph_2 - fields->A_1*dz_ph_2;
   AcrossGradPh_am[2][2] = fields->A_1*dy_ph_2 - fields->A_2*dx_ph_2;
   
   AcrossGradPh_am[3][0] = fields->A_2*dz_ph_3 - fields->A_3*dy_ph_3;
   AcrossGradPh_am[3][1] = fields->A_3*dx_ph_3 - fields->A_1*dz_ph_3;
   AcrossGradPh_am[3][2] = fields->A_1*dy_ph_3 - fields->A_2*dx_ph_3;

   double rho_BR2 = ( pow(dy_A_3 - dz_A_2, 2) + pow(dz_A_1 - dx_A_3, 2) + pow(dx_A_2 - dy_A_1, 2) ) / twoaR4;

   double rho_rH =  (  h_00*(  AcrossGradPh_am[0][0]*AcrossGradPh_am[0][0]
                             + AcrossGradPh_am[0][1]*AcrossGradPh_am[0][1]
                             + AcrossGradPh_am[0][2]*AcrossGradPh_am[0][2])
                     + h_11*(  AcrossGradPh_am[1][0]*AcrossGradPh_am[1][0]
                             + AcrossGradPh_am[1][1]*AcrossGradPh_am[1][1]
                             + AcrossGradPh_am[1][2]*AcrossGradPh_am[1][2])
                     + h_22*(  AcrossGradPh_am[2][0]*AcrossGradPh_am[2][0]
                             + AcrossGradPh_am[2][1]*AcrossGradPh_am[2][1]
                             + AcrossGradPh_am[2][2]*AcrossGradPh_am[2][2])
                     + h_33*(  AcrossGradPh_am[3][0]*AcrossGradPh_am[3][0]
                             + AcrossGradPh_am[3][1]*AcrossGradPh_am[3][1]
                             + AcrossGradPh_am[3][2]*AcrossGradPh_am[3][2])
                     + 2*(  h_01*(  AcrossGradPh_am[0][0]*AcrossGradPh_am[1][0]
                                  + AcrossGradPh_am[0][1]*AcrossGradPh_am[1][1]
                                  + AcrossGradPh_am[0][2]*AcrossGradPh_am[1][2])
                          + h_02*(  AcrossGradPh_am[0][0]*AcrossGradPh_am[2][0]
                                  + AcrossGradPh_am[0][1]*AcrossGradPh_am[2][1]
                                  + AcrossGradPh_am[0][2]*AcrossGradPh_am[2][2])
                          + h_03*(  AcrossGradPh_am[0][0]*AcrossGradPh_am[3][0]
                                  + AcrossGradPh_am[0][1]*AcrossGradPh_am[3][1]
                                  + AcrossGradPh_am[0][2]*AcrossGradPh_am[3][2])
                          + h_12*(  AcrossGradPh_am[1][0]*AcrossGradPh_am[2][0]
                                  + AcrossGradPh_am[1][1]*AcrossGradPh_am[2][1]
                                  + AcrossGradPh_am[1][2]*AcrossGradPh_am[2][2])
                          + h_13*(  AcrossGradPh_am[1][0]*AcrossGradPh_am[3][0]
                                  + AcrossGradPh_am[1][1]*AcrossGradPh_am[3][1]
                                  + AcrossGradPh_am[1][2]*AcrossGradPh_am[3][2])
                          + h_23*(  AcrossGradPh_am[2][0]*AcrossGradPh_am[3][0]
                                  + AcrossGradPh_am[2][1]*AcrossGradPh_am[3][1]
                                  + AcrossGradPh_am[2][2]*AcrossGradPh_am[3][2])
                         )
                    );

   if (0 > rho_rH) rho_rH = 0; else rho_rH /= twoaR4;

   double rho_rG = (  dx_ph_0*dx_ph_0 + dy_ph_0*dy_ph_0 + dz_ph_0*dz_ph_0
                    + dx_ph_1*dx_ph_1 + dy_ph_1*dy_ph_1 + dz_ph_1*dz_ph_1
                    + dx_ph_2*dx_ph_2 + dy_ph_2*dy_ph_2 + dz_ph_2*dz_ph_2
                    + dx_ph_3*dx_ph_3 + dy_ph_3*dy_ph_3 + dz_ph_3*dz_ph_3
                   ) / twoaR2;

   if (includeLambdaTerm)
   {
      double rho_rL = pow(lambda * (  pow(fields->ph_0, 2) + pow(fields->ph_1, 2) 
                                    + pow(fields->ph_2, 2) + pow(fields->ph_3, 2) - 1 ), 2) / 2;

      return(rho_BR2 + rho_rH + rho_rG + rho_rL);
   }

   return(rho_BR2 + rho_rH + rho_rG);
}

double TRunLattice::compute_rho(int i, int j, int k, double aR2, double aR4,
                                TFields fields[MAX_NSIDE][MAX_NSIDE][MAX_NSIDE], 
                                TFields momenta[MAX_NSIDE][MAX_NSIDE][MAX_NSIDE],
                                bool includeLambdaTerm)
/* IMPORTANT: Expects derivatives to be ready. */
{
   double dx_ph_0, dx_ph_1, dx_ph_2, dx_ph_3, dx_A_2, dx_A_3, 
          dy_ph_0, dy_ph_1, dy_ph_2, dy_ph_3, dy_A_1, dy_A_3, 
          dz_ph_0, dz_ph_1, dz_ph_2, dz_ph_3, dz_A_1, dz_A_2;
           
   aux[i][j][k].get_derivatives(dx_ph_0, dx_ph_1, dx_ph_2, dx_ph_3, dx_A_2, dx_A_3, 
                                dy_ph_0, dy_ph_1, dy_ph_2, dy_ph_3, dy_A_1, dy_A_3, 
                                dz_ph_0, dz_ph_1, dz_ph_2, dz_ph_3, dz_A_1, dz_A_2);

   double rho_E = ( sqr(momenta[i][j][k].A_1) + sqr(momenta[i][j][k].A_2) + sqr(momenta[i][j][k].A_3) ) / (2 * aR2);

   double AmAmOaR2 = ( sqr(fields[i][j][k].A_1) + sqr(fields[i][j][k].A_2) + sqr(fields[i][j][k].A_3) ) / aR2;

   double h_00, h_01, h_02, h_11, h_12, h_22,
          invM_00, invM_01, invM_02, invM_03, invM_11, invM_12, invM_13, invM_22, invM_23, invM_33;

   compute_invMh(fields[i][j][k].ph_0, fields[i][j][k].ph_1, fields[i][j][k].ph_2, fields[i][j][k].ph_3, AmAmOaR2,
                 invM_00, invM_01, invM_02, invM_03, invM_11, invM_12, invM_13, invM_22, invM_23, invM_33,
                 h_00, h_01, h_02, h_11, h_12, h_22);

   double rho_static = compute_rho_static_core(&(fields[i][j][k]), 2 * aR2, 2 * aR4,
                                               h_00, h_01, h_02, h_11, h_12, h_22,
                                               dx_ph_0, dx_ph_1, dx_ph_2, dx_ph_3, dx_A_2, dx_A_3, 
                                               dy_ph_0, dy_ph_1, dy_ph_2, dy_ph_3, dy_A_1, dy_A_3, 
                                               dz_ph_0, dz_ph_1, dz_ph_2, dz_ph_3, dz_A_1, dz_A_2,
                                               includeLambdaTerm);

   double rho_M = (  momenta[i][j][k].ph_0 * (momenta[i][j][k].ph_3*invM_03 + momenta[i][j][k].ph_2*invM_02 + momenta[i][j][k].ph_1*invM_01 + momenta[i][j][k].ph_0*invM_00)
                   + momenta[i][j][k].ph_1 * (momenta[i][j][k].ph_3*invM_13 + momenta[i][j][k].ph_2*invM_12 + momenta[i][j][k].ph_1*invM_11 + momenta[i][j][k].ph_0*invM_01)
                   + momenta[i][j][k].ph_2 * (momenta[i][j][k].ph_3*invM_23 + momenta[i][j][k].ph_2*invM_22 + momenta[i][j][k].ph_1*invM_12 + momenta[i][j][k].ph_0*invM_02)
                   + momenta[i][j][k].ph_3 * (momenta[i][j][k].ph_3*invM_33 + momenta[i][j][k].ph_2*invM_23 + momenta[i][j][k].ph_1*invM_13 + momenta[i][j][k].ph_0*invM_03) 
                  ) / (2 * aR4);

   return(rho_static + rho_E + rho_M);
}

double TRunLattice::compute_rho(TFields fields[MAX_NSIDE][MAX_NSIDE][MAX_NSIDE], 
                                TFields momenta[MAX_NSIDE][MAX_NSIDE][MAX_NSIDE],
                                bool includeLambdaTerm)
/* IMPORTANT: Expects derivatives to be ready. */
{
   double rho = 0;
   
   double aR2 = a * a;
   double aR4 = aR2 * aR2;
   
#ifdef _OPENMP
#pragma omp parallel reduction(+ : rho)
{
#pragma omp for
#endif
   for (int i = 0; i <= nSideM1; i++)
       for (int j = 0; j <= nSideM1; j++)
           for (int k = 0; k <= nSideM1; k++)
           {
              rho += compute_rho(i, j, k, aR2, aR4, fields, momenta, includeLambdaTerm);
           }
#ifdef _OPENMP
#pragma omp barrier
}
#endif

   return(rho / ((nSideM1+1)*(nSideM1+1)*(nSideM1+1)));
}

void TRunLattice::init_a(double h0, double o_sim, double o_d, double o_v, double rho_tol)
{
   // Compute initial scale factor from optimized Minkowski rho, omega_sim and H_0.

   a = 1;

   rho_s = RHO_SIM_TO_MKS * compute_rho(fields[active], momenta[active], false);
    
   a = sqrt(G_FACTOR * rho_s / o_sim) * MPC_TO_KM / h0;
   a_0 = a;
    
   // Compute energy densities
    
   rho_s /= a*a;
    
   rho = rho_s / o_sim;

   rho_r = rho * (1 - (o_sim + o_d + o_v)); rho_r_0 = rho_r; // Assuming flatness, Omega_radiation makes up the difference
   rho_d = rho * o_d; rho_d_0 = rho_d;
   rho_v = rho * o_v;
    
   // Rescale A and pph so rho optimization remains valid after spatial scaling
    
#ifdef _OPENMP
#pragma omp parallel
{
#pragma omp for
#endif
   for (int i = 0; i <= nSideM1; i++)
       for (int j = 0; j <= nSideM1; j++)
           for (int k = 0; k <= nSideM1; k++)
           {
              fields[active][i][j][k].A_1 *= a;
              fields[active][i][j][k].A_2 *= a;
              fields[active][i][j][k].A_3 *= a;

              momenta[active][i][j][k].ph_0 *= a;
              momenta[active][i][j][k].ph_1 *= a;
              momenta[active][i][j][k].ph_2 *= a;
              momenta[active][i][j][k].ph_3 *= a;
           }
#ifdef _OPENMP
#pragma omp barrier
}
#endif
}

void TRunLattice::evolve_a(double dtn)
/* IMPORTANT: Expects derivatives to be ready. */
{
   rho_s = RHO_SIM_TO_MKS * compute_rho(fields[active], momenta[active], false);
   rho_r = rho_r_0 * pow(a_0 / a, 4);
   rho_d = rho_d_0 * pow(a_0 / a, 3);
   
   rho = rho_s + rho_r + rho_d + rho_v;
   
   // Forward Euler
   
   a += a*a*sqrt(G_FACTOR * rho) * dtn;
}

int main(int argc, char* argv[])
{
   cout.setf(ios::left, ios::adjustfield);
   cout.precision(OUTPUT_PRECISION);

   int nMaxIterations         = DEFAULT_MAX_ITERATIONS;
   int nSaveSteps             = DEFAULT_SAVE_STEPS;
   int nReportSteps           = DEFAULT_REPORT_STEPS;
   int maxSteps               = -1;       // If > -1, terminate after maxSteps time steps
   double maxEta              = -1.0;     // If > 0, terminate when abs(eta)  >= maxEta
   double daoa0               = DEFAULT_DAOA_0;
   double epsilon             = DEFAULT_EPSILON;
   double lambda_minkowski    = DEFAULT_LAMBDA;
   double h0                  = DEFAULT_H_0;
   double o_sim               = 1;
   double o_d                 = 0;
   double o_v                 = 0;
   double eta                 = 0;
   double detaFactor          = 1;

   string iFile("cartesian");             // Default initial VTK filename (.VTK extension implied)

   string eFile("evolved");               // Default evolved VTK filename (.VTK extension implied)
   string eTitle("Evolved");              // Default evolved VTK title

   // Parse command line

   long int i = 1;
   while (i < argc)
   {
      if (!strcmp("-if", argv[i])) iFile.assign(argv[i + 1]); else
      if (!strcmp("-of", argv[i])) eFile.assign(argv[i + 1]); else
      if (!strcmp("-ot", argv[i])) eTitle.assign(argv[i + 1]); else

      if (!strcmp("-da", argv[i])) daoa0 = atof(argv[i + 1]); else
      if (!strcmp("-la", argv[i])) lambda_minkowski = atof(argv[i + 1]); else
#ifndef EULER
#ifndef RK4
      if (!strcmp("-eo", argv[i])) epsilon = fabs(atof(argv[i + 1])); else
      if (!strcmp("-mi", argv[i])) nMaxIterations = atol(argv[i + 1]); else
#endif
#endif       

      if (!strcmp("-mt", argv[i])) maxEta = atof(argv[i + 1]); else
      if (!strcmp("-ms", argv[i])) maxSteps = atol(argv[i + 1]); else
      if (!strcmp("-sp", argv[i])) nSaveSteps = atol(argv[i + 1]); else
      if (!strcmp("-rp", argv[i])) nReportSteps = atol(argv[i + 1]); else
       
      if (!strcmp("-t0", argv[i])) eta = atof(argv[i + 1]); else
      if (!strcmp("-df", argv[i])) detaFactor = atof(argv[i + 1]); else
      if (!strcmp("-h0", argv[i])) h0 = atof(argv[i + 1]); else
       
      if (!strcmp("-os", argv[i])) o_sim = atof(argv[i + 1]); else
      if (!strcmp("-od", argv[i])) o_d = atof(argv[i + 1]); else
      if (!strcmp("-ov", argv[i])) o_v = atof(argv[i + 1]); else
      {
         cout << "Unrecognized command line argument: " << argv[i] << endl;
         cout << endl;
         cout << "Supported options [defaults in square brackets]:" << endl;
         cout << endl;
         cout << "Filename for initial data input:        -if [\"" << iFile << "\"]" << endl;
         cout << "Filename for evolved data output:       -of [\"" << eFile << "\"]" << endl;
         cout << "Title of evolved data file:             -ot [\"" << eTitle << "\"]" << endl;
         cout << endl;
         cout << "Initial da/a:                           -da [" << daoa0 << "]" << endl;
         cout << "Lambda (potential constant):            -la [" << lambda_minkowski << "]" << endl;
#ifndef EULER
#ifndef RK4
         cout << "Epsilon (negligible error):             -eo [" << epsilon << "]" << endl;
         cout << "Max iterations in solver:               -mi [" << nMaxIterations << "]" << endl;
#endif
#endif
         cout << endl;

         cout << "Max simulated time (conformal):         -mt [" << maxEta << "]" << endl;
         cout << "Max number of time steps:               -ms [" << maxSteps << "]" << endl;
         cout << "Save period (steps):                    -sp [" << nSaveSteps << "]" << endl;
         cout << "Report period (steps):                  -rp [" << nReportSteps << "]" << endl;
         cout << endl;

         cout << "Start time (conformal):                 -t0 [" << eta << "]" << endl;
         cout << "Time step factor for a(t):              -df [" << detaFactor << "]" << endl;
         cout << "Initial Hubble parameter (km/s/Mpc):    -h0 [" << h0 << "]" << endl;

         cout << "Omega_simulation:                       -os [" << o_sim << "]" << endl;
         cout << "Omega_dust:                             -od [" << o_d << "]" << endl;
         cout << "Omega_vacuum:                           -ov [" << o_v << "]" << endl;

         exit(1);
      }

      i += 2;
   }

   cout << "daoa: " << daoa0 << endl;
   cout << "df: " << detaFactor << endl;
#ifndef EULER
#ifndef RK4
   cout << "epsilon: " << epsilon << endl;
   cout << "mi: " << nMaxIterations << endl;
#endif
#endif
   cout << "lambda: " << lambda_minkowski << endl;
   cout << "Integrator: fixed step ";
#ifdef EULER
   cout << "euler";
#else
#ifdef RK4
   cout << "RK4";
#else
   cout << "SV";
#endif
#endif

   cout << endl << flush;

   // Initialize lattice

   TRunLattice lattice;
   lattice.set_lambda(lambda_minkowski);
   lattice.set_a(-1);
    
   // Load initial data

   cout << "Loading initial data." << endl << flush;
    
   double t = 0;
    
   string iTitle;
   lattice.read_vtk(iFile + ".vtk", iTitle, t);
    
   cout << "Done loading." << endl << flush;
   cout << "nSide: " << lattice.get_nSide() << endl;

   lattice.compute_derivatives(fields[lattice.get_active()]); // Initialize derivatives

   if (0 > lattice.get_a())
   {
      t = 0;

      cout << "Using h0 = " << h0 << endl;

      lattice.init_a(h0, o_sim, o_d, o_v, epsilon);

      cout << "Setting a0 = " << lattice.get_a() << endl;
   }

   cout << "rho_c = " << sqr(h0 / MPC_TO_KM) / G_FACTOR << endl;
    
   lattice.set_lambda(lambda_minkowski / lattice.get_a());

   // Prepare for main loop

   char cBuff[128]; // Character buffer for integer to string conversion
   string sBuff("");

   i = 0;
    
   double deta = daoa0 / (lattice.get_a() * sqrt(G_FACTOR * lattice.get_rho()));
    
   // Save initial data

   cout << "Step " << i << "; eta = " << eta << "; a = " << lattice.get_a() << "; t = " << t
        << "; rho_s = " << lattice.get_rho_s() << "; rho_r = " << lattice.get_rho_r()
        << "; rho_d = " << lattice.get_rho_d() << "; rho = " << lattice.get_rho();
   cout << endl << flush;
    
   sprintf(cBuff, "%lu", i);
   sBuff = eFile + "_" + cBuff;

   cout << "Saving initial data to \"" << sBuff << ".vtk\"." << endl << flush;

   lattice.write_vtk(sBuff + ".vtk", eTitle, t);

   cout << "Done saving." << endl << flush;
    
   // Main loop
    
   while ((0 > maxSteps) || (i < maxSteps))
   {
      int err = 0;

      lattice.evolve_a(detaFactor * deta / 2);                              // Half step, so we use mid-point a
      lattice.set_lambda(lambda_minkowski / lattice.get_a());
      t += deta * lattice.get_a() / 2;                                      // Time at half step
          
#ifdef EULER
      lattice.euler(deta * SEC_TO_SIM);
#else
#ifdef RK4
      lattice.rk4(deta * SEC_TO_SIM);
#else
      err = lattice.sv(deta * SEC_TO_SIM, epsilon, nMaxIterations);
#endif
#endif

      lattice.compute_derivatives(fields[lattice.get_active()]);            // Update derivatives ahead of next iteration
      lattice.evolve_a(detaFactor * deta / 2);                              // Second half of step in a evolution
      lattice.set_lambda(lambda_minkowski / lattice.get_a());
      t += deta * lattice.get_a() / 2;                                      // Time at end of step

      i++;

      eta += deta;

      bool doSave = ((i == maxSteps) || (lattice.get_a() < 0) || isnan(lattice.get_a()) ||
                     (lattice.get_rho_s() < 0) || isnan(lattice.get_rho_s()) ||
                     ((i / nSaveSteps) != ((i - 1) / nSaveSteps)) || 
                     ((maxEta > 0.0) && (abs(eta) >= maxEta)) || (0.0 == deta));
      bool doReport = ((err) || (doSave) || (0 == nReportSteps) || 
                       ((i / nReportSteps) != ((i - 1) / nReportSteps)));

      if (doReport)
      {
         cout << "Step " << i << "; eta = " << eta << "; a = " << lattice.get_a() << "; t = " << t
              << "; rho_s = " << lattice.get_rho_s() << "; rho_r = " << lattice.get_rho_r()
              << "; rho_d = " << lattice.get_rho_d() << "; rho = " << lattice.get_rho();

         if (err & 1) cout << " Max number of iterations exceeded in field step.";
         if (err & 2) cout << " Max number of iterations exceeded in momentum step.";

         cout << endl << flush;
      }

      if (doSave)
      {
         sprintf(cBuff, "%lu", i);
         sBuff = eFile + "_" + cBuff;

         cout << "Saving data to \"" << sBuff << ".vtk\"." << endl << flush;

         lattice.write_vtk(sBuff + ".vtk", eTitle, t);

         cout << "Done saving." << endl << flush;
      }

      if ((i == maxSteps) || (lattice.get_rho_s() < 0) || isnan(lattice.get_rho_s()) || 
                             (lattice.get_a() < 0) || isnan(lattice.get_a()) || 
                             ((maxEta > 0) && (abs(eta) >= maxEta))) break;
   }

   return(0);
}
