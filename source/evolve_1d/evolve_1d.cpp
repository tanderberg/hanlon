/* evolve_1d.cpp : Euler, RK4 and symplectic (default) integrators (static 1D).
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
#include "core_1d.cpp"

const double DEFAULT_DTN                    = 0.0001;    /* Default time step, in lattice units */
const double DEFAULT_EPSILON                = 1e-6;      /* Default "negligible" error */
const double DEFAULT_REPORT_PERIOD          = 0.0;       /* Approximate simulated time between successive reports */
const double DEFAULT_SAVE_PERIOD            = 1.0;       /* Approximate simulated time between successive saves */
const int DEFAULT_MAX_ITERATIONS            = 100;       /* Default max number of iterations in SV solver */

#ifdef RK4
TFields dotFields[3][MAX_NSIDE];                         /* Max allowed number of cells always allocated. Ugly but fast. */
TFields dotMomenta[3][MAX_NSIDE];                        /* Ditto */
#else
TFields dotMomenta[MAX_NSIDE];                           /* Ditto */
TFields dotFields[MAX_NSIDE];                            /* Ditto */
#endif

class TRunLattice : public TLattice
{
public:

#ifdef RK4
   void rk4(double dtn);
#else
   int sv(double dtn, double epsilon, int nMaxIterations);
   void euler(double dtn);
#endif
};

#ifdef RK4

void TRunLattice::rk4(double dtn)
{
   int substep_1 = (active + 1) % 3;
   int substep_2 = (substep_1 + 1) % 3;

   double dtHalf = dtn / 2;
   double dtThird = dtn / 3;
   double dtSixth = dtThird / 2;

   // First Euler substep
   
   compute_derivatives(fields[active]);

#ifdef _OPENMP
#pragma omp parallel
{
#pragma omp for
#endif
   for (int i = 0; i <= nSideM1; i++)
   {
      compute_dots(i, fields[active], &(momenta[active][i]), dotFields[active][i], dotMomenta[active][i]);
           
      fields[substep_1][i].ph_0 = fields[active][i].ph_0 + dtHalf * dotFields[active][i].ph_0;
      fields[substep_1][i].ph_1 = fields[active][i].ph_1 + dtHalf * dotFields[active][i].ph_1;
      fields[substep_1][i].ph_2 = fields[active][i].ph_2 + dtHalf * dotFields[active][i].ph_2;
      fields[substep_1][i].ph_3 = fields[active][i].ph_3 + dtHalf * dotFields[active][i].ph_3;

      fields[substep_1][i].A_1 = fields[active][i].A_1 + dtHalf * dotFields[active][i].A_1;
      fields[substep_1][i].A_2 = fields[active][i].A_2 + dtHalf * dotFields[active][i].A_2;
      fields[substep_1][i].A_3 = fields[active][i].A_3 + dtHalf * dotFields[active][i].A_3;

      momenta[substep_1][i].ph_0 = momenta[active][i].ph_0 + dtHalf * dotMomenta[active][i].ph_0;
      momenta[substep_1][i].ph_1 = momenta[active][i].ph_1 + dtHalf * dotMomenta[active][i].ph_1;
      momenta[substep_1][i].ph_2 = momenta[active][i].ph_2 + dtHalf * dotMomenta[active][i].ph_2;
      momenta[substep_1][i].ph_3 = momenta[active][i].ph_3 + dtHalf * dotMomenta[active][i].ph_3;

      momenta[substep_1][i].A_1 = momenta[active][i].A_1 + dtHalf * dotMomenta[active][i].A_1;
      momenta[substep_1][i].A_2 = momenta[active][i].A_2 + dtHalf * dotMomenta[active][i].A_2;
      momenta[substep_1][i].A_3 = momenta[active][i].A_3 + dtHalf * dotMomenta[active][i].A_3;
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
   {
      compute_dots(i, fields[substep_1], &(momenta[substep_1][i]), dotFields[substep_1][i], dotMomenta[substep_1][i]);
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
   {
      fields[substep_1][i].ph_0 = fields[active][i].ph_0 + dtHalf * dotFields[substep_1][i].ph_0;
      fields[substep_1][i].ph_1 = fields[active][i].ph_1 + dtHalf * dotFields[substep_1][i].ph_1;
      fields[substep_1][i].ph_2 = fields[active][i].ph_2 + dtHalf * dotFields[substep_1][i].ph_2;
      fields[substep_1][i].ph_3 = fields[active][i].ph_3 + dtHalf * dotFields[substep_1][i].ph_3;

      fields[substep_1][i].A_1 = fields[active][i].A_1 + dtHalf * dotFields[substep_1][i].A_1;
      fields[substep_1][i].A_2 = fields[active][i].A_2 + dtHalf * dotFields[substep_1][i].A_2;
      fields[substep_1][i].A_3 = fields[active][i].A_3 + dtHalf * dotFields[substep_1][i].A_3;

      momenta[substep_1][i].ph_0 = momenta[active][i].ph_0 + dtHalf * dotMomenta[substep_1][i].ph_0;
      momenta[substep_1][i].ph_1 = momenta[active][i].ph_1 + dtHalf * dotMomenta[substep_1][i].ph_1;
      momenta[substep_1][i].ph_2 = momenta[active][i].ph_2 + dtHalf * dotMomenta[substep_1][i].ph_2;
      momenta[substep_1][i].ph_3 = momenta[active][i].ph_3 + dtHalf * dotMomenta[substep_1][i].ph_3;

      momenta[substep_1][i].A_1 = momenta[active][i].A_1 + dtHalf * dotMomenta[substep_1][i].A_1;
      momenta[substep_1][i].A_2 = momenta[active][i].A_2 + dtHalf * dotMomenta[substep_1][i].A_2;
      momenta[substep_1][i].A_3 = momenta[active][i].A_3 + dtHalf * dotMomenta[substep_1][i].A_3;
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
   {
      compute_dots(i, fields[substep_1], &(momenta[substep_1][i]), dotFields[substep_2][i], dotMomenta[substep_2][i]);
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
   {
      fields[substep_1][i].ph_0 = fields[active][i].ph_0 + dtHalf * dotFields[substep_2][i].ph_0;
      fields[substep_1][i].ph_1 = fields[active][i].ph_1 + dtHalf * dotFields[substep_2][i].ph_1;
      fields[substep_1][i].ph_2 = fields[active][i].ph_2 + dtHalf * dotFields[substep_2][i].ph_2;
      fields[substep_1][i].ph_3 = fields[active][i].ph_3 + dtHalf * dotFields[substep_2][i].ph_3;

      fields[substep_1][i].A_1 = fields[active][i].A_1 + dtHalf * dotFields[substep_2][i].A_1;
      fields[substep_1][i].A_2 = fields[active][i].A_2 + dtHalf * dotFields[substep_2][i].A_2;
      fields[substep_1][i].A_3 = fields[active][i].A_3 + dtHalf * dotFields[substep_2][i].A_3;

      momenta[substep_1][i].ph_0 = momenta[active][i].ph_0 + dtHalf * dotMomenta[substep_2][i].ph_0;
      momenta[substep_1][i].ph_1 = momenta[active][i].ph_1 + dtHalf * dotMomenta[substep_2][i].ph_1;
      momenta[substep_1][i].ph_2 = momenta[active][i].ph_2 + dtHalf * dotMomenta[substep_2][i].ph_2;
      momenta[substep_1][i].ph_3 = momenta[active][i].ph_3 + dtHalf * dotMomenta[substep_2][i].ph_3;

      momenta[substep_1][i].A_1 = momenta[active][i].A_1 + dtHalf * dotMomenta[substep_2][i].A_1;
      momenta[substep_1][i].A_2 = momenta[active][i].A_2 + dtHalf * dotMomenta[substep_2][i].A_2;
      momenta[substep_1][i].A_3 = momenta[active][i].A_3 + dtHalf * dotMomenta[substep_2][i].A_3;
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
   {
      TFields df, dm;
      compute_dots(i, j, k, fields[substep_1], &(momenta[substep_1][i]), df, dm);
           
      fields[active][i].ph_0 += dtSixth * (dotFields[active][i].ph_0 + df.ph_0) +
                                           dtThird * (dotFields[substep_1][i].ph_0 + dotFields[substep_2][i].ph_0);
      fields[active][i].ph_1 += dtSixth * (dotFields[active][i].ph_1 + df.ph_1) +
                                           dtThird * (dotFields[substep_1][i].ph_1 + dotFields[substep_2][i].ph_1);
      fields[active][i].ph_2 += dtSixth * (dotFields[active][i].ph_2 + df.ph_2) +
                                           dtThird * (dotFields[substep_1][i].ph_2 + dotFields[substep_2][i].ph_2);
      fields[active][i].ph_3 += dtSixth * (dotFields[active][i].ph_3 + df.ph_3) +
                                           dtThird * (dotFields[substep_1][i].ph_3 + dotFields[substep_2][i].ph_3);

      fields[active][i].A_1 += dtSixth * (dotFields[active][i].A_1 + df.A_1) +
                                          dtThird * (dotFields[substep_1][i].A_1 + dotFields[substep_2][i].A_1);
      fields[active][i].A_2 += dtSixth * (dotFields[active][i].A_2 + df.A_2) +
                                          dtThird * (dotFields[substep_1][i].A_2 + dotFields[substep_2][i].A_2);
      fields[active][i].A_3 += dtSixth * (dotFields[active][i].A_3 + df.A_3) +
                                          dtThird * (dotFields[substep_1][i].A_3 + dotFields[substep_2][i].A_3);

      momenta[active][i].ph_0 += dtSixth * (dotMomenta[active][i].ph_0 + dm.ph_0) +
                                            dtThird * (dotMomenta[substep_1][i].ph_0 + dotMomenta[substep_2][i].ph_0);
      momenta[active][i].ph_1 += dtSixth * (dotMomenta[active][i].ph_1 + dm.ph_1) +
                                            dtThird * (dotMomenta[substep_1][i].ph_1 + dotMomenta[substep_2][i].ph_1);
      momenta[active][i].ph_2 += dtSixth * (dotMomenta[active][i].ph_2 + dm.ph_2) +
                                            dtThird * (dotMomenta[substep_1][i].ph_2 + dotMomenta[substep_2][i].ph_2);
      momenta[active][i].ph_3 += dtSixth * (dotMomenta[active][i].ph_3 + dm.ph_3) +
                                            dtThird * (dotMomenta[substep_1][i].ph_3 + dotMomenta[substep_2][i].ph_3);

      momenta[active][i].A_1 += dtSixth * (dotMomenta[active][i].A_1 + dm.A_1) +
                                           dtThird * (dotMomenta[substep_1][i].A_1 + dotMomenta[substep_2][i].A_1);
      momenta[active][i].A_2 += dtSixth * (dotMomenta[active][i].A_2 + dm.A_2) +
                                           dtThird * (dotMomenta[substep_1][i].A_2 + dotMomenta[substep_2][i].A_2);
      momenta[active][i].A_3 += dtSixth * (dotMomenta[active][i].A_3 + dm.A_3) +
                                           dtThird * (dotMomenta[substep_1][i].A_3 + dotMomenta[substep_2][i].A_3);
   }
#ifdef _OPENMP
#pragma omp barrier
}
#endif

   // Done!
}

#else

int TRunLattice::sv(double dtn, double epsilon, int nMaxIterations)
{
   int scratch = (active + 1) % 3;
   int inActive = (scratch + 1) % 3;

   int err = 0;

   // Precompute auxiliary data

   double dtnHalf = dtn / 2;
   
   compute_derivatives(fields[active]);

   // Initial scratch fields guess: forward Euler

#ifdef _OPENMP
#pragma omp parallel
{
#pragma omp for
#endif
   for (int i = 0; i <= nSideM1; i++)
   {
      TFields df;
      compute_dotFields(&(fields[active][i]), &(momenta[active][i]), df);

      fields[scratch][i].ph_0 = fields[active][i].ph_0 + dtnHalf * df.ph_0;
      fields[scratch][i].ph_1 = fields[active][i].ph_1 + dtnHalf * df.ph_1;
      fields[scratch][i].ph_2 = fields[active][i].ph_2 + dtnHalf * df.ph_2;
      fields[scratch][i].ph_3 = fields[active][i].ph_3 + dtnHalf * df.ph_3;
              
      fields[scratch][i].A_1 = fields[active][i].A_1 + dtnHalf * df.A_1;
      fields[scratch][i].A_2 = fields[active][i].A_2 + dtnHalf * df.A_2;
      fields[scratch][i].A_3 = fields[active][i].A_3 + dtnHalf * df.A_3;
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
      {
         compute_dotFields_and_dDotFields_dFields(&(fields[scratch][i]), &(momenta[active][i]), df, A);

         A[0][0] *= -dtnHalf; A[0][0] += 1;
         A[0][1] *= -dtnHalf;
         A[0][2] *= -dtnHalf;
         A[0][3] *= -dtnHalf;
         A[0][4] *= -dtnHalf;
         A[0][5] *= -dtnHalf;
         A[0][6] *= -dtnHalf;

         A[1][0] *= -dtnHalf;
         A[1][1] *= -dtnHalf; A[1][1] += 1;
         A[1][2] *= -dtnHalf;
         A[1][3] *= -dtnHalf;
         A[1][4] *= -dtnHalf;
         A[1][5] *= -dtnHalf;
         A[1][6] *= -dtnHalf;

         A[2][0] *= -dtnHalf;
         A[2][1] *= -dtnHalf;
         A[2][2] *= -dtnHalf; A[2][2] += 1;
         A[2][3] *= -dtnHalf;
         A[2][4] *= -dtnHalf;
         A[2][5] *= -dtnHalf;
         A[2][6] *= -dtnHalf;

         A[3][0] *= -dtnHalf;
         A[3][1] *= -dtnHalf;
         A[3][2] *= -dtnHalf;
         A[3][3] *= -dtnHalf; A[3][3] += 1;
         A[3][4] *= -dtnHalf;
         A[3][5] *= -dtnHalf;
         A[3][6] *= -dtnHalf;

         A[4][4] = 1;
         A[5][5] = 1;
         A[6][6] = 1;
                   
         b[0] = (fields[scratch][i].ph_0 - fields[active][i].ph_0) - dtnHalf * df.ph_0;
         b[1] = (fields[scratch][i].ph_1 - fields[active][i].ph_1) - dtnHalf * df.ph_1;
         b[2] = (fields[scratch][i].ph_2 - fields[active][i].ph_2) - dtnHalf * df.ph_2;
         b[3] = (fields[scratch][i].ph_3 - fields[active][i].ph_3) - dtnHalf * df.ph_3;
                   
         b[4] = (fields[scratch][i].A_1 - fields[active][i].A_1) - dtnHalf * df.A_1;
         b[5] = (fields[scratch][i].A_2 - fields[active][i].A_2) - dtnHalf * df.A_2;
         b[6] = (fields[scratch][i].A_3 - fields[active][i].A_3) - dtnHalf * df.A_3;

         solve7special(&(A[0][0]), &(indx[0]), &(b[0]));

         f.ph_0 = fields[scratch][i].ph_0 - b[0];
         f.ph_1 = fields[scratch][i].ph_1 - b[1];
         f.ph_2 = fields[scratch][i].ph_2 - b[2];
         f.ph_3 = fields[scratch][i].ph_3 - b[3];

         f.A_1 = fields[scratch][i].A_1 - b[4];
         f.A_2 = fields[scratch][i].A_2 - b[5];
         f.A_3 = fields[scratch][i].A_3 - b[6];

         deviation = fabs(b[0]) / (fabs(f.ph_0) + 1); if (local_max_rel_dy < deviation) local_max_rel_dy = deviation;
         deviation = fabs(b[1]) / (fabs(f.ph_1) + 1); if (local_max_rel_dy < deviation) local_max_rel_dy = deviation;
         deviation = fabs(b[2]) / (fabs(f.ph_2) + 1); if (local_max_rel_dy < deviation) local_max_rel_dy = deviation;
         deviation = fabs(b[3]) / (fabs(f.ph_3) + 1); if (local_max_rel_dy < deviation) local_max_rel_dy = deviation;

         deviation = fabs(b[4]) / (fabs(f.A_1) + 1); if (local_max_rel_dy < deviation) local_max_rel_dy = deviation;
         deviation = fabs(b[5]) / (fabs(f.A_2) + 1); if (local_max_rel_dy < deviation) local_max_rel_dy = deviation;
         deviation = fabs(b[6]) / (fabs(f.A_3) + 1); if (local_max_rel_dy < deviation) local_max_rel_dy = deviation;

         fields[inActive][i] = f;
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
   {
      compute_dotMomenta(i, fields[scratch], &(momenta[active][i]), dotMomenta[i]);

      momenta[scratch][i].ph_0 = momenta[active][i].ph_0 + dtn * dotMomenta[i].ph_0;
      momenta[scratch][i].ph_1 = momenta[active][i].ph_1 + dtn * dotMomenta[i].ph_1;
      momenta[scratch][i].ph_2 = momenta[active][i].ph_2 + dtn * dotMomenta[i].ph_2;
      momenta[scratch][i].ph_3 = momenta[active][i].ph_3 + dtn * dotMomenta[i].ph_3;
              
      momenta[scratch][i].A_1 = momenta[active][i].A_1 + dtn * dotMomenta[i].A_1;
      momenta[scratch][i].A_2 = momenta[active][i].A_2 + dtn * dotMomenta[i].A_2;
      momenta[scratch][i].A_3 = momenta[active][i].A_3 + dtn * dotMomenta[i].A_3;
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
       {
          compute_dotMomenta_and_dDotFields_dFields(i, fields[scratch], &(momenta[from_p][i]), dm, A);

          A[0][0] *= dtnHalf; A[0][0] += 1;
          A[0][1] *= dtnHalf;
          A[0][2] *= dtnHalf;
          A[0][3] *= dtnHalf;
          A[0][4] *= dtnHalf;
          A[0][5] *= dtnHalf;
          A[0][6] *= dtnHalf;

          A[1][0] *= dtnHalf;
          A[1][1] *= dtnHalf; A[1][1] += 1;
          A[1][2] *= dtnHalf;
          A[1][3] *= dtnHalf;
          A[1][4] *= dtnHalf;
          A[1][5] *= dtnHalf;
          A[1][6] *= dtnHalf;

          A[2][0] *= dtnHalf;
          A[2][1] *= dtnHalf;
          A[2][2] *= dtnHalf; A[2][2] += 1;
          A[2][3] *= dtnHalf;
          A[2][4] *= dtnHalf;
          A[2][5] *= dtnHalf;
          A[2][6] *= dtnHalf;

          A[3][0] *= dtnHalf;
          A[3][1] *= dtnHalf;
          A[3][2] *= dtnHalf;
          A[3][3] *= dtnHalf; A[3][3] += 1;
          A[3][4] *= dtnHalf;
          A[3][5] *= dtnHalf;
          A[3][6] *= dtnHalf;

          A[4][4] = 1;
          A[5][5] = 1;
          A[6][6] = 1;
 
          b[0] = (momenta[from_p][i].ph_0 - momenta[active][i].ph_0) - dtnHalf * (dotMomenta[i].ph_0 + dm.ph_0);
          b[1] = (momenta[from_p][i].ph_1 - momenta[active][i].ph_1) - dtnHalf * (dotMomenta[i].ph_1 + dm.ph_1);
          b[2] = (momenta[from_p][i].ph_2 - momenta[active][i].ph_2) - dtnHalf * (dotMomenta[i].ph_2 + dm.ph_2);
          b[3] = (momenta[from_p][i].ph_3 - momenta[active][i].ph_3) - dtnHalf * (dotMomenta[i].ph_3 + dm.ph_3);
                   
          b[4] = (momenta[from_p][i].A_1 - momenta[active][i].A_1) - dtnHalf * (dotMomenta[i].A_1 + dm.A_1);
          b[5] = (momenta[from_p][i].A_2 - momenta[active][i].A_2) - dtnHalf * (dotMomenta[i].A_2 + dm.A_2);
          b[6] = (momenta[from_p][i].A_3 - momenta[active][i].A_3) - dtnHalf * (dotMomenta[i].A_3 + dm.A_3);

          solve7special(&(A[0][0]), &(indx[0]), &(b[0]));

          m.ph_0 = momenta[from_p][i].ph_0 - b[0];
          m.ph_1 = momenta[from_p][i].ph_1 - b[1];
          m.ph_2 = momenta[from_p][i].ph_2 - b[2];
          m.ph_3 = momenta[from_p][i].ph_3 - b[3];

          m.A_1 = momenta[from_p][i].A_1 - b[4];
          m.A_2 = momenta[from_p][i].A_2 - b[5];
          m.A_3 = momenta[from_p][i].A_3 - b[6];

          deviation = fabs(b[0]) / (fabs(m.ph_0) + 1); if (local_max_rel_dy < deviation) local_max_rel_dy = deviation;
          deviation = fabs(b[1]) / (fabs(m.ph_1) + 1); if (local_max_rel_dy < deviation) local_max_rel_dy = deviation;
          deviation = fabs(b[2]) / (fabs(m.ph_2) + 1); if (local_max_rel_dy < deviation) local_max_rel_dy = deviation;
          deviation = fabs(b[3]) / (fabs(m.ph_3) + 1); if (local_max_rel_dy < deviation) local_max_rel_dy = deviation;

          deviation = fabs(b[4]) / (fabs(m.A_1) + 1); if (local_max_rel_dy < deviation) local_max_rel_dy = deviation;
          deviation = fabs(b[5]) / (fabs(m.A_2) + 1); if (local_max_rel_dy < deviation) local_max_rel_dy = deviation;
          deviation = fabs(b[6]) / (fabs(m.A_3) + 1); if (local_max_rel_dy < deviation) local_max_rel_dy = deviation;

          momenta[to_p][i] = m;
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
      memcpy(&(momenta[scratch][0]), &(momenta[from_p][0]), (nSideM1+1) * sizeof(double) * NCOMPS);
   }

   // Compute final field values in scratch array

#ifdef _OPENMP
#pragma omp parallel
{
#pragma omp for
#endif
   for (int i = 0; i <= nSideM1; i++)
   {
      TFields df;
      compute_dotFields(&(fields[scratch][i]), &(momenta[scratch][i]), df);

      fields[scratch][i].ph_0 += dtnHalf * df.ph_0;
      fields[scratch][i].ph_1 += dtnHalf * df.ph_1;
      fields[scratch][i].ph_2 += dtnHalf * df.ph_2;
      fields[scratch][i].ph_3 += dtnHalf * df.ph_3;

      fields[scratch][i].A_1 += dtnHalf * df.A_1;
      fields[scratch][i].A_2 += dtnHalf * df.A_2;
      fields[scratch][i].A_3 += dtnHalf * df.A_3;
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
{
   int inActive;

   if (0 == active) inActive = 1; else inActive = 0;

   // Forward Euler
   
   compute_derivatives(fields[active]);

#ifdef _OPENMP
#pragma omp parallel
{
#pragma omp for
#endif
   for (int i = 0; i <= nSideM1; i++)
   {
      TFields df, dm;
      compute_dots(i, fields[active], &(momenta[active][i]), df, dm);

      fields[inActive][i].ph_0 = fields[active][i].ph_0 + dtn * df.ph_0;
      fields[inActive][i].ph_1 = fields[active][i].ph_1 + dtn * df.ph_1;
      fields[inActive][i].ph_2 = fields[active][i].ph_2 + dtn * df.ph_2;
      fields[inActive][i].ph_3 = fields[active][i].ph_3 + dtn * df.ph_3;

      fields[inActive][i].A_1 = fields[active][i].A_1 + dtn * df.A_1;
      fields[inActive][i].A_2 = fields[active][i].A_2 + dtn * df.A_2;
      fields[inActive][i].A_3 = fields[active][i].A_3 + dtn * df.A_3;

      momenta[inActive][i].ph_0 = momenta[active][i].ph_0 + dtn * dm.ph_0;
      momenta[inActive][i].ph_1 = momenta[active][i].ph_1 + dtn * dm.ph_1;
      momenta[inActive][i].ph_2 = momenta[active][i].ph_2 + dtn * dm.ph_2;
      momenta[inActive][i].ph_3 = momenta[active][i].ph_3 + dtn * dm.ph_3;

      momenta[inActive][i].A_1 = momenta[active][i].A_1 + dtn * dm.A_1;
      momenta[inActive][i].A_2 = momenta[active][i].A_2 + dtn * dm.A_2;
      momenta[inActive][i].A_3 = momenta[active][i].A_3 + dtn * dm.A_3;
   }
#ifdef _OPENMP
#pragma omp barrier
}
#endif

   active = inActive;
}

#endif

int main(int argc, char* argv[])
{
   cout.setf(ios::left, ios::adjustfield);
   cout.precision(OUTPUT_PRECISION);

   int nMaxIterations         = DEFAULT_MAX_ITERATIONS;
   int maxSteps               = -1;       // If > -1, terminate after maxSteps time steps
   double maxTime             = -1.0;     // If > 0, terminate when time >= maxTime
   double savePeriod          = DEFAULT_SAVE_PERIOD;
   double reportPeriod        = DEFAULT_REPORT_PERIOD;
   double dtn                 = DEFAULT_DTN;
   double epsilon             = DEFAULT_EPSILON;
   double lambda              = DEFAULT_LAMBDA;

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

      if (!strcmp("-dt", argv[i])) dtn = atof(argv[i + 1]); else
      if (!strcmp("-la", argv[i])) lambda = atof(argv[i + 1]); else
#ifndef EULER
#ifndef RK4
      if (!strcmp("-eo", argv[i])) epsilon = fabs(atof(argv[i + 1])); else
      if (!strcmp("-mi", argv[i])) nMaxIterations = atol(argv[i + 1]); else
#endif
#endif       

      if (!strcmp("-mt", argv[i])) maxTime = atof(argv[i + 1]); else
      if (!strcmp("-ms", argv[i])) maxSteps = atol(argv[i + 1]); else
      if (!strcmp("-sp", argv[i])) savePeriod = atof(argv[i + 1]); else
      if (!strcmp("-rp", argv[i])) reportPeriod = atof(argv[i + 1]); else
      {
         cout << "Unrecognized command line argument: " << argv[i] << endl;
         cout << endl;
         cout << "Supported options [defaults in square brackets]:" << endl;
         cout << endl;
         cout << "Filename for initial data input:        -if [\"" << iFile << "\"]" << endl;
         cout << "Filename for evolved data output:       -of [\"" << eFile << "\"]" << endl;
         cout << "Title of evolved data file:             -ot [\"" << eTitle << "\"]" << endl;
         cout << endl;
         cout << "Time step (lattice units):              -dt [" << dtn << "]" << endl;
         cout << "Lambda (potential constant):            -la [" << lambda << "]" << endl;
#ifndef EULER
#ifndef RK4
         cout << "Epsilon (negligible error):             -eo [" << epsilon << "]" << endl;
         cout << "Max iterations in solver:               -mi [" << nMaxIterations << "]" << endl;
#endif
#endif
         cout << endl;
         cout << "Max simulated time (lattice units):     -mt [" << maxTime << "]" << endl;
         cout << "Max number of time steps:               -ms [" << maxSteps << "]" << endl;
         cout << "Save period (lattice units):            -sp [" << savePeriod << "]" << endl;
         cout << "Report period (lattice units):          -rp [" << reportPeriod << "]" << endl;
         cout << endl;

         exit(1);
      }

      i += 2;
   }

   cout << "dt: " << dtn << endl;
#ifndef EULER
#ifndef RK4
   cout << "epsilon: " << epsilon << endl;
   cout << "mi: " << nMaxIterations << endl;
#endif
#endif
   cout << "lambda: " << lambda << endl;
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
   lattice.set_lambda(lambda);
    
   // Load initial data

   cout << "Loading initial data." << endl << flush;

   string iTitle;
   lattice.read_vtk(iFile + ".vtk", iTitle);

   cout << "Done loading." << endl << flush;
   cout << "nSide: " << lattice.get_nSide() << endl;

   // Main loop

   char cBuff[128]; // Character buffer for integer to string conversion
   string sBuff("");

   double t = 0;
   double prev_t;

   i = 0;

   while ((0 > maxSteps) || (i < maxSteps))
   {
      int err = 0;
#ifdef EULER
      lattice.euler(dtn);
#else
#ifdef RK4
      lattice.rk4(dtn);
#else
      err = lattice.sv(dtn, epsilon, nMaxIterations);
#endif
#endif

      i++;

      prev_t = t;
      t += dtn;

      bool doSave = ((i == maxSteps) || 
                     ((long)((t + dtn / 2) / savePeriod) > (long)((prev_t + dtn / 2) / savePeriod)) || 
                     ((maxTime > 0.0) && (t >= maxTime)) || (0.0 == dtn));
      bool doReport = ((err) || (doSave) || (0.0 == reportPeriod) || 
                       ((long)((t + dtn / 2) / reportPeriod) > (long)((prev_t + dtn / 2) / reportPeriod)));

      if (doReport)
      {
         cout << "Step " << i << "; t = " << t;

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

      if ((i == maxSteps) || ((maxTime > 0) && (t >= maxTime))) break;
   }

   return(0);
}
