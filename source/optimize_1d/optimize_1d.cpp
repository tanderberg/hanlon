/* optimize_1d.cpp : Minimizes local energy density (scalars in polar form, 1D).
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
#include "fields_aux_1d.cpp"

#ifndef N_THREADS
const int N_THREADS              = 4;                 /* Must divide lattice side evenly, ideally = number of CPU cores. */
#endif

const double DEFAULT_RHO_REL_TOL = 0.0001;            /* Default relative tolerance in rho */

const double DEFAULT_STEP_SIZE   = 0.1;
const double DEFAULT_TOL         = 1e-3;
const double DEFAULT_GRAD_TOL    = 1e-6;
const size_t DEFAULT_MAX_ITER    = 1000;

TFields fields[2][MAX_NSIDE];                         /* Max allowed number of cells always allocated. Ugly but fast. */
TAux aux[MAX_NSIDE];                                  /* Ditto */

class TOptimizeLattice;

struct TParams
{
   int i;
   TOptimizeLattice *lattice;

   gsl_function f_rho_th1;
   gsl_function f_rho_th2;
   gsl_function f_rho_th3;
   gsl_function f_rho_A1;
   gsl_function f_rho_A2;
   gsl_function f_rho_A3;
};

#include "lattice_1d.cpp"

using namespace std;

class TOptimizeLattice : public TLattice
{
   TParams par;
   gsl_multimin_function_fdf f_rho;
   gsl_vector *v;
   const gsl_multimin_fdfminimizer_type *T;
   gsl_multimin_fdfminimizer *s;
   
public:

   TOptimizeLattice();
   ~TOptimizeLattice();

   void minimize(const double step_size, const double tol, const double grad_tol, const size_t max_iter, 
                 const int from_i, const int to_i);
};

double rho_th1(double th1, void * params)
{
   TParams *par = (TParams*)params;

   int i = par->i;
   
   par->lattice->copy2private(i);
   par->lattice->privateFields[4].th_1 = wrap(th1);

   par->lattice->compute_affected_private_derivatives(0);
   par->lattice->privateAux[4].compute_hG(par->lattice->privateFields[4].th_1, 
                                          par->lattice->privateFields[4].th_2, 
                                          par->lattice->privateFields[4].th_3);

   return par->lattice->affected_private_static_rho();
}

double rho_th2(double th2, void * params)
{
   TParams *par = (TParams*)params;

   int i = par->i;
   
   par->lattice->copy2private(i);
   par->lattice->privateFields[4].th_2 = wrap(th2);

   par->lattice->compute_affected_private_derivatives(1);
   par->lattice->privateAux[4].compute_hG(par->lattice->privateFields[4].th_1, 
                                          par->lattice->privateFields[4].th_2, 
                                          par->lattice->privateFields[4].th_3);

   return par->lattice->affected_private_static_rho();
}

double rho_th3(double th3, void * params)
{
   TParams *par = (TParams*)params;

   int i = par->i;
   
   par->lattice->copy2private(i);
   par->lattice->privateFields[4].th_3 = wrap(th3);

   par->lattice->compute_affected_private_derivatives(2);
   par->lattice->privateAux[4].compute_hG(par->lattice->privateFields[4].th_1, 
                                          par->lattice->privateFields[4].th_2, 
                                          par->lattice->privateFields[4].th_3);

   return par->lattice->affected_private_static_rho();
}

double rho_A1(double A1, void * params)
{
   TParams *par = (TParams*)params;

   int i = par->i;
   
   par->lattice->copy2private(i);
   par->lattice->privateFields[4].A_1 = A1;

   par->lattice->compute_affected_private_derivatives(3);

   return par->lattice->affected_private_static_rho();
}

double rho_A2(double A2, void * params)
{
   TParams *par = (TParams*)params;

   int i = par->i;
   
   par->lattice->copy2private(i);
   par->lattice->privateFields[4].A_2 = A2;

   par->lattice->compute_affected_private_derivatives(4);

   return par->lattice->affected_private_static_rho();
}

double rho_A3(double A3, void * params)
{
   TParams *par = (TParams*)params;

   int i = par->i;
   
   par->lattice->copy2private(i);
   par->lattice->privateFields[4].A_3 = A3;

   par->lattice->compute_affected_private_derivatives(5);

   return par->lattice->affected_private_static_rho();
}

double rho(const gsl_vector *v, void *params)
{
   TParams *par = (TParams*)params;

   int i = par->i;
   
   par->lattice->copy2private(i);

   par->lattice->privateFields[4].th_1 = wrap(gsl_vector_get(v, 0));
   par->lattice->privateFields[4].th_2 = wrap(gsl_vector_get(v, 1));
   par->lattice->privateFields[4].th_3 = wrap(gsl_vector_get(v, 2));

   par->lattice->privateFields[4].A_1 = gsl_vector_get(v, 3);
   par->lattice->privateFields[4].A_2 = gsl_vector_get(v, 4);
   par->lattice->privateFields[4].A_3 = gsl_vector_get(v, 5);

   par->lattice->compute_affected_private_derivatives(0);
   par->lattice->compute_affected_private_derivatives(1);
   par->lattice->compute_affected_private_derivatives(2);
   par->lattice->compute_affected_private_derivatives(3);
   par->lattice->compute_affected_private_derivatives(4);
   par->lattice->compute_affected_private_derivatives(5);
   
   par->lattice->privateAux[4].compute_hG(par->lattice->privateFields[4].th_1, 
                                          par->lattice->privateFields[4].th_2, 
                                          par->lattice->privateFields[4].th_3);

   return par->lattice->affected_private_static_rho();
}

void drho(const gsl_vector *v, void *params, gsl_vector *result)
{
   TParams *par = (TParams*)params;

   int i = par->i;
   
   par->lattice->copy2private(i);

   par->lattice->privateFields[4].th_1 = wrap(gsl_vector_get(v, 0));
   par->lattice->privateFields[4].th_2 = wrap(gsl_vector_get(v, 1));
   par->lattice->privateFields[4].th_3 = wrap(gsl_vector_get(v, 2));

   par->lattice->privateFields[4].A_1 = gsl_vector_get(v, 3);
   par->lattice->privateFields[4].A_2 = gsl_vector_get(v, 4);
   par->lattice->privateFields[4].A_3 = gsl_vector_get(v, 5);

   par->lattice->compute_affected_private_derivatives(0);
   par->lattice->compute_affected_private_derivatives(1);
   par->lattice->compute_affected_private_derivatives(2);
   par->lattice->compute_affected_private_derivatives(3);
   par->lattice->compute_affected_private_derivatives(4);
   par->lattice->compute_affected_private_derivatives(5);
   
   par->lattice->privateAux[4].compute_hG(par->lattice->privateFields[4].th_1, 
                                          par->lattice->privateFields[4].th_2, 
                                          par->lattice->privateFields[4].th_3);
      
   TFields dotMomenta;

   par->lattice->compute_private_staticDotMomenta(dotMomenta);

   // dot p = - dH/dq

   gsl_vector_set(result, 0, -dotMomenta.th_1);
   gsl_vector_set(result, 1, -dotMomenta.th_2);
   gsl_vector_set(result, 2, -dotMomenta.th_3);
   gsl_vector_set(result, 3, -dotMomenta.A_1);
   gsl_vector_set(result, 4, -dotMomenta.A_2);
   gsl_vector_set(result, 5, -dotMomenta.A_3);
}

void rhodrho(const gsl_vector *v, void *params, double *rho, gsl_vector *drho) 
{
   TParams *par = (TParams*)params;

   int i = par->i;
   
   par->lattice->copy2private(i);

   par->lattice->privateFields[4].th_1 = wrap(gsl_vector_get(v, 0));
   par->lattice->privateFields[4].th_2 = wrap(gsl_vector_get(v, 1));
   par->lattice->privateFields[4].th_3 = wrap(gsl_vector_get(v, 2));

   par->lattice->privateFields[4].A_1 = gsl_vector_get(v, 3);
   par->lattice->privateFields[4].A_2 = gsl_vector_get(v, 4);
   par->lattice->privateFields[4].A_3 = gsl_vector_get(v, 5);

   par->lattice->compute_affected_private_derivatives(0);
   par->lattice->compute_affected_private_derivatives(1);
   par->lattice->compute_affected_private_derivatives(2);
   par->lattice->compute_affected_private_derivatives(3);
   par->lattice->compute_affected_private_derivatives(4);
   par->lattice->compute_affected_private_derivatives(5);
   
   par->lattice->privateAux[4].compute_hG(par->lattice->privateFields[4].th_1, 
                                          par->lattice->privateFields[4].th_2, 
                                          par->lattice->privateFields[4].th_3);

   *rho = par->lattice->affected_private_static_rho();

   TFields dotMomenta;

   par->lattice->compute_private_staticDotMomenta(dotMomenta);

   // dot p = - dH/dq

   gsl_vector_set(drho, 0, -dotMomenta.th_1);
   gsl_vector_set(drho, 1, -dotMomenta.th_2);
   gsl_vector_set(drho, 2, -dotMomenta.th_3);
   gsl_vector_set(drho, 3, -dotMomenta.A_1);
   gsl_vector_set(drho, 4, -dotMomenta.A_2);
   gsl_vector_set(drho, 5, -dotMomenta.A_3);
}

TOptimizeLattice::TOptimizeLattice()
{   
   par.lattice = this;

   f_rho.f = &rho;
   f_rho.df = &drho;
   f_rho.fdf = &rhodrho;
   f_rho.n = 6;
   f_rho.params = &par;

   par.f_rho_th1.function = &rho_th1; par.f_rho_th1.params = &par;
   par.f_rho_th2.function = &rho_th2; par.f_rho_th2.params = &par;
   par.f_rho_th3.function = &rho_th3; par.f_rho_th3.params = &par;
   par.f_rho_A1.function = &rho_A1;   par.f_rho_A1.params = &par;
   par.f_rho_A2.function = &rho_A2;   par.f_rho_A2.params = &par;
   par.f_rho_A3.function = &rho_A3;   par.f_rho_A3.params = &par;

   v = gsl_vector_alloc(6);

   T = gsl_multimin_fdfminimizer_vector_bfgs2;
   s = gsl_multimin_fdfminimizer_alloc(T, 6);
}

TOptimizeLattice::~TOptimizeLattice()
{
   gsl_multimin_fdfminimizer_free(s);
   gsl_vector_free(v);
}

void TOptimizeLattice::minimize(const double step_size, const double tol, const double grad_tol, const size_t max_iter,
                                const int from_i, const int to_i)
{
   int inActive;
   if (0 == active) inActive = 1; else inActive = 0;

   for (int i = from_i; i <= to_i; i++)
   {
       par.i = i;

       gsl_vector_set(v, 0, fields[active][i].th_1);
       gsl_vector_set(v, 1, fields[active][i].th_2);
       gsl_vector_set(v, 2, fields[active][i].th_3);
       gsl_vector_set(v, 3, fields[active][i].A_1);
       gsl_vector_set(v, 4, fields[active][i].A_2);
       gsl_vector_set(v, 5, fields[active][i].A_3);

       size_t iter = 0;
       int status;

       gsl_multimin_fdfminimizer_set(s, &f_rho, v, step_size, tol);

       do
       {
          iter++;
          status = gsl_multimin_fdfminimizer_iterate(s);

          if (status) break;

          status = gsl_multimin_test_gradient(s->gradient, grad_tol);
       }
       while ((status == GSL_CONTINUE) && (iter < max_iter));
               
       if (iter == max_iter) cout << "Iteration limit hit at " << i << endl << flush;
               
       fields[inActive][i].th_1 = wrap(gsl_vector_get(s->x, 0));
       fields[inActive][i].th_2 = wrap(gsl_vector_get(s->x, 1));
       fields[inActive][i].th_3 = wrap(gsl_vector_get(s->x, 2));
       fields[inActive][i].A_1 = gsl_vector_get(s->x, 3);
       fields[inActive][i].A_2 = gsl_vector_get(s->x, 4);
       fields[inActive][i].A_3 = gsl_vector_get(s->x, 5);
   }

   // Done; swap grids
   
   active = inActive;
}

int main(int argc, char* argv[])
{
   cout.setf(ios::left, ios::adjustfield);
   cout.precision(OUTPUT_PRECISION);

   // Set default parameter values

   double rho_rel_tol = DEFAULT_RHO_REL_TOL;

   double step_size   = DEFAULT_STEP_SIZE;
   double tol         = DEFAULT_TOL;
   double grad_tol    = DEFAULT_GRAD_TOL;
   size_t max_iter    = DEFAULT_MAX_ITER;

   char cBuff[128];               // Character buffer for integer to string conversion

   string sFile("signal.txt");    // Default signal filename
   string sBuff("");              // Used to hold input from signal file

   string iFile("relaxed");       // Default initial VTK filename
   string rFile("optimized");     // Default minimized VTK filename
   string rTitle("Optimized");    // Default minimized VTK title

   // Parse command line

   int i = 1;
   while (i < argc)
   {
      if (!strcmp("-sf", argv[i])) sFile.assign(argv[i + 1]); else
      if (!strcmp("-if", argv[i])) iFile.assign(argv[i + 1]); else
      if (!strcmp("-of", argv[i])) rFile.assign(argv[i + 1]); else

      if (!strcmp("-ot", argv[i])) rTitle.assign(argv[i + 1]); else

      if (!strcmp("-mi", argv[i])) max_iter = atol(argv[i + 1]); else
      if (!strcmp("-ss", argv[i])) step_size = atof(argv[i + 1]); else
      if (!strcmp("-at", argv[i])) tol = atof(argv[i + 1]); else
      if (!strcmp("-gt", argv[i])) grad_tol = atof(argv[i + 1]); else
      if (!strcmp("-rr", argv[i])) rho_rel_tol = atof(argv[i + 1]); else
      
      {
         cout << "Unrecognized command line argument: " << argv[i] << endl;
         cout << endl;
         cout << "Supported options [defaults in square brackets]:" << endl;
         cout << endl;
         cout << "Filename for signal file:                -sf [\"" << sFile << "\"]" << endl;
         cout << "Filename for initial data input:         -if [\"" << iFile << "\"]" << endl;
         cout << "Filename for optimized data output:      -of [\"" << rFile << "\"]" << endl;
         cout << endl;
         cout << "Title of optimized data file:            -ot [\"" << rTitle << "\"]" << endl;
         cout << endl;
         cout << "Max number of optimizer iterations:      -mi [" << max_iter << "]" << endl;
         cout << "Step size:                               -ss [" << step_size << "]" << endl;
         cout << "Absolute rho tolerance (local):          -at [" << tol << "]" << endl;
         cout << "Absolute rho gradient tolerance (local): -gt [" << grad_tol << "]" << endl;
         cout << "Relative rho tolerance (global):         -rr [" << rho_rel_tol << "]" << endl;

         exit(1);
      }

      i += 2;
   }

   cout << "mi: " << max_iter << endl;
   cout << "ss: " << step_size << endl;
   cout << "at: " << tol << endl;
   cout << "gt: " << grad_tol << endl; 
   cout << "rr: " << rho_rel_tol << endl << endl << flush;
   
   // Initialize lattice

   TOptimizeLattice lattice[N_THREADS];
   lattice[0].clear(1);

   // Load data

   cout << "Loading initial data." << endl << flush;

   string iTitle;
   lattice[0].read_vtk(iFile + ".vtk", iTitle);

   cout << "Done loading." << endl << endl << flush;
   
   // Initialize
   
   int i_stride = lattice[0].get_nSide() / N_THREADS;

   if (lattice[0].get_nSide() != i_stride * N_THREADS)
   {
      cout << "Can not divide lattice side " << lattice[0].get_nSide() << " by " 
           << N_THREADS << " threads." << endl << flush;

       exit(2);
   }

   for (i = 1; i < N_THREADS; i++) 
   {
      lattice[i].set_nSide(lattice[0].get_nSide(), false);
   }

   lattice[0].init_aux_static();
      
   double rho_rel = 1e6;
   double prev_sr, sr = lattice[0].sum_rho_static();

   long saveAfterSweepsDone = 0;
   long lastSavedSweep = 0;
   long nSweepsDone = 0;

   cout << "Sweep, energy, relative energy:" << endl;
   cout << nSweepsDone << ", " << sr << ", 1.0" << endl << flush;

   do
   {
      prev_sr = sr;

#ifdef _OPENMP
#pragma omp parallel
{
#pragma omp for
#endif
      for (i = 0; i < N_THREADS; i++)
      {
         lattice[i].minimize(step_size, tol, grad_tol, max_iter, i * i_stride, (i + 1) * i_stride - 1);
      }
#ifdef _OPENMP
#pragma omp barrier
}
#endif

      lattice[0].init_aux_static();
      sr = lattice[0].sum_rho_static();

      nSweepsDone++;

      cout << nSweepsDone << ", " << sr << ", " << sr/prev_sr << endl << flush;

      if ("" != sFile)
      {
         ifstream fin(sFile.c_str());

         if (fin)
         {
            sBuff = "";

            fin >> sBuff;

            if ("abort" == sBuff) { fin.close(); return(1); }
            if ("quit" == sBuff)  { fin.close(); break; }
            if ("save" == sBuff)
            {
               fin >> sBuff;
               saveAfterSweepsDone = atol(sBuff.c_str());

               if (0 == saveAfterSweepsDone) saveAfterSweepsDone = nSweepsDone;
            }
         }

         fin.close();
      }

      if ((lastSavedSweep < saveAfterSweepsDone) && (saveAfterSweepsDone <= nSweepsDone))
      {
         saveAfterSweepsDone = 0;

         sprintf(cBuff, "%lu", nSweepsDone);
         sBuff = rFile + "_" + cBuff;

         cout << "Saving optimized data to \"" << sBuff << ".vtk\"." << endl << flush;

         lattice[0].write_vtk(sBuff + ".vtk", rTitle, nSweepsDone);

         lastSavedSweep = nSweepsDone;

         cout << "Done saving." << endl << flush;
      }

      if (sr <= prev_sr) rho_rel = 1.0 - sr/prev_sr; else break;
   }
   while (rho_rel > rho_rel_tol);
   
   cout << endl << flush;

   if (sr > prev_sr)
   {
      if (0 == lattice[0].get_active()) lattice[0].set_active(1); else lattice[0].set_active(0);

      nSweepsDone--;

      cout << "Backtracking." << endl << endl << flush;      
   }

   cout << "Saving optimized data to \"" << rFile << ".vtk\"." << endl << flush;

   lattice[0].write_vtk(rFile + ".vtk", rTitle, nSweepsDone);

   cout << "Done saving." << endl << flush;

   return(0);
}
