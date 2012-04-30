/* relax_1D.cpp : Reduces local energy density by low pass filtering (1D).
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
const int N_THREADS                         = 4;      /* Must divide lattice side evenly, ideally = number of CPU cores. */
#endif

const int DEFAULT_NSWEEPS                   = 10;     /* Default number of sweeps per relax call */                                                                                                                                                                
const double DEFAULT_RHO_REL_TOL            = 1e-6;   /* Default relative tolerance in rho */
const double DEFAULT_WEIGHT_THRESHOLD_TOL   = 1e-4;   /* Default relative error in rho under which central weight is increased */
const double DEFAULT_CENTRAL_WEIGHT_STEP    = 0.01;   /* Amount central weight is increased */

TFields fields[2][MAX_NSIDE];                         /* Max allowed number of cells always allocated. Ugly but fast. */
double affectedRho[MAX_NSIDE];                        /* Ditto */
TAux aux[MAX_NSIDE];                                  /* Ditto */

#include "lattice_1d.cpp"

using namespace std;

class TRelaxLattice : public TLattice
{
public:

   void relax_lowpass(int from_i, int to_i);
   void init_affectedRho();
};

void TRelaxLattice::relax_lowpass(int from_i, int to_i)
/* IMPORTANT: Assumes all h, G matrices and derivatives are ready for use. */
{
   int inActive;
   
   if (0 == active) inActive = 1; else inActive = 0;

   for (int i = from_i; i <= to_i; i++)
   {
      copy2private(i);
           
      double min_rho = affectedRho[i];
      enum { conf_original, conf_theta, conf_a, conf_both } min_conf = conf_original;

      // Get lowpassed field values

      TFields lowpassed;
                   
      lowpass.compute(privateFields, &lowpassed);
                 
      // Wiggle field values, compute resulting energy

      // Case #1: only theta lowpassed

      privateFields[4].th_1 = lowpassed.th_1;
      privateFields[4].th_2 = lowpassed.th_2;
      privateFields[4].th_3 = lowpassed.th_3;

      compute_affected_private_derivatives(0);
      compute_affected_private_derivatives(1);
      compute_affected_private_derivatives(2);
                  
      privateAux[4].compute_hG(lowpassed.th_1, lowpassed.th_2, lowpassed.th_3);

      double wiggled_rho = affected_private_static_rho();

      if (wiggled_rho < min_rho)
      {
         min_rho = wiggled_rho;
         min_conf = conf_theta;
      }

      // Case #2: theta and A both lowpassed

      privateFields[4].A_1 = lowpassed.A_1;
      privateFields[4].A_2 = lowpassed.A_2;
      privateFields[4].A_3 = lowpassed.A_3;

      compute_affected_private_derivatives(3);
      compute_affected_private_derivatives(4);
      compute_affected_private_derivatives(5);

      wiggled_rho = affected_private_static_rho();

      if (wiggled_rho < min_rho)
      {
         min_rho = wiggled_rho;
         min_conf = conf_both;
      }

      // Case #3: only A lowpassed

      privateFields[4].th_1 = fields[active][i].th_1;
      privateFields[4].th_2 = fields[active][i].th_2;
      privateFields[4].th_3 = fields[active][i].th_3;
                  
      copy2affected_private_derivatives(i, 0);
      copy2affected_private_derivatives(i, 1);
      copy2affected_private_derivatives(i, 2);

      privateAux[4] = aux[i];

      wiggled_rho = affected_private_static_rho();

      if (wiggled_rho < min_rho)
      {
         min_rho = wiggled_rho;
         min_conf = conf_a;
      }

      // Copy lowest energy configuration to inActive grid

      switch (min_conf)
      {
         case conf_theta: {
                             fields[inActive][i].th_1 = lowpassed.th_1;
                             fields[inActive][i].th_2 = lowpassed.th_2;
                             fields[inActive][i].th_3 = lowpassed.th_3;

                             fields[inActive][i].A_1 = fields[active][i].A_1;
                             fields[inActive][i].A_2 = fields[active][i].A_2;
                             fields[inActive][i].A_3 = fields[active][i].A_3;

                             affectedRho[i] = min_rho;

                             break;
                          }
         case conf_a:     {
                             fields[inActive][i].th_1 = fields[active][i].th_1;
                             fields[inActive][i].th_2 = fields[active][i].th_2;
                             fields[inActive][i].th_3 = fields[active][i].th_3;

                             fields[inActive][i].A_1 = lowpassed.A_1;
                             fields[inActive][i].A_2 = lowpassed.A_2;
                             fields[inActive][i].A_3 = lowpassed.A_3;

                             affectedRho[i] = min_rho;

                             break;
                          }
         case conf_both:  {
                             fields[inActive][i].th_1 = lowpassed.th_1;
                             fields[inActive][i].th_2 = lowpassed.th_2;
                             fields[inActive][i].th_3 = lowpassed.th_3;

                             fields[inActive][i].A_1 = lowpassed.A_1;
                             fields[inActive][i].A_2 = lowpassed.A_2;
                             fields[inActive][i].A_3 = lowpassed.A_3;

                             affectedRho[i] = min_rho;

                             break;
                          }
         default :        {
                             fields[inActive][i] = fields[active][i];
                          }
      }
   }

   // Sweep done, swap grids

   active = inActive;
}

void TRelaxLattice::init_affectedRho()
/* IMPORTANT: Assumes all h, G matrices and derivatives are ready for use. */
{
#ifdef _OPENMP
#pragma omp parallel default(shared)
{
#pragma omp for
#endif
       for (int i = 0; i <= nSideM1; i++)
           affectedRho[i] = affected_static_rho(i);
#ifdef _OPENMP
#pragma omp barrier
}
#endif
}

int main(int argc, char *argv[])
{
   cout.setf(ios::left, ios::adjustfield);
   cout.precision(OUTPUT_PRECISION);

   // Set default parameter values

   long nSweeps = DEFAULT_NSWEEPS;
   
   double rho_rel_tol = DEFAULT_RHO_REL_TOL;
   double central_weight_step = DEFAULT_CENTRAL_WEIGHT_STEP;
   double weight_threshold_tol = DEFAULT_WEIGHT_THRESHOLD_TOL;

   char cBuff[128];               // Character buffer for integer to string conversion

   string sFile("signal.txt");    // Default signal filename
   string sBuff("");              // Used to hold input from signal file

   string iFile("initial");       // Default initial VTK filename
   string rFile("relaxed");       // Default relaxed VTK filename
   string rTitle("Relaxed");      // Default relaxed VTK title

   // Parse command line

   int i = 1;
   while (i < argc)
   {
      if (!strcmp("-sf", argv[i])) sFile.assign(argv[i + 1]); else
      if (!strcmp("-if", argv[i])) iFile.assign(argv[i + 1]); else
      if (!strcmp("-of", argv[i])) rFile.assign(argv[i + 1]); else

      if (!strcmp("-ot", argv[i])) rTitle.assign(argv[i + 1]); else

      if (!strcmp("-ns", argv[i])) nSweeps = atol(argv[i + 1]); else
      if (!strcmp("-rr", argv[i])) rho_rel_tol = atof(argv[i + 1]); else
      if (!strcmp("-wt", argv[i])) weight_threshold_tol = atof(argv[i + 1]); else
      if (!strcmp("-ws", argv[i])) central_weight_step = atof(argv[i + 1]); else
      {
         cout << "Unrecognized command line argument: " << argv[i] << endl;
         cout << endl;
         cout << "Supported options [defaults in square brackets]:" << endl;
         cout << endl;
         cout << "Filename for signal file:                    -sf [\"" << sFile << "\"]" << endl;
         cout << "Filename for initial data output:            -if [\"" << iFile << "\"]" << endl;
         cout << "Filename for relaxed data output:            -of [\"" << rFile << "\"]" << endl;
         cout << endl;
         cout << "Title of relaxed data file:                  -ot [\"" << rTitle << "\"]" << endl;
         cout << endl;
         cout << "Number of sweeps between convergence checks: -ns [" << nSweeps << "]" << endl;
         cout << endl;
         cout << "Relative energy tolerance:                   -rr [" << rho_rel_tol << "]" << endl;
         cout << "Relative central weight stepping threshold:  -wt [" << weight_threshold_tol << "]" << endl;
         cout << "Central weight step:                         -ws [" << central_weight_step << "]" << endl;

         exit(1);
      }

      i += 2;
   }

   cout << "Relative energy tolerance:                  -rr = " << rho_rel_tol << endl;
   cout << "Relative central weight stepping threshold: -wt = " << weight_threshold_tol << endl;
   cout << "Central weight step:                        -ws = " << central_weight_step << endl;
   cout << endl;
   
   // Initialize lattice

   TRelaxLattice lattice[N_THREADS];
   lattice[0].clear(1);

   // Load seed data

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
   lattice[0].init_affectedRho();
   
   double rho_rel = 1e6;
   double prev_sr, sr = lattice[0].sum_rho_static();
  
   long saveAfterSweepsDone = 0;
   long lastSavedSweep = 0;
   long nSweepsDone = 0;

   long i_central_weight = 0;
   double central_weight = i_central_weight * central_weight_step;

   cout << "Sweep, central weight, energy, relative energy:" << endl;
   cout << nSweepsDone << ", "  << central_weight << ", " << sr << ", 1.0" << endl << flush;
   
   do
   {
      prev_sr = sr;

      for (int sweep = 0; sweep < nSweeps; sweep++)
      {
#ifdef _OPENMP
#pragma omp parallel
{
#pragma omp for
#endif
         for (int i = 0; i < N_THREADS; i++)
         {
            lattice[i].relax_lowpass(i * i_stride, (i + 1) * i_stride - 1);
         }
#ifdef _OPENMP
#pragma omp barrier
}
#endif

         lattice[0].init_aux_static();
      }

      sr = lattice[0].sum_rho_static();

      nSweepsDone += nSweeps;

      cout << nSweepsDone << ", "  << central_weight << ", " << sr << ", " << sr/prev_sr << endl << flush;

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

         cout << "Saving relaxed data to \"" << sBuff << ".vtk\"." << endl << flush;

         lattice[0].write_vtk(sBuff + ".vtk", rTitle, nSweepsDone);

         lastSavedSweep = nSweepsDone;

         cout << "Done saving." << endl << flush;
      }

      if (sr < prev_sr) 
      {
         rho_rel = 1.0 - sr/prev_sr;
            
         if ((rho_rel < weight_threshold_tol) && (central_weight < 1.0))
         {
            i_central_weight++;
            central_weight = i_central_weight * central_weight_step;

           for (int i = 0; i < N_THREADS; i++)
               lattice[i].set_central_weight(central_weight);
         }
      }
      else break;
   }
   while (rho_rel > rho_rel_tol);
   
   cout << endl << flush;

   cout << "Saving relaxed data to \"" << rFile << ".vtk\"." << endl << flush;

   lattice[0].write_vtk(rFile + ".vtk", rTitle, nSweepsDone);

   cout << "Done saving." << endl << flush;

   return(0);
}
