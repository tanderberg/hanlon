/* seed_1d.cpp : Creates random field configurations (1D).
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
#include "math_aux.cpp"

#ifdef USE_MKL
#include "mkl_vsl.h"
#endif

#ifdef _OPENMP
#include <omp.h>
#endif

#ifndef MAX_NSIDE
const int MAX_NSIDE                 = 8196;          /* Max number of lattice sites per side.
                                                        Power of 2 for faster indexing (maybe). */
#endif

const int DEFAULT_NSIDE             = 1024;

#ifndef EPSILON
const double EPSILON                = 1.0e-10;
#endif

const double DEFAULT_BOUNDARY_WIDTH = 20;            /* Default initial width (in lattice units) of cell boundaries */
const double MAX_INITIAL_A          = 10;            /* Max initial |A| */

const int OUTPUT_PRECISION          = 15;

const int DEFAULT_NCELLS            = 10;            /* Default initial number of Voronoi cells in grid */
const int NCOMPS                    = 7;             /* Number of field components used INTERNALLY */

typedef double TCube[NCOMPS][MAX_NSIDE];             /* First index = field (n_0..n_3, A_1..A_3).
                                                                                                                                                   Second index: x position of cell */
TCube nA;

#ifdef USE_MKL
VSLStreamStatePtr pVslStream = NULL;
#endif

struct TRndm
{
   double v;

   void add(double x) { v += x; }
   void sub(double x) { v -= x; }
   void mul(double x) { v *= x; }
   void div(double x) { v /= x; }

   double sqr()     { return(v*v); }
   double norm()    { return(sqrt(sqr())); }
   void normalize() { double r = norm(); v /= r; }

   void boxin(double lo, double hi, double side)
   {
      while (v > hi) v -= side;
      while (v < lo) v += side;
   }

   void halfline(double side) 
   {
#ifdef USE_MKL
      vdRngUniform(VSL_METHOD_DUNIFORM_STD, pVslStream, 1, &v, ZERO, side);
#else
      v = rand() * side / RAND_MAX;
#endif
   }
   
   TRndm() { v = 0; }
};

struct TRndm3
{
   double v[3];

   void add(double x) { v[0] += x; v[1] += x; v[2] += x; }
   void sub(double x) { v[0] -= x; v[1] -= x; v[2] -= x; }
   void mul(double x) { v[0] *= x; v[1] *= x; v[2] *= x; }
   void div(double x) { v[0] /= x; v[1] /= x; v[2] /= x; }

   double sqr()     { return(v[0]*v[0] + v[1]*v[1] + v[2]*v[2]); }
   double norm()    { return(sqrt(sqr())); }
   void normalize() { double r = norm(); v[0] /= r; v[1] /= r; v[2] /= r; }

   void boxin(double lo, double hi, double side)
   {
      while (v[0] > hi) v[0] -= side;
      while (v[0] < lo) v[0] += side;
      while (v[1] > hi) v[1] -= side;
      while (v[1] < lo) v[1] += side;
      while (v[2] > hi) v[2] -= side;
      while (v[2] < lo) v[2] += side;
   }

   void octant(double side) 
   {
#ifdef USE_MKL
      vdRngUniform(VSL_METHOD_DUNIFORM_STD, pVslStream, 3, v, ZERO, side);
#else
      v[0] = rand() * side / RAND_MAX;
      v[1] = rand() * side / RAND_MAX;
      v[2] = rand() * side / RAND_MAX;
#endif
   }

   void ball(double radius, double diameter, double radiusSquared)
   {
      do
      {
#ifdef USE_MKL

         vdRngUniform(VSL_METHOD_DUNIFORM_STD, pVslStream, 3, v, -radius, radius);
#else
         v[0] = rand() * diameter / RAND_MAX - radius;
         v[1] = rand() * diameter / RAND_MAX - radius;
         v[2] = rand() * diameter / RAND_MAX - radius;
#endif
      }
      while (sqr() > radiusSquared);
   }

   void sphere(double radius, double diameter, double radiusSquared)
   {
      double r;

      do
      {
#ifdef USE_MKL
         vdRngUniform(VSL_METHOD_DUNIFORM_STD, pVslStream, 3, v, -radius, radius);
#else
         v[0] = rand() * diameter / RAND_MAX - radius;
         v[1] = rand() * diameter / RAND_MAX - radius;
         v[2] = rand() * diameter / RAND_MAX - radius;
#endif
         r = sqr();
      }
      while (r > radiusSquared);

      r = sqrt(radiusSquared / r);

      v[0] *= r;
      v[1] *= r;
      v[2] *= r;
   }
};

struct TRndm4
{
   double v[4];

   void add(double x) { v[0] += x; v[1] += x; v[2] += x; v[3] += x; }
   void sub(double x) { v[0] -= x; v[1] -= x; v[2] -= x; v[3] -= x; }
   void mul(double x) { v[0] *= x; v[1] *= x; v[2] *= x; v[3] *= x; }
   void div(double x) { v[0] /= x; v[1] /= x; v[2] /= x; v[3] /= x; }

   double sqr()     { return(v[0]*v[0] + v[1]*v[1] + v[2]*v[2] + v[3]*v[3]); }
   double norm()    { return(sqrt(sqr())); }
   void normalize() { double r = norm(); v[0] /= r; v[1] /= r; v[2] /= r; v[3] /= r; }

   void boxin(double lo, double hi, double side)
   {
      while (v[0] > hi) v[0] -= side;
      while (v[0] < lo) v[0] += side;
      while (v[1] > hi) v[1] -= side;
      while (v[1] < lo) v[1] += side;
      while (v[2] > hi) v[2] -= side;
      while (v[2] < lo) v[2] += side;
      while (v[3] > hi) v[3] -= side;
      while (v[3] < lo) v[3] += side;
   }

   void ball(double radius, double diameter, double radiusSquared)
   {
      do
      {
#ifdef USE_MKL
         vdRngUniform(VSL_METHOD_DUNIFORM_STD, pVslStream, 4, v, -radius, radius);
#else
         v[0] = rand() * diameter / RAND_MAX - radius;
         v[1] = rand() * diameter / RAND_MAX - radius;
         v[2] = rand() * diameter / RAND_MAX - radius;
         v[3] = rand() * diameter / RAND_MAX - radius;
#endif
      }
      while (sqr() > radiusSquared);
   }

   void sphere(double radius, double diameter, double radiusSquared)
   {
      double r;

      do
      {
#ifdef USE_MKL
         vdRngUniform(VSL_METHOD_DUNIFORM_STD, pVslStream, 4, v, -radius, radius);
#else
         v[0] = rand() * diameter / RAND_MAX - radius;
         v[1] = rand() * diameter / RAND_MAX - radius;
         v[2] = rand() * diameter / RAND_MAX - radius;
         v[3] = rand() * diameter / RAND_MAX - radius;
#endif
         r = sqr();
      }
      while (r > radiusSquared);

      r = sqrt(radiusSquared / r);

      v[0] *= r; 
      v[1] *= r; 
      v[2] *= r;
      v[3] *= r;
   }
};

struct TCell
{
   double y[NCOMPS];
   double i;
};

void voronoi(int nSide, int nCells, double boundaryWidth, double aFactorFactor)
/*
   boundaryWidth is expressed as a multiple of lattice spacing, not in physical units
*/
{
   // Create list of random cell centers with random n_a and magnetic potential

   TRndm rndm;
   TRndm3 rndm3;
   TRndm4 rndm4;

   double aFactor = TWO * MAX_INITIAL_A * aFactorFactor;

   TCell* cell = new TCell[nCells];

   // Do this outside parallel section so same random seed always yields same initial configuration

   for (int i = 0; i < nCells; i++)
   {
      rndm.halfline(nSide);

      cell[i].i = rndm.v;

      rndm4.sphere(ONE, TWO, ONE);

      cell[i].y[0] = rndm4.v[0];
      cell[i].y[1] = rndm4.v[1];
      cell[i].y[2] = rndm4.v[2];
      cell[i].y[3] = rndm4.v[3];
   }

   // Create random A after n so that same random seed yields same n cells whether A is quenched or not

   for (int i = 0; i < nCells; i++)
   {
      rndm3.ball(HALF, ONE, QUARTER);
      rndm3.mul(aFactor);

      cell[i].y[4] = rndm3.v[0];
      cell[i].y[5] = rndm3.v[1];
      cell[i].y[6] = rndm3.v[2];
   }

#ifdef _OPENMP
#pragma omp parallel default(shared) private(rndm, rndm4)
{
   // Loop through grid and set its field values as spline interpolation between two closest cell centers

#pragma omp for
#endif
   for (int i = 0; i < nSide; i++)
   {
      // Determine two closest cell centers

      double mrR2 = TWO*nSide*nSide; int m =-1; // Closest
      double nrR2 = TWO*nSide*nSide; int n =-1; // Next closest

      for (int l = 0; l < nCells; l++)
      {
         rndm.v = fabs(cell[l].i - i); rndm.v = min(rndm.v, nSide - rndm.v);

         double rR2 = rndm.sqr();

         if (rR2 < mrR2)
         {
            nrR2 = mrR2;
            n = m;

            mrR2 = rR2;
            m = l;
         }
         else if (rR2 < nrR2)
         {
            nrR2 = rR2;
            n = l;
         }
      }

      // Distance parameter for spline is computed from midpoint between cell centers

      double r = spline(sqrt(mrR2) - sqrt(nrR2), boundaryWidth);

      // Interpolate n_a

      rndm4.v[0] = (cell[m].y[0] - cell[n].y[0])*r + cell[n].y[0];
      rndm4.v[1] = (cell[m].y[1] - cell[n].y[1])*r + cell[n].y[1];
      rndm4.v[2] = (cell[m].y[2] - cell[n].y[2])*r + cell[n].y[2];
      rndm4.v[3] = (cell[m].y[3] - cell[n].y[3])*r + cell[n].y[3];

      rndm4.normalize();

      nA[0][i] = rndm4.v[0];
      nA[1][i] = rndm4.v[1];
      nA[2][i] = rndm4.v[2];
      nA[3][i] = rndm4.v[3];

      // Interpolate A

      nA[4][i] = (cell[m].y[4] - cell[n].y[4])*r + cell[n].y[4];
      nA[5][i] = (cell[m].y[5] - cell[n].y[5])*r + cell[n].y[5];
      nA[6][i] = (cell[m].y[6] - cell[n].y[6])*r + cell[n].y[6];
   }
           
#ifdef _OPENMP
#pragma omp barrier
}
#endif

   delete [] cell;
}

void pulses(int nSide, int nCells, double boundaryWidth, double aFactorFactor)
/*
   boundaryWidth is expressed as a multiple of lattice spacing, not in physical units
*/
{
   // Create list of random cell centers with random n_a and magnetic potential

   TRndm rndm;
   TRndm3 rndm3;
   TRndm4 rndm4;

   double aFactor = TWO * MAX_INITIAL_A * aFactorFactor;

   TCell* cell = new TCell[nCells];

   // Do this outside parallel section so same random seed always yields same initial configuration

   for (int i = 0; i < nCells; i++)
   {
      rndm.halfline(nSide);

      if (1 == nCells)
      {
         cell[i].i = nSide / 2;
      }
      else
      {
         cell[i].i = rndm.v;
      }

      rndm4.sphere(ONE, TWO, ONE);

      cell[i].y[0] = rndm4.v[0];
      cell[i].y[1] = rndm4.v[1];
      cell[i].y[2] = rndm4.v[2];
      cell[i].y[3] = rndm4.v[3];
   }

   // Create random A after n so that same random seed yields same n cells whether A is quenched or not

   for (int i = 0; i < nCells; i++)
   {
      rndm3.ball(HALF, ONE, QUARTER);
      rndm3.mul(aFactor);

      cell[i].y[4] = rndm3.v[0];
      cell[i].y[5] = rndm3.v[1];
      cell[i].y[6] = rndm3.v[2];
   }

#ifdef _OPENMP
#pragma omp parallel default(shared) private(rndm3, rndm4)
{
   // Loop through grid and set its field values as spline interpolation between closest cell center and 0

#pragma omp for
#endif
   for (int i = 0; i < nSide; i++)
   {
      // Determine closest cell center

      double mrR2 = TWO*nSide*nSide;
      double rR2;

      int m = -1; // Closest

      for (int l = 0; l < nCells; l++)
      {
         rndm.v = fabs(cell[l].i - i); rndm.v = min(rndm.v, nSide - rndm.v);

         rR2 = rndm.sqr();

         if (rR2 < mrR2)
         {
            mrR2 = rR2;
            m = l;
         }
      }

      // Spline interpolation from cell center to 0

      double r = spline(sqrt(mrR2) - HALF * boundaryWidth, boundaryWidth);

      // Interpolate n_a

      rndm4.v[0] = cell[m].y[0] * r + (ONE - r);
      rndm4.v[1] = cell[m].y[1] * r;
      rndm4.v[2] = cell[m].y[2] * r;
      rndm4.v[3] = cell[m].y[3] * r;

      rndm4.normalize();

      nA[0][i] = rndm4.v[0];
      nA[1][i] = rndm4.v[1];
      nA[2][i] = rndm4.v[2];
      nA[3][i] = rndm4.v[3];

      // Interpolate A

      nA[4][i] = cell[m].y[4] * r;
      nA[5][i] = cell[m].y[5] * r;
      nA[6][i] = cell[m].y[6] * r;
   }

#ifdef _OPENMP
#pragma omp barrier
}
#endif

   delete [] cell;
}

void periodic(int nSide, int nx, double aFactorFactor)
{
   // Create 4 cells (central + one for each axis) with random n_a and magnetic potential

   TRndm3 rndm3;
   TRndm4 rndm4;

   double aFactor = TWO * MAX_INITIAL_A * aFactorFactor;

   TCell* cell = new TCell[2];

   // Do this outside parallel section so same random seed always yields same initial configuration

   for (int i = 0; i < 2; i++)
   {
      cell[i].i = nSide / 2;

      if (0 == i)
      {
         // Central cell. Store theta values rather than n_a for better interpolation.

         rndm4.sphere(ONE, TWO, ONE);

         double thetaHalf = acos(rndm4.v[0]);
         double thetaOverSinThetaHalf = TWO;

         if (fabs(thetaHalf) > EPSILON) thetaOverSinThetaHalf *= thetaHalf / sin(thetaHalf);

         cell[0].y[0] = thetaOverSinThetaHalf;
         cell[0].y[1] = thetaOverSinThetaHalf * rndm4.v[1];
         cell[0].y[2] = thetaOverSinThetaHalf * rndm4.v[2];
         cell[0].y[3] = thetaOverSinThetaHalf * rndm4.v[3];
      }
      else
      {
         // Set off-center n_a values antipodal to central cell's
      
         cell[i].y[0] = cell[0].y[0];
         cell[i].y[1] = wrap(cell[0].y[1] + 2*PI);
         cell[i].y[2] = wrap(cell[0].y[2] + 2*PI);
         cell[i].y[3] = wrap(cell[0].y[3] + 2*PI);
      }

      rndm3.ball(HALF, ONE, QUARTER);
      rndm3.mul(aFactor);

      cell[i].y[4] = rndm3.v[0];
      cell[i].y[5] = rndm3.v[1];
      cell[i].y[6] = rndm3.v[2];
   }
   
   double kx = (2.0 * nx) / nSide;

#ifdef _OPENMP
#pragma omp parallel default(shared) private(rndm4)
{
   // Loop through grid and interpolate between central cell and the three axis cells

#pragma omp for
#endif
   for (int i = 0; i < nSide; i++)
   {
       double r_ii;
       double r_i = (i - nSide / 2) * kx;
       r_i = fabs(modf(r_i, &r_ii));
	            
        // Interpolate theta

        rndm4.v[1] = wrap((1 - r_i) * cell[0].y[1] + r_i * cell[1].y[1]);
        rndm4.v[2] = wrap((1 - r_i) * cell[0].y[2] + r_i * cell[1].y[2]);
        rndm4.v[3] = wrap((1 - r_i) * cell[0].y[3] + r_i * cell[1].y[3]);

        // Turn it back into n_a
              
        double thR2 = sqr(rndm4.v[1]) + sqr(rndm4.v[2]) + sqr(rndm4.v[3]);

        if (thR2 < EPSILON) // O(th_i^3) series expansion 
        {
           rndm4.v[0] = 1 - thR2 / 8;
           rndm4.v[1] *= (24 - thR2) / 48;
           rndm4.v[2] *= (24 - thR2) / 48;
           rndm4.v[3] *= (24 - thR2) / 48;
        }
        else
        {
           double th = sqrt(thR2);
   
           rndm4.v[0] = cos(th / 2);
           rndm4.v[1] *= sin(th / 2) / th;
           rndm4.v[2] *= sin(th / 2) / th;
           rndm4.v[3] *= sin(th / 2) / th;
        }

        rndm4.normalize();

        nA[0][i] = rndm4.v[0];
        nA[1][i] = rndm4.v[1];
        nA[2][i] = rndm4.v[2];
        nA[3][i] = rndm4.v[3];

        // Interpolate A

        nA[4][i] = (1 - r_i) * cell[0].y[4] + r_i * cell[1].y[4];
        nA[5][i] = (1 - r_i) * cell[0].y[5] + r_i * cell[1].y[5];
        nA[6][i] = (1 - r_i) * cell[0].y[6] + r_i * cell[1].y[6];
   }

#ifdef _OPENMP
#pragma omp barrier
}
#endif

   delete [] cell;
}

bool write_vtk(std::string fileName, std::string title, int nSide, double dxn, double t)
{
   using namespace std;

   ofstream vtk;

   vtk.open(fileName.c_str());
   if (!vtk)
   {
      cout << "ERROR: Could not open VTK file " << fileName << endl << flush;
      return(false);
   }

   vtk.setf(ios::left, ios::adjustfield);
   vtk.precision(OUTPUT_PRECISION);

   vtk << "# vtk DataFile Version 3.0" << endl;
   vtk << title << "; t = " << t << endl;
   vtk << "ASCII" << endl;
   vtk << "DATASET STRUCTURED_POINTS" << endl;
   vtk << "DIMENSIONS " << nSide << endl;
   vtk << "ORIGIN 0" << endl;
   vtk << "SPACING " << dxn << endl;
   vtk << "POINT_DATA " << nSide << endl;

   vtk << "VECTORS T DOUBLE" << endl;

   for (int i = 0; i < nSide; i++)
   {
      double thetaHalf = acos(nA[0][i]);
      double thetaOverSinThetaHalf = TWO;

      if (fabs(thetaHalf) > EPSILON) thetaOverSinThetaHalf *= thetaHalf / sin(thetaHalf);

      vtk << thetaOverSinThetaHalf * nA[1][i] << ' ' 
          << thetaOverSinThetaHalf * nA[2][i] << ' ' 
          << thetaOverSinThetaHalf * nA[3][i] << endl;
   }

   vtk << "VECTORS A DOUBLE" << endl;

   for (int i = 0; i < nSide; i++)
   {
      vtk << nA[4][i] << ' ' 
          << nA[5][i] << ' ' 
          << nA[6][i] << endl;
   }

   vtk << "VECTORS pT DOUBLE" << endl;

   for (int i = 0; i < nSide; i++)
   {
      vtk << 0 << endl;
   }

   vtk << "VECTORS pA DOUBLE" << endl;

   for (int i = 0; i < nSide; i++)
   {
      vtk << 0 << endl;
   }

   vtk.close();

   return(true);
}

int main(int argc, char* argv[])
{
    using namespace std;

    cout.setf(ios::left, ios::adjustfield);
    cout.precision(OUTPUT_PRECISION);

    long rs              = 0;            // Random generator seed

    int nSide            = DEFAULT_NSIDE;
    int nCells           = DEFAULT_NCELLS;
    int nx               = 1;
    int mode             = 0;
    double boundaryWidth = DEFAULT_BOUNDARY_WIDTH;
    double aFactorFactor = ONE;

    string iFile("initial");             // Default initial VTK filename (.VTK extension implied)
    string iTitle("Initial");            // Default initial VTK title

    // Parse command line

    long int i = 1;
    while (i < argc)
    {
       if (!strcmp("-rs", argv[i])) rs = atol(argv[i + 1]); else
       if (!strcmp("-of", argv[i])) iFile.assign(argv[i + 1]); else
       if (!strcmp("-ot", argv[i])) iTitle.assign(argv[i + 1]); else
       if (!strcmp("-ns", argv[i])) nSide = atol(argv[i + 1]); else
       if (!strcmp("-om", argv[i])) mode = atol(argv[i + 1]); else
       if (!strcmp("-nc", argv[i])) nCells = atol(argv[i + 1]); else
       if (!strcmp("-bw", argv[i])) boundaryWidth = atof(argv[i + 1]); else
       if (!strcmp("-af", argv[i])) aFactorFactor = atof(argv[i + 1]); else
       if (!strcmp("-nx", argv[i])) nx = atol(argv[i + 1]); else
       {
          cout << "Unrecognized command line argument: " << argv[i] << endl;
          cout << endl;
          cout << "Supported options [defaults in square brackets]:" << endl;
          cout << endl;
          cout << "Random seed (long integer):                           -rs [" << rs << "]" << endl;
          cout << endl;
          cout << "Filename for initial data output:                     -of [\"" << iFile << "\"]" << endl;
          cout << "Title of initial data file:                           -ot [\"" << iTitle << "\"]" << endl;
          cout << endl;
          cout << "Operating mode (0 = cells, 1 = pulses, 2 = periodic): -om [" << mode << "]" << endl;
          cout << "Number of lattice sites (2 to " << MAX_NSIDE << "):                  -ns [" << nSide << "]" << endl;
          cout << "Number of control points:                             -nc [" << nCells << "]" << endl;
          cout << "Boundary width (lattice units):                       -bw [" << boundaryWidth << "]" << endl;
          cout << "Initial A factor:                                     -af [" << aFactorFactor << "]" << endl;
          cout << "Periods along x axis (mode 2):                        -nx [" << nx << "]" << endl;

          exit(1);
       }

       i += 2;
    }

    cout << "nSide: " << nSide << endl;

    if ((2 > nSide) || (nSide > MAX_NSIDE))
    {
       cout << "Number of lattice sites must be in range 2 to " << MAX_NSIDE << endl;
       exit(2);
    }

    // Initialize random generator

    cout << "Initializing random number generator." << endl << flush;

#ifdef USE_MKL
    vslNewStream(&pVslStream, VSL_BRNG_WH, rs);
#else
    srand(rs);
#endif

    // Initialize data

    cout << "Generating initial data." << endl << flush;

    switch (mode)
    {
       case 1:  pulses(nSide, nCells, boundaryWidth, aFactorFactor);
                break;
       case 2:  periodic(nSide, nx, aFactorFactor);
                break;
       default: voronoi(nSide, nCells, boundaryWidth, aFactorFactor);
    }

    cout << "Saving initial data to \"" << iFile << ".vtk\"." << endl << flush;

    write_vtk(iFile + ".vtk", iTitle, nSide, ONE, ZERO);

    cout << "Done saving." << endl << flush;
    
    return(0);
}
