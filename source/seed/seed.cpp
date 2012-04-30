/* seed.cpp : Creates random field configurations.
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
const int MAX_NSIDE                 = 64;            /* Number of lattice sites per side.
                                                        Power of 2 for faster indexing (maybe). */
#endif

const int DEFAULT_NSIDE             = 64;

#ifndef EPSILON
const double EPSILON                = 1.0e-10;
#endif

const double DEFAULT_BOUNDARY_WIDTH = 5;             /* Default initial width (in lattice units) of cell boundaries */
const double MAX_INITIAL_A          = 10;            /* Max initial |A| */

const int OUTPUT_PRECISION          = 15;

const int DEFAULT_NCELLS            = 10;            /* Default initial number of Voronoi cells in grid */
const int NCOMPS                    = 7;             /* Number of field components used INTERNALLY */

typedef double TCube[NCOMPS][MAX_NSIDE][MAX_NSIDE][MAX_NSIDE];   /* First index = field (n_0..n_3, A_1..A_3).
                                                                                                                                                   Last 3 indices: x, y, z position of cell */
TCube nA;

#ifdef USE_MKL
VSLStreamStatePtr pVslStream = NULL;
#endif

struct TRndm3
{
   double v[3];
   
   TRndm3() { v[0] = 0; v[1] = 0; v[2] = 0; }

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
   
   TRndm4() { v[0] = 0; v[1] = 0; v[2] = 0; v[3] = 0; }

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
   double i, j, k;
};

void voronoi(int nSide, int nCells, double boundaryWidth, double aFactorFactor)
/*
   boundaryWidth is expressed as a multiple of lattice spacing, not in physical units
*/
{
   // Create list of random cell centers with random n_a and magnetic potential

   TRndm3 rndm3;
   TRndm4 rndm4;

   double aFactor = TWO * MAX_INITIAL_A * aFactorFactor;

   TCell* cell = new TCell[nCells];

   // Do this outside parallel section so same random seed always yields same initial configuration

   for (int i = 0; i < nCells; i++)
   {
      rndm3.octant(nSide);

      cell[i].i = rndm3.v[0];
      cell[i].j = rndm3.v[1];
      cell[i].k = rndm3.v[2];

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
   // Loop through grid and set its field values as spline interpolation between two closest cell centers

#pragma omp for
#endif
   for (int i = 0; i < nSide; i++)
       for (int j = 0; j < nSide; j++)
           for (int k = 0; k < nSide; k++)
           {
              // Determine two closest cell centers

              double mrR2 = TWO*nSide*nSide; int m =-1; // Closest
              double nrR2 = TWO*nSide*nSide; int n =-1; // Next closest

              for (int l = 0; l < nCells; l++)
              {
                 rndm3.v[0] = fabs(cell[l].i - i); rndm3.v[0] = min(rndm3.v[0], nSide - rndm3.v[0]);
                 rndm3.v[1] = fabs(cell[l].j - j); rndm3.v[1] = min(rndm3.v[1], nSide - rndm3.v[1]);
                 rndm3.v[2] = fabs(cell[l].k - k); rndm3.v[2] = min(rndm3.v[2], nSide - rndm3.v[2]);

                 double rR2 = rndm3.sqr();

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

              // Project vector from nearest center to current position onto vector
              // from nearest center to next nearest center for distance along spline

              rndm3.v[0] = fabs(cell[n].i - cell[m].i); rndm3.v[0] = min(rndm3.v[0], nSide - rndm3.v[0]);
              rndm3.v[1] = fabs(cell[n].j - cell[m].j); rndm3.v[1] = min(rndm3.v[1], nSide - rndm3.v[1]);
              rndm3.v[2] = fabs(cell[n].k - cell[m].k); rndm3.v[2] = min(rndm3.v[2], nSide - rndm3.v[2]);

              double r = rndm3.norm();
              rndm3.div(r);

              double v1 = fabs(i - cell[m].i); v1 = min(v1, nSide - v1);
              double v2 = fabs(j - cell[m].j); v2 = min(v2, nSide - v2);
              double v3 = fabs(k - cell[m].k); v3 = min(v3, nSide - v3);

              // Projection of current position vector (from closest center) on line between centers

              v1 *= rndm3.v[0];
              v2 *= rndm3.v[1];
              v3 *= rndm3.v[2];

              // Distance parameter for spline is computed from midpoint of line between centers

              r = spline(sqrt(v1*v1 + v2*v2 + v3*v3) - HALF*r, boundaryWidth);

              // Interpolate n_a

              rndm4.v[0] = (cell[m].y[0] - cell[n].y[0])*r + cell[n].y[0];
              rndm4.v[1] = (cell[m].y[1] - cell[n].y[1])*r + cell[n].y[1];
              rndm4.v[2] = (cell[m].y[2] - cell[n].y[2])*r + cell[n].y[2];
              rndm4.v[3] = (cell[m].y[3] - cell[n].y[3])*r + cell[n].y[3];

              rndm4.normalize();

              nA[0][i][j][k] = rndm4.v[0];
              nA[1][i][j][k] = rndm4.v[1];
              nA[2][i][j][k] = rndm4.v[2];
              nA[3][i][j][k] = rndm4.v[3];

              // Interpolate A

              nA[4][i][j][k] = (cell[m].y[4] - cell[n].y[4])*r + cell[n].y[4];
              nA[5][i][j][k] = (cell[m].y[5] - cell[n].y[5])*r + cell[n].y[5];
              nA[6][i][j][k] = (cell[m].y[6] - cell[n].y[6])*r + cell[n].y[6];
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

   TRndm3 rndm3;
   TRndm4 rndm4;

   double aFactor = TWO * MAX_INITIAL_A * aFactorFactor;

   TCell* cell = new TCell[nCells];

   // Do this outside parallel section so same random seed always yields same initial configuration

   for (int i = 0; i < nCells; i++)
   {
      rndm3.octant(nSide);

      if (1 == nCells)
      {
         cell[i].i = nSide / 2;
         cell[i].j = nSide / 2;
         cell[i].k = nSide / 2;
      }
      else
      {
         cell[i].i = rndm3.v[0];
         cell[i].j = rndm3.v[1];
         cell[i].k = rndm3.v[2];
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
       for (int j = 0; j < nSide; j++)
           for (int k = 0; k < nSide; k++)
           {
              // Determine closest cell center

              double mrR2 = TWO*nSide*nSide;
              double rR2;

              int m =-1; // Closest

              for (int l = 0; l < nCells; l++)
              {
                 rndm3.v[0] = fabs(cell[l].i - i); rndm3.v[0] = min(rndm3.v[0], nSide - rndm3.v[0]);
                 rndm3.v[1] = fabs(cell[l].j - j); rndm3.v[1] = min(rndm3.v[1], nSide - rndm3.v[1]);
                 rndm3.v[2] = fabs(cell[l].k - k); rndm3.v[2] = min(rndm3.v[2], nSide - rndm3.v[2]);

                 rR2 = rndm3.sqr();

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

              nA[0][i][j][k] = rndm4.v[0];
              nA[1][i][j][k] = rndm4.v[1];
              nA[2][i][j][k] = rndm4.v[2];
              nA[3][i][j][k] = rndm4.v[3];

              // Interpolate A

              nA[4][i][j][k] = cell[m].y[4] * r;
              nA[5][i][j][k] = cell[m].y[5] * r;
              nA[6][i][j][k] = cell[m].y[6] * r;
           }

#ifdef _OPENMP
#pragma omp barrier
}
#endif

   delete [] cell;
}

void tubes(int nSide, int nTubes, int nPoints, double boundaryWidth, double xRadius, double yRadius, 
           int nzTwists, int npTwists, double aFactorFactor)
{
#ifdef _OPENMP
#pragma omp parallel default(shared)
{
   // Loop through grid, make sure all cells are initialized to 0.

#pragma omp for
#endif
   for (int i = 0; i < nSide; i++)
       for (int j = 0; j < nSide; j++)
           for (int k = 0; k < nSide; k++)
           {
              nA[0][i][j][k] = 0;
              nA[1][i][j][k] = 0;
              nA[2][i][j][k] = 0;
              nA[3][i][j][k] = 0;

              nA[4][i][j][k] = 0;
              nA[5][i][j][k] = 0;
              nA[6][i][j][k] = 0;
           }

#ifdef _OPENMP
#pragma omp barrier
}
#endif

   TRndm3 rndm3;
   TRndm4 rndm4;

   double aFactor = TWO * MAX_INITIAL_A * aFactorFactor;
   double z_phase_step = nzTwists * TWO * PI / nPoints;

   TCell* point = new TCell[nPoints];
   
   // Since there is no guarantee that paths will not overlap, this is not parallellized; otherwise, multiple threads might end
   // up writing to the same location at the same time.

   for (int tubeNr = 0; tubeNr < nTubes; tubeNr++)
   {
      // Create list of random control point coordinates. Do this outside parallel section so same random seed always yields same initial configuration.

      for (int i = 0; i < nPoints; i++)
      {
         rndm3.octant(nSide);

         point[i].i = rndm3.v[0];
         point[i].j = rndm3.v[1];
         point[i].k = rndm3.v[2];
      }
      
      // Start point. Store theta values rather than n_a for better interpolation.

      rndm4.sphere(ONE, TWO, ONE);

      double thetaHalf = acos(rndm4.v[0]);
      double thetaOverSinThetaHalf = TWO;

      if (fabs(thetaHalf) > EPSILON) thetaOverSinThetaHalf *= thetaHalf / sin(thetaHalf);

      point[0].y[0] = thetaOverSinThetaHalf;
      point[0].y[1] = thetaOverSinThetaHalf * rndm4.v[1];
      point[0].y[2] = thetaOverSinThetaHalf * rndm4.v[2];
      point[0].y[3] = thetaOverSinThetaHalf * rndm4.v[3];      
      
      rndm3.ball(HALF, ONE, QUARTER);
      rndm3.mul(aFactor);

      point[0].y[4] = rndm3.v[0];
      point[0].y[5] = rndm3.v[1];
      point[0].y[6] = rndm3.v[2];
      
      // Step along Catmull-Rom splines between control points, generate field values: A along path, theta on ellipses about it.
      
      int h = nPoints - 1;
   
      for (int i = 0; i < nPoints; i++)
      {
         int j = (i + 1) % nPoints;
         int k = (i + 2) % nPoints;

         double mi_x = point[j].i - point[h].i;
         double mi_y = point[j].j - point[h].j;
         double mi_z = point[j].k - point[h].k;

         double mj_x = point[k].i - point[i].i;
         double mj_y = point[k].j - point[i].j;
         double mj_z = point[k].k - point[i].k;

         double z_phase = i * z_phase_step;
         
         // Step along spline between the two control points point[i] and point[j].

         for (int ii = 0; ii < 2 * nSide; ii++)
         {
            double t = ii / (TWO * nSide);
            double local_z_phase = z_phase + t * z_phase_step;

            double di = point[j].i - point[i].i;
            double dj = point[j].j - point[i].j;
            double dk = point[j].k - point[i].k;

            double pos_x = point[i].i + t*(mi_x + t*((3*di - mj_x - 2*mi_x) + t*(mi_x + mj_x - 2*di)));
            double pos_y = point[i].j + t*(mi_y + t*((3*dj - mj_y - 2*mi_y) + t*(mi_y + mj_y - 2*dj)));
            double pos_z = point[i].k + t*(mi_z + t*((3*dk - mj_z - 2*mi_z) + t*(mi_z + mj_z - 2*dk)));

            double z_axis_x = mi_x + t*(2*(- 2*mi_x - mj_x + 3*di) - t*3*(- mi_x - mj_x + 2*di));
            double z_axis_y = mi_y + t*(2*(- 2*mi_y - mj_y + 3*dj) - t*3*(- mi_y - mj_y + 2*dj));
            double z_axis_z = mi_z + t*(2*(- 2*mi_z - mj_z + 3*dk) - t*3*(- mi_z - mj_z + 2*dk));

            double norm = sqrt(z_axis_x*z_axis_x + z_axis_y*z_axis_y + z_axis_z*z_axis_z);

            z_axis_x /= norm;
            z_axis_y /= norm;
            z_axis_z /= norm;

            int ix = int(pos_x); while (0 > ix) ix += nSide; while (ix >= nSide) ix -= nSide;
            int iy = int(pos_y); while (0 > iy) iy += nSide; while (iy >= nSide) iy -= nSide;
            int iz = int(pos_z); while (0 > iz) iz += nSide; while (iz >= nSide) iz -= nSide;

            // Create ellipse centered on current position and orthogonal to direction of spline.
            
            double s = sqrt(1 - z_axis_z*z_axis_z);
            double tt = 1 - z_axis_z;
   
            // Axis of rotation
   
            double nrm = sqrt(z_axis_x * z_axis_x + z_axis_y * z_axis_y);
            double x = z_axis_y / nrm;
            double y = -z_axis_x / nrm;
   
            // Rotation matrix.
   
            double R00 = tt*x*x + z_axis_z;
            double R02 = y*s;       
            double R10 = tt*x*y;     
            double R11 = tt*y*y + z_axis_z;
            double R12 = -x*s;

            // Write rotated value of A vector to grid site at center of tube.

            double cos_local_z_phase = cos(local_z_phase);
            double sin_local_z_phase = sin(local_z_phase);

            // First rotate A about z axis in its own frame by local z phase...
 
            double xp = point[0].y[4] * cos_local_z_phase - point[0].y[5] * sin_local_z_phase;
            double yp = point[0].y[4] * sin_local_z_phase + point[0].y[5] * cos_local_z_phase;

            // ...then rotate local frame so its z axis coincides with z_axis (direction of spline).

            nA[4][ix][iy][iz] = R00 * xp + R10 * yp - R02 * point[0].y[6];
            nA[5][ix][iy][iz] = R10 * xp + R11 * yp - R12 * point[0].y[6];
            nA[6][ix][iy][iz] = R02 * xp + R12 * yp + z_axis_z * point[0].y[6];
            
            nA[0][ix][iy][iz] = 1; // Mark this grid site as set
            
            // Create ellipse in (x,y) plane, centered on (0,0); rotate it to plane with normal z_axis; offset to position.
   
            for (int iii = 0; iii < 4 * nSide; iii++)
            {
               double theta = iii * FOUR * PI / (4 * nSide);

               xp = xRadius * cos(theta);
               yp = yRadius * sin(theta);

               ix = int(R00 * xp + R10 * yp + pos_x); while (0 > ix) ix += nSide; while (ix >= nSide) ix -= nSide;
               iy = int(R10 * xp + R11 * yp + pos_y); while (0 > iy) iy += nSide; while (iy >= nSide) iy -= nSide;
               iz = int(R02 * xp + R12 * yp + pos_z); while (0 > iz) iz += nSide; while (iz >= nSide) iz -= nSide;
      
               double phase = local_z_phase + npTwists * theta;
      
               nA[1][ix][iy][iz] = wrap(point[0].y[1] + phase);
               nA[2][ix][iy][iz] = wrap(point[0].y[2] + phase);
               nA[3][ix][iy][iz] = wrap(point[0].y[3] + phase);
               
               nA[0][ix][iy][iz] = 1; // Mark this grid site as set
            }
         }
     
         h = i;
      }
   
      // Scan the grid, count sites which have been set.
   
      int nCells = 0;
   
#ifdef _OPENMP
#pragma omp parallel reduction(+ : nCells)
{
#pragma omp for
#endif
      for (int i = 0; i < nSide; i++)
          for (int j = 0; j < nSide; j++)
              for (int k = 0; k < nSide; k++)
              {
                 if (nA[0][i][j][k]) nCells++;
              }
#ifdef _OPENMP
#pragma omp barrier
}
#endif

      // Assemble list of sites which have been set. Those become our cell centers for interpolation.

      TCell* cell = new TCell[nCells];
   
      int idx = 0;

      for (int i = 0; i < nSide; i++)
          for (int j = 0; j < nSide; j++)
              for (int k = 0; k < nSide; k++)
              {
                 if (nA[0][i][j][k])
                 {
                    // Turn theta back into n_a

                    double thR2 = sqr(nA[1][i][j][k]) + sqr(nA[2][i][j][k]) + sqr(nA[3][i][j][k]);

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

                    cell[idx].y[0] = rndm4.v[0];
                    cell[idx].y[1] = rndm4.v[1];
                    cell[idx].y[2] = rndm4.v[2];
                    cell[idx].y[3] = rndm4.v[3];

                    cell[idx].y[4] = nA[4][i][j][k];
                    cell[idx].y[5] = nA[5][i][j][k];
                    cell[idx].y[6] = nA[6][i][j][k];

                    cell[idx].i = i;
                    cell[idx].j = j;
                    cell[idx].k = k;

                    idx++;
                 }
              }

#ifdef _OPENMP
#pragma omp parallel default(shared) private(rndm3, rndm4)
{
      // Loop through grid and set its field values as spline interpolation between closest cell center and 0.

#pragma omp for
#endif
      for (int i = 0; i < nSide; i++)
          for (int j = 0; j < nSide; j++)
              for (int k = 0; k < nSide; k++)
              {
                 // Determine closest cell center

                 double mrR2 = TWO*nSide*nSide;
                 double rR2;

                 int m =-1; // Closest

                 for (int l = 0; l < nCells; l++)
                 {
                    rndm3.v[0] = fabs(cell[l].i - i); rndm3.v[0] = min(rndm3.v[0], nSide - rndm3.v[0]);
                    rndm3.v[1] = fabs(cell[l].j - j); rndm3.v[1] = min(rndm3.v[1], nSide - rndm3.v[1]);
                    rndm3.v[2] = fabs(cell[l].k - k); rndm3.v[2] = min(rndm3.v[2], nSide - rndm3.v[2]);

                    rR2 = rndm3.sqr();

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

                 nA[0][i][j][k] = rndm4.v[0];
                 nA[1][i][j][k] = rndm4.v[1];
                 nA[2][i][j][k] = rndm4.v[2];
                 nA[3][i][j][k] = rndm4.v[3];

                 // Interpolate A

                 nA[4][i][j][k] = cell[m].y[4] * r;
                 nA[5][i][j][k] = cell[m].y[5] * r;
                 nA[6][i][j][k] = cell[m].y[6] * r;
              }

#ifdef _OPENMP
#pragma omp barrier
}
#endif

      delete [] cell;
   }
   
   delete [] point;
}

void periodic(int nSide, int nx, int ny, int nz, double aFactorFactor)
{
   // Create 4 cells (central + one for each axis) with random n_a and magnetic potential

   TRndm3 rndm3;
   TRndm4 rndm4;

   double aFactor = TWO * MAX_INITIAL_A * aFactorFactor;

   TCell* cell = new TCell[4];

   // Do this outside parallel section so same random seed always yields same initial configuration

   for (int i = 0; i < 4; i++)
   {
      cell[i].i = nSide / 2;
      cell[i].j = nSide / 2;
      cell[i].k = nSide / 2;

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
   double ky = (2.0 * ny) / nSide;
   double kz = (2.0 * nz) / nSide;
   
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
       
       for (int j = 0; j < nSide; j++)
       {
           double r_jj;
           double r_j = (j - nSide / 2) * ky;
           r_j = fabs(modf(r_j, &r_jj));
       
           for (int k = 0; k < nSide; k++)
           {
              double r_kk;
              double r_k = (k - nSide / 2) * kz;
              r_k = fabs(modf(r_k, &r_kk));
            
              // Interpolate theta

              rndm4.v[1] = wrap(((3 - r_i - r_j - r_k) * cell[0].y[1] + r_i * cell[1].y[1] + r_j * cell[2].y[1] + r_k * cell[3].y[1]) / 3);
              rndm4.v[2] = wrap(((3 - r_i - r_j - r_k) * cell[0].y[2] + r_i * cell[1].y[2] + r_j * cell[2].y[2] + r_k * cell[3].y[2]) / 3);
              rndm4.v[3] = wrap(((3 - r_i - r_j - r_k) * cell[0].y[3] + r_i * cell[1].y[3] + r_j * cell[2].y[3] + r_k * cell[3].y[3]) / 3);

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

              nA[0][i][j][k] = rndm4.v[0];
              nA[1][i][j][k] = rndm4.v[1];
              nA[2][i][j][k] = rndm4.v[2];
              nA[3][i][j][k] = rndm4.v[3];

              // Interpolate A

              nA[4][i][j][k] = ((3 - r_i - r_j - r_k) * cell[0].y[4] + r_i * cell[1].y[4] + r_j * cell[2].y[4] + r_k * cell[3].y[4]) / 3;
              nA[5][i][j][k] = ((3 - r_i - r_j - r_k) * cell[0].y[5] + r_i * cell[1].y[5] + r_j * cell[2].y[5] + r_k * cell[3].y[5]) / 3;
              nA[6][i][j][k] = ((3 - r_i - r_j - r_k) * cell[0].y[6] + r_i * cell[1].y[6] + r_j * cell[2].y[6] + r_k * cell[3].y[6]) / 3;
           }
       }
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
   vtk << "DIMENSIONS " << nSide << ' ' << nSide << ' ' << nSide << endl;
   vtk << "ORIGIN 0 0 0" << endl;
   vtk << "SPACING " << dxn << ' ' << dxn << ' ' << dxn << endl;
   vtk << "POINT_DATA " << nSide*nSide*nSide << endl;

   vtk << "VECTORS T DOUBLE" << endl;

   for (int i = 0; i < nSide; i++)
       for (int j = 0; j < nSide; j++)
           for (int k = 0; k < nSide; k++)
           {
              double thetaHalf = acos(nA[0][i][j][k]);
              double thetaOverSinThetaHalf = TWO;

              if (fabs(thetaHalf) > EPSILON) thetaOverSinThetaHalf *= thetaHalf / sin(thetaHalf);

              vtk << thetaOverSinThetaHalf * nA[1][i][j][k] << ' ' 
                  << thetaOverSinThetaHalf * nA[2][i][j][k] << ' ' 
                  << thetaOverSinThetaHalf * nA[3][i][j][k] << endl;
           }

   vtk << "VECTORS A DOUBLE" << endl;

   for (int i = 0; i < nSide; i++)
       for (int j = 0; j < nSide; j++)
           for (int k = 0; k < nSide; k++)
           {
              vtk << nA[4][i][j][k] << ' ' 
                  << nA[5][i][j][k] << ' ' 
                  << nA[6][i][j][k] << endl;
           }

   vtk << "VECTORS pT DOUBLE" << endl;

   for (int i = 0; i < nSide; i++)
       for (int j = 0; j < nSide; j++)
           for (int k = 0; k < nSide; k++)
           {
              vtk << 0 << ' ' << 0 << ' ' << 0 << endl;
           }

   vtk << "VECTORS pA DOUBLE" << endl;

   for (int i = 0; i < nSide; i++)
       for (int j = 0; j < nSide; j++)
           for (int k = 0; k < nSide; k++)
           {
              vtk << 0 << ' ' << 0 << ' ' << 0 << endl;
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
    int ny               = 1;
    int nz               = 1;
    int np               = 1;
    int nt               = 1;
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
       if (!strcmp("-nc", argv[i])) nCells = atol(argv[i + 1]); else
       if (!strcmp("-bw", argv[i])) boundaryWidth = atof(argv[i + 1]); else
       if (!strcmp("-af", argv[i])) aFactorFactor = atof(argv[i + 1]); else
       if (!strcmp("-om", argv[i])) mode = atol(argv[i + 1]); else
       if (!strcmp("-nx", argv[i])) nx = atol(argv[i + 1]); else
       if (!strcmp("-ny", argv[i])) ny = atol(argv[i + 1]); else
       if (!strcmp("-nz", argv[i])) nz = atol(argv[i + 1]); else
       if (!strcmp("-np", argv[i])) np = atol(argv[i + 1]); else
       if (!strcmp("-nt", argv[i])) nt = atol(argv[i + 1]); else
       {
          cout << "Unrecognized command line argument: " << argv[i] << endl;
          cout << endl;
          cout << "Supported options [defaults in square brackets]:" << endl;
          cout << endl;
          cout << "Random seed (long integer):                                             -rs [" << rs << "]" << endl;
          cout << endl;
          cout << "Filename for initial data output:                                       -of [\"" << iFile << "\"]" << endl;
          cout << "Title of initial data file:                                             -ot [\"" << iTitle << "\"]" << endl;
          cout << endl;
          cout << "Lattice sites per side (2 to " << MAX_NSIDE << "):                                       -ns [" << nSide << "]" << endl;
          cout << "Operating mode (0 = cells, 1 = pulses, 2 = periodic, 3 = tubes):        -om [" << mode << "]" << endl;
          cout << "Number of control points:                                               -nc [" << nCells << "]" << endl;
          cout << "Number of tubes (mode 3):                                               -nt [" << nt << "]" << endl;
          cout << "Boundary width (lattice units):                                         -bw [" << boundaryWidth << "]" << endl;
          cout << "Initial A factor:                                                       -af [" << aFactorFactor << "]" << endl;
          cout << "Periods along x axis (mode 2) or tube x radius (mode 3, lattice units): -nx [" << nx << "]" << endl;
          cout << "Periods along y axis (mode 2) or tube y radius (mode 3, lattice units): -ny [" << ny << "]" << endl;
          cout << "Periods along z direction (modes 2 and 3):                              -nz [" << nz << "]" << endl;
          cout << "Periods along poloidal direction (mode 3):                              -np [" << np << "]" << endl;

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
    
    if (3 == mode)
    {
       if (3 > nCells)
       {
          cout << "Number of control points must be > 2 for tubes" << endl;
          exit(3);
       }
       
       if (1 > nt)
       {
          cout << "Number of tubes must be > 0" << endl;
          exit(4); 
       }
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
       case 2:  periodic(nSide, nx, ny, nz, aFactorFactor);
                break;
       case 3:  tubes(nSide, nt, nCells, boundaryWidth, nx, ny, nz, np, aFactorFactor);
                break;
       default: voronoi(nSide, nCells, boundaryWidth, aFactorFactor);
    }

    cout << "Saving initial data to \"" << iFile << ".vtk\"." << endl << flush;

    write_vtk(iFile + ".vtk", iTitle, nSide, ONE, ZERO);

    cout << "Done saving." << endl << flush;
    
    return(0);
}
