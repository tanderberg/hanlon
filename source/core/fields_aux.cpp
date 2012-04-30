/* fields_aux.cpp : Auxiliary functionality for scalars in polar form.
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

#ifndef MAX_NSIDE
const int MAX_NSIDE              = 128;   /* Max number of cells per side of the lattice cube.
                                             Power of 2 for faster array indexing (maybe). */
#endif

#define WRAP_DTH                 = TRUE;

const int DEFAULT_NSIDE          = 64;    /* Default number of cells per side of the lattice cube */
const double SMALL_THETA_SQUARED = 1e-14; /* Below which series expansion is used to avoid 1/0 */

#include "lu.cpp"
#include "math_aux.cpp"

struct TFields
{
   double th_1, th_2, th_3,
          A_1, A_2, A_3;

   void clear();
   void add_fields(TFields *otherFields);
   void sub_fields(TFields *otherFields);
   void mul_fields(double factor);
};

void TFields::clear()
{
   th_1 = 0.0;
   th_2 = 0.0;
   th_3 = 0.0;

   A_1 = 0.0;
   A_2 = 0.0;
   A_3 = 0.0;
}

void TFields::add_fields(TFields *otherFields)
{
   th_1 += otherFields->th_1;
   th_2 += otherFields->th_2;
   th_3 += otherFields->th_3;

   A_1 += otherFields->A_1;
   A_2 += otherFields->A_2;
   A_3 += otherFields->A_3;
}

void TFields::sub_fields(TFields *otherFields)
{
   th_1 -= otherFields->th_1;
   th_2 -= otherFields->th_2;
   th_3 -= otherFields->th_3;

   A_1 -= otherFields->A_1;
   A_2 -= otherFields->A_2;
   A_3 -= otherFields->A_3;
}

void TFields::mul_fields(double factor)
{
   th_1 *= factor;
   th_2 *= factor;
   th_3 *= factor;

   A_1 *= factor;
   A_2 *= factor;
   A_3 *= factor;
}

struct TRho
{
   double BR2, rH, rG;
};

class TAux
{
protected:

   double h_00, h_01, h_02, h_11, h_12, h_22, // gR2 * H_ab
          G_00, G_01, G_02, G_11, G_12, G_22,
          
          dx_th_1, dx_th_2, dx_th_3, // IMPORTANT! Derivatives stored here are not normalized by lattice spacing
          dy_th_1, dy_th_2, dy_th_3,
          dz_th_1, dz_th_2, dz_th_3,

          dx_A_2, dx_A_3,
          dy_A_1, dy_A_3,
          dz_A_1, dz_A_2;

   double get_rho_static_core(TFields *fields, TRho &rho, 
                              double h_00, double h_01, double h_02, double h_11, double h_12, double h_22, 
                              double G_00, double G_01, double G_02, double G_11, double G_12, double G_22);
                               
public: 

   void compute_hG_core(double th_1, double th_2, double th_3, 
                        double &h_00, double &h_01, double &h_02, double &h_11, double &h_12, double &h_22, 
                        double &G_00, double &G_01, double &G_02, double &G_11, double &G_12, double &G_22);
   void compute_hG(double th_1, double th_2, double th_3);
   double get_derivative(int axis, int component);
   void set_derivative(int axis, int component, double value);
   void compute_derivative(TFields fields[9][9][9], int i, int j, int k, int axis, int component);
   void compute_derivative(TFields fields[MAX_NSIDE][MAX_NSIDE][MAX_NSIDE], int i, int j, int k, int nSideM1, int axis, int component);
   void compute_derivatives(TFields fields[MAX_NSIDE][MAX_NSIDE][MAX_NSIDE], int i, int j, int k, int nSideM1);
   double get_rho_static(TFields *fields, TRho &rho);
   void subtract_drho_static_dA(TFields *fields, double &drdA1, double &drdA2, double &drdA3);
   void subtract_offset_drho_dA(int di, int dj, int dk, double &drdA1, double &drdA2, double &drdA3);
   void subtract_drho_static_dth(TFields *fields, double &drdth1, double &drdth2, double &drdth3);
   void subtract_offset_drho_dth(int di, int dj, int dk, TFields *fields, double &drdth1, double &drdth2, double &drdth3);
};

void TAux::compute_hG_core(double th_1, double th_2, double th_3, 
                           double &h_00, double &h_01, double &h_02, double &h_11, double &h_12, double &h_22, 
                           double &G_00, double &G_01, double &G_02, double &G_11, double &G_12, double &G_22)
{
   double th1R2 = th_1*th_1;
   double th2R2 = th_2*th_2;
   double th3R2 = th_3*th_3;

   double thTR2 = th1R2 + th2R2;
   double thR2 = thTR2 + th3R2;

   if (thR2 < SMALL_THETA_SQUARED) // O(th_i^3) series expansion
   {
      h_00 = 1-(2*th_1*th_2*th_3 + 4*th2R2 + th3R2)/12;
      h_01 = (4*th_1*th_2 + th_3*(th2R2 - th1R2))/12;
      h_02 = -(2*th_1*th_3 + th_2*(3*thTR2 + th3R2 - 12))/24;
      h_11 = 1-(-2*th_1*th_2*th_3 + 4*th1R2 + th3R2)/12;
      h_12 = (-2*th_2*th_3 + th_1*(3*thTR2 + th3R2 - 12))/24;
      h_22 = thTR2/4;

      G_00 = 1.0/4.0 - (th_2*th_2 + th_3*th_3)/48.0;
      G_01 = th_2*th_1/48.0;
      G_02 = th_3*th_1/48.0;
      G_11 = 1.0/4.0 - (th_1*th_1 + th_3*th_3)/48.0;
      G_12 = th_3*th_2/48.0;
      G_22 = 1.0/4.0 - thTR2/48.0;
   }   
   else
   {
      double th = sqrt(thR2);
      double thR3 = thR2*th;
      double thR4 = thR2*thR2;
      double thR5 = thR4*th;
      double thR6 = thR4*thR2;
      double thR8 = thR6*thR2;
      
      double th2R3 = th2R2 * th_2;
      double th2R4 = th2R2 * th2R2;
      
      double th3R4 = th3R2 * th3R2;

      double mess = th2R4 + th1R2*(th1R2 + th2R2*2.0);
      
      double th12 = th_1*th_2;
      double th13 = th_1*th_3;
      double th23 = th_2*th_3;
      
      double th123 = th12*th_3;

      double cth = cos(th);
      double sth = sin(th);
      
      double cthR2 = cth*cth;
      double sthR2 = sth*sth;
      
      double th_sth = th * sth;
      double cth2Pth_sthR2 = pow(cth*2.0 + th_sth, 2.0);
      
      double cthM1 = cth - 1.0;
      double cthM1R2 = cthM1*cthM1;

      double hF = gR2 / thR8; /* Common factor */

      h_00 = hF*(th3R2*(4.0*th1R2*th3R2 + thR2*(thR2*cthM1R2 - 2.0*th1R2*(cthR2*2.0 + th_sth*(cthM1) + 2.0))) + th1R2*(cth2Pth_sthR2*mess + th3R2*thTR2*(4.0*cth*(cth + th_sth) + thR2*sthR2 + 4.0)) + thR2*(th*(thR3*sthR2 + 2.0*th123*(th*cth*(cthM1) + sth*(cth + th_sth - 1.0))) + th1R2*(-2.0*thR2*(sthR2 - cth*(cth*2.0 + th_sth)) + thTR2*(sthR2 + cth*(-6.0*th_sth + cth*(thR2 - 8.0))))));

      h_01 = hF*(th_2*(th_1*cth2Pth_sthR2*mess + th_3*(thR3*th_2*(th*cth*(cthM1) + sth*(cth + th_sth - 1.0)) + th13*thTR2*(4.0*cth*(cth + th_sth) + thR2*sthR2 + 4.0))) + th_1*(-(thR3*th13*(th*cth*(cthM1) + sth*(cth + th_sth - 1.0))) + th_2*(2.0*th3R2*(th3R2*2.0 - thR2*(cthR2*2.0 + th_sth*(cthM1) + 2.0)) + thR2*(-2.0*thR2*(sthR2 - cth*(cth*2.0 + th_sth)) + thTR2*(sthR2 + cth*(-6.0*th_sth + cth*(thR2 - 8.0)))))));

      h_02 = -hF*(-4.0*th13*th3R4 - th13*(th3R2*thTR2*(thR2 + cth*(4.0*th_sth + cth*(-thR2 + 4.0)) + 4.0) + mess*(thR2 + cth*(4.0*th_sth + cth*(-thR2 + 4.0)))) + thR2*(th*(-(th*th13*cth*(cth*2.0 + th_sth + 2.0)) + th_2*(cthM1)*(thR2*sth - th3R2*(-th + sth))) + th13*(th3R2*(cthR2*2.0 + th_sth*(cthM1) + 6.0) + thTR2*(th_sth*(cth*5.0 + 1.0) + cthR2*(-thR2 + 7.0) + 1.0))));

      h_11 = hF*(th3R2*(4.0*th2R2*th3R2 + thR2*(thR2*cthM1R2 - 2.0*th2R2*(cthR2*2.0 + th_sth*(cthM1) + 2.0))) + th2R2*(cth2Pth_sthR2*mess + th3R2*thTR2*(4.0*cth*(cth + th_sth) + thR2*sthR2 + 4.0)) + thR2*(th*(thR3*sthR2 - 2.0*th123*(th*cth*(cthM1) + sth*(cth + th_sth - 1.0))) + th2R2*(-2.0*thR2*(sthR2 - cth*(cth*2.0 + th_sth)) + thTR2*(sthR2 + cth*(-6.0*th_sth + cth*(thR2 - 8.0))))));

      h_12 = -hF*(-4.0*th23*th3R4 - thR5*th_1*sth*(cthM1) + th_3*(-(th_2*(th3R2*thTR2*(thR2 + cth*(4.0*th_sth + cth*(-thR2 + 4.0)) + 4.0) + mess*(thR2 + cth*(4.0*th_sth + cth*(-thR2 + 4.0))))) + thR2*(th*(th13*(cthM1)*(-th + sth) - th*th_2*cth*(cth*2.0 + th_sth + 2.0)) + th2R3*(th_sth*(cth*5.0 + 1.0) + cthR2*(-thR2 + 7.0) + 1.0) + th_2*(th3R2*(cthR2*2.0 + th_sth*(cthM1) + 6.0) + th1R2*(th_sth*(cth*5.0 + 1.0) + cthR2*(-thR2 + 7.0) + 1.0)))));

      h_22 = hF*(4.0*th3R2*(thR4 + th3R2*(thR2*-2.0 + th3R2)) + th3R2*(cth2Pth_sthR2*mess + th3R2*thTR2*(4.0*cth*(cth + th_sth) + thR2*sthR2 + 4.0)) + thR2*(thR2*cthM1R2*thTR2 + th3R2*thTR2*(cthR2*-4.0 + sthR2 + th*(th*cthR2 - 2.0*sth*(cth*2.0 + 1.0)) - 4.0)));
   
      double GF1 = 0.25 / thR2;
      double GF2 = -cthM1 / (2.0 * thR4);
      double GFD = GF1 - GF2;

      G_00 = th1R2 * GF1 + (th2R2 + th3R2) * GF2;
      G_01 = th12 * GFD;
      G_02 = th13 * GFD;

      G_11 = th2R2 * GF1 + (th1R2 + th3R2) * GF2;
      G_12 = th23 * GFD;

      G_22 = th3R2 * GF1 + thTR2 * GF2;
   }
}

void TAux::compute_hG(double th_1, double th_2, double th_3)
{
   compute_hG_core(th_1, th_2, th_3, 
                   h_00, h_01, h_02, h_11, h_12, h_22, 
                   G_00, G_01, G_02, G_11, G_12, G_22);
}

double TAux::get_derivative(int axis, int component)
{
   switch (axis)
   {
      case 0: // x derivative
         switch (component)
         {
             case 0: return(dx_th_1);
             case 1: return(dx_th_2);
             case 2: return(dx_th_3);
             case 4: return(dx_A_2);
             case 5: return(dx_A_3);

             default: return(0.0);
         }
         break;
         
      case 1: // y derivative
         switch (component)
         {
             case 0: return(dy_th_1);
             case 1: return(dy_th_2);
             case 2: return(dy_th_3);
             case 3: return(dy_A_1);
             case 5: return(dy_A_3);

             default: return(0.0);
         }
         break;
         
      case 2: // z derivative
         switch (component)
         {
             case 0: return(dz_th_1);
             case 1: return(dz_th_2);
             case 2: return(dz_th_3);
             case 3: return(dz_A_1);
             case 4: return(dz_A_2);

             default: return(0.0);
         }
         break;

      default: return(0.0);
   }
}

void TAux::set_derivative(int axis, int component, double value)
{
   switch (axis)
   {
      case 0: // x derivative
         switch (component)
         {
             case 0: dx_th_1 = value; break;
             case 1: dx_th_2 = value; break;
             case 2: dx_th_3 = value; break;
             case 4: dx_A_2 = value; break;
             case 5: dx_A_3 = value; break;
         }
         break;
         
      case 1: // y derivative
         switch (component)
         {
             case 0: dy_th_1 = value; break;
             case 1: dy_th_2 = value; break;
             case 2: dy_th_3 = value; break;
             case 3: dy_A_1 = value; break;
             case 5: dy_A_3 = value; break;
         }
         break;
         
      case 2: // z derivative
         switch (component)
         {
             case 0: dz_th_1 = value; break;
             case 1: dz_th_2 = value; break;
             case 2: dz_th_3 = value; break;
             case 3: dz_A_1 = value; break;
             case 4: dz_A_2 = value; break;
         }
         break;
   }
}

void TAux::compute_derivative(TFields fields[9][9][9], int i, int j, int k, int axis, int component)
/* IMPORTANT: Does NOT divide by lattice spacing, so if absolute (rather than relative) magnitude of derivative is needed, user must do that. */
{
   int iM1 = i - 1;
   int jM1 = j - 1;
   int kM1 = k - 1;

   int iM2 = iM1 - 1;
   int jM2 = jM1 - 1;
   int kM2 = kM1 - 1;

   int iP1 = i + 1;
   int jP1 = j + 1;
   int kP1 = k + 1;

   int iP2 = iP1 + 1;
   int jP2 = jP1 + 1;
   int kP2 = kP1 + 1;

   switch (axis)
   {
      case 0: // x derivative

         switch (component)
         {
             case 0: // th_1
                  dx_th_1 = 
                    (  fields[iM1][jM1][kM2].th_1 + fields[iM1][jM1][kP2].th_1 + fields[iM1][jP1][kM2].th_1 + fields[iM1][jP1][kP2].th_1
                     + fields[iM1][jM2][kM1].th_1 + fields[iM1][jM2][kP1].th_1 + fields[iM1][jP2][kM1].th_1 + fields[iM1][jP2][kP1].th_1
                     - fields[iP1][jM1][kM2].th_1 - fields[iP1][jM1][kP2].th_1 - fields[iP1][jP1][kM2].th_1 - fields[iP1][jP1][kP2].th_1
                     - fields[iP1][jM2][kM1].th_1 - fields[iP1][jM2][kP1].th_1 - fields[iP1][jP2][kM1].th_1 - fields[iP1][jP2][kP1].th_1)/120
                  + (  fields[iM2][j][kM1].th_1 + fields[iM2][j][kP1].th_1 + fields[iM2][jM1][k].th_1 + fields[iM2][jP1][k].th_1
                     - fields[iP2][j][kM1].th_1 - fields[iP2][j][kP1].th_1 - fields[iP2][jM1][k].th_1 - fields[iP2][jP1][k].th_1
                     - fields[iM1][jM1][kM1].th_1 - fields[iM1][jM1][kP1].th_1 - fields[iM1][jP1][kM1].th_1 - fields[iM1][jP1][kP1].th_1
                     + fields[iP1][jM1][kM1].th_1 + fields[iP1][jM1][kP1].th_1 + fields[iP1][jP1][kM1].th_1 + fields[iP1][jP1][kP1].th_1)/30
                  + (- fields[iM2][j][k].th_1 + fields[iP2][j][k].th_1)/20
                  - (  fields[iM1][j][kM1].th_1 + fields[iM1][j][kP1].th_1 + fields[iM1][jM1][k].th_1 + fields[iM1][jP1][k].th_1
                     - fields[iP1][j][kM1].th_1 - fields[iP1][j][kP1].th_1 - fields[iP1][jM1][k].th_1 - fields[iP1][jP1][k].th_1)/12
                  + (- fields[iM1][j][k].th_1 + fields[iP1][j][k].th_1)*4.0/15;
#ifdef WRAP_DTH
                  dx_th_1 = wrap(dx_th_1);
#endif
                break;
             case 1: // th_2
                  dx_th_2 = 
                    (  fields[iM1][jM1][kM2].th_2 + fields[iM1][jM1][kP2].th_2 + fields[iM1][jP1][kM2].th_2 + fields[iM1][jP1][kP2].th_2
                     + fields[iM1][jM2][kM1].th_2 + fields[iM1][jM2][kP1].th_2 + fields[iM1][jP2][kM1].th_2 + fields[iM1][jP2][kP1].th_2
                     - fields[iP1][jM1][kM2].th_2 - fields[iP1][jM1][kP2].th_2 - fields[iP1][jP1][kM2].th_2 - fields[iP1][jP1][kP2].th_2
                     - fields[iP1][jM2][kM1].th_2 - fields[iP1][jM2][kP1].th_2 - fields[iP1][jP2][kM1].th_2 - fields[iP1][jP2][kP1].th_2)/120
                  + (  fields[iM2][j][kM1].th_2 + fields[iM2][j][kP1].th_2 + fields[iM2][jM1][k].th_2 + fields[iM2][jP1][k].th_2
                     - fields[iP2][j][kM1].th_2 - fields[iP2][j][kP1].th_2 - fields[iP2][jM1][k].th_2 - fields[iP2][jP1][k].th_2
                     - fields[iM1][jM1][kM1].th_2 - fields[iM1][jM1][kP1].th_2 - fields[iM1][jP1][kM1].th_2 - fields[iM1][jP1][kP1].th_2
                     + fields[iP1][jM1][kM1].th_2 + fields[iP1][jM1][kP1].th_2 + fields[iP1][jP1][kM1].th_2 + fields[iP1][jP1][kP1].th_2)/30
                  - (  fields[iM2][j][k].th_2 - fields[iP2][j][k].th_2)/20
                  + (- fields[iM1][j][kM1].th_2 - fields[iM1][j][kP1].th_2 - fields[iM1][jM1][k].th_2 - fields[iM1][jP1][k].th_2
                     + fields[iP1][j][kM1].th_2 + fields[iP1][j][kP1].th_2 + fields[iP1][jM1][k].th_2 + fields[iP1][jP1][k].th_2)/12
                  - (  fields[iM1][j][k].th_2 - fields[iP1][j][k].th_2)*4.0/15;
#ifdef WRAP_DTH
                  dx_th_2 = wrap(dx_th_2);
#endif
                break;
             case 2: // th_3
                  dx_th_3 = 
                    (  fields[iM1][jM1][kM2].th_3 + fields[iM1][jM1][kP2].th_3 + fields[iM1][jP1][kM2].th_3 + fields[iM1][jP1][kP2].th_3
                     + fields[iM1][jM2][kM1].th_3 + fields[iM1][jM2][kP1].th_3 + fields[iM1][jP2][kM1].th_3 + fields[iM1][jP2][kP1].th_3
                     - fields[iP1][jM1][kM2].th_3 - fields[iP1][jM1][kP2].th_3 - fields[iP1][jP1][kM2].th_3 - fields[iP1][jP1][kP2].th_3
                     - fields[iP1][jM2][kM1].th_3 - fields[iP1][jM2][kP1].th_3 - fields[iP1][jP2][kM1].th_3 - fields[iP1][jP2][kP1].th_3)/120
                  + (  fields[iM2][j][kM1].th_3 + fields[iM2][j][kP1].th_3 + fields[iM2][jM1][k].th_3 + fields[iM2][jP1][k].th_3
                     - fields[iP2][j][kM1].th_3 - fields[iP2][j][kP1].th_3 - fields[iP2][jM1][k].th_3 - fields[iP2][jP1][k].th_3
                     - fields[iM1][jM1][kM1].th_3 - fields[iM1][jM1][kP1].th_3 - fields[iM1][jP1][kM1].th_3 - fields[iM1][jP1][kP1].th_3
                     + fields[iP1][jM1][kM1].th_3 + fields[iP1][jM1][kP1].th_3 + fields[iP1][jP1][kM1].th_3 + fields[iP1][jP1][kP1].th_3)/30
                  + (- fields[iM2][j][k].th_3 + fields[iP2][j][k].th_3)/20
                  + (- fields[iM1][j][kM1].th_3 - fields[iM1][j][kP1].th_3 - fields[iM1][jM1][k].th_3 - fields[iM1][jP1][k].th_3
                     + fields[iP1][j][kM1].th_3 + fields[iP1][j][kP1].th_3 + fields[iP1][jM1][k].th_3 + fields[iP1][jP1][k].th_3)/12
                  + (- fields[iM1][j][k].th_3 + fields[iP1][j][k].th_3)*4.0/15;
#ifdef WRAP_DTH
                  dx_th_3 = wrap(dx_th_3);
#endif
                break;
             case 4: // A_2
                  dx_A_2 = 
                   (  fields[iM1][jM1][kM2].A_2 + fields[iM1][jM1][kP2].A_2 + fields[iM1][jP1][kM2].A_2 + fields[iM1][jP1][kP2].A_2
                    + fields[iM1][jM2][kM1].A_2 + fields[iM1][jM2][kP1].A_2 + fields[iM1][jP2][kM1].A_2 + fields[iM1][jP2][kP1].A_2
                    - fields[iP1][jM1][kM2].A_2 - fields[iP1][jM1][kP2].A_2 - fields[iP1][jP1][kM2].A_2 - fields[iP1][jP1][kP2].A_2
                    - fields[iP1][jM2][kM1].A_2 - fields[iP1][jM2][kP1].A_2 - fields[iP1][jP2][kM1].A_2 - fields[iP1][jP2][kP1].A_2)/120
                 + (  fields[iM2][j][kM1].A_2 + fields[iM2][j][kP1].A_2 + fields[iM2][jM1][k].A_2 + fields[iM2][jP1][k].A_2
                    - fields[iP2][j][kM1].A_2 - fields[iP2][j][kP1].A_2 - fields[iP2][jM1][k].A_2 - fields[iP2][jP1][k].A_2
                    - fields[iM1][jM1][kM1].A_2 - fields[iM1][jM1][kP1].A_2 - fields[iM1][jP1][kM1].A_2 - fields[iM1][jP1][kP1].A_2
                    + fields[iP1][jM1][kM1].A_2 + fields[iP1][jM1][kP1].A_2 + fields[iP1][jP1][kM1].A_2 + fields[iP1][jP1][kP1].A_2)/30
                 + (- fields[iM2][j][k].A_2 + fields[iP2][j][k].A_2)/20
                 + (- fields[iM1][j][kM1].A_2 - fields[iM1][j][kP1].A_2 - fields[iM1][jM1][k].A_2 - fields[iM1][jP1][k].A_2
                    + fields[iP1][j][kM1].A_2 + fields[iP1][j][kP1].A_2 + fields[iP1][jM1][k].A_2 + fields[iP1][jP1][k].A_2)/12
                 + (- fields[iM1][j][k].A_2 + fields[iP1][j][k].A_2)*4.0/15;
                break;
             case 5: // A_3
                  dx_A_3 = 
                   (  fields[iM1][jM1][kM2].A_3 + fields[iM1][jM1][kP2].A_3 + fields[iM1][jP1][kM2].A_3 + fields[iM1][jP1][kP2].A_3
                    + fields[iM1][jM2][kM1].A_3 + fields[iM1][jM2][kP1].A_3 + fields[iM1][jP2][kM1].A_3 + fields[iM1][jP2][kP1].A_3
                    - fields[iP1][jM1][kM2].A_3 - fields[iP1][jM1][kP2].A_3 - fields[iP1][jP1][kM2].A_3 - fields[iP1][jP1][kP2].A_3
                    - fields[iP1][jM2][kM1].A_3 - fields[iP1][jM2][kP1].A_3 - fields[iP1][jP2][kM1].A_3 - fields[iP1][jP2][kP1].A_3)/120
                 + (  fields[iM2][j][kM1].A_3 + fields[iM2][j][kP1].A_3 + fields[iM2][jM1][k].A_3 + fields[iM2][jP1][k].A_3
                    - fields[iP2][j][kM1].A_3 - fields[iP2][j][kP1].A_3 - fields[iP2][jM1][k].A_3 - fields[iP2][jP1][k].A_3
                    - fields[iM1][jM1][kM1].A_3 - fields[iM1][jM1][kP1].A_3 - fields[iM1][jP1][kM1].A_3 - fields[iM1][jP1][kP1].A_3
                    + fields[iP1][jM1][kM1].A_3 + fields[iP1][jM1][kP1].A_3 + fields[iP1][jP1][kM1].A_3 + fields[iP1][jP1][kP1].A_3)/30
                 + (- fields[iM2][j][k].A_3 + fields[iP2][j][k].A_3)/20
                 - (  fields[iM1][j][kM1].A_3 + fields[iM1][j][kP1].A_3 + fields[iM1][jM1][k].A_3 + fields[iM1][jP1][k].A_3
                    - fields[iP1][j][kM1].A_3 - fields[iP1][j][kP1].A_3 - fields[iP1][jM1][k].A_3 - fields[iP1][jP1][k].A_3)/12
                 + (- fields[iM1][j][k].A_3 + fields[iP1][j][k].A_3)*4.0/15;
               break;
         }
         break;
         
      case 1: // y derivative

         switch (component)
         {
             case 0: // th_1
                  dy_th_1 = 
                    (  fields[iM1][jM1][kM2].th_1 + fields[iM1][jM1][kP2].th_1 - fields[iM1][jP1][kM2].th_1 - fields[iM1][jP1][kP2].th_1
                     + fields[iP1][jM1][kM2].th_1 + fields[iP1][jM1][kP2].th_1 - fields[iP1][jP1][kM2].th_1 - fields[iP1][jP1][kP2].th_1
                     + fields[iM2][jM1][kM1].th_1 + fields[iM2][jM1][kP1].th_1 - fields[iM2][jP1][kM1].th_1 - fields[iM2][jP1][kP1].th_1
                     + fields[iP2][jM1][kM1].th_1 + fields[iP2][jM1][kP1].th_1 - fields[iP2][jP1][kM1].th_1 - fields[iP2][jP1][kP1].th_1)/120
                  + (  fields[i][jM2][kM1].th_1 + fields[i][jM2][kP1].th_1 - fields[i][jP2][kM1].th_1 - fields[i][jP2][kP1].th_1
                     + fields[iM1][jM2][k].th_1 - fields[iM1][jP2][k].th_1 + fields[iP1][jM2][k].th_1 - fields[iP1][jP2][k].th_1
                     - fields[iM1][jM1][kM1].th_1 - fields[iM1][jM1][kP1].th_1 + fields[iM1][jP1][kM1].th_1 + fields[iM1][jP1][kP1].th_1
                     - fields[iP1][jM1][kM1].th_1 - fields[iP1][jM1][kP1].th_1 + fields[iP1][jP1][kM1].th_1 + fields[iP1][jP1][kP1].th_1)/30
                  + (- fields[i][jM2][k].th_1 + fields[i][jP2][k].th_1)/20
                  - (  fields[i][jM1][kM1].th_1 + fields[i][jM1][kP1].th_1 - fields[i][jP1][kM1].th_1 - fields[i][jP1][kP1].th_1
                     + fields[iM1][jM1][k].th_1 - fields[iM1][jP1][k].th_1 + fields[iP1][jM1][k].th_1 - fields[iP1][jP1][k].th_1)/12
                  + (- fields[i][jM1][k].th_1 + fields[i][jP1][k].th_1)*4.0/15;
#ifdef WRAP_DTH
                  dy_th_1 = wrap(dy_th_1);
#endif
                break;
             case 1: // th_2
                  dy_th_2 = 
                    (  fields[iM1][jM1][kM2].th_2 + fields[iM1][jM1][kP2].th_2 - fields[iM1][jP1][kM2].th_2 - fields[iM1][jP1][kP2].th_2
                     + fields[iP1][jM1][kM2].th_2 + fields[iP1][jM1][kP2].th_2 - fields[iP1][jP1][kM2].th_2 - fields[iP1][jP1][kP2].th_2
                     + fields[iM2][jM1][kM1].th_2 + fields[iM2][jM1][kP1].th_2 - fields[iM2][jP1][kM1].th_2 - fields[iM2][jP1][kP1].th_2
                     + fields[iP2][jM1][kM1].th_2 + fields[iP2][jM1][kP1].th_2 - fields[iP2][jP1][kM1].th_2 - fields[iP2][jP1][kP1].th_2)/120
                  + (  fields[i][jM2][kM1].th_2 + fields[i][jM2][kP1].th_2 - fields[i][jP2][kM1].th_2 - fields[i][jP2][kP1].th_2
                     + fields[iM1][jM2][k].th_2 - fields[iM1][jP2][k].th_2 + fields[iP1][jM2][k].th_2 - fields[iP1][jP2][k].th_2
                     - fields[iM1][jM1][kM1].th_2 - fields[iM1][jM1][kP1].th_2 + fields[iM1][jP1][kM1].th_2 + fields[iM1][jP1][kP1].th_2
                     - fields[iP1][jM1][kM1].th_2 - fields[iP1][jM1][kP1].th_2 + fields[iP1][jP1][kM1].th_2 + fields[iP1][jP1][kP1].th_2)/30
                  + (- fields[i][jM2][k].th_2 + fields[i][jP2][k].th_2)/20
                  + (- fields[i][jM1][kM1].th_2 - fields[i][jM1][kP1].th_2 + fields[i][jP1][kM1].th_2 + fields[i][jP1][kP1].th_2
                     - fields[iM1][jM1][k].th_2 + fields[iM1][jP1][k].th_2 - fields[iP1][jM1][k].th_2 + fields[iP1][jP1][k].th_2)/12
                  + (- fields[i][jM1][k].th_2 + fields[i][jP1][k].th_2)*4.0/15;
#ifdef WRAP_DTH
                  dy_th_2 = wrap(dy_th_2);
#endif
                break;
             case 2: // th_3
                  dy_th_3 = 
                    (  fields[iM1][jM1][kM2].th_3 + fields[iM1][jM1][kP2].th_3 - fields[iM1][jP1][kM2].th_3 - fields[iM1][jP1][kP2].th_3
                     + fields[iP1][jM1][kM2].th_3 + fields[iP1][jM1][kP2].th_3 - fields[iP1][jP1][kM2].th_3 - fields[iP1][jP1][kP2].th_3
                     + fields[iM2][jM1][kM1].th_3 + fields[iM2][jM1][kP1].th_3 - fields[iM2][jP1][kM1].th_3 - fields[iM2][jP1][kP1].th_3
                     + fields[iP2][jM1][kM1].th_3 + fields[iP2][jM1][kP1].th_3 - fields[iP2][jP1][kM1].th_3 - fields[iP2][jP1][kP1].th_3)/120
                  + (  fields[i][jM2][kM1].th_3 + fields[i][jM2][kP1].th_3 - fields[i][jP2][kM1].th_3 - fields[i][jP2][kP1].th_3
                     + fields[iM1][jM2][k].th_3 - fields[iM1][jP2][k].th_3 + fields[iP1][jM2][k].th_3 - fields[iP1][jP2][k].th_3
                     - fields[iM1][jM1][kM1].th_3 - fields[iM1][jM1][kP1].th_3 + fields[iM1][jP1][kM1].th_3 + fields[iM1][jP1][kP1].th_3
                     - fields[iP1][jM1][kM1].th_3 - fields[iP1][jM1][kP1].th_3 + fields[iP1][jP1][kM1].th_3 + fields[iP1][jP1][kP1].th_3)/30
                  - (  fields[i][jM2][k].th_3 - fields[i][jP2][k].th_3)/20
                  + (- fields[i][jM1][kM1].th_3 - fields[i][jM1][kP1].th_3 + fields[i][jP1][kM1].th_3 + fields[i][jP1][kP1].th_3
                     - fields[iM1][jM1][k].th_3 + fields[iM1][jP1][k].th_3 - fields[iP1][jM1][k].th_3 + fields[iP1][jP1][k].th_3)/12
                  - (  fields[i][jM1][k].th_3 - fields[i][jP1][k].th_3)*4.0/15;
#ifdef WRAP_DTH
                  dy_th_3 = wrap(dy_th_3);
#endif
                break;
             case 3: // A_1
                 dy_A_1 = 
                   (  fields[iM1][jM1][kM2].A_1 + fields[iM1][jM1][kP2].A_1 - fields[iM1][jP1][kM2].A_1 - fields[iM1][jP1][kP2].A_1
                    + fields[iP1][jM1][kM2].A_1 + fields[iP1][jM1][kP2].A_1 - fields[iP1][jP1][kM2].A_1 - fields[iP1][jP1][kP2].A_1
                    + fields[iM2][jM1][kM1].A_1 + fields[iM2][jM1][kP1].A_1 - fields[iM2][jP1][kM1].A_1 - fields[iM2][jP1][kP1].A_1
                    + fields[iP2][jM1][kM1].A_1 + fields[iP2][jM1][kP1].A_1 - fields[iP2][jP1][kM1].A_1 - fields[iP2][jP1][kP1].A_1)/120
                 + (  fields[i][jM2][kM1].A_1 + fields[i][jM2][kP1].A_1 - fields[i][jP2][kM1].A_1 - fields[i][jP2][kP1].A_1
                    + fields[iM1][jM2][k].A_1 - fields[iM1][jP2][k].A_1 + fields[iP1][jM2][k].A_1 - fields[iP1][jP2][k].A_1
                    - fields[iM1][jM1][kM1].A_1 - fields[iM1][jM1][kP1].A_1 + fields[iM1][jP1][kM1].A_1 + fields[iM1][jP1][kP1].A_1
                    - fields[iP1][jM1][kM1].A_1 - fields[iP1][jM1][kP1].A_1 + fields[iP1][jP1][kM1].A_1 + fields[iP1][jP1][kP1].A_1)/30
                 + (- fields[i][jM2][k].A_1 + fields[i][jP2][k].A_1)/20
                 - (  fields[i][jM1][kM1].A_1 + fields[i][jM1][kP1].A_1 - fields[i][jP1][kM1].A_1 - fields[i][jP1][kP1].A_1
                    + fields[iM1][jM1][k].A_1 - fields[iM1][jP1][k].A_1 + fields[iP1][jM1][k].A_1 - fields[iP1][jP1][k].A_1)/12
                 + (- fields[i][jM1][k].A_1 + fields[i][jP1][k].A_1)*4.0/15;
                break;
             case 5: // A_3
                 dy_A_3 = 
                   (  fields[iM1][jM1][kM2].A_3 + fields[iM1][jM1][kP2].A_3 - fields[iM1][jP1][kM2].A_3 - fields[iM1][jP1][kP2].A_3
                    + fields[iP1][jM1][kM2].A_3 + fields[iP1][jM1][kP2].A_3 - fields[iP1][jP1][kM2].A_3 - fields[iP1][jP1][kP2].A_3
                    + fields[iM2][jM1][kM1].A_3 + fields[iM2][jM1][kP1].A_3 - fields[iM2][jP1][kM1].A_3 - fields[iM2][jP1][kP1].A_3
                    + fields[iP2][jM1][kM1].A_3 + fields[iP2][jM1][kP1].A_3 - fields[iP2][jP1][kM1].A_3 - fields[iP2][jP1][kP1].A_3)/120
                 + (  fields[i][jM2][kM1].A_3 + fields[i][jM2][kP1].A_3 - fields[i][jP2][kM1].A_3 - fields[i][jP2][kP1].A_3
                    + fields[iM1][jM2][k].A_3 - fields[iM1][jP2][k].A_3 + fields[iP1][jM2][k].A_3 - fields[iP1][jP2][k].A_3
                    - fields[iM1][jM1][kM1].A_3 - fields[iM1][jM1][kP1].A_3 + fields[iM1][jP1][kM1].A_3 + fields[iM1][jP1][kP1].A_3
                    - fields[iP1][jM1][kM1].A_3 - fields[iP1][jM1][kP1].A_3 + fields[iP1][jP1][kM1].A_3 + fields[iP1][jP1][kP1].A_3)/30
                 - (  fields[i][jM2][k].A_3 - fields[i][jP2][k].A_3)/20
                 + (- fields[i][jM1][kM1].A_3 - fields[i][jM1][kP1].A_3 + fields[i][jP1][kM1].A_3 + fields[i][jP1][kP1].A_3
                    - fields[iM1][jM1][k].A_3 + fields[iM1][jP1][k].A_3 - fields[iP1][jM1][k].A_3 + fields[iP1][jP1][k].A_3)/12
                 - (  fields[i][jM1][k].A_3 - fields[i][jP1][k].A_3)*4.0/15;
      break;
         }
         break;
         
      case 2: // z derivative

         switch (component)
         {
             case 0: // th_1
                  dz_th_1 = 
                    (  fields[iM1][jM2][kM1].th_1 - fields[iM1][jM2][kP1].th_1 + fields[iM1][jP2][kM1].th_1 - fields[iM1][jP2][kP1].th_1
                     + fields[iP1][jM2][kM1].th_1 - fields[iP1][jM2][kP1].th_1 + fields[iP1][jP2][kM1].th_1 - fields[iP1][jP2][kP1].th_1
                     + fields[iM2][jM1][kM1].th_1 - fields[iM2][jM1][kP1].th_1 + fields[iM2][jP1][kM1].th_1 - fields[iM2][jP1][kP1].th_1
                     + fields[iP2][jM1][kM1].th_1 - fields[iP2][jM1][kP1].th_1 + fields[iP2][jP1][kM1].th_1 - fields[iP2][jP1][kP1].th_1)/120
                  + (  fields[i][jM1][kM2].th_1 - fields[i][jM1][kP2].th_1 + fields[i][jP1][kM2].th_1 - fields[i][jP1][kP2].th_1
                     + fields[iM1][j][kM2].th_1 - fields[iM1][j][kP2].th_1 + fields[iP1][j][kM2].th_1 - fields[iP1][j][kP2].th_1
                     - fields[iM1][jM1][kM1].th_1 + fields[iM1][jM1][kP1].th_1 - fields[iM1][jP1][kM1].th_1 + fields[iM1][jP1][kP1].th_1
                     - fields[iP1][jM1][kM1].th_1 + fields[iP1][jM1][kP1].th_1 - fields[iP1][jP1][kM1].th_1 + fields[iP1][jP1][kP1].th_1)/30
                  + (- fields[i][j][kM2].th_1 + fields[i][j][kP2].th_1)/20
                  - (  fields[i][jM1][kM1].th_1 - fields[i][jM1][kP1].th_1 + fields[i][jP1][kM1].th_1 - fields[i][jP1][kP1].th_1
                     + fields[iM1][j][kM1].th_1 - fields[iM1][j][kP1].th_1 + fields[iP1][j][kM1].th_1 - fields[iP1][j][kP1].th_1)/12
                  + (- fields[i][j][kM1].th_1 + fields[i][j][kP1].th_1)*4.0/15;
#ifdef WRAP_DTH
                  dz_th_1 = wrap(dz_th_1);
#endif
                break;
             case 1: // th_2
                  dz_th_2 = 
                    (  fields[iM1][jM2][kM1].th_2 - fields[iM1][jM2][kP1].th_2 + fields[iM1][jP2][kM1].th_2 - fields[iM1][jP2][kP1].th_2
                     + fields[iP1][jM2][kM1].th_2 - fields[iP1][jM2][kP1].th_2 + fields[iP1][jP2][kM1].th_2 - fields[iP1][jP2][kP1].th_2
                     + fields[iM2][jM1][kM1].th_2 - fields[iM2][jM1][kP1].th_2 + fields[iM2][jP1][kM1].th_2 - fields[iM2][jP1][kP1].th_2
                     + fields[iP2][jM1][kM1].th_2 - fields[iP2][jM1][kP1].th_2 + fields[iP2][jP1][kM1].th_2 - fields[iP2][jP1][kP1].th_2)/120
                  + (  fields[i][jM1][kM2].th_2 - fields[i][jM1][kP2].th_2 + fields[i][jP1][kM2].th_2 - fields[i][jP1][kP2].th_2
                     + fields[iM1][j][kM2].th_2 - fields[iM1][j][kP2].th_2 + fields[iP1][j][kM2].th_2 - fields[iP1][j][kP2].th_2
                     - fields[iM1][jM1][kM1].th_2 + fields[iM1][jM1][kP1].th_2 - fields[iM1][jP1][kM1].th_2 + fields[iM1][jP1][kP1].th_2
                     - fields[iP1][jM1][kM1].th_2 + fields[iP1][jM1][kP1].th_2 - fields[iP1][jP1][kM1].th_2 + fields[iP1][jP1][kP1].th_2)/30
                  + (- fields[i][j][kM2].th_2 + fields[i][j][kP2].th_2)/20
                  - (  fields[i][jM1][kM1].th_2 - fields[i][jM1][kP1].th_2 + fields[i][jP1][kM1].th_2 - fields[i][jP1][kP1].th_2
                     + fields[iM1][j][kM1].th_2 - fields[iM1][j][kP1].th_2 + fields[iP1][j][kM1].th_2 - fields[iP1][j][kP1].th_2)/12
                  + (- fields[i][j][kM1].th_2 + fields[i][j][kP1].th_2)*4.0/15;
#ifdef WRAP_DTH
                  dz_th_2 = wrap(dz_th_2);
#endif
                break;
             case 2: // th_3
                  dz_th_3 = 
                    (  fields[iM1][jM2][kM1].th_3 - fields[iM1][jM2][kP1].th_3 + fields[iM1][jP2][kM1].th_3 - fields[iM1][jP2][kP1].th_3
                     + fields[iP1][jM2][kM1].th_3 - fields[iP1][jM2][kP1].th_3 + fields[iP1][jP2][kM1].th_3 - fields[iP1][jP2][kP1].th_3
                     + fields[iM2][jM1][kM1].th_3 - fields[iM2][jM1][kP1].th_3 + fields[iM2][jP1][kM1].th_3 - fields[iM2][jP1][kP1].th_3
                     + fields[iP2][jM1][kM1].th_3 - fields[iP2][jM1][kP1].th_3 + fields[iP2][jP1][kM1].th_3 - fields[iP2][jP1][kP1].th_3)/120
                  + (  fields[i][jM1][kM2].th_3 - fields[i][jM1][kP2].th_3 + fields[i][jP1][kM2].th_3 - fields[i][jP1][kP2].th_3
                     + fields[iM1][j][kM2].th_3 - fields[iM1][j][kP2].th_3 + fields[iP1][j][kM2].th_3 - fields[iP1][j][kP2].th_3
                     - fields[iM1][jM1][kM1].th_3 + fields[iM1][jM1][kP1].th_3 - fields[iM1][jP1][kM1].th_3 + fields[iM1][jP1][kP1].th_3
                     - fields[iP1][jM1][kM1].th_3 + fields[iP1][jM1][kP1].th_3 - fields[iP1][jP1][kM1].th_3 + fields[iP1][jP1][kP1].th_3)/30
                  + (- fields[i][j][kM2].th_3 + fields[i][j][kP2].th_3)/20
                  - (  fields[i][jM1][kM1].th_3 - fields[i][jM1][kP1].th_3 + fields[i][jP1][kM1].th_3 - fields[i][jP1][kP1].th_3
                     + fields[iM1][j][kM1].th_3 - fields[iM1][j][kP1].th_3 + fields[iP1][j][kM1].th_3 - fields[iP1][j][kP1].th_3)/12
                  + (- fields[i][j][kM1].th_3 + fields[i][j][kP1].th_3)*4.0/15;
#ifdef WRAP_DTH
                  dz_th_3 = wrap(dz_th_3);
#endif
                break;
             case 3: // A_1
                 dz_A_1 =
                   (  fields[iM1][jM2][kM1].A_1 - fields[iM1][jM2][kP1].A_1 + fields[iM1][jP2][kM1].A_1 - fields[iM1][jP2][kP1].A_1
                    + fields[iP1][jM2][kM1].A_1 - fields[iP1][jM2][kP1].A_1 + fields[iP1][jP2][kM1].A_1 - fields[iP1][jP2][kP1].A_1
                    + fields[iM2][jM1][kM1].A_1 - fields[iM2][jM1][kP1].A_1 + fields[iM2][jP1][kM1].A_1 - fields[iM2][jP1][kP1].A_1
                    + fields[iP2][jM1][kM1].A_1 - fields[iP2][jM1][kP1].A_1 + fields[iP2][jP1][kM1].A_1 - fields[iP2][jP1][kP1].A_1)/120
                 + (  fields[i][jM1][kM2].A_1 - fields[i][jM1][kP2].A_1 + fields[i][jP1][kM2].A_1 - fields[i][jP1][kP2].A_1
                    + fields[iM1][j][kM2].A_1 - fields[iM1][j][kP2].A_1 + fields[iP1][j][kM2].A_1 - fields[iP1][j][kP2].A_1
                    - fields[iM1][jM1][kM1].A_1 + fields[iM1][jM1][kP1].A_1 - fields[iM1][jP1][kM1].A_1 + fields[iM1][jP1][kP1].A_1
                    - fields[iP1][jM1][kM1].A_1 + fields[iP1][jM1][kP1].A_1 - fields[iP1][jP1][kM1].A_1 + fields[iP1][jP1][kP1].A_1)/30
                 + (- fields[i][j][kM2].A_1 + fields[i][j][kP2].A_1)/20
                 - (  fields[i][jM1][kM1].A_1 - fields[i][jM1][kP1].A_1 + fields[i][jP1][kM1].A_1 - fields[i][jP1][kP1].A_1
                    + fields[iM1][j][kM1].A_1 - fields[iM1][j][kP1].A_1 + fields[iP1][j][kM1].A_1 - fields[iP1][j][kP1].A_1)/12
                 + (- fields[i][j][kM1].A_1 + fields[i][j][kP1].A_1)*4.0/15;
                break;
             case 4: // A_2
                 dz_A_2 =
                   (  fields[iM1][jM2][kM1].A_2 - fields[iM1][jM2][kP1].A_2 + fields[iM1][jP2][kM1].A_2 - fields[iM1][jP2][kP1].A_2
                    + fields[iP1][jM2][kM1].A_2 - fields[iP1][jM2][kP1].A_2 + fields[iP1][jP2][kM1].A_2 - fields[iP1][jP2][kP1].A_2
                    + fields[iM2][jM1][kM1].A_2 - fields[iM2][jM1][kP1].A_2 + fields[iM2][jP1][kM1].A_2 - fields[iM2][jP1][kP1].A_2
                    + fields[iP2][jM1][kM1].A_2 - fields[iP2][jM1][kP1].A_2 + fields[iP2][jP1][kM1].A_2 - fields[iP2][jP1][kP1].A_2)/120
                 + (  fields[i][jM1][kM2].A_2 - fields[i][jM1][kP2].A_2 + fields[i][jP1][kM2].A_2 - fields[i][jP1][kP2].A_2
                    + fields[iM1][j][kM2].A_2 - fields[iM1][j][kP2].A_2 + fields[iP1][j][kM2].A_2 - fields[iP1][j][kP2].A_2
                    - fields[iM1][jM1][kM1].A_2 + fields[iM1][jM1][kP1].A_2 - fields[iM1][jP1][kM1].A_2 + fields[iM1][jP1][kP1].A_2
                    - fields[iP1][jM1][kM1].A_2 + fields[iP1][jM1][kP1].A_2 - fields[iP1][jP1][kM1].A_2 + fields[iP1][jP1][kP1].A_2)/30
                 + (- fields[i][j][kM2].A_2 + fields[i][j][kP2].A_2)/20
                 - (  fields[i][jM1][kM1].A_2 - fields[i][jM1][kP1].A_2 + fields[i][jP1][kM1].A_2 - fields[i][jP1][kP1].A_2
                    + fields[iM1][j][kM1].A_2 - fields[iM1][j][kP1].A_2 + fields[iP1][j][kM1].A_2 - fields[iP1][j][kP1].A_2)/12
                 + (- fields[i][j][kM1].A_2 + fields[i][j][kP1].A_2)*4.0/15;
                break;
         }
         break;
   }
}

void TAux::compute_derivative(TFields fields[MAX_NSIDE][MAX_NSIDE][MAX_NSIDE], int i, int j, int k, int nSideM1,
                              int axis, int component)
/* IMPORTANT: Does NOT divide by lattice spacing, so if absolute (rather than relative) magnitude of derivative is needed, user must do that. */
{
   int iM1, jM1, kM1, iM2, jM2, kM2, iP1, jP1, kP1, iP2, jP2, kP2;

   if (0 == i) iM1 = nSideM1; else iM1 = i - 1;
   if (0 == j) jM1 = nSideM1; else jM1 = j - 1;
   if (0 == k) kM1 = nSideM1; else kM1 = k - 1;

   if (0 == iM1) iM2 = nSideM1; else iM2 = iM1 - 1;
   if (0 == jM1) jM2 = nSideM1; else jM2 = jM1 - 1;
   if (0 == kM1) kM2 = nSideM1; else kM2 = kM1 - 1;

   if (nSideM1 == i) iP1 = 0; else iP1 = i + 1;
   if (nSideM1 == j) jP1 = 0; else jP1 = j + 1;
   if (nSideM1 == k) kP1 = 0; else kP1 = k + 1;

   if (nSideM1 == iP1) iP2 = 0; else iP2 = iP1 + 1;
   if (nSideM1 == jP1) jP2 = 0; else jP2 = jP1 + 1;
   if (nSideM1 == kP1) kP2 = 0; else kP2 = kP1 + 1;

   switch (axis)
   {
      case 0: // x derivative

         switch (component)
         {
             case 0: // th_1
                  dx_th_1 = 
                    (  fields[iM1][jM1][kM2].th_1 + fields[iM1][jM1][kP2].th_1 + fields[iM1][jP1][kM2].th_1 + fields[iM1][jP1][kP2].th_1
                     + fields[iM1][jM2][kM1].th_1 + fields[iM1][jM2][kP1].th_1 + fields[iM1][jP2][kM1].th_1 + fields[iM1][jP2][kP1].th_1
                     - fields[iP1][jM1][kM2].th_1 - fields[iP1][jM1][kP2].th_1 - fields[iP1][jP1][kM2].th_1 - fields[iP1][jP1][kP2].th_1
                     - fields[iP1][jM2][kM1].th_1 - fields[iP1][jM2][kP1].th_1 - fields[iP1][jP2][kM1].th_1 - fields[iP1][jP2][kP1].th_1)/120
                  + (  fields[iM2][j][kM1].th_1 + fields[iM2][j][kP1].th_1 + fields[iM2][jM1][k].th_1 + fields[iM2][jP1][k].th_1
                     - fields[iP2][j][kM1].th_1 - fields[iP2][j][kP1].th_1 - fields[iP2][jM1][k].th_1 - fields[iP2][jP1][k].th_1
                     - fields[iM1][jM1][kM1].th_1 - fields[iM1][jM1][kP1].th_1 - fields[iM1][jP1][kM1].th_1 - fields[iM1][jP1][kP1].th_1
                     + fields[iP1][jM1][kM1].th_1 + fields[iP1][jM1][kP1].th_1 + fields[iP1][jP1][kM1].th_1 + fields[iP1][jP1][kP1].th_1)/30
                  + (- fields[iM2][j][k].th_1 + fields[iP2][j][k].th_1)/20
                  - (  fields[iM1][j][kM1].th_1 + fields[iM1][j][kP1].th_1 + fields[iM1][jM1][k].th_1 + fields[iM1][jP1][k].th_1
                     - fields[iP1][j][kM1].th_1 - fields[iP1][j][kP1].th_1 - fields[iP1][jM1][k].th_1 - fields[iP1][jP1][k].th_1)/12
                  + (- fields[iM1][j][k].th_1 + fields[iP1][j][k].th_1)*4.0/15;
#ifdef WRAP_DTH
                  dx_th_1 = wrap(dx_th_1);
#endif
                break;
             case 1: // th_2
                  dx_th_2 = 
                    (  fields[iM1][jM1][kM2].th_2 + fields[iM1][jM1][kP2].th_2 + fields[iM1][jP1][kM2].th_2 + fields[iM1][jP1][kP2].th_2
                     + fields[iM1][jM2][kM1].th_2 + fields[iM1][jM2][kP1].th_2 + fields[iM1][jP2][kM1].th_2 + fields[iM1][jP2][kP1].th_2
                     - fields[iP1][jM1][kM2].th_2 - fields[iP1][jM1][kP2].th_2 - fields[iP1][jP1][kM2].th_2 - fields[iP1][jP1][kP2].th_2
                     - fields[iP1][jM2][kM1].th_2 - fields[iP1][jM2][kP1].th_2 - fields[iP1][jP2][kM1].th_2 - fields[iP1][jP2][kP1].th_2)/120
                  + (  fields[iM2][j][kM1].th_2 + fields[iM2][j][kP1].th_2 + fields[iM2][jM1][k].th_2 + fields[iM2][jP1][k].th_2
                     - fields[iP2][j][kM1].th_2 - fields[iP2][j][kP1].th_2 - fields[iP2][jM1][k].th_2 - fields[iP2][jP1][k].th_2
                     - fields[iM1][jM1][kM1].th_2 - fields[iM1][jM1][kP1].th_2 - fields[iM1][jP1][kM1].th_2 - fields[iM1][jP1][kP1].th_2
                     + fields[iP1][jM1][kM1].th_2 + fields[iP1][jM1][kP1].th_2 + fields[iP1][jP1][kM1].th_2 + fields[iP1][jP1][kP1].th_2)/30
                  - (  fields[iM2][j][k].th_2 - fields[iP2][j][k].th_2)/20
                  + (- fields[iM1][j][kM1].th_2 - fields[iM1][j][kP1].th_2 - fields[iM1][jM1][k].th_2 - fields[iM1][jP1][k].th_2
                     + fields[iP1][j][kM1].th_2 + fields[iP1][j][kP1].th_2 + fields[iP1][jM1][k].th_2 + fields[iP1][jP1][k].th_2)/12
                  - (  fields[iM1][j][k].th_2 - fields[iP1][j][k].th_2)*4.0/15;
#ifdef WRAP_DTH
                  dx_th_2 = wrap(dx_th_2);
#endif
                break;
             case 2: // th_3
                  dx_th_3 = 
                    (  fields[iM1][jM1][kM2].th_3 + fields[iM1][jM1][kP2].th_3 + fields[iM1][jP1][kM2].th_3 + fields[iM1][jP1][kP2].th_3
                     + fields[iM1][jM2][kM1].th_3 + fields[iM1][jM2][kP1].th_3 + fields[iM1][jP2][kM1].th_3 + fields[iM1][jP2][kP1].th_3
                     - fields[iP1][jM1][kM2].th_3 - fields[iP1][jM1][kP2].th_3 - fields[iP1][jP1][kM2].th_3 - fields[iP1][jP1][kP2].th_3
                     - fields[iP1][jM2][kM1].th_3 - fields[iP1][jM2][kP1].th_3 - fields[iP1][jP2][kM1].th_3 - fields[iP1][jP2][kP1].th_3)/120
                  + (  fields[iM2][j][kM1].th_3 + fields[iM2][j][kP1].th_3 + fields[iM2][jM1][k].th_3 + fields[iM2][jP1][k].th_3
                     - fields[iP2][j][kM1].th_3 - fields[iP2][j][kP1].th_3 - fields[iP2][jM1][k].th_3 - fields[iP2][jP1][k].th_3
                     - fields[iM1][jM1][kM1].th_3 - fields[iM1][jM1][kP1].th_3 - fields[iM1][jP1][kM1].th_3 - fields[iM1][jP1][kP1].th_3
                     + fields[iP1][jM1][kM1].th_3 + fields[iP1][jM1][kP1].th_3 + fields[iP1][jP1][kM1].th_3 + fields[iP1][jP1][kP1].th_3)/30
                  + (- fields[iM2][j][k].th_3 + fields[iP2][j][k].th_3)/20
                  + (- fields[iM1][j][kM1].th_3 - fields[iM1][j][kP1].th_3 - fields[iM1][jM1][k].th_3 - fields[iM1][jP1][k].th_3
                     + fields[iP1][j][kM1].th_3 + fields[iP1][j][kP1].th_3 + fields[iP1][jM1][k].th_3 + fields[iP1][jP1][k].th_3)/12
                  + (- fields[iM1][j][k].th_3 + fields[iP1][j][k].th_3)*4.0/15;
#ifdef WRAP_DTH
                  dx_th_3 = wrap(dx_th_3);
#endif
                break;
             case 4: // A_2
                  dx_A_2 = 
                   (  fields[iM1][jM1][kM2].A_2 + fields[iM1][jM1][kP2].A_2 + fields[iM1][jP1][kM2].A_2 + fields[iM1][jP1][kP2].A_2
                    + fields[iM1][jM2][kM1].A_2 + fields[iM1][jM2][kP1].A_2 + fields[iM1][jP2][kM1].A_2 + fields[iM1][jP2][kP1].A_2
                    - fields[iP1][jM1][kM2].A_2 - fields[iP1][jM1][kP2].A_2 - fields[iP1][jP1][kM2].A_2 - fields[iP1][jP1][kP2].A_2
                    - fields[iP1][jM2][kM1].A_2 - fields[iP1][jM2][kP1].A_2 - fields[iP1][jP2][kM1].A_2 - fields[iP1][jP2][kP1].A_2)/120
                 + (  fields[iM2][j][kM1].A_2 + fields[iM2][j][kP1].A_2 + fields[iM2][jM1][k].A_2 + fields[iM2][jP1][k].A_2
                    - fields[iP2][j][kM1].A_2 - fields[iP2][j][kP1].A_2 - fields[iP2][jM1][k].A_2 - fields[iP2][jP1][k].A_2
                    - fields[iM1][jM1][kM1].A_2 - fields[iM1][jM1][kP1].A_2 - fields[iM1][jP1][kM1].A_2 - fields[iM1][jP1][kP1].A_2
                    + fields[iP1][jM1][kM1].A_2 + fields[iP1][jM1][kP1].A_2 + fields[iP1][jP1][kM1].A_2 + fields[iP1][jP1][kP1].A_2)/30
                 + (- fields[iM2][j][k].A_2 + fields[iP2][j][k].A_2)/20
                 + (- fields[iM1][j][kM1].A_2 - fields[iM1][j][kP1].A_2 - fields[iM1][jM1][k].A_2 - fields[iM1][jP1][k].A_2
                    + fields[iP1][j][kM1].A_2 + fields[iP1][j][kP1].A_2 + fields[iP1][jM1][k].A_2 + fields[iP1][jP1][k].A_2)/12
                 + (- fields[iM1][j][k].A_2 + fields[iP1][j][k].A_2)*4.0/15;
                break;
             case 5: // A_3
                  dx_A_3 = 
                   (  fields[iM1][jM1][kM2].A_3 + fields[iM1][jM1][kP2].A_3 + fields[iM1][jP1][kM2].A_3 + fields[iM1][jP1][kP2].A_3
                    + fields[iM1][jM2][kM1].A_3 + fields[iM1][jM2][kP1].A_3 + fields[iM1][jP2][kM1].A_3 + fields[iM1][jP2][kP1].A_3
                    - fields[iP1][jM1][kM2].A_3 - fields[iP1][jM1][kP2].A_3 - fields[iP1][jP1][kM2].A_3 - fields[iP1][jP1][kP2].A_3
                    - fields[iP1][jM2][kM1].A_3 - fields[iP1][jM2][kP1].A_3 - fields[iP1][jP2][kM1].A_3 - fields[iP1][jP2][kP1].A_3)/120
                 + (  fields[iM2][j][kM1].A_3 + fields[iM2][j][kP1].A_3 + fields[iM2][jM1][k].A_3 + fields[iM2][jP1][k].A_3
                    - fields[iP2][j][kM1].A_3 - fields[iP2][j][kP1].A_3 - fields[iP2][jM1][k].A_3 - fields[iP2][jP1][k].A_3
                    - fields[iM1][jM1][kM1].A_3 - fields[iM1][jM1][kP1].A_3 - fields[iM1][jP1][kM1].A_3 - fields[iM1][jP1][kP1].A_3
                    + fields[iP1][jM1][kM1].A_3 + fields[iP1][jM1][kP1].A_3 + fields[iP1][jP1][kM1].A_3 + fields[iP1][jP1][kP1].A_3)/30
                 + (- fields[iM2][j][k].A_3 + fields[iP2][j][k].A_3)/20
                 - (  fields[iM1][j][kM1].A_3 + fields[iM1][j][kP1].A_3 + fields[iM1][jM1][k].A_3 + fields[iM1][jP1][k].A_3
                    - fields[iP1][j][kM1].A_3 - fields[iP1][j][kP1].A_3 - fields[iP1][jM1][k].A_3 - fields[iP1][jP1][k].A_3)/12
                 + (- fields[iM1][j][k].A_3 + fields[iP1][j][k].A_3)*4.0/15;
               break;
         }
         break;
         
      case 1: // y derivative

         switch (component)
         {
             case 0: // th_1
                  dy_th_1 = 
                    (  fields[iM1][jM1][kM2].th_1 + fields[iM1][jM1][kP2].th_1 - fields[iM1][jP1][kM2].th_1 - fields[iM1][jP1][kP2].th_1
                     + fields[iP1][jM1][kM2].th_1 + fields[iP1][jM1][kP2].th_1 - fields[iP1][jP1][kM2].th_1 - fields[iP1][jP1][kP2].th_1
                     + fields[iM2][jM1][kM1].th_1 + fields[iM2][jM1][kP1].th_1 - fields[iM2][jP1][kM1].th_1 - fields[iM2][jP1][kP1].th_1
                     + fields[iP2][jM1][kM1].th_1 + fields[iP2][jM1][kP1].th_1 - fields[iP2][jP1][kM1].th_1 - fields[iP2][jP1][kP1].th_1)/120
                  + (  fields[i][jM2][kM1].th_1 + fields[i][jM2][kP1].th_1 - fields[i][jP2][kM1].th_1 - fields[i][jP2][kP1].th_1
                     + fields[iM1][jM2][k].th_1 - fields[iM1][jP2][k].th_1 + fields[iP1][jM2][k].th_1 - fields[iP1][jP2][k].th_1
                     - fields[iM1][jM1][kM1].th_1 - fields[iM1][jM1][kP1].th_1 + fields[iM1][jP1][kM1].th_1 + fields[iM1][jP1][kP1].th_1
                     - fields[iP1][jM1][kM1].th_1 - fields[iP1][jM1][kP1].th_1 + fields[iP1][jP1][kM1].th_1 + fields[iP1][jP1][kP1].th_1)/30
                  + (- fields[i][jM2][k].th_1 + fields[i][jP2][k].th_1)/20
                  - (  fields[i][jM1][kM1].th_1 + fields[i][jM1][kP1].th_1 - fields[i][jP1][kM1].th_1 - fields[i][jP1][kP1].th_1
                     + fields[iM1][jM1][k].th_1 - fields[iM1][jP1][k].th_1 + fields[iP1][jM1][k].th_1 - fields[iP1][jP1][k].th_1)/12
                  + (- fields[i][jM1][k].th_1 + fields[i][jP1][k].th_1)*4.0/15;
#ifdef WRAP_DTH
                  dy_th_1 = wrap(dy_th_1);
#endif
                break;
             case 1: // th_2
                  dy_th_2 = 
                    (  fields[iM1][jM1][kM2].th_2 + fields[iM1][jM1][kP2].th_2 - fields[iM1][jP1][kM2].th_2 - fields[iM1][jP1][kP2].th_2
                     + fields[iP1][jM1][kM2].th_2 + fields[iP1][jM1][kP2].th_2 - fields[iP1][jP1][kM2].th_2 - fields[iP1][jP1][kP2].th_2
                     + fields[iM2][jM1][kM1].th_2 + fields[iM2][jM1][kP1].th_2 - fields[iM2][jP1][kM1].th_2 - fields[iM2][jP1][kP1].th_2
                     + fields[iP2][jM1][kM1].th_2 + fields[iP2][jM1][kP1].th_2 - fields[iP2][jP1][kM1].th_2 - fields[iP2][jP1][kP1].th_2)/120
                  + (  fields[i][jM2][kM1].th_2 + fields[i][jM2][kP1].th_2 - fields[i][jP2][kM1].th_2 - fields[i][jP2][kP1].th_2
                     + fields[iM1][jM2][k].th_2 - fields[iM1][jP2][k].th_2 + fields[iP1][jM2][k].th_2 - fields[iP1][jP2][k].th_2
                     - fields[iM1][jM1][kM1].th_2 - fields[iM1][jM1][kP1].th_2 + fields[iM1][jP1][kM1].th_2 + fields[iM1][jP1][kP1].th_2
                     - fields[iP1][jM1][kM1].th_2 - fields[iP1][jM1][kP1].th_2 + fields[iP1][jP1][kM1].th_2 + fields[iP1][jP1][kP1].th_2)/30
                  + (- fields[i][jM2][k].th_2 + fields[i][jP2][k].th_2)/20
                  + (- fields[i][jM1][kM1].th_2 - fields[i][jM1][kP1].th_2 + fields[i][jP1][kM1].th_2 + fields[i][jP1][kP1].th_2
                     - fields[iM1][jM1][k].th_2 + fields[iM1][jP1][k].th_2 - fields[iP1][jM1][k].th_2 + fields[iP1][jP1][k].th_2)/12
                  + (- fields[i][jM1][k].th_2 + fields[i][jP1][k].th_2)*4.0/15;
#ifdef WRAP_DTH
                  dy_th_2 = wrap(dy_th_2);
#endif
                break;
             case 2: // th_3
                  dy_th_3 = 
                    (  fields[iM1][jM1][kM2].th_3 + fields[iM1][jM1][kP2].th_3 - fields[iM1][jP1][kM2].th_3 - fields[iM1][jP1][kP2].th_3
                     + fields[iP1][jM1][kM2].th_3 + fields[iP1][jM1][kP2].th_3 - fields[iP1][jP1][kM2].th_3 - fields[iP1][jP1][kP2].th_3
                     + fields[iM2][jM1][kM1].th_3 + fields[iM2][jM1][kP1].th_3 - fields[iM2][jP1][kM1].th_3 - fields[iM2][jP1][kP1].th_3
                     + fields[iP2][jM1][kM1].th_3 + fields[iP2][jM1][kP1].th_3 - fields[iP2][jP1][kM1].th_3 - fields[iP2][jP1][kP1].th_3)/120
                  + (  fields[i][jM2][kM1].th_3 + fields[i][jM2][kP1].th_3 - fields[i][jP2][kM1].th_3 - fields[i][jP2][kP1].th_3
                     + fields[iM1][jM2][k].th_3 - fields[iM1][jP2][k].th_3 + fields[iP1][jM2][k].th_3 - fields[iP1][jP2][k].th_3
                     - fields[iM1][jM1][kM1].th_3 - fields[iM1][jM1][kP1].th_3 + fields[iM1][jP1][kM1].th_3 + fields[iM1][jP1][kP1].th_3
                     - fields[iP1][jM1][kM1].th_3 - fields[iP1][jM1][kP1].th_3 + fields[iP1][jP1][kM1].th_3 + fields[iP1][jP1][kP1].th_3)/30
                  - (  fields[i][jM2][k].th_3 - fields[i][jP2][k].th_3)/20
                  + (- fields[i][jM1][kM1].th_3 - fields[i][jM1][kP1].th_3 + fields[i][jP1][kM1].th_3 + fields[i][jP1][kP1].th_3
                     - fields[iM1][jM1][k].th_3 + fields[iM1][jP1][k].th_3 - fields[iP1][jM1][k].th_3 + fields[iP1][jP1][k].th_3)/12
                  - (  fields[i][jM1][k].th_3 - fields[i][jP1][k].th_3)*4.0/15;
#ifdef WRAP_DTH
                  dy_th_3 = wrap(dy_th_3);
#endif
                break;
             case 3: // A_1
                 dy_A_1 = 
                   (  fields[iM1][jM1][kM2].A_1 + fields[iM1][jM1][kP2].A_1 - fields[iM1][jP1][kM2].A_1 - fields[iM1][jP1][kP2].A_1
                    + fields[iP1][jM1][kM2].A_1 + fields[iP1][jM1][kP2].A_1 - fields[iP1][jP1][kM2].A_1 - fields[iP1][jP1][kP2].A_1
                    + fields[iM2][jM1][kM1].A_1 + fields[iM2][jM1][kP1].A_1 - fields[iM2][jP1][kM1].A_1 - fields[iM2][jP1][kP1].A_1
                    + fields[iP2][jM1][kM1].A_1 + fields[iP2][jM1][kP1].A_1 - fields[iP2][jP1][kM1].A_1 - fields[iP2][jP1][kP1].A_1)/120
                 + (  fields[i][jM2][kM1].A_1 + fields[i][jM2][kP1].A_1 - fields[i][jP2][kM1].A_1 - fields[i][jP2][kP1].A_1
                    + fields[iM1][jM2][k].A_1 - fields[iM1][jP2][k].A_1 + fields[iP1][jM2][k].A_1 - fields[iP1][jP2][k].A_1
                    - fields[iM1][jM1][kM1].A_1 - fields[iM1][jM1][kP1].A_1 + fields[iM1][jP1][kM1].A_1 + fields[iM1][jP1][kP1].A_1
                    - fields[iP1][jM1][kM1].A_1 - fields[iP1][jM1][kP1].A_1 + fields[iP1][jP1][kM1].A_1 + fields[iP1][jP1][kP1].A_1)/30
                 + (- fields[i][jM2][k].A_1 + fields[i][jP2][k].A_1)/20
                 - (  fields[i][jM1][kM1].A_1 + fields[i][jM1][kP1].A_1 - fields[i][jP1][kM1].A_1 - fields[i][jP1][kP1].A_1
                    + fields[iM1][jM1][k].A_1 - fields[iM1][jP1][k].A_1 + fields[iP1][jM1][k].A_1 - fields[iP1][jP1][k].A_1)/12
                 + (- fields[i][jM1][k].A_1 + fields[i][jP1][k].A_1)*4.0/15;
                break;
             case 5: // A_3
                 dy_A_3 = 
                   (  fields[iM1][jM1][kM2].A_3 + fields[iM1][jM1][kP2].A_3 - fields[iM1][jP1][kM2].A_3 - fields[iM1][jP1][kP2].A_3
                    + fields[iP1][jM1][kM2].A_3 + fields[iP1][jM1][kP2].A_3 - fields[iP1][jP1][kM2].A_3 - fields[iP1][jP1][kP2].A_3
                    + fields[iM2][jM1][kM1].A_3 + fields[iM2][jM1][kP1].A_3 - fields[iM2][jP1][kM1].A_3 - fields[iM2][jP1][kP1].A_3
                    + fields[iP2][jM1][kM1].A_3 + fields[iP2][jM1][kP1].A_3 - fields[iP2][jP1][kM1].A_3 - fields[iP2][jP1][kP1].A_3)/120
                 + (  fields[i][jM2][kM1].A_3 + fields[i][jM2][kP1].A_3 - fields[i][jP2][kM1].A_3 - fields[i][jP2][kP1].A_3
                    + fields[iM1][jM2][k].A_3 - fields[iM1][jP2][k].A_3 + fields[iP1][jM2][k].A_3 - fields[iP1][jP2][k].A_3
                    - fields[iM1][jM1][kM1].A_3 - fields[iM1][jM1][kP1].A_3 + fields[iM1][jP1][kM1].A_3 + fields[iM1][jP1][kP1].A_3
                    - fields[iP1][jM1][kM1].A_3 - fields[iP1][jM1][kP1].A_3 + fields[iP1][jP1][kM1].A_3 + fields[iP1][jP1][kP1].A_3)/30
                 - (  fields[i][jM2][k].A_3 - fields[i][jP2][k].A_3)/20
                 + (- fields[i][jM1][kM1].A_3 - fields[i][jM1][kP1].A_3 + fields[i][jP1][kM1].A_3 + fields[i][jP1][kP1].A_3
                    - fields[iM1][jM1][k].A_3 + fields[iM1][jP1][k].A_3 - fields[iP1][jM1][k].A_3 + fields[iP1][jP1][k].A_3)/12
                 - (  fields[i][jM1][k].A_3 - fields[i][jP1][k].A_3)*4.0/15;
      break;
         }
         break;
         
      case 2: // z derivative

         switch (component)
         {
             case 0: // th_1
                  dz_th_1 = 
                    (  fields[iM1][jM2][kM1].th_1 - fields[iM1][jM2][kP1].th_1 + fields[iM1][jP2][kM1].th_1 - fields[iM1][jP2][kP1].th_1
                     + fields[iP1][jM2][kM1].th_1 - fields[iP1][jM2][kP1].th_1 + fields[iP1][jP2][kM1].th_1 - fields[iP1][jP2][kP1].th_1
                     + fields[iM2][jM1][kM1].th_1 - fields[iM2][jM1][kP1].th_1 + fields[iM2][jP1][kM1].th_1 - fields[iM2][jP1][kP1].th_1
                     + fields[iP2][jM1][kM1].th_1 - fields[iP2][jM1][kP1].th_1 + fields[iP2][jP1][kM1].th_1 - fields[iP2][jP1][kP1].th_1)/120
                  + (  fields[i][jM1][kM2].th_1 - fields[i][jM1][kP2].th_1 + fields[i][jP1][kM2].th_1 - fields[i][jP1][kP2].th_1
                     + fields[iM1][j][kM2].th_1 - fields[iM1][j][kP2].th_1 + fields[iP1][j][kM2].th_1 - fields[iP1][j][kP2].th_1
                     - fields[iM1][jM1][kM1].th_1 + fields[iM1][jM1][kP1].th_1 - fields[iM1][jP1][kM1].th_1 + fields[iM1][jP1][kP1].th_1
                     - fields[iP1][jM1][kM1].th_1 + fields[iP1][jM1][kP1].th_1 - fields[iP1][jP1][kM1].th_1 + fields[iP1][jP1][kP1].th_1)/30
                  + (- fields[i][j][kM2].th_1 + fields[i][j][kP2].th_1)/20
                  - (  fields[i][jM1][kM1].th_1 - fields[i][jM1][kP1].th_1 + fields[i][jP1][kM1].th_1 - fields[i][jP1][kP1].th_1
                     + fields[iM1][j][kM1].th_1 - fields[iM1][j][kP1].th_1 + fields[iP1][j][kM1].th_1 - fields[iP1][j][kP1].th_1)/12
                  + (- fields[i][j][kM1].th_1 + fields[i][j][kP1].th_1)*4.0/15;
#ifdef WRAP_DTH
                  dz_th_1 = wrap(dz_th_1);
#endif
                break;
             case 1: // th_2
                  dz_th_2 = 
                    (  fields[iM1][jM2][kM1].th_2 - fields[iM1][jM2][kP1].th_2 + fields[iM1][jP2][kM1].th_2 - fields[iM1][jP2][kP1].th_2
                     + fields[iP1][jM2][kM1].th_2 - fields[iP1][jM2][kP1].th_2 + fields[iP1][jP2][kM1].th_2 - fields[iP1][jP2][kP1].th_2
                     + fields[iM2][jM1][kM1].th_2 - fields[iM2][jM1][kP1].th_2 + fields[iM2][jP1][kM1].th_2 - fields[iM2][jP1][kP1].th_2
                     + fields[iP2][jM1][kM1].th_2 - fields[iP2][jM1][kP1].th_2 + fields[iP2][jP1][kM1].th_2 - fields[iP2][jP1][kP1].th_2)/120
                  + (  fields[i][jM1][kM2].th_2 - fields[i][jM1][kP2].th_2 + fields[i][jP1][kM2].th_2 - fields[i][jP1][kP2].th_2
                     + fields[iM1][j][kM2].th_2 - fields[iM1][j][kP2].th_2 + fields[iP1][j][kM2].th_2 - fields[iP1][j][kP2].th_2
                     - fields[iM1][jM1][kM1].th_2 + fields[iM1][jM1][kP1].th_2 - fields[iM1][jP1][kM1].th_2 + fields[iM1][jP1][kP1].th_2
                     - fields[iP1][jM1][kM1].th_2 + fields[iP1][jM1][kP1].th_2 - fields[iP1][jP1][kM1].th_2 + fields[iP1][jP1][kP1].th_2)/30
                  + (- fields[i][j][kM2].th_2 + fields[i][j][kP2].th_2)/20
                  - (  fields[i][jM1][kM1].th_2 - fields[i][jM1][kP1].th_2 + fields[i][jP1][kM1].th_2 - fields[i][jP1][kP1].th_2
                     + fields[iM1][j][kM1].th_2 - fields[iM1][j][kP1].th_2 + fields[iP1][j][kM1].th_2 - fields[iP1][j][kP1].th_2)/12
                  + (- fields[i][j][kM1].th_2 + fields[i][j][kP1].th_2)*4.0/15;
#ifdef WRAP_DTH
                  dz_th_2 = wrap(dz_th_2);
#endif
                break;
             case 2: // th_3
                  dz_th_3 = 
                    (  fields[iM1][jM2][kM1].th_3 - fields[iM1][jM2][kP1].th_3 + fields[iM1][jP2][kM1].th_3 - fields[iM1][jP2][kP1].th_3
                     + fields[iP1][jM2][kM1].th_3 - fields[iP1][jM2][kP1].th_3 + fields[iP1][jP2][kM1].th_3 - fields[iP1][jP2][kP1].th_3
                     + fields[iM2][jM1][kM1].th_3 - fields[iM2][jM1][kP1].th_3 + fields[iM2][jP1][kM1].th_3 - fields[iM2][jP1][kP1].th_3
                     + fields[iP2][jM1][kM1].th_3 - fields[iP2][jM1][kP1].th_3 + fields[iP2][jP1][kM1].th_3 - fields[iP2][jP1][kP1].th_3)/120
                  + (  fields[i][jM1][kM2].th_3 - fields[i][jM1][kP2].th_3 + fields[i][jP1][kM2].th_3 - fields[i][jP1][kP2].th_3
                     + fields[iM1][j][kM2].th_3 - fields[iM1][j][kP2].th_3 + fields[iP1][j][kM2].th_3 - fields[iP1][j][kP2].th_3
                     - fields[iM1][jM1][kM1].th_3 + fields[iM1][jM1][kP1].th_3 - fields[iM1][jP1][kM1].th_3 + fields[iM1][jP1][kP1].th_3
                     - fields[iP1][jM1][kM1].th_3 + fields[iP1][jM1][kP1].th_3 - fields[iP1][jP1][kM1].th_3 + fields[iP1][jP1][kP1].th_3)/30
                  + (- fields[i][j][kM2].th_3 + fields[i][j][kP2].th_3)/20
                  - (  fields[i][jM1][kM1].th_3 - fields[i][jM1][kP1].th_3 + fields[i][jP1][kM1].th_3 - fields[i][jP1][kP1].th_3
                     + fields[iM1][j][kM1].th_3 - fields[iM1][j][kP1].th_3 + fields[iP1][j][kM1].th_3 - fields[iP1][j][kP1].th_3)/12
                  + (- fields[i][j][kM1].th_3 + fields[i][j][kP1].th_3)*4.0/15;
#ifdef WRAP_DTH
                  dz_th_3 = wrap(dz_th_3);
#endif
                break;
             case 3: // A_1
                 dz_A_1 =
                   (  fields[iM1][jM2][kM1].A_1 - fields[iM1][jM2][kP1].A_1 + fields[iM1][jP2][kM1].A_1 - fields[iM1][jP2][kP1].A_1
                    + fields[iP1][jM2][kM1].A_1 - fields[iP1][jM2][kP1].A_1 + fields[iP1][jP2][kM1].A_1 - fields[iP1][jP2][kP1].A_1
                    + fields[iM2][jM1][kM1].A_1 - fields[iM2][jM1][kP1].A_1 + fields[iM2][jP1][kM1].A_1 - fields[iM2][jP1][kP1].A_1
                    + fields[iP2][jM1][kM1].A_1 - fields[iP2][jM1][kP1].A_1 + fields[iP2][jP1][kM1].A_1 - fields[iP2][jP1][kP1].A_1)/120
                 + (  fields[i][jM1][kM2].A_1 - fields[i][jM1][kP2].A_1 + fields[i][jP1][kM2].A_1 - fields[i][jP1][kP2].A_1
                    + fields[iM1][j][kM2].A_1 - fields[iM1][j][kP2].A_1 + fields[iP1][j][kM2].A_1 - fields[iP1][j][kP2].A_1
                    - fields[iM1][jM1][kM1].A_1 + fields[iM1][jM1][kP1].A_1 - fields[iM1][jP1][kM1].A_1 + fields[iM1][jP1][kP1].A_1
                    - fields[iP1][jM1][kM1].A_1 + fields[iP1][jM1][kP1].A_1 - fields[iP1][jP1][kM1].A_1 + fields[iP1][jP1][kP1].A_1)/30
                 + (- fields[i][j][kM2].A_1 + fields[i][j][kP2].A_1)/20
                 - (  fields[i][jM1][kM1].A_1 - fields[i][jM1][kP1].A_1 + fields[i][jP1][kM1].A_1 - fields[i][jP1][kP1].A_1
                    + fields[iM1][j][kM1].A_1 - fields[iM1][j][kP1].A_1 + fields[iP1][j][kM1].A_1 - fields[iP1][j][kP1].A_1)/12
                 + (- fields[i][j][kM1].A_1 + fields[i][j][kP1].A_1)*4.0/15;
                break;
             case 4: // A_2
                 dz_A_2 =
                   (  fields[iM1][jM2][kM1].A_2 - fields[iM1][jM2][kP1].A_2 + fields[iM1][jP2][kM1].A_2 - fields[iM1][jP2][kP1].A_2
                    + fields[iP1][jM2][kM1].A_2 - fields[iP1][jM2][kP1].A_2 + fields[iP1][jP2][kM1].A_2 - fields[iP1][jP2][kP1].A_2
                    + fields[iM2][jM1][kM1].A_2 - fields[iM2][jM1][kP1].A_2 + fields[iM2][jP1][kM1].A_2 - fields[iM2][jP1][kP1].A_2
                    + fields[iP2][jM1][kM1].A_2 - fields[iP2][jM1][kP1].A_2 + fields[iP2][jP1][kM1].A_2 - fields[iP2][jP1][kP1].A_2)/120
                 + (  fields[i][jM1][kM2].A_2 - fields[i][jM1][kP2].A_2 + fields[i][jP1][kM2].A_2 - fields[i][jP1][kP2].A_2
                    + fields[iM1][j][kM2].A_2 - fields[iM1][j][kP2].A_2 + fields[iP1][j][kM2].A_2 - fields[iP1][j][kP2].A_2
                    - fields[iM1][jM1][kM1].A_2 + fields[iM1][jM1][kP1].A_2 - fields[iM1][jP1][kM1].A_2 + fields[iM1][jP1][kP1].A_2
                    - fields[iP1][jM1][kM1].A_2 + fields[iP1][jM1][kP1].A_2 - fields[iP1][jP1][kM1].A_2 + fields[iP1][jP1][kP1].A_2)/30
                 + (- fields[i][j][kM2].A_2 + fields[i][j][kP2].A_2)/20
                 - (  fields[i][jM1][kM1].A_2 - fields[i][jM1][kP1].A_2 + fields[i][jP1][kM1].A_2 - fields[i][jP1][kP1].A_2
                    + fields[iM1][j][kM1].A_2 - fields[iM1][j][kP1].A_2 + fields[iP1][j][kM1].A_2 - fields[iP1][j][kP1].A_2)/12
                 + (- fields[i][j][kM1].A_2 + fields[i][j][kP1].A_2)*4.0/15;
                break;
         }
         break;
   }
}

void TAux::compute_derivatives(TFields fields[MAX_NSIDE][MAX_NSIDE][MAX_NSIDE], int i, int j, int k, int nSideM1)
/* IMPORTANT: Does NOT divide by lattice spacing, so if absolute (rather than relative) magnitude of derivative is needed, user must do that. */
{
   int iM1, jM1, kM1, iM2, jM2, kM2, iP1, jP1, kP1, iP2, jP2, kP2;

   if (0 == i) iM1 = nSideM1; else iM1 = i - 1;
   if (0 == j) jM1 = nSideM1; else jM1 = j - 1;
   if (0 == k) kM1 = nSideM1; else kM1 = k - 1;

   if (0 == iM1) iM2 = nSideM1; else iM2 = iM1 - 1;
   if (0 == jM1) jM2 = nSideM1; else jM2 = jM1 - 1;
   if (0 == kM1) kM2 = nSideM1; else kM2 = kM1 - 1;

   if (nSideM1 == i) iP1 = 0; else iP1 = i + 1;
   if (nSideM1 == j) jP1 = 0; else jP1 = j + 1;
   if (nSideM1 == k) kP1 = 0; else kP1 = k + 1;

   if (nSideM1 == iP1) iP2 = 0; else iP2 = iP1 + 1;
   if (nSideM1 == jP1) jP2 = 0; else jP2 = jP1 + 1;
   if (nSideM1 == kP1) kP2 = 0; else kP2 = kP1 + 1;

   dx_th_1 = 
                    (  fields[iM1][jM1][kM2].th_1 + fields[iM1][jM1][kP2].th_1 + fields[iM1][jP1][kM2].th_1 + fields[iM1][jP1][kP2].th_1
                     + fields[iM1][jM2][kM1].th_1 + fields[iM1][jM2][kP1].th_1 + fields[iM1][jP2][kM1].th_1 + fields[iM1][jP2][kP1].th_1
                     - fields[iP1][jM1][kM2].th_1 - fields[iP1][jM1][kP2].th_1 - fields[iP1][jP1][kM2].th_1 - fields[iP1][jP1][kP2].th_1
                     - fields[iP1][jM2][kM1].th_1 - fields[iP1][jM2][kP1].th_1 - fields[iP1][jP2][kM1].th_1 - fields[iP1][jP2][kP1].th_1)/120
                  + (  fields[iM2][j][kM1].th_1 + fields[iM2][j][kP1].th_1 + fields[iM2][jM1][k].th_1 + fields[iM2][jP1][k].th_1
                     - fields[iP2][j][kM1].th_1 - fields[iP2][j][kP1].th_1 - fields[iP2][jM1][k].th_1 - fields[iP2][jP1][k].th_1
                     - fields[iM1][jM1][kM1].th_1 - fields[iM1][jM1][kP1].th_1 - fields[iM1][jP1][kM1].th_1 - fields[iM1][jP1][kP1].th_1
                     + fields[iP1][jM1][kM1].th_1 + fields[iP1][jM1][kP1].th_1 + fields[iP1][jP1][kM1].th_1 + fields[iP1][jP1][kP1].th_1)/30
                  + (- fields[iM2][j][k].th_1 + fields[iP2][j][k].th_1)/20
                  - (  fields[iM1][j][kM1].th_1 + fields[iM1][j][kP1].th_1 + fields[iM1][jM1][k].th_1 + fields[iM1][jP1][k].th_1
                     - fields[iP1][j][kM1].th_1 - fields[iP1][j][kP1].th_1 - fields[iP1][jM1][k].th_1 - fields[iP1][jP1][k].th_1)/12
                  + (- fields[iM1][j][k].th_1 + fields[iP1][j][k].th_1)*4.0/15;
#ifdef WRAP_DTH
                  dx_th_1 = wrap(dx_th_1);
#endif
   dx_th_2 = 
                    (  fields[iM1][jM1][kM2].th_2 + fields[iM1][jM1][kP2].th_2 + fields[iM1][jP1][kM2].th_2 + fields[iM1][jP1][kP2].th_2
                     + fields[iM1][jM2][kM1].th_2 + fields[iM1][jM2][kP1].th_2 + fields[iM1][jP2][kM1].th_2 + fields[iM1][jP2][kP1].th_2
                     - fields[iP1][jM1][kM2].th_2 - fields[iP1][jM1][kP2].th_2 - fields[iP1][jP1][kM2].th_2 - fields[iP1][jP1][kP2].th_2
                     - fields[iP1][jM2][kM1].th_2 - fields[iP1][jM2][kP1].th_2 - fields[iP1][jP2][kM1].th_2 - fields[iP1][jP2][kP1].th_2)/120
                  + (  fields[iM2][j][kM1].th_2 + fields[iM2][j][kP1].th_2 + fields[iM2][jM1][k].th_2 + fields[iM2][jP1][k].th_2
                     - fields[iP2][j][kM1].th_2 - fields[iP2][j][kP1].th_2 - fields[iP2][jM1][k].th_2 - fields[iP2][jP1][k].th_2
                     - fields[iM1][jM1][kM1].th_2 - fields[iM1][jM1][kP1].th_2 - fields[iM1][jP1][kM1].th_2 - fields[iM1][jP1][kP1].th_2
                     + fields[iP1][jM1][kM1].th_2 + fields[iP1][jM1][kP1].th_2 + fields[iP1][jP1][kM1].th_2 + fields[iP1][jP1][kP1].th_2)/30
                  - (  fields[iM2][j][k].th_2 - fields[iP2][j][k].th_2)/20
                  + (- fields[iM1][j][kM1].th_2 - fields[iM1][j][kP1].th_2 - fields[iM1][jM1][k].th_2 - fields[iM1][jP1][k].th_2
                     + fields[iP1][j][kM1].th_2 + fields[iP1][j][kP1].th_2 + fields[iP1][jM1][k].th_2 + fields[iP1][jP1][k].th_2)/12
                  - (  fields[iM1][j][k].th_2 - fields[iP1][j][k].th_2)*4.0/15;
#ifdef WRAP_DTH
                  dx_th_2 = wrap(dx_th_2);
#endif
   dx_th_3 = 
                    (  fields[iM1][jM1][kM2].th_3 + fields[iM1][jM1][kP2].th_3 + fields[iM1][jP1][kM2].th_3 + fields[iM1][jP1][kP2].th_3
                     + fields[iM1][jM2][kM1].th_3 + fields[iM1][jM2][kP1].th_3 + fields[iM1][jP2][kM1].th_3 + fields[iM1][jP2][kP1].th_3
                     - fields[iP1][jM1][kM2].th_3 - fields[iP1][jM1][kP2].th_3 - fields[iP1][jP1][kM2].th_3 - fields[iP1][jP1][kP2].th_3
                     - fields[iP1][jM2][kM1].th_3 - fields[iP1][jM2][kP1].th_3 - fields[iP1][jP2][kM1].th_3 - fields[iP1][jP2][kP1].th_3)/120
                  + (  fields[iM2][j][kM1].th_3 + fields[iM2][j][kP1].th_3 + fields[iM2][jM1][k].th_3 + fields[iM2][jP1][k].th_3
                     - fields[iP2][j][kM1].th_3 - fields[iP2][j][kP1].th_3 - fields[iP2][jM1][k].th_3 - fields[iP2][jP1][k].th_3
                     - fields[iM1][jM1][kM1].th_3 - fields[iM1][jM1][kP1].th_3 - fields[iM1][jP1][kM1].th_3 - fields[iM1][jP1][kP1].th_3
                     + fields[iP1][jM1][kM1].th_3 + fields[iP1][jM1][kP1].th_3 + fields[iP1][jP1][kM1].th_3 + fields[iP1][jP1][kP1].th_3)/30
                  + (- fields[iM2][j][k].th_3 + fields[iP2][j][k].th_3)/20
                  + (- fields[iM1][j][kM1].th_3 - fields[iM1][j][kP1].th_3 - fields[iM1][jM1][k].th_3 - fields[iM1][jP1][k].th_3
                     + fields[iP1][j][kM1].th_3 + fields[iP1][j][kP1].th_3 + fields[iP1][jM1][k].th_3 + fields[iP1][jP1][k].th_3)/12
                  + (- fields[iM1][j][k].th_3 + fields[iP1][j][k].th_3)*4.0/15;
#ifdef WRAP_DTH
                  dx_th_3 = wrap(dx_th_3);
#endif
   dx_A_2 = 
                   (  fields[iM1][jM1][kM2].A_2 + fields[iM1][jM1][kP2].A_2 + fields[iM1][jP1][kM2].A_2 + fields[iM1][jP1][kP2].A_2
                    + fields[iM1][jM2][kM1].A_2 + fields[iM1][jM2][kP1].A_2 + fields[iM1][jP2][kM1].A_2 + fields[iM1][jP2][kP1].A_2
                    - fields[iP1][jM1][kM2].A_2 - fields[iP1][jM1][kP2].A_2 - fields[iP1][jP1][kM2].A_2 - fields[iP1][jP1][kP2].A_2
                    - fields[iP1][jM2][kM1].A_2 - fields[iP1][jM2][kP1].A_2 - fields[iP1][jP2][kM1].A_2 - fields[iP1][jP2][kP1].A_2)/120
                 + (  fields[iM2][j][kM1].A_2 + fields[iM2][j][kP1].A_2 + fields[iM2][jM1][k].A_2 + fields[iM2][jP1][k].A_2
                    - fields[iP2][j][kM1].A_2 - fields[iP2][j][kP1].A_2 - fields[iP2][jM1][k].A_2 - fields[iP2][jP1][k].A_2
                    - fields[iM1][jM1][kM1].A_2 - fields[iM1][jM1][kP1].A_2 - fields[iM1][jP1][kM1].A_2 - fields[iM1][jP1][kP1].A_2
                    + fields[iP1][jM1][kM1].A_2 + fields[iP1][jM1][kP1].A_2 + fields[iP1][jP1][kM1].A_2 + fields[iP1][jP1][kP1].A_2)/30
                 + (- fields[iM2][j][k].A_2 + fields[iP2][j][k].A_2)/20
                 + (- fields[iM1][j][kM1].A_2 - fields[iM1][j][kP1].A_2 - fields[iM1][jM1][k].A_2 - fields[iM1][jP1][k].A_2
                    + fields[iP1][j][kM1].A_2 + fields[iP1][j][kP1].A_2 + fields[iP1][jM1][k].A_2 + fields[iP1][jP1][k].A_2)/12
                 + (- fields[iM1][j][k].A_2 + fields[iP1][j][k].A_2)*4.0/15;
   dx_A_3 = 
                   (  fields[iM1][jM1][kM2].A_3 + fields[iM1][jM1][kP2].A_3 + fields[iM1][jP1][kM2].A_3 + fields[iM1][jP1][kP2].A_3
                    + fields[iM1][jM2][kM1].A_3 + fields[iM1][jM2][kP1].A_3 + fields[iM1][jP2][kM1].A_3 + fields[iM1][jP2][kP1].A_3
                    - fields[iP1][jM1][kM2].A_3 - fields[iP1][jM1][kP2].A_3 - fields[iP1][jP1][kM2].A_3 - fields[iP1][jP1][kP2].A_3
                    - fields[iP1][jM2][kM1].A_3 - fields[iP1][jM2][kP1].A_3 - fields[iP1][jP2][kM1].A_3 - fields[iP1][jP2][kP1].A_3)/120
                 + (  fields[iM2][j][kM1].A_3 + fields[iM2][j][kP1].A_3 + fields[iM2][jM1][k].A_3 + fields[iM2][jP1][k].A_3
                    - fields[iP2][j][kM1].A_3 - fields[iP2][j][kP1].A_3 - fields[iP2][jM1][k].A_3 - fields[iP2][jP1][k].A_3
                    - fields[iM1][jM1][kM1].A_3 - fields[iM1][jM1][kP1].A_3 - fields[iM1][jP1][kM1].A_3 - fields[iM1][jP1][kP1].A_3
                    + fields[iP1][jM1][kM1].A_3 + fields[iP1][jM1][kP1].A_3 + fields[iP1][jP1][kM1].A_3 + fields[iP1][jP1][kP1].A_3)/30
                 + (- fields[iM2][j][k].A_3 + fields[iP2][j][k].A_3)/20
                 - (  fields[iM1][j][kM1].A_3 + fields[iM1][j][kP1].A_3 + fields[iM1][jM1][k].A_3 + fields[iM1][jP1][k].A_3
                    - fields[iP1][j][kM1].A_3 - fields[iP1][j][kP1].A_3 - fields[iP1][jM1][k].A_3 - fields[iP1][jP1][k].A_3)/12
                 + (- fields[iM1][j][k].A_3 + fields[iP1][j][k].A_3)*4.0/15;
   dy_th_1 = 
                    (  fields[iM1][jM1][kM2].th_1 + fields[iM1][jM1][kP2].th_1 - fields[iM1][jP1][kM2].th_1 - fields[iM1][jP1][kP2].th_1
                     + fields[iP1][jM1][kM2].th_1 + fields[iP1][jM1][kP2].th_1 - fields[iP1][jP1][kM2].th_1 - fields[iP1][jP1][kP2].th_1
                     + fields[iM2][jM1][kM1].th_1 + fields[iM2][jM1][kP1].th_1 - fields[iM2][jP1][kM1].th_1 - fields[iM2][jP1][kP1].th_1
                     + fields[iP2][jM1][kM1].th_1 + fields[iP2][jM1][kP1].th_1 - fields[iP2][jP1][kM1].th_1 - fields[iP2][jP1][kP1].th_1)/120
                  + (  fields[i][jM2][kM1].th_1 + fields[i][jM2][kP1].th_1 - fields[i][jP2][kM1].th_1 - fields[i][jP2][kP1].th_1
                     + fields[iM1][jM2][k].th_1 - fields[iM1][jP2][k].th_1 + fields[iP1][jM2][k].th_1 - fields[iP1][jP2][k].th_1
                     - fields[iM1][jM1][kM1].th_1 - fields[iM1][jM1][kP1].th_1 + fields[iM1][jP1][kM1].th_1 + fields[iM1][jP1][kP1].th_1
                     - fields[iP1][jM1][kM1].th_1 - fields[iP1][jM1][kP1].th_1 + fields[iP1][jP1][kM1].th_1 + fields[iP1][jP1][kP1].th_1)/30
                  + (- fields[i][jM2][k].th_1 + fields[i][jP2][k].th_1)/20
                  - (  fields[i][jM1][kM1].th_1 + fields[i][jM1][kP1].th_1 - fields[i][jP1][kM1].th_1 - fields[i][jP1][kP1].th_1
                     + fields[iM1][jM1][k].th_1 - fields[iM1][jP1][k].th_1 + fields[iP1][jM1][k].th_1 - fields[iP1][jP1][k].th_1)/12
                  + (- fields[i][jM1][k].th_1 + fields[i][jP1][k].th_1)*4.0/15;
#ifdef WRAP_DTH
                  dy_th_1 = wrap(dy_th_1);
#endif
   dy_th_2 = 
                    (  fields[iM1][jM1][kM2].th_2 + fields[iM1][jM1][kP2].th_2 - fields[iM1][jP1][kM2].th_2 - fields[iM1][jP1][kP2].th_2
                     + fields[iP1][jM1][kM2].th_2 + fields[iP1][jM1][kP2].th_2 - fields[iP1][jP1][kM2].th_2 - fields[iP1][jP1][kP2].th_2
                     + fields[iM2][jM1][kM1].th_2 + fields[iM2][jM1][kP1].th_2 - fields[iM2][jP1][kM1].th_2 - fields[iM2][jP1][kP1].th_2
                     + fields[iP2][jM1][kM1].th_2 + fields[iP2][jM1][kP1].th_2 - fields[iP2][jP1][kM1].th_2 - fields[iP2][jP1][kP1].th_2)/120
                  + (  fields[i][jM2][kM1].th_2 + fields[i][jM2][kP1].th_2 - fields[i][jP2][kM1].th_2 - fields[i][jP2][kP1].th_2
                     + fields[iM1][jM2][k].th_2 - fields[iM1][jP2][k].th_2 + fields[iP1][jM2][k].th_2 - fields[iP1][jP2][k].th_2
                     - fields[iM1][jM1][kM1].th_2 - fields[iM1][jM1][kP1].th_2 + fields[iM1][jP1][kM1].th_2 + fields[iM1][jP1][kP1].th_2
                     - fields[iP1][jM1][kM1].th_2 - fields[iP1][jM1][kP1].th_2 + fields[iP1][jP1][kM1].th_2 + fields[iP1][jP1][kP1].th_2)/30
                  + (- fields[i][jM2][k].th_2 + fields[i][jP2][k].th_2)/20
                  + (- fields[i][jM1][kM1].th_2 - fields[i][jM1][kP1].th_2 + fields[i][jP1][kM1].th_2 + fields[i][jP1][kP1].th_2
                     - fields[iM1][jM1][k].th_2 + fields[iM1][jP1][k].th_2 - fields[iP1][jM1][k].th_2 + fields[iP1][jP1][k].th_2)/12
                  + (- fields[i][jM1][k].th_2 + fields[i][jP1][k].th_2)*4.0/15;
#ifdef WRAP_DTH
                  dy_th_2 = wrap(dy_th_2);
#endif
   dy_th_3 = 
                    (  fields[iM1][jM1][kM2].th_3 + fields[iM1][jM1][kP2].th_3 - fields[iM1][jP1][kM2].th_3 - fields[iM1][jP1][kP2].th_3
                     + fields[iP1][jM1][kM2].th_3 + fields[iP1][jM1][kP2].th_3 - fields[iP1][jP1][kM2].th_3 - fields[iP1][jP1][kP2].th_3
                     + fields[iM2][jM1][kM1].th_3 + fields[iM2][jM1][kP1].th_3 - fields[iM2][jP1][kM1].th_3 - fields[iM2][jP1][kP1].th_3
                     + fields[iP2][jM1][kM1].th_3 + fields[iP2][jM1][kP1].th_3 - fields[iP2][jP1][kM1].th_3 - fields[iP2][jP1][kP1].th_3)/120
                  + (  fields[i][jM2][kM1].th_3 + fields[i][jM2][kP1].th_3 - fields[i][jP2][kM1].th_3 - fields[i][jP2][kP1].th_3
                     + fields[iM1][jM2][k].th_3 - fields[iM1][jP2][k].th_3 + fields[iP1][jM2][k].th_3 - fields[iP1][jP2][k].th_3
                     - fields[iM1][jM1][kM1].th_3 - fields[iM1][jM1][kP1].th_3 + fields[iM1][jP1][kM1].th_3 + fields[iM1][jP1][kP1].th_3
                     - fields[iP1][jM1][kM1].th_3 - fields[iP1][jM1][kP1].th_3 + fields[iP1][jP1][kM1].th_3 + fields[iP1][jP1][kP1].th_3)/30
                  - (  fields[i][jM2][k].th_3 - fields[i][jP2][k].th_3)/20
                  + (- fields[i][jM1][kM1].th_3 - fields[i][jM1][kP1].th_3 + fields[i][jP1][kM1].th_3 + fields[i][jP1][kP1].th_3
                     - fields[iM1][jM1][k].th_3 + fields[iM1][jP1][k].th_3 - fields[iP1][jM1][k].th_3 + fields[iP1][jP1][k].th_3)/12
                  - (  fields[i][jM1][k].th_3 - fields[i][jP1][k].th_3)*4.0/15;
#ifdef WRAP_DTH
                  dy_th_3 = wrap(dy_th_3);
#endif
   dy_A_1 = 
                   (  fields[iM1][jM1][kM2].A_1 + fields[iM1][jM1][kP2].A_1 - fields[iM1][jP1][kM2].A_1 - fields[iM1][jP1][kP2].A_1
                    + fields[iP1][jM1][kM2].A_1 + fields[iP1][jM1][kP2].A_1 - fields[iP1][jP1][kM2].A_1 - fields[iP1][jP1][kP2].A_1
                    + fields[iM2][jM1][kM1].A_1 + fields[iM2][jM1][kP1].A_1 - fields[iM2][jP1][kM1].A_1 - fields[iM2][jP1][kP1].A_1
                    + fields[iP2][jM1][kM1].A_1 + fields[iP2][jM1][kP1].A_1 - fields[iP2][jP1][kM1].A_1 - fields[iP2][jP1][kP1].A_1)/120
                 + (  fields[i][jM2][kM1].A_1 + fields[i][jM2][kP1].A_1 - fields[i][jP2][kM1].A_1 - fields[i][jP2][kP1].A_1
                    + fields[iM1][jM2][k].A_1 - fields[iM1][jP2][k].A_1 + fields[iP1][jM2][k].A_1 - fields[iP1][jP2][k].A_1
                    - fields[iM1][jM1][kM1].A_1 - fields[iM1][jM1][kP1].A_1 + fields[iM1][jP1][kM1].A_1 + fields[iM1][jP1][kP1].A_1
                    - fields[iP1][jM1][kM1].A_1 - fields[iP1][jM1][kP1].A_1 + fields[iP1][jP1][kM1].A_1 + fields[iP1][jP1][kP1].A_1)/30
                 + (- fields[i][jM2][k].A_1 + fields[i][jP2][k].A_1)/20
                 - (  fields[i][jM1][kM1].A_1 + fields[i][jM1][kP1].A_1 - fields[i][jP1][kM1].A_1 - fields[i][jP1][kP1].A_1
                    + fields[iM1][jM1][k].A_1 - fields[iM1][jP1][k].A_1 + fields[iP1][jM1][k].A_1 - fields[iP1][jP1][k].A_1)/12
                 + (- fields[i][jM1][k].A_1 + fields[i][jP1][k].A_1)*4.0/15;
   dy_A_3 = 
                   (  fields[iM1][jM1][kM2].A_3 + fields[iM1][jM1][kP2].A_3 - fields[iM1][jP1][kM2].A_3 - fields[iM1][jP1][kP2].A_3
                    + fields[iP1][jM1][kM2].A_3 + fields[iP1][jM1][kP2].A_3 - fields[iP1][jP1][kM2].A_3 - fields[iP1][jP1][kP2].A_3
                    + fields[iM2][jM1][kM1].A_3 + fields[iM2][jM1][kP1].A_3 - fields[iM2][jP1][kM1].A_3 - fields[iM2][jP1][kP1].A_3
                    + fields[iP2][jM1][kM1].A_3 + fields[iP2][jM1][kP1].A_3 - fields[iP2][jP1][kM1].A_3 - fields[iP2][jP1][kP1].A_3)/120
                 + (  fields[i][jM2][kM1].A_3 + fields[i][jM2][kP1].A_3 - fields[i][jP2][kM1].A_3 - fields[i][jP2][kP1].A_3
                    + fields[iM1][jM2][k].A_3 - fields[iM1][jP2][k].A_3 + fields[iP1][jM2][k].A_3 - fields[iP1][jP2][k].A_3
                    - fields[iM1][jM1][kM1].A_3 - fields[iM1][jM1][kP1].A_3 + fields[iM1][jP1][kM1].A_3 + fields[iM1][jP1][kP1].A_3
                    - fields[iP1][jM1][kM1].A_3 - fields[iP1][jM1][kP1].A_3 + fields[iP1][jP1][kM1].A_3 + fields[iP1][jP1][kP1].A_3)/30
                 - (  fields[i][jM2][k].A_3 - fields[i][jP2][k].A_3)/20
                 + (- fields[i][jM1][kM1].A_3 - fields[i][jM1][kP1].A_3 + fields[i][jP1][kM1].A_3 + fields[i][jP1][kP1].A_3
                    - fields[iM1][jM1][k].A_3 + fields[iM1][jP1][k].A_3 - fields[iP1][jM1][k].A_3 + fields[iP1][jP1][k].A_3)/12
                 - (  fields[i][jM1][k].A_3 - fields[i][jP1][k].A_3)*4.0/15;
   dz_th_1 = 
                    (  fields[iM1][jM2][kM1].th_1 - fields[iM1][jM2][kP1].th_1 + fields[iM1][jP2][kM1].th_1 - fields[iM1][jP2][kP1].th_1
                     + fields[iP1][jM2][kM1].th_1 - fields[iP1][jM2][kP1].th_1 + fields[iP1][jP2][kM1].th_1 - fields[iP1][jP2][kP1].th_1
                     + fields[iM2][jM1][kM1].th_1 - fields[iM2][jM1][kP1].th_1 + fields[iM2][jP1][kM1].th_1 - fields[iM2][jP1][kP1].th_1
                     + fields[iP2][jM1][kM1].th_1 - fields[iP2][jM1][kP1].th_1 + fields[iP2][jP1][kM1].th_1 - fields[iP2][jP1][kP1].th_1)/120
                  + (  fields[i][jM1][kM2].th_1 - fields[i][jM1][kP2].th_1 + fields[i][jP1][kM2].th_1 - fields[i][jP1][kP2].th_1
                     + fields[iM1][j][kM2].th_1 - fields[iM1][j][kP2].th_1 + fields[iP1][j][kM2].th_1 - fields[iP1][j][kP2].th_1
                     - fields[iM1][jM1][kM1].th_1 + fields[iM1][jM1][kP1].th_1 - fields[iM1][jP1][kM1].th_1 + fields[iM1][jP1][kP1].th_1
                     - fields[iP1][jM1][kM1].th_1 + fields[iP1][jM1][kP1].th_1 - fields[iP1][jP1][kM1].th_1 + fields[iP1][jP1][kP1].th_1)/30
                  + (- fields[i][j][kM2].th_1 + fields[i][j][kP2].th_1)/20
                  - (  fields[i][jM1][kM1].th_1 - fields[i][jM1][kP1].th_1 + fields[i][jP1][kM1].th_1 - fields[i][jP1][kP1].th_1
                     + fields[iM1][j][kM1].th_1 - fields[iM1][j][kP1].th_1 + fields[iP1][j][kM1].th_1 - fields[iP1][j][kP1].th_1)/12
                  + (- fields[i][j][kM1].th_1 + fields[i][j][kP1].th_1)*4.0/15;
#ifdef WRAP_DTH
                  dz_th_1 = wrap(dz_th_1);
#endif
   dz_th_2 = 
                    (  fields[iM1][jM2][kM1].th_2 - fields[iM1][jM2][kP1].th_2 + fields[iM1][jP2][kM1].th_2 - fields[iM1][jP2][kP1].th_2
                     + fields[iP1][jM2][kM1].th_2 - fields[iP1][jM2][kP1].th_2 + fields[iP1][jP2][kM1].th_2 - fields[iP1][jP2][kP1].th_2
                     + fields[iM2][jM1][kM1].th_2 - fields[iM2][jM1][kP1].th_2 + fields[iM2][jP1][kM1].th_2 - fields[iM2][jP1][kP1].th_2
                     + fields[iP2][jM1][kM1].th_2 - fields[iP2][jM1][kP1].th_2 + fields[iP2][jP1][kM1].th_2 - fields[iP2][jP1][kP1].th_2)/120
                  + (  fields[i][jM1][kM2].th_2 - fields[i][jM1][kP2].th_2 + fields[i][jP1][kM2].th_2 - fields[i][jP1][kP2].th_2
                     + fields[iM1][j][kM2].th_2 - fields[iM1][j][kP2].th_2 + fields[iP1][j][kM2].th_2 - fields[iP1][j][kP2].th_2
                     - fields[iM1][jM1][kM1].th_2 + fields[iM1][jM1][kP1].th_2 - fields[iM1][jP1][kM1].th_2 + fields[iM1][jP1][kP1].th_2
                     - fields[iP1][jM1][kM1].th_2 + fields[iP1][jM1][kP1].th_2 - fields[iP1][jP1][kM1].th_2 + fields[iP1][jP1][kP1].th_2)/30
                  + (- fields[i][j][kM2].th_2 + fields[i][j][kP2].th_2)/20
                  - (  fields[i][jM1][kM1].th_2 - fields[i][jM1][kP1].th_2 + fields[i][jP1][kM1].th_2 - fields[i][jP1][kP1].th_2
                     + fields[iM1][j][kM1].th_2 - fields[iM1][j][kP1].th_2 + fields[iP1][j][kM1].th_2 - fields[iP1][j][kP1].th_2)/12
                  + (- fields[i][j][kM1].th_2 + fields[i][j][kP1].th_2)*4.0/15;
#ifdef WRAP_DTH
                  dz_th_2 = wrap(dz_th_2);
#endif
   dz_th_3 = 
                    (  fields[iM1][jM2][kM1].th_3 - fields[iM1][jM2][kP1].th_3 + fields[iM1][jP2][kM1].th_3 - fields[iM1][jP2][kP1].th_3
                     + fields[iP1][jM2][kM1].th_3 - fields[iP1][jM2][kP1].th_3 + fields[iP1][jP2][kM1].th_3 - fields[iP1][jP2][kP1].th_3
                     + fields[iM2][jM1][kM1].th_3 - fields[iM2][jM1][kP1].th_3 + fields[iM2][jP1][kM1].th_3 - fields[iM2][jP1][kP1].th_3
                     + fields[iP2][jM1][kM1].th_3 - fields[iP2][jM1][kP1].th_3 + fields[iP2][jP1][kM1].th_3 - fields[iP2][jP1][kP1].th_3)/120
                  + (  fields[i][jM1][kM2].th_3 - fields[i][jM1][kP2].th_3 + fields[i][jP1][kM2].th_3 - fields[i][jP1][kP2].th_3
                     + fields[iM1][j][kM2].th_3 - fields[iM1][j][kP2].th_3 + fields[iP1][j][kM2].th_3 - fields[iP1][j][kP2].th_3
                     - fields[iM1][jM1][kM1].th_3 + fields[iM1][jM1][kP1].th_3 - fields[iM1][jP1][kM1].th_3 + fields[iM1][jP1][kP1].th_3
                     - fields[iP1][jM1][kM1].th_3 + fields[iP1][jM1][kP1].th_3 - fields[iP1][jP1][kM1].th_3 + fields[iP1][jP1][kP1].th_3)/30
                  + (- fields[i][j][kM2].th_3 + fields[i][j][kP2].th_3)/20
                  - (  fields[i][jM1][kM1].th_3 - fields[i][jM1][kP1].th_3 + fields[i][jP1][kM1].th_3 - fields[i][jP1][kP1].th_3
                     + fields[iM1][j][kM1].th_3 - fields[iM1][j][kP1].th_3 + fields[iP1][j][kM1].th_3 - fields[iP1][j][kP1].th_3)/12
                  + (- fields[i][j][kM1].th_3 + fields[i][j][kP1].th_3)*4.0/15;
#ifdef WRAP_DTH
                  dz_th_3 = wrap(dz_th_3);
#endif
   dz_A_1 =
                   (  fields[iM1][jM2][kM1].A_1 - fields[iM1][jM2][kP1].A_1 + fields[iM1][jP2][kM1].A_1 - fields[iM1][jP2][kP1].A_1
                    + fields[iP1][jM2][kM1].A_1 - fields[iP1][jM2][kP1].A_1 + fields[iP1][jP2][kM1].A_1 - fields[iP1][jP2][kP1].A_1
                    + fields[iM2][jM1][kM1].A_1 - fields[iM2][jM1][kP1].A_1 + fields[iM2][jP1][kM1].A_1 - fields[iM2][jP1][kP1].A_1
                    + fields[iP2][jM1][kM1].A_1 - fields[iP2][jM1][kP1].A_1 + fields[iP2][jP1][kM1].A_1 - fields[iP2][jP1][kP1].A_1)/120
                 + (  fields[i][jM1][kM2].A_1 - fields[i][jM1][kP2].A_1 + fields[i][jP1][kM2].A_1 - fields[i][jP1][kP2].A_1
                    + fields[iM1][j][kM2].A_1 - fields[iM1][j][kP2].A_1 + fields[iP1][j][kM2].A_1 - fields[iP1][j][kP2].A_1
                    - fields[iM1][jM1][kM1].A_1 + fields[iM1][jM1][kP1].A_1 - fields[iM1][jP1][kM1].A_1 + fields[iM1][jP1][kP1].A_1
                    - fields[iP1][jM1][kM1].A_1 + fields[iP1][jM1][kP1].A_1 - fields[iP1][jP1][kM1].A_1 + fields[iP1][jP1][kP1].A_1)/30
                 + (- fields[i][j][kM2].A_1 + fields[i][j][kP2].A_1)/20
                 - (  fields[i][jM1][kM1].A_1 - fields[i][jM1][kP1].A_1 + fields[i][jP1][kM1].A_1 - fields[i][jP1][kP1].A_1
                    + fields[iM1][j][kM1].A_1 - fields[iM1][j][kP1].A_1 + fields[iP1][j][kM1].A_1 - fields[iP1][j][kP1].A_1)/12
                 + (- fields[i][j][kM1].A_1 + fields[i][j][kP1].A_1)*4.0/15;
   dz_A_2 =
                   (  fields[iM1][jM2][kM1].A_2 - fields[iM1][jM2][kP1].A_2 + fields[iM1][jP2][kM1].A_2 - fields[iM1][jP2][kP1].A_2
                    + fields[iP1][jM2][kM1].A_2 - fields[iP1][jM2][kP1].A_2 + fields[iP1][jP2][kM1].A_2 - fields[iP1][jP2][kP1].A_2
                    + fields[iM2][jM1][kM1].A_2 - fields[iM2][jM1][kP1].A_2 + fields[iM2][jP1][kM1].A_2 - fields[iM2][jP1][kP1].A_2
                    + fields[iP2][jM1][kM1].A_2 - fields[iP2][jM1][kP1].A_2 + fields[iP2][jP1][kM1].A_2 - fields[iP2][jP1][kP1].A_2)/120
                 + (  fields[i][jM1][kM2].A_2 - fields[i][jM1][kP2].A_2 + fields[i][jP1][kM2].A_2 - fields[i][jP1][kP2].A_2
                    + fields[iM1][j][kM2].A_2 - fields[iM1][j][kP2].A_2 + fields[iP1][j][kM2].A_2 - fields[iP1][j][kP2].A_2
                    - fields[iM1][jM1][kM1].A_2 + fields[iM1][jM1][kP1].A_2 - fields[iM1][jP1][kM1].A_2 + fields[iM1][jP1][kP1].A_2
                    - fields[iP1][jM1][kM1].A_2 + fields[iP1][jM1][kP1].A_2 - fields[iP1][jP1][kM1].A_2 + fields[iP1][jP1][kP1].A_2)/30
                 + (- fields[i][j][kM2].A_2 + fields[i][j][kP2].A_2)/20
                 - (  fields[i][jM1][kM1].A_2 - fields[i][jM1][kP1].A_2 + fields[i][jP1][kM1].A_2 - fields[i][jP1][kP1].A_2
                    + fields[iM1][j][kM1].A_2 - fields[iM1][j][kP1].A_2 + fields[iP1][j][kM1].A_2 - fields[iP1][j][kP1].A_2)/12
                 + (- fields[i][j][kM1].A_2 + fields[i][j][kP1].A_2)*4.0/15;
}

double TAux::get_rho_static_core(TFields *fields, TRho &rho, 
                                 double h_00, double h_01, double h_02, double h_11, double h_12, double h_22, 
                                 double G_00, double G_01, double G_02, double G_11, double G_12, double G_22)
/* Individual terms are returned in the TRho struct.  */
{
   double AcrossGradTh_am[3][3];

   AcrossGradTh_am[0][0] = fields->A_2*dz_th_1 - fields->A_3*dy_th_1;
   AcrossGradTh_am[0][1] = fields->A_3*dx_th_1 - fields->A_1*dz_th_1;
   AcrossGradTh_am[0][2] = fields->A_1*dy_th_1 - fields->A_2*dx_th_1;
   AcrossGradTh_am[1][0] = fields->A_2*dz_th_2 - fields->A_3*dy_th_2;
   AcrossGradTh_am[1][1] = fields->A_3*dx_th_2 - fields->A_1*dz_th_2;
   AcrossGradTh_am[1][2] = fields->A_1*dy_th_2 - fields->A_2*dx_th_2;
   AcrossGradTh_am[2][0] = fields->A_2*dz_th_3 - fields->A_3*dy_th_3;
   AcrossGradTh_am[2][1] = fields->A_3*dx_th_3 - fields->A_1*dz_th_3;
   AcrossGradTh_am[2][2] = fields->A_1*dy_th_3 - fields->A_2*dx_th_3;

   rho.BR2 = ( pow(dy_A_3 - dz_A_2, 2) + pow(dz_A_1 - dx_A_3, 2) + pow(dx_A_2 - dy_A_1, 2) );

   rho.rH =  (  h_00*(  AcrossGradTh_am[0][0]*AcrossGradTh_am[0][0]
                      + AcrossGradTh_am[0][1]*AcrossGradTh_am[0][1]
                      + AcrossGradTh_am[0][2]*AcrossGradTh_am[0][2])
              + h_11*(  AcrossGradTh_am[1][0]*AcrossGradTh_am[1][0]
                      + AcrossGradTh_am[1][1]*AcrossGradTh_am[1][1]
                      + AcrossGradTh_am[1][2]*AcrossGradTh_am[1][2])
              + h_22*(  AcrossGradTh_am[2][0]*AcrossGradTh_am[2][0]
                      + AcrossGradTh_am[2][1]*AcrossGradTh_am[2][1]
                      + AcrossGradTh_am[2][2]*AcrossGradTh_am[2][2])
              + 2*(  h_01*(  AcrossGradTh_am[0][0]*AcrossGradTh_am[1][0]
                           + AcrossGradTh_am[0][1]*AcrossGradTh_am[1][1]
                           + AcrossGradTh_am[0][2]*AcrossGradTh_am[1][2])
                   + h_02*(  AcrossGradTh_am[0][0]*AcrossGradTh_am[2][0]
                           + AcrossGradTh_am[0][1]*AcrossGradTh_am[2][1]
                           + AcrossGradTh_am[0][2]*AcrossGradTh_am[2][2])
                   + h_12*(  AcrossGradTh_am[1][0]*AcrossGradTh_am[2][0]
                           + AcrossGradTh_am[1][1]*AcrossGradTh_am[2][1]
                           + AcrossGradTh_am[1][2]*AcrossGradTh_am[2][2])
                  )
             );

   if (0.0 > rho.rH) rho.rH = 0.0;

   rho.rG = (  G_00*( dx_th_1*dx_th_1 + dy_th_1*dy_th_1 + dz_th_1*dz_th_1 )
             + G_11*( dx_th_2*dx_th_2 + dy_th_2*dy_th_2 + dz_th_2*dz_th_2 )
             + G_22*( dx_th_3*dx_th_3 + dy_th_3*dy_th_3 + dz_th_3*dz_th_3 )
             + 2*(  G_01*( dx_th_1*dx_th_2 + dy_th_1*dy_th_2 + dz_th_1*dz_th_2 )
                  + G_02*( dx_th_1*dx_th_3 + dy_th_1*dy_th_3 + dz_th_1*dz_th_3 )
                  + G_12*( dx_th_2*dx_th_3 + dy_th_2*dy_th_3 + dz_th_2*dz_th_3 ) )
            );

   if (0.0 > rho.rG) rho.rG = 0.0;

   return(rho.BR2 + rho.rH + rho.rG);
}

double TAux::get_rho_static(TFields *fields, TRho &rho)
/* Individual terms are returned in the TRho struct.
     IMPORTANT: assumes that h, G and field derivatives are ready for use! */
{
   return(get_rho_static_core(fields, rho, h_00, h_01, h_02, h_11, h_12, h_22, G_00, G_01, G_02, G_11, G_12, G_22));
}

void TAux::subtract_drho_static_dA(TFields *fields, double &drdA1, double &drdA2, double &drdA3)
/* Effect on local cell only, NOT all affected cells.
     IMPORTANT: assumes h, G, M, theta field derivatives are ready for use! */
{
   double AcrossGradTh_am[3][3];

   AcrossGradTh_am[0][0] = fields->A_2*dz_th_1 - fields->A_3*dy_th_1;
   AcrossGradTh_am[0][1] = fields->A_3*dx_th_1 - fields->A_1*dz_th_1;
   AcrossGradTh_am[0][2] = fields->A_1*dy_th_1 - fields->A_2*dx_th_1;
   AcrossGradTh_am[1][0] = fields->A_2*dz_th_2 - fields->A_3*dy_th_2;
   AcrossGradTh_am[1][1] = fields->A_3*dx_th_2 - fields->A_1*dz_th_2;
   AcrossGradTh_am[1][2] = fields->A_1*dy_th_2 - fields->A_2*dx_th_2;
   AcrossGradTh_am[2][0] = fields->A_2*dz_th_3 - fields->A_3*dy_th_3;
   AcrossGradTh_am[2][1] = fields->A_3*dx_th_3 - fields->A_1*dz_th_3;
   AcrossGradTh_am[2][2] = fields->A_1*dy_th_3 - fields->A_2*dx_th_3;
   
   drdA1 -= 2 * (  h_00 * (dy_th_1 * AcrossGradTh_am[0][2] - dz_th_1 * AcrossGradTh_am[0][1])
                 + h_11 * (dy_th_2 * AcrossGradTh_am[1][2] - dz_th_2 * AcrossGradTh_am[1][1])
                 + h_22 * (dy_th_3 * AcrossGradTh_am[2][2] - dz_th_3 * AcrossGradTh_am[2][1])
                 + h_01 * (  (dy_th_1 * AcrossGradTh_am[1][2] - dz_th_1 * AcrossGradTh_am[1][1])
                           + (dy_th_2 * AcrossGradTh_am[0][2] - dz_th_2 * AcrossGradTh_am[0][1]) )
                 + h_02 * (  (dy_th_1 * AcrossGradTh_am[2][2] - dz_th_1 * AcrossGradTh_am[2][1])
                           + (dy_th_3 * AcrossGradTh_am[0][2] - dz_th_3 * AcrossGradTh_am[0][1]) )
                 + h_12 * (  (dy_th_2 * AcrossGradTh_am[2][2] - dz_th_2 * AcrossGradTh_am[2][1])
                           + (dy_th_3 * AcrossGradTh_am[1][2] - dz_th_3 * AcrossGradTh_am[1][1]) ) );

   drdA2 -= 2 * (  h_00 * (dz_th_1 * AcrossGradTh_am[0][0] - dx_th_1 * AcrossGradTh_am[0][2])
                 + h_11 * (dz_th_2 * AcrossGradTh_am[1][0] - dx_th_2 * AcrossGradTh_am[1][2])
                 + h_22 * (dz_th_3 * AcrossGradTh_am[2][0] - dx_th_3 * AcrossGradTh_am[2][2])
                 + h_01 * (  (dz_th_1 * AcrossGradTh_am[1][0] - dx_th_1 * AcrossGradTh_am[1][2])
                           + (dz_th_2 * AcrossGradTh_am[0][0] - dx_th_2 * AcrossGradTh_am[0][2]) )
                 + h_02 * (  (dz_th_1 * AcrossGradTh_am[2][0] - dx_th_1 * AcrossGradTh_am[2][2])
                           + (dz_th_3 * AcrossGradTh_am[0][0] - dx_th_3 * AcrossGradTh_am[0][2]) )
                 + h_12 * (  (dz_th_2 * AcrossGradTh_am[2][0] - dx_th_2 * AcrossGradTh_am[2][2])
                           + (dz_th_3 * AcrossGradTh_am[1][0] - dx_th_3 * AcrossGradTh_am[1][2]) ) );

   drdA3 -= 2 * (  h_00 * (dx_th_1 * AcrossGradTh_am[0][1] - dy_th_1 * AcrossGradTh_am[0][0])
                 + h_11 * (dx_th_2 * AcrossGradTh_am[1][1] - dy_th_2 * AcrossGradTh_am[1][0])
                 + h_22 * (dx_th_3 * AcrossGradTh_am[2][1] - dy_th_3 * AcrossGradTh_am[2][0])
                 + h_01 * (  (dx_th_1 * AcrossGradTh_am[1][1] - dy_th_1 * AcrossGradTh_am[1][0])
                           + (dx_th_2 * AcrossGradTh_am[0][1] - dy_th_2 * AcrossGradTh_am[0][0]) )
                 + h_02 * (  (dx_th_1 * AcrossGradTh_am[2][1] - dy_th_1 * AcrossGradTh_am[2][0])
                           + (dx_th_3 * AcrossGradTh_am[0][1] - dy_th_3 * AcrossGradTh_am[0][0]) )
                 + h_12 * (  (dx_th_2 * AcrossGradTh_am[2][1] - dy_th_2 * AcrossGradTh_am[2][0])
                           + (dx_th_3 * AcrossGradTh_am[1][1] - dy_th_3 * AcrossGradTh_am[1][0]) ) );
}

void TAux::subtract_offset_drho_dA(int di, int dj, int dk, double &drdA1, double &drdA2, double &drdA3)
/* (di, dj, dk)  is position in lattice units of the cell where A is varied, relative to the cell managed by this TAux object.
     Variations in the fields of the latter, central cell are NOT handled by this method!
     IMPORTANT: assumes that A field derivatives are ready for use! */
{
   double curlA_x = dy_A_3 - dz_A_2;
   double curlA_y = dz_A_1 - dx_A_3;
   double curlA_z = dx_A_2 - dy_A_1;

   switch (di)
   {
     case -2: switch (dj)
              {                         
                case -1: switch (dk)
                         {
                           case -1: drdA1 -= -(curlA_z-curlA_y)/60;
                                    drdA2 -= -curlA_x/60;
                                    drdA3 -=  curlA_x/60;
                                    break;
                           case 0:  drdA2 -= curlA_z/15;
                                    drdA3 -= -curlA_y/15;
                                    break;
                           case 1:  drdA1 -= -(curlA_z+curlA_y)/60;
                                    drdA2 -=  curlA_x/60;
                                    drdA3 -=  curlA_x/60;
                                    break;
                         }
                         break;
                case 0:  switch (dk)
                         {
                           case -1: drdA2 -= curlA_z/15;
                                    drdA3 -= -curlA_y/15;
                                    break;
                           case 0:  drdA2 -= -curlA_z/10;
                                    drdA3 -= curlA_y/10;
                                    break;
                           case 1:  drdA2 -= curlA_z/15;
                                    drdA3 -= -curlA_y/15;
                                    break;
                         }
                         break;
                case 1:  switch (dk)
                         {
                           case -1: drdA1 -= (curlA_z+curlA_y)/60;
                                    drdA2 -= -curlA_x/60;
                                    drdA3 -= -curlA_x/60;
                                    break;
                           case 0:  drdA2 -= curlA_z/15;
                                    drdA3 -= -curlA_y/15;
                                    break;
                           case 1:  drdA1 -= (curlA_z-curlA_y)/60;
                                    drdA2 -= curlA_x/60;
                                    drdA3 -= -curlA_x/60;
                                    break;
                         }
                         break;
              }
              break;
     case -1: switch (dj)
              {                         
                case -2: switch (dk)
                         {
                           case -1: drdA1 -= curlA_y/60;
                                    drdA2 -= (curlA_z-curlA_x)/60;
                                    drdA3 -= -curlA_y/60;
                                    break;
                           case 0:  drdA1 -= -curlA_z/15;
                                    drdA3 -= curlA_x/15;
                                    break;
                           case 1:  drdA1 -= -curlA_y/60;
                                    drdA2 -= (curlA_z+curlA_x)/60;
                                    drdA3 -= -curlA_y/60;
                                    break;
                         }
                         break;
                case -1: switch (dk)
                         {
                           case -2: drdA1 -= -curlA_z/60;
                                    drdA2 -= curlA_z/60;
                                    drdA3 -= -(curlA_y-curlA_x)/60;
                                    break;
                           case -1: drdA1 -= (curlA_z-curlA_y)/15;
                                    drdA2 -= -(curlA_z-curlA_x)/15;
                                    drdA3 -= (curlA_y-curlA_x)/15;
                                    break;
                           case 0:  drdA1 -= curlA_z/6;
                                    drdA2 -= -curlA_z/6;
                                    drdA3 -= (curlA_y-curlA_x)/6;
                                    break;
                           case 1:  drdA1 -= (curlA_z+curlA_y)/15;
                                    drdA2 -= -(curlA_z+curlA_x)/15;
                                    drdA3 -= (curlA_y-curlA_x)/15;
                                    break;
                           case 2:  drdA1 -= -curlA_z/60;
                                    drdA2 -= curlA_z/60;
                                    drdA3 -= -(curlA_y-curlA_x)/60;
                                    break;
                         }
                         break;
                 case 0: switch (dk)
                         {
                           case -2: drdA1 -= curlA_y/15;
                                    drdA2 -= -curlA_x/15;
                                    break;
                           case -1: drdA1 -= -curlA_y/6;
                                    drdA2 -= -(curlA_z-curlA_x)/6;
                                    drdA3 -= curlA_y/6;
                                    break;
                           case 0:  drdA2 -= -8*curlA_z/15;
                                    drdA3 -=  8*curlA_y/15;
                                    break;
                           case 1:  drdA1 -= curlA_y/6;
                                    drdA2 -= -(curlA_z+curlA_x)/6;
                                    drdA3 -= curlA_y/6;
                                    break;
                           case 2:  drdA1 -= -curlA_y/15;
                                    drdA2 -=  curlA_x/15;
                                    break;
                         }
                         break;
                 case 1: switch (dk)
                         {
                           case -2: drdA1 -= curlA_z/60;
                                    drdA2 -= curlA_z/60;
                                    drdA3 -= -(curlA_y+curlA_x)/60;
                                    break;
                           case -1: drdA1 -= -(curlA_z+curlA_y)/15;
                                    drdA2 -= -(curlA_z-curlA_x)/15;
                                    drdA3 -=  (curlA_y+curlA_x)/15;
                                    break;
                           case 0:  drdA1 -= -curlA_z/6;
                                    drdA2 -= -curlA_z/6;
                                    drdA3 -= (curlA_y+curlA_x)/6;
                                    break;
                           case 1:  drdA1 -= -(curlA_z-curlA_y)/15;
                                    drdA2 -= -(curlA_z+curlA_x)/15;
                                    drdA3 -=  (curlA_y+curlA_x)/15;
                                    break;
                           case 2:  drdA1 -= curlA_z/60;
                                    drdA2 -= curlA_z/60;
                                    drdA3 -= -(curlA_y+curlA_x)/60;
                                    break;
                         }
                         break;
                 case 2: switch (dk)
                         {
                           case -1: drdA1 -= curlA_y/60;
                                    drdA2 -= (curlA_z-curlA_x)/60;
                                    drdA3 -= -curlA_y/60;
                                    break;
                           case 0:  drdA1 -=  curlA_z/15;
                                    drdA3 -= -curlA_x/15;
                                    break;
                           case 1:  drdA1 -= -curlA_y/60;
                                    drdA2 -= (curlA_z+curlA_x)/60;
                                    drdA3 -= -curlA_y/60;
                                    break;
                         }
                         break;
              }
              break;
     case 0:  switch (dj)
              {                         
                case -2: switch (dk)
                         {
                           case -1: drdA1 -= -curlA_z/15;
                                    drdA3 -=  curlA_x/15;
                                    break;
                           case 0:  drdA1 -=  curlA_z/10;
                                    drdA3 -= -curlA_x/10;
                                    break;
                           case 1:  drdA1 -= -curlA_z/15;
                                    drdA3 -=  curlA_x/15;
                                    break;
                         }
                         break;
                case -1: switch (dk)
                         {
                           case -2: drdA1 -=  curlA_y/15;
                                    drdA2 -= -curlA_x/15;
                                    break;
                           case -1: drdA1 -= (curlA_z-curlA_y)/6;
                                    drdA2 -= curlA_x/6;
                                    drdA3 -= -curlA_x/6;
                                    break;
                           case 0:  drdA1 -= 8*curlA_z/15;
                                    drdA3 -= -8*curlA_x/15;
                                    break;
                           case 1:  drdA1 -= (curlA_z+curlA_y)/6;
                                    drdA2 -= -curlA_x/6;
                                    drdA3 -= -curlA_x/6;
                                    break;
                           case 2:  drdA1 -= -curlA_y/15;
                                    drdA2 -=  curlA_x/15;
                                    break;
                         }
                         break;
                 case 0: switch (dk)
                         { 
                           case -2: drdA1 -= -curlA_y/10;
                                    drdA2 -=  curlA_x/10;
                                    break;
                           case -1: drdA1 -= -8*curlA_y/15;
                                    drdA2 -=  8*curlA_x/15;
                                    break;
                           case 1:  drdA1 -=  8*curlA_y/15;
                                    drdA2 -= -8*curlA_x/15;
                                    break;
                           case 2:  drdA1 -=  curlA_y/10;
                                    drdA2 -= -curlA_x/10;
                                    break;
                         }
                         break;
                 case 1: switch (dk)
                         {
                           case -2: drdA1 -=  curlA_y/15;
                                    drdA2 -= -curlA_x/15;
                                    break;
                           case -1: drdA1 -= -(curlA_z+curlA_y)/6;
                                    drdA2 -= curlA_x/6;
                                    drdA3 -= curlA_x/6;
                                    break;
                           case 0:  drdA1 -= -8*curlA_z/15;
                                    drdA3 -=  8*curlA_x/15;
                                    break;
                            case 1: drdA1 -= -(curlA_z-curlA_y)/6;
                                    drdA2 -= -curlA_x/6;
                                    drdA3 -=  curlA_x/6;
                                    break;
                            case 2: drdA1 -= -curlA_y/15;
                                    drdA2 -=  curlA_x/15;
                                    break;
                         }
                         break;
                 case 2: switch (dk)
                         {
                           case -1: drdA1 -=  curlA_z/15;
                                    drdA3 -= -curlA_x/15;
                                    break;
                           case 0:  drdA1 -= -curlA_z/10;
                                    drdA3 -=  curlA_x/10;
                                    break;
                           case 1:  drdA1 -=  curlA_z/15;
                                    drdA3 -= -curlA_x/15;
                                    break;
                         }
                         break;
              }
              break;
     case 1:  switch (dj)
              {                         
                case -2: switch (dk)
                         {
                           case -1: drdA1 -= curlA_y/60;
                                    drdA2 -= -(curlA_z+curlA_x)/60;
                                    drdA3 -= curlA_y/60;
                                    break;
                           case 0:  drdA1 -= -curlA_z/15;
                                    drdA3 -=  curlA_x/15;
                                    break;
                           case 1:  drdA1 -= -curlA_y/60;
                                    drdA2 -= -(curlA_z-curlA_x)/60;
                                    drdA3 -= curlA_y/60;
                                    break;
                         }
                         break;
                case -1: switch (dk)
                         {
                           case -2: drdA1 -= -curlA_z/60;
                                    drdA2 -= -curlA_z/60;
                                    drdA3 -= (curlA_y+curlA_x)/60;
                                    break;
                           case -1: drdA1 -= (curlA_z-curlA_y)/15;
                                    drdA2 -= (curlA_z+curlA_x)/15;
                                    drdA3 -= -(curlA_y+curlA_x)/15;
                                    break;
                           case 0:  drdA1 -= curlA_z/6;
                                    drdA2 -= curlA_z/6;
                                    drdA3 -= -(curlA_y+curlA_x)/6;
                                    break;
                           case 1:  drdA1 -= (curlA_z+curlA_y)/15;
                                    drdA2 -= (curlA_z-curlA_x)/15;
                                    drdA3 -= -(curlA_y+curlA_x)/15;
                                    break;
                           case 2:  drdA1 -= -curlA_z/60;
                                    drdA2 -= -curlA_z/60;
                                    drdA3 -= (curlA_y+curlA_x)/60;
                                    break;
                         }
                         break;
                 case 0: switch (dk)
                         { 
                           case -2: drdA1 -=  curlA_y/15;
                                    drdA2 -= -curlA_x/15;
                                    break;
                           case -1: drdA1 -= -curlA_y/6;
                                    drdA2 -= (curlA_z+curlA_x)/6;
                                    drdA3 -= -curlA_y/6;
                                    break;
                           case 0:  drdA2 -=  8*curlA_z/15;
                                    drdA3 -= -8*curlA_y/15;
                                    break;
                           case 1:  drdA1 -= curlA_y/6;
                                    drdA2 -= (curlA_z-curlA_x)/6;
                                    drdA3 -= -curlA_y/6;
                                    break;
                           case 2:  drdA1 -= -curlA_y/15;
                                    drdA2 -=  curlA_x/15;
                                    break;
                         }
                         break;
                 case 1: switch (dk)
                         { 
                           case -2: drdA1 -=  curlA_z/60;
                                    drdA2 -= -curlA_z/60;
                                    drdA3 -= (curlA_y-curlA_x)/60;
                                    break;
                           case -1: drdA1 -= -(curlA_z+curlA_y)/15;
                                    drdA2 -=  (curlA_z+curlA_x)/15;
                                    drdA3 -= -(curlA_y-curlA_x)/15;
                                    break;
                           case 0:  drdA1 -= -curlA_z/6;
                                    drdA2 -=  curlA_z/6;
                                    drdA3 -= -(curlA_y-curlA_x)/6;
                                    break;
                           case 1:  drdA1 -= -(curlA_z-curlA_y)/15;
                                    drdA2 -=  (curlA_z-curlA_x)/15;
                                    drdA3 -= -(curlA_y-curlA_x)/15;
                                    break;
                           case 2:  drdA1 -=  curlA_z/60;
                                    drdA2 -= -curlA_z/60;
                                    drdA3 -= (curlA_y-curlA_x)/60;
                                    break;
                         }
                         break;
                 case 2: switch (dk)
                         {
                           case -1: drdA1 -= curlA_y/60;
                                    drdA2 -= -(curlA_z+curlA_x)/60;
                                    drdA3 -= curlA_y/60;
                                    break;
                           case 0:  drdA1 -=  curlA_z/15;
                                    drdA3 -= -curlA_x/15;
                                    break;
                           case 1:  drdA1 -= -curlA_y/60;
                                    drdA2 -= -(curlA_z-curlA_x)/60;
                                    drdA3 -= curlA_y/60;
                                    break;
                         }
                         break;
              }
              break;
     case 2:  switch (dj)
              {                         
                case -1: switch (dk)
                         {
                           case -1: drdA1 -= -(curlA_z-curlA_y)/60;
                                    drdA2 -= -curlA_x/60;
                                    drdA3 -=  curlA_x/60;
                                    break;
                           case 0:  drdA2 -= -curlA_z/15;
                                    drdA3 -=  curlA_y/15;
                                    break;
                           case 1:  drdA1 -= -(curlA_z+curlA_y)/60;
                                    drdA2 -= curlA_x/60;
                                    drdA3 -= curlA_x/60;
                                    break;
                         }
                         break;
                 case 0: switch (dk)
                         { 
                           case -1: drdA2 -= -curlA_z/15;
                                    drdA3 -=  curlA_y/15;
                                    break;
                           case 0:  drdA2 -=  curlA_z/10;
                                    drdA3 -= -curlA_y/10;
                                    break;
                           case 1:  drdA2 -= -curlA_z/15;
                                    drdA3 -=  curlA_y/15;
                                    break;
                         }
                         break;
                 case 1: switch (dk)
                         { 
                           case -1: drdA1 -= (curlA_z+curlA_y)/60;
                                    drdA2 -= -curlA_x/60;
                                    drdA3 -= -curlA_x/60;
                                    break;
                           case 0:  drdA2 -= -curlA_z/15;
                                    drdA3 -=  curlA_y/15;
                                    break;
                           case 1:  drdA1 -= (curlA_z-curlA_y)/60;
                                    drdA2 -=  curlA_x/60;
                                    drdA3 -= -curlA_x/60;
                                    break;
                         }
                         break;
              }
              break;
    }
}

void TAux::subtract_drho_static_dth(TFields *fields, double &drdth1, double &drdth2, double &drdth3)
/* Effect on local cell only, NOT all affected cells.
     IMPORTANT: assumes th derivatives are ready for use!  */
{
   double A_1 = fields->A_1;
   double A_2 = fields->A_2;
   double A_3 = fields->A_3;

   double th_1 = fields->th_1;
   double th_2 = fields->th_2;
   double th_3 = fields->th_3;

   double dh00dth, dh01dth, dh02dth, dh11dth, dh12dth, dh22dth, dg00dth, dg01dth, dg02dth, dg11dth, dg12dth, dg22dth;
   
   double th_1R2 = th_1 * th_1;
   double th_2R2 = th_2 * th_2;
   double th_3R2 = th_3 * th_3;

   double thR2 = th_1R2 + th_2R2 + th_3R2;
   double th = sqrt(thR2);

   double th_1R3 = th_1 * th_1R2;
   double th_2R3 = th_2 * th_2R2;
   double th_3R3 = th_3 * th_3R2;

   double th_1R4 = th_1R2 * th_1R2;
   double th_2R4 = th_2R2 * th_2R2;
   double th_3R4 = th_3R2 * th_3R2;

   double th_1R5 = th_1R3 * th_1R2;
   double th_2R5 = th_2R3 * th_2R2;
   double th_3R5 = th_3R3 * th_3R2;

   double th_3R6 = th_3R4 * th_3R2;

   double cth, sth, cthR2, sthR2, thR3, thR4, thR5, thR6, thR7, sththP2cth, sththP2cthR2, cthM1, cthM1R2;

   // d rho / d th_1
   
   if (thR2 < SMALL_THETA_SQUARED)
   {
      dh00dth = th_3*th_2/6+(16*th_2R2-9*th_3R2)*th_1/180;
      dh01dth = th_2/3-th_3*th_1/6-(24*th_2*th_1R2+8*th_2R3+13*th_3R2*th_2)/180;
      dh02dth = -th_3/12-th_2*th_1/4+(18*th_3*th_1R2+6*th_3*th_2R2+th_3R3)/180;
      dh11dth = -2*th_1/3-th_3*th_2/6+(32*th_1R3+(16*th_2R2+17*th_3R2)*th_1)/180;
      dh12dth = -1/2+(9*th_1R2+3*th_2R2+th_3R2)/24+th_3*th_2*th_1/15;
      dh22dth = th_1/2-(6*th_1R3+(6*th_2R2+th_3R2)*th_1)/36;

      dg00dth = (th_2R2+th_3R2)*th_1/720;
      dg01dth = th_2/48-(3*th_2*th_1R2+th_2R3+th_3R2*th_2)/1440;
      dg02dth = th_3/48-(3*th_3*th_1R2+th_3*th_2R2+th_3R3)/1440;
      dg11dth = -th_1/24+(2*th_1R3+(th_2R2+2*th_3R2)*th_1)/720;
      dg12dth = -th_3*th_2*th_1/720;
      dg22dth = -th_1/24+(2*th_1R3+(2*th_2R2+th_3R2)*th_1)/720;
   }
   else
   {
      sth = sin(th);
      cth = cos(th);

      cthR2 = cth * cth;
      sthR2 = sth * sth;
      
      sththP2cth = sth*th+2*cth;
      sththP2cthR2 = sththP2cth*sththP2cth;

      cthM1 = cth - 1.0;
      cthM1R2 = cthM1*cthM1;

      thR3 = thR2 * th;
      thR4 = thR2 * thR2;
      thR5 = thR2 * thR3;
      thR6 = thR2 * thR4;
      thR7 = thR6 * th;
      
      double norm = 1 / (thR3 * thR2 * thR3 * thR2);

      dh00dth = norm * ( 2*(th_1*(thR2*(th_2R4*sththP2cthR2+thR2*(th_2R2 *(cthR2*thR2+sth*(sth-6*cth*th) -8*cthR2) -th_1R2 *(sthR2*thR2 +3 *(cth*th*(7*sth-cth*th) +2*(4*cthR2-sthR2)))) +th_2R2*th_3R2 *(sthR2*thR2+4*(cth*(sth*th+cth)+1))) +th_1R2*(th_2R4+th_1R2*(2*th_2R2+th_1R2))*sththP2cth *(cth*thR2-(5*sth*th+8*cth)) +thR2*(th_1R2*(th_3R2*(th*(cth*((-cthM1)*th+17*sth)+3*sthR2*th) +5*(4*(cthR2+1)-sth*th)) -th_2R2*(thR2*(sth*(cth*th-7*sth)+5*cthR2) +sth*(3*sth-40*cth*th)-40*cthR2)) +th*(th_1*th_2*th_3*(th*(sth*(th-5*sth)+3*(-cthM1)*cth)+5*(-cthM1)*sth) +thR3*(3*sth*(cth*th-sth)+4*cthR2)) -th_1R4*(thR2*(sth*(cth*th-6*sth)+5*cthR2) +3*(sth*(sth-12*cth*th)-12*cthR2))+4*th_3R4)) +th_3*(th_1*th_3 *(th_1R4*(sth*thR2*(cth*th-5*sth) +2*(cthR2*thR2-(9*cth*sth*th+8*(cthR2+1)))) -(thR4*(3*cthM1*sth*th+2*(cth*(3*cth-2)+3))+16*th_1R2*th_3R2)) +th_2*(th_1R3*th_2*th_3 *(sth*thR2*(cth*th-5*sth) +2*(cthR2*thR2-(9*cth*sth*th+8*(cthR2+1)))) +thR5*((sthR2+cthM1*cth)*th+cthM1*sth)))) );
      dh01dth = norm * ( th_2*(th_2R4*thR2*sththP2cthR2+th_1R2*(2*(th_2R2+th_1R2)*th_3R2 *(sth*thR2*(cth*th-5*sth) +2 *(cthR2*thR2 -(9*cth*sth*th+8*(cthR2+1)))) -th_1R2*thR2 *(sth*thR2*(2*cth*th-11*sth) +2 *(5*cthR2*thR2 +sth*(3*sth-34*cth*th) -34*cthR2)))) +thR2*(th_2*(th_2R2*th_3R2*(sthR2*thR2+4*(cth*(sth*th+cth)+1)) -th_1R2*thR2 *((2*sthR2-5*cthR2)*thR2+sth*(36*cth*th-11*sth)+40*cthR2)) +th*(th_2R3*th*(cthR2*thR2+sth*(sth-6*cth*th)-8*cthR2) +th_1*(th_2R2-th_1R2)*th_3 *(th*(sth*(th-5*sth)+3*(-cthM1)*cth)+5*(-cthM1)*sth)) +th_1R2*th_2 *(th_3R2*(5*sthR2*thR2+2*(th*((-cthM1)*cth*th+5*(3*cth-1)*sth) +18*(cthR2+1))) -2*th_2R2 *(thR2*(sth*(cth*th-6*sth)+5*cthR2) +3*(sth*(sth-12*cth*th)-12*cthR2)))) +2*(th_2*(th_3R2*(2*th_3R2*thR2-(thR4*(cthM1*sth*th+2*(cthR2+1))+16*th_1R2*th_3R2)) +thR6*(sth*(cth*th-sth)+2*cthR2)) +th_1*(th_1*th_2*(th_2R4+th_1R2*(2*th_2R2+th_1R2))*sththP2cth *(cth*thR2-(5*sth*th+8*cth)) -th_3*thR5*((sthR2+cthM1*cth)*th+cthM1*sth))) );
      dh02dth = norm * ( th_3*(th_2R4*thR2*sththP2cthR2+th_1*(th_1*(2*(th_2R4+th_1R2*(2*th_2R2+th_1R2)) *sththP2cth *(cth*thR2-(5*sth*th+8*cth)) -thR4 *((sthR2-4*cthR2)*thR2 +sth*((22*cth+5)*th-7*sth) +2*(cth*(11*cth+4)+5))) +th_2*th_3*thR3 *(th*(sth*(th-5*sth)+3*(-cthM1)*cth) +5*(-cthM1)*sth)) +thR2*(th_2R2 *(th_3R2*(sthR2*thR2+4*(cth*(sth*th+cth)+1)) +thR2 *(cthR2*thR2+sth*(sth-(5*cth+1)*th) -2*(3*cthR2+1))) +th_1R2*th_3R2 *(th *((4*sthR2+(-cthM1)*cth)*th +(21*cth-5)*sth) +24*(cthR2+2))) +th_1R2*(2*(th_2R2+th_1R2)*th_3R2 *(sth*thR2*(cth*th-5*sth) +2 *(cthR2*thR2 -(9*cth*sth*th+8*(cthR2+1)))) -thR2 *(th_2R2 *(th *(th *(sth*(2*cth*th-11*sth)+cth*(9*cth+1)) -(63*cth+5)*sth) +6*(sthR2-2*(5*cthR2+1))) +th_1R2 *(th *(th *(2*sth*(cth*th-5*sth)+cth*(9*cth+1)) -(59*cth+5)*sth) +2*(3*sthR2-2*(14*cthR2+3)))))) +thR2*(th_1*th_2*thR3*((sthR2+(-cthM1)*cth)*th+3*cthM1*sth)+4*th_3R5) +th_3*(thR6*(sth*(cth*th-sth)+cth*(cth+2)+1) -th_3R2*(thR4*(cthM1*sth*th+2*(cthR2+3))+32*th_1R2*th_3R2)) );
      dh11dth = norm * ( -2*(th_2*(th_1*th_2 *(thR4*((sthR2-2*cthR2)*thR2+5*sth*(3*cth*th-sth)+16*cthR2) -(th_2R4+th_1R2*(2*th_2R2+th_1R2))*sththP2cth *(cth*thR2-(5*sth*th+8*cth))) +th_3*(th_1*(th_1*thR3*(th*(sth*(th-5*sth)+3*(-cthM1)*cth)+5*(-cthM1)*sth) -th_2*(th_2R2+th_1R2)*th_3 *(sth*thR2*(cth*th-5*sth) +2*(cthR2*thR2-(9*cth*sth*th+8*(cthR2+1))))) +thR5*((sthR2+cthM1*cth)*th+cthM1*sth)) +th_1*th_2*thR2 *((th_2R2+th_1R2)*(thR2*(cth*sth*th+5*(cthR2-sthR2)) +sth*(3*sth-32*cth*th)-32*cthR2) -th_3R2*(th*((2*sthR2+(-cthM1)*cth)*th+(13*cth-5)*sth)+16*(cthR2+1)))) +th_1*(thR4*(cthM1*th_3R2*(sth*th+2*cthM1)-sth*thR2*(cth*th-sth)) +16*th_2R2*th_3R4)) );
      dh12dth = norm * ( cthM1*sth*thR7+th_1*th_3 *(th_2*(2*(th_2R4+th_1R2*(2*th_2R2+th_1R2))*sththP2cth *(cth*thR2-(5*sth*th+8*cth)) -th_1R2*thR2 *(th*(th*(sth*(2*cth*th-9*sth)+cth*(9*cth+1)) -5*(11*cth+1)*sth) +2*(3*sthR2-2*(13*cthR2+3))) +2*(th_2R2+th_1R2)*th_3R2 *(sth*thR2*(cth*th-5*sth) +2*(cthR2*thR2-(9*cth*sth*th+8*(cthR2+1))))) +thR2*(th_2*th_3R2 *(th*((3*sthR2+(-cthM1)*cth)*th+(17*cth-5)*sth) +4*(5*cthR2+11)) -(th_2R3*(th*(th*(sth*(2*cth*th-9*sth)+cth*(9*cth+1)) -5*(11*cth+1)*sth) +2*(3*sthR2-2*(13*cthR2+3))) +th*(th_2*th *(th*(cth*(17*sth-3*cth*th)+sthR2*th) +2*(sth*(2*th-3*sth)+4*(cth*(2*cth+1)+1))) +th_1*th_3 *(th*(sth*(th-5*sth)+3*(-cthM1)*cth) +5*(-cthM1)*sth))))) -(th_1*(th_1*thR5*((sthR2+(-cthM1)*cth)*th+3*cthM1*sth) +32*th_2*th_3R5) +th_3R2*thR5*((sthR2+cthM1*cth)*th+cthM1*sth)) );
      dh22dth = norm * ( 2*th_1 *(cthM1R2*thR6+th_3R2*(thR2*(th_3R2*(sthR2*thR2+4*(cth*(sth*th+cth)+7)) +thR2*(cthR2*thR2+sth*(sth-2*(2*cth+1)*th) -4*(cthR2+3))) +(th_2R4+th_1R2*(2*th_2R2+th_1R2))*sththP2cth *(cth*thR2-(5*sth*th+8*cth)) +(th_2R2+th_1R2)*(th_3R2 *(sth*thR2*(cth*th-5*sth) +2*(cthR2*thR2 -(9*cth*sth*th+8*(cthR2+1)))) -thR2 *(thR2*(sth*(cth*th-4*sth)+cth*(4*cth+1)) +sth*(3*sth-(23*cth+5)*th) -4*(5*cthR2+3)))) -(cthM1*(th_2R2+th_1R2)*thR4*(sth*th+2*cthM1)+16*th_3R6)) );

      norm = HALF / (thR3 * thR3);

      dg00dth = norm * ( th_1*(th_3R2+th_2R2)*(thR2+(sth*th+4*cthM1)) );
      dg01dth = norm * ( th_2*(thR2*(thR2+2*cthM1)-2*th_1R2*(thR2+sth*th+4*cthM1))/2 );
      dg02dth = norm * ( th_3*(thR2*(thR2+2*cthM1)-2*th_1R2*(thR2+sth*th+4*cthM1))/2 );
      dg11dth = norm * ( th_1*(-(th_2R2+2*cthM1)*thR2+(th_3R2+th_1R2)*(sth*th+4*cthM1)) );
      dg12dth = norm * ( -th_1*th_2*th_3*((th+sth)*th+4*cthM1) );
      dg22dth = norm * ( th_2*(-(th_1R2+2*cthM1)*thR2 + (th_3R2+th_2R2)*(sth*th+4*cthM1)) );
   }

   double r =   (  dx_th_1*(dx_th_1*dg00dth + dx_th_2*dg01dth + dx_th_3*dg02dth)
                 + dx_th_2*(dx_th_1*dg01dth + dx_th_2*dg11dth + dx_th_3*dg12dth)
                 + dx_th_3*(dx_th_1*dg02dth + dx_th_2*dg12dth + dx_th_3*dg22dth) )

              + (  (A_2*dx_th_3 - A_3*dx_th_2)*(dh02dth*(A_1*dx_th_2 - A_2*dx_th_1) + dh01dth*(A_3*dx_th_1 - A_1*dx_th_3) + dh00dth*(A_2*dx_th_3 - A_3*dx_th_2))
                 + (A_3*dx_th_1 - A_1*dx_th_3)*(dh12dth*(A_1*dx_th_2 - A_2*dx_th_1) + dh11dth*(A_3*dx_th_1 - A_1*dx_th_3) + dh01dth*(A_2*dx_th_3 - A_3*dx_th_2))
                 + (A_1*dx_th_2 - A_2*dx_th_1)*(dh22dth*(A_1*dx_th_2 - A_2*dx_th_1) + dh12dth*(A_3*dx_th_1 - A_1*dx_th_3) + dh02dth*(A_2*dx_th_3 - A_3*dx_th_2)) );

   r +=   (  dy_th_1*(dy_th_1*dg00dth + dy_th_2*dg01dth + dy_th_3*dg02dth)
           + dy_th_2*(dy_th_1*dg01dth + dy_th_2*dg11dth + dy_th_3*dg12dth)
           + dy_th_3*(dy_th_1*dg02dth + dy_th_2*dg12dth + dy_th_3*dg22dth) )

        + (  (A_2*dy_th_3 - A_3*dy_th_2)*(dh02dth*(A_1*dy_th_2 - A_2*dy_th_1) + dh01dth*(A_3*dy_th_1 - A_1*dy_th_3) + dh00dth*(A_2*dy_th_3 - A_3*dy_th_2))
           + (A_3*dy_th_1 - A_1*dy_th_3)*(dh12dth*(A_1*dy_th_2 - A_2*dy_th_1) + dh11dth*(A_3*dy_th_1 - A_1*dy_th_3) + dh01dth*(A_2*dy_th_3 - A_3*dy_th_2))
           + (A_1*dy_th_2 - A_2*dy_th_1)*(dh22dth*(A_1*dy_th_2 - A_2*dy_th_1) + dh12dth*(A_3*dy_th_1 - A_1*dy_th_3) + dh02dth*(A_2*dy_th_3 - A_3*dy_th_2)) );

   r +=   (  dz_th_1*(dz_th_1*dg00dth + dz_th_2*dg01dth + dz_th_3*dg02dth)
           + dz_th_2*(dz_th_1*dg01dth + dz_th_2*dg11dth + dz_th_3*dg12dth)
           + dz_th_3*(dz_th_1*dg02dth + dz_th_2*dg12dth + dz_th_3*dg22dth) )

        + (  (A_2*dz_th_3 - A_3*dz_th_2)*(dh02dth*(A_1*dz_th_2 - A_2*dz_th_1) + dh01dth*(A_3*dz_th_1 - A_1*dz_th_3) + dh00dth*(A_2*dz_th_3 - A_3*dz_th_2))
           + (A_3*dz_th_1 - A_1*dz_th_3)*(dh12dth*(A_1*dz_th_2 - A_2*dz_th_1) + dh11dth*(A_3*dz_th_1 - A_1*dz_th_3) + dh01dth*(A_2*dz_th_3 - A_3*dz_th_2))
           + (A_1*dz_th_2 - A_2*dz_th_1)*(dh22dth*(A_1*dz_th_2 - A_2*dz_th_1) + dh12dth*(A_3*dz_th_1 - A_1*dz_th_3) + dh02dth*(A_2*dz_th_3 - A_3*dz_th_2)) );

   drdth1 -= r;
   
   // d rho / d th_2

   if (thR2 < SMALL_THETA_SQUARED)
   {
      dh00dth = -2*th_2/3+th_3*th_1/6+(16*th_2*th_1R2+32*th_2R3+17*th_3R2*th_2)/180;
      dh01dth = th_1/3+th_3*th_2/6-(8*th_1R3+(24*th_2R2+13*th_3R2)*th_1)/180;
      dh02dth = 1/2-(3*th_1R2+9*th_2R2+th_3R2)/24+th_3*th_2*th_1/15;
      dh11dth = -th_3*th_1/6+(16*th_2*th_1R2-9*th_3R2*th_2)/180;
      dh12dth = -th_3/12+th_2*th_1/4+(6*th_3*th_1R2+18*th_3*th_2R2+th_3R3)/180;
      dh22dth = th_2/2-(6*th_2*th_1R2+6*th_2R3+th_3R2*th_2)/36;

      dg00dth = -th_2/24+(th_2*th_1R2+2*th_2R3+2*th_3R2*th_2)/720;
      dg01dth = th_1/48-(th_1R3+(3*th_2R2+th_3R2)*th_1)/1440;
      dg02dth = -th_3*th_2*th_1/720;
      dg11dth = (th_2*th_1R2+th_3R2*th_2)/720;
      dg12dth = th_3/48-(th_3*th_1R2+3*th_3*th_2R2+th_3R3)/1440;
      dg22dth = -th_2/24+(2*th_2*th_1R2+2*th_2R3+th_3R2*th_2)/720;
   }
   else
   {
      double norm = 1 / (thR3 * thR2 * thR3 * thR2);

      dh00dth = norm * ( -2*(th_1*(th_1*th_2 *(thR4*((sthR2-2*cthR2)*thR2+5*sth*(3*cth*th-sth)+16*cthR2) -2*th_1R2*th_2R2*sththP2cth*(cth*thR2-(5*sth*th+8*cth)) +thR2*((th_2R2+th_1R2)*(thR2*(cth*sth*th+5*(cthR2-sthR2)) +sth*(3*sth-32*cth*th)-32*cthR2) -th_3R2*(th*((2*sthR2+(-cthM1)*cth)*th+(13*cth-5)*sth) +16*(cthR2+1)))) -(th_1*th_2*(th_2R4+th_1R4)*sththP2cth*(cth*thR2-(5*sth*th+8*cth)) +th_3*(th_2*(th_2*thR3*(th*(sth*(th-5*sth)+3*(-cthM1)*cth)+5*(-cthM1)*sth) +th_1*(th_2R2+th_1R2)*th_3 *(sth*thR2*(cth*th-5*sth) +2*(cthR2*thR2-(9*cth*sth*th+8*(cthR2+1))))) +thR5*((sthR2+cthM1*cth)*th+cthM1*sth)))) +th_2*(thR4*(cthM1*th_3R2*(sth*th+2*cthM1)-sth*thR2*(cth*th-sth)) +16*th_1R2*th_3R4)) );
      dh01dth = norm * ( th_1*(th_1R4*thR2*sththP2cthR2+th_2R2*(2*(th_2R2+th_1R2)*th_3R2 *(sth*thR2*(cth*th-5*sth) +2 *(cthR2*thR2 -(9*cth*sth*th+8*(cthR2+1)))) -th_2R2*thR2 *(sth*thR2*(2*cth*th-11*sth) +2 *(5*cthR2*thR2 +sth*(3*sth-34*cth*th) -34*cthR2)))) +thR2*(th*(th_1R2*(th_1*th*(cthR2*thR2+sth*(sth-6*cth*th)-8*cthR2) -th_2*th_3*(th*(sth*(th-5*sth)+3*(-cthM1)*cth)+5*(-cthM1)*sth)) +th_2R3*th_3*(th*(sth*(th-5*sth)+3*(-cthM1)*cth)+5*(-cthM1)*sth)) +th_1*(th_3R2*(th_2R2*(5*sthR2*thR2+2*(th*((-cthM1)*cth*th+5*(3*cth-1)*sth) +18*(cthR2+1))) +th_1R2*(sthR2*thR2+4*(cth*(sth*th+cth)+1))) -th_2R2*(thR2*((2*sthR2-5*cthR2)*thR2+sth*(36*cth*th-11*sth) +40*cthR2) +2*th_1R2 *(thR2*(sth*(cth*th-6*sth)+5*cthR2) +3*(sth*(sth-12*cth*th)-12*cthR2))))) +2*(th_1*(th_3R2*(2*th_3R2*thR2-(thR4*(cthM1*sth*th+2*(cthR2+1))+16*th_2R2*th_3R2)) +thR6*(sth*(cth*th-sth)+2*cthR2)) +th_2*(th_1*th_2*(th_2R4+th_1R2*(2*th_2R2+th_1R2))*sththP2cth *(cth*thR2-(5*sth*th+8*cth)) +th_3*thR5*((sthR2+cthM1*cth)*th+cthM1*sth))) );
      dh02dth = norm * ( -(cthM1*sth*thR7+th_3*(th_1*th_2 *(th_2R2*(thR2*(th *(th *(sth*(2*cth*th-9*sth)+cth*(9*cth+1)) -5*(11*cth+1)*sth) +2*(3*sthR2-2*(13*cthR2+3))) -2*th_2R2*sththP2cth *(cth*thR2-(5*sth*th+8*cth))) -2*th_1R4*sththP2cth*(cth*thR2-(5*sth*th+8*cth)) +thR2*(th_1R2*(th *(th *(sth*(2*cth*th-9*sth)+cth*(9*cth+1)) -5*(11*cth+1)*sth) +2*(3*sthR2-2*(13*cthR2+3))) +thR2*(th*(cth*(17*sth-3*cth*th)+sthR2*th) +2 *(sth*(2*th-3*sth) +4*(cth*(2*cth+1)+1))))) -(4*th_1R3*th_2R3*sththP2cth*(cth*thR2-(5*sth*th+8*cth)) +th_3*(th_1*th_2*th_3 *(thR2*(th *((3*sthR2+(-cthM1)*cth)*th +(17*cth-5)*sth) +4*(5*cthR2+11)) +2*(th_2R2+th_1R2) *(sth*thR2*(cth*th-5*sth) +2*(cthR2*thR2-(9*cth*sth*th+8*(cthR2+1))))) +thR3*(th_2R2*(th*(sth*(th-5*sth)+3*(-cthM1)*cth) +5*(-cthM1)*sth) +thR2*((sthR2+cthM1*cth)*th +cthM1*sth))))) +th_2*(32*th_1*th_3R5-th_2*thR5 *((sthR2+(-cthM1)*cth)*th+3*cthM1*sth))) );
      dh11dth = norm * ( 2*(th_2*(thR2*(th_1R4*sththP2cthR2+thR2*(th_1R2 *(cthR2*thR2+sth*(sth-6*cth*th) -8*cthR2) -th_2R2 *(sthR2*thR2 +3 *(cth*th*(7*sth-cth*th) +2*(4*cthR2-sthR2)))) +th_1R2*th_3R2 *(sthR2*thR2+4*(cth*(sth*th+cth)+1))) +th_2R2*(th_2R4+th_1R2*(2*th_2R2+th_1R2))*sththP2cth *(cth*thR2-(5*sth*th+8*cth)) +thR2*(th_2R2*(th_3R2*(th*(cth*((-cthM1)*th+17*sth)+3*sthR2*th) +5*(4*(cthR2+1)-sth*th)) -th_1R2*(thR2*(sth*(cth*th-7*sth)+5*cthR2) +sth*(3*sth-40*cth*th)-40*cthR2)) +th*(thR3*(3*sth*(cth*th-sth)+4*cthR2) -th_1*th_2*th_3*(th*(sth*(th-5*sth)+3*(-cthM1)*cth)+5*(-cthM1)*sth)) -th_2R4*(thR2*(sth*(cth*th-6*sth)+5*cthR2) +3*(sth*(sth-12*cth*th)-12*cthR2))+4*th_3R4)) -th_3*(th_2*th_3 *(-th_2R2*(th_2R2+th_1R2) *(sth*thR2*(cth*th-5*sth) +2*(cthR2*thR2-(9*cth*sth*th+8*(cthR2+1)))) +thR4*(3*cthM1*sth*th+2*(cth*(3*cth-2)+3))+16*th_2R2*th_3R2) +th_1*thR5*((sthR2+cthM1*cth)*th+cthM1*sth))) );
      dh12dth = norm * ( th_3*(th_1R4*thR2*sththP2cthR2+th_2*(th_2*(2*(th_2R4+th_1R2*(2*th_2R2+th_1R2)) *sththP2cth *(cth*thR2-(5*sth*th+8*cth)) -thR4 *((sthR2-4*cthR2)*thR2 +sth*((22*cth+5)*th-7*sth) +2*(cth*(11*cth+4)+5))) -th_1*th_3*thR3 *(th*(sth*(th-5*sth)+3*(-cthM1)*cth) +5*(-cthM1)*sth)) +thR2*(th_1R2 *(th_3R2*(sthR2*thR2+4*(cth*(sth*th+cth)+1)) +thR2 *(cthR2*thR2+sth*(sth-(5*cth+1)*th) -2*(3*cthR2+1))) +th_2R2*th_3R2 *(th *((4*sthR2+(-cthM1)*cth)*th +(21*cth-5)*sth) +24*(cthR2+2))) +th_2R2*(2*(th_2R2+th_1R2)*th_3R2 *(sth*thR2*(cth*th-5*sth) +2 *(cthR2*thR2 -(9*cth*sth*th+8*(cthR2+1)))) -thR2 *(th_1R2 *(th *(th *(sth*(2*cth*th-11*sth)+cth*(9*cth+1)) -(63*cth+5)*sth) +6*(sthR2-2*(5*cthR2+1))) +th_2R2 *(th *(th *(2*sth*(cth*th-5*sth)+cth*(9*cth+1)) -(59*cth+5)*sth) +2*(3*sthR2-2*(14*cthR2+3)))))) +thR2*(4*th_3R5-th_1*th_2*thR3*((sthR2+(-cthM1)*cth)*th+3*cthM1*sth)) +th_3*(thR6*(sth*(cth*th-sth)+cth*(cth+2)+1) -th_3R2*(thR4*(cthM1*sth*th+2*(cthR2+3))+32*th_2R2*th_3R2)) );
      dh22dth = norm * ( 2*th_2 *(cthM1R2*thR6+th_3R2*(thR2*(th_3R2*(sthR2*thR2+4*(cth*(sth*th+cth)+7)) +thR2*(cthR2*thR2+sth*(sth-2*(2*cth+1)*th) -4*(cthR2+3))) +(th_2R4+th_1R2*(2*th_2R2+th_1R2))*sththP2cth *(cth*thR2-(5*sth*th+8*cth)) +(th_2R2+th_1R2)*(th_3R2 *(sth*thR2*(cth*th-5*sth) +2*(cthR2*thR2 -(9*cth*sth*th+8*(cthR2+1)))) -thR2 *(thR2*(sth*(cth*th-4*sth)+cth*(4*cth+1)) +sth*(3*sth-(23*cth+5)*th) -4*(5*cthR2+3)))) -(cthM1*(th_2R2+th_1R2)*thR4*(sth*th+2*cthM1)+16*th_3R6)) );

      norm = HALF / (thR3 * thR3);

      dg00dth = norm * ( th_2*(-(th_1R2+2*cthM1)*thR2+(th_3R2+th_2R2)*(sth*th+4*cthM1)) );
      dg01dth = norm * ( th_1*(thR2*(thR2+2*cthM1)-2*th_2R2*(thR2+sth*th+4*cthM1))/2 );
      dg02dth = norm * ( -th_1*th_2*th_3*(thR2+sth*th+4*cthM1) );
      dg11dth = norm * ( th_2*(th_3R2+th_1R2)*(thR2+(sth*th+4*cthM1)) );
      dg12dth = norm * ( th_3*(thR2*(thR2+2*cthM1)-2*th_2R2*(th*(th+sth)+4*cthM1))/2 );
      dg22dth = norm * ( th_2*(-(th_3R2+2*cthM1)*thR2+(th_2R2+th_1R2)*(sth*th+4*cthM1)) );
   }
   
   r =  (  dx_th_1*(dx_th_1*dg00dth + dx_th_2*dg01dth + dx_th_3*dg02dth)
         + dx_th_2*(dx_th_1*dg01dth + dx_th_2*dg11dth + dx_th_3*dg12dth)
         + dx_th_3*(dx_th_1*dg02dth + dx_th_2*dg12dth + dx_th_3*dg22dth) )

      + (  (A_2*dx_th_3 - A_3*dx_th_2)*(dh02dth*(A_1*dx_th_2 - A_2*dx_th_1) + dh01dth*(A_3*dx_th_1 - A_1*dx_th_3) + dh00dth*(A_2*dx_th_3 - A_3*dx_th_2))
         + (A_3*dx_th_1 - A_1*dx_th_3)*(dh12dth*(A_1*dx_th_2 - A_2*dx_th_1) + dh11dth*(A_3*dx_th_1 - A_1*dx_th_3) + dh01dth*(A_2*dx_th_3 - A_3*dx_th_2))
         + (A_1*dx_th_2 - A_2*dx_th_1)*(dh22dth*(A_1*dx_th_2 - A_2*dx_th_1) + dh12dth*(A_3*dx_th_1 - A_1*dx_th_3) + dh02dth*(A_2*dx_th_3 - A_3*dx_th_2)) );

   r +=  (  dy_th_1*(dy_th_1*dg00dth + dy_th_2*dg01dth + dy_th_3*dg02dth)
          + dy_th_2*(dy_th_1*dg01dth + dy_th_2*dg11dth + dy_th_3*dg12dth)
          + dy_th_3*(dy_th_1*dg02dth + dy_th_2*dg12dth + dy_th_3*dg22dth) )

       + (  (A_2*dy_th_3 - A_3*dy_th_2)*(dh02dth*(A_1*dy_th_2 - A_2*dy_th_1) + dh01dth*(A_3*dy_th_1 - A_1*dy_th_3) + dh00dth*(A_2*dy_th_3 - A_3*dy_th_2))
          + (A_3*dy_th_1 - A_1*dy_th_3)*(dh12dth*(A_1*dy_th_2 - A_2*dy_th_1) + dh11dth*(A_3*dy_th_1 - A_1*dy_th_3) + dh01dth*(A_2*dy_th_3 - A_3*dy_th_2))
          + (A_1*dy_th_2 - A_2*dy_th_1)*(dh22dth*(A_1*dy_th_2 - A_2*dy_th_1) + dh12dth*(A_3*dy_th_1 - A_1*dy_th_3) + dh02dth*(A_2*dy_th_3 - A_3*dy_th_2)) );

   r +=  (  dz_th_1*(dz_th_1*dg00dth + dz_th_2*dg01dth + dz_th_3*dg02dth)
          + dz_th_2*(dz_th_1*dg01dth + dz_th_2*dg11dth + dz_th_3*dg12dth)
          + dz_th_3*(dz_th_1*dg02dth + dz_th_2*dg12dth + dz_th_3*dg22dth) )

       + (  (A_2*dz_th_3 - A_3*dz_th_2)*(dh02dth*(A_1*dz_th_2 - A_2*dz_th_1) + dh01dth*(A_3*dz_th_1 - A_1*dz_th_3) + dh00dth*(A_2*dz_th_3 - A_3*dz_th_2))
          + (A_3*dz_th_1 - A_1*dz_th_3)*(dh12dth*(A_1*dz_th_2 - A_2*dz_th_1) + dh11dth*(A_3*dz_th_1 - A_1*dz_th_3) + dh01dth*(A_2*dz_th_3 - A_3*dz_th_2))
          + (A_1*dz_th_2 - A_2*dz_th_1)*(dh22dth*(A_1*dz_th_2 - A_2*dz_th_1) + dh12dth*(A_3*dz_th_1 - A_1*dz_th_3) + dh02dth*(A_2*dz_th_3 - A_3*dz_th_2)) );

   drdth2 -= r;

   // d rho / d th_3

   if (thR2 < SMALL_THETA_SQUARED)
   {
      dh00dth = -th_3/6+th_2*th_1/6-(9*th_3*th_1R2-17*th_3*th_2R2-2*th_3R3)/180;
      dh01dth = -(th_1R2-th_2R2)/12-13*th_3*th_2*th_1/90;
      dh02dth = -th_1/12-th_3*th_2/12+(2*th_1R3+(2*th_2R2+th_3R2)*th_1)/60;
      dh11dth = -th_3/6-th_2*th_1/6+(17*th_3*th_1R2-9*th_3*th_2R2+2*th_3R3)/180;
      dh12dth = -th_2/12+th_3*th_1/12+(2*th_2*th_1R2+2*th_2R3+th_3R2*th_2)/60;
      dh22dth = +(-(th_3*th_1R2+th_3*th_2R2)/36);

      dg00dth = -th_3/24+(th_3*th_1R2+2*th_3*th_2R2+2*th_3R3)/720;
      dg01dth = -th_3*th_2*th_1/720;
      dg02dth = th_1/48-(th_1R3+(th_2R2+3*th_3R2)*th_1)/1440;
      dg11dth = -th_3/24+(2*th_3*th_1R2+th_3*th_2R2+2*th_3R3)/720;
      dg12dth = th_2/48-(th_2*th_1R2+th_2R3+3*th_3R2*th_2)/1440;
      dg22dth = (th_3*th_1R2+th_3*th_2R2)/720;
   }
   else
   {
      double norm = 1 / (thR3 * thR2 * thR3 * thR2);

      dh00dth = norm * ( -2*(th_3R3*(cthM1*thR4*(sth*th+2*cthM1)+16*th_1R2*th_3R2) -(th_1*(th_3*(th_1R5*sththP2cth*(cth*thR2-(5*sth*th+8*cth)) +th_3*(th_2*thR3*(th*(sth*(th-5*sth)+3*(-cthM1)*cth)+5*(-cthM1)*sth) +th_1*(th_2R2+th_1R2)*th_3 *(sth*thR2*(cth*th-5*sth) +2*(cthR2*thR2-(9*cth*sth*th+8*(cthR2+1))))) -th_1*th_2R2*thR2 *(thR2*(sth*(cth*th-4*sth)+5*cthR2) +sth*(3*sth-28*cth*th)-4*(7*cthR2+1))) +thR2*(th_2*thR3*((sthR2+cthM1*cth)*th+cthM1*sth) -th_1R3*th_3 *(thR2*(sth*(cth*th-4*sth)+5*cthR2) +sth*(3*sth-28*cth*th)-4*(7*cthR2+1)))) +th_3*(th_1R2*th_2R2*(th_2R2+2*th_1R2)*sththP2cth*(cth*thR2-(5*sth*th+8*cth)) +thR2*(th_1R2*(th_3R2*(th*((sthR2+(-cthM1)*cth)*th+(9*cth-5)*sth) +4*(3*cthR2+5)) -thR2*(th*(cth*(11*sth-cth*th)+sthR2*th) +2*(2*(-sthR2+3*cthR2+1)-sth*th))) +thR4*(sth*(cth*th-sth)+(cth-2)*cth+1))))) );
      dh01dth = norm * ( -(th_3*(thR2*(th*(2*th_1*th_2*th *(th*(cth*(11*sth-cth*th)+sthR2*th) +2*(2*(-sthR2+3*cthR2+1)-sth*th)) +(th_1R2-th_2R2)*th_3 *(th*(sth*(th-5*sth)+3*(-cthM1)*cth)+5*(-cthM1)*sth)) +2*th_1*th_2*(th_2R2+th_1R2) *(thR2*(sth*(cth*th-4*sth)+5*cthR2) +sth*(3*sth-28*cth*th)-4*(7*cthR2+1))) -2*th_1*th_2*(th_2R4+th_1R2*(2*th_2R2+th_1R2))*sththP2cth *(cth*thR2-(5*sth*th+8*cth))) -th_2*(2*th_1*th_3R3 *(thR2*(th*((sthR2+(-cthM1)*cth)*th+(9*cth-5)*sth)+4*(3*cthR2+5)) +(th_2R2+th_1R2)*(sth*thR2*(cth*th-5*sth) +2*(cthR2*thR2-(9*cth*sth*th+8*(cthR2+1))))) +th_2*thR5*((sthR2+cthM1*cth)*th+cthM1*sth)) +th_1*(th_1*thR5*((sthR2+cthM1*cth)*th+cthM1*sth)+32*th_2*th_3R5)) );
      dh02dth = norm * ( th_1*((th_2R4+th_1R2*(2*th_2R2+th_1R2))*thR2*sththP2cthR2 +th_3R2*(2*(th_2R2+th_1R2)*th_3R2 *(sth*thR2*(cth*th-5*sth)+2*(cthR2*thR2 -(9*cth*sth*th+8*(cthR2+1)))) -(th_2R2+th_1R2)*thR2 *(th*(th*(2*sth*(cth*th-4*sth)+cth*(9*cth+1)) -(51*cth+5)*sth) +6*(sthR2-4*(2*cthR2+1))))) +thR2*(th*(th_1*(th_2R2+th_1R2)*th *(cthR2*thR2+sth*(sth-(5*cth+1)*th)-2*(3*cthR2+1)) +th_3R2*(th_2*th_3*(th*(sth*(th-5*sth)+3*(-cthM1)*cth)+5*(-cthM1)*sth) -th_1*th *(th*((sthR2-cthR2)*th+(10*cth-1)*sth) +2*(-2*sthR2+cth*(5*cth+4)+11)))) +th_1*th_3R4*(th*((sthR2+(-cthM1)*cth)*th+(9*cth-5)*sth)+4*(3*cthR2+14))) +th_3*(2*th_1*(th_2R4+th_1R2*(2*th_2R2+th_1R2))*th_3*sththP2cth *(cth*thR2-(5*sth*th+8*cth)) +th_2*thR5*((3*sthR2+cthM1*cth)*th+5*cthM1*sth)) +th_1*(thR6*(sth*(cth*th-sth)+cth*(cth+2)+1)-32*th_3R6) );
      dh11dth = norm * ( -2*(th_3*(thR2*(th_2*th *(th_2*th *(th*(cth*(11*sth-cth*th)+sthR2*th) +2*(2*(-sthR2+3*cthR2+1)-sth*th)) +th_1*th_3*(th*(sth*(th-5*sth)+3*(-cthM1)*cth)+5*(-cthM1)*sth)) -(th_2R2*th_3R2 *(th*((sthR2+(-cthM1)*cth)*th+(9*cth-5)*sth)+4*(3*cthR2+5)) +thR4*(sth*(cth*th-sth)+(cth-2)*cth+1)) +th_2R2*(th_2R2+th_1R2) *(thR2*(sth*(cth*th-4*sth)+5*cthR2) +sth*(3*sth-28*cth*th)-4*(7*cthR2+1))) -th_2R2*(th_2R4+th_1R2*(2*th_2R2+th_1R2))*sththP2cth *(cth*thR2-(5*sth*th+8*cth))) +th_2*(th_1*(thR5*((sthR2+cthM1*cth)*th+cthM1*sth) -th_1*th_2*th_3R3 *(sth*thR2*(cth*th-5*sth)+2 *(cthR2*thR2 -(9*cth*sth*th+8*(cthR2+1))))) -th_2R3*th_3R3 *(sth*thR2*(cth*th-5*sth)+2*(cthR2*thR2 -(9*cth*sth*th+8*(cthR2+1))))) +th_3R3*(cthM1*thR4*(sth*th+2*cthM1)+16*th_2R2*th_3R2)) );
      dh12dth = norm * ( th_2*((th_2R4+th_1R2*(2*th_2R2+th_1R2))*thR2*sththP2cthR2 +th_3R2*(2*(th_2R2+th_1R2)*th_3R2 *(sth*thR2*(cth*th-5*sth)+2*(cthR2*thR2 -(9*cth*sth*th+8*(cthR2+1)))) -(th_2R2+th_1R2)*thR2 *(th*(th*(2*sth*(cth*th-4*sth)+cth*(9*cth+1)) -(51*cth+5)*sth) +6*(sthR2-4*(2*cthR2+1))))) +thR2*(th*(th_2*(th_2R2+th_1R2)*th *(cthR2*thR2+sth*(sth-(5*cth+1)*th)-2*(3*cthR2+1)) -th_3R2*(th_1*th_3*(th*(sth*(th-5*sth)+3*(-cthM1)*cth)+5*(-cthM1)*sth) +th_2*th *(th*((sthR2-cthR2)*th+(10*cth-1)*sth) +2*(-2*sthR2+cth*(5*cth+4)+11)))) +th_2*th_3R4*(th*((sthR2+(-cthM1)*cth)*th+(9*cth-5)*sth)+4*(3*cthR2+14))) +th_3*(th_1*(2*th_1*th_2*(2*th_2R2+th_1R2)*th_3*sththP2cth*(cth*thR2-(5*sth*th+8*cth)) -thR5*((3*sthR2+cthM1*cth)*th+5*cthM1*sth)) +2*th_2R5*th_3*sththP2cth*(cth*thR2-(5*sth*th+8*cth))) +th_2*(thR6*(sth*(cth*th-sth)+cth*(cth+2)+1)-32*th_3R6) );
      dh22dth = norm * ( 2*th_3 *(thR2*((th_2R4+th_1R2*(2*th_2R2+th_1R2))*sththP2cthR2 +thR2*((th_2R2+th_1R2)*(cthR2*thR2+sth*(sth-(5*cth+1)*th) +2*((2-3*cth)*cth-3)) +4*(thR2-6*th_3R2))) +th_3R2*((th_2R4+th_1R2*(2*th_2R2+th_1R2))*sththP2cth*(cth*thR2-(5*sth*th+8*cth)) +thR2*(36*th_3R2-(th_2R2+th_1R2)*(thR2*(sth*(cth*th-4*sth)+cth*(4*cth+1)) +sth*(3*sth-(23*cth+5)*th)-20*(cthR2+1))) +th_3R2*((th_2R2+th_1R2)*(sth*thR2*(cth*th-5*sth) +2*(cthR2*thR2-(9*cth*sth*th+8*(cthR2+1)))) -16*th_3R2))) );

      norm = HALF / (thR3 * thR3);

      dg00dth = norm * ( th_3*(-(th_1R2+2*cthM1)*thR2+(th_3R2+th_2R2)*(sth*th+4*cthM1)) );
      dg01dth = norm * ( -th_1*th_2*th_3*((th+sth)*th+4*cthM1) );
      dg02dth = norm * ( th_1*(thR2*(thR2+2*cthM1)-2*th_3R2*((th+sth)*th+4*cthM1))/2 );
      dg11dth = norm * ( th_3*(-(th_2R2+2*cthM1)*thR2+(th_3R2+th_1R2)*(sth*th+4*cthM1)) );
      dg12dth = norm * ( th_2*(thR2*(thR2+2*cthM1)-2*th_3R2*(thR2+sth*th+4*cthM1))/2 );
      dg22dth = norm * ( th_3*(th_2R2+th_1R2)*(thR2+(sth*th+4*cthM1)) );
   }
   
   r =   (  dx_th_1*(dx_th_1*dg00dth + dx_th_2*dg01dth + dx_th_3*dg02dth)
          + dx_th_2*(dx_th_1*dg01dth + dx_th_2*dg11dth + dx_th_3*dg12dth)
          + dx_th_3*(dx_th_1*dg02dth + dx_th_2*dg12dth + dx_th_3*dg22dth) )

       + (  (A_2*dx_th_3 - A_3*dx_th_2)*(dh02dth*(A_1*dx_th_2 - A_2*dx_th_1) + dh01dth*(A_3*dx_th_1 - A_1*dx_th_3) + dh00dth*(A_2*dx_th_3 - A_3*dx_th_2))
          + (A_3*dx_th_1 - A_1*dx_th_3)*(dh12dth*(A_1*dx_th_2 - A_2*dx_th_1) + dh11dth*(A_3*dx_th_1 - A_1*dx_th_3) + dh01dth*(A_2*dx_th_3 - A_3*dx_th_2))
          + (A_1*dx_th_2 - A_2*dx_th_1)*(dh22dth*(A_1*dx_th_2 - A_2*dx_th_1) + dh12dth*(A_3*dx_th_1 - A_1*dx_th_3) + dh02dth*(A_2*dx_th_3 - A_3*dx_th_2)) );

   r +=   (  dy_th_1*(dy_th_1*dg00dth + dy_th_2*dg01dth + dy_th_3*dg02dth)
           + dy_th_2*(dy_th_1*dg01dth + dy_th_2*dg11dth + dy_th_3*dg12dth)
           + dy_th_3*(dy_th_1*dg02dth + dy_th_2*dg12dth + dy_th_3*dg22dth) )

        + (  (A_2*dy_th_3 - A_3*dy_th_2)*(dh02dth*(A_1*dy_th_2 - A_2*dy_th_1) + dh01dth*(A_3*dy_th_1 - A_1*dy_th_3) + dh00dth*(A_2*dy_th_3 - A_3*dy_th_2))
           + (A_3*dy_th_1 - A_1*dy_th_3)*(dh12dth*(A_1*dy_th_2 - A_2*dy_th_1) + dh11dth*(A_3*dy_th_1 - A_1*dy_th_3) + dh01dth*(A_2*dy_th_3 - A_3*dy_th_2))
           + (A_1*dy_th_2 - A_2*dy_th_1)*(dh22dth*(A_1*dy_th_2 - A_2*dy_th_1) + dh12dth*(A_3*dy_th_1 - A_1*dy_th_3) + dh02dth*(A_2*dy_th_3 - A_3*dy_th_2)) );

   r +=   (  dz_th_1*(dz_th_1*dg00dth + dz_th_2*dg01dth + dz_th_3*dg02dth)
           + dz_th_2*(dz_th_1*dg01dth + dz_th_2*dg11dth + dz_th_3*dg12dth)
           + dz_th_3*(dz_th_1*dg02dth + dz_th_2*dg12dth + dz_th_3*dg22dth) )

        + (  (A_2*dz_th_3 - A_3*dz_th_2)*(dh02dth*(A_1*dz_th_2 - A_2*dz_th_1) + dh01dth*(A_3*dz_th_1 - A_1*dz_th_3) + dh00dth*(A_2*dz_th_3 - A_3*dz_th_2))
           + (A_3*dz_th_1 - A_1*dz_th_3)*(dh12dth*(A_1*dz_th_2 - A_2*dz_th_1) + dh11dth*(A_3*dz_th_1 - A_1*dz_th_3) + dh01dth*(A_2*dz_th_3 - A_3*dz_th_2))
           + (A_1*dz_th_2 - A_2*dz_th_1)*(dh22dth*(A_1*dz_th_2 - A_2*dz_th_1) + dh12dth*(A_3*dz_th_1 - A_1*dz_th_3) + dh02dth*(A_2*dz_th_3 - A_3*dz_th_2)) );

   drdth3 -= r;
}

void TAux::subtract_offset_drho_dth(int di, int dj, int dk, TFields *fields, double &drdth1, double &drdth2, double &drdth3)
/* (di, dj, dk)  is position in lattice units of the cell where th is varied, relative to the cell managed by this TAux object and the TFields argument. 
     Variations in the fields of the latter, central cell are NOT handled by this method!
     IMPORTANT: assumes that h, G and theta field derivatives are ready for use! */
{
   switch (di)
   {
     case -2: switch (dj)
              {                         
                case -1: switch (dk)
                         {
                             case -1: drdth1 -= (dz_th_3*G_02+dy_th_3*G_02+dz_th_2*G_01+dy_th_2*G_01+dz_th_1*G_00+dy_th_1*G_00+fields->A_3*((dy_th_3*h_02+dy_th_2*h_01+dy_th_1*h_00)*fields->A_3-(dz_th_3*h_02+dy_th_3*h_02+dz_th_2*h_01+dy_th_2*h_01+dz_th_1*h_00+dy_th_1*h_00)*fields->A_2-(dx_th_3*h_02+dx_th_2*h_01+dx_th_1*h_00)*fields->A_1)+(dz_th_3*h_02+dz_th_2*h_01+dz_th_1*h_00)*fields->A_2*fields->A_2+fields->A_1*((dz_th_3*h_02+dy_th_3*h_02+dz_th_2*h_01+dy_th_2*h_01+dz_th_1*h_00+dy_th_1*h_00)*fields->A_1-(dx_th_3*h_02+dx_th_2*h_01+dx_th_1*h_00)*fields->A_2))/120;
                                      drdth2 -= (dz_th_3*G_12+dy_th_3*G_12+dz_th_2*(G_11+G_01)+dy_th_2*(G_11+G_01)+fields->A_3*((dy_th_3*h_12+dy_th_2*(h_11+h_01))*fields->A_3-(dz_th_3*h_12+dy_th_3*h_12+dz_th_2*(h_11+h_01)+dy_th_2*(h_11+h_01))*fields->A_2-(dx_th_3*h_12+dx_th_2*(h_11+h_01))*fields->A_1)+(dz_th_3*h_12+dz_th_2*(h_11+h_01))*fields->A_2*fields->A_2+fields->A_1*((dz_th_3*h_12+dy_th_3*h_12+dz_th_2*(h_11+h_01)+dy_th_2*(h_11+h_01))*fields->A_1-(dx_th_3*h_12+dx_th_2*(h_11+h_01))*fields->A_2))/120;
                                      drdth3 -= ((dz_th_3+dy_th_3)*(G_22+G_12+G_02+(h_22+h_12+h_02)*fields->A_1*fields->A_1)+(h_22+h_12+h_02)*(fields->A_3*(dy_th_3*fields->A_3-((dz_th_3+dy_th_3)*fields->A_2+dx_th_3*fields->A_1))+fields->A_2*(dz_th_3*fields->A_2-dx_th_3*fields->A_1)))/120;
                                      break;
                             case 0:  drdth1 -= (dx_th_3*G_02+dx_th_2*G_01+dx_th_1*G_00+fields->A_3*((dx_th_3*h_02+dx_th_2*h_01+dx_th_1*h_00)*fields->A_3-(dz_th_3*h_02+dz_th_2*h_01+dz_th_1*h_00)*fields->A_1)+fields->A_2*((dx_th_3*h_02+dx_th_2*h_01+dx_th_1*h_00)*fields->A_2-(dy_th_3*h_02+dy_th_2*h_01+dy_th_1*h_00)*fields->A_1))/30;
                                      drdth2 -= (dx_th_3*G_12+dx_th_2*(G_11+G_01)+fields->A_3*((dx_th_3*h_12+dx_th_2*(h_11+h_01))*fields->A_3-(dz_th_3*h_12+dz_th_2*(h_11+h_01))*fields->A_1)+fields->A_2*((dx_th_3*h_12+dx_th_2*(h_11+h_01))*fields->A_2-(dy_th_3*h_12+dy_th_2*(h_11+h_01))*fields->A_1))/30;
                                      drdth3 -= (dx_th_3*(G_22+G_12+G_02)+(h_22+h_12+h_02)*(fields->A_3*(dx_th_3*fields->A_3-dz_th_3*fields->A_1)+fields->A_2*(dx_th_3*fields->A_2-dy_th_3*fields->A_1)))/30;
                                      break;
                             case 1:  drdth1 -= (-dz_th_3*G_02+dy_th_3*G_02-dz_th_2*G_01+dy_th_2*G_01-dz_th_1*G_00+dy_th_1*G_00+fields->A_3*((dy_th_3*h_02+dy_th_2*h_01+dy_th_1*h_00)*fields->A_3-(dz_th_3*h_02-dy_th_3*h_02+dz_th_2*h_01-dy_th_2*h_01+dz_th_1*h_00-dy_th_1*h_00)*fields->A_2+(dx_th_3*h_02+dx_th_2*h_01+dx_th_1*h_00)*fields->A_1)-(dz_th_3*h_02+dz_th_2*h_01+dz_th_1*h_00)*fields->A_2*fields->A_2-fields->A_1*((dx_th_3*h_02+dx_th_2*h_01+dx_th_1*h_00)*fields->A_2+(dz_th_3*h_02-dy_th_3*h_02+dz_th_2*h_01-dy_th_2*h_01+dz_th_1*h_00-dy_th_1*h_00)*fields->A_1))/120;
                                      drdth2 -= (-dz_th_3*G_12+dy_th_3*G_12-dz_th_2*(G_11+G_01)+dy_th_2*(G_11+G_01)+fields->A_3*((dy_th_3*h_12+dy_th_2*(h_11+h_01))*fields->A_3-(dz_th_3*h_12-dy_th_3*h_12+dz_th_2*(h_11+h_01)-dy_th_2*(h_11+h_01))*fields->A_2+(dx_th_3*h_12+dx_th_2*(h_11+h_01))*fields->A_1)-(dz_th_3*h_12+dz_th_2*(h_11+h_01))*fields->A_2*fields->A_2-fields->A_1*((dx_th_3*h_12+dx_th_2*(h_11+h_01))*fields->A_2+(dz_th_3*h_12-dy_th_3*h_12+dz_th_2*(h_11+h_01)-dy_th_2*(h_11+h_01))*fields->A_1))/120;
                                      drdth3 -= ((dy_th_3-dz_th_3)*(G_22+G_12+G_02+(h_22+h_12+h_02)*fields->A_1*fields->A_1)+(h_22+h_12+h_02)*(fields->A_3*(dy_th_3*fields->A_3+(dy_th_3-dz_th_3)*fields->A_2+dx_th_3*fields->A_1)-fields->A_2*(dz_th_3*fields->A_2+dx_th_3*fields->A_1)))/120;
                                      break;
                         }
                         break;
                case 0:  switch(dk)
                         {
                             case -1: drdth1 -= (dx_th_3*G_02+dx_th_2*G_01+dx_th_1*G_00+fields->A_3*((dx_th_3*h_02+dx_th_2*h_01+dx_th_1*h_00)*fields->A_3-(dz_th_3*h_02+dz_th_2*h_01+dz_th_1*h_00)*fields->A_1)+fields->A_2*((dx_th_3*h_02+dx_th_2*h_01+dx_th_1*h_00)*fields->A_2-(dy_th_3*h_02+dy_th_2*h_01+dy_th_1*h_00)*fields->A_1))/30;
                                      drdth2 -= (dx_th_3*G_12+dx_th_2*(G_11+G_01)+fields->A_3*((dx_th_3*h_12+dx_th_2*(h_11+h_01))*fields->A_3-(dz_th_3*h_12+dz_th_2*(h_11+h_01))*fields->A_1)+fields->A_2*((dx_th_3*h_12+dx_th_2*(h_11+h_01))*fields->A_2-(dy_th_3*h_12+dy_th_2*(h_11+h_01))*fields->A_1))/30;
                                      drdth3 -= (dx_th_3*(G_22+G_12+G_02)+(h_22+h_12+h_02)*(fields->A_3*(dx_th_3*fields->A_3-dz_th_3*fields->A_1)+fields->A_2*(dx_th_3*fields->A_2-dy_th_3*fields->A_1)))/30;
                                      break;
                             case 0:  drdth1 -= -(dx_th_3*G_02+dx_th_2*G_01+dx_th_1*G_00+fields->A_3*((dx_th_3*h_02+dx_th_2*h_01+dx_th_1*h_00)*fields->A_3-(dz_th_3*h_02+dz_th_2*h_01+dz_th_1*h_00)*fields->A_1)+fields->A_2*((dx_th_3*h_02+dx_th_2*h_01+dx_th_1*h_00)*fields->A_2-(dy_th_3*h_02+dy_th_2*h_01+dy_th_1*h_00)*fields->A_1))/20;
                                      drdth2 -= -(dx_th_3*G_12+dx_th_2*(G_11+G_01)+fields->A_3*((dx_th_3*h_12+dx_th_2*(h_11+h_01))*fields->A_3-(dz_th_3*h_12+dz_th_2*(h_11+h_01))*fields->A_1)+fields->A_2*((dx_th_3*h_12+dx_th_2*(h_11+h_01))*fields->A_2-(dy_th_3*h_12+dy_th_2*(h_11+h_01))*fields->A_1))/20;
                                      drdth3 -= -(dx_th_3*(G_22+G_12+G_02)+(h_22+h_12+h_02)*(fields->A_3*(dx_th_3*fields->A_3-dz_th_3*fields->A_1)+fields->A_2*(dx_th_3*fields->A_2-dy_th_3*fields->A_1)))/20;
                                      break;
                             case 1:  drdth1 -= (dx_th_3*G_02+dx_th_2*G_01+dx_th_1*G_00+fields->A_3*((dx_th_3*h_02+dx_th_2*h_01+dx_th_1*h_00)*fields->A_3-(dz_th_3*h_02+dz_th_2*h_01+dz_th_1*h_00)*fields->A_1)+fields->A_2*((dx_th_3*h_02+dx_th_2*h_01+dx_th_1*h_00)*fields->A_2-(dy_th_3*h_02+dy_th_2*h_01+dy_th_1*h_00)*fields->A_1))/30;
                                      drdth2 -= (dx_th_3*G_12+dx_th_2*(G_11+G_01)+fields->A_3*((dx_th_3*h_12+dx_th_2*(h_11+h_01))*fields->A_3-(dz_th_3*h_12+dz_th_2*(h_11+h_01))*fields->A_1)+fields->A_2*((dx_th_3*h_12+dx_th_2*(h_11+h_01))*fields->A_2-(dy_th_3*h_12+dy_th_2*(h_11+h_01))*fields->A_1))/30;
                                      drdth3 -= (dx_th_3*(G_22+G_12+G_02)+(h_22+h_12+h_02)*(fields->A_3*(dx_th_3*fields->A_3-dz_th_3*fields->A_1)+fields->A_2*(dx_th_3*fields->A_2-dy_th_3*fields->A_1)))/30;
                                      break;
                         }
                         break;
                case 1:  switch(dk)
                         {
                             case -1: drdth1 -= -(-dz_th_3*G_02+dy_th_3*G_02-dz_th_2*G_01+dy_th_2*G_01-dz_th_1*G_00+dy_th_1*G_00+fields->A_3*((dy_th_3*h_02+dy_th_2*h_01+dy_th_1*h_00)*fields->A_3-(dz_th_3*h_02-dy_th_3*h_02+dz_th_2*h_01-dy_th_2*h_01+dz_th_1*h_00-dy_th_1*h_00)*fields->A_2+(dx_th_3*h_02+dx_th_2*h_01+dx_th_1*h_00)*fields->A_1)-(dz_th_3*h_02+dz_th_2*h_01+dz_th_1*h_00)*fields->A_2*fields->A_2-fields->A_1*((dx_th_3*h_02+dx_th_2*h_01+dx_th_1*h_00)*fields->A_2+(dz_th_3*h_02-dy_th_3*h_02+dz_th_2*h_01-dy_th_2*h_01+dz_th_1*h_00-dy_th_1*h_00)*fields->A_1))/120;
                                      drdth2 -= -(-dz_th_3*G_12+dy_th_3*G_12-dz_th_2*(G_11+G_01)+dy_th_2*(G_11+G_01)+fields->A_3*((dy_th_3*h_12+dy_th_2*(h_11+h_01))*fields->A_3-(dz_th_3*h_12-dy_th_3*h_12+dz_th_2*(h_11+h_01)-dy_th_2*(h_11+h_01))*fields->A_2+(dx_th_3*h_12+dx_th_2*(h_11+h_01))*fields->A_1)-(dz_th_3*h_12+dz_th_2*(h_11+h_01))*fields->A_2*fields->A_2-fields->A_1*((dx_th_3*h_12+dx_th_2*(h_11+h_01))*fields->A_2+(dz_th_3*h_12-dy_th_3*h_12+dz_th_2*(h_11+h_01)-dy_th_2*(h_11+h_01))*fields->A_1))/120;
                                      drdth3 -= -((dy_th_3-dz_th_3)*(G_22+G_12+G_02+(h_22+h_12+h_02)*fields->A_1*fields->A_1)+(h_22+h_12+h_02)*(fields->A_3*(dy_th_3*fields->A_3+(dy_th_3-dz_th_3)*fields->A_2+dx_th_3*fields->A_1)-fields->A_2*(dz_th_3*fields->A_2+dx_th_3*fields->A_1)))/120;
                                      break;
                             case 0:  drdth1 -= (dx_th_3*G_02+dx_th_2*G_01+dx_th_1*G_00+fields->A_3*((dx_th_3*h_02+dx_th_2*h_01+dx_th_1*h_00)*fields->A_3-(dz_th_3*h_02+dz_th_2*h_01+dz_th_1*h_00)*fields->A_1)+fields->A_2*((dx_th_3*h_02+dx_th_2*h_01+dx_th_1*h_00)*fields->A_2-(dy_th_3*h_02+dy_th_2*h_01+dy_th_1*h_00)*fields->A_1))/30;
                                      drdth2 -= (dx_th_3*G_12+dx_th_2*(G_11+G_01)+fields->A_3*((dx_th_3*h_12+dx_th_2*(h_11+h_01))*fields->A_3-(dz_th_3*h_12+dz_th_2*(h_11+h_01))*fields->A_1)+fields->A_2*((dx_th_3*h_12+dx_th_2*(h_11+h_01))*fields->A_2-(dy_th_3*h_12+dy_th_2*(h_11+h_01))*fields->A_1))/30;
                                      drdth3 -= (dx_th_3*(G_22+G_12+G_02)+(h_22+h_12+h_02)*(fields->A_3*(dx_th_3*fields->A_3-dz_th_3*fields->A_1)+fields->A_2*(dx_th_3*fields->A_2-dy_th_3*fields->A_1)))/30;
                                      break;
                             case 1:  drdth1 -= -(dz_th_3*G_02+dy_th_3*G_02+dz_th_2*G_01+dy_th_2*G_01+dz_th_1*G_00+dy_th_1*G_00+fields->A_3*((dy_th_3*h_02+dy_th_2*h_01+dy_th_1*h_00)*fields->A_3-(dz_th_3*h_02+dy_th_3*h_02+dz_th_2*h_01+dy_th_2*h_01+dz_th_1*h_00+dy_th_1*h_00)*fields->A_2-(dx_th_3*h_02+dx_th_2*h_01+dx_th_1*h_00)*fields->A_1)+(dz_th_3*h_02+dz_th_2*h_01+dz_th_1*h_00)*fields->A_2*fields->A_2+fields->A_1*((dz_th_3*h_02+dy_th_3*h_02+dz_th_2*h_01+dy_th_2*h_01+dz_th_1*h_00+dy_th_1*h_00)*fields->A_1-(dx_th_3*h_02+dx_th_2*h_01+dx_th_1*h_00)*fields->A_2))/120;
                                      drdth2 -= -(dz_th_3*G_12+dy_th_3*G_12+dz_th_2*(G_11+G_01)+dy_th_2*(G_11+G_01)+fields->A_3*((dy_th_3*h_12+dy_th_2*(h_11+h_01))*fields->A_3-(dz_th_3*h_12+dy_th_3*h_12+dz_th_2*(h_11+h_01)+dy_th_2*(h_11+h_01))*fields->A_2-(dx_th_3*h_12+dx_th_2*(h_11+h_01))*fields->A_1)+(dz_th_3*h_12+dz_th_2*(h_11+h_01))*fields->A_2*fields->A_2+fields->A_1*((dz_th_3*h_12+dy_th_3*h_12+dz_th_2*(h_11+h_01)+dy_th_2*(h_11+h_01))*fields->A_1-(dx_th_3*h_12+dx_th_2*(h_11+h_01))*fields->A_2))/120;
                                      drdth3 -= -((dz_th_3+dy_th_3)*(G_22+G_12+G_02+(h_22+h_12+h_02)*fields->A_1*fields->A_1)+(h_22+h_12+h_02)*(fields->A_3*(dy_th_3*fields->A_3-((dz_th_3+dy_th_3)*fields->A_2+dx_th_3*fields->A_1))+fields->A_2*(dz_th_3*fields->A_2-dx_th_3*fields->A_1)))/120;
                                      break;
                         }
                         break;
              }
              break;
     case -1: switch(dj)
              {
                case -2: switch(dk)
                         {
                             case -1: drdth1 -= (dz_th_3*G_02+dx_th_3*G_02+dz_th_2*G_01+dx_th_2*G_01+dz_th_1*G_00+dx_th_1*G_00+fields->A_3*((dx_th_3*h_02+dx_th_2*h_01+dx_th_1*h_00)*fields->A_3-(dy_th_3*h_02+dy_th_2*h_01+dy_th_1*h_00)*fields->A_2-(dz_th_3*h_02+dx_th_3*h_02+dz_th_2*h_01+dx_th_2*h_01+dz_th_1*h_00+dx_th_1*h_00)*fields->A_1)+(dz_th_3*h_02+dx_th_3*h_02+dz_th_2*h_01+dx_th_2*h_01+dz_th_1*h_00+dx_th_1*h_00)*fields->A_2*fields->A_2+fields->A_1*((dz_th_3*h_02+dz_th_2*h_01+dz_th_1*h_00)*fields->A_1-(dy_th_3*h_02+dy_th_2*h_01+dy_th_1*h_00)*fields->A_2))/120;
                                      drdth2 -= (dz_th_3*G_12+dx_th_3*G_12+dz_th_2*(G_11+G_01)+dx_th_2*(G_11+G_01)+fields->A_3*((dx_th_3*h_12+dx_th_2*(h_11+h_01))*fields->A_3-(dy_th_3*h_12+dy_th_2*(h_11+h_01))*fields->A_2-(dz_th_3*h_12+dx_th_3*h_12+dz_th_2*(h_11+h_01)+dx_th_2*(h_11+h_01))*fields->A_1)+(dz_th_3*h_12+dx_th_3*h_12+dz_th_2*(h_11+h_01)+dx_th_2*(h_11+h_01))*fields->A_2*fields->A_2+fields->A_1*((dz_th_3*h_12+dz_th_2*(h_11+h_01))*fields->A_1-(dy_th_3*h_12+dy_th_2*(h_11+h_01))*fields->A_2))/120;
                                      drdth3 -= ((dz_th_3+dx_th_3)*(G_22+G_12+G_02)+(h_22+h_12+h_02)*(fields->A_3*(dx_th_3*fields->A_3-(dy_th_3*fields->A_2+(dz_th_3+dx_th_3)*fields->A_1))+(dz_th_3+dx_th_3)*fields->A_2*fields->A_2+fields->A_1*(dz_th_3*fields->A_1-dy_th_3*fields->A_2)))/120;
                                      break;
                             case 0:  drdth1 -= (dy_th_3*G_02+dy_th_2*G_01+dy_th_1*G_00+fields->A_3*((dy_th_3*h_02+dy_th_2*h_01+dy_th_1*h_00)*fields->A_3-(dz_th_3*h_02+dz_th_2*h_01+dz_th_1*h_00)*fields->A_2)+fields->A_1*((dy_th_3*h_02+dy_th_2*h_01+dy_th_1*h_00)*fields->A_1-(dx_th_3*h_02+dx_th_2*h_01+dx_th_1*h_00)*fields->A_2))/30;
                                      drdth2 -= (dy_th_3*G_12+dy_th_2*(G_11+G_01)+fields->A_3*((dy_th_3*h_12+dy_th_2*(h_11+h_01))*fields->A_3-(dz_th_3*h_12+dz_th_2*(h_11+h_01))*fields->A_2)+fields->A_1*((dy_th_3*h_12+dy_th_2*(h_11+h_01))*fields->A_1-(dx_th_3*h_12+dx_th_2*(h_11+h_01))*fields->A_2))/30;
                                      drdth3 -= (dy_th_3*(G_22+G_12+G_02+(h_22+h_12+h_02)*fields->A_1*fields->A_1)+(h_22+h_12+h_02)*(dy_th_3*fields->A_3*fields->A_3-fields->A_2*(dz_th_3*fields->A_3+dx_th_3*fields->A_1)))/30;
                                      break;
                             case 1:  drdth1 -= (-dz_th_3*G_02+dx_th_3*G_02-dz_th_2*G_01+dx_th_2*G_01-dz_th_1*G_00+dx_th_1*G_00+fields->A_3*((dx_th_3*h_02+dx_th_2*h_01+dx_th_1*h_00)*fields->A_3+(dy_th_3*h_02+dy_th_2*h_01+dy_th_1*h_00)*fields->A_2-(dz_th_3*h_02-dx_th_3*h_02+dz_th_2*h_01-dx_th_2*h_01+dz_th_1*h_00-dx_th_1*h_00)*fields->A_1)-(dz_th_3*h_02-dx_th_3*h_02+dz_th_2*h_01-dx_th_2*h_01+dz_th_1*h_00-dx_th_1*h_00)*fields->A_2*fields->A_2-fields->A_1*((dy_th_3*h_02+dy_th_2*h_01+dy_th_1*h_00)*fields->A_2+(dz_th_3*h_02+dz_th_2*h_01+dz_th_1*h_00)*fields->A_1))/120;
                                      drdth2 -= (-dz_th_3*G_12+dx_th_3*G_12-dz_th_2*(G_11+G_01)+dx_th_2*(G_11+G_01)+fields->A_3*((dx_th_3*h_12+dx_th_2*(h_11+h_01))*fields->A_3+(dy_th_3*h_12+dy_th_2*(h_11+h_01))*fields->A_2-(dz_th_3*h_12-dx_th_3*h_12+dz_th_2*(h_11+h_01)-dx_th_2*(h_11+h_01))*fields->A_1)-(dz_th_3*h_12-dx_th_3*h_12+dz_th_2*(h_11+h_01)-dx_th_2*(h_11+h_01))*fields->A_2*fields->A_2-fields->A_1*((dy_th_3*h_12+dy_th_2*(h_11+h_01))*fields->A_2+(dz_th_3*h_12+dz_th_2*(h_11+h_01))*fields->A_1))/120;
                                      drdth3 -= ((dx_th_3-dz_th_3)*(G_22+G_12+G_02)+(h_22+h_12+h_02)*(fields->A_3*(dx_th_3*fields->A_3+dy_th_3*fields->A_2+(dx_th_3-dz_th_3)*fields->A_1)+(dx_th_3-dz_th_3)*fields->A_2*fields->A_2-fields->A_1*(dy_th_3*fields->A_2+dz_th_3*fields->A_1)))/120;
                                      break;
                         }
                         break;
                case -1: switch(dk)
                         {
                             case -2: drdth1 -= (dy_th_3*G_02+dx_th_3*G_02+dy_th_2*G_01+dx_th_2*G_01+dy_th_1*G_00+dx_th_1*G_00+fields->A_3*((dy_th_3*h_02+dx_th_3*h_02+dy_th_2*h_01+dx_th_2*h_01+dy_th_1*h_00+dx_th_1*h_00)*fields->A_3-(dz_th_3*h_02+dz_th_2*h_01+dz_th_1*h_00)*(fields->A_2+fields->A_1))+(dx_th_3*h_02+dx_th_2*h_01+dx_th_1*h_00)*fields->A_2*fields->A_2+fields->A_1*((dy_th_3*h_02+dy_th_2*h_01+dy_th_1*h_00)*fields->A_1-(dy_th_3*h_02+dx_th_3*h_02+dy_th_2*h_01+dx_th_2*h_01+dy_th_1*h_00+dx_th_1*h_00)*fields->A_2))/120;
                                      drdth2 -= (dy_th_3*G_12+dx_th_3*G_12+dy_th_2*(G_11+G_01)+dx_th_2*(G_11+G_01)+fields->A_3*((dy_th_3*h_12+dx_th_3*h_12+dy_th_2*(h_11+h_01)+dx_th_2*(h_11+h_01))*fields->A_3-(dz_th_3*h_12+dz_th_2*(h_11+h_01))*(fields->A_2+fields->A_1))+(dx_th_3*h_12+dx_th_2*(h_11+h_01))*fields->A_2*fields->A_2+fields->A_1*((dy_th_3*h_12+dy_th_2*(h_11+h_01))*fields->A_1-(dy_th_3*h_12+dx_th_3*h_12+dy_th_2*(h_11+h_01)+dx_th_2*(h_11+h_01))*fields->A_2))/120;
                                      drdth3 -= ((dy_th_3+dx_th_3)*(G_22+G_12+G_02)+(h_22+h_12+h_02)*(fields->A_3*((dy_th_3+dx_th_3)*fields->A_3-dz_th_3*(fields->A_2+fields->A_1))+dx_th_3*fields->A_2*fields->A_2+fields->A_1*(dy_th_3*fields->A_1-(dy_th_3+dx_th_3)*fields->A_2)))/120;
                                      break;
                             case -1: drdth1 -= -(dz_th_3*G_02+dy_th_3*G_02+dx_th_3*G_02+dz_th_2*G_01+dy_th_2*G_01+dx_th_2*G_01+dz_th_1*G_00+dy_th_1*G_00+dx_th_1*G_00+fields->A_3*((dy_th_3*h_02+dx_th_3*h_02+dy_th_2*h_01+dx_th_2*h_01+dy_th_1*h_00+dx_th_1*h_00)*fields->A_3-(dz_th_3*h_02+dy_th_3*h_02+dz_th_2*h_01+dy_th_2*h_01+dz_th_1*h_00+dy_th_1*h_00)*fields->A_2-(dz_th_3*h_02+dx_th_3*h_02+dz_th_2*h_01+dx_th_2*h_01+dz_th_1*h_00+dx_th_1*h_00)*fields->A_1)+(dz_th_3*h_02+dx_th_3*h_02+dz_th_2*h_01+dx_th_2*h_01+dz_th_1*h_00+dx_th_1*h_00)*fields->A_2*fields->A_2+fields->A_1*((dz_th_3*h_02+dy_th_3*h_02+dz_th_2*h_01+dy_th_2*h_01+dz_th_1*h_00+dy_th_1*h_00)*fields->A_1-(dy_th_3*h_02+dx_th_3*h_02+dy_th_2*h_01+dx_th_2*h_01+dy_th_1*h_00+dx_th_1*h_00)*fields->A_2))/30;
                                      drdth2 -= -(dz_th_3*G_12+dy_th_3*G_12+dx_th_3*G_12+dz_th_2*(G_11+G_01)+dy_th_2*(G_11+G_01)+dx_th_2*(G_11+G_01)+fields->A_3*((dy_th_3*h_12+dx_th_3*h_12+dy_th_2*(h_11+h_01)+dx_th_2*(h_11+h_01))*fields->A_3-(dz_th_3*h_12+dy_th_3*h_12+dz_th_2*(h_11+h_01)+dy_th_2*(h_11+h_01))*fields->A_2-(dz_th_3*h_12+dx_th_3*h_12+dz_th_2*(h_11+h_01)+dx_th_2*(h_11+h_01))*fields->A_1)+(dz_th_3*h_12+dx_th_3*h_12+dz_th_2*(h_11+h_01)+dx_th_2*(h_11+h_01))*fields->A_2*fields->A_2+fields->A_1*((dz_th_3*h_12+dy_th_3*h_12+dz_th_2*(h_11+h_01)+dy_th_2*(h_11+h_01))*fields->A_1-(dy_th_3*h_12+dx_th_3*h_12+dy_th_2*(h_11+h_01)+dx_th_2*(h_11+h_01))*fields->A_2))/30;
                                      drdth3 -= -((dz_th_3+dy_th_3+dx_th_3)*(G_22+G_12+G_02)+(h_22+h_12+h_02)*(fields->A_3*((dy_th_3+dx_th_3)*fields->A_3-((dz_th_3+dy_th_3)*fields->A_2+(dz_th_3+dx_th_3)*fields->A_1))+(dz_th_3+dx_th_3)*fields->A_2*fields->A_2+fields->A_1*((dz_th_3+dy_th_3)*fields->A_1-(dy_th_3+dx_th_3)*fields->A_2)))/30;
                                      break;
                             case 0:  drdth1 -= -(dy_th_3*G_02+dx_th_3*G_02+dy_th_2*G_01+dx_th_2*G_01+dy_th_1*G_00+dx_th_1*G_00+fields->A_3*((dy_th_3*h_02+dx_th_3*h_02+dy_th_2*h_01+dx_th_2*h_01+dy_th_1*h_00+dx_th_1*h_00)*fields->A_3-(dz_th_3*h_02+dz_th_2*h_01+dz_th_1*h_00)*(fields->A_2+fields->A_1))+(dx_th_3*h_02+dx_th_2*h_01+dx_th_1*h_00)*fields->A_2*fields->A_2+fields->A_1*((dy_th_3*h_02+dy_th_2*h_01+dy_th_1*h_00)*fields->A_1-(dy_th_3*h_02+dx_th_3*h_02+dy_th_2*h_01+dx_th_2*h_01+dy_th_1*h_00+dx_th_1*h_00)*fields->A_2))/12;
                                      drdth2 -= -(dy_th_3*G_12+dx_th_3*G_12+dy_th_2*(G_11+G_01)+dx_th_2*(G_11+G_01)+fields->A_3*((dy_th_3*h_12+dx_th_3*h_12+dy_th_2*(h_11+h_01)+dx_th_2*(h_11+h_01))*fields->A_3-(dz_th_3*h_12+dz_th_2*(h_11+h_01))*(fields->A_2+fields->A_1))+(dx_th_3*h_12+dx_th_2*(h_11+h_01))*fields->A_2*fields->A_2+fields->A_1*((dy_th_3*h_12+dy_th_2*(h_11+h_01))*fields->A_1-(dy_th_3*h_12+dx_th_3*h_12+dy_th_2*(h_11+h_01)+dx_th_2*(h_11+h_01))*fields->A_2))/12;
                                      drdth3 -= -((dy_th_3+dx_th_3)*(G_22+G_12+G_02)+(h_22+h_12+h_02)*(fields->A_3*((dy_th_3+dx_th_3)*fields->A_3-dz_th_3*(fields->A_2+fields->A_1))+dx_th_3*fields->A_2*fields->A_2+fields->A_1*(dy_th_3*fields->A_1-(dy_th_3+dx_th_3)*fields->A_2)))/12;
                                      break;
                             case 1:  drdth1 -= -(-dz_th_3*G_02+dy_th_3*G_02+dx_th_3*G_02-dz_th_2*G_01+dy_th_2*G_01+dx_th_2*G_01-dz_th_1*G_00+dy_th_1*G_00+dx_th_1*G_00+(dy_th_3*h_02+dx_th_3*h_02+dy_th_2*h_01+dx_th_2*h_01+dy_th_1*h_00+dx_th_1*h_00)*fields->A_3*fields->A_3-((dz_th_3*h_02-dy_th_3*h_02+dz_th_2*h_01-dy_th_2*h_01+dz_th_1*h_00-dy_th_1*h_00)*fields->A_2+(dz_th_3*h_02-dx_th_3*h_02+dz_th_2*h_01-dx_th_2*h_01+dz_th_1*h_00-dx_th_1*h_00)*fields->A_1)*fields->A_3-(dz_th_3*h_02-dx_th_3*h_02+dz_th_2*h_01-dx_th_2*h_01+dz_th_1*h_00-dx_th_1*h_00)*fields->A_2*fields->A_2-fields->A_1*((dy_th_3*h_02+dx_th_3*h_02+dy_th_2*h_01+dx_th_2*h_01+dy_th_1*h_00+dx_th_1*h_00)*fields->A_2+(dz_th_3*h_02-dy_th_3*h_02+dz_th_2*h_01-dy_th_2*h_01+dz_th_1*h_00-dy_th_1*h_00)*fields->A_1))/30;
                                      drdth2 -= -(-dz_th_3*G_12+dy_th_3*G_12+dx_th_3*G_12-dz_th_2*(G_11+G_01)+dy_th_2*(G_11+G_01)+dx_th_2*(G_11+G_01)+(dy_th_3*h_12+dx_th_3*h_12+dy_th_2*(h_11+h_01)+dx_th_2*(h_11+h_01))*fields->A_3*fields->A_3-((dz_th_3*h_12-dy_th_3*h_12+dz_th_2*(h_11+h_01)-dy_th_2*(h_11+h_01))*fields->A_2+(dz_th_3*h_12-dx_th_3*h_12+dz_th_2*(h_11+h_01)-dx_th_2*(h_11+h_01))*fields->A_1)*fields->A_3-(dz_th_3*h_12-dx_th_3*h_12+dz_th_2*(h_11+h_01)-dx_th_2*(h_11+h_01))*fields->A_2*fields->A_2-fields->A_1*((dy_th_3*h_12+dx_th_3*h_12+dy_th_2*(h_11+h_01)+dx_th_2*(h_11+h_01))*fields->A_2+(dz_th_3*h_12-dy_th_3*h_12+dz_th_2*(h_11+h_01)-dy_th_2*(h_11+h_01))*fields->A_1))/30;
                                      drdth3 -= -((-dz_th_3+dy_th_3+dx_th_3)*(G_22+G_12+G_02)+(h_22+h_12+h_02)*(fields->A_3*((dy_th_3+dx_th_3)*fields->A_3+(dy_th_3-dz_th_3)*fields->A_2+(dx_th_3-dz_th_3)*fields->A_1)+(dx_th_3-dz_th_3)*fields->A_2*fields->A_2+fields->A_1*((dy_th_3-dz_th_3)*fields->A_1-(dy_th_3+dx_th_3)*fields->A_2)))/30;
                                      break;
                             case 2:  drdth1 -= (dy_th_3*G_02+dx_th_3*G_02+dy_th_2*G_01+dx_th_2*G_01+dy_th_1*G_00+dx_th_1*G_00+fields->A_3*((dy_th_3*h_02+dx_th_3*h_02+dy_th_2*h_01+dx_th_2*h_01+dy_th_1*h_00+dx_th_1*h_00)*fields->A_3-(dz_th_3*h_02+dz_th_2*h_01+dz_th_1*h_00)*(fields->A_2+fields->A_1))+(dx_th_3*h_02+dx_th_2*h_01+dx_th_1*h_00)*fields->A_2*fields->A_2+fields->A_1*((dy_th_3*h_02+dy_th_2*h_01+dy_th_1*h_00)*fields->A_1-(dy_th_3*h_02+dx_th_3*h_02+dy_th_2*h_01+dx_th_2*h_01+dy_th_1*h_00+dx_th_1*h_00)*fields->A_2))/120;
                                      drdth2 -= (dy_th_3*G_12+dx_th_3*G_12+dy_th_2*(G_11+G_01)+dx_th_2*(G_11+G_01)+fields->A_3*((dy_th_3*h_12+dx_th_3*h_12+dy_th_2*(h_11+h_01)+dx_th_2*(h_11+h_01))*fields->A_3-(dz_th_3*h_12+dz_th_2*(h_11+h_01))*(fields->A_2+fields->A_1))+(dx_th_3*h_12+dx_th_2*(h_11+h_01))*fields->A_2*fields->A_2+fields->A_1*((dy_th_3*h_12+dy_th_2*(h_11+h_01))*fields->A_1-(dy_th_3*h_12+dx_th_3*h_12+dy_th_2*(h_11+h_01)+dx_th_2*(h_11+h_01))*fields->A_2))/120;
                                      drdth3 -= ((dy_th_3+dx_th_3)*(G_22+G_12+G_02)+(h_22+h_12+h_02)*(fields->A_3*((dy_th_3+dx_th_3)*fields->A_3-dz_th_3*(fields->A_2+fields->A_1))+dx_th_3*fields->A_2*fields->A_2+fields->A_1*(dy_th_3*fields->A_1-(dy_th_3+dx_th_3)*fields->A_2)))/120;
                                      break;
                         }
                         break;
                case 0:  switch(dk)
                         {
                             case-2 : drdth1 -= -(-dz_th_3*G_02-dz_th_2*G_01-dz_th_1*G_00+((dy_th_3*h_02+dy_th_2*h_01+dy_th_1*h_00)*fields->A_2+(dx_th_3*h_02+dx_th_2*h_01+dx_th_1*h_00)*fields->A_1)*fields->A_3-(dz_th_3*h_02+dz_th_2*h_01+dz_th_1*h_00)*(fields->A_2*fields->A_2+fields->A_1*fields->A_1))/30;
                                      drdth2 -= -(-dz_th_3*G_12-dz_th_2*(G_11+G_01)+((dy_th_3*h_12+dy_th_2*(h_11+h_01))*fields->A_2+(dx_th_3*h_12+dx_th_2*(h_11+h_01))*fields->A_1)*fields->A_3-(dz_th_3*h_12+dz_th_2*(h_11+h_01))*(fields->A_2*fields->A_2+fields->A_1*fields->A_1))/30;
                                      drdth3 -= -((h_22+h_12+h_02)*(dy_th_3*fields->A_2+dx_th_3*fields->A_1)*fields->A_3-dz_th_3*(G_22+G_12+G_02+(h_22+h_12+h_02)*(fields->A_2*fields->A_2+fields->A_1*fields->A_1)))/30;
                                      break;
                             case -1: drdth1 -= -(dz_th_3*G_02+dx_th_3*G_02+dz_th_2*G_01+dx_th_2*G_01+dz_th_1*G_00+dx_th_1*G_00+fields->A_3*((dx_th_3*h_02+dx_th_2*h_01+dx_th_1*h_00)*fields->A_3-(dy_th_3*h_02+dy_th_2*h_01+dy_th_1*h_00)*fields->A_2-(dz_th_3*h_02+dx_th_3*h_02+dz_th_2*h_01+dx_th_2*h_01+dz_th_1*h_00+dx_th_1*h_00)*fields->A_1)+(dz_th_3*h_02+dx_th_3*h_02+dz_th_2*h_01+dx_th_2*h_01+dz_th_1*h_00+dx_th_1*h_00)*fields->A_2*fields->A_2+fields->A_1*((dz_th_3*h_02+dz_th_2*h_01+dz_th_1*h_00)*fields->A_1-(dy_th_3*h_02+dy_th_2*h_01+dy_th_1*h_00)*fields->A_2))/12;
                                      drdth2 -= -(dz_th_3*G_12+dx_th_3*G_12+dz_th_2*(G_11+G_01)+dx_th_2*(G_11+G_01)+fields->A_3*((dx_th_3*h_12+dx_th_2*(h_11+h_01))*fields->A_3-(dy_th_3*h_12+dy_th_2*(h_11+h_01))*fields->A_2-(dz_th_3*h_12+dx_th_3*h_12+dz_th_2*(h_11+h_01)+dx_th_2*(h_11+h_01))*fields->A_1)+(dz_th_3*h_12+dx_th_3*h_12+dz_th_2*(h_11+h_01)+dx_th_2*(h_11+h_01))*fields->A_2*fields->A_2+fields->A_1*((dz_th_3*h_12+dz_th_2*(h_11+h_01))*fields->A_1-(dy_th_3*h_12+dy_th_2*(h_11+h_01))*fields->A_2))/12;
                                      drdth3 -= -((dz_th_3+dx_th_3)*(G_22+G_12+G_02)+(h_22+h_12+h_02)*(fields->A_3*(dx_th_3*fields->A_3-(dy_th_3*fields->A_2+(dz_th_3+dx_th_3)*fields->A_1))+(dz_th_3+dx_th_3)*fields->A_2*fields->A_2+fields->A_1*(dz_th_3*fields->A_1-dy_th_3*fields->A_2)))/12;
                                      break;
                             case 0:  drdth1 -= -4*(dx_th_3*G_02+dx_th_2*G_01+dx_th_1*G_00+fields->A_3*((dx_th_3*h_02+dx_th_2*h_01+dx_th_1*h_00)*fields->A_3-(dz_th_3*h_02+dz_th_2*h_01+dz_th_1*h_00)*fields->A_1)+fields->A_2*((dx_th_3*h_02+dx_th_2*h_01+dx_th_1*h_00)*fields->A_2-(dy_th_3*h_02+dy_th_2*h_01+dy_th_1*h_00)*fields->A_1))/15;
                                      drdth2 -= -4*(dx_th_3*G_12+dx_th_2*(G_11+G_01)+fields->A_3*((dx_th_3*h_12+dx_th_2*(h_11+h_01))*fields->A_3-(dz_th_3*h_12+dz_th_2*(h_11+h_01))*fields->A_1)+fields->A_2*((dx_th_3*h_12+dx_th_2*(h_11+h_01))*fields->A_2-(dy_th_3*h_12+dy_th_2*(h_11+h_01))*fields->A_1))/15;
                                      drdth3 -= -4*(dx_th_3*(G_22+G_12+G_02)+(h_22+h_12+h_02)*(fields->A_3*(dx_th_3*fields->A_3-dz_th_3*fields->A_1)+fields->A_2*(dx_th_3*fields->A_2-dy_th_3*fields->A_1)))/15;
                                      break;
                             case 1:  drdth1 -= -(-dz_th_3*G_02+dx_th_3*G_02-dz_th_2*G_01+dx_th_2*G_01-dz_th_1*G_00+dx_th_1*G_00+fields->A_3*((dx_th_3*h_02+dx_th_2*h_01+dx_th_1*h_00)*fields->A_3+(dy_th_3*h_02+dy_th_2*h_01+dy_th_1*h_00)*fields->A_2-(dz_th_3*h_02-dx_th_3*h_02+dz_th_2*h_01-dx_th_2*h_01+dz_th_1*h_00-dx_th_1*h_00)*fields->A_1)-(dz_th_3*h_02-dx_th_3*h_02+dz_th_2*h_01-dx_th_2*h_01+dz_th_1*h_00-dx_th_1*h_00)*fields->A_2*fields->A_2-fields->A_1*((dy_th_3*h_02+dy_th_2*h_01+dy_th_1*h_00)*fields->A_2+(dz_th_3*h_02+dz_th_2*h_01+dz_th_1*h_00)*fields->A_1))/12;
                                      drdth2 -= -(-dz_th_3*G_12+dx_th_3*G_12-dz_th_2*(G_11+G_01)+dx_th_2*(G_11+G_01)+fields->A_3*((dx_th_3*h_12+dx_th_2*(h_11+h_01))*fields->A_3+(dy_th_3*h_12+dy_th_2*(h_11+h_01))*fields->A_2-(dz_th_3*h_12-dx_th_3*h_12+dz_th_2*(h_11+h_01)-dx_th_2*(h_11+h_01))*fields->A_1)-(dz_th_3*h_12-dx_th_3*h_12+dz_th_2*(h_11+h_01)-dx_th_2*(h_11+h_01))*fields->A_2*fields->A_2-fields->A_1*((dy_th_3*h_12+dy_th_2*(h_11+h_01))*fields->A_2+(dz_th_3*h_12+dz_th_2*(h_11+h_01))*fields->A_1))/12;
                                      drdth3 -= -((dx_th_3-dz_th_3)*(G_22+G_12+G_02)+(h_22+h_12+h_02)*(fields->A_3*(dx_th_3*fields->A_3+dy_th_3*fields->A_2+(dx_th_3-dz_th_3)*fields->A_1)+(dx_th_3-dz_th_3)*fields->A_2*fields->A_2-fields->A_1*(dy_th_3*fields->A_2+dz_th_3*fields->A_1)))/12;
                                      break;
                             case 2:  drdth1 -= (-dz_th_3*G_02-dz_th_2*G_01-dz_th_1*G_00+((dy_th_3*h_02+dy_th_2*h_01+dy_th_1*h_00)*fields->A_2+(dx_th_3*h_02+dx_th_2*h_01+dx_th_1*h_00)*fields->A_1)*fields->A_3-(dz_th_3*h_02+dz_th_2*h_01+dz_th_1*h_00)*(fields->A_2*fields->A_2+fields->A_1*fields->A_1))/30;
                                      drdth2 -= (-dz_th_3*G_12-dz_th_2*(G_11+G_01)+((dy_th_3*h_12+dy_th_2*(h_11+h_01))*fields->A_2+(dx_th_3*h_12+dx_th_2*(h_11+h_01))*fields->A_1)*fields->A_3-(dz_th_3*h_12+dz_th_2*(h_11+h_01))*(fields->A_2*fields->A_2+fields->A_1*fields->A_1))/30;
                                      drdth3 -= ((h_22+h_12+h_02)*(dy_th_3*fields->A_2+dx_th_3*fields->A_1)*fields->A_3-dz_th_3*(G_22+G_12+G_02+(h_22+h_12+h_02)*(fields->A_2*fields->A_2+fields->A_1*fields->A_1)))/30;
                                      break;
                         }
                         break;
                case 1:  switch(dk)
                         {
                             case -2: drdth1 -= -(dy_th_3*G_02-dx_th_3*G_02+dy_th_2*G_01-dx_th_2*G_01+dy_th_1*G_00-dx_th_1*G_00+fields->A_3*((dy_th_3*h_02-dx_th_3*h_02+dy_th_2*h_01-dx_th_2*h_01+dy_th_1*h_00-dx_th_1*h_00)*fields->A_3+(dz_th_3*h_02+dz_th_2*h_01+dz_th_1*h_00)*(fields->A_1-fields->A_2))-(dx_th_3*h_02+dx_th_2*h_01+dx_th_1*h_00)*fields->A_2*fields->A_2+fields->A_1*((dy_th_3*h_02-dx_th_3*h_02+dy_th_2*h_01-dx_th_2*h_01+dy_th_1*h_00-dx_th_1*h_00)*fields->A_2+(dy_th_3*h_02+dy_th_2*h_01+dy_th_1*h_00)*fields->A_1))/120;
                                      drdth2 -= -(dy_th_3*G_12-dx_th_3*G_12+dy_th_2*(G_11+G_01)-dx_th_2*(G_11+G_01)+fields->A_3*((dy_th_3*h_12-dx_th_3*h_12+dy_th_2*(h_11+h_01)-dx_th_2*(h_11+h_01))*fields->A_3+(dz_th_3*h_12+dz_th_2*(h_11+h_01))*(fields->A_1-fields->A_2))-(dx_th_3*h_12+dx_th_2*(h_11+h_01))*fields->A_2*fields->A_2+fields->A_1*((dy_th_3*h_12-dx_th_3*h_12+dy_th_2*(h_11+h_01)-dx_th_2*(h_11+h_01))*fields->A_2+(dy_th_3*h_12+dy_th_2*(h_11+h_01))*fields->A_1))/120;
                                      drdth3 -= -((dy_th_3-dx_th_3)*(G_22+G_12+G_02)+(h_22+h_12+h_02)*(fields->A_3*((dy_th_3-dx_th_3)*fields->A_3+dz_th_3*(fields->A_1-fields->A_2))-dx_th_3*fields->A_2*fields->A_2+fields->A_1*((dy_th_3-dx_th_3)*fields->A_2+dy_th_3*fields->A_1)))/120;
                                      break;
                             case -1: drdth1 -= (-dz_th_3*G_02+dy_th_3*G_02-dx_th_3*G_02-dz_th_2*G_01+dy_th_2*G_01-dx_th_2*G_01-dz_th_1*G_00+dy_th_1*G_00-dx_th_1*G_00+fields->A_3*((dy_th_3*h_02-dx_th_3*h_02+dy_th_2*h_01-dx_th_2*h_01+dy_th_1*h_00-dx_th_1*h_00)*fields->A_3-(dz_th_3*h_02-dy_th_3*h_02+dz_th_2*h_01-dy_th_2*h_01+dz_th_1*h_00-dy_th_1*h_00)*fields->A_2+(dz_th_3*h_02+dx_th_3*h_02+dz_th_2*h_01+dx_th_2*h_01+dz_th_1*h_00+dx_th_1*h_00)*fields->A_1)+fields->A_2*((dy_th_3*h_02-dx_th_3*h_02+dy_th_2*h_01-dx_th_2*h_01+dy_th_1*h_00-dx_th_1*h_00)*fields->A_1-(dz_th_3*h_02+dx_th_3*h_02+dz_th_2*h_01+dx_th_2*h_01+dz_th_1*h_00+dx_th_1*h_00)*fields->A_2)-(dz_th_3*h_02-dy_th_3*h_02+dz_th_2*h_01-dy_th_2*h_01+dz_th_1*h_00-dy_th_1*h_00)*fields->A_1*fields->A_1)/30;
                                      drdth2 -= (-dz_th_3*G_12+dy_th_3*G_12-dx_th_3*G_12-dz_th_2*(G_11+G_01)+dy_th_2*(G_11+G_01)-dx_th_2*(G_11+G_01)+fields->A_3*((dy_th_3*h_12-dx_th_3*h_12+dy_th_2*(h_11+h_01)-dx_th_2*(h_11+h_01))*fields->A_3-(dz_th_3*h_12-dy_th_3*h_12+dz_th_2*(h_11+h_01)-dy_th_2*(h_11+h_01))*fields->A_2+(dz_th_3*h_12+dx_th_3*h_12+dz_th_2*(h_11+h_01)+dx_th_2*(h_11+h_01))*fields->A_1)+fields->A_2*((dy_th_3*h_12-dx_th_3*h_12+dy_th_2*(h_11+h_01)-dx_th_2*(h_11+h_01))*fields->A_1-(dz_th_3*h_12+dx_th_3*h_12+dz_th_2*(h_11+h_01)+dx_th_2*(h_11+h_01))*fields->A_2)-(dz_th_3*h_12-dy_th_3*h_12+dz_th_2*(h_11+h_01)-dy_th_2*(h_11+h_01))*fields->A_1*fields->A_1)/30;
                                      drdth3 -= ((-dz_th_3+dy_th_3-dx_th_3)*(G_22+G_12+G_02)+(h_22+h_12+h_02)*(fields->A_3*((dy_th_3-dx_th_3)*fields->A_3+(dy_th_3-dz_th_3)*fields->A_2+(dz_th_3+dx_th_3)*fields->A_1)-(dz_th_3+dx_th_3)*fields->A_2*fields->A_2+fields->A_1*((dy_th_3-dx_th_3)*fields->A_2+(dy_th_3-dz_th_3)*fields->A_1)))/30;
                                      break;
                             case 0:  drdth1 -= (dy_th_3*G_02-dx_th_3*G_02+dy_th_2*G_01-dx_th_2*G_01+dy_th_1*G_00-dx_th_1*G_00+fields->A_3*((dy_th_3*h_02-dx_th_3*h_02+dy_th_2*h_01-dx_th_2*h_01+dy_th_1*h_00-dx_th_1*h_00)*fields->A_3+(dz_th_3*h_02+dz_th_2*h_01+dz_th_1*h_00)*(fields->A_1-fields->A_2))-(dx_th_3*h_02+dx_th_2*h_01+dx_th_1*h_00)*fields->A_2*fields->A_2+fields->A_1*((dy_th_3*h_02-dx_th_3*h_02+dy_th_2*h_01-dx_th_2*h_01+dy_th_1*h_00-dx_th_1*h_00)*fields->A_2+(dy_th_3*h_02+dy_th_2*h_01+dy_th_1*h_00)*fields->A_1))/12;
                                      drdth2 -= (dy_th_3*G_12-dx_th_3*G_12+dy_th_2*(G_11+G_01)-dx_th_2*(G_11+G_01)+fields->A_3*((dy_th_3*h_12-dx_th_3*h_12+dy_th_2*(h_11+h_01)-dx_th_2*(h_11+h_01))*fields->A_3+(dz_th_3*h_12+dz_th_2*(h_11+h_01))*(fields->A_1-fields->A_2))-(dx_th_3*h_12+dx_th_2*(h_11+h_01))*fields->A_2*fields->A_2+fields->A_1*((dy_th_3*h_12-dx_th_3*h_12+dy_th_2*(h_11+h_01)-dx_th_2*(h_11+h_01))*fields->A_2+(dy_th_3*h_12+dy_th_2*(h_11+h_01))*fields->A_1))/12;
                                      drdth3 -= ((dy_th_3-dx_th_3)*(G_22+G_12+G_02)+(h_22+h_12+h_02)*(fields->A_3*((dy_th_3-dx_th_3)*fields->A_3+dz_th_3*(fields->A_1-fields->A_2))-dx_th_3*fields->A_2*fields->A_2+fields->A_1*((dy_th_3-dx_th_3)*fields->A_2+dy_th_3*fields->A_1)))/12;
                                      break;
                             case 1:  drdth1 -= (dz_th_3*G_02+dy_th_3*G_02-dx_th_3*G_02+dz_th_2*G_01+dy_th_2*G_01-dx_th_2*G_01+dz_th_1*G_00+dy_th_1*G_00-dx_th_1*G_00+fields->A_3*((dy_th_3*h_02-dx_th_3*h_02+dy_th_2*h_01-dx_th_2*h_01+dy_th_1*h_00-dx_th_1*h_00)*fields->A_3-(dz_th_3*h_02+dy_th_3*h_02+dz_th_2*h_01+dy_th_2*h_01+dz_th_1*h_00+dy_th_1*h_00)*fields->A_2+(dz_th_3*h_02-dx_th_3*h_02+dz_th_2*h_01-dx_th_2*h_01+dz_th_1*h_00-dx_th_1*h_00)*fields->A_1)+(dz_th_3*h_02-dx_th_3*h_02+dz_th_2*h_01-dx_th_2*h_01+dz_th_1*h_00-dx_th_1*h_00)*fields->A_2*fields->A_2+fields->A_1*((dy_th_3*h_02-dx_th_3*h_02+dy_th_2*h_01-dx_th_2*h_01+dy_th_1*h_00-dx_th_1*h_00)*fields->A_2+(dz_th_3*h_02+dy_th_3*h_02+dz_th_2*h_01+dy_th_2*h_01+dz_th_1*h_00+dy_th_1*h_00)*fields->A_1))/30;
                                      drdth2 -= (dz_th_3*G_12+dy_th_3*G_12-dx_th_3*G_12+dz_th_2*(G_11+G_01)+dy_th_2*(G_11+G_01)-dx_th_2*(G_11+G_01)+fields->A_3*((dy_th_3*h_12-dx_th_3*h_12+dy_th_2*(h_11+h_01)-dx_th_2*(h_11+h_01))*fields->A_3-(dz_th_3*h_12+dy_th_3*h_12+dz_th_2*(h_11+h_01)+dy_th_2*(h_11+h_01))*fields->A_2+(dz_th_3*h_12-dx_th_3*h_12+dz_th_2*(h_11+h_01)-dx_th_2*(h_11+h_01))*fields->A_1)+(dz_th_3*h_12-dx_th_3*h_12+dz_th_2*(h_11+h_01)-dx_th_2*(h_11+h_01))*fields->A_2*fields->A_2+fields->A_1*((dy_th_3*h_12-dx_th_3*h_12+dy_th_2*(h_11+h_01)-dx_th_2*(h_11+h_01))*fields->A_2+(dz_th_3*h_12+dy_th_3*h_12+dz_th_2*(h_11+h_01)+dy_th_2*(h_11+h_01))*fields->A_1))/30;
                                      drdth3 -= ((dz_th_3+dy_th_3-dx_th_3)*(G_22+G_12+G_02)+(h_22+h_12+h_02)*(fields->A_3*((dy_th_3-dx_th_3)*fields->A_3-(dz_th_3+dy_th_3)*fields->A_2+(dz_th_3-dx_th_3)*fields->A_1)+(dz_th_3-dx_th_3)*fields->A_2*fields->A_2+fields->A_1*((dy_th_3-dx_th_3)*fields->A_2+(dz_th_3+dy_th_3)*fields->A_1)))/30;
                                      break;
                             case 2:  drdth1 -= -(dy_th_3*G_02-dx_th_3*G_02+dy_th_2*G_01-dx_th_2*G_01+dy_th_1*G_00-dx_th_1*G_00+fields->A_3*((dy_th_3*h_02-dx_th_3*h_02+dy_th_2*h_01-dx_th_2*h_01+dy_th_1*h_00-dx_th_1*h_00)*fields->A_3+(dz_th_3*h_02+dz_th_2*h_01+dz_th_1*h_00)*(fields->A_1-fields->A_2))-(dx_th_3*h_02+dx_th_2*h_01+dx_th_1*h_00)*fields->A_2*fields->A_2+fields->A_1*((dy_th_3*h_02-dx_th_3*h_02+dy_th_2*h_01-dx_th_2*h_01+dy_th_1*h_00-dx_th_1*h_00)*fields->A_2+(dy_th_3*h_02+dy_th_2*h_01+dy_th_1*h_00)*fields->A_1))/120;
                                      drdth2 -= -(dy_th_3*G_12-dx_th_3*G_12+dy_th_2*(G_11+G_01)-dx_th_2*(G_11+G_01)+fields->A_3*((dy_th_3*h_12-dx_th_3*h_12+dy_th_2*(h_11+h_01)-dx_th_2*(h_11+h_01))*fields->A_3+(dz_th_3*h_12+dz_th_2*(h_11+h_01))*(fields->A_1-fields->A_2))-(dx_th_3*h_12+dx_th_2*(h_11+h_01))*fields->A_2*fields->A_2+fields->A_1*((dy_th_3*h_12-dx_th_3*h_12+dy_th_2*(h_11+h_01)-dx_th_2*(h_11+h_01))*fields->A_2+(dy_th_3*h_12+dy_th_2*(h_11+h_01))*fields->A_1))/120;
                                      drdth3 -= -((dy_th_3-dx_th_3)*(G_22+G_12+G_02)+(h_22+h_12+h_02)*(fields->A_3*((dy_th_3-dx_th_3)*fields->A_3+dz_th_3*(fields->A_1-fields->A_2))-dx_th_3*fields->A_2*fields->A_2+fields->A_1*((dy_th_3-dx_th_3)*fields->A_2+dy_th_3*fields->A_1)))/120;
                                      break;
                         }
                         break;
                case 2:  switch(dk)
                         {
                             case -1: drdth1 -= (dz_th_3*G_02+dx_th_3*G_02+dz_th_2*G_01+dx_th_2*G_01+dz_th_1*G_00+dx_th_1*G_00+fields->A_3*((dx_th_3*h_02+dx_th_2*h_01+dx_th_1*h_00)*fields->A_3-(dy_th_3*h_02+dy_th_2*h_01+dy_th_1*h_00)*fields->A_2-(dz_th_3*h_02+dx_th_3*h_02+dz_th_2*h_01+dx_th_2*h_01+dz_th_1*h_00+dx_th_1*h_00)*fields->A_1)+(dz_th_3*h_02+dx_th_3*h_02+dz_th_2*h_01+dx_th_2*h_01+dz_th_1*h_00+dx_th_1*h_00)*fields->A_2*fields->A_2+fields->A_1*((dz_th_3*h_02+dz_th_2*h_01+dz_th_1*h_00)*fields->A_1-(dy_th_3*h_02+dy_th_2*h_01+dy_th_1*h_00)*fields->A_2))/120;
                                      drdth2 -= (dz_th_3*G_12+dx_th_3*G_12+dz_th_2*(G_11+G_01)+dx_th_2*(G_11+G_01)+fields->A_3*((dx_th_3*h_12+dx_th_2*(h_11+h_01))*fields->A_3-(dy_th_3*h_12+dy_th_2*(h_11+h_01))*fields->A_2-(dz_th_3*h_12+dx_th_3*h_12+dz_th_2*(h_11+h_01)+dx_th_2*(h_11+h_01))*fields->A_1)+(dz_th_3*h_12+dx_th_3*h_12+dz_th_2*(h_11+h_01)+dx_th_2*(h_11+h_01))*fields->A_2*fields->A_2+fields->A_1*((dz_th_3*h_12+dz_th_2*(h_11+h_01))*fields->A_1-(dy_th_3*h_12+dy_th_2*(h_11+h_01))*fields->A_2))/120;
                                      drdth3 -= ((dz_th_3+dx_th_3)*(G_22+G_12+G_02)+(h_22+h_12+h_02)*(fields->A_3*(dx_th_3*fields->A_3-(dy_th_3*fields->A_2+(dz_th_3+dx_th_3)*fields->A_1))+(dz_th_3+dx_th_3)*fields->A_2*fields->A_2+fields->A_1*(dz_th_3*fields->A_1-dy_th_3*fields->A_2)))/120;
                                      break;
                             case 0:  drdth1 -= -(dy_th_3*G_02+dy_th_2*G_01+dy_th_1*G_00+fields->A_3*((dy_th_3*h_02+dy_th_2*h_01+dy_th_1*h_00)*fields->A_3-(dz_th_3*h_02+dz_th_2*h_01+dz_th_1*h_00)*fields->A_2)+fields->A_1*((dy_th_3*h_02+dy_th_2*h_01+dy_th_1*h_00)*fields->A_1-(dx_th_3*h_02+dx_th_2*h_01+dx_th_1*h_00)*fields->A_2))/30;
                                      drdth2 -= -(dy_th_3*G_12+dy_th_2*(G_11+G_01)+fields->A_3*((dy_th_3*h_12+dy_th_2*(h_11+h_01))*fields->A_3-(dz_th_3*h_12+dz_th_2*(h_11+h_01))*fields->A_2)+fields->A_1*((dy_th_3*h_12+dy_th_2*(h_11+h_01))*fields->A_1-(dx_th_3*h_12+dx_th_2*(h_11+h_01))*fields->A_2))/30;
                                      drdth3 -= -(dy_th_3*(G_22+G_12+G_02+(h_22+h_12+h_02)*fields->A_1*fields->A_1)+(h_22+h_12+h_02)*(dy_th_3*fields->A_3*fields->A_3-fields->A_2*(dz_th_3*fields->A_3+dx_th_3*fields->A_1)))/30;
                                      break;
                             case 1:  drdth1 -= (-dz_th_3*G_02+dx_th_3*G_02-dz_th_2*G_01+dx_th_2*G_01-dz_th_1*G_00+dx_th_1*G_00+fields->A_3*((dx_th_3*h_02+dx_th_2*h_01+dx_th_1*h_00)*fields->A_3+(dy_th_3*h_02+dy_th_2*h_01+dy_th_1*h_00)*fields->A_2-(dz_th_3*h_02-dx_th_3*h_02+dz_th_2*h_01-dx_th_2*h_01+dz_th_1*h_00-dx_th_1*h_00)*fields->A_1)-(dz_th_3*h_02-dx_th_3*h_02+dz_th_2*h_01-dx_th_2*h_01+dz_th_1*h_00-dx_th_1*h_00)*fields->A_2*fields->A_2-fields->A_1*((dy_th_3*h_02+dy_th_2*h_01+dy_th_1*h_00)*fields->A_2+(dz_th_3*h_02+dz_th_2*h_01+dz_th_1*h_00)*fields->A_1))/120;
                                      drdth2 -= (-dz_th_3*G_12+dx_th_3*G_12-dz_th_2*(G_11+G_01)+dx_th_2*(G_11+G_01)+fields->A_3*((dx_th_3*h_12+dx_th_2*(h_11+h_01))*fields->A_3+(dy_th_3*h_12+dy_th_2*(h_11+h_01))*fields->A_2-(dz_th_3*h_12-dx_th_3*h_12+dz_th_2*(h_11+h_01)-dx_th_2*(h_11+h_01))*fields->A_1)-(dz_th_3*h_12-dx_th_3*h_12+dz_th_2*(h_11+h_01)-dx_th_2*(h_11+h_01))*fields->A_2*fields->A_2-fields->A_1*((dy_th_3*h_12+dy_th_2*(h_11+h_01))*fields->A_2+(dz_th_3*h_12+dz_th_2*(h_11+h_01))*fields->A_1))/120;
                                      drdth3 -= ((dx_th_3-dz_th_3)*(G_22+G_12+G_02)+(h_22+h_12+h_02)*(fields->A_3*(dx_th_3*fields->A_3+dy_th_3*fields->A_2+(dx_th_3-dz_th_3)*fields->A_1)+(dx_th_3-dz_th_3)*fields->A_2*fields->A_2-fields->A_1*(dy_th_3*fields->A_2+dz_th_3*fields->A_1)))/120;
                                      break;
                         }
                         break;
              }
              break;
     case 0:  switch(dj)
              {
                case -2: switch(dk)
                         {
                             case -1: drdth1 -= (dy_th_3*G_02+dy_th_2*G_01+dy_th_1*G_00+fields->A_3*((dy_th_3*h_02+dy_th_2*h_01+dy_th_1*h_00)*fields->A_3-(dz_th_3*h_02+dz_th_2*h_01+dz_th_1*h_00)*fields->A_2)+fields->A_1*((dy_th_3*h_02+dy_th_2*h_01+dy_th_1*h_00)*fields->A_1-(dx_th_3*h_02+dx_th_2*h_01+dx_th_1*h_00)*fields->A_2))/30;
                                      drdth2 -= (dy_th_3*G_12+dy_th_2*(G_11+G_01)+fields->A_3*((dy_th_3*h_12+dy_th_2*(h_11+h_01))*fields->A_3-(dz_th_3*h_12+dz_th_2*(h_11+h_01))*fields->A_2)+fields->A_1*((dy_th_3*h_12+dy_th_2*(h_11+h_01))*fields->A_1-(dx_th_3*h_12+dx_th_2*(h_11+h_01))*fields->A_2))/30;
                                      drdth3 -= (dy_th_3*(G_22+G_12+G_02+(h_22+h_12+h_02)*fields->A_1*fields->A_1)+(h_22+h_12+h_02)*(dy_th_3*fields->A_3*fields->A_3-fields->A_2*(dz_th_3*fields->A_3+dx_th_3*fields->A_1)))/30;
                                      break;
                             case 0:  drdth1 -= -(dy_th_3*G_02+dy_th_2*G_01+dy_th_1*G_00+fields->A_3*((dy_th_3*h_02+dy_th_2*h_01+dy_th_1*h_00)*fields->A_3-(dz_th_3*h_02+dz_th_2*h_01+dz_th_1*h_00)*fields->A_2)+fields->A_1*((dy_th_3*h_02+dy_th_2*h_01+dy_th_1*h_00)*fields->A_1-(dx_th_3*h_02+dx_th_2*h_01+dx_th_1*h_00)*fields->A_2))/20;
                                      drdth2 -= -(dy_th_3*G_12+dy_th_2*(G_11+G_01)+fields->A_3*((dy_th_3*h_12+dy_th_2*(h_11+h_01))*fields->A_3-(dz_th_3*h_12+dz_th_2*(h_11+h_01))*fields->A_2)+fields->A_1*((dy_th_3*h_12+dy_th_2*(h_11+h_01))*fields->A_1-(dx_th_3*h_12+dx_th_2*(h_11+h_01))*fields->A_2))/20;
                                      drdth3 -= -(dy_th_3*(G_22+G_12+G_02+(h_22+h_12+h_02)*fields->A_1*fields->A_1)+(h_22+h_12+h_02)*(dy_th_3*fields->A_3*fields->A_3-fields->A_2*(dz_th_3*fields->A_3+dx_th_3*fields->A_1)))/20;
                                      break;
                             case 1:  drdth1 -= (dy_th_3*G_02+dy_th_2*G_01+dy_th_1*G_00+fields->A_3*((dy_th_3*h_02+dy_th_2*h_01+dy_th_1*h_00)*fields->A_3-(dz_th_3*h_02+dz_th_2*h_01+dz_th_1*h_00)*fields->A_2)+fields->A_1*((dy_th_3*h_02+dy_th_2*h_01+dy_th_1*h_00)*fields->A_1-(dx_th_3*h_02+dx_th_2*h_01+dx_th_1*h_00)*fields->A_2))/30;
                                      drdth2 -= (dy_th_3*G_12+dy_th_2*(G_11+G_01)+fields->A_3*((dy_th_3*h_12+dy_th_2*(h_11+h_01))*fields->A_3-(dz_th_3*h_12+dz_th_2*(h_11+h_01))*fields->A_2)+fields->A_1*((dy_th_3*h_12+dy_th_2*(h_11+h_01))*fields->A_1-(dx_th_3*h_12+dx_th_2*(h_11+h_01))*fields->A_2))/30;
                                      drdth3 -= (dy_th_3*(G_22+G_12+G_02+(h_22+h_12+h_02)*fields->A_1*fields->A_1)+(h_22+h_12+h_02)*(dy_th_3*fields->A_3*fields->A_3-fields->A_2*(dz_th_3*fields->A_3+dx_th_3*fields->A_1)))/30;
                                      break;
                         }
                         break;
                case -1: switch(dk)
                         {
                             case -2: drdth1 -= -(-dz_th_3*G_02-dz_th_2*G_01-dz_th_1*G_00+((dy_th_3*h_02+dy_th_2*h_01+dy_th_1*h_00)*fields->A_2+(dx_th_3*h_02+dx_th_2*h_01+dx_th_1*h_00)*fields->A_1)*fields->A_3-(dz_th_3*h_02+dz_th_2*h_01+dz_th_1*h_00)*(fields->A_2*fields->A_2+fields->A_1*fields->A_1))/30;
                                      drdth2 -= -(-dz_th_3*G_12-dz_th_2*(G_11+G_01)+((dy_th_3*h_12+dy_th_2*(h_11+h_01))*fields->A_2+(dx_th_3*h_12+dx_th_2*(h_11+h_01))*fields->A_1)*fields->A_3-(dz_th_3*h_12+dz_th_2*(h_11+h_01))*(fields->A_2*fields->A_2+fields->A_1*fields->A_1))/30;
                                      drdth3 -= -((h_22+h_12+h_02)*(dy_th_3*fields->A_2+dx_th_3*fields->A_1)*fields->A_3-dz_th_3*(G_22+G_12+G_02+(h_22+h_12+h_02)*(fields->A_2*fields->A_2+fields->A_1*fields->A_1)))/30;
                                      break;
                             case -1: drdth1 -= -(dz_th_3*G_02+dy_th_3*G_02+dz_th_2*G_01+dy_th_2*G_01+dz_th_1*G_00+dy_th_1*G_00+fields->A_3*((dy_th_3*h_02+dy_th_2*h_01+dy_th_1*h_00)*fields->A_3-(dz_th_3*h_02+dy_th_3*h_02+dz_th_2*h_01+dy_th_2*h_01+dz_th_1*h_00+dy_th_1*h_00)*fields->A_2-(dx_th_3*h_02+dx_th_2*h_01+dx_th_1*h_00)*fields->A_1)+(dz_th_3*h_02+dz_th_2*h_01+dz_th_1*h_00)*fields->A_2*fields->A_2+fields->A_1*((dz_th_3*h_02+dy_th_3*h_02+dz_th_2*h_01+dy_th_2*h_01+dz_th_1*h_00+dy_th_1*h_00)*fields->A_1-(dx_th_3*h_02+dx_th_2*h_01+dx_th_1*h_00)*fields->A_2))/12;
                                      drdth2 -= -(dz_th_3*G_12+dy_th_3*G_12+dz_th_2*(G_11+G_01)+dy_th_2*(G_11+G_01)+fields->A_3*((dy_th_3*h_12+dy_th_2*(h_11+h_01))*fields->A_3-(dz_th_3*h_12+dy_th_3*h_12+dz_th_2*(h_11+h_01)+dy_th_2*(h_11+h_01))*fields->A_2-(dx_th_3*h_12+dx_th_2*(h_11+h_01))*fields->A_1)+(dz_th_3*h_12+dz_th_2*(h_11+h_01))*fields->A_2*fields->A_2+fields->A_1*((dz_th_3*h_12+dy_th_3*h_12+dz_th_2*(h_11+h_01)+dy_th_2*(h_11+h_01))*fields->A_1-(dx_th_3*h_12+dx_th_2*(h_11+h_01))*fields->A_2))/12;
                                      drdth3 -= -((dz_th_3+dy_th_3)*(G_22+G_12+G_02+(h_22+h_12+h_02)*fields->A_1*fields->A_1)+(h_22+h_12+h_02)*(fields->A_3*(dy_th_3*fields->A_3-((dz_th_3+dy_th_3)*fields->A_2+dx_th_3*fields->A_1))+fields->A_2*(dz_th_3*fields->A_2-dx_th_3*fields->A_1)))/12;
                                      break;
                             case 0:  drdth1 -= -4*(dy_th_3*G_02+dy_th_2*G_01+dy_th_1*G_00+fields->A_3*((dy_th_3*h_02+dy_th_2*h_01+dy_th_1*h_00)*fields->A_3-(dz_th_3*h_02+dz_th_2*h_01+dz_th_1*h_00)*fields->A_2)+fields->A_1*((dy_th_3*h_02+dy_th_2*h_01+dy_th_1*h_00)*fields->A_1-(dx_th_3*h_02+dx_th_2*h_01+dx_th_1*h_00)*fields->A_2))/15;
                                      drdth2 -= -4*(dy_th_3*G_12+dy_th_2*(G_11+G_01)+fields->A_3*((dy_th_3*h_12+dy_th_2*(h_11+h_01))*fields->A_3-(dz_th_3*h_12+dz_th_2*(h_11+h_01))*fields->A_2)+fields->A_1*((dy_th_3*h_12+dy_th_2*(h_11+h_01))*fields->A_1-(dx_th_3*h_12+dx_th_2*(h_11+h_01))*fields->A_2))/15;
                                      drdth3 -= -4*(dy_th_3*(G_22+G_12+G_02+(h_22+h_12+h_02)*fields->A_1*fields->A_1)+(h_22+h_12+h_02)*(dy_th_3*fields->A_3*fields->A_3-fields->A_2*(dz_th_3*fields->A_3+dx_th_3*fields->A_1)))/15;
                                      break;
                             case 1:  drdth1 -= -(-dz_th_3*G_02+dy_th_3*G_02-dz_th_2*G_01+dy_th_2*G_01-dz_th_1*G_00+dy_th_1*G_00+fields->A_3*((dy_th_3*h_02+dy_th_2*h_01+dy_th_1*h_00)*fields->A_3-(dz_th_3*h_02-dy_th_3*h_02+dz_th_2*h_01-dy_th_2*h_01+dz_th_1*h_00-dy_th_1*h_00)*fields->A_2+(dx_th_3*h_02+dx_th_2*h_01+dx_th_1*h_00)*fields->A_1)-(dz_th_3*h_02+dz_th_2*h_01+dz_th_1*h_00)*fields->A_2*fields->A_2-fields->A_1*((dx_th_3*h_02+dx_th_2*h_01+dx_th_1*h_00)*fields->A_2+(dz_th_3*h_02-dy_th_3*h_02+dz_th_2*h_01-dy_th_2*h_01+dz_th_1*h_00-dy_th_1*h_00)*fields->A_1))/12;
                                      drdth2 -= -(-dz_th_3*G_12+dy_th_3*G_12-dz_th_2*(G_11+G_01)+dy_th_2*(G_11+G_01)+fields->A_3*((dy_th_3*h_12+dy_th_2*(h_11+h_01))*fields->A_3-(dz_th_3*h_12-dy_th_3*h_12+dz_th_2*(h_11+h_01)-dy_th_2*(h_11+h_01))*fields->A_2+(dx_th_3*h_12+dx_th_2*(h_11+h_01))*fields->A_1)-(dz_th_3*h_12+dz_th_2*(h_11+h_01))*fields->A_2*fields->A_2-fields->A_1*((dx_th_3*h_12+dx_th_2*(h_11+h_01))*fields->A_2+(dz_th_3*h_12-dy_th_3*h_12+dz_th_2*(h_11+h_01)-dy_th_2*(h_11+h_01))*fields->A_1))/12;
                                      drdth3 -= -((dy_th_3-dz_th_3)*(G_22+G_12+G_02+(h_22+h_12+h_02)*fields->A_1*fields->A_1)+(h_22+h_12+h_02)*(fields->A_3*(dy_th_3*fields->A_3+(dy_th_3-dz_th_3)*fields->A_2+dx_th_3*fields->A_1)-fields->A_2*(dz_th_3*fields->A_2+dx_th_3*fields->A_1)))/12;
                                      break;
                             case 2:  drdth1 -= (-dz_th_3*G_02-dz_th_2*G_01-dz_th_1*G_00+((dy_th_3*h_02+dy_th_2*h_01+dy_th_1*h_00)*fields->A_2+(dx_th_3*h_02+dx_th_2*h_01+dx_th_1*h_00)*fields->A_1)*fields->A_3-(dz_th_3*h_02+dz_th_2*h_01+dz_th_1*h_00)*(fields->A_2*fields->A_2+fields->A_1*fields->A_1))/30;
                                      drdth2 -= (-dz_th_3*G_12-dz_th_2*(G_11+G_01)+((dy_th_3*h_12+dy_th_2*(h_11+h_01))*fields->A_2+(dx_th_3*h_12+dx_th_2*(h_11+h_01))*fields->A_1)*fields->A_3-(dz_th_3*h_12+dz_th_2*(h_11+h_01))*(fields->A_2*fields->A_2+fields->A_1*fields->A_1))/30;
                                      drdth3 -= ((h_22+h_12+h_02)*(dy_th_3*fields->A_2+dx_th_3*fields->A_1)*fields->A_3-dz_th_3*(G_22+G_12+G_02+(h_22+h_12+h_02)*(fields->A_2*fields->A_2+fields->A_1*fields->A_1)))/30;
                                      break;
                         }
                         break;
                case 0:  switch(dk)
                         {
                             case -2: drdth1 -= (-dz_th_3*G_02-dz_th_2*G_01-dz_th_1*G_00+((dy_th_3*h_02+dy_th_2*h_01+dy_th_1*h_00)*fields->A_2+(dx_th_3*h_02+dx_th_2*h_01+dx_th_1*h_00)*fields->A_1)*fields->A_3-(dz_th_3*h_02+dz_th_2*h_01+dz_th_1*h_00)*(fields->A_2*fields->A_2+fields->A_1*fields->A_1))/20;
                                      drdth2 -= (-dz_th_3*G_12-dz_th_2*(G_11+G_01)+((dy_th_3*h_12+dy_th_2*(h_11+h_01))*fields->A_2+(dx_th_3*h_12+dx_th_2*(h_11+h_01))*fields->A_1)*fields->A_3-(dz_th_3*h_12+dz_th_2*(h_11+h_01))*(fields->A_2*fields->A_2+fields->A_1*fields->A_1))/20;
                                      drdth3 -= ((h_22+h_12+h_02)*(dy_th_3*fields->A_2+dx_th_3*fields->A_1)*fields->A_3-dz_th_3*(G_22+G_12+G_02+(h_22+h_12+h_02)*(fields->A_2*fields->A_2+fields->A_1*fields->A_1)))/20;
                                      break;
                             case -1: drdth1 -= 4*(-dz_th_3*G_02-dz_th_2*G_01-dz_th_1*G_00+((dy_th_3*h_02+dy_th_2*h_01+dy_th_1*h_00)*fields->A_2+(dx_th_3*h_02+dx_th_2*h_01+dx_th_1*h_00)*fields->A_1)*fields->A_3-(dz_th_3*h_02+dz_th_2*h_01+dz_th_1*h_00)*(fields->A_2*fields->A_2+fields->A_1*fields->A_1))/15;
                                      drdth2 -= 4*(-dz_th_3*G_12-dz_th_2*(G_11+G_01)+((dy_th_3*h_12+dy_th_2*(h_11+h_01))*fields->A_2+(dx_th_3*h_12+dx_th_2*(h_11+h_01))*fields->A_1)*fields->A_3-(dz_th_3*h_12+dz_th_2*(h_11+h_01))*(fields->A_2*fields->A_2+fields->A_1*fields->A_1))/15;
                                      drdth3 -= 4*((h_22+h_12+h_02)*(dy_th_3*fields->A_2+dx_th_3*fields->A_1)*fields->A_3-dz_th_3*(G_22+G_12+G_02+(h_22+h_12+h_02)*(fields->A_2*fields->A_2+fields->A_1*fields->A_1)))/15;
                                      break;
                             case 1:  drdth1 -= -4*(-dz_th_3*G_02-dz_th_2*G_01-dz_th_1*G_00+((dy_th_3*h_02+dy_th_2*h_01+dy_th_1*h_00)*fields->A_2+(dx_th_3*h_02+dx_th_2*h_01+dx_th_1*h_00)*fields->A_1)*fields->A_3-(dz_th_3*h_02+dz_th_2*h_01+dz_th_1*h_00)*(fields->A_2*fields->A_2+fields->A_1*fields->A_1))/15;
                                      drdth2 -= -4*(-dz_th_3*G_12-dz_th_2*(G_11+G_01)+((dy_th_3*h_12+dy_th_2*(h_11+h_01))*fields->A_2+(dx_th_3*h_12+dx_th_2*(h_11+h_01))*fields->A_1)*fields->A_3-(dz_th_3*h_12+dz_th_2*(h_11+h_01))*(fields->A_2*fields->A_2+fields->A_1*fields->A_1))/15;
                                      drdth3 -= -4*((h_22+h_12+h_02)*(dy_th_3*fields->A_2+dx_th_3*fields->A_1)*fields->A_3-dz_th_3*(G_22+G_12+G_02+(h_22+h_12+h_02)*(fields->A_2*fields->A_2+fields->A_1*fields->A_1)))/15;
                                      break;
                             case 2:  drdth1 -= -(-dz_th_3*G_02-dz_th_2*G_01-dz_th_1*G_00+((dy_th_3*h_02+dy_th_2*h_01+dy_th_1*h_00)*fields->A_2+(dx_th_3*h_02+dx_th_2*h_01+dx_th_1*h_00)*fields->A_1)*fields->A_3-(dz_th_3*h_02+dz_th_2*h_01+dz_th_1*h_00)*(fields->A_2*fields->A_2+fields->A_1*fields->A_1))/20;
                                      drdth2 -= -(-dz_th_3*G_12-dz_th_2*(G_11+G_01)+((dy_th_3*h_12+dy_th_2*(h_11+h_01))*fields->A_2+(dx_th_3*h_12+dx_th_2*(h_11+h_01))*fields->A_1)*fields->A_3-(dz_th_3*h_12+dz_th_2*(h_11+h_01))*(fields->A_2*fields->A_2+fields->A_1*fields->A_1))/20;
                                      drdth3 -= -((h_22+h_12+h_02)*(dy_th_3*fields->A_2+dx_th_3*fields->A_1)*fields->A_3-dz_th_3*(G_22+G_12+G_02+(h_22+h_12+h_02)*(fields->A_2*fields->A_2+fields->A_1*fields->A_1)))/20;
                                      break;
                         }
                         break;
                case 1:  switch(dk)
                         {
                             case -2: drdth1 -= -(-dz_th_3*G_02-dz_th_2*G_01-dz_th_1*G_00+((dy_th_3*h_02+dy_th_2*h_01+dy_th_1*h_00)*fields->A_2+(dx_th_3*h_02+dx_th_2*h_01+dx_th_1*h_00)*fields->A_1)*fields->A_3-(dz_th_3*h_02+dz_th_2*h_01+dz_th_1*h_00)*(fields->A_2*fields->A_2+fields->A_1*fields->A_1))/30;
                                      drdth2 -= -(-dz_th_3*G_12-dz_th_2*(G_11+G_01)+((dy_th_3*h_12+dy_th_2*(h_11+h_01))*fields->A_2+(dx_th_3*h_12+dx_th_2*(h_11+h_01))*fields->A_1)*fields->A_3-(dz_th_3*h_12+dz_th_2*(h_11+h_01))*(fields->A_2*fields->A_2+fields->A_1*fields->A_1))/30;
                                      drdth3 -= -((h_22+h_12+h_02)*(dy_th_3*fields->A_2+dx_th_3*fields->A_1)*fields->A_3-dz_th_3*(G_22+G_12+G_02+(h_22+h_12+h_02)*(fields->A_2*fields->A_2+fields->A_1*fields->A_1)))/30;
                                      break;
                             case -1: drdth1 -= (-dz_th_3*G_02+dy_th_3*G_02-dz_th_2*G_01+dy_th_2*G_01-dz_th_1*G_00 +dy_th_1*G_00 +fields->A_3*((dy_th_3*h_02+dy_th_2*h_01+dy_th_1*h_00)*fields->A_3 -(dz_th_3*h_02-dy_th_3*h_02+dz_th_2*h_01-dy_th_2*h_01 +dz_th_1*h_00-dy_th_1*h_00) *fields->A_2+(dx_th_3*h_02+dx_th_2*h_01+dx_th_1*h_00)*fields->A_1) -(dz_th_3*h_02+dz_th_2*h_01+dz_th_1*h_00)*fields->A_2*fields->A_2 -fields->A_1*((dx_th_3*h_02+dx_th_2*h_01+dx_th_1*h_00)*fields->A_2 +(dz_th_3*h_02-dy_th_3*h_02+dz_th_2*h_01-dy_th_2*h_01 +dz_th_1*h_00-dy_th_1*h_00) *fields->A_1)) /12;
                                      drdth2 -= (-dz_th_3*G_12+dy_th_3*G_12-dz_th_2*(G_11+G_01)+dy_th_2*(G_11+G_01) +fields->A_3*((dy_th_3*h_12+dy_th_2*(h_11+h_01))*fields->A_3 -(dz_th_3*h_12-dy_th_3*h_12+dz_th_2*(h_11+h_01) -dy_th_2*(h_11+h_01)) *fields->A_2+(dx_th_3*h_12+dx_th_2*(h_11+h_01))*fields->A_1) -(dz_th_3*h_12+dz_th_2*(h_11+h_01))*fields->A_2*fields->A_2 -fields->A_1*((dx_th_3*h_12+dx_th_2*(h_11+h_01))*fields->A_2 +(dz_th_3*h_12-dy_th_3*h_12+dz_th_2*(h_11+h_01) -dy_th_2*(h_11+h_01)) *fields->A_1)) /12;
                                      drdth3 -= ((dy_th_3-dz_th_3)*(G_22+G_12+G_02+(h_22+h_12+h_02)*fields->A_1*fields->A_1) +(h_22+h_12+h_02)*(fields->A_3*(dy_th_3*fields->A_3+(dy_th_3-dz_th_3)*fields->A_2+dx_th_3*fields->A_1) -fields->A_2*(dz_th_3*fields->A_2+dx_th_3*fields->A_1))) /12;
                                      break;
                             case 0:  drdth1 -= 4*(dy_th_3*G_02+dy_th_2*G_01+dy_th_1*G_00 +fields->A_3*((dy_th_3*h_02+dy_th_2*h_01+dy_th_1*h_00)*fields->A_3 -(dz_th_3*h_02+dz_th_2*h_01+dz_th_1*h_00)*fields->A_2) +fields->A_1*((dy_th_3*h_02+dy_th_2*h_01+dy_th_1*h_00)*fields->A_1 -(dx_th_3*h_02+dx_th_2*h_01+dx_th_1*h_00)*fields->A_2)) /15;
                                      drdth2 -= 4*(dy_th_3*G_12+dy_th_2*(G_11+G_01) +fields->A_3*((dy_th_3*h_12+dy_th_2*(h_11+h_01))*fields->A_3 -(dz_th_3*h_12+dz_th_2*(h_11+h_01))*fields->A_2) +fields->A_1*((dy_th_3*h_12+dy_th_2*(h_11+h_01))*fields->A_1 -(dx_th_3*h_12+dx_th_2*(h_11+h_01))*fields->A_2)) /15;
                                      drdth3 -= 4*(dy_th_3*(G_22+G_12+G_02+(h_22+h_12+h_02)*fields->A_1*fields->A_1) +(h_22+h_12+h_02)*(dy_th_3*fields->A_3*fields->A_3-fields->A_2*(dz_th_3*fields->A_3+dx_th_3*fields->A_1))) /15;
                                      break;
                             case 1:  drdth1 -= (dz_th_3*G_02+dy_th_3*G_02+dz_th_2*G_01+dy_th_2*G_01+dz_th_1*G_00+dy_th_1*G_00 +fields->A_3*((dy_th_3*h_02+dy_th_2*h_01+dy_th_1*h_00)*fields->A_3 -(dz_th_3*h_02+dy_th_3*h_02+dz_th_2*h_01+dy_th_2*h_01 +dz_th_1*h_00+dy_th_1*h_00) *fields->A_2-(dx_th_3*h_02+dx_th_2*h_01+dx_th_1*h_00)*fields->A_1) +(dz_th_3*h_02+dz_th_2*h_01+dz_th_1*h_00)*fields->A_2*fields->A_2 +fields->A_1*((dz_th_3*h_02+dy_th_3*h_02+dz_th_2*h_01+dy_th_2*h_01 +dz_th_1*h_00+dy_th_1*h_00) *fields->A_1 -(dx_th_3*h_02+dx_th_2*h_01+dx_th_1*h_00)*fields->A_2)) /12;
                                      drdth2 -= (dz_th_3*G_12+dy_th_3*G_12+dz_th_2*(G_11+G_01)+dy_th_2*(G_11+G_01) +fields->A_3*((dy_th_3*h_12+dy_th_2*(h_11+h_01))*fields->A_3 -(dz_th_3*h_12+dy_th_3*h_12+dz_th_2*(h_11+h_01) +dy_th_2*(h_11+h_01)) *fields->A_2-(dx_th_3*h_12+dx_th_2*(h_11+h_01))*fields->A_1) +(dz_th_3*h_12+dz_th_2*(h_11+h_01))*fields->A_2*fields->A_2 +fields->A_1*((dz_th_3*h_12+dy_th_3*h_12+dz_th_2*(h_11+h_01) +dy_th_2*(h_11+h_01)) *fields->A_1 -(dx_th_3*h_12+dx_th_2*(h_11+h_01))*fields->A_2)) /12;
                                      drdth3 -= ((dz_th_3+dy_th_3)*(G_22+G_12+G_02+(h_22+h_12+h_02)*fields->A_1*fields->A_1) +(h_22+h_12+h_02)*(fields->A_3*(dy_th_3*fields->A_3-((dz_th_3+dy_th_3)*fields->A_2+dx_th_3*fields->A_1)) +fields->A_2*(dz_th_3*fields->A_2-dx_th_3*fields->A_1))) /12;
                                      break;
                             case 2:  drdth1 -= (-dz_th_3*G_02-dz_th_2*G_01-dz_th_1*G_00 +((dy_th_3*h_02+dy_th_2*h_01+dy_th_1*h_00)*fields->A_2 +(dx_th_3*h_02+dx_th_2*h_01+dx_th_1*h_00)*fields->A_1) *fields->A_3-(dz_th_3*h_02+dz_th_2*h_01+dz_th_1*h_00)*(fields->A_2*fields->A_2+fields->A_1*fields->A_1)) /30;
                                      drdth2 -= (-dz_th_3*G_12-dz_th_2*(G_11+G_01) +((dy_th_3*h_12+dy_th_2*(h_11+h_01))*fields->A_2 +(dx_th_3*h_12+dx_th_2*(h_11+h_01))*fields->A_1) *fields->A_3-(dz_th_3*h_12+dz_th_2*(h_11+h_01))*(fields->A_2*fields->A_2+fields->A_1*fields->A_1)) /30;
                                      drdth3 -= ((h_22+h_12+h_02)*(dy_th_3*fields->A_2+dx_th_3*fields->A_1)*fields->A_3 -dz_th_3*(G_22+G_12+G_02+(h_22+h_12+h_02)*(fields->A_2*fields->A_2+fields->A_1*fields->A_1))) /30;
                                      break;
                         }
                         break;
                case 2:  switch(dk)
                         {
                             case -1: drdth1 -= -(dy_th_3*G_02+dy_th_2*G_01+dy_th_1*G_00 +fields->A_3*((dy_th_3*h_02+dy_th_2*h_01+dy_th_1*h_00)*fields->A_3 -(dz_th_3*h_02+dz_th_2*h_01+dz_th_1*h_00)*fields->A_2) +fields->A_1*((dy_th_3*h_02+dy_th_2*h_01+dy_th_1*h_00)*fields->A_1 -(dx_th_3*h_02+dx_th_2*h_01+dx_th_1*h_00)*fields->A_2)) /30;
                                      drdth2 -= -(dy_th_3*G_12+dy_th_2*(G_11+G_01) +fields->A_3*((dy_th_3*h_12+dy_th_2*(h_11+h_01))*fields->A_3 -(dz_th_3*h_12+dz_th_2*(h_11+h_01))*fields->A_2) +fields->A_1*((dy_th_3*h_12+dy_th_2*(h_11+h_01))*fields->A_1 -(dx_th_3*h_12+dx_th_2*(h_11+h_01))*fields->A_2)) /30;
                                      drdth3 -= -(dy_th_3*(G_22+G_12+G_02+(h_22+h_12+h_02)*fields->A_1*fields->A_1) +(h_22+h_12+h_02)*(dy_th_3*fields->A_3*fields->A_3-fields->A_2*(dz_th_3*fields->A_3+dx_th_3*fields->A_1))) /30;
                                      break;
                             case 0:  drdth1 -= (dy_th_3*G_02+dy_th_2*G_01+dy_th_1*G_00 +fields->A_3*((dy_th_3*h_02+dy_th_2*h_01+dy_th_1*h_00)*fields->A_3 -(dz_th_3*h_02+dz_th_2*h_01+dz_th_1*h_00)*fields->A_2) +fields->A_1*((dy_th_3*h_02+dy_th_2*h_01+dy_th_1*h_00)*fields->A_1 -(dx_th_3*h_02+dx_th_2*h_01+dx_th_1*h_00)*fields->A_2)) /20;
                                      drdth2 -= (dy_th_3*G_12+dy_th_2*(G_11+G_01) +fields->A_3*((dy_th_3*h_12+dy_th_2*(h_11+h_01))*fields->A_3 -(dz_th_3*h_12+dz_th_2*(h_11+h_01))*fields->A_2) +fields->A_1*((dy_th_3*h_12+dy_th_2*(h_11+h_01))*fields->A_1 -(dx_th_3*h_12+dx_th_2*(h_11+h_01))*fields->A_2)) /20;
                                      drdth3 -= (dy_th_3*(G_22+G_12+G_02+(h_22+h_12+h_02)*fields->A_1*fields->A_1) +(h_22+h_12+h_02)*(dy_th_3*fields->A_3*fields->A_3-fields->A_2*(dz_th_3*fields->A_3+dx_th_3*fields->A_1))) /20;
                                      break;
                             case 1:  drdth1 -= -(dy_th_3*G_02+dy_th_2*G_01+dy_th_1*G_00 +fields->A_3*((dy_th_3*h_02+dy_th_2*h_01+dy_th_1*h_00)*fields->A_3 -(dz_th_3*h_02+dz_th_2*h_01+dz_th_1*h_00)*fields->A_2) +fields->A_1*((dy_th_3*h_02+dy_th_2*h_01+dy_th_1*h_00)*fields->A_1 -(dx_th_3*h_02+dx_th_2*h_01+dx_th_1*h_00)*fields->A_2)) /30;
                                      drdth2 -= -(dy_th_3*G_12+dy_th_2*(G_11+G_01) +fields->A_3*((dy_th_3*h_12+dy_th_2*(h_11+h_01))*fields->A_3 -(dz_th_3*h_12+dz_th_2*(h_11+h_01))*fields->A_2) +fields->A_1*((dy_th_3*h_12+dy_th_2*(h_11+h_01))*fields->A_1 -(dx_th_3*h_12+dx_th_2*(h_11+h_01))*fields->A_2)) /30;
                                      drdth3 -= -(dy_th_3*(G_22+G_12+G_02+(h_22+h_12+h_02)*fields->A_1*fields->A_1) +(h_22+h_12+h_02)*(dy_th_3*fields->A_3*fields->A_3-fields->A_2*(dz_th_3*fields->A_3+dx_th_3*fields->A_1))) /30;
                                      break;
                         }
                         break;
              }
              break;
     case 1:  switch(dj)
              {
                case -2: switch(dk)
                         {
                             case -1: drdth1 -= -(-dz_th_3*G_02+dx_th_3*G_02-dz_th_2*G_01+dx_th_2*G_01-dz_th_1*G_00 +dx_th_1*G_00 +fields->A_3*((dx_th_3*h_02+dx_th_2*h_01+dx_th_1*h_00)*fields->A_3 +(dy_th_3*h_02+dy_th_2*h_01+dy_th_1*h_00)*fields->A_2 -(dz_th_3*h_02-dx_th_3*h_02+dz_th_2*h_01-dx_th_2*h_01 +dz_th_1*h_00-dx_th_1*h_00) *fields->A_1) -(dz_th_3*h_02-dx_th_3*h_02+dz_th_2*h_01-dx_th_2*h_01 +dz_th_1*h_00-dx_th_1*h_00) *fields->A_2*fields->A_2 -fields->A_1*((dy_th_3*h_02+dy_th_2*h_01+dy_th_1*h_00)*fields->A_2 +(dz_th_3*h_02+dz_th_2*h_01+dz_th_1*h_00)*fields->A_1)) /120;
                                      drdth2 -= -(-dz_th_3*G_12+dx_th_3*G_12-dz_th_2*(G_11+G_01)+dx_th_2*(G_11+G_01) +fields->A_3*((dx_th_3*h_12+dx_th_2*(h_11+h_01))*fields->A_3 +(dy_th_3*h_12+dy_th_2*(h_11+h_01))*fields->A_2 -(dz_th_3*h_12-dx_th_3*h_12+dz_th_2*(h_11+h_01) -dx_th_2*(h_11+h_01)) *fields->A_1) -(dz_th_3*h_12-dx_th_3*h_12+dz_th_2*(h_11+h_01) -dx_th_2*(h_11+h_01)) *fields->A_2*fields->A_2 -fields->A_1*((dy_th_3*h_12+dy_th_2*(h_11+h_01))*fields->A_2 +(dz_th_3*h_12+dz_th_2*(h_11+h_01))*fields->A_1)) /120;
                                      drdth3 -= -((dx_th_3-dz_th_3)*(G_22+G_12+G_02)+(h_22+h_12+h_02) *(fields->A_3 *(dx_th_3*fields->A_3+dy_th_3*fields->A_2 +(dx_th_3-dz_th_3)*fields->A_1) +(dx_th_3-dz_th_3)*fields->A_2*fields->A_2 -fields->A_1*(dy_th_3*fields->A_2+dz_th_3*fields->A_1))) /120;
                                      break;
                             case 0:  drdth1 -= (dy_th_3*G_02+dy_th_2*G_01+dy_th_1*G_00 +fields->A_3*((dy_th_3*h_02+dy_th_2*h_01+dy_th_1*h_00)*fields->A_3 -(dz_th_3*h_02+dz_th_2*h_01+dz_th_1*h_00)*fields->A_2) +fields->A_1*((dy_th_3*h_02+dy_th_2*h_01+dy_th_1*h_00)*fields->A_1 -(dx_th_3*h_02+dx_th_2*h_01+dx_th_1*h_00)*fields->A_2)) /30;
                                      drdth2 -= (dy_th_3*G_12+dy_th_2*(G_11+G_01) +fields->A_3*((dy_th_3*h_12+dy_th_2*(h_11+h_01))*fields->A_3 -(dz_th_3*h_12+dz_th_2*(h_11+h_01))*fields->A_2) +fields->A_1*((dy_th_3*h_12+dy_th_2*(h_11+h_01))*fields->A_1 -(dx_th_3*h_12+dx_th_2*(h_11+h_01))*fields->A_2)) /30;
                                      drdth3 -= (dy_th_3*(G_22+G_12+G_02+(h_22+h_12+h_02)*fields->A_1*fields->A_1) +(h_22+h_12+h_02)*(dy_th_3*fields->A_3*fields->A_3-fields->A_2*(dz_th_3*fields->A_3+dx_th_3*fields->A_1))) /30;
                                      break;
                             case 1:  drdth1 -= -(dz_th_3*G_02+dx_th_3*G_02+dz_th_2*G_01+dx_th_2*G_01+dz_th_1*G_00 +dx_th_1*G_00 +fields->A_3*((dx_th_3*h_02+dx_th_2*h_01+dx_th_1*h_00)*fields->A_3 -(dy_th_3*h_02+dy_th_2*h_01+dy_th_1*h_00)*fields->A_2 -(dz_th_3*h_02+dx_th_3*h_02+dz_th_2*h_01+dx_th_2*h_01 +dz_th_1*h_00+dx_th_1*h_00) *fields->A_1) +(dz_th_3*h_02+dx_th_3*h_02+dz_th_2*h_01+dx_th_2*h_01 +dz_th_1*h_00+dx_th_1*h_00) *fields->A_2*fields->A_2 +fields->A_1*((dz_th_3*h_02+dz_th_2*h_01+dz_th_1*h_00)*fields->A_1 -(dy_th_3*h_02+dy_th_2*h_01+dy_th_1*h_00)*fields->A_2)) /120;
                                      drdth2 -= -(dz_th_3*G_12+dx_th_3*G_12+dz_th_2*(G_11+G_01)+dx_th_2*(G_11+G_01) +fields->A_3*((dx_th_3*h_12+dx_th_2*(h_11+h_01))*fields->A_3 -(dy_th_3*h_12+dy_th_2*(h_11+h_01))*fields->A_2 -(dz_th_3*h_12+dx_th_3*h_12+dz_th_2*(h_11+h_01) +dx_th_2*(h_11+h_01)) *fields->A_1) +(dz_th_3*h_12+dx_th_3*h_12+dz_th_2*(h_11+h_01) +dx_th_2*(h_11+h_01)) *fields->A_2*fields->A_2 +fields->A_1*((dz_th_3*h_12+dz_th_2*(h_11+h_01))*fields->A_1 -(dy_th_3*h_12+dy_th_2*(h_11+h_01))*fields->A_2)) /120;
                                      drdth3 -= -((dz_th_3+dx_th_3)*(G_22+G_12+G_02)+(h_22+h_12+h_02) *(fields->A_3 *(dx_th_3*fields->A_3 -(dy_th_3*fields->A_2+(dz_th_3+dx_th_3)*fields->A_1)) +(dz_th_3+dx_th_3)*fields->A_2*fields->A_2 +fields->A_1*(dz_th_3*fields->A_1-dy_th_3*fields->A_2))) /120;
                                      break;
                         }
                         break;
                case -1: switch(dk)
                         {
                             case -2: drdth1 -= (dy_th_3*G_02-dx_th_3*G_02+dy_th_2*G_01-dx_th_2*G_01+dy_th_1*G_00-dx_th_1*G_00 +fields->A_3*((dy_th_3*h_02-dx_th_3*h_02+dy_th_2*h_01-dx_th_2*h_01 +dy_th_1*h_00-dx_th_1*h_00) *fields->A_3 +(dz_th_3*h_02+dz_th_2*h_01+dz_th_1*h_00)*(fields->A_1-fields->A_2)) -(dx_th_3*h_02+dx_th_2*h_01+dx_th_1*h_00)*fields->A_2*fields->A_2 +fields->A_1*((dy_th_3*h_02-dx_th_3*h_02+dy_th_2*h_01-dx_th_2*h_01 +dy_th_1*h_00-dx_th_1*h_00) *fields->A_2 +(dy_th_3*h_02+dy_th_2*h_01+dy_th_1*h_00)*fields->A_1)) /120;
                                      drdth2 -= (dy_th_3*G_12-dx_th_3*G_12+dy_th_2*(G_11+G_01)-dx_th_2*(G_11+G_01) +fields->A_3*((dy_th_3*h_12-dx_th_3*h_12+dy_th_2*(h_11+h_01) -dx_th_2*(h_11+h_01)) *fields->A_3 +(dz_th_3*h_12+dz_th_2*(h_11+h_01))*(fields->A_1-fields->A_2)) -(dx_th_3*h_12+dx_th_2*(h_11+h_01))*fields->A_2*fields->A_2 +fields->A_1*((dy_th_3*h_12-dx_th_3*h_12+dy_th_2*(h_11+h_01) -dx_th_2*(h_11+h_01)) *fields->A_2 +(dy_th_3*h_12+dy_th_2*(h_11+h_01))*fields->A_1)) /120;
                                      drdth3 -= ((dy_th_3-dx_th_3)*(G_22+G_12+G_02)+(h_22+h_12+h_02) *(fields->A_3 *((dy_th_3-dx_th_3)*fields->A_3 +dz_th_3*(fields->A_1-fields->A_2)) -dx_th_3*fields->A_2*fields->A_2 +fields->A_1 *((dy_th_3-dx_th_3)*fields->A_2+dy_th_3*fields->A_1))) /120;
                                      break;
                             case -1: drdth1 -= -(dz_th_3*G_02+dy_th_3*G_02-dx_th_3*G_02+dz_th_2*G_01+dy_th_2*G_01 -dx_th_2*G_01+dz_th_1*G_00+dy_th_1*G_00-dx_th_1*G_00 +fields->A_3*((dy_th_3*h_02-dx_th_3*h_02+dy_th_2*h_01-dx_th_2*h_01 +dy_th_1*h_00-dx_th_1*h_00) *fields->A_3 -(dz_th_3*h_02+dy_th_3*h_02+dz_th_2*h_01+dy_th_2*h_01 +dz_th_1*h_00+dy_th_1*h_00) *fields->A_2 +(dz_th_3*h_02-dx_th_3*h_02+dz_th_2*h_01-dx_th_2*h_01 +dz_th_1*h_00-dx_th_1*h_00) *fields->A_1) +(dz_th_3*h_02-dx_th_3*h_02+dz_th_2*h_01-dx_th_2*h_01 +dz_th_1*h_00-dx_th_1*h_00) *fields->A_2*fields->A_2 +fields->A_1*((dy_th_3*h_02-dx_th_3*h_02+dy_th_2*h_01-dx_th_2*h_01 +dy_th_1*h_00-dx_th_1*h_00) *fields->A_2 +(dz_th_3*h_02+dy_th_3*h_02+dz_th_2*h_01+dy_th_2*h_01 +dz_th_1*h_00+dy_th_1*h_00) *fields->A_1)) /30;
                                      drdth2 -= -(dz_th_3*G_12+dy_th_3*G_12-dx_th_3*G_12+dz_th_2*(G_11+G_01) +dy_th_2*(G_11+G_01)-dx_th_2*(G_11+G_01) +fields->A_3*((dy_th_3*h_12-dx_th_3*h_12+dy_th_2*(h_11+h_01) -dx_th_2*(h_11+h_01)) *fields->A_3 -(dz_th_3*h_12+dy_th_3*h_12+dz_th_2*(h_11+h_01) +dy_th_2*(h_11+h_01)) *fields->A_2 +(dz_th_3*h_12-dx_th_3*h_12+dz_th_2*(h_11+h_01) -dx_th_2*(h_11+h_01)) *fields->A_1) +(dz_th_3*h_12-dx_th_3*h_12+dz_th_2*(h_11+h_01) -dx_th_2*(h_11+h_01)) *fields->A_2*fields->A_2 +fields->A_1*((dy_th_3*h_12-dx_th_3*h_12+dy_th_2*(h_11+h_01) -dx_th_2*(h_11+h_01)) *fields->A_2 +(dz_th_3*h_12+dy_th_3*h_12+dz_th_2*(h_11+h_01) +dy_th_2*(h_11+h_01)) *fields->A_1)) /30;
                                      drdth3 -= -((dz_th_3+dy_th_3-dx_th_3)*(G_22+G_12+G_02) +(h_22+h_12+h_02)*(fields->A_3*((dy_th_3-dx_th_3)*fields->A_3-(dz_th_3+dy_th_3)*fields->A_2 +(dz_th_3-dx_th_3)*fields->A_1) +(dz_th_3-dx_th_3)*fields->A_2*fields->A_2 +fields->A_1*((dy_th_3-dx_th_3)*fields->A_2+(dz_th_3+dy_th_3)*fields->A_1))) /30;
                                      break;
                             case 0:  drdth1 -= -(dy_th_3*G_02-dx_th_3*G_02+dy_th_2*G_01-dx_th_2*G_01+dy_th_1*G_00 -dx_th_1*G_00 +fields->A_3*((dy_th_3*h_02-dx_th_3*h_02+dy_th_2*h_01-dx_th_2*h_01 +dy_th_1*h_00-dx_th_1*h_00) *fields->A_3 +(dz_th_3*h_02+dz_th_2*h_01+dz_th_1*h_00)*(fields->A_1-fields->A_2)) -(dx_th_3*h_02+dx_th_2*h_01+dx_th_1*h_00)*fields->A_2*fields->A_2 +fields->A_1*((dy_th_3*h_02-dx_th_3*h_02+dy_th_2*h_01-dx_th_2*h_01 +dy_th_1*h_00-dx_th_1*h_00) *fields->A_2 +(dy_th_3*h_02+dy_th_2*h_01+dy_th_1*h_00)*fields->A_1)) /12;
                                      drdth2 -= -(dy_th_3*G_12-dx_th_3*G_12+dy_th_2*(G_11+G_01)-dx_th_2*(G_11+G_01) +fields->A_3*((dy_th_3*h_12-dx_th_3*h_12+dy_th_2*(h_11+h_01) -dx_th_2*(h_11+h_01)) *fields->A_3 +(dz_th_3*h_12+dz_th_2*(h_11+h_01))*(fields->A_1-fields->A_2)) -(dx_th_3*h_12+dx_th_2*(h_11+h_01))*fields->A_2*fields->A_2 +fields->A_1*((dy_th_3*h_12-dx_th_3*h_12+dy_th_2*(h_11+h_01) -dx_th_2*(h_11+h_01)) *fields->A_2 +(dy_th_3*h_12+dy_th_2*(h_11+h_01))*fields->A_1)) /12;
                                      drdth3 -= -((dy_th_3-dx_th_3)*(G_22+G_12+G_02)+(h_22+h_12+h_02) *(fields->A_3 *((dy_th_3-dx_th_3)*fields->A_3 +dz_th_3*(fields->A_1-fields->A_2)) -dx_th_3*fields->A_2*fields->A_2 +fields->A_1 *((dy_th_3-dx_th_3)*fields->A_2+dy_th_3*fields->A_1))) /12;
                                      break;
                             case 1:  drdth1 -= -(-dz_th_3*G_02+dy_th_3*G_02-dx_th_3*G_02-dz_th_2*G_01+dy_th_2*G_01 -dx_th_2*G_01-dz_th_1*G_00+dy_th_1*G_00-dx_th_1*G_00 +fields->A_3*((dy_th_3*h_02-dx_th_3*h_02+dy_th_2*h_01-dx_th_2*h_01 +dy_th_1*h_00-dx_th_1*h_00) *fields->A_3 -(dz_th_3*h_02-dy_th_3*h_02+dz_th_2*h_01-dy_th_2*h_01 +dz_th_1*h_00-dy_th_1*h_00) *fields->A_2 +(dz_th_3*h_02+dx_th_3*h_02+dz_th_2*h_01+dx_th_2*h_01 +dz_th_1*h_00+dx_th_1*h_00) *fields->A_1) +fields->A_2*((dy_th_3*h_02-dx_th_3*h_02+dy_th_2*h_01-dx_th_2*h_01 +dy_th_1*h_00-dx_th_1*h_00) *fields->A_1 -(dz_th_3*h_02+dx_th_3*h_02+dz_th_2*h_01+dx_th_2*h_01 +dz_th_1*h_00+dx_th_1*h_00) *fields->A_2) -(dz_th_3*h_02-dy_th_3*h_02+dz_th_2*h_01-dy_th_2*h_01 +dz_th_1*h_00-dy_th_1*h_00) *fields->A_1*fields->A_1) /30;
                                      drdth2 -= -(-dz_th_3*G_12+dy_th_3*G_12-dx_th_3*G_12-dz_th_2*(G_11+G_01) +dy_th_2*(G_11+G_01)-dx_th_2*(G_11+G_01) +fields->A_3*((dy_th_3*h_12-dx_th_3*h_12+dy_th_2*(h_11+h_01) -dx_th_2*(h_11+h_01)) *fields->A_3 -(dz_th_3*h_12-dy_th_3*h_12+dz_th_2*(h_11+h_01) -dy_th_2*(h_11+h_01)) *fields->A_2 +(dz_th_3*h_12+dx_th_3*h_12+dz_th_2*(h_11+h_01) +dx_th_2*(h_11+h_01)) *fields->A_1) +fields->A_2*((dy_th_3*h_12-dx_th_3*h_12+dy_th_2*(h_11+h_01) -dx_th_2*(h_11+h_01)) *fields->A_1 -(dz_th_3*h_12+dx_th_3*h_12+dz_th_2*(h_11+h_01) +dx_th_2*(h_11+h_01)) *fields->A_2) -(dz_th_3*h_12-dy_th_3*h_12+dz_th_2*(h_11+h_01) -dy_th_2*(h_11+h_01)) *fields->A_1*fields->A_1) /30;
                                      drdth3 -= -((-dz_th_3+dy_th_3-dx_th_3)*(G_22+G_12+G_02) +(h_22+h_12+h_02)*(fields->A_3*((dy_th_3-dx_th_3)*fields->A_3+(dy_th_3-dz_th_3)*fields->A_2 +(dz_th_3+dx_th_3)*fields->A_1) -(dz_th_3+dx_th_3)*fields->A_2*fields->A_2 +fields->A_1*((dy_th_3-dx_th_3)*fields->A_2+(dy_th_3-dz_th_3)*fields->A_1))) /30;
                                      break;
                             case 2:  drdth1 -= (dy_th_3*G_02-dx_th_3*G_02+dy_th_2*G_01-dx_th_2*G_01+dy_th_1*G_00-dx_th_1*G_00 +fields->A_3*((dy_th_3*h_02-dx_th_3*h_02+dy_th_2*h_01-dx_th_2*h_01 +dy_th_1*h_00-dx_th_1*h_00) *fields->A_3 +(dz_th_3*h_02+dz_th_2*h_01+dz_th_1*h_00)*(fields->A_1-fields->A_2)) -(dx_th_3*h_02+dx_th_2*h_01+dx_th_1*h_00)*fields->A_2*fields->A_2 +fields->A_1*((dy_th_3*h_02-dx_th_3*h_02+dy_th_2*h_01-dx_th_2*h_01 +dy_th_1*h_00-dx_th_1*h_00) *fields->A_2 +(dy_th_3*h_02+dy_th_2*h_01+dy_th_1*h_00)*fields->A_1)) /120;
                                      drdth2 -= (dy_th_3*G_12-dx_th_3*G_12+dy_th_2*(G_11+G_01)-dx_th_2*(G_11+G_01) +fields->A_3*((dy_th_3*h_12-dx_th_3*h_12+dy_th_2*(h_11+h_01) -dx_th_2*(h_11+h_01)) *fields->A_3 +(dz_th_3*h_12+dz_th_2*(h_11+h_01))*(fields->A_1-fields->A_2)) -(dx_th_3*h_12+dx_th_2*(h_11+h_01))*fields->A_2*fields->A_2 +fields->A_1*((dy_th_3*h_12-dx_th_3*h_12+dy_th_2*(h_11+h_01) -dx_th_2*(h_11+h_01)) *fields->A_2 +(dy_th_3*h_12+dy_th_2*(h_11+h_01))*fields->A_1)) /120;
                                      drdth3 -= ((dy_th_3-dx_th_3)*(G_22+G_12+G_02)+(h_22+h_12+h_02) *(fields->A_3 *((dy_th_3-dx_th_3)*fields->A_3 +dz_th_3*(fields->A_1-fields->A_2)) -dx_th_3*fields->A_2*fields->A_2 +fields->A_1 *((dy_th_3-dx_th_3)*fields->A_2+dy_th_3*fields->A_1))) /120;
                                      break;
                         }
                         break;
                case 0:  switch(dk)
                         {
                             case -2: drdth1 -= -(-dz_th_3*G_02-dz_th_2*G_01-dz_th_1*G_00 +((dy_th_3*h_02+dy_th_2*h_01+dy_th_1*h_00)*fields->A_2 +(dx_th_3*h_02+dx_th_2*h_01+dx_th_1*h_00)*fields->A_1) *fields->A_3-(dz_th_3*h_02+dz_th_2*h_01+dz_th_1*h_00)*(fields->A_2*fields->A_2+fields->A_1*fields->A_1)) /30;
                                      drdth2 -= -(-dz_th_3*G_12-dz_th_2*(G_11+G_01) +((dy_th_3*h_12+dy_th_2*(h_11+h_01))*fields->A_2 +(dx_th_3*h_12+dx_th_2*(h_11+h_01))*fields->A_1) *fields->A_3-(dz_th_3*h_12+dz_th_2*(h_11+h_01))*(fields->A_2*fields->A_2+fields->A_1*fields->A_1)) /30;
                                      drdth3 -= -((h_22+h_12+h_02)*(dy_th_3*fields->A_2+dx_th_3*fields->A_1)*fields->A_3 -dz_th_3*(G_22+G_12+G_02+(h_22+h_12+h_02)*(fields->A_2*fields->A_2+fields->A_1*fields->A_1))) /30;
                                      break;
                             case -1: drdth1 -= (-dz_th_3*G_02+dx_th_3*G_02-dz_th_2*G_01+dx_th_2*G_01-dz_th_1*G_00 +dx_th_1*G_00 +fields->A_3*((dx_th_3*h_02+dx_th_2*h_01+dx_th_1*h_00)*fields->A_3 +(dy_th_3*h_02+dy_th_2*h_01+dy_th_1*h_00)*fields->A_2 -(dz_th_3*h_02-dx_th_3*h_02+dz_th_2*h_01-dx_th_2*h_01 +dz_th_1*h_00-dx_th_1*h_00) *fields->A_1) -(dz_th_3*h_02-dx_th_3*h_02+dz_th_2*h_01-dx_th_2*h_01 +dz_th_1*h_00-dx_th_1*h_00) *fields->A_2*fields->A_2 -fields->A_1*((dy_th_3*h_02+dy_th_2*h_01+dy_th_1*h_00)*fields->A_2 +(dz_th_3*h_02+dz_th_2*h_01+dz_th_1*h_00)*fields->A_1)) /12;
                                      drdth2 -= (-dz_th_3*G_12+dx_th_3*G_12-dz_th_2*(G_11+G_01)+dx_th_2*(G_11+G_01) +fields->A_3*((dx_th_3*h_12+dx_th_2*(h_11+h_01))*fields->A_3 +(dy_th_3*h_12+dy_th_2*(h_11+h_01))*fields->A_2 -(dz_th_3*h_12-dx_th_3*h_12+dz_th_2*(h_11+h_01) -dx_th_2*(h_11+h_01)) *fields->A_1) -(dz_th_3*h_12-dx_th_3*h_12+dz_th_2*(h_11+h_01) -dx_th_2*(h_11+h_01)) *fields->A_2*fields->A_2 -fields->A_1*((dy_th_3*h_12+dy_th_2*(h_11+h_01))*fields->A_2 +(dz_th_3*h_12+dz_th_2*(h_11+h_01))*fields->A_1)) /12;
                                      drdth3 -= ((dx_th_3-dz_th_3)*(G_22+G_12+G_02)+(h_22+h_12+h_02) *(fields->A_3 *(dx_th_3*fields->A_3+dy_th_3*fields->A_2 +(dx_th_3-dz_th_3)*fields->A_1) +(dx_th_3-dz_th_3)*fields->A_2*fields->A_2 -fields->A_1*(dy_th_3*fields->A_2+dz_th_3*fields->A_1))) /12;
                                      break;
                             case 0:  drdth1 -= 4*(dx_th_3*G_02+dx_th_2*G_01+dx_th_1*G_00 +fields->A_3*((dx_th_3*h_02+dx_th_2*h_01+dx_th_1*h_00)*fields->A_3 -(dz_th_3*h_02+dz_th_2*h_01+dz_th_1*h_00)*fields->A_1) +fields->A_2*((dx_th_3*h_02+dx_th_2*h_01+dx_th_1*h_00)*fields->A_2 -(dy_th_3*h_02+dy_th_2*h_01+dy_th_1*h_00)*fields->A_1)) /15;
                                      drdth2 -= 4*(dx_th_3*G_12+dx_th_2*(G_11+G_01) +fields->A_3*((dx_th_3*h_12+dx_th_2*(h_11+h_01))*fields->A_3 -(dz_th_3*h_12+dz_th_2*(h_11+h_01))*fields->A_1) +fields->A_2*((dx_th_3*h_12+dx_th_2*(h_11+h_01))*fields->A_2 -(dy_th_3*h_12+dy_th_2*(h_11+h_01))*fields->A_1)) /15;
                                      drdth3 -= 4*(dx_th_3*(G_22+G_12+G_02)+(h_22+h_12+h_02) *(fields->A_3*(dx_th_3*fields->A_3-dz_th_3*fields->A_1) +fields->A_2*(dx_th_3*fields->A_2-dy_th_3*fields->A_1))) /15;
                                      break;
                             case 1:  drdth1 -= (dz_th_3*G_02+dx_th_3*G_02+dz_th_2*G_01+dx_th_2*G_01+dz_th_1*G_00+dx_th_1*G_00 +fields->A_3*((dx_th_3*h_02+dx_th_2*h_01+dx_th_1*h_00)*fields->A_3 -(dy_th_3*h_02+dy_th_2*h_01+dy_th_1*h_00)*fields->A_2 -(dz_th_3*h_02+dx_th_3*h_02+dz_th_2*h_01+dx_th_2*h_01 +dz_th_1*h_00+dx_th_1*h_00) *fields->A_1) +(dz_th_3*h_02+dx_th_3*h_02+dz_th_2*h_01+dx_th_2*h_01 +dz_th_1*h_00+dx_th_1*h_00) *fields->A_2*fields->A_2 +fields->A_1*((dz_th_3*h_02+dz_th_2*h_01+dz_th_1*h_00)*fields->A_1 -(dy_th_3*h_02+dy_th_2*h_01+dy_th_1*h_00)*fields->A_2)) /12;
                                      drdth2 -= (dz_th_3*G_12+dx_th_3*G_12+dz_th_2*(G_11+G_01)+dx_th_2*(G_11+G_01) +fields->A_3*((dx_th_3*h_12+dx_th_2*(h_11+h_01))*fields->A_3 -(dy_th_3*h_12+dy_th_2*(h_11+h_01))*fields->A_2 -(dz_th_3*h_12+dx_th_3*h_12+dz_th_2*(h_11+h_01) +dx_th_2*(h_11+h_01)) *fields->A_1) +(dz_th_3*h_12+dx_th_3*h_12+dz_th_2*(h_11+h_01) +dx_th_2*(h_11+h_01)) *fields->A_2*fields->A_2 +fields->A_1*((dz_th_3*h_12+dz_th_2*(h_11+h_01))*fields->A_1 -(dy_th_3*h_12+dy_th_2*(h_11+h_01))*fields->A_2)) /12;
                                      drdth3 -= ((dz_th_3+dx_th_3)*(G_22+G_12+G_02)+(h_22+h_12+h_02) *(fields->A_3 *(dx_th_3*fields->A_3 -(dy_th_3*fields->A_2+(dz_th_3+dx_th_3)*fields->A_1)) +(dz_th_3+dx_th_3)*fields->A_2*fields->A_2 +fields->A_1*(dz_th_3*fields->A_1-dy_th_3*fields->A_2))) /12;
                                      break;
                             case 2:  drdth1 -= (-dz_th_3*G_02-dz_th_2*G_01-dz_th_1*G_00 +((dy_th_3*h_02+dy_th_2*h_01+dy_th_1*h_00)*fields->A_2 +(dx_th_3*h_02+dx_th_2*h_01+dx_th_1*h_00)*fields->A_1) *fields->A_3-(dz_th_3*h_02+dz_th_2*h_01+dz_th_1*h_00)*(fields->A_2*fields->A_2+fields->A_1*fields->A_1)) /30;
                                      drdth2 -= (-dz_th_3*G_12-dz_th_2*(G_11+G_01) +((dy_th_3*h_12+dy_th_2*(h_11+h_01))*fields->A_2 +(dx_th_3*h_12+dx_th_2*(h_11+h_01))*fields->A_1) *fields->A_3-(dz_th_3*h_12+dz_th_2*(h_11+h_01))*(fields->A_2*fields->A_2+fields->A_1*fields->A_1)) /30;
                                      drdth3 -= ((h_22+h_12+h_02)*(dy_th_3*fields->A_2+dx_th_3*fields->A_1)*fields->A_3 -dz_th_3*(G_22+G_12+G_02+(h_22+h_12+h_02)*(fields->A_2*fields->A_2+fields->A_1*fields->A_1))) /30;
                                      break;
                         }
                         break;
                case 1:  switch(dk)
                         {
                             case -2: drdth1 -= -(dy_th_3*G_02+dx_th_3*G_02+dy_th_2*G_01+dx_th_2*G_01+dy_th_1*G_00 +dx_th_1*G_00 +fields->A_3*((dy_th_3*h_02+dx_th_3*h_02+dy_th_2*h_01+dx_th_2*h_01 +dy_th_1*h_00+dx_th_1*h_00) *fields->A_3 -(dz_th_3*h_02+dz_th_2*h_01+dz_th_1*h_00)*(fields->A_2+fields->A_1)) +(dx_th_3*h_02+dx_th_2*h_01+dx_th_1*h_00)*fields->A_2*fields->A_2 +fields->A_1*((dy_th_3*h_02+dy_th_2*h_01+dy_th_1*h_00)*fields->A_1 -(dy_th_3*h_02+dx_th_3*h_02+dy_th_2*h_01+dx_th_2*h_01 +dy_th_1*h_00+dx_th_1*h_00) *fields->A_2)) /120;
                                      drdth2 -= -(dy_th_3*G_12+dx_th_3*G_12+dy_th_2*(G_11+G_01)+dx_th_2*(G_11+G_01) +fields->A_3*((dy_th_3*h_12+dx_th_3*h_12+dy_th_2*(h_11+h_01) +dx_th_2*(h_11+h_01)) *fields->A_3 -(dz_th_3*h_12+dz_th_2*(h_11+h_01))*(fields->A_2+fields->A_1)) +(dx_th_3*h_12+dx_th_2*(h_11+h_01))*fields->A_2*fields->A_2 +fields->A_1*((dy_th_3*h_12+dy_th_2*(h_11+h_01))*fields->A_1 -(dy_th_3*h_12+dx_th_3*h_12+dy_th_2*(h_11+h_01) +dx_th_2*(h_11+h_01)) *fields->A_2)) /120;
                                      drdth3 -= -((dy_th_3+dx_th_3)*(G_22+G_12+G_02)+(h_22+h_12+h_02) *(fields->A_3 *((dy_th_3+dx_th_3)*fields->A_3 -dz_th_3*(fields->A_2+fields->A_1)) +dx_th_3*fields->A_2*fields->A_2 +fields->A_1 *(dy_th_3*fields->A_1-(dy_th_3+dx_th_3)*fields->A_2))) /120;
                                      break;
                             case -1: drdth1 -= (-dz_th_3*G_02+dy_th_3*G_02+dx_th_3*G_02-dz_th_2*G_01+dy_th_2*G_01 +dx_th_2*G_01-dz_th_1*G_00+dy_th_1*G_00+dx_th_1*G_00 +(dy_th_3*h_02+dx_th_3*h_02+dy_th_2*h_01+dx_th_2*h_01 +dy_th_1*h_00+dx_th_1*h_00) *fields->A_3*fields->A_3 -((dz_th_3*h_02-dy_th_3*h_02+dz_th_2*h_01-dy_th_2*h_01 +dz_th_1*h_00-dy_th_1*h_00) *fields->A_2 +(dz_th_3*h_02-dx_th_3*h_02+dz_th_2*h_01-dx_th_2*h_01 +dz_th_1*h_00-dx_th_1*h_00) *fields->A_1) *fields->A_3 -(dz_th_3*h_02-dx_th_3*h_02+dz_th_2*h_01-dx_th_2*h_01 +dz_th_1*h_00-dx_th_1*h_00) *fields->A_2*fields->A_2 -fields->A_1*((dy_th_3*h_02+dx_th_3*h_02+dy_th_2*h_01+dx_th_2*h_01 +dy_th_1*h_00+dx_th_1*h_00) *fields->A_2 +(dz_th_3*h_02-dy_th_3*h_02+dz_th_2*h_01-dy_th_2*h_01 +dz_th_1*h_00-dy_th_1*h_00) *fields->A_1)) /30;
                                      drdth2 -= (-dz_th_3*G_12+dy_th_3*G_12+dx_th_3*G_12-dz_th_2*(G_11+G_01) +dy_th_2*(G_11+G_01)+dx_th_2*(G_11+G_01) +(dy_th_3*h_12+dx_th_3*h_12+dy_th_2*(h_11+h_01) +dx_th_2*(h_11+h_01)) *fields->A_3*fields->A_3 -((dz_th_3*h_12-dy_th_3*h_12+dz_th_2*(h_11+h_01) -dy_th_2*(h_11+h_01)) *fields->A_2 +(dz_th_3*h_12-dx_th_3*h_12+dz_th_2*(h_11+h_01) -dx_th_2*(h_11+h_01)) *fields->A_1) *fields->A_3 -(dz_th_3*h_12-dx_th_3*h_12+dz_th_2*(h_11+h_01) -dx_th_2*(h_11+h_01)) *fields->A_2*fields->A_2 -fields->A_1*((dy_th_3*h_12+dx_th_3*h_12+dy_th_2*(h_11+h_01) +dx_th_2*(h_11+h_01)) *fields->A_2 +(dz_th_3*h_12-dy_th_3*h_12+dz_th_2*(h_11+h_01) -dy_th_2*(h_11+h_01)) *fields->A_1)) /30;
                                      drdth3 -= ((-dz_th_3+dy_th_3+dx_th_3)*(G_22+G_12+G_02) +(h_22+h_12+h_02)*(fields->A_3*((dy_th_3+dx_th_3)*fields->A_3+(dy_th_3-dz_th_3)*fields->A_2 +(dx_th_3-dz_th_3)*fields->A_1) +(dx_th_3-dz_th_3)*fields->A_2*fields->A_2 +fields->A_1*((dy_th_3-dz_th_3)*fields->A_1-(dy_th_3+dx_th_3)*fields->A_2))) /30;
                                      break;
                             case 0:  drdth1 -= (dy_th_3*G_02+dx_th_3*G_02+dy_th_2*G_01+dx_th_2*G_01+dy_th_1*G_00+dx_th_1*G_00 +fields->A_3*((dy_th_3*h_02+dx_th_3*h_02+dy_th_2*h_01+dx_th_2*h_01 +dy_th_1*h_00+dx_th_1*h_00) *fields->A_3 -(dz_th_3*h_02+dz_th_2*h_01+dz_th_1*h_00)*(fields->A_2+fields->A_1)) +(dx_th_3*h_02+dx_th_2*h_01+dx_th_1*h_00)*fields->A_2*fields->A_2 +fields->A_1*((dy_th_3*h_02+dy_th_2*h_01+dy_th_1*h_00)*fields->A_1 -(dy_th_3*h_02+dx_th_3*h_02+dy_th_2*h_01+dx_th_2*h_01 +dy_th_1*h_00+dx_th_1*h_00) *fields->A_2)) /12;
                                      drdth2 -= (dy_th_3*G_12+dx_th_3*G_12+dy_th_2*(G_11+G_01)+dx_th_2*(G_11+G_01) +fields->A_3*((dy_th_3*h_12+dx_th_3*h_12+dy_th_2*(h_11+h_01) +dx_th_2*(h_11+h_01)) *fields->A_3 -(dz_th_3*h_12+dz_th_2*(h_11+h_01))*(fields->A_2+fields->A_1)) +(dx_th_3*h_12+dx_th_2*(h_11+h_01))*fields->A_2*fields->A_2 +fields->A_1*((dy_th_3*h_12+dy_th_2*(h_11+h_01))*fields->A_1 -(dy_th_3*h_12+dx_th_3*h_12+dy_th_2*(h_11+h_01) +dx_th_2*(h_11+h_01)) *fields->A_2)) /12;
                                      drdth3 -= ((dy_th_3+dx_th_3)*(G_22+G_12+G_02)+(h_22+h_12+h_02) *(fields->A_3 *((dy_th_3+dx_th_3)*fields->A_3 -dz_th_3*(fields->A_2+fields->A_1)) +dx_th_3*fields->A_2*fields->A_2 +fields->A_1 *(dy_th_3*fields->A_1-(dy_th_3+dx_th_3)*fields->A_2))) /12;
                                      break;
                             case 1:  drdth1 -= (dz_th_3*G_02+dy_th_3*G_02+dx_th_3*G_02+dz_th_2*G_01+dy_th_2*G_01+dx_th_2*G_01 +dz_th_1*G_00+dy_th_1*G_00+dx_th_1*G_00 +fields->A_3*((dy_th_3*h_02+dx_th_3*h_02+dy_th_2*h_01+dx_th_2*h_01 +dy_th_1*h_00+dx_th_1*h_00) *fields->A_3 -(dz_th_3*h_02+dy_th_3*h_02+dz_th_2*h_01+dy_th_2*h_01 +dz_th_1*h_00+dy_th_1*h_00) *fields->A_2 -(dz_th_3*h_02+dx_th_3*h_02+dz_th_2*h_01+dx_th_2*h_01 +dz_th_1*h_00+dx_th_1*h_00) *fields->A_1) +(dz_th_3*h_02+dx_th_3*h_02+dz_th_2*h_01+dx_th_2*h_01 +dz_th_1*h_00+dx_th_1*h_00) *fields->A_2*fields->A_2 +fields->A_1*((dz_th_3*h_02+dy_th_3*h_02+dz_th_2*h_01+dy_th_2*h_01 +dz_th_1*h_00+dy_th_1*h_00) *fields->A_1 -(dy_th_3*h_02+dx_th_3*h_02+dy_th_2*h_01+dx_th_2*h_01 +dy_th_1*h_00+dx_th_1*h_00) *fields->A_2)) /30;
                                      drdth2 -= (dz_th_3*G_12+dy_th_3*G_12+dx_th_3*G_12+dz_th_2*(G_11+G_01) +dy_th_2*(G_11+G_01)+dx_th_2*(G_11+G_01) +fields->A_3*((dy_th_3*h_12+dx_th_3*h_12+dy_th_2*(h_11+h_01) +dx_th_2*(h_11+h_01)) *fields->A_3 -(dz_th_3*h_12+dy_th_3*h_12+dz_th_2*(h_11+h_01) +dy_th_2*(h_11+h_01)) *fields->A_2 -(dz_th_3*h_12+dx_th_3*h_12+dz_th_2*(h_11+h_01) +dx_th_2*(h_11+h_01)) *fields->A_1) +(dz_th_3*h_12+dx_th_3*h_12+dz_th_2*(h_11+h_01) +dx_th_2*(h_11+h_01)) *fields->A_2*fields->A_2 +fields->A_1*((dz_th_3*h_12+dy_th_3*h_12+dz_th_2*(h_11+h_01) +dy_th_2*(h_11+h_01)) *fields->A_1 -(dy_th_3*h_12+dx_th_3*h_12+dy_th_2*(h_11+h_01) +dx_th_2*(h_11+h_01)) *fields->A_2)) /30;
                                      drdth3 -= ((dz_th_3+dy_th_3+dx_th_3)*(G_22+G_12+G_02) +(h_22+h_12+h_02)*(fields->A_3*((dy_th_3+dx_th_3)*fields->A_3 -((dz_th_3+dy_th_3)*fields->A_2+(dz_th_3+dx_th_3)*fields->A_1)) +(dz_th_3+dx_th_3)*fields->A_2*fields->A_2 +fields->A_1*((dz_th_3+dy_th_3)*fields->A_1-(dy_th_3+dx_th_3)*fields->A_2))) /30;
                                      break;
                             case 2:  drdth1 -= -(dy_th_3*G_02+dx_th_3*G_02+dy_th_2*G_01+dx_th_2*G_01+dy_th_1*G_00 +dx_th_1*G_00 +fields->A_3*((dy_th_3*h_02+dx_th_3*h_02+dy_th_2*h_01+dx_th_2*h_01 +dy_th_1*h_00+dx_th_1*h_00) *fields->A_3 -(dz_th_3*h_02+dz_th_2*h_01+dz_th_1*h_00)*(fields->A_2+fields->A_1)) +(dx_th_3*h_02+dx_th_2*h_01+dx_th_1*h_00)*fields->A_2*fields->A_2 +fields->A_1*((dy_th_3*h_02+dy_th_2*h_01+dy_th_1*h_00)*fields->A_1 -(dy_th_3*h_02+dx_th_3*h_02+dy_th_2*h_01+dx_th_2*h_01 +dy_th_1*h_00+dx_th_1*h_00) *fields->A_2)) /120;
                                      drdth2 -= -(dy_th_3*G_12+dx_th_3*G_12+dy_th_2*(G_11+G_01)+dx_th_2*(G_11+G_01) +fields->A_3*((dy_th_3*h_12+dx_th_3*h_12+dy_th_2*(h_11+h_01) +dx_th_2*(h_11+h_01)) *fields->A_3 -(dz_th_3*h_12+dz_th_2*(h_11+h_01))*(fields->A_2+fields->A_1)) +(dx_th_3*h_12+dx_th_2*(h_11+h_01))*fields->A_2*fields->A_2 +fields->A_1*((dy_th_3*h_12+dy_th_2*(h_11+h_01))*fields->A_1 -(dy_th_3*h_12+dx_th_3*h_12+dy_th_2*(h_11+h_01) +dx_th_2*(h_11+h_01)) *fields->A_2)) /120;
                                      drdth3 -= -((dy_th_3+dx_th_3)*(G_22+G_12+G_02)+(h_22+h_12+h_02) *(fields->A_3 *((dy_th_3+dx_th_3)*fields->A_3 -dz_th_3*(fields->A_2+fields->A_1)) +dx_th_3*fields->A_2*fields->A_2 +fields->A_1 *(dy_th_3*fields->A_1-(dy_th_3+dx_th_3)*fields->A_2))) /120;
                                      break;
                         }
                         break;
                case 2:  switch(dk)
                         {
                             case -1: drdth1 -= -(-dz_th_3*G_02+dx_th_3*G_02-dz_th_2*G_01+dx_th_2*G_01-dz_th_1*G_00 +dx_th_1*G_00 +fields->A_3*((dx_th_3*h_02+dx_th_2*h_01+dx_th_1*h_00)*fields->A_3 +(dy_th_3*h_02+dy_th_2*h_01+dy_th_1*h_00)*fields->A_2 -(dz_th_3*h_02-dx_th_3*h_02+dz_th_2*h_01-dx_th_2*h_01 +dz_th_1*h_00-dx_th_1*h_00) *fields->A_1) -(dz_th_3*h_02-dx_th_3*h_02+dz_th_2*h_01-dx_th_2*h_01 +dz_th_1*h_00-dx_th_1*h_00) *fields->A_2*fields->A_2 -fields->A_1*((dy_th_3*h_02+dy_th_2*h_01+dy_th_1*h_00)*fields->A_2 +(dz_th_3*h_02+dz_th_2*h_01+dz_th_1*h_00)*fields->A_1)) /120;
                                      drdth2 -= -(-dz_th_3*G_12+dx_th_3*G_12-dz_th_2*(G_11+G_01)+dx_th_2*(G_11+G_01) +fields->A_3*((dx_th_3*h_12+dx_th_2*(h_11+h_01))*fields->A_3 +(dy_th_3*h_12+dy_th_2*(h_11+h_01))*fields->A_2 -(dz_th_3*h_12-dx_th_3*h_12+dz_th_2*(h_11+h_01) -dx_th_2*(h_11+h_01)) *fields->A_1) -(dz_th_3*h_12-dx_th_3*h_12+dz_th_2*(h_11+h_01) -dx_th_2*(h_11+h_01)) *fields->A_2*fields->A_2 -fields->A_1*((dy_th_3*h_12+dy_th_2*(h_11+h_01))*fields->A_2 +(dz_th_3*h_12+dz_th_2*(h_11+h_01))*fields->A_1)) /120;
                                      drdth3 -= -((dx_th_3-dz_th_3)*(G_22+G_12+G_02)+(h_22+h_12+h_02) *(fields->A_3 *(dx_th_3*fields->A_3+dy_th_3*fields->A_2 +(dx_th_3-dz_th_3)*fields->A_1) +(dx_th_3-dz_th_3)*fields->A_2*fields->A_2 -fields->A_1*(dy_th_3*fields->A_2+dz_th_3*fields->A_1))) /120;
                                      break;
                             case 0:  drdth1 -= -(dy_th_3*G_02+dy_th_2*G_01+dy_th_1*G_00 +fields->A_3*((dy_th_3*h_02+dy_th_2*h_01+dy_th_1*h_00)*fields->A_3 -(dz_th_3*h_02+dz_th_2*h_01+dz_th_1*h_00)*fields->A_2) +fields->A_1*((dy_th_3*h_02+dy_th_2*h_01+dy_th_1*h_00)*fields->A_1 -(dx_th_3*h_02+dx_th_2*h_01+dx_th_1*h_00)*fields->A_2)) /30;
                                      drdth2 -= -(dy_th_3*G_12+dy_th_2*(G_11+G_01) +fields->A_3*((dy_th_3*h_12+dy_th_2*(h_11+h_01))*fields->A_3 -(dz_th_3*h_12+dz_th_2*(h_11+h_01))*fields->A_2) +fields->A_1*((dy_th_3*h_12+dy_th_2*(h_11+h_01))*fields->A_1 -(dx_th_3*h_12+dx_th_2*(h_11+h_01))*fields->A_2)) /30;
                                      drdth3 -= -(dy_th_3*(G_22+G_12+G_02+(h_22+h_12+h_02)*fields->A_1*fields->A_1) +(h_22+h_12+h_02)*(dy_th_3*fields->A_3*fields->A_3-fields->A_2*(dz_th_3*fields->A_3+dx_th_3*fields->A_1))) /30;
                                      break;
                             case 1:  drdth1 -= -(dz_th_3*G_02+dx_th_3*G_02+dz_th_2*G_01+dx_th_2*G_01+dz_th_1*G_00 +dx_th_1*G_00 +fields->A_3*((dx_th_3*h_02+dx_th_2*h_01+dx_th_1*h_00)*fields->A_3 -(dy_th_3*h_02+dy_th_2*h_01+dy_th_1*h_00)*fields->A_2 -(dz_th_3*h_02+dx_th_3*h_02+dz_th_2*h_01+dx_th_2*h_01 +dz_th_1*h_00+dx_th_1*h_00) *fields->A_1) +(dz_th_3*h_02+dx_th_3*h_02+dz_th_2*h_01+dx_th_2*h_01 +dz_th_1*h_00+dx_th_1*h_00) *fields->A_2*fields->A_2 +fields->A_1*((dz_th_3*h_02+dz_th_2*h_01+dz_th_1*h_00)*fields->A_1 -(dy_th_3*h_02+dy_th_2*h_01+dy_th_1*h_00)*fields->A_2)) /120;
                                      drdth2 -= -(dz_th_3*G_12+dx_th_3*G_12+dz_th_2*(G_11+G_01)+dx_th_2*(G_11+G_01) +fields->A_3*((dx_th_3*h_12+dx_th_2*(h_11+h_01))*fields->A_3 -(dy_th_3*h_12+dy_th_2*(h_11+h_01))*fields->A_2 -(dz_th_3*h_12+dx_th_3*h_12+dz_th_2*(h_11+h_01) +dx_th_2*(h_11+h_01)) *fields->A_1) +(dz_th_3*h_12+dx_th_3*h_12+dz_th_2*(h_11+h_01) +dx_th_2*(h_11+h_01)) *fields->A_2*fields->A_2 +fields->A_1*((dz_th_3*h_12+dz_th_2*(h_11+h_01))*fields->A_1 -(dy_th_3*h_12+dy_th_2*(h_11+h_01))*fields->A_2)) /120;
                                      drdth3 -= -((dz_th_3+dx_th_3)*(G_22+G_12+G_02)+(h_22+h_12+h_02) *(fields->A_3 *(dx_th_3*fields->A_3 -(dy_th_3*fields->A_2+(dz_th_3+dx_th_3)*fields->A_1)) +(dz_th_3+dx_th_3)*fields->A_2*fields->A_2 +fields->A_1*(dz_th_3*fields->A_1-dy_th_3*fields->A_2))) /120;
                                      break;
                         }
                         break;
              }
              break;
     case 2:  switch(dj)
              {
                case -1: switch(dk)
                         {
                             case -1: drdth1 -= (dz_th_3*G_02+dy_th_3*G_02+dz_th_2*G_01+dy_th_2*G_01+dz_th_1*G_00+dy_th_1*G_00 +fields->A_3*((dy_th_3*h_02+dy_th_2*h_01+dy_th_1*h_00)*fields->A_3 -(dz_th_3*h_02+dy_th_3*h_02+dz_th_2*h_01+dy_th_2*h_01 +dz_th_1*h_00+dy_th_1*h_00) *fields->A_2-(dx_th_3*h_02+dx_th_2*h_01+dx_th_1*h_00)*fields->A_1) +(dz_th_3*h_02+dz_th_2*h_01+dz_th_1*h_00)*fields->A_2*fields->A_2 +fields->A_1*((dz_th_3*h_02+dy_th_3*h_02+dz_th_2*h_01+dy_th_2*h_01 +dz_th_1*h_00+dy_th_1*h_00) *fields->A_1 -(dx_th_3*h_02+dx_th_2*h_01+dx_th_1*h_00)*fields->A_2)) /120;
                                      drdth2 -= (dz_th_3*G_12+dy_th_3*G_12+dz_th_2*(G_11+G_01)+dy_th_2*(G_11+G_01) +fields->A_3*((dy_th_3*h_12+dy_th_2*(h_11+h_01))*fields->A_3 -(dz_th_3*h_12+dy_th_3*h_12+dz_th_2*(h_11+h_01) +dy_th_2*(h_11+h_01)) *fields->A_2-(dx_th_3*h_12+dx_th_2*(h_11+h_01))*fields->A_1) +(dz_th_3*h_12+dz_th_2*(h_11+h_01))*fields->A_2*fields->A_2 +fields->A_1*((dz_th_3*h_12+dy_th_3*h_12+dz_th_2*(h_11+h_01) +dy_th_2*(h_11+h_01)) *fields->A_1 -(dx_th_3*h_12+dx_th_2*(h_11+h_01))*fields->A_2)) /120;
                                      drdth3 -= ((dz_th_3+dy_th_3)*(G_22+G_12+G_02+(h_22+h_12+h_02)*fields->A_1*fields->A_1) +(h_22+h_12+h_02)*(fields->A_3*(dy_th_3*fields->A_3-((dz_th_3+dy_th_3)*fields->A_2+dx_th_3*fields->A_1)) +fields->A_2*(dz_th_3*fields->A_2-dx_th_3*fields->A_1))) /120;
                                      break;
                             case 0:  drdth1 -= -(dx_th_3*G_02+dx_th_2*G_01+dx_th_1*G_00 +fields->A_3*((dx_th_3*h_02+dx_th_2*h_01+dx_th_1*h_00)*fields->A_3 -(dz_th_3*h_02+dz_th_2*h_01+dz_th_1*h_00)*fields->A_1) +fields->A_2*((dx_th_3*h_02+dx_th_2*h_01+dx_th_1*h_00)*fields->A_2 -(dy_th_3*h_02+dy_th_2*h_01+dy_th_1*h_00)*fields->A_1)) /30;
                                      drdth2 -= -(dx_th_3*G_12+dx_th_2*(G_11+G_01) +fields->A_3*((dx_th_3*h_12+dx_th_2*(h_11+h_01))*fields->A_3 -(dz_th_3*h_12+dz_th_2*(h_11+h_01))*fields->A_1) +fields->A_2*((dx_th_3*h_12+dx_th_2*(h_11+h_01))*fields->A_2 -(dy_th_3*h_12+dy_th_2*(h_11+h_01))*fields->A_1)) /30;
                                      drdth3 -= -(dx_th_3*(G_22+G_12+G_02)+(h_22+h_12+h_02) *(fields->A_3*(dx_th_3*fields->A_3-dz_th_3*fields->A_1) +fields->A_2*(dx_th_3*fields->A_2-dy_th_3*fields->A_1))) /30;
                                      break;
                             case 1:  drdth1 -= (-dz_th_3*G_02+dy_th_3*G_02-dz_th_2*G_01+dy_th_2*G_01-dz_th_1*G_00 +dy_th_1*G_00 +fields->A_3*((dy_th_3*h_02+dy_th_2*h_01+dy_th_1*h_00)*fields->A_3 -(dz_th_3*h_02-dy_th_3*h_02+dz_th_2*h_01-dy_th_2*h_01 +dz_th_1*h_00-dy_th_1*h_00) *fields->A_2+(dx_th_3*h_02+dx_th_2*h_01+dx_th_1*h_00)*fields->A_1) -(dz_th_3*h_02+dz_th_2*h_01+dz_th_1*h_00)*fields->A_2*fields->A_2 -fields->A_1*((dx_th_3*h_02+dx_th_2*h_01+dx_th_1*h_00)*fields->A_2 +(dz_th_3*h_02-dy_th_3*h_02+dz_th_2*h_01-dy_th_2*h_01 +dz_th_1*h_00-dy_th_1*h_00) *fields->A_1)) /120;
                                      drdth2 -= (-dz_th_3*G_12+dy_th_3*G_12-dz_th_2*(G_11+G_01)+dy_th_2*(G_11+G_01) +fields->A_3*((dy_th_3*h_12+dy_th_2*(h_11+h_01))*fields->A_3 -(dz_th_3*h_12-dy_th_3*h_12+dz_th_2*(h_11+h_01) -dy_th_2*(h_11+h_01)) *fields->A_2+(dx_th_3*h_12+dx_th_2*(h_11+h_01))*fields->A_1) -(dz_th_3*h_12+dz_th_2*(h_11+h_01))*fields->A_2*fields->A_2 -fields->A_1*((dx_th_3*h_12+dx_th_2*(h_11+h_01))*fields->A_2 +(dz_th_3*h_12-dy_th_3*h_12+dz_th_2*(h_11+h_01) -dy_th_2*(h_11+h_01)) *fields->A_1)) /120;
                                      drdth3 -= ((dy_th_3-dz_th_3)*(G_22+G_12+G_02+(h_22+h_12+h_02)*fields->A_1*fields->A_1) +(h_22+h_12+h_02)*(fields->A_3*(dy_th_3*fields->A_3+(dy_th_3-dz_th_3)*fields->A_2+dx_th_3*fields->A_1) -fields->A_2*(dz_th_3*fields->A_2+dx_th_3*fields->A_1))) /120;
                                      break;
                         }
                         break;
                case 0:  switch(dk)
                         {
                             case -1: drdth1 -= -(dx_th_3*G_02+dx_th_2*G_01+dx_th_1*G_00 +fields->A_3*((dx_th_3*h_02+dx_th_2*h_01+dx_th_1*h_00)*fields->A_3 -(dz_th_3*h_02+dz_th_2*h_01+dz_th_1*h_00)*fields->A_1) +fields->A_2*((dx_th_3*h_02+dx_th_2*h_01+dx_th_1*h_00)*fields->A_2 -(dy_th_3*h_02+dy_th_2*h_01+dy_th_1*h_00)*fields->A_1)) /30;
                                      drdth2 -= -(dx_th_3*G_12+dx_th_2*(G_11+G_01) +fields->A_3*((dx_th_3*h_12+dx_th_2*(h_11+h_01))*fields->A_3 -(dz_th_3*h_12+dz_th_2*(h_11+h_01))*fields->A_1) +fields->A_2*((dx_th_3*h_12+dx_th_2*(h_11+h_01))*fields->A_2 -(dy_th_3*h_12+dy_th_2*(h_11+h_01))*fields->A_1)) /30;
                                      drdth3 -= -(dx_th_3*(G_22+G_12+G_02)+(h_22+h_12+h_02) *(fields->A_3*(dx_th_3*fields->A_3-dz_th_3*fields->A_1) +fields->A_2*(dx_th_3*fields->A_2-dy_th_3*fields->A_1))) /30;
                                      break;
                             case 0:  drdth1 -= (dx_th_3*G_02+dx_th_2*G_01+dx_th_1*G_00 +fields->A_3*((dx_th_3*h_02+dx_th_2*h_01+dx_th_1*h_00)*fields->A_3 -(dz_th_3*h_02+dz_th_2*h_01+dz_th_1*h_00)*fields->A_1) +fields->A_2*((dx_th_3*h_02+dx_th_2*h_01+dx_th_1*h_00)*fields->A_2 -(dy_th_3*h_02+dy_th_2*h_01+dy_th_1*h_00)*fields->A_1)) /20;
                                      drdth2 -= (dx_th_3*G_12+dx_th_2*(G_11+G_01) +fields->A_3*((dx_th_3*h_12+dx_th_2*(h_11+h_01))*fields->A_3 -(dz_th_3*h_12+dz_th_2*(h_11+h_01))*fields->A_1) +fields->A_2*((dx_th_3*h_12+dx_th_2*(h_11+h_01))*fields->A_2 -(dy_th_3*h_12+dy_th_2*(h_11+h_01))*fields->A_1)) /20;
                                      drdth3 -= (dx_th_3*(G_22+G_12+G_02)+(h_22+h_12+h_02) *(fields->A_3*(dx_th_3*fields->A_3-dz_th_3*fields->A_1) +fields->A_2*(dx_th_3*fields->A_2-dy_th_3*fields->A_1))) /20;
                                      break;
                             case 1:  drdth1 -= -(dx_th_3*G_02+dx_th_2*G_01+dx_th_1*G_00 +fields->A_3*((dx_th_3*h_02+dx_th_2*h_01+dx_th_1*h_00)*fields->A_3 -(dz_th_3*h_02+dz_th_2*h_01+dz_th_1*h_00)*fields->A_1) +fields->A_2*((dx_th_3*h_02+dx_th_2*h_01+dx_th_1*h_00)*fields->A_2 -(dy_th_3*h_02+dy_th_2*h_01+dy_th_1*h_00)*fields->A_1)) /30;
                                      drdth2 -= -(dx_th_3*G_12+dx_th_2*(G_11+G_01) +fields->A_3*((dx_th_3*h_12+dx_th_2*(h_11+h_01))*fields->A_3 -(dz_th_3*h_12+dz_th_2*(h_11+h_01))*fields->A_1) +fields->A_2*((dx_th_3*h_12+dx_th_2*(h_11+h_01))*fields->A_2 -(dy_th_3*h_12+dy_th_2*(h_11+h_01))*fields->A_1)) /30;
                                      drdth3 -= -(dx_th_3*(G_22+G_12+G_02)+(h_22+h_12+h_02) *(fields->A_3*(dx_th_3*fields->A_3-dz_th_3*fields->A_1) +fields->A_2*(dx_th_3*fields->A_2-dy_th_3*fields->A_1))) /30;
                                      break;
                         }
                         break;
                case 1:  switch(dk)
                         {
                             case -1: drdth1 -= -(-dz_th_3*G_02+dy_th_3*G_02-dz_th_2*G_01+dy_th_2*G_01-dz_th_1*G_00 +dy_th_1*G_00 +fields->A_3*((dy_th_3*h_02+dy_th_2*h_01+dy_th_1*h_00)*fields->A_3 -(dz_th_3*h_02-dy_th_3*h_02+dz_th_2*h_01-dy_th_2*h_01 +dz_th_1*h_00-dy_th_1*h_00) *fields->A_2+(dx_th_3*h_02+dx_th_2*h_01+dx_th_1*h_00)*fields->A_1) -(dz_th_3*h_02+dz_th_2*h_01+dz_th_1*h_00)*fields->A_2*fields->A_2 -fields->A_1*((dx_th_3*h_02+dx_th_2*h_01+dx_th_1*h_00)*fields->A_2 +(dz_th_3*h_02-dy_th_3*h_02+dz_th_2*h_01-dy_th_2*h_01 +dz_th_1*h_00-dy_th_1*h_00) *fields->A_1)) /120;
                                      drdth2 -= -(-dz_th_3*G_12+dy_th_3*G_12-dz_th_2*(G_11+G_01)+dy_th_2*(G_11+G_01) +fields->A_3*((dy_th_3*h_12+dy_th_2*(h_11+h_01))*fields->A_3 -(dz_th_3*h_12-dy_th_3*h_12+dz_th_2*(h_11+h_01) -dy_th_2*(h_11+h_01)) *fields->A_2+(dx_th_3*h_12+dx_th_2*(h_11+h_01))*fields->A_1) -(dz_th_3*h_12+dz_th_2*(h_11+h_01))*fields->A_2*fields->A_2 -fields->A_1*((dx_th_3*h_12+dx_th_2*(h_11+h_01))*fields->A_2 +(dz_th_3*h_12-dy_th_3*h_12+dz_th_2*(h_11+h_01) -dy_th_2*(h_11+h_01)) *fields->A_1)) /120;
                                      drdth3 -= -((dy_th_3-dz_th_3)*(G_22+G_12+G_02+(h_22+h_12+h_02)*fields->A_1*fields->A_1) +(h_22+h_12+h_02)*(fields->A_3*(dy_th_3*fields->A_3+(dy_th_3-dz_th_3)*fields->A_2+dx_th_3*fields->A_1) -fields->A_2*(dz_th_3*fields->A_2+dx_th_3*fields->A_1))) /120;
                                      break;
                             case 0:  drdth1 -= -(dx_th_3*G_02+dx_th_2*G_01+dx_th_1*G_00 +fields->A_3*((dx_th_3*h_02+dx_th_2*h_01+dx_th_1*h_00)*fields->A_3 -(dz_th_3*h_02+dz_th_2*h_01+dz_th_1*h_00)*fields->A_1) +fields->A_2*((dx_th_3*h_02+dx_th_2*h_01+dx_th_1*h_00)*fields->A_2 -(dy_th_3*h_02+dy_th_2*h_01+dy_th_1*h_00)*fields->A_1)) /30;
                                      drdth2 -= -(dx_th_3*G_12+dx_th_2*(G_11+G_01) +fields->A_3*((dx_th_3*h_12+dx_th_2*(h_11+h_01))*fields->A_3 -(dz_th_3*h_12+dz_th_2*(h_11+h_01))*fields->A_1) +fields->A_2*((dx_th_3*h_12+dx_th_2*(h_11+h_01))*fields->A_2 -(dy_th_3*h_12+dy_th_2*(h_11+h_01))*fields->A_1)) /30;
                                      drdth3 -= -(dx_th_3*(G_22+G_12+G_02)+(h_22+h_12+h_02) *(fields->A_3*(dx_th_3*fields->A_3-dz_th_3*fields->A_1) +fields->A_2*(dx_th_3*fields->A_2-dy_th_3*fields->A_1))) /30;
                                      break;
                             case 1:  drdth1 -= -(dz_th_3*G_02+dy_th_3*G_02+dz_th_2*G_01+dy_th_2*G_01+dz_th_1*G_00 +dy_th_1*G_00 +fields->A_3*((dy_th_3*h_02+dy_th_2*h_01+dy_th_1*h_00)*fields->A_3 -(dz_th_3*h_02+dy_th_3*h_02+dz_th_2*h_01+dy_th_2*h_01 +dz_th_1*h_00+dy_th_1*h_00) *fields->A_2-(dx_th_3*h_02+dx_th_2*h_01+dx_th_1*h_00)*fields->A_1) +(dz_th_3*h_02+dz_th_2*h_01+dz_th_1*h_00)*fields->A_2*fields->A_2 +fields->A_1*((dz_th_3*h_02+dy_th_3*h_02+dz_th_2*h_01+dy_th_2*h_01 +dz_th_1*h_00+dy_th_1*h_00) *fields->A_1 -(dx_th_3*h_02+dx_th_2*h_01+dx_th_1*h_00)*fields->A_2)) /120;
                                      drdth2 -= -(dz_th_3*G_12+dy_th_3*G_12+dz_th_2*(G_11+G_01)+dy_th_2*(G_11+G_01) +fields->A_3*((dy_th_3*h_12+dy_th_2*(h_11+h_01))*fields->A_3 -(dz_th_3*h_12+dy_th_3*h_12+dz_th_2*(h_11+h_01) +dy_th_2*(h_11+h_01)) *fields->A_2-(dx_th_3*h_12+dx_th_2*(h_11+h_01))*fields->A_1) +(dz_th_3*h_12+dz_th_2*(h_11+h_01))*fields->A_2*fields->A_2 +fields->A_1*((dz_th_3*h_12+dy_th_3*h_12+dz_th_2*(h_11+h_01) +dy_th_2*(h_11+h_01)) *fields->A_1 -(dx_th_3*h_12+dx_th_2*(h_11+h_01))*fields->A_2)) /120;
                                      drdth3 -= -((dz_th_3+dy_th_3)*(G_22+G_12+G_02+(h_22+h_12+h_02)*fields->A_1*fields->A_1) +(h_22+h_12+h_02)*(fields->A_3*(dy_th_3*fields->A_3-((dz_th_3+dy_th_3)*fields->A_2+dx_th_3*fields->A_1)) +fields->A_2*(dz_th_3*fields->A_2-dx_th_3*fields->A_1))) /120;
                                      break;
                         }
                         break;
              }
              break;
   }
}
