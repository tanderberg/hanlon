/* fields_aux_c.cpp : Auxiliary functionality for scalars in cartesian form.
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
const int MAX_NSIDE         = 64;     /* Max number of lattice sites per side.
                                         Power of 2 for faster array indexing (maybe). */
#endif
const int DEFAULT_NSIDE     = 64;     /* Default number lattice sites per side. */
const int NCOMPS            = 7;      /* Number of field components. */

const int OUTPUT_PRECISION  = 15;

const double DEFAULT_LAMBDA = 1000.0; /* Large constant for fake harmonic potential */

#include "math_aux.cpp"

struct TFields
{
   double ph_0, ph_1, ph_2, ph_3,
          A_1, A_2, A_3;

   void clear();
   void add_fields(TFields *otherFields);
   void sub_fields(TFields *otherFields);
   void mul_fields(double factor);
   double scalarproduct_ph(TFields *otherFields);
   double norm_ph();
   void normalize_ph();
};

void TFields::clear()
{
   ph_0 = 0;
   ph_1 = 0;
   ph_2 = 0;
   ph_3 = 0;

   A_1 = 0;
   A_2 = 0;
   A_3 = 0;
}

void TFields::add_fields(TFields *otherFields)
{
   ph_0 += otherFields->ph_0;
   ph_1 += otherFields->ph_1;
   ph_2 += otherFields->ph_2;
   ph_3 += otherFields->ph_3;

   A_1 += otherFields->A_1;
   A_2 += otherFields->A_2;
   A_3 += otherFields->A_3;
}

void TFields::sub_fields(TFields *otherFields)
{
   ph_0 -= otherFields->ph_0;
   ph_1 -= otherFields->ph_1;
   ph_2 -= otherFields->ph_2;
   ph_3 -= otherFields->ph_3;

   A_1 -= otherFields->A_1;
   A_2 -= otherFields->A_2;
   A_3 -= otherFields->A_3;
}

void TFields::mul_fields(double factor)
{
   ph_0 *= factor;
   ph_1 *= factor;
   ph_2 *= factor;
   ph_3 *= factor;

   A_1 *= factor;
   A_2 *= factor;
   A_3 *= factor;
}

double TFields::scalarproduct_ph(TFields *otherFields)
{
   return ph_0 * otherFields->ph_0 + 
          ph_1 * otherFields->ph_1 + 
          ph_2 * otherFields->ph_2 + 
          ph_3 * otherFields->ph_3;
}

double TFields::norm_ph()
{
   return scalarproduct_ph(this);
}

void TFields::normalize_ph()
{
   double norm = sqrt(norm_ph());
   
   ph_0 /= norm;
   ph_1 /= norm;
   ph_2 /= norm;
   ph_3 /= norm;
}

struct TAux
{
   double dx_ph_0, dx_ph_1, dx_ph_2, dx_ph_3, // IMPORTANT! Derivatives stored here are not normalized by lattice spacing
          dy_ph_0, dy_ph_1, dy_ph_2, dy_ph_3,
          dz_ph_0, dz_ph_1, dz_ph_2, dz_ph_3,

          dx_A_2, dx_A_3,
          dy_A_1, dy_A_3,
          dz_A_1, dz_A_2;
   
   TAux();
   
   void clear();
   
   void compute_ph_derivatives_core(int i, int j, int k, int iM1, int jM1, int kM1, int iM2, int jM2, int kM2, 
                                    int iP1, int jP1, int kP1, int iP2, int jP2, int kP2,
                                    TFields fields[MAX_NSIDE][MAX_NSIDE][MAX_NSIDE]);
   void compute_ph_derivatives(int i, int j, int k, int nSideM1, TFields fields[MAX_NSIDE][MAX_NSIDE][MAX_NSIDE]);
   void compute_ph_derivatives(int i, int j, int k, int nSideM1, TFields fields[MAX_NSIDE][MAX_NSIDE][MAX_NSIDE],
                               double &out_dx_ph_0, double &out_dx_ph_1, double &out_dx_ph_2, double &out_dx_ph_3,  
                               double &out_dy_ph_0, double &out_dy_ph_1, double &out_dy_ph_2, double &out_dy_ph_3, 
                               double &out_dz_ph_0, double &out_dz_ph_1, double &out_dz_ph_2, double &out_dz_ph_3);
   void get_ph_derivatives(double &out_dx_ph_0, double &out_dx_ph_1, double &out_dx_ph_2, double &out_dx_ph_3,  
                           double &out_dy_ph_0, double &out_dy_ph_1, double &out_dy_ph_2, double &out_dy_ph_3, 
                           double &out_dz_ph_0, double &out_dz_ph_1, double &out_dz_ph_2, double &out_dz_ph_3);
   void compute_A_derivatives_core(int i, int j, int k, 
                                   int iM1, int jM1, int kM1, int iM2, int jM2, int kM2, 
                                   int iP1, int jP1, int kP1, int iP2, int jP2, int kP2, 
                                   TFields fields[MAX_NSIDE][MAX_NSIDE][MAX_NSIDE]);
   void compute_A_derivatives(int i, int j, int k, int nSideM1, TFields fields[MAX_NSIDE][MAX_NSIDE][MAX_NSIDE]);
   void compute_A_derivatives(int i, int j, int k, int nSideM1, TFields fields[MAX_NSIDE][MAX_NSIDE][MAX_NSIDE],
                              double &out_dx_A_2, double &out_dx_A_3, double &out_dy_A_1, 
                              double &out_dy_A_3, double &out_dz_A_1, double &out_dz_A_2);
   void get_A_derivatives(double &out_dx_A_2, double &out_dx_A_3, double &out_dy_A_1, 
                          double &out_dy_A_3, double &out_dz_A_1, double &out_dz_A_2);
                          
   double compute_divA(int i, int j, int k, int nSideM1, TFields fields[MAX_NSIDE][MAX_NSIDE][MAX_NSIDE]);

   void compute_derivatives(int i, int j, int k, int nSideM1, TFields fields[MAX_NSIDE][MAX_NSIDE][MAX_NSIDE]);
   void compute_derivatives(int i, int j, int k, int nSideM1, TFields fields[MAX_NSIDE][MAX_NSIDE][MAX_NSIDE],
                            double &out_dx_ph_0, double &out_dx_ph_1, double &out_dx_ph_2, double &out_dx_ph_3, 
                            double &out_dx_A_2, double &out_dx_A_3, 
                            double &out_dy_ph_0, double &out_dy_ph_1, double &out_dy_ph_2, double &out_dy_ph_3, 
                            double &out_dy_A_1, double &out_dy_A_3, 
                            double &out_dz_ph_0, double &out_dz_ph_1, double &out_dz_ph_2, double &out_dz_ph_3, 
                            double &out_dz_A_1, double &out_dz_A_2);
   void get_derivatives(double &out_dx_ph_0, double &out_dx_ph_1, double &out_dx_ph_2, double &out_dx_ph_3, 
                        double &out_dx_A_2, double &out_dx_A_3, 
                        double &out_dy_ph_0, double &out_dy_ph_1, double &out_dy_ph_2, double &out_dy_ph_3, 
                        double &out_dy_A_1, double &out_dy_A_3, 
                        double &out_dz_ph_0, double &out_dz_ph_1, double &out_dz_ph_2, double &out_dz_ph_3, 
                        double &out_dz_A_1, double &out_dz_A_2);
};

TAux::TAux()
{
   clear();
}

void TAux::clear()
{
   dx_ph_0 = 0;
   dx_ph_1 = 0;
   dx_ph_2 = 0;
   dx_ph_3 = 0;
   dy_ph_0 = 0;
   dy_ph_1 = 0;
   dy_ph_2 = 0;
   dy_ph_3 = 0;
   dz_ph_0 = 0;
   dz_ph_1 = 0;
   dz_ph_2 = 0;
   dz_ph_3 = 0;

   dx_A_2 = 0;
   dx_A_3 = 0;
   dy_A_1 = 0;
   dy_A_3 = 0;
   dz_A_1 = 0;
   dz_A_2 = 0;
}

void TAux::compute_ph_derivatives_core(int i, int j, int k, int iM1, int jM1, int kM1, int iM2, int jM2, int kM2, 
                                       int iP1, int jP1, int kP1, int iP2, int jP2, int kP2,
                                       TFields fields[MAX_NSIDE][MAX_NSIDE][MAX_NSIDE])
/* IMPORTANT: Does NOT divide by lattice spacing, i.e. if lattice spacing != 1, caller must divide results by it. */
{
   dx_ph_0 = 
                    (  fields[iM1][jM1][kM2].ph_0 + fields[iM1][jM1][kP2].ph_0 + fields[iM1][jP1][kM2].ph_0 + fields[iM1][jP1][kP2].ph_0
                     + fields[iM1][jM2][kM1].ph_0 + fields[iM1][jM2][kP1].ph_0 + fields[iM1][jP2][kM1].ph_0 + fields[iM1][jP2][kP1].ph_0
                     - fields[iP1][jM1][kM2].ph_0 - fields[iP1][jM1][kP2].ph_0 - fields[iP1][jP1][kM2].ph_0 - fields[iP1][jP1][kP2].ph_0
                     - fields[iP1][jM2][kM1].ph_0 - fields[iP1][jM2][kP1].ph_0 - fields[iP1][jP2][kM1].ph_0 - fields[iP1][jP2][kP1].ph_0)/120
                  + (  fields[iM2][j][kM1].ph_0 + fields[iM2][j][kP1].ph_0 + fields[iM2][jM1][k].ph_0 + fields[iM2][jP1][k].ph_0
                     - fields[iP2][j][kM1].ph_0 - fields[iP2][j][kP1].ph_0 - fields[iP2][jM1][k].ph_0 - fields[iP2][jP1][k].ph_0
                     - fields[iM1][jM1][kM1].ph_0 - fields[iM1][jM1][kP1].ph_0 - fields[iM1][jP1][kM1].ph_0 - fields[iM1][jP1][kP1].ph_0
                     + fields[iP1][jM1][kM1].ph_0 + fields[iP1][jM1][kP1].ph_0 + fields[iP1][jP1][kM1].ph_0 + fields[iP1][jP1][kP1].ph_0)/30
                  + (- fields[iM2][j][k].ph_0 + fields[iP2][j][k].ph_0)/20
                  - (  fields[iM1][j][kM1].ph_0 + fields[iM1][j][kP1].ph_0 + fields[iM1][jM1][k].ph_0 + fields[iM1][jP1][k].ph_0
                     - fields[iP1][j][kM1].ph_0 - fields[iP1][j][kP1].ph_0 - fields[iP1][jM1][k].ph_0 - fields[iP1][jP1][k].ph_0)/12
                  + (- fields[iM1][j][k].ph_0 + fields[iP1][j][k].ph_0)*4.0/15;
   dx_ph_1 = 
                    (  fields[iM1][jM1][kM2].ph_1 + fields[iM1][jM1][kP2].ph_1 + fields[iM1][jP1][kM2].ph_1 + fields[iM1][jP1][kP2].ph_1
                     + fields[iM1][jM2][kM1].ph_1 + fields[iM1][jM2][kP1].ph_1 + fields[iM1][jP2][kM1].ph_1 + fields[iM1][jP2][kP1].ph_1
                     - fields[iP1][jM1][kM2].ph_1 - fields[iP1][jM1][kP2].ph_1 - fields[iP1][jP1][kM2].ph_1 - fields[iP1][jP1][kP2].ph_1
                     - fields[iP1][jM2][kM1].ph_1 - fields[iP1][jM2][kP1].ph_1 - fields[iP1][jP2][kM1].ph_1 - fields[iP1][jP2][kP1].ph_1)/120
                  + (  fields[iM2][j][kM1].ph_1 + fields[iM2][j][kP1].ph_1 + fields[iM2][jM1][k].ph_1 + fields[iM2][jP1][k].ph_1
                     - fields[iP2][j][kM1].ph_1 - fields[iP2][j][kP1].ph_1 - fields[iP2][jM1][k].ph_1 - fields[iP2][jP1][k].ph_1
                     - fields[iM1][jM1][kM1].ph_1 - fields[iM1][jM1][kP1].ph_1 - fields[iM1][jP1][kM1].ph_1 - fields[iM1][jP1][kP1].ph_1
                     + fields[iP1][jM1][kM1].ph_1 + fields[iP1][jM1][kP1].ph_1 + fields[iP1][jP1][kM1].ph_1 + fields[iP1][jP1][kP1].ph_1)/30
                  + (- fields[iM2][j][k].ph_1 + fields[iP2][j][k].ph_1)/20
                  - (  fields[iM1][j][kM1].ph_1 + fields[iM1][j][kP1].ph_1 + fields[iM1][jM1][k].ph_1 + fields[iM1][jP1][k].ph_1
                     - fields[iP1][j][kM1].ph_1 - fields[iP1][j][kP1].ph_1 - fields[iP1][jM1][k].ph_1 - fields[iP1][jP1][k].ph_1)/12
                  + (- fields[iM1][j][k].ph_1 + fields[iP1][j][k].ph_1)*4.0/15;
   dx_ph_2 = 
                    (  fields[iM1][jM1][kM2].ph_2 + fields[iM1][jM1][kP2].ph_2 + fields[iM1][jP1][kM2].ph_2 + fields[iM1][jP1][kP2].ph_2
                     + fields[iM1][jM2][kM1].ph_2 + fields[iM1][jM2][kP1].ph_2 + fields[iM1][jP2][kM1].ph_2 + fields[iM1][jP2][kP1].ph_2
                     - fields[iP1][jM1][kM2].ph_2 - fields[iP1][jM1][kP2].ph_2 - fields[iP1][jP1][kM2].ph_2 - fields[iP1][jP1][kP2].ph_2
                     - fields[iP1][jM2][kM1].ph_2 - fields[iP1][jM2][kP1].ph_2 - fields[iP1][jP2][kM1].ph_2 - fields[iP1][jP2][kP1].ph_2)/120
                  + (  fields[iM2][j][kM1].ph_2 + fields[iM2][j][kP1].ph_2 + fields[iM2][jM1][k].ph_2 + fields[iM2][jP1][k].ph_2
                     - fields[iP2][j][kM1].ph_2 - fields[iP2][j][kP1].ph_2 - fields[iP2][jM1][k].ph_2 - fields[iP2][jP1][k].ph_2
                     - fields[iM1][jM1][kM1].ph_2 - fields[iM1][jM1][kP1].ph_2 - fields[iM1][jP1][kM1].ph_2 - fields[iM1][jP1][kP1].ph_2
                     + fields[iP1][jM1][kM1].ph_2 + fields[iP1][jM1][kP1].ph_2 + fields[iP1][jP1][kM1].ph_2 + fields[iP1][jP1][kP1].ph_2)/30
                  + (- fields[iM2][j][k].ph_2 + fields[iP2][j][k].ph_2)/20
                  - (  fields[iM1][j][kM1].ph_2 + fields[iM1][j][kP1].ph_2 + fields[iM1][jM1][k].ph_2 + fields[iM1][jP1][k].ph_2
                     - fields[iP1][j][kM1].ph_2 - fields[iP1][j][kP1].ph_2 - fields[iP1][jM1][k].ph_2 - fields[iP1][jP1][k].ph_2)/12
                  + (- fields[iM1][j][k].ph_2 + fields[iP1][j][k].ph_2)*4.0/15;
   dx_ph_3 = 
                    (  fields[iM1][jM1][kM2].ph_3 + fields[iM1][jM1][kP2].ph_3 + fields[iM1][jP1][kM2].ph_3 + fields[iM1][jP1][kP2].ph_3
                     + fields[iM1][jM2][kM1].ph_3 + fields[iM1][jM2][kP1].ph_3 + fields[iM1][jP2][kM1].ph_3 + fields[iM1][jP2][kP1].ph_3
                     - fields[iP1][jM1][kM2].ph_3 - fields[iP1][jM1][kP2].ph_3 - fields[iP1][jP1][kM2].ph_3 - fields[iP1][jP1][kP2].ph_3
                     - fields[iP1][jM2][kM1].ph_3 - fields[iP1][jM2][kP1].ph_3 - fields[iP1][jP2][kM1].ph_3 - fields[iP1][jP2][kP1].ph_3)/120
                  + (  fields[iM2][j][kM1].ph_3 + fields[iM2][j][kP1].ph_3 + fields[iM2][jM1][k].ph_3 + fields[iM2][jP1][k].ph_3
                     - fields[iP2][j][kM1].ph_3 - fields[iP2][j][kP1].ph_3 - fields[iP2][jM1][k].ph_3 - fields[iP2][jP1][k].ph_3
                     - fields[iM1][jM1][kM1].ph_3 - fields[iM1][jM1][kP1].ph_3 - fields[iM1][jP1][kM1].ph_3 - fields[iM1][jP1][kP1].ph_3
                     + fields[iP1][jM1][kM1].ph_3 + fields[iP1][jM1][kP1].ph_3 + fields[iP1][jP1][kM1].ph_3 + fields[iP1][jP1][kP1].ph_3)/30
                  + (- fields[iM2][j][k].ph_3 + fields[iP2][j][k].ph_3)/20
                  - (  fields[iM1][j][kM1].ph_3 + fields[iM1][j][kP1].ph_3 + fields[iM1][jM1][k].ph_3 + fields[iM1][jP1][k].ph_3
                     - fields[iP1][j][kM1].ph_3 - fields[iP1][j][kP1].ph_3 - fields[iP1][jM1][k].ph_3 - fields[iP1][jP1][k].ph_3)/12
                  + (- fields[iM1][j][k].ph_3 + fields[iP1][j][k].ph_3)*4.0/15;
   dy_ph_0 = 
                    (  fields[iM1][jM1][kM2].ph_0 + fields[iM1][jM1][kP2].ph_0 - fields[iM1][jP1][kM2].ph_0 - fields[iM1][jP1][kP2].ph_0
                     + fields[iP1][jM1][kM2].ph_0 + fields[iP1][jM1][kP2].ph_0 - fields[iP1][jP1][kM2].ph_0 - fields[iP1][jP1][kP2].ph_0
                     + fields[iM2][jM1][kM1].ph_0 + fields[iM2][jM1][kP1].ph_0 - fields[iM2][jP1][kM1].ph_0 - fields[iM2][jP1][kP1].ph_0
                     + fields[iP2][jM1][kM1].ph_0 + fields[iP2][jM1][kP1].ph_0 - fields[iP2][jP1][kM1].ph_0 - fields[iP2][jP1][kP1].ph_0)/120
                  + (  fields[i][jM2][kM1].ph_0 + fields[i][jM2][kP1].ph_0 - fields[i][jP2][kM1].ph_0 - fields[i][jP2][kP1].ph_0
                     + fields[iM1][jM2][k].ph_0 - fields[iM1][jP2][k].ph_0 + fields[iP1][jM2][k].ph_0 - fields[iP1][jP2][k].ph_0
                     - fields[iM1][jM1][kM1].ph_0 - fields[iM1][jM1][kP1].ph_0 + fields[iM1][jP1][kM1].ph_0 + fields[iM1][jP1][kP1].ph_0
                     - fields[iP1][jM1][kM1].ph_0 - fields[iP1][jM1][kP1].ph_0 + fields[iP1][jP1][kM1].ph_0 + fields[iP1][jP1][kP1].ph_0)/30
                  + (- fields[i][jM2][k].ph_0 + fields[i][jP2][k].ph_0)/20
                  - (  fields[i][jM1][kM1].ph_0 + fields[i][jM1][kP1].ph_0 - fields[i][jP1][kM1].ph_0 - fields[i][jP1][kP1].ph_0
                     + fields[iM1][jM1][k].ph_0 - fields[iM1][jP1][k].ph_0 + fields[iP1][jM1][k].ph_0 - fields[iP1][jP1][k].ph_0)/12
                  + (- fields[i][jM1][k].ph_0 + fields[i][jP1][k].ph_0)*4.0/15;
   dy_ph_1 = 
                    (  fields[iM1][jM1][kM2].ph_1 + fields[iM1][jM1][kP2].ph_1 - fields[iM1][jP1][kM2].ph_1 - fields[iM1][jP1][kP2].ph_1
                     + fields[iP1][jM1][kM2].ph_1 + fields[iP1][jM1][kP2].ph_1 - fields[iP1][jP1][kM2].ph_1 - fields[iP1][jP1][kP2].ph_1
                     + fields[iM2][jM1][kM1].ph_1 + fields[iM2][jM1][kP1].ph_1 - fields[iM2][jP1][kM1].ph_1 - fields[iM2][jP1][kP1].ph_1
                     + fields[iP2][jM1][kM1].ph_1 + fields[iP2][jM1][kP1].ph_1 - fields[iP2][jP1][kM1].ph_1 - fields[iP2][jP1][kP1].ph_1)/120
                  + (  fields[i][jM2][kM1].ph_1 + fields[i][jM2][kP1].ph_1 - fields[i][jP2][kM1].ph_1 - fields[i][jP2][kP1].ph_1
                     + fields[iM1][jM2][k].ph_1 - fields[iM1][jP2][k].ph_1 + fields[iP1][jM2][k].ph_1 - fields[iP1][jP2][k].ph_1
                     - fields[iM1][jM1][kM1].ph_1 - fields[iM1][jM1][kP1].ph_1 + fields[iM1][jP1][kM1].ph_1 + fields[iM1][jP1][kP1].ph_1
                     - fields[iP1][jM1][kM1].ph_1 - fields[iP1][jM1][kP1].ph_1 + fields[iP1][jP1][kM1].ph_1 + fields[iP1][jP1][kP1].ph_1)/30
                  + (- fields[i][jM2][k].ph_1 + fields[i][jP2][k].ph_1)/20
                  - (  fields[i][jM1][kM1].ph_1 + fields[i][jM1][kP1].ph_1 - fields[i][jP1][kM1].ph_1 - fields[i][jP1][kP1].ph_1
                     + fields[iM1][jM1][k].ph_1 - fields[iM1][jP1][k].ph_1 + fields[iP1][jM1][k].ph_1 - fields[iP1][jP1][k].ph_1)/12
                  + (- fields[i][jM1][k].ph_1 + fields[i][jP1][k].ph_1)*4.0/15;
   dy_ph_2 = 
                    (  fields[iM1][jM1][kM2].ph_2 + fields[iM1][jM1][kP2].ph_2 - fields[iM1][jP1][kM2].ph_2 - fields[iM1][jP1][kP2].ph_2
                     + fields[iP1][jM1][kM2].ph_2 + fields[iP1][jM1][kP2].ph_2 - fields[iP1][jP1][kM2].ph_2 - fields[iP1][jP1][kP2].ph_2
                     + fields[iM2][jM1][kM1].ph_2 + fields[iM2][jM1][kP1].ph_2 - fields[iM2][jP1][kM1].ph_2 - fields[iM2][jP1][kP1].ph_2
                     + fields[iP2][jM1][kM1].ph_2 + fields[iP2][jM1][kP1].ph_2 - fields[iP2][jP1][kM1].ph_2 - fields[iP2][jP1][kP1].ph_2)/120
                  + (  fields[i][jM2][kM1].ph_2 + fields[i][jM2][kP1].ph_2 - fields[i][jP2][kM1].ph_2 - fields[i][jP2][kP1].ph_2
                     + fields[iM1][jM2][k].ph_2 - fields[iM1][jP2][k].ph_2 + fields[iP1][jM2][k].ph_2 - fields[iP1][jP2][k].ph_2
                     - fields[iM1][jM1][kM1].ph_2 - fields[iM1][jM1][kP1].ph_2 + fields[iM1][jP1][kM1].ph_2 + fields[iM1][jP1][kP1].ph_2
                     - fields[iP1][jM1][kM1].ph_2 - fields[iP1][jM1][kP1].ph_2 + fields[iP1][jP1][kM1].ph_2 + fields[iP1][jP1][kP1].ph_2)/30
                  + (- fields[i][jM2][k].ph_2 + fields[i][jP2][k].ph_2)/20
                  + (- fields[i][jM1][kM1].ph_2 - fields[i][jM1][kP1].ph_2 + fields[i][jP1][kM1].ph_2 + fields[i][jP1][kP1].ph_2
                     - fields[iM1][jM1][k].ph_2 + fields[iM1][jP1][k].ph_2 - fields[iP1][jM1][k].ph_2 + fields[iP1][jP1][k].ph_2)/12
                  + (- fields[i][jM1][k].ph_2 + fields[i][jP1][k].ph_2)*4.0/15;
   dy_ph_3 = 
                    (  fields[iM1][jM1][kM2].ph_3 + fields[iM1][jM1][kP2].ph_3 - fields[iM1][jP1][kM2].ph_3 - fields[iM1][jP1][kP2].ph_3
                     + fields[iP1][jM1][kM2].ph_3 + fields[iP1][jM1][kP2].ph_3 - fields[iP1][jP1][kM2].ph_3 - fields[iP1][jP1][kP2].ph_3
                     + fields[iM2][jM1][kM1].ph_3 + fields[iM2][jM1][kP1].ph_3 - fields[iM2][jP1][kM1].ph_3 - fields[iM2][jP1][kP1].ph_3
                     + fields[iP2][jM1][kM1].ph_3 + fields[iP2][jM1][kP1].ph_3 - fields[iP2][jP1][kM1].ph_3 - fields[iP2][jP1][kP1].ph_3)/120
                  + (  fields[i][jM2][kM1].ph_3 + fields[i][jM2][kP1].ph_3 - fields[i][jP2][kM1].ph_3 - fields[i][jP2][kP1].ph_3
                     + fields[iM1][jM2][k].ph_3 - fields[iM1][jP2][k].ph_3 + fields[iP1][jM2][k].ph_3 - fields[iP1][jP2][k].ph_3
                     - fields[iM1][jM1][kM1].ph_3 - fields[iM1][jM1][kP1].ph_3 + fields[iM1][jP1][kM1].ph_3 + fields[iM1][jP1][kP1].ph_3
                     - fields[iP1][jM1][kM1].ph_3 - fields[iP1][jM1][kP1].ph_3 + fields[iP1][jP1][kM1].ph_3 + fields[iP1][jP1][kP1].ph_3)/30
                  - (  fields[i][jM2][k].ph_3 - fields[i][jP2][k].ph_3)/20
                  + (- fields[i][jM1][kM1].ph_3 - fields[i][jM1][kP1].ph_3 + fields[i][jP1][kM1].ph_3 + fields[i][jP1][kP1].ph_3
                     - fields[iM1][jM1][k].ph_3 + fields[iM1][jP1][k].ph_3 - fields[iP1][jM1][k].ph_3 + fields[iP1][jP1][k].ph_3)/12
                  - (  fields[i][jM1][k].ph_3 - fields[i][jP1][k].ph_3)*4.0/15;
   dz_ph_0 = 
                    (  fields[iM1][jM2][kM1].ph_0 - fields[iM1][jM2][kP1].ph_0 + fields[iM1][jP2][kM1].ph_0 - fields[iM1][jP2][kP1].ph_0
                     + fields[iP1][jM2][kM1].ph_0 - fields[iP1][jM2][kP1].ph_0 + fields[iP1][jP2][kM1].ph_0 - fields[iP1][jP2][kP1].ph_0
                     + fields[iM2][jM1][kM1].ph_0 - fields[iM2][jM1][kP1].ph_0 + fields[iM2][jP1][kM1].ph_0 - fields[iM2][jP1][kP1].ph_0
                     + fields[iP2][jM1][kM1].ph_0 - fields[iP2][jM1][kP1].ph_0 + fields[iP2][jP1][kM1].ph_0 - fields[iP2][jP1][kP1].ph_0)/120
                  + (  fields[i][jM1][kM2].ph_0 - fields[i][jM1][kP2].ph_0 + fields[i][jP1][kM2].ph_0 - fields[i][jP1][kP2].ph_0
                     + fields[iM1][j][kM2].ph_0 - fields[iM1][j][kP2].ph_0 + fields[iP1][j][kM2].ph_0 - fields[iP1][j][kP2].ph_0
                     - fields[iM1][jM1][kM1].ph_0 + fields[iM1][jM1][kP1].ph_0 - fields[iM1][jP1][kM1].ph_0 + fields[iM1][jP1][kP1].ph_0
                     - fields[iP1][jM1][kM1].ph_0 + fields[iP1][jM1][kP1].ph_0 - fields[iP1][jP1][kM1].ph_0 + fields[iP1][jP1][kP1].ph_0)/30
                  + (- fields[i][j][kM2].ph_0 + fields[i][j][kP2].ph_0)/20
                  - (  fields[i][jM1][kM1].ph_0 - fields[i][jM1][kP1].ph_0 + fields[i][jP1][kM1].ph_0 - fields[i][jP1][kP1].ph_0
                     + fields[iM1][j][kM1].ph_0 - fields[iM1][j][kP1].ph_0 + fields[iP1][j][kM1].ph_0 - fields[iP1][j][kP1].ph_0)/12
                  + (- fields[i][j][kM1].ph_0 + fields[i][j][kP1].ph_0)*4.0/15;
   dz_ph_1 = 
                    (  fields[iM1][jM2][kM1].ph_1 - fields[iM1][jM2][kP1].ph_1 + fields[iM1][jP2][kM1].ph_1 - fields[iM1][jP2][kP1].ph_1
                     + fields[iP1][jM2][kM1].ph_1 - fields[iP1][jM2][kP1].ph_1 + fields[iP1][jP2][kM1].ph_1 - fields[iP1][jP2][kP1].ph_1
                     + fields[iM2][jM1][kM1].ph_1 - fields[iM2][jM1][kP1].ph_1 + fields[iM2][jP1][kM1].ph_1 - fields[iM2][jP1][kP1].ph_1
                     + fields[iP2][jM1][kM1].ph_1 - fields[iP2][jM1][kP1].ph_1 + fields[iP2][jP1][kM1].ph_1 - fields[iP2][jP1][kP1].ph_1)/120
                  + (  fields[i][jM1][kM2].ph_1 - fields[i][jM1][kP2].ph_1 + fields[i][jP1][kM2].ph_1 - fields[i][jP1][kP2].ph_1
                     + fields[iM1][j][kM2].ph_1 - fields[iM1][j][kP2].ph_1 + fields[iP1][j][kM2].ph_1 - fields[iP1][j][kP2].ph_1
                     - fields[iM1][jM1][kM1].ph_1 + fields[iM1][jM1][kP1].ph_1 - fields[iM1][jP1][kM1].ph_1 + fields[iM1][jP1][kP1].ph_1
                     - fields[iP1][jM1][kM1].ph_1 + fields[iP1][jM1][kP1].ph_1 - fields[iP1][jP1][kM1].ph_1 + fields[iP1][jP1][kP1].ph_1)/30
                  + (- fields[i][j][kM2].ph_1 + fields[i][j][kP2].ph_1)/20
                  - (  fields[i][jM1][kM1].ph_1 - fields[i][jM1][kP1].ph_1 + fields[i][jP1][kM1].ph_1 - fields[i][jP1][kP1].ph_1
                     + fields[iM1][j][kM1].ph_1 - fields[iM1][j][kP1].ph_1 + fields[iP1][j][kM1].ph_1 - fields[iP1][j][kP1].ph_1)/12
                  + (- fields[i][j][kM1].ph_1 + fields[i][j][kP1].ph_1)*4.0/15;
   dz_ph_2 = 
                    (  fields[iM1][jM2][kM1].ph_2 - fields[iM1][jM2][kP1].ph_2 + fields[iM1][jP2][kM1].ph_2 - fields[iM1][jP2][kP1].ph_2
                     + fields[iP1][jM2][kM1].ph_2 - fields[iP1][jM2][kP1].ph_2 + fields[iP1][jP2][kM1].ph_2 - fields[iP1][jP2][kP1].ph_2
                     + fields[iM2][jM1][kM1].ph_2 - fields[iM2][jM1][kP1].ph_2 + fields[iM2][jP1][kM1].ph_2 - fields[iM2][jP1][kP1].ph_2
                     + fields[iP2][jM1][kM1].ph_2 - fields[iP2][jM1][kP1].ph_2 + fields[iP2][jP1][kM1].ph_2 - fields[iP2][jP1][kP1].ph_2)/120
                  + (  fields[i][jM1][kM2].ph_2 - fields[i][jM1][kP2].ph_2 + fields[i][jP1][kM2].ph_2 - fields[i][jP1][kP2].ph_2
                     + fields[iM1][j][kM2].ph_2 - fields[iM1][j][kP2].ph_2 + fields[iP1][j][kM2].ph_2 - fields[iP1][j][kP2].ph_2
                     - fields[iM1][jM1][kM1].ph_2 + fields[iM1][jM1][kP1].ph_2 - fields[iM1][jP1][kM1].ph_2 + fields[iM1][jP1][kP1].ph_2
                     - fields[iP1][jM1][kM1].ph_2 + fields[iP1][jM1][kP1].ph_2 - fields[iP1][jP1][kM1].ph_2 + fields[iP1][jP1][kP1].ph_2)/30
                  + (- fields[i][j][kM2].ph_2 + fields[i][j][kP2].ph_2)/20
                  - (  fields[i][jM1][kM1].ph_2 - fields[i][jM1][kP1].ph_2 + fields[i][jP1][kM1].ph_2 - fields[i][jP1][kP1].ph_2
                     + fields[iM1][j][kM1].ph_2 - fields[iM1][j][kP1].ph_2 + fields[iP1][j][kM1].ph_2 - fields[iP1][j][kP1].ph_2)/12
                  + (- fields[i][j][kM1].ph_2 + fields[i][j][kP1].ph_2)*4.0/15;
   dz_ph_3 = 
                    (  fields[iM1][jM2][kM1].ph_3 - fields[iM1][jM2][kP1].ph_3 + fields[iM1][jP2][kM1].ph_3 - fields[iM1][jP2][kP1].ph_3
                     + fields[iP1][jM2][kM1].ph_3 - fields[iP1][jM2][kP1].ph_3 + fields[iP1][jP2][kM1].ph_3 - fields[iP1][jP2][kP1].ph_3
                     + fields[iM2][jM1][kM1].ph_3 - fields[iM2][jM1][kP1].ph_3 + fields[iM2][jP1][kM1].ph_3 - fields[iM2][jP1][kP1].ph_3
                     + fields[iP2][jM1][kM1].ph_3 - fields[iP2][jM1][kP1].ph_3 + fields[iP2][jP1][kM1].ph_3 - fields[iP2][jP1][kP1].ph_3)/120
                  + (  fields[i][jM1][kM2].ph_3 - fields[i][jM1][kP2].ph_3 + fields[i][jP1][kM2].ph_3 - fields[i][jP1][kP2].ph_3
                     + fields[iM1][j][kM2].ph_3 - fields[iM1][j][kP2].ph_3 + fields[iP1][j][kM2].ph_3 - fields[iP1][j][kP2].ph_3
                     - fields[iM1][jM1][kM1].ph_3 + fields[iM1][jM1][kP1].ph_3 - fields[iM1][jP1][kM1].ph_3 + fields[iM1][jP1][kP1].ph_3
                     - fields[iP1][jM1][kM1].ph_3 + fields[iP1][jM1][kP1].ph_3 - fields[iP1][jP1][kM1].ph_3 + fields[iP1][jP1][kP1].ph_3)/30
                  + (- fields[i][j][kM2].ph_3 + fields[i][j][kP2].ph_3)/20
                  - (  fields[i][jM1][kM1].ph_3 - fields[i][jM1][kP1].ph_3 + fields[i][jP1][kM1].ph_3 - fields[i][jP1][kP1].ph_3
                     + fields[iM1][j][kM1].ph_3 - fields[iM1][j][kP1].ph_3 + fields[iP1][j][kM1].ph_3 - fields[iP1][j][kP1].ph_3)/12
                  + (- fields[i][j][kM1].ph_3 + fields[i][j][kP1].ph_3)*4.0/15;
}

void TAux::compute_ph_derivatives(int i, int j, int k, int nSideM1, TFields fields[MAX_NSIDE][MAX_NSIDE][MAX_NSIDE])
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

   compute_ph_derivatives_core(i, j, k, iM1, jM1, kM1, iM2, jM2, kM2, iP1, jP1, kP1, iP2, jP2, kP2, fields);
}

void TAux::compute_ph_derivatives(int i, int j, int k, int nSideM1, TFields fields[MAX_NSIDE][MAX_NSIDE][MAX_NSIDE],
                                  double &out_dx_ph_0, double &out_dx_ph_1, double &out_dx_ph_2, double &out_dx_ph_3,  
                                  double &out_dy_ph_0, double &out_dy_ph_1, double &out_dy_ph_2, double &out_dy_ph_3, 
                                  double &out_dz_ph_0, double &out_dz_ph_1, double &out_dz_ph_2, double &out_dz_ph_3)
{
   compute_ph_derivatives(i, j, k, nSideM1, fields);
   get_ph_derivatives(out_dx_ph_0, out_dx_ph_1, out_dx_ph_2, out_dx_ph_3,  
                      out_dy_ph_0, out_dy_ph_1, out_dy_ph_2, out_dy_ph_3, 
                      out_dz_ph_0, out_dz_ph_1, out_dz_ph_2, out_dz_ph_3);
}

void TAux::get_ph_derivatives(double &out_dx_ph_0, double &out_dx_ph_1, double &out_dx_ph_2, double &out_dx_ph_3,  
                              double &out_dy_ph_0, double &out_dy_ph_1, double &out_dy_ph_2, double &out_dy_ph_3, 
                              double &out_dz_ph_0, double &out_dz_ph_1, double &out_dz_ph_2, double &out_dz_ph_3)
{
   out_dx_ph_0 = dx_ph_0;
   out_dx_ph_1 = dx_ph_1;
   out_dx_ph_2 = dx_ph_2;
   out_dx_ph_3 = dx_ph_3; 
   
   out_dy_ph_0 = dy_ph_0;
   out_dy_ph_1 = dy_ph_1;
   out_dy_ph_2 = dy_ph_2;
   out_dy_ph_3 = dy_ph_3;
   
   out_dz_ph_0 = dz_ph_0;
   out_dz_ph_1 = dz_ph_1;
   out_dz_ph_2 = dz_ph_2;
   out_dz_ph_3 = dz_ph_3;
}

void TAux::compute_A_derivatives_core(int i, int j, int k, 
                                      int iM1, int jM1, int kM1, int iM2, int jM2, int kM2, 
                                      int iP1, int jP1, int kP1, int iP2, int jP2, int kP2, 
                                      TFields fields[MAX_NSIDE][MAX_NSIDE][MAX_NSIDE])
/* IMPORTANT: Does NOT divide by lattice spacing, i.e. if lattice spacing != 1, caller must divide results by it. */
{
/*
   dx_A_1 = 
                   (  fields[iM1][jM1][kM2].A_1 + fields[iM1][jM1][kP2].A_1 + fields[iM1][jP1][kM2].A_1 + fields[iM1][jP1][kP2].A_1
                    + fields[iM1][jM2][kM1].A_1 + fields[iM1][jM2][kP1].A_1 + fields[iM1][jP2][kM1].A_1 + fields[iM1][jP2][kP1].A_1
                    - fields[iP1][jM1][kM2].A_1 - fields[iP1][jM1][kP2].A_1 - fields[iP1][jP1][kM2].A_1 - fields[iP1][jP1][kP2].A_1
                    - fields[iP1][jM2][kM1].A_1 - fields[iP1][jM2][kP1].A_1 - fields[iP1][jP2][kM1].A_1 - fields[iP1][jP2][kP1].A_1)/120
                 + (  fields[iM2][j][kM1].A_1 + fields[iM2][j][kP1].A_1 + fields[iM2][jM1][k].A_1 + fields[iM2][jP1][k].A_1
                    - fields[iP2][j][kM1].A_1 - fields[iP2][j][kP1].A_1 - fields[iP2][jM1][k].A_1 - fields[iP2][jP1][k].A_1
                    - fields[iM1][jM1][kM1].A_1 - fields[iM1][jM1][kP1].A_1 - fields[iM1][jP1][kM1].A_1 - fields[iM1][jP1][kP1].A_1
                    + fields[iP1][jM1][kM1].A_1 + fields[iP1][jM1][kP1].A_1 + fields[iP1][jP1][kM1].A_1 + fields[iP1][jP1][kP1].A_1)/30
                 + (- fields[iM2][j][k].A_1 + fields[iP2][j][k].A_1)/20
                 + (- fields[iM1][j][kM1].A_1 - fields[iM1][j][kP1].A_1 - fields[iM1][jM1][k].A_1 - fields[iM1][jP1][k].A_1
                    + fields[iP1][j][kM1].A_1 + fields[iP1][j][kP1].A_1 + fields[iP1][jM1][k].A_1 + fields[iP1][jP1][k].A_1)/12
                 + (- fields[iM1][j][k].A_1 + fields[iP1][j][k].A_1)*4.0/15;
*/
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
/*
   dy_A_2 = 
                   (  fields[iM1][jM1][kM2].A_2 + fields[iM1][jM1][kP2].A_2 - fields[iM1][jP1][kM2].A_2 - fields[iM1][jP1][kP2].A_2
                    + fields[iP1][jM1][kM2].A_2 + fields[iP1][jM1][kP2].A_2 - fields[iP1][jP1][kM2].A_2 - fields[iP1][jP1][kP2].A_2
                    + fields[iM2][jM1][kM1].A_2 + fields[iM2][jM1][kP1].A_2 - fields[iM2][jP1][kM1].A_2 - fields[iM2][jP1][kP1].A_2
                    + fields[iP2][jM1][kM1].A_2 + fields[iP2][jM1][kP1].A_2 - fields[iP2][jP1][kM1].A_2 - fields[iP2][jP1][kP1].A_2)/120
                 + (  fields[i][jM2][kM1].A_2 + fields[i][jM2][kP1].A_2 - fields[i][jP2][kM1].A_2 - fields[i][jP2][kP1].A_2
                    + fields[iM1][jM2][k].A_2 - fields[iM1][jP2][k].A_2 + fields[iP1][jM2][k].A_2 - fields[iP1][jP2][k].A_2
                    - fields[iM1][jM1][kM1].A_2 - fields[iM1][jM1][kP1].A_2 + fields[iM1][jP1][kM1].A_2 + fields[iM1][jP1][kP1].A_2
                    - fields[iP1][jM1][kM1].A_2 - fields[iP1][jM1][kP1].A_2 + fields[iP1][jP1][kM1].A_2 + fields[iP1][jP1][kP1].A_2)/30
                 + (- fields[i][jM2][k].A_2 + fields[i][jP2][k].A_2)/20
                 - (  fields[i][jM1][kM1].A_2 + fields[i][jM1][kP1].A_2 - fields[i][jP1][kM1].A_2 - fields[i][jP1][kP1].A_2
                    + fields[iM1][jM1][k].A_2 - fields[iM1][jP1][k].A_2 + fields[iP1][jM1][k].A_2 - fields[iP1][jP1][k].A_2)/12
                 + (- fields[i][jM1][k].A_2 + fields[i][jP1][k].A_2)*4.0/15;
*/
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
/*
   dz_A_3 =
                   (  fields[iM1][jM2][kM1].A_3 - fields[iM1][jM2][kP1].A_3 + fields[iM1][jP2][kM1].A_3 - fields[iM1][jP2][kP1].A_3
                    + fields[iP1][jM2][kM1].A_3 - fields[iP1][jM2][kP1].A_3 + fields[iP1][jP2][kM1].A_3 - fields[iP1][jP2][kP1].A_3
                    + fields[iM2][jM1][kM1].A_3 - fields[iM2][jM1][kP1].A_3 + fields[iM2][jP1][kM1].A_3 - fields[iM2][jP1][kP1].A_3
                    + fields[iP2][jM1][kM1].A_3 - fields[iP2][jM1][kP1].A_3 + fields[iP2][jP1][kM1].A_3 - fields[iP2][jP1][kP1].A_3)/120
                 + (  fields[i][jM1][kM2].A_3 - fields[i][jM1][kP2].A_3 + fields[i][jP1][kM2].A_3 - fields[i][jP1][kP2].A_3
                    + fields[iM1][j][kM2].A_3 - fields[iM1][j][kP2].A_3 + fields[iP1][j][kM2].A_3 - fields[iP1][j][kP2].A_3
                    - fields[iM1][jM1][kM1].A_3 + fields[iM1][jM1][kP1].A_3 - fields[iM1][jP1][kM1].A_3 + fields[iM1][jP1][kP1].A_3
                    - fields[iP1][jM1][kM1].A_3 + fields[iP1][jM1][kP1].A_3 - fields[iP1][jP1][kM1].A_3 + fields[iP1][jP1][kP1].A_3)/30
                 + (- fields[i][j][kM2].A_3 + fields[i][j][kP2].A_3)/20
                 - (  fields[i][jM1][kM1].A_3 - fields[i][jM1][kP1].A_3 + fields[i][jP1][kM1].A_3 - fields[i][jP1][kP1].A_3
                    + fields[iM1][j][kM1].A_3 - fields[iM1][j][kP1].A_3 + fields[iP1][j][kM1].A_3 - fields[iP1][j][kP1].A_3)/12
                 + (- fields[i][j][kM1].A_3 + fields[i][j][kP1].A_3)*4.0/15;
*/
}

void TAux::compute_A_derivatives(int i, int j, int k, int nSideM1, TFields fields[MAX_NSIDE][MAX_NSIDE][MAX_NSIDE])
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

   compute_A_derivatives_core(i, j, k, iM1, jM1, kM1, iM2, jM2, kM2, iP1, jP1, kP1, iP2, jP2, kP2, fields);
}

void TAux::compute_A_derivatives(int i, int j, int k, int nSideM1, TFields fields[MAX_NSIDE][MAX_NSIDE][MAX_NSIDE],
                                 double &out_dx_A_2, double &out_dx_A_3, double &out_dy_A_1, 
                                 double &out_dy_A_3, double &out_dz_A_1, double &out_dz_A_2)
{
   compute_A_derivatives(i, j, k, nSideM1, fields);
   get_A_derivatives(out_dx_A_2, out_dx_A_3, out_dy_A_1, 
                     out_dy_A_3, out_dz_A_1, out_dz_A_2);
}

void TAux::get_A_derivatives(double &out_dx_A_2, double &out_dx_A_3, double &out_dy_A_1, 
                             double &out_dy_A_3, double &out_dz_A_1, double &out_dz_A_2)
{
   out_dx_A_2 = dx_A_2;
   out_dx_A_3 = dx_A_3;
   out_dy_A_1 = dy_A_1;
   out_dy_A_3 = dy_A_3;
   out_dz_A_1 = dz_A_1;
   out_dz_A_2 = dz_A_2;
}

double TAux::compute_divA(int i, int j, int k, int nSideM1, TFields fields[MAX_NSIDE][MAX_NSIDE][MAX_NSIDE])
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

   double dA1dx =  (  fields[iM1][jM1][kM2].A_1 + fields[iM1][jM1][kP2].A_1 + fields[iM1][jP1][kM2].A_1 + fields[iM1][jP1][kP2].A_1
                    + fields[iM1][jM2][kM1].A_1 + fields[iM1][jM2][kP1].A_1 + fields[iM1][jP2][kM1].A_1 + fields[iM1][jP2][kP1].A_1
                    - fields[iP1][jM1][kM2].A_1 - fields[iP1][jM1][kP2].A_1 - fields[iP1][jP1][kM2].A_1 - fields[iP1][jP1][kP2].A_1
                    - fields[iP1][jM2][kM1].A_1 - fields[iP1][jM2][kP1].A_1 - fields[iP1][jP2][kM1].A_1 - fields[iP1][jP2][kP1].A_1)/120
                 + (  fields[iM2][j][kM1].A_1 + fields[iM2][j][kP1].A_1 + fields[iM2][jM1][k].A_1 + fields[iM2][jP1][k].A_1
                    - fields[iP2][j][kM1].A_1 - fields[iP2][j][kP1].A_1 - fields[iP2][jM1][k].A_1 - fields[iP2][jP1][k].A_1
                    - fields[iM1][jM1][kM1].A_1 - fields[iM1][jM1][kP1].A_1 - fields[iM1][jP1][kM1].A_1 - fields[iM1][jP1][kP1].A_1
                    + fields[iP1][jM1][kM1].A_1 + fields[iP1][jM1][kP1].A_1 + fields[iP1][jP1][kM1].A_1 + fields[iP1][jP1][kP1].A_1)/30
                 + (- fields[iM2][j][k].A_1 + fields[iP2][j][k].A_1)/20
                 + (- fields[iM1][j][kM1].A_1 - fields[iM1][j][kP1].A_1 - fields[iM1][jM1][k].A_1 - fields[iM1][jP1][k].A_1
                    + fields[iP1][j][kM1].A_1 + fields[iP1][j][kP1].A_1 + fields[iP1][jM1][k].A_1 + fields[iP1][jP1][k].A_1)/12
                 + (- fields[iM1][j][k].A_1 + fields[iP1][j][k].A_1)*4.0/15;

   double dA2dy =  (  fields[iM1][jM1][kM2].A_2 + fields[iM1][jM1][kP2].A_2 - fields[iM1][jP1][kM2].A_2 - fields[iM1][jP1][kP2].A_2
                    + fields[iP1][jM1][kM2].A_2 + fields[iP1][jM1][kP2].A_2 - fields[iP1][jP1][kM2].A_2 - fields[iP1][jP1][kP2].A_2
                    + fields[iM2][jM1][kM1].A_2 + fields[iM2][jM1][kP1].A_2 - fields[iM2][jP1][kM1].A_2 - fields[iM2][jP1][kP1].A_2
                    + fields[iP2][jM1][kM1].A_2 + fields[iP2][jM1][kP1].A_2 - fields[iP2][jP1][kM1].A_2 - fields[iP2][jP1][kP1].A_2)/120
                 + (  fields[i][jM2][kM1].A_2 + fields[i][jM2][kP1].A_2 - fields[i][jP2][kM1].A_2 - fields[i][jP2][kP1].A_2
                    + fields[iM1][jM2][k].A_2 - fields[iM1][jP2][k].A_2 + fields[iP1][jM2][k].A_2 - fields[iP1][jP2][k].A_2
                    - fields[iM1][jM1][kM1].A_2 - fields[iM1][jM1][kP1].A_2 + fields[iM1][jP1][kM1].A_2 + fields[iM1][jP1][kP1].A_2
                    - fields[iP1][jM1][kM1].A_2 - fields[iP1][jM1][kP1].A_2 + fields[iP1][jP1][kM1].A_2 + fields[iP1][jP1][kP1].A_2)/30
                 + (- fields[i][jM2][k].A_2 + fields[i][jP2][k].A_2)/20
                 - (  fields[i][jM1][kM1].A_2 + fields[i][jM1][kP1].A_2 - fields[i][jP1][kM1].A_2 - fields[i][jP1][kP1].A_2
                    + fields[iM1][jM1][k].A_2 - fields[iM1][jP1][k].A_2 + fields[iP1][jM1][k].A_2 - fields[iP1][jP1][k].A_2)/12
                 + (- fields[i][jM1][k].A_2 + fields[i][jP1][k].A_2)*4.0/15;

   double dA3dz =  (  fields[iM1][jM2][kM1].A_3 - fields[iM1][jM2][kP1].A_3 + fields[iM1][jP2][kM1].A_3 - fields[iM1][jP2][kP1].A_3
                    + fields[iP1][jM2][kM1].A_3 - fields[iP1][jM2][kP1].A_3 + fields[iP1][jP2][kM1].A_3 - fields[iP1][jP2][kP1].A_3
                    + fields[iM2][jM1][kM1].A_3 - fields[iM2][jM1][kP1].A_3 + fields[iM2][jP1][kM1].A_3 - fields[iM2][jP1][kP1].A_3
                    + fields[iP2][jM1][kM1].A_3 - fields[iP2][jM1][kP1].A_3 + fields[iP2][jP1][kM1].A_3 - fields[iP2][jP1][kP1].A_3)/120
                 + (  fields[i][jM1][kM2].A_3 - fields[i][jM1][kP2].A_3 + fields[i][jP1][kM2].A_3 - fields[i][jP1][kP2].A_3
                    + fields[iM1][j][kM2].A_3 - fields[iM1][j][kP2].A_3 + fields[iP1][j][kM2].A_3 - fields[iP1][j][kP2].A_3
                    - fields[iM1][jM1][kM1].A_3 + fields[iM1][jM1][kP1].A_3 - fields[iM1][jP1][kM1].A_3 + fields[iM1][jP1][kP1].A_3
                    - fields[iP1][jM1][kM1].A_3 + fields[iP1][jM1][kP1].A_3 - fields[iP1][jP1][kM1].A_3 + fields[iP1][jP1][kP1].A_3)/30
                 + (- fields[i][j][kM2].A_3 + fields[i][j][kP2].A_3)/20
                 - (  fields[i][jM1][kM1].A_3 - fields[i][jM1][kP1].A_3 + fields[i][jP1][kM1].A_3 - fields[i][jP1][kP1].A_3
                    + fields[iM1][j][kM1].A_3 - fields[iM1][j][kP1].A_3 + fields[iP1][j][kM1].A_3 - fields[iP1][j][kP1].A_3)/12
                 + (- fields[i][j][kM1].A_3 + fields[i][j][kP1].A_3)*4.0/15;
                    
   return(dA1dx + dA2dy + dA3dz);
}

void TAux::compute_derivatives(int i, int j, int k, int nSideM1, TFields fields[MAX_NSIDE][MAX_NSIDE][MAX_NSIDE])
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
   
   compute_A_derivatives_core(i, j, k, iM1, jM1, kM1, iM2, jM2, kM2, iP1, jP1, kP1, iP2, jP2, kP2, fields);
   compute_ph_derivatives_core(i, j, k, iM1, jM1, kM1, iM2, jM2, kM2, iP1, jP1, kP1, iP2, jP2, kP2, fields);
}

void TAux::compute_derivatives(int i, int j, int k, int nSideM1, TFields fields[MAX_NSIDE][MAX_NSIDE][MAX_NSIDE],
                               double &out_dx_ph_0, double &out_dx_ph_1, double &out_dx_ph_2, double &out_dx_ph_3, 
                               double &out_dx_A_2, double &out_dx_A_3, 
                               double &out_dy_ph_0, double &out_dy_ph_1, double &out_dy_ph_2, double &out_dy_ph_3, 
                               double &out_dy_A_1, double &out_dy_A_3, 
                               double &out_dz_ph_0, double &out_dz_ph_1, double &out_dz_ph_2, double &out_dz_ph_3, 
                               double &out_dz_A_1, double &out_dz_A_2)
{
   compute_derivatives(i, j, k, nSideM1, fields);
   
   get_A_derivatives(out_dx_A_2, out_dx_A_3, out_dy_A_1, 
                     out_dy_A_3, out_dz_A_1, out_dz_A_2);
   get_ph_derivatives(out_dx_ph_0, out_dx_ph_1, out_dx_ph_2, out_dx_ph_3,  
                      out_dy_ph_0, out_dy_ph_1, out_dy_ph_2, out_dy_ph_3, 
                      out_dz_ph_0, out_dz_ph_1, out_dz_ph_2, out_dz_ph_3);
}

void TAux::get_derivatives(double &out_dx_ph_0, double &out_dx_ph_1, double &out_dx_ph_2, double &out_dx_ph_3, 
                           double &out_dx_A_2, double &out_dx_A_3, 
                           double &out_dy_ph_0, double &out_dy_ph_1, double &out_dy_ph_2, double &out_dy_ph_3, 
                           double &out_dy_A_1, double &out_dy_A_3, 
                           double &out_dz_ph_0, double &out_dz_ph_1, double &out_dz_ph_2, double &out_dz_ph_3, 
                           double &out_dz_A_1, double &out_dz_A_2)
{
   get_A_derivatives(out_dx_A_2, out_dx_A_3, out_dy_A_1, 
                     out_dy_A_3, out_dz_A_1, out_dz_A_2);
   get_ph_derivatives(out_dx_ph_0, out_dx_ph_1, out_dx_ph_2, out_dx_ph_3,  
                      out_dy_ph_0, out_dy_ph_1, out_dy_ph_2, out_dy_ph_3, 
                      out_dz_ph_0, out_dz_ph_1, out_dz_ph_2, out_dz_ph_3);
}

TFields fields[3][MAX_NSIDE][MAX_NSIDE][MAX_NSIDE];      /* Max allowed number of cells always allocated. Ugly but fast. */
TFields momenta[3][MAX_NSIDE][MAX_NSIDE][MAX_NSIDE];     /* Ditto */
TAux aux[MAX_NSIDE][MAX_NSIDE][MAX_NSIDE];               /* Ditto */
 