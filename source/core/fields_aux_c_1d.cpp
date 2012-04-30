/* fields_aux_c_1d.cpp : Auxiliary functionality for (cartesian) scalars (1D).
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
const int MAX_NSIDE         = 8196; /* Max number of lattice sites per side.
                                       Power of 2 for faster array indexing (maybe). */
#endif
const int DEFAULT_NSIDE     = 8196; /* Default number lattice sites per side. */
const int NCOMPS            = 7;    /* Number of field components. */

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

          dx_A_2, dx_A_3;
   
   TAux();
   
   void clear();
   
   void compute_ph_derivatives_core(int i, int iM1, int iM2, int iP1, int iP2, TFields fields[MAX_NSIDE]);
   void compute_ph_derivatives(int i, int nSideM1, TFields fields[MAX_NSIDE]);
   void compute_ph_derivatives(int i, int nSideM1, TFields fields[MAX_NSIDE],
                               double &out_dx_ph_0, double &out_dx_ph_1, double &out_dx_ph_2, double &out_dx_ph_3);
   void get_ph_derivatives(double &out_dx_ph_0, double &out_dx_ph_1, double &out_dx_ph_2, double &out_dx_ph_3);
   void compute_A_derivatives_core(int i, int iM1, int iM2, int iP1, int iP2, TFields fields[MAX_NSIDE]);
   void compute_A_derivatives(int i, int nSideM1, TFields fields[MAX_NSIDE]);
   void compute_A_derivatives(int i, int nSideM1, TFields fields[MAX_NSIDE], double &out_dx_A_2, double &out_dx_A_3);
   void get_A_derivatives(double &out_dx_A_2, double &out_dx_A_3);
                          
   double compute_divA(int i, int nSideM1, TFields fields[MAX_NSIDE]);

   void compute_derivatives(int i, int nSideM1, TFields fields[MAX_NSIDE]);
   void compute_derivatives(int i, int nSideM1, TFields fields[MAX_NSIDE],
                            double &out_dx_ph_0, double &out_dx_ph_1, double &out_dx_ph_2, double &out_dx_ph_3, 
                            double &out_dx_A_2, double &out_dx_A_3);
   void get_derivatives(double &out_dx_ph_0, double &out_dx_ph_1, double &out_dx_ph_2, double &out_dx_ph_3, 
                        double &out_dx_A_2, double &out_dx_A_3);
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

   dx_A_2 = 0;
   dx_A_3 = 0;
}

void TAux::compute_ph_derivatives_core(int i, int iM1, int iM2, int iP1, int iP2, TFields fields[MAX_NSIDE])
/* IMPORTANT: Does NOT divide by lattice spacing, i.e. if lattice spacing != 1, caller must divide results by it. */
{
   dx_ph_0 = (fields[iM2].ph_0 - fields[iP2].ph_0)/12 + (fields[iP1].ph_0 - fields[iM1].ph_0)*(2.0/3.0);
   dx_ph_1 = (fields[iM2].ph_1 - fields[iP2].ph_1)/12 + (fields[iP1].ph_1 - fields[iM1].ph_1)*(2.0/3.0);
   dx_ph_2 = (fields[iM2].ph_2 - fields[iP2].ph_2)/12 + (fields[iP1].ph_2 - fields[iM1].ph_2)*(2.0/3.0);
   dx_ph_3 = (fields[iM2].ph_3 - fields[iP2].ph_3)/12 + (fields[iP1].ph_3 - fields[iM1].ph_3)*(2.0/3.0);
}

void TAux::compute_ph_derivatives(int i, int nSideM1, TFields fields[MAX_NSIDE])
{
   int iM1, iM2, iP1, iP2;

   if (0 == i) iM1 = nSideM1; else iM1 = i - 1;

   if (0 == iM1) iM2 = nSideM1; else iM2 = iM1 - 1;

   if (nSideM1 == i) iP1 = 0; else iP1 = i + 1;

   if (nSideM1 == iP1) iP2 = 0; else iP2 = iP1 + 1;

   compute_ph_derivatives_core(i, iM1, iM2, iP1, iP2, fields);
}

void TAux::compute_ph_derivatives(int i, int nSideM1, TFields fields[MAX_NSIDE],
                                  double &out_dx_ph_0, double &out_dx_ph_1, double &out_dx_ph_2, double &out_dx_ph_3)
{
   compute_ph_derivatives(i, nSideM1, fields);
   get_ph_derivatives(out_dx_ph_0, out_dx_ph_1, out_dx_ph_2, out_dx_ph_3);
}

void TAux::get_ph_derivatives(double &out_dx_ph_0, double &out_dx_ph_1, double &out_dx_ph_2, double &out_dx_ph_3)
{
   out_dx_ph_0 = dx_ph_0;
   out_dx_ph_1 = dx_ph_1;
   out_dx_ph_2 = dx_ph_2;
   out_dx_ph_3 = dx_ph_3; 
}

void TAux::compute_A_derivatives_core(int i, int iM1, int iM2, int iP1, int iP2, TFields fields[MAX_NSIDE])
/* IMPORTANT: Does NOT divide by lattice spacing, i.e. if lattice spacing != 1, caller must divide results by it. */
{
   dx_A_2 = (fields[iM2].A_2 - fields[iP2].A_2)/12 + (fields[iP1].A_2 - fields[iM1].A_2)*(2.0/3.0);
   dx_A_3 = (fields[iM2].A_3 - fields[iP2].A_3)/12 + (fields[iP1].A_3 - fields[iM1].A_3)*(2.0/3.0);
}

void TAux::compute_A_derivatives(int i, int nSideM1, TFields fields[MAX_NSIDE])
{
   int iM1, iM2, iP1, iP2;

   if (0 == i) iM1 = nSideM1; else iM1 = i - 1;

   if (0 == iM1) iM2 = nSideM1; else iM2 = iM1 - 1;

   if (nSideM1 == i) iP1 = 0; else iP1 = i + 1;

   if (nSideM1 == iP1) iP2 = 0; else iP2 = iP1 + 1;

   compute_A_derivatives_core(i, iM1, iM2, iP1, iP2, fields);
}

void TAux::compute_A_derivatives(int i, int nSideM1, TFields fields[MAX_NSIDE], double &out_dx_A_2, double &out_dx_A_3)
{
   compute_A_derivatives(i, nSideM1, fields);
   get_A_derivatives(out_dx_A_2, out_dx_A_3);
}

void TAux::get_A_derivatives(double &out_dx_A_2, double &out_dx_A_3)
{
   out_dx_A_2 = dx_A_2;
   out_dx_A_3 = dx_A_3;
}

double TAux::compute_divA(int i, int nSideM1, TFields fields[MAX_NSIDE])
{
   int iM1, iM2, iP1, iP2;

   if (0 == i) iM1 = nSideM1; else iM1 = i - 1;

   if (0 == iM1) iM2 = nSideM1; else iM2 = iM1 - 1;

   if (nSideM1 == i) iP1 = 0; else iP1 = i + 1;

   if (nSideM1 == iP1) iP2 = 0; else iP2 = iP1 + 1;

   return( (- fields[iM2].A_1 + fields[iP2].A_1)/12 + (- fields[iM1].A_1 + fields[iP1].A_1)*2.0/3.0 );
}

void TAux::compute_derivatives(int i, int nSideM1, TFields fields[MAX_NSIDE])
{
   int iM1, iM2, iP1, iP2;

   if (0 == i) iM1 = nSideM1; else iM1 = i - 1;

   if (0 == iM1) iM2 = nSideM1; else iM2 = iM1 - 1;

   if (nSideM1 == i) iP1 = 0; else iP1 = i + 1;

   if (nSideM1 == iP1) iP2 = 0; else iP2 = iP1 + 1;
   
   compute_A_derivatives_core(i, iM1, iM2, iP1, iP2, fields);
   compute_ph_derivatives_core(i, iM1, iM2, iP1, iP2, fields);
}

void TAux::compute_derivatives(int i, int nSideM1, TFields fields[MAX_NSIDE],
                               double &out_dx_ph_0, double &out_dx_ph_1, double &out_dx_ph_2, double &out_dx_ph_3, 
                               double &out_dx_A_2, double &out_dx_A_3)
{
   compute_derivatives(i, nSideM1, fields);
   
   get_A_derivatives(out_dx_A_2, out_dx_A_3);
   get_ph_derivatives(out_dx_ph_0, out_dx_ph_1, out_dx_ph_2, out_dx_ph_3);
}

void TAux::get_derivatives(double &out_dx_ph_0, double &out_dx_ph_1, double &out_dx_ph_2, double &out_dx_ph_3, 
                           double &out_dx_A_2, double &out_dx_A_3)
{
   get_A_derivatives(out_dx_A_2, out_dx_A_3);
   get_ph_derivatives(out_dx_ph_0, out_dx_ph_1, out_dx_ph_2, out_dx_ph_3);
}

TFields fields[3][MAX_NSIDE];      /* Max allowed number of cells always allocated. Ugly but fast. */
TFields momenta[3][MAX_NSIDE];     /* Ditto */
TAux aux[MAX_NSIDE];               /* Ditto */
 