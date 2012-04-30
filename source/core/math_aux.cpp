/* math_aux.cpp : Auxiliary mathematical definitions and functions.
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

const double PI      = 3.1415926535897931;
const double gR2     = 0.232;              /* sin(theta_Weinberg)^2 = gB^2/(gB^2 + gW^2) */

const double ZERO    = 0.0;
const double ONE     = 1.0;
const double TWO     = 2.0;
const double THREE   = 3.0;
const double FOUR    = 4.0;
const double QUARTER = 0.25;
const double HALF    = 0.5;

#define min(x, y) ((x) < (y) ? (x) : (y))
#define max(x, y) ((x) > (y) ? (x) : (y))

#define sign(x) ((x) < (0) ? (-1) : (1))

#ifndef isnan
#define isnan(x) ((x) != (x))
#endif

double sqr(double x) { return(x*x); }

double spline(double x, double width)
{
   double r = x / width;

   if (r < -HALF) return(ONE);
   if (r > HALF) return(ZERO);

   return(HALF + r * (TWO * r * r - THREE/TWO));
}

double wrap(double angle)
{
   if (fabs(angle) > TWO*PI)
   {
      double intpart;

      return(TWO*PI * (TWO * modf((angle/(FOUR*PI) + HALF), &intpart) - ONE));
   }

   return(angle);   
}

void swap(double &x, double &y)
{
   double z = x;
   x = y;
   y = z;
}

void compute_eigenvalues(double M_11, double M_12, double M_13, double M_22, double M_23, double M_33,
                         double &e_1, double &e_2, double &e_3)
/* Assumes M_ab is symmetric. */
{
   double c0 = M_11*M_22*M_33 + 2*M_12*M_13*M_23 - M_11*M_23*M_23 - M_22*M_13*M_13 - M_33*M_12*M_12;
   double c1 = M_11*M_22 - M_12*M_12 + M_11*M_33 - M_13*M_13 + M_22*M_33 - M_23*M_23;
   double c2 = M_11 + M_22 + M_33;

   double c2Div3 = c2 / 3;

   double aDiv3 = (c1 - c2 * c2Div3) / 3;
   if (aDiv3 > 0) aDiv3 = 0.0;
   double mbDiv2 = (c0 + c2Div3*(2 * c2Div3 * c2Div3 - c1)) / 2;

   double q = mbDiv2 * mbDiv2 + aDiv3 * aDiv3 * aDiv3;
   if (q > 0) q = 0;

   double magnitude = sqrt(-aDiv3);
   double angle = atan2(sqrt(-q), mbDiv2) / 3;
   double cs = cos(angle);
   double sn = sin(angle);

   e_1 = c2Div3 + 2 * magnitude * cs;
   e_2 = c2Div3 - magnitude * (cs + sqrt(3.0) * sn);
   e_3 = c2Div3 - magnitude * (cs - sqrt(3.0) * sn);

   if (e_1 > e_2) swap(e_1, e_2);
   if (e_2 > e_3) 
   {
      swap(e_2, e_3);
      if (e_1 > e_2) swap(e_1, e_2);
   }
}
