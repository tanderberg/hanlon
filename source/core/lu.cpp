/* lu.cpp : Functions for LU decomposition and backsubstitution.
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

const double TINY = 1e-20;

void lubksb_3(double* a, int *indx, double b[])
/* LU back substitution, n = 3. */
{
   int ii = 0, ip;
   double sum;

   ip = indx[0]; 
   sum = b[ip]; 
   b[ip] = b[0];

   if (sum) ii = 1;

   b[0] = sum;

   ip = indx[1];
   sum = b[ip];
   b[ip] = b[1];

   if (ii)
   {
      sum -= a[3] * b[0];
   }
   else
   {
      if (sum) ii = 2;
   }

   b[1] = sum;

   sum = b[2];
   
   if (ii > 0)
   {
	  sum -= a[7] * b[1];
   }

   if (1 == ii)
   {
      sum -= a[6] * b[0];
   }

   b[2] = sum;

   b[2] = b[2] / a[8];
   b[1] = (b[1] - a[5] * b[2]) / a[4];
   b[0] = (b[0] - a[1] * b[1] - a[2] * b[2]) / a[0];
}

void ludcmp_3(double* a, int *indx)
/* LU decomposition, n = 3. */
{
   int imax;
   double big, dum;
   double vv[3];

   big = fabs(a[0]);
   if ((dum = fabs(a[1])) > big) big = dum;
   if ((dum = fabs(a[2])) > big) big = dum;

   vv[0] = 1.0 / big;

   big = fabs(a[3]);
   if ((dum = fabs(a[4])) > big) big = dum;
   if ((dum = fabs(a[5])) > big) big = dum;

   vv[1] = 1.0 / big;

   big = fabs(a[6]);
   if ((dum = fabs(a[7])) > big) big = dum;
   if ((dum = fabs(a[8])) > big) big = dum;

   vv[2] = 1.0 / big;

   big = vv[0] * fabs(a[0]);
   imax = 0;

   if ((dum = vv[1] * fabs(a[3])) >= big)
   {
      big = dum;
      imax = 1;
   }
      
   if ((dum = vv[2] * fabs(a[6])) >= big)
   {
      imax = 2;

      dum = a[6]; a[6] = a[0]; a[0] = dum;
      dum = a[7]; a[7] = a[1]; a[1] = dum;
      dum = a[8]; a[8] = a[2]; a[2] = dum;

      vv[2] = vv[0];
   }
   else if (1 == imax)
   {
      dum = a[3]; a[3] = a[0]; a[0] = dum;
      dum = a[4]; a[4] = a[1]; a[1] = dum;
      dum = a[5]; a[5] = a[2]; a[2] = dum;

      vv[1] = vv[0];
   }

   indx[0] = imax;

   if (0.0 == a[0]) a[0] = TINY;

   a[3] /= a[0];
   a[6] /= a[0];
   
   a[4] -= a[3] * a[1];

   big = vv[1] * fabs(a[4]);
   imax = 1;

   a[7] -= a[6] * a[1];

   if ((dum = vv[2] * fabs(a[7])) >= big)
   {
      imax = 2;

      dum = a[6]; a[6] = a[3]; a[3] = dum;
      dum = a[7]; a[7] = a[4]; a[4] = dum;
      dum = a[8]; a[8] = a[5]; a[5] = dum;
   }

   indx[1] = imax;

   if (0.0 == a[4]) a[4] = TINY;

   a[7] /= a[4];
   
   a[5] -= a[3] * a[2];

   a[8] -= a[6] * a[2] + a[7] * a[5];

   if (0.0 == a[8]) a[8] = TINY;
}

void lubksb4(double* a, int *indx, double b[])
/* LU back substitution, n = 4. */
{
   int ii = 0, ip;
   double sum;

   ip = indx[0]; sum = b[ip]; b[ip] = b[0];
   if (sum) ii = 1;
   b[0] = sum;

   ip = indx[1]; sum = b[ip]; b[ip] = b[1];
   if (ii)
   {
      sum -= a[4] * b[0];
   }
   else
   {
      if (sum) ii = 2;
   }
   b[1] = sum;

   ip = indx[2]; sum = b[ip]; b[ip] = b[2];
   if (ii)
   {
      if (2 > ii) sum -= a[8] * b[0];
      sum -= a[9] * b[1];
   }
   else
   {
      if (sum) ii = 3;
   }
   b[2] = sum;

   sum = b[3];
   if (ii)
   {
      if (2 > ii) sum -= a[12] * b[0];
      if (3 > ii) sum -= a[13] * b[1];
      sum -= a[14] * b[2];
   }
   b[3] = sum / a[15];

   b[2] = (b[2] - a[11] * b[3]) / a[10];
   b[1] = (b[1] - a[6] * b[2] - a[7] * b[3]) / a[5];
   b[0] = (b[0] - a[1] * b[1] - a[2] * b[2] - a[3] * b[3]) / a[0];
}

void ludcmp4(double* a, int *indx)
/* LU decomposition, n = 4. */
{
   int imax;
   double big, dum;
   double vv[4];

   big = fabs(a[0]);
   if ((dum = fabs(a[1])) > big) big = dum;
   if ((dum = fabs(a[2])) > big) big = dum;
   if ((dum = fabs(a[3])) > big) big = dum;
   vv[0] = 1.0 / big;

   big = fabs(a[4]);
   if ((dum = fabs(a[5])) > big) big = dum;
   if ((dum = fabs(a[6])) > big) big = dum;
   if ((dum = fabs(a[7])) > big) big = dum;
   vv[1] = 1.0 / big;

   big = fabs(a[8]);
   if ((dum = fabs(a[9]))  > big) big = dum;
   if ((dum = fabs(a[10])) > big) big = dum;
   if ((dum = fabs(a[11])) > big) big = dum;
   vv[2] = 1.0 / big;

   big = fabs(a[12]);
   if ((dum = fabs(a[13])) > big) big = dum;
   if ((dum = fabs(a[14])) > big) big = dum;
   if ((dum = fabs(a[15])) > big) big = dum;
   vv[3] = 1.0 / big;

   big = vv[0] * fabs(a[0]); imax = 0;
   if ((dum = vv[1] * fabs(a[4]))  >= big) { big = dum; imax = 1; }
   if ((dum = vv[2] * fabs(a[8]))  >= big) { big = dum; imax = 2; }
   if ((dum = vv[3] * fabs(a[12])) >= big) { big = dum; imax = 3; }

   if (imax)
   {
      dum = a[imax*4];     a[imax*4]     = a[0]; a[0] = dum;
      dum = a[imax*4 + 1]; a[imax*4 + 1] = a[1]; a[1] = dum;
      dum = a[imax*4 + 2]; a[imax*4 + 2] = a[2]; a[2] = dum;
      dum = a[imax*4 + 3]; a[imax*4 + 3] = a[3]; a[3] = dum;

      vv[imax] = vv[0];
   }
   indx[0] = imax;

   if (0.0 == a[0]) a[0] = TINY;
      
   a[4]  /= a[0];
   a[8]  /= a[0];
   a[12] /= a[0];

   a[5]  -= a[4]  * a[1]; big = vv[1] * fabs(a[5]); imax = 1;
   a[9]  -= a[8]  * a[1]; if ((dum = vv[2] * fabs(a[9])) >= big)  { big = dum; imax = 2; }
   a[13] -= a[12] * a[1]; if ((dum = vv[3] * fabs(a[13])) >= big) { big = dum; imax = 3; }

   if (1 != imax)
   {
      dum = a[imax*4];     a[imax*4]     = a[4]; a[4] = dum;
      dum = a[imax*4 + 1]; a[imax*4 + 1] = a[5]; a[5] = dum;
      dum = a[imax*4 + 2]; a[imax*4 + 2] = a[6]; a[6] = dum;
      dum = a[imax*4 + 3]; a[imax*4 + 3] = a[7]; a[7] = dum;

      vv[imax] = vv[1];
   }
   indx[1] = imax;

   if (0.0 == a[5]) a[5] = TINY;

   a[9]  /= a[5];
   a[13] /= a[5];

   a[6] -= a[4] * a[2];

   a[10] -= a[8]  * a[2] + a[9] * a[6]; big = vv[2] * fabs(a[10]); imax = 2;
   a[14] -= a[12] * a[2] + a[13] * a[6]; if ((dum = vv[3] * fabs(a[14])) >= big) { big = dum; imax = 3; }

   if (3 == imax)
   {
      dum = a[12]; a[12] = a[8];  a[8]  = dum;
      dum = a[13]; a[13] = a[9];  a[9]  = dum;
      dum = a[14]; a[14] = a[10]; a[10] = dum;
      dum = a[15]; a[15] = a[11]; a[11] = dum;

      vv[3] = vv[2];
   }
   indx[2] = imax;

   if (0.0 == a[10]) a[10] = TINY;

   a[14] /= a[10];

   a[7]  -= a[4]  * a[3];
   a[11] -= a[8]  * a[3] + a[9]  * a[7];
   a[15] -= a[12] * a[3] + a[13] * a[7] + a[14] * a[11];

   if (0.0 == a[15]) a[15] = TINY;
}

void ludcmp7special(double* a, int *indx)
/* LU decomposition, 7x7 input matrix assumed to be diagonal in last three rows. */
{
   int imax;
   double big, dum;
   double vv[4];

   big = fabs(a[0]);
   if ((dum = fabs(a[1])) > big) big = dum;
   if ((dum = fabs(a[2])) > big) big = dum;
   if ((dum = fabs(a[3])) > big) big = dum;
   if ((dum = fabs(a[4])) > big) big = dum;
   if ((dum = fabs(a[5])) > big) big = dum;
   if ((dum = fabs(a[6])) > big) big = dum;
   vv[0] = 1.0 / big;

   big = fabs(a[7]);
   if ((dum = fabs(a[8]))  > big) big = dum;
   if ((dum = fabs(a[9]))  > big) big = dum;
   if ((dum = fabs(a[10])) > big) big = dum;
   if ((dum = fabs(a[11])) > big) big = dum;
   if ((dum = fabs(a[12])) > big) big = dum;
   if ((dum = fabs(a[13])) > big) big = dum;
   vv[1] = 1.0 / big;

   big = fabs(a[14]);
   if ((dum = fabs(a[15])) > big) big = dum;
   if ((dum = fabs(a[16])) > big) big = dum;
   if ((dum = fabs(a[17])) > big) big = dum;
   if ((dum = fabs(a[18])) > big) big = dum;
   if ((dum = fabs(a[19])) > big) big = dum;
   if ((dum = fabs(a[20])) > big) big = dum;
   vv[2] = 1.0 / big;

   big = fabs(a[21]);
   if ((dum = fabs(a[22])) > big) big = dum;
   if ((dum = fabs(a[23])) > big) big = dum;
   if ((dum = fabs(a[24])) > big) big = dum;
   if ((dum = fabs(a[25])) > big) big = dum;
   if ((dum = fabs(a[26])) > big) big = dum;
   if ((dum = fabs(a[27])) > big) big = dum;
   vv[3] = 1.0 / big;

   big = vv[0] * fabs(a[0]); imax = 0;
   if ((dum = vv[1] * fabs(a[7]))  >= big) { big = dum; imax = 1; }
   if ((dum = vv[2] * fabs(a[14])) >= big) { big = dum; imax = 2; }
   if ((dum = vv[3] * fabs(a[21])) >= big) { big = dum; imax = 3; }

   if (imax)
   {
      dum = a[imax*7];     a[imax*7]     = a[0]; a[0] = dum;
      dum = a[imax*7 + 1]; a[imax*7 + 1] = a[1]; a[1] = dum;
      dum = a[imax*7 + 2]; a[imax*7 + 2] = a[2]; a[2] = dum;
      dum = a[imax*7 + 3]; a[imax*7 + 3] = a[3]; a[3] = dum;
      dum = a[imax*7 + 4]; a[imax*7 + 4] = a[4]; a[4] = dum;
      dum = a[imax*7 + 5]; a[imax*7 + 5] = a[5]; a[5] = dum;
      dum = a[imax*7 + 6]; a[imax*7 + 6] = a[6]; a[6] = dum;

      vv[imax] = vv[0];
   }

   indx[0] = imax;

   if (0.0 == a[0]) a[0] = TINY;

   a[7]  /= a[0];
   a[14] /= a[0];   
   a[21] /= a[0];  
 
   a[8]  -= a[7]  * a[1]; big = vv[1] * fabs(a[8]); imax = 1;
   a[15] -= a[14] * a[1]; if ((dum = vv[2] * fabs(a[15])) >= big) { big = dum; imax = 2; }
   a[22] -= a[21] * a[1]; if ((dum = vv[3] * fabs(a[22])) >= big) { big = dum; imax = 3; }

   if (1 != imax)
   {
      dum = a[imax*7];     a[imax*7]     = a[7];  a[7]  = dum;
      dum = a[imax*7 + 1]; a[imax*7 + 1] = a[8];  a[8]  = dum;
      dum = a[imax*7 + 2]; a[imax*7 + 2] = a[9];  a[9]  = dum;
      dum = a[imax*7 + 3]; a[imax*7 + 3] = a[10]; a[10] = dum;
      dum = a[imax*7 + 4]; a[imax*7 + 4] = a[11]; a[11] = dum;
      dum = a[imax*7 + 5]; a[imax*7 + 5] = a[12]; a[12] = dum;
      dum = a[imax*7 + 6]; a[imax*7 + 6] = a[13]; a[13] = dum;

      vv[imax] = vv[1];
   }

   indx[1] = imax;

   if (0.0 == a[8]) a[8] = TINY;

   a[15] /= a[8];
   a[22] /= a[8];
   
   a[9] -= a[7] * a[2];

   a[16] -= a[14] * a[2] + a[15] * a[9]; big = vv[2] * fabs(a[16]); imax = 2;
   a[23] -= a[21] * a[2] + a[22] * a[9]; if ((dum = vv[3] * fabs(a[23])) >= big) { big = dum; imax = 3; }

   if (2 != imax)
   {
      dum = a[imax*7];     a[imax*7]     = a[14]; a[14] = dum;
      dum = a[imax*7 + 1]; a[imax*7 + 1] = a[15]; a[15] = dum;
      dum = a[imax*7 + 2]; a[imax*7 + 2] = a[16]; a[16] = dum;
      dum = a[imax*7 + 3]; a[imax*7 + 3] = a[17]; a[17] = dum;
      dum = a[imax*7 + 4]; a[imax*7 + 4] = a[18]; a[18] = dum;
      dum = a[imax*7 + 5]; a[imax*7 + 5] = a[19]; a[19] = dum;
      dum = a[imax*7 + 6]; a[imax*7 + 6] = a[20]; a[20] = dum;

      vv[imax] = vv[2];
   }

   indx[2] = imax;

   if (0.0 == a[16]) a[16] = TINY;

   a[23] /= a[16];
   
   a[10] -= a[7] * a[3];
   a[17] -= a[14] * a[3] + a[15] * a[10];

   a[24] -= a[21] * a[3] + a[22] * a[10] + a[23] * a[17];

   if (0.0 == a[24]) a[24] = TINY;
   
   a[11] -= a[7] * a[4];
   a[18] -= a[14] * a[4] + a[15] * a[11];
   a[25] -= a[21] * a[4] + a[22] * a[11] + a[23] * a[18];
      
   a[12] -= a[7]  * a[5];
   a[19] -= a[14] * a[5] + a[15] * a[12];
   a[26] -= a[21] * a[5] + a[22] * a[12] + a[23] * a[19];
   
   a[13] -= a[7]  * a[6];
   a[20] -= a[14] * a[6] + a[15] * a[13];
   a[27] -= a[21] * a[6] + a[22] * a[13] + a[23] * a[20];
}

void lubksb7special(double* a, int *indx, double b[])
/* LU back substitution, indx assumed to be trivial for last four rows, with a[32] = a[40] = a[48] = 1. */
{
   int ii = 0, ip;
   double sum;

   ip = indx[0]; sum = b[ip]; b[ip] = b[0];
   if (sum) ii = 1;
   b[0] = sum;
   
   ip = indx[1]; sum = b[ip]; b[ip] = b[1];
   if (ii) sum -= a[7] * b[0]; else if (sum) ii = 2;
   b[1] = sum;

   ip = indx[2]; sum = b[ip]; b[ip] = b[2];
   if (ii)
   {
      if (2 > ii) sum -= a[14] * b[0];
      sum -= a[15] * b[1];
   }
   else
   {
      if (sum) ii = 3;
   }
   b[2] = sum;

   if (ii)
   {
      if (2 > ii) b[3] -= a[21] * b[0];
      if (3 > ii) b[3] -= a[22] * b[1];
      b[3] -= a[23] * b[2];
   }

   b[3] = (b[3] - a[25] * b[4] - a[26] * b[5] - a[27] * b[6]) / a[24];
   b[2] = (b[2] - a[17] * b[3] - a[18] * b[4] - a[19] * b[5] - a[20] * b[6]) / a[16];
   b[1] = (b[1] -  a[9] * b[2] - a[10] * b[3] - a[11] * b[4] - a[12] * b[5] - a[13] * b[6]) / a[8];
   b[0] = (b[0] -  a[1] * b[1] -  a[2] * b[2] -  a[3] * b[3] -  a[4] * b[4] -  a[5] * b[5] - a[6] * b[6]) / a[0];
}

void solve7special(double* A, int *indx, double b[])
/* Solves *special* system of 7 equations A * x = b for x by LU decomposition and backsubstitution.
   Answer x is returned in b, original content of A matrix is overwritten. */
{
   ludcmp7special(A, indx);
   lubksb7special(A, indx, b);
}