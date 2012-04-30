/* analyze_frw.cpp : Analyzes raw (cartesian) field data on FRW metric.
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
#include "core_frw.cpp"

TFields dotFields[MAX_NSIDE][MAX_NSIDE][MAX_NSIDE];   /* Max allowed number of fields always allocated. Ugly but fast. */

typedef struct
{
   double E_x, E_y, E_z, 
          B_x, B_y, B_z, 
          theta_x, theta_y, theta_z, 
          rho_E, rho_B, rho_theta, rho_mixed, rho_lambda, rho, 
          p_E, p_B, p_theta, p_mixed, p, 
          w,
          divE, q,
          m_hi, m_mid, m_lo;
} TObservable;

typedef struct
{
   double A_x, A_y, A_z,
          B_x, B_y, B_z, 
          E_x, E_y, E_z, 
          ph_0, ph_1, ph_2, ph_3;
} TCov;

struct TRho
{
   double BR2, rH, rG, rL;
};

class TAnalyzeLattice : public TLattice
{
   void writeVTKheader(ostream& vtk, string title);
   
public:

   void writeVTKs(string dirPath, string title, int analysis_level);

   void compute_photon_msqr(TFields *dotFields,
                            double h_00, double h_01, double h_02, double h_11, double h_12, double h_22,
                            double dx_ph_0, double dx_ph_1, double dx_ph_2, double dx_ph_3,  
                            double dy_ph_0, double dy_ph_1, double dy_ph_2, double dy_ph_3, 
                            double dz_ph_0, double dz_ph_1, double dz_ph_2, double dz_ph_3,
                            double &M_11, double &M_12, double &M_13, double &M_22, double &M_23, double &M_33);
   double get_rho_static_core(TFields *fields, double twoaR2, double twoaR4, 
                              double h_00, double h_01, double h_02, double h_11, double h_12, double h_22, 
                              double dx_ph_0, double dx_ph_1, double dx_ph_2, double dx_ph_3, 
                              double dx_A_2, double dx_A_3, 
                              double dy_ph_0, double dy_ph_1, double dy_ph_2, double dy_ph_3, 
                              double dy_A_1, double dy_A_3, 
                              double dz_ph_0, double dz_ph_1, double dz_ph_2, double dz_ph_3, 
                              double dz_A_1, double dz_A_2,
                              TRho &rho);
   double compute_q(TFields *fields, TFields *dotFields,
                    double h_00, double h_01, double h_02, double h_11, double h_12, double h_22, 
                    double dx_ph_0, double dx_ph_1, double dx_ph_2, double dx_ph_3,
                    double dy_ph_0, double dy_ph_1, double dy_ph_2, double dy_ph_3,
                    double dz_ph_0, double dz_ph_1, double dz_ph_2, double dz_ph_3);
   void compute_observables(int i, int j, int k, TFields fields[MAX_NSIDE][MAX_NSIDE][MAX_NSIDE], 
                            TFields dotFields[MAX_NSIDE][MAX_NSIDE][MAX_NSIDE], 
                            TObservable& o);
};

void TAnalyzeLattice::compute_photon_msqr(TFields *dotFields,
                                          double h_00, double h_01, double h_02, double h_11, double h_12, double h_22,
                                          double dx_ph_0, double dx_ph_1, double dx_ph_2, double dx_ph_3,  
                                          double dy_ph_0, double dy_ph_1, double dy_ph_2, double dy_ph_3, 
                                          double dz_ph_0, double dz_ph_1, double dz_ph_2, double dz_ph_3,
                                          double &M_11, double &M_12, double &M_13, double &M_22, double &M_23, double &M_33)
/* Returns matrix with squared photon masses. */
{
   double h[4][4];
   
   h[0][0] = h_00;
   h[0][1] = h_01;
   h[0][2] = h_02;
   h[1][1] = h_11;
   h[1][2] = h_12;
   h[2][2] = h_22;
   
   h[0][3] = 0;
   h[1][0] = h[0][1];
   h[1][3] = h[0][2];
   h[2][0] = h[0][2];
   h[2][1] = h[1][2];
   h[2][3] = -h[0][1];
   h[3][0] = h[0][3];
   h[3][1] = h[1][3];
   h[3][2] = h[2][3];
   h[3][3] = h[0][0];

   double dth[4][4]; // First index is Lorentz (0 = time), second is group

   dth[0][0] = dotFields->ph_0 / a;
   dth[0][1] = dotFields->ph_1 / a;
   dth[0][2] = dotFields->ph_2 / a;
   dth[0][3] = dotFields->ph_3 / a;
   
   dth[1][0] = dx_ph_0 / a;
   dth[1][1] = dx_ph_1 / a;
   dth[1][2] = dx_ph_2 / a;
   dth[1][3] = dx_ph_3 / a;
   dth[2][0] = dy_ph_0 / a;
   dth[2][1] = dy_ph_1 / a;
   dth[2][2] = dy_ph_2 / a;
   dth[2][3] = dy_ph_3 / a;
   dth[3][0] = dz_ph_0 / a;
   dth[3][1] = dz_ph_1 / a;
   dth[3][2] = dz_ph_2 / a;
   dth[3][3] = dz_ph_3 / a;

   M_11 = 0;
   M_12 = 0;
   M_13 = 0;
   M_22 = 0;
   M_23 = 0;
   M_33 = 0;
   
   for (int a = 0; a < 4; a++)
       for (int b = 0; b < 4; b++)
       {
          M_11 += h[a][b]*(dth[0][a]*dth[0][b] + dth[2][a]*dth[2][b] + dth[3][a]*dth[3][b]);
          M_12 -= h[a][b]*dth[2][a]*dth[1][b];
          M_13 -= h[a][b]*dth[3][a]*dth[1][b];
          M_22 += h[a][b]*(dth[0][a]*dth[0][b] + dth[1][a]*dth[1][b] + dth[3][a]*dth[3][b]);
          M_23 -= h[a][b]*dth[3][a]*dth[2][b];
          M_33 += h[a][b]*(dth[0][a]*dth[0][b] + dth[1][a]*dth[1][b] + dth[2][a]*dth[2][b]);
       }
}

double TAnalyzeLattice::get_rho_static_core(TFields *fields, double twoaR2, double twoaR4,
                                            double h_00, double h_01, double h_02, double h_11, double h_12, double h_22, 
                                            double dx_ph_0, double dx_ph_1, double dx_ph_2, double dx_ph_3, 
                                            double dx_A_2, double dx_A_3, 
                                            double dy_ph_0, double dy_ph_1, double dy_ph_2, double dy_ph_3, 
                                            double dy_A_1, double dy_A_3, 
                                            double dz_ph_0, double dz_ph_1, double dz_ph_2, double dz_ph_3, 
                                            double dz_A_1, double dz_A_2,
                                            TRho &rho)
/* Individual terms are returned in the TRho struct.  */
{
   double AcrossGradPh_am[4][3];

   AcrossGradPh_am[0][0] = fields->A_2*dz_ph_0 - fields->A_3*dy_ph_0;
   AcrossGradPh_am[0][1] = fields->A_3*dx_ph_0 - fields->A_1*dz_ph_0;
   AcrossGradPh_am[0][2] = fields->A_1*dy_ph_0 - fields->A_2*dx_ph_0;
   
   AcrossGradPh_am[1][0] = fields->A_2*dz_ph_1 - fields->A_3*dy_ph_1;
   AcrossGradPh_am[1][1] = fields->A_3*dx_ph_1 - fields->A_1*dz_ph_1;
   AcrossGradPh_am[1][2] = fields->A_1*dy_ph_1 - fields->A_2*dx_ph_1;
   
   AcrossGradPh_am[2][0] = fields->A_2*dz_ph_2 - fields->A_3*dy_ph_2;
   AcrossGradPh_am[2][1] = fields->A_3*dx_ph_2 - fields->A_1*dz_ph_2;
   AcrossGradPh_am[2][2] = fields->A_1*dy_ph_2 - fields->A_2*dx_ph_2;
   
   AcrossGradPh_am[3][0] = fields->A_2*dz_ph_3 - fields->A_3*dy_ph_3;
   AcrossGradPh_am[3][1] = fields->A_3*dx_ph_3 - fields->A_1*dz_ph_3;
   AcrossGradPh_am[3][2] = fields->A_1*dy_ph_3 - fields->A_2*dx_ph_3;

   rho.BR2 = ( pow(dy_A_3 - dz_A_2, 2) + pow(dz_A_1 - dx_A_3, 2) + pow(dx_A_2 - dy_A_1, 2) ) / twoaR4;

   rho.rH =  (  h_00*(  AcrossGradPh_am[0][0]*AcrossGradPh_am[0][0]
                      + AcrossGradPh_am[0][1]*AcrossGradPh_am[0][1]
                      + AcrossGradPh_am[0][2]*AcrossGradPh_am[0][2])
              + h_11*(  AcrossGradPh_am[1][0]*AcrossGradPh_am[1][0]
                      + AcrossGradPh_am[1][1]*AcrossGradPh_am[1][1]
                      + AcrossGradPh_am[1][2]*AcrossGradPh_am[1][2])
              + h_22*(  AcrossGradPh_am[2][0]*AcrossGradPh_am[2][0]
                      + AcrossGradPh_am[2][1]*AcrossGradPh_am[2][1]
                      + AcrossGradPh_am[2][2]*AcrossGradPh_am[2][2])
              + h_33*(  AcrossGradPh_am[3][0]*AcrossGradPh_am[3][0]
                      + AcrossGradPh_am[3][1]*AcrossGradPh_am[3][1]
                      + AcrossGradPh_am[3][2]*AcrossGradPh_am[3][2])
              + 2*(  h_01*(  AcrossGradPh_am[0][0]*AcrossGradPh_am[1][0]
                           + AcrossGradPh_am[0][1]*AcrossGradPh_am[1][1]
                           + AcrossGradPh_am[0][2]*AcrossGradPh_am[1][2])
                   + h_02*(  AcrossGradPh_am[0][0]*AcrossGradPh_am[2][0]
                           + AcrossGradPh_am[0][1]*AcrossGradPh_am[2][1]
                           + AcrossGradPh_am[0][2]*AcrossGradPh_am[2][2])
                   + h_03*(  AcrossGradPh_am[0][0]*AcrossGradPh_am[3][0]
                           + AcrossGradPh_am[0][1]*AcrossGradPh_am[3][1]
                           + AcrossGradPh_am[0][2]*AcrossGradPh_am[3][2])
                   + h_12*(  AcrossGradPh_am[1][0]*AcrossGradPh_am[2][0]
                           + AcrossGradPh_am[1][1]*AcrossGradPh_am[2][1]
                           + AcrossGradPh_am[1][2]*AcrossGradPh_am[2][2])
                   + h_13*(  AcrossGradPh_am[1][0]*AcrossGradPh_am[3][0]
                           + AcrossGradPh_am[1][1]*AcrossGradPh_am[3][1]
                           + AcrossGradPh_am[1][2]*AcrossGradPh_am[3][2])
                   + h_23*(  AcrossGradPh_am[2][0]*AcrossGradPh_am[3][0]
                           + AcrossGradPh_am[2][1]*AcrossGradPh_am[3][1]
                           + AcrossGradPh_am[2][2]*AcrossGradPh_am[3][2])
                  )
             );

   if (0 > rho.rH) rho.rH = 0; else rho.rH /= twoaR4;

   rho.rG = (  dx_ph_0*dx_ph_0 + dy_ph_0*dy_ph_0 + dz_ph_0*dz_ph_0
             + dx_ph_1*dx_ph_1 + dy_ph_1*dy_ph_1 + dz_ph_1*dz_ph_1
             + dx_ph_2*dx_ph_2 + dy_ph_2*dy_ph_2 + dz_ph_2*dz_ph_2
             + dx_ph_3*dx_ph_3 + dy_ph_3*dy_ph_3 + dz_ph_3*dz_ph_3
            ) / twoaR2;

   rho.rL = pow((lambda / a) * (  pow(fields->ph_0, 2) + pow(fields->ph_1, 2) 
                                + pow(fields->ph_2, 2) + pow(fields->ph_3, 2) - 1 ), 2) / 2;

   return(rho.BR2 + rho.rH + rho.rG + rho.rL);
}

double TAnalyzeLattice::compute_q(TFields *fields, TFields *dotFields, 
                                  double h_00, double h_01, double h_02, double h_11, double h_12, double h_22,
                                  double dx_ph_0, double dx_ph_1, double dx_ph_2, double dx_ph_3,
                                  double dy_ph_0, double dy_ph_1, double dy_ph_2, double dy_ph_3,
                                  double dz_ph_0, double dz_ph_1, double dz_ph_2, double dz_ph_3)
{
   double Ad_ph_0 = fields->A_1 * dx_ph_0 + fields->A_2 * dy_ph_0 + fields->A_3 * dz_ph_0;
   double Ad_ph_1 = fields->A_1 * dx_ph_1 + fields->A_2 * dy_ph_1 + fields->A_3 * dz_ph_1;
   double Ad_ph_2 = fields->A_1 * dx_ph_2 + fields->A_2 * dy_ph_2 + fields->A_3 * dz_ph_2;
   double Ad_ph_3 = fields->A_1 * dx_ph_3 + fields->A_2 * dy_ph_3 + fields->A_3 * dz_ph_3;

   return( (-dotFields->ph_0 * (h_00 * Ad_ph_0 + h_01 * Ad_ph_1 + h_02 * Ad_ph_2 + h_03 * Ad_ph_3) +
            -dotFields->ph_1 * (h_01 * Ad_ph_0 + h_11 * Ad_ph_1 + h_12 * Ad_ph_2 + h_13 * Ad_ph_3) +
            -dotFields->ph_2 * (h_02 * Ad_ph_0 + h_12 * Ad_ph_1 + h_22 * Ad_ph_2 + h_23 * Ad_ph_3) +
            -dotFields->ph_3 * (h_03 * Ad_ph_0 + h_13 * Ad_ph_1 + h_23 * Ad_ph_2 + h_33 * Ad_ph_3)) / (a*a*a) );
}

void TAnalyzeLattice::compute_observables(int i, int j, int k, 
                                          TFields fields[MAX_NSIDE][MAX_NSIDE][MAX_NSIDE], 
                                          TFields dotFields[MAX_NSIDE][MAX_NSIDE][MAX_NSIDE],
                                          TObservable& o)
/* IMPORTANT: Expects correctly normalized derivatives to be ready. */
{
   o.E_x = -dotFields[i][j][k].A_1 / a;
   o.E_y = -dotFields[i][j][k].A_2 / a;
   o.E_z = -dotFields[i][j][k].A_3 / a;

   double dx_ph_0, dx_ph_1, dx_ph_2, dx_ph_3, dx_A_2, dx_A_3, 
          dy_ph_0, dy_ph_1, dy_ph_2, dy_ph_3, dy_A_1, dy_A_3, 
          dz_ph_0, dz_ph_1, dz_ph_2, dz_ph_3, dz_A_1, dz_A_2;
           
   aux[i][j][k].get_derivatives(dx_ph_0, dx_ph_1, dx_ph_2, dx_ph_3, dx_A_2, dx_A_3, 
                                dy_ph_0, dy_ph_1, dy_ph_2, dy_ph_3, dy_A_1, dy_A_3, 
                                dz_ph_0, dz_ph_1, dz_ph_2, dz_ph_3, dz_A_1, dz_A_2);

   double aR2 = a * a;
   double aR4 = aR2 * aR2;
   double twoaR2 = 2 * aR2;
   double twoaR4 = 2 * aR4;
    
   o.B_x = (dy_A_3 - dz_A_2) / aR2;
   o.B_y = (dz_A_1 - dx_A_3) / aR2;
   o.B_z = (dx_A_2 - dy_A_1) / aR2;

   double ph = 1 - sqr(fields[i][j][k].ph_0);
   double trig = 2*acos(fields[i][j][k].ph_0)/sqrt(ph);

   if (isnan(trig))
   {
      trig = 0;
   
      o.theta_x = 0;
      o.theta_y = 0;
      o.theta_z = 0;
      
     // ph_0 ~ 0; ph_0 is the radial direction here, so dotph_0 must be 0 too
   }
   else
   {
      o.theta_x = trig * fields[i][j][k].ph_1;
      o.theta_y = trig * fields[i][j][k].ph_2;
      o.theta_z = trig * fields[i][j][k].ph_3;
   }

   double EU2 = o.E_x * o.E_x + o.E_y * o.E_y + o.E_z * o.E_z;

   double AU2 =   fields[i][j][k].A_1 * fields[i][j][k].A_1
                + fields[i][j][k].A_2 * fields[i][j][k].A_2
                + fields[i][j][k].A_3 * fields[i][j][k].A_3;

   o.rho_E = EU2 / twoaR2;

   double h_00, h_01, h_02, h_11, h_12, h_22;
    
   compute_h(fields[i][j][k].ph_0, fields[i][j][k].ph_1, fields[i][j][k].ph_2, fields[i][j][k].ph_3, 
             h_00, h_01, h_02, h_11, h_12, h_22);

   TRho rho_static;
   get_rho_static_core(&(fields[i][j][k]), twoaR2, twoaR4,
                       h_00, h_01, h_02, h_11, h_12, h_22,
                       dx_ph_0, dx_ph_1, dx_ph_2, dx_ph_3, dx_A_2, dx_A_3, 
                       dy_ph_0, dy_ph_1, dy_ph_2, dy_ph_3, dy_A_1, dy_A_3, 
                       dz_ph_0, dz_ph_1, dz_ph_2, dz_ph_3, dz_A_1, dz_A_2,
                      rho_static);
    
   o.rho_B = rho_static.BR2;

   double rHt = AU2 * (  h_00 * dotFields[i][j][k].ph_0 * dotFields[i][j][k].ph_0
                       + h_11 * dotFields[i][j][k].ph_1 * dotFields[i][j][k].ph_1
                       + h_22 * dotFields[i][j][k].ph_2 * dotFields[i][j][k].ph_2
                       + h_33 * dotFields[i][j][k].ph_3 * dotFields[i][j][k].ph_3
                    
                       + 2 * (  h_01 * dotFields[i][j][k].ph_0 * dotFields[i][j][k].ph_1
                              + h_02 * dotFields[i][j][k].ph_0 * dotFields[i][j][k].ph_2
                              + h_03 * dotFields[i][j][k].ph_0 * dotFields[i][j][k].ph_3
                              + h_12 * dotFields[i][j][k].ph_1 * dotFields[i][j][k].ph_2
                              + h_13 * dotFields[i][j][k].ph_1 * dotFields[i][j][k].ph_3
                              + h_23 * dotFields[i][j][k].ph_2 * dotFields[i][j][k].ph_3 )
                      );

   if (0 > rHt) rHt = 0; else rHt /=  twoaR4;

   double rGt = (  dotFields[i][j][k].ph_0 * dotFields[i][j][k].ph_0
                 + dotFields[i][j][k].ph_1 * dotFields[i][j][k].ph_1
                 + dotFields[i][j][k].ph_2 * dotFields[i][j][k].ph_2
                 + dotFields[i][j][k].ph_3 * dotFields[i][j][k].ph_3 ) / twoaR2;

   o.rho_theta  = rGt + rho_static.rG + rho_static.rL;
   o.rho_lambda = rho_static.rL;
   o.rho_mixed  = rHt + rho_static.rH;
   o.rho = o.rho_E + o.rho_B + o.rho_theta + o.rho_mixed;

   o.p_E = o.rho_E / 3;
   o.p_B = o.rho_B / 3;

   o.p_theta = rGt - rho_static.rL - rho_static.rG / 3;
   o.p_mixed = o.rho_mixed / 3;
   o.p = o.p_E + o.p_B + o.p_theta + o.p_mixed;

   if (0 != o.rho) o.w = o.p / o.rho; else o.w = 0;

   o.divE = -aux[i][j][k].compute_divA(i, j, k, nSideM1, dotFields) / aR2;
   o.q = compute_q(&(fields[i][j][k]), &(dotFields[i][j][k]), 
                   h_00, h_01, h_02, h_11, h_12, h_22, 
                   dx_ph_0, dx_ph_1, dx_ph_2, dx_ph_3,
                   dy_ph_0, dy_ph_1, dy_ph_2, dy_ph_3,
                   dz_ph_0, dz_ph_1, dz_ph_2, dz_ph_3);

   double M_11, M_12, M_13, M_22, M_23, M_33;   
   compute_photon_msqr(&(dotFields[i][j][k]),
                       h_00, h_01, h_02, h_11, h_12, h_22,
                       dx_ph_0, dx_ph_1, dx_ph_2, dx_ph_3,
                       dy_ph_0, dy_ph_1, dy_ph_2, dy_ph_3,
                       dz_ph_0, dz_ph_1, dz_ph_2, dz_ph_3,
                       M_11, M_12, M_13, M_22, M_23, M_33);
   
   compute_eigenvalues(M_11, M_12, M_13, M_22, M_23, M_33, o.m_lo, o.m_mid, o.m_hi);
   
   o.m_lo = sqrt(max(0, o.m_lo));
   o.m_mid = sqrt(max(0, o.m_mid));
   o.m_hi = sqrt(max(0, o.m_hi));   
}

void TAnalyzeLattice::writeVTKheader(ostream& vtk, string title)
{
   vtk.setf(ios::left, ios::adjustfield);
   vtk.precision(OUTPUT_PRECISION);

   vtk << "# vtk DataFile Version 3.0" << endl;
   vtk << title << endl;
   vtk << "ASCII" << endl;
   vtk << "DATASET STRUCTURED_POINTS" << endl;
   vtk << "DIMENSIONS " << get_nSide() << ' ' << get_nSide() << ' ' << get_nSide() << endl;
   vtk << "ORIGIN 0 0 0" << endl;
   vtk << "SPACING 1 1 1" << endl;
   vtk << "POINT_DATA " << get_nSide()*get_nSide()*get_nSide() << endl;
}

void TAnalyzeLattice::writeVTKs(string dirPath, string title, int analysis_level)
{
   ofstream vtk_A, vtk_E, vtk_B, vtk_theta, 
            vtk_rho_E, vtk_rho_B, vtk_rho_theta, vtk_rho_mixed, vtk_rho_lambda, vtk_rho, 
            vtk_p_E, vtk_p_B, vtk_p_theta, vtk_p_mixed, vtk_p, 
            vtk_w, 
            vtk_divE, vtk_q, vtk_divEMq, 
            vtk_m_hi, vtk_m_mid, vtk_m_lo;

   string path_A         = dirPath + "/" + "A.vtk";
   string path_E         = dirPath + "/" + "E.vtk";
   string path_B         = dirPath + "/" + "B.vtk";
   string path_theta     = dirPath + "/" + "theta.vtk";

   string rho_E          = dirPath + "/" + "rho_E.vtk";
   string rho_B          = dirPath + "/" + "rho_B.vtk";
   string rho_theta      = dirPath + "/" + "rho_theta.vtk";
   string rho_mixed      = dirPath + "/" + "rho_mixed.vtk";
   string rho_lambda     = dirPath + "/" + "rho_lambda.vtk";
   string path_rho       = dirPath + "/" + "rho.vtk";
    
   string p_E            = dirPath + "/" + "p_E.vtk";
   string p_B            = dirPath + "/" + "p_B.vtk";
   string p_theta        = dirPath + "/" + "p_theta.vtk";
   string p_mixed        = dirPath + "/" + "p_mixed.vtk";
   string path_p         = dirPath + "/" + "p.vtk";

   string path_w         = dirPath + "/" + "w.vtk";

   string path_divE      = dirPath + "/" + "divE.vtk";
   string path_q         = dirPath + "/" + "q.vtk";
   string path_divEMq    = dirPath + "/" + "divEMq.vtk";    
    
   string path_m_hi      = dirPath + "/" + "m_hi.vtk";
   string path_m_mid     = dirPath + "/" + "m_mid.vtk";
   string path_m_lo      = dirPath + "/" + "m_lo.vtk";

   vtk_A.open(path_A.c_str());
   vtk_E.open(path_E.c_str());
   vtk_B.open(path_B.c_str());
   vtk_theta.open(path_theta.c_str());

   vtk_rho_E.open(rho_E.c_str());
   vtk_rho_B.open(rho_B.c_str());
   vtk_rho_theta.open(rho_theta.c_str());
   vtk_rho_mixed.open(rho_mixed.c_str());
   vtk_rho_lambda.open(rho_lambda.c_str());
   vtk_rho.open(path_rho.c_str());
    
   vtk_p_E.open(p_E.c_str());
   vtk_p_B.open(p_B.c_str());
   vtk_p_theta.open(p_theta.c_str());
   vtk_p_mixed.open(p_mixed.c_str());
   vtk_p.open(path_p.c_str());
    
   vtk_w.open(path_w.c_str());

   vtk_divE.open(path_divE.c_str());
   vtk_q.open(path_q.c_str());
   vtk_divEMq.open(path_divEMq.c_str());

   vtk_m_hi.open(path_m_hi.c_str());
   vtk_m_mid.open(path_m_mid.c_str());
   vtk_m_lo.open(path_m_lo.c_str());

   writeVTKheader(vtk_A, title + "; A"); vtk_A << "VECTORS A DOUBLE" << endl;
   writeVTKheader(vtk_E, title + "; E"); vtk_E << "VECTORS E DOUBLE" << endl;
   writeVTKheader(vtk_B, title + "; B"); vtk_B << "VECTORS B DOUBLE" << endl;
   writeVTKheader(vtk_theta, title + "; theta"); vtk_theta << "VECTORS theta DOUBLE" << endl;

   writeVTKheader(vtk_rho_E, title + "; rho_E"); vtk_rho_E << "SCALARS rho_E DOUBLE" << endl << "LOOKUP_TABLE default" << endl;
   writeVTKheader(vtk_rho_B, title + "; rho_B"); vtk_rho_B << "SCALARS rho_B DOUBLE" << endl << "LOOKUP_TABLE default" << endl;
   writeVTKheader(vtk_rho_theta, title + "; rho_theta"); vtk_rho_theta << "SCALARS rho_theta DOUBLE" << endl << "LOOKUP_TABLE default" << endl;
   writeVTKheader(vtk_rho_lambda, title + "; rho_lambda"); vtk_rho_lambda << "SCALARS rho_lambda DOUBLE" << endl << "LOOKUP_TABLE default" << endl;
   writeVTKheader(vtk_rho_mixed, title + "; rho_mixed"); vtk_rho_mixed << "SCALARS rho_mixed DOUBLE" << endl << "LOOKUP_TABLE default" << endl;
   writeVTKheader(vtk_rho, title + "; rho"); vtk_rho << "SCALARS rho DOUBLE" << endl << "LOOKUP_TABLE default" << endl;

   writeVTKheader(vtk_p_E, title + "; p_E"); vtk_p_E << "SCALARS p_E DOUBLE" << endl << "LOOKUP_TABLE default" << endl;
   writeVTKheader(vtk_p_B, title + "; p_B"); vtk_p_B << "SCALARS p_B DOUBLE" << endl << "LOOKUP_TABLE default" << endl;
   writeVTKheader(vtk_p_theta, title + "; p_theta"); vtk_p_theta << "SCALARS p_theta DOUBLE" << endl << "LOOKUP_TABLE default" << endl;
   writeVTKheader(vtk_p_mixed, title + "; p_mixed"); vtk_p_mixed << "SCALARS p_mixed DOUBLE" << endl << "LOOKUP_TABLE default" << endl;
   writeVTKheader(vtk_p, title + "; p"); vtk_p << "SCALARS p DOUBLE" << endl << "LOOKUP_TABLE default" << endl;

   writeVTKheader(vtk_w, title + "; w"); vtk_w << "SCALARS w DOUBLE" << endl << "LOOKUP_TABLE default" << endl;

   writeVTKheader(vtk_divE, title + "; divE"); vtk_divE << "SCALARS divE DOUBLE" << endl << "LOOKUP_TABLE default" << endl;
   writeVTKheader(vtk_q, title + "; q"); vtk_q << "SCALARS q DOUBLE" << endl << "LOOKUP_TABLE default" << endl;
   writeVTKheader(vtk_divEMq, title + "; div EMq"); vtk_divEMq << "SCALARS divEMq DOUBLE" << endl << "LOOKUP_TABLE default" << endl;

   writeVTKheader(vtk_m_hi, title + "; m_hi"); vtk_m_hi << "SCALARS m_hi DOUBLE" << endl << "LOOKUP_TABLE default" << endl;
   writeVTKheader(vtk_m_mid, title + "; m_mid"); vtk_m_mid << "SCALARS m_mid DOUBLE" << endl << "LOOKUP_TABLE default" << endl;
   writeVTKheader(vtk_m_lo, title + "; m_lo"); vtk_m_lo << "SCALARS m_lo DOUBLE" << endl << "LOOKUP_TABLE default" << endl;
    
   TObservable avg = {};
   TCov mean = {};

   double max_dt_phi = 0;
   int max_dt_phi_i = -1;
   int max_dt_phi_j = -1;
   int max_dt_phi_k = -1;

   double max_dt_A  = 0;
   int max_dt_A_i  = -1;
   int max_dt_A_j  = -1;
   int max_dt_A_k  = -1;

   double max_rel_dt_phi = 0;
   int max_rel_dt_phi_i = -1;
   int max_rel_dt_phi_j = -1;
   int max_rel_dt_phi_k = -1;

   double max_rel_dt_A   = 0;
   int max_rel_dt_A_i  = -1;
   int max_rel_dt_A_j  = -1;
   int max_rel_dt_A_k  = -1;
    
   int nSide = get_nSide();
    
   compute_derivatives(fields[active]);
    
#ifdef _OPENMP
#pragma omp parallel
{
#pragma omp for
#endif
   for (int i = 0; i < nSide; i++)
       for (int j = 0; j < nSide; j++)
           for (int k = 0; k < nSide; k++)
           {
              compute_dotFields(&(fields[active][i][j][k]), &(momenta[active][i][j][k]), dotFields[i][j][k]);
           }

#ifdef _OPENMP
#pragma omp barrier
}
#endif

   for (int i = 0; i < nSide; i++)
       for (int j = 0; j < nSide; j++)
           for (int k = 0; k < nSide; k++)
           {
              TObservable o;
              compute_observables(i, j, k, fields[active], dotFields, o);

              vtk_A     << fields[active][i][j][k].A_1 << ' ' 
                        << fields[active][i][j][k].A_2 << ' ' 
                        << fields[active][i][j][k].A_3 << endl;

              vtk_E     << o.E_x << ' ' << o.E_y << ' ' << o.E_z << endl;
              vtk_B     << o.B_x << ' ' << o.B_y << ' ' << o.B_z << endl;
              
              vtk_theta << o.theta_x << ' ' << o.theta_y << ' ' << o.theta_z << endl;
                        
              vtk_rho_E << o.rho_E << endl;
              vtk_rho_B << o.rho_B << endl;
              vtk_rho_theta << o.rho_theta << endl;
              vtk_rho_lambda << o.rho_lambda << endl;
              vtk_rho_mixed << o.rho_mixed << endl;
              vtk_rho << o.rho << endl;
              
              vtk_p_E << o.p_E << endl;
              vtk_p_B << o.p_B << endl;
              vtk_p_theta << o.p_theta << endl;
              vtk_p_mixed << o.p_mixed << endl;
              vtk_p << o.p << endl;
              
              vtk_w << o.w << endl;
              
              vtk_q << o.q << endl;
              vtk_divE << o.divE << endl;
              vtk_divEMq << o.divE - o.q << endl;
              
              vtk_m_hi << o.m_hi << endl;
              vtk_m_mid << o.m_mid << endl;
              vtk_m_lo << o.m_lo << endl;
              
              double r = max(max(max(abs(dotFields[i][j][k].ph_0),
                                     abs(dotFields[i][j][k].ph_1)), 
                                     abs(dotFields[i][j][k].ph_2)), 
                                     abs(dotFields[i][j][k].ph_3));

              if (r > max_dt_phi)
              {
                 max_dt_phi = r;

                 max_dt_phi_i = i;
                 max_dt_phi_j = j;
                 max_dt_phi_k = k;
              }

              r = max(max(abs(dotFields[i][j][k].A_1), 
                          abs(dotFields[i][j][k].A_2)), 
                          abs(dotFields[i][j][k].A_3));

              if (r > max_dt_A)
              {
                 max_dt_A = r;

                 max_dt_A_i = i;
                 max_dt_A_j = j;
                 max_dt_A_k = k;
              }

              r = max(max(max(abs(dotFields[i][j][k].ph_0 / fields[active][i][j][k].ph_0), 
                              abs(dotFields[i][j][k].ph_1 / fields[active][i][j][k].ph_1)), 
                              abs(dotFields[i][j][k].ph_2 / fields[active][i][j][k].ph_2)), 
                              abs(dotFields[i][j][k].ph_3 / fields[active][i][j][k].ph_3));

              if (r > max_rel_dt_phi)
              {
                 max_rel_dt_phi = r;

                 max_rel_dt_phi_i = i;
                 max_rel_dt_phi_j = j;
                 max_rel_dt_phi_k = k;
              }

              r = max(max(abs(dotFields[i][j][k].A_1 / fields[active][i][j][k].A_1), 
                          abs(dotFields[i][j][k].A_2 / fields[active][i][j][k].A_2)), 
                          abs(dotFields[i][j][k].A_3 / fields[active][i][j][k].A_3));

              if (r > max_rel_dt_A)
              {
                 max_rel_dt_A = r;

                 max_rel_dt_A_i = i;
                 max_rel_dt_A_j = j;
                 max_rel_dt_A_k = k;
              }

              avg.E_x         += o.E_x;
              avg.E_y         += o.E_y;
              avg.E_z         += o.E_z; 
              avg.B_x         += o.B_x; 
              avg.B_y         += o.B_y; 
              avg.B_z         += o.B_z; 
              avg.theta_x     += o.theta_x; 
              avg.theta_y     += o.theta_y; 
              avg.theta_z     += o.theta_z; 
              avg.rho_E       += o.rho_E; 
              avg.rho_B       += o.rho_B; 
              avg.rho_theta   += o.rho_theta; 
              avg.rho_lambda  += o.rho_lambda; 
              avg.rho_mixed   += o.rho_mixed; 
              avg.rho         += o.rho; 
              avg.p_E         += o.p_E; 
              avg.p_B         += o.p_B; 
              avg.p_theta     += o.p_theta; 
              avg.p_mixed     += o.p_mixed; 
              avg.p           += o.p; 
              avg.w           += o.w; 
              avg.divE        += o.divE; 
              avg.q           += o.q;
              avg.m_hi        += o.m_hi;
              avg.m_mid       += o.m_mid;
              avg.m_lo        += o.m_lo;

              mean.A_x        += fields[active][i][j][k].A_1;
              mean.A_y        += fields[active][i][j][k].A_2;
              mean.A_z        += fields[active][i][j][k].A_3;

              mean.ph_0       += fields[active][i][j][k].ph_0;
              mean.ph_1       += fields[active][i][j][k].ph_1;
              mean.ph_2       += fields[active][i][j][k].ph_2;
              mean.ph_3       += fields[active][i][j][k].ph_3;
           }

   vtk_m_lo.close();
   vtk_m_mid.close();
   vtk_m_hi.close();
   vtk_divEMq.close();
   vtk_q.close();
   vtk_divE.close();
   vtk_w.close();
   vtk_p.close();
   vtk_p_mixed.close();
   vtk_p_theta.close();
   vtk_p_B.close();
   vtk_p_E.close();
   vtk_rho.close();    
   vtk_rho_mixed.close();
   vtk_rho_theta.close();
   vtk_rho_lambda.close();
   vtk_rho_B.close();
   vtk_rho_E.close();
   vtk_theta.close();
   vtk_B.close();
   vtk_E.close();
   vtk_A.close();
   
   double nLatticeSites = nSide * nSide * nSide;

   avg.E_x         /= nLatticeSites;
   avg.E_y         /= nLatticeSites;
   avg.E_z         /= nLatticeSites;
   avg.B_x         /= nLatticeSites;
   avg.B_y         /= nLatticeSites;
   avg.B_z         /= nLatticeSites;
   avg.theta_x     /= nLatticeSites;
   avg.theta_y     /= nLatticeSites;
   avg.theta_z     /= nLatticeSites;
   avg.rho_E       /= nLatticeSites;
   avg.rho_B       /= nLatticeSites;
   avg.rho_theta   /= nLatticeSites;
   avg.rho_lambda  /= nLatticeSites;
   avg.rho_mixed   /= nLatticeSites;
   avg.rho         /= nLatticeSites;
   avg.p_E         /= nLatticeSites;
   avg.p_B         /= nLatticeSites;
   avg.p_theta     /= nLatticeSites;
   avg.p_mixed     /= nLatticeSites;
   avg.p           /= nLatticeSites;
   avg.w           /= nLatticeSites;
   avg.divE        /= nLatticeSites;
   avg.q           /= nLatticeSites;
   avg.m_hi        /= nLatticeSites;
   avg.m_mid       /= nLatticeSites;
   avg.m_lo        /= nLatticeSites;

   mean.A_x        /= nLatticeSites;
   mean.A_y        /= nLatticeSites;
   mean.A_z        /= nLatticeSites;

   mean.ph_0       /= nLatticeSites;
   mean.ph_1       /= nLatticeSites;
   mean.ph_2       /= nLatticeSites;
   mean.ph_3       /= nLatticeSites;

   mean.B_x        = avg.B_x;
   mean.B_y		= avg.B_y;
   mean.B_z		= avg.B_z;

   mean.E_x		= avg.E_x;
   mean.E_y		= avg.E_y;
   mean.E_z		= avg.E_z;
    
   ofstream averages;
   
   string path_avg = dirPath + "/" + "averages.txt";
   averages.open(path_avg.c_str());

   averages << "A: " << mean.A_x << ' ' << mean.A_y << ' ' << mean.A_z << endl;
   averages << "E: " << avg.E_x << ' ' << avg.E_y << ' ' << avg.E_z << endl;
   averages << "B: " << avg.B_x << ' ' << avg.B_y << ' ' << avg.B_z << endl;
   averages << "theta: " << avg.theta_x << ' ' << avg.theta_y << ' ' << avg.theta_z << endl;
   averages << "phi: " << mean.ph_0 << ' ' << mean.ph_1 << ' ' << mean.ph_2 << ' ' << mean.ph_3 << endl;
   averages << "rho_E: " << avg.rho_E << endl;
   averages << "rho_B: " << avg.rho_B << endl;
   averages << "rho_theta: " << avg.rho_theta << endl;
   averages << "rho_lambda: " << avg.rho_lambda << endl;
   averages << "rho_mixed: " << avg.rho_mixed << endl;
   averages << "rho: " << avg.rho << endl;
   averages << "p_E: " << avg.p_E << endl;
   averages << "p_B: " << avg.p_B << endl;
   averages << "p_theta: " << avg.p_theta << endl;
   averages << "p_mixed: " << avg.p_mixed << endl;
   averages << "p: " << avg.p << endl;
   averages << "w: " << avg.w << endl;
   averages << "div E: " << avg.divE << endl;
   averages << "q: " << avg.q << endl;
   averages << "m_hi: " << avg.m_hi << endl;
   averages << "m_mid: " << avg.m_mid << endl;
   averages << "m_lo: " << avg.m_lo << endl;
   averages << "max deta_phi_a: " << max_dt_phi << " (" << max_dt_phi_i << ", " << max_dt_phi_j << ", " << max_dt_phi_k << ")" << endl;
   averages << "max deta_A_a: " << max_dt_A << " (" << max_dt_A_i << ", " << max_dt_A_j << ", " << max_dt_A_k << ")" << endl;
   averages << "max deta_phi_a/phi_a: " << max_rel_dt_phi << " (" << max_rel_dt_phi_i << ", " << max_rel_dt_phi_j << ", " << max_rel_dt_phi_k << ")" << endl;
   averages << "max deta_A_a/A_a: " << max_rel_dt_A << " (" << max_rel_dt_A_i << ", " << max_rel_dt_A_j << ", " << max_rel_dt_A_k << ")" << endl;

   averages.close();

   if (1 > analysis_level) return;
   
   int nSideO2 = nSide / 2;

   TObservable var = {};
   TCov var2 = {};

#ifdef _OPENMP
#pragma omp parallel
{
#pragma omp for
#endif
   for (int i = 0; i < nSide; i++)
       for (int j = 0; j < nSide; j++)
           for (int k = 0; k < nSide; k++)
           {
              TObservable o;
              
              compute_observables(i, j, k, fields[active], dotFields, o);

              var.E_x         += sqr(o.E_x - avg.E_x);
              var.E_y         += sqr(o.E_y - avg.E_y);
              var.E_z         += sqr(o.E_z - avg.E_z); 
              var.B_x         += sqr(o.B_x - avg.B_x); 
              var.B_y         += sqr(o.B_y - avg.B_y); 
              var.B_z         += sqr(o.B_z - avg.B_z); 
              var.theta_x     += sqr(o.theta_x - avg.theta_x); 
              var.theta_y     += sqr(o.theta_y - avg.theta_y); 
              var.theta_z     += sqr(o.theta_z - avg.theta_z); 
              var.rho_E       += sqr(o.rho_E - avg.rho_E); 
              var.rho_B       += sqr(o.rho_B - avg.rho_B); 
              var.rho_theta   += sqr(o.rho_theta - avg.rho_theta); 
              var.rho_lambda  += sqr(o.rho_lambda - avg.rho_lambda); 
              var.rho_mixed   += sqr(o.rho_mixed - avg.rho_mixed); 
              var.rho         += sqr(o.rho - avg.rho); 
              var.p_E         += sqr(o.p_E - avg.p_E); 
              var.p_B         += sqr(o.p_B - avg.p_B); 
              var.p_theta     += sqr(o.p_theta - avg.p_theta); 
              var.p_mixed     += sqr(o.p_mixed - avg.p_mixed); 
              var.p           += sqr(o.p - avg.p);
              var.w           += sqr(o.w - avg.w);
              var.divE        += sqr(o.divE - avg.divE);
              var.q           += sqr(o.q - avg.q);
              var.m_hi        += sqr(o.m_hi - avg.m_hi);
              var.m_mid       += sqr(o.m_mid - avg.m_mid);
              var.m_lo        += sqr(o.m_lo - avg.m_lo);

              var2.A_x        += sqr(fields[active][i][j][k].A_1 - mean.A_x);
              var2.A_y        += sqr(fields[active][i][j][k].A_2 - mean.A_y);
              var2.A_z        += sqr(fields[active][i][j][k].A_3 - mean.A_z);

              var2.ph_0       += sqr(fields[active][i][j][k].ph_0 - mean.ph_0);
              var2.ph_1       += sqr(fields[active][i][j][k].ph_1 - mean.ph_1);
              var2.ph_2       += sqr(fields[active][i][j][k].ph_2 - mean.ph_2);
              var2.ph_3       += sqr(fields[active][i][j][k].ph_3 - mean.ph_3);
           }

#ifdef _OPENMP
#pragma omp barrier
}
#endif

   var.E_x        /= nLatticeSites;
   var.E_y        /= nLatticeSites;
   var.E_z        /= nLatticeSites;
   var.B_x        /= nLatticeSites;
   var.B_y        /= nLatticeSites;
   var.B_z        /= nLatticeSites;
   var.theta_x    /= nLatticeSites;
   var.theta_y    /= nLatticeSites;
   var.theta_z    /= nLatticeSites;
   var.rho_E      /= nLatticeSites;
   var.rho_B      /= nLatticeSites;
   var.rho_theta  /= nLatticeSites;
   var.rho_lambda /= nLatticeSites;
   var.rho_mixed  /= nLatticeSites;
   var.rho        /= nLatticeSites;
   var.p_E        /= nLatticeSites;
   var.p_B        /= nLatticeSites;
   var.p_theta    /= nLatticeSites;
   var.p_mixed    /= nLatticeSites;
   var.p          /= nLatticeSites;
   var.w          /= nLatticeSites;
   var.divE       /= nLatticeSites;
   var.q          /= nLatticeSites;
   var.m_hi       /= nLatticeSites;
   var.m_mid      /= nLatticeSites;
   var.m_lo       /= nLatticeSites;

   var2.A_x       /= nLatticeSites;
   var2.A_y       /= nLatticeSites;
   var2.A_z       /= nLatticeSites;

   var2.ph_0      /= nLatticeSites;
   var2.ph_1      /= nLatticeSites;
   var2.ph_2      /= nLatticeSites;
   var2.ph_3      /= nLatticeSites;

   var2.B_x       = var.B_x;
   var2.B_y       = var.B_y;
   var2.B_z       = var.B_z;

   var2.E_x       = var.E_x;
   var2.E_y       = var.E_y;
   var2.E_z       = var.E_z;

   ofstream variances;
   string path_var = dirPath + "/" + "variances.txt";
   variances.open(path_var.c_str());

   variances << "A: " << var2.A_x << ' ' << var2.A_y << ' ' << var2.A_z << endl;
   variances << "E: " << var.E_x << ' ' << var.E_y << ' ' << var.E_z << endl;
   variances << "B: " << var.B_x << ' ' << var.B_y << ' ' << var.B_z << endl;
   variances << "theta: " << var.theta_x << ' ' << var.theta_y << ' ' << var.theta_z << endl;
   variances << "phi: " << var2.ph_0 << ' ' << var2.ph_1 << ' ' << var2.ph_2 << ' ' << var2.ph_3 << endl;
   variances << "rho_E: " << var.rho_E << endl;
   variances << "rho_B: " << var.rho_B << endl;
   variances << "rho_theta: " << var.rho_theta << endl;
   variances << "rho_lambda: " << var.rho_lambda << endl;
   variances << "rho_mixed: " << var.rho_mixed << endl;
   variances << "rho: " << var.rho << endl;
   variances << "p_E: " << var.p_E << endl;
   variances << "p_B: " << var.p_B << endl;
   variances << "p_theta: " << var.p_theta << endl;
   variances << "p_mixed: " << var.p_mixed << endl;
   variances << "p: " << var.p << endl;
   variances << "w: " << var.w << endl;
   variances << "div E: " << var.divE << endl;
   variances << "q: " << var.q << endl;
   variances << "m_hi: " << var.m_hi << endl;
   variances << "m_mid: " << var.m_mid << endl;
   variances << "m_lo: " << var.m_lo << endl;
    
   variances.close();

   if (2 > analysis_level) return;

   TCov cov[MAX_NSIDE] = {};
   int n[MAX_NSIDE] = {};

   double aR2 = a * a;
    
#ifdef _OPENMP
#pragma omp parallel
{
#pragma omp for
#endif
   for (int i_a = 0; i_a < nSide; i_a++)
   {
       int i_a_seq = i_a * nSide;
       
       for (int j_a = 0; j_a < nSide; j_a++)
       {
           int j_a_seq = (j_a + i_a_seq) * nSide;
       
           for (int k_a = 0; k_a < nSide; k_a++)
           {
               int a_seq = k_a + j_a_seq;
           
               double E_x_a = -dotFields[i_a][j_a][k_a].A_1 / a;
               double E_y_a = -dotFields[i_a][j_a][k_a].A_2 / a;
               double E_z_a = -dotFields[i_a][j_a][k_a].A_3 / a;

               double B_x_a = (aux[i_a][j_a][k_a].dy_A_3 - aux[i_a][j_a][k_a].dz_A_2) / aR2;
               double B_y_a = (aux[i_a][j_a][k_a].dz_A_1 - aux[i_a][j_a][k_a].dx_A_3) / aR2;
               double B_z_a = (aux[i_a][j_a][k_a].dx_A_2 - aux[i_a][j_a][k_a].dy_A_1) / aR2;

               for (int i_b = 0; i_b < nSide; i_b++)
               {
                   int i_b_seq = i_b * nSide;
                   int i_distance = (abs(i_b - i_a)) % nSideO2;

                   for (int j_b = 0; j_b < nSide; j_b++)
                   {
                       int j_b_seq = (j_b + i_b_seq) * nSide;
                       int j_distance = (abs(j_b - j_a)) % nSideO2;

                       for (int k_b = 0; k_b < nSide; k_b++)
                       {
                           if ((k_b + j_b_seq) <= a_seq) continue;
                       
                           double distance = sqrt(  sqr(i_distance) + sqr(j_distance) 
                                                  + sqr((abs(k_b - k_a)) % nSideO2) );

                           int distance_bin = int(distance + 0.5);

                           double E_x_b = -dotFields[i_b][j_b][k_b].A_1 / a;
                           double E_y_b = -dotFields[i_b][j_b][k_b].A_2 / a;
                           double E_z_b = -dotFields[i_b][j_b][k_b].A_3 / a;

                           double B_x_b = (aux[i_b][j_b][k_b].dy_A_3 - aux[i_b][j_b][k_b].dz_A_2) / aR2;
                           double B_y_b = (aux[i_b][j_b][k_b].dz_A_1 - aux[i_b][j_b][k_b].dx_A_3) / aR2;
                           double B_z_b = (aux[i_b][j_b][k_b].dx_A_2 - aux[i_b][j_b][k_b].dy_A_1) / aR2;

                           cov[distance_bin].E_x += E_x_a * E_x_b;
                           cov[distance_bin].E_y += E_y_a * E_y_b;
                           cov[distance_bin].E_z += E_z_a * E_z_b;

                           cov[distance_bin].B_x += B_x_a * B_x_b;
                           cov[distance_bin].B_y += B_y_a * B_y_b;
                           cov[distance_bin].B_z += B_z_a * B_z_b;

                           cov[distance_bin].ph_0 += fields[active][i_a][j_a][k_a].ph_0 * fields[active][i_b][j_b][k_b].ph_0;
                           cov[distance_bin].ph_1 += fields[active][i_a][j_a][k_a].ph_1 * fields[active][i_b][j_b][k_b].ph_1;
                           cov[distance_bin].ph_2 += fields[active][i_a][j_a][k_a].ph_2 * fields[active][i_b][j_b][k_b].ph_2;
                           cov[distance_bin].ph_3 += fields[active][i_a][j_a][k_a].ph_3 * fields[active][i_b][j_b][k_b].ph_3;

                           cov[distance_bin].A_x += fields[active][i_a][j_a][k_a].A_1 * fields[active][i_b][j_b][k_b].A_1;
                           cov[distance_bin].A_y += fields[active][i_a][j_a][k_a].A_2 * fields[active][i_b][j_b][k_b].A_2;
                           cov[distance_bin].A_z += fields[active][i_a][j_a][k_a].A_3 * fields[active][i_b][j_b][k_b].A_3;

                           n[distance_bin]++;
                       }
                   }
               }
           }
       }
   }

#ifdef _OPENMP
#pragma omp barrier
}
#endif

   for (int i = 0; i < nSide; i++)
   {
      cov[i].A_x /= n[i]; cov[i].A_x -= sqr(mean.A_x);
      cov[i].A_y /= n[i]; cov[i].A_y -= sqr(mean.A_y);
      cov[i].A_z /= n[i]; cov[i].A_z -= sqr(mean.A_z);

      cov[i].B_x /= n[i]; cov[i].B_x -= sqr(mean.B_x);
      cov[i].B_y /= n[i]; cov[i].B_y -= sqr(mean.B_y);
      cov[i].B_z /= n[i]; cov[i].B_z -= sqr(mean.B_z);

      cov[i].E_x /= n[i]; cov[i].E_x -= sqr(mean.E_x);
      cov[i].E_y /= n[i]; cov[i].E_y -= sqr(mean.E_y);
      cov[i].E_z /= n[i]; cov[i].E_z -= sqr(mean.E_z);

      cov[i].ph_0 /= n[i]; cov[i].ph_0 -= sqr(mean.ph_0);
      cov[i].ph_1 /= n[i]; cov[i].ph_1 -= sqr(mean.ph_1);
      cov[i].ph_2 /= n[i]; cov[i].ph_2 -= sqr(mean.ph_2);
      cov[i].ph_3 /= n[i]; cov[i].ph_3 -= sqr(mean.ph_3);
   }

   ofstream covariance;
   string path_covar = dirPath + "/" + "covariance.txt";
   covariance.open(path_covar.c_str());

   covariance << "Distance, Cov(A), Cov(B), Cov(E), Cov(phi)" << endl;

   for (int i = 0; i < int(nSideO2 * sqrt(3.0) + 0.5); i++)
   {
      covariance << i << ", "
                 << cov[i].A_x + cov[i].A_y + cov[i].A_z << ", "
                 << cov[i].B_x + cov[i].B_y + cov[i].B_z << ", "
                 << cov[i].E_x + cov[i].E_y + cov[i].E_z << ", "
                 << cov[i].ph_0 + cov[i].ph_1 + cov[i].ph_2 + cov[i].ph_3 << endl;
   }

   covariance.close();
}

int main(int argc, char* argv[])
{
   cout.setf(ios::left, ios::adjustfield);
   cout.precision(OUTPUT_PRECISION);

   double lambda = 0;
   int analysis_level = 0;

   string iFile("raw"); // Default input VTK filename
   string oDir = "";    // Output VTK directory
    
   // Parse command line

   int i = 1;
   while (i < argc)
   {
      if (!strcmp("-if", argv[i])) iFile.assign(argv[i + 1]); else
      if (!strcmp("-od", argv[i])) oDir.assign(argv[i + 1]); else

      if (!strcmp("-la", argv[i])) lambda = atof(argv[i + 1]); else
      if (!strcmp("-al", argv[i])) analysis_level = atol(argv[i + 1]); else
      {
         cout << "Unrecognized command line argument: " << argv[i] << endl;
         cout << endl;
         cout << "Supported options [defaults in square brackets]:" << endl;
         cout << endl;
         cout << "Input filename:                      -if [\"" << iFile << "\"]" << endl;
         cout << "Output directory:                    -od [\"" << oDir << "\"]" << endl;
         cout << endl;
         cout << "Lambda (potential constant):         -la [" << lambda << "]" << endl;
         cout << "Analysis level [0..2]:               -al [" << analysis_level << "]" << endl;

         cout << flush;

         exit(1);
      }

      i += 2;
   }

   if ("" == oDir) oDir = stripExtension(iFile);

   // Initialize lattice

   TAnalyzeLattice lattice;

   // Load initial data

   cout << "Reading input file \"" << iFile << ".vtk\"" << endl << flush;

   string title = "";
   double t;
   lattice.read_vtk(iFile + ".vtk", title, t);

   cout << "Done loading." << endl << flush;
   cout << "nSide: " << lattice.get_nSide() << endl;
    
   // Create output files

   cout << "Creating output directory " << oDir << endl << flush;

#ifdef _mkdir
   _mkdir(oDir.c_str());
#else
   string s = "mkdir \"" + oDir + "\"";
   system(s.c_str());
#endif

   cout << "Writing output files..." << endl << flush;

   lattice.set_lambda(lambda);
   lattice.writeVTKs(oDir, title, analysis_level);

   return(0);
}
