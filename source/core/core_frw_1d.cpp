/* core_frw_1d.cpp : Basic 1D lattice for symplectic integration on FRW metric.
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

#include "fields_aux_c_1d.cpp"
#include "lu.cpp"

using namespace std;

#include "string_aux.cpp"

class TLattice
{
   void compute_invM(double ph_0, double ph_1, double ph_2, double ph_3, double AmAmOaR2,
                     double &invM_00, double &invM_01, double &invM_02, double &invM_03,
                     double &invM_11, double &invM_12, double &invM_13,
                     double &invM_22, double &invM_23,
                     double &invM_33);

   void compute_dotMomenta_core(int i, TFields fields[MAX_NSIDE], 
                                TFields* momenta,
                                double invM_00, double invM_01, double invM_02, double invM_03,
                                double invM_11, double invM_12, double invM_13,
                                double invM_22, double invM_23,
                                double invM_33,
                                TFields &dotMomenta);
                                
   void dDotFields_dFields_aR3_core(TFields *fields, TFields *momenta, 
                                    double invM_00, double invM_01, double invM_02, double invM_03,
                                    double invM_11, double invM_12, double invM_13,
                                    double invM_22, double invM_23,
                                    double invM_33,
                                    double dDotFields[NCOMPS][NCOMPS]);

protected:

   int nSideM1, active;
   double lambda, lambdaR2, 
          a_0, rho_r_0, rho_d_0, rho_v, 
          a, rho_r, rho_d, rho_s, rho;

   void compute_h(double ph_0, double ph_1, double ph_2, double ph_3, 
                  double &h_00, double &h_01, double &h_02, double &h_11, double &h_12, double &h_22);

   void compute_invMh(double ph_0, double ph_1, double ph_2, double ph_3, double AmAmOaR2,
                      double &invM_00, double &invM_01, double &invM_02, double &invM_03,
                      double &invM_11, double &invM_12, double &invM_13,
                      double &invM_22, double &invM_23,
                      double &invM_33,
                      double &h_00, double &h_01, double &h_02, double &h_11, double &h_12, double &h_22);

public:

   TLattice();

   void set_nSide(int new_nSide = MAX_NSIDE, bool clear = true);
   int get_nSide() { return(nSideM1 + 1); }
   int get_active() { return(active); }

   void set_lambda(double new_lambda = DEFAULT_LAMBDA) { lambda = new_lambda; lambdaR2 = sqr(lambda); }
   double get_lambda() { return(lambda); }
   
   void set_a(double new_a) { a = new_a; }
   double get_a() { return(a); }
   
   double get_rho_s() { return(rho_s); }
   double get_rho_r() { return(rho_r); }
   double get_rho_d() { return(rho_d); }
   double get_rho() { return(rho); }

   bool write_vtk(string fileName, string title, double t);
   void read_vtk(string fileName, string& title, double& t);
   
   void compute_derivatives(TFields fields[MAX_NSIDE]);

   void compute_dotMomenta(int i, TFields fields[MAX_NSIDE], TFields* momenta, TFields &dotMomenta);
   void compute_dotFields(TFields *fields, TFields *momenta, TFields &dotFields);
   void compute_dots(int i, TFields fields[MAX_NSIDE], TFields* momenta, TFields &dotFields, TFields &dotMomenta);
   void compute_dotFields_and_dDotFields_dFields_aR3(TFields *fields, TFields *momenta, 
                                                     TFields &dotFields, double dDotFields[NCOMPS][NCOMPS]);
   void compute_dotMomenta_and_dDotFields_dFields_aR3(int i, TFields fields[MAX_NSIDE], TFields *momenta, 
                                                      TFields &dotMomenta, double dDotFields[NCOMPS][NCOMPS]);
   void compute_staticDotMomenta(int i, TFields fields[MAX_NSIDE], TFields &dotMomenta);
};

TLattice::TLattice()
{
   active = 0;
   set_nSide();
   set_lambda();
   
   a_0 = 1;
   rho_r_0 = 0;
   rho_d_0 = 0;
   rho_v = 0; 
   
   a = a_0;
   rho_r = rho_r_0;
   rho_d = rho_d_0;
   rho_s = 0;
   
   rho = rho_s + rho_r + rho_d + rho_v;
}

void TLattice::set_nSide(int new_nSide, bool clear)
{
   if (new_nSide < 5) nSideM1 = 4; else
      if (new_nSide > MAX_NSIDE) nSideM1 = MAX_NSIDE - 1; else
         nSideM1 = new_nSide - 1;

   if (clear)
   {
#ifdef _OPENMP
#pragma omp parallel
{
#pragma omp for
#endif
      for (int i = 0; i <= nSideM1; i++)
      {
         fields[active][i].clear();
         momenta[active][i].clear();
      }
#ifdef _OPENMP
#pragma omp barrier
}
#endif
   }
}

bool TLattice::write_vtk(string fileName, string title, double t)
{
   int i;

   ofstream vtk;

   vtk.open(fileName.c_str());
   if (!vtk)
   {
      cout << "ERROR: Could not open VTK file " << fileName << endl << flush;
      return(false);
   }

   vtk.setf(ios::left, ios::adjustfield);
   vtk.precision(OUTPUT_PRECISION);

   vtk << "# vtk DataFile Version 3.0" << endl ;
   vtk << title << "; t = " << t << "; a = " << a 
       << "; rho_r = " << rho_r << "; rho_d = " << rho_d << "; rho_v = " << rho_v << endl;
   vtk << "ASCII" << endl ;
   vtk << "DATASET STRUCTURED_POINTS" << endl ;
   vtk << "DIMENSIONS " << nSideM1+1 << endl;
   vtk << "ORIGIN 0" << endl ;
   vtk << "SPACING " << 1 << endl;
   vtk << "POINT_DATA " << nSideM1+1 << endl;

   vtk << "SCALARS PHI_0 DOUBLE" << endl;
   vtk << "LOOKUP_TABLE default" << endl;

   for (i = 0; i <= nSideM1; i++)
   {
      vtk << fields[active][i].ph_0 << endl;
   }

   vtk << "SCALARS PHI_1 DOUBLE" << endl;
   vtk << "LOOKUP_TABLE default" << endl;

   for (i = 0; i <= nSideM1; i++)
   {
      vtk << fields[active][i].ph_1 << endl;
   }

   vtk << "SCALARS PHI_2 DOUBLE" << endl;
   vtk << "LOOKUP_TABLE default" << endl;

   for (i = 0; i <= nSideM1; i++)
   {
      vtk << fields[active][i].ph_2 << endl;
   }

   vtk << "SCALARS PHI_3 DOUBLE" << endl;
   vtk << "LOOKUP_TABLE default" << endl;

   for (i = 0; i <= nSideM1; i++)
   {
      vtk << fields[active][i].ph_3 << endl;
   }

   vtk << "VECTORS A DOUBLE" << endl;

   for (i = 0; i <= nSideM1; i++)
   {
      vtk << fields[active][i].A_1 << ' ' << fields[active][i].A_2 << ' ' << fields[active][i].A_3 << endl;
   }

   vtk << "SCALARS pPHI_0 DOUBLE" << endl;
   vtk << "LOOKUP_TABLE default" << endl;

   for (i = 0; i <= nSideM1; i++)
   {
      vtk << momenta[active][i].ph_0 << endl;
   }

   vtk << "SCALARS pPHI_1 DOUBLE" << endl;
   vtk << "LOOKUP_TABLE default" << endl;

   for (i = 0; i <= nSideM1; i++)
   {
      vtk << momenta[active][i].ph_1 << endl;
   }

   vtk << "SCALARS pPHI_2 DOUBLE" << endl;
   vtk << "LOOKUP_TABLE default" << endl;

   for (i = 0; i <= nSideM1; i++)
   {
      vtk << momenta[active][i].ph_2 << endl;
   }

   vtk << "SCALARS pPHI_3 DOUBLE" << endl;
   vtk << "LOOKUP_TABLE default" << endl;

   for (i = 0; i <= nSideM1; i++)
   {
      vtk << momenta[active][i].ph_3 << endl;
   }

   vtk << "VECTORS pA DOUBLE" << endl;

   for (i = 0; i <= nSideM1; i++)
   {
      vtk << momenta[active][i].A_1 << ' ' << momenta[active][i].A_2 << ' ' << momenta[active][i].A_3 << endl;
   }

   vtk.close();

   return(true);
}

void TLattice::read_vtk(string fileName, string& title, double& t)
{
   string line;
   ifstream inFile(fileName.c_str());

   int iLine = 0;

   if (inFile.is_open())
   {
      while (!inFile.eof())
      {
         getline(inFile, line);

         iLine++;

         if (2 == iLine)
         {
            title = line;
                
            // Extract t, a, omega_sim, omega_m, omega_l from title line

            int i = line.find("; t = ");                
            if (0 <= i) { i += 5; t = doubleToken(line, i, " ;\n\r"); }

            i = line.find("; a = ");
            if (0 <= i) { i += 5; a = doubleToken(line, i, " ;\n\r"); a_0 = a; }
                
            i = line.find("; rho_r = ");
            if (0 <= i) { i += 9; rho_r = doubleToken(line, i, " ;\n\r"); rho_r_0 = rho_r; }

            i = line.find("; rho_d = ");
            if (0 <= i) { i += 9; rho_d = doubleToken(line, i, " ;\n\r"); rho_d_0 = rho_d; }

            i = line.find("; rho_v = ");
            if (0 <= i) { i += 9; rho_v = doubleToken(line, i, " ;\n\r"); }

            continue;
         }

         line = trim(line);

         if (0 == line.find("DIMENSIONS "))
         {
            string s;

            int i = 10;
            s = token(line, i);

            int j = atoi(s.c_str());

            if ((5 > j) || (j > MAX_NSIDE))
            {
               cout << "Bad x dimension " << s << " (should be 5 to " << MAX_NSIDE << ")" << endl << flush;
               exit(2);
            }
            
            set_nSide(j);

            const string sNSide = s;
         }
         else if ("SCALARS PHI_0 DOUBLE" == line)
         {
            getline(inFile, line); // Skip table line
                
            for (int i = 0; i <= nSideM1; i++)
            {
               string line;
               getline(inFile, line);

               int fromTo = 0;

               fields[active][i].ph_0 = doubleToken(line, fromTo);
            }
         }
         else if ("SCALARS PHI_1 DOUBLE" == line)
         {
            getline(inFile, line); // Skip table line
                
            for (int i = 0; i <= nSideM1; i++)
            {
               string line;
               getline(inFile, line);

               int fromTo = 0;

               fields[active][i].ph_1 = doubleToken(line, fromTo);
            }
         }
         else if ("SCALARS PHI_2 DOUBLE" == line)
         {
            getline(inFile, line); // Skip table line
                
            for (int i = 0; i <= nSideM1; i++)
            {
               string line;
               getline(inFile, line);

               int fromTo = 0;

               fields[active][i].ph_2 = doubleToken(line, fromTo);
            }
         }
         else if ("SCALARS PHI_3 DOUBLE" == line)
         {
            getline(inFile, line); // Skip table line
                
            for (int i = 0; i <= nSideM1; i++)
            {
               string line;
               getline(inFile, line);

               int fromTo = 0;

               fields[active][i].ph_3 = doubleToken(line, fromTo);
            }
         }
         else if ("VECTORS A DOUBLE" == line)
         {
            for (int i = 0; i <= nSideM1; i++)
            {
               string line;
               getline(inFile, line);

               int fromTo = 0;

               fields[active][i].A_1 = doubleToken(line, fromTo);
               fields[active][i].A_2 = doubleToken(line, fromTo);
               fields[active][i].A_3 = doubleToken(line, fromTo);
            }
         }
         else if ("SCALARS pPHI_0 DOUBLE" == line)
         {
            getline(inFile, line); // Skip table line
                
            for (int i = 0; i <= nSideM1; i++)
            {
               string line;
               getline(inFile, line);

               int fromTo = 0;

               momenta[active][i].ph_0 = doubleToken(line, fromTo);
            }
         }
         else if ("SCALARS pPHI_1 DOUBLE" == line)
         {
            getline(inFile, line); // Skip table line
                
            for (int i = 0; i <= nSideM1; i++)
            {
               string line;
               getline(inFile, line);

               int fromTo = 0;

               momenta[active][i].ph_1 = doubleToken(line, fromTo);
            }
         }
         else if ("SCALARS pPHI_2 DOUBLE" == line)
         {
            getline(inFile, line); // Skip table line
                
            for (int i = 0; i <= nSideM1; i++)
            {
               string line;
               getline(inFile, line);

               int fromTo = 0;

               momenta[active][i].ph_2 = doubleToken(line, fromTo);
            }
         }
         else if ("SCALARS pPHI_3 DOUBLE" == line)
         {
            getline(inFile, line); // Skip table line
                
            for (int i = 0; i <= nSideM1; i++)
            {
               string line;
               getline(inFile, line);

               int fromTo = 0;

               momenta[active][i].ph_3 = doubleToken(line, fromTo);
            }
         }
         else if ("VECTORS pA DOUBLE" == line)
         {
            for (int i = 0; i <= nSideM1; i++)
            {
               string line;
               getline(inFile, line);

               int fromTo = 0;

               momenta[active][i].A_1 = doubleToken(line, fromTo);
               momenta[active][i].A_2 = doubleToken(line, fromTo);
               momenta[active][i].A_3 = doubleToken(line, fromTo);
            }
         }
      }

      inFile.close();
   }
}

#define h_03 0
#define h_13 h_02
#define h_23 (-h_01)
#define h_33 h_00

void TLattice::compute_h(double ph_0, double ph_1, double ph_2, double ph_3, 
                         double &h_00, double &h_01, double &h_02, double &h_11, double &h_12, double &h_22)
{
   h_00 = 4 * gR2  * (ph_1*ph_1 + ph_2*ph_2);
   h_01 = 4 * gR2  * (ph_0*ph_1 - ph_2*ph_3);
   h_02 = 4 * gR2  * (ph_0*ph_2 + ph_1*ph_3);

   h_11 = 4 * gR2  * (ph_0*ph_0 + ph_3*ph_3);
   h_12 = 16 * gR2 * ph_1*ph_2;

   h_22 = h_11 + 16 * gR2 * ph_2*ph_2;
   h_11 += 16 * gR2 * ph_1*ph_1;
}                            

void TLattice::compute_invMh(double ph_0, double ph_1, double ph_2, double ph_3, double AmAmOaR2,
                             double &invM_00, double &invM_01, double &invM_02, double &invM_03,
                             double &invM_11, double &invM_12, double &invM_13,
                             double &invM_22, double &invM_23,
                             double &invM_33,
                             double &h_00, double &h_01, double &h_02, double &h_11, double &h_12, double &h_22)
{
   compute_h(ph_0, ph_1, ph_2, ph_3, h_00, h_01, h_02, h_11, h_12, h_22);

   double A[16];

   A[0]  = 1 + AmAmOaR2 * h_00;
   A[1]  = AmAmOaR2 * h_01;
   A[2]  = AmAmOaR2 * h_02;
   A[3]  = AmAmOaR2 * h_03;
   
   A[4]  = A[1];
   A[5]  = 1 + AmAmOaR2 * h_11;
   A[6]  = AmAmOaR2 * h_12;
   A[7]  = AmAmOaR2 * h_13;

   A[8]  = A[2];
   A[9]  = A[6];
   A[10] = 1 + AmAmOaR2 * h_22;
   A[11] = AmAmOaR2 * h_23;
   
   A[12] = A[3];
   A[13] = A[7];
   A[14] = A[11];
   A[15] = 1 + AmAmOaR2 * h_33;
   
   int indx[4];
   
   ludcmp4(A, indx);

   double col[4];
   
   col[0] = 1; col[1] = 0; col[2] = 0; col[3] = 0;
   lubksb4(A, indx, col);
   invM_00 = col[0]; invM_01 = col[1]; invM_02 = col[2]; invM_03 = col[3];

   col[0] = 0; col[1] = 1; col[2] = 0; col[3] = 0;
   lubksb4(A, indx, col);
   invM_11 = col[1]; invM_12 = col[2]; invM_13 = col[3];

   col[0] = 0; col[1] = 0; col[2] = 1; col[3] = 0;
   lubksb4(A, indx, col);
   invM_22 = col[2]; invM_23 = col[3];

   col[0] = 0; col[1] = 0; col[2] = 0; col[3] = 1;
   lubksb4(A, indx, col);
   invM_33 = col[3];
}

void TLattice::compute_invM(double ph_0, double ph_1, double ph_2, double ph_3, double AmAmOaR2,
                            double &invM_00, double &invM_01, double &invM_02, double &invM_03,
                            double &invM_11, double &invM_12, double &invM_13,
                            double &invM_22, double &invM_23,
                            double &invM_33)
{
   double h_00, h_01, h_02, h_11, h_12, h_22;
   
   compute_invMh(ph_0, ph_1, ph_2, ph_3, AmAmOaR2,
                 invM_00, invM_01, invM_02, invM_03, invM_11, invM_12, invM_13, invM_22, invM_23, invM_33,
                 h_00, h_01, h_02, h_11, h_12, h_22);
}

void TLattice::compute_derivatives(TFields fields[MAX_NSIDE])
{
#ifdef _OPENMP
#pragma omp parallel
{
#pragma omp for
#endif
   for (int i = 0; i <= nSideM1; i++)
       aux[i].compute_derivatives(i, nSideM1, fields);
#ifdef _OPENMP
#pragma omp barrier
}
#endif
}

void TLattice::compute_dotMomenta(int i, TFields fields[MAX_NSIDE], TFields* momenta, TFields &dotMomenta)
{
   double invM_00, invM_01, invM_02, invM_03, invM_11, invM_12, invM_13, invM_22, invM_23, invM_33;

   compute_invM(fields[i].ph_0, fields[i].ph_1, fields[i].ph_2, fields[i].ph_3, 
                (sqr(fields[i].A_1) + sqr(fields[i].A_2) + sqr(fields[i].A_3)) / (a * a),
                invM_00, invM_01, invM_02, invM_03, invM_11, invM_12, invM_13, invM_22, invM_23, invM_33);

   compute_dotMomenta_core(i, fields, momenta,
                           invM_00, invM_01, invM_02, invM_03, invM_11, invM_12, invM_13, invM_22, invM_23, invM_33,
                           dotMomenta);
}

void TLattice::compute_dotFields(TFields *fields, TFields *momenta, TFields &dotFields)
{
   // dot A_m = dH/d pA_m
   
   dotFields.A_1 = a * momenta->A_1;
   dotFields.A_2 = a * momenta->A_2;
   dotFields.A_3 = a * momenta->A_3;
   
   // Prepare to compute dot ph

   double invM_00, invM_01, invM_02, invM_03, invM_11, invM_12, invM_13, invM_22, invM_23, invM_33;

   compute_invM(fields->ph_0, fields->ph_1, fields->ph_2, fields->ph_3, 
                (sqr(fields->A_1) + sqr(fields->A_2) + sqr(fields->A_3)) / (a * a),
                invM_00, invM_01, invM_02, invM_03, invM_11, invM_12, invM_13, invM_22, invM_23, invM_33);

   // dot ph_a = dH/d pph_a

   dotFields.ph_0 = (momenta->ph_3*invM_03 + momenta->ph_2*invM_02 + momenta->ph_1*invM_01 + momenta->ph_0*invM_00) / a;
   dotFields.ph_1 = (momenta->ph_3*invM_13 + momenta->ph_2*invM_12 + momenta->ph_1*invM_11 + momenta->ph_0*invM_01) / a;
   dotFields.ph_2 = (momenta->ph_3*invM_23 + momenta->ph_2*invM_22 + momenta->ph_1*invM_12 + momenta->ph_0*invM_02) / a;
   dotFields.ph_3 = (momenta->ph_3*invM_33 + momenta->ph_2*invM_23 + momenta->ph_1*invM_13 + momenta->ph_0*invM_03) / a;
}

void TLattice::compute_dots(int i, TFields fields[MAX_NSIDE], TFields* momenta, TFields &dotFields, TFields &dotMomenta)
{
   // dot A_m = dH/d pA_m
   
   dotFields.A_1 = a * momenta->A_1;
   dotFields.A_2 = a * momenta->A_2;
   dotFields.A_3 = a * momenta->A_3;
   
   // Prepare to compute dot ph

   double invM_00, invM_01, invM_02, invM_03, invM_11, invM_12, invM_13, invM_22, invM_23, invM_33;

   compute_invM(fields[i].ph_0, fields[i].ph_1, fields[i].ph_2, fields[i].ph_3, 
                (sqr(fields[i].A_1) + sqr(fields[i].A_2) + sqr(fields[i].A_3)) / (a * a),
                invM_00, invM_01, invM_02, invM_03, invM_11, invM_12, invM_13, invM_22, invM_23, invM_33);

   // dot ph_a = dH/d pph_a

   dotFields.ph_0 = (momenta->ph_3*invM_03 + momenta->ph_2*invM_02 + momenta->ph_1*invM_01 + momenta->ph_0*invM_00) / a;
   dotFields.ph_1 = (momenta->ph_3*invM_13 + momenta->ph_2*invM_12 + momenta->ph_1*invM_11 + momenta->ph_0*invM_01) / a;
   dotFields.ph_2 = (momenta->ph_3*invM_23 + momenta->ph_2*invM_22 + momenta->ph_1*invM_12 + momenta->ph_0*invM_02) / a;
   dotFields.ph_3 = (momenta->ph_3*invM_33 + momenta->ph_2*invM_23 + momenta->ph_1*invM_13 + momenta->ph_0*invM_03) / a;

   // Compute dot pph, dot pA
   
   compute_dotMomenta_core(i, fields, momenta,
                           invM_00, invM_01, invM_02, invM_03, invM_11, invM_12, invM_13, invM_22, invM_23, invM_33,
                           dotMomenta);
}

void TLattice::compute_dotFields_and_dDotFields_dFields_aR3(TFields *fields, TFields *momenta, 
                                                            TFields &dotFields, double dDotFields[NCOMPS][NCOMPS])
{
   // dot A_m = dH/d pA_m
   
   dotFields.A_1 = a * momenta->A_1;
   dotFields.A_2 = a * momenta->A_2;
   dotFields.A_3 = a * momenta->A_3;
   
   // Prepare to compute dot ph

   double invM_00, invM_01, invM_02, invM_03, invM_11, invM_12, invM_13, invM_22, invM_23, invM_33;

   compute_invM(fields->ph_0, fields->ph_1, fields->ph_2, fields->ph_3, 
                (sqr(fields->A_1) + sqr(fields->A_2) + sqr(fields->A_3)) / (a * a),
                invM_00, invM_01, invM_02, invM_03, invM_11, invM_12, invM_13, invM_22, invM_23, invM_33);

   // dot ph_a = dH/d pph_a

   dotFields.ph_0 = (momenta->ph_3*invM_03 + momenta->ph_2*invM_02 + momenta->ph_1*invM_01 + momenta->ph_0*invM_00) / a;
   dotFields.ph_1 = (momenta->ph_3*invM_13 + momenta->ph_2*invM_12 + momenta->ph_1*invM_11 + momenta->ph_0*invM_01) / a;
   dotFields.ph_2 = (momenta->ph_3*invM_23 + momenta->ph_2*invM_22 + momenta->ph_1*invM_12 + momenta->ph_0*invM_02) / a;
   dotFields.ph_3 = (momenta->ph_3*invM_33 + momenta->ph_2*invM_23 + momenta->ph_1*invM_13 + momenta->ph_0*invM_03) / a;
   
   // Compute d dotFields / dFields

   dDotFields_dFields_aR3_core(fields, momenta,
                               invM_00, invM_01, invM_02, invM_03, invM_11, invM_12, invM_13, invM_22, invM_23, invM_33,
                               dDotFields);
}

void TLattice::compute_dotMomenta_and_dDotFields_dFields_aR3(int i, TFields fields[MAX_NSIDE], TFields *momenta, 
                                                             TFields &dotMomenta, double dDotFields[NCOMPS][NCOMPS])
{
   // Prepare to compute dotMomenta

   double invM_00, invM_01, invM_02, invM_03, invM_11, invM_12, invM_13, invM_22, invM_23, invM_33;

   compute_invM(fields[i].ph_0, fields[i].ph_1, fields[i].ph_2, fields[i].ph_3, 
                (sqr(fields[i].A_1) + sqr(fields[i].A_2) + sqr(fields[i].A_3)) / (a * a),
                invM_00, invM_01, invM_02, invM_03, invM_11, invM_12, invM_13, invM_22, invM_23, invM_33);
   
   // Compute dotMomenta

   compute_dotMomenta_core(i, fields, momenta, 
                           invM_00, invM_01, invM_02, invM_03, invM_11, invM_12, invM_13, invM_22, invM_23, invM_33,
                           dotMomenta);

   // Compute d dotFields / dFields
   
   dDotFields_dFields_aR3_core(&(fields[i]), momenta,
                               invM_00, invM_01, invM_02, invM_03, invM_11, invM_12, invM_13, invM_22, invM_23, invM_33,
                               dDotFields);
}

void TLattice::compute_dotMomenta_core(int i, TFields fields[MAX_NSIDE], 
                                       TFields* momenta, 
                                       double invM_00, double invM_01, double invM_02, double invM_03,
                                       double invM_11, double invM_12, double invM_13,
                                       double invM_22, double invM_23,
                                       double invM_33,
                                       TFields &dotMomenta)
{
   // Get static contribution
   
   compute_staticDotMomenta(i, fields, dotMomenta);
   
   // Prepare to compute kinetic contribution
   
   double ph_0 = fields[i].ph_0;
   double ph_1 = fields[i].ph_1;
   double ph_2 = fields[i].ph_2;
   double ph_3 = fields[i].ph_3;
   
   double h_00, h_01, h_02, h_11, h_12, h_22;
   compute_h(ph_0, ph_1, ph_2, ph_3, h_00, h_01, h_02, h_11, h_12, h_22);
   
   double pph0 = momenta->ph_0;
   double pph1 = momenta->ph_1;
   double pph2 = momenta->ph_2;
   double pph3 = momenta->ph_3;
   
   double poM0 = pph0*invM_00 + pph1*invM_01 + pph2*invM_02 + pph3*invM_03;
   double poM1 = pph0*invM_01 + pph1*invM_11 + pph2*invM_12 + pph3*invM_13;
   double poM2 = pph0*invM_02 + pph1*invM_12 + pph2*invM_22 + pph3*invM_23;
   double poM3 = pph0*invM_03 + pph1*invM_13 + pph2*invM_23 + pph3*invM_33;

   double aR3 = a * a * a;
   
   double r = (  poM0*(h_00*poM0 + h_01*poM1 + h_02*poM2 + h_03*poM3)
               + poM1*(h_01*poM0 + h_11*poM1 + h_12*poM2 + h_13*poM3)
               + poM2*(h_02*poM0 + h_12*poM1 + h_22*poM2 + h_23*poM3)
               + poM3*(h_03*poM0 + h_13*poM1 + h_23*poM2 + h_33*poM3)
              ) / aR3;

   double A_1 = fields[i].A_1;
   double A_2 = fields[i].A_2;
   double A_3 = fields[i].A_3;
   
   // Add kinetic contribution

   dotMomenta.A_1 += r*A_1;
   dotMomenta.A_2 += r*A_2;
   dotMomenta.A_3 += r*A_3;
   
   r = 2 * gR2 * (A_1*A_1 + A_2*A_2 + A_3*A_3) / aR3;

   dotMomenta.ph_0 += r * (  poM0*(                 ph_1*poM1 +     ph_2*poM2              )
                           + poM1*( ph_1*poM0   + 2*ph_0*poM1                   + ph_2*poM3)
                           + poM2*( ph_2*poM0                   + 2*ph_0*poM2   - ph_1*poM3)
                           + poM3*(                 ph_2*poM1     - ph_1*poM2              ) );

   dotMomenta.ph_1 += r * (  poM0*( 2*ph_1*poM0   + ph_0*poM1     + ph_3*poM2              )
                           + poM1*(   ph_0*poM0 + 8*ph_1*poM1   + 4*ph_2*poM2   + ph_3*poM3)
                           + poM2*(   ph_3*poM0 + 4*ph_2*poM1                   - ph_0*poM3)
                           + poM3*(                 ph_3*poM1     - ph_0*poM2 + 2*ph_1*poM3) );

   dotMomenta.ph_2 += r * (  poM0*( 2*ph_2*poM0   - ph_3*poM1     + ph_0*poM2              )
                           + poM1*(  -ph_3*poM0                 + 4*ph_1*poM2   + ph_0*poM3)
                           + poM2*(   ph_0*poM0 + 4*ph_1*poM1   + 8*ph_2*poM2   + ph_3*poM3)
                           + poM3*(                 ph_0*poM1     + ph_3*poM2 + 2*ph_2*poM3) );

   dotMomenta.ph_3 += r * (  poM0*(                -ph_2*poM1 +     ph_1*poM2              )
                           + poM1*(-ph_2*poM0   + 2*ph_3*poM1                   + ph_1*poM3)
                           + poM2*( ph_1*poM0                   + 2*ph_3*poM2   + ph_2*poM3)
                           + poM3*(                 ph_1*poM1     + ph_2*poM2              ) );
}

void TLattice::dDotFields_dFields_aR3_core(TFields *fields, TFields *momenta,
                                           double invM_00, double invM_01, double invM_02, double invM_03,
                                           double invM_11, double invM_12, double invM_13,
                                           double invM_22, double invM_23,
                                           double invM_33,
                                           double dDotFields[NCOMPS][NCOMPS])
// [ [d dotph_0 dph_0, d dot ph_0 dph_1,  d dot ph_0 dph_2,  d dot ph_0 dph_3, d dot ph_0 dA_1, d dot ph_0 dA_2, d dot ph_0 dA_3],
//    [d dotph_1 dph_0, d dot ph_1 dph_1,  d dot ph_1 dph_2,  d dot ph_1 dph_3, d dot ph_1 dA_1, d dot ph_1 dA_2, d dot ph_1 dA_3],
//    [d dotph_2 dph_0, d dot ph_2 dph_1,  d dot ph_2 dph_2,  d dot ph_2 dph_3, d dot ph_2 dA_1, d dot ph_2 dA_2, d dot ph_2 dA_3], 
//    [d dotph_3 dph_0, d dot ph_3 dph_1,  d dot ph_3 dph_2,  d dot ph_3 dph_3, d dot ph_3 dA_1, d dot ph_3 dA_2, d dot ph_3 dA_3],
//    [d dot A_1 dph_0, d dot A_1 dph_1,  d dot A_1 dph_2,  d dot A_1 dph_3, d dot A_1 dA_1, d dot A_1 dA_2, d dot A_1 dA_3],
//    [d dot A_2 dph_0, d dot A_2 dph_1,  d dot A_2 dph_2,  d dot A_2 dph_3, d dot A_2 dA_1, d dot A_2 dA_2, d dot A_2 dA_3],
//    [d dot A_3 dph_0, d dot A_3 dph_1,  d dot A_3 dph_2,  d dot A_3 dph_3, d dot A_3 dA_1, d dot A_3 dA_2, d dot A_3 dA_3] ] 
//
//  The only non-zero entries are pph_i [(M^-1) (dM / dq) (M^-1)]_ia with q = field in column (ph_0, ph_1, ph_2, ph_3, A_1, A_2, A_3)
//   i.e. the three last lines are all zeros.
//
{
   double ph_0 = fields->ph_0;
   double ph_1 = fields->ph_1;
   double ph_2 = fields->ph_2;
   double ph_3 = fields->ph_3;
   
   double A_1 = fields->A_1;
   double A_2 = fields->A_2;
   double A_3 = fields->A_3;
   
   double pph0 = momenta->ph_0;
   double pph1 = momenta->ph_1;
   double pph2 = momenta->ph_2;
   double pph3 = momenta->ph_3;
   
   double r = -4*gR2*(A_1*A_1 + A_2*A_2 + A_3*A_3);

   memset(dDotFields, 0, sizeof(dDotFields));
   
   // ddotph0_dph0 
   dDotFields[0][0] = r*(  2*ph_0*(    pph0*(invM_01*invM_01 + invM_02*invM_02)
                                     + pph1*(invM_01*invM_11 + invM_02*invM_12)
                                     + pph2*(invM_01*invM_12 + invM_02*invM_22)
                                     + pph3*(invM_01*invM_13 + invM_02*invM_23))
                           + ph_2*(  2*pph0*(invM_00*invM_02 + invM_01*invM_03)
                                     + pph1*(invM_00*invM_12 + invM_01*invM_02 + invM_01*invM_13 + invM_11*invM_03)
                                     + pph2*(invM_00*invM_22 + invM_01*invM_23 + invM_03*invM_12 + invM_02*invM_02) 
                                     + pph3*(invM_00*invM_23 + invM_02*invM_03 + invM_01*invM_33 + invM_03*invM_13)) 
                           + ph_1*(  2*pph0*(invM_00*invM_01 - invM_02*invM_03)
                                     + pph1*(invM_00*invM_11 - invM_02*invM_13 - invM_03*invM_12 + invM_01*invM_01)
                                     + pph2*(invM_00*invM_12 + invM_01*invM_02 - invM_02*invM_23 - invM_03*invM_22)
                                     + pph3*(invM_00*invM_13 + invM_01*invM_03 - invM_02*invM_33 - invM_03*invM_23)));

   // ddotph0_dph1 
   dDotFields[0][1] = r*(    ph_0*(  2*pph0*(invM_00*invM_01 - invM_02*invM_03)
                                     + pph1*(invM_00*invM_11 - invM_02*invM_13 - invM_03*invM_12 + invM_01*invM_01)
                                     + pph2*(invM_00*invM_12 + invM_01*invM_02 - invM_02*invM_23 - invM_03*invM_22)
                                     + pph3*(invM_00*invM_13 + invM_01*invM_03 - invM_02*invM_33 - invM_03*invM_23))
                         + 2*ph_1*(    pph0*(invM_00*invM_00 + 4*invM_01*invM_01 + invM_03*invM_03)
                                     + pph1*(invM_00*invM_01 + 4*invM_01*invM_11 + invM_03*invM_13)
                                     + pph2*(invM_00*invM_02 + 4*invM_01*invM_12 + invM_03*invM_23)
                                     + pph3*(invM_00*invM_03 + 4*invM_01*invM_13 + invM_03*invM_33))
                         + 4*ph_2*(  2*pph0*invM_01*invM_02
                                     + pph1*(invM_01*invM_12 + invM_02*invM_11)
                                     + pph2*(invM_01*invM_22 + invM_02*invM_12)
                                     + pph3*(invM_01*invM_23 + invM_02*invM_13))
                           + ph_3*(  2*pph0*(invM_00*invM_02 + invM_01*invM_03)
                                     + pph1*(invM_00*invM_12 + invM_01*invM_02 + invM_01*invM_13 + invM_11*invM_03)
                                     + pph2*(invM_00*invM_22 + invM_01*invM_23 + invM_03*invM_12 + invM_02*invM_02)
                                     + pph3*(invM_00*invM_23 + invM_02*invM_03 + invM_01*invM_33 + invM_03*invM_13)));

   // ddotph0_dph2 
   dDotFields[0][2] = r*(    ph_0*(  2*pph0*(invM_00*invM_02 + invM_01*invM_03)
                                     + pph1*(invM_00*invM_12 + invM_01*invM_02 + invM_01*invM_13 + invM_11*invM_03)
                                     + pph2*(invM_00*invM_22 + invM_01*invM_23 + invM_03*invM_12 + invM_02*invM_02)
                                     + pph3*(invM_00*invM_23 + invM_02*invM_03 + invM_01*invM_33 + invM_03*invM_13))
                         + 4*ph_1*(  2*pph0*invM_01*invM_02
                                     + pph1*(invM_01*invM_12 + invM_02*invM_11)
                                     + pph2*(invM_01*invM_22 + invM_02*invM_12)
                                     + pph3*(invM_01*invM_23 + invM_02*invM_13))
                         + 2*ph_2*(    pph0*(invM_00*invM_00 + 4*invM_02*invM_02 + invM_03*invM_03)
                                     + pph1*(invM_00*invM_01 + 4*invM_02*invM_12 + invM_03*invM_13)
                                     + pph2*(invM_00*invM_02 + 4*invM_02*invM_22 + invM_03*invM_23)
                                     + pph3*(invM_00*invM_03 + 4*invM_02*invM_23 + invM_03*invM_33))
                           + ph_3*(  2*pph0*(invM_02*invM_03 - invM_00*invM_01)
                                     + pph1*(invM_02*invM_13 - invM_00*invM_11 + invM_03*invM_12 - invM_01*invM_01)
                                     + pph2*(invM_02*invM_23 - invM_01*invM_02 - invM_00*invM_12 + invM_03*invM_22)
                                     + pph3*(invM_02*invM_33 - invM_01*invM_03 - invM_00*invM_13 + invM_03*invM_23)));

   // ddotph0_dph3 
   dDotFields[0][3] = r*(    ph_1*(  2*pph0*(invM_00*invM_02 + invM_01*invM_03)
                                     + pph1*(invM_00*invM_12 + invM_01*invM_02 + invM_01*invM_13 + invM_11*invM_03)
                                     + pph2*(invM_00*invM_22 + invM_01*invM_23 + invM_03*invM_12 + invM_02*invM_02)
                                     + pph3*(invM_00*invM_23 + invM_02*invM_03 + invM_01*invM_33 + invM_03*invM_13))
                           + ph_2*(  2*pph0*(invM_02*invM_03 - invM_00*invM_01)
                                     + pph1*(invM_02*invM_13 - invM_00*invM_11 + invM_03*invM_12 - invM_01*invM_01)
                                     + pph2*(invM_02*invM_23 - invM_01*invM_02 - invM_00*invM_12 + invM_03*invM_22)
                                     + pph3*(invM_02*invM_33 - invM_01*invM_03 - invM_00*invM_13 + invM_03*invM_23))
                         + 2*ph_3*(    pph0*(invM_01*invM_01 + invM_02*invM_02)
                                     + pph1*(invM_01*invM_11 + invM_02*invM_12)
                                     + pph2*(invM_01*invM_12 + invM_02*invM_22)
                                     + pph3*(invM_01*invM_13 + invM_02*invM_23)));

   // ddotph1_dph0 
   dDotFields[1][0] = r*(  2*ph_0*(    pph0*(invM_01*invM_11 + invM_02*invM_12)
                                     + pph1*(invM_11*invM_11 + invM_12*invM_12)
                                     + pph2*(invM_11*invM_12 + invM_12*invM_22)
                                     + pph3*(invM_11*invM_13 + invM_12*invM_23))
                           + ph_1*(    pph0*(invM_00*invM_11 - invM_02*invM_13 - invM_03*invM_12 + invM_01*invM_01)
                                   + 2*pph1*(invM_01*invM_11 - invM_12*invM_13)
                                     + pph2*(invM_01*invM_12 + invM_02*invM_11 - invM_12*invM_23 - invM_13*invM_22)
                                     + pph3*(invM_01*invM_13 + invM_11*invM_03 - invM_12*invM_33 - invM_13*invM_23))
                           + ph_2*(    pph0*(invM_00*invM_12 + invM_01*invM_02 + invM_01*invM_13 + invM_11*invM_03)
                                   + 2*pph1*(invM_01*invM_12 + invM_11*invM_13)
                                     + pph2*(invM_01*invM_22 + invM_02*invM_12 + invM_11*invM_23 + invM_12*invM_13)
                                     + pph3*(invM_01*invM_23 + invM_03*invM_12 + invM_11*invM_33 + invM_13*invM_13)));

   // ddotph1_dph1 
   dDotFields[1][1] = r*(    ph_0*(    pph0*(invM_00*invM_11 - invM_02*invM_13 - invM_03*invM_12 + invM_01*invM_01)
                                   + 2*pph1*(invM_01*invM_11 - invM_12*invM_13)
                                     + pph2*(invM_01*invM_12 + invM_02*invM_11 - invM_12*invM_23 - invM_13*invM_22)
                                     + pph3*(invM_01*invM_13 + invM_11*invM_03 - invM_12*invM_33 - invM_13*invM_23))
                         + 2*ph_1*(    pph0*(invM_00*invM_01 + 4*invM_01*invM_11 + invM_03*invM_13)
                                     + pph1*(invM_01*invM_01 + 4*invM_11*invM_11 + invM_13*invM_13)
                                     + pph2*(invM_01*invM_02 + 4*invM_11*invM_12 + invM_13*invM_23)
                                     + pph3*(invM_01*invM_03 + 4*invM_11*invM_13 + invM_13*invM_33))
                         + 4*ph_2*(    pph0*(invM_01*invM_12 + invM_02*invM_11)
                                   + 2*pph1*invM_11*invM_12
                                     + pph2*(invM_11*invM_22 + invM_12*invM_12)
                                     + pph3*(invM_11*invM_23 + invM_12*invM_13))
                           + ph_3*(    pph0*(invM_00*invM_12 + invM_01*invM_02 + invM_01*invM_13 + invM_11*invM_03)
                                   + 2*pph1*(invM_01*invM_12 + invM_11*invM_13)
                                     + pph2*(invM_01*invM_22 + invM_02*invM_12 + invM_11*invM_23 + invM_12*invM_13)
                                     + pph3*(invM_01*invM_23 + invM_03*invM_12 + invM_11*invM_33 + invM_13*invM_13)));

   // ddotph1_dph2 
   dDotFields[1][2] = r*(    ph_0*(    pph0*(invM_00*invM_12 + invM_01*invM_02 + invM_01*invM_13 + invM_11*invM_03)
                                   + 2*pph1*(invM_01*invM_12 + invM_11*invM_13)
                                     + pph2*(invM_01*invM_22 + invM_02*invM_12 + invM_11*invM_23 + invM_12*invM_13)
                                     + pph3*(invM_01*invM_23 + invM_03*invM_12 + invM_11*invM_33 + invM_13*invM_13))
                         + 4*ph_1*(    pph0*(invM_01*invM_12 + invM_02*invM_11)
                                   + 2*pph1*invM_11*invM_12
                                     + pph2*(invM_11*invM_22 + invM_12*invM_12)
                                     + pph3*(invM_11*invM_23 + invM_12*invM_13))
                         + 2*ph_2*(    pph0*(invM_00*invM_01 + 4*invM_02*invM_12 + invM_03*invM_13)
                                     + pph1*(invM_01*invM_01 + 4*invM_12*invM_12 + invM_13*invM_13)
                                     + pph2*(invM_01*invM_02 + 4*invM_12*invM_22 + invM_13*invM_23)
                                     + pph3*(invM_01*invM_03 + 4*invM_12*invM_23 + invM_13*invM_33))
                           + ph_3*(    pph0*(invM_02*invM_13 - invM_00*invM_11 + invM_03*invM_12 - invM_01*invM_01)
                                   + 2*pph1*(invM_12*invM_13 - invM_01*invM_11)
                                     + pph2*(invM_12*invM_23 - invM_02*invM_11 - invM_01*invM_12 + invM_13*invM_22)
                                     + pph3*(invM_12*invM_33 - invM_11*invM_03 - invM_01*invM_13 + invM_13*invM_23)));

   // ddotph1_dph3 
   dDotFields[1][3] = r*(    ph_1*(    pph0*(invM_00*invM_12 + invM_01*invM_02 + invM_01*invM_13 + invM_11*invM_03)
                                   + 2*pph1*(invM_01*invM_12 + invM_11*invM_13)
                                     + pph2*(invM_01*invM_22 + invM_02*invM_12 + invM_11*invM_23 + invM_12*invM_13)
                                     + pph3*(invM_01*invM_23 + invM_03*invM_12 + invM_11*invM_33 + invM_13*invM_13))
                           + ph_2*(    pph0*(invM_02*invM_13 - invM_00*invM_11 + invM_03*invM_12 - invM_01*invM_01)
                                   + 2*pph1*(invM_12*invM_13 - invM_01*invM_11)
                                     + pph2*(invM_12*invM_23 - invM_02*invM_11 - invM_01*invM_12 + invM_13*invM_22)
                                     + pph3*(invM_12*invM_33 - invM_11*invM_03 - invM_01*invM_13 + invM_13*invM_23))
                         + 2*ph_3*(    pph0*(invM_01*invM_11 + invM_02*invM_12)
                                     + pph1*(invM_11*invM_11 + invM_12*invM_12)
                                     + pph2*(invM_11*invM_12 + invM_12*invM_22)
                                     + pph3*(invM_11*invM_13 + invM_12*invM_23)));

   // ddotph2_dph0
   dDotFields[2][0] = r*(  2*ph_0*(    pph0*(invM_01*invM_12 + invM_02*invM_22)
                                     + pph1*(invM_11*invM_12 + invM_12*invM_22)
                                     + pph2*(invM_12*invM_12 + invM_22*invM_22)
                                     + pph3*(invM_12*invM_13 + invM_22*invM_23))
                           + ph_1*(    pph0*(invM_00*invM_12 + invM_01*invM_02 - invM_02*invM_23 - invM_03*invM_22)
                                     + pph1*(invM_01*invM_12 + invM_02*invM_11 - invM_12*invM_23 - invM_13*invM_22)
                                   + 2*pph2*(invM_02*invM_12 - invM_22*invM_23)
                                     + pph3*(invM_02*invM_13 + invM_03*invM_12 - invM_22*invM_33 - invM_23*invM_23))
                           + ph_2*(    pph0*(invM_00*invM_22 + invM_01*invM_23 + invM_03*invM_12 + invM_02*invM_02)
                                     + pph1*(invM_01*invM_22 + invM_02*invM_12 + invM_11*invM_23 + invM_12*invM_13)
                                   + 2*pph2*(invM_02*invM_22 + invM_12*invM_23)
                                     + pph3*(invM_02*invM_23 + invM_03*invM_22 + invM_12*invM_33 + invM_13*invM_23)));

   // ddotph2_dph1
   dDotFields[2][1] = r*(    ph_0*(    pph0*(invM_00*invM_12 + invM_01*invM_02 - invM_02*invM_23 - invM_03*invM_22)
                                     + pph1*(invM_01*invM_12 + invM_02*invM_11 - invM_12*invM_23 - invM_13*invM_22)
                                   + 2*pph2*(invM_02*invM_12 - invM_22*invM_23)
                                     + pph3*(invM_02*invM_13 + invM_03*invM_12 - invM_22*invM_33 - invM_23*invM_23))
                         + 2*ph_1*(    pph0*(invM_00*invM_02 + 4*invM_01*invM_12 + invM_03*invM_23)
                                     + pph1*(invM_01*invM_02 + 4*invM_11*invM_12 + invM_13*invM_23)
                                     + pph2*(invM_02*invM_02 + 4*invM_12*invM_12 + invM_23*invM_23)
                                     + pph3*(invM_02*invM_03 + 4*invM_12*invM_13 + invM_23*invM_33))
                         + 4*ph_2*(    pph0*(invM_01*invM_22 + invM_02*invM_12)
                                     + pph1*(invM_11*invM_22 + invM_12*invM_12)
                                   + 2*pph2*invM_12*invM_22
                                     + pph3*(invM_12*invM_23 + invM_13*invM_22))
                           + ph_3*(    pph0*(invM_00*invM_22 + invM_01*invM_23 + invM_03*invM_12 + invM_02*invM_02)
                                     + pph1*(invM_01*invM_22 + invM_02*invM_12 + invM_11*invM_23 + invM_12*invM_13)
                                   + 2*pph2*(invM_02*invM_22 + invM_12*invM_23)
                                     + pph3*(invM_02*invM_23 + invM_03*invM_22 + invM_12*invM_33 + invM_13*invM_23)));

   // ddotph2_dph2 
   dDotFields[2][2] = r*(    ph_0*(    pph0*(invM_00*invM_22 + invM_01*invM_23 + invM_03*invM_12 + invM_02*invM_02)
                                     + pph1*(invM_01*invM_22 + invM_02*invM_12 + invM_11*invM_23 + invM_12*invM_13)
                                   + 2*pph2*(invM_02*invM_22 + invM_12*invM_23)
                                     + pph3*(invM_02*invM_23 + invM_03*invM_22 + invM_12*invM_33 + invM_13*invM_23))
                         + 4*ph_1*(    pph0*(invM_01*invM_22 + invM_02*invM_12)
                                     + pph1*(invM_11*invM_22 + invM_12*invM_12)
                                   + 2*pph2*invM_12*invM_22
                                     + pph3*(invM_12*invM_23 + invM_13*invM_22))
                         + 2*ph_2*(    pph0*(invM_00*invM_02 + 4*invM_02*invM_22 + invM_03*invM_23)
                                     + pph1*(invM_01*invM_02 + 4*invM_12*invM_22 + invM_13*invM_23)
                                     + pph2*(invM_02*invM_02 + 4*invM_22*invM_22 + invM_23*invM_23)
                                     + pph3*(invM_02*invM_03 + 4*invM_22*invM_23 + invM_23*invM_33))
                           + ph_3*(    pph0*(invM_02*invM_23 - invM_01*invM_02 - invM_00*invM_12 + invM_03*invM_22)
                                     + pph1*(invM_12*invM_23 - invM_02*invM_11 - invM_01*invM_12 + invM_13*invM_22)
                                   + 2*pph2*(invM_22*invM_23 - invM_02*invM_12)
                                     + pph3*(invM_22*invM_33 - invM_03*invM_12 - invM_02*invM_13 + invM_23*invM_23)));

   // ddotph2_dph3
   dDotFields[2][3] = r*(    ph_1*(    pph0*(invM_00*invM_22 + invM_01*invM_23 + invM_03*invM_12 + invM_02*invM_02)
                                   +   pph1*(invM_01*invM_22 + invM_02*invM_12 + invM_11*invM_23 + invM_12*invM_13)
                                   + 2*pph2*(invM_02*invM_22 + invM_12*invM_23)
                                   +   pph3*(invM_02*invM_23 + invM_03*invM_22 + invM_12*invM_33 + invM_13*invM_23))
                           + ph_2*(    pph0*(invM_02*invM_23 - invM_01*invM_02 - invM_00*invM_12 + invM_03*invM_22)
                                     + pph1*(invM_12*invM_23 - invM_02*invM_11 - invM_01*invM_12 + invM_13*invM_22)
                                   + 2*pph2*(invM_22*invM_23 - invM_02*invM_12)
                                     + pph3*(invM_22*invM_33 - invM_03*invM_12 - invM_02*invM_13 + invM_23*invM_23))
                         + 2*ph_3*(    pph0*(invM_01*invM_12 + invM_02*invM_22)
                                     + pph1*(invM_11*invM_12 + invM_12*invM_22)
                                     + pph2*(invM_12*invM_12 + invM_22*invM_22)
                                     + pph3*(invM_12*invM_13 + invM_22*invM_23)));

   // ddotph3_dph0
   dDotFields[3][0] = r*(  2*ph_0*(    pph0*(invM_01*invM_13 + invM_02*invM_23)
                                     + pph1*(invM_11*invM_13 + invM_12*invM_23)
                                     + pph2*(invM_12*invM_13 + invM_22*invM_23)
                                     + pph3*(invM_13*invM_13 + invM_23*invM_23))
                           + ph_1*(    pph0*(invM_00*invM_13 + invM_01*invM_03 - invM_02*invM_33 - invM_03*invM_23)
                                     + pph1*(invM_01*invM_13 + invM_11*invM_03 - invM_12*invM_33 - invM_13*invM_23)
                                     + pph2*(invM_02*invM_13 + invM_03*invM_12 - invM_22*invM_33 - invM_23*invM_23)
                                   + 2*pph3*(invM_03*invM_13 - invM_23*invM_33))
                           + ph_2*(    pph0*(invM_00*invM_23 + invM_02*invM_03 + invM_01*invM_33 + invM_03*invM_13)
                                     + pph1*(invM_01*invM_23 + invM_03*invM_12 + invM_11*invM_33 + invM_13*invM_13)
                                     + pph2*(invM_02*invM_23 + invM_03*invM_22 + invM_12*invM_33 + invM_13*invM_23)
                                   + 2*pph3*(invM_03*invM_23 + invM_13*invM_33)));

   // ddotph3_dph1
   dDotFields[3][1] = r*(    ph_3*(    pph0*(invM_00*invM_23 + invM_02*invM_03 + invM_01*invM_33 + invM_03*invM_13)
                                     + pph1*(invM_01*invM_23 + invM_03*invM_12 + invM_11*invM_33 + invM_13*invM_13)
                                     + pph2*(invM_02*invM_23 + invM_03*invM_22 + invM_12*invM_33 + invM_13*invM_23)
                                   + 2*pph3*(invM_03*invM_23 + invM_13*invM_33))
                         + 4*ph_2*(    pph0*(invM_01*invM_23 + invM_02*invM_13)
                                     + pph1*(invM_11*invM_23 + invM_12*invM_13)
                                     + pph2*(invM_12*invM_23 + invM_13*invM_22)
                                   + 2*pph3*invM_13*invM_23)
                           + ph_0*(    pph0*(invM_00*invM_13 + invM_01*invM_03 - invM_02*invM_33 - invM_03*invM_23)
                                     + pph1*(invM_01*invM_13 + invM_11*invM_03 - invM_12*invM_33 - invM_13*invM_23)
                                     + pph2*(invM_02*invM_13 + invM_03*invM_12 - invM_22*invM_33 - invM_23*invM_23)
                                   + 2*pph3*(invM_03*invM_13 - invM_23*invM_33))
                         + 2*ph_1*(    pph0*(invM_00*invM_03 + 4*invM_01*invM_13 + invM_03*invM_33)
                                     + pph1*(invM_01*invM_03 + 4*invM_11*invM_13 + invM_13*invM_33)
                                     + pph2*(invM_02*invM_03 + 4*invM_12*invM_13 + invM_23*invM_33)
                                     + pph3*(invM_03*invM_03 + 4*invM_13*invM_13 + invM_33*invM_33)));

   // ddotph3_dph2
   dDotFields[3][2] = r*(    ph_0*(    pph0*(invM_00*invM_23 + invM_02*invM_03 + invM_01*invM_33 + invM_03*invM_13)
                                     + pph1*(invM_01*invM_23 + invM_03*invM_12 + invM_11*invM_33 + invM_13*invM_13)
                                     + pph2*(invM_02*invM_23 + invM_03*invM_22 + invM_12*invM_33 + invM_13*invM_23)
                                     + pph3*(2*invM_03*invM_23 + 2*invM_13*invM_33))
                         + 4*ph_1*(    pph0*(invM_01*invM_23 + invM_02*invM_13)
                                     + pph1*(invM_11*invM_23 + invM_12*invM_13)
                                     + pph2*(invM_12*invM_23 + invM_13*invM_22)
                                   + 2*pph3*invM_13*invM_23)
                         + 2*ph_2*(    pph0*(invM_00*invM_03 + 4*invM_02*invM_23 + invM_03*invM_33)
                                     + pph1*(invM_01*invM_03 + 4*invM_12*invM_23 + invM_13*invM_33)
                                     + pph2*(invM_02*invM_03 + 4*invM_22*invM_23 + invM_23*invM_33)
                                     + pph3*(invM_03*invM_03 + 4*invM_23*invM_23 + invM_33*invM_33))
                           + ph_3*(    pph0*(invM_02*invM_33 - invM_01*invM_03 - invM_00*invM_13 + invM_03*invM_23)
                                     + pph1*(invM_12*invM_33 - invM_11*invM_03 - invM_01*invM_13 + invM_13*invM_23)
                                     + pph2*(invM_22*invM_33 - invM_03*invM_12 - invM_02*invM_13 + invM_23*invM_23)
                                   + 2*pph3*(invM_23*invM_33 - invM_03*invM_13)));

   // ddotph3_dph3
   dDotFields[3][3] = r*(    ph_1*(    pph0*(invM_00*invM_23 + invM_02*invM_03 + invM_01*invM_33 + invM_03*invM_13)
                                     + pph1*(invM_01*invM_23 + invM_03*invM_12 + invM_11*invM_33 + invM_13*invM_13)
                                     + pph2*(invM_02*invM_23 + invM_03*invM_22 + invM_12*invM_33 + invM_13*invM_23)
                                   + 2*pph3*(invM_03*invM_23 + invM_13*invM_33))
                           + ph_2*(    pph0*(invM_02*invM_33 - invM_01*invM_03 - invM_00*invM_13 + invM_03*invM_23)
                                     + pph1*(invM_12*invM_33 - invM_11*invM_03 - invM_01*invM_13 + invM_13*invM_23)
                                     + pph2*(invM_22*invM_33 - invM_03*invM_12 - invM_02*invM_13 + invM_23*invM_23)
                                   + 2*pph3*(invM_23*invM_33 - invM_03*invM_13))
                         + 2*ph_3*(    pph0*(invM_01*invM_13 + invM_02*invM_23)
                                     + pph1*(invM_11*invM_13 + invM_12*invM_23)
                                     + pph2*(invM_12*invM_13 + invM_22*invM_23)
                                     + pph3*(invM_13*invM_13 + invM_23*invM_23)));

   r = -8*gR2*(  4*ph_1*ph_2*(2*pph0*invM_01*invM_02 + pph1*(invM_01*invM_12 + invM_02*invM_11) + pph2*(invM_01*invM_22 + invM_02*invM_12) + pph3*(invM_01*invM_23 + invM_02*invM_13))
               + (ph_0*ph_0+ph_3*ph_3)*(pph1*(invM_01*invM_11 + invM_02*invM_12) + pph2*(invM_01*invM_12 + invM_02*invM_22) + pph3*(invM_01*invM_13 + invM_02*invM_23) + pph0*(invM_01*invM_01 + invM_02*invM_02))
               + (ph_2*ph_3-ph_0*ph_1)*(2*pph0*(invM_02*invM_03 - invM_00*invM_01) + pph2*(invM_02*invM_23 - invM_01*invM_02 - invM_00*invM_12 + invM_03*invM_22) + pph3*(invM_02*invM_33 - invM_01*invM_03 - invM_00*invM_13 + invM_03*invM_23) + pph1*(invM_02*invM_13 - invM_00*invM_11 + invM_03*invM_12 - invM_01*invM_01))
               + (ph_0*ph_2+ph_1*ph_3)*(pph1*(invM_00*invM_12 + invM_01*invM_02 + invM_01*invM_13 + invM_11*invM_03) + pph3*(invM_00*invM_23 + invM_02*invM_03 + invM_01*invM_33 + invM_03*invM_13) + pph2*(invM_00*invM_22 + invM_01*invM_23 + invM_03*invM_12 + invM_02*invM_02) + pph0*(2*invM_00*invM_02 + 2*invM_01*invM_03))
               + ph_1*ph_1*(pph1*(invM_00*invM_01 + 4*invM_01*invM_11 + invM_03*invM_13) + pph2*(invM_00*invM_02 + 4*invM_01*invM_12 + invM_03*invM_23) + pph3*(invM_00*invM_03 + 4*invM_01*invM_13 + invM_03*invM_33) + pph0*(invM_00*invM_00 + 4*invM_01*invM_01 + invM_03*invM_03))
               + ph_2*ph_2*(pph1*(invM_00*invM_01 + 4*invM_02*invM_12 + invM_03*invM_13) + pph2*(invM_00*invM_02 + 4*invM_02*invM_22 + invM_03*invM_23) + pph3*(invM_00*invM_03 + 4*invM_02*invM_23 + invM_03*invM_33) + pph0*(invM_00*invM_00 + 4*invM_02*invM_02 + invM_03*invM_03)));

   // ddotph0_dA1 
   dDotFields[0][4] = A_1*r;

   // ddotph0_dA2 
   dDotFields[0][5] = A_2*r;

   // ddotph0_dA3 
   dDotFields[0][6] = A_3*r;
   
   r = -8*gR2*(  (ph_0*ph_0 + ph_3*ph_3)*(pph0*(invM_01*invM_11 + invM_02*invM_12) + pph2*(invM_11*invM_12 + invM_12*invM_22) + pph3*(invM_11*invM_13 + invM_12*invM_23) + pph1*(invM_11*invM_11 + invM_12*invM_12))
               + (ph_0*ph_2 + ph_1*ph_3)*(pph0*(invM_00*invM_12 + invM_01*invM_02 + invM_01*invM_13 + invM_11*invM_03) + pph2*(invM_01*invM_22 + invM_02*invM_12 + invM_11*invM_23 + invM_12*invM_13) + pph3*(invM_01*invM_23 + invM_03*invM_12 + invM_11*invM_33 + invM_13*invM_13) + pph1*(2*invM_01*invM_12 + 2*invM_11*invM_13))
               + 4*ph_1*ph_2*(2*pph1*invM_11*invM_12 + pph0*(invM_01*invM_12 + invM_02*invM_11) + pph3*(invM_11*invM_23 + invM_12*invM_13) + pph2*(invM_11*invM_22 + invM_12*invM_12))
               + ph_1*ph_1*(pph0*(invM_00*invM_01 + 4*invM_01*invM_11 + invM_03*invM_13) + pph2*(invM_01*invM_02 + 4*invM_11*invM_12 + invM_13*invM_23) + pph3*(invM_01*invM_03 + 4*invM_11*invM_13 + invM_13*invM_33) + pph1*(invM_01*invM_01 + 4*invM_11*invM_11 + invM_13*invM_13))
               + ph_2*ph_2*(pph0*(invM_00*invM_01 + 4*invM_02*invM_12 + invM_03*invM_13) + pph2*(invM_01*invM_02 + 4*invM_12*invM_22 + invM_13*invM_23) + pph3*(invM_01*invM_03 + 4*invM_12*invM_23 + invM_13*invM_33) + pph1*(invM_01*invM_01 + 4*invM_12*invM_12 + invM_13*invM_13)) 
               + ph_0*ph_1*(pph1*(2*invM_01*invM_11 - 2*invM_12*invM_13) + pph2*(invM_01*invM_12 + invM_02*invM_11 - invM_12*invM_23 - invM_13*invM_22) + pph3*(invM_01*invM_13 + invM_11*invM_03 - invM_12*invM_33 - invM_13*invM_23) + pph0*(invM_00*invM_11 - invM_02*invM_13 - invM_03*invM_12 + invM_01*invM_01)) 
               + ph_2*ph_3*(pph1*(2*invM_12*invM_13 - 2*invM_01*invM_11) + pph2*(invM_12*invM_23 - invM_02*invM_11 - invM_01*invM_12 + invM_13*invM_22) + pph3*(invM_12*invM_33 - invM_11*invM_03 - invM_01*invM_13 + invM_13*invM_23) + pph0*(invM_02*invM_13 - invM_00*invM_11 + invM_03*invM_12 - invM_01*invM_01)));

   // ddotph1_dA1 
   dDotFields[1][4] = A_1*r;

   // ddotph1_dA2 
   dDotFields[1][5] = A_2*r;

   // ddotph1_dA3 
   dDotFields[1][6] = A_3*r;
   
   r = -8*gR2*(  (ph_0*ph_0 + ph_3*ph_3)*(pph0*(invM_01*invM_12 + invM_02*invM_22) + pph1*(invM_11*invM_12 + invM_12*invM_22) + pph3*(invM_12*invM_13 + invM_22*invM_23) + pph2*(invM_12*invM_12 + invM_22*invM_22)) 
               + (ph_0*ph_2 + ph_1*ph_3)*(pph1*(invM_01*invM_22 + invM_02*invM_12 + invM_11*invM_23 + invM_12*invM_13) + pph3*(invM_02*invM_23 + invM_03*invM_22 + invM_12*invM_33 + invM_13*invM_23) + pph0*(invM_00*invM_22 + invM_01*invM_23 + invM_03*invM_12 + invM_02*invM_02) + pph2*(2*invM_02*invM_22 + 2*invM_12*invM_23))
               + 4*ph_1*ph_2*(2*pph2*invM_12*invM_22 + pph0*(invM_01*invM_22 + invM_02*invM_12) + pph3*(invM_12*invM_23 + invM_13*invM_22) + pph1*(invM_11*invM_22 + invM_12*invM_12))
               + ph_1*ph_1*(pph0*(invM_00*invM_02 + 4*invM_01*invM_12 + invM_03*invM_23) + pph1*(invM_01*invM_02 + 4*invM_11*invM_12 + invM_13*invM_23) + pph3*(invM_02*invM_03 + 4*invM_12*invM_13 + invM_23*invM_33) + pph2*(invM_02*invM_02 + 4*invM_12*invM_12 + invM_23*invM_23))
               + ph_2*ph_2*(pph0*(invM_00*invM_02 + 4*invM_02*invM_22 + invM_03*invM_23) + pph1*(invM_01*invM_02 + 4*invM_12*invM_22 + invM_13*invM_23) + pph3*(invM_02*invM_03 + 4*invM_22*invM_23 + invM_23*invM_33) + pph2*(invM_02*invM_02 + 4*invM_22*invM_22 + invM_23*invM_23)) 
               + ph_0*ph_1*(pph2*(2*invM_02*invM_12 - 2*invM_22*invM_23) + pph0*(invM_00*invM_12 + invM_01*invM_02 - invM_02*invM_23 - invM_03*invM_22) + pph1*(invM_01*invM_12 + invM_02*invM_11 - invM_12*invM_23 - invM_13*invM_22) + pph3*(invM_02*invM_13 + invM_03*invM_12 - invM_22*invM_33 - invM_23*invM_23)) 
               + ph_2*ph_3*(pph2*(2*invM_22*invM_23 - 2*invM_02*invM_12) + pph0*(invM_02*invM_23 - invM_01*invM_02 - invM_00*invM_12 + invM_03*invM_22) + pph1*(invM_12*invM_23 - invM_02*invM_11 - invM_01*invM_12 + invM_13*invM_22) + pph3*(invM_22*invM_33 - invM_03*invM_12 - invM_02*invM_13 + invM_23*invM_23)));

   // ddotph2_dA1 
   dDotFields[2][4] = A_1*r;

   // ddotph2_dA2
   dDotFields[2][5] = A_2*r;

   // ddotph2_dA3
   dDotFields[2][6] = A_3*r;
   
   r = -8*gR2*(  (ph_0*ph_0+ ph_3*ph_3)*(pph0*(invM_01*invM_13 + invM_02*invM_23) + pph1*(invM_11*invM_13 + invM_12*invM_23) + pph2*(invM_12*invM_13 + invM_22*invM_23) + pph3*(invM_13*invM_13 + invM_23*invM_23))
               + (ph_0*ph_2 + ph_1*ph_3)*(pph0*(invM_00*invM_23 + invM_02*invM_03 + invM_01*invM_33 + invM_03*invM_13) + pph2*(invM_02*invM_23 + invM_03*invM_22 + invM_12*invM_33 + invM_13*invM_23) + pph1*(invM_01*invM_23 + invM_03*invM_12 + invM_11*invM_33 + invM_13*invM_13) + 2*pph3*(invM_03*invM_23 + invM_13*invM_33))
               + 4*ph_1*ph_2*(pph0*invM_01*invM_23 + pph0*invM_02*invM_13 + pph1*invM_11*invM_23 + pph1*invM_12*invM_13 + pph2*invM_12*invM_23 + pph2*invM_13*invM_22 + 2*pph3*invM_13*invM_23)
               + ph_1*ph_1*(pph0*(invM_00*invM_03 + 4*invM_01*invM_13 + invM_03*invM_33) + pph1*(invM_01*invM_03 + 4*invM_11*invM_13 + invM_13*invM_33) + pph2*(invM_02*invM_03 + 4*invM_12*invM_13 + invM_23*invM_33) + pph3*(invM_03*invM_03 + 4*invM_13*invM_13 + invM_33*invM_33))
               + ph_2*ph_2*(pph0*(invM_00*invM_03 + 4*invM_02*invM_23 + invM_03*invM_33) + pph1*(invM_01*invM_03 + 4*invM_12*invM_23 + invM_13*invM_33) + pph2*(invM_02*invM_03 + 4*invM_22*invM_23 + invM_23*invM_33) + pph3*(invM_03*invM_03 + 4*invM_23*invM_23 + invM_33*invM_33))
               + ph_0*ph_1*(pph3*(2*invM_03*invM_13 - 2*invM_23*invM_33) + pph0*(invM_00*invM_13 + invM_01*invM_03 - invM_02*invM_33 - invM_03*invM_23) + pph1*(invM_01*invM_13 + invM_11*invM_03 - invM_12*invM_33 - invM_13*invM_23) + pph2*(invM_02*invM_13 + invM_03*invM_12 - invM_22*invM_33 - invM_23*invM_23))
               + ph_2*ph_3*(pph3*(2*invM_23*invM_33 - 2*invM_03*invM_13) + pph0*(invM_02*invM_33 - invM_01*invM_03 - invM_00*invM_13 + invM_03*invM_23) + pph1*(invM_12*invM_33 - invM_11*invM_03 - invM_01*invM_13 + invM_13*invM_23) + pph2*(invM_22*invM_33 - invM_03*invM_12 - invM_02*invM_13 + invM_23*invM_23)));

   // ddotph3_dA1
   dDotFields[3][4] = A_1*r;

   // ddotph3_dA2
   dDotFields[3][5] = A_2*r;

   // ddotph3_dA3
   dDotFields[3][6] = A_3*r;
}

#define ph_0(x) (fields[x].ph_0)
#define ph_1(x) (fields[x].ph_1)
#define ph_2(x) (fields[x].ph_2)
#define ph_3(x) (fields[x].ph_3)
#define A_1(x) (fields[x].A_1)
#define A_2(x) (fields[x].A_2)
#define A_3(x) (fields[x].A_3)

#define d1ph_0(x) (aux[x].dx_ph_0)
#define d1ph_1(x) (aux[x].dx_ph_1)
#define d1ph_2(x) (aux[x].dx_ph_2)
#define d1ph_3(x) (aux[x].dx_ph_3)
#define d1A_2(x) (aux[x].dx_A_2)
#define d1A_3(x) (aux[x].dx_A_3)

void TLattice::compute_staticDotMomenta(int i, TFields fields[MAX_NSIDE], TFields &dotMomenta)
/* IMPORTANT: Requires field derivatives to be ready for use! */
{
   int iM1, iM2, iP1, iP2;

   if (0 == i) iM1 = nSideM1; else iM1 = i - 1;

   if (0 == iM1) iM2 = nSideM1; else iM2 = iM1 - 1;

   if (nSideM1 == i) iP1 = 0; else iP1 = i + 1;

   if (nSideM1 == iP1) iP2 = 0; else iP2 = iP1 + 1;

   double aR2 = a*a;
   
   dotMomenta.ph_0 = a*(d1ph_0(iM1)*(-2.0/3.0) + d1ph_0(iP1)*(2.0/3.0) + d1ph_0(iM2)*(1.0/12.0) - d1ph_0(iP2)*(1.0/12.0) + ((-
    1.0/3.0)*gR2*(d1ph_1(i)*d1ph_1(i))*(12.0*(A_2(i)*A_2(i))*ph_0(i) + 12.0*(A_3(i)*A_3(i))*ph_0(i)) - (1.0/3.0)*gR2*(d1ph_2(i)*
    d1ph_2(i))*(12.0*(A_2(i)*A_2(i))*ph_0(i) + 12.0*(A_3(i)*A_3(i))*ph_0(i)) - (1.0/3.0)*gR2*d1ph_0(i)*d1ph_1(i)*(12.0*(A_2(i)*
    A_2(i))*ph_1(i) + 12.0*(A_3(i)*A_3(i))*ph_1(i)) - (1.0/3.0)*gR2*d1ph_0(i)*d1ph_2(i)*(12.0*(A_2(i)*A_2(i))*ph_2(i) + 12.0*(A_3
    (i)*A_3(i))*ph_2(i)) - (1.0/3.0)*gR2*d1ph_2(i)*d1ph_3(i)*(-12.0*(A_2(i)*A_2(i))*ph_1(i) - 12.0*(A_3(i)*A_3(i))*ph_1(i)) - (
    1.0/3.0)*gR2*d1ph_1(i)*d1ph_3(i)*(12.0*(A_2(i)*A_2(i))*ph_2(i) + 12.0*(A_3(i)*A_3(i))*ph_2(i)) - (1.0/3.0)*gR2*d1ph_0(iP2)*((
    A_2(iP2)*A_2(iP2))*(ph_1(iP2)*ph_1(iP2)) + (A_2(iP2)*A_2(iP2))*(ph_2(iP2)*ph_2(iP2)) + (A_3(iP2)*A_3(iP2))*(ph_1(iP2)*ph_1(
    iP2)) + (A_3(iP2)*A_3(iP2))*(ph_2(iP2)*ph_2(iP2))) - (1.0/3.0)*gR2*d1ph_0(iM2)*(-((A_2(iM2)*A_2(iM2))*(ph_1(iM2)*ph_1(iM2))) 
    - (A_2(iM2)*A_2(iM2))*(ph_2(iM2)*ph_2(iM2)) - (A_3(iM2)*A_3(iM2))*(ph_1(iM2)*ph_1(iM2)) - (A_3(iM2)*A_3(iM2))*(ph_2(iM2)*ph_2
    (iM2))) - (1.0/3.0)*gR2*d1ph_0(iM1)*(8.0*(A_2(iM1)*A_2(iM1))*(ph_1(iM1)*ph_1(iM1)) + 8.0*(A_2(iM1)*A_2(iM1))*(ph_2(iM1)*ph_2(
    iM1)) + 8.0*(A_3(iM1)*A_3(iM1))*(ph_1(iM1)*ph_1(iM1)) + 8.0*(A_3(iM1)*A_3(iM1))*(ph_2(iM1)*ph_2(iM1))) - (1.0/3.0)*gR2*d1ph_0
    (iP1)*(-8.0*(A_2(iP1)*A_2(iP1))*(ph_1(iP1)*ph_1(iP1)) - 8.0*(A_2(iP1)*A_2(iP1))*(ph_2(iP1)*ph_2(iP1)) - 8.0*(A_3(iP1)*A_3(iP1
    ))*(ph_1(iP1)*ph_1(iP1)) - 8.0*(A_3(iP1)*A_3(iP1))*(ph_2(iP1)*ph_2(iP1))) - (1.0/3.0)*gR2*d1ph_2(iP2)*((A_2(iP2)*A_2(iP2))*
    ph_0(iP2)*ph_2(iP2) + (A_3(iP2)*A_3(iP2))*ph_0(iP2)*ph_2(iP2) + (A_2(iP2)*A_2(iP2))*ph_1(iP2)*ph_3(iP2) + (A_3(iP2)*A_3(iP2
    ))*ph_1(iP2)*ph_3(iP2)) - (1.0/3.0)*gR2*d1ph_1(iM2)*(-((A_2(iM2)*A_2(iM2))*ph_0(iM2)*ph_1(iM2)) - (A_3(iM2)*A_3(iM2))*ph_0(
    iM2)*ph_1(iM2) + (A_2(iM2)*A_2(iM2))*ph_2(iM2)*ph_3(iM2) + (A_3(iM2)*A_3(iM2))*ph_2(iM2)*ph_3(iM2)) - (1.0/3.0)*gR2*d1ph_1(
    iP2)*((A_2(iP2)*A_2(iP2))*ph_0(iP2)*ph_1(iP2) + (A_3(iP2)*A_3(iP2))*ph_0(iP2)*ph_1(iP2) - (A_2(iP2)*A_2(iP2))*ph_2(iP2)*ph_3(
    iP2) - (A_3(iP2)*A_3(iP2))*ph_2(iP2)*ph_3(iP2)) - (1.0/3.0)*gR2*d1ph_2(iM2)*(-((A_2(iM2)*A_2(iM2))*ph_0(iM2)*ph_2(iM2)) - (
    A_3(iM2)*A_3(iM2))*ph_0(iM2)*ph_2(iM2) - (A_2(iM2)*A_2(iM2))*ph_1(iM2)*ph_3(iM2) - (A_3(iM2)*A_3(iM2))*ph_1(iM2)*ph_3(iM2)) -
     (1.0/3.0)*gR2*d1ph_1(iM1)*(8.0*(A_2(iM1)*A_2(iM1))*ph_0(iM1)*ph_1(iM1) + 8.0*(A_3(iM1)*A_3(iM1))*ph_0(iM1)*ph_1(iM1) - 8.0*(
    A_2(iM1)*A_2(iM1))*ph_2(iM1)*ph_3(iM1) - 8.0*(A_3(iM1)*A_3(iM1))*ph_2(iM1)*ph_3(iM1)) - (1.0/3.0)*gR2*d1ph_1(iP1)*(-8.0*(A_2(
    iP1)*A_2(iP1))*ph_0(iP1)*ph_1(iP1) - 8.0*(A_3(iP1)*A_3(iP1))*ph_0(iP1)*ph_1(iP1) + 8.0*(A_2(iP1)*A_2(iP1))*ph_2(iP1)*ph_3(iP1
    ) + 8.0*(A_3(iP1)*A_3(iP1))*ph_2(iP1)*ph_3(iP1)) - (1.0/3.0)*gR2*d1ph_2(iM1)*(8.0*(A_2(iM1)*A_2(iM1))*ph_0(iM1)*ph_2(iM1) + 
    8.0*(A_3(iM1)*A_3(iM1))*ph_0(iM1)*ph_2(iM1) + 8.0*(A_2(iM1)*A_2(iM1))*ph_1(iM1)*ph_3(iM1) + 8.0*(A_3(iM1)*A_3(iM1))*ph_1(iM1
    )*ph_3(iM1)) - (1.0/3.0)*gR2*d1ph_2(iP1)*(-8.0*(A_2(iP1)*A_2(iP1))*ph_0(iP1)*ph_2(iP1) - 8.0*(A_3(iP1)*A_3(iP1))*ph_0(iP1)*
    ph_2(iP1) - 8.0*(A_2(iP1)*A_2(iP1))*ph_1(iP1)*ph_3(iP1) - 8.0*(A_3(iP1)*A_3(iP1))*ph_1(iP1)*ph_3(iP1)))/aR2 - 2.0*(a*a)*
    lambdaR2*ph_0(i)*(ph_0(i)*ph_0(i) + ph_1(i)*ph_1(i) + ph_2(i)*ph_2(i) + ph_3(i)*ph_3(i) - 1.0));

   dotMomenta.ph_1 = a*(d1ph_1(iM1)*(-2.0/3.0) + d1ph_1(iP1)*(2.0/3.0) + d1ph_1(iM2)*(1.0/12.0) - d1ph_1(iP2)*(1.0/12.0) + ((-
    1.0/3.0)*gR2*(d1ph_0(i)*d1ph_0(i))*(12.0*(A_2(i)*A_2(i))*ph_1(i) + 12.0*(A_3(i)*A_3(i))*ph_1(i)) - (1.0/3.0)*gR2*(d1ph_3(i)*
    d1ph_3(i))*(12.0*(A_2(i)*A_2(i))*ph_1(i) + 12.0*(A_3(i)*A_3(i))*ph_1(i)) - (1.0/3.0)*gR2*(d1ph_1(i)*d1ph_1(i))*(48.0*(A_2(i)*
    A_2(i))*ph_1(i) + 48.0*(A_3(i)*A_3(i))*ph_1(i)) - (1.0/3.0)*gR2*d1ph_0(i)*d1ph_1(i)*(12.0*(A_2(i)*A_2(i))*ph_0(i) + 12.0*(A_3
    (i)*A_3(i))*ph_0(i)) - (1.0/3.0)*gR2*d1ph_2(i)*d1ph_3(i)*(-12.0*(A_2(i)*A_2(i))*ph_0(i) - 12.0*(A_3(i)*A_3(i))*ph_0(i)) - (
    1.0/3.0)*gR2*d1ph_0(i)*d1ph_2(i)*(12.0*(A_2(i)*A_2(i))*ph_3(i) + 12.0*(A_3(i)*A_3(i))*ph_3(i)) - (1.0/3.0)*gR2*d1ph_1(i)*
    d1ph_3(i)*(12.0*(A_2(i)*A_2(i))*ph_3(i) + 12.0*(A_3(i)*A_3(i))*ph_3(i)) - (1.0/3.0)*gR2*d1ph_1(i)*d1ph_2(i)*(48.0*(A_2(i)*A_2
    (i))*ph_2(i) + 48.0*(A_3(i)*A_3(i))*ph_2(i)) - (1.0/3.0)*gR2*d1ph_2(iM2)*(-4.0*(A_2(iM2)*A_2(iM2))*ph_1(iM2)*ph_2(iM2) - 4.0
    *(A_3(iM2)*A_3(iM2))*ph_1(iM2)*ph_2(iM2)) - (1.0/3.0)*gR2*d1ph_2(iP2)*(4.0*(A_2(iP2)*A_2(iP2))*ph_1(iP2)*ph_2(iP2) + 4.0*(A_3
    (iP2)*A_3(iP2))*ph_1(iP2)*ph_2(iP2)) - (1.0/3.0)*gR2*d1ph_2(iM1)*(32.0*(A_2(iM1)*A_2(iM1))*ph_1(iM1)*ph_2(iM1) + 32.0*(A_3(
    iM1)*A_3(iM1))*ph_1(iM1)*ph_2(iM1)) - (1.0/3.0)*gR2*d1ph_2(iP1)*(-32.0*(A_2(iP1)*A_2(iP1))*ph_1(iP1)*ph_2(iP1) - 32.0*(A_3(
    iP1)*A_3(iP1))*ph_1(iP1)*ph_2(iP1)) - (1.0/3.0)*gR2*d1ph_3(iP2)*((A_2(iP2)*A_2(iP2))*ph_0(iP2)*ph_2(iP2) + (A_3(iP2)*A_3(iP2
    ))*ph_0(iP2)*ph_2(iP2) + (A_2(iP2)*A_2(iP2))*ph_1(iP2)*ph_3(iP2) + (A_3(iP2)*A_3(iP2))*ph_1(iP2)*ph_3(iP2)) - (1.0/3.0)*gR2*
    d1ph_0(iM2)*(-((A_2(iM2)*A_2(iM2))*ph_0(iM2)*ph_1(iM2)) - (A_3(iM2)*A_3(iM2))*ph_0(iM2)*ph_1(iM2) + (A_2(iM2)*A_2(iM2))*ph_2(
    iM2)*ph_3(iM2) + (A_3(iM2)*A_3(iM2))*ph_2(iM2)*ph_3(iM2)) - (1.0/3.0)*gR2*d1ph_0(iP2)*((A_2(iP2)*A_2(iP2))*ph_0(iP2)*ph_1(iP2
    ) + (A_3(iP2)*A_3(iP2))*ph_0(iP2)*ph_1(iP2) - (A_2(iP2)*A_2(iP2))*ph_2(iP2)*ph_3(iP2) - (A_3(iP2)*A_3(iP2))*ph_2(iP2)*ph_3(
    iP2)) - (1.0/3.0)*gR2*d1ph_3(iM2)*(-((A_2(iM2)*A_2(iM2))*ph_0(iM2)*ph_2(iM2)) - (A_3(iM2)*A_3(iM2))*ph_0(iM2)*ph_2(iM2) - (
    A_2(iM2)*A_2(iM2))*ph_1(iM2)*ph_3(iM2) - (A_3(iM2)*A_3(iM2))*ph_1(iM2)*ph_3(iM2)) - (1.0/3.0)*gR2*d1ph_0(iM1)*(8.0*(A_2(iM1)*
    A_2(iM1))*ph_0(iM1)*ph_1(iM1) + 8.0*(A_3(iM1)*A_3(iM1))*ph_0(iM1)*ph_1(iM1) - 8.0*(A_2(iM1)*A_2(iM1))*ph_2(iM1)*ph_3(iM1) - 
    8.0*(A_3(iM1)*A_3(iM1))*ph_2(iM1)*ph_3(iM1)) - (1.0/3.0)*gR2*d1ph_0(iP1)*(-8.0*(A_2(iP1)*A_2(iP1))*ph_0(iP1)*ph_1(iP1) - 8.0
    *(A_3(iP1)*A_3(iP1))*ph_0(iP1)*ph_1(iP1) + 8.0*(A_2(iP1)*A_2(iP1))*ph_2(iP1)*ph_3(iP1) + 8.0*(A_3(iP1)*A_3(iP1))*ph_2(iP1)*
    ph_3(iP1)) - (1.0/3.0)*gR2*d1ph_3(iM1)*(8.0*(A_2(iM1)*A_2(iM1))*ph_0(iM1)*ph_2(iM1) + 8.0*(A_3(iM1)*A_3(iM1))*ph_0(iM1)*ph_2(
    iM1) + 8.0*(A_2(iM1)*A_2(iM1))*ph_1(iM1)*ph_3(iM1) + 8.0*(A_3(iM1)*A_3(iM1))*ph_1(iM1)*ph_3(iM1)) - (1.0/3.0)*gR2*d1ph_3(iP1
    )*(-8.0*(A_2(iP1)*A_2(iP1))*ph_0(iP1)*ph_2(iP1) - 8.0*(A_3(iP1)*A_3(iP1))*ph_0(iP1)*ph_2(iP1) - 8.0*(A_2(iP1)*A_2(iP1))*ph_1(
    iP1)*ph_3(iP1) - 8.0*(A_3(iP1)*A_3(iP1))*ph_1(iP1)*ph_3(iP1)) - (1.0/3.0)*gR2*d1ph_1(iP2)*((A_2(iP2)*A_2(iP2))*(ph_0(iP2)*
    ph_0(iP2)) + 4.0*(A_2(iP2)*A_2(iP2))*(ph_1(iP2)*ph_1(iP2)) + (A_3(iP2)*A_3(iP2))*(ph_0(iP2)*ph_0(iP2)) + 4.0*(A_3(iP2)*A_3(
    iP2))*(ph_1(iP2)*ph_1(iP2)) + (A_2(iP2)*A_2(iP2))*(ph_3(iP2)*ph_3(iP2)) + (A_3(iP2)*A_3(iP2))*(ph_3(iP2)*ph_3(iP2))) - (1.0/
    3.0)*gR2*d1ph_1(iM2)*(-((A_2(iM2)*A_2(iM2))*(ph_0(iM2)*ph_0(iM2))) - 4.0*(A_2(iM2)*A_2(iM2))*(ph_1(iM2)*ph_1(iM2)) - (A_3(iM2
    )*A_3(iM2))*(ph_0(iM2)*ph_0(iM2)) - 4.0*(A_3(iM2)*A_3(iM2))*(ph_1(iM2)*ph_1(iM2)) - (A_2(iM2)*A_2(iM2))*(ph_3(iM2)*ph_3(iM2))
     - (A_3(iM2)*A_3(iM2))*(ph_3(iM2)*ph_3(iM2))) - (1.0/3.0)*gR2*d1ph_1(iM1)*(8.0*(A_2(iM1)*A_2(iM1))*(ph_0(iM1)*ph_0(iM1)) + 
    32.0*(A_2(iM1)*A_2(iM1))*(ph_1(iM1)*ph_1(iM1)) + 8.0*(A_3(iM1)*A_3(iM1))*(ph_0(iM1)*ph_0(iM1)) + 32.0*(A_3(iM1)*A_3(iM1))*(
    ph_1(iM1)*ph_1(iM1)) + 8.0*(A_2(iM1)*A_2(iM1))*(ph_3(iM1)*ph_3(iM1)) + 8.0*(A_3(iM1)*A_3(iM1))*(ph_3(iM1)*ph_3(iM1))) - (1.0/
    3.0)*gR2*d1ph_1(iP1)*(-8.0*(A_2(iP1)*A_2(iP1))*(ph_0(iP1)*ph_0(iP1)) - 32.0*(A_2(iP1)*A_2(iP1))*(ph_1(iP1)*ph_1(iP1)) - 8.0*(
    A_3(iP1)*A_3(iP1))*(ph_0(iP1)*ph_0(iP1)) - 32.0*(A_3(iP1)*A_3(iP1))*(ph_1(iP1)*ph_1(iP1)) - 8.0*(A_2(iP1)*A_2(iP1))*(ph_3(iP1
    )*ph_3(iP1)) - 8.0*(A_3(iP1)*A_3(iP1))*(ph_3(iP1)*ph_3(iP1))))/aR2 - 2.0*(a*a)*lambdaR2*ph_1(i)*(ph_0(i)*ph_0(i) + ph_1(i)*
    ph_1(i) + ph_2(i)*ph_2(i) + ph_3(i)*ph_3(i) - 1.0));

   dotMomenta.ph_2 = a*(d1ph_2(iM1)*(-2.0/3.0) + d1ph_2(iP1)*(2.0/3.0) + d1ph_2(iM2)*(1.0/12.0) - d1ph_2(iP2)*(1.0/12.0) + ((-
    1.0/3.0)*gR2*(d1ph_0(i)*d1ph_0(i))*(12.0*(A_2(i)*A_2(i))*ph_2(i) + 12.0*(A_3(i)*A_3(i))*ph_2(i)) - (1.0/3.0)*gR2*(d1ph_3(i)*
    d1ph_3(i))*(12.0*(A_2(i)*A_2(i))*ph_2(i) + 12.0*(A_3(i)*A_3(i))*ph_2(i)) - (1.0/3.0)*gR2*(d1ph_2(i)*d1ph_2(i))*(48.0*(A_2(i)*
    A_2(i))*ph_2(i) + 48.0*(A_3(i)*A_3(i))*ph_2(i)) - (1.0/3.0)*gR2*d1ph_0(i)*d1ph_2(i)*(12.0*(A_2(i)*A_2(i))*ph_0(i) + 12.0*(A_3
    (i)*A_3(i))*ph_0(i)) - (1.0/3.0)*gR2*d1ph_1(i)*d1ph_3(i)*(12.0*(A_2(i)*A_2(i))*ph_0(i) + 12.0*(A_3(i)*A_3(i))*ph_0(i)) - (1.0
    /3.0)*gR2*d1ph_0(i)*d1ph_1(i)*(-12.0*(A_2(i)*A_2(i))*ph_3(i) - 12.0*(A_3(i)*A_3(i))*ph_3(i)) - (1.0/3.0)*gR2*d1ph_2(i)*d1ph_3
    (i)*(12.0*(A_2(i)*A_2(i))*ph_3(i) + 12.0*(A_3(i)*A_3(i))*ph_3(i)) - (1.0/3.0)*gR2*d1ph_1(i)*d1ph_2(i)*(48.0*(A_2(i)*A_2(i))*
    ph_1(i) + 48.0*(A_3(i)*A_3(i))*ph_1(i)) - (1.0/3.0)*gR2*d1ph_1(iM2)*(-4.0*(A_2(iM2)*A_2(iM2))*ph_1(iM2)*ph_2(iM2) - 4.0*(A_3(
    iM2)*A_3(iM2))*ph_1(iM2)*ph_2(iM2)) - (1.0/3.0)*gR2*d1ph_1(iP2)*(4.0*(A_2(iP2)*A_2(iP2))*ph_1(iP2)*ph_2(iP2) + 4.0*(A_3(iP2)*
    A_3(iP2))*ph_1(iP2)*ph_2(iP2)) - (1.0/3.0)*gR2*d1ph_1(iM1)*(32.0*(A_2(iM1)*A_2(iM1))*ph_1(iM1)*ph_2(iM1) + 32.0*(A_3(iM1)*A_3
    (iM1))*ph_1(iM1)*ph_2(iM1)) - (1.0/3.0)*gR2*d1ph_1(iP1)*(-32.0*(A_2(iP1)*A_2(iP1))*ph_1(iP1)*ph_2(iP1) - 32.0*(A_3(iP1)*A_3(
    iP1))*ph_1(iP1)*ph_2(iP1)) - (1.0/3.0)*gR2*d1ph_0(iP2)*((A_2(iP2)*A_2(iP2))*ph_0(iP2)*ph_2(iP2) + (A_3(iP2)*A_3(iP2))*ph_0(
    iP2)*ph_2(iP2) + (A_2(iP2)*A_2(iP2))*ph_1(iP2)*ph_3(iP2) + (A_3(iP2)*A_3(iP2))*ph_1(iP2)*ph_3(iP2)) - (1.0/3.0)*gR2*d1ph_3(
    iM2)*((A_2(iM2)*A_2(iM2))*ph_0(iM2)*ph_1(iM2) + (A_3(iM2)*A_3(iM2))*ph_0(iM2)*ph_1(iM2) - (A_2(iM2)*A_2(iM2))*ph_2(iM2)*ph_3(
    iM2) - (A_3(iM2)*A_3(iM2))*ph_2(iM2)*ph_3(iM2)) - (1.0/3.0)*gR2*d1ph_3(iP2)*(-((A_2(iP2)*A_2(iP2))*ph_0(iP2)*ph_1(iP2)) - (
    A_3(iP2)*A_3(iP2))*ph_0(iP2)*ph_1(iP2) + (A_2(iP2)*A_2(iP2))*ph_2(iP2)*ph_3(iP2) + (A_3(iP2)*A_3(iP2))*ph_2(iP2)*ph_3(iP2)) -
     (1.0/3.0)*gR2*d1ph_0(iM2)*(-((A_2(iM2)*A_2(iM2))*ph_0(iM2)*ph_2(iM2)) - (A_3(iM2)*A_3(iM2))*ph_0(iM2)*ph_2(iM2) - (A_2(iM2)*
    A_2(iM2))*ph_1(iM2)*ph_3(iM2) - (A_3(iM2)*A_3(iM2))*ph_1(iM2)*ph_3(iM2)) - (1.0/3.0)*gR2*d1ph_0(iM1)*(8.0*(A_2(iM1)*A_2(iM1
    ))*ph_0(iM1)*ph_2(iM1) + 8.0*(A_3(iM1)*A_3(iM1))*ph_0(iM1)*ph_2(iM1) + 8.0*(A_2(iM1)*A_2(iM1))*ph_1(iM1)*ph_3(iM1) + 8.0*(A_3
    (iM1)*A_3(iM1))*ph_1(iM1)*ph_3(iM1)) - (1.0/3.0)*gR2*d1ph_0(iP1)*(-8.0*(A_2(iP1)*A_2(iP1))*ph_0(iP1)*ph_2(iP1) - 8.0*(A_3(iP1
    )*A_3(iP1))*ph_0(iP1)*ph_2(iP1) - 8.0*(A_2(iP1)*A_2(iP1))*ph_1(iP1)*ph_3(iP1) - 8.0*(A_3(iP1)*A_3(iP1))*ph_1(iP1)*ph_3(iP1)) 
    - (1.0/3.0)*gR2*d1ph_3(iM1)*(-8.0*(A_2(iM1)*A_2(iM1))*ph_0(iM1)*ph_1(iM1) - 8.0*(A_3(iM1)*A_3(iM1))*ph_0(iM1)*ph_1(iM1) + 8.0
    *(A_2(iM1)*A_2(iM1))*ph_2(iM1)*ph_3(iM1) + 8.0*(A_3(iM1)*A_3(iM1))*ph_2(iM1)*ph_3(iM1)) - (1.0/3.0)*gR2*d1ph_3(iP1)*(8.0*(A_2
    (iP1)*A_2(iP1))*ph_0(iP1)*ph_1(iP1) + 8.0*(A_3(iP1)*A_3(iP1))*ph_0(iP1)*ph_1(iP1) - 8.0*(A_2(iP1)*A_2(iP1))*ph_2(iP1)*ph_3(
    iP1) - 8.0*(A_3(iP1)*A_3(iP1))*ph_2(iP1)*ph_3(iP1)) - (1.0/3.0)*gR2*d1ph_2(iP2)*((A_2(iP2)*A_2(iP2))*(ph_0(iP2)*ph_0(iP2)) + 
    (A_3(iP2)*A_3(iP2))*(ph_0(iP2)*ph_0(iP2)) + 4.0*(A_2(iP2)*A_2(iP2))*(ph_2(iP2)*ph_2(iP2)) + (A_2(iP2)*A_2(iP2))*(ph_3(iP2)*
    ph_3(iP2)) + 4.0*(A_3(iP2)*A_3(iP2))*(ph_2(iP2)*ph_2(iP2)) + (A_3(iP2)*A_3(iP2))*(ph_3(iP2)*ph_3(iP2))) - (1.0/3.0)*gR2*
    d1ph_2(iM2)*(-((A_2(iM2)*A_2(iM2))*(ph_0(iM2)*ph_0(iM2))) - (A_3(iM2)*A_3(iM2))*(ph_0(iM2)*ph_0(iM2)) - 4.0*(A_2(iM2)*A_2(iM2
    ))*(ph_2(iM2)*ph_2(iM2)) - (A_2(iM2)*A_2(iM2))*(ph_3(iM2)*ph_3(iM2)) - 4.0*(A_3(iM2)*A_3(iM2))*(ph_2(iM2)*ph_2(iM2)) - (A_3(
    iM2)*A_3(iM2))*(ph_3(iM2)*ph_3(iM2))) - (1.0/3.0)*gR2*d1ph_2(iM1)*(8.0*(A_2(iM1)*A_2(iM1))*(ph_0(iM1)*ph_0(iM1)) + 8.0*(A_3(
    iM1)*A_3(iM1))*(ph_0(iM1)*ph_0(iM1)) + 32.0*(A_2(iM1)*A_2(iM1))*(ph_2(iM1)*ph_2(iM1)) + 8.0*(A_2(iM1)*A_2(iM1))*(ph_3(iM1)*
    ph_3(iM1)) + 32.0*(A_3(iM1)*A_3(iM1))*(ph_2(iM1)*ph_2(iM1)) + 8.0*(A_3(iM1)*A_3(iM1))*(ph_3(iM1)*ph_3(iM1))) - (1.0/3.0)*gR2*
    d1ph_2(iP1)*(-8.0*(A_2(iP1)*A_2(iP1))*(ph_0(iP1)*ph_0(iP1)) - 8.0*(A_3(iP1)*A_3(iP1))*(ph_0(iP1)*ph_0(iP1)) - 32.0*(A_2(iP1)*
    A_2(iP1))*(ph_2(iP1)*ph_2(iP1)) - 8.0*(A_2(iP1)*A_2(iP1))*(ph_3(iP1)*ph_3(iP1)) - 32.0*(A_3(iP1)*A_3(iP1))*(ph_2(iP1)*ph_2(
    iP1)) - 8.0*(A_3(iP1)*A_3(iP1))*(ph_3(iP1)*ph_3(iP1))))/aR2 - 2.0*(a*a)*lambdaR2*ph_2(i)*(ph_0(i)*ph_0(i) + ph_1(i)*ph_1(i) +
     ph_2(i)*ph_2(i) + ph_3(i)*ph_3(i) - 1.0));

   dotMomenta.ph_3 = a*(d1ph_3(iM1)*(-2.0/3.0) + d1ph_3(iP1)*(2.0/3.0) + d1ph_3(iM2)*(1.0/12.0) - d1ph_3(iP2)*(1.0/12.0) + ((-
    1.0/3.0)*gR2*(d1ph_1(i)*d1ph_1(i))*(12.0*(A_2(i)*A_2(i))*ph_3(i) + 12.0*(A_3(i)*A_3(i))*ph_3(i)) - (1.0/3.0)*gR2*(d1ph_2(i)*
    d1ph_2(i))*(12.0*(A_2(i)*A_2(i))*ph_3(i) + 12.0*(A_3(i)*A_3(i))*ph_3(i)) - (1.0/3.0)*gR2*d1ph_0(i)*d1ph_2(i)*(12.0*(A_2(i)*
    A_2(i))*ph_1(i) + 12.0*(A_3(i)*A_3(i))*ph_1(i)) - (1.0/3.0)*gR2*d1ph_0(i)*d1ph_1(i)*(-12.0*(A_2(i)*A_2(i))*ph_2(i) - 12.0*(
    A_3(i)*A_3(i))*ph_2(i)) - (1.0/3.0)*gR2*d1ph_1(i)*d1ph_3(i)*(12.0*(A_2(i)*A_2(i))*ph_1(i) + 12.0*(A_3(i)*A_3(i))*ph_1(i)) - (
    1.0/3.0)*gR2*d1ph_2(i)*d1ph_3(i)*(12.0*(A_2(i)*A_2(i))*ph_2(i) + 12.0*(A_3(i)*A_3(i))*ph_2(i)) - (1.0/3.0)*gR2*d1ph_3(iP2)*((
    A_2(iP2)*A_2(iP2))*(ph_1(iP2)*ph_1(iP2)) + (A_2(iP2)*A_2(iP2))*(ph_2(iP2)*ph_2(iP2)) + (A_3(iP2)*A_3(iP2))*(ph_1(iP2)*ph_1(
    iP2)) + (A_3(iP2)*A_3(iP2))*(ph_2(iP2)*ph_2(iP2))) - (1.0/3.0)*gR2*d1ph_3(iM2)*(-((A_2(iM2)*A_2(iM2))*(ph_1(iM2)*ph_1(iM2))) 
    - (A_2(iM2)*A_2(iM2))*(ph_2(iM2)*ph_2(iM2)) - (A_3(iM2)*A_3(iM2))*(ph_1(iM2)*ph_1(iM2)) - (A_3(iM2)*A_3(iM2))*(ph_2(iM2)*ph_2
    (iM2))) - (1.0/3.0)*gR2*d1ph_3(iM1)*(8.0*(A_2(iM1)*A_2(iM1))*(ph_1(iM1)*ph_1(iM1)) + 8.0*(A_2(iM1)*A_2(iM1))*(ph_2(iM1)*ph_2(
    iM1)) + 8.0*(A_3(iM1)*A_3(iM1))*(ph_1(iM1)*ph_1(iM1)) + 8.0*(A_3(iM1)*A_3(iM1))*(ph_2(iM1)*ph_2(iM1))) - (1.0/3.0)*gR2*d1ph_3
    (iP1)*(-8.0*(A_2(iP1)*A_2(iP1))*(ph_1(iP1)*ph_1(iP1)) - 8.0*(A_2(iP1)*A_2(iP1))*(ph_2(iP1)*ph_2(iP1)) - 8.0*(A_3(iP1)*A_3(iP1
    ))*(ph_1(iP1)*ph_1(iP1)) - 8.0*(A_3(iP1)*A_3(iP1))*(ph_2(iP1)*ph_2(iP1))) - (1.0/3.0)*gR2*d1ph_1(iP2)*((A_2(iP2)*A_2(iP2))*
    ph_0(iP2)*ph_2(iP2) + (A_3(iP2)*A_3(iP2))*ph_0(iP2)*ph_2(iP2) + (A_2(iP2)*A_2(iP2))*ph_1(iP2)*ph_3(iP2) + (A_3(iP2)*A_3(iP2
    ))*ph_1(iP2)*ph_3(iP2)) - (1.0/3.0)*gR2*d1ph_2(iM2)*((A_2(iM2)*A_2(iM2))*ph_0(iM2)*ph_1(iM2) + (A_3(iM2)*A_3(iM2))*ph_0(iM2)*
    ph_1(iM2) - (A_2(iM2)*A_2(iM2))*ph_2(iM2)*ph_3(iM2) - (A_3(iM2)*A_3(iM2))*ph_2(iM2)*ph_3(iM2)) - (1.0/3.0)*gR2*d1ph_2(iP2
    )*(-((A_2(iP2)*A_2(iP2))*ph_0(iP2)*ph_1(iP2)) - (A_3(iP2)*A_3(iP2))*ph_0(iP2)*ph_1(iP2) + (A_2(iP2)*A_2(iP2))*ph_2(iP2)*ph_3(
    iP2) + (A_3(iP2)*A_3(iP2))*ph_2(iP2)*ph_3(iP2)) - (1.0/3.0)*gR2*d1ph_1(iM2)*(-((A_2(iM2)*A_2(iM2))*ph_0(iM2)*ph_2(iM2)) - (
    A_3(iM2)*A_3(iM2))*ph_0(iM2)*ph_2(iM2) - (A_2(iM2)*A_2(iM2))*ph_1(iM2)*ph_3(iM2) - (A_3(iM2)*A_3(iM2))*ph_1(iM2)*ph_3(iM2)) -
     (1.0/3.0)*gR2*d1ph_1(iM1)*(8.0*(A_2(iM1)*A_2(iM1))*ph_0(iM1)*ph_2(iM1) + 8.0*(A_3(iM1)*A_3(iM1))*ph_0(iM1)*ph_2(iM1) + 8.0*(
    A_2(iM1)*A_2(iM1))*ph_1(iM1)*ph_3(iM1) + 8.0*(A_3(iM1)*A_3(iM1))*ph_1(iM1)*ph_3(iM1)) - (1.0/3.0)*gR2*d1ph_1(iP1)*(-8.0*(A_2(
    iP1)*A_2(iP1))*ph_0(iP1)*ph_2(iP1) - 8.0*(A_3(iP1)*A_3(iP1))*ph_0(iP1)*ph_2(iP1) - 8.0*(A_2(iP1)*A_2(iP1))*ph_1(iP1)*ph_3(iP1
    ) - 8.0*(A_3(iP1)*A_3(iP1))*ph_1(iP1)*ph_3(iP1)) - (1.0/3.0)*gR2*d1ph_2(iM1)*(-8.0*(A_2(iM1)*A_2(iM1))*ph_0(iM1)*ph_1(iM1) - 
    8.0*(A_3(iM1)*A_3(iM1))*ph_0(iM1)*ph_1(iM1) + 8.0*(A_2(iM1)*A_2(iM1))*ph_2(iM1)*ph_3(iM1) + 8.0*(A_3(iM1)*A_3(iM1))*ph_2(iM1
    )*ph_3(iM1)) - (1.0/3.0)*gR2*d1ph_2(iP1)*(8.0*(A_2(iP1)*A_2(iP1))*ph_0(iP1)*ph_1(iP1) + 8.0*(A_3(iP1)*A_3(iP1))*ph_0(iP1)*
    ph_1(iP1) - 8.0*(A_2(iP1)*A_2(iP1))*ph_2(iP1)*ph_3(iP1) - 8.0*(A_3(iP1)*A_3(iP1))*ph_2(iP1)*ph_3(iP1)))/aR2 - 2.0*(a*a)*
    lambdaR2*ph_3(i)*(ph_0(i)*ph_0(i) + ph_1(i)*ph_1(i) + ph_2(i)*ph_2(i) + ph_3(i)*ph_3(i) - 1.0));

   dotMomenta.A_1 = 0.0;

   dotMomenta.A_2 = (d1A_2(iM1)*(-2.0/3.0) + d1A_2(iP1)*(2.0/3.0) + d1A_2(iM2)*(1.0/12.0) - d1A_2(iP2)*(1.0/12.0) - 32.0*gR2*A_2
    (i)*ph_1(i)*ph_2(i)*d1ph_1(i)*d1ph_2(i) + (d1ph_0(i)*d1ph_0(i))*(-4.0*gR2*A_2(i)*(ph_1(i)*ph_1(i)) - 4.0*gR2*A_2(i)*(ph_2(i)*
    ph_2(i))) + (d1ph_3(i)*d1ph_3(i))*(-4.0*gR2*A_2(i)*(ph_1(i)*ph_1(i)) - 4.0*gR2*A_2(i)*(ph_2(i)*ph_2(i))) + d1ph_0(i)*d1ph_1(i
    )*(-8.0*gR2*A_2(i)*ph_0(i)*ph_1(i) + 8.0*gR2*A_2(i)*ph_2(i)*ph_3(i)) + d1ph_0(i)*d1ph_2(i)*(-8.0*gR2*A_2(i)*ph_0(i)*ph_2(i) -
     8.0*gR2*A_2(i)*ph_1(i)*ph_3(i)) + d1ph_1(i)*d1ph_3(i)*(-8.0*gR2*A_2(i)*ph_0(i)*ph_2(i) - 8.0*gR2*A_2(i)*ph_1(i)*ph_3(i)) + 
    d1ph_2(i)*d1ph_3(i)*(8.0*gR2*A_2(i)*ph_0(i)*ph_1(i) - 8.0*gR2*A_2(i)*ph_2(i)*ph_3(i)) + (d1ph_1(i)*d1ph_1(i))*(-4.0*gR2*A_2(i
    )*(ph_0(i)*ph_0(i)) - 16.0*gR2*A_2(i)*(ph_1(i)*ph_1(i)) - 4.0*gR2*A_2(i)*(ph_3(i)*ph_3(i))) + (d1ph_2(i)*d1ph_2(i))*(-4.0*gR2
    *A_2(i)*(ph_0(i)*ph_0(i)) - 16.0*gR2*A_2(i)*(ph_2(i)*ph_2(i)) - 4.0*gR2*A_2(i)*(ph_3(i)*ph_3(i))))/a;

   dotMomenta.A_3 = (d1A_3(iM1)*(-2.0/3.0) + d1A_3(iP1)*(2.0/3.0) + d1A_3(iM2)*(1.0/12.0) - d1A_3(iP2)*(1.0/12.0) - 32.0*gR2*A_3
    (i)*ph_1(i)*ph_2(i)*d1ph_1(i)*d1ph_2(i) + (d1ph_0(i)*d1ph_0(i))*(-4.0*gR2*A_3(i)*(ph_1(i)*ph_1(i)) - 4.0*gR2*A_3(i)*(ph_2(i)*
    ph_2(i))) + (d1ph_3(i)*d1ph_3(i))*(-4.0*gR2*A_3(i)*(ph_1(i)*ph_1(i)) - 4.0*gR2*A_3(i)*(ph_2(i)*ph_2(i))) + d1ph_0(i)*d1ph_1(i
    )*(-8.0*gR2*A_3(i)*ph_0(i)*ph_1(i) + 8.0*gR2*A_3(i)*ph_2(i)*ph_3(i)) + d1ph_0(i)*d1ph_2(i)*(-8.0*gR2*A_3(i)*ph_0(i)*ph_2(i) -
     8.0*gR2*A_3(i)*ph_1(i)*ph_3(i)) + d1ph_1(i)*d1ph_3(i)*(-8.0*gR2*A_3(i)*ph_0(i)*ph_2(i) - 8.0*gR2*A_3(i)*ph_1(i)*ph_3(i)) + 
    d1ph_2(i)*d1ph_3(i)*(8.0*gR2*A_3(i)*ph_0(i)*ph_1(i) - 8.0*gR2*A_3(i)*ph_2(i)*ph_3(i)) + (d1ph_1(i)*d1ph_1(i))*(-4.0*gR2*A_3(i
    )*(ph_0(i)*ph_0(i)) - 16.0*gR2*A_3(i)*(ph_1(i)*ph_1(i)) - 4.0*gR2*A_3(i)*(ph_3(i)*ph_3(i))) + (d1ph_2(i)*d1ph_2(i))*(-4.0*gR2
    *A_3(i)*(ph_0(i)*ph_0(i)) - 16.0*gR2*A_3(i)*(ph_2(i)*ph_2(i)) - 4.0*gR2*A_3(i)*(ph_3(i)*ph_3(i))))/a;
}
