/* lattice.cpp : Basic lattice for scalar fields in polar form.
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

const int OUTPUT_PRECISION = 15;

typedef TFields TPrivateFields[9][9][9];
typedef TAux TPrivateAux[9][9][9];

class TLowpass
{
   double c_000, c_001, c_002, c_011, c_111, c_112, c_122, c_222;

public:

   void set_central_weight(double central_weight);

   TLowpass() { set_central_weight(0.0); };

   void compute(TFields fields[MAX_NSIDE][MAX_NSIDE][MAX_NSIDE], int i, int j, int k, int nSideM1, TFields *result);
   void compute(TPrivateFields fields, TFields *result);
};

void TLowpass::set_central_weight(double central_weight)
{
   c_000 = central_weight;

   c_001 = (-1368*c_000+305)/3240;
   c_002 = (63*c_000+2)/1620;
   c_011 = 2*(9*c_000+2)/135;
   c_111 = 37.0/3888;
   c_112 = -2*c_000/135;
   c_122 = (72*c_000-7)/38880;
   c_222 = (27*c_000+1)/19440;
}

void TLowpass::compute(TFields fields[MAX_NSIDE][MAX_NSIDE][MAX_NSIDE], int i, int j, int k, int nSideM1, TFields *result)
{
   *result = fields[i][j][k];
   result->mul_fields(c_000);

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

   TFields partial_sum;

   partial_sum = fields[iM2][jM2][kM2];
   partial_sum.add_fields(&(fields[iM2][jM2][kP2]));
   partial_sum.add_fields(&(fields[iM2][jP2][kM2]));
   partial_sum.add_fields(&(fields[iM2][jP2][kP2]));
   partial_sum.add_fields(&(fields[iP2][jM2][kM2]));
   partial_sum.add_fields(&(fields[iP2][jM2][kP2]));
   partial_sum.add_fields(&(fields[iP2][jP2][kM2]));
   partial_sum.add_fields(&(fields[iP2][jP2][kP2]));

   partial_sum.mul_fields(c_222);
   result->add_fields(&partial_sum);

   partial_sum = fields[iM2][jM2][kM1];
   partial_sum.add_fields(&(fields[iM2][jM2][kP1]));
   partial_sum.add_fields(&(fields[iM2][jM1][kM2]));
   partial_sum.add_fields(&(fields[iM2][jM1][kP2]));
   partial_sum.add_fields(&(fields[iM2][jP1][kM2]));
   partial_sum.add_fields(&(fields[iM2][jP1][kP2]));
   partial_sum.add_fields(&(fields[iM2][jP2][kM1]));
   partial_sum.add_fields(&(fields[iM2][jP2][kP1]));

   partial_sum.add_fields(&(fields[iM1][jM2][kM2]));
   partial_sum.add_fields(&(fields[iM1][jM2][kP2]));
   partial_sum.add_fields(&(fields[iM1][jP2][kM2]));
   partial_sum.add_fields(&(fields[iM1][jP2][kP2]));

   partial_sum.add_fields(&(fields[iP1][jM2][kM2]));
   partial_sum.add_fields(&(fields[iP1][jM2][kP2]));
   partial_sum.add_fields(&(fields[iP1][jP2][kM2]));
   partial_sum.add_fields(&(fields[iP1][jP2][kP2]));
   
   partial_sum.add_fields(&(fields[iP2][jM2][kM1]));
   partial_sum.add_fields(&(fields[iP2][jM2][kP1]));
   partial_sum.add_fields(&(fields[iP2][jM1][kM2]));
   partial_sum.add_fields(&(fields[iP2][jM1][kP2]));
   partial_sum.add_fields(&(fields[iP2][jP1][kM2]));
   partial_sum.add_fields(&(fields[iP2][jP1][kP2]));
   partial_sum.add_fields(&(fields[iP2][jP2][kM1]));
   partial_sum.add_fields(&(fields[iP2][jP2][kP1]));
   
   partial_sum.mul_fields(c_122);
   result->add_fields(&partial_sum);

   partial_sum = fields[iM2][jM1][kM1];
   partial_sum.add_fields(&(fields[iM2][jM1][kP1]));
   partial_sum.add_fields(&(fields[iM2][jP1][kM1]));
   partial_sum.add_fields(&(fields[iM2][jP1][kP1]));

   partial_sum.add_fields(&(fields[iM1][jM2][kM1]));
   partial_sum.add_fields(&(fields[iM1][jM2][kP1]));
   partial_sum.add_fields(&(fields[iM1][jM1][kM2]));
   partial_sum.add_fields(&(fields[iM1][jM1][kP2]));
   partial_sum.add_fields(&(fields[iM1][jP1][kM2]));
   partial_sum.add_fields(&(fields[iM1][jP1][kP2]));
   partial_sum.add_fields(&(fields[iM1][jP2][kM1]));
   partial_sum.add_fields(&(fields[iM1][jP2][kP1]));

   partial_sum.add_fields(&(fields[iP1][jM2][kM1]));
   partial_sum.add_fields(&(fields[iP1][jM2][kP1]));
   partial_sum.add_fields(&(fields[iP1][jM1][kM2]));
   partial_sum.add_fields(&(fields[iP1][jM1][kP2]));
   partial_sum.add_fields(&(fields[iP1][jP1][kM2]));
   partial_sum.add_fields(&(fields[iP1][jP1][kP2]));
   partial_sum.add_fields(&(fields[iP1][jP2][kM1]));
   partial_sum.add_fields(&(fields[iP1][jP2][kP1]));

   partial_sum.add_fields(&(fields[iP2][jM1][kM1]));
   partial_sum.add_fields(&(fields[iP2][jM1][kP1]));
   partial_sum.add_fields(&(fields[iP2][jP1][kM1]));
   partial_sum.add_fields(&(fields[iP2][jP1][kP1]));

   partial_sum.mul_fields(c_112);
   result->add_fields(&partial_sum);
   
   partial_sum = fields[iM1][jM1][kM1];
   partial_sum.add_fields(&(fields[iM1][jM1][kP1]));
   partial_sum.add_fields(&(fields[iM1][jP1][kM1]));
   partial_sum.add_fields(&(fields[iM1][jP1][kP1]));
   partial_sum.add_fields(&(fields[iP1][jM1][kM1]));
   partial_sum.add_fields(&(fields[iP1][jM1][kP1]));
   partial_sum.add_fields(&(fields[iP1][jP1][kM1]));
   partial_sum.add_fields(&(fields[iP1][jP1][kP1]));

   partial_sum.mul_fields(c_111);
   result->add_fields(&partial_sum);

   partial_sum = fields[iM1][jM1][k];
   partial_sum.add_fields(&(fields[iM1][j][kM1]));
   partial_sum.add_fields(&(fields[iM1][j][kP1]));
   partial_sum.add_fields(&(fields[iM1][jP1][k]));

   partial_sum.add_fields(&(fields[i][jM1][kM1]));
   partial_sum.add_fields(&(fields[i][jM1][kP1]));
   partial_sum.add_fields(&(fields[i][jP1][kM1]));
   partial_sum.add_fields(&(fields[i][jP1][kP1]));

   partial_sum.add_fields(&(fields[iP1][jM1][k]));
   partial_sum.add_fields(&(fields[iP1][j][kM1]));
   partial_sum.add_fields(&(fields[iP1][j][kP1]));
   partial_sum.add_fields(&(fields[iP1][jP1][k]));

   partial_sum.mul_fields(c_011);
   result->add_fields(&partial_sum);

   partial_sum = fields[iM2][j][k];
   partial_sum.add_fields(&(fields[i][jM2][k]));
   partial_sum.add_fields(&(fields[i][j][kM2]));
   partial_sum.add_fields(&(fields[i][j][kP2]));
   partial_sum.add_fields(&(fields[i][jP2][k]));
   partial_sum.add_fields(&(fields[iP2][j][k]));

   partial_sum.mul_fields(c_002);
   result->add_fields(&partial_sum);

   partial_sum = fields[iM1][j][k];
   partial_sum.add_fields(&(fields[i][jM1][k]));
   partial_sum.add_fields(&(fields[i][j][kM1]));
   partial_sum.add_fields(&(fields[i][j][kP1]));
   partial_sum.add_fields(&(fields[i][jP1][k]));
   partial_sum.add_fields(&(fields[iP1][j][k]));

   partial_sum.mul_fields(c_001);
   result->add_fields(&partial_sum);
}

void TLowpass::compute(TPrivateFields fields, TFields *result)
{
   *result = fields[4][4][4];
   result->mul_fields(c_000);

   TFields partial_sum;

   partial_sum = fields[2][2][2];
   partial_sum.add_fields(&(fields[2][2][6]));
   partial_sum.add_fields(&(fields[2][6][2]));
   partial_sum.add_fields(&(fields[2][6][6]));
   partial_sum.add_fields(&(fields[6][2][2]));
   partial_sum.add_fields(&(fields[6][2][6]));
   partial_sum.add_fields(&(fields[6][6][2]));
   partial_sum.add_fields(&(fields[6][6][6]));

   partial_sum.mul_fields(c_222);
   result->add_fields(&partial_sum);

   partial_sum = fields[2][2][3];
   partial_sum.add_fields(&(fields[2][2][5]));
   partial_sum.add_fields(&(fields[2][3][2]));
   partial_sum.add_fields(&(fields[2][3][6]));
   partial_sum.add_fields(&(fields[2][5][2]));
   partial_sum.add_fields(&(fields[2][5][6]));
   partial_sum.add_fields(&(fields[2][6][3]));
   partial_sum.add_fields(&(fields[2][6][5]));

   partial_sum.add_fields(&(fields[3][2][2]));
   partial_sum.add_fields(&(fields[3][2][6]));
   partial_sum.add_fields(&(fields[3][6][2]));
   partial_sum.add_fields(&(fields[3][6][6]));

   partial_sum.add_fields(&(fields[5][2][2]));
   partial_sum.add_fields(&(fields[5][2][6]));
   partial_sum.add_fields(&(fields[5][6][2]));
   partial_sum.add_fields(&(fields[5][6][6]));
   
   partial_sum.add_fields(&(fields[6][2][3]));
   partial_sum.add_fields(&(fields[6][2][5]));
   partial_sum.add_fields(&(fields[6][3][2]));
   partial_sum.add_fields(&(fields[6][3][6]));
   partial_sum.add_fields(&(fields[6][5][2]));
   partial_sum.add_fields(&(fields[6][5][6]));
   partial_sum.add_fields(&(fields[6][6][3]));
   partial_sum.add_fields(&(fields[6][6][5]));
   
   partial_sum.mul_fields(c_122);
   result->add_fields(&partial_sum);

   partial_sum = fields[2][3][3];
   partial_sum.add_fields(&(fields[2][3][5]));
   partial_sum.add_fields(&(fields[2][5][3]));
   partial_sum.add_fields(&(fields[2][5][5]));

   partial_sum.add_fields(&(fields[3][2][3]));
   partial_sum.add_fields(&(fields[3][2][5]));
   partial_sum.add_fields(&(fields[3][3][2]));
   partial_sum.add_fields(&(fields[3][3][6]));
   partial_sum.add_fields(&(fields[3][5][2]));
   partial_sum.add_fields(&(fields[3][5][6]));
   partial_sum.add_fields(&(fields[3][6][3]));
   partial_sum.add_fields(&(fields[3][6][5]));

   partial_sum.add_fields(&(fields[5][2][3]));
   partial_sum.add_fields(&(fields[5][2][5]));
   partial_sum.add_fields(&(fields[5][3][2]));
   partial_sum.add_fields(&(fields[5][3][6]));
   partial_sum.add_fields(&(fields[5][5][2]));
   partial_sum.add_fields(&(fields[5][5][6]));
   partial_sum.add_fields(&(fields[5][6][3]));
   partial_sum.add_fields(&(fields[5][6][5]));

   partial_sum.add_fields(&(fields[6][3][3]));
   partial_sum.add_fields(&(fields[6][3][5]));
   partial_sum.add_fields(&(fields[6][5][3]));
   partial_sum.add_fields(&(fields[6][5][5]));

   partial_sum.mul_fields(c_112);
   result->add_fields(&partial_sum);
   
   partial_sum = fields[3][3][3];
   partial_sum.add_fields(&(fields[3][3][5]));
   partial_sum.add_fields(&(fields[3][5][3]));
   partial_sum.add_fields(&(fields[3][5][5]));
   partial_sum.add_fields(&(fields[5][3][3]));
   partial_sum.add_fields(&(fields[5][3][5]));
   partial_sum.add_fields(&(fields[5][5][3]));
   partial_sum.add_fields(&(fields[5][5][5]));

   partial_sum.mul_fields(c_111);
   result->add_fields(&partial_sum);

   partial_sum = fields[3][3][4];
   partial_sum.add_fields(&(fields[3][4][3]));
   partial_sum.add_fields(&(fields[3][4][5]));
   partial_sum.add_fields(&(fields[3][5][4]));

   partial_sum.add_fields(&(fields[4][3][3]));
   partial_sum.add_fields(&(fields[4][3][5]));
   partial_sum.add_fields(&(fields[4][5][3]));
   partial_sum.add_fields(&(fields[4][5][5]));

   partial_sum.add_fields(&(fields[5][3][4]));
   partial_sum.add_fields(&(fields[5][4][3]));
   partial_sum.add_fields(&(fields[5][4][5]));
   partial_sum.add_fields(&(fields[5][5][4]));

   partial_sum.mul_fields(c_011);
   result->add_fields(&partial_sum);

   partial_sum = fields[2][4][4];
   partial_sum.add_fields(&(fields[4][2][4]));
   partial_sum.add_fields(&(fields[4][4][2]));
   partial_sum.add_fields(&(fields[4][4][6]));
   partial_sum.add_fields(&(fields[4][6][4]));
   partial_sum.add_fields(&(fields[6][4][4]));

   partial_sum.mul_fields(c_002);
   result->add_fields(&partial_sum);

   partial_sum = fields[3][4][4];
   partial_sum.add_fields(&(fields[4][3][4]));
   partial_sum.add_fields(&(fields[4][4][3]));
   partial_sum.add_fields(&(fields[4][4][5]));
   partial_sum.add_fields(&(fields[4][5][4]));
   partial_sum.add_fields(&(fields[5][4][4]));

   partial_sum.mul_fields(c_001);
   result->add_fields(&partial_sum);
}

using namespace std;

#include "string_aux.cpp"

class TLattice
{
protected:

   int nSideM1, active;
   
   TLowpass lowpass;

   void dDotFields_dFields(int indx, int momIndx, int i, int j, int k, double* dDotFields);

public:

   TPrivateFields privateFields;
   TPrivateAux privateAux;

   TLattice();

   void clear(const int idx);
   void set_nSide(const int new_nSide = MAX_NSIDE, const bool clear = true);
   int get_nSide() { return(nSideM1 + 1); }
   int get_active() { return(active); }
   void set_active(const int newActive) { active = newActive;}

   void compute_private_staticDotMomenta(TPrivateFields privateFields, TPrivateAux privateAux, TFields &dotMomenta);
   void compute_private_staticDotMomenta(TFields &dotMomenta) { compute_private_staticDotMomenta(privateFields, privateAux, dotMomenta); }
   void compute_staticDotMomenta(int indx, int i, int j, int k, TFields &dotMomenta);

   void copy2private(int indx, int i, int j, int k, TPrivateFields privateFields, TPrivateAux privateAux);
   void copy2private(int i, int j, int k) { copy2private(active, i, j, k, privateFields, privateAux); }
   void copy2affected_private_derivatives(int i, int j, int k, int component, TPrivateAux privateAux);
   void copy2affected_private_derivatives(int i, int j, int k, int component)
        { copy2affected_private_derivatives(i, j, k, component, privateAux); }
   void compute_affected_private_derivatives(int component, TPrivateFields privateFields, TPrivateAux privateAux);
   void compute_affected_private_derivatives(int component) 
        { compute_affected_private_derivatives(component, privateFields, privateAux); }
   
   void init_aux_static();
   double sum_rho_static();
   double affected_static_rho(int i, int j, int k);
   double affected_private_static_rho(TPrivateFields privateFields, TPrivateAux privateAux);
   double affected_private_static_rho() { return affected_private_static_rho(privateFields, privateAux); }

   void set_central_weight(double central_weight) { lowpass.set_central_weight(central_weight); };
   void lowpass_fields();
   
   void wrap_th(const int indx);

   bool write_vtk(string fileName, string title, double t);
   void read_vtk(string fileName, string& title);
};

TLattice::TLattice()
{
   active = 0;
   set_nSide();
}

void TLattice::clear(const int idx)
{
#ifdef _OPENMP
#pragma omp parallel
{
#pragma omp for
#endif
      for (int i = 0; i <= nSideM1; i++)
          for (int j = 0; j <= nSideM1; j++)
              for (int k = 0; k <= nSideM1; k++)
              {
                  fields[idx][i][j][k].clear();
              }
#ifdef _OPENMP
#pragma omp barrier
}
#endif   
}

void TLattice::set_nSide(const int new_nSide, const bool doClear)
{
   if (new_nSide < 5) nSideM1 = 4; else
      if (new_nSide > MAX_NSIDE) nSideM1 = MAX_NSIDE - 1; else
         nSideM1 = new_nSide - 1;

   if (doClear) clear(active);
}

bool TLattice::write_vtk(string fileName, string title, double t)
{
   int i, j, k;

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
   vtk << title << "; t = " << t << endl;
   vtk << "ASCII" << endl ;
   vtk << "DATASET STRUCTURED_POINTS" << endl ;
   vtk << "DIMENSIONS " << nSideM1+1 << ' ' << nSideM1+1 << ' ' << nSideM1+1 << endl;
   vtk << "ORIGIN 0 0 0" << endl ;
   vtk << "SPACING " << 1 << ' ' << 1 << ' ' << 1 << endl;
   vtk << "POINT_DATA " << (nSideM1+1)*(nSideM1+1)*(nSideM1+1) << endl;

   vtk << "VECTORS T DOUBLE" << endl;

   for (i = 0; i <= nSideM1; i++)
       for (j = 0; j <= nSideM1; j++)
           for (k = 0; k <= nSideM1; k++)
           {
              vtk << fields[active][i][j][k].th_1 << ' ' << fields[active][i][j][k].th_2 << ' ' << fields[active][i][j][k].th_3 << endl;
           }

   vtk << "VECTORS A DOUBLE" << endl;

   for (i = 0; i <= nSideM1; i++)
       for (j = 0; j <= nSideM1; j++)
           for (k = 0; k <= nSideM1; k++)
           {
              vtk << fields[active][i][j][k].A_1 << ' ' << fields[active][i][j][k].A_2 << ' ' << fields[active][i][j][k].A_3 << endl;
           }

   vtk << "VECTORS pT DOUBLE" << endl;

   for (i = 0; i <= nSideM1; i++)
       for (j = 0; j <= nSideM1; j++)
           for (k = 0; k <= nSideM1; k++)
           {
              vtk << 0 << ' ' << 0 << ' ' << 0 << endl;
           }

   vtk << "VECTORS pA DOUBLE" << endl;

   for (i = 0; i <= nSideM1; i++)
       for (j = 0; j <= nSideM1; j++)
           for (k = 0; k <= nSideM1; k++)
           {
              vtk << 0 << ' ' << 0 << ' ' << 0 << endl;
           }

   vtk.close();

   return(true);
}

void TLattice::read_vtk(string fileName, string& title)
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

            s = token(line, i);

            if (sNSide != s)
            {
               cout << "Wrong y dimension " << s << " (should be " << sNSide << ")" << endl << flush;
               exit(3);
            }

            s = token(line, i);

            if (sNSide != s)
            {
               cout << "Wrong z dimension " << s << " (should be " << sNSide << ")" << endl << flush;
               exit(4);
            }
         }
         else if ("VECTORS T DOUBLE" == line)
         {
            for (int i = 0; i <= nSideM1; i++)
                for (int j = 0; j <= nSideM1; j++)
                    for (int k = 0; k <= nSideM1; k++)
                    {
                       string line;
                       getline(inFile, line);

                       int fromTo = 0;

                       fields[active][i][j][k].th_1 = doubleToken(line, fromTo);
                       fields[active][i][j][k].th_2 = doubleToken(line, fromTo);
                       fields[active][i][j][k].th_3 = doubleToken(line, fromTo);
                    }
         }
         else if ("VECTORS A DOUBLE" == line)
         {
            for (int i = 0; i <= nSideM1; i++)
                for (int j = 0; j <= nSideM1; j++)
                    for (int k = 0; k <= nSideM1; k++)
                    {
                       string line;
                       getline(inFile, line);

                       int fromTo = 0;

                       fields[active][i][j][k].A_1 = doubleToken(line, fromTo);
                       fields[active][i][j][k].A_2 = doubleToken(line, fromTo);
                       fields[active][i][j][k].A_3 = doubleToken(line, fromTo);
                    }
         }
         else if ("VECTORS pT DOUBLE" == line)
         {
            for (int i = 0; i <= nSideM1; i++)
                for (int j = 0; j <= nSideM1; j++)
                    for (int k = 0; k <= nSideM1; k++)
                    {
                       string line;
                       getline(inFile, line);
                    }
         }
         else if ("VECTORS pA DOUBLE" == line)
         {
            for (int i = 0; i <= nSideM1; i++)
                for (int j = 0; j <= nSideM1; j++)
                    for (int k = 0; k <= nSideM1; k++)
                    {
                       string line;
                       getline(inFile, line);
                    }
         }
      }

      inFile.close();
   }
}

void TLattice::compute_private_staticDotMomenta(TPrivateFields privateFields, TPrivateAux privateAux, TFields &dotMomenta)
/* Computes time derivatives of  conjugate momenta by differentiating the Hamiltonian w.r.t. fields.
   IMPORTANT: assumes privateAux data (h, G, M, field derivatives) for cell and affected neighbors are ready for use! */
{
   // dot pth_a = - dH/d th_a

   dotMomenta.th_1 = 0.0;
   dotMomenta.th_2 = 0.0;
   dotMomenta.th_3 = 0.0;

   privateAux[6][5][5].subtract_offset_drho_dth(-2,-1,-1, &(privateFields[6][5][5]), dotMomenta.th_1, dotMomenta.th_2, dotMomenta.th_3);
   privateAux[6][5][4].subtract_offset_drho_dth(-2,-1,0, &(privateFields[6][5][4]), dotMomenta.th_1, dotMomenta.th_2, dotMomenta.th_3);
   privateAux[6][5][3].subtract_offset_drho_dth(-2,-1,1, &(privateFields[6][5][3]), dotMomenta.th_1, dotMomenta.th_2, dotMomenta.th_3);
   privateAux[6][4][5].subtract_offset_drho_dth(-2,0,-1, &(privateFields[6][4][5]), dotMomenta.th_1, dotMomenta.th_2, dotMomenta.th_3);
   privateAux[6][4][4].subtract_offset_drho_dth(-2,0,0, &(privateFields[6][4][4]), dotMomenta.th_1, dotMomenta.th_2, dotMomenta.th_3);
   privateAux[6][4][3].subtract_offset_drho_dth(-2,0,1, &(privateFields[6][4][3]), dotMomenta.th_1, dotMomenta.th_2, dotMomenta.th_3);
   privateAux[6][3][5].subtract_offset_drho_dth(-2,1,-1, &(privateFields[6][3][5]), dotMomenta.th_1, dotMomenta.th_2, dotMomenta.th_3);
   privateAux[6][3][4].subtract_offset_drho_dth(-2,1,0, &(privateFields[6][3][4]), dotMomenta.th_1, dotMomenta.th_2, dotMomenta.th_3);
   privateAux[6][3][3].subtract_offset_drho_dth(-2,1,1, &(privateFields[6][3][3]), dotMomenta.th_1, dotMomenta.th_2, dotMomenta.th_3);

   privateAux[5][6][5].subtract_offset_drho_dth(-1,-2,-1, &(privateFields[5][6][5]), dotMomenta.th_1, dotMomenta.th_2, dotMomenta.th_3);
   privateAux[5][6][4].subtract_offset_drho_dth(-1,-2,0, &(privateFields[5][6][4]), dotMomenta.th_1, dotMomenta.th_2, dotMomenta.th_3);
   privateAux[5][6][3].subtract_offset_drho_dth(-1,-2,1, &(privateFields[5][6][3]), dotMomenta.th_1, dotMomenta.th_2, dotMomenta.th_3);
   privateAux[5][5][6].subtract_offset_drho_dth(-1,-1,-2, &(privateFields[5][5][6]), dotMomenta.th_1, dotMomenta.th_2, dotMomenta.th_3);
   privateAux[5][5][5].subtract_offset_drho_dth(-1,-1,-1, &(privateFields[5][5][5]), dotMomenta.th_1, dotMomenta.th_2, dotMomenta.th_3);
   privateAux[5][5][4].subtract_offset_drho_dth(-1,-1,0, &(privateFields[5][5][4]), dotMomenta.th_1, dotMomenta.th_2, dotMomenta.th_3);
   privateAux[5][5][3].subtract_offset_drho_dth(-1,-1,1, &(privateFields[5][5][3]), dotMomenta.th_1, dotMomenta.th_2, dotMomenta.th_3);
   privateAux[5][5][2].subtract_offset_drho_dth(-1,-1,2, &(privateFields[5][5][2]), dotMomenta.th_1, dotMomenta.th_2, dotMomenta.th_3);
   privateAux[5][4][6].subtract_offset_drho_dth(-1,0,-2, &(privateFields[5][4][6]), dotMomenta.th_1, dotMomenta.th_2, dotMomenta.th_3);
   privateAux[5][4][5].subtract_offset_drho_dth(-1,0,-1, &(privateFields[5][4][5]), dotMomenta.th_1, dotMomenta.th_2, dotMomenta.th_3);
   privateAux[5][4][4].subtract_offset_drho_dth(-1,0,0, &(privateFields[5][4][4]), dotMomenta.th_1, dotMomenta.th_2, dotMomenta.th_3);
   privateAux[5][4][3].subtract_offset_drho_dth(-1,0,1, &(privateFields[5][4][3]), dotMomenta.th_1, dotMomenta.th_2, dotMomenta.th_3);
   privateAux[5][4][2].subtract_offset_drho_dth(-1,0,2, &(privateFields[5][4][2]), dotMomenta.th_1, dotMomenta.th_2, dotMomenta.th_3);
   privateAux[5][3][6].subtract_offset_drho_dth(-1,1,-2, &(privateFields[5][3][6]), dotMomenta.th_1, dotMomenta.th_2, dotMomenta.th_3);
   privateAux[5][3][5].subtract_offset_drho_dth(-1,1,-1, &(privateFields[5][3][5]), dotMomenta.th_1, dotMomenta.th_2, dotMomenta.th_3);
   privateAux[5][3][4].subtract_offset_drho_dth(-1,1,0, &(privateFields[5][3][4]), dotMomenta.th_1, dotMomenta.th_2, dotMomenta.th_3);
   privateAux[5][3][3].subtract_offset_drho_dth(-1,1,1, &(privateFields[5][3][3]), dotMomenta.th_1, dotMomenta.th_2, dotMomenta.th_3);
   privateAux[5][3][2].subtract_offset_drho_dth(-1,1,2, &(privateFields[5][3][2]), dotMomenta.th_1, dotMomenta.th_2, dotMomenta.th_3);
   privateAux[5][2][5].subtract_offset_drho_dth(-1,2,-1, &(privateFields[5][2][5]), dotMomenta.th_1, dotMomenta.th_2, dotMomenta.th_3);
   privateAux[5][2][4].subtract_offset_drho_dth(-1,2,0, &(privateFields[5][2][4]), dotMomenta.th_1, dotMomenta.th_2, dotMomenta.th_3);
   privateAux[5][2][3].subtract_offset_drho_dth(-1,2,1, &(privateFields[5][2][3]), dotMomenta.th_1, dotMomenta.th_2, dotMomenta.th_3);

   privateAux[4][6][5].subtract_offset_drho_dth(0,-2,-1, &(privateFields[4][6][5]), dotMomenta.th_1, dotMomenta.th_2, dotMomenta.th_3);
   privateAux[4][6][4].subtract_offset_drho_dth(0,-2,0, &(privateFields[4][6][4]), dotMomenta.th_1, dotMomenta.th_2, dotMomenta.th_3);
   privateAux[4][6][3].subtract_offset_drho_dth(0,-2,1, &(privateFields[4][6][3]), dotMomenta.th_1, dotMomenta.th_2, dotMomenta.th_3);
   privateAux[4][5][6].subtract_offset_drho_dth(0,-1,-2, &(privateFields[4][5][6]), dotMomenta.th_1, dotMomenta.th_2, dotMomenta.th_3);
   privateAux[4][5][5].subtract_offset_drho_dth(0,-1,-1, &(privateFields[4][5][5]), dotMomenta.th_1, dotMomenta.th_2, dotMomenta.th_3);
   privateAux[4][5][4].subtract_offset_drho_dth(0,-1,0, &(privateFields[4][5][4]), dotMomenta.th_1, dotMomenta.th_2, dotMomenta.th_3);
   privateAux[4][5][3].subtract_offset_drho_dth(0,-1,1, &(privateFields[4][5][3]), dotMomenta.th_1, dotMomenta.th_2, dotMomenta.th_3);
   privateAux[4][5][2].subtract_offset_drho_dth(0,-1,2, &(privateFields[4][5][2]), dotMomenta.th_1, dotMomenta.th_2, dotMomenta.th_3);
   privateAux[4][4][6].subtract_offset_drho_dth(0,0,-2, &(privateFields[4][4][6]), dotMomenta.th_1, dotMomenta.th_2, dotMomenta.th_3);
   privateAux[4][4][5].subtract_offset_drho_dth(0,0,-1, &(privateFields[4][4][5]), dotMomenta.th_1, dotMomenta.th_2, dotMomenta.th_3);
   privateAux[4][4][3].subtract_offset_drho_dth(0,0,1, &(privateFields[4][4][3]), dotMomenta.th_1, dotMomenta.th_2, dotMomenta.th_3);
   privateAux[4][4][2].subtract_offset_drho_dth(0,0,2, &(privateFields[4][4][2]), dotMomenta.th_1, dotMomenta.th_2, dotMomenta.th_3);
   privateAux[4][3][6].subtract_offset_drho_dth(0,1,-2, &(privateFields[4][3][6]), dotMomenta.th_1, dotMomenta.th_2, dotMomenta.th_3);
   privateAux[4][3][5].subtract_offset_drho_dth(0,1,-1, &(privateFields[4][3][5]), dotMomenta.th_1, dotMomenta.th_2, dotMomenta.th_3);
   privateAux[4][3][4].subtract_offset_drho_dth(0,1,0, &(privateFields[4][3][4]), dotMomenta.th_1, dotMomenta.th_2, dotMomenta.th_3);
   privateAux[4][3][3].subtract_offset_drho_dth(0,1,1, &(privateFields[4][3][3]), dotMomenta.th_1, dotMomenta.th_2, dotMomenta.th_3);
   privateAux[4][3][2].subtract_offset_drho_dth(0,1,2, &(privateFields[4][3][2]), dotMomenta.th_1, dotMomenta.th_2, dotMomenta.th_3);
   privateAux[4][2][5].subtract_offset_drho_dth(0,2,-1, &(privateFields[4][2][5]), dotMomenta.th_1, dotMomenta.th_2, dotMomenta.th_3);
   privateAux[4][2][4].subtract_offset_drho_dth(0,2,0, &(privateFields[4][2][4]), dotMomenta.th_1, dotMomenta.th_2, dotMomenta.th_3);
   privateAux[4][2][3].subtract_offset_drho_dth(0,2,1, &(privateFields[4][2][3]), dotMomenta.th_1, dotMomenta.th_2, dotMomenta.th_3);

   privateAux[3][6][5].subtract_offset_drho_dth(1,-2,-1, &(privateFields[3][6][5]), dotMomenta.th_1, dotMomenta.th_2, dotMomenta.th_3);
   privateAux[3][6][4].subtract_offset_drho_dth(1,-2,0, &(privateFields[3][6][4]), dotMomenta.th_1, dotMomenta.th_2, dotMomenta.th_3);
   privateAux[3][6][3].subtract_offset_drho_dth(1,-2,1, &(privateFields[3][6][3]), dotMomenta.th_1, dotMomenta.th_2, dotMomenta.th_3);
   privateAux[3][5][6].subtract_offset_drho_dth(1,-1,-2, &(privateFields[3][5][6]), dotMomenta.th_1, dotMomenta.th_2, dotMomenta.th_3);
   privateAux[3][5][5].subtract_offset_drho_dth(1,-1,-1, &(privateFields[3][5][5]), dotMomenta.th_1, dotMomenta.th_2, dotMomenta.th_3);
   privateAux[3][5][4].subtract_offset_drho_dth(1,-1,0, &(privateFields[3][5][4]), dotMomenta.th_1, dotMomenta.th_2, dotMomenta.th_3);
   privateAux[3][5][3].subtract_offset_drho_dth(1,-1,1, &(privateFields[3][5][3]), dotMomenta.th_1, dotMomenta.th_2, dotMomenta.th_3);
   privateAux[3][5][2].subtract_offset_drho_dth(1,-1,2, &(privateFields[3][5][2]), dotMomenta.th_1, dotMomenta.th_2, dotMomenta.th_3);
   privateAux[3][4][6].subtract_offset_drho_dth(1,0,-2, &(privateFields[3][4][6]), dotMomenta.th_1, dotMomenta.th_2, dotMomenta.th_3);
   privateAux[3][4][5].subtract_offset_drho_dth(1,0,-1, &(privateFields[3][4][5]), dotMomenta.th_1, dotMomenta.th_2, dotMomenta.th_3);
   privateAux[3][4][4].subtract_offset_drho_dth(1,0,0, &(privateFields[3][4][4]), dotMomenta.th_1, dotMomenta.th_2, dotMomenta.th_3);
   privateAux[3][4][3].subtract_offset_drho_dth(1,0,1, &(privateFields[3][4][3]), dotMomenta.th_1, dotMomenta.th_2, dotMomenta.th_3);
   privateAux[3][4][2].subtract_offset_drho_dth(1,0,2, &(privateFields[3][4][2]), dotMomenta.th_1, dotMomenta.th_2, dotMomenta.th_3);
   privateAux[3][3][6].subtract_offset_drho_dth(1,1,-2, &(privateFields[3][3][6]), dotMomenta.th_1, dotMomenta.th_2, dotMomenta.th_3);
   privateAux[3][3][5].subtract_offset_drho_dth(1,1,-1, &(privateFields[3][3][5]), dotMomenta.th_1, dotMomenta.th_2, dotMomenta.th_3);
   privateAux[3][3][4].subtract_offset_drho_dth(1,1,0, &(privateFields[3][3][4]), dotMomenta.th_1, dotMomenta.th_2, dotMomenta.th_3);
   privateAux[3][3][3].subtract_offset_drho_dth(1,1,1, &(privateFields[3][3][3]), dotMomenta.th_1, dotMomenta.th_2, dotMomenta.th_3);
   privateAux[3][3][2].subtract_offset_drho_dth(1,1,2, &(privateFields[3][3][2]), dotMomenta.th_1, dotMomenta.th_2, dotMomenta.th_3);
   privateAux[3][2][5].subtract_offset_drho_dth(1,2,-1, &(privateFields[3][2][5]), dotMomenta.th_1, dotMomenta.th_2, dotMomenta.th_3);
   privateAux[3][2][4].subtract_offset_drho_dth(1,2,0, &(privateFields[3][2][4]), dotMomenta.th_1, dotMomenta.th_2, dotMomenta.th_3);
   privateAux[3][2][3].subtract_offset_drho_dth(1,2,1, &(privateFields[3][2][3]), dotMomenta.th_1, dotMomenta.th_2, dotMomenta.th_3);

   privateAux[2][5][5].subtract_offset_drho_dth(2,-1,-1, &(privateFields[2][5][5]), dotMomenta.th_1, dotMomenta.th_2, dotMomenta.th_3);
   privateAux[2][5][4].subtract_offset_drho_dth(2,-1,0, &(privateFields[2][5][4]), dotMomenta.th_1, dotMomenta.th_2, dotMomenta.th_3);
   privateAux[2][5][3].subtract_offset_drho_dth(2,-1,1, &(privateFields[2][5][3]), dotMomenta.th_1, dotMomenta.th_2, dotMomenta.th_3);
   privateAux[2][4][5].subtract_offset_drho_dth(2,0,-1, &(privateFields[2][4][5]), dotMomenta.th_1, dotMomenta.th_2, dotMomenta.th_3);
   privateAux[2][4][4].subtract_offset_drho_dth(2,0,0, &(privateFields[2][4][4]), dotMomenta.th_1, dotMomenta.th_2, dotMomenta.th_3);
   privateAux[2][4][3].subtract_offset_drho_dth(2,0,1, &(privateFields[2][4][3]), dotMomenta.th_1, dotMomenta.th_2, dotMomenta.th_3);
   privateAux[2][3][5].subtract_offset_drho_dth(2,1,-1, &(privateFields[2][3][5]), dotMomenta.th_1, dotMomenta.th_2, dotMomenta.th_3);
   privateAux[2][3][4].subtract_offset_drho_dth(2,1,0, &(privateFields[2][3][4]), dotMomenta.th_1, dotMomenta.th_2, dotMomenta.th_3);
   privateAux[2][3][3].subtract_offset_drho_dth(2,1,1, &(privateFields[2][3][3]), dotMomenta.th_1, dotMomenta.th_2, dotMomenta.th_3);

   privateAux[4][4][4].subtract_drho_static_dth(&(privateFields[4][4][4]), dotMomenta.th_1, dotMomenta.th_2, dotMomenta.th_3);

   // dot pA_m =-dH/d A_m

   dotMomenta.A_1 = 0.0;
   dotMomenta.A_2 = 0.0;
   dotMomenta.A_3 = 0.0;

   privateAux[6][5][5].subtract_offset_drho_dA(-2,-1,-1, dotMomenta.A_1, dotMomenta.A_2, dotMomenta.A_3);
   privateAux[6][5][4].subtract_offset_drho_dA(-2,-1,0, dotMomenta.A_1, dotMomenta.A_2, dotMomenta.A_3);
   privateAux[6][5][3].subtract_offset_drho_dA(-2,-1,1, dotMomenta.A_1, dotMomenta.A_2, dotMomenta.A_3);
   privateAux[6][4][5].subtract_offset_drho_dA(-2,0,-1, dotMomenta.A_1, dotMomenta.A_2, dotMomenta.A_3);
   privateAux[6][4][4].subtract_offset_drho_dA(-2,0,0, dotMomenta.A_1, dotMomenta.A_2, dotMomenta.A_3);
   privateAux[6][4][3].subtract_offset_drho_dA(-2,0,1, dotMomenta.A_1, dotMomenta.A_2, dotMomenta.A_3);
   privateAux[6][3][5].subtract_offset_drho_dA(-2,1,-1, dotMomenta.A_1, dotMomenta.A_2, dotMomenta.A_3);
   privateAux[6][3][4].subtract_offset_drho_dA(-2,1,0, dotMomenta.A_1, dotMomenta.A_2, dotMomenta.A_3);
   privateAux[6][3][3].subtract_offset_drho_dA(-2,1,1, dotMomenta.A_1, dotMomenta.A_2, dotMomenta.A_3);

   privateAux[5][6][5].subtract_offset_drho_dA(-1,-2,-1, dotMomenta.A_1, dotMomenta.A_2, dotMomenta.A_3);
   privateAux[5][6][4].subtract_offset_drho_dA(-1,-2,0, dotMomenta.A_1, dotMomenta.A_2, dotMomenta.A_3);
   privateAux[5][6][3].subtract_offset_drho_dA(-1,-2,1, dotMomenta.A_1, dotMomenta.A_2, dotMomenta.A_3);
   privateAux[5][5][6].subtract_offset_drho_dA(-1,-1,-2, dotMomenta.A_1, dotMomenta.A_2, dotMomenta.A_3);
   privateAux[5][5][5].subtract_offset_drho_dA(-1,-1,-1, dotMomenta.A_1, dotMomenta.A_2, dotMomenta.A_3);
   privateAux[5][5][4].subtract_offset_drho_dA(-1,-1,0, dotMomenta.A_1, dotMomenta.A_2, dotMomenta.A_3);
   privateAux[5][5][3].subtract_offset_drho_dA(-1,-1,1, dotMomenta.A_1, dotMomenta.A_2, dotMomenta.A_3);
   privateAux[5][5][2].subtract_offset_drho_dA(-1,-1,2, dotMomenta.A_1, dotMomenta.A_2, dotMomenta.A_3);
   privateAux[5][4][6].subtract_offset_drho_dA(-1,0,-2, dotMomenta.A_1, dotMomenta.A_2, dotMomenta.A_3);
   privateAux[5][4][5].subtract_offset_drho_dA(-1,0,-1, dotMomenta.A_1, dotMomenta.A_2, dotMomenta.A_3);
   privateAux[5][4][4].subtract_offset_drho_dA(-1,0,0, dotMomenta.A_1, dotMomenta.A_2, dotMomenta.A_3);
   privateAux[5][4][3].subtract_offset_drho_dA(-1,0,1, dotMomenta.A_1, dotMomenta.A_2, dotMomenta.A_3);
   privateAux[5][4][2].subtract_offset_drho_dA(-1,0,2, dotMomenta.A_1, dotMomenta.A_2, dotMomenta.A_3);
   privateAux[5][3][6].subtract_offset_drho_dA(-1,1,-2, dotMomenta.A_1, dotMomenta.A_2, dotMomenta.A_3);
   privateAux[5][3][5].subtract_offset_drho_dA(-1,1,-1, dotMomenta.A_1, dotMomenta.A_2, dotMomenta.A_3);
   privateAux[5][3][4].subtract_offset_drho_dA(-1,1,0, dotMomenta.A_1, dotMomenta.A_2, dotMomenta.A_3);
   privateAux[5][3][3].subtract_offset_drho_dA(-1,1,1, dotMomenta.A_1, dotMomenta.A_2, dotMomenta.A_3);
   privateAux[5][3][2].subtract_offset_drho_dA(-1,1,2, dotMomenta.A_1, dotMomenta.A_2, dotMomenta.A_3);
   privateAux[5][2][5].subtract_offset_drho_dA(-1,2,-1, dotMomenta.A_1, dotMomenta.A_2, dotMomenta.A_3);
   privateAux[5][2][4].subtract_offset_drho_dA(-1,2,0, dotMomenta.A_1, dotMomenta.A_2, dotMomenta.A_3);
   privateAux[5][2][3].subtract_offset_drho_dA(-1,2,1, dotMomenta.A_1, dotMomenta.A_2, dotMomenta.A_3);

   privateAux[4][6][5].subtract_offset_drho_dA(0,-2,-1, dotMomenta.A_1, dotMomenta.A_2, dotMomenta.A_3);
   privateAux[4][6][4].subtract_offset_drho_dA(0,-2,0, dotMomenta.A_1, dotMomenta.A_2, dotMomenta.A_3);
   privateAux[4][6][3].subtract_offset_drho_dA(0,-2,1, dotMomenta.A_1, dotMomenta.A_2, dotMomenta.A_3);
   privateAux[4][5][6].subtract_offset_drho_dA(0,-1,-2, dotMomenta.A_1, dotMomenta.A_2, dotMomenta.A_3);
   privateAux[4][5][5].subtract_offset_drho_dA(0,-1,-1, dotMomenta.A_1, dotMomenta.A_2, dotMomenta.A_3);
   privateAux[4][5][4].subtract_offset_drho_dA(0,-1,0, dotMomenta.A_1, dotMomenta.A_2, dotMomenta.A_3);
   privateAux[4][5][3].subtract_offset_drho_dA(0,-1,1, dotMomenta.A_1, dotMomenta.A_2, dotMomenta.A_3);
   privateAux[4][5][2].subtract_offset_drho_dA(0,-1,2, dotMomenta.A_1, dotMomenta.A_2, dotMomenta.A_3);
   privateAux[4][4][6].subtract_offset_drho_dA(0,0,-2, dotMomenta.A_1, dotMomenta.A_2, dotMomenta.A_3);
   privateAux[4][4][5].subtract_offset_drho_dA(0,0,-1, dotMomenta.A_1, dotMomenta.A_2, dotMomenta.A_3);
   privateAux[4][4][3].subtract_offset_drho_dA(0,0,1, dotMomenta.A_1, dotMomenta.A_2, dotMomenta.A_3);
   privateAux[4][4][2].subtract_offset_drho_dA(0,0,2, dotMomenta.A_1, dotMomenta.A_2, dotMomenta.A_3);
   privateAux[4][3][6].subtract_offset_drho_dA(0,1,-2, dotMomenta.A_1, dotMomenta.A_2, dotMomenta.A_3);
   privateAux[4][3][5].subtract_offset_drho_dA(0,1,-1, dotMomenta.A_1, dotMomenta.A_2, dotMomenta.A_3);
   privateAux[4][3][4].subtract_offset_drho_dA(0,1,0, dotMomenta.A_1, dotMomenta.A_2, dotMomenta.A_3);
   privateAux[4][3][3].subtract_offset_drho_dA(0,1,1, dotMomenta.A_1, dotMomenta.A_2, dotMomenta.A_3);
   privateAux[4][3][2].subtract_offset_drho_dA(0,1,2, dotMomenta.A_1, dotMomenta.A_2, dotMomenta.A_3);
   privateAux[4][2][5].subtract_offset_drho_dA(0,2,-1, dotMomenta.A_1, dotMomenta.A_2, dotMomenta.A_3);
   privateAux[4][2][4].subtract_offset_drho_dA(0,2,0, dotMomenta.A_1, dotMomenta.A_2, dotMomenta.A_3);
   privateAux[4][2][3].subtract_offset_drho_dA(0,2,1, dotMomenta.A_1, dotMomenta.A_2, dotMomenta.A_3);

   privateAux[3][6][5].subtract_offset_drho_dA(1,-2,-1, dotMomenta.A_1, dotMomenta.A_2, dotMomenta.A_3);
   privateAux[3][6][4].subtract_offset_drho_dA(1,-2,0, dotMomenta.A_1, dotMomenta.A_2, dotMomenta.A_3);
   privateAux[3][6][3].subtract_offset_drho_dA(1,-2,1, dotMomenta.A_1, dotMomenta.A_2, dotMomenta.A_3);
   privateAux[3][5][6].subtract_offset_drho_dA(1,-1,-2, dotMomenta.A_1, dotMomenta.A_2, dotMomenta.A_3);
   privateAux[3][5][5].subtract_offset_drho_dA(1,-1,-1, dotMomenta.A_1, dotMomenta.A_2, dotMomenta.A_3);
   privateAux[3][5][4].subtract_offset_drho_dA(1,-1,0, dotMomenta.A_1, dotMomenta.A_2, dotMomenta.A_3);
   privateAux[3][5][3].subtract_offset_drho_dA(1,-1,1, dotMomenta.A_1, dotMomenta.A_2, dotMomenta.A_3);
   privateAux[3][5][2].subtract_offset_drho_dA(1,-1,2, dotMomenta.A_1, dotMomenta.A_2, dotMomenta.A_3);
   privateAux[3][4][6].subtract_offset_drho_dA(1,0,-2, dotMomenta.A_1, dotMomenta.A_2, dotMomenta.A_3);
   privateAux[3][4][5].subtract_offset_drho_dA(1,0,-1, dotMomenta.A_1, dotMomenta.A_2, dotMomenta.A_3);
   privateAux[3][4][4].subtract_offset_drho_dA(1,0,0, dotMomenta.A_1, dotMomenta.A_2, dotMomenta.A_3);
   privateAux[3][4][3].subtract_offset_drho_dA(1,0,1, dotMomenta.A_1, dotMomenta.A_2, dotMomenta.A_3);
   privateAux[3][4][2].subtract_offset_drho_dA(1,0,2, dotMomenta.A_1, dotMomenta.A_2, dotMomenta.A_3);
   privateAux[3][3][6].subtract_offset_drho_dA(1,1,-2, dotMomenta.A_1, dotMomenta.A_2, dotMomenta.A_3);
   privateAux[3][3][5].subtract_offset_drho_dA(1,1,-1, dotMomenta.A_1, dotMomenta.A_2, dotMomenta.A_3);
   privateAux[3][3][4].subtract_offset_drho_dA(1,1,0, dotMomenta.A_1, dotMomenta.A_2, dotMomenta.A_3);
   privateAux[3][3][3].subtract_offset_drho_dA(1,1,1, dotMomenta.A_1, dotMomenta.A_2, dotMomenta.A_3);
   privateAux[3][3][2].subtract_offset_drho_dA(1,1,2, dotMomenta.A_1, dotMomenta.A_2, dotMomenta.A_3);
   privateAux[3][2][5].subtract_offset_drho_dA(1,2,-1, dotMomenta.A_1, dotMomenta.A_2, dotMomenta.A_3);
   privateAux[3][2][4].subtract_offset_drho_dA(1,2,0, dotMomenta.A_1, dotMomenta.A_2, dotMomenta.A_3);
   privateAux[3][2][3].subtract_offset_drho_dA(1,2,1, dotMomenta.A_1, dotMomenta.A_2, dotMomenta.A_3);

   privateAux[2][5][5].subtract_offset_drho_dA(2,-1,-1, dotMomenta.A_1, dotMomenta.A_2, dotMomenta.A_3);
   privateAux[2][5][4].subtract_offset_drho_dA(2,-1,0, dotMomenta.A_1, dotMomenta.A_2, dotMomenta.A_3);
   privateAux[2][5][3].subtract_offset_drho_dA(2,-1,1, dotMomenta.A_1, dotMomenta.A_2, dotMomenta.A_3);
   privateAux[2][4][5].subtract_offset_drho_dA(2,0,-1, dotMomenta.A_1, dotMomenta.A_2, dotMomenta.A_3);
   privateAux[2][4][4].subtract_offset_drho_dA(2,0,0, dotMomenta.A_1, dotMomenta.A_2, dotMomenta.A_3);
   privateAux[2][4][3].subtract_offset_drho_dA(2,0,1, dotMomenta.A_1, dotMomenta.A_2, dotMomenta.A_3);
   privateAux[2][3][5].subtract_offset_drho_dA(2,1,-1, dotMomenta.A_1, dotMomenta.A_2, dotMomenta.A_3);
   privateAux[2][3][4].subtract_offset_drho_dA(2,1,0, dotMomenta.A_1, dotMomenta.A_2, dotMomenta.A_3);
   privateAux[2][3][3].subtract_offset_drho_dA(2,1,1, dotMomenta.A_1, dotMomenta.A_2, dotMomenta.A_3);

   privateAux[4][4][4].subtract_drho_static_dA(&(privateFields[4][4][4]), dotMomenta.A_1, dotMomenta.A_2, dotMomenta.A_3);
}

void TLattice::compute_staticDotMomenta(int indx, int i, int j, int k, TFields &dotMomenta)
/* Computes time derivatives of  conjugate momenta by differentiating the Hamiltonian w.r.t. fields.
   IMPORTANT: assumes aux data (h, G, M, field derivatives) for cell and affected neighbors are ready for use! */
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

   // dot pth_a = - dH/d th_a

   dotMomenta.th_1 = 0.0;
   dotMomenta.th_2 = 0.0;
   dotMomenta.th_3 = 0.0;

   aux[iP2][jP1][kP1].subtract_offset_drho_dth(-2,-1,-1, &(fields[indx][iP2][jP1][kP1]), dotMomenta.th_1, dotMomenta.th_2, dotMomenta.th_3);
   aux[iP2][jP1][k].subtract_offset_drho_dth(-2,-1,0, &(fields[indx][iP2][jP1][k]), dotMomenta.th_1, dotMomenta.th_2, dotMomenta.th_3);
   aux[iP2][jP1][kM1].subtract_offset_drho_dth(-2,-1,1, &(fields[indx][iP2][jP1][kM1]), dotMomenta.th_1, dotMomenta.th_2, dotMomenta.th_3);
   aux[iP2][j][kP1].subtract_offset_drho_dth(-2,0,-1, &(fields[indx][iP2][j][kP1]), dotMomenta.th_1, dotMomenta.th_2, dotMomenta.th_3);
   aux[iP2][j][k].subtract_offset_drho_dth(-2,0,0, &(fields[indx][iP2][j][k]), dotMomenta.th_1, dotMomenta.th_2, dotMomenta.th_3);
   aux[iP2][j][kM1].subtract_offset_drho_dth(-2,0,1, &(fields[indx][iP2][j][kM1]), dotMomenta.th_1, dotMomenta.th_2, dotMomenta.th_3);
   aux[iP2][jM1][kP1].subtract_offset_drho_dth(-2,1,-1, &(fields[indx][iP2][jM1][kP1]), dotMomenta.th_1, dotMomenta.th_2, dotMomenta.th_3);
   aux[iP2][jM1][k].subtract_offset_drho_dth(-2,1,0, &(fields[indx][iP2][jM1][k]), dotMomenta.th_1, dotMomenta.th_2, dotMomenta.th_3);
   aux[iP2][jM1][kM1].subtract_offset_drho_dth(-2,1,1, &(fields[indx][iP2][jM1][kM1]), dotMomenta.th_1, dotMomenta.th_2, dotMomenta.th_3);

   aux[iP1][jP2][kP1].subtract_offset_drho_dth(-1,-2,-1, &(fields[indx][iP1][jP2][kP1]), dotMomenta.th_1, dotMomenta.th_2, dotMomenta.th_3);
   aux[iP1][jP2][k].subtract_offset_drho_dth(-1,-2,0, &(fields[indx][iP1][jP2][k]), dotMomenta.th_1, dotMomenta.th_2, dotMomenta.th_3);
   aux[iP1][jP2][kM1].subtract_offset_drho_dth(-1,-2,1, &(fields[indx][iP1][jP2][kM1]), dotMomenta.th_1, dotMomenta.th_2, dotMomenta.th_3);
   aux[iP1][jP1][kP2].subtract_offset_drho_dth(-1,-1,-2, &(fields[indx][iP1][jP1][kP2]), dotMomenta.th_1, dotMomenta.th_2, dotMomenta.th_3);
   aux[iP1][jP1][kP1].subtract_offset_drho_dth(-1,-1,-1, &(fields[indx][iP1][jP1][kP1]), dotMomenta.th_1, dotMomenta.th_2, dotMomenta.th_3);
   aux[iP1][jP1][k].subtract_offset_drho_dth(-1,-1,0, &(fields[indx][iP1][jP1][k]), dotMomenta.th_1, dotMomenta.th_2, dotMomenta.th_3);
   aux[iP1][jP1][kM1].subtract_offset_drho_dth(-1,-1,1, &(fields[indx][iP1][jP1][kM1]), dotMomenta.th_1, dotMomenta.th_2, dotMomenta.th_3);
   aux[iP1][jP1][kM2].subtract_offset_drho_dth(-1,-1,2, &(fields[indx][iP1][jP1][kM2]), dotMomenta.th_1, dotMomenta.th_2, dotMomenta.th_3);
   aux[iP1][j][kP2].subtract_offset_drho_dth(-1,0,-2, &(fields[indx][iP1][j][kP2]), dotMomenta.th_1, dotMomenta.th_2, dotMomenta.th_3);
   aux[iP1][j][kP1].subtract_offset_drho_dth(-1,0,-1, &(fields[indx][iP1][j][kP1]), dotMomenta.th_1, dotMomenta.th_2, dotMomenta.th_3);
   aux[iP1][j][k].subtract_offset_drho_dth(-1,0,0, &(fields[indx][iP1][j][k]), dotMomenta.th_1, dotMomenta.th_2, dotMomenta.th_3);
   aux[iP1][j][kM1].subtract_offset_drho_dth(-1,0,1, &(fields[indx][iP1][j][kM1]), dotMomenta.th_1, dotMomenta.th_2, dotMomenta.th_3);
   aux[iP1][j][kM2].subtract_offset_drho_dth(-1,0,2, &(fields[indx][iP1][j][kM2]), dotMomenta.th_1, dotMomenta.th_2, dotMomenta.th_3);
   aux[iP1][jM1][kP2].subtract_offset_drho_dth(-1,1,-2, &(fields[indx][iP1][jM1][kP2]), dotMomenta.th_1, dotMomenta.th_2, dotMomenta.th_3);
   aux[iP1][jM1][kP1].subtract_offset_drho_dth(-1,1,-1, &(fields[indx][iP1][jM1][kP1]), dotMomenta.th_1, dotMomenta.th_2, dotMomenta.th_3);
   aux[iP1][jM1][k].subtract_offset_drho_dth(-1,1,0, &(fields[indx][iP1][jM1][k]), dotMomenta.th_1, dotMomenta.th_2, dotMomenta.th_3);
   aux[iP1][jM1][kM1].subtract_offset_drho_dth(-1,1,1, &(fields[indx][iP1][jM1][kM1]), dotMomenta.th_1, dotMomenta.th_2, dotMomenta.th_3);
   aux[iP1][jM1][kM2].subtract_offset_drho_dth(-1,1,2, &(fields[indx][iP1][jM1][kM2]), dotMomenta.th_1, dotMomenta.th_2, dotMomenta.th_3);
   aux[iP1][jM2][kP1].subtract_offset_drho_dth(-1,2,-1, &(fields[indx][iP1][jM2][kP1]), dotMomenta.th_1, dotMomenta.th_2, dotMomenta.th_3);
   aux[iP1][jM2][k].subtract_offset_drho_dth(-1,2,0, &(fields[indx][iP1][jM2][k]), dotMomenta.th_1, dotMomenta.th_2, dotMomenta.th_3);
   aux[iP1][jM2][kM1].subtract_offset_drho_dth(-1,2,1, &(fields[indx][iP1][jM2][kM1]), dotMomenta.th_1, dotMomenta.th_2, dotMomenta.th_3);

   aux[i][jP2][kP1].subtract_offset_drho_dth(0,-2,-1, &(fields[indx][i][jP2][kP1]), dotMomenta.th_1, dotMomenta.th_2, dotMomenta.th_3);
   aux[i][jP2][k].subtract_offset_drho_dth(0,-2,0, &(fields[indx][i][jP2][k]), dotMomenta.th_1, dotMomenta.th_2, dotMomenta.th_3);
   aux[i][jP2][kM1].subtract_offset_drho_dth(0,-2,1, &(fields[indx][i][jP2][kM1]), dotMomenta.th_1, dotMomenta.th_2, dotMomenta.th_3);
   aux[i][jP1][kP2].subtract_offset_drho_dth(0,-1,-2, &(fields[indx][i][jP1][kP2]), dotMomenta.th_1, dotMomenta.th_2, dotMomenta.th_3);
   aux[i][jP1][kP1].subtract_offset_drho_dth(0,-1,-1, &(fields[indx][i][jP1][kP1]), dotMomenta.th_1, dotMomenta.th_2, dotMomenta.th_3);
   aux[i][jP1][k].subtract_offset_drho_dth(0,-1,0, &(fields[indx][i][jP1][k]), dotMomenta.th_1, dotMomenta.th_2, dotMomenta.th_3);
   aux[i][jP1][kM1].subtract_offset_drho_dth(0,-1,1, &(fields[indx][i][jP1][kM1]), dotMomenta.th_1, dotMomenta.th_2, dotMomenta.th_3);
   aux[i][jP1][kM2].subtract_offset_drho_dth(0,-1,2, &(fields[indx][i][jP1][kM2]), dotMomenta.th_1, dotMomenta.th_2, dotMomenta.th_3);
   aux[i][j][kP2].subtract_offset_drho_dth(0,0,-2, &(fields[indx][i][j][kP2]), dotMomenta.th_1, dotMomenta.th_2, dotMomenta.th_3);
   aux[i][j][kP1].subtract_offset_drho_dth(0,0,-1, &(fields[indx][i][j][kP1]), dotMomenta.th_1, dotMomenta.th_2, dotMomenta.th_3);
   aux[i][j][kM1].subtract_offset_drho_dth(0,0,1, &(fields[indx][i][j][kM1]), dotMomenta.th_1, dotMomenta.th_2, dotMomenta.th_3);
   aux[i][j][kM2].subtract_offset_drho_dth(0,0,2, &(fields[indx][i][j][kM2]), dotMomenta.th_1, dotMomenta.th_2, dotMomenta.th_3);
   aux[i][jM1][kP2].subtract_offset_drho_dth(0,1,-2, &(fields[indx][i][jM1][kP2]), dotMomenta.th_1, dotMomenta.th_2, dotMomenta.th_3);
   aux[i][jM1][kP1].subtract_offset_drho_dth(0,1,-1, &(fields[indx][i][jM1][kP1]), dotMomenta.th_1, dotMomenta.th_2, dotMomenta.th_3);
   aux[i][jM1][k].subtract_offset_drho_dth(0,1,0, &(fields[indx][i][jM1][k]), dotMomenta.th_1, dotMomenta.th_2, dotMomenta.th_3);
   aux[i][jM1][kM1].subtract_offset_drho_dth(0,1,1, &(fields[indx][i][jM1][kM1]), dotMomenta.th_1, dotMomenta.th_2, dotMomenta.th_3);
   aux[i][jM1][kM2].subtract_offset_drho_dth(0,1,2, &(fields[indx][i][jM1][kM2]), dotMomenta.th_1, dotMomenta.th_2, dotMomenta.th_3);
   aux[i][jM2][kP1].subtract_offset_drho_dth(0,2,-1, &(fields[indx][i][jM2][kP1]), dotMomenta.th_1, dotMomenta.th_2, dotMomenta.th_3);
   aux[i][jM2][k].subtract_offset_drho_dth(0,2,0, &(fields[indx][i][jM2][k]), dotMomenta.th_1, dotMomenta.th_2, dotMomenta.th_3);
   aux[i][jM2][kM1].subtract_offset_drho_dth(0,2,1, &(fields[indx][i][jM2][kM1]), dotMomenta.th_1, dotMomenta.th_2, dotMomenta.th_3);

   aux[iM1][jP2][kP1].subtract_offset_drho_dth(1,-2,-1, &(fields[indx][iM1][jP2][kP1]), dotMomenta.th_1, dotMomenta.th_2, dotMomenta.th_3);
   aux[iM1][jP2][k].subtract_offset_drho_dth(1,-2,0, &(fields[indx][iM1][jP2][k]), dotMomenta.th_1, dotMomenta.th_2, dotMomenta.th_3);
   aux[iM1][jP2][kM1].subtract_offset_drho_dth(1,-2,1, &(fields[indx][iM1][jP2][kM1]), dotMomenta.th_1, dotMomenta.th_2, dotMomenta.th_3);
   aux[iM1][jP1][kP2].subtract_offset_drho_dth(1,-1,-2, &(fields[indx][iM1][jP1][kP2]), dotMomenta.th_1, dotMomenta.th_2, dotMomenta.th_3);
   aux[iM1][jP1][kP1].subtract_offset_drho_dth(1,-1,-1, &(fields[indx][iM1][jP1][kP1]), dotMomenta.th_1, dotMomenta.th_2, dotMomenta.th_3);
   aux[iM1][jP1][k].subtract_offset_drho_dth(1,-1,0, &(fields[indx][iM1][jP1][k]), dotMomenta.th_1, dotMomenta.th_2, dotMomenta.th_3);
   aux[iM1][jP1][kM1].subtract_offset_drho_dth(1,-1,1, &(fields[indx][iM1][jP1][kM1]), dotMomenta.th_1, dotMomenta.th_2, dotMomenta.th_3);
   aux[iM1][jP1][kM2].subtract_offset_drho_dth(1,-1,2, &(fields[indx][iM1][jP1][kM2]), dotMomenta.th_1, dotMomenta.th_2, dotMomenta.th_3);
   aux[iM1][j][kP2].subtract_offset_drho_dth(1,0,-2, &(fields[indx][iM1][j][kP2]), dotMomenta.th_1, dotMomenta.th_2, dotMomenta.th_3);
   aux[iM1][j][kP1].subtract_offset_drho_dth(1,0,-1, &(fields[indx][iM1][j][kP1]), dotMomenta.th_1, dotMomenta.th_2, dotMomenta.th_3);
   aux[iM1][j][k].subtract_offset_drho_dth(1,0,0, &(fields[indx][iM1][j][k]), dotMomenta.th_1, dotMomenta.th_2, dotMomenta.th_3);
   aux[iM1][j][kM1].subtract_offset_drho_dth(1,0,1, &(fields[indx][iM1][j][kM1]), dotMomenta.th_1, dotMomenta.th_2, dotMomenta.th_3);
   aux[iM1][j][kM2].subtract_offset_drho_dth(1,0,2, &(fields[indx][iM1][j][kM2]), dotMomenta.th_1, dotMomenta.th_2, dotMomenta.th_3);
   aux[iM1][jM1][kP2].subtract_offset_drho_dth(1,1,-2, &(fields[indx][iM1][jM1][kP2]), dotMomenta.th_1, dotMomenta.th_2, dotMomenta.th_3);
   aux[iM1][jM1][kP1].subtract_offset_drho_dth(1,1,-1, &(fields[indx][iM1][jM1][kP1]), dotMomenta.th_1, dotMomenta.th_2, dotMomenta.th_3);
   aux[iM1][jM1][k].subtract_offset_drho_dth(1,1,0, &(fields[indx][iM1][jM1][k]), dotMomenta.th_1, dotMomenta.th_2, dotMomenta.th_3);
   aux[iM1][jM1][kM1].subtract_offset_drho_dth(1,1,1, &(fields[indx][iM1][jM1][kM1]), dotMomenta.th_1, dotMomenta.th_2, dotMomenta.th_3);
   aux[iM1][jM1][kM2].subtract_offset_drho_dth(1,1,2, &(fields[indx][iM1][jM1][kM2]), dotMomenta.th_1, dotMomenta.th_2, dotMomenta.th_3);
   aux[iM1][jM2][kP1].subtract_offset_drho_dth(1,2,-1, &(fields[indx][iM1][jM2][kP1]), dotMomenta.th_1, dotMomenta.th_2, dotMomenta.th_3);
   aux[iM1][jM2][k].subtract_offset_drho_dth(1,2,0, &(fields[indx][iM1][jM2][k]), dotMomenta.th_1, dotMomenta.th_2, dotMomenta.th_3);
   aux[iM1][jM2][kM1].subtract_offset_drho_dth(1,2,1, &(fields[indx][iM1][jM2][kM1]), dotMomenta.th_1, dotMomenta.th_2, dotMomenta.th_3);

   aux[iM2][jP1][kP1].subtract_offset_drho_dth(2,-1,-1, &(fields[indx][iM2][jP1][kP1]), dotMomenta.th_1, dotMomenta.th_2, dotMomenta.th_3);
   aux[iM2][jP1][k].subtract_offset_drho_dth(2,-1,0, &(fields[indx][iM2][jP1][k]), dotMomenta.th_1, dotMomenta.th_2, dotMomenta.th_3);
   aux[iM2][jP1][kM1].subtract_offset_drho_dth(2,-1,1, &(fields[indx][iM2][jP1][kM1]), dotMomenta.th_1, dotMomenta.th_2, dotMomenta.th_3);
   aux[iM2][j][kP1].subtract_offset_drho_dth(2,0,-1, &(fields[indx][iM2][j][kP1]), dotMomenta.th_1, dotMomenta.th_2, dotMomenta.th_3);
   aux[iM2][j][k].subtract_offset_drho_dth(2,0,0, &(fields[indx][iM2][j][k]), dotMomenta.th_1, dotMomenta.th_2, dotMomenta.th_3);
   aux[iM2][j][kM1].subtract_offset_drho_dth(2,0,1, &(fields[indx][iM2][j][kM1]), dotMomenta.th_1, dotMomenta.th_2, dotMomenta.th_3);
   aux[iM2][jM1][kP1].subtract_offset_drho_dth(2,1,-1, &(fields[indx][iM2][jM1][kP1]), dotMomenta.th_1, dotMomenta.th_2, dotMomenta.th_3);
   aux[iM2][jM1][k].subtract_offset_drho_dth(2,1,0, &(fields[indx][iM2][jM1][k]), dotMomenta.th_1, dotMomenta.th_2, dotMomenta.th_3);
   aux[iM2][jM1][kM1].subtract_offset_drho_dth(2,1,1, &(fields[indx][iM2][jM1][kM1]), dotMomenta.th_1, dotMomenta.th_2, dotMomenta.th_3);

   aux[i][j][k].subtract_drho_static_dth(&(fields[indx][i][j][k]), dotMomenta.th_1, dotMomenta.th_2, dotMomenta.th_3);

   // dot pA_m =-dH/d A_m

   dotMomenta.A_1 = 0.0;
   dotMomenta.A_2 = 0.0;
   dotMomenta.A_3 = 0.0;

   aux[iP2][jP1][kP1].subtract_offset_drho_dA(-2,-1,-1, dotMomenta.A_1, dotMomenta.A_2, dotMomenta.A_3);
   aux[iP2][jP1][k].subtract_offset_drho_dA(-2,-1,0, dotMomenta.A_1, dotMomenta.A_2, dotMomenta.A_3);
   aux[iP2][jP1][kM1].subtract_offset_drho_dA(-2,-1,1, dotMomenta.A_1, dotMomenta.A_2, dotMomenta.A_3);
   aux[iP2][j][kP1].subtract_offset_drho_dA(-2,0,-1, dotMomenta.A_1, dotMomenta.A_2, dotMomenta.A_3);
   aux[iP2][j][k].subtract_offset_drho_dA(-2,0,0, dotMomenta.A_1, dotMomenta.A_2, dotMomenta.A_3);
   aux[iP2][j][kM1].subtract_offset_drho_dA(-2,0,1, dotMomenta.A_1, dotMomenta.A_2, dotMomenta.A_3);
   aux[iP2][jM1][kP1].subtract_offset_drho_dA(-2,1,-1, dotMomenta.A_1, dotMomenta.A_2, dotMomenta.A_3);
   aux[iP2][jM1][k].subtract_offset_drho_dA(-2,1,0, dotMomenta.A_1, dotMomenta.A_2, dotMomenta.A_3);
   aux[iP2][jM1][kM1].subtract_offset_drho_dA(-2,1,1, dotMomenta.A_1, dotMomenta.A_2, dotMomenta.A_3);

   aux[iP1][jP2][kP1].subtract_offset_drho_dA(-1,-2,-1, dotMomenta.A_1, dotMomenta.A_2, dotMomenta.A_3);
   aux[iP1][jP2][k].subtract_offset_drho_dA(-1,-2,0, dotMomenta.A_1, dotMomenta.A_2, dotMomenta.A_3);
   aux[iP1][jP2][kM1].subtract_offset_drho_dA(-1,-2,1, dotMomenta.A_1, dotMomenta.A_2, dotMomenta.A_3);
   aux[iP1][jP1][kP2].subtract_offset_drho_dA(-1,-1,-2, dotMomenta.A_1, dotMomenta.A_2, dotMomenta.A_3);
   aux[iP1][jP1][kP1].subtract_offset_drho_dA(-1,-1,-1, dotMomenta.A_1, dotMomenta.A_2, dotMomenta.A_3);
   aux[iP1][jP1][k].subtract_offset_drho_dA(-1,-1,0, dotMomenta.A_1, dotMomenta.A_2, dotMomenta.A_3);
   aux[iP1][jP1][kM1].subtract_offset_drho_dA(-1,-1,1, dotMomenta.A_1, dotMomenta.A_2, dotMomenta.A_3);
   aux[iP1][jP1][kM2].subtract_offset_drho_dA(-1,-1,2, dotMomenta.A_1, dotMomenta.A_2, dotMomenta.A_3);
   aux[iP1][j][kP2].subtract_offset_drho_dA(-1,0,-2, dotMomenta.A_1, dotMomenta.A_2, dotMomenta.A_3);
   aux[iP1][j][kP1].subtract_offset_drho_dA(-1,0,-1, dotMomenta.A_1, dotMomenta.A_2, dotMomenta.A_3);
   aux[iP1][j][k].subtract_offset_drho_dA(-1,0,0, dotMomenta.A_1, dotMomenta.A_2, dotMomenta.A_3);
   aux[iP1][j][kM1].subtract_offset_drho_dA(-1,0,1, dotMomenta.A_1, dotMomenta.A_2, dotMomenta.A_3);
   aux[iP1][j][kM2].subtract_offset_drho_dA(-1,0,2, dotMomenta.A_1, dotMomenta.A_2, dotMomenta.A_3);
   aux[iP1][jM1][kP2].subtract_offset_drho_dA(-1,1,-2, dotMomenta.A_1, dotMomenta.A_2, dotMomenta.A_3);
   aux[iP1][jM1][kP1].subtract_offset_drho_dA(-1,1,-1, dotMomenta.A_1, dotMomenta.A_2, dotMomenta.A_3);
   aux[iP1][jM1][k].subtract_offset_drho_dA(-1,1,0, dotMomenta.A_1, dotMomenta.A_2, dotMomenta.A_3);
   aux[iP1][jM1][kM1].subtract_offset_drho_dA(-1,1,1, dotMomenta.A_1, dotMomenta.A_2, dotMomenta.A_3);
   aux[iP1][jM1][kM2].subtract_offset_drho_dA(-1,1,2, dotMomenta.A_1, dotMomenta.A_2, dotMomenta.A_3);
   aux[iP1][jM2][kP1].subtract_offset_drho_dA(-1,2,-1, dotMomenta.A_1, dotMomenta.A_2, dotMomenta.A_3);
   aux[iP1][jM2][k].subtract_offset_drho_dA(-1,2,0, dotMomenta.A_1, dotMomenta.A_2, dotMomenta.A_3);
   aux[iP1][jM2][kM1].subtract_offset_drho_dA(-1,2,1, dotMomenta.A_1, dotMomenta.A_2, dotMomenta.A_3);

   aux[i][jP2][kP1].subtract_offset_drho_dA(0,-2,-1, dotMomenta.A_1, dotMomenta.A_2, dotMomenta.A_3);
   aux[i][jP2][k].subtract_offset_drho_dA(0,-2,0, dotMomenta.A_1, dotMomenta.A_2, dotMomenta.A_3);
   aux[i][jP2][kM1].subtract_offset_drho_dA(0,-2,1, dotMomenta.A_1, dotMomenta.A_2, dotMomenta.A_3);
   aux[i][jP1][kP2].subtract_offset_drho_dA(0,-1,-2, dotMomenta.A_1, dotMomenta.A_2, dotMomenta.A_3);
   aux[i][jP1][kP1].subtract_offset_drho_dA(0,-1,-1, dotMomenta.A_1, dotMomenta.A_2, dotMomenta.A_3);
   aux[i][jP1][k].subtract_offset_drho_dA(0,-1,0, dotMomenta.A_1, dotMomenta.A_2, dotMomenta.A_3);
   aux[i][jP1][kM1].subtract_offset_drho_dA(0,-1,1, dotMomenta.A_1, dotMomenta.A_2, dotMomenta.A_3);
   aux[i][jP1][kM2].subtract_offset_drho_dA(0,-1,2, dotMomenta.A_1, dotMomenta.A_2, dotMomenta.A_3);
   aux[i][j][kP2].subtract_offset_drho_dA(0,0,-2, dotMomenta.A_1, dotMomenta.A_2, dotMomenta.A_3);
   aux[i][j][kP1].subtract_offset_drho_dA(0,0,-1, dotMomenta.A_1, dotMomenta.A_2, dotMomenta.A_3);
   aux[i][j][kM1].subtract_offset_drho_dA(0,0,1, dotMomenta.A_1, dotMomenta.A_2, dotMomenta.A_3);
   aux[i][j][kM2].subtract_offset_drho_dA(0,0,2, dotMomenta.A_1, dotMomenta.A_2, dotMomenta.A_3);
   aux[i][jM1][kP2].subtract_offset_drho_dA(0,1,-2, dotMomenta.A_1, dotMomenta.A_2, dotMomenta.A_3);
   aux[i][jM1][kP1].subtract_offset_drho_dA(0,1,-1, dotMomenta.A_1, dotMomenta.A_2, dotMomenta.A_3);
   aux[i][jM1][k].subtract_offset_drho_dA(0,1,0, dotMomenta.A_1, dotMomenta.A_2, dotMomenta.A_3);
   aux[i][jM1][kM1].subtract_offset_drho_dA(0,1,1, dotMomenta.A_1, dotMomenta.A_2, dotMomenta.A_3);
   aux[i][jM1][kM2].subtract_offset_drho_dA(0,1,2, dotMomenta.A_1, dotMomenta.A_2, dotMomenta.A_3);
   aux[i][jM2][kP1].subtract_offset_drho_dA(0,2,-1, dotMomenta.A_1, dotMomenta.A_2, dotMomenta.A_3);
   aux[i][jM2][k].subtract_offset_drho_dA(0,2,0, dotMomenta.A_1, dotMomenta.A_2, dotMomenta.A_3);
   aux[i][jM2][kM1].subtract_offset_drho_dA(0,2,1, dotMomenta.A_1, dotMomenta.A_2, dotMomenta.A_3);

   aux[iM1][jP2][kP1].subtract_offset_drho_dA(1,-2,-1, dotMomenta.A_1, dotMomenta.A_2, dotMomenta.A_3);
   aux[iM1][jP2][k].subtract_offset_drho_dA(1,-2,0, dotMomenta.A_1, dotMomenta.A_2, dotMomenta.A_3);
   aux[iM1][jP2][kM1].subtract_offset_drho_dA(1,-2,1, dotMomenta.A_1, dotMomenta.A_2, dotMomenta.A_3);
   aux[iM1][jP1][kP2].subtract_offset_drho_dA(1,-1,-2, dotMomenta.A_1, dotMomenta.A_2, dotMomenta.A_3);
   aux[iM1][jP1][kP1].subtract_offset_drho_dA(1,-1,-1, dotMomenta.A_1, dotMomenta.A_2, dotMomenta.A_3);
   aux[iM1][jP1][k].subtract_offset_drho_dA(1,-1,0, dotMomenta.A_1, dotMomenta.A_2, dotMomenta.A_3);
   aux[iM1][jP1][kM1].subtract_offset_drho_dA(1,-1,1, dotMomenta.A_1, dotMomenta.A_2, dotMomenta.A_3);
   aux[iM1][jP1][kM2].subtract_offset_drho_dA(1,-1,2, dotMomenta.A_1, dotMomenta.A_2, dotMomenta.A_3);
   aux[iM1][j][kP2].subtract_offset_drho_dA(1,0,-2, dotMomenta.A_1, dotMomenta.A_2, dotMomenta.A_3);
   aux[iM1][j][kP1].subtract_offset_drho_dA(1,0,-1, dotMomenta.A_1, dotMomenta.A_2, dotMomenta.A_3);
   aux[iM1][j][k].subtract_offset_drho_dA(1,0,0, dotMomenta.A_1, dotMomenta.A_2, dotMomenta.A_3);
   aux[iM1][j][kM1].subtract_offset_drho_dA(1,0,1, dotMomenta.A_1, dotMomenta.A_2, dotMomenta.A_3);
   aux[iM1][j][kM2].subtract_offset_drho_dA(1,0,2, dotMomenta.A_1, dotMomenta.A_2, dotMomenta.A_3);
   aux[iM1][jM1][kP2].subtract_offset_drho_dA(1,1,-2, dotMomenta.A_1, dotMomenta.A_2, dotMomenta.A_3);
   aux[iM1][jM1][kP1].subtract_offset_drho_dA(1,1,-1, dotMomenta.A_1, dotMomenta.A_2, dotMomenta.A_3);
   aux[iM1][jM1][k].subtract_offset_drho_dA(1,1,0, dotMomenta.A_1, dotMomenta.A_2, dotMomenta.A_3);
   aux[iM1][jM1][kM1].subtract_offset_drho_dA(1,1,1, dotMomenta.A_1, dotMomenta.A_2, dotMomenta.A_3);
   aux[iM1][jM1][kM2].subtract_offset_drho_dA(1,1,2, dotMomenta.A_1, dotMomenta.A_2, dotMomenta.A_3);
   aux[iM1][jM2][kP1].subtract_offset_drho_dA(1,2,-1, dotMomenta.A_1, dotMomenta.A_2, dotMomenta.A_3);
   aux[iM1][jM2][k].subtract_offset_drho_dA(1,2,0, dotMomenta.A_1, dotMomenta.A_2, dotMomenta.A_3);
   aux[iM1][jM2][kM1].subtract_offset_drho_dA(1,2,1, dotMomenta.A_1, dotMomenta.A_2, dotMomenta.A_3);

   aux[iM2][jP1][kP1].subtract_offset_drho_dA(2,-1,-1, dotMomenta.A_1, dotMomenta.A_2, dotMomenta.A_3);
   aux[iM2][jP1][k].subtract_offset_drho_dA(2,-1,0, dotMomenta.A_1, dotMomenta.A_2, dotMomenta.A_3);
   aux[iM2][jP1][kM1].subtract_offset_drho_dA(2,-1,1, dotMomenta.A_1, dotMomenta.A_2, dotMomenta.A_3);
   aux[iM2][j][kP1].subtract_offset_drho_dA(2,0,-1, dotMomenta.A_1, dotMomenta.A_2, dotMomenta.A_3);
   aux[iM2][j][k].subtract_offset_drho_dA(2,0,0, dotMomenta.A_1, dotMomenta.A_2, dotMomenta.A_3);
   aux[iM2][j][kM1].subtract_offset_drho_dA(2,0,1, dotMomenta.A_1, dotMomenta.A_2, dotMomenta.A_3);
   aux[iM2][jM1][kP1].subtract_offset_drho_dA(2,1,-1, dotMomenta.A_1, dotMomenta.A_2, dotMomenta.A_3);
   aux[iM2][jM1][k].subtract_offset_drho_dA(2,1,0, dotMomenta.A_1, dotMomenta.A_2, dotMomenta.A_3);
   aux[iM2][jM1][kM1].subtract_offset_drho_dA(2,1,1, dotMomenta.A_1, dotMomenta.A_2, dotMomenta.A_3);

   aux[i][j][k].subtract_drho_static_dA(&(fields[indx][i][j][k]), dotMomenta.A_1, dotMomenta.A_2, dotMomenta.A_3);
}

void TLattice::copy2private(int indx, int i, int j, int k, TPrivateFields privateFields, TPrivateAux privateAux)
/* Copies contents of global fields[indx] and aux matrices around (i,j,k) to privateFields and privateAux. 
     Overkill since only the inner [2..6][2..6][2..6] cube of privateAux is needed, and the corners are not needed in either privateFields or privateAux.
*/
{
   int ii[9], jj[9], kk[9];
   
   ii[4] = i;
   jj[4] = j;
   kk[4] = k;

   if (0 == i) ii[3] = nSideM1; else ii[3] = i - 1;
   if (0 == j) jj[3] = nSideM1; else jj[3] = j - 1;
   if (0 == k) kk[3] = nSideM1; else kk[3] = k - 1;

   if (0 == ii[3]) ii[2] = nSideM1; else ii[2] = ii[3] - 1;
   if (0 == jj[3]) jj[2] = nSideM1; else jj[2] = jj[3] - 1;
   if (0 == kk[3]) kk[2] = nSideM1; else kk[2] = kk[3] - 1;

   if (0 == ii[2]) ii[1] = nSideM1; else ii[1] = ii[2] - 1;
   if (0 == jj[2]) jj[1] = nSideM1; else jj[1] = jj[2] - 1;
   if (0 == kk[2]) kk[1] = nSideM1; else kk[1] = kk[2] - 1;

   if (0 == ii[1]) ii[0] = nSideM1; else ii[0] = ii[1] - 1;
   if (0 == jj[1]) jj[0] = nSideM1; else jj[0] = jj[1] - 1;
   if (0 == kk[1]) kk[0] = nSideM1; else kk[0] = kk[1] - 1;

   if (nSideM1 == i) ii[5] = 0; else ii[5] = i + 1;
   if (nSideM1 == j) jj[5] = 0; else jj[5] = j + 1;
   if (nSideM1 == k) kk[5] = 0; else kk[5] = k + 1;

   if (nSideM1 == ii[5]) ii[6] = 0; else ii[6] = ii[5] + 1;
   if (nSideM1 == jj[5]) jj[6] = 0; else jj[6] = jj[5] + 1;
   if (nSideM1 == kk[5]) kk[6] = 0; else kk[6] = kk[5] + 1;

   if (nSideM1 == ii[6]) ii[7] = 0; else ii[7] = ii[6] + 1;
   if (nSideM1 == jj[6]) jj[7] = 0; else jj[7] = jj[6] + 1;
   if (nSideM1 == kk[6]) kk[7] = 0; else kk[7] = kk[6] + 1;

   if (nSideM1 == ii[7]) ii[8] = 0; else ii[8] = ii[7] + 1;
   if (nSideM1 == jj[7]) jj[8] = 0; else jj[8] = jj[7] + 1;
   if (nSideM1 == kk[7]) kk[8] = 0; else kk[8] = kk[7] + 1;
   
   for (int toi = 0; toi < 9; toi++)
   {
       int fromi = ii[toi];
       
       for (int toj = 0; toj < 9; toj++)
       {
           int fromj = jj[toj];
       
           for (int tok = 0; tok < 9; tok++)
           {
              int fromk = kk[tok];
           
              privateFields[toi][toj][tok] = fields[indx][fromi][fromj][fromk];
              privateAux[toi][toj][tok] = aux[fromi][fromj][fromk];
           }
       }
   }
}

void TLattice::copy2affected_private_derivatives(int i, int j, int k, int component, TPrivateAux privateAux)
/* Copies contents of derivatives of component from global aux matrix around (i,j,k) to privateAux. */
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

   privateAux[2][3][4].set_derivative(0, component, aux[iM2][jM1][k].get_derivative(0, component));
   privateAux[2][4][3].set_derivative(0, component, aux[iM2][j][kM1].get_derivative(0, component));
   privateAux[2][4][4].set_derivative(0, component, aux[iM2][j][k].get_derivative(0, component));
   privateAux[2][4][5].set_derivative(0, component, aux[iM2][j][kP1].get_derivative(0, component));
   privateAux[2][5][4].set_derivative(0, component, aux[iM2][jP1][k].get_derivative(0, component));
   privateAux[3][2][3].set_derivative(0, component, aux[iM1][jM2][kM1].get_derivative(0, component));
   privateAux[3][2][5].set_derivative(0, component, aux[iM1][jM2][kP1].get_derivative(0, component));
   privateAux[3][3][2].set_derivative(0, component, aux[iM1][jM1][kM2].get_derivative(0, component));
   privateAux[3][3][3].set_derivative(0, component, aux[iM1][jM1][kM1].get_derivative(0, component));
   privateAux[3][3][4].set_derivative(0, component, aux[iM1][jM1][k].get_derivative(0, component));
   privateAux[3][3][5].set_derivative(0, component, aux[iM1][jM1][kP1].get_derivative(0, component));
   privateAux[3][3][6].set_derivative(0, component, aux[iM1][jM1][kP2].get_derivative(0, component));
   privateAux[3][4][3].set_derivative(0, component, aux[iM1][j][kM1].get_derivative(0, component));
   privateAux[3][4][4].set_derivative(0, component, aux[iM1][j][k].get_derivative(0, component));
   privateAux[3][4][5].set_derivative(0, component, aux[iM1][j][kP1].get_derivative(0, component));
   privateAux[3][5][2].set_derivative(0, component, aux[iM1][jP1][kM2].get_derivative(0, component));
   privateAux[3][5][3].set_derivative(0, component, aux[iM1][jP1][kM1].get_derivative(0, component));
   privateAux[3][5][4].set_derivative(0, component, aux[iM1][jP1][k].get_derivative(0, component));
   privateAux[3][5][5].set_derivative(0, component, aux[iM1][jP1][kP1].get_derivative(0, component));
   privateAux[3][5][6].set_derivative(0, component, aux[iM1][jP1][kP2].get_derivative(0, component));
   privateAux[3][6][3].set_derivative(0, component, aux[iM1][jP2][kM1].get_derivative(0, component));
   privateAux[3][6][5].set_derivative(0, component, aux[iM1][jP2][kP1].get_derivative(0, component));
   privateAux[5][2][3].set_derivative(0, component, aux[iP1][jM2][kM1].get_derivative(0, component));
   privateAux[5][2][5].set_derivative(0, component, aux[iP1][jM2][kP1].get_derivative(0, component));
   privateAux[5][3][2].set_derivative(0, component, aux[iP1][jM1][kM2].get_derivative(0, component));
   privateAux[5][3][3].set_derivative(0, component, aux[iP1][jM1][kM1].get_derivative(0, component));
   privateAux[5][3][4].set_derivative(0, component, aux[iP1][jM1][k].get_derivative(0, component));
   privateAux[5][3][5].set_derivative(0, component, aux[iP1][jM1][kP1].get_derivative(0, component));
   privateAux[5][3][6].set_derivative(0, component, aux[iP1][jM1][kP2].get_derivative(0, component));
   privateAux[5][4][3].set_derivative(0, component, aux[iP1][j][kM1].get_derivative(0, component));
   privateAux[5][4][4].set_derivative(0, component, aux[iP1][j][k].get_derivative(0, component));
   privateAux[5][4][5].set_derivative(0, component, aux[iP1][j][kP1].get_derivative(0, component));
   privateAux[5][5][2].set_derivative(0, component, aux[iP1][jP1][kM2].get_derivative(0, component));
   privateAux[5][5][3].set_derivative(0, component, aux[iP1][jP1][kM1].get_derivative(0, component));
   privateAux[5][5][4].set_derivative(0, component, aux[iP1][jP1][k].get_derivative(0, component));
   privateAux[5][5][5].set_derivative(0, component, aux[iP1][jP1][kP1].get_derivative(0, component));
   privateAux[5][5][6].set_derivative(0, component, aux[iP1][jP1][kP2].get_derivative(0, component));
   privateAux[5][6][3].set_derivative(0, component, aux[iP1][jP2][kM1].get_derivative(0, component));
   privateAux[5][6][5].set_derivative(0, component, aux[iP1][jP2][kP1].get_derivative(0, component));
   privateAux[6][3][4].set_derivative(0, component, aux[iP2][jM1][k].get_derivative(0, component));
   privateAux[6][4][3].set_derivative(0, component, aux[iP2][j][kM1].get_derivative(0, component));
   privateAux[6][4][4].set_derivative(0, component, aux[iP2][j][k].get_derivative(0, component));
   privateAux[6][4][5].set_derivative(0, component, aux[iP2][j][kP1].get_derivative(0, component));
   privateAux[6][5][4].set_derivative(0, component, aux[iP2][jP1][k].get_derivative(0, component));

   privateAux[2][3][3].set_derivative(1, component, aux[iM2][jM1][kM1].get_derivative(1, component));
   privateAux[2][3][5].set_derivative(1, component, aux[iM2][jM1][kP1].get_derivative(1, component));
   privateAux[2][5][3].set_derivative(1, component, aux[iM2][jP1][kM1].get_derivative(1, component));
   privateAux[2][5][5].set_derivative(1, component, aux[iM2][jP1][kP1].get_derivative(1, component));
   privateAux[3][2][4].set_derivative(1, component, aux[iM1][jM2][k].get_derivative(1, component));
   privateAux[3][3][2].set_derivative(1, component, aux[iM1][jM1][kM2].get_derivative(1, component));
   privateAux[3][3][3].set_derivative(1, component, aux[iM1][jM1][kM1].get_derivative(1, component));
   privateAux[3][3][4].set_derivative(1, component, aux[iM1][jM1][k].get_derivative(1, component));
   privateAux[3][3][5].set_derivative(1, component, aux[iM1][jM1][kP1].get_derivative(1, component));
   privateAux[3][3][6].set_derivative(1, component, aux[iM1][jM1][kP2].get_derivative(1, component));
   privateAux[3][5][2].set_derivative(1, component, aux[iM1][jP1][kM2].get_derivative(1, component));
   privateAux[3][5][3].set_derivative(1, component, aux[iM1][jP1][kM1].get_derivative(1, component));
   privateAux[3][5][4].set_derivative(1, component, aux[iM1][jP1][k].get_derivative(1, component));
   privateAux[3][5][5].set_derivative(1, component, aux[iM1][jP1][kP1].get_derivative(1, component));
   privateAux[3][5][6].set_derivative(1, component, aux[iM1][jP1][kP2].get_derivative(1, component));
   privateAux[3][6][4].set_derivative(1, component, aux[iM1][jP2][k].get_derivative(1, component));
   privateAux[4][2][3].set_derivative(1, component, aux[i][jM2][kM1].get_derivative(1, component));
   privateAux[4][2][4].set_derivative(1, component, aux[i][jM2][k].get_derivative(1, component));
   privateAux[4][2][5].set_derivative(1, component, aux[i][jM2][kP1].get_derivative(1, component));
   privateAux[4][3][3].set_derivative(1, component, aux[i][jM1][kM1].get_derivative(1, component));
   privateAux[4][3][4].set_derivative(1, component, aux[i][jM1][k].get_derivative(1, component));
   privateAux[4][3][5].set_derivative(1, component, aux[i][jM1][kP1].get_derivative(1, component));
   privateAux[4][5][3].set_derivative(1, component, aux[i][jP1][kM1].get_derivative(1, component));
   privateAux[4][5][4].set_derivative(1, component, aux[i][jP1][k].get_derivative(1, component));
   privateAux[4][5][5].set_derivative(1, component, aux[i][jP1][kP1].get_derivative(1, component));
   privateAux[4][6][3].set_derivative(1, component, aux[i][jP2][kM1].get_derivative(1, component));
   privateAux[4][6][4].set_derivative(1, component, aux[i][jP2][k].get_derivative(1, component));
   privateAux[4][6][5].set_derivative(1, component, aux[i][jP2][kP1].get_derivative(1, component));
   privateAux[5][2][4].set_derivative(1, component, aux[iP1][jM2][k].get_derivative(1, component));
   privateAux[5][3][2].set_derivative(1, component, aux[iP1][jM1][kM2].get_derivative(1, component));
   privateAux[5][3][3].set_derivative(1, component, aux[iP1][jM1][kM1].get_derivative(1, component));
   privateAux[5][3][4].set_derivative(1, component, aux[iP1][jM1][k].get_derivative(1, component));
   privateAux[5][3][5].set_derivative(1, component, aux[iP1][jM1][kP1].get_derivative(1, component));
   privateAux[5][3][6].set_derivative(1, component, aux[iP1][jM1][kP2].get_derivative(1, component));
   privateAux[5][5][2].set_derivative(1, component, aux[iP1][jP1][kM2].get_derivative(1, component));
   privateAux[5][5][3].set_derivative(1, component, aux[iP1][jP1][kM1].get_derivative(1, component));
   privateAux[5][5][4].set_derivative(1, component, aux[iP1][jP1][k].get_derivative(1, component));
   privateAux[5][5][5].set_derivative(1, component, aux[iP1][jP1][kP1].get_derivative(1, component));
   privateAux[5][5][6].set_derivative(1, component, aux[iP1][jP1][kP2].get_derivative(1, component));
   privateAux[5][6][4].set_derivative(1, component, aux[iP1][jP2][k].get_derivative(1, component));
   privateAux[6][3][3].set_derivative(1, component, aux[iP2][jM1][kM1].get_derivative(1, component));
   privateAux[6][3][5].set_derivative(1, component, aux[iP2][jM1][kP1].get_derivative(1, component));
   privateAux[6][5][3].set_derivative(1, component, aux[iP2][jP1][kM1].get_derivative(1, component));
   privateAux[6][5][5].set_derivative(1, component, aux[iP2][jP1][kP1].get_derivative(1, component));

   privateAux[2][3][3].set_derivative(2, component, aux[iM2][jM1][kM1].get_derivative(2, component));
   privateAux[2][3][5].set_derivative(2, component, aux[iM2][jM1][kP1].get_derivative(2, component));
   privateAux[2][5][3].set_derivative(2, component, aux[iM2][jP1][kM1].get_derivative(2, component));
   privateAux[2][5][5].set_derivative(2, component, aux[iM2][jP1][kP1].get_derivative(2, component));
   privateAux[3][2][3].set_derivative(2, component, aux[iM1][jM2][kM1].get_derivative(2, component));
   privateAux[3][2][5].set_derivative(2, component, aux[iM1][jM2][kP1].get_derivative(2, component));
   privateAux[3][3][3].set_derivative(2, component, aux[iM1][jM1][kM1].get_derivative(2, component));
   privateAux[3][3][5].set_derivative(2, component, aux[iM1][jM1][kP1].get_derivative(2, component));
   privateAux[3][4][2].set_derivative(2, component, aux[iM1][j][kM2].get_derivative(2, component));
   privateAux[3][4][3].set_derivative(2, component, aux[iM1][j][kM1].get_derivative(2, component));
   privateAux[3][4][5].set_derivative(2, component, aux[iM1][j][kP1].get_derivative(2, component));
   privateAux[3][4][6].set_derivative(2, component, aux[iM1][j][kP2].get_derivative(2, component));
   privateAux[3][5][3].set_derivative(2, component, aux[iM1][jP1][kM1].get_derivative(2, component));
   privateAux[3][5][5].set_derivative(2, component, aux[iM1][jP1][kP1].get_derivative(2, component));
   privateAux[3][6][3].set_derivative(2, component, aux[iM1][jP2][kM1].get_derivative(2, component));
   privateAux[3][6][5].set_derivative(2, component, aux[iM1][jP2][kP1].get_derivative(2, component));
   privateAux[4][3][2].set_derivative(2, component, aux[i][jM1][kM2].get_derivative(2, component));
   privateAux[4][3][3].set_derivative(2, component, aux[i][jM1][kM1].get_derivative(2, component));
   privateAux[4][3][5].set_derivative(2, component, aux[i][jM1][kP1].get_derivative(2, component));
   privateAux[4][3][6].set_derivative(2, component, aux[i][jM1][kP2].get_derivative(2, component));
   privateAux[4][4][2].set_derivative(2, component, aux[i][j][kM2].get_derivative(2, component));
   privateAux[4][4][3].set_derivative(2, component, aux[i][j][kM1].get_derivative(2, component));
   privateAux[4][4][5].set_derivative(2, component, aux[i][j][kP1].get_derivative(2, component));
   privateAux[4][4][6].set_derivative(2, component, aux[i][j][kP2].get_derivative(2, component));
   privateAux[4][5][2].set_derivative(2, component, aux[i][jP1][kM2].get_derivative(2, component));
   privateAux[4][5][3].set_derivative(2, component, aux[i][jP1][kM1].get_derivative(2, component));
   privateAux[4][5][5].set_derivative(2, component, aux[i][jP1][kP1].get_derivative(2, component));
   privateAux[4][5][6].set_derivative(2, component, aux[i][jP1][kP2].get_derivative(2, component));
   privateAux[5][2][3].set_derivative(2, component, aux[iP1][jM2][kM1].get_derivative(2, component));
   privateAux[5][2][5].set_derivative(2, component, aux[iP1][jM2][kP1].get_derivative(2, component));
   privateAux[5][3][3].set_derivative(2, component, aux[iP1][jM1][kM1].get_derivative(2, component));
   privateAux[5][3][5].set_derivative(2, component, aux[iP1][jM1][kP1].get_derivative(2, component));
   privateAux[5][4][2].set_derivative(2, component, aux[iP1][j][kM2].get_derivative(2, component));
   privateAux[5][4][3].set_derivative(2, component, aux[iP1][j][kM1].get_derivative(2, component));
   privateAux[5][4][5].set_derivative(2, component, aux[iP1][j][kP1].get_derivative(2, component));
   privateAux[5][4][6].set_derivative(2, component, aux[iP1][j][kP2].get_derivative(2, component));
   privateAux[5][5][3].set_derivative(2, component, aux[iP1][jP1][kM1].get_derivative(2, component));
   privateAux[5][5][5].set_derivative(2, component, aux[iP1][jP1][kP1].get_derivative(2, component));
   privateAux[5][6][3].set_derivative(2, component, aux[iP1][jP2][kM1].get_derivative(2, component));
   privateAux[5][6][5].set_derivative(2, component, aux[iP1][jP2][kP1].get_derivative(2, component));
   privateAux[6][3][3].set_derivative(2, component, aux[iP2][jM1][kM1].get_derivative(2, component));
   privateAux[6][3][5].set_derivative(2, component, aux[iP2][jM1][kP1].get_derivative(2, component));
   privateAux[6][5][3].set_derivative(2, component, aux[iP2][jP1][kM1].get_derivative(2, component));
   privateAux[6][5][5].set_derivative(2, component, aux[iP2][jP1][kP1].get_derivative(2, component));
}

void TLattice::compute_affected_private_derivatives(int component, TPrivateFields privateFields, TPrivateAux privateAux)
/* Computes spatial derivatives of privateFields[9][9][9] affected by changing the value of field component at (4,4,4) and stores them in privateAux [2..6][2..6][2..6]. 
     Components 0 thru 2 are th_1 thru th_3;
     components 3 thru 5 are A_1 thru A_3. */
{
   privateAux[2][3][4].compute_derivative(privateFields, 2, 3, 4, 0, component);
   privateAux[2][4][3].compute_derivative(privateFields, 2, 4, 3, 0, component);
   privateAux[2][4][4].compute_derivative(privateFields, 2, 4, 4, 0, component);
   privateAux[2][4][5].compute_derivative(privateFields, 2, 4, 5, 0, component);
   privateAux[2][5][4].compute_derivative(privateFields, 2, 5, 4, 0, component);
   privateAux[3][2][3].compute_derivative(privateFields, 3, 2, 3, 0, component);
   privateAux[3][2][5].compute_derivative(privateFields, 3, 2, 5, 0, component);
   privateAux[3][3][2].compute_derivative(privateFields, 3, 3, 2, 0, component);
   privateAux[3][3][3].compute_derivative(privateFields, 3, 3, 3, 0, component);
   privateAux[3][3][4].compute_derivative(privateFields, 3, 3, 4, 0, component);
   privateAux[3][3][5].compute_derivative(privateFields, 3, 3, 5, 0, component);
   privateAux[3][3][6].compute_derivative(privateFields, 3, 3, 6, 0, component);
   privateAux[3][4][3].compute_derivative(privateFields, 3, 4, 3, 0, component);
   privateAux[3][4][4].compute_derivative(privateFields, 3, 4, 4, 0, component);
   privateAux[3][4][5].compute_derivative(privateFields, 3, 4, 5, 0, component);
   privateAux[3][5][2].compute_derivative(privateFields, 3, 5, 2, 0, component);
   privateAux[3][5][3].compute_derivative(privateFields, 3, 5, 3, 0, component);
   privateAux[3][5][4].compute_derivative(privateFields, 3, 5, 4, 0, component);
   privateAux[3][5][5].compute_derivative(privateFields, 3, 5, 5, 0, component);
   privateAux[3][5][6].compute_derivative(privateFields, 3, 5, 6, 0, component);
   privateAux[3][6][3].compute_derivative(privateFields, 3, 6, 3, 0, component);
   privateAux[3][6][5].compute_derivative(privateFields, 3, 6, 5, 0, component);
   privateAux[5][2][3].compute_derivative(privateFields, 5, 2, 3, 0, component);
   privateAux[5][2][5].compute_derivative(privateFields, 5, 2, 5, 0, component);
   privateAux[5][3][2].compute_derivative(privateFields, 5, 3, 2, 0, component);
   privateAux[5][3][3].compute_derivative(privateFields, 5, 3, 3, 0, component);
   privateAux[5][3][4].compute_derivative(privateFields, 5, 3, 4, 0, component);
   privateAux[5][3][5].compute_derivative(privateFields, 5, 3, 5, 0, component);
   privateAux[5][3][6].compute_derivative(privateFields, 5, 3, 6, 0, component);
   privateAux[5][4][3].compute_derivative(privateFields, 5, 4, 3, 0, component);
   privateAux[5][4][4].compute_derivative(privateFields, 5, 4, 4, 0, component);
   privateAux[5][4][5].compute_derivative(privateFields, 5, 4, 5, 0, component);
   privateAux[5][5][2].compute_derivative(privateFields, 5, 5, 2, 0, component);
   privateAux[5][5][3].compute_derivative(privateFields, 5, 5, 3, 0, component);
   privateAux[5][5][4].compute_derivative(privateFields, 5, 5, 4, 0, component);
   privateAux[5][5][5].compute_derivative(privateFields, 5, 5, 5, 0, component);
   privateAux[5][5][6].compute_derivative(privateFields, 5, 5, 6, 0, component);
   privateAux[5][6][3].compute_derivative(privateFields, 5, 6, 3, 0, component);
   privateAux[5][6][5].compute_derivative(privateFields, 5, 6, 5, 0, component);
   privateAux[6][3][4].compute_derivative(privateFields, 6, 3, 4, 0, component);
   privateAux[6][4][3].compute_derivative(privateFields, 6, 4, 3, 0, component);
   privateAux[6][4][4].compute_derivative(privateFields, 6, 4, 4, 0, component);
   privateAux[6][4][5].compute_derivative(privateFields, 6, 4, 5, 0, component);
   privateAux[6][5][4].compute_derivative(privateFields, 6, 5, 4, 0, component);

   privateAux[2][3][3].compute_derivative(privateFields, 2, 3, 3, 1, component);
   privateAux[2][3][5].compute_derivative(privateFields, 2, 3, 5, 1, component);
   privateAux[2][5][3].compute_derivative(privateFields, 2, 5, 3, 1, component);
   privateAux[2][5][5].compute_derivative(privateFields, 2, 5, 5, 1, component);
   privateAux[3][2][4].compute_derivative(privateFields, 3, 2, 4, 1, component);
   privateAux[3][3][2].compute_derivative(privateFields, 3, 3, 2, 1, component);
   privateAux[3][3][3].compute_derivative(privateFields, 3, 3, 3, 1, component);
   privateAux[3][3][4].compute_derivative(privateFields, 3, 3, 4, 1, component);
   privateAux[3][3][5].compute_derivative(privateFields, 3, 3, 5, 1, component);
   privateAux[3][3][6].compute_derivative(privateFields, 3, 3, 6, 1, component);
   privateAux[3][5][2].compute_derivative(privateFields, 3, 5, 2, 1, component);
   privateAux[3][5][3].compute_derivative(privateFields, 3, 5, 3, 1, component);
   privateAux[3][5][4].compute_derivative(privateFields, 3, 5, 4, 1, component);
   privateAux[3][5][5].compute_derivative(privateFields, 3, 5, 5, 1, component);
   privateAux[3][5][6].compute_derivative(privateFields, 3, 5, 6, 1, component);
   privateAux[3][6][4].compute_derivative(privateFields, 3, 6, 4, 1, component);
   privateAux[4][2][3].compute_derivative(privateFields, 4, 2, 3, 1, component);
   privateAux[4][2][4].compute_derivative(privateFields, 4, 2, 4, 1, component);
   privateAux[4][2][5].compute_derivative(privateFields, 4, 2, 5, 1, component);
   privateAux[4][3][3].compute_derivative(privateFields, 4, 3, 3, 1, component);
   privateAux[4][3][4].compute_derivative(privateFields, 4, 3, 4, 1, component);
   privateAux[4][3][5].compute_derivative(privateFields, 4, 3, 5, 1, component);
   privateAux[4][5][3].compute_derivative(privateFields, 4, 5, 3, 1, component);
   privateAux[4][5][4].compute_derivative(privateFields, 4, 5, 4, 1, component);
   privateAux[4][5][5].compute_derivative(privateFields, 4, 5, 5, 1, component);
   privateAux[4][6][3].compute_derivative(privateFields, 4, 6, 3, 1, component);
   privateAux[4][6][4].compute_derivative(privateFields, 4, 6, 4, 1, component);
   privateAux[4][6][5].compute_derivative(privateFields, 4, 6, 5, 1, component);
   privateAux[5][2][4].compute_derivative(privateFields, 5, 2, 4, 1, component);
   privateAux[5][3][2].compute_derivative(privateFields, 5, 3, 2, 1, component);
   privateAux[5][3][3].compute_derivative(privateFields, 5, 3, 3, 1, component);
   privateAux[5][3][4].compute_derivative(privateFields, 5, 3, 4, 1, component);
   privateAux[5][3][5].compute_derivative(privateFields, 5, 3, 5, 1, component);
   privateAux[5][3][6].compute_derivative(privateFields, 5, 3, 6, 1, component);
   privateAux[5][5][2].compute_derivative(privateFields, 5, 5, 2, 1, component);
   privateAux[5][5][3].compute_derivative(privateFields, 5, 5, 3, 1, component);
   privateAux[5][5][4].compute_derivative(privateFields, 5, 5, 4, 1, component);
   privateAux[5][5][5].compute_derivative(privateFields, 5, 5, 5, 1, component);
   privateAux[5][5][6].compute_derivative(privateFields, 5, 5, 6, 1, component);
   privateAux[5][6][4].compute_derivative(privateFields, 5, 6, 4, 1, component);
   privateAux[6][3][3].compute_derivative(privateFields, 6, 3, 3, 1, component);
   privateAux[6][3][5].compute_derivative(privateFields, 6, 3, 5, 1, component);
   privateAux[6][5][3].compute_derivative(privateFields, 6, 5, 3, 1, component);
   privateAux[6][5][5].compute_derivative(privateFields, 6, 5, 5, 1, component);

   privateAux[2][3][3].compute_derivative(privateFields, 2, 3, 3, 2, component);
   privateAux[2][3][5].compute_derivative(privateFields, 2, 3, 5, 2, component);
   privateAux[2][5][3].compute_derivative(privateFields, 2, 5, 3, 2, component);
   privateAux[2][5][5].compute_derivative(privateFields, 2, 5, 5, 2, component);
   privateAux[3][2][3].compute_derivative(privateFields, 3, 2, 3, 2, component);
   privateAux[3][2][5].compute_derivative(privateFields, 3, 2, 5, 2, component);
   privateAux[3][3][3].compute_derivative(privateFields, 3, 3, 3, 2, component);
   privateAux[3][3][5].compute_derivative(privateFields, 3, 3, 5, 2, component);
   privateAux[3][4][2].compute_derivative(privateFields, 3, 4, 2, 2, component);
   privateAux[3][4][3].compute_derivative(privateFields, 3, 4, 3, 2, component);
   privateAux[3][4][5].compute_derivative(privateFields, 3, 4, 5, 2, component);
   privateAux[3][4][6].compute_derivative(privateFields, 3, 4, 6, 2, component);
   privateAux[3][5][3].compute_derivative(privateFields, 3, 5, 3, 2, component);
   privateAux[3][5][5].compute_derivative(privateFields, 3, 5, 5, 2, component);
   privateAux[3][6][3].compute_derivative(privateFields, 3, 6, 3, 2, component);
   privateAux[3][6][5].compute_derivative(privateFields, 3, 6, 5, 2, component);
   privateAux[4][3][2].compute_derivative(privateFields, 4, 3, 2, 2, component);
   privateAux[4][3][3].compute_derivative(privateFields, 4, 3, 3, 2, component);
   privateAux[4][3][5].compute_derivative(privateFields, 4, 3, 5, 2, component);
   privateAux[4][3][6].compute_derivative(privateFields, 4, 3, 6, 2, component);
   privateAux[4][4][2].compute_derivative(privateFields, 4, 4, 2, 2, component);
   privateAux[4][4][3].compute_derivative(privateFields, 4, 4, 3, 2, component);
   privateAux[4][4][5].compute_derivative(privateFields, 4, 4, 5, 2, component);
   privateAux[4][4][6].compute_derivative(privateFields, 4, 4, 6, 2, component);
   privateAux[4][5][2].compute_derivative(privateFields, 4, 5, 2, 2, component);
   privateAux[4][5][3].compute_derivative(privateFields, 4, 5, 3, 2, component);
   privateAux[4][5][5].compute_derivative(privateFields, 4, 5, 5, 2, component);
   privateAux[4][5][6].compute_derivative(privateFields, 4, 5, 6, 2, component);
   privateAux[5][2][3].compute_derivative(privateFields, 5, 2, 3, 2, component);
   privateAux[5][2][5].compute_derivative(privateFields, 5, 2, 5, 2, component);
   privateAux[5][3][3].compute_derivative(privateFields, 5, 3, 3, 2, component);
   privateAux[5][3][5].compute_derivative(privateFields, 5, 3, 5, 2, component);
   privateAux[5][4][2].compute_derivative(privateFields, 5, 4, 2, 2, component);
   privateAux[5][4][3].compute_derivative(privateFields, 5, 4, 3, 2, component);
   privateAux[5][4][5].compute_derivative(privateFields, 5, 4, 5, 2, component);
   privateAux[5][4][6].compute_derivative(privateFields, 5, 4, 6, 2, component);
   privateAux[5][5][3].compute_derivative(privateFields, 5, 5, 3, 2, component);
   privateAux[5][5][5].compute_derivative(privateFields, 5, 5, 5, 2, component);
   privateAux[5][6][3].compute_derivative(privateFields, 5, 6, 3, 2, component);
   privateAux[5][6][5].compute_derivative(privateFields, 5, 6, 5, 2, component);
   privateAux[6][3][3].compute_derivative(privateFields, 6, 3, 3, 2, component);
   privateAux[6][3][5].compute_derivative(privateFields, 6, 3, 5, 2, component);
   privateAux[6][5][3].compute_derivative(privateFields, 6, 5, 3, 2, component);
   privateAux[6][5][5].compute_derivative(privateFields, 6, 5, 5, 2, component);
}

void TLattice::init_aux_static()
/* Almost init_aux(), except invM is not computed (since it's not needed in the static case) */
{
#ifdef _OPENMP
#pragma omp parallel
{
#pragma omp for
#endif
   for (int i = 0; i <= nSideM1; i++)
       for (int j = 0; j <= nSideM1; j++)
           for (int k = 0; k <= nSideM1; k++)
           {           
              aux[i][j][k].compute_hG(fields[active][i][j][k].th_1,
                                      fields[active][i][j][k].th_2,
                                      fields[active][i][j][k].th_3);
              aux[i][j][k].compute_derivatives(fields[active], i, j, k, nSideM1);
           }
#ifdef _OPENMP
#pragma omp barrier
}
#endif
}

double TLattice::sum_rho_static()
/*  IMPORTANT: Assumes all h, G matrices and derivatives are ready for use.
     Multiply by dx^3/dx^2 = dx for total energy */
{
   double result = 0.0;

#ifdef _OPENMP
#pragma omp parallel default(shared) reduction(+: result)
{
#pragma omp for
#endif
   for (int i = 0; i <= nSideM1; i++)
       for (int j = 0; j <= nSideM1; j++)
           for (int k = 0; k <= nSideM1; k++)
           {
              TRho rho;
              result += aux[i][j][k].get_rho_static(&(fields[active][i][j][k]), rho);
           }
#ifdef _OPENMP
#pragma omp barrier
}
#endif

   return(result);
}

double TLattice::affected_private_static_rho(TPrivateFields privateFields, TPrivateAux privateAux)
/* Sums part of total static rho affected by varying privateFields at (4, 4, 4).
     IMPORTANT: Assumes all (center and off-center) h, G matrices and derivatives of off-center cells are ready for use in privateAux.
     Also, does NOT divide by number of fields summed over, so if absolute (rather than relative) magnitude of rho is needed, caller must do that;
     and does NOT divide by dx^2 from derivatives  */
{
   double rhoSum;
   TRho rho;

   // Sum up static rho of off-center fields affected by varying fields at (4,4,4).
   // IMPORTANT: requires derivatives and h, G matrices of off-center cells to be ready in privateAux!

   rhoSum  = privateAux[6][5][5].get_rho_static(&(privateFields[6][5][5]), rho);
   rhoSum += privateAux[6][5][4].get_rho_static(&(privateFields[6][5][4]), rho);
   rhoSum += privateAux[6][5][3].get_rho_static(&(privateFields[6][5][3]), rho);
   rhoSum += privateAux[6][4][5].get_rho_static(&(privateFields[6][4][5]), rho);
   rhoSum += privateAux[6][4][4].get_rho_static(&(privateFields[6][4][4]), rho);
   rhoSum += privateAux[6][4][3].get_rho_static(&(privateFields[6][4][3]), rho);
   rhoSum += privateAux[6][3][5].get_rho_static(&(privateFields[6][3][5]), rho);
   rhoSum += privateAux[6][3][4].get_rho_static(&(privateFields[6][3][4]), rho);
   rhoSum += privateAux[6][3][3].get_rho_static(&(privateFields[6][3][3]), rho);

   rhoSum += privateAux[5][6][5].get_rho_static(&(privateFields[5][6][5]), rho);
   rhoSum += privateAux[5][6][4].get_rho_static(&(privateFields[5][6][4]), rho);
   rhoSum += privateAux[5][6][3].get_rho_static(&(privateFields[5][6][3]), rho);
   rhoSum += privateAux[5][5][6].get_rho_static(&(privateFields[5][5][6]), rho);
   rhoSum += privateAux[5][5][5].get_rho_static(&(privateFields[5][5][5]), rho);
   rhoSum += privateAux[5][5][4].get_rho_static(&(privateFields[5][5][4]), rho);
   rhoSum += privateAux[5][5][3].get_rho_static(&(privateFields[5][5][3]), rho);
   rhoSum += privateAux[5][5][2].get_rho_static(&(privateFields[5][5][2]), rho);
   rhoSum += privateAux[5][4][6].get_rho_static(&(privateFields[5][4][6]), rho);
   rhoSum += privateAux[5][4][5].get_rho_static(&(privateFields[5][4][5]), rho);
   rhoSum += privateAux[5][4][4].get_rho_static(&(privateFields[5][4][4]), rho);
   rhoSum += privateAux[5][4][3].get_rho_static(&(privateFields[5][4][3]), rho);
   rhoSum += privateAux[5][4][2].get_rho_static(&(privateFields[5][4][2]), rho);
   rhoSum += privateAux[5][3][6].get_rho_static(&(privateFields[5][3][6]), rho);
   rhoSum += privateAux[5][3][5].get_rho_static(&(privateFields[5][3][5]), rho);
   rhoSum += privateAux[5][3][4].get_rho_static(&(privateFields[5][3][4]), rho);
   rhoSum += privateAux[5][3][3].get_rho_static(&(privateFields[5][3][3]), rho);
   rhoSum += privateAux[5][3][2].get_rho_static(&(privateFields[5][3][2]), rho);
   rhoSum += privateAux[5][2][5].get_rho_static(&(privateFields[5][2][5]), rho);
   rhoSum += privateAux[5][2][4].get_rho_static(&(privateFields[5][2][4]), rho);
   rhoSum += privateAux[5][2][3].get_rho_static(&(privateFields[5][2][3]), rho);

   rhoSum += privateAux[4][6][5].get_rho_static(&(privateFields[4][6][5]), rho);
   rhoSum += privateAux[4][6][4].get_rho_static(&(privateFields[4][6][4]), rho);
   rhoSum += privateAux[4][6][3].get_rho_static(&(privateFields[4][6][3]), rho);
   rhoSum += privateAux[4][5][6].get_rho_static(&(privateFields[4][5][6]), rho);
   rhoSum += privateAux[4][5][5].get_rho_static(&(privateFields[4][5][5]), rho);
   rhoSum += privateAux[4][5][4].get_rho_static(&(privateFields[4][5][4]), rho);
   rhoSum += privateAux[4][5][3].get_rho_static(&(privateFields[4][5][3]), rho);
   rhoSum += privateAux[4][5][2].get_rho_static(&(privateFields[4][5][2]), rho);
   rhoSum += privateAux[4][4][6].get_rho_static(&(privateFields[4][4][6]), rho);
   rhoSum += privateAux[4][4][5].get_rho_static(&(privateFields[4][4][5]), rho);
   rhoSum += privateAux[4][4][3].get_rho_static(&(privateFields[4][4][3]), rho);
   rhoSum += privateAux[4][4][2].get_rho_static(&(privateFields[4][4][2]), rho);
   rhoSum += privateAux[4][3][6].get_rho_static(&(privateFields[4][3][6]), rho);
   rhoSum += privateAux[4][3][5].get_rho_static(&(privateFields[4][3][5]), rho);
   rhoSum += privateAux[4][3][4].get_rho_static(&(privateFields[4][3][4]), rho);
   rhoSum += privateAux[4][3][3].get_rho_static(&(privateFields[4][3][3]), rho);
   rhoSum += privateAux[4][3][2].get_rho_static(&(privateFields[4][3][2]), rho);
   rhoSum += privateAux[4][2][5].get_rho_static(&(privateFields[4][2][5]), rho);
   rhoSum += privateAux[4][2][4].get_rho_static(&(privateFields[4][2][4]), rho);
   rhoSum += privateAux[4][2][3].get_rho_static(&(privateFields[4][2][3]), rho);

   rhoSum += privateAux[3][6][5].get_rho_static(&(privateFields[3][6][5]), rho);
   rhoSum += privateAux[3][6][4].get_rho_static(&(privateFields[3][6][4]), rho);
   rhoSum += privateAux[3][6][3].get_rho_static(&(privateFields[3][6][3]), rho);
   rhoSum += privateAux[3][5][6].get_rho_static(&(privateFields[3][5][6]), rho);
   rhoSum += privateAux[3][5][5].get_rho_static(&(privateFields[3][5][5]), rho);
   rhoSum += privateAux[3][5][4].get_rho_static(&(privateFields[3][5][4]), rho);
   rhoSum += privateAux[3][5][3].get_rho_static(&(privateFields[3][5][3]), rho);
   rhoSum += privateAux[3][5][2].get_rho_static(&(privateFields[3][5][2]), rho);
   rhoSum += privateAux[3][4][6].get_rho_static(&(privateFields[3][4][6]), rho);
   rhoSum += privateAux[3][4][5].get_rho_static(&(privateFields[3][4][5]), rho);
   rhoSum += privateAux[3][4][4].get_rho_static(&(privateFields[3][4][4]), rho);
   rhoSum += privateAux[3][4][3].get_rho_static(&(privateFields[3][4][3]), rho);
   rhoSum += privateAux[3][4][2].get_rho_static(&(privateFields[3][4][2]), rho);
   rhoSum += privateAux[3][3][6].get_rho_static(&(privateFields[3][3][6]), rho);
   rhoSum += privateAux[3][3][5].get_rho_static(&(privateFields[3][3][5]), rho);
   rhoSum += privateAux[3][3][4].get_rho_static(&(privateFields[3][3][4]), rho);
   rhoSum += privateAux[3][3][3].get_rho_static(&(privateFields[3][3][3]), rho);
   rhoSum += privateAux[3][3][2].get_rho_static(&(privateFields[3][3][2]), rho);
   rhoSum += privateAux[3][2][5].get_rho_static(&(privateFields[3][2][5]), rho);
   rhoSum += privateAux[3][2][4].get_rho_static(&(privateFields[3][2][4]), rho);
   rhoSum += privateAux[3][2][3].get_rho_static(&(privateFields[3][2][3]), rho);

   rhoSum += privateAux[2][5][5].get_rho_static(&(privateFields[2][5][5]), rho);
   rhoSum += privateAux[2][5][4].get_rho_static(&(privateFields[2][5][4]), rho);
   rhoSum += privateAux[2][5][3].get_rho_static(&(privateFields[2][5][3]), rho);
   rhoSum += privateAux[2][4][5].get_rho_static(&(privateFields[2][4][5]), rho);
   rhoSum += privateAux[2][4][4].get_rho_static(&(privateFields[2][4][4]), rho);
   rhoSum += privateAux[2][4][3].get_rho_static(&(privateFields[2][4][3]), rho);
   rhoSum += privateAux[2][3][5].get_rho_static(&(privateFields[2][3][5]), rho);
   rhoSum += privateAux[2][3][4].get_rho_static(&(privateFields[2][3][4]), rho);
   rhoSum += privateAux[2][3][3].get_rho_static(&(privateFields[2][3][3]), rho);

   // Add static rho of central site. IMPORTANT: requires h, G matrices to be ready.

   rhoSum += privateAux[4][4][4].get_rho_static(&(privateFields[4][4][4]), rho);

   return(rhoSum);
}

double TLattice::affected_static_rho(int i, int j, int k)
/* Sums part of total static rho affected by varying fields at (i, j, k).
     IMPORTANT: Assumes all (center and off-center) h, G matrices and derivatives of off-center cells are ready for use.
     Also, does NOT divide by number of fields summed over, so if absolute (rather than relative) magnitude of rho is needed, caller must do that;
     and does NOT divide by dx^2 from derivatives  */
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

   double rhoSum;
   TRho rho;

   // Sum up static rho of off-center fields affected by varying fields at (i,j,k).
   // IMPORTANT: requires derivatives and h, G matrices of off-center cells to be ready!

   rhoSum  = aux[iP2][jP1][kP1].get_rho_static(&(fields[active][iP2][jP1][kP1]), rho);
   rhoSum += aux[iP2][jP1][k].get_rho_static(&(fields[active][iP2][jP1][k]), rho);
   rhoSum += aux[iP2][jP1][kM1].get_rho_static(&(fields[active][iP2][jP1][kM1]), rho);
   rhoSum += aux[iP2][j][kP1].get_rho_static(&(fields[active][iP2][j][kP1]), rho);
   rhoSum += aux[iP2][j][k].get_rho_static(&(fields[active][iP2][j][k]), rho);
   rhoSum += aux[iP2][j][kM1].get_rho_static(&(fields[active][iP2][j][kM1]), rho);
   rhoSum += aux[iP2][jM1][kP1].get_rho_static(&(fields[active][iP2][jM1][kP1]), rho);
   rhoSum += aux[iP2][jM1][k].get_rho_static(&(fields[active][iP2][jM1][k]), rho);
   rhoSum += aux[iP2][jM1][kM1].get_rho_static(&(fields[active][iP2][jM1][kM1]), rho);

   rhoSum += aux[iP1][jP2][kP1].get_rho_static(&(fields[active][iP1][jP2][kP1]), rho);
   rhoSum += aux[iP1][jP2][k].get_rho_static(&(fields[active][iP1][jP2][k]), rho);
   rhoSum += aux[iP1][jP2][kM1].get_rho_static(&(fields[active][iP1][jP2][kM1]), rho);
   rhoSum += aux[iP1][jP1][kP2].get_rho_static(&(fields[active][iP1][jP1][kP2]), rho);
   rhoSum += aux[iP1][jP1][kP1].get_rho_static(&(fields[active][iP1][jP1][kP1]), rho);
   rhoSum += aux[iP1][jP1][k].get_rho_static(&(fields[active][iP1][jP1][k]), rho);
   rhoSum += aux[iP1][jP1][kM1].get_rho_static(&(fields[active][iP1][jP1][kM1]), rho);
   rhoSum += aux[iP1][jP1][kM2].get_rho_static(&(fields[active][iP1][jP1][kM2]), rho);
   rhoSum += aux[iP1][j][kP2].get_rho_static(&(fields[active][iP1][j][kP2]), rho);
   rhoSum += aux[iP1][j][kP1].get_rho_static(&(fields[active][iP1][j][kP1]), rho);
   rhoSum += aux[iP1][j][k].get_rho_static(&(fields[active][iP1][j][k]), rho);
   rhoSum += aux[iP1][j][kM1].get_rho_static(&(fields[active][iP1][j][kM1]), rho);
   rhoSum += aux[iP1][j][kM2].get_rho_static(&(fields[active][iP1][j][kM2]), rho);
   rhoSum += aux[iP1][jM1][kP2].get_rho_static(&(fields[active][iP1][jM1][kP2]), rho);
   rhoSum += aux[iP1][jM1][kP1].get_rho_static(&(fields[active][iP1][jM1][kP1]), rho);
   rhoSum += aux[iP1][jM1][k].get_rho_static(&(fields[active][iP1][jM1][k]), rho);
   rhoSum += aux[iP1][jM1][kM1].get_rho_static(&(fields[active][iP1][jM1][kM1]), rho);
   rhoSum += aux[iP1][jM1][kM2].get_rho_static(&(fields[active][iP1][jM1][kM2]), rho);
   rhoSum += aux[iP1][jM2][kP1].get_rho_static(&(fields[active][iP1][jM2][kP1]), rho);
   rhoSum += aux[iP1][jM2][k].get_rho_static(&(fields[active][iP1][jM2][k]), rho);
   rhoSum += aux[iP1][jM2][kM1].get_rho_static(&(fields[active][iP1][jM2][kM1]), rho);

   rhoSum += aux[i][jP2][kP1].get_rho_static(&(fields[active][i][jP2][kP1]), rho);
   rhoSum += aux[i][jP2][k].get_rho_static(&(fields[active][i][jP2][k]), rho);
   rhoSum += aux[i][jP2][kM1].get_rho_static(&(fields[active][i][jP2][kM1]), rho);
   rhoSum += aux[i][jP1][kP2].get_rho_static(&(fields[active][i][jP1][kP2]), rho);
   rhoSum += aux[i][jP1][kP1].get_rho_static(&(fields[active][i][jP1][kP1]), rho);
   rhoSum += aux[i][jP1][k].get_rho_static(&(fields[active][i][jP1][k]), rho);
   rhoSum += aux[i][jP1][kM1].get_rho_static(&(fields[active][i][jP1][kM1]), rho);
   rhoSum += aux[i][jP1][kM2].get_rho_static(&(fields[active][i][jP1][kM2]), rho);
   rhoSum += aux[i][j][kP2].get_rho_static(&(fields[active][i][j][kP2]), rho);
   rhoSum += aux[i][j][kP1].get_rho_static(&(fields[active][i][j][kP1]), rho);
   rhoSum += aux[i][j][kM1].get_rho_static(&(fields[active][i][j][kM1]), rho);
   rhoSum += aux[i][j][kM2].get_rho_static(&(fields[active][i][j][kM2]), rho);
   rhoSum += aux[i][jM1][kP2].get_rho_static(&(fields[active][i][jM1][kP2]), rho);
   rhoSum += aux[i][jM1][kP1].get_rho_static(&(fields[active][i][jM1][kP1]), rho);
   rhoSum += aux[i][jM1][k].get_rho_static(&(fields[active][i][jM1][k]), rho);
   rhoSum += aux[i][jM1][kM1].get_rho_static(&(fields[active][i][jM1][kM1]), rho);
   rhoSum += aux[i][jM1][kM2].get_rho_static(&(fields[active][i][jM1][kM2]), rho);
   rhoSum += aux[i][jM2][kP1].get_rho_static(&(fields[active][i][jM2][kP1]), rho);
   rhoSum += aux[i][jM2][k].get_rho_static(&(fields[active][i][jM2][k]), rho);
   rhoSum += aux[i][jM2][kM1].get_rho_static(&(fields[active][i][jM2][kM1]), rho);

   rhoSum += aux[iM1][jP2][kP1].get_rho_static(&(fields[active][iM1][jP2][kP1]), rho);
   rhoSum += aux[iM1][jP2][k].get_rho_static(&(fields[active][iM1][jP2][k]), rho);
   rhoSum += aux[iM1][jP2][kM1].get_rho_static(&(fields[active][iM1][jP2][kM1]), rho);
   rhoSum += aux[iM1][jP1][kP2].get_rho_static(&(fields[active][iM1][jP1][kP2]), rho);
   rhoSum += aux[iM1][jP1][kP1].get_rho_static(&(fields[active][iM1][jP1][kP1]), rho);
   rhoSum += aux[iM1][jP1][k].get_rho_static(&(fields[active][iM1][jP1][k]), rho);
   rhoSum += aux[iM1][jP1][kM1].get_rho_static(&(fields[active][iM1][jP1][kM1]), rho);
   rhoSum += aux[iM1][jP1][kM2].get_rho_static(&(fields[active][iM1][jP1][kM2]), rho);
   rhoSum += aux[iM1][j][kP2].get_rho_static(&(fields[active][iM1][j][kP2]), rho);
   rhoSum += aux[iM1][j][kP1].get_rho_static(&(fields[active][iM1][j][kP1]), rho);
   rhoSum += aux[iM1][j][k].get_rho_static(&(fields[active][iM1][j][k]), rho);
   rhoSum += aux[iM1][j][kM1].get_rho_static(&(fields[active][iM1][j][kM1]), rho);
   rhoSum += aux[iM1][j][kM2].get_rho_static(&(fields[active][iM1][j][kM2]), rho);
   rhoSum += aux[iM1][jM1][kP2].get_rho_static(&(fields[active][iM1][jM1][kP2]), rho);
   rhoSum += aux[iM1][jM1][kP1].get_rho_static(&(fields[active][iM1][jM1][kP1]), rho);
   rhoSum += aux[iM1][jM1][k].get_rho_static(&(fields[active][iM1][jM1][k]), rho);
   rhoSum += aux[iM1][jM1][kM1].get_rho_static(&(fields[active][iM1][jM1][kM1]), rho);
   rhoSum += aux[iM1][jM1][kM2].get_rho_static(&(fields[active][iM1][jM1][kM2]), rho);
   rhoSum += aux[iM1][jM2][kP1].get_rho_static(&(fields[active][iM1][jM2][kP1]), rho);
   rhoSum += aux[iM1][jM2][k].get_rho_static(&(fields[active][iM1][jM2][k]), rho);
   rhoSum += aux[iM1][jM2][kM1].get_rho_static(&(fields[active][iM1][jM2][kM1]), rho);

   rhoSum += aux[iM2][jP1][kP1].get_rho_static(&(fields[active][iM2][jP1][kP1]), rho);
   rhoSum += aux[iM2][jP1][k].get_rho_static(&(fields[active][iM2][jP1][k]), rho);
   rhoSum += aux[iM2][jP1][kM1].get_rho_static(&(fields[active][iM2][jP1][kM1]), rho);
   rhoSum += aux[iM2][j][kP1].get_rho_static(&(fields[active][iM2][j][kP1]), rho);
   rhoSum += aux[iM2][j][k].get_rho_static(&(fields[active][iM2][j][k]), rho);
   rhoSum += aux[iM2][j][kM1].get_rho_static(&(fields[active][iM2][j][kM1]), rho);
   rhoSum += aux[iM2][jM1][kP1].get_rho_static(&(fields[active][iM2][jM1][kP1]), rho);
   rhoSum += aux[iM2][jM1][k].get_rho_static(&(fields[active][iM2][jM1][k]), rho);
   rhoSum += aux[iM2][jM1][kM1].get_rho_static(&(fields[active][iM2][jM1][kM1]), rho);

   // Add static rho of central site. IMPORTANT: requires h, G matrices to be ready.

   rhoSum += aux[i][j][k].get_rho_static(&(fields[active][i][j][k]), rho);

   return(rhoSum);
}

void TLattice::lowpass_fields()
{
   int scratch = (active + 1) % 2;

#ifdef _OPENMP
#pragma omp parallel default(shared)
{
#pragma omp for
#endif
   for (int i = 0; i <= nSideM1; i++)
       for (int j = 0; j <= nSideM1; j++)
           for (int k = 0; k <= nSideM1; k++)
           {
               lowpass.compute(fields[active], i, j, k, nSideM1, &(fields[scratch][i][j][k]));
           }

#ifdef _OPENMP
#pragma omp barrier

#pragma omp for
#endif
   for (int i = 0; i <= nSideM1; i++)
       for (int j = 0; j <= nSideM1; j++)
           for (int k = 0; k <= nSideM1; k++)
           {
               fields[active][i][j][k] = fields[scratch][i][j][k];

               fields[active][i][j][k].th_1 = wrap(fields[active][i][j][k].th_1);
               fields[active][i][j][k].th_2 = wrap(fields[active][i][j][k].th_2);
               fields[active][i][j][k].th_3 = wrap(fields[active][i][j][k].th_3);
           }
#ifdef _OPENMP
#pragma omp barrier
}
#endif
}

void TLattice::wrap_th(const int indx)
{
#ifdef _OPENMP
#pragma omp parallel default(shared)
{
#pragma omp for
#endif
   for (int i = 0; i <= nSideM1; i++)
       for (int j = 0; j <= nSideM1; j++)
           for (int k = 0; k <= nSideM1; k++)
           {
              fields[indx][i][j][k].th_1 = wrap(fields[indx][i][j][k].th_1);
              fields[indx][i][j][k].th_2 = wrap(fields[indx][i][j][k].th_2);
              fields[indx][i][j][k].th_3 = wrap(fields[indx][i][j][k].th_3);
           }
#ifdef _OPENMP
#pragma omp barrier
}
#endif
}
