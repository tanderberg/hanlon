/* lattice_1d.cpp : Basic 1D lattice for scalar fields in polar form.
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

typedef TFields TPrivateFields[9];
typedef TAux TPrivateAux[9];

class TLowpass
{
   double c_000, c_001;

public:

   void set_central_weight(double central_weight);

   TLowpass() { set_central_weight(0.0); };

   void compute(TFields fields[MAX_NSIDE], int i, int nSideM1, TFields *result);
   void compute(TPrivateFields fields, TFields *result);
};

void TLowpass::set_central_weight(double central_weight)
{
   c_000 = central_weight;

   c_001 = (1 - central_weight) / 2;
}

void TLowpass::compute(TFields fields[MAX_NSIDE], int i, int nSideM1, TFields *result)
{
   *result = fields[i];
   result->mul_fields(c_000);

   int iM1, iP1;

   if (0 == i) iM1 = nSideM1; else iM1 = i - 1;
   if (nSideM1 == i) iP1 = 0; else iP1 = i + 1;

   TFields partial_sum;

   partial_sum = fields[iM1];
   partial_sum.add_fields(&(fields[iP1]));

   partial_sum.mul_fields(c_001);
   result->add_fields(&partial_sum);
}

void TLowpass::compute(TPrivateFields fields, TFields *result)
{
   *result = fields[4];
   result->mul_fields(c_000);

   TFields partial_sum;

   partial_sum = fields[3];
   partial_sum.add_fields(&(fields[5]));

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

   void dDotFields_dFields(int indx, int momIndx, int i, double* dDotFields);

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
   void compute_staticDotMomenta(int indx, int i, TFields &dotMomenta);

   void copy2private(int indx, int i, TPrivateFields privateFields, TPrivateAux privateAux);
   void copy2private(int i) { copy2private(active, i, privateFields, privateAux); }
   void copy2affected_private_derivatives(int i, int component, TPrivateAux privateAux);
   void copy2affected_private_derivatives(int i, int component)
        { copy2affected_private_derivatives(i, component, privateAux); }
   void compute_affected_private_derivatives(int component, TPrivateFields privateFields, TPrivateAux privateAux);
   void compute_affected_private_derivatives(int component) 
        { compute_affected_private_derivatives(component, privateFields, privateAux); }
   
   void init_aux_static();
   double sum_rho_static();
   double affected_static_rho(int i);
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
      {
         fields[idx][i].clear();
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
   vtk << title << "; t = " << t << endl;
   vtk << "ASCII" << endl ;
   vtk << "DATASET STRUCTURED_POINTS" << endl ;
   vtk << "DIMENSIONS " << nSideM1+1  << endl;
   vtk << "ORIGIN 0" << endl ;
   vtk << "SPACING " << 1 << endl;
   vtk << "POINT_DATA " << nSideM1+1 << endl;

   vtk << "VECTORS T DOUBLE" << endl;

   for (i = 0; i <= nSideM1; i++)
   {
      vtk << fields[active][i].th_1 << ' ' << fields[active][i].th_2 << ' ' << fields[active][i].th_3 << endl;
   }

   vtk << "VECTORS A DOUBLE" << endl;

   for (i = 0; i <= nSideM1; i++)
   {
      vtk << fields[active][i].A_1 << ' ' << fields[active][i].A_2 << ' ' << fields[active][i].A_3 << endl;
   }

   vtk << "VECTORS pT DOUBLE" << endl;

   for (i = 0; i <= nSideM1; i++)
   {
      vtk << 0 << ' ' << 0 << ' ' << 0 << endl;
   }

   vtk << "VECTORS pA DOUBLE" << endl;

   for (i = 0; i <= nSideM1; i++)
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
         }
         else if ("VECTORS T DOUBLE" == line)
         {
            for (int i = 0; i <= nSideM1; i++)
            {
               string line;
               getline(inFile, line);

               int fromTo = 0;

               fields[active][i].th_1 = doubleToken(line, fromTo);
               fields[active][i].th_2 = doubleToken(line, fromTo);
               fields[active][i].th_3 = doubleToken(line, fromTo);
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
         else if ("VECTORS pT DOUBLE" == line)
         {
            for (int i = 0; i <= nSideM1; i++)
            {
               string line;
               getline(inFile, line);
            }
         }
         else if ("VECTORS pA DOUBLE" == line)
         {
            for (int i = 0; i <= nSideM1; i++)
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

   privateAux[6].subtract_offset_drho_dth(-2, &(privateFields[6]), dotMomenta.th_1, dotMomenta.th_2, dotMomenta.th_3);
   privateAux[5].subtract_offset_drho_dth(-1, &(privateFields[5]), dotMomenta.th_1, dotMomenta.th_2, dotMomenta.th_3);
   privateAux[3].subtract_offset_drho_dth(1, &(privateFields[3]), dotMomenta.th_1, dotMomenta.th_2, dotMomenta.th_3);
   privateAux[2].subtract_offset_drho_dth(2, &(privateFields[2]), dotMomenta.th_1, dotMomenta.th_2, dotMomenta.th_3);

   privateAux[4].subtract_drho_static_dth(&(privateFields[4]), dotMomenta.th_1, dotMomenta.th_2, dotMomenta.th_3);

   // dot pA_m =-dH/d A_m

   dotMomenta.A_1 = 0.0;
   dotMomenta.A_2 = 0.0;
   dotMomenta.A_3 = 0.0;

   privateAux[6].subtract_offset_drho_dA(-2, dotMomenta.A_1, dotMomenta.A_2, dotMomenta.A_3);
   privateAux[5].subtract_offset_drho_dA(-1, dotMomenta.A_1, dotMomenta.A_2, dotMomenta.A_3);
   privateAux[3].subtract_offset_drho_dA(1, dotMomenta.A_1, dotMomenta.A_2, dotMomenta.A_3);
   privateAux[2].subtract_offset_drho_dA(2, dotMomenta.A_1, dotMomenta.A_2, dotMomenta.A_3);

   privateAux[4].subtract_drho_static_dA(&(privateFields[4]), dotMomenta.A_1, dotMomenta.A_2, dotMomenta.A_3);
}

void TLattice::compute_staticDotMomenta(int indx, int i, TFields &dotMomenta)
/* Computes time derivatives of  conjugate momenta by differentiating the Hamiltonian w.r.t. fields.
   IMPORTANT: assumes aux data (h, G, M, field derivatives) for cell and affected neighbors are ready for use! */
{
   int iM1, iM2, iP1, iP2;

   if (0 == i) iM1 = nSideM1; else iM1 = i - 1;
   if (0 == iM1) iM2 = nSideM1; else iM2 = iM1 - 1;
   if (nSideM1 == i) iP1 = 0; else iP1 = i + 1;
   if (nSideM1 == iP1) iP2 = 0; else iP2 = iP1 + 1;

   // dot pth_a = - dH/d th_a

   dotMomenta.th_1 = 0.0;
   dotMomenta.th_2 = 0.0;
   dotMomenta.th_3 = 0.0;

   aux[iP2].subtract_offset_drho_dth(-2, &(fields[indx][iP2]), dotMomenta.th_1, dotMomenta.th_2, dotMomenta.th_3);
   aux[iP1].subtract_offset_drho_dth(-1, &(fields[indx][iP1]), dotMomenta.th_1, dotMomenta.th_2, dotMomenta.th_3);
   aux[iM1].subtract_offset_drho_dth(1, &(fields[indx][iM1]), dotMomenta.th_1, dotMomenta.th_2, dotMomenta.th_3);
   aux[iM2].subtract_offset_drho_dth(2, &(fields[indx][iM2]), dotMomenta.th_1, dotMomenta.th_2, dotMomenta.th_3);

   aux[i].subtract_drho_static_dth(&(fields[indx][i]), dotMomenta.th_1, dotMomenta.th_2, dotMomenta.th_3);

   // dot pA_m =-dH/d A_m

   dotMomenta.A_1 = 0.0;
   dotMomenta.A_2 = 0.0;
   dotMomenta.A_3 = 0.0;

   aux[iP2].subtract_offset_drho_dA(-2, dotMomenta.A_1, dotMomenta.A_2, dotMomenta.A_3);
   aux[iP1].subtract_offset_drho_dA(-1, dotMomenta.A_1, dotMomenta.A_2, dotMomenta.A_3);
   aux[iM1].subtract_offset_drho_dA(1, dotMomenta.A_1, dotMomenta.A_2, dotMomenta.A_3);
   aux[iM2].subtract_offset_drho_dA(2, dotMomenta.A_1, dotMomenta.A_2, dotMomenta.A_3);

   aux[i].subtract_drho_static_dA(&(fields[indx][i]), dotMomenta.A_1, dotMomenta.A_2, dotMomenta.A_3);
}

void TLattice::copy2private(int indx, int i, TPrivateFields privateFields, TPrivateAux privateAux)
/* Copies contents of global fields[indx] and aux matrices around i to privateFields and privateAux. 
     Overkill since only the inner [2..6] part of privateAux is needed, and the corners are not needed in either privateFields or privateAux.
*/
{
   int ii[9];
   
   ii[4] = i;

   if (0 == i) ii[3] = nSideM1; else ii[3] = i - 1;
   if (0 == ii[3]) ii[2] = nSideM1; else ii[2] = ii[3] - 1;
   if (0 == ii[2]) ii[1] = nSideM1; else ii[1] = ii[2] - 1;
   if (0 == ii[1]) ii[0] = nSideM1; else ii[0] = ii[1] - 1;

   if (nSideM1 == i) ii[5] = 0; else ii[5] = i + 1;
   if (nSideM1 == ii[5]) ii[6] = 0; else ii[6] = ii[5] + 1;
   if (nSideM1 == ii[6]) ii[7] = 0; else ii[7] = ii[6] + 1;
   if (nSideM1 == ii[7]) ii[8] = 0; else ii[8] = ii[7] + 1;
   
   for (int toi = 0; toi < 9; toi++)
   {
       int fromi = ii[toi];

       privateFields[toi] = fields[indx][fromi];
       privateAux[toi] = aux[fromi];
   }
}

void TLattice::copy2affected_private_derivatives(int i, int component, TPrivateAux privateAux)
/* Copies contents of derivatives of component from global aux matrix around i to privateAux. */
{
   int iM1, iM2, iP1, iP2;

   if (0 == i) iM1 = nSideM1; else iM1 = i - 1;
   if (0 == iM1) iM2 = nSideM1; else iM2 = iM1 - 1;
   if (nSideM1 == i) iP1 = 0; else iP1 = i + 1;
   if (nSideM1 == iP1) iP2 = 0; else iP2 = iP1 + 1;

   privateAux[2].set_derivative(component, aux[iM2].get_derivative(component));
   privateAux[3].set_derivative(component, aux[iM1].get_derivative(component));
   privateAux[5].set_derivative(component, aux[iP1].get_derivative(component));
   privateAux[6].set_derivative(component, aux[iP2].get_derivative(component));
}

void TLattice::compute_affected_private_derivatives(int component, TPrivateFields privateFields, TPrivateAux privateAux)
/* Computes spatial derivatives of privateFields[9] affected by changing the value of field component at 4 and stores them in privateAux [2..6]. 
     Components 0 thru 2 are th_1 thru th_3;
     components 3 thru 5 are A_1 thru A_3. */
{
   privateAux[2].compute_derivative(privateFields, 2, 0, component);
   privateAux[3].compute_derivative(privateFields, 3, 0, component);
   privateAux[5].compute_derivative(privateFields, 5, 0, component);
   privateAux[6].compute_derivative(privateFields, 6, 0, component);
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
   {           
      aux[i].compute_hG(fields[active][i].th_1, fields[active][i].th_2, fields[active][i].th_3);
      aux[i].compute_derivatives(fields[active], i, nSideM1);
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
   {
      TRho rho;
      result += aux[i].get_rho_static(&(fields[active][i]), rho);
   }
#ifdef _OPENMP
#pragma omp barrier
}
#endif

   return(result);
}

double TLattice::affected_private_static_rho(TPrivateFields privateFields, TPrivateAux privateAux)
/* Sums part of total static rho affected by varying privateFields at 4.
     IMPORTANT: Assumes all (center and off-center) h, G matrices and derivatives of off-center cells are ready for use in privateAux.
     Also, does NOT divide by number of fields summed over, so if absolute (rather than relative) magnitude of rho is needed, caller must do that;
     and does NOT divide by dx^2 from derivatives  */
{
   double rhoSum;
   TRho rho;

   // Sum up static rho of off-center fields affected by varying fields at 4.
   // IMPORTANT: requires derivatives and h, G matrices of off-center cells to be ready in privateAux!

   rhoSum  = privateAux[6].get_rho_static(&(privateFields[6]), rho);
   rhoSum += privateAux[5].get_rho_static(&(privateFields[5]), rho);
   rhoSum += privateAux[3].get_rho_static(&(privateFields[3]), rho);
   rhoSum += privateAux[2].get_rho_static(&(privateFields[2]), rho);

   // Add static rho of central site. IMPORTANT: requires h, G matrices to be ready.

   rhoSum += privateAux[4].get_rho_static(&(privateFields[4]), rho);

   return(rhoSum);
}

double TLattice::affected_static_rho(int i)
/* Sums part of total static rho affected by varying fields at i.
     IMPORTANT: Assumes all (center and off-center) h, G matrices and derivatives of off-center cells are ready for use.
     Also, does NOT divide by number of fields summed over, so if absolute (rather than relative) magnitude of rho is needed, caller must do that;
     and does NOT divide by dx^2 from derivatives  */
{
   int iM1, iM2, iP1, iP2;

   if (0 == i) iM1 = nSideM1; else iM1 = i - 1;
   if (0 == iM1) iM2 = nSideM1; else iM2 = iM1 - 1;
   if (nSideM1 == i) iP1 = 0; else iP1 = i + 1;
   if (nSideM1 == iP1) iP2 = 0; else iP2 = iP1 + 1;

   double rhoSum;
   TRho rho;

   // Sum up static rho of off-center fields affected by varying fields at i.
   // IMPORTANT: requires derivatives and h, G matrices of off-center cells to be ready!

   rhoSum  = aux[iP2].get_rho_static(&(fields[active][iP2]), rho);
   rhoSum += aux[iP1].get_rho_static(&(fields[active][iP1]), rho);
   rhoSum += aux[iM1].get_rho_static(&(fields[active][iM1]), rho);
   rhoSum += aux[iM2].get_rho_static(&(fields[active][iM2]), rho);

   // Add static rho of central site. IMPORTANT: requires h, G matrices to be ready.

   rhoSum += aux[i].get_rho_static(&(fields[active][i]), rho);

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
   {
      lowpass.compute(fields[active], i, nSideM1, &(fields[scratch][i]));
   }

#ifdef _OPENMP
#pragma omp barrier

#pragma omp for
#endif
   for (int i = 0; i <= nSideM1; i++)
   {
      fields[active][i] = fields[scratch][i];

      fields[active][i].th_1 = wrap(fields[active][i].th_1);
      fields[active][i].th_2 = wrap(fields[active][i].th_2);
      fields[active][i].th_3 = wrap(fields[active][i].th_3);
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
   {
      fields[indx][i].th_1 = wrap(fields[indx][i].th_1);
      fields[indx][i].th_2 = wrap(fields[indx][i].th_2);
      fields[indx][i].th_3 = wrap(fields[indx][i].th_3);
   }
#ifdef _OPENMP
#pragma omp barrier
}
#endif
}
