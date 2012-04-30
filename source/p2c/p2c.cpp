/* p2c.cpp : Converts vtk in polar form to cartesian, discards momenta.
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
#include "fields_aux.cpp"

TFields fields[2][MAX_NSIDE][MAX_NSIDE][MAX_NSIDE];   /* Max allowed number of fieldss always allocated. Ugly but fast. */
TAux aux[MAX_NSIDE][MAX_NSIDE][MAX_NSIDE];            /* Ditto */

#include "lattice.cpp"

class TP2CLattice : public TLattice
{
public:
   
   void writeCartesianVTK(string fileName, string title);
};

void TP2CLattice::writeCartesianVTK(string fileName, string title)
{
   int i, j, k;

   ofstream vtk;

   vtk.open(fileName.c_str());
   if (!vtk)
   {
      cout << "ERROR: Could not open VTK file " << fileName << endl << flush;
      return;
   }

   vtk.setf(ios::left, ios::adjustfield);
   vtk.precision(OUTPUT_PRECISION);

   vtk << "# vtk DataFile Version 3.0" << endl ;
   vtk << title << endl;
   vtk << "ASCII" << endl ;
   vtk << "DATASET STRUCTURED_POINTS" << endl ;
   vtk << "DIMENSIONS " << nSideM1+1 << ' ' << nSideM1+1 << ' ' << nSideM1+1 << endl;
   vtk << "ORIGIN 0 0 0" << endl ;
   vtk << "SPACING " << 1 << ' ' << 1 << ' ' << 1 << endl;
   vtk << "POINT_DATA " << (nSideM1+1)*(nSideM1+1)*(nSideM1+1) << endl;

   vtk << "SCALARS PHI_0 DOUBLE" << endl;
   vtk << "LOOKUP_TABLE default" << endl;

   for (i = 0; i <= nSideM1; i++)
       for (j = 0; j <= nSideM1; j++)
           for (k = 0; k <= nSideM1; k++)
           {
              double theta = sqrt(pow(fields[active][i][j][k].th_1, 2) + 
                                  pow(fields[active][i][j][k].th_2, 2) + 
                                  pow(fields[active][i][j][k].th_3, 2)); 

              vtk << cos(theta / 2) << endl;
           }

   vtk << "SCALARS PHI_1 DOUBLE" << endl;
   vtk << "LOOKUP_TABLE default" << endl;

   for (i = 0; i <= nSideM1; i++)
       for (j = 0; j <= nSideM1; j++)
           for (k = 0; k <= nSideM1; k++)
           {
              double thR2 = pow(fields[active][i][j][k].th_1, 2) + 
                            pow(fields[active][i][j][k].th_2, 2) + 
                            pow(fields[active][i][j][k].th_3, 2);

              if (thR2 < SMALL_THETA_SQUARED) // O(th_i^3) series expansion 
              {
                 vtk << fields[active][i][j][k].th_1 * (24 - thR2) / 48 << endl;
              }
              else
              {
                 double theta = sqrt(thR2);
                 vtk << (fields[active][i][j][k].th_1 / theta) * sin(theta / 2) << endl;
              }
           }

   vtk << "SCALARS PHI_2 DOUBLE" << endl;
   vtk << "LOOKUP_TABLE default" << endl;

   for (i = 0; i <= nSideM1; i++)
       for (j = 0; j <= nSideM1; j++)
           for (k = 0; k <= nSideM1; k++)
           {
              double thR2 = pow(fields[active][i][j][k].th_1, 2) + 
                            pow(fields[active][i][j][k].th_2, 2) + 
                            pow(fields[active][i][j][k].th_3, 2);

              if (thR2 < SMALL_THETA_SQUARED) // O(th_i^3) series expansion 
              {
                 vtk << fields[active][i][j][k].th_2 * (24 - thR2) / 48 << endl;
              }
              else
              {
                 double theta = sqrt(thR2);
                 vtk << (fields[active][i][j][k].th_2 / theta) * sin(theta / 2) << endl;
              }
           }

   vtk << "SCALARS PHI_3 DOUBLE" << endl;
   vtk << "LOOKUP_TABLE default" << endl;

   for (i = 0; i <= nSideM1; i++)
       for (j = 0; j <= nSideM1; j++)
           for (k = 0; k <= nSideM1; k++)
           {
              double thR2 = pow(fields[active][i][j][k].th_1, 2) + 
                            pow(fields[active][i][j][k].th_2, 2) + 
                            pow(fields[active][i][j][k].th_3, 2);

              if (thR2 < SMALL_THETA_SQUARED) // O(th_i^3) series expansion 
              {
                 vtk << fields[active][i][j][k].th_3 * (24 - thR2) / 48 << endl;
              }
              else
              {
                 double theta = sqrt(thR2);
                 vtk << (fields[active][i][j][k].th_3 / theta) * sin(theta / 2) << endl;
              }
           }

   vtk << "VECTORS A DOUBLE" << endl;

   for (i = 0; i <= nSideM1; i++)
       for (j = 0; j <= nSideM1; j++)
           for (k = 0; k <= nSideM1; k++)
           {
              vtk << fields[active][i][j][k].A_1 << ' ' << fields[active][i][j][k].A_2 << ' ' << fields[active][i][j][k].A_3 << endl;
           }

   vtk << "SCALARS pPHI_0 DOUBLE" << endl;
   vtk << "LOOKUP_TABLE default" << endl;

   for (i = 0; i <= nSideM1; i++)
       for (j = 0; j <= nSideM1; j++)
           for (k = 0; k <= nSideM1; k++)
           {
              vtk << 0 << endl;
           }

   vtk << "SCALARS pPHI_1 DOUBLE" << endl;
   vtk << "LOOKUP_TABLE default" << endl;

   for (i = 0; i <= nSideM1; i++)
       for (j = 0; j <= nSideM1; j++)
           for (k = 0; k <= nSideM1; k++)
           {
              vtk << 0 << endl;
           }

   vtk << "SCALARS pPHI_2 DOUBLE" << endl;
   vtk << "LOOKUP_TABLE default" << endl;

   for (i = 0; i <= nSideM1; i++)
       for (j = 0; j <= nSideM1; j++)
           for (k = 0; k <= nSideM1; k++)
           {
              vtk << 0 << endl;
           }

   vtk << "SCALARS pPHI_3 DOUBLE" << endl;
   vtk << "LOOKUP_TABLE default" << endl;

   for (i = 0; i <= nSideM1; i++)
       for (j = 0; j <= nSideM1; j++)
           for (k = 0; k <= nSideM1; k++)
           {
              vtk << 0 << endl;
           }

   vtk << "VECTORS pA DOUBLE" << endl;

   for (i = 0; i <= nSideM1; i++)
       for (j = 0; j <= nSideM1; j++)
           for (k = 0; k <= nSideM1; k++)
           {
              vtk << 0 << ' ' << 0 << ' ' << 0 << endl;
           }

   vtk.close();
}

int main(int argc, char* argv[])
{
    cout.setf(ios::left, ios::adjustfield);
    cout.precision(OUTPUT_PRECISION);

    string iFile("optimized");       // Default input VTK filename
    string oFile("cartesian");       // Default output VTK filename
    
    // Parse command line

    int i = 1;
    while (i < argc)
    {
       if (!strcmp("-if", argv[i])) iFile.assign(argv[i + 1]); else
       if (!strcmp("-of", argv[i])) oFile.assign(argv[i + 1]); else
       {
          cout << "Unrecognized command line argument: " << argv[i] << endl;
          cout << endl;
          cout << "Supported options [defaults in square brackets]:" << endl;
          cout << endl;
          cout << "Input filename:  -if [\"" << iFile << "\"]" << endl;
          cout << "Output filename: -of [\"" << oFile << "\"]" << endl;

          cout << flush;

          exit(1);
       }

       i += 2;
    }

    // Initialize lattice

    TP2CLattice lattice;

    // Load initial data

    cout << "Reading input file \"" << iFile << ".vtk\"" << endl << flush;

    string title = "";
    lattice.read_vtk(iFile + ".vtk", title);

    cout << "Done loading." << endl << flush;
    cout << "nSide: " << lattice.get_nSide() << endl;
    
    // Create output file

    cout << "Writing output file \"" << oFile << ".vtk\"" << endl << flush;

    lattice.writeCartesianVTK(oFile + ".vtk", title);

    return(0);
}
