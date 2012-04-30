/* string_aux.cpp : Auxiliary string functions.
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

string trim(string line)
{
   int i = line.find_first_not_of(" \n\r", 0);

   if (0 > i) return("");

   int j = line.find_first_of("\n\r", i);

   if (0 > j) j = line.length();

   if (j <= i) return("");

   string s = line.substr(i, j - i);

   return(s);
}

string token(string line, int& fromTo, string separators = " \n\r")
{
   int i = line.find_first_not_of(separators, fromTo);

   if (0 > i) return("");

   fromTo = line.find_first_of(separators, i);

   if (0 > fromTo) fromTo = line.length();

   if (fromTo <= i) return("");

   string s = line.substr(i, fromTo - i);

   return(s);
}

double doubleToken(string line, int& fromTo, string separators = " \n\r")
{
   string s = token(line, fromTo, separators);
   if ("" == s) return(0);
   return(atof(s.c_str()));
}

string stripExtension(string s)
{
   int i = s.find_last_of(".");
   if (0 > i) return(s);
   return(s.substr(0, i));
}

