/**
 *  @file      bispec.cpp
 *  @brief     The covariant BSSN evolution for spherically symmetric bimetric spacetimes, with the pseudospectral method.
 *  @authors   Francesco Torsello
 *  @copyright GNU General Public License (GPLv3).
 */

#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <cmath>
#include <cstdio>
#include <chrono>
#include <vector>
#include <bitset>

#ifndef OBSERVER
    #define OBSERVER 1
#endif // OBSERVER

#ifndef _EVOLVE_DSIG
    #define _EVOLVE_DSIG 0
#endif // _EVOLVE_DSIG

#ifndef _DETECT_NAN
    #define _DETECT_NAN 1
#endif // _DETECT_NAN

using namespace std;

void dumpbytes(const vector<char>& v)
{
    for (int i=0; i<v.size(); ++i)
    {
        printf("%u ", (unsigned char)v[i]);
        if ((i+1) % 16 == 0)
            printf("\n");
    }
    printf("\n");
}

long fromBin(long n)
{
    long factor = 1;
    long total = 0;

    while (n != 0)
    {
        total += (n%10) * factor;
        n /= 10;
        factor *= 2;
    }

    return total;
}

/////////////////////////////////////////////////////////////////////////////////////////
/// The main entry point of `bispec`.
///
int main( int argc, char* argv[] )
{
    string line;
    ifstream chebyshev;
    chebyshev.open( "include/chebyshev-values/test.dat", ios::in | ios::binary );

    if( chebyshev.is_open() )
    {
        // get the starting position
      streampos startpos = chebyshev.tellg();

      // go to the end
      chebyshev.seekg(0, std::ios::end);

      // get the ending position
      streampos endpos = chebyshev.tellg();

      // go back to the start
      chebyshev.seekg(0, std::ios::beg);

        // create a vector to hold the data that
      // is resized to the total size of the file
        std::vector<char> cheby_values;
        cheby_values.resize(static_cast<size_t>(endpos - startpos));

          // create a vector to hold the data that
      // is resized to the total size of the file
        std::bitset<endpos - startpos> cheby_values_bin;
        //cheby_values_bin.resize(static_cast<size_t>(endpos - startpos));

        // read it in
      chebyshev.read(&cheby_values[0], cheby_values.size());

        // print it out (for clarity)
      //for(const char& c : cheby_values)
      //{
      //   cout << (int) c << endl;
      //}

      //dumpbytes( cheby_values );

        std::vector<char> cheby_values_dec;
        cheby_values_dec.resize(cheby_values.size());

        int i;

        for( i = 0; i <= cheby_values.size(); i++ )
        {
            cheby_values_dec[i] = fromBin( cheby_values[i] );
        }

      cout.write(cheby_values_dec.data(), cheby_values_dec.size());

    } else {

        cout << "Unable to open chebyshev file.";

    }

    return 0;

}
