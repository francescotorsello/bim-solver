/**
 *  @file      gridInitialData.h
 *  @brief     Reads initial data into the grid.
 *  @author    Mikica Kocic
 *  @copyright GNU General Public License (GPLv3).
 */

#ifndef _GRID_INITIAL_DATA_H_INCLUDED
#define _GRID_INITIAL_DATA_H_INCLUDED

#include <string>
#include <vector>
#include <cstdio>

/** Reads the initial data from the given file.
 */
class GridInitialData : GridUser
{
    bool  binaryFormat;   //!< Read binary data from the input file.
    std::string fileName; //!< Input filename.
    std::vector<Int> input; //! The list of grid functions

public:

    /** Constructs the initial data reader as specified in the parameter file.
     */
    GridInitialData( Parameters& params, UniformGrid& ug )
        : GridUser(ug), binaryFormat(false)
    {
        input.reserve( 256 );
        input.push_back( fld::r );

        params.get( "input.file",   fileName, "stdin" );
        params.get( "input.binary", binaryFormat, true );

        slog << "Initial Data:" << std::endl << std::endl
             << "    file = " << fileName << ",  binary = " << binaryFormat
             << std::endl << std::endl;
    }

    /** The given grid functions will be read as the initial data.
     */
    GridInitialData& addGFs( const std::vector<Int>& gfs )
    {
        input.insert( input.end(), gfs.begin(), gfs.end() );
        return *this;
    }

    /** Loads the initial data.
     */
    bool load ()
    {
        FILE* inf = fileName == "stdin" ? stdin
                  : fopen( fileName.c_str(), binaryFormat ? "rb" : "r" );
        if( ! inf ) {
            std::cerr << "Error: Could not open " << fileName << std::endl;
            return false;
        }

        if ( mpiSize () > 1 ) {
            Int offset = gridDriver->get_nOffset() * input.size() * sizeof(Real);
            fseek( inf, offset, SEEK_SET );
        }

        for( Int n = nGhost; n < nGhost + nLen; ++n )
        {
            // GF( fld::mpiRank, 0, n  ) = mpiRank ();
            if( ! GridPoint( *gridDriver, 0, n ).read( inf, binaryFormat, input ) ) {
                return false;
            }
        }

        if ( inf != stdin ) {
            fclose( inf );
        }

        return true;
    }
};

#endif // _GRID_INITIAL_DATA_H_INCLUDED
