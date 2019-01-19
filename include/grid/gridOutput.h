/**
 *  @file      gridOutput.h
 *  @brief     Writes checkpoints from a grid to an output data file.
 *  @author    Mikica Kocic
 *  @copyright GNU General Public License (GPLv3).
 */

#ifndef _GRID_OUTPUT_H_INCLUDED
#define _GRID_OUTPUT_H_INCLUDED

#include <string>
#include <vector>
#include <deque>
#include <cstdio>

/** Writes checkpoints from the grid to the output data file.
 */
class GridOutputWriter : GridUser
{
    bool  binaryFormat;          //!< Write binary data to output file.
    FILE* outf;                  //!< Output stream for the results.
    std::string fileName;        //!< Output filename
    std::string prefixFileName;  //!< Filename that is copied out (reversed) before start
    std::vector<GF_Descriptor> output; //!< List of grid functions to write out

    /** skipInfo declares a time point at which to set a new `mSkip` and `delta_t` for
     *  the integration.
     */
    struct skipInfo { Real t; Int mSkip; Real delta_t; };

    std::deque<skipInfo> skipSetter;  //!< Apply new "skips" at the given time point

    Int mSkip;     //!< Output every `mSkip` rows
    Int nOut;      //!< Number of grid points to get in the output
    Int nSkip;     //!< Output every `nSkip` points
    Int chunkSize; //!< Size (in bytes) of the output chunk

    Int checkNaNs_nFrom;  //!< Check for NaNs from this `n`
    Int checkNaNs_nTo;    //!< Check for NaNs to this `n`

    /** Writes the results for a particular grid row to output.
     */
    void write( Int m, Int nFrom, Int nTo, Int nSkip = 1 )
    {
        for( Int n = nFrom; n < nTo; n += nSkip ) {
            GridPoint( *gridDriver, m, n ).write( outf, binaryFormat, output );
        }
    }

public:

    Int get_mSkip () const { return mSkip; }
    Int get_nOut  () const { return nOut;  }
    Int get_nSkip () const { return nSkip; }

    Int get_nFrom () const { return checkNaNs_nFrom;  }
    Int get_nTo   () const { return checkNaNs_nTo;    }

    /** Retrieve the size (in bytes) of the output record.
     */
    Int recordSizeInBytes () const { return chunkSize; }

    /** Output size (number of grid functions).
     */
    Int GFCount () const { return output.size(); }

    /** Constructs the output writer as specified in the parameter file.
     */
    GridOutputWriter( Parameters& params, UniformGrid& ug )
        : GridUser( ug ), binaryFormat(false), outf(NULL)
    {
        chunkSize = 0; // Still not defined

        output.reserve( 256 );
        output.push_back( { fld::t, "t", "t" } );
        output.push_back( { fld::r, "r", "r" } );

        params.get( "output.file",   fileName, "stdout" );
        params.get( "output.binary", binaryFormat, true );
        params.get( "output.prefix", prefixFileName, "" );

        slog << "Output Writer:" << std::endl << std::endl
             << "    file = " << fileName << ",  binary = " << binaryFormat
             << std::endl;

        params.get( "output.nOut",  nOut,  10 );
        params.get( "output.nSkip", nSkip,  1 );
        params.get( "output.mSkip", mSkip,  1 );

        slog << "    nOut = " << nOut << ",  nSkip = " << nSkip << ",  mSkip = " << mSkip
             << std::endl;

        params.get( "checkNaNs.nFrom", checkNaNs_nFrom, nOut*9/10 );
        params.get( "checkNaNs.nTo",   checkNaNs_nTo,   nOut      );

        if( checkNaNs_nTo >= 0 ) {
             slog << "    Check for NaNs: nFrom = " << checkNaNs_nFrom
                  << ",  nTo = " << checkNaNs_nTo << std::endl;
        }

        slog << std::endl;

        // Get skip setters
        //
        for( int i = 1; i < 10; ++i )
        {
            std::string tag = "at" + std::to_string( i );
            skipInfo info;
            params.get( tag + ".t",       info.t,       NAN );
            params.get( tag + ".mSkip",   info.mSkip,   1   );
            params.get( tag + ".delta_t", info.delta_t, NAN );
            if ( ! /*std::*/ISNAN( info.t ) ) {
                skipSetter.push_back( info );
            }
        }

        if( skipSetter.size() > 0 )
        {
            slog << "Skip Setter:" << std::endl << std::endl;
            for( auto i: skipSetter ) {
                slog << "    t = " << i.t << ", mSkip = " << i.mSkip
                    << ", delta_t = " << i.delta_t << std::endl;
            }
            slog << std::endl;
        }
    }

    ~GridOutputWriter ()
    {
        close ();
    }

    /** The given grid functions will be written to the output.
     */
    void gridFunctions( const std::vector<GF_Descriptor>& gfs )
    {
        output.insert( output.end(), gfs.begin(), gfs.end() );
        // slog << "    Added " << gfs.size () << " GFs to the output" << std::endl;
    }

    /** Opens output data file for storing the results.
     */
    bool open ()
    {
        // Each record has two ghost regions in addition to the middle nOut/nSkip part
        //
        chunkSize = ( 2 * nGhost + nOut / nSkip ) * output.size() * sizeof(Real);

        if( fileName == "stdout" ) {
            binaryFormat = false;
            outf = stdout;
            return false;
        }

        // When using MPI, insert our rank number in the filename
        //
        if ( mpiSize() > 1 )
        {
            size_t lastdot = fileName.find_last_of( "." );
            fileName = fileName.substr( 0, lastdot ) + "-"
                       + std::to_string( 1000 + mpiRank() ).substr( 1 )
                       + fileName.substr( lastdot );
        }

        outf = mpiRank() >= 4 ? NULL
                : fopen( fileName.c_str(), binaryFormat ? "wb+" : "w+" );
/*
        if( binaryFormat && prefixFileName != "" )
        {
            std::ifstream pref( prefixFileName, std::ios::in | std::ios::binary );
            if( pref )
            {
                slog << "*** Appending the reversed-t: " << prefixFileName << std::endl;

                pref.seekg( 0, std::ios::end );  // Go to end of file
                auto len = pref.tellg ();   // Get the length of file (in bytes)

                if ( ( len % chunkSize ) != 0 ) {
                    std::cerr << "Error: Invalid grid alignment." << std::endl;
                    return false;
                }

                // While len > 0 means skip the last record (at t_0)
                for( len -= chunkSize; len > 0; len -= chunkSize )
                {
                    pref.seekg( len, std::ios::beg );
                    pref.read( (char*)&grid[0][0], chunkSize );
                    fwrite( &grid[0][0], chunkSize, 1, outf );
                }
                pref.close ();
            }
        }
*/
        return true;
    }

    /** Closes the output file.
     */
    void close ()
    {
        if( outf != NULL ) {
            fclose( outf );
            outf = NULL;
        }
    }

    /** Writes the header to output file with the information about grid functions.
     */
    void writeHeader ()
    {
        if( outf != NULL )
        {
            for( auto f: output ) {
                fprintf( outf, "%s -> \"%s\"\n", f.name.c_str(), f.tex.c_str() );
            }
            fprintf( outf, "*** DATA ***\n");
        }
    }

    void write( Int m, Real cur_t, Real& delta_t )
    {
        if ( outf != NULL )
        {
            write( m, 0, nGhost );
            write( m, nGhost, nGhost + std::min( nLen, nOut ), nSkip );
            write( m, nGhost + nLen, nTotal );
        }

        if( skipSetter.size() > 0 && cur_t >= skipSetter.front().t )
        {
            mSkip = skipSetter.front().mSkip;
            Real new_dt = skipSetter.front().delta_t;
            if( ! /*std::*/ISNAN( new_dt ) )
            {
                delta_t = new_dt;
                slog << std::endl << std::endl << "*** at t = " << cur_t
                     << " new mSkip = " << mSkip << ", delta_t = " << delta_t
                     << std::endl << std::endl;
            }
            else
            {
                slog << std::endl << std::endl << "*** at t = " << cur_t
                     << " new mSkip = " << mSkip << std::endl << std::endl;
            }
            skipSetter.pop_front ();
        }
    }
};

#endif // _GRID_OUTPUT_H_INCLUDED
