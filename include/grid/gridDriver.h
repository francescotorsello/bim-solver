/**
 *  @file      gridDriver.h
 *  @brief     Implements a spatially cell-centered uniform grid.
 *  @author    Mikica Kocic
 *  @copyright GNU General Public License (GPLv3).
 *
 *  @page fdm Finite Difference Approximation
 *
 *  At the core of finite difference approximation is a discretization of the
 *  spacetime, or a numerical grid (shown below).\ Then, a function `f(t,r)` is
 *  represented by values at a discrete set of points `(m,n)`.
 *
 *  <pre>
 *                p in M          <-->     (t,r)      <-->        (m,n)
 *         (point in a manifold)        (in a chart)        (discretized point)
 *                                           |                      |
 *                                           V                      V
 *                                         f(t,r)                 f(m,n)    </pre>
 *
 *  The points `(m,n)` are called the grid points or nodes. The grid points are chosen
 *  to reside at the center of grid cells (the other approach would be to associate the
 *  grid points with the vertices of the grid cells).&nbsp;
 *
 *  The data structure @ref GridPoint is a wrapper for the set of fields `{f,...}` that
 *  arekept in the memory for the point `(m,n)`.
 *
 *  The grid consists of the *interior points*, and the *virtual points* at the boundary.
 *  The virtual points are used to impose the boundary conditions.&nbsp;
 *
 *  A cell-centered uniform grid with `nLen` grid points in each grid row is shown below
 *  (see also @cite Baumgarte:2010nr, p.192).
 *
 *  <pre>
 *         left virtual                grid cell        grid point     right virtual
 *            point(s)                 .---^---.            |             point(s)
 *      .        |               n-1  '    n    '  n+1      |                |
 *      :   +----|----+-------+-------+---------+-------+---|---+-------+----|----+
 *     m+1  |    O    |       |       | (m+1,n) |       |   X   |       |    O    |
 *          +---------+-------+-------.=========.-------+-------+-------+---------+
 *      m   |         |       |(m,n-1)|  (m,n)  |(m,n+1)|       |       |         |
 *          +---------+-------+-------'========='-------+-------+-------+---------+
 *     m-1  |         |       |       |         |       |       |       |         |
 *      .   +---------'-------+-------+---------+-------+-------+-------'---------+
 *      :        0    |   1       2      . . .                    nLen  |  nLen+1
 *                    ^                                                 ^
 *                   r_min                                            r_max       </pre>
 *
 *  The cell-centered uniform grid implies:
 *   @code
 *     r = r_min + ( n - 1/2 ) delta_r
 *     delta_r = ( r_max - r_min ) / nLen
 *   @endcode
 *   @li if `nGhost = 1`, `n` ranges between 1 and `nLen`;
 *       the boundary virtual points are `grid[m][0]` and `grid[m][nLen+1]`
 *
 *   @li `m` ranges between 1 and mLen, and wraps around;
 *       the slices `grid[mLen+1]`... are used for intermediate integration steps.
 *
 *  In parallel simulations, the grid is split across processes using MPI,
 *  allowing storage of large grid functions and parallel work on each part.
 *
 */

#ifndef _GRID_DRIVER_H_INCLUDED
#define _GRID_DRIVER_H_INCLUDED

/////////////////////////////////////////////////////////////////////////////////////////
/** @defgroup g4 Grid driver                                                           */
/////////////////////////////////////////////////////////////////////////////////////////
                                                                                   /*@{*/
/** Uniform grid of grid functions.
 */
class UniformGrid
    : public MPIWorld
{
    /** Allocates memory for the numerical grid.
     *
     *  We allocate `mExtra` extra m grid rows for the intermediate steps.
     *  We also allocate 2*nGhost n-cells as the L/R spatial boundary.
     */
    void createGrid ()
    {
        size_t mDim = mLen + mExtra;
        grid = new Real[ mDim * nTotal * GFCNT ];
        std::fill_n( grid, mDim * nTotal * GFCNT, NAN );

        slog << "    Allocated " << mDim << " * " << nTotal
             << " * " << GFCNT
             << " * Real" << 8 * sizeof(Real)
             << " = " << mDim * nTotal * GFCNT
             << " bytes" << std::endl
             << "    Tracking " << GFCNT << " grid functions at a point"
             << std::endl << std::endl;
    }

    /** Releases memory allocated for the numerical grid.
     */
    void deleteGrid ()
    {
        if ( grid != NULL ) {
            delete[] grid;
            grid = NULL;
        }
    }

protected:

    Real* grid;      //!< The grid storage
                     //   Spatial discretization:
    Int   nOffset;   //!< Offset in the MPI world
    Int   nWorldLen; //!< Number of cells in the MPI world
    Int   nLen;      //!< Number of our cells in `r`-direction
    Int   nGhost;    //!< Number of virtual (ghost) cells on the `r`-boundary
    Int   nTotal;    //!< Total number of `r`-cells (`nTotal = nLen + 2 * nGhost`)
    Real  delta_r;   //!< The grid spacing across `r`
                     //   Temporal discretization:
    Int   mLen;      //!< Number of cached cells in `t`-direction
    Int   mExtra;    //!< Number of extra cells (for integration substeps)

public:

    Real* get_grid () { return grid; }

    Int  get_nOffset   () const { return nOffset;   }
    Int  get_nWorldLen () const { return nWorldLen; }
    Int  get_nLen      () const { return nLen;      }
    Int  get_nGhost    () const { return nGhost;    }
    Int  get_nTotal    () const { return nTotal;    }
    Int  get_mLen      () const { return mLen;      }
    Real get_delta_r   () const { return delta_r;   }

    Int  mpiSize       () const { return size ();   }
    Int  mpiRank       () const { return rank ();   }

    /** Constructs the uniform grid as specified in the parameter file.
     */
    UniformGrid( Parameters& params )
    {
        if ( ! isFirstInRank () ) {
            slog.disabled = true; // only the first in rank can access the output
        }

        params.get( "grid.nLen",    nWorldLen, 10  );
        params.get( "grid.nGhost",  nGhost,    10  );
        params.get( "grid.delta_r", delta_r,   0.1 );
        params.get( "grid.mLen",    mLen,      5   );
        params.get( "grid.mExtra",  mExtra,    9   );

        /// - Split the data equally among all the ranks in the MPI world
        ///
        nLen = nWorldLen / mpiSize();
        nOffset = mpiRank() * nLen; // Seek offset in the MPI world

        /// - The last in rank gets 'all the rest'
        ///
        if ( mpiRank() == mpiSize() - 1 && mpiSize() > 1 ) {
            nLen += ( nWorldLen % mpiSize() );
        }

        nTotal = nLen + 2 * nGhost;

        slog << "Uniform Grid Driver:" << std::endl << std::endl
             << "    nLen = " << nLen << ",  nGhost = " << nGhost
             << ",  delta_r = " << delta_r
             << ",  mLen = " << mLen << ",  mExtra = " << mExtra << std::endl
             << "    Compile-time: CFDS_ORDER = " << CFDS_ORDER << std::endl;

        #if _OPENMP
            if ( mpiSize() > 1 ) {
                omp_set_num_threads( 1 ); /// @todo make omp.nthreads configurable
            }
            slog << "    Using OpenMP: #procs = "
                 << omp_get_num_procs() << ", max #threads = "
                 << omp_get_max_threads() << std::endl;
        #endif

        if ( mpiSize() > 1 ) {
             slog << "    Using MPI: world size = " << mpiSize()
                  << ",  sum of world's nLen " << nWorldLen << ",  nOffset = " << nOffset
                  << std::endl;
        }

        /// - Sanity check that nGhost > CFDS_ORDER/2
        ///
        if( nGhost < CFDS_ORDER/2 ) {
            slog << "Error: Ghost region does not fit the FD order." << std::endl;
            quit( -1 );
        }

        /// - Allocate memory for the grid
        ///
        createGrid ();

        /// - Configure MPI ghost chunk datatype
        ///
        /// @warning  Using fld::mpiBoundary to split GFs works only if the grid index
        ///           is organized as `(<m> * GFCNT + <gf> ) * nTotal + <n>`.
        ///           Replace fld::mpiBoundary by GFCNT to exchange all the GFs.
        ///
        defineGhostChunk( fld::mpiBoundary + 1, nTotal, nGhost );
    }

    /** Access the given grid function data.
     */
    inline Real& GF( Int gf, Int m, Int n ) {
        return  grid[ (m * GFCNT + gf ) * nTotal + n ];
    }

    /** Exchange the boundaries at the given time slice.
     */
    bool exchangeBoundaries( Int m, bool receiveEdges = true )
    {
        Real* timeSlice = &GF( fld::t, m, 0 );
        return MPIWorld::exchangeBoundaries(
            timeSlice,                           // Left ghost  (received)
            timeSlice + nGhost,                  // Left edge   (sent)
            timeSlice + nGhost + nLen - nGhost,  // Right_edge  (sent)
            timeSlice + nGhost + nLen,           // Right ghost (received)
            receiveEdges  // True if to receive the edges (otherwise just send ghosts)
        );
    }

    /** Nonblocking signal to all the ranks that we are quitting integration.
     */
    void abortExchange ()
    {
        MPIWorld::abortExchange( /* dummy data */ &GF( fld::t, 0, 0 ) );
    }
};

/////////////////////////////////////////////////////////////////////////////////////////

/** GridUser contains cached variables from the grid-driver.
 *  @note Finite differences are defined for the grid users.
 */
class GridUser
{
private:
    Real* grid;  //!< Private access to the grid storage (not shared with descendants)

protected:
    UniformGrid* gridDriver;  //!< The attached driver (available to descendants)

    Int nOffset;        //!< Offset in the MPI world
    Int nWorldLen;      //!< Number of cells in the MPI world
    Int nLen;           //!< Number of cells in `r`-direction
    Int nGhost;         //!< Number of virtual (ghost) cells on the `r`-boundary
    Int nTotal;         //!< Total number of cells in `r`-direction including ghosts
    Int mLen;           //!< Number of cells in `t`-direction

    Real delta_r;       //!< The grid spacing accross `r`
    Real inv_delta_r;   //!< `1 / delta_r`
    Real inv_delta_rr;  //!< `1 / delta_r^2`

public:
    /////////////////////////////////////////////////////////////////////////////////////

    Int mpiSize () const { return gridDriver->mpiSize (); }
    Int mpiRank () const { return gridDriver->mpiRank (); }

    GridUser( UniformGrid& ug )
        : gridDriver( &ug )
    {
        grid      = ug.get_grid ();
        nOffset   = ug.get_nOffset ();
        nWorldLen = ug.get_nWorldLen ();
        nLen      = ug.get_nLen ();
        mLen      = ug.get_mLen ();
        nGhost    = ug.get_nGhost ();
        nTotal    = ug.get_nTotal ();
        delta_r   = ug.get_delta_r ();

        // Calculate various constants involving 1/delta_r
        //
        inv_delta_r  = 1 / ug.get_delta_r ();
        inv_delta_rr = 1 / ( ug.get_delta_r () * ug.get_delta_r () );
    }

    /** Access the given grid function data.
     */
    inline Real& GF( Int gf, Int m, Int n ) {
        return  grid[ (m * GFCNT + gf ) * nTotal + n ];
    }

    /** `emitField` is a macro that emits the method which encapsulates a grid function.
     *  Such method can be later used to access data in a grid at a point given by (m,n).
     *  For example, `emitField(r)` defines `r(m,n)`.
     */
    #define emitField(gf) \
        inline Real& gf( Int m, Int n ) { return GF( fld::gf, m, n ); }

#if 1
    #define emitDerivative_r(gf)   emitField(gf##_r)
    #define emitDerivative_rr(gf)  emitField(gf##_rr)
#else
    #define emitDerivative_r(gf) \
        inline Real gf##_r( Int m, Int n ) { return GF_r( gf, m, n ); }
    #define emitDerivative_rr(gf) \
        inline Real gf##_rr( Int m, Int n ) { return GF_rr( gf, m, n ); }
#endif

    // The system GFs are 'time' and 'space', which are defined in "gridFunctions.h"
    //
    emitField( t )
    emitField( r )

    /////////////////////////////////////////////////////////////////////////////////////

    /** Cubic spline smoother of a grid function.
     */
    void cubicSplineSmooth( Int m, Int gf, Int lin2n, Int cub2n )
    {
        // The accuracy is very low at low r. Assume that a few first cells are linear.
        //
        if( lin2n > 0 )
        {
            Real dydx = delta_r * GF( gf, m, nGhost + lin2n ) / r( m, nGhost + lin2n );
            for ( Int i = 0; i < lin2n; ++i ) {
                GF( gf , m, nGhost + i ) = ( i + 0.5 ) * dydx;
            }
        }

        // A six-point cubic spline to smooth the region around r = 0.
        //
        if( cub2n > 0 )
        {
            static CubicSpline spline( 6 );

            Real pts[6][2] = {
                { r( m, nGhost           ),   GF( gf, m, nGhost           ) },
                { r( m, nGhost + 3       ),   GF( gf, m, nGhost + 3       ) },
                { r( m, nGhost + cub2n/2 ),   GF( gf, m, nGhost + cub2n/2 ) },
                { r( m, nGhost + cub2n-6 ),   GF( gf, m, nGhost + cub2n-6 ) },
                { r( m, nGhost + cub2n-3 ),   GF( gf, m, nGhost + cub2n-3 ) },
                { r( m, nGhost + cub2n   ),   GF( gf, m, nGhost + cub2n   ) }
            };
            spline.initialize( pts );

            for( Int i = 0; i < cub2n; ++i ) {
                GF( gf, m, nGhost + i ) = spline( r( m, nGhost + i ) );
            }
        }
    }

    /////////////////////////////////////////////////////////////////////////////////////

    void applyBoundaryConditions( Int m, Int gf, Int parity );
    void smoothenGF( Int m, Int outgf, Int tmpgf, Int ingf, Int parity );
    void smoothenGF2( Int m, Int outgf, Int tmpgf, Int ingf, Int parity );
    void smoothenGF0( Int m, Int nCopyFrom, Int nCopyTo, Int sgRadius, Int outgf, Int tmpgf, Int ingf, Int parity );
};

/////////////////////////////////////////////////////////////////////////////////////////

/** A wrapper for a grid point so that the all of the grid functions
 *  can be treated as a compound object when reading/writing.
 */
class GridPoint : GridUser
{
    Int m, n;
public:

    GridPoint( UniformGrid& ug, Int mPos, Int nPos )
        : GridUser( ug ), m( mPos ), n( nPos )
    {}

    /** Resets all the variables to a specific value (by default NAN).
     */
    void clear( Real value = NAN )
    {
        for( Int gf = 0; gf < GFCNT; ++gf ) {
            GF( gf, m, n ) = value;
        }
    }

    /** Reads the variable values from a file stream.
     */
    bool read( FILE* inf, bool isBinary, const std::vector<Int>& input )
    {
        auto len = input.size();
        Real* data = (Real*)alloca( len * sizeof(Real) );

        if( isBinary )
        {
            if( len != fread( data, sizeof(data[0]), len, inf ) ) {
                std::cerr << "Error: Premature end of grid." << std::endl;
                return false;
            }
        }
        else
        {
            char* line = (char*)alloca( input.size() * 80 * sizeof(char) );
            line[0] = '\0';
            do {
                if( ! fgets( line, input.size() * 80 - 1, inf ) ) {
                    std::cerr << "Error: Premature end of grid." << std::endl;
                    return false;
                }
            } while( line[0] == '*' ); // Ignore comments (starting with '*')

            size_t count = sscanf_Real( line, data, len );
            if( count != len ) {
                std::cerr << "Error: Wrong number of variables." << std::endl;
                return false;
            }
        }

        for( size_t gf = 0; gf < len; ++gf ) {
            GF( input[gf], m, n ) = data[gf];
        }

        return true;
    }

    /** Writes the variable values to a file stream.
     *  @warning The order of fields in data[] (output data record)
     *           should match the Mathematica notebook!
     */
    void write( FILE* outf, bool isBinary, const std::vector<GF_Descriptor> output )
    {
        auto len = output.size(); /// @todo fld::output should be dynamic in ID class
        Real* data = (Real*)alloca( len * sizeof(Real) );

        for( size_t gf = 0; gf < len; ++gf ) {
            data[gf] = GF( output[gf].gf, m, n );
        }

        if( isBinary ) {
            fwrite( data, sizeof(data[0]), len, outf );
        }
        else {
            for( size_t i = 0; i < len; ++i ) {
                if ( i > 0 ) fputs( "\t", outf );
                fputReal( outf,  data[i] );
            }
            fputs( "\n", outf );
        }
    }
};

/////////////////////////////////////////////////////////////////////////////////////////

void GridUser::applyBoundaryConditions( Int m, Int gf, Int parity )
{
    // Left ghost region
    for( Int i = 0; i < nGhost; ++i ) {
        GF( gf, m, nGhost - i - 1 ) = parity * GF( gf, m, nGhost + i );
    }

    // Right ghost region
    for( Int n = nGhost + nLen; n < nTotal; ++n ) {
        extrapolate_R( gf, m, n );
    }
}

void GridUser::smoothenGF2( Int m, Int copy2gf, Int outgf, Int ingf, Int parity )
{
    smoothenGF( m, -1, fld::gdbg, ingf, parity );
    smoothenGF( m, copy2gf, outgf, fld::gdbg, parity );
}

/** Smooth data using a local polynomial regression. The convolution is done with
 *  respect to the Savitzky-Golay matrix that corresponds to a smoothing kernel of
 *  radius `sgRadius` and a polynomial regression of degree `polyDeg`.
 *  @pre `outf` must not be equal to `ingf`
 *  @post `outf` is copied to `copy2gf` if `copy2gf >= 0`
 *  @warning Does not work with MPI!
 *  @see Mathematica notebook `Savitzky-Golay smoothing.nb`, which tests the algorithm.
 */
void GridUser::smoothenGF( Int m, Int copy2gf, Int outgf, Int ingf, Int parity )
{
    // Smoothing parameters:
    Int sgRadius = 32;  // Default kernel radius
    Int order    = 2;   // Polynomial order is twice of this
    Int guard    = 3;   // guard buffer between NaN and the convolution window

    // Find the first not-NaN from the left
    //
    Int nFrom = nGhost;
    while( nFrom < nTotal && std::isnan( GF( ingf, m, nFrom ) ) ) {
        ++nFrom;
    }

    // If there were no NaN's, move nFrom as far as possible to the left
    //
    if( nFrom == nGhost ) {
        nFrom = -nLen;
    }

    // Find the last not-NaN from the right (do not average with the extrapolated region)
    //
    Int nTo = nGhost + nLen;
    while( nTo >= nGhost && std::isnan( GF( ingf, m, nTo - 1 ) ) ) {
        --nTo;
    }

    OMP_parallel_for( Int n = nFrom < 0 ? 0 : nFrom ; n < nTo; ++n )
    {
        // Find the minimum radius (between 0 and sgRadius)
        //
        Int left = n - nFrom;
        Int right = nTo - guard - 1 - n;
        Int rad = right < left ? right : left;
        if ( rad < 0 ) {
            rad = 0;
        }
        else if( rad > sgRadius ) {
            rad = sgRadius;
        }

        // Convolve with the coefficients.
        // Parity pair on the staggered grid where r(m,nGhost-1) = -r(m,nGhost) reads:
        //
        //    gf( m, nGhost + k ) == gf( m, nGhost-1 - k ) * parity,  where k >= 0
        //
        // Let k = n + i - nGhost; then:
        //
        //    if( k >= 0 )
        //    then: use gf( m, nGhost + k )
        //    else: use gf( m, nGhost-1 - k ) * parity
        //          which is the same as gf( m, nGhost-1 - (n + i - nGhost) ) * parity
        //          or gf( m, 2*nGhost - n - i - 1 ) * parity
        //
        Real sum = 0;
        for( Int i = -rad; i <= rad; ++i ) {
            sum += SG_coeff[rad][order][ rad + i ]
                   * ( n + i >= nGhost ? GF( ingf, m, n + i )
                                       : parity * GF( ingf, m, 2 * nGhost - n - i - 1 ) );
        }
        GF( outgf, m, n ) = sum;
    }

    // Extrapolate to the right
    //
    for( Int n = nTo; n < nTotal; ++n )
    {
        extrapolate_R( outgf, m, n );  /// @todo should be a choice of extrapolation
        //extrapolate_D( outgf, m, n );
    }

    // Optionally copy the output to a new location
    //
    if( copy2gf >= 0 )
    {
        OMP_parallel_for( Int n = nFrom; n < nTotal; ++n ) {
            GF( copy2gf, m, n ) = GF( outgf, m, n );
        }
    }
}

void GridUser::smoothenGF0( Int m, Int nCopyFrom, Int nCopyTo, Int sgRadius, Int copy2gf, Int outgf, Int ingf, Int parity )
{
    // Smoothing parameters:
    //Int sgRadius = 32;  // Default kernel radius
    Int order    = 2;   // Polynomial order is twice of this
    Int guard    = 3;   // guard buffer between NaN and the convolution window

    // Find the first not-NaN from the left
    //
    Int nFrom = nGhost;
    while( nFrom < nTotal && std::isnan( GF( ingf, m, nFrom ) ) ) {
        ++nFrom;
    }

    // If there were no NaN's, move nFrom as far as possible to the left
    //
    if( nFrom == nGhost ) {
        nFrom = -nLen;
    }

    // Smoothen only within a region close to r=0
    //
    Int nTo = nGhost + nLen;
    while( nTo >= nGhost && std::isnan( GF( ingf, m, nTo - 1 ) ) ) {
        --nTo;
    }

    OMP_parallel_for( Int n = nFrom < 0 ? 0 : nFrom ; n < nTo; ++n )
    {
        // Find the minimum radius (between 0 and sgRadius)
        //
        Int left = n - nFrom;
        Int right = nTo - guard - 1 - n;
        Int rad = right < left ? right : left;
        if ( rad < 0 ) {
            rad = 0;
        }
        else if( rad > sgRadius ) {
            rad = sgRadius;
        }

        // Convolve with the coefficients.
        // Parity pair on the staggered grid where r(m,nGhost-1) = -r(m,nGhost) reads:
        //
        //    gf( m, nGhost + k ) == gf( m, nGhost-1 - k ) * parity,  where k >= 0
        //
        // Let k = n + i - nGhost; then:
        //
        //    if( k >= 0 )
        //    then: use gf( m, nGhost + k )
        //    else: use gf( m, nGhost-1 - k ) * parity
        //          which is the same as gf( m, nGhost-1 - (n + i - nGhost) ) * parity
        //          or gf( m, 2*nGhost - n - i - 1 ) * parity
        //
        Real sum = 0;
        for( Int i = -rad; i <= rad; ++i ) {
            sum += SG_coeff[rad][order][ rad + i ]
                   * ( n + i >= nGhost ? GF( ingf, m, n + i )
                                       : parity * GF( ingf, m, 2 * nGhost - n - i - 1 ) );
        }
        GF( outgf, m, n ) = sum;
    }

    // Optionally copy the output to a new location, up to the grid point nCopyTo (this is to avoid issues with the right boundary conditions)
    //
    if( copy2gf >= 0 )
    {
        OMP_parallel_for( Int n = nCopyFrom/*nFrom < 0 ? 0 : nFrom*/; n < nCopyTo; ++n ) {
            GF( copy2gf, m, n ) = GF( outgf, m, n );
        }
    }
}
                                                                                   /*@}*/
/////////////////////////////////////////////////////////////////////////////////////////

//#include "numMethods/SavitzkyGolayFilter.h"

#endif // _GRID_DRIVER_H_INCLUDED
