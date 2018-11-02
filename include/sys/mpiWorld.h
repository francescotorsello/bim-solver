/**
 *  @file      mpiWorld.h
 *  @brief     The MPI implementation that splits a grid among the ranks.
 *  @author    Mikica Kocic
 *  @copyright GNU General Public License (GPLv3).
 */

#ifndef _MPI_WORLD_H_INCLUDED
#define _MPI_WORLD_H_INCLUDED

#include <mpi.h>

/** Implements the message passing interface (MPI) for high-performance computing (HPC).
 *
 *  In the MPI environment, the calculations on the grid are split among the ranks
 *  so the edges are exchanged between each two neighbors in the following way:
 *  <pre>
 *        .---- left -----.-------.---- right ----.
 *        | ghost | edge  |       | edge  | ghost |
 *        +-------+-----------------------+-------+
 *   ...  |       |     grid chunk i      |       |
 *        +-------+-----------------------+-------+------------+-------+
 *                                |       |   grid chunk i+1   |       |  ...
 *                                +-------+--------------------+-------+
 *                                :     left      :
 *                                : ghost | edge  :
 *                                `-------v-------´
 *                         exchanged between two neighbors              </pre>
 *
 *  The overlapping edges of one rank become other's ghosts (of the same size).
 *  The ghost of the left-most rank is determined from the inner boundary conditions and
 *  the ghost of the right-most rank is determined from the outer boundary conditions.
 *
 *  When compiling, use `mpic++` with `-D_USEMPI`. To run, use `mpiexec` or `mpirun`.
 *
 *  @warning The maximal slicing is not compliant with MPI since the boundary value
 *  problem for the slicing differential equation requires access to the whole grid!
 *
 *  @todo Implement the collective I/O with MPI_File_open and MPI_File_write.
 */
class MPIWorld
{
    int worldSize;           //!< The size of the MPI world
    int myRank;              //!< Our rank in the MPI world
    int leftRank;            //!< The rank of the left neighbor
    int rightRank;           //!< The rank of the right neighbor
    MPI_Request  waitLeft;   //!< Outstanding MPI_Isend request to the left rank
    MPI_Request  waitRight;  //!< Outstanding MPI_Isend request to the right rank
    MPI_Datatype ghostChunk; //!< Describes a ghost-size block at each slice

public:

    /** Initializes the MPI environment and finds our rank in the MPI world.
     */
    MPIWorld ()
    {
        // Initialize the MPI environment
        //
        MPI_Init( NULL, NULL );

        // Find out the rank and the size of the MPI
        //
        MPI_Comm_rank( MPI_COMM_WORLD, &myRank );
        MPI_Comm_size( MPI_COMM_WORLD, &worldSize );

        leftRank  = isFirstInRank() ? -1 : myRank - 1;
        rightRank = isLastInRank()  ? -1 : myRank + 1;

        waitLeft = waitRight = MPI_REQUEST_NULL;
    }

    /** Wait all outstanding MPI requests to finish, then finalize MPI.
     */
    void cleanup ()
    {
        if( waitLeft != MPI_REQUEST_NULL ) {
            MPI_Cancel( &waitLeft );
        }
        if( waitRight != MPI_REQUEST_NULL ) {
            MPI_Cancel( &waitRight );
        }
        MPI_Barrier( MPI_COMM_WORLD );
        MPI_Finalize ();
    }

    ~MPIWorld ()
    {
        cleanup ();
    }

    void quit( int rc )
    {
        cleanup ();
        std::exit( rc );
    }

    int size () const {
        return worldSize;
    }

    int rank () const {
        return myRank;
    }

    bool isFirstInRank () const {
        return myRank == 0;
    }

    bool isLastInRank () const {
        return myRank == worldSize - 1;
    }

    /** Defines a MPI datatype for the ghost chunk, which consists of nGFs blocks
     *  where each block size is nGhost Reals, repeating every nTotal strides.
     *  We have:
     *  `grid index (gf,m,n) = (m) * nGFs * nTotal + (gf) * nTotal + (n)`, hence:
     *
     *  Note that nGFs does not need to be GFCNT (we can transfer a truncated list
     *  of grid functions).
     *
     *  @code
     *  MPI_Type_vector( count, blocklen, stride, oldtype, newtype )
     *
     *     count    = nGFs,         number of blocks
     *     blocklen = nGhost,       number of elements in each block
     *     stride   = nTotal,       number of elements between start of each block
     *     oldtype  = MPI_DOUBLE,   old datatype (handle)
     *     newtype  = &ghostChunk,  new datatype (handle)
     *  @endcode
     *  <pre>
     *   |<- blocklen ->|     |<- blocklen ->|    ...    |<- blocklen ->|
     *   :                    :                                         :
     *   |<----- stride ----->|                                         :
     *   :                                                              :
     *    `------------------------------.-----------------------------´
     *                                 count
     *   Total size = count * blocklen                             </pre>
     */
    void defineGhostChunk( Int nGFs, Int nTotal, Int nGhost )
    {
        MPI_Datatype data_type = sizeof(Real) == sizeof(double)
                               ? MPI_DOUBLE : MPI_LONG_DOUBLE;
        MPI_Type_vector( nGFs, nGhost, nTotal, data_type, &ghostChunk );
        MPI_Type_commit( &ghostChunk );
    }

    /** Exchange the grid regions at boundaries with the other MPI ranks.
     */
    bool exchangeBoundaries
    (
        Real* leftGhost,  //!< Where to put the right edge of the left neighbor
        Real* leftEdge,   //!< Pointer to the left edge of our data
        Real* rightEdge,  //!< Pointer to the right edge of our data
        Real* rightGhost, //!< Where to put the left edge of the right neighbor
        bool receiveEdges = true //!< Whether to receive the edges (or just send ghosts)
        )
    {
        waitLeft = waitRight = MPI_REQUEST_NULL;

        /// - Send our edges (they become other's ghosts)
        ///
        if( leftRank >= 0 ) {
            if( MPI_SUCCESS != MPI_Isend( leftEdge, 1, ghostChunk, leftRank, 0,
                                          MPI_COMM_WORLD, &waitLeft ) ) {
                return false;
            }
        }
        if( rightRank >= 0 ) {
            if( MPI_SUCCESS != MPI_Isend( rightEdge, 1, ghostChunk, rightRank, 0,
                       MPI_COMM_WORLD, &waitRight ) ) {
                return false;
            }
        }

        /// - Optionally receive other's edges so they become our ghosts
        ///
        if( ! receiveEdges ) {
            return true;
        }

        if( leftRank >= 0 )
        {
            MPI_Status status;
            if( MPI_SUCCESS != MPI_Recv( leftGhost, 1, ghostChunk,
                                         leftRank, MPI_ANY_TAG,
                                         MPI_COMM_WORLD, &status )
                || status.MPI_TAG != 0 )
            {
                // std::cerr << "[MPI #" << rank() << " got abort from left]";
                abortExchange( leftGhost );
                return false;
            }
        }
        if( rightRank >= 0 )
        {
            MPI_Status status;
            if( MPI_SUCCESS != MPI_Recv( rightGhost, 1, ghostChunk,
                                         rightRank, MPI_ANY_TAG,
                                         MPI_COMM_WORLD, &status )
                ||status.MPI_TAG != 0 )
            {
                // std::cerr << "[MPI #" << rank() << " got abort from right]";
                abortExchange( leftGhost );
                return false;
            }
        }

        return true;
    }

    /** Wait the other ranks to have received our edges.
     */
    bool waitExchange ()
    {
        if( waitLeft != MPI_REQUEST_NULL ) {
            MPI_Status status;
            if ( MPI_SUCCESS != MPI_Wait( &waitLeft, &status ) ) {
                return false;
            }
        }
        if( waitRight != MPI_REQUEST_NULL ) {
            MPI_Status status;
            if( MPI_SUCCESS != MPI_Wait( &waitRight, &status ) ) {
                return false;
            }
        }
        return true;
    }

    /** Nonblocking signal to all the ranks that we are quitting the integration.
     */
    void abortExchange( Real* dummyData )
    {
        if( leftRank >= 0 ) {
            MPI_Request request;
            MPI_Isend( dummyData, 1, ghostChunk, leftRank, 1, MPI_COMM_WORLD, &request );
            MPI_Request_free( &request );
        }
        if( rightRank >= 0 ) {
            MPI_Request request;
            MPI_Isend( dummyData, 1, ghostChunk, rightRank, 1, MPI_COMM_WORLD, &request );
            MPI_Request_free( &request );
        }
    }
};

#endif // _MPI_WORLD_H_INCLUDED
