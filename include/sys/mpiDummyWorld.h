/**
 *  @file      mpiDummyWorld.h
 *  @brief     A dummy message passing interface (MPI).
 *  @author    Mikica Kocic
 *  @copyright GNU General Public License (GPLv3).
 */

#ifndef _MPI_DUMMY_WORLD_H_INCLUDED
#define _MPI_DUMMY_WORLD_H_INCLUDED

/** Implements a dummy message passing interface (MPI).
 *  To enable the real MPI code, use `mpic++` with `-D_USEMPI`.
 *
 *  @warning The maximal slicing is not compliant with MPI since the boundary value
 *  problem for the slicing differential equation requires access to the whole grid!
 */
class DummyMPIWorld
{
public:
    void quit( int rc )
    {
        exit( rc );
    }

    int size () const { return 1; } //!< The size of the MPI world
    int rank () const { return 0; } //!< Our rank in the MPI world

    bool isFirstInRank () const { return true; }
    bool isLastInRank () const { return true; }

    /** Defines a MPI datatype for the ghost chunk, which consists of nGFs blocks
     *  where each block size is nGhost Reals, repeating every nTotal strides.
     */
    void defineGhostChunk( Int nGFs, Int nTotal, Int nGhost ) {}

    /** Exchange the boundary data with the neighboring ranks.
     */
    bool exchangeBoundaries
    (
        void* left_ghost,  //!< Where to put the right edge of the left neighbor
        void* left_edge,   //!< Pointer to the left edge of our data
        void* right_edge,  //!< Pointer to the right edge of our data
        void* right_ghost, //!< Where to put the left edge of the right neighbor
        bool receiveEdges  //!< If to receive the edges (otherwise just send the ghosts)
        )
    {
        return true;
    }

    /** Wait the other ranks to have received our edges.
     */
    bool waitExchange ()
    {
        return true;
    }

    /** Nonblocking signal to all the ranks that we are quitting integration.
     */
    void abortExchange( Real* dummyData )
    {
    }
};

typedef DummyMPIWorld MPIWorld;

#endif // _MPI_DUMMY_WORLD_H_INCLUDED
