/*
 * cl_MPI_XTK.cpp
 *
 *  Created on: Jan 27, 2017
 *      Author: doble
 */
#include "cl_Communication_Tools.hpp" // COM/src

namespace moris
{
    //-----------------------------------------------------

    moris_id par_rank()
    {

        int tProcRank;
        MPI_Comm_rank(gMorisComm.get_global_comm(), &tProcRank);
        return ( moris_id ) tProcRank;
    }

    //-----------------------------------------------------

    moris_id par_size()
    {
        int tProcSize;
        MPI_Comm_size(gMorisComm.get_global_comm(), &tProcSize);
        return ( moris_id ) tProcSize;
    }

    //-----------------------------------------------------

    void barrier()
    {
        MPI_Barrier(gMorisComm.get_global_comm());
    }

    //-----------------------------------------------------

    void sum_all_int(
            const moris::uint & aLocalInput,
            moris::uint       & aGlobalSum)
    {
        MPI_Allreduce(&aLocalInput,&aGlobalSum,1,MPI_INT,MPI_SUM,MPI_COMM_WORLD);
    }

    void
    sum_all_real(
            const moris::real & aLocalInput,
            moris::real       & aGlobalSum)
    {
        MPI_Allreduce(&aLocalInput,&aGlobalSum,1,get_comm_datatype(aLocalInput),MPI_SUM,MPI_COMM_WORLD);
    }

    //-----------------------------------------------------

    void broadcast(uint & aMessage)
    {
        if(par_size() > 1)
        {
            MPI_Bcast(&aMessage,1,MPI_UNSIGNED,0,gMorisComm.get_global_comm());
        }
    }


    //-----------------------------------------------------

    void broadcast(std::string & aMessage)
    {
        if(par_size() > 1)
        {
            char* cstr = strdup(aMessage.c_str());
            MPI_Bcast(cstr,aMessage.length() ,MPI_CHAR,0,gMorisComm.get_global_comm());
            aMessage.assign(cstr);
            delete [] cstr;
        }
    }

    //---------------------------------------------

    void
    begin_log_communication(enum CommunicationType    CommType)
    {


        // Get count integer
        //    short_uint tCount = gMorisComm.get_global_comm()Counter((short_uint) CommType,0);

        // Get integer corresponding to enum
        //    uint tEnumInt  = (uint)CommType;

        //Print horizontal line
        //MORIS_LOG_INFO<< "========================================";
        //MORIS_LOG_INFO<< "Begin Communication Type: "<< tEnumInt <<" Count:"<< tCount;
    }

    //-----------------------------------------------------

    void
    end_log_communication(enum CommunicationType    CommType)
    {
        // Get count integer
        //    short_uint tCount = gMorisComm.get_global_comm()Counter((short_uint) CommType,0);

        // Get integer corresponding to enum
        //    uint tEnumInt  = (uint)CommType;
        //MORIS_LOG_INFO<< "End Communication Type: "<< tEnumInt <<" Count:"<< tCount;
        //MORIS_LOG_INFO<< "======================================== \n";

        // Advance counter for this type of communication
//        gMorisComm.get_global_comm()Counter((short_uint)CommType,0)++;
    }

    //-----------------------------------------------------


    //-----------------------------------------------------

    uint gather_value_and_bcast_max( uint aMessage )
    {
        int tProcSize = (int)par_size();
        int tProcRank = (int)par_rank();
        if(tProcSize >1)
        {
            int aMessageInt = (int)aMessage;
            int *recvbuf;
            recvbuf = (int *)malloc(tProcSize*1*sizeof(int));
            int tMaxValue = 0;
            if(tProcRank == 0)
            {
                MPI_Gather(&aMessageInt,1,MPI_INT,recvbuf,1,MPI_INT,0,gMorisComm.get_global_comm());
                tMaxValue = *std::max_element(recvbuf, recvbuf+tProcSize);
                aMessage = (uint)tMaxValue;
            }
            if(tProcRank > 0)
            {
                MPI_Gather(&aMessageInt,1,MPI_INT,recvbuf,1,MPI_INT,0,gMorisComm.get_global_comm());
            }

            broadcast(aMessage);
        }

        return aMessage;
    }

    //-----------------------------------------------------
    void
    create_proc_cart(
            const uint        & aNumberOfDimensions,
            Matrix < DDUMat >       & aProcDims,
            Matrix < DDUMat >       & aProcCoords,
            Matrix < IdMat >        & aProcNeighbors )
    {

        /* Table for aProcNeighbors
         *
         *  neighbors for 1D case
         *
         *  .--------------.
         *  |  0 | me |  2 |   o --> i
         *  '--------------'
         *
         *  neighbors for 2D case
         *
         *  .--------------.
         *  |  6 |  7 |  8 |   j
         *  |--------------|   ^
         *  |  3 | me |  5 |   |
         *  |--------------|   o --> i
         *  |  0 |  1 |  2 |
         *  '--------------'
         *
         *  neighbors for 3D case
         *
         *  .--------------.
         *  |  6 |  7 |  8 |   j
         *  |--------------|   ^
         *  |  3 |  4 |  5 |   |
         *  |--------------|   o --> i
         *  |  0 |  1 |  2 |   k = -1
         *  '--------------'
         *
         *  .--------------.
         *  | 15 | 16 | 17 |   j
         *  |--------------|   ^
         *  | 12 | me | 14 |   |
         *  |--------------|   o --> i
         *  |  9 | 10 | 11 |   k =  0
         *  '--------------'
         *
         *  .--------------.
         *  | 24 | 25 | 26 |   j
         *  |--------------|   ^
         *  | 21 | 22 | 23 |   |
         *  |--------------|   o --> i
         *  | 18 | 19 | 20 |   k =  1
         *  '--------------'
         */

        // Number of processors in i, j and k-direction
        int tDims[3] = { 0, 0, 0 };

        // Creates a grid of processors
        MPI_Dims_create(
                par_size(),
                aNumberOfDimensions,
                tDims );

        // No periodic boundary conditions are needed
        int tPeriods[3]  = { 0, 0, 0 };

        // Ranks of the processors stay unaltered in the new communicator "new_comm"
        int tReorder = 1;
        int tMyRank  = ( int ) par_rank();
        int tCoords[3]  = { 0, 0, 0 };

        // create new communicator
        // New communicator for cart and coordinates
        MPI_Comm tNewComm;

        MPI_Cart_create(
                gMorisComm.get_global_comm(),
                aNumberOfDimensions,
                tDims,
                tPeriods,
                tReorder,
                &tNewComm );

        //Define coordinates for the current processor
        MPI_Cart_coords(
                tNewComm,
                tMyRank,
                3, // must be equal to length of coords
                tCoords );

        int tNeighborCoords[ 3 ] = { 0, 0, 0 };

        // counter for neighbor
        uint tCount = 0;

        // copy proc coordinates
        aProcCoords.set_size( aNumberOfDimensions, 1 );
        for( uint k=0; k<aNumberOfDimensions; ++k )
        {
            aProcCoords( k ) = ( uint ) tCoords[ k ];
        }

        // copy dims
        aProcDims.set_size( aNumberOfDimensions, 1 );
        for( uint k=0; k<aNumberOfDimensions; ++k )
        {
            aProcDims( k ) = ( uint ) tDims[ k ];
        }

        if ( aNumberOfDimensions == 1 )
        {
            // assign memory for neighbors
            aProcNeighbors.resize( 3, 1 );
            aProcNeighbors.fill( gNoProcID );
            // loop in i-direction
            for ( int i=-1; i<=1; ++i )
            {
                // calculate coordinate of neighbor
                tNeighborCoords[ 0 ] = tCoords [ 0 ] + i ;

                int tRank;

                // error from MPI
                int tError = 1;

                if (       tNeighborCoords[ 0 ] >= 0
                        && tNeighborCoords[ 0 ] < tDims[0] )
                {
                    // get rank of neighbor proc
                    tError = MPI_Cart_rank(
                            tNewComm,
                            tNeighborCoords,
                            &tRank );
                }
                // if getting the rank was successful
                if ( tError == 0 )
                {
                    // add rank to neighbor list
                    aProcNeighbors( tCount ) = ( moris_id ) tRank;
                }

                // increment counter
                ++tCount;
            }
        }
        else if ( aNumberOfDimensions == 2 )
        {
            // assign memory for neighbors
            aProcNeighbors.resize( 9, 1 );
            aProcNeighbors.fill( gNoProcID );
            // loop over all j
            for ( int j=-1; j<= 1; ++j )
            {
                for ( int i=-1; i<=1; ++i )
                {
                    // calculate coordinate of neighbor
                    tNeighborCoords[ 0 ] = tCoords [ 0 ] + i ;
                    tNeighborCoords[ 1 ] = tCoords [ 1 ] + j ;

                    int tRank;

                    // error from MPI
                    int tError = 1;

                    // test if coordinate exists
                    if (       tNeighborCoords[ 0 ] >= 0
                            && tNeighborCoords[ 1 ] >= 0
                            && tNeighborCoords[ 0 ] < tDims[0]
                            && tNeighborCoords[ 1 ] < tDims[1] )
                    {
                        // get rank of neighbor proc
                        tError = MPI_Cart_rank(
                                tNewComm,
                                tNeighborCoords,
                                &tRank );
                    }

                    // if getting the rank was successful
                    if ( tError == 0 )
                    {
                        // add rank to neighbor list
                        aProcNeighbors( tCount ) = ( moris_id ) tRank;
                    }

                    // increment counter
                    ++tCount;
                }
            }
        }
        else if ( aNumberOfDimensions == 3 )
        {
            // assign memory for neighbors
            aProcNeighbors.resize( 27, 1 );
            aProcNeighbors.fill( gNoProcID );

            for( int k=-1; k<=1; ++k )
            {
                for ( int j=-1; j<= 1; ++j )
                {
                    for ( int i=-1; i<=1; ++i )
                    {
                        // calculate coordinate of neighbor
                        tNeighborCoords[ 0 ] = tCoords [ 0 ] + i ;
                        tNeighborCoords[ 1 ] = tCoords [ 1 ] + j ;
                        tNeighborCoords[ 2 ] = tCoords [ 2 ] + k ;

                        int tRank;

                        // error from MPI
                        int tError = 1;

                        // test if coordinate can exist
                        if (       tNeighborCoords[ 0 ] >= 0
                                && tNeighborCoords[ 1 ] >= 0
                                && tNeighborCoords[ 2 ] >= 0
                                && tNeighborCoords[ 0 ] < tDims[0]
                                && tNeighborCoords[ 1 ] < tDims[1]
                                && tNeighborCoords[ 2 ] < tDims[2] )
                        {
                            // get rank of neighbor proc
                            tError = MPI_Cart_rank(
                                    tNewComm,
                                    tNeighborCoords,
                                    &tRank );
                        }

                        // if getting the rank was successful
                        if ( tError == 0 )
                        {
                            // add rank to neighbor list
                            aProcNeighbors( tCount ) = ( uint ) tRank;
                        }

                        // increment counter
                        ++tCount;
                    }
                }
            }
        }
        else
        {
            MORIS_ASSERT(
                    false ,
                    "create_proc_cart: invalid number of dimensions ");
        }
    }

//------------------------------------------------------------------------------

    int
    create_comm_tag ( const int & aSource, const int & aTarget )
    {
        int tMin = ( aSource < aTarget ) ? ( aSource ) : ( aTarget );
        int tMax = ( aSource > aTarget ) ? ( aSource ) : ( aTarget );

        if ( tMin == aTarget )
        {
            return 4*( tMax*par_size() + tMin );
        }
        else
        {
            return 4*( tMax*par_size() + tMin )  + 2;
        }
        return tMax*par_size() + tMin ;
    }

//------------------------------------------------------------------------------
}
