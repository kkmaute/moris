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

    size_t par_rank()
    {

        int tProcRank;
        MPI_Comm_rank(gMorisComm.get_global_comm(), &tProcRank);
        return (size_t)tProcRank;
    }

    //-----------------------------------------------------

    size_t par_size()
    {
        int tProcSize;
        MPI_Comm_size(gMorisComm.get_global_comm(), &tProcSize);
        return (size_t)tProcSize;
    }

    //-----------------------------------------------------

    void barrier()
    {
        MPI_Barrier(gMorisComm.get_global_comm());
    }

    //-----------------------------------------------------

    void Sum_All_Local_Int(
            const moris::uint & aLocalInput,
            moris::uint       & aGlobalSum)
    {
        MPI_Allreduce(&aLocalInput,&aGlobalSum,1,MPI_INT,MPI_SUM,MPI_COMM_WORLD);
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
    /*
    std::vector<uint>
    convert_mat_to_vect(Mat<uint> aMat)
    {
        uint n_cols = aMat.n_cols();
        uint n_rows = aMat.n_rows();
        uint* memPtr= mem_pointer(aMat);
        std::vector<uint> buffer(n_cols*n_rows+2,0);

        buffer[0] = n_cols;
        buffer[1] = n_rows;

        for(uint i = 0; i<n_cols*n_rows; i++)
        {
            buffer[i+2] = memPtr[i];
        }
        return buffer;
    }


    std::vector<int>
    convert_mat_to_vect(Mat<int> aMat)
    {
        uint n_cols = aMat.n_cols();
        uint n_rows = aMat.n_rows();
        int* memPtr= mem_pointer(aMat);
        std::vector<int> buffer(n_cols*n_rows+2,0);

        buffer[0] = n_cols;
        buffer[1] = n_rows;

        for(uint i = 0; i<n_cols*n_rows; i++)
        {
            buffer[i+2] = memPtr[i];
        }
        return buffer;
    }

    //-----------------------------------------------------

    std::vector<real>
    convert_mat_to_vect(Mat<real> aMat)
    {
        real n_cols = aMat.n_cols();
        real n_rows = aMat.n_rows();
        real* memPtr= mem_pointer(aMat);
        std::vector<real> buffer(n_cols*n_rows+2,0);

        buffer[0] = n_cols;
        buffer[1] = n_rows;
        for(auto i = 0; i<n_cols*n_rows; i++)
        {
            buffer[i+2] = memPtr[i];
        }

        return buffer;
    }
    //-----------------------------------------------------

    Mat<real>
    convert_vect_to_mat(std::vector<real> aBuffer)
    {
        real n_cols = aBuffer[0];
        real n_rows = aBuffer[1];
        Mat<real> tMat(n_rows,n_cols,0);

        real i = 2;
        for(real c =0; c<n_cols; c++)
        {
            for(real r =0; r<n_rows; r++)
            {
                tMat(r,c) = aBuffer[i];
                i++;
            }
        }
        return tMat;
    }

    //-----------------------------------------------------

    Mat<uint>
    convert_vect_to_mat(std::vector<uint> aBuffer)
    {
        uint n_cols = aBuffer[0];
        uint n_rows = aBuffer[1];
        Mat<uint> tMat(n_rows,n_cols,0);

        uint i = 2;
        for(uint c =0; c<n_cols; c++)
        {
            for(uint r =0; r<n_rows; r++)
            {
                tMat(r,c) = aBuffer[i];
                i++;
            }
        }
        return tMat;
    }

    Mat<int>
    convert_vect_to_mat(std::vector<int> aBuffer)
    {
        uint n_cols = aBuffer[0];
        uint n_rows = aBuffer[1];
        Mat<int> tMat(n_rows,n_cols,0);

        uint i = 2;
        for(uint c =0; c<n_cols; c++)
        {
            for(uint r =0; r<n_rows; r++)
            {
                tMat(r,c) = aBuffer[i];
                i++;
            }
        }
        return tMat;
    }


    //-----------------------------------------------------

    uint*
    convert_mat_to_array(Mat<uint> aMat)
    {
        uint n_cols  = aMat.n_cols();
        uint n_rows  = aMat.n_rows();
        uint* memPtr = mem_pointer(aMat);
        uint  tSize  = aMat.n_cols()*aMat.n_rows()+2;
        uint* buffer;
        //
        buffer = new uint[tSize];
        buffer[0] = n_cols;
        buffer[1] = n_rows;
        for(uint i = 0; i<n_cols*n_rows; i++)
        {
            buffer[i+2] = memPtr[i];
        }
        return buffer;
    }

    //-----------------------------------------------------

    Mat<uint>
    convert_array_to_mat(uint* aBuffer)
    {
        uint n_cols = aBuffer[0];
        uint n_rows = aBuffer[1];
        Mat<uint> tMat(n_rows,n_cols,0);

        uint i = 2;
        for(uint c =0; c<n_cols; c++)
        {
            for(uint r =0; r<n_rows; r++)
            {
                tMat(r,c) = aBuffer[i];
                i++;
            }
        }
        return tMat;
    }
    */

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
            Mat< uint >       & aProcDims,
            Mat< uint >       & aProcCoords,
            Mat< uint >       & aProcNeighbors )
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
            aProcNeighbors.set_size( 3, 1, MORIS_UINT_MAX );
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
                    aProcNeighbors( tCount ) = ( uint ) tRank;
                }

                // increment counter
                ++tCount;
            }
        }
        else if ( aNumberOfDimensions == 2 )
        {
            // assign memory for neighbors
            aProcNeighbors.set_size( 9, 1, MORIS_UINT_MAX );

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
                        aProcNeighbors( tCount ) = ( uint ) tRank;
                    }

                    // increment counter
                    ++tCount;
                }
            }
        }
        else if ( aNumberOfDimensions == 3 )
        {
            // assign memory for neighbors
            aProcNeighbors.set_size( 27, 1, MORIS_UINT_MAX );
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
