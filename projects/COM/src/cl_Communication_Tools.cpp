/*
 * cl_MPI_XTK.cpp
 *
 *  Created on: Jan 27, 2017
 *      Author: doble
 */
#include "cl_Communication_Tools.hpp" // COM/src

//#include <thread>         // std::this_thread::sleep_for
//#include <chrono>         // std::chrono::seconds

namespace moris
{
    MPI_Comm
    get_comm()
    {
        return gMorisComm.get_comm();
    }

    //------------------------------------------------------------------------------------------------------------------

    void comm_split(int color, int key, const std::string& aCommName)
    {
        MPI_Comm aNewComm;
        MPI_Comm_split(gMorisComm.get_comm(), color, key, &aNewComm);
        gMorisComm.add_communicator(aNewComm, aCommName);
    }

    //------------------------------------------------------------------------------------------------------------------

    void comm_join()
    {
        MORIS_ERROR(gMorisComm.mActiveCommunicator > 0, "Cannot comm_join(), active communicator is already the global one.");
        gMorisComm.remove_communicator(gMorisComm.mActiveCommunicator);
        gMorisComm.mActiveCommunicator -= 1;
    }

    //------------------------------------------------------------------------------------------------------------------

    moris_id par_rank()
    {

        int tProcRank;
        MPI_Comm_rank(gMorisComm.get_comm(), &tProcRank);
        return ( moris_id ) tProcRank;
    }

    //------------------------------------------------------------------------------------------------------------------

    moris_id par_size()
    {
        int tProcSize;
        MPI_Comm_size(gMorisComm.get_comm(), &tProcSize);
        return ( moris_id ) tProcSize;
    }

    //------------------------------------------------------------------------------------------------------------------

    void barrier(std::string aBarrierName)
    {
        MPI_Barrier(gMorisComm.get_comm());
    }

    //------------------------------------------------------------------------------------------------------------------
    void broadcast(uint & aMessage)
    {
        if(par_size() > 1)
        {
            MPI_Bcast(&aMessage,1,MPI_UNSIGNED,0,gMorisComm.get_comm());
        }
    }


    //------------------------------------------------------------------------------------------------------------------

    void broadcast(std::string & aMessage)
    {
        if(par_size() > 1)
        {
            char* cstr = strdup(aMessage.c_str());
            MPI_Bcast(cstr,aMessage.length() ,MPI_CHAR,0,gMorisComm.get_comm());
            aMessage.assign(cstr);
            delete [] cstr;
        }
    }

    //------------------------------------------------------------------------------------------------------------------

    void
    begin_log_communication(enum CommunicationType    CommType)
    {


        // Get count integer
        //    short_uint tCount = gMorisComm.get_comm()Counter((short_uint) CommType,0);

        // Get integer corresponding to enum
        //    uint tEnumInt  = (uint)CommType;

        //Print horizontal line
        //MORIS_LOG_INFO<< "========================================";
        //MORIS_LOG_INFO<< "Begin Communication Type: "<< tEnumInt <<" Count:"<< tCount;
    }

    //------------------------------------------------------------------------------------------------------------------

    void
    end_log_communication(enum CommunicationType    CommType)
    {
        // Get count integer
        //    short_uint tCount = gMorisComm.get_comm()Counter((short_uint) CommType,0);

        // Get integer corresponding to enum
        //    uint tEnumInt  = (uint)CommType;
        //MORIS_LOG_INFO<< "End Communication Type: "<< tEnumInt <<" Count:"<< tCount;
        //MORIS_LOG_INFO<< "======================================== \n";

        // Advance counter for this type of communication
//        gMorisComm.get_comm()Counter((short_uint)CommType,0)++;
    }

    //------------------------------------------------------------------------------------------------------------------


    //------------------------------------------------------------------------------------------------------------------

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
                MPI_Gather(&aMessageInt,1,MPI_INT,recvbuf,1,MPI_INT,0,gMorisComm.get_comm());
                tMaxValue = *std::max_element(recvbuf, recvbuf+tProcSize);
                aMessage = (uint)tMaxValue;
            }
            if(tProcRank > 0)
            {
                MPI_Gather(&aMessageInt,1,MPI_INT,recvbuf,1,MPI_INT,0,gMorisComm.get_comm());
            }

            broadcast(aMessage);
        }

        return aMessage;
    }

    //------------------------------------------------------------------------------------------------------------------
    void create_proc_cart(
            const uint        & aNumberOfDimensions,
            Matrix < DDUMat > & aProcDims,
            Matrix < DDUMat > & aProcCoords,
            Matrix < IdMat >  & aProcNeighbors )
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
        MPI_Dims_create( par_size(),
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

        MPI_Cart_create( gMorisComm.get_comm(),
                         aNumberOfDimensions,
                         tDims,
                         tPeriods,
                         tReorder,
                         &tNewComm );

        //Define coordinates for the current processor
        MPI_Cart_coords( tNewComm,
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

                if ( tNeighborCoords[ 0 ] >= 0 && tNeighborCoords[ 0 ] < tDims[0] )
                {
                    // get rank of neighbor proc
                    tError = MPI_Cart_rank( tNewComm,
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
                        tError = MPI_Cart_rank( tNewComm,
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
                            tError = MPI_Cart_rank( tNewComm,
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
            MORIS_ERROR( false, "create_proc_cart: invalid number of dimensions" );
        }
    }

    //------------------------------------------------------------------------------------------------------------------
    void create_proc_dims(
            const uint           & aNumberOfDimensions,
            Matrix < DDUMat >    & aMeshDims,
            Matrix < DDUMat >    & aProcDimsTemp)

    {

        //This function determines the processor dimensions based on minimizing mesh grid interfaces.
        //Eg. If the mesh is 100x10 elements, and we have 10 processors, the processor layout will be 10x1.

        aProcDimsTemp.set_size( 3, 1 );

        real tInterfaceCount=0;



        //1D Processor Grid
        if (aNumberOfDimensions==1)
        {
            aProcDimsTemp(0) = par_size();
            aProcDimsTemp(1) = 1;
            aProcDimsTemp(2) = 1;
        }

        //2D Processor Grid
        else if (aNumberOfDimensions==2)
        {
            //Iterating through "i" (x) direction processors
            for (uint i=1; i <= (uint) par_size(); i++)
            {

                //initiallizing interface count value
                if (i==1)
                {
                    tInterfaceCount = aMeshDims(1)*(i-1) + aMeshDims(0)*(par_size()/i-1);
                    //Assigning processor grid dimensions
                    aProcDimsTemp(0) = i;
                    aProcDimsTemp(1) = par_size()/i;
                    aProcDimsTemp(2) = 1;
                    //std::cout<<"interface"<<tInterfaceCount<<std::endl;
                    //std::cout<<"i"<<aProcDimsTemp(0)<<std::endl;
                    //std::cout<<"j"<<aProcDimsTemp(1)<<std::endl;
                }

                //Otherwise, is i a factor of the number of processors?
                else if ((uint) par_size()%i == 0)
                {
                    //Is the last calculated interface count larger than the current iterations?
                    if (tInterfaceCount > (aMeshDims(1)*(i-1) + aMeshDims(0)*(par_size()/i-1)))
                    {
                        tInterfaceCount = aMeshDims(1)*(i-1) + aMeshDims(0)*(par_size()/i-1);
                        //Assigning processor grid dimensions
                        aProcDimsTemp(0) = i;
                        aProcDimsTemp(1) = par_size()/i;
                        aProcDimsTemp(2) = 1;
                        //std::cout<<"interface"<<tInterfaceCount<<std::endl;
                        //std::cout<<"i"<<aProcDimsTemp(0)<<std::endl;
                        //std::cout<<"j"<<aProcDimsTemp(1)<<std::endl;
                    }
                    else
                    {
                        break;
                    }
                }
            }
        }

        //3D Processor Grid
        else if (aNumberOfDimensions == 3)
        {
            // Iterating through i-direction processor possibilities
            for (uint i = 1; i <= (uint) par_size(); i++)
            {
                // Is "i" a factor of total processors used?
                if((uint) par_size()%i == 0)
                {
                    // Iterating through j-direction processors possibilities
                    for (uint j = 1; j <= (uint) par_size()/i; j++)
                    {
                        // Initializing interface count value
                        if (i == 1 && j == 1)
                        {
                            tInterfaceCount = aMeshDims(0)*aMeshDims(1)*(par_size()/(i*j)-1) +
                                    aMeshDims(0)*aMeshDims(2)*(j-1) + aMeshDims(1)*aMeshDims(2)*(i-1);

                            //Repace processor grid dimensions
                            aProcDimsTemp(0) = i;
                            aProcDimsTemp(1) = j;
                            aProcDimsTemp(2) = par_size()/(i*j);
                        }

                        // Is i*j a factor of processor count?
                        else if ((uint) par_size()%(i*j) == 0)
                        {

                            //If the old interface count is higher than the current iteration
                            if (tInterfaceCount > (aMeshDims(0)*aMeshDims(1)*(par_size()/(i*j)-1) +
                                    aMeshDims(0)*aMeshDims(2)*(j-1) + aMeshDims(1)*aMeshDims(2)*(i-1)))
                            {

                                //Replace interface count
                                tInterfaceCount = aMeshDims(0)*aMeshDims(1)*(par_size()/(i*j)-1) +
                                        aMeshDims(0)*aMeshDims(2)*(j-1) + aMeshDims(1)*aMeshDims(2)*(i-1);

                                //Repace processor grid dimensions
                                aProcDimsTemp(0) = i;
                                aProcDimsTemp(1) = j;
                                aProcDimsTemp(2) = par_size()/(i*j);
                            }
                        }

                    }
                }

            }
        }
    }

    //------------------------------------------------------------------------------------------------------------------
    void create_proc_cart_2(
            const uint        & aDecompMethod,
            const uint        & aNumberOfDimensions,
            Matrix < DDUMat > & aMeshDims,
            Matrix < DDUMat > & aProcDims,
            Matrix < DDUMat > & aProcCoords,
            Matrix < IdMat >  & aProcNeighbors )
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

        //determining an appropriate processor layout to minimize grid interface
        int tDims[3] = { 0, 0, 0 };

        //Original Decomp Method. Minimizes processor interfaces

        switch (aDecompMethod)
        {

            case 1 : // MPI Default Decomp Method
            {
                MPI_Dims_create( par_size(),
                        aNumberOfDimensions,
                        tDims );
                break;
            }

            case 2 :// Decomposition method to minimize mesh interface
            {
                Matrix < DDUMat >     tProcDimsTemp;
                create_proc_dims(aNumberOfDimensions,
                        aMeshDims,
                        tProcDimsTemp);

                tDims[0]=(int)tProcDimsTemp(0);
                tDims[1]=(int)tProcDimsTemp(1);
                tDims[2]=(int)tProcDimsTemp(2);
                break;
            }

            default:
            {
                MORIS_ERROR( false, "create_proc_cart: Undefined decomposition method" );
            }
        }

        // No periodic boundary conditions are needed
        int tPeriods[3]  = { 0, 0, 0 };

        // Ranks of the processors stay unaltered in the new communicator "new_comm"
        int tReorder = 1;
        int tMyRank  = ( int ) par_rank();
        int tCoords[3]  = { 0, 0, 0 };

        // create new communicator
        // New communicator for cart and coordinates
        MPI_Comm tNewComm;

        MPI_Cart_create( gMorisComm.get_comm(),
                aNumberOfDimensions,
                tDims,
                tPeriods,
                tReorder,
                &tNewComm );

        //Define coordinates for the current processor
        MPI_Cart_coords( tNewComm,
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

                if ( tNeighborCoords[ 0 ] >= 0 && tNeighborCoords[ 0 ] < tDims[0] )
                {
                    // get rank of neighbor proc
                    tError = MPI_Cart_rank( tNewComm,
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
                        tError = MPI_Cart_rank( tNewComm,
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
                            tError = MPI_Cart_rank( tNewComm,
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
            MORIS_ERROR( false, "create_proc_cart: invalid number of dimensions" );
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

    void
    all_gather_cell_of_str(
        Cell<std::string> const & aCellToGather,
        Cell<Cell<std::string>> & aGatheredCells,
        moris_index aTag,
        moris_index aBaseProc)
    {
        
        MPI_Request tRequest;


    std::vector<char> cstrings;
    cstrings.reserve(aCellToGather.size());
    for(std::string s: aCellToGather.data())
    {
        for(size_t i = 0; i < strlen(s.c_str()); ++i)
        {
            cstrings.push_back(s.c_str()[i]);
        }
        // terminate str
        cstrings.push_back('\0');
    }

    MPI_Isend(cstrings.data(), cstrings.size(), MPI_CHAR, aBaseProc, aTag, moris::get_comm(), &tRequest);

    barrier();

    // on base rank go ahead and receive the data
    if(par_rank() == aBaseProc)
    {
        aGatheredCells.resize(par_size());
        for(int i = 0; i < par_size(); i++)
        {
            MPI_Status tStatus;
            MPI_Probe(i, aTag, moris::get_comm(), &tStatus);

            //    MORIS_ERROR(tExists,"Trying to receive a message that does not exists");

            int tLength = 0;
            MPI_Get_count(
                    &tStatus,
                    MPI_CHAR,
                    &tLength);
            
            char* tChars = new char[tLength];

            MPI_Recv(
                tChars,
                tLength,
                MPI_CHAR,
                i,
                aTag,
                moris::get_comm(),
                &tStatus);

            moris::uint tCellIndex = 0;
            aGatheredCells(i).push_back("");
            for (int  j = 0; j < tLength; j++)
            {
                if(tChars[j] == '\0')
                {
                    tCellIndex++;
                    aGatheredCells(i).push_back("");
                }
                else
                {
                    aGatheredCells(i)(tCellIndex).push_back(tChars[j]);
                }
                
            }
            aGatheredCells(i).pop_back();
        }
    }

    barrier();
    }

//------------------------------------------------------------------------------
}
