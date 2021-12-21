/*
 * cl_mpi_xtk.hpp
 *
 *  Created on: Jan 27, 2017
 *      Author: doble
 */
#ifndef SRC_COMM_CL_COMMUNICATION_TOOLS_HPP_
#define SRC_COMM_CL_COMMUNICATION_TOOLS_HPP_

// Standard headers
#include <iostream>
#include <string>

#include "cl_Communication_Enums.hpp" // COM/src
#include "cl_Communication_Manager.hpp" // COM/src
#include "cl_Matrix.hpp"
#include "linalg_typedefs.hpp"
#include "fn_isvector.hpp"

//#include "cl_Bitset.hpp" // CON/src
#include "cl_Cell.hpp" // CON/src

// Mesh specific headers
//#include "cl_Mesh_Enums.hpp" // MTK/src

// Externally Defined Global Communicator
extern moris::Comm_Manager gMorisComm;

namespace moris
{

    // value to be used if a proc has no neighbor at the
    // specified position
    const moris_id   gNoProcID    =  MORIS_SINT_MAX;

    //------------------------------------------------
    //  GENERAL MPI FUNCTIONS

    /**
     * Gets the current global communicator
     *
     * @return MPI_Comm, current global communicator
     */
    MPI_Comm
    get_comm();

    /**
     * Splits the current communicator
     *
     * @param color control of subset assignment
     * @param key control of rank assignment
     */
    void comm_split(int color, int key, const std::string& aCommName);

    /**
     * Returns to global communicator
     */
    void comm_join();

    /**
     * Returns the current processor rank
     *
     * @return moris_id, current processor rank
     */
    moris::moris_id
    par_rank();

    /**
     * Returns the size of the processor pool
     *
     * @return moris_id, size of the processor pool
     */
    moris::moris_id
    par_size();

    /**
     * Holds all processors
     *
     * @param aBarrierName string naming this barrier
     */
    void barrier(std::string aBarrierName = " ");

    /**
     * Broadcast a Message to all Procs
     *
     * @param[in] aMessage      A uint
     *
     * @param[out] aMessage     The message of Proc 0
     *
     */
    void broadcast(uint & aMessage);

    /**
     * Broadcast a Message to all Procs
     *
     * @param[in] aMessage      A string to broadcast
     *
     * @param[out] aMessage     The message of Proc 0
     *
     */
    void broadcast(std::string & aMessage);

    /**
     * Broadcast a Message to all Procs
     *
     * @param[in] aMessage      A bitset
     *
     * @param[out] aMessage     The message of Proc 0
     *
     */
    //    void broadcast(BoostBitset & aMessage);

    /*
     *  Scatters information from processor 0 to all processors
     */
    moris::size_t
    scatter(std::vector<moris::size_t>  & aMessage);

    moris::uint gather_value_and_bcast_max( moris::uint aMessage );

    /*
     * Prints communication header to log
     * @param[in]  CommType - Enum for communication type to be documented
     */
    void begin_log_communication(enum CommunicationType    CommType);

    /*
     * Prints communication header to log
     * @param[in]  CommType - Enum for communication type to be documented
     */
    void end_log_communication(enum CommunicationType    CommType);

    /*
     * setup log sink should only be done once by proc 0
     * @param[in]  CommType - Enum for communication type to be documented
     */
    void setup_log_sink(std::string    aCommLogName);

    /*
     * Generates filename which is used as a filter for the specific log sink
     *
     * The filename is morisroot/log/comm_log/EnumInt_Count.log
     */
    std::string
    generate_filename(enum CommunicationType    CommType);

    //------------------------------------------------------------------------------
    /**
     *
     * @brief creates a proc cart similar to MPI_CART, depending
     *        on the number of dimensions specified.
     *
     * @param[in]  aDecompMethod         Decomposition method used. 0=User Defined.
     *                   1=Min Proc Interface (MPI Cart). 2=Min Mesh Interface.
     * @param[in]  aNumberOfDimensions   Must be 1, 2 or 3.
     * @param[in,out] aProcDims          Dimensions of generated cart. Input used
     *                   for aDecompMethod = 0 or 2. Output used for aDecompMethod = 1.
     * @param[out] aProcCoords           Coordinates of current processor.
     * @param[out] aProcNeighbors        Neighbors of current proc.
     *                                   Contains UINT_MAX if no neighbor exists.
     *
     *  neighbor pattern as returned to aProcNeighbors:
     *
     *  <table>
     *  <caption> 1D case </caption>
     *  <tr><th> 0 </th><th> me </th><th> 2 </th></th></tr>
     *  </table>
     *  ( cols: i-direction )
     *
     *
     *  <table>
     *  <caption> 2D case</caption>
     *  <tr><th> 6 </th><th> 7 </th><th> 8 </th></th></tr>
     *  <tr><th> 3 </th><th> me </th><th> 5 </th></th></tr>
     *  <tr><th> 0 </th><th> 1 </th><th> 2 </th></th></tr>
     *  </table>
     *  ( cols: i-direction, rows: j-direction )
     *
     *
     *  <table>
     *  <caption> 3D case</caption>
     *  <tr>
     *  <th>
     *  <table>
     *  <caption> layer k=-1</caption>
     *  <tr><th> 6 </th><th> 7 </th><th> 8 </th></th></tr>
     *  <tr><th> 3 </th><th> 4 </th><th> 5 </th></th></tr>
     *  <tr><th> 0 </th><th> 1 </th><th> 2 </th></th></tr>
     *  </table>
     *  </th>
     *  <th>
     *  <table>
     *  <caption> layer k=0</caption>
     *  <tr><th> 15 </th><th> 16 </th><th> 17 </th></th></tr>
     *  <tr><th> 12 </th><th> me </th><th> 14 </th></th></tr>
     *  <tr><th> 9 </th><th> 10 </th><th> 11 </th></th></tr>
     *  </table>
     *  </th>
     *  <th>
     *  <table>
     *  <caption> layer k=+1</caption>
     *  <tr><th> 24 </th><th> 25 </th><th> 26 </th></th></tr>
     *  <tr><th> 21 </th><th> 22 </th><th> 23 </th></th></tr>
     *  <tr><th> 18 </th><th> 19 </th><th> 20 </th></th></tr>
     *  </table>
     *  </th>
     *  </tr>
     *  </table>
     *  ( cols: i-direction, rows: j-direction )
     *
     */

    void create_proc_cart(
            const uint              & aDecompMethod,
            const uint              & aNumberOfDimensions,
            Matrix < DDUMat >       & aProcDims,
            Matrix < DDUMat >       & aProcCoords,
            Matrix < IdMat >        & aProcNeighbors );


    //------------------------------------------------------------------------------

    /**
     *
     * @brief                 returns an MPI enum defining the
     *                        data type that is to be communicated.
     *
     * @param[in] aSample     primitive data type with arbitrary value
     *
     * see also http://mpitutorial.com/tutorials/mpi-send-and-receive/
     */
    template < typename T > MPI_Datatype
    get_comm_datatype( const T & aSample );

    // moris::lint (32-bit)

    inline MPI_Datatype
    get_comm_datatype( const int & aSample )
    {
        return MPI_INT;
    }

    // moris::lint (64-bit)

    inline MPI_Datatype
    get_comm_datatype( const long int & aSample )
    {
        return MPI_LONG;
    }

    // moris::uint (32-bit)

    inline MPI_Datatype
    get_comm_datatype( const unsigned int & aSample )
    {
        return MPI_UNSIGNED;
    }

    // moris::uint (64-bit)

    inline MPI_Datatype
    get_comm_datatype( const long unsigned int & aSample )
    {
        return MPI_UNSIGNED_LONG;
    }

    // moris::real (32-bit)

    inline MPI_Datatype
    get_comm_datatype( const double & aSample )
    {
        return MPI_DOUBLE;
    }

    // moris::real (64-bit)

    inline MPI_Datatype
    get_comm_datatype( const long double & aSample )
    {
        return MPI_LONG_DOUBLE;
    }

    // MPI_CXX_DOUBLE_COMPLEX is supported since MPI-3.0
#ifdef MPI_CXX_DOUBLE_COMPLEX
    // moris::cplx (32-bit)

    inline MPI_Datatype
    get_comm_datatype( const std::complex<double> & aSample )
    {
        return MPI_CXX_DOUBLE_COMPLEX;
    }
#endif

    // MPI_CXX_LONG_DOUBLE_COMPLEX is supported since MPI-3.0
#ifdef MPI_CXX_LONG_DOUBLE_COMPLEX
    // moris::cplx (64-bit)

    inline MPI_Datatype
    get_comm_datatype( const std::complex<long double> & aSample )
    {
        return MPI_CXX_LONG_DOUBLE_COMPLEX;
    }
#endif

    //------------------------------------------------------------------------------

    /*
     * @brief Determines sum of all local values on all processors
     *
     * @param[in] aLocalInput     local value
     *
     * returns global sum
     */
    template <typename Data>
    Data
    sum_all( const Data & aLocalInput )
    {
        Data tGlobalSum;
        MPI_Allreduce(&aLocalInput,&tGlobalSum,1,get_comm_datatype(aLocalInput),MPI_SUM,MPI_COMM_WORLD);
        return tGlobalSum;
    }

    //------------------------------------------------------------------------------

    /*
     * @brief Determines sum of all components of a matrix on all processors
     *
     * @param[in] aLocalInput     local matrix
     *
     * returns matrix with global sums
     */

    template <typename E>
    E
    sum_all_matrix( const E & aLocalMatrix )
    {
        E tGlobalMatrix(aLocalMatrix.n_rows(),aLocalMatrix.n_cols());

        void* aLocalPtr  = (void*)aLocalMatrix.data();
        void* aGlobalPtr = (void*)tGlobalMatrix.data();

        int tNumElems = aLocalMatrix.numel();

        MPI_Datatype tType = get_comm_datatype ( ( typename E::Data_Type ) 0 );

        MPI_Allreduce(aLocalPtr,aGlobalPtr,tNumElems,tType,MPI_SUM,MPI_COMM_WORLD);

        return tGlobalMatrix;
    }

    //------------------------------------------------------------------------------

    /*
     * @brief Determines maximal value of all local values
     *
     * @param[in] aLocalInput     local value
     * @param[in] aGlobalMax      reference to maximal global value
     */
    template <typename Data>
    Data
    max_all( const Data & aLocalInput )
    {
        Data aGlobalMax;
        MPI_Allreduce(&aLocalInput,&aGlobalMax,1,get_comm_datatype(aLocalInput),MPI_MAX,MPI_COMM_WORLD);
        return aGlobalMax;
    }

    //------------------------------------------------------------------------------

    /*
     * @brief Determines maximal value of all local values
     *
     * @param[in] aLocalInput     local value
     * @param[in] aGlobalMax      reference to maximal global value
     */
    template <typename Data>
    Data
    min_all( const Data & aLocalInput )
    {
        Data aGlobalMin;
        MPI_Allreduce(&aLocalInput,&aGlobalMin,1,get_comm_datatype(aLocalInput),MPI_MIN,MPI_COMM_WORLD);
        return aGlobalMin;
    }

    /*
    * MPI_LAND is a predefined reduction operation that will calculate the logical and of the values reduced.
    * https://www.rookiehpc.com/mpi/docs/mpi_land.php
    */
    
    bool
    all_land(bool aMyBool);

    //------------------------------------------------------------------------------
    /**
     *
     * @brief                 Communication tags are used within MPI to label
     *                        messages that are exchanged between to procs.
     *                        If this is not wanted, one can use MPI_ANY_TAG
     *                        instead.
     *
     * @param[in] aSource     rank of sending proc
     * @param[in] aTraget     rank of receiving proc
     *
     */
    int
    create_comm_tag ( const int & aSource, const int & aTarget ) ;

    //------------------------------------------------------------------------------

    /**
     *
     * @brief                 sends a value to each proc and receives values
     *                        from each proc
     *
     * @param[in] aSend       values to send to each individual proc
     * @param[in] aReceive    values to receive from each individual proc
     */
    template <typename T> void
    communicate_scalars( 
            const Matrix < DDUMat > & aCommunicationList,
            const Matrix< T >      & aScalarsToSend,
            Matrix< T >            & aScalarsToReceive )
    {
        // only call this when we are in parallel mode
        if ( par_size() > 1 )
        {
            // get number of procs to communicate with
            uint tNumberOfProcs = aCommunicationList.length();

            // get my ID
            uint tMyRank = par_rank();

            aScalarsToReceive.set_size( tNumberOfProcs, 1, 0 );

            // Allocate memory for status/request vector
            // These vectors will be used to determine if the exchange has been completed across all processors
            MPI_Status*  tSendStatus  = ( MPI_Status*  ) alloca( sizeof( MPI_Status  ) * tNumberOfProcs );
            MPI_Status*  tRecvStatus  = ( MPI_Status*  ) alloca( sizeof( MPI_Status  ) * tNumberOfProcs );
            MPI_Request* tSendRequest = ( MPI_Request* ) alloca( sizeof( MPI_Request ) * tNumberOfProcs );
            MPI_Request* tRecvRequest = ( MPI_Request* ) alloca( sizeof( MPI_Request ) * tNumberOfProcs );

            // determine MPI datatype
            MPI_Datatype tType = get_comm_datatype ( ( typename Matrix< T >::Data_Type ) 0 );

            // loop over all procs
            for( uint k=0; k<tNumberOfProcs; ++k )
            {
                // only communicate if proc neighbor exists and is not me
                if (    ( aCommunicationList( k ) < MORIS_UINT_MAX ) &&
                        ( aCommunicationList( k ) != tMyRank ) )
                {
                    // create tag
                    int tSendTag = create_comm_tag ( tMyRank, aCommunicationList( k ) );

                    // send values
                    MPI_Isend(
                            &aScalarsToSend( k ),
                            1,
                            tType,
                            aCommunicationList( k ),
                            tSendTag,
                            gMorisComm.get_global_comm(),
                            &tSendRequest[ k ] );

                    // create tag
                    int tRecvTag = create_comm_tag ( aCommunicationList( k ), tMyRank );

                    // receive values
                    MPI_Irecv(
                            &aScalarsToReceive( k ),
                            1,
                            tType,
                            aCommunicationList( k ),
                            tRecvTag,
                            gMorisComm.get_global_comm(),
                            & tRecvRequest[ k ] );

                    MPI_Wait( &tSendRequest[ k ], &tSendStatus[ k ] );
                    MPI_Wait( &tRecvRequest[ k ], &tRecvStatus[ k ] );
                }
            }
        }
    }

    //------------------------------------------------------------------------------
    /**
     *
     * @brief                 Receives a moris::Mat<T> from the proc aSource.
     *                        Sending proc must call send_mat_to_proc.
     *                        Does nothing if MORIS is run in serial.
     *
     * @param[in] aMatrix     matrix to be communicated
     * @param[in] aSource     rank of sending proc
     *
     */
    template <typename T> void
    communicate_mats(
            const Matrix < IdMat >     & aCommunicationList,
            const Cell< Matrix< T > >  & aMatsToSend,
            Cell< Matrix< T > >        & aMatsToReceive )
    {
        moris_id tParSize = par_size();
        // only call this when we are in parallel mode
        if ( tParSize > 1 )
        {
            // get number of procs to communicate with
            moris_id tNumberOfProcs = aCommunicationList.length();

            // get my ID
            moris_id tMyRank = par_rank();

            // Allocate memory for status/request vector
            // These vectors will be used to determine if the exchange has been completed across all processors
            MPI_Status*  tSendStatus  = ( MPI_Status*  ) alloca( 2 * sizeof( MPI_Status  ) * tNumberOfProcs );
            MPI_Status*  tRecvStatus  = ( MPI_Status*  ) alloca( 2 * sizeof( MPI_Status  ) * tNumberOfProcs );
            MPI_Request* tSendRequest = ( MPI_Request* ) alloca( 2 * sizeof( MPI_Request ) * tNumberOfProcs );
            MPI_Request* tRecvRequest = ( MPI_Request* ) alloca( 2 * sizeof( MPI_Request ) * tNumberOfProcs );

            // ncows and ncols of mats to be sent
            Matrix< DDUMat >tSendRowCols( 2, tNumberOfProcs, 0 );

            // ncows and ncols of mats to be received
            Matrix< DDUMat >tRecvRowCols( 2, tNumberOfProcs, 0 );

            // loop over all procs
            for( moris_id k = 0; k < tNumberOfProcs; ++k )
            {
                // only communicate if proc neighbor exists and is not me
                if ( ( aCommunicationList( k ) < tParSize ) && ( aCommunicationList( k ) != tMyRank ) )
                {
                    tSendRowCols( 0, k ) = aMatsToSend( k ).n_rows();
                    tSendRowCols( 1, k ) = aMatsToSend( k ).n_cols();

                    // make sure that MPI can send this data set
                    MORIS_ASSERT( aMatsToSend( k ).numel() < INT_MAX, "send_mat_to_proc: matrix too big" );
                }
            }

            // determine MPI datatype
            MPI_Datatype tRowsColsType = get_comm_datatype ( ( uint ) 0 );

            // loop over all procs
            for( moris_id k=0; k<tNumberOfProcs; ++k )
            {
                // only communicate if proc neighbor exists and is not me
                if ( ( aCommunicationList( k ) < tParSize ) && ( aCommunicationList( k ) != tMyRank ) )
                {
                    // create send data
                    uint tSendArray[ 2 ];
                    tSendArray[ 0 ] = tSendRowCols( 0, k );
                    tSendArray[ 1 ] = tSendRowCols( 1, k );

                    // create tag
                    int tSendTag = create_comm_tag ( tMyRank, aCommunicationList( k ) );

                    // send array size
                    MPI_Isend( &tSendArray,
                            2,
                            tRowsColsType,
                            aCommunicationList( k ),
                            tSendTag,
                            gMorisComm.get_global_comm(),
                            &tSendRequest[ k ] );

                    // create receive data
                    uint tRecvArray[ 2 ];

                    // create tag
                    int tRecvTag = create_comm_tag ( aCommunicationList( k ), tMyRank );

                    // receive array size
                    MPI_Irecv( &tRecvArray,
                            2,
                            tRowsColsType,
                            aCommunicationList( k ),
                            tRecvTag,
                            gMorisComm.get_global_comm(),
                            &tRecvRequest[ k ] );

                    MPI_Wait( &tSendRequest[ k ], &tSendStatus[ k ] );
                    MPI_Wait( &tRecvRequest[ k ], &tRecvStatus[ k ] );

                    // store received data
                    tRecvRowCols( 0, k ) = tRecvArray[ 0 ];
                    tRecvRowCols( 1, k ) = tRecvArray[ 1 ];
                }
            }

            // clear output matrix
            Matrix< T > tEmpty;
            aMatsToReceive.clear();
            aMatsToReceive.resize( tNumberOfProcs, tEmpty );

            // send and receive matrices
            for( moris_id k=0; k<tNumberOfProcs; ++k )
            {
                // only communicate if proc neighbor exists and is not me
                if ( ( aCommunicationList( k ) < tParSize ) && ( aCommunicationList( k ) != tMyRank ) )
                {
                    uint l = k + tNumberOfProcs;

                    // get data type
                    MPI_Datatype tDataType = get_comm_datatype ( ( typename Matrix< T >::Data_Type ) 0 );

                    // calculate length of array to send
                    int tSendLength = tSendRowCols( 0, k )*tSendRowCols( 1, k );

                    if ( tSendLength > 0 )
                    {
                        // create tag
                        int tSendTag = create_comm_tag ( tMyRank, aCommunicationList( k ) ) + 1;

                        // send array
                        MPI_Isend( aMatsToSend( k ).data(),
                                tSendLength,
                                tDataType,
                                aCommunicationList( k ),
                                tSendTag,
                                gMorisComm.get_global_comm(),
                                &tSendRequest[ l ] );
                    }

                    // length of array to receive
                    int tRecvLength = tRecvRowCols( 0, k )*tRecvRowCols( 1, k );

                    if ( tRecvLength > 0 )
                    {
                        // assign memory for matrix
                        aMatsToReceive( k ).set_size( tRecvRowCols( 0, k ), tRecvRowCols( 1, k ) );

                        // create tag
                        int tRecvTag = create_comm_tag ( aCommunicationList( k ), tMyRank ) + 1;

                        // receive array size
                        MPI_Irecv( aMatsToReceive( k ).data(),  //tRecvArray,
                                tRecvLength,
                                tDataType,
                                aCommunicationList( k ),
                                tRecvTag,
                                gMorisComm.get_global_comm(),
                                &tRecvRequest[ l ] );
                    }

                    // wait until both messages are complete
                    if ( tSendLength > 0 )
                    {
                        MPI_Wait( &tSendRequest[ l ], &tSendStatus[ l ] );
                    }

                    if ( tRecvLength > 0 )
                    {
                        MPI_Wait( &tRecvRequest[ l ], &tRecvStatus[ l ] );
                    }
                }
            }
        }
    }

    //------------------------------------------------------------------------------

    // FIXME: the following rountine needs to be simplified using all_gather
    /**
     * This function collects the value aValue from all procs and puts
     * the result in a moris::Mat
     *
     * @param[in] aValue    scalar value from current proc
     *
     * @return  Mat<T>      matrix containing values from all procs
     *
     */
    template <typename T>
    void
    comm_gather_and_broadcast( typename Matrix< T >::Data_Type aValue, Matrix< T > & aMatrix )
    {
        moris_id tParSize = par_size();

        aMatrix.set_size( tParSize, 1 );

        if( tParSize > 1 )
        {
            // get data type
            MPI_Datatype tDataType = get_comm_datatype ( ( typename Matrix< T >::Data_Type) 0 );

            // create send array
            typename Matrix< T >::Data_Type tSendArray[ 1 ] = { aValue };

            // create receive buffer
            typename Matrix< T >::Data_Type* tRecvArray =  new typename Matrix< T >::Data_Type[ tParSize ];

            // gather value over all procs
            MPI_Gather(
                    tSendArray,
                    1,
                    tDataType,
                    tRecvArray,
                    1,
                    tDataType,
                    0,
                    gMorisComm.get_global_comm() );

            // broadcast values over all procs
            MPI_Bcast(
                    tRecvArray,
                    tParSize,
                    tDataType,
                    0,
                    gMorisComm.get_global_comm() );

            // write values of buffer into output matrix
            for( moris_id k=0; k<tParSize; ++k )
            {
                aMatrix( k ) = tRecvArray[ k ];
            }

            // delete array
            delete [] tRecvArray;
        }
        else
        {
            aMatrix( 0 ) = aValue;
        }
    }

    //------------------------------------------------------------------------------

    inline
    std::string parallelize_path( const std::string & aFilePath )
    {
        if( par_size() == 1 || aFilePath.size() == 0 )
        {
            // leave path untouched
            return aFilePath;
        }
        else
        {
            return        aFilePath.substr(0,aFilePath.find_last_of(".")) // base path
                    + "." + std::to_string( par_size() ) // rank of this processor
                    + "." + std::to_string( par_rank() ) // number of procs
                    +  aFilePath.substr( aFilePath.find_last_of("."), aFilePath.length() ); // file extension
        }
    }

    //------------------------------------------------------------------------------
    // print dots for nice output

    inline
    std::string proc_string()
    {
        std::string tString = "              ";

        if( par_size() > 1 )
        {
            uint tMyRank = par_rank();
            tString = "  proc " + std::to_string( tMyRank );

            if ( tMyRank < 10 )
            {
                tString +=" ... :" ;
            }
            else if ( tMyRank < 100 )
            {
                tString +=" .. :" ;
            }
            else if ( tMyRank < 1000 )
            {
                tString +=" . :" ;
            }
            else if ( tMyRank < 10000 )
            {
                tString +="  :" ;
            }
            else
            {
                tString +=" :" ;
            }
        }

        return tString;
    }

    //------------------------------------------------------------------------------

    /*!
     * Gathers matrix from procs on root proc
     * Only aGatheredMats is populated on base proc
     */
    template <typename MatrixType>
    void
    all_gather_vector(
            Matrix<MatrixType> const & aMatToGather,
            Cell<Matrix<MatrixType>> & aGatheredMats,
            moris_index                aTag,
            moris_index                aFixedDim,
            moris_index                aBaseProc = 0)
    {
        // if rows are fixed
        moris_index tNumRow = 0;
        moris_index tNumCol = 0;

        typename Matrix< MatrixType >::Data_Type tTyper = 1;

        if(aFixedDim == 0)
        {
            
            tNumRow = aMatToGather.n_rows();
        }
        else
        {
            tNumCol = aMatToGather.n_cols();
        }
  
        MPI_Request tRequest;

        // send the matrix
        MPI_Isend(
                aMatToGather.data(),
                aMatToGather.numel(),
                moris::get_comm_datatype(tTyper),
                aBaseProc,
                aTag,
                moris::get_comm(),
                &tRequest);

        barrier();

        // on base rank go ahead and receive the data
        if(par_rank() == aBaseProc)
        {
            aGatheredMats.resize(par_size());

            for(int i = 0; i < par_size(); i++)
            {
                MPI_Status tStatus;
                MPI_Probe(i, aTag, moris::get_comm(), &tStatus);

                //    MORIS_ERROR(tExists,"Trying to receive a message that does not exists");
                 
                int tLength = 0;
                MPI_Get_count(
                        &tStatus,
                        moris::get_comm_datatype(tTyper),
                        &tLength);

                if(tLength != 0 )
                {
                    if(aFixedDim == 0)
                    {
                    tNumCol = tLength/tNumRow;
                    }
                    else
                    {
                    tNumRow = tLength/tNumCol;
                    }
                }
                // Resize the matrix
                aGatheredMats(i).resize(tNumRow, tNumCol);

                MPI_Recv(
                        aGatheredMats(i).data(),
                        tLength,
                        moris::get_comm_datatype(tTyper),
                        i,
                        aTag,
                        moris::get_comm(),
                        &tStatus);
            }
        }

        barrier();
    }

    //------------------------------------------------------------------------------

    /*!
     * Gather cell of strings
     */
    void
    all_gather_cell_of_str(
            Cell<std::string> const & aCellToGather,
            Cell<Cell<std::string>> & aGatheredCells,
            moris_index aTag,
            moris_index aBaseProc = 0);

    //------------------------------------------------------------------------------

    /* computes offset for each processor given a local value
     *
     * @param aLocalInput  local value
     *
     * @ return            offset value (same type as input value)
     */
    template< typename T>
    T
    get_processor_offset(const T & aLocalInput)
    {
        // get number of processors
        moris_id tNumProcs = par_size();

        // get local rank
        moris_id tMyPID = par_rank();

        // create matrix of reals of size tNumProcs
        Matrix<DDRMat> tInputMatrix(tNumProcs,1,0.0);

        // fill in local value
        tInputMatrix(tMyPID) = (real)aLocalInput;

        // sum matrices across all procs
        Matrix<DDRMat> tOutputMtrax = sum_all_matrix(tInputMatrix);

        // compute local offset
        T tOffset=0;

        for (moris_id ip=0;ip<tMyPID;++ip)
        {
            tOffset += (T) tOutputMtrax(ip);
        }

        return tOffset;
    }
}

#endif /* SRC_COMM_CL_COMMUNICATION_TOOLS_HPP_ */
