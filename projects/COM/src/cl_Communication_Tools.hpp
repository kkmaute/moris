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
#include "cl_Matrix.hpp" // LNA/src
#include "linalg_typedefs.hpp"


//#include "fn_mem_pointer.hpp" // LNA/src
//#include "cl_Bitset.hpp" // CON/src
#include "cl_Cell.hpp" // CON/src

// Mesh specific headers
//#include "cl_Mesh_Enums.hpp" // MTK/src

// Externally Defined Global Communicator
extern moris::Comm_Manager gMorisComm;

namespace moris
{
    //------------------------------------------------
    //  GENERAL MPI FUNCTIONS


    /*
     * Returns the current processor rank
     */
    moris::moris_id
    par_rank();

    /*
     * Returns the size of the processor pool
     */
    moris::moris_id
    par_size();

    /*
     * Holds all processors
     */
    void barrier();


    /*
     * Wrapper around MPI_Allreduce to sum up integers
     */
    void Sum_All_Local_Int(
            const moris::uint & aLocalInput,
            moris::uint       & aGlobalSum);


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

    /*
     * Communicate Id of node created on an edge shared by a face
     * @param[in] tNID - ID of node created on a face
     * @param[in]
     */
    void
    communicate_info(Matrix< IdMat >       & aSendProcs,
            Matrix< IdMat >                & aRecvProcs,
            moris::Cell<moris::uint>       & aSendTags,
            moris::Cell<moris::uint>       & aRecvTags,
            moris::Cell<Matrix< DDUMat >>  & aSendMessage,
            moris::Cell<Matrix< DDUMat >>  & aRecvMessage);



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
     * @param[in]  aNumberOfDimensions   must be 1, 2 or 3
     * @param[out] aProcDims             dimensions of generated cart
     * @param[out] aProcCoords           coordinate of current proc
     * @param[out] aProcNeighbors        neighbors of current proc,
     *                                   contains UINT_MAX if no neighbor exists
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
    void
    create_proc_cart(
            const uint              & aNumberOfDimensions,
            Matrix < DDUMat >       & aProcDims,
            Matrix < DDUMat >       & aProcCoords,
            Matrix < IdMat >       & aProcNeighbors );

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
    get_comm_datatype( const T & aSample )
    {
        MORIS_ASSERT( false , "get_comm_datatype: unknown data type.");
        return MPI_DATATYPE_NULL;
    }
    
    // moris::lint (32-bit)
    template <> MPI_Datatype
    get_comm_datatype( const int & aSample )
    {
        return MPI_INT;
    }

    // moris::lint (64-bit)
    template <> MPI_Datatype
    get_comm_datatype( const long int & aSample )
    {
        return MPI_LONG;
    }

    // moris::uint (32-bit)
    template <> MPI_Datatype
    get_comm_datatype( const unsigned int & aSample )
    {
        return MPI_UNSIGNED;
    }

    // moris::uint (64-bit)
    template <> MPI_Datatype
    get_comm_datatype( const long unsigned int & aSample )
    {
        return MPI_UNSIGNED_LONG;
    }

    // moris::real (32-bit)
    template <> MPI_Datatype
    get_comm_datatype( const double & aSample )
    {
        return MPI_DOUBLE;
    }

    // moris::real (64-bit)
    template <> MPI_Datatype
    get_comm_datatype( const long double & aSample )
    {
        return MPI_LONG_DOUBLE;
    }
    
// MPI_CXX_DOUBLE_COMPLEX is supported since MPI-3.0
#ifdef MPI_CXX_DOUBLE_COMPLEX
    // moris::cplx (32-bit)
    template <> MPI_Datatype
    get_comm_datatype( const std::complex<double> & aSample )
    {
        return MPI_CXX_DOUBLE_COMPLEX;
    }
#endif

// MPI_CXX_LONG_DOUBLE_COMPLEX is supported since MPI-3.0
#ifdef MPI_CXX_LONG_DOUBLE_COMPLEX
    // moris::cplx (64-bit)
    template <> MPI_Datatype
    get_comm_datatype( const std::complex<long double> & aSample )
    {
        return MPI_CXX_LONG_DOUBLE_COMPLEX;
    }
#endif

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
     * @brief                 Sends a moris::Mat<T> to the proc aTarget.
     *                        Receiving proc must call recv_mat_from_proc.
     *                        Does nothing if MORIS is run in serial.
     *
     * @param[in] aMatrix     matrix to be communicated
     * @param[in] aTarget     rank of receiving proc
     *
     */
/*    template <typename T> void
    send_mat_to_proc( const Matrix< T > & aMatrix,
                      const int      & aTarget )
    {
		
        if ( par_size() > 1 )
        {
			
			
            // get dimensions of matrix
            uint tRowsCols[ 2 ] ;
            tRowsCols[ 0 ] = aMatrix.n_rows();
            tRowsCols[ 1 ] = aMatrix.n_cols();

            // determine length of matrix to send
            uint tLength = tRowsCols[ 0 ]*tRowsCols[ 1 ] ;

            // make sure that MPI can send this data set
            MORIS_ASSERT( tLength < INT_MAX,
                          "send_mat_to_proc: matrix too big" );
			
            // determine MPI datatype
            MPI_Datatype tRowsColsType = get_comm_datatype ( tLength );
			
            // create message tag
            int tTag = create_comm_tag ( par_rank(), aTarget );

            // send array size
            MPI_Send(
                    tRowsCols,
                    2,
                    tRowsColsType,
                    aTarget,
                    tTag,
                    gMorisComm.get_global_comm() );
                    
            // create temporary buffer array
            typename Matrix< T >::Data_Type* tArray  = new typename Matrix< T >::Data_Type[ tLength ];

            // counter for flattening
            uint tCount = 0;

            // flatten matrix
            for ( uint j=0; j<tRowsCols[ 1 ] ; ++j )
            {
                for ( uint i=0; i<tRowsCols[ 0 ] ; ++i )
                {
                    tArray[ tCount++ ] = aMatrix( i, j );
                }
            }

            // get data type
            typename Matrix< T >::Data_Type tSample = 0;
            MPI_Datatype tType = get_comm_datatype ( tSample );

            // increment tag
            ++tTag;
			
            // send data
            MPI_Send(
                    tArray,
                    tLength,
                    tType,
                    aTarget,
                    tTag,
                    gMorisComm.get_global_comm() );
			
			
            // delete buffer
            delete [] tArray;

        }
    } */
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
    /*template <typename T> void
    recv_mat_from_proc( Matrix< T >  & aMatrix,
                        const int & aSource )
        {
            if ( par_size() > 1 )
            {
                // status
                MPI_Status tStatus;

                // get dimensions of matrix
                uint tRowsCols[ 2 ] ;

                uint tLength = 0;

                // determine MPI data type
                MPI_Datatype tRowsColsType = get_comm_datatype ( tLength );

                // create message tag
                int tTag = create_comm_tag ( aSource, par_rank() );

                // receive array size
                MPI_Recv(
                        tRowsCols,
                        2,
                        tRowsColsType,
                        aSource,
                        tTag,
                        gMorisComm.get_global_comm(),
                        &tStatus );

                // calculate length of array to expect
                tLength = tRowsCols[ 0 ]*tRowsCols[ 1 ];

                // create temporary buffer array
                typename Matrix< T >::Data_Type* tArray  = new typename Matrix< T >::Data_Type[ tLength ];

                // get data type
                typename Matrix< T >::Data_Type tSample = 0;
                MPI_Datatype tType = get_comm_datatype ( tSample );

                // increment tag
                ++tTag;

                // receive data
                MPI_Recv(
                        tArray,
                        tLength,
                        tType,
                        aSource,
                        tTag,
                        gMorisComm.get_global_comm(),
                        &tStatus );

                // assign memory for matrix
                uint tRows = tRowsCols[ 0 ];
                uint tCols = tRowsCols[ 1 ];
                aMatrix.set_size( tRows, tCols );

                // counter for unflattening
                uint tCount = 0;

                // unflatten matrix
                for ( uint j=0; j<tCols ; ++j )
                {
                    for ( uint i=0; i<tRows ; ++i )
                    {
                        aMatrix( i, j ) = tArray[ tCount++ ];
                    }
                }

                // delete buffer
                delete [] tArray;
            }
        } */

//------------------------------------------------------------------------------
    /**
     *
     * @brief                 sends a value to each proc and receives values
     *                        from each proc
     *
     *
     * @param[in] aSend       values to send to each individyal proc
     * @param[in] aReceive    values to receive from each individual proc
     */
    template <typename T> void
    communicate_scalars( const Matrix < IdMat > & aCommunicationList,
                         const Matrix< T >      & aScalarsToSend,
                               Matrix< T >      & aScalarsToReceive )
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
        communicate_mats( const Matrix < IdMat >     & aCommunicationList,
                          const Cell< Matrix< T > >  & aMatsToSend,
                                Cell< Matrix< T > >  & aMatsToReceive )
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
                for( moris_id k=0; k<tNumberOfProcs; ++k )
                {
                    // only communicate if proc neighbor exists and is not me
                    if ( ( aCommunicationList( k ) < tParSize ) && ( aCommunicationList( k ) != tMyRank ) )
                    {
                        tSendRowCols( 0, k ) = aMatsToSend( k ).n_rows();
                        tSendRowCols( 1, k ) = aMatsToSend( k ).n_cols();

                        // make sure that MPI can send this data set
                        MORIS_ASSERT( tSendRowCols( 0, k )*tSendRowCols( 1, k ) < INT_MAX,
                                "send_mat_to_proc: matrix too big" );
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
                        MPI_Isend(
                                &tSendArray,
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
                        MPI_Irecv(
                                &tRecvArray,
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

                        // pointer for sending data
                        typename Matrix< T >::Data_Type* tSendArray = nullptr;

                        // pointer for receiving data
                        typename Matrix< T >::Data_Type* tRecvArray = nullptr;

                        // get data type
                        MPI_Datatype tDataType = get_comm_datatype ( ( typename Matrix< T >::Data_Type ) 0 );

                        // calculate length of array to send
                        int tSendLength = tSendRowCols( 0, k )*tSendRowCols( 1, k );

                        if ( tSendLength > 0 )
                        {
                            // create temporary buffer array for sending
                            tSendArray  = new typename Matrix< T >::Data_Type[ tSendLength ];

                            // counter for flattening
                            uint tCount = 0;

                            // flatten matrix
                            for ( uint j=0; j<tSendRowCols( 1, k ) ; ++j )
                            {
                                for ( uint i=0; i<tSendRowCols( 0, k ) ; ++i )
                                {
                                    tSendArray[ tCount++ ] = aMatsToSend( k )( i, j );
                                }
                            }

                            // create tag
                            int tSendTag = create_comm_tag ( tMyRank, aCommunicationList( k ) ) + 1;

                            // send array
                            MPI_Isend(
                                    tSendArray,
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
                            // create temporary buffer array for receiving
                            tRecvArray  = new typename Matrix< T >::Data_Type[ tRecvLength ];

                            // create tag
                            int tRecvTag = create_comm_tag ( aCommunicationList( k ), tMyRank ) + 1;

                            // receive array size
                            MPI_Irecv(
                                    tRecvArray,
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
                            delete [] tSendArray;
                        }

                        if ( tRecvLength > 0 )
                        {
                            MPI_Wait( &tRecvRequest[ l ], &tRecvStatus[ l ] );

                            // assign memory for matrix
                            aMatsToReceive( k ).set_size( tRecvRowCols( 0, k ), tRecvRowCols( 1, k ) );

                            // counter for unflattening
                            uint tCount = 0;

                            // unflatten matrix
                            for ( uint j=0; j<tRecvRowCols( 1, k ) ; ++j )
                            {
                                for ( uint i=0; i<tRecvRowCols( 0, k ) ; ++i )
                                {
                                    aMatsToReceive( k )( i, j ) = tRecvArray[ tCount++ ];
                                }
                            }
                            delete [] tRecvArray;
                        }
                    }
                }
            }
        }

//------------------------------------------------------------------------------

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
        comm_gather_and_broadcast( typename Matrix< T >::Data_Type aValue, Matrix< T > aMatrix )
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
}


#endif /* SRC_COMM_CL_COMMUNICATION_TOOLS_HPP_ */
