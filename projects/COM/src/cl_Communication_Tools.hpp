/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_Communication_Tools.hpp
 *
 */

#ifndef SRC_COMM_CL_COMMUNICATION_TOOLS_HPP_
#define SRC_COMM_CL_COMMUNICATION_TOOLS_HPP_

// Standard headers
#include <iostream>
#include <string>

#include "cl_Communication_Enums.hpp"      // COM/src
#include "cl_Communication_Manager.hpp"    // COM/src
#include "cl_Matrix.hpp"
#include "linalg_typedefs.hpp"
#include "fn_isvector.hpp"

// #include "cl_Bitset.hpp" // CON/src
#include "cl_Vector.hpp"    // CON/src

// Mesh specific headers
// #include "cl_MTK_Enums.hpp" // MTK/src

// Externally Defined Global Communicator
extern moris::Comm_Manager gMorisComm;

namespace moris
{

    // value to be used if a proc has no neighbor at the
    // specified position
    const moris_id gNoProcID = MORIS_SINT_MAX;

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
    void comm_split( int color, int key, const std::string& aCommName );

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
    void barrier( std::string aBarrierName = " " );

    /**
     * Broadcast a Message to all Procs
     *
     * @param[in] aMessage      A uint
     *
     * @param[out] aMessage     The message of Proc 0
     *
     */
    void broadcast( uint& aMessage );

    /**
     * Broadcast a Message to all Procs
     *
     * @param[in] aMessage      A string to broadcast
     *
     * @param[out] aMessage     The message of Proc 0
     *
     */
    void broadcast( std::string& aMessage );

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
    scatter( std::vector< moris::size_t >& aMessage );

    moris::uint gather_value_and_bcast_max( moris::uint aMessage );

    /*
     * Prints communication header to log
     * @param[in]  CommType - Enum for communication type to be documented
     */
    void begin_log_communication( enum CommunicationType CommType );

    /*
     * Prints communication header to log
     * @param[in]  CommType - Enum for communication type to be documented
     */
    void end_log_communication( enum CommunicationType CommType );

    /*
     * setup log sink should only be done once by proc 0
     * @param[in]  CommType - Enum for communication type to be documented
     */
    void setup_log_sink( std::string aCommLogName );

    /*
     * Generates filename which is used as a filter for the specific log sink
     *
     * The filename is morisroot/log/comm_log/EnumInt_Count.log
     */
    std::string
    generate_filename( enum CommunicationType CommType );

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
            const uint&       aDecompMethod,
            const uint&       aNumberOfDimensions,
            Matrix< DDUMat >& aProcDims,
            Matrix< DDUMat >& aProcCoords,
            Matrix< IdMat >&  aProcNeighbors );

    //------------------------------------------------------------------------------

    /**
     *
     * @brief                 returns an MPI enum defining the
     *                        data type that is to be communicated.
     *
     * @param[in]      primitive data type with arbitrary value
     *
     * see also http://mpitutorial.com/tutorials/mpi-send-and-receive/
     */
    template< typename T >
    MPI_Datatype
    get_comm_datatype( const T& );

    // moris::lint (32-bit)

    inline MPI_Datatype
    get_comm_datatype( const int& )
    {
        return MPI_INT;
    }

    // moris::lint (64-bit)

    inline MPI_Datatype
    get_comm_datatype( const long int& )
    {
        return MPI_LONG;
    }

    // moris::uint (32-bit)

    inline MPI_Datatype
    get_comm_datatype( const unsigned int& )
    {
        return MPI_UNSIGNED;
    }

    // moris::uint (64-bit)

    inline MPI_Datatype
    get_comm_datatype( const long unsigned int& )
    {
        return MPI_UNSIGNED_LONG;
    }

    // moris::real (32-bit)

    inline MPI_Datatype
    get_comm_datatype( const double& )
    {
        return MPI_DOUBLE;
    }

    // moris::real (64-bit)

    inline MPI_Datatype
    get_comm_datatype( const long double& )
    {
        return MPI_LONG_DOUBLE;
    }

    // MPI_CXX_DOUBLE_COMPLEX is supported since MPI-3.0
#ifdef MPI_CXX_DOUBLE_COMPLEX
    // moris::cplx (32-bit)

    inline MPI_Datatype
    get_comm_datatype( const std::complex< double >& )
    {
        return MPI_CXX_DOUBLE_COMPLEX;
    }
#endif

    // MPI_CXX_LONG_DOUBLE_COMPLEX is supported since MPI-3.0
#ifdef MPI_CXX_LONG_DOUBLE_COMPLEX
    // moris::cplx (64-bit)

    inline MPI_Datatype
    get_comm_datatype( const std::complex< long double >& )
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
    template< typename Data >
    Data
    sum_all( const Data& aLocalInput )
    {
        Data tGlobalSum;
        MPI_Allreduce( &aLocalInput, &tGlobalSum, 1, get_comm_datatype( aLocalInput ), MPI_SUM, MPI_COMM_WORLD );
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

    template< typename E >
    E
    sum_all_matrix( const E& aLocalMatrix )
    {
        E tGlobalMatrix( aLocalMatrix.n_rows(), aLocalMatrix.n_cols() );

        void* aLocalPtr  = (void*)aLocalMatrix.data();
        void* aGlobalPtr = (void*)tGlobalMatrix.data();

        int tNumElems = aLocalMatrix.numel();

        MPI_Datatype tType = get_comm_datatype( (typename E::Data_Type)0 );

        MPI_Allreduce( aLocalPtr, aGlobalPtr, tNumElems, tType, MPI_SUM, MPI_COMM_WORLD );

        return tGlobalMatrix;
    }

    //------------------------------------------------------------------------------

    /*
     * @brief Determines minimum of all components of a matrix on all processors
     *
     * @param[in] aLocalInput     local matrix
     *
     * returns matrix with global minimums
     */

    template< typename E >
    E
    min_all_matrix( const E& aLocalMatrix )
    {
        E tGlobalMatrix( aLocalMatrix.n_rows(), aLocalMatrix.n_cols() );

        void* aLocalPtr  = (void*)aLocalMatrix.data();
        void* aGlobalPtr = (void*)tGlobalMatrix.data();

        int tNumElems = aLocalMatrix.numel();

        MPI_Datatype tType = get_comm_datatype( (typename E::Data_Type)0 );

        MPI_Allreduce( aLocalPtr, aGlobalPtr, tNumElems, tType, MPI_MIN, MPI_COMM_WORLD );

        return tGlobalMatrix;
    }

    //------------------------------------------------------------------------------

    /*
     * @brief Determines maximal value of all local values
     *
     * @param[in] aLocalInput     local value
     * @param[in] aGlobalMax      reference to maximal global value
     */
    template< typename Data >
    Data
    max_all( const Data& aLocalInput )
    {
        Data aGlobalMax;
        MPI_Allreduce( &aLocalInput, &aGlobalMax, 1, get_comm_datatype( aLocalInput ), MPI_MAX, MPI_COMM_WORLD );
        return aGlobalMax;
    }

    //------------------------------------------------------------------------------

    /*
     * @brief Determines maximal value of all local values
     *
     * @param[in] aLocalInput     local value
     * @param[in] aGlobalMax      reference to maximal global value
     */
    template< typename Data >
    Data
    min_all( const Data& aLocalInput )
    {
        Data aGlobalMin;
        MPI_Allreduce( &aLocalInput, &aGlobalMin, 1, get_comm_datatype( aLocalInput ), MPI_MIN, MPI_COMM_WORLD );
        return aGlobalMin;
    }

    //------------------------------------------------------------------------------

    /*
     * MPI_LAND is a predefined reduction operation that will calculate the logical and of the values reduced.
     * https://www.rookiehpc.com/mpi/docs/mpi_land.php
     */

    bool
    all_land( bool aMyBool );

    //------------------------------------------------------------------------------

    /**
     * @brief logical or statement reduction for all processes
     *
     * @param aMyBool
     * @return true
     * @return false
     */

    bool
    all_lor( bool aMyBool );

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
    create_comm_tag( const int& aSource, const int& aTarget );

    //------------------------------------------------------------------------------

    /**
     * Builds a communication table map that
     * @param aCommunicationTable List of processors to communicate with
     * @return Map going from processor to position in the communication table
     */
    Vector< moris_id > build_communication_table_map( const Matrix< IdMat >& aCommunicationTable );

    //------------------------------------------------------------------------------

    /**
     *
     * @brief sends and receives matrix to and from each processor in communication list
     *
     * @param[in] aCommunicationList List of ranks of communicating processors
     * @param[in] aSend              Matrix of scalars to be sent to each processor
     * @param[in] aReceive           Matrix to scalars be received from each processor
     */
    template< typename T >
    void
    communicate_scalars(
            const Vector< moris_index >& aCommunicationList,
            const Vector< T >&           aScalarsToSend,
            Vector< T >&                 aScalarsToReceive )
    {
        // only call this when we are in parallel mode
        if ( par_size() > 1 )
        {
            // get number of procs to communicate with
            uint tNumberOfProcs = aCommunicationList.size();

            // get my ID
            moris_index tMyRank = par_rank();

            // size matrix for storing receiving values
            aScalarsToReceive.resize( tNumberOfProcs, 0 );

            // Allocate memory for request vector
            MPI_Request* tSendRequest = (MPI_Request*)alloca( sizeof( MPI_Request ) * tNumberOfProcs );
            MPI_Request* tRecvRequest = (MPI_Request*)alloca( sizeof( MPI_Request ) * tNumberOfProcs );

            // determine MPI datatype
            MPI_Datatype tType = get_comm_datatype( static_cast< T >( 0 ) );

            // loop over all procs
            for ( uint k = 0; k < tNumberOfProcs; ++k )
            {
                // only communicate if proc neighbor exists and is not me
                if ( ( aCommunicationList( k ) < MORIS_SINT_MAX ) &&    //
                        ( aCommunicationList( k ) != tMyRank ) )
                {
                    // create tag
                    int tSendTag = create_comm_tag( tMyRank, aCommunicationList( k ) );
                    int tRecvTag = create_comm_tag( aCommunicationList( k ), tMyRank );

                    // non-blocking request to receive values
                    MPI_Irecv(
                            &aScalarsToReceive( k ),
                            1,
                            tType,
                            aCommunicationList( k ),
                            tRecvTag,
                            gMorisComm.get_global_comm(),
                            &tRecvRequest[ k ] );

                    // non-blocking request to send values
                    MPI_Isend(
                            &aScalarsToSend( k ),
                            1,
                            tType,
                            aCommunicationList( k ),
                            tSendTag,
                            gMorisComm.get_global_comm(),
                            &tSendRequest[ k ] );
                }
            }

            // loop over all procs
            for ( uint k = 0; k < tNumberOfProcs; ++k )
            {
                // only communicate if proc neighbor exists and is not me
                if ( ( aCommunicationList( k ) < MORIS_SINT_MAX ) &&    //
                        ( aCommunicationList( k ) != tMyRank ) )
                {
                    // wait until send and receive requests have been processed
                    MPI_Wait( &tSendRequest[ k ], MPI_STATUS_IGNORE );
                    MPI_Wait( &tRecvRequest[ k ], MPI_STATUS_IGNORE );
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
    template< typename T >
    void
    communicate_mats(
            const Matrix< IdMat >&     aCommunicationList,
            const Vector< Matrix< T > >& aMatsToSend,
            Vector< Matrix< T > >&       aMatsToReceive )
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
            MPI_Request* tSendRequest = (MPI_Request*)alloca( 2 * sizeof( MPI_Request ) * tNumberOfProcs );
            MPI_Request* tRecvRequest = (MPI_Request*)alloca( 2 * sizeof( MPI_Request ) * tNumberOfProcs );

            // nrows and ncols of mats to be sent
            Matrix< DDUMat > tSendRowCols( 2 * tNumberOfProcs, 1, MORIS_UINT_MAX );

            // nrows and ncols of mats to be received
            Matrix< DDUMat > tRecvRowCols( 2 * tNumberOfProcs, 1, MORIS_UINT_MAX );

            // loop over all procs
            for ( moris_id k = 0; k < tNumberOfProcs; ++k )
            {
                // FIXME: function should only be called with proper communication list such that following assert can be activated
                // MORIS_ASSERT( aCommunicationList( k ) < tParSize,
                //         "communicate_mats - rank of communicating processor exceeds limit.");

                // FIXME: function should only be called with proper communication list such that following assert can be activated
                // MORIS_ASSERT( aCommunicationList( k ) != tMyRank,
                //        "communicate_mats - processor wants to send matrices to itself.");

                // only communicate if proc neighbor exists and is not me
                if ( ( aCommunicationList( k ) < tParSize ) && ( aCommunicationList( k ) != tMyRank ) )
                {
                    tSendRowCols( k * 2 )     = aMatsToSend( k ).n_rows();
                    tSendRowCols( k * 2 + 1 ) = aMatsToSend( k ).n_cols();

                    // make sure that MPI can send this data set
                    MORIS_ASSERT( aMatsToSend( k ).numel() < INT_MAX,
                            "communicate_mats - matrix to be sent is too big." );
                }
            }

            // determine MPI datatype
            MPI_Datatype tRowsColsType = get_comm_datatype( (uint)0 );

            // loop over all procs
            for ( moris_id k = 0; k < tNumberOfProcs; ++k )
            {
                // only communicate if proc neighbor exists and is not me
                if ( ( aCommunicationList( k ) < tParSize ) && ( aCommunicationList( k ) != tMyRank ) )
                {
                    // create tag
                    int tSendTag = create_comm_tag( tMyRank, aCommunicationList( k ) );
                    int tRecvTag = create_comm_tag( aCommunicationList( k ), tMyRank );

                    // receive array size
                    MPI_Irecv(
                            tRecvRowCols.data() + k * 2,
                            2,
                            tRowsColsType,
                            aCommunicationList( k ),
                            tRecvTag,
                            gMorisComm.get_global_comm(),
                            &tRecvRequest[ k ] );

                    // send array size
                    MPI_Isend(
                            tSendRowCols.data() + k * 2,
                            2,
                            tRowsColsType,
                            aCommunicationList( k ),
                            tSendTag,
                            gMorisComm.get_global_comm(),
                            &tSendRequest[ k ] );
                }
            }

            // wait for all send and receive requests have been fulfilled
            for ( moris_id k = 0; k < tNumberOfProcs; ++k )
            {
                // only communicate if proc neighbor exists and is not me
                if ( ( aCommunicationList( k ) < tParSize ) && ( aCommunicationList( k ) != tMyRank ) )
                {
                    MPI_Wait( &tSendRequest[ k ], MPI_STATUS_IGNORE );
                    MPI_Wait( &tRecvRequest[ k ], MPI_STATUS_IGNORE );
                }
            }

            // clear output matrix
            Matrix< T > tEmpty;
            aMatsToReceive.clear();
            aMatsToReceive.resize( tNumberOfProcs, tEmpty );

            // get data type
            MPI_Datatype tDataType = get_comm_datatype( (typename Matrix< T >::Data_Type)0 );

            // send and receive matrices
            for ( moris_id k = 0; k < tNumberOfProcs; ++k )
            {
                // only communicate if proc neighbor exists and is not me
                if ( ( aCommunicationList( k ) < tParSize ) && ( aCommunicationList( k ) != tMyRank ) )
                {
                    // increment request counter
                    uint kk = k + tNumberOfProcs;

                    // calculate length of array to send
                    int tSendLength = tSendRowCols( k * 2 ) * tSendRowCols( k * 2 + 1 );

                    if ( tSendLength > 0 )
                    {
                        // create tag
                        int tSendTag = create_comm_tag( tMyRank, aCommunicationList( k ) ) + 1;

                        // send array
                        MPI_Isend( aMatsToSend( k ).data(),
                                tSendLength,
                                tDataType,
                                aCommunicationList( k ),
                                tSendTag,
                                gMorisComm.get_global_comm(),
                                &tSendRequest[ kk ] );
                    }

                    // length of array to receive
                    int tRecvLength = tRecvRowCols( k * 2 ) * tRecvRowCols( k * 2 + 1 );

                    if ( tRecvLength > 0 )
                    {
                        // assign memory for matrix
                        aMatsToReceive( k ).set_size( tRecvRowCols( k * 2 ), tRecvRowCols( k * 2 + 1 ) );

                        // create tag
                        int tRecvTag = create_comm_tag( aCommunicationList( k ), tMyRank ) + 1;

                        // receive array size
                        MPI_Irecv( aMatsToReceive( k ).data(),    // tRecvArray,
                                tRecvLength,
                                tDataType,
                                aCommunicationList( k ),
                                tRecvTag,
                                gMorisComm.get_global_comm(),
                                &tRecvRequest[ kk ] );
                    }
                }
            }

            // wait for all send and receive requests have been fulfilled
            for ( moris_id k = 0; k < tNumberOfProcs; ++k )
            {
                // only communicate if proc neighbor exists and is not me
                if ( ( aCommunicationList( k ) < tParSize ) && ( aCommunicationList( k ) != tMyRank ) )
                {
                    // increment request counter
                    uint kk = k + tNumberOfProcs;

                    // calculate length of array to send
                    int tSendLength = tSendRowCols( k * 2 ) * tSendRowCols( k * 2 + 1 );

                    // wait until both messages are complete
                    if ( tSendLength > 0 )
                    {
                        MPI_Wait( &tSendRequest[ kk ], MPI_STATUS_IGNORE );
                    }

                    // length of array to receive
                    int tRecvLength = tRecvRowCols( k * 2 ) * tRecvRowCols( k * 2 + 1 );

                    if ( tRecvLength > 0 )
                    {
                        MPI_Wait( &tRecvRequest[ kk ], MPI_STATUS_IGNORE );
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
     * @param[in] aValue    scalar value from current processor
     *
     * @return  Mat<T>      matrix containing values from all processors
     *
     */
    template< typename T >
    void
    comm_gather_and_broadcast(
            typename Matrix< T >::Data_Type aValue,
            Matrix< T >&                    aMatrix )
    {
        // get number of processors
        moris_id tParSize = par_size();

        // set size of matrix to be populated
        aMatrix.set_size( tParSize, 1 );

        // communicate value
        if ( tParSize > 1 )
        {
            // get data type
            MPI_Datatype tDataType = get_comm_datatype( (typename Matrix< T >::Data_Type)0 );

            // gather value over all procs
            MPI_Allgather(
                    &aValue,
                    1,
                    tDataType,
                    aMatrix.data(),
                    1,
                    tDataType,
                    gMorisComm.get_global_comm() );
        }
        else
        {
            aMatrix( 0 ) = aValue;
        }
    }

    //------------------------------------------------------------------------------

    inline std::string
    parallelize_path( const std::string& aFilePath )
    {
        if ( par_size() == 1 || aFilePath.size() == 0 )
        {
            // leave path untouched
            return aFilePath;
        }
        else
        {
            return aFilePath.substr( 0, aFilePath.find_last_of( "." ) )                      // base path
                 + "." + std::to_string( par_size() )                                        // rank of this processor
                 + "." + std::to_string( par_rank() )                                        // number of procs
                 + aFilePath.substr( aFilePath.find_last_of( "." ), aFilePath.length() );    // file extension
        }
    }

    //------------------------------------------------------------------------------
    // print dots for nice output

    inline std::string
    proc_string()
    {
        std::string tString = "              ";

        if ( par_size() > 1 )
        {
            uint tMyRank = par_rank();
            tString      = "  proc " + std::to_string( tMyRank );

            if ( tMyRank < 10 )
            {
                tString += " ... :";
            }
            else if ( tMyRank < 100 )
            {
                tString += " .. :";
            }
            else if ( tMyRank < 1000 )
            {
                tString += " . :";
            }
            else if ( tMyRank < 10000 )
            {
                tString += "  :";
            }
            else
            {
                tString += " :";
            }
        }

        return tString;
    }

    //------------------------------------------------------------------------------

    /**
     * @brief Gathers matrix from all processors on base processor
     * Only aGatheredMats is populated on base processor
     *
     * @param aMatToGather
     * @param aGatheredMats
     * @param aTag
     * @param aFixedDim
     * @param aBaseProc
     * @return * template< typename MatrixType >
     */
    template< typename MatrixType >
    void
    all_gather_vector(
            Matrix< MatrixType > const &  aMatToGather,
            Vector< Matrix< MatrixType > >& aGatheredMats,
            moris_index                   aTag,
            moris_index                   aFixedDim,
            moris_index                   aBaseProc = 0 )
    {
        // if rows are fixed
        moris_index tNumRow = 0;
        moris_index tNumCol = 0;

        typename Matrix< MatrixType >::Data_Type tType = 1;

        if ( aFixedDim == 0 )
        {

            tNumRow = aMatToGather.n_rows();
        }
        else
        {
            tNumCol = aMatToGather.n_cols();
        }

        MPI_Request tRequest;

        // send the matrix via nonblocking request
        MPI_Isend(
                aMatToGather.data(),
                aMatToGather.numel(),
                moris::get_comm_datatype( tType ),
                aBaseProc,
                aTag,
                moris::get_comm(),
                &tRequest );

        // FIXME: should not be needed as base processor issues blocking probe and receive requests
        barrier();

        // on base rank go ahead and receive the data
        if ( par_rank() == aBaseProc )
        {
            aGatheredMats.resize( par_size() );

            for ( int i = 0; i < par_size(); i++ )
            {
                // check and wait until message from processor "i" is ready
                MPI_Status tStatus;
                MPI_Probe( i, aTag, moris::get_comm(), &tStatus );

                // get length of message
                int tLength = 0;
                MPI_Get_count(
                        &tStatus,
                        moris::get_comm_datatype( tType ),
                        &tLength );

                // compute number of rows or columns
                if ( tLength != 0 )
                {
                    if ( aFixedDim == 0 )
                    {
                        tNumCol = tLength / tNumRow;
                    }
                    else
                    {
                        tNumRow = tLength / tNumCol;
                    }
                }

                // check for correct rows or columns as incorrect rounding may happen
                MORIS_ASSERT( tNumCol * tNumRow == tLength,
                        "all_gather_vector - incorrect number of rows or columns." );

                // Resize the matrix
                aGatheredMats( i ).resize( tNumRow, tNumCol );

                // receive matrix
                MPI_Recv(
                        aGatheredMats( i ).data(),
                        tLength,
                        moris::get_comm_datatype( tType ),
                        i,
                        aTag,
                        moris::get_comm(),
                        &tStatus );
            }
        }

        // wait until send request has been processed
        MPI_Wait( &tRequest, MPI_STATUS_IGNORE );

        // FIXME: should not be needed as this point can only be reached if all
        //        send and receive request have been processed
        barrier();
    }

    //------------------------------------------------------------------------------

    /*!
     * Gather cell of strings
     */
    void
    all_gather_cell_of_str(
            Vector< std::string > const &  aCellToGather,
            Vector< Vector< std::string > >& aGatheredCells,
            moris_index                  aTag,
            moris_index                  aBaseProc = 0 );

    //------------------------------------------------------------------------------

    /* computes offset for each processor given a local value
     *
     * @param aLocalInput  local value
     *
     * @ return            offset value (same type as input value)
     */
    template< typename T >
    T
    get_processor_offset( const T& aLocalInput )
    {
        // get local rank
        moris_id tMyPID = par_rank();

        Matrix< DDRMat > tOutputMatrix;
        comm_gather_and_broadcast( aLocalInput, tOutputMatrix );

        // compute local offset
        T tOffset = 0;

        for ( moris_id ip = 0; ip < tMyPID; ++ip )
        {
            tOffset += (T)tOutputMatrix( ip );
        }

        return tOffset;
    }

    //------------------------------------------------------------------------------

    /**
     * @brief print the matrix in a sequential format
     *
     * @tparam T data type
     * @param aMat
     * @param aString
     * @param aPrintFunc
     */

    // How to use this function: because of the issue in the overload resolution template instantiation fails, thus type need to be specified explicitly
    // void (*func)(Matrix<T>, std::string) = print/ print_as_row_vector / print_fancy ; print_mat<T>(aMat, aString,func );
    // for example: void (*func)(Matrix<DDRMat>, std::string) = print_as_row_vector ; print_mat<DDRMat>(aMat, aString,func );
    template< typename T >
    void
    print_mat( const Matrix< T >& aMat, std::string aString, std::function< void( Matrix< T >, std::string ) > aPrintFunc = print )
    {
        // define tag
        int tRowCommTag  = 10;
        int tDataCommTag = 110;

        // get the processr size and rank
        moris_id tParSize = par_size();
        moris_id tMyRank  = par_rank();

        // get data type
        MPI_Datatype tDataType = get_comm_datatype( (typename Matrix< T >::Data_Type)0 );

        // send everything to processor zero
        if ( tMyRank == 0 )
        {
            // create the
            moris::Vector< Matrix< T > > tMats( tParSize );
            moris::Vector< moris::uint > tNumRows( tParSize );

            tNumRows( tMyRank ) = aMat.n_rows();
            tMats( tMyRank )    = aMat;

            aPrintFunc( tMats( tMyRank ), aString + "_" + std::to_string( tMyRank ) );

            // receive the info form all other processors
            for ( int iProc = 1; iProc < tParSize; iProc++ )
            {
                // get number of rows
                MPI_Recv( tNumRows.memptr() + iProc, 1, tDataType, iProc, tRowCommTag, moris::get_comm(), MPI_STATUS_IGNORE );

                // check and wait until message from processor "iProc" is ready
                MPI_Status tStatus;
                MPI_Probe( iProc, tDataCommTag, moris::get_comm(), &tStatus );

                // get length of message to determine columns
                int tNumEl = 0;
                MPI_Get_count(
                        &tStatus,
                        tDataType,
                        &tNumEl );

                // determine the columns
                int tCol = tNumEl / tNumRows( iProc );

                // set the size
                tMats( iProc ).set_size( tNumRows( iProc ), tCol );

                // Receive the next value from the next processor
                MPI_Recv( tMats( iProc ).data(), tMats( iProc ).numel(), tDataType, iProc, tDataCommTag, moris::get_comm(), &tStatus );

                // call the print function
                aPrintFunc( tMats( iProc ), aString + "_" + std::to_string( iProc ) );
            }
        }
        else
        {
            size_t tRows = aMat.n_rows();
            // send number of rows
            MPI_Send( &tRows, 1, tDataType, 0, tRowCommTag, moris::get_comm() );

            // send data
            MPI_Send( aMat.data(), (int)aMat.numel(), tDataType, 0, tDataCommTag, moris::get_comm() );
        }
    }

    //------------------------------------------------------------------------------

    /**
     * @brief  the print cell version of the above function
     *
     * @tparam T
     * @param aCell
     * @param aString
     * @param aPrintFunc
     */


    template< typename T >
    void
    print_cell( const Vector< T >& aCell, std::string aString, std::function< void( Vector< T > const &, std::string ) > aPrintFunc = print )
    {
        // define tag
        int tRowCommTag  = 10;
        int tDataCommTag = 110;

        // get the processr size and rank
        moris_id tParSize = par_size();
        moris_id tMyRank  = par_rank();

        // get data type
        MPI_Datatype tDataType = get_comm_datatype( (typename Vector< T >::value_type)0 );

        // send everything to processor zero
        if ( tMyRank == 0 )
        {
            // create the
            moris::Vector< moris::Vector< T > > tCells( tParSize );
            moris::Vector< moris::uint >      tCellSizes( tParSize );

            tCellSizes( tMyRank ) = aCell.size();
            tCells( tMyRank )     = aCell;

            aPrintFunc( tCells( tMyRank ), aString + "_" + std::to_string( tMyRank ) );

            // recive the info form all other processors
            for ( int iProc = 1; iProc < tParSize; iProc++ )
            {
                // get number of rows
                MPI_Recv( tCellSizes.memptr() + iProc, 1, tDataType, iProc, tRowCommTag, moris::get_comm(), MPI_STATUS_IGNORE );

                tCells( iProc ).resize( tCellSizes( iProc ) );

                // Receive the next value from the next processor
                MPI_Recv( tCells( iProc ).memptr(), tCells( iProc ).size(), tDataType, iProc, tDataCommTag, moris::get_comm(), MPI_STATUS_IGNORE );

                aPrintFunc( tCells( iProc ), aString + "_" + std::to_string( iProc ) );
            }
        }
        else
        {
            size_t tSize = aCell.size();
            // send number of rows
            MPI_Send( &tSize, 1, tDataType, 0, tRowCommTag, moris::get_comm() );

            // send data
            MPI_Send( aCell.memptr(), (int)aCell.size(), tDataType, 0, tDataCommTag, moris::get_comm() );
        }
    }


    //------------------------------------------------------------------------------

    /**
     *
     * @brief                 Communicate cell version of the  matrix version above
     *
     * @param[in] aCell     matrix to be communicated
     * @param[in] aSource     rank of sending proc
     *
     */

    // There is minimal change from the matrix version to diagnose and debug since matrix version is
    // thoroughly tested
    template< typename T >
    void
    communicate_cells(
            const Vector< moris_index >& aCommunicationList,
            const Vector< Vector< T > >&   aCellsToSend,
            Vector< Vector< T > >&         aCellsToReceive )
    {
        moris_id tParSize = par_size();
        // only call this when we are in parallel mode
        if ( tParSize > 1 )
        {
            // get number of procs to communicate with
            moris_id tNumberOfProcs = aCommunicationList.size();

            // get my ID
            moris_id tMyRank = par_rank();

            // Allocate memory for status/request vector
            // These vectors will be used to determine if the exchange has been completed across all processors
            MPI_Request* tSendRequest = (MPI_Request*)alloca( 2 * sizeof( MPI_Request ) * tNumberOfProcs );
            MPI_Request* tRecvRequest = (MPI_Request*)alloca( 2 * sizeof( MPI_Request ) * tNumberOfProcs );

            // nrows and ncols of mats to be sent
            Matrix< DDUMat > tSendRowCols( tNumberOfProcs, 1, MORIS_UINT_MAX );

            // nrows and ncols of mats to be received
            Matrix< DDUMat > tRecvRowCols( tNumberOfProcs, 1, MORIS_UINT_MAX );

            // loop over all procs
            for ( moris_id k = 0; k < tNumberOfProcs; ++k )
            {
                // FIXME: function should only be called with proper communication list such that following assert can be activated
                // MORIS_ASSERT( aCommunicationList( k ) < tParSize,
                //         "communicate_mats - rank of communicating processor exceeds limit.");

                // FIXME: function should only be called with proper communication list such that following assert can be activated
                // MORIS_ASSERT( aCommunicationList( k ) != tMyRank,
                //        "communicate_mats - processor wants to send matrices to itself.");

                // only communicate if proc neighbor exists and is not me
                if ( ( aCommunicationList( k ) < tParSize ) && ( aCommunicationList( k ) != tMyRank ) )
                {
                    tSendRowCols( k ) = aCellsToSend( k ).size();

                    // make sure that MPI can send this data set
                    MORIS_ASSERT( aCellsToSend( k ).size() < INT_MAX,
                            "communicate_mats - matrix to be sent is too big." );
                }
            }

            // determine MPI datatype
            MPI_Datatype tRowsColsType = get_comm_datatype( (uint)0 );

            // loop over all procs
            for ( moris_id k = 0; k < tNumberOfProcs; ++k )
            {
                // only communicate if proc neighbor exists and is not me
                if ( ( aCommunicationList( k ) < tParSize ) && ( aCommunicationList( k ) != tMyRank ) )
                {
                    // create tag
                    int tSendTag = create_comm_tag( tMyRank, aCommunicationList( k ) );
                    int tRecvTag = create_comm_tag( aCommunicationList( k ), tMyRank );

                    // receive array size
                    MPI_Irecv(
                            tRecvRowCols.data() + k,
                            1,
                            tRowsColsType,
                            aCommunicationList( k ),
                            tRecvTag,
                            gMorisComm.get_global_comm(),
                            &tRecvRequest[ k ] );

                    // send array size
                    MPI_Isend(
                            tSendRowCols.data() + k,
                            1,
                            tRowsColsType,
                            aCommunicationList( k ),
                            tSendTag,
                            gMorisComm.get_global_comm(),
                            &tSendRequest[ k ] );
                }
            }

            // wait for all send and receive requests have been fulfilled
            for ( moris_id k = 0; k < tNumberOfProcs; ++k )
            {
                // only communicate if proc neighbor exists and is not me
                if ( ( aCommunicationList( k ) < tParSize ) && ( aCommunicationList( k ) != tMyRank ) )
                {
                    MPI_Wait( &tSendRequest[ k ], MPI_STATUS_IGNORE );
                    MPI_Wait( &tRecvRequest[ k ], MPI_STATUS_IGNORE );
                }
            }

            // clear output matrix
            Vector< T > tEmpty;
            aCellsToReceive.clear();
            aCellsToReceive.resize( tNumberOfProcs, tEmpty );

            // get data type
            MPI_Datatype tDataType = get_comm_datatype( (typename Vector< T >::value_type)0 );

            // send and receive matrices
            for ( moris_id k = 0; k < tNumberOfProcs; ++k )
            {
                // only communicate if proc neighbor exists and is not me
                if ( ( aCommunicationList( k ) < tParSize ) && ( aCommunicationList( k ) != tMyRank ) )
                {
                    // increment request counter
                    uint kk = k + tNumberOfProcs;

                    // calculate length of array to send
                    int tSendLength = tSendRowCols( k );

                    if ( tSendLength > 0 )
                    {
                        // create tag
                        int tSendTag = create_comm_tag( tMyRank, aCommunicationList( k ) ) + 1;

                        // send array
                        MPI_Isend( aCellsToSend( k ).memptr(),
                                tSendLength,
                                tDataType,
                                aCommunicationList( k ),
                                tSendTag,
                                gMorisComm.get_global_comm(),
                                &tSendRequest[ kk ] );
                    }

                    // length of array to receive
                    int tRecvLength = tRecvRowCols( k );

                    if ( tRecvLength > 0 )
                    {
                        // assign memory for matrix
                        aCellsToReceive( k ).resize( tRecvRowCols( k ) );

                        // create tag
                        int tRecvTag = create_comm_tag( aCommunicationList( k ), tMyRank ) + 1;

                        // receive array size
                        MPI_Irecv( aCellsToReceive( k ).memptr(),    // tRecvArray,
                                tRecvLength,
                                tDataType,
                                aCommunicationList( k ),
                                tRecvTag,
                                gMorisComm.get_global_comm(),
                                &tRecvRequest[ kk ] );
                    }
                }
            }

            // wait for all send and receive requests have been fulfilled
            for ( moris_id k = 0; k < tNumberOfProcs; ++k )
            {
                // only communicate if proc neighbor exists and is not me
                if ( ( aCommunicationList( k ) < tParSize ) && ( aCommunicationList( k ) != tMyRank ) )
                {
                    // increment request counter
                    uint kk = k + tNumberOfProcs;

                    // calculate length of array to send
                    int tSendLength = tSendRowCols( k );

                    // wait until both messages are complete
                    if ( tSendLength > 0 )
                    {
                        MPI_Wait( &tSendRequest[ kk ], MPI_STATUS_IGNORE );
                    }

                    // length of array to receive
                    int tRecvLength = tRecvRowCols( k );

                    if ( tRecvLength > 0 )
                    {
                        MPI_Wait( &tRecvRequest[ kk ], MPI_STATUS_IGNORE );
                    }
                }
            }
        }
    }
}    // namespace moris

#endif /* SRC_COMM_CL_COMMUNICATION_TOOLS_HPP_ */
