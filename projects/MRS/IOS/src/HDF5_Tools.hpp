/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * HDF5_Tools.hpp
 *
 */

#ifndef PROJECTS_MRS_IOS_SRC_HDF5_TOOLS_HPP_
#define PROJECTS_MRS_IOS_SRC_HDF5_TOOLS_HPP_

#include <cstdio>
#include <string>
#include <fstream>

// HD5 c-interface
#include "hdf5.h"

// communicator
#include "cl_Communication_Manager.hpp"
#include "cl_Communication_Tools.hpp"

// #include "assert.hpp"

#include "moris_typedefs.hpp"     //COR/src
#include "cl_Matrix.hpp"          //LINALG/src
#include "linalg_typedefs.hpp"    //LINALG/src

namespace moris
{
    //------------------------------------------------------------------------------

    /**
     * this function takes a path and makes it parallel
     */
    inline std::string
    make_path_parallel(
            const std::string& aPath,
            bool               aAddParExt = true )
    {
        // test if running in parallel mode
        if ( par_size() > 1 && aAddParExt )
        {
            // get file extension
            auto tFileExt = aPath.substr(
                    aPath.find_last_of( '.' ),
                    aPath.length() );

            // get base path
            auto tBasePath = aPath.substr(
                    0,
                    aPath.find_last_of( '.' ) );

            // add proc number to path
            std::string aParallelPath = tBasePath + "_"
                                      + std::to_string( par_size() ) + "."
                                      + std::to_string( par_rank() ) + tFileExt;
            return aParallelPath;
        }

        // do not modify path
        return aPath;
    }

    //------------------------------------------------------------------------------

    /**
     * create a new HDF5 file
     */
    inline hid_t
    create_hdf5_file(
            const std::string& aPath,
            bool               aAddParExt = true )
    {
        MORIS_ERROR( aPath.size() > 0, "No file path given." );

        // create parallel path
        std::string tPath = make_path_parallel( aPath, aAddParExt );

        // create file
        hid_t tFileID = H5Fcreate(
                tPath.c_str(),
                H5F_ACC_TRUNC,
                H5P_DEFAULT,
                H5P_DEFAULT );

        MORIS_ERROR( tFileID > 0, "create_hdf5_file - could not create hdf5 file %s", tPath.c_str() );

        return tFileID;
    }

    //------------------------------------------------------------------------------

    /**
     * open an existing hdf5 file
     */
    inline hid_t
    open_hdf5_file(
            const std::string& aPath,
            bool               aAddParExt = true,
            bool               aReadOnly  = false )
    {
        MORIS_ERROR( aPath.size() > 0, "No file path given." );

        // check that if file is opened with write permission, i.e. aReadOnly is not true, and
        // it is a parallel process, parallel extension is enabled
        MORIS_ASSERT( par_size() > 1 && !aReadOnly ? aAddParExt : true,
                "open_hdf5_file - same file cannot be opened simultaneously by multiple parallel processes \n"
                "with write permission; use set parallel name extension flag to true for file: %s",
                aPath.c_str() );

        // create parallel path
        std::string tPath = make_path_parallel( aPath, aAddParExt );

        // test if file exists
        std::ifstream tFile( tPath );

        // throw error if file does not exist
        MORIS_ERROR( tFile, "open_hdf5_file - could not open HDF5 file %s", tPath.c_str() );

        // close file
        tFile.close();

        // open file as HDF5 handler - read only access
        if ( aReadOnly )
        {
            return H5Fopen(
                    tPath.c_str(),
                    H5F_ACC_RDONLY,
                    H5P_DEFAULT );
        }

        // open file as HDF5 handler - read / write access
        return H5Fopen(
                tPath.c_str(),
                H5F_ACC_RDWR,
                H5P_DEFAULT );
    }

    //------------------------------------------------------------------------------

    /**
     * close an open hdf5 file
     */
    inline herr_t
    close_hdf5_file( hid_t aFileID )
    {
        return H5Fclose( aFileID );
    }

    //------------------------------------------------------------------------------

    /**
     * test if a data set exists
     */
    inline bool
    dataset_exists( hid_t aFileID, const std::string& aLabel )
    {
        hid_t tDataSet = 0;
        return H5Lexists( aFileID, aLabel.c_str(), tDataSet );
    }

    //------------------------------------------------------------------------------
    /**
     *
     * @brief                 returns a HDF5 enum defining the
     *                        data type that is to be communicated.
     *
     * @param[in] aSample     primitive data type with arbitrary value
     *
     * see also https://support.hdfgroup.org/HDF5/doc/H5.user/Datatypes.html
     */
    template< typename T >
    inline hid_t
    get_hdf5_datatype( const T& aSample )
    {
        MORIS_ERROR( false, "get_hdf5_datatype: unknown data type." );
        return 0;
    }

    //------------------------------------------------------------------------------

    template<>
    inline hid_t
    get_hdf5_datatype( const int& )
    {
        return H5T_NATIVE_INT;
    }

    //------------------------------------------------------------------------------

    template<>
    inline hid_t
    get_hdf5_datatype( const long int& )
    {
        return H5T_NATIVE_LONG;
    }

    //------------------------------------------------------------------------------

    // moris::uint
    template<>
    inline hid_t
    get_hdf5_datatype( const unsigned int& )
    {
        return H5T_NATIVE_UINT;
    }

    //------------------------------------------------------------------------------

    // moris::luint
    template<>
    inline hid_t
    get_hdf5_datatype( const long unsigned int& )
    {
        return H5T_NATIVE_ULONG;
    }

    //------------------------------------------------------------------------------

    template<>
    inline hid_t
    get_hdf5_datatype( const double& )
    {
        return H5T_NATIVE_DOUBLE;
    }

    //------------------------------------------------------------------------------

    template<>
    inline hid_t
    get_hdf5_datatype( const long double& )
    {
        return H5T_NATIVE_LDOUBLE;
    }

    //------------------------------------------------------------------------------

    template<>
    inline hid_t
    get_hdf5_datatype( const bool& )
    {
        return H5T_NATIVE_HBOOL;
    }

    //------------------------------------------------------------------------------

    /**
     * this function returns true of both the HDF5 datatype
     * and the passed datatype have the same size
     */
    template< typename T >
    inline bool
    test_size_of_datatype( const T aValue )
    {
        return H5Tget_size( get_hdf5_datatype( aValue ) )
            == sizeof( T );
    }

    //------------------------------------------------------------------------------

    /**
     * unpacks a moris::Mat and stores it into the file
     * file must be open
     *
     * @param[ inout ] aFileID  handler to hdf5 file
     * @param[ in ]    aLabel   label of matrix to save
     * @param[ in ]    aMatrix  moris mat that is to be stored
     * @param[ in ]    aStatus  error handler
     *
     * see also
     * https://support.hdfgroup.org/ftp/HDF5/examples/misc-examples/
     */
    template< typename T >
    inline void
    save_matrix_to_hdf5_file(
            hid_t&             aFileID,
            const std::string& aLabel,
            const Matrix< T >& aMatrix,
            herr_t&            aStatus )
    {
        // check datatype
        MORIS_ASSERT( test_size_of_datatype( (typename Matrix< T >::Data_Type)0 ),
                "Sizes of MORIS datatype and HDF5 datatype do not match." );

        //        MORIS_ASSERT( aMatrix.numel() != 0, "save_matrix_to_hdf5_file(), Output matrix has zero elements");

        // test if data set exists
        if ( dataset_exists( aFileID, aLabel ) )
        {
            // create message
            std::string tMessage = "The data set " + aLabel + " can not be created because it already exist.";

            // throw error
            // MORIS_ERROR( false, tMessage.c_str() );
        }

        // matrix dimensions
        hsize_t tDims[ 2 ];

        // get dimensions from matrix
        tDims[ 0 ] = aMatrix.n_rows();
        tDims[ 1 ] = aMatrix.n_cols();

        // create data space
        hid_t tDataSpace = H5Screate_simple( 2, tDims, nullptr );

        // select data type for matrix to save
        hid_t tDataType = H5Tcopy( get_hdf5_datatype( (typename Matrix< T >::Data_Type)0 ) );

        // set data type to little endian
        aStatus        = H5Tset_order( tDataType, H5T_ORDER_LE );
        hid_t tDataSet = 0;

        // create new dataset
        tDataSet = H5Dcreate(
                aFileID,
                aLabel.c_str(),
                tDataType,
                tDataSpace,
                H5P_DEFAULT,
                H5P_DEFAULT,
                H5P_DEFAULT );

        // test if matrix is not empty
        if ( tDims[ 0 ] * tDims[ 1 ] > 0 )
        {
            // allocate top level array which contains rows
            typename Matrix< T >::Data_Type** tData = (          //
                    typename Matrix< T >::Data_Type**)malloc(    //
                    tDims[ 0 ] * sizeof( typename Matrix< T >::Data_Type* ) );

            // allocate memory for data
            tData[ 0 ] =
                    (typename Matrix< T >::Data_Type*)malloc(
                            tDims[ 0 ] * tDims[ 1 ] * sizeof( typename Matrix< T >::Data_Type ) );

            // loop over all rows and allocate colums
            for ( hsize_t i = 0; i < tDims[ 0 ]; ++i )
            {
                tData[ i ] = tData[ 0 ] + i * tDims[ 1 ];
            }

            // convert moris::Mat to data
            for ( hsize_t i = 0; i < tDims[ 0 ]; ++i )
            {
                for ( hsize_t j = 0; j < tDims[ 1 ]; ++j )
                {
                    tData[ i ][ j ] = aMatrix( i, j );
                }
            }

            // write data into data set
            aStatus = H5Dwrite(
                    tDataSet,
                    tDataType,
                    H5S_ALL,
                    H5S_ALL,
                    H5P_DEFAULT,
                    &tData[ 0 ][ 0 ] );

            // tidy up memory
            free( tData[ 0 ] );
            free( tData );
        }

        // close open hids
        H5Sclose( tDataSpace );
        H5Tclose( tDataType );
        H5Dclose( tDataSet );

        // check for error
        MORIS_ASSERT( aStatus == 0, "Error in HDF5 save_matrix_to_hdf5_file()" );
    }

    //------------------------------------------------------------------------------

    /**
     * unpacks a moris::Mat and stores it into the file
     * file must be open
     *
     * @param[ inout ] aFileID  handler to hdf5 file
     * @param[ in ]    aLabel   label of matrix to save
     * @param[ out ]    aMatrix  moris mat that is to be loaded
     * @param[ in ]    aStatus  error handler
     *
     * see also
     * https://support.hdfgroup.org/ftp/HDF5/examples/misc-examples/
     */
    template< typename T >
    inline void
    load_matrix_from_hdf5_file(
            hid_t&             aFileID,
            const std::string& aLabel,
            Matrix< T >&       aMatrix,
            herr_t&            aStatus )
    {
        // check datatype
        MORIS_ASSERT( test_size_of_datatype( (typename Matrix< T >::Data_Type)0 ),
                "Sizes of MORIS data type and HDF5 data type do not match." );

        // test if data set exists
        MORIS_ERROR( dataset_exists( aFileID, aLabel ),
                "The data set %s cannot be opened because it does not exist.\n",
                aLabel.c_str() );

        // open the data set
        hid_t tDataSet = H5Dopen1( aFileID, aLabel.c_str() );

        // get the data type of the set
        hid_t tDataType = H5Dget_type( tDataSet );

        // make sure that data type fits to type of matrix
        MORIS_ERROR( H5Tget_class( tDataType ) == H5Tget_class( get_hdf5_datatype( (typename Matrix< T >::Data_Type)0 ) ),
                "ERROR in reading from file: matrix %s has the wrong data type.",
                aLabel.c_str() );

        // get handler to data space
        hid_t tDataSpace = H5Dget_space( tDataSet );

        // matrix dimensions
        hsize_t tDims[ 2 ];

        // ask hdf for dimensions
        aStatus = H5Sget_simple_extent_dims( tDataSpace, tDims, nullptr );

        // allocate memory for output
        aMatrix.set_size( tDims[ 0 ], tDims[ 1 ] );

        // test if matrix is not empty
        if ( tDims[ 0 ] * tDims[ 1 ] > 0 )
        {
            // allocate top level array which contains rows
            typename Matrix< T >::Data_Type** tData =
                    (typename Matrix< T >::Data_Type**)malloc(
                            tDims[ 0 ] * sizeof( typename Matrix< T >::Data_Type* ) );

            // allocate memory for data
            tData[ 0 ] = (typename Matrix< T >::Data_Type*)malloc(
                    tDims[ 0 ] * tDims[ 1 ] * sizeof( typename Matrix< T >::Data_Type ) );

            // loop over all rows and allocate columns
            for ( hsize_t i = 1; i < tDims[ 0 ]; ++i )
            {
                tData[ i ] = tData[ 0 ] + i * tDims[ 1 ];
            }

            // read data from file
            aStatus = H5Dread(
                    tDataSet,
                    tDataType,
                    H5S_ALL,
                    H5S_ALL,
                    H5P_DEFAULT,
                    &tData[ 0 ][ 0 ] );

            // write values into matrix
            for ( hsize_t j = 0; j < tDims[ 1 ]; ++j )
            {
                for ( hsize_t i = 0; i < tDims[ 0 ]; ++i )
                {
                    aMatrix( i, j ) = tData[ i ][ j ];
                }
            }

            // tidy up memory
            free( tData[ 0 ] );
            free( tData );
        }
        else if ( aStatus == 2 )
        {
            // all good, reset status
            aStatus = 0;
        }
        // Close/release resources
        H5Tclose( tDataType );
        H5Dclose( tDataSet );
        H5Sclose( tDataSpace );

        // check for error
        MORIS_ASSERT( aStatus == 0, "Error in HDF5 load_matrix_from_hdf5_file()" );
    }

    //------------------------------------------------------------------------------

    /**
     * unpacks a std::vector and stores it into the file
     * file must be open - used for plato integration do not remove
     *
     * @param[ inout ] aFileID  handler to hdf5 file
     * @param[ in ]    aLabel   label of vector to save
     * @param[ in ]    aVector  vector that is to be stored
     * @param[ in ]    aStatus  error handler
     *
     * see also
     * https://support.hdfgroup.org/ftp/HDF5/examples/misc-examples/
     */
    template< typename T >
    void
    save_vector_to_hdf5_file(
            hid_t&                  aFileID,
            const std::string&      aLabel,
            const std::vector< T >& aVector,
            herr_t&                 aStatus )
    {
        // check datatype
        MORIS_ASSERT( test_size_of_datatype( (T)0 ), "Sizes of MORIS data type and HDF5 data type do not match." );

        // test if dataset exists
        if ( dataset_exists( aFileID, aLabel ) )
        {
            // create message
            std::string tMessage = "The dataset " + aLabel + " can not be created because it already exist.";

            // MORIS_ERROR( false, tMessage.c_str() );
        }

        // get dimensions from vector
        hsize_t tSize = aVector.size();

        // create data space
        hid_t tDataSpace = H5Screate_simple( 1, &tSize, nullptr );

        // select data type for vector to save
        hid_t tDataType = H5Tcopy( get_hdf5_datatype( (T)0 ) );

        // set data type to little endian
        aStatus        = H5Tset_order( tDataType, H5T_ORDER_LE );
        hid_t tDataSet = 0;

        // create new dataset
        tDataSet = H5Dcreate(
                aFileID,
                aLabel.c_str(),
                tDataType,
                tDataSpace,
                H5P_DEFAULT,
                H5P_DEFAULT,
                H5P_DEFAULT );

        // test if vector is not empty
        if ( tSize > 0 )
        {
            // write data into dataset
            aStatus = H5Dwrite(
                    tDataSet,
                    tDataType,
                    H5S_ALL,
                    H5S_ALL,
                    H5P_DEFAULT,
                    aVector.data() );
        }

        // close open hids
        H5Sclose( tDataSpace );
        H5Tclose( tDataType );
        H5Dclose( tDataSet );

        // check for error
        MORIS_ASSERT( aStatus == 0, "Error in HDF5 save_vector_to_hdf5_file()" );
    }

    //------------------------------------------------------------------------------

    /**
     * unpacks a std::vector and stores it into the file
     * file must be open
     *
     * @param[ inout ] aFileID  handler to hdf5 file
     * @param[ in ]    aLabel   label of vector to save
     * @param[ out ]    aVector   vector that is to be loaded
     * @param[ in ]    aStatus  error handler
     *
     * see also
     * https://support.hdfgroup.org/ftp/HDF5/examples/misc-examples/
     */
    template< typename T >
    void
    load_vector_from_hdf5_file(
            hid_t&             aFileID,
            const std::string& aLabel,
            std::vector< T >&  aVector,
            herr_t&            aStatus )
    {
        // check datatype
        MORIS_ASSERT( test_size_of_datatype( (T)0 ), "Sizes of MORIS data type and HDF5 data type do not match." );

        // test if dataset exists
        if ( !dataset_exists( aFileID, aLabel ) )
        {
            std::ostringstream tErrorMessage;
            tErrorMessage << "\n\n************** ERROR IN FILE: " << __FILE__ << ", FUNCTION: " << __PRETTY_FUNCTION__
                          << ", LINE: " << __LINE__
                          << ", MESSAGE: The dataset " + aLabel + " can not be opened because it does not exist.";
            throw std::runtime_error( tErrorMessage.str().c_str() );
        }

        // open the data set
        hid_t tDataSet = H5Dopen1( aFileID, aLabel.c_str() );

        // get the data type of the set
        hid_t tDataType = H5Dget_type( tDataSet );

        // make sure that datatype fits to type of vector
        MORIS_ERROR( H5Tget_class( tDataType ) == H5Tget_class( get_hdf5_datatype( (T)0 ) ),
                "ERROR in reading from file: vector %s has the wrong datatype.\n",
                aLabel.c_str() );

        // get handler to dataspace
        hid_t tDataSpace = H5Dget_space( tDataSet );

        // vector dimensions
        hsize_t tDims[ 1 ];

        // ask hdf for dimensions
        aStatus = H5Sget_simple_extent_dims( tDataSpace, tDims, nullptr );

        // allocate memory for output
        aVector.resize( tDims[ 0 ] );

        // test if vector is not empty
        if ( tDims[ 0 ] > 0 )
        {
            // read data from file
            aStatus = H5Dread(
                    tDataSet,
                    tDataType,
                    H5S_ALL,
                    H5S_ALL,
                    H5P_DEFAULT,
                    aVector.data() );
        }
        else if ( aStatus == 2 )
        {
            // all good, reset status
            aStatus = 0;
        }
        // Close/release resources
        H5Tclose( tDataType );
        H5Dclose( tDataSet );
        H5Sclose( tDataSpace );

        // check for error
        MORIS_ASSERT( aStatus == 0, "Error in HDF5 load_vector_from_hdf5_file()" );
    }

    //------------------------------------------------------------------------------

    /**
     * saves a scalar value to a file
     * file must be open
     *
     * @param[ inout ] aFileID  handler to hdf5 file
     * @param[ in ]    aLabel   label of matrix to save
     * @param[ in ]    aValue   value that is to be stored
     * @param[ in ]    aStatus  error handler
     */
    template< typename T >
    inline void
    save_scalar_to_hdf5_file(
            hid_t&             aFileID,
            const std::string& aLabel,
            const T&           aValue,
            herr_t&            aStatus )
    {
        // test if data set exists
        MORIS_ERROR( !dataset_exists( aFileID, aLabel ),
                "The data set %s can not be created because it already exist.\n",
                aLabel.c_str() );

        // check data type
        MORIS_ASSERT( test_size_of_datatype( (T)0 ),
                "Sizes of MORIS data type and HDF5 data type do not match." );

        // select data type for matrix to save
        hid_t tDataType = H5Tcopy( get_hdf5_datatype( (T)0 ) );

        // set data type to little Indian
        aStatus = H5Tset_order( tDataType, H5T_ORDER_LE );

        // matrix dimensions
        hsize_t tDims[ 1 ] = { 1 };

        // create data space
        hid_t tDataSpace = H5Screate_simple( 1, tDims, nullptr );

        // create new data set
        hid_t tDataSet = H5Dcreate(
                aFileID,
                aLabel.c_str(),
                tDataType,
                tDataSpace,
                H5P_DEFAULT,
                H5P_DEFAULT,
                H5P_DEFAULT );

        // write data into data set
        aStatus = H5Dwrite(
                tDataSet,
                tDataType,
                H5S_ALL,
                H5S_ALL,
                H5P_DEFAULT,
                &aValue );

        // close open hids
        H5Sclose( tDataSpace );
        H5Tclose( tDataType );
        H5Dclose( tDataSet );

        // check for error
        MORIS_ASSERT( aStatus == 0, "Error in HDF5 save_scalar_to_hdf5_file()" );
    }

    //------------------------------------------------------------------------------

    template<>
    inline void
    save_scalar_to_hdf5_file(
            hid_t&             aFileID,
            const std::string& aLabel,
            const bool&        aValue,
            herr_t&            aStatus )
    {
        // test if data set exists
        MORIS_ERROR( !dataset_exists( aFileID, aLabel ),
                "The data set %s can not be created because it already exist.\n",
                aLabel.c_str() );

        // select data type for matrix to save
        hid_t tDataType = H5Tcopy( H5T_NATIVE_HBOOL );

        // matrix dimensions
        hsize_t tDims[ 1 ] = { 1 };

        // create data space
        hid_t tDataSpace = H5Screate_simple( 1, tDims, nullptr );

        // create new data set
        hid_t tDataSet = H5Dcreate(
                aFileID,
                aLabel.c_str(),
                tDataType,
                tDataSpace,
                H5P_DEFAULT,
                H5P_DEFAULT,
                H5P_DEFAULT );

        // value to cast bool to
        hbool_t tValue = (hbool_t)aValue;

        // write data into data set
        aStatus = H5Dwrite(
                tDataSet,
                tDataType,
                H5S_ALL,
                H5S_ALL,
                H5P_DEFAULT,
                &tValue );

        // close open hids
        H5Sclose( tDataSpace );
        H5Tclose( tDataType );
        H5Dclose( tDataSet );

        // check for error
        MORIS_ASSERT( aStatus == 0, "Error in HDF5 save_scalar_to_hdf5_file()" );
    }

    //------------------------------------------------------------------------------

    /**
     * loads a scalar value from a file
     * file must be open
     *
     * @param[ inout ] aFileID  handler to hdf5 file
     * @param[ in ]    aLabel   label of matrix to save
     * @param[ out ]   aValue  moris mat that is to be loaded
     * @param[ in ]    aStatus  error handler
     *
     * see also
     * https://support.hdfgroup.org/ftp/HDF5/examples/misc-examples/
     */
    template< typename T >
    inline void
    load_scalar_from_hdf5_file(
            hid_t&             aFileID,
            const std::string& aLabel,
            T&                 aValue,
            herr_t&            aStatus )
    {
        // test if data set exists
        MORIS_ERROR( dataset_exists( aFileID, aLabel ),
                "The data set %s can not be opened because it does not exist.\n",
                aLabel.c_str() );

        // check datatype
        MORIS_ASSERT( test_size_of_datatype( (T)0 ),
                "Sizes of MORIS data type and HDF5 data type do not match." );

        // open the data set
        hid_t tDataSet = H5Dopen1( aFileID, aLabel.c_str() );

        // get the data type of the set
        hid_t tDataType = H5Dget_type( tDataSet );

        // make sure that datatype i scorrect
        /*if ( H5Tget_class( tDataType ) != H5Tget_class( ( T ) 0 ) )
        {
            std::string tMessage = "ERROR in reading from file: scalar "
                                   + aLabel + " has the wrong datatype.\n";

            MORIS_ERROR( false, tMessage.c_str() );
        } */

        // get handler to dataspace
        hid_t tDataSpace = H5Dget_space( tDataSet );

        // read data from file
        aStatus = H5Dread(
                tDataSet,
                tDataType,
                H5S_ALL,
                H5S_ALL,
                H5P_DEFAULT,
                &aValue );

        // Close/release resources
        H5Tclose( tDataType );
        H5Dclose( tDataSet );
        H5Sclose( tDataSpace );

        // check for error
        MORIS_ASSERT( aStatus == 0, "Error in HDF5 load_scalar_from_hdf5_file()" );
    }

    //------------------------------------------------------------------------------

    template<>
    inline void
    load_scalar_from_hdf5_file(
            hid_t&             aFileID,
            const std::string& aLabel,
            bool&              aValue,
            herr_t&            aStatus )
    {
        // test if data set exists
        MORIS_ERROR( dataset_exists( aFileID, aLabel ),
                "The data set %s can not be opened because it does not exist.\n",
                aLabel.c_str() );

        // open the data set
        hid_t tDataSet = H5Dopen1( aFileID, aLabel.c_str() );

        // get the data type of the set
        hid_t tDataType = H5Dget_type( tDataSet );

        // make sure that datatype fits to type of matrix
        /*if (       H5Tget_class( tDataType )
                != H5Tget_class( H5T_NATIVE_HBOOL ) )
        {
            std::string tMessage = "ERROR in reading from file: scalar "
                                   + aLabel + " has the wrong datatype.\n";

            MORIS_ERROR( false, tMessage.c_str() );
        }*/

        // get handler to data space
        hid_t tDataSpace = H5Dget_space( tDataSet );

        // value to cast bool to
        hbool_t tValue;

        // read data from file
        aStatus = H5Dread(
                tDataSet,
                tDataType,
                H5S_ALL,
                H5S_ALL,
                H5P_DEFAULT,
                &tValue );

        // cast output value
        aValue = (bool)tValue;

        // Close/release resources
        H5Tclose( tDataType );
        H5Dclose( tDataSet );
        H5Sclose( tDataSpace );

        // check for error
        MORIS_ASSERT( aStatus == 0, "Error in HDF5 load_scalar_from_hdf5_file()" );
    }

    //------------------------------------------------------------------------------
    inline void
    save_string_to_hdf5_file(
            hid_t&             aFileID,
            const std::string& aLabel,
            const std::string& aValue,
            herr_t&            aStatus )
    {
        // test if data set exists
        MORIS_ERROR( !dataset_exists( aFileID, aLabel ),
                "The data set %s can not be created because it already exist.\n",
                aLabel.c_str() );

        // select data type for string
        hid_t tDataType = H5Tcopy( H5T_C_S1 );

        // set size of output type
        aStatus = H5Tset_size( tDataType, aValue.length() );

        // matrix dimensions
        hsize_t tDims[ 1 ] = { 1 };

        // create data space
        hid_t tDataSpace = H5Screate_simple( 1, tDims, nullptr );

        // create new dataset
        hid_t tDataSet = H5Dcreate(
                aFileID,
                aLabel.c_str(),
                tDataType,
                tDataSpace,
                H5P_DEFAULT,
                H5P_DEFAULT,
                H5P_DEFAULT );

        // write data into data set
        aStatus = H5Dwrite(
                tDataSet,
                tDataType,
                H5S_ALL,
                H5S_ALL,
                H5P_DEFAULT,
                aValue.c_str() );

        // close open hids
        H5Sclose( tDataSpace );
        H5Tclose( tDataType );
        H5Dclose( tDataSet );

        // check for error
        MORIS_ASSERT( aStatus == 0, "Error in HDF5 save_string_to_hdf5_file()" );
    }

    //------------------------------------------------------------------------------
    inline void
    load_string_from_hdf5_file(
            hid_t&             aFileID,
            const std::string& aLabel,
            std::string&       aValue,
            herr_t&            aStatus )
    {
        // test if data set exists
        MORIS_ERROR( dataset_exists( aFileID, aLabel ),
                "The data set %scan not be opened because it does not exist.\n",
                aLabel.c_str() );

        // open the data set
        hid_t tDataSet = H5Dopen1( aFileID, aLabel.c_str() );

        // get handler to data space
        hid_t tDataSpace = H5Dget_space( tDataSet );

        // get the data type of the set
        hid_t tDataType = H5Dget_type( tDataSet );

        // get length of string
        hsize_t tSize = H5Dget_storage_size( tDataSet );

        // allocate buffer for string
        char* tBuffer = (char*)malloc( tSize * sizeof( char ) );

        // load string from hdf5
        aStatus = H5Dread(
                tDataSet,
                tDataType,
                H5S_ALL,
                H5S_ALL,
                H5P_DEFAULT,
                tBuffer );

        // create string from buffer
        aValue.assign( tBuffer, tSize );

        // delete buffer
        free( tBuffer );

        // close open hids
        H5Sclose( tDataSpace );
        H5Tclose( tDataType );
        H5Dclose( tDataSet );

        // check for error
        MORIS_ASSERT( aStatus == 0, "Error in HDF5 load_string_from_hdf5_file()" );
    }

    //------------------------------------------------------------------------------
}    // namespace moris

#endif /* PROJECTS_MRS_IOS_SRC_HDF5_TOOLS_HPP_ */
