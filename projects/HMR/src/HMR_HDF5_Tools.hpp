/*
 * HMR_HDF5_Tools.hpp
 *
 *  Created on: Jul 1, 2018
 *      Author: messe
 */

#ifndef SRC_HMR_HMR_HDF5_TOOLS_HPP_
#define SRC_HMR_HMR_HDF5_TOOLS_HPP_

// HD5 c-interface
#include "hdf5.h"

#include "typedefs.hpp" //COR/src
#include "cl_Mat.hpp" //LNA/src

namespace moris
{
    namespace hmr
    {
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
        template < typename T > hid_t
        get_hdf5_datatype( const T & aSample )
        {
            MORIS_ASSERT( false , "get_hdf5_datatype: unknown data type.");
            return NULL;
        }

//------------------------------------------------------------------------------

        template <> hid_t
        get_hdf5_datatype( const int & aSample )
        {
            return H5T_NATIVE_INT;
        }

//------------------------------------------------------------------------------

        template <> hid_t
        get_hdf5_datatype( const long int & aSample )
        {
            return H5T_NATIVE_LONG;
        }

//------------------------------------------------------------------------------

        // moris::uint
        template <> hid_t
        get_hdf5_datatype( const unsigned int & aSample )
        {
            return H5T_NATIVE_UINT;
        }

//------------------------------------------------------------------------------

        // moris::luint
        template <> hid_t
        get_hdf5_datatype( const long unsigned int & aSample )
        {
            return H5T_NATIVE_ULONG;
        }

//------------------------------------------------------------------------------

        template <> hid_t
        get_hdf5_datatype( const double & aSample )
        {
            return H5T_NATIVE_DOUBLE;
        }

//------------------------------------------------------------------------------

        template <> hid_t
        get_hdf5_datatype( const long double & aSample )
        {
            return H5T_NATIVE_LDOUBLE;
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
        template < typename T >
        void
        save_matrix_to_hdf5_file(
                hid_t               & aFileID,
                const std::string   & aLabel,
                const Mat< T >      & aMatrix,
                herr_t              & aStatus
                )
        {
            // matrix dimensions
            hsize_t tDims[ 2 ];

            // get dimensions from matrix
            tDims[ 0 ] = aMatrix.n_rows();
            tDims[ 1 ] = aMatrix.n_cols();

            // allocate top level array which contains rows
            T** tData = (T**) malloc( tDims[ 0 ]*sizeof( T* ) );

            // allocate memory for data
            tData[ 0 ] = ( T* ) malloc( tDims[ 0 ]*  tDims[ 1 ] * sizeof( T ) );

            // loop over all rows and allocate colums
            for( hsize_t i=0; i<tDims[ 0 ]; ++i )
            {
                tData[ i ] = tData[0]+ i*tDims[ 1 ];
            }

            // convert moris::Mat to data
            for ( hsize_t i = 0; i < tDims[ 0 ]; ++i )
            {
                for ( hsize_t j = 0; j < tDims[ 1 ]; ++j )
                {
                    tData[ i ][ j ] = aMatrix( i, j );
                }
            }

            // create data space
            hid_t  tDataSpace
                = H5Screate_simple( 2, tDims, NULL);

            // select data type for matrix to save
            hid_t tDataType = H5Tcopy( get_hdf5_datatype( ( T ) 0 ) );

            // set data type to little endian
            aStatus = H5Tset_order( tDataType, H5T_ORDER_LE );

            // create new dataset
            hid_t tDataSet = H5Dcreate(
                    aFileID,
                    aLabel.c_str(),
                    tDataType,
                    tDataSpace,
                    H5P_DEFAULT,
                    H5P_DEFAULT,
                    H5P_DEFAULT );

            // write data into dataset
            aStatus = H5Dwrite(
                    tDataSet,
                    tDataType,
                    H5S_ALL,
                    H5S_ALL,
                    H5P_DEFAULT,
                    &tData[ 0 ][ 0 ]);

            // tidy up memory
            free( tData[ 0 ] );
            free( tData );

            // close open hids
            H5Sclose( tDataSpace );
            H5Tclose( tDataType );
            H5Dclose( tDataSet );

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
        template < typename T >
        void
        load_matrix_from_hdf5_file(
                hid_t               & aFileID,
                const std::string   & aLabel,
                Mat< T >            & aMatrix,
                herr_t              & aStatus
        )
        {


            // open the data set
            hid_t tDataSet = H5Dopen1( aFileID, aLabel.c_str() );

            // get the data type of the set
            hid_t tDataType = H5Dget_type( tDataSet );

            // make sure that datatype fits to type of matrix
            if (       H5Tget_class( tDataType )
                    != H5Tget_class( get_hdf5_datatype( ( T ) 0 ) ) )
            {
                std::fprintf( stdout,"ERROR in reading from file: field %s has the wrong datatype.\n",
                        aLabel.c_str() );
                exit( - 1 );
            }

            // get handler to dataspace
            hid_t tDataSpace = H5Dget_space( tDataSet );

            // matrix dimensions
            hsize_t tDims[ 2 ];

            // ask hdf for dimensions
            aStatus  = H5Sget_simple_extent_dims( tDataSpace, tDims, NULL);

            // allocate top level array which contains rows
            T** tData = (T**) malloc( tDims[ 0 ]*sizeof( T* ) );

            // allocate memory for data
            tData[ 0 ] = ( T* ) malloc( tDims[ 0 ]*  tDims[ 1 ] * sizeof( T ) );

            // loop over all rows and allocate colums
            for( hsize_t i=1; i<tDims[ 0 ]; ++i )
            {
                tData[ i ] = tData[ 0 ]+ i*tDims[ 1 ];
            }


            // read data from file
            aStatus = H5Dread(
                    tDataSet,
                    tDataType,
                    H5S_ALL,
                    H5S_ALL,
                    H5P_DEFAULT,
                    &tData[ 0 ][ 0 ] );

            // Close/release resources
            H5Tclose( tDataType );
            H5Dclose( tDataSet );
            H5Sclose( tDataSpace );

            // allocate memory for output
            aMatrix.set_size( tDims[ 0 ], tDims[ 1 ] );

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
        template < typename T >
        void
        save_scalar_to_hdf5_file(
                hid_t               & aFileID,
                const std::string   & aLabel,
                const T             & aValue,
                herr_t              & aStatus
        )
        {
            // select data type for matrix to save
            hid_t tDataType = H5Tcopy( get_hdf5_datatype( ( T ) 0 ) );

            // set data type to little endian
            aStatus = H5Tset_order( tDataType, H5T_ORDER_LE );

            // matrix dimensions
            hsize_t tDims[ 1 ] = { 1 };

            // create data space
            hid_t tDataSpace
                = H5Screate_simple( 1, tDims, NULL );

            // create new dataset
            hid_t tDataSet = H5Dcreate(
                    aFileID,
                    aLabel.c_str(),
                    tDataType,
                    tDataSpace,
                    H5P_DEFAULT,
                    H5P_DEFAULT,
                    H5P_DEFAULT );

            // write data into dataset
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
        }

//------------------------------------------------------------------------------

        template <>
        void
        save_scalar_to_hdf5_file(
                hid_t               & aFileID,
                const std::string   & aLabel,
                const bool          & aValue,
                herr_t              & aStatus
        )
        {
            // select data type for matrix to save
            hid_t tDataType = H5Tcopy( H5T_NATIVE_HBOOL );

            // matrix dimensions
            hsize_t tDims[ 1 ] = { 1 };

            // create data space
            hid_t tDataSpace
            = H5Screate_simple( 1, tDims, NULL );

            // create new dataset
            hid_t tDataSet = H5Dcreate(
                    aFileID,
                    aLabel.c_str(),
                    tDataType,
                    tDataSpace,
                    H5P_DEFAULT,
                    H5P_DEFAULT,
                    H5P_DEFAULT );

            // value to cast bool to
            hbool_t tValue = ( hbool_t ) aValue;

            // write data into dataset
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
        template < typename T >
        void
        load_scalar_from_hdf5_file(
                hid_t               & aFileID,
                const std::string   & aLabel,
                T                   & aValue,
                herr_t              & aStatus
        )
        {
            // open the data set
            hid_t tDataSet = H5Dopen1( aFileID, aLabel.c_str() );

            // get the data type of the set
            hid_t tDataType = H5Dget_type( tDataSet );

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
        }

//------------------------------------------------------------------------------

        template <>
        void
        load_scalar_from_hdf5_file(
                hid_t               & aFileID,
                const std::string   & aLabel,
                bool                & aValue,
                herr_t              & aStatus
        )
        {
            // open the data set
            hid_t tDataSet = H5Dopen1( aFileID, aLabel.c_str() );

            // get the data type of the set
            hid_t tDataType = H5Dget_type( tDataSet );

            // make sure that datatype fits to type of matrix
            if (       H5Tget_class( tDataType )
                    != H5Tget_class( H5T_NATIVE_HBOOL ) )
            {
                std::fprintf( stdout,"ERROR in reading from file: field %s has the wrong datatype.\n",
                        aLabel.c_str() );
                exit( -1 );
            }

            // get handler to dataspace
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
            aValue = ( bool ) tValue;

            // Close/release resources
            H5Tclose( tDataType );
            H5Dclose( tDataSet );
            H5Sclose( tDataSpace );
        }

//------------------------------------------------------------------------------
    }
}



#endif /* SRC_HMR_HMR_HDF5_TOOLS_HPP_ */
