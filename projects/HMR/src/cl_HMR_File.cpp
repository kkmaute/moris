/*
 * cl_HMR_File.cpp
 *
 *  Created on: Jun 29, 2018
 *      Author: messe
 */

#include "HMR_HDF5_Tools.hpp" //HMR/src
#include "cl_HMR_Factory.hpp" //HMR/src
#include "cl_HMR_File.hpp" //HMR/src

namespace moris
{
    namespace hmr
    {

//------------------------------------------------------------------------------

    void
    File::create( const std::string & aPath )
    {
        // Create a new file using default properties
        mFileID = H5Fcreate(
                this->parralize_filename( aPath ).c_str(),
                H5F_ACC_TRUNC,
                H5P_DEFAULT,
                H5P_DEFAULT);
    }

//------------------------------------------------------------------------------

    void
    File::open( const std::string & aPath )
    {
        // opens an existing file with read and write access
        mFileID = H5Fopen( this->parralize_filename( aPath ).c_str(),
                H5F_ACC_RDWR,
                H5P_DEFAULT);
    }

//------------------------------------------------------------------------------

    void
    File::close()
    {
        // close the hdf file
        mStatus = H5Fclose( mFileID );
    }

//------------------------------------------------------------------------------

        void
        File::save_settings( const Parameters * aParameters )
        {
            // save dimensions of field
            save_matrix_to_hdf5_file(
                    mFileID,
                    "DomainDimensions",
                    aParameters->get_domain_dimensions(),
                    mStatus );

            // save domain offset
            save_matrix_to_hdf5_file(
                    mFileID,
                    "DomainOffset",
                    aParameters->get_domain_offset(),
                    mStatus );

            // save number of elements on coarsest mesh
            save_matrix_to_hdf5_file(
                    mFileID,
                    "CoarsestElements",
                    aParameters->get_number_of_elements_per_dimension(),
                    mStatus );

            // save buffer size
            save_scalar_to_hdf5_file(
                    mFileID,
                    "BufferSize",
                    aParameters->get_buffer_size(),
                    mStatus);

            // save max polynomial
            /*save_scalar_to_hdf5_file(
                    mFileID,
                    "MaxPolynomial",
                    aParameters->get_max_polynomial(),
                    mStatus ); */

            // save verbosity flag
            save_scalar_to_hdf5_file(
                    mFileID,
                    "VerboseFlag",
                    aParameters->is_verbose(),
                    mStatus );

            // save truncation flag
            save_scalar_to_hdf5_file(
                    mFileID,
                    "BSplineTruncationFlag",
                    aParameters->truncate_bsplines(),
                    mStatus );

            // save max volume level
            save_scalar_to_hdf5_file(
                    mFileID,
                    "MaxVolumeLevel",
                    aParameters->get_max_volume_level(),
                    mStatus );

            // save max surface level
            save_scalar_to_hdf5_file(
                    mFileID,
                    "MaxSurfaceLevel",
                    aParameters->get_max_surface_level(),
                    mStatus );

            // save mesh scaling factor for gmsh
            save_scalar_to_hdf5_file(
                    mFileID,
                    "GmshScale",
                    aParameters->get_gmsh_scale(),
                    mStatus );

            // save Lagrange mesh associations
            save_matrix_to_hdf5_file(
                    mFileID,
                    "LagrangeOrders",
                    aParameters->get_lagrange_orders(),
                    mStatus );

            // save Lagrange mesh associations
            save_matrix_to_hdf5_file(
                    mFileID,
                    "LagrangePatterns",
                    aParameters->get_lagrange_patterns(),
                    mStatus );

            // save bspline mesh associations
            save_matrix_to_hdf5_file(
                    mFileID,
                    "BSplineOrders",
                    aParameters->get_bspline_orders(),
                    mStatus );

            // save bspline mesh associations
            save_matrix_to_hdf5_file(
                    mFileID,
                    "BSplinePatterns",
                    aParameters->get_bspline_patterns(),
                    mStatus );

            // save linking flags
            save_matrix_to_hdf5_file(
                    mFileID,
                    "LagrangeToBSpline",
                    aParameters->get_lagrange_to_bspline(),
                    mStatus );
        }

//------------------------------------------------------------------------------

        void
        File::load_settings( Parameters * aParameters )
        {
            // placeholders for data read from file
            Matrix< DDRMat >  tMatReal;
            Matrix< DDLUMat > tMatLuint;
            Matrix< DDUMat > tMatUint;
            real         tValReal;
            uint         tValUint;
            luint        tValLuint;
            bool         tValBool;

            // load dimensions from field
            load_matrix_from_hdf5_file(
                    mFileID,
                    "DomainDimensions",
                    tMatReal,
                    mStatus );

            // set domain dimensions
            aParameters->set_domain_dimensions( tMatReal );

            // load domain offset
            load_matrix_from_hdf5_file(
                    mFileID,
                    "DomainOffset",
                    tMatReal,
                    mStatus );

            // set domain offset
            aParameters->set_domain_offset( tMatReal );

            // load number of elements on coarsest mesh
            load_matrix_from_hdf5_file(
                    mFileID,
                    "CoarsestElements",
                    tMatLuint,
                    mStatus );

            // set number of elements
            aParameters->set_number_of_elements_per_dimension( tMatLuint );

            // load buffer size
            load_scalar_from_hdf5_file(
                    mFileID,
                    "BufferSize",
                    tValLuint,
                    mStatus);

            // set buffer size
            aParameters->set_buffer_size( tValLuint );

            // load max polynomial
            /*load_scalar_from_hdf5_file(
                    mFileID,
                    "MaxPolynomial",
                    tValLuint,
                    mStatus ); */

            // set max polynomial
            //aParameters->set_max_polynomial( tValLuint );

            // load truncation flag
            load_scalar_from_hdf5_file(
                    mFileID,
                    "BSplineTruncationFlag",
                    tValBool,
                    mStatus );

            // set truncation flag
            aParameters->set_bspline_truncation( tValBool );

            // load verbosity flag
            load_scalar_from_hdf5_file(
                    mFileID,
                    "VerboseFlag",
                    tValBool,
                    mStatus );

            // set verbose flag
            aParameters->set_verbose( tValBool );

            // load max volume level
            load_scalar_from_hdf5_file(
                    mFileID,
                    "MaxVolumeLevel",
                    tValUint,
                    mStatus );

            // set max volume level
            aParameters->set_max_volume_level( tValUint );

            // load max surface level
            load_scalar_from_hdf5_file(
                    mFileID,
                    "MaxSurfaceLevel",
                    tValUint,
                    mStatus );

            // set max surface level
            aParameters->set_max_surface_level( tValUint );

            // load scaling factor for gmsh
            load_scalar_from_hdf5_file(
                    mFileID,
                    "GmshScale",
                    tValReal,
                    mStatus );

            // set scaling factor for gmsh
            aParameters->set_gmsh_scale( tValReal );

            // load orders of meshes
            load_matrix_from_hdf5_file(
                    mFileID,
                    "LagrangeOrders",
                    tMatUint,
                    mStatus );

            aParameters->set_lagrange_orders( tMatUint );

            // load Lagrange mesh associations
            load_matrix_from_hdf5_file(
                    mFileID,
                    "LagrangePatterns",
                    tMatUint,
                    mStatus );

            aParameters->set_lagrange_patterns( tMatUint );

            // load orders of meshes
            load_matrix_from_hdf5_file(
                    mFileID,
                    "BSplineOrders",
                    tMatUint,
                    mStatus );

            aParameters->set_bspline_orders( tMatUint );

            // load bspline mesh associations
            load_matrix_from_hdf5_file(
                    mFileID,
                    "BSplinePatterns",
                    tMatUint,
                    mStatus );

            aParameters->set_bspline_patterns( tMatUint );

            // load lagrange to bspline links
            // save linking flags
            load_matrix_from_hdf5_file(
                    mFileID,
                    "LagrangeToBSpline",
                    tMatUint,
                    mStatus );

            aParameters->set_lagrange_to_bspline( tMatUint );
        }

//------------------------------------------------------------------------------

        void
        File::save_refinement_pattern(
                Background_Mesh_Base * aMesh,
                const uint           & aPattern )
        {
            // step 1: count how many elements need are refined on each level
            uint tMaxLevel = aMesh->get_max_level();

            // element counter
            Matrix< DDLUMat > tElementCounter ( tMaxLevel+1, 1, 0 );

            // activate this pattern
            aMesh->set_activation_pattern( aPattern );

            // collect all elements that are flagged for refinement
            for( uint l=0; l<tMaxLevel; ++l )
            {
                // cell which contains elements
                Cell< Background_Element_Base* > tElements;

                // collect elements from this level
                aMesh->collect_elements_on_level_within_proc_domain( l, tElements );

                // loop over all elements
                for( auto tElement : tElements )
                {
                    // test if element is refined
                    if( tElement->is_refined( aPattern ) )
                    {
                        // increment counter
                        ++tElementCounter ( l );
                    }
                }
            }

            // allocate pattern
            luint* tPattern = new luint[ sum( tElementCounter ) ];

            // write values into pattern
            // collect all elements that are flagged for refinement
            hsize_t tCount = 0;
            for( uint l=0; l<tMaxLevel; ++l )
            {
                // cell which contains elements
                Cell< Background_Element_Base* > tElements;

                // collect elements from this level
                aMesh->collect_elements_on_level_within_proc_domain( l, tElements );

                // reset element counter
                luint tElementCount = 0;

                // loop over all elements
                for( auto tElement : tElements )
                {
                    // test if element is refined
                    if( tElement->is_refined( aPattern ) )
                    {
                        tPattern[ tCount++ ] = tElementCount;
                    }

                    // deactive elements must not be counted.
                    if( ! tElement->is_deactive( aPattern ) )
                    {
                        ++tElementCount;
                    }
                }
            }

            // create dataspace
            hid_t  tDataSpace = H5Screate_simple( 1, &tCount, NULL );

            // select data type for matrix to save
            hid_t tDataType = H5Tcopy( get_hdf5_datatype( ( luint ) 0 ) );

            // set data type to little endian
            mStatus = H5Tset_order( tDataType, H5T_ORDER_LE );

            // count number of active elements
            save_scalar_to_hdf5_file(
                    mFileID,
                    "ActiveElements",
                    aMesh->get_number_of_active_elements_on_proc(),
                    mStatus );

            std::string tLabel = "RefinementPattern";

            // create new dataset
            hid_t tDataSet = H5Dcreate(
                    mFileID,
                    tLabel.c_str(),
                    tDataType,
                    tDataSpace,
                    H5P_DEFAULT,
                    H5P_DEFAULT,
                    H5P_DEFAULT );

            // write data into dataset
            mStatus = H5Dwrite(
                    tDataSet,
                    tDataType,
                    H5S_ALL,
                    H5S_ALL,
                    H5P_DEFAULT,
                    &tPattern[ 0 ] );

            // tidy up memory
            delete[] tPattern;

            // close open hids
            H5Sclose( tDataSpace );
            H5Tclose( tDataType );
            H5Dclose( tDataSet );

            // create name
            //std::string tCounterlabel
            //    = "RefinementCounter_" + std::to_string( aMesh->get_activation_pattern() );

            std::string tCounterlabel = "RefinementCounter";

            // save counter
            save_matrix_to_hdf5_file(
                    mFileID,
                    tCounterlabel.c_str(),
                    tElementCounter,
                    mStatus );

        }

//------------------------------------------------------------------------------

        void
        File::load_refinement_pattern(
                Background_Mesh_Base * aMesh,
                const uint           & aPattern )
        {

            // matrix containing counter
            Matrix< DDLUMat > tElementCounter;

            // activate this pattern
            aMesh->set_activation_pattern( aPattern );


           // std::string tCounterlabel
           //     = "RefinementCounter_" + std::to_string( tActivePattern );

            std::string tCounterlabel  = "RefinementCounter";

            // load counter
            load_matrix_from_hdf5_file(
                    mFileID,
                    tCounterlabel.c_str(),
                    tElementCounter,
                    mStatus );

            // get number of levels
            uint tNumberOfLevels = tElementCounter.length();


            // allocate pattern
            luint* tPattern = new luint[ sum( tElementCounter ) ];

            // create name
            //std::string tLabel
            //    = "RefinementPattern_" + std::to_string( aMesh->get_activation_pattern() );
            std::string tLabel  = "RefinementPattern";

            // open the data set
            hid_t tDataSet = H5Dopen1( mFileID, tLabel.c_str() );

            // get the data type of the set
            hid_t tDataType = H5Dget_type( tDataSet );

            // make sure that datatype fits to type of matrix
            if (       H5Tget_class( tDataType )
                    != H5Tget_class( get_hdf5_datatype( ( luint ) 0 ) ) )
            {
                std::fprintf( stdout,"ERROR in reading from file: field RefinementPattern has the wrong datatype.\n" );
                exit( -1 );
            }



            // get handler to dataspace
            hid_t tDataSpace = H5Dget_space( tDataSet );

            // read data
            mStatus = H5Dread(
                    tDataSet,
                    tDataType,
                    H5S_ALL,
                    H5S_ALL,
                    H5P_DEFAULT,
                    &tPattern[ 0 ] );

            // Close/release resources
            H5Tclose( tDataType );
            H5Dclose( tDataSet );
            H5Sclose( tDataSpace );

            // counter for array
            hsize_t tCount = 0;

            // loop over all levels
            for( uint l=0; l<tNumberOfLevels; ++l )
            {
                // cell which contains elements
                Cell< Background_Element_Base* > tElements;

                // collect elements from this level
                aMesh->collect_elements_on_level_within_proc_domain( l, tElements );



                luint tNumberOfElements = tElementCounter( l );

                for( luint k=0; k<tNumberOfElements; ++k )
                {
                    tElements( tPattern[ tCount++ ] )->put_on_refinement_queue();
                }

                // refine mesh
                aMesh->perform_refinement();
            }

            aMesh->update_database();


            // get number of active elements
            luint tNumberOfElements;
            load_scalar_from_hdf5_file(
                    mFileID,
                    "ActiveElements",
                    tNumberOfElements,
                    mStatus );

            MORIS_ERROR(
                    aMesh->get_number_of_active_elements_on_proc() == tNumberOfElements,
                    "Error in loading HDF5 file" );

            // tidy up memory
            delete [] tPattern;
        }

//-------------------------------------------------------------------------------

        std::string
        File::parralize_filename( const std::string & aPath )
        {

            // test if running in parallel mode
            if ( par_size() > 1 )
            {
                // get file extesion
                auto tFileExt = aPath.substr(
                        aPath.find_last_of("."),
                        aPath.length() );

                // get base path
                auto tBasePath = aPath.substr(
                        0,
                        aPath.find_last_of(".") );

                // add proc number to path
                std::string aParallelPath = tBasePath + "_"
                                           +  std::to_string( par_rank() )
                                           + tFileExt;
                return aParallelPath;
            }
            else
            {
                // do not modify path
                return aPath;
            }
        }

//-------------------------------------------------------------------------------

        /**
         * free function needed by loading constructor
         */
        Parameters *
        create_hmr_parameters_from_hdf5_file( const std::string & aPath )
        {
            // create file object
            File tHDF5;

            // open file on disk
            tHDF5.open( aPath );

            // create new parameter pointer
            Parameters * aParameters = new Parameters;

            // load settings
            tHDF5.load_settings( aParameters );

            // close file
            tHDF5.close();

            // return pointer
            return aParameters;
        }

    }
} /* namespace moris */
