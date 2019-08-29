/*
 * cl_HMR_File.cpp
 *
 *  Created on: Jun 29, 2018
 *      Author: messe
 */

#include "cl_HMR_File.hpp" //HMR/src

#include "cl_HMR_Factory.hpp" //HMR/src
#include "HDF5_Tools.hpp" //HMR/src
#include "cl_Map.hpp"

namespace moris
{
    namespace hmr
    {

//------------------------------------------------------------------------------

    void File::create( const std::string & aPath )
    {
        // Create a new file using default properties
        mFileID = create_hdf5_file( aPath );
    }

//------------------------------------------------------------------------------

    void File::open( const std::string & aPath )
    {
        // opens an existing file with read and write access
        mFileID = open_hdf5_file( aPath );
    }

//------------------------------------------------------------------------------

    void File::close()
    {
        // close the hdf file
        mStatus = close_hdf5_file( mFileID );
    }

//------------------------------------------------------------------------------

        void File::save_settings( const Parameters * aParameters )
        {
            // save dimensions of field
            save_matrix_to_hdf5_file( mFileID,
                                      "DomainDimensions",
                                      aParameters->get_domain_dimensions(),
                                      mStatus );

            // save domain offset
            save_matrix_to_hdf5_file( mFileID,
                                      "DomainOffset",
                                      aParameters->get_domain_offset(),
                                      mStatus );

            // save number of elements on coarsest mesh
            save_matrix_to_hdf5_file( mFileID,
                                      "CoarsestElements",
                                      aParameters->get_number_of_elements_per_dimension(),
                                      mStatus );

            // save buffer size
            save_scalar_to_hdf5_file( mFileID,
                                      "RefinementBuffer",
                                      aParameters->get_refinement_buffer(),
                                      mStatus);

            // save max polynomial
            /*save_scalar_to_hdf5_file(
                    mFileID,
                    "MaxPolynomial",
                    aParameters->get_max_polynomial(),
                    mStatus ); */

            // save verbosity flag
            save_scalar_to_hdf5_file( mFileID,
                                      "SeverityFlag",
                                      gLogger.get_severity_level(),
                                      mStatus );

            // save multigrid flag
            save_scalar_to_hdf5_file( mFileID,
                                      "MultigridFlag",
                                      aParameters->use_multigrid(),
                                      mStatus );

            // save truncation flag
            save_scalar_to_hdf5_file( mFileID,
                                      "BSplineTruncationFlag",
                                      aParameters->truncate_bsplines(),
                                      mStatus );

            // save initial refinement
            save_scalar_to_hdf5_file( mFileID,
                                      "InitialBSplineRefinement",
                                      aParameters->get_initial_refinement(),
                                      mStatus );

            // save initial refinement
            save_scalar_to_hdf5_file( mFileID,
                                      "AdditionalLagrangeRefinement",
                                      aParameters->get_additional_lagrange_refinement(),
                                      mStatus );

            // save maximal refinement level
            save_scalar_to_hdf5_file( mFileID,
                                      "MaxRefinementLevel",
                                      aParameters->get_max_refinement_level(),
                                      mStatus );

            // save mesh scaling factor for gmsh
            save_scalar_to_hdf5_file( mFileID,
                                      "GmshScale",
                                      aParameters->get_gmsh_scale(),
                                      mStatus );

            // save Lagrange mesh associations
            save_matrix_to_hdf5_file( mFileID,
                                      "LagrangeOrders",
                                      aParameters->get_lagrange_orders(),
                                      mStatus );

            // save Lagrange mesh associations
            save_matrix_to_hdf5_file( mFileID,
                                      "LagrangePatterns",
                                      aParameters->get_lagrange_patterns(),
                                      mStatus );

            // save bspline mesh associations
            save_matrix_to_hdf5_file( mFileID,
                                      "BSplineOrders",
                                      aParameters->get_bspline_orders(),
                                      mStatus );

            // save bspline mesh associations
            save_matrix_to_hdf5_file( mFileID,
                                      "BSplinePatterns",
                                      aParameters->get_bspline_patterns(),
                                      mStatus );

            // save B-Spline Maps
            save_matrix_to_hdf5_file( mFileID,
                                      "BSplineInputMap",
                                      aParameters->get_bspline_input_map(),
                                      mStatus );

            save_matrix_to_hdf5_file( mFileID,
                                      "BSplineOutputMap",
                                      aParameters->get_bspline_output_map(),
                                      mStatus );

            // save Sidesets
            Matrix< DDUMat > tSideSets = aParameters->get_side_sets();
            if ( tSideSets.length() == 0 )
            {
                tSideSets.set_size( 1, 1, 0 );
            }

            save_matrix_to_hdf5_file( mFileID,
                                      "SideSets",
                                      tSideSets,
                                      mStatus );
        }

//------------------------------------------------------------------------------

        void File::load_settings( Parameters * aParameters )
        {
            // placeholders for data read from file
            Matrix< DDRMat >  tMatReal;
            Matrix< DDLUMat > tMatLuint;
            Matrix< DDUMat >  tMatUint;
            real              tValReal;
            uint              tValUint;
            sint              tValSint;
            luint             tValLuint;
            bool              tValBool;

            // load dimensions from field
            load_matrix_from_hdf5_file( mFileID,
                                        "DomainDimensions",
                                        tMatReal,
                                        mStatus );

            // set domain dimensions
            aParameters->set_domain_dimensions( tMatReal );

            // load domain offset
            load_matrix_from_hdf5_file( mFileID,
                                        "DomainOffset",
                                        tMatReal,
                                        mStatus );

            // set domain offset
            aParameters->set_domain_offset( tMatReal );

            // load number of elements on coarsest mesh
            load_matrix_from_hdf5_file( mFileID,
                                        "CoarsestElements",
                                        tMatLuint,
                                        mStatus );

            // set number of elements
            aParameters->set_number_of_elements_per_dimension( tMatLuint );

            // load buffer size
            load_scalar_from_hdf5_file( mFileID,
                                        "RefinementBuffer",
                                        tValLuint,
                                        mStatus);

            // set buffer size
            aParameters->set_refinement_buffer( tValLuint );

            // load max polynomial
            /*load_scalar_from_hdf5_file( mFileID,
                                          "MaxPolynomial",
                                          tValLuint,
                                          mStatus ); */

            // set max polynomial
            //aParameters->set_max_polynomial( tValLuint );

            // load truncation flag
            load_scalar_from_hdf5_file( mFileID,
                                        "BSplineTruncationFlag",
                                        tValBool,
                                        mStatus );

            // set truncation flag
            aParameters->set_bspline_truncation( tValBool );

            // load verbosity flag
            load_scalar_from_hdf5_file( mFileID,
                                        "SeverityFlag",
                                        tValSint,
                                        mStatus );

            // set verbose flag
            aParameters->set_severity_level( tValSint );

            // load multigrid flag
            load_scalar_from_hdf5_file( mFileID,
                                        "MultigridFlag",
                                        tValBool,
                                        mStatus );

            // set Multigrid flag
            aParameters->set_multigrid( tValBool );

            // load initial refinement
            load_scalar_from_hdf5_file( mFileID,
                                        "InitialBSplineRefinement",
                                        tValUint,
                                        mStatus );

            aParameters->set_initial_refinement( tValUint );

            // load initial refinement
            load_scalar_from_hdf5_file( mFileID,
                                        "AdditionalLagrangeRefinement",
                                        tValUint,
                                        mStatus );

            aParameters->set_additional_lagrange_refinement( tValUint );

            // loadmaximal refinement level
            load_scalar_from_hdf5_file( mFileID,
                                        "MaxRefinementLevel",
                                        tValUint,
                                        mStatus );
            aParameters->set_max_refinement_level( tValUint );

            // load scaling factor for gmsh
            load_scalar_from_hdf5_file( mFileID,
                                        "GmshScale",
                                        tValReal,
                                        mStatus );

            // set scaling factor for gmsh
            aParameters->set_gmsh_scale( tValReal );

            // load orders of meshes
            load_matrix_from_hdf5_file( mFileID,
                                        "LagrangeOrders",
                                        tMatUint,
                                        mStatus );

            aParameters->set_lagrange_orders( tMatUint );

            // load Lagrange mesh associations
            load_matrix_from_hdf5_file( mFileID,
                                        "LagrangePatterns",
                                        tMatUint,
                                        mStatus );

            aParameters->set_lagrange_patterns( tMatUint );

            // load orders of meshes
            load_matrix_from_hdf5_file( mFileID,
                                        "BSplineOrders",
                                        tMatUint,
                                        mStatus );

            aParameters->set_bspline_orders( tMatUint );

            // load bspline mesh associations
            load_matrix_from_hdf5_file( mFileID,
                                        "BSplinePatterns",
                                        tMatUint,
                                        mStatus );

            aParameters->set_bspline_patterns( tMatUint );

            // load B-Spline input maps
            load_matrix_from_hdf5_file( mFileID,
                                        "BSplineInputMap",
                                        tMatUint,
                                        mStatus );

            aParameters->set_bspline_input_map( tMatUint );

            // load B-Spline output maps
            load_matrix_from_hdf5_file( mFileID,
                                        "BSplineOutputMap",
                                        tMatUint,
                                        mStatus );
            aParameters->set_bspline_output_map( tMatUint );

            // load side sets
            load_matrix_from_hdf5_file( mFileID,
                                        "SideSets",
                                        tMatUint,
                                        mStatus );

            // test if matrix has values
            if( tMatUint.length() > 0 )
            {
                if( tMatUint( 0 ) != 0 )
                {
                    // reset matrix
                    aParameters->set_side_sets( tMatUint );
                }
            }
        }

//------------------------------------------------------------------------------

        void File::save_refinement_pattern( Background_Mesh_Base * aMesh  )
        {
            // step 1: count how many elements need are refined on each level
            uint tMaxLevel = aMesh->get_max_level();

            // element counter
            Matrix< DDLUMat > tElementCounter ( tMaxLevel+1, 2, 0 );

            uint tBSplinePattern  = aMesh->get_parameters()->get_bspline_output_pattern();
            uint tLagrangePattern = aMesh->get_parameters()->get_lagrange_output_pattern();

            // collect all elements that are flagged for refinement
            for( uint l = 0; l < tMaxLevel; ++l )
            {
                // cell which contains elements
                Cell< Background_Element_Base* > tElements;

                // collect elements from this level
                aMesh->collect_elements_on_level_within_proc_domain( l, tElements );

                // loop over all elements
                for( auto tElement : tElements )
                {
                    // test if B-Spline Element is refined
                    if( tElement->is_refined( tBSplinePattern ) )
                    {
                        // increment counter
                        ++tElementCounter ( l, 0 );
                    }
                    else if( tElement->is_refined( tLagrangePattern ) )
                    {
                        // increment counter
                        ++tElementCounter ( l, 1 );
                    }
                }
            }

            // allocate patterns
            Matrix< DDLUMat > tBSplineElements  ( sum( tElementCounter.get_column( 0 ) ), 1 );
            Matrix< DDLUMat > tLagrangeElements ( sum( tElementCounter.get_column( 1 ) ), 1 );

            // reset counters
            hsize_t tBSplineCount = 0;
            hsize_t tLagrangeCount = 0;

            for( uint l=0; l<tMaxLevel; ++l )
            {
                // cell which contains elements
                Cell< Background_Element_Base* > tElements;

                // collect elements from this level
                aMesh->collect_elements_on_level_within_proc_domain( l, tElements );

                // loop over all elements
                for( Background_Element_Base * tElement : tElements )
                {
                    // test if element is refined
                    if( tElement->is_refined( tBSplinePattern ) )
                    {
                        tBSplineElements( tBSplineCount++ ) = tElement->get_hmr_id();
                    }
                    else if ( tElement->is_refined( tLagrangePattern ) )
                    {
                        tLagrangeElements( tLagrangeCount++ ) = tElement->get_hmr_id();
                    }
                }
            }

            save_matrix_to_hdf5_file( mFileID,
                                      "ElementCounter",
                                      tElementCounter,
                                      mStatus );

            save_matrix_to_hdf5_file( mFileID,
                                      "BSplineElements",
                                      tBSplineElements,
                                      mStatus );

            save_matrix_to_hdf5_file( mFileID,
                                      "LagrangeElements",
                                      tLagrangeElements,
                                      mStatus );
        }

//------------------------------------------------------------------------------

        void File::load_refinement_pattern( Background_Mesh_Base * aMesh,
                                            const bool             aMode )
        {
            uint tBSplinePattern;
            uint tLagrangePattern;

            // check if we are loading input or output
            if( aMode )
            {
                tBSplinePattern  = aMesh->get_parameters()->get_bspline_output_pattern();
                tLagrangePattern = aMesh->get_parameters()->get_lagrange_output_pattern();
            }
            else
            {
                tBSplinePattern  = aMesh->get_parameters()->get_bspline_input_pattern();
                tLagrangePattern = aMesh->get_parameters()->get_lagrange_input_pattern();
            }
            // matrix containing counter
            Matrix< DDLUMat > tElementCounter;

            // load counter
            load_matrix_from_hdf5_file( mFileID,
                                        "ElementCounter",
                                        tElementCounter,
                                        mStatus );

            // allocate pattern
            Matrix< DDLUMat > tBSplineElements;
            load_matrix_from_hdf5_file( mFileID,
                                        "BSplineElements",
                                        tBSplineElements,
                                        mStatus );

            Matrix< DDLUMat > tLagrangeElements;
            load_matrix_from_hdf5_file( mFileID,
                                        "LagrangeElements",
                                        tLagrangeElements,
                                        mStatus );

            // get number of levels
            uint tNumberOfLevels = tElementCounter.n_rows();

            // select B-Spline pattern
            aMesh->set_activation_pattern( tBSplinePattern );

            // reset counter
            luint tCount = 0;

            // loop over all levels
            for( uint l=0; l<tNumberOfLevels; ++l )
            {
                // cell which contains elements
                Cell< Background_Element_Base* > tElements;

                // collect elements from this level
                aMesh->collect_elements_on_level_within_proc_domain( l, tElements );

                // create a map with ids
                map< moris_id, luint > tMap;

                luint j = 0;
                for( Background_Element_Base* tElement : tElements )
                {
                    tMap[ tElement->get_hmr_id() ] = j++;
                }

                luint tNumberOfElements = tElementCounter( l, 0 );

                for( luint k=0; k<tNumberOfElements; ++k )
                {
                    tElements( tMap.find( tBSplineElements( tCount++ ) ) )->put_on_refinement_queue();
                }

                // refine mesh
                aMesh->perform_refinement( tBSplinePattern );
            }

            // clone B-Spline to Lagrange
            aMesh->copy_pattern( tBSplinePattern, tLagrangePattern );

            // select Lagrange pattern
            aMesh->set_activation_pattern( tLagrangePattern );

            // reset counter
            tCount = 0;

            // loop over all levels
            for( uint l=0; l<tNumberOfLevels; ++l )
            {
                // cell which contains elements
                Cell< Background_Element_Base* > tElements;

                // collect elements from this level
                aMesh->collect_elements_on_level_within_proc_domain( l, tElements );

                // create a map with ids
                map< moris_id, luint > tMap;

                luint j = 0;
                for( Background_Element_Base* tElement : tElements )
                {
                    tMap[ tElement->get_hmr_id() ] = j++;
                }

                luint tNumberOfElements = tElementCounter( l, 1 );

                for( luint k=0; k<tNumberOfElements; ++k )
                {
                    tElements( tMap.find( tLagrangeElements( tCount++ ) ) )->put_on_refinement_queue();
                }

                // refine mesh
                aMesh->perform_refinement( tLagrangePattern );
            }

            aMesh->update_database();
        }

//-------------------------------------------------------------------------------

        std::string File::parralize_filename( const std::string & aPath )
        {
            // test if running in parallel mode
            if ( par_size() > 1 )
            {
                // get file extesion
                auto tFileExt = aPath.substr( aPath.find_last_of("."),
                                              aPath.length() );

                // get base path
                auto tBasePath = aPath.substr( 0, aPath.find_last_of(".") );

                // add proc number to path
                std::string aParallelPath = tBasePath + "_" +  std::to_string( par_rank() ) + tFileExt;
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
        Parameters * create_hmr_parameters_from_hdf5_file( const std::string & aPath )
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
