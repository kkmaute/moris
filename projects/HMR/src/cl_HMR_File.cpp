/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_HMR_File.cpp
 *
 */

#include "cl_HMR_File.hpp" //HMR/src

#include "cl_HMR_Factory.hpp" //HMR/src
#include "HDF5_Tools.hpp" //HMR/src
#include "cl_Map.hpp"

namespace moris::hmr
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
        save_matrix_to_hdf5_file( mFileID,
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

        load_scalar_from_hdf5_file( mFileID,
                "InitialBSplineRefinement",
                tValUint,
                mStatus );

        aParameters->set_initial_refinement( {{tValUint}} );

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
                "OrderToLagrangeMeshList",
                tMatUint,
                mStatus );

        aParameters->set_lagrange_orders( tMatUint );

        // load Lagrange mesh associations
        load_matrix_from_hdf5_file( mFileID,
                "PatternToLagrangeMeshList",
                tMatUint,
                mStatus );

        aParameters->set_lagrange_patterns( tMatUint );

        // load orders of meshes
        load_matrix_from_hdf5_file( mFileID,
                "OrderToBspMeshList",
                tMatUint,
                mStatus );

        aParameters->set_bspline_orders( tMatUint );

        // load bspline mesh associations
        load_matrix_from_hdf5_file( mFileID,
                "PatternToBspMeshList",
                tMatUint,
                mStatus );

        aParameters->set_bspline_patterns( tMatUint );

        // set lagrange to bpline mesh dependecies. since we read one lag mesh from file all bsplines belong to this mesh
        moris::Vector< Matrix< DDSMat > >tMatBspToLag( 1 );
        tMatBspToLag( 0 ).set_size(tMatUint.numel(), 1);

        for( uint Ik = 0; Ik < tMatUint.numel(); Ik++ )
        {
            tMatBspToLag( 0 )( Ik ) = Ik;
        }

        aParameters->set_lagrange_to_bspline_mesh( tMatBspToLag );

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

    void File::save_refinement_pattern( Lagrange_Mesh_Base * aLagrangeMesh )
    {
        Background_Mesh_Base * aBackgroundMesh = aLagrangeMesh->get_background_mesh();
        // step 1: count how many elements need are refined on each level
        uint tMaxLevel = aBackgroundMesh->get_max_level();

        uint tNumBSplineMeshes = aLagrangeMesh->get_number_of_bspline_meshes();

        // Initialize Cells and Mat for Pattern and order list. Both to do the unique on Cell and to write the Mat to hdf5
        moris::Vector< moris::uint > tPatternList( 1 + tNumBSplineMeshes, MORIS_UINT_MAX );
        moris::Vector< moris::uint > tOrderList( 1 + tNumBSplineMeshes, MORIS_UINT_MAX );

        moris::Matrix< DDUMat > tPatternLagMat( 1 , 1, MORIS_UINT_MAX );
        moris::Matrix< DDUMat > tPatternBspMat( tNumBSplineMeshes, 1, MORIS_UINT_MAX );
        moris::Matrix< DDUMat > tOrderLagMat( 1 , 1, MORIS_UINT_MAX );
        moris::Matrix< DDUMat > tOrderBspMat( tNumBSplineMeshes, 1, MORIS_UINT_MAX );

        // Fill Cells and Mats with pattern index and order
        tPatternList( 0 ) = aLagrangeMesh->get_activation_pattern();
        tPatternLagMat( 0 ) =aLagrangeMesh->get_activation_pattern();
        tOrderList( 0 ) = aLagrangeMesh->get_order();
        tOrderLagMat( 0 ) = aLagrangeMesh->get_order();

        for( uint Ik = 0; Ik < tNumBSplineMeshes; Ik++ )
        {
            tPatternList( Ik+1 ) = aLagrangeMesh->get_bspline_mesh( Ik )->get_activation_pattern();
            tPatternBspMat( Ik ) = aLagrangeMesh->get_bspline_mesh( Ik )->get_activation_pattern();
            tOrderList( Ik+1 )   = aLagrangeMesh->get_bspline_mesh( Ik )->get_min_order();
            tOrderBspMat( Ik )   = aLagrangeMesh->get_bspline_mesh( Ik )->get_min_order();
        }

        MORIS_ERROR( tPatternLagMat.max() != MORIS_UINT_MAX, "File::save_refinement_pattern(); the pattern list is not initialized correctly and has a MORIS_UINT_MAX entry" );
        MORIS_ERROR( tPatternBspMat.max() != MORIS_UINT_MAX, "File::save_refinement_pattern(); the order list is not initialized correctly and has a MORIS_UINT_MAX entry" );

        save_matrix_to_hdf5_file( mFileID,
                "PatternToLagrangeMeshList",
                tPatternLagMat,
                mStatus );

        save_matrix_to_hdf5_file( mFileID,
                "PatternToBspMeshList",
                tPatternBspMat,
                mStatus );

        save_matrix_to_hdf5_file( mFileID,
                "OrderToLagrangeMeshList",
                tOrderLagMat,
                mStatus );

        save_matrix_to_hdf5_file( mFileID,
                "OrderToBspMeshList",
                tOrderBspMat,
                mStatus );

        // Sort this created list
        std::sort( ( tPatternList.data() ).data(), ( tPatternList.data() ).data() + tPatternList.size() );

        // use std::unique and std::distance to create list containing all used dof types. This list is unique
        auto last = std::unique( ( tPatternList.data() ).data(), ( tPatternList.data() ).data() + tPatternList.size() );
        auto pos  = std::distance( ( tPatternList.data() ).data(), last );

        tPatternList.resize( pos );
        uint tNumUniquePattern = tPatternList.size();
        moris::Matrix< DDUMat >tPatternListUniqueMat( tNumUniquePattern, 1, MORIS_UINT_MAX );

        // Copy unique list in Mat
        for( uint Ik = 0; Ik < tPatternList.size(); Ik++ )
        {
            tPatternListUniqueMat( Ik ) = tPatternList( Ik );
        }

        save_matrix_to_hdf5_file( mFileID,
                "PatternInd",
                tPatternListUniqueMat,
                mStatus );

        // element counter
        Matrix< DDLUMat > tElementCounter ( tMaxLevel+1, tNumUniquePattern, 0 );

        // collect all elements that are flagged for refinement
        for( uint l = 0; l < tMaxLevel; ++l )
        {
            // cell which contains elements
            Vector< Background_Element_Base* > tElements;

            // collect elements from this level
            aBackgroundMesh->collect_elements_on_level_within_proc_domain( l, tElements );

            // loop over all elements
            for( auto tElement : tElements )
            {
                for( uint Ik = 0; Ik < tNumUniquePattern; ++Ik )
                {
                    // test if B-Spline Element is refined
                    if( tElement->is_refined( tPatternList( Ik ) ) )
                    {
                        // increment counter
                        ++tElementCounter ( l, Ik );
                    }
                }
            }
        }

        moris::Vector< Matrix< DDLUMat > > tPatternElement( tNumUniquePattern );
        moris::Vector< hsize_t > tElementPerPatternCount( tNumUniquePattern, 0 );

        for( uint Ik = 0; Ik < tNumUniquePattern; ++Ik )
        {
            tPatternElement( Ik ).set_size( sum( tElementCounter.get_column( Ik ) ) , 1 );
        }

        for( uint l = 0; l < tMaxLevel; ++l )
        {
            // cell which contains elements
            Vector< Background_Element_Base* > tElements;

            // collect elements from this level
            aBackgroundMesh->collect_elements_on_level_within_proc_domain( l, tElements );

            // loop over all elements
            for( Background_Element_Base * tElement : tElements )
            {
                for( uint Ik = 0; Ik < tNumUniquePattern; ++Ik )
                {
                    // test if element is refined
                    if( tElement->is_refined( tPatternList( Ik ) ) )
                    {
                        tPatternElement( Ik )( tElementPerPatternCount( Ik )++ ) = tElement->get_hmr_id();
                    }
                }
            }
        }

        save_matrix_to_hdf5_file( mFileID,
                "ElementCounter",
                tElementCounter,
                mStatus );

        for( uint Ik = 0; Ik < tNumUniquePattern; ++Ik )
        {
            std::string tSubsectionStr = "Pattern_" + std::to_string( tPatternList( Ik ) ) + "_Elements";

            save_matrix_to_hdf5_file( mFileID,
                    tSubsectionStr,
                    tPatternElement( Ik ),
                    mStatus );
        }
    }

    //------------------------------------------------------------------------------

    void File::save_refinement_pattern( Background_Mesh_Base          * aBackgroundMesh,
                                        const moris::Matrix< DDUMat > & tPatternToSave  )
    {
        // step 1: count how many elements need are refined on each level
        uint tMaxLevel = aBackgroundMesh->get_max_level();

        uint tNumPattern = tPatternToSave.numel();

        save_matrix_to_hdf5_file( mFileID,
                "PatternInd",
                tPatternToSave,
                mStatus );

        // element counter
        Matrix< DDLUMat > tElementCounter ( tMaxLevel+1, tNumPattern, 0 );

        // collect all elements that are flagged for refinement
        for( uint l = 0; l < tMaxLevel; ++l )
        {
            // cell which contains elements
            Vector< Background_Element_Base* > tElements;

            // collect elements from this level
            aBackgroundMesh->collect_elements_on_level_within_proc_domain( l, tElements );

            // loop over all elements
            for( auto tElement : tElements )
            {
                for( uint Ik = 0; Ik < tNumPattern; ++Ik )
                {
                    // test if B-Spline Element is refined
                    if( tElement->is_refined( tPatternToSave( Ik ) ) )
                    {
                        // increment counter
                        ++tElementCounter ( l, Ik );
                    }
                }
            }
        }

        moris::Vector< Matrix< DDLUMat > > tPatternElement( tNumPattern );
        moris::Vector< hsize_t > tElementPerPatternCount( tNumPattern, 0 );

        for( uint Ik = 0; Ik < tNumPattern; ++Ik )
        {
            tPatternElement( Ik ).set_size( sum( tElementCounter.get_column( Ik ) ) , 1 );
        }

        for( uint l = 0; l < tMaxLevel; ++l )
        {
            // cell which contains elements
            Vector< Background_Element_Base* > tElements;

            // collect elements from this level
            aBackgroundMesh->collect_elements_on_level_within_proc_domain( l, tElements );

            // loop over all elements
            for( Background_Element_Base * tElement : tElements )
            {
                for( uint Ik = 0; Ik < tNumPattern; ++Ik )
                {
                    // test if element is refined
                    if( tElement->is_refined( tPatternToSave( Ik ) ) )
                    {
                        tPatternElement( Ik )( tElementPerPatternCount( Ik )++ ) = tElement->get_hmr_id();
                    }
                }
            }
        }

        save_matrix_to_hdf5_file( mFileID,
                "ElementCounter",
                tElementCounter,
                mStatus );

        for( uint Ik = 0; Ik < tNumPattern; ++Ik )
        {
            std::string tSubsectionStr = "Pattern_" + std::to_string( tPatternToSave( Ik ) ) + "_Elements";

            save_matrix_to_hdf5_file( mFileID,
                    tSubsectionStr,
                    tPatternElement( Ik ),
                    mStatus );
        }
    }

    //------------------------------------------------------------------------------

    void File::save_refinement_pattern(
            Background_Mesh_Base             * aBackgroundMesh,
            const moris::Matrix< DDUMat >    & tPatternToSave,
            Matrix< DDLUMat >                & aElementCounterPerLevelAndPattern,
            moris::Vector< Matrix< DDLUMat > > & aElementPerPattern)
    {
        uint tMaxLevel = aBackgroundMesh->get_max_level();

        uint tNumPattern = tPatternToSave.numel();

        moris::Matrix< DDUMat > tPatterns(tNumPattern, 1 );

        // element counter
        aElementCounterPerLevelAndPattern.set_size( tMaxLevel+1, tNumPattern, 0 );

        // collect all elements that are flagged for refinement
        for( uint l = 0; l < tMaxLevel; ++l )
        {
            // cell which contains elements
            Vector< Background_Element_Base* > tElements;

            // collect elements from this level
            aBackgroundMesh->collect_elements_on_level_within_proc_domain( l, tElements );

            // loop over all elements
            for( auto tElement : tElements )
            {
                for( uint Ik = 0; Ik < tNumPattern; ++Ik )
                {
                    // test if B-Spline Element is refined
                    if( tElement->is_refined( tPatternToSave( Ik ) ) )
                    {
                        // increment counter
                        ++aElementCounterPerLevelAndPattern ( l, Ik );
                    }
                }
            }
        }

        aElementPerPattern.resize( tNumPattern );
        moris::Vector< luint > tElementPerPatternCount( tNumPattern, 0 );

        for( uint Ik = 0; Ik < tNumPattern; ++Ik )
        {
            aElementPerPattern( Ik ).set_size( sum( aElementCounterPerLevelAndPattern.get_column( Ik ) ) , 1 );
        }

        for( uint l = 0; l < tMaxLevel; ++l )
        {
            // cell which contains elements
            Vector< Background_Element_Base* > tElements;

            // collect elements from this level
            aBackgroundMesh->collect_elements_on_level_within_proc_domain( l, tElements );

            // loop over all elements
            for( Background_Element_Base * tElement : tElements )
            {
                for( uint Ik = 0; Ik < tNumPattern; ++Ik )
                {
                    // test if element is refined
                    if( tElement->is_refined( tPatternToSave( Ik ) ) )
                    {
                        aElementPerPattern( Ik )( tElementPerPatternCount( Ik )++ ) = tElement->get_hmr_id();
                    }
                }
            }
        }
    }

    //------------------------------------------------------------------------------

    void File::save_refinement_pattern(
            Lagrange_Mesh_Base               * aLagrangeMesh,
            const uint                         aDiscretizationMeshIndex,
            Matrix< DDLUMat >                & aElementCounterPerLevelAndPattern,
            moris::Vector< Matrix< DDLUMat > > & aElementPerPattern )
    {
        Background_Mesh_Base * aBackgroundMesh = aLagrangeMesh->get_background_mesh();
        // step 1: count how many elements need are refined on each level
        uint tMaxLevel = aBackgroundMesh->get_max_level();

        moris::Matrix< DDUMat > tPatterns( 2, 1 );

        uint tSourceLagrangePattern = aLagrangeMesh->get_activation_pattern();
        uint tSourceBSplinePattern = aLagrangeMesh ->get_bspline_pattern( aDiscretizationMeshIndex );

        moris::Matrix< DDUMat > tPatternList;

        if( tSourceLagrangePattern == tSourceBSplinePattern )
        {
            tPatternList = { {tSourceLagrangePattern} };
        }
        else
        {
            tPatternList = { {tSourceLagrangePattern}, { tSourceBSplinePattern } };
        }

        uint tNumPattern = tPatternList.numel();

        // element counter
        aElementCounterPerLevelAndPattern.set_size( tMaxLevel+1, tNumPattern, 0 );

        // collect all elements that are flagged for refinement
        for( uint l = 0; l < tMaxLevel; ++l )
        {
            // cell which contains elements
            Vector< Background_Element_Base* > tElements;

            // collect elements from this level
            aBackgroundMesh->collect_elements_on_level_within_proc_domain( l, tElements );

            // loop over all elements
            for( auto tElement : tElements )
            {
                for( uint Ik = 0; Ik < tNumPattern; ++Ik )
                {
                    // test if B-Spline Element is refined
                    if( tElement->is_refined( tPatternList( Ik ) ) )
                    {
                        // increment counter
                        ++aElementCounterPerLevelAndPattern ( l, Ik );
                    }
                }
            }
        }

        aElementPerPattern.resize( tNumPattern );
        moris::Vector< luint > tElementPerPatternCount( tNumPattern, 0 );

        for( uint Ik = 0; Ik < tNumPattern; ++Ik )
        {
            aElementPerPattern( Ik ).set_size( sum( aElementCounterPerLevelAndPattern.get_column( Ik ) ) , 1 );
        }

        for( uint l = 0; l < tMaxLevel; ++l )
        {
            // cell which contains elements
            Vector< Background_Element_Base* > tElements;

            // collect elements from this level
            aBackgroundMesh->collect_elements_on_level_within_proc_domain( l, tElements );

            // loop over all elements
            for( Background_Element_Base * tElement : tElements )
            {
                for( uint Ik = 0; Ik < tNumPattern; ++Ik )
                {
                    // test if element is refined
                    if( tElement->is_refined( tPatternList( Ik ) ) )
                    {
                        aElementPerPattern( Ik )( tElementPerPatternCount( Ik )++ ) = tElement->get_hmr_id();
                    }
                }
            }
        }
    }

    //------------------------------------------------------------------------------

    void File::load_refinement_pattern(
            Background_Mesh_Base             * aMesh,
            Matrix< DDLUMat >                & aElementCounterPerLevelAndPattern,
            moris::Vector< Matrix< DDLUMat > > & aElementPerPattern)
    {
        uint tNumPattern = aElementPerPattern.size();

        Matrix< DDUMat > tPatternListUniqueMat( tNumPattern, 1 );
        tPatternListUniqueMat( 0 ) = 5;
        if( tNumPattern == 2 )
        {
            tPatternListUniqueMat( 1 ) = 6;
        }

        // get number of levels
        uint tNumberOfLevels = aElementCounterPerLevelAndPattern.n_rows();

        for(uint Ik = 0; Ik<tNumPattern; Ik++)
        {
            // reset counter
            luint tCount = 0;

            // select B-Spline pattern
            aMesh->set_activation_pattern( tPatternListUniqueMat( Ik ) );

            // loop over all levels
            for( uint l=0; l<tNumberOfLevels; ++l )
            {
                // cell which contains elements
                Vector< Background_Element_Base* > tElements;

                // collect elements from this level
                aMesh->collect_elements_on_level_within_proc_domain( l, tElements );

                // create a map with ids
                map< moris_id, luint > tMap;

                luint j = 0;
                for( Background_Element_Base* tElement : tElements )
                {
                    tMap[ tElement->get_hmr_id() ] = j++;
                }

                luint tNumberOfElements = aElementCounterPerLevelAndPattern( l, Ik );

                for( luint k=0; k<tNumberOfElements; ++k )
                {
                    tElements( tMap.find( aElementPerPattern( Ik )( tCount++ ) ) )->put_on_refinement_queue();
                }

                // refine mesh
                aMesh->perform_refinement( tPatternListUniqueMat( Ik ) );
            }
        }

        aMesh->update_database();
    }

    //------------------------------------------------------------------------------

    void File::load_refinement_pattern(
            Background_Mesh_Base * aMesh,
            const bool             aMode )
    {
        Matrix< DDUMat > tPatternListUniqueMat;
        load_matrix_from_hdf5_file( mFileID,
                "PatternInd",
                tPatternListUniqueMat,
                mStatus );

        uint tMaxPatternInd = tPatternListUniqueMat.max();
        Matrix< DDUMat > tPatternMap( tMaxPatternInd +1 , 1, MORIS_UINT_MAX );

        uint tNumUniquePattern = tPatternListUniqueMat.numel();

        for(uint Ik = 0; Ik < tNumUniquePattern; Ik++)
        {
            tPatternMap( tPatternListUniqueMat( Ik ) ) = Ik;
        }

        // matrix containing counter
        Matrix< DDLUMat > tElementCounter;
        load_matrix_from_hdf5_file( mFileID,
                "ElementCounter",
                tElementCounter,
                mStatus );

        moris::Vector< Matrix< DDLUMat > > tPatternElement( tNumUniquePattern );

        for(uint Ik = 0; Ik<tNumUniquePattern; Ik++)
        {
            std::string tSubsectionStr = "Pattern_" + std::to_string( tPatternListUniqueMat( Ik ) ) + "_Elements";

            // allocate pattern
            load_matrix_from_hdf5_file( mFileID,
                    tSubsectionStr,
                    tPatternElement( Ik ),
                    mStatus );
        }

        // get number of levels
        uint tNumberOfLevels = tElementCounter.n_rows();

        for(uint Ik = 0; Ik<tNumUniquePattern; Ik++)
        {
            // reset counter
            luint tCount = 0;

            // select B-Spline pattern
            aMesh->set_activation_pattern( tPatternListUniqueMat( Ik ) );

            // loop over all levels
            for( uint l=0; l<tNumberOfLevels; ++l )
            {
                // cell which contains elements
                Vector< Background_Element_Base* > tElements;

                // collect elements from this level
                aMesh->collect_elements_on_level_within_proc_domain( l, tElements );

                // create a map with ids
                map< moris_id, luint > tMap;

                luint j = 0;
                for( Background_Element_Base* tElement : tElements )
                {
                    tMap[ tElement->get_hmr_id() ] = j++;
                }

                luint tNumberOfElements = tElementCounter( l, Ik );

                for( luint k=0; k<tNumberOfElements; ++k )
                {
                    tElements( tMap.find( tPatternElement( Ik )( tCount++ ) ) )->put_on_refinement_queue();
                }

                // refine mesh
                aMesh->perform_refinement( tPatternListUniqueMat( Ik ) );
            }
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

} /* namespace moris */
