/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * ut_HMR_BSpline_Mesh_Private.cpp
 *
 */

#include <catch.hpp>

#include "cl_Communication_Manager.hpp" //COM/src
#include "cl_Communication_Tools.hpp" //COM/src
#include "moris_typedefs.hpp" //COR/src
#include "paths.hpp"
#include "cl_Matrix.hpp" //LINALG/src
#include "linalg_typedefs.hpp"
#include "op_equal_equal.hpp"
#include "fn_all_true.hpp"
#include "HDF5_Tools.hpp"

#define protected public
#define private   public
#include "cl_HMR_BSpline_Mesh_Base.hpp" //HMR/src
#include "cl_HMR_Background_Mesh.hpp" //HMR/src
#include "cl_HMR_Database.hpp"
#include "cl_HMR_Element.hpp" //HMR/src
#include "cl_HMR_Factory.hpp" //HMR/src
#include "cl_HMR_Parameters.hpp" //HMR/src
#include "cl_HMR.hpp" //HMR/src
#undef protected
#undef private

namespace moris::hmr
{

    TEST_CASE( "HMR_Bspline_Mesh_Private", "[moris],[mesh],[hmr],[BSplineMesh_private],[BsplineMesh]" )
    {
        //-------------------------------------------------------------------------------

        if ( par_size() == 1 || par_size() == 2 || par_size() == 4 )
        {
            //-------------------------------------------------------------------------------

            SECTION( "B-Spline Mesh test basis numbering" )
            {
                std::string tMorisRoot = get_base_moris_dir();

                // do this for first and second dimension
                for ( uint tDimension = 2; tDimension <= 3; ++tDimension )
                {
                    // do this for first second and third order
                    for ( uint tOrder = 1; tOrder < 3; ++tOrder )
                    {
                        for ( uint tMultigrid = 0; tMultigrid < 2; ++tMultigrid )
                        {
                            // The parameter object controls the behavior of HMR.
                            Parameters tParameters;

                            tParameters.set_multigrid( tMultigrid == 1 );

                            Vector< uint > tNumberOfElements( tDimension, 2 * tOrder );
                            tParameters.set_number_of_elements_per_dimension( tNumberOfElements );

                            Vector< real > tDomainOffset( tDimension, 0.0 );
                            tParameters.set_domain_offset( tDomainOffset );

                            tParameters.set_lagrange_orders( { 1 } );
                            tParameters.set_lagrange_patterns( { 0 } );

                            tParameters.set_bspline_orders( { tOrder } );
                            tParameters.set_bspline_patterns( { 0 } );

                            // set buffer
                            tParameters.set_refinement_buffer( tOrder );
                            tParameters.set_staircase_buffer( tOrder );

                            // create HMR Object
                            HMR tHMR( tParameters );

                            // refine the mesh three times
                            for ( uint k = 0; k < 3; ++k )
                            {
                                tHMR.flag_element( 0 );
                                tHMR.perform_refinement_based_on_working_pattern( 0 );
                            }

                            // finish mesh
                            tHMR.finalize();

                            // create matrix with IDs
                            Matrix< IdMat > tActiveBasis;
                            Matrix< IdMat > tRefinedBasis;

                            // reset counter
                            uint tCount = 0;

                            BSpline_Mesh_Base* tMesh = tHMR.mDatabase->mBSplineMeshes( 0 );

                            Vector< Basis* >& mActiveBasisOnProc = tMesh->mActiveBasisOnProc;

                            tActiveBasis.set_size( mActiveBasisOnProc.size(), 1 );

                            // loop over all active basis
                            for ( Basis* tBasis: mActiveBasisOnProc )
                            {
                                tActiveBasis( tCount++ ) = tBasis->get_id();
                            }

                            if ( tParameters.use_multigrid() )
                            {
                                Vector< Basis* >& mRefinedBasisOnProc = tMesh->mRefinedBasisOnProc;

                                tCount = 0;
                                tRefinedBasis.set_size( mRefinedBasisOnProc.size(), 1 );

                                for ( Basis* tBasis: mRefinedBasisOnProc )
                                {
                                    tRefinedBasis( tCount++ ) = tBasis->get_id();
                                }
                            }

                            std::string tFilePath = tMorisRoot + "/projects/HMR/test/data/hmr_bspline_ids_" + std::to_string(
                                    tDimension ) + std::to_string( tOrder ) + std::to_string( tMultigrid ) + ".hdf5";

                            /*
                            hid_t tFile = create_hdf5_file( tFilePath  );
                            herr_t tStatus = 0;
                            save_matrix_to_hdf5_file( tFile, "ActiveBasis", tActiveBasis, tStatus );
    
                            if( tParameters.use_multigrid() )
                            {
                                save_matrix_to_hdf5_file( tFile, "RefinedBasis", tRefinedBasis, tStatus );
                            }
                            close_hdf5_file( tFile ); */

                            Matrix< IdMat > tActiveBasisExpect;
                            Matrix< IdMat > tRefinedBasisExpect;

                            hid_t tFile = open_hdf5_file( tFilePath );
                            herr_t tStatus = 0;
                            load_matrix_from_hdf5_file( tFile, "ActiveBasis", tActiveBasisExpect, tStatus );
                            REQUIRE( all_true( tActiveBasis == tActiveBasisExpect ) );
                            if ( tParameters.use_multigrid() )
                            {
                                load_matrix_from_hdf5_file( tFile, "RefinedBasis", tRefinedBasisExpect,
                                                                   tStatus );
                                REQUIRE( all_true( tRefinedBasis == tRefinedBasisExpect ) );
                            }

                            close_hdf5_file( tFile );
                        }
                    }
                }
            }
        }
    }
}
