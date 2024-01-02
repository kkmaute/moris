/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_MSI_Dof_Manager_Test.cpp
 *
 */

#include "catch.hpp"
#include "moris_typedefs.hpp"
#include "cl_Map.hpp"
#include "cl_Matrix.hpp"
#include "linalg_typedefs.hpp"
#include "fn_equal_to.hpp"

#include "cl_Communication_Tools.hpp"
#include "cl_Communication_Manager.hpp"

#define protected public
#define private public
#include "cl_MSI_Adof.hpp"
#include "cl_MSI_Pdof_Host.hpp"
#include "cl_MSI_Equation_Object.hpp"
#include "cl_MSI_Dof_Manager.hpp"
#include "cl_MSI_Model_Solver_Interface.hpp"
#include "cl_MSI_Node_Proxy.hpp"
#include "cl_MSI_Solver_Interface_Proxy.hpp"
#include "cl_DLA_Solver_Interface.hpp"
#undef protected
#undef private

#include "cl_SOL_Warehouse.hpp"

#include "fn_PRM_MSI_Parameters.hpp"
namespace moris::MSI
{
    TEST_CASE( "Sparsity_Pattern", "[MSI],[Compute_Sparsity_Pattern]" )
    {
        if ( par_size() == 4 )
        {
            // initialize the solver interface
            MSI_Solver_Interface_Proxy* tSolverInterfaceProxy = new MSI_Solver_Interface_Proxy();
            MSI_Solver_Interface*       tSolverInterface      = tSolverInterfaceProxy;

            // get the adof map and owning processors
            Matrix< DDSMat > tAdofLocalGlobalOverLappingMap      = tSolverInterfaceProxy->get_my_local_global_overlapping_map();
            Matrix< DDSMat > tAdofLocalGlobalOverLappingMapOwner = tSolverInterfaceProxy->get_my_local_global_overlapping_map_owners();

            // create dummy variables to be set
            tSolverInterface->mDofMgn          = new Dof_Manager();
            tSolverInterface->mSolverWarehouse = std::make_shared< sol::SOL_Warehouse >();

            // set the comm table for the dof manager
            tSolverInterface->mDofMgn->mCommTable = tSolverInterfaceProxy->mCommTable;

            // create hypothetical adofs and assign id and owner based on the
            // matrix info
            for ( size_t i = 0; i < tAdofLocalGlobalOverLappingMap.numel(); i++ )
            {
                tSolverInterface->mDofMgn->mAdofList.push_back( new Adof() );

                tSolverInterface->mDofMgn->mAdofList( i )->set_adof_id( tAdofLocalGlobalOverLappingMap( i ) );
                tSolverInterface->mDofMgn->mAdofList( i )->set_adof_owning_processor( tAdofLocalGlobalOverLappingMapOwner( i ) );
            }

            // call the sparsity pattern
            tSolverInterface->compute_sparsity_pattern();

            // cast to solver interface such that we can access member data that was
            // computed in the compute_sparsity_pattern
            moris::Solver_Interface* tDLASolInterface = static_cast< moris::Solver_Interface* >( tSolverInterface );

            // print the data
            void ( *tPrintFunc )( Vector< moris_id > const &, std::string ) = print_as_row_vector;
            print_cell< moris_id >( tDLASolInterface->mNonZeroDigonal, "mNonZeroDigonal", tPrintFunc );
            print_cell< moris_id >( tDLASolInterface->mNonZeroOffDigonal, "mNonZeroOffDigonal", tPrintFunc );

            // set what data should be for each processor
            Vector< moris_id > tNonZeroDigonal;
            Vector< moris_id > tNonZeroOffDigonal;
            switch ( par_rank() )
            {
                case 0:
                    tNonZeroDigonal    = { 4, 4, 4, 4, 4, 4, 4, 4 };
                    tNonZeroOffDigonal = { 0, 0, 2, 2, 2, 6, 0, 0 };
                    break;
                case 1:
                    tNonZeroDigonal    = { 2, 2, 2, 2 };
                    tNonZeroOffDigonal = { 2, 2, 6, 5 };
                    break;
                case 2:
                    tNonZeroDigonal    = { 2, 2 };
                    tNonZeroOffDigonal = { 5, 5 };
                    break;
                case 3:
                    tNonZeroDigonal    = { 2, 2, 2, 2 };
                    tNonZeroOffDigonal = { 5, 6, 0, 0 };
                    break;

                default:
                    break;
            }

            // check the the data
            REQUIRE( std::equal( tNonZeroDigonal.begin(), tNonZeroDigonal.end(), tDLASolInterface->mNonZeroDigonal.begin() ) );
            REQUIRE( std::equal( tNonZeroOffDigonal.begin(), tNonZeroOffDigonal.end(), tDLASolInterface->mNonZeroOffDigonal.begin() ) );
        }
    }

}    // namespace moris::MSI
