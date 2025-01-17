/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 * ------------------------------------------------------------------------------------
 *
 * cl_XTK_Decomposition_Manager.hpp
 *
 */

#include "cl_XTK_Decomposition_Manager.hpp"
#include "cl_XTK_Decomposition_Algorithm_Factory.hpp"

namespace moris::xtk
{
    //--------------------------------------------------------------------------------------------------

    Decomposition_Manager::Decomposition_Manager(
            const Vector< Subdivision_Method >& aSubdivisionMethods,
            moris::Parameter_List&              aParameterList,
            mtk::Mesh*                          aMesh )
            : mAlgorithms( aSubdivisionMethods.size() )
            , mMesh( aMesh )
    {
        // Get number of mtk cells
        uint tNumCells = mMesh->get_num_elems();
        
        // Create map for the mesh
        for( uint iCell = 0; iCell < tNumCells: iCell++ )
        {
            mDecompositionMap[ mMesh->get_mtk_cell( iCell ) ] = Vector< Decomposition_Algorithm_Type >();
        }
        
        // Loop through each of the subdivision methods
        for ( uint iSubdivision = 0; iSubdivision < aSubdivisionMethods.size(); iSubdivision++ )
        {
            // Create the decomposition algorithm
            mAlgorithms( iSubdivision ) = create_decomposition_algorithm( aSubdivisionMethods( iSubdivision ), aParameterList );
        }
    }

    //--------------------------------------------------------------------------------------------------

    void Decomposition_Manager::perform()
    {
        // Update the map for the subdivision methods
        this->update_eligible_elements();

        // Loop through each of the decomposition algorithms
        for ( auto tAlgorithm : mAlgorithms )
        {
            // Generate the decomposition struct
            Integration_Mesh_Generation_Data tDecompositionStruct = this->get_decomposition_struct( tAlgorithm );

            // Perform the decomposition algorithm
            tAlgorithm->perform(); // Brendan need args here
        }
    }

    //--------------------------------------------------------------------------------------------------

    void
    Decomposition_Manager::update_eligible_elements()
    {
        // Get the number of cells in the background mesh
        uint tNumCells = mMesh->get_num_elems();

        // Loop through each decomposition algorithm
        for ( auto tAlgorithm : mDecompositionMap )
        {
            // Loop through each element
            for( auto& iElementContext : mDecompositionMap )
            {
                // Check if this element should be decomposed by this algorithm
                if( tAlgorithm->is_eligible( iElementContext ) )
                {
                    // Add this algorithm to the list of algorithms that can operate on this element
                    iElementContext.second.push_back( tAlgorithm->get_algorithm_type() );
            }
        }
    }

    //--------------------------------------------------------------------------------------------------

    Integration_Mesh_Generation_Data
    get_decomposition_struct( std::unique_ptr< Decomposition_Algorithm > const & aAlgorithm )
    {
        // Loop through all the elements in the map and get the ones that have this algorithm
    }
    //--------------------------------------------------------------------------------------------------
}    // namespace moris::xtk