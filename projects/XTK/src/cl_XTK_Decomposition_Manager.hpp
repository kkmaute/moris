/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 * ------------------------------------------------------------------------------------
 *
 * cl_XTK_Decomposition_Manager.hpp
 *
 */

#pragma once

#include "cl_XTK_Node_Hierarchy_Interface.hpp"
#include "cl_XTK_Regular_Subdivision_Interface.hpp"
#include "cl_XTK_Integration_Mesh_Generator_Data.hpp"

namespace moris::xtk
{
    class Decomposition_Manager
    {
        std::unordered_map< mtk::Cell*, Vector< Decomposition_Algorithm_Type > > mDecompositionMap;    // Determines which subdivision methods have been applied to each element
        Vector< std::unique_ptr< Decomposition_Algorithm > >                     mAlgorithms;          // List of decomposition algorithms to be performed
        mtk::Mesh*                                                               mMesh;                // Mesh to be decomposed

      public:
        Decomposition_Manager(
                Vector< Subdivision_Method > const & aSubdivisionMethods,
                moris::Parameter_List&               aParameterList,
                mtk::Mesh*                           aMesh );

        /**
         * Performs each decomposition algorithm on the mesh, allowing for algorithms to blacklist certain elements from other algorithms
         */
        void perform();

      private:
        /**
         * Checks with each decomposition algorithm and updates the list of elements it can operate on
         */
        void update_eligible_elements();

        /**
         * Generates the data struct needed to perform the decomposition
         */
        Integration_Mesh_Generation_Data get_decomposition_struct( std::unique_ptr< Decomposition_Algorithm > const & aAlgorithm );
    };
}    // namespace moris::xtk