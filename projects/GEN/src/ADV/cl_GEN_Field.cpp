/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_GEN_Field.cpp
 *
 */

#include "cl_GEN_Field.hpp"

namespace moris::ge
{
    //--------------------------------------------------------------------------------------------------------------

    Field::Field( Matrix< DDRMat > aConstants )
            : ADV_Manager( aConstants )
    {
    }
    
    //--------------------------------------------------------------------------------------------------------------
    
    Field::Field(
            const Matrix< DDUMat >& aFieldVariableIndices,
            const Matrix< DDSMat >& aSharedADVIds )
            : ADV_Manager( aFieldVariableIndices, aSharedADVIds )
    {
    }

    //--------------------------------------------------------------------------------------------------------------

    void Field::set_dependencies( Cell< std::shared_ptr< Field > > aDependencyFields )
    {
    }

    //--------------------------------------------------------------------------------------------------------------

    void
    Field::add_child_node( uint aNodeIndex, std::shared_ptr< Child_Node > aChildNode )
    {
    }

    //--------------------------------------------------------------------------------------------------------------

    void
    Field::add_nodal_data( mtk::Interpolation_Mesh* aMesh )
    {
    }

    //--------------------------------------------------------------------------------------------------------------

    void
    Field::reset_nodal_data()
    {
    }

    //--------------------------------------------------------------------------------------------------------------

    void
    Field::set_num_original_nodes( uint aNumOriginalNodes )
    {
        mNumOriginalNodes = aNumOriginalNodes;
    }

    //--------------------------------------------------------------------------------------------------------------
}
