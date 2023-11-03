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

    Field::Field(
            Matrix< DDRMat >& aADVs,
            Matrix< DDUMat >  aFieldVariableIndices,
            Matrix< DDUMat >  aADVIndices,
            Matrix< DDRMat >  aConstants,
            std::string       aName )
            : mADVManager( aADVs, aFieldVariableIndices, aADVIndices, aConstants )
            , mSensitivities( 1, aFieldVariableIndices.length() + aConstants.length() )
            , mName( std::move( aName ) )
    {
    }

    //--------------------------------------------------------------------------------------------------------------

    Field::Field( Matrix< DDRMat > aConstants )
            : mADVManager( aConstants )
            , mSensitivities( 1, aConstants.length(), 0.0 )
    {
    }
    
    //--------------------------------------------------------------------------------------------------------------
    
    Field::Field( const Matrix< DDSMat >& aSharedADVIds )
            : mADVManager( aSharedADVIds )
            , mSensitivities( 1, aSharedADVIds.length() )
    {
    }

    //--------------------------------------------------------------------------------------------------------------

    Field::Field( const Field& aCopy,
            const Cell< uint >& aReplaceVariables,
            const Cell< real >& aNewConstants )
            : mADVManager( aCopy.mADVManager, aReplaceVariables, aNewConstants )
            , mSensitivities( aCopy.mSensitivities )
    {
    }

    //--------------------------------------------------------------------------------------------------------------

    std::shared_ptr< Field > Field::copy(
            const Cell< uint >& aReplaceVariables,
            const Cell< real >& aNewConstants )
    {
        return nullptr;
    }

    //--------------------------------------------------------------------------------------------------------------

    void Field::import_advs( sol::Dist_Vector* aOwnedADVs )
    {
        mADVManager.import_advs( aOwnedADVs );
    }

    //--------------------------------------------------------------------------------------------------------------

    uint Field::get_number_of_reference_coordinates()
    {
        return 0;
    }

    //--------------------------------------------------------------------------------------------------------------

    Matrix< DDSMat > Field::get_determining_adv_ids(
            uint                    aNodeIndex,
            const Matrix< DDRMat >& aCoordinates )
    {
        return mADVManager.get_determining_adv_ids();
    }

    //--------------------------------------------------------------------------------------------------------------

    bool Field::has_advs()
    {
        return mADVManager.has_advs();
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

    std::string Field::get_name()
    {
        return mName;
    }

    //--------------------------------------------------------------------------------------------------------------
}
