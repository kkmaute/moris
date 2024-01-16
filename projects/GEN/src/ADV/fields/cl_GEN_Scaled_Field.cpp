/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_GEN_Scaled_Field.cpp
 *
 */

#include "cl_GEN_Scaled_Field.hpp"

namespace moris::ge
{
    //--------------------------------------------------------------------------------------------------------------

    real Scaled_Field::get_field_value(
            uint                    aNodeIndex,
            const Matrix< DDRMat >& aCoordinates )
    {
        return mADVManager.get_variable( 0 ) * mField->get_field_value( aNodeIndex, aCoordinates );
    }

    //--------------------------------------------------------------------------------------------------------------

    real Scaled_Field::get_field_value( Derived_Node* aDerivedNode, const Node_Manager& aNodeManager )
    {
        return mADVManager.get_variable( 0 ) * mField->get_field_value( aDerivedNode, aNodeManager );
    }

    //--------------------------------------------------------------------------------------------------------------

    const Matrix<DDRMat>& Scaled_Field::get_dfield_dadvs(
            uint                    aNodeIndex,
            const Matrix< DDRMat >& aCoordinates )
    {
        mSensitivities = mADVManager.get_variable( 0 ) * mField->get_dfield_dadvs( aNodeIndex, aCoordinates );
        return mSensitivities;
    }
    //--------------------------------------------------------------------------------------------------------------

    void Scaled_Field::get_dfield_dadvs( Matrix< DDRMat >& aSensitivities, Derived_Node* aDerivedNode, const Node_Manager& aNodeManager )
    {
        mField->get_dfield_dadvs( aSensitivities, aDerivedNode, aNodeManager );
        aSensitivities = aSensitivities * mADVManager.get_variable( 0 );
    }


    //--------------------------------------------------------------------------------------------------------------

    Matrix<DDSMat> Scaled_Field::get_determining_adv_ids(
            uint                    aNodeIndex,
            const Matrix< DDRMat >& aCoordinates )
    {
        return mField->get_determining_adv_ids( aNodeIndex, aCoordinates );
    }

    //--------------------------------------------------------------------------------------------------------------

    void Scaled_Field::get_determining_adv_ids( Matrix< DDSMat >& aDeterminingADVIDs, Derived_Node* aDerivedNode, const Node_Manager& aNodeManager )
    {
        return mField->get_determining_adv_ids( aDeterminingADVIDs, aDerivedNode, aNodeManager );
    }

    //--------------------------------------------------------------------------------------------------------------

    void Scaled_Field::get_dfield_dcoordinates(
            uint                    aNodeIndex,
            const Matrix< DDRMat >& aCoordinates,
            Matrix< DDRMat >&       aSensitivities )
    {
        mField->get_dfield_dcoordinates( aNodeIndex, aCoordinates, aSensitivities );
        aSensitivities = aSensitivities * mADVManager.get_variable( 0 );
    }

    //--------------------------------------------------------------------------------------------------------------

    void Scaled_Field::set_dependencies( Cell< std::shared_ptr< Field > > aDependencyFields )
    {
        MORIS_ERROR( aDependencyFields.size() == 1, "A scaled field only depends on one field." );
        mField = aDependencyFields( 0 );
    }

    //--------------------------------------------------------------------------------------------------------------

    std::shared_ptr< mtk::Field > Scaled_Field::get_mtk_field()
    {
        return mField->get_mtk_field();
    }

    //--------------------------------------------------------------------------------------------------------------

}
