/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_GEN_Property.cpp
 *
 */

#include "cl_GEN_Property.hpp"

namespace moris::ge
{

    //--------------------------------------------------------------------------------------------------------------

    Property_Parameters::Property_Parameters( const ParameterList& aParameterList )
            : Field_Parameters( aParameterList )
            , mDependencyNames( aParameterList.get_cell< std::string >( "dependencies" ) )
            , mPDVType( get_pdv_type_map()[ aParameterList.get< std::string >( "pdv_type" ) ] )
            , mInterpolationPDV( aParameterList.get< std::string >( "pdv_mesh_type" ) == "interpolation" )
            , mPDVMeshSetIndices( aParameterList.get_cell< uint >( "pdv_mesh_set_indices" ) )
            , mPDVMeshSetNames( aParameterList.get_cell< std::string >( "pdv_mesh_set_names" ) )
    {
    }

    //--------------------------------------------------------------------------------------------------------------

    Property::Property(
            std::shared_ptr< Field > aField,
            Property_Parameters      aParameters )
            : Design_Field( aField, aParameters )
            , mParameters( aParameters )
    {
    }

    //--------------------------------------------------------------------------------------------------------------

    void Property::update_dependencies( Cell< std::shared_ptr< Design_Field > > aAllUpdatedFields )
    {
        // Set up dependency fields
        uint tNumDependencies = mParameters.mDependencyNames.size();
        Cell< std::shared_ptr< Field > > tDependencyFields( tNumDependencies );

        // Grab dependencies
        for ( uint tDependencyIndex = 0; tDependencyIndex < tNumDependencies; tDependencyIndex++ )
        {
            for ( uint tFieldIndex = 0; tFieldIndex < aAllUpdatedFields.size(); tFieldIndex++ )
            {
                if ( aAllUpdatedFields( tFieldIndex )->get_name() == mParameters.mDependencyNames( tDependencyIndex ) )
                {
                    tDependencyFields( tDependencyIndex ) = aAllUpdatedFields( tFieldIndex )->get_field();
                }
            }
        }

        // Set dependencies
        this->get_field()->set_dependencies( tDependencyFields );
    }

    //--------------------------------------------------------------------------------------------------------------

    PDV_Type Property::get_pdv_type()
    {
        return mParameters.mPDVType;
    }

    //--------------------------------------------------------------------------------------------------------------

    bool Property::is_interpolation_pdv()
    {
        return mParameters.mInterpolationPDV;
    }

    //--------------------------------------------------------------------------------------------------------------

    Cell< uint > Property::get_pdv_mesh_set_indices()
    {
        return mParameters.mPDVMeshSetIndices;
    }

    //--------------------------------------------------------------------------------------------------------------

    Cell< std::string > Property::get_pdv_mesh_set_names()
    {
        return mParameters.mPDVMeshSetNames;
    }

    //--------------------------------------------------------------------------------------------------------------

}
