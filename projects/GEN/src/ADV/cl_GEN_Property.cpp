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

#include <utility>
#include "cl_MTK_Integration_Mesh.hpp"

namespace moris::gen
{

    //--------------------------------------------------------------------------------------------------------------

    Property_Parameters::Property_Parameters( const Parameter_List& aParameterList )
            : Field_Parameters( aParameterList )
            , Design_Parameters( aParameterList )
            , mPDVType( get_pdv_type_map()[ aParameterList.get< std::string >( "pdv_type" ) ] )
            , mInterpolationPDV( aParameterList.get< std::string >( "pdv_mesh_type" ) == "interpolation" )
            , mPDVMeshSetIndices( aParameterList.get_cell< uint >( "pdv_mesh_set_indices" ) )
            , mPDVMeshSetNames( aParameterList.get_cell< std::string >( "pdv_mesh_set_names" ) )
    {
    }

    //--------------------------------------------------------------------------------------------------------------

    Property::Property(
            std::shared_ptr< Field >   aField,
            const Property_Parameters& aParameters,
            Node_Manager&              aNodeManager )
            : Design_Field( std::move( aField ), aParameters, aNodeManager )
            , Design( aParameters )
            , mParameters( aParameters )
    {
    }

    //--------------------------------------------------------------------------------------------------------------

    void Property::set_node_manager( Node_Manager& aNodeManager )
    {
        mNodeManager = &aNodeManager;
    }

    //--------------------------------------------------------------------------------------------------------------

    void Property::update_dependencies( const Vector< std::shared_ptr< Design > >& aAllUpdatedDesigns )
    {
        // Get fields from designs
        Vector< std::shared_ptr< Field > > tUpdatedFields( aAllUpdatedDesigns.size() );
        for ( uint iFieldIndex = 0; iFieldIndex < tUpdatedFields.size(); iFieldIndex++ )
        {
            tUpdatedFields( iFieldIndex ) = aAllUpdatedDesigns( iFieldIndex )->get_field();
        }

        // Update fields
        Design_Field::update_dependencies( tUpdatedFields );
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

    Vector< uint > Property::get_pdv_mesh_set_indices( mtk::Integration_Mesh* aMesh )
    {
        // Number of mesh set names
        uint tNumMeshSetNames = mParameters.mPDVMeshSetNames.size();

        // Create new cell of indices from names
        Vector< uint > tMeshSetIndicesFromNames( tNumMeshSetNames );

        // Set for each property index the list of mesh set indices
        for ( uint iSetNameIndex = 0; iSetNameIndex < tNumMeshSetNames; iSetNameIndex++ )
        {
            tMeshSetIndicesFromNames( iSetNameIndex ) =
                    aMesh->get_set_index_by_name( mParameters.mPDVMeshSetNames( iSetNameIndex ) );
        }

        // Append given indices
        tMeshSetIndicesFromNames.append( mParameters.mPDVMeshSetIndices );

        return tMeshSetIndicesFromNames;
    }

    //--------------------------------------------------------------------------------------------------------------

    void Property::get_design_info(
            uint                    aNodeIndex,
            const Matrix< DDRMat >& aCoordinates,
            Vector< real >&         aOutputDesignInfo )
    {
        aOutputDesignInfo.resize( 1 );
        aOutputDesignInfo( 0 ) = Design_Field::get_field_value( aNodeIndex, aCoordinates );
    }

    std::string
    Property::get_name()
    {
        return Design_Field::get_name();
    }

    bool Property::intended_discretization()
    {
        return ( mParameters.mDiscretizationIndex >= 0 );
    }

    //--------------------------------------------------------------------------------------------------------------

    moris_index
    Property::get_discretization_mesh_index()
    {
        MORIS_ASSERT( mParameters.mDiscretizationIndex >= 0,
                "A discretization is not intended for this field. Check this with intended_discretization() first." );

        return mParameters.mDiscretizationIndex;
    }

    //--------------------------------------------------------------------------------------------------------------

    real
    Property::get_discretization_lower_bound()
    {
        return mParameters.mDiscretizationLowerBound;
    }

    //--------------------------------------------------------------------------------------------------------------

    real
    Property::get_discretization_upper_bound()
    {
        return mParameters.mDiscretizationUpperBound;
    }

    //--------------------------------------------------------------------------------------------------------------

}    // namespace moris::gen
