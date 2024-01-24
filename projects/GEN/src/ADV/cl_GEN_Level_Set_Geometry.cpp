/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_GEN_Level_Set_Geometry.cpp
 *
 */

#include "cl_GEN_Level_Set_Geometry.hpp"
#include "cl_GEN_Intersection_Node_Linear.hpp"
#include "cl_GEN_Intersection_Node_Bilinear.hpp"
#include "cl_GEN_Derived_Node.hpp"
#include "cl_GEN_Basis_Node.hpp"
#include "fn_dot.hpp"

namespace moris::ge
{

    //--------------------------------------------------------------------------------------------------------------

    Level_Set_Parameters::Level_Set_Parameters( const ParameterList& aParameterList )
            : Field_Parameters( aParameterList )
            , Design_Parameters( aParameterList )
            , mIsocontourThreshold( aParameterList.get< real >( "isocontour_threshold" ) )
            , mIsocontourTolerance( aParameterList.get< real >( "isocontour_tolerance" ) )
            , mIntersectionTolerance( aParameterList.get< real >( "intersection_tolerance" ) )
    {
    }

    //--------------------------------------------------------------------------------------------------------------

    Level_Set_Geometry::Level_Set_Geometry(
            std::shared_ptr< Field > aField,
            Level_Set_Parameters     aParameters,
            Node_Manager&            aNodeManager )
            : Geometry( aParameters )
            , Design_Field( aField, aParameters, aNodeManager )
            , mParameters( aParameters )
    {
    }

    //--------------------------------------------------------------------------------------------------------------

    void Level_Set_Geometry::set_node_manager( Node_Manager& aNodeManager )
    {
        mNodeManager = &aNodeManager;
    }

    //--------------------------------------------------------------------------------------------------------------

    real
    Level_Set_Geometry::get_isocontour_threshold() const
    {
        return mParameters.mIsocontourThreshold;
    }

    //--------------------------------------------------------------------------------------------------------------

    bool Level_Set_Geometry::depends_on_advs() const
    {
        return mField->has_advs();
    }

    //--------------------------------------------------------------------------------------------------------------

    Geometric_Region Level_Set_Geometry::get_geometric_region(
            uint                    aNodeIndex,
            const Matrix< DDRMat >& aNodeCoordinates )
    {
        // If it's a background node, can get geometric region from field value
        if ( mNodeManager->is_background_node( aNodeIndex ) )
        {
            return this->determine_geometric_region( this->get_field_value( aNodeIndex, aNodeCoordinates ) );
        }
        else
        {
            // Get derived node
            const Derived_Node& tDerivedNode = mNodeManager->get_derived_node( aNodeIndex );

            // If derived node knows it is on this interface, can return interface
            if ( tDerivedNode.is_on_interface( *this ) )
            {
                return Geometric_Region::INTERFACE;
            }
            else
            {
                // Get locators
                const Cell< Basis_Node >& tLocators = tDerivedNode.get_locator_nodes();

                // If we only have 2 locators, can use special logic if at least one node is on the interface
                if ( tLocators.size() == 2 )
                {
                    if ( this->get_geometric_region( tLocators( 0 ).get_index(), tLocators( 0 ).get_global_coordinates() ) == Geometric_Region::INTERFACE )
                    {
                        // Take geometric region of second node
                        return this->get_geometric_region( tLocators( 1 ).get_index(), tLocators( 1 ).get_global_coordinates() );
                    }
                    else if ( this->get_geometric_region( tLocators( 1 ).get_index(), tLocators( 1 ).get_global_coordinates() ) == Geometric_Region::INTERFACE )
                    {
                        // Take geometric region of first node
                        return this->get_geometric_region( tLocators( 0 ).get_index(), tLocators( 0 ).get_global_coordinates() );
                    }
                    else
                    {
                        // Resort to field value
                        return this->determine_geometric_region( this->get_field_value( aNodeIndex, aNodeCoordinates ) );
                    }
                }
                else
                {
                    // Must use field value for more than 2 locators
                    return this->determine_geometric_region( this->get_field_value( aNodeIndex, aNodeCoordinates ) );
                }
            }
        }
    }

    //--------------------------------------------------------------------------------------------------------------

    Intersection_Node* Level_Set_Geometry::create_intersection_node(
            uint                     aNodeIndex,
            const Cell< Node* >&     aBackgroundNodes,
            const Parent_Node&       aFirstParentNode,
            const Parent_Node&       aSecondParentNode,
            mtk::Geometry_Type       aBackgroundGeometryType,
            mtk::Interpolation_Order aBackgroundInterpolationOrder )
    {
        if ( this->use_multilinear_interpolation() )
        {
            // Create multilinear intersection node
            return new Intersection_Node_Bilinear(
                    aNodeIndex,
                    aBackgroundNodes,
                    aFirstParentNode,
                    aSecondParentNode,
                    aBackgroundGeometryType,
                    aBackgroundInterpolationOrder,
                    *this );
        }
        else
        {
            // Create linear intersection node
            return new Intersection_Node_Linear(
                    aNodeIndex,
                    aBackgroundNodes,
                    aFirstParentNode,
                    aSecondParentNode,
                    aBackgroundGeometryType,
                    aBackgroundInterpolationOrder,
                    *this );
        }
    }

    //--------------------------------------------------------------------------------------------------------------

    void Level_Set_Geometry::get_dfield_dcoordinates(
            const Basis_Node& aParentNode,
            Matrix< DDRMat >& aSensitivities ) const
    {
        if ( this->use_multilinear_interpolation() )
        {
            // TODO
            MORIS_ERROR( false, "Multilinear interpolation was detected for an intersection on another intersection. This needs to be implemented." );
        }
        else
        {
            // Get parents
            const Cell< Basis_Node >& tParentNodes = aParentNode.get_locator_nodes();

            // Geometry values
            real tDeltaPhi = this->get_field_value( tParentNodes( 1 ).get_index(), tParentNodes( 1 ).get_global_coordinates() )
                           - this->get_field_value( tParentNodes( 0 ).get_index(), tParentNodes( 0 ).get_global_coordinates() );

            // Compute parent vector
            Matrix< DDRMat > tParentVector = tParentNodes( 1 ).get_global_coordinates() - tParentNodes( 0 ).get_global_coordinates();

            // get number of spatial dimensions;
            uint tNumDim = tParentVector.length();

            // Compute square of length of parent vector
            real tParentLengthSquared = dot( tParentVector, tParentVector );

            // Sensitivities: dPhi/dx_i  = delta(Phi) / L_i where L_i = PaerentVectorLenth^2 / (ParentVector * e_i)
            for ( uint tCoordinateIndex = 0; tCoordinateIndex < tNumDim; tCoordinateIndex++ )
            {
                aSensitivities( tCoordinateIndex ) = tDeltaPhi * tParentVector( tCoordinateIndex ) / ( tParentLengthSquared + MORIS_REAL_EPS );
            }
        }
    }

    //--------------------------------------------------------------------------------------------------------------

    Cell< std::shared_ptr< mtk::Field > > Level_Set_Geometry::get_mtk_fields()
    {
        // Get MTK field
        std::shared_ptr< mtk::Field > tMTKField = this->get_mtk_field();

        // Add to cell if it exists
        if ( tMTKField )
        {
            return { this->get_mtk_field() };
        }
        else
        {
            return {};
        }
    }

    //--------------------------------------------------------------------------------------------------------------

    Geometric_Region Level_Set_Geometry::determine_geometric_region( real aLevelSetValue ) const
    {
        // Determine if value indicates that this point is on the interface
        if ( std::abs( aLevelSetValue - mParameters.mIsocontourThreshold ) <= mParameters.mIsocontourTolerance )
        {
            return Geometric_Region::INTERFACE;
        }

        // Otherwise, give the region relative to the isocontour threshold
        else if ( aLevelSetValue - mParameters.mIsocontourThreshold < 0 )
        {
            return Geometric_Region::NEGATIVE;
        }
        else
        {
            return Geometric_Region::POSITIVE;
        }
    }

    //--------------------------------------------------------------------------------------------------------------

    std::string
    Level_Set_Geometry::get_name()
    {
        return Design_Field::get_name();
    }

    //--------------------------------------------------------------------------------------------------------------

    void
    Level_Set_Geometry::import_advs( sol::Dist_Vector* aOwnedADVs )
    {
        return Design_Field::import_advs( aOwnedADVs );
    }

    //--------------------------------------------------------------------------------------------------------------

    void
    Level_Set_Geometry::reset_nodal_data( mtk::Interpolation_Mesh* aInterpolationMesh )
    {
        Design_Field::reset_nodal_data( aInterpolationMesh );
    }

    //--------------------------------------------------------------------------------------------------------------

    void
    Level_Set_Geometry::discretize(
            mtk::Mesh_Pair          aMeshPair,
            sol::Dist_Vector*       aOwnedADVs,
            const Matrix< DDSMat >& aSharedADVIds,
            uint                    aADVOffsetID )
    {
        Design_Field::discretize( aMeshPair, aOwnedADVs, aSharedADVIds, aADVOffsetID );
    }

    //--------------------------------------------------------------------------------------------------------------

    void
    Level_Set_Geometry::discretize(
            std::shared_ptr< mtk::Field > aMTKField,
            mtk::Mesh_Pair                aMeshPair,
            sol::Dist_Vector*             aOwnedADVs,
            const Matrix< DDSMat >&       aSharedADVIds,
            uint                          aADVOffsetID )
    {
        Design_Field::discretize( aMTKField, aMeshPair, aOwnedADVs, aSharedADVIds, aADVOffsetID );
    }

    //--------------------------------------------------------------------------------------------------------------

    void Level_Set_Geometry::get_design_info(
            uint                    aNodeIndex,
            const Matrix< DDRMat >& aCoordinates,
            Cell< real >&           aOutputDesignInfo )
    {
        aOutputDesignInfo.resize( 1 );
        aOutputDesignInfo( 0 ) = Design_Field::get_field_value( aNodeIndex, aCoordinates );
    }

    bool Level_Set_Geometry::intended_discretization()
    {
        return ( mParameters.mDiscretizationIndex >= 0 );
    }

    //--------------------------------------------------------------------------------------------------------------

    moris_index
    Level_Set_Geometry::get_discretization_mesh_index()
    {
        MORIS_ASSERT( mParameters.mDiscretizationIndex >= 0,
                "A discretization is not intended for this field. Check this with intended_discretization() first." );

        return mParameters.mDiscretizationIndex;
    }

    //--------------------------------------------------------------------------------------------------------------

    real
    Level_Set_Geometry::get_discretization_lower_bound()
    {
        return mParameters.mDiscretizationLowerBound;
    }

    //--------------------------------------------------------------------------------------------------------------

    real
    Level_Set_Geometry::get_discretization_upper_bound()
    {
        return mParameters.mDiscretizationUpperBound;
    }

}    // namespace moris::ge
