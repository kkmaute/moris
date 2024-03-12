/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_GEN_Voxel_Geometry.cpp
 *
 */

#include "cl_GEN_Voxel_Geometry.hpp"
#include "cl_GEN_Voxel_Input.hpp"
#include "cl_GEN_Intersection_Node_Voxel.hpp"
#include "cl_GEN_Parent_Node.hpp"
#include "cl_GEN_Basis_Node.hpp"

namespace moris::gen
{
    //--------------------------------------------------------------------------------------------------------------

    Voxel_Geometry::Voxel_Geometry(
            std::shared_ptr< Voxel_Input > aVoxelInput,
            uint                           aIndex )
            : Geometry( Design_Parameters(), 1e-12 )
            , mVoxelInput( aVoxelInput )
            , mIndex( aIndex )
    {
    }

    //--------------------------------------------------------------------------------------------------------------

    bool Voxel_Geometry::depends_on_advs() const
    {
        return false;
    }

    //--------------------------------------------------------------------------------------------------------------

    Geometric_Region Voxel_Geometry::get_geometric_region(
            uint                    aNodeIndex,
            const Matrix< DDRMat >& aNodeCoordinates )
    {
        // Get node manager
        const Node_Manager& tNodeManager = mVoxelInput->get_node_manager();

        // Test for background node
        if ( mVoxelInput->get_node_manager().is_background_node( aNodeIndex ) )
        {
            // Background node, so we can check the voxel ID from the voxel input
            if ( mVoxelInput->get_voxel_ID( aNodeCoordinates ) == mIndex )
            {
                return Geometric_Region::POSITIVE;
            }
            else
            {
                return Geometric_Region::NEGATIVE;
            }
        }
        else
        {
            // Get derived node
            const Derived_Node& tDerivedNode = tNodeManager.get_derived_node( aNodeIndex );

            // If derived node knows it is on this interface, can return interface
            if ( tDerivedNode.is_on_interface( *this ) )
            {
                return Geometric_Region::INTERFACE;
            }
            else
            {
                // Get locators
                const Vector< Basis_Node >& tLocators = tDerivedNode.get_locator_nodes();

                // Start voting
                uint tNegativeVotes = 0;
                uint tPositiveVotes = 0;

                // Vote using locators
                for ( auto iLocator : tLocators )
                {
                    // Get geometric region
                    Geometric_Region tLocatorRegion = this->get_geometric_region( iLocator.get_index(), iLocator.get_global_coordinates() );

                    // Votes only matter if basis is nonzero
                    if ( iLocator.get_basis() > 1E-8 )
                    {
                        // Register geometric region vote
                        if ( tLocatorRegion == Geometric_Region::NEGATIVE )
                        {
                            tNegativeVotes++;
                        }
                        else if ( tLocatorRegion == Geometric_Region::POSITIVE )
                        {
                            tPositiveVotes++;
                        }
                    }
                }

                // Return based on voting results
                if ( tPositiveVotes == tNegativeVotes )
                {
                    return Geometric_Region::INTERFACE;
                }
                else if ( tPositiveVotes > tNegativeVotes )
                {
                    return Geometric_Region::POSITIVE;
                }
                else
                {
                    return Geometric_Region::NEGATIVE;
                }
            }
        }
    }

    //--------------------------------------------------------------------------------------------------------------

    Intersection_Node* Voxel_Geometry::create_intersection_node(
            uint                              aNodeIndex,
            const Vector< Background_Node* >& aBackgroundNodes,
            const Parent_Node&                aFirstParentNode,
            const Parent_Node&                aSecondParentNode,
            mtk::Geometry_Type                aBackgroundGeometryType,
            mtk::Interpolation_Order          aBackgroundInterpolationOrder )
    {
        return new Intersection_Node_Voxel(
                aNodeIndex,
                aBackgroundNodes,
                aFirstParentNode,
                aSecondParentNode,
                aBackgroundGeometryType,
                aBackgroundInterpolationOrder,
                *this );
    }

    //--------------------------------------------------------------------------------------------------------------

    real Voxel_Geometry::compute_intersection_local_coordinate(
            const Vector< Background_Node* >& aBackgroundNodes,
            const Parent_Node&                aFirstParentNode,
            const Parent_Node&                aSecondParentNode )
    {
        // Get parent geometric regions
        Geometric_Region tFirstParentGeometricRegion  = this->get_geometric_region( aFirstParentNode.get_index(), aFirstParentNode.get_global_coordinates() );
        Geometric_Region tSecondParentGeometricRegion = this->get_geometric_region( aSecondParentNode.get_index(), aSecondParentNode.get_global_coordinates() );

        // Return local coordinate based on geometric regions of the parent nodes
        if ( tFirstParentGeometricRegion == tSecondParentGeometricRegion and tFirstParentGeometricRegion not_eq Geometric_Region::INTERFACE )
        {
            // Geometric regions are the same; no intersection
            return MORIS_REAL_MAX;
        }
        else if ( tFirstParentGeometricRegion == Geometric_Region::INTERFACE )
        {
            // First parent on interface
            return -1.0;
        }
        else if ( tSecondParentGeometricRegion == Geometric_Region::INTERFACE )
        {
            // Second parent on interface
            return 1.0;
        }
        else
        {
            // Interface is midway between the parent nodes
            return 0.0;
        }
    }

    //--------------------------------------------------------------------------------------------------------------

    Vector< std::shared_ptr< mtk::Field > > Voxel_Geometry::get_mtk_fields()
    {
        return {};
    }

    //--------------------------------------------------------------------------------------------------------------

    uint Voxel_Geometry::get_num_fields()
    {
        return 0;
    }

    //--------------------------------------------------------------------------------------------------------------

    std::string Voxel_Geometry::get_name()
    {
        return mVoxelInput->get_file_name() + "_" + std::to_string( mIndex );
    }

    //--------------------------------------------------------------------------------------------------------------

    Vector< std::string > Voxel_Geometry::get_field_names()
    {
        return {};
    }

    //--------------------------------------------------------------------------------------------------------------

    bool Voxel_Geometry::intended_discretization()
    {
        return false;
    }

    //--------------------------------------------------------------------------------------------------------------

    moris_index Voxel_Geometry::get_discretization_mesh_index()
    {
        return -1;
    }

    //--------------------------------------------------------------------------------------------------------------

    real Voxel_Geometry::get_discretization_lower_bound()
    {
        return -1.0;
    }

    //--------------------------------------------------------------------------------------------------------------

    real Voxel_Geometry::get_discretization_upper_bound()
    {
        return 1.0;
    }

    //--------------------------------------------------------------------------------------------------------------

    Vector< std::shared_ptr< Field > > Voxel_Geometry::get_fields()
    {
        return {};
    }

    //--------------------------------------------------------------------------------------------------------------

    void Voxel_Geometry::set_advs( sol::Dist_Vector* aADVs )
    {
    }

    //--------------------------------------------------------------------------------------------------------------

    void Voxel_Geometry::import_advs( sol::Dist_Vector* aOwnedADVs )
    {
    }

    //--------------------------------------------------------------------------------------------------------------

    void Voxel_Geometry::reset_nodal_data( mtk::Interpolation_Mesh* aInterpolationMesh )
    {
    }

    //--------------------------------------------------------------------------------------------------------------

    void Voxel_Geometry::discretize(
            mtk::Mesh_Pair    aMeshPair,
            sol::Dist_Vector* aOwnedADVs )
    {
    }

    //--------------------------------------------------------------------------------------------------------------

    void Voxel_Geometry::discretize(
            std::shared_ptr< mtk::Field > aMTKField,
            mtk::Mesh_Pair                aMeshPair,
            sol::Dist_Vector*             aOwnedADVs )
    {
    }

    //--------------------------------------------------------------------------------------------------------------

    void Voxel_Geometry::get_design_info(
            uint                    aNodeIndex,
            const Matrix< DDRMat >& aCoordinates,
            Vector< real >&         aOutputDesignInfo )
    {
    }

    //--------------------------------------------------------------------------------------------------------------

}    // namespace moris::gen
