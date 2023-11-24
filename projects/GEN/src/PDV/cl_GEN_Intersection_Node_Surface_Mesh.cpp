/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_GEN_Intersection_Node_Surface_Mesh.cpp
 *
 */

#include <limits>
#include "cl_GEN_Intersection_Node_Surface_Mesh.hpp"

namespace moris::ge
{
    Intersection_Node_Surface_Mesh::Intersection_Node_Surface_Mesh(
        std::shared_ptr< Intersection_Node > aFirstParentNode,
        std::shared_ptr< Intersection_Node > aSecondParentNode,
        uint                                 aFirstParentNodeIndex,
        uint                                 aSecondParentNodeIndex,
        const Matrix< DDRMat >&              aFirstParentNodeCoordinates,
        const Matrix< DDRMat >&              aSecondParentNodeCoordinates,
        const Matrix< DDUMat >&              aAncestorNodeIndices,
        const Cell< Matrix< DDRMat > >&      aAncestorNodeCoordinates )
        : Intersection_Node(
            compute_local_coordinate(
                aFirstParentNodeCoordinates,
                aSecondParentNodeCoordinates,
                aAncestorNodeCoordinates
            ),
            aFirstParentNode,
            aSecondParentNode,
            aFirstParentNodeIndex,
            aSecondParentNodeIndex,
            { { -1 } },
            { { 1 } },
            aAncestorNodeIndices,
            aAncestorNodeCoordinates,
            Element_Intersection_Type::Linear_1D )
    {
        // Ensure that only two ancestor nodes are provided
        MORIS_ERROR( aAncestorNodeCoordinates.size() != 2, "GEN - Intersection_Node_Surface_Mesh - Exactly two surface mesh nodes must be supplied." );
    }
    
    Matrix< DDRMat >
    Intersection_Node_Surface_Mesh::compute_global_coordinates()
    {
        // TODO
        MORIS_ERROR( false, "Intersection_Node_Surface_Mesh - compute_global_coordinates() not implemented yet.");
        return { { std::numeric_limits<double>::quiet_NaN() } };
    }

    bool
    Intersection_Node_Surface_Mesh::determine_is_intersected(
        const Element_Intersection_Type aAncestorBasisFunction,
        const Matrix< DDRMat >&         aFirstParentNodeLocalCoordinates,
        const Matrix< DDRMat >&         aSecondParentNodeLocalCoordinates
    )
    {
        // TODO
        MORIS_ERROR( false, "Intersection_Node_Surface_Mesh - determine_is_intersected() not implemented yet.");
        return std::numeric_limits<double>::quiet_NaN();
    }

    real
    Intersection_Node_Surface_Mesh::compute_local_coordinate( 
        const Matrix< DDRMat >& aFirstParentNodeCoordinates, 
        const Matrix< DDRMat >& aSecondParentNodeCoordinates, 
        const Cell< Matrix< DDRMat > >& aAncestorNodeCoordinates )
    {
        switch( mAncestorNodeCoordinates( 0 ).numel() )
        {
            case 2:
                {
                    // assign the x and y coordinates of the parent nodes
                    real tXp1 = aFirstParentNodeCoordinates( 0, 0);
                    real tXp2 = aSecondParentNodeCoordinates( 0, 0 );
                    real tYp1 = aFirstParentNodeCoordinates( 1, 0 );
                    real tYp2 = aSecondParentNodeCoordinates( 1, 0 );
                    real tXa1 = aAncestorNodeCoordinates( 0 )( 0, 0 );
                    real tXa2 = aAncestorNodeCoordinates( 1 )( 0, 0 );
                    real tYa1 = aAncestorNodeCoordinates( 0 )( 1, 0 );
                    real tYa2 = aAncestorNodeCoordinates( 1 )( 1, 0 );

                    real tLocalCoordinate = ( 0.5*( tYp1 + tYp2 ) - tYa1 + ( ( tYa1 - tYa2 ) * ( tXp1 - 0.5 * ( tXp1 + tXp2 ) ) ) / ( tXa1 - tXa2 )) 
                        / ( 0.5 * ( tYp2 - tYp1 ) + 0.5 * ( tXp1 - tXp2 ) / ( tXa1 - tXa2 ) );

                    return tLocalCoordinate;
                }
            case 3:
                {
                    // TODO: facet intersection
                    MORIS_ERROR( false, "GEN - Intersection_Node_Surface_Mesh - compute_local_coordinate: 3D Surface mesh intersections not yet implemented." );
                    return std::numeric_limits<double>::quiet_NaN();
                }
            default:
                {
                    MORIS_ERROR( false, "GEN - Intersection_Node_Surface_Mesh - compute_local_coordinate: Surface Mesh intersection not implemented for %lu dimensional geometries", aFirstParentNodeCoordinates.numel() );
                    return std::numeric_limits<double>::quiet_NaN();
                }
        }
        return std::numeric_limits<double>::quiet_NaN();
    }
    
    void 
    Intersection_Node_Surface_Mesh::get_dcoordinate_dadv( Matrix< DDRMat >& aCoordinateSensitivities, const Matrix< DDRMat >& aSensitivityFactor )
    {
        // TODO
        MORIS_ERROR( false, "Intersection_Node_Surface_Mesh - get_dcoordinate_dadv() not implemented yet.");

        return;
    };
    
    Matrix< DDSMat > 
    Intersection_Node_Surface_Mesh::get_coordinate_determining_adv_ids()
    {
        // TODO
        MORIS_ERROR( false, "Intersection_Node_Surface_Mesh - get_coordinate_determining_adv_ids() not implemented yet.");
        return { { -1 } };
    }

    real 
    Intersection_Node_Surface_Mesh::get_dxi_dfield_from_ancestor( uint aAncestorIndex )
    {
        // TODO
        MORIS_ERROR( false, "Intersection_Node_Surface_Mesh - get_dxi_dfield_from_ancestor() not implemented yet.");

        return std::numeric_limits<double>::quiet_NaN();
    }

    Matrix< DDRMat > 
    Intersection_Node_Surface_Mesh::get_dxi_dcoordinate_first_parent()
    {
        // TODO
        MORIS_ERROR( false, "Intersection_Node_Surface_Mesh - get_dxi_dcoordinate_first_parent() not implemented yet.");
        return { { std::numeric_limits<double>::quiet_NaN() } };
    }

    Matrix< DDRMat > 
    Intersection_Node_Surface_Mesh::get_dxi_dcoordinate_second_parent()
    {
        // TODO
        MORIS_ERROR( false, "Intersection_Node_Surface_Mesh - get_dxi_dcoordinate_second_parent() not implemented yet.");
        return { { std::numeric_limits<double>::quiet_NaN() } };
    }
} // namespace moris::ge
