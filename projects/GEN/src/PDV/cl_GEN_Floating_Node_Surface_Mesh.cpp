/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_GEN_Floating_Node_Surface_Mesh.cpp
 *
 */

#include "cl_GEN_Floating_Node_Surface_Mesh.hpp"
#include "cl_GEN_Parent_Node.hpp"
#include "cl_GEN_Surface_Mesh_Geometry.hpp"

namespace moris::gen
{
    Floating_Node_Surface_Mesh::Floating_Node_Surface_Mesh(
            uint                              aNodeIndex,
            const Vector< Background_Node* >& aBackgroundNodes,
            const Matrix< DDRMat >&           aParametricCoordinates,
            uint                              aParentVertex,
            mtk::Geometry_Type                aBackgroundGeometryType,
            mtk::Interpolation_Order          aBackgroundInterpolationOrder,
            Surface_Mesh_Geometry&            aInterfaceGeometry )
            : Floating_Node(
                      aNodeIndex,
                      aBackgroundNodes,
                      aParametricCoordinates,
                      aBackgroundGeometryType,
                      aBackgroundInterpolationOrder )
            , mParentVertex( aParentVertex )
            , mInterfaceGeometry( aInterfaceGeometry )
    {
    }

    //--------------------------------------------------------------------------------------------------------------

    Geometry& Floating_Node_Surface_Mesh::get_interface_geometry()
    {
        return mInterfaceGeometry;
    }

    //--------------------------------------------------------------------------------------------------------------

    const Geometry& Floating_Node_Surface_Mesh::get_interface_geometry() const
    {
        return mInterfaceGeometry;
    }

    //--------------------------------------------------------------------------------------------------------------

    bool Floating_Node_Surface_Mesh::depends_on_advs() const
    {
        return mInterfaceGeometry.facet_vertex_depends_on_advs( mParentVertex );
    }

    //--------------------------------------------------------------------------------------------------------------

    void Floating_Node_Surface_Mesh::append_dcoordinate_dadv(
            Matrix< DDRMat >&       aCoordinateSensitivities,
            const Matrix< DDRMat >& aSensitivityFactor ) const
    {
        // Since the floating node lies on the surface mesh vertex, its sensitivities are identical to the vertex
        Matrix< DDRMat > tSensitivitiesToAdd = mInterfaceGeometry.facet_vertex_depends_on_advs( mParentVertex ) ? aSensitivityFactor * mInterfaceGeometry.get_dvertex_dadv( mParentVertex )
                                                                                                                : Matrix< DDRMat >();

        // Resize sensitivities
        uint tJoinedSensitivityLength = aCoordinateSensitivities.n_cols();
        aCoordinateSensitivities.resize( tSensitivitiesToAdd.n_rows(),
                tJoinedSensitivityLength + tSensitivitiesToAdd.n_cols() );

        // Join sensitivities
        for ( uint iCoordinateIndex = 0; iCoordinateIndex < tSensitivitiesToAdd.n_rows(); iCoordinateIndex++ )
        {
            for ( uint iAddedSensitivity = 0; iAddedSensitivity < tSensitivitiesToAdd.n_cols(); iAddedSensitivity++ )
            {
                aCoordinateSensitivities( iCoordinateIndex, tJoinedSensitivityLength + iAddedSensitivity ) =
                        tSensitivitiesToAdd( iCoordinateIndex, iAddedSensitivity );
            }
        }
    }

    //--------------------------------------------------------------------------------------------------------------

    Vector< sint >
    Floating_Node_Surface_Mesh::get_coordinate_determining_adv_ids() const
    {
        // Since the floating node is exactly its parent vertex, its sensitivities are exactly the sensitivities of the parent vertex
        return mInterfaceGeometry.get_vertex_adv_ids( mParentVertex );
    }

    //--------------------------------------------------------------------------------------------------------------

}    // namespace moris::gen
