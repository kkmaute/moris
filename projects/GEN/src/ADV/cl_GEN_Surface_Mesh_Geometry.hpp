/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_GEN_Surface_Mesh_Geometry.hpp
 *
 */

#pragma once

#include "fn_SDF_Raycast.hpp"

#include "cl_GEN_Design_Field.hpp"
#include "cl_GEN_Geometry.hpp"
#include "GEN_Data_Types.hpp"

namespace moris::ge
{
    /**
     * This is a struct used to simplify \ref moris::ge::Surface_Mesh_Geometry constructors. It contains all field and level-set parameters.
     */
    struct Surface_Mesh_Parameters
    {
        Matrix< DDRMat > mOffsets;
        std::string      mFilePath;

        /**
         * Constructor with a given parameter list
         *
         * @param aParameterList Parameter list with level set geometry parameters
         */
        explicit Surface_Mesh_Parameters( const ParameterList& aParameterList = prm::create_surface_mesh_geometry_parameter_list() );
    };

    class Surface_Mesh_Geometry : public sdf::Object
            , public Geometry
            , public std::enable_shared_from_this< Surface_Mesh_Geometry >    // TODO make it so we don't need enable_shared_from_this, should be possible in intersection node
    {
      private:
        Surface_Mesh_Parameters mParameters;

      public:
        /**
         * Constructor taking in a field pointer and a set of parameters.
         *
         * @param aField Field for computing nodal values
         * @param aParameters Field parameters
         */
        explicit Surface_Mesh_Geometry( Surface_Mesh_Parameters aParameters = Surface_Mesh_Parameters() );

        /**
         * Gets the geometric region of a node, based on this geometry.
         *
         * @param aNodeIndex Node index
         * @param aNodeCoordinates Node coordinates
         * @return Geometric region enum
         */
        Geometric_Region get_geometric_region(
                uint                    aNodeIndex,
                const Matrix< DDRMat >& aNodeCoordinates ) override;

        /**
         * Creates an intersection node based on the given information. The intersection node may or may not represent an intersection;
         * that is, its position may lie outside of the edge definition based on the given nodal coordinates. This information can be
         * requested from the created intersection node.
         *
         * @param aNodeIndex Node index of the new intersection node
         * @param aBaseNodes Base nodes of the element where the intersection lies
         * @param aFirstParentNode Node marking the starting point of the intersection edge
         * @param aSecondParentNode Node marking the ending point of the intersection edge
         * @param aBaseGeometryType Geometry type of the base node element
         * @return New intersection node
         */
        Intersection_Node* create_intersection_node(
                uint                 aNodeIndex,
                const Cell< Node* >& aBaseNodes,
                const Parent_Node&   aFirstParentNode,
                const Parent_Node&   aSecondParentNode,
                mtk::Geometry_Type   aBaseGeometryType ) override;

        /**
         * Gets an MTK field, if this geometry uses one that needs to be remapped to a new mesh
         *
         * @return MTK field
         */
        Cell< std::shared_ptr< mtk::Field > > get_mtk_fields() override;

      private:
        void find_candidate_ancestors()
        {
        }
    };
}    // namespace moris::ge
