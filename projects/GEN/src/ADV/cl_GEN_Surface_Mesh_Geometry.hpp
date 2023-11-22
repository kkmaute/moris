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

#include "cl_GEN_Design_Field.hpp"
#include "cl_GEN_Geometry.hpp"
#include "GEN_Data_Types.hpp"

namespace moris::ge
{
    enum class Int_Interpolation
    {
        LINEAR,
        MULTILINEAR
    };

    /**
     * This is a struct used to simplify \ref moris::ge::Surface_Mesh_Geometry constructors. It contains all field and level-set parameters.
     */
    struct Surface_Mesh_Parameters : public Field_Parameters
    {
        Int_Interpolation mIntersectionInterpolation; // The type of interpolation used to determine intersection location

        /**
         * Constructor with a given parameter list
         *
         * @param aParameterList Parameter list with level set geometry parameters
         */
        explicit Surface_Mesh_Parameters( const ParameterList& aParameterList = prm::create_Surface_Mesh_Geometry_parameter_list() );
    };

    class Surface_Mesh_Geometry : public sdf::Object, public Geometry, public std::enable_shared_from_this< Surface_Mesh_Geometry > // TODO make it so we don't need enable_shared_from_this, should be possible in intersection node
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
        explicit Surface_Mesh_Geometry(
              std::shared_ptr< Field > aField,
              Level_Set_Parameters     aParameters = Level_Set_Parameters() );

        /**
         * Gets the intersection interpolation type for this geometry.
         *
         * @return Intersection interpolation
         */
        Int_Interpolation get_intersection_interpolation();

        /**
         * Gets the mode of intersection used for this geometry
         *
         * @return Intersection_Mode enum
         */
        Intersection_Mode get_intersection_mode();

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
         * @param aEdgeFirstNodeIndex First node index on the intersection edge
         * @param aEdgeSecondNodeIndex Second node index on the intersection edge
         * @param aEdgeFirstIntersectionNode First intersection node on the intersection edge, if it is also an intersection
         * @param aEdgeSecondIntersectionNode Second intersection node on the intersection edge, if it is also an intersection
         * @param aEdgeFirstNodeLocalCoordinates Local coordinates of the first node inside the background element
         * @param aEdgeSecondNodeLocalCoordinates Local coordinates of the second node inside the background element
         * @param aEdgeFirstNodeGlobalCoordinates Global coordinates of the first node
         * @param aEdgeSecondNodeGlobalCoordinates Global coordinates of the second node
         * @param aBackgroundElementNodeIndices Node indices of the background element
         * @param aBackgroundElementNodeCoordinates Node coordinates of the background element
         * @return Created intersection node
         */
        std::shared_ptr< Intersection_Node > create_intersection_node(
                uint                                 aEdgeFirstNodeIndex,
                uint                                 aEdgeSecondNodeIndex,
                std::shared_ptr< Intersection_Node > aEdgeFirstIntersectionNode,
                std::shared_ptr< Intersection_Node > aEdgeSecondIntersectionNode,
                const Matrix< DDRMat >&              aEdgeFirstNodeLocalCoordinates,
                const Matrix< DDRMat >&              aEdgeSecondNodeLocalCoordinates,
                const Matrix< DDRMat >&              aEdgeFirstNodeGlobalCoordinates,
                const Matrix< DDRMat >&              aEdgeSecondNodeGlobalCoordinates,
                const Matrix< DDUMat >&              aBackgroundElementNodeIndices,
                const Cell< Matrix< DDRMat > >&      aBackgroundElementNodeCoordinates ) override;

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
}
