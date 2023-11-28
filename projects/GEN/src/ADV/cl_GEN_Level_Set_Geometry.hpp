/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_GEN_Level_Set_Geometry.hpp
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
     * This is a struct used to simplify \ref moris::ge::Level_Set_Geometry constructors. It contains all field and level-set parameters.
     */
    struct Level_Set_Parameters : public Field_Parameters
    {
        Int_Interpolation mIntersectionInterpolation; // The type of interpolation used to determine intersection location
        real              mIsocontourThreshold;       // Level set isocontour level
        real              mIsocontourTolerance;       // Interface tolerance based on geometry value
        real              mIntersectionTolerance;     // Interface tolerance based on intersecction distance

        /**
         * Constructor with a given parameter list
         *
         * @param aParameterList Parameter list with level set geometry parameters
         */
        explicit Level_Set_Parameters( const ParameterList& aParameterList = prm::create_level_set_geometry_parameter_list() );
    };

    class Level_Set_Geometry : public Design_Field, public Geometry, public std::enable_shared_from_this< Level_Set_Geometry > // TODO make it so we don't need enable_shared_from_this, should be possible in intersection node
    {
    private:
        Level_Set_Parameters mParameters;

    public:

        /**
         * Constructor taking in a field pointer and a set of parameters.
         *
         * @param aField Field for computing nodal values
         * @param aParameters Field parameters
         */
        explicit Level_Set_Geometry(
              std::shared_ptr< Field > aField,
              Level_Set_Parameters     aParameters = Level_Set_Parameters(),
              Node_Manager&            aNodeManager = Node_Manager::get_trivial_instance() );

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
         * Accesses the isocontour level that determines the interface for this geometry
         *
         * @return the isocontour level that determines the geometry interface
         */
        real get_isocontour_threshold();

        /**
         * Acccesses the isocontour tolerance for this geometry
         *
         * @return isocontour tolerance
         */
        real get_isocontour_tolerance();

        /**
         * Accesses the intersection tolerance for this geometry
         *
         * @return The real value of the intersection tolerance
         */
        real get_intersection_tolerance();

        /**
         * Gets the geometric region of a node, based on this geometry.
         *
         * @param aNodeIndex Node index
         * @param aNodeCoordinates Node coordinates
         * @return Geometric region enum
         */
        Geometric_Region get_base_geometric_region(
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
        /**
         * Determines the geometric region of a point based on a level set value
         *
         * @param aLevelSetValue Value of the level set function
         * @return Geometric region enum
         */
        Geometric_Region determine_geometric_region( real aLevelSetValue );

    };
}
