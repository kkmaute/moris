/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_Geometry.hpp
 *
 */

#ifndef INCLUDE_GEOMETRY_CL_GEOMETRY_HPP_
#define INCLUDE_GEOMETRY_CL_GEOMETRY_HPP_

#include "assert/fn_xtk_assert.hpp"
#include "linalg/cl_XTK_Matrix.hpp"
#include "cl_MTK_Enums.hpp"

using namespace moris;
namespace xtk
{

    /*
     * Any geometry must inherit either the discrete or analytic functions interfaces provided in include/geometry
     * For examples on creating analytic or discrete geometries see examples/cl_Geometry_Examples.cpp.
     */

    class Geometry
    {
      public:
        ~Geometry()
        {
        }

        /*
         * **************************************************************************************
         *              REQUIRED INTERFACE FOR BOTH ANALYTIC AND DISCRETE GEOMETRY
         * **************************************************************************************
         */

        /**
         * This tells the user whether the given geometry is an analytical expression
         * @param[out] True if analytic geometry (i.e. Sphere from an equation)
         */
        virtual bool is_analytic() const = 0;

        /**
         *
         */
        virtual void get_dphi_dp_size( moris::size_t &aNumRows, moris::size_t &aNumCols ) const = 0;

        /*
         * **************************************************************************************
         *                             ANALYTIC GEOMETRY INTERFACE
         * **************************************************************************************
         */

        /*
         * Given a node coordinate, the geometry needs to return the distance to the nearest function.
         */
        virtual moris::real evaluate_field_value_with_coordinate( moris::size_t const &aRowIndex,
                moris::Matrix< DDRMat > const                                         &aCoordinates ) const
        {
            XTK_ERROR << "evaluate_field_value_with_coordinate not implemented. This could be due to a geometry not being based on an analytic expression.";
            return 0;
        }

        /*
         * Given a node coordinate @param[in] aCoordinates, the function returns a matrix of relevant node coordinates
         * Where each cell represents a design variable x,y,z vector
         */
        virtual moris::Matrix< DDRMat > evaluate_sensitivity_dphi_dp_with_coordinate( moris::size_t const &aRowIndex,
                moris::Matrix< DDRMat > const                                                             &aCoordinates ) const
        {
            XTK_ERROR << "evaluate_sensitivity_dx_dp not implemented. This could be due to a geometry not being based on an analytic expression.";
            return moris::Matrix< DDRMat >();
        }

        /*
         * **************************************************************************************
         *                             DISCRETE GEOMETRY INTERFACE
         * **************************************************************************************
         */

        /*
         * Given a node index and rank, provide the field value or minimum distance to a geometry feature
         * NOTE: THIS CURRENTLY REQUIRES THE BACKGROUND MESH AND GEOMETRY ARE COINCIDENT (NEEDS AN OBJECT THAT MAPS BETWEEN THE TWO)
         */
        virtual moris::real access_field_value_with_entity_index(
                moris::moris_index aEntityIndex,
                mtk::EntityRank    aEntityRank ) const
        {
            XTK_ERROR << " access_field_value_with_entity_index not implemented. This could be due to a geometry not being based on a mesh.";
            return 0;
        }

        /*
         * Given a local coordinate, node index and entity rank, the function returns a matrix of relevant node coordinates
         *  NOTE: THIS CURRENTLY REQUIRES THE BACKGROUND MESH AND GEOMETRY ARE COINCIDENT (NEEDS AN OBJECT THAT MAPS BETWEEN THE TWO)
         */
        virtual moris::Matrix< DDRMat > evaluate_sensitivity_dphi_dp( moris::Matrix< DDRMat > const &aLocalCoordinate, moris::size_t aEntityIndex, mtk::EntityRank aEntityRank )
        {
            //        XTK_ERROR<<" evaluate_sensitivity_dx_dp not implemented. This could be due to a geometry not being based on a mesh.";
            return moris::Matrix< DDRMat >( 1, 1, 0 );
        }

        /*
         * Returns the advs that the provided node indices are dependent on
         */
        virtual moris::Matrix< IndexMat > get_node_adv_indices( moris::Matrix< IndexMat > const &aNodeIndices )
        {
            XTK_ERROR << " get_node_adv_indices not implemented, This could be due to a geometry not being based on a mesh.";
            return moris::Matrix< IndexMat >( 1, 1, 0 );
        }
    };
}    // namespace xtk

#endif /* INCLUDE_GEOMETRY_CL_GEOMETRY_HPP_ */
