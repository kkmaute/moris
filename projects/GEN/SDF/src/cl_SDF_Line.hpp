/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_SDF_Line.hpp
 *
 */

#pragma once

#include "moris_typedefs.hpp"
#include "cl_Matrix.hpp"
#include "linalg_typedefs.hpp"
#include "op_minus.hpp"
#include "op_times.hpp"

#include "cl_Vector.hpp"
#include "cl_MTK_Cell.hpp"
#include "cl_SDF_Facet_Vertex.hpp"
#include "cl_SDF_Facet.hpp"

namespace moris::sdf
{
    //-------------------------------------------------------------------------------

    /**
     * the line class for the sdf generator
     */
    class Line : public Facet
    {
        bool mFlag = false;
        //-------------------------------------------------------------------------------

      public:
        //-------------------------------------------------------------------------------

        Line(
                moris_index                                aIndex,
                Vector< std::shared_ptr< Facet_Vertex > >& aVertices );

        //-------------------------------------------------------------------------------

        ~Line(){};

        //-------------------------------------------------------------------------------

        void
        update_data() override;

        /**
         * Computes the Hesse normal distance and the normal vector of the line. Values stored in mHesse and mNormal
         *
         */
        void
        calculate_hesse_normal_form();

        //-------------------------------------------------------------------------------
        // SDF functions
        //-------------------------------------------------------------------------------

        /**
         * @brief Unused function for Line. Errors out when called.
         */
        bool
        check_edge(
                const uint              aEdge,
                const uint              aAxis,
                const Matrix< DDRMat >& aPoint ) override;

        //-------------------------------------------------------------------------------

        /**
         *
         * @brief Returns the distance of a Point to the line.
         * @param[in]  aPoint  point to be considered
         *
         */
        real
        get_distance_to_point(
                const Matrix< DDRMat >& aPoint ) override;

        //-------------------------------------------------------------------------------

        /**
         * Computes the sensitivity of the local coordinate of a parent edge with respect to the facet vertices
         *
         * @param aXi Local coordinate of an intersection along a parent edge, NOT the local coordinate of the facet.
         * @param aRotationMatrix Rotation matrix of the parent edge
         * @return Vector< real > Sensitivities of the local coordinate with respect to the facet vertices. Same size as the number of vertices
         */
        Matrix< DDRMat >
        compute_dxi_dvertices( moris::gen::Intersection_Node_Surface_Mesh& aIntersectionNode ) override;

      private:
    };    // class Line

    //-------------------------------------------------------------------------------

}    // namespace moris::sdf
