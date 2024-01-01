/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_MTK_Mapper.hpp
 *
 */

#ifndef PROJECTS_MTK_MAP_SRC_CL_MTK_MAPPER_HPP_
#define PROJECTS_MTK_MAP_SRC_CL_MTK_MAPPER_HPP_

#include "moris_typedefs.hpp"
#include "cl_Vector.hpp"

#include "cl_Matrix.hpp"
#include "linalg_typedefs.hpp"

namespace moris
{
    //------------------------------------------------------------------------------

    namespace MSI
    {
        class Equation_Object;
    }

    //------------------------------------------------------------------------------

    namespace fem
    {
        class IWG_L2;
        class Node_Base;
    }

    //------------------------------------------------------------------------------

    namespace mdl
    {
        class Model;
    }

    //------------------------------------------------------------------------------

    namespace mtk
    {
        class Mesh;
        class Mesh_Manager;
        class Field;
    }

    //------------------------------------------------------------------------------
    namespace mtk
    {
        //------------------------------------------------------------------------------

        class Mapper
        {
                mdl::Model   * mModel = nullptr;

                bool mHaveIwgAndModel = false;

                //------------------------------------------------------------------------------
            public:
                //------------------------------------------------------------------------------

                //------------------------------------------------------------------------------

                Mapper();

                //------------------------------------------------------------------------------

                /**
                 * destructor
                 */
                ~Mapper();

                //------------------------------------------------------------------------------

                void map_input_field_to_output_field(
                        mtk::Field * aFieldSource,
                        mtk::Field * aFieldTarget );

                void map_input_field_to_output_field_2( mtk::Field * aFieldSource);

                //------------------------------------------------------------------------------

                void perform_mapping(
                        mtk::Field * aField,
                        const enum EntityRank aSourceEntityRank,
                        const enum EntityRank aTargetEntityRank );

                //------------------------------------------------------------------------------

                /*
                 void perform_filter
                         const std::string      & aSourceLabel,
                         const real             & aFilterRadius,
                         Matrix< DDRMat > & aValues );
                 */

                //------------------------------------------------------------------------------

                /*void
                project_coeffs_from_node_data(
                 const Matrix< DDRMat > & aNodeValues,
                 const uint             & aBSplineOrder,
                 std::shared_ptr< Mesh >  aMesh,
                 Matrix< DDRMat >       & aCoeffs ); */

                //------------------------------------------------------------------------------
            private:
                //------------------------------------------------------------------------------

                //------------------------------------------------------------------------------

                void interpolate_field(
                        mtk::Field * aFieldSource,
                        mtk::Field * aFieldTarget);

                //------------------------------------------------------------------------------

                void change_field_order(
                        mtk::Field * aFieldSource,
                        mtk::Field * aFieldTarget );

                //------------------------------------------------------------------------------

                void map_node_to_bspline(  mtk::Field * aField );

                //------------------------------------------------------------------------------

                void map_node_to_bspline_from_field( mtk::Field * aField );

                ////------------------------------------------------------------------------------
                //
                //         void map_node_to_element_same_mesh(
                //                          const moris_index   aSourceIndex,
                //                          const moris_index   aTargetIndex );

                //------------------------------------------------------------------------------

                void
                map_bspline_to_node_same_mesh( mtk::Field * aField );

                //------------------------------------------------------------------------------

                void create_iwg_and_model(
                        mtk::Field * aField,
                        const real aAlpha = 0.0 );

                //------------------------------------------------------------------------------

                void create_nodes_for_filter();

                //------------------------------------------------------------------------------

                void calculate_filter_weights( const real & aFilterRadius );

                //------------------------------------------------------------------------------
        };

        //------------------------------------------------------------------------------
    } /* namespace mtk */
} /* namespace moris */

#endif /* PROJECTS_MTK_MAP_SRC_CL_MTK_MAPPER_HPP_ */

