/*
 * cl_MTK_Mapper.hpp
 *
 *  Created on: Oct 30, 2018
 *      Author: messe
 */

#ifndef PROJECTS_MTK_MAP_SRC_CL_MTK_MAPPER_HPP_
#define PROJECTS_MTK_MAP_SRC_CL_MTK_MAPPER_HPP_

#include "typedefs.hpp"
#include "cl_Cell.hpp"

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
    namespace mapper
    {
        //------------------------------------------------------------------------------

        class Mapper
        {
                // Source meshes
                moris::mtk::Mesh* mSourceMesh;

                // Target meshes
                moris::mtk::Mesh* mTargetMesh;

                // Mesh manager- needed for FEM Model
                moris::moris_index             mMeshPairIndex_In;
                moris::moris_index             mMeshPairIndex_Out;
                std::shared_ptr<mtk::Mesh_Manager> mMeshManager;
                //         fem::IWG_L2                  * mIWG;
                mdl::Model                   * mModel = nullptr;
                uint                           mBSplineMeshIndex;

                mtk::Field * mFieldIn = nullptr;
                mtk::Field * mFieldOut = nullptr;

                //         moris::Cell< Node* >                  mNodes;

                bool mHaveIwgAndModel = false;
                //         bool mHaveNodes       = false;
                //         real mFilterRadius = 0;

                //------------------------------------------------------------------------------
            public:
                //------------------------------------------------------------------------------

                /**
                 * constructor with only one mesh
                 */
                Mapper( std::shared_ptr<mtk::Mesh_Manager> aMeshManager,
                        const moris_index aMeshPairIndex,
                        const uint aBSplineMeshIndex = 0 );

                //------------------------------------------------------------------------------

                Mapper( mtk::Field                         * aFieldIn,
                        mtk::Field                         * aFieldOut );

                //------------------------------------------------------------------------------

                /**
                 * destructor
                 */
                ~Mapper();

                //------------------------------------------------------------------------------

                void map_input_field_to_output_field();

                //------------------------------------------------------------------------------

                void interpolate_field(
                        mtk::Field * aFieldSource,
                        mtk::Field * aFieldTarget);

                //------------------------------------------------------------------------------

                void perform_mapping( const std::string      & aSourceLabel,
                        const enum EntityRank    aSourceEntityRank,
                        const std::string      & aTargetLabel,
                        const enum EntityRank    aTargetEntityRank );

                void perform_mapping( const Matrix<DDRMat>& aSourceField,
                        const enum EntityRank aSourceEntityRank,
                        Matrix<DDRMat>&       aTargetField,
                        const enum EntityRank aTargetEntityRank );

                void perform_mapping(
                        mtk::Field * aField,
                        const enum EntityRank aSourceEntityRank,
                        const enum EntityRank aTargetEntityRank );

                //------------------------------------------------------------------------------

                void perform_filter( const std::string      & aSourceLabel,
                        const real             & aFilterRadius,
                        Matrix< DDRMat > & aValues );

                //------------------------------------------------------------------------------

                /*void
         project_coeffs_from_node_data(
                 const Matrix< DDRMat > & aNodeValues,
                 const uint             & aBSplineOrder,
                 std::shared_ptr< Mesh >  aMesh,
                 Matrix< DDRMat >       & aCoeffs ); */

                //------------------------------------------------------------------------------

                /*
                 * set the parameter for the L2 projection
                 */
                void set_l2_alpha( const real & aAlpha );

                //------------------------------------------------------------------------------
            private:
                //------------------------------------------------------------------------------

                void map_node_to_bspline( Matrix<DDRMat>& aSolution );

                void map_node_to_bspline_same_mesh(
                        const moris_index        aSourceIndex,
                        const moris_index        aTargetIndex,
                        const enum EntityRank    aBSplineRank );

                void map_node_to_bspline_from_field( const Matrix<DDRMat>& aSourceField,
                        Matrix<DDRMat>&       aTargetField,
                        const enum EntityRank aBSplineRank );

                //------------------------------------------------------------------------------

                void map_node_to_bspline_from_field( mtk::Field * aField );

                ////------------------------------------------------------------------------------
                //
                //         void map_node_to_element_same_mesh(
                //                          const moris_index   aSourceIndex,
                //                          const moris_index   aTargetIndex );

                //------------------------------------------------------------------------------

                void map_bspline_to_node_same_mesh(
                        const moris_index     aSourceIndex,
                        const enum EntityRank aBSplineRank,
                        const moris_index     aTargetIndex );

                //------------------------------------------------------------------------------

                void create_iwg_and_model( const real aAlpha = 0.0 );

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
