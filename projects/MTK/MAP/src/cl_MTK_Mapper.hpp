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
        class Discrete_Field;
        class BSpline_Field;

        class Mapper
        {
            // Source meshes
            Mesh* mSourceMesh;

            // Target meshes
            Mesh* mTargetMesh; // FIXME

            // Mesh manager- needed for FEM Model
            std::shared_ptr<Mesh_Manager> mInputMeshManager;
            uint mMeshPairIndex;

            //         fem::IWG_L2                  * mIWG;
            mdl::Model                   * mModel = nullptr;
            uint                           mBSplineMeshIndex; // FIXME

            //         Cell< Node* >                  mNodes;

            bool mHaveIwgAndModel = false;
            //         bool mHaveNodes       = false;
            //         real mFilterRadius = 0;

        public:

            /**
             * Constructor
             *
             * @param aInputMeshManager Mesh manager for input mesh
             * @param aMeshPairIndex Mesh pair index
             * @param aBSplineMeshIndex B-spline mesh index
             */
            Mapper( std::shared_ptr<Mesh_Manager> aInputMeshManager,
                    uint aMeshPairIndex,
                    uint aBSplineMeshIndex = 0 );

            /**
             * Destructor
             */
            ~Mapper();

            //------------------------------------------------------------------------------

            void map_input_field_to_output_field(
                    std::shared_ptr<Field>       aFieldSource,
                    std::shared_ptr<BSpline_Field> aFieldTarget );

            //------------------------------------------------------------------------------

            void interpolate_field(
                    std::shared_ptr<Field>          aFieldSource,
                    std::shared_ptr<Discrete_Field> aFieldTarget);

            //------------------------------------------------------------------------------

            void change_field_order(
                    std::shared_ptr<Field>          aFieldSource,
                    std::shared_ptr<Discrete_Field> aFieldTarget );

            //------------------------------------------------------------------------------

            void perform_mapping(
                    std::string aSourceLabel,
                    EntityRank  aSourceEntityRank,
                    std::string aTargetLabel,
                    EntityRank  aTargetEntityRank );

            void perform_mapping(
                    const Matrix<DDRMat>& aSourceField,
                    EntityRank            aSourceEntityRank,
                    Matrix<DDRMat>&       aTargetField,
                    EntityRank            aTargetEntityRank );

            //------------------------------------------------------------------------------

//            void perform_filter(
//                    std::string      aSourceLabel,
//                    real             aFilterRadius,
//                    Matrix< DDRMat > & aValues );

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
            void set_l2_alpha( real aAlpha );

        private:

            void map_node_to_bspline( Matrix<DDRMat>& aSolution );

            void map_node_to_bspline_same_mesh(
                    moris_index aSourceIndex,
                    moris_index aTargetIndex,
                    EntityRank  aBSplineRank );

            void map_node_to_bspline_from_field(
                    const Matrix<DDRMat>& aSourceField,
                    Matrix<DDRMat>&       aTargetField,
                    EntityRank            aBSplineRank );

            //------------------------------------------------------------------------------

            void map_node_to_bspline_from_field(
                    std::shared_ptr<Field>         aSourceField,
                    std::shared_ptr<BSpline_Field> aTargetField,
                    std::shared_ptr<Mesh_Manager> aMeshManager,
                    uint aMeshIndex);

            //------------------------------------------------------------------------------

            void map_bspline_to_node_same_mesh(
                    moris_index     aSourceIndex,
                    EntityRank      aBSplineRank,
                    moris_index     aTargetIndex );

            //------------------------------------------------------------------------------

//            void map_bspline_to_node_same_mesh( mtk::std::shared_ptr<BSpline_Field> aField );

            void create_iwg_and_model(real aAlpha = 0.0);

            //------------------------------------------------------------------------------

            void create_iwg_and_model(
                    std::shared_ptr<Mesh_Manager> aMeshManager,
                    uint aMeshIndex,
                    uint aDiscretizationMeshIndex,
                    real aAlpha = 0.0);

            //------------------------------------------------------------------------------

            void create_nodes_for_filter();

            //------------------------------------------------------------------------------

            void calculate_filter_weights( real aFilterRadius );

            //------------------------------------------------------------------------------
        };
    } /* namespace mtk */
} /* namespace moris */


#endif /* PROJECTS_MTK_MAP_SRC_CL_MTK_MAPPER_HPP_ */
