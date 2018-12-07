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
// #include "cl_MTK_Field.hpp"
#include "cl_MTK_Mesh.hpp"
#include "cl_MTK_Mapper_Node.hpp"

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
    namespace mapper
    {
//------------------------------------------------------------------------------

     class Mapper
     {
         std::shared_ptr< mtk::Mesh >       mSourceMesh;
         std::shared_ptr< mtk::Mesh >       mTargetMesh;
         fem::IWG_L2                        * mIWG;
         mdl::Model                         * mModel;
         const uint                           mBSplineOrder;

         moris::Cell< Node* >                  mNodes;

         bool mHaveIwgAndModel = false;
         bool mHaveNodes       = false;
         real mFilterRadius = 0;

//------------------------------------------------------------------------------
     public:
//------------------------------------------------------------------------------

         /**
          * constructor with only one mesh
          */
         Mapper( std::shared_ptr< mtk::Mesh > aMesh, const uint aBSplineOrder=0 );

//------------------------------------------------------------------------------

         /**
          * destructor
          */
         ~Mapper();

//------------------------------------------------------------------------------

         void
         perform_mapping(
                 const std::string      & aSourceLabel,
                 const enum EntityRank    aSourceEntityRank,
                 const std::string      & aTargetLabel,
                 const enum EntityRank    aTargetEntityRank );

//------------------------------------------------------------------------------

         void
         perform_filter(
                 const std::string & aSourceLabel,
                 const real        & aFilterRadius,
                 Matrix< DDRMat >  & aValues );

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
         void
         set_l2_alpha( const real & aAlpha );

//------------------------------------------------------------------------------
     private:
//------------------------------------------------------------------------------

         void
         map_node_to_bspline_same_mesh(
                 const moris_index        aSourceIndex,
                 const moris_index        aTargetIndex,
                 const enum EntityRank    aBSplineRank );

//------------------------------------------------------------------------------

         void
         map_node_to_element_same_mesh(
                          const moris_index   aSourceIndex,
                          const moris_index   aTargetIndex );

//------------------------------------------------------------------------------

         void
         map_bspline_to_node_same_mesh(
                 const moris_index     aSourceIndex,
                 const enum EntityRank aBSplineRank,
                 const moris_index     aTargetIndex );

//------------------------------------------------------------------------------

         void
         create_iwg_and_model( const real aAlpha = 0.0 );

//------------------------------------------------------------------------------

         void
         create_nodes_for_filter();

//------------------------------------------------------------------------------

         void
         calculate_filter_weights( const real & aFilterRadius );

//------------------------------------------------------------------------------
     };

//------------------------------------------------------------------------------
    } /* namespace mtk */
} /* namespace moris */


#endif /* PROJECTS_MTK_MAP_SRC_CL_MTK_MAPPER_HPP_ */
