/*
 * cl_GE_Discrete.hpp
 *
 *  Created on: Mar 21, 2019
 *      Author: sonne
 */

#ifndef PROJECTS_GEN_SRC_CL_GE_DISCRETE_HPP_
#define PROJECTS_GEN_SRC_CL_GE_DISCRETE_HPP_

// GE includes
#include "cl_GE_Geometry.hpp"

// MTK includes
//#include "cl_MTK_Mesh.hpp"
//#include "cl_Mesh_Enums.hpp"

// LINALG includes
#include "cl_Matrix.hpp"


namespace moris
{
namespace ge
{
    class Discrete : public Geometry
    {
        /*
         * note: this is taken from XTK's discrete geometry class, currently merging into a general Geometry Engine
         */
    public:
        Discrete()
        {
            mMySpaceInterpType  = fem::Interpolation_Type::LAGRANGE;
            mMySpaceInterpOrder = mtk::Interpolation_Order::LINEAR;
            mMyTimeInterpType   = fem::Interpolation_Type::CONSTANT;
            mMyTimeInterpOrder  = mtk::Interpolation_Order::CONSTANT;
        };

        ~Discrete(){};
        //------------------------------------------------------------------------------
        /*
         * @brief function to report geometry representation type
         */
        enum GeomType
        get_geom_type() const
        {
            return GeomType::DISCRETE;
        }
        //------------------------------------------------------------------------------
        /*
         * @brief set associated mesh
         */
        void
        set_my_mesh(mtk::Mesh_Manager* aMyMesh)
        {
            MORIS_ASSERT( mMyMesh != nullptr, "ge::Geometry::get_my_mesh(): the associated mesh has not been set" );
            mMyMesh = aMyMesh;
        }
        //------------------------------------------------------------------------------
        void
        set_my_interpolation_rules( fem::Interpolation_Type  aSpaceInterpType,
                                    mtk::Interpolation_Order aSpaceInterpOrder,
                                    fem::Interpolation_Type  aTimeInterpType  = fem::Interpolation_Type::CONSTANT,
                                    mtk::Interpolation_Order aTimeInterpOrder = mtk::Interpolation_Order::CONSTANT )
        {
            mMySpaceInterpType  = aSpaceInterpType;
            mMySpaceInterpOrder = aSpaceInterpOrder;
            mMyTimeInterpType   = aTimeInterpType;
            mMyTimeInterpOrder  = aTimeInterpOrder;
        }
        //------------------------------------------------------------------------------
        void
        set_my_target_field( std::shared_ptr< hmr::Field > &aField )
        {
            MORIS_ASSERT(mMyMesh->get_interpolation_mesh(0)->get_mesh_type() == MeshType::HMR, "ge::Discrete::set_my_target_field(): currently only set for use with an hmr field");
            mMyTargetField = aField;
        }
        //------------------------------------------------------------------------------
        void
        set_my_output_field( std::shared_ptr< hmr::Field > &aField )
        {
            MORIS_ASSERT(mMyMesh->get_interpolation_mesh(0)->get_mesh_type() == MeshType::HMR, "ge::Discrete::set_my_target_field(): currently only set for use with an hmr field");
            mMyOutputField = aField;
        }
        //------------------------------------------------------------------------------
        std::shared_ptr< hmr::Field >
        get_my_target_field()
        {
            MORIS_ASSERT( mMyTargetField != nullptr, "ge::Geometry::get_my_target_field(): target field not set, make sure to set the field before using geometry rep" );
            return mMyTargetField;
        }
        //------------------------------------------------------------------------------
        std::shared_ptr< hmr::Field >
        get_my_output_field()
        {
            MORIS_ASSERT( mMyOutputField != nullptr, "ge::Geometry::get_my_output_field(): output field not set, make sure to set the field before using geometry rep" );
            return mMyOutputField;
        }
        //------------------------------------------------------------------------------
        mtk::Mesh_Manager*
        get_my_mesh()
        {
            return mMyMesh;
        }
        //------------------------------------------------------------------------------
        fem::Interpolation_Type
        get_my_space_interpolation_type()
        {
            return mMySpaceInterpType;
        }
        //------------------------------------------------------------------------------
        mtk::Interpolation_Order
        get_my_space_interpolation_order()
        {
            return mMySpaceInterpOrder;
        }
        //------------------------------------------------------------------------------
        fem::Interpolation_Type
        get_my_time_interpolation_type()
        {
            return mMyTimeInterpType;
        }
        //------------------------------------------------------------------------------
        mtk::Interpolation_Order
        get_my_time_interpolation_order()
        {
            return mMyTimeInterpOrder;
        }
        //------------------------------------------------------------------------------
        /*
         * @brief function to set the member variables
         *
         * param[in] - aMeshWithLevelSetFields      mtk mesh pointer with fields
         * param[in] - aFieldNames                  cell of field names
         */
        void
        set_member_variables(moris::mtk::Mesh_Manager*         aMeshWithLevelSetFields,
                             Cell<std::string> const & aFieldNames)
        {
            mNumLevelSets        = aFieldNames.size();
            mActiveLevelSetIndex = 0;
            mLevelSetFieldNames  = aFieldNames;
            mLevelSetMesh        = aMeshWithLevelSetFields->get_interpolation_mesh( mActiveLevelSetIndex );
        }

        //------------------------------------------------------------------------------
        /**
         * This assumes you are working with the active level set mesh
         */
        moris::real
        access_field_value_with_entity_index(moris::moris_index aEntityIndex,
                                             enum EntityRank    aEntityRank) const
        {
            MORIS_ASSERT(aEntityRank==EntityRank::NODE,"Only nodal levelset values are supported");
            std::string const & tActiveFieldName = get_active_level_set_field_name();
            return mLevelSetMesh->get_entity_field_value_real_scalar({{aEntityIndex}}, tActiveFieldName, (moris::EntityRank)aEntityRank)(0,0);
        }

        //------------------------------------------------------------------------------
        moris::Matrix< moris::DDRMat >
        evaluate_sensitivity_dphi_dp(moris::Matrix< moris::DDRMat > const & aLocalCoordinate, moris::size_t aEntityIndex, enum EntityRank aEntityRank)
        {
            //TODO: Implement this function
            moris::Matrix< moris::DDRMat > tSensitivityDxDp(1,1,0);
            std::cout<<"evaluate_sensitivity_dx_dp function is not implemented in level set mesh";
            return tSensitivityDxDp;
        }

        //------------------------------------------------------------------------------
        std::string const & get_active_level_set_field_name() const
        {
            return mLevelSetFieldNames(mActiveLevelSetIndex);
        }

        //------------------------------------------------------------------------------
        moris::mtk::Mesh*  get_level_set_mesh()
        {
            return mLevelSetMesh;
        }

        //------------------------------------------------------------------------------
    private:
        moris::Cell< real > mMyConstants = 0.0;
        mtk::Mesh_Manager*  mMyMesh = nullptr;
        std::shared_ptr< hmr::Field > mMyTargetField = nullptr;   // currently only set up for hmr::field, need generalized field class?
        std::shared_ptr< hmr::Field > mMyOutputField = nullptr;   // reiterate above; note that this is the field resulting from the map of the target field

        fem::Interpolation_Type  mMySpaceInterpType;
        mtk::Interpolation_Order mMySpaceInterpOrder;
        fem::Interpolation_Type  mMyTimeInterpType;
        mtk::Interpolation_Order mMyTimeInterpOrder;

        // don't think these are necessary anymore?
        moris::size_t     mNumLevelSets;
        moris::size_t     mActiveLevelSetIndex;
        Cell<std::string> mLevelSetFieldNames;
        moris::mtk::Mesh* mLevelSetMesh;

        //------------------------------------------------------------------------------
    protected:

        //------------------------------------------------------------------------------
    };
} /* namespace ge */
} /* namespace moris */



#endif /* PROJECTS_GEN_SRC_CL_GE_DISCRETE_HPP_ */
