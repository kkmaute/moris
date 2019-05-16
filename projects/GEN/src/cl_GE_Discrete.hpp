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
        Discrete(){};

        ~Discrete(){};
        //------------------------------------------------------------------------------
        /*
         * @brief function to report if geometry representation type
         */
        enum type get_geom_type() const
        {
            return type::DISCRETE;
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
