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

        Discrete(moris::mtk::Mesh*         aMeshWithLevelSetFields,
                 Cell<std::string> const & aFieldNames) :
                    mNumLevelSets(aFieldNames.size()),
                    mActiveLevelSetIndex(0),
                    mLevelSetFieldNames(aFieldNames),
                    mLevelSetMesh(aMeshWithLevelSetFields)
                {
                };

        ~Discrete(){};
        //------------------------------------------------------------------------------
        /*
         * @brief function to report if geometry is analytic or not
         */
        bool is_analytic() const
        {
            return false;
        }
        //------------------------------------------------------------------------------
        /*
         * @brief dummy function to set the member variables,
         *              going to merge into a single constructor through the factory
         *
         * param[in] - aMeshWithLevelSetFields      mtk mesh pointer with fields
         * param[in] - aFieldNames                  cell of field names
         */
        void
        set_member_variables(moris::mtk::Mesh*         aMeshWithLevelSetFields,
                             Cell<std::string> const & aFieldNames)
        {
            mNumLevelSets        = aFieldNames.size();
            mActiveLevelSetIndex = 0;
            mLevelSetFieldNames  = aFieldNames;
            mLevelSetMesh        = aMeshWithLevelSetFields;
        }

        //------------------------------------------------------------------------------
        moris::Matrix< moris::IndexMat >
        get_node_adv_indices(moris::Matrix< moris::IndexMat > const & aNodeIndices)
        {
            moris::size_t tNumADVS = 2;
            moris::Matrix< moris::IndexMat > tADVIndices(1,tNumADVS);

            for(moris::size_t i = 0; i<tNumADVS; i++)
            {
                tADVIndices(0,i) = mLevelSetMesh->get_glb_entity_id_from_entity_loc_index(aNodeIndices(0,i),moris::EntityRank::NODE);
            }

            return tADVIndices;
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
        moris::size_t get_num_levelset() const
        {
            return mNumLevelSets;
        }

        //------------------------------------------------------------------------------
        bool advance_to_next_level_set()
        {
            bool tAnotherLevelset;
            if(mActiveLevelSetIndex == get_num_levelset()-1)
            {
                mActiveLevelSetIndex++;
                tAnotherLevelset = false;
            }
            else
            {
                mActiveLevelSetIndex++;
                tAnotherLevelset = true;
            }

            return tAnotherLevelset;
        }

        //------------------------------------------------------------------------------
        std::string const & get_level_set_field_name(moris::size_t aLevelSetIndex) const
        {
            MORIS_ASSERT(aLevelSetIndex < mNumLevelSets, "Requested level set field name is outside of bounds");
            return mLevelSetFieldNames(aLevelSetIndex);
        }

        //------------------------------------------------------------------------------
        Cell<std::string> const & get_level_set_field_name() const
        {
            return mLevelSetFieldNames;
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
