#ifndef MORIS_CL_GEN_DISCRETE_LEVEL_SET_HPP
#define MORIS_CL_GEN_DISCRETE_LEVEL_SET_HPP

#include "cl_GEN_Geometry_Discrete.hpp"
#include "cl_Matrix.hpp"
#include "cl_Cell.hpp"

namespace moris
{
    namespace ge
    {
        class Discrete_Level_Set: public Geometry_Discrete
        {
        private:
            size_t mActiveFieldIndex = 0;
            Cell<std::string> mFieldNames;
            Interpolation_Mesh* mMesh;
            EntityRank mEntityRank;

        public:
            /**
             * Constructor
             *
             * @param aMeshWithLevelSetFields Mesh with the level set fields
             * @param aFieldNames Names of the fields
             */
            Discrete_Level_Set(moris::mtk::Interpolation_Mesh* aMeshWithLevelSetFields,
                    moris::Cell<std::string> const & aFieldNames, aEntityRank = EntityRank::NODE);

            /**
             * Given an index, the discrete geometry needs to return a field value.
             *
             * @param aEntityIndex the index of the field value
             * @return field value at the specified index
             */
            real evaluate_field_value(moris_index aEntityIndex);

            /**
             * Given an index, the discrete geometry needs to return sensitivites with respect to the field value
             *
             * @param aEntityIndex the index of the field value
             * @return matrix of sensitivities
             */
            Matrix<DDRMat> evaluate_sensitivity(moris_index aEntityIndex);

            /**
             * Sets the active field to the given value
             *
             * @param aActiveFieldIndex The index of the desired field to become active
             */
            void set_active_field(size_t aActiveFieldIndex);


//            moris::Matrix< moris::IndexMat >
//            get_node_adv_indices(moris::Matrix< moris::IndexMat > const & aNodeIndices)
//            {
//                moris::size_t tNumADVS = 2;
//                moris::Matrix< moris::IndexMat > tADVIndices(1,tNumADVS);
//
//                for(moris::size_t i = 0; i<tNumADVS; i++)
//                {
//                    tADVIndices(0,i) = mLevelSetMesh->get_glb_entity_id_from_entity_loc_index(aNodeIndices(0,i),moris::EntityRank::NODE);
//                }
//
//                return tADVIndices;
//            }
//
//
//            /**
//             * This assumes you are working with the active level set mesh
//             */
//            moris::real
//            access_field_value_with_entity_index(moris::moris_index aEntityIndex,
//                                                 enum moris::EntityRank    aEntityRank) const
//            {
//                MORIS_ASSERT(aEntityRank==moris::EntityRank::NODE,"Only nodal levelset values are supported");
//                std::string const & tActiveFieldName = get_active_level_set_field_name();
//                return mLevelSetMesh->get_entity_field_value_real_scalar({{aEntityIndex}}, tActiveFieldName, (moris::EntityRank)aEntityRank)(0,0);
//            }
//
//            moris::Matrix< moris::DDRMat > evaluate_sensitivity_dx_dp(moris::Matrix< moris::DDRMat > const & aLocalCoordinate, moris::size_t aEntityIndex, enum moris::EntityRank aEntityRank)
//            {
//                //TODO: Implement this function
//                moris::Matrix< moris::DDRMat > tSensitivityDxDp(1,1,0);
//                std::cout<<"evaluate_sensitivity_dx_dp function is not implemented in level set mesh";
//                return tSensitivityDxDp;
//            }
//
//
//
//            std::string const & get_active_level_set_field_name() const
//            {
//                return mLevelSetFieldNames(mActiveLevelSetIndex);
//            }
//
//            moris::size_t get_num_levelset() const
//            {
//                return mNumLevelSets;
//            }
//
//            bool advance_to_next_level_set()
//            {
//                bool tAnotherLevelset;
//                if(mActiveLevelSetIndex == get_num_levelset()-1)
//                {
//                    mActiveLevelSetIndex++;
//                    tAnotherLevelset = false;
//                }
//                else
//                {
//                    mActiveLevelSetIndex++;
//                    tAnotherLevelset = true;
//                }
//
//                return tAnotherLevelset;
//            }
//
//            std::string const & get_level_set_field_name(moris::size_t aLevelSetIndex) const
//            {
//                MORIS_ASSERT(aLevelSetIndex < mNumLevelSets, "Requested level set field name is outside of bounds");
//                return mLevelSetFieldNames(aLevelSetIndex);
//            }
//
//            moris::Cell<std::string> const & get_level_set_field_name() const
//            {
//                return mLevelSetFieldNames;
//            }
//
//            moris::mtk::Interpolation_Mesh*  get_level_set_mesh()
//            {
//                return mLevelSetMesh;
//            }
        };
    }
}

#endif /* MORIS_CL_GEN_DISCRETE_LEVEL_SET_HPP */
