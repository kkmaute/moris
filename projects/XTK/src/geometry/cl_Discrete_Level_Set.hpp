/*
 * cl_Mesh_Level_Set_Mesh.hpp
 *
 *  Created on: Jul 4, 2017
 *      Author: ktdoble
 */

#ifndef SRC_MESH_CL_MESH_LEVEL_SET_MESH_HPP_
#define SRC_MESH_CL_MESH_LEVEL_SET_MESH_HPP_

// Standard includes
#include <string>

#include "cl_XTK_Matrix_Base_Utilities.hpp" // For print
#include "cl_Matrix.hpp"
#include "cl_Geometry.hpp"

// XTKL: Linear Algebra Includes



// XTKL: Containers
#include"cl_Cell.hpp"

//XTKL: Mesh Interface Includes
#include "cl_MTK_Mesh.hpp"
#include "cl_Mesh_Enums.hpp"

//XTKL: Topology
#include "cl_XTK_Topology.hpp"
#include "cl_XTK_Basis_Function.hpp"

namespace xtk
{
class Discrete_Level_Set: public Geometry
{
public:
    Discrete_Level_Set()
    {

    }


    Discrete_Level_Set( moris::mtk::Mesh*         aMeshWithLevelSetFields,
                        Cell<std::string> const & aFieldNames) :
                            mNumLevelSets(aFieldNames.size()),
                            mActiveLevelSetIndex(0),
                            mLevelSetFieldNames(aFieldNames),
                            mLevelSetMesh(aMeshWithLevelSetFields)
    {

    }


    bool is_analytic() const
    {
        return false;
    }

    void get_dphi_dp_size(moris::size_t & aNumRows, moris::size_t & aNumCols) const
    {
        aNumRows = 2;
        aNumCols = 3;
    }


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

    moris::Matrix< moris::DDRMat > evaluate_sensitivity_dx_dp(moris::Matrix< moris::DDRMat > const & aLocalCoordinate, moris::size_t aEntityIndex, enum EntityRank aEntityRank)
    {
        //TODO: Implement this function
        moris::Matrix< moris::DDRMat > tSensitivityDxDp(1,1,0);
        std::cout<<"evaluate_sensitivity_dx_dp function is not implemented in level set mesh";
        return tSensitivityDxDp;
    }



    std::string const & get_active_level_set_field_name() const
    {
        return mLevelSetFieldNames(mActiveLevelSetIndex);
    }

    moris::size_t get_num_levelset() const
    {
        return mNumLevelSets;
    }

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

    std::string const & get_level_set_field_name(moris::size_t aLevelSetIndex) const
    {
        MORIS_ASSERT(aLevelSetIndex < mNumLevelSets, "Requested level set field name is outside of bounds");
        return mLevelSetFieldNames(aLevelSetIndex);
    }

    Cell<std::string> const & get_level_set_field_name() const
    {
        return mLevelSetFieldNames;
    }

    moris::mtk::Mesh*  get_level_set_mesh()
    {
        return mLevelSetMesh;
    }

private:
    moris::size_t     mNumLevelSets;
    moris::size_t     mActiveLevelSetIndex;
    Cell<std::string> mLevelSetFieldNames;
    // TODO: test multiple level set meshes
    // Right now keeping it simple and assuming all the level set data is applied as field on the same mesh
    // Cell<std::shared_ptr<mesh::Mesh_Data<moris::real,moris::size_t>>> mLevelSetMesh;
    moris::mtk::Mesh* mLevelSetMesh;



};
}

#endif /* SRC_MESH_CL_MESH_LEVEL_SET_MESH_HPP_ */
