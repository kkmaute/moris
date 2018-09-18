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

#include "../linalg/cl_XTK_Matrix_Base_Utilities.hpp" // For print
#include "linalg/cl_XTK_Matrix.hpp"
#include "geometry/cl_Geometry.hpp"

// XTKL: Linear Algebra Includes

#include "assert/fn_xtk_assert.hpp"

// XTKL: Containers
#include"containers/cl_XTK_Cell.hpp"

//XTKL: Mesh Interface Includes
#include "mesh/cl_Mesh_Data.hpp"
#include "mesh/cl_Mesh_Enums.hpp"
#include "mesh/cl_Mesh_Builder.hpp"

//XTKL: Topology
#include "topology/cl_XTK_Topology.hpp"
#include "topology/cl_XTK_Basis_Function.hpp"

namespace xtk
{
template<typename Real, typename Integer, typename Real_Matrix, typename Integer_Matrix>
class Discrete_Level_Set: public Geometry<Real, Integer, Real_Matrix, Integer_Matrix>
{
public:
    Discrete_Level_Set()
    {

    }
    Discrete_Level_Set(Cell<Geometry<Real, Integer, Real_Matrix, Integer_Matrix>*> & aLevelSetFunctionsToDiscretize,
                   std::string const & aMeshFile,
                   Cell<std::string> const & aFieldNames,
                   mesh::Mesh_Builder<Real, Integer, Real_Matrix, Integer_Matrix> & aMeshBuilder) :
    mNumLevelSets(aLevelSetFunctionsToDiscretize.size()), mActiveLevelSetIndex(0), mLevelSetFieldNames(aFieldNames)
    {
        mLevelSetMesh = aMeshBuilder.build_mesh_from_string(aMeshFile, aFieldNames, false);
        this->discretize_level_set_functions(aLevelSetFunctionsToDiscretize);
    }

    Discrete_Level_Set( std::shared_ptr<mesh::Mesh_Data<Real, Integer, Real_Matrix, Integer_Matrix>> & aMeshWithLevelSetFields,
                        Cell<std::string> const & aFieldNames) :
    mNumLevelSets(aFieldNames.size()), mActiveLevelSetIndex(0), mLevelSetFieldNames(aFieldNames),  mLevelSetMesh(aMeshWithLevelSetFields)
    {

    }


    bool is_analytic() const
    {
        return false;
    }

    void get_dphi_dp_size(Integer & aNumRows, Integer & aNumCols) const
    {
        aNumRows = 2;
        aNumCols = 3;
    }


    moris::Matrix< Integer_Matrix >
    get_node_adv_indices(moris::Matrix< Integer_Matrix > const & aNodeIndices)
    {
        Integer tNumADVS = 2;
        moris::Matrix< Integer_Matrix > tADVIndices(1,tNumADVS);

        for(Integer i = 0; i<tNumADVS; i++)
        {
            tADVIndices(0,i) = mLevelSetMesh->get_glb_entity_id_from_entity_loc_index(aNodeIndices(0,i),EntityRank::NODE);
        }

        return tADVIndices;
    }


    /**
     * This assumes you are working with the active level set mesh
     */
    Real access_field_value_with_entity_index(Integer aEntityIndex, enum EntityRank aEntityRank) const
    {
        XTK_ASSERT(aEntityRank==EntityRank::NODE,"Only nodal levelset values are supported");
        std::string const & tActiveFieldName = get_active_level_set_field_name();
        return mLevelSetMesh->get_entity_field_value(aEntityIndex, tActiveFieldName, aEntityRank);
    }

    moris::Matrix< Real_Matrix > evaluate_sensitivity_dx_dp(moris::Matrix< Real_Matrix > const & aLocalCoordinate, Integer aEntityIndex, enum EntityRank aEntityRank)
    {
        //TODO: Implement this function
        moris::Matrix< Real_Matrix > tSensitivityDxDp(1,1,0);
        XTK_ERROR<<"evaluate_sensitivity_dx_dp function is not implemented in level set mesh";
        return tSensitivityDxDp;
    }



    std::string const & get_active_level_set_field_name() const
    {
        return mLevelSetFieldNames(mActiveLevelSetIndex);
    }

    Integer get_num_levelset() const
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

    std::string const & get_level_set_field_name(Integer aLevelSetIndex) const
    {
        XTK_ASSERT(aLevelSetIndex < mNumLevelSets, "Requested level set field name is outside of bounds");
        return mLevelSetFieldNames(aLevelSetIndex);
    }

    Cell<std::string> const & get_level_set_field_name() const
    {
        return mLevelSetFieldNames;
    }

    std::shared_ptr<mesh::Mesh_Data<Real, Integer, Real_Matrix, Integer_Matrix>>  get_level_set_mesh()
    {
        return mLevelSetMesh;
    }

private:
    Integer mNumLevelSets;
    Integer mActiveLevelSetIndex;
    Cell<std::string> mLevelSetFieldNames;
    //TODO: Accomodate multiple level set meshes
    // Right now keeping it simple and assuming all the level set data is applied as field on the same mesh
//    Cell<std::shared_ptr<mesh::Mesh_Data<Real,Integer>>> mLevelSetMesh;
    std::shared_ptr<mesh::Mesh_Data<Real, Integer, Real_Matrix, Integer_Matrix>> mLevelSetMesh;


private:
    void discretize_level_set_functions(Cell<Geometry<Real, Integer, Real_Matrix, Integer_Matrix>*> & aLevelSetFunctionsToDiscretize)
    {
        // Get information about number of nodes and their coordinates
        // Split into two loops to avoid rewriting add_mesh_field_data function and to collect all field data first then apply to mesh
        Integer tNumNodes = mLevelSetMesh->get_num_entities(EntityRank::NODE);
        moris::Matrix< Real_Matrix > tCoordinates = mLevelSetMesh->get_all_node_coordinates_loc_inds();
        moris::Matrix< Real_Matrix > tNodeCoordinates(1, 3);
        Cell < Cell < Real >> tFieldData(mNumLevelSets, tNumNodes);

        for (Integer i = 0; i < mNumLevelSets; i++)
        {
            XTK_ASSERT(aLevelSetFunctionsToDiscretize(i)->is_analytic(),"Cannot discretize a non-analytic level set field");
            for (Integer n = 0; n < tNumNodes; n++)
            {
                tFieldData(i)(n) = aLevelSetFunctionsToDiscretize(i)->evaluate_field_value_with_coordinate(n,tCoordinates);
            }
        }

        for (Integer i = 0; i < mNumLevelSets; i++)
        {
            mLevelSetMesh->add_mesh_field_data_loc_indices(mLevelSetFieldNames(i), EntityRank::NODE, tFieldData(i));
        }
    }

};
}

#endif /* SRC_MESH_CL_MESH_LEVEL_SET_MESH_HPP_ */
