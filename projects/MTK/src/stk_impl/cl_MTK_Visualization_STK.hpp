/*
 * cl_MTK__VISUALIZATION_STK.hpp
 *
 *  Created on: Jun 6, 2019
 *      Author: doble
 */

#ifndef PROJECTS_MTK_SRC_STK_IMPL_CL_MTK__VISUALIZATION_STK_HPP_
#define PROJECTS_MTK_SRC_STK_IMPL_CL_MTK__VISUALIZATION_STK_HPP_

#include "cl_Matrix.hpp"
#include "cl_Cell.hpp"
#include "cl_MTK_Scalar_Field_Info.hpp"
#include "typedefs.hpp"


namespace moris
{
namespace mtk
{

class Visualization_STK
{
public:
    Visualization_STK()
    {

    }

    /*!
     * Declare fields for cell ownership, and cells on proc to show up in exodus file.
     * The outputted mtk fields info must be provided to the create mesh call after setup.
     */
    moris::mtk::MtkFieldsInfo*
    setup_parallel_cell_fields_for_declaration()
    {

        // declare cell owner field
        this->setup_cell_owner_field_for_declaration();

        // declare cells on proc field
        this->setup_cells_on_proc_field_for_declaration();

        return &mFieldMtkFieldsInfo;

    }

    void
    setup_cell_owner_field_for_declaration()
    {
        // declare cell owner field
        std::string tOwnerCellField = "cell_owner";
        mCellOwnersFields.set_field_name(tOwnerCellField);
        mCellOwnersFields.set_field_entity_rank(EntityRank::ELEMENT);
        add_field_for_mesh_input(&mCellOwnersFields,mFieldMtkFieldsInfo);
    }

    void
    setup_cells_on_proc_field_for_declaration()
    {
        // allocate cells on procs fields
        moris::uint tParSize = par_size();
        mCellsOnProcFields.resize(tParSize);

        for(moris::uint  i = 0; i < tParSize; i++)
        {
            std::string tProcsOnFieldName = "cells_on_proc_"+ std::to_string(i);
            mCellsOnProcFields(i).set_field_name(tProcsOnFieldName);
            mCellsOnProcFields(i).set_field_entity_rank(EntityRank::ELEMENT);
            add_field_for_mesh_input(&mCellsOnProcFields(i),mFieldMtkFieldsInfo);
        }
    }

    /*
     * After mesh has been setup with the fields declared in setup_parallel_cell_fields_for_declaration
     * then populate the information on the mesh with this call
     */
    void
    populate_parallel_cell_fields_on_mesh(moris::mtk::Mesh* aMesh)
    {
        // populate the cell owner field
        this->populate_cell_owner_field_on_mesh(aMesh);

        // populate the cells on procs field
        this->populate_cells_on_proc_field_on_mesh(aMesh);
    }

    void
    populate_cell_owner_field_on_mesh(moris::mtk::Mesh* aMesh)
    {
        moris::uint tNumCells = aMesh->get_num_entities(EntityRank::ELEMENT);
        moris::Matrix<moris::DDRMat> tCellOwnerField(tNumCells,1);

        // iterate through cells and get owner
        for(moris::uint  i = 0; i < tNumCells; i++)
        {
            tCellOwnerField(i) = aMesh->get_entity_owner(i,EntityRank::ELEMENT);
        }

        // add to mesh
        aMesh->add_mesh_field_real_scalar_data_loc_inds(mCellOwnersFields.get_field_name(),EntityRank::ELEMENT,tCellOwnerField);
    }

    void
    populate_cells_on_proc_field_on_mesh(moris::mtk::Mesh* aMesh)
    {
        moris::uint tNumCells = aMesh->get_num_entities(EntityRank::ELEMENT);

        moris::Cell<moris::Cell<moris::real>> tCellsOnProcs(par_size());

        moris::uint tMyRank = par_rank();

        for(moris::uint i = 0; i < tNumCells; i++)
        {
            moris::Matrix<moris::IndexMat> tProcsSharing;
            aMesh->get_processors_whom_share_entity(i,EntityRank::ELEMENT,tProcsSharing);

            // always add cell to cells on the current proc
            tCellsOnProcs(tMyRank).push_back(i);

            // iterate through sharing procs
            for(moris::uint j = 0; j < tProcsSharing.numel(); j++ )
            {
                if(tProcsSharing(j) != (moris_id)tMyRank)
                {
                    // add cell index as being on other processor
                    tCellsOnProcs(tProcsSharing(j)).push_back(i);
                }
            }
        }

        // add field to mesh
        moris::Matrix<moris::DDRMat> tFlagCellsOnProc(tNumCells,1);

        for(moris::uint  i = 0; i <(uint) par_size(); i++)
        {
            if(tCellsOnProcs.size() > 0)
            {
                // reset the flag matrix
                tFlagCellsOnProc.fill(0);

                //iterate through cells in tCellsOnProc for this procs
                for(moris::uint j = 0; j < tCellsOnProcs(i).size(); j++)
                {
                    tFlagCellsOnProc(tCellsOnProcs(i)(j)) = 1;
                }

                aMesh->add_mesh_field_real_scalar_data_loc_inds(mCellsOnProcFields(i).get_field_name(),EntityRank::ELEMENT,tFlagCellsOnProc);

            }
        }

    }


private:
    moris::mtk::MtkFieldsInfo mFieldMtkFieldsInfo;
    Scalar_Field_Info<DDRMat> mCellOwnersFields;
    moris::Cell< Scalar_Field_Info<DDRMat>> mCellsOnProcFields;

};
}
}

#endif /* PROJECTS_MTK_SRC_STK_IMPL_CL_MTK__VISUALIZATION_STK_HPP_ */
