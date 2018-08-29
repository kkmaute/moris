/*
 * fn_write_element_ownership_as_field.hpp
 *
 *  Created on: May 24, 2018
 *      Author: ktdoble
 */

#ifndef SRC_MESH_FN_WRITE_ELEMENT_OWNERSHIP_AS_FIELD_HPP_
#define SRC_MESH_FN_WRITE_ELEMENT_OWNERSHIP_AS_FIELD_HPP_

#include "mesh/cl_Mesh_Data.hpp"
#include "xtk/cl_XTK_Mesh.hpp"
#include "xtk/cl_XTK_Cut_Mesh.hpp"
#include "containers/cl_XTK_Cell.hpp"
namespace xtk
{

/*
 * Write the elemental owner data as a field in the output mesh, this function writes this data on the output mesh taking into account the
 * children elements created in the XTK decomposition. For visualization purposes the elemental field needs to be a real type number rather
 * than an integer. The field with the name "aOwnerFieldName" needs to be declared in the output options prior to
 * call to get_output_mesh in the XTK model.
 */

template<typename Real, typename Integer, typename Real_Matrix, typename Integer_Matrix>
void
write_element_ownership_as_field(std::string aOwnerFieldName,
                                 XTK_Mesh<Real,Integer,Real_Matrix,Integer_Matrix> & aXTKMesh,
                                 Cut_Mesh<Real,Integer,Real_Matrix,Integer_Matrix> & aCutMesh,
                                 mesh::Mesh_Data<Real,Integer,Real_Matrix,Integer_Matrix> & aOutputMesh
                                 )
{
    // Get the background mesh
    mesh::Mesh_Data<Real, Integer, Real_Matrix, Integer_Matrix> const & tBackgroundMesh = aXTKMesh.get_mesh_data();

    // Initialize Data (using output mesh because this knows the total number of elements
    Integer tNumElementsOutput = aOutputMesh.get_num_entities(EntityRank::ELEMENT);
    Cell<Real> tOwnerData(tNumElementsOutput);

    // Iterate through background mesh elements (i here corresponds to elemental index)
    for(Integer i = 0; i<tBackgroundMesh.get_num_entities(EntityRank::ELEMENT); i++)
    {
        // Process that owns the element
        Integer tElementOwner = tBackgroundMesh.get_entity_parallel_owner_rank(i, EntityRank::ELEMENT);

        // Get the element Id (needed to translate between background and output mesh)
        Integer tElementId = tBackgroundMesh.get_glb_entity_id_from_entity_loc_index(i,EntityRank::ELEMENT);

        // Check to see if this element has any children
        if(aXTKMesh.entity_has_children(i,EntityRank::ELEMENT))
        {
            // The location of the child mesh in the cut mesh
            Integer tChildMeshIndex = aXTKMesh.child_mesh_index(i,EntityRank::ELEMENT);

            // Retrieve all the element Ids of the children
            moris::Mat_New<Integer, Integer_Matrix> const & tElementIds = aCutMesh.get_element_ids(tChildMeshIndex);

            //Iterate through children elements and ask the output mesh for the indices using ids.
            // The index is then used to place the data in the correct location of tOwnerData.
            for(Integer j = 0; j<tElementIds.n_cols(); j++ )
            {
                // Get element index in output mesh using element Id

                Integer tElementIndex = aOutputMesh.get_loc_entity_index_from_entity_glb_id(tElementIds(0,j), EntityRank::ELEMENT);

                // Add to owner data
                tOwnerData(tElementIndex) = (Real) tElementOwner;

            }
        }

        // No children elements case
        else
        {

            // Get element index in output mesh using element Id
            Integer tElementIndex = aOutputMesh.get_loc_entity_index_from_entity_glb_id(tElementId, EntityRank::ELEMENT);

            // Add to owner data
            tOwnerData(tElementIndex) = (Real) tElementOwner;

        }

    }

    // Write the data to the mesh
    aOutputMesh.add_mesh_field_data_loc_indices(aOwnerFieldName, EntityRank::ELEMENT, tOwnerData);
}

}




#endif /* SRC_MESH_FN_WRITE_ELEMENT_OWNERSHIP_AS_FIELD_HPP_ */
