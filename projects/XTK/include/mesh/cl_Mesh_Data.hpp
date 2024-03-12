/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_Mesh_Data.hpp
 *
 */

#ifndef INCLUDE_MESH_CL_MESH_DATA_HPP_
#define INCLUDE_MESH_CL_MESH_DATA_HPP_

#include "mesh/cl_MTK_Enums.hpp"
#include "containers/cl_XTK_Cell.hpp"
#include "xtk/cl_XTK_Pending_Node.hpp"
#include "linalg/cl_XTK_Matrix.hpp"

// TODO: Move these with get_coordinate_field and get field
#include <stk_mesh/base/Field.hpp>                // for Field
#include "stk_mesh/base/FieldBase.hpp"            // for field_data, etc
#include <stk_mesh/base/CoordinateSystems.hpp>    // for Cartesian

using namespace xtk;
namespace mesh
{
    template< typename Real, typename Integer, typename Real_Matrix, typename Integer_Matrix >
    class Mesh_Data
    {
      public:
        virtual ~Mesh_Data()
        {
        }

        virtual moris::Matrix< Integer_Matrix >
        get_entity_connected_to_entity_loc_inds( Integer aEntityIndex,
                mtk::EntityRank                          aInputEntityRank,
                mtk::EntityRank                          aOutputEntityRank ) const = 0;

        // Since some mesh implementations do not support connectivities aInputEntityRank == aOutputEntityRank, XTK only needs Element to element
        // connectivity therefore this is implemented here

        /**
         *
         * @param[in]  aElementIndex - element index to get the neighbors for
         * @param[out] Matrix of size 2xNum Neighbors where the first row is
         *             the neighbor and the second is the face index shared
         *             by the elements
         */
        virtual moris::Matrix< Integer_Matrix >
        get_elements_connected_to_element_loc_inds( Integer aElementIndex ) const = 0;

        /*
         * For outputtingg
         */
        virtual moris::Matrix< Integer_Matrix >
        get_entity_connected_to_entity_glb_ids( Integer aEntityIndex,
                mtk::EntityRank                         aInputEntityRank,
                mtk::EntityRank                         aOutputEntityRank ) const = 0;

        virtual Integer get_element_face_ordinal_loc_inds( Integer const &aElementIndex,
                Integer const                                            &aFaceIndex ) const                             = 0;
        virtual Integer get_num_entities( mtk::EntityRank aEntityRank ) const = 0;

        virtual moris::Matrix< Real_Matrix > get_selected_node_coordinates_loc_inds( moris::Matrix< Integer_Matrix > const &aNodeIndices ) const = 0;

        virtual moris::Matrix< Real_Matrix > get_all_node_coordinates_loc_inds() const = 0;

        //    virtual moris::Matrix< Integer_Matrix > get_all_entity_indices(mtk::EntityRank aEntityRank) const  = 0;

        virtual Integer get_glb_entity_id_from_entity_loc_index( Integer aEntityIndex, mtk::EntityRank aEntityRank ) const = 0;

        virtual Integer get_loc_entity_index_from_entity_glb_id( Integer aEntityId, mtk::EntityRank aEntityRank ) const = 0;

        virtual moris::Matrix< Integer_Matrix > const &get_local_to_global_map( mtk::EntityRank aEntityRank ) const = 0;

        virtual moris::Matrix< Real_Matrix > get_entity_field_value( moris::Matrix< Integer_Matrix > const &aEntityIndex, std::string const &aFieldName, mtk::EntityRank aFieldEntityRank ) const = 0;

        virtual Real get_entity_field_value( Integer const &aEntityIndex, std::string const &aFieldName, mtk::EntityRank aFieldEntityRank ) = 0;

        virtual Integer get_entity_field_value_integer( Integer const &aEntityIndex, std::string const &aFieldName, mtk::EntityRank aFieldEntityRank ) = 0;

        // Parallel Implementation
        virtual Integer get_entity_parallel_owner_rank( Integer aEntityIndex, mtk::EntityRank aEntityRank ) const = 0;

        virtual Integer get_num_of_entities_shared_with_processor( Integer aProcessorRank, mtk::EntityRank aEntityRank, bool aSendFlag ) const = 0;

        virtual void get_processors_whom_share_entity( Integer aEntityIndex,
                mtk::EntityRank                                aEntityRank,
                moris::Matrix< Integer_Matrix >               &aProcsWhomShareEntity ) const = 0;

        // NOTE: this should not change the underlying meshes data (i.e. do not use stk's internal function because it alters certain bits of data internally while I want it to be constant)
        virtual Integer allocate_entity_ids( Integer aNumIdstoAllocate,
                mtk::EntityRank                      aEntityRank ) const = 0;

        virtual void add_mesh_field_data_loc_indices( std::string const &aFieldName,
                mtk::EntityRank const                                   &aFieldEntityRank,
                xtk::Vector< Real > const                                 &aFieldData ) = 0;

        virtual void add_mesh_field_data_loc_indices_integer( std::string const &aFieldName,
                mtk::EntityRank const                                           &aFieldEntityRank,
                xtk::Vector< Integer > const                                      &aFieldData ) = 0;

        virtual void write_output_mesh( std::string const &aMeshName,
                xtk::Vector< std::string > const            &aRealNodeFieldsToOutput       = {},
                xtk::Vector< std::string > const            &aIntNodeFieldsToOutput        = {},
                xtk::Vector< std::string > const            &aRealElementFieldsToOutput    = {},
                xtk::Vector< std::string > const            &aIntElementFieldsToOutput     = {},
                xtk::Vector< std::string > const            &aRealVectorNodeFieldsToOutput = {},
                Integer                                    aTime                         = 0 ) const = 0;

        /*
         * EVERYTHING BELOW THIS POINT IS RELATED TO PACKAGING MESH DATA FOR OUTPUTTING OF THE CONFORMAL MESH WHICH HAS THE SAME BLOCKS AS THE PREVIOUS MESH
         */

        virtual void get_all_part_names( mtk::EntityRank const &aPrimaryRank,
                xtk::Vector< std::string >                       &aPartNames ) const = 0;

        /*
         * Get buckets, where a requirement  of a bucket is have the same entity rank and the same parts.
         */

        virtual Integer get_num_buckets( mtk::EntityRank aEntityRank ) const = 0;

        virtual moris::Matrix< Integer_Matrix > get_entities_in_bucket_loc_index( Integer aBucketOrdinal, mtk::EntityRank aEntityRank ) const = 0;

        /*
         * This should only return not internal part membership
         * For example using STK, there should not be universal, globally shared, locally owned, or aura parts included in the returned part ordinal
         */
        virtual void get_entity_part_membership_ordinals( Integer const &aEntityIndex, mtk::EntityRank aEntityRank, xtk::Vector< Integer > &aPartOrdinal, bool const &aInduced = false ) const = 0;

        /*
         * This function returns the string associated with a part which are accessed using part ordinal
         */
        virtual void get_part_name_from_part_ordinals( xtk::Vector< Integer > const &aPartOrdinals, xtk::Vector< std::string > &aPartNames ) const = 0;

        // Basis function accessing functions --------------------------------------------
        // For Lagrangian type meshes where nodes and basis functions coincide these are the same as the accessing node functions above

        /*
         * Returns the number of basis functions in the mesh
         */
        virtual Integer get_num_basis_functions() const = 0;

        /*
         *  Returns the elements in the support of a the basis function with index (aBasisIndex)
         */
        virtual moris::Matrix< Integer_Matrix > get_elements_in_basis_support( Integer aBasisIndex ) const = 0;

        // TODO: Figure out a clean way to access this information
        /*
         * DO NOT USE THESE IN CODE THEY ARE HERE RIGHT NOW PROVIDE EASY ACCESS TO WORK ON RAW STK DATA IN TEST CASES
         */
        //    virtual stk::mesh::Field<Real, stk::mesh::Cartesian3d> * get_coordinate_field() {};
        //    virtual stk::mesh::FieldBase * get_field(mtk::EntityRank aFieldRank,
        //                                             std::string aFieldName) {};
        //    virtual stk::mesh::BulkData & mesh_bulk_data()=0;
        //    virtual stk::mesh::MetaData & mesh_meta_data()=0;
    };
}    // namespace mesh

#endif /* INCLUDE_MESH_CL_MESH_DATA_HPP_ */
