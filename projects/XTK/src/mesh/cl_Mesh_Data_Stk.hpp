/*
 * cl_Mesh_Data_Stk.hpp
 *
 *  Created on: Jun 20, 2017
 *      Author: ktdoble
 */

#ifndef SRC_MESH_CL_MESH_DATA_STK_HPP_
#define SRC_MESH_CL_MESH_DATA_STK_HPP_

// Std: Standard Includes and MPI
#include <memory>
#include <mpi.h>

// XTKL: Linear Algebra Includes
#include "containers/cl_XTK_Cell.hpp"

// XTKL: Mesh implementations depending on STK
#include "mesh/cl_Mesh_External_Entity_Data_Stk.hpp"

// XTKL: Logging and Assertion Functions
#include "ios/cl_Logger.hpp"
#include "assert/fn_xtk_assert.hpp"

// TPL: STK Includes
#include <stk_io/StkMeshIoBroker.hpp>    // for StkMeshIoBroker
#include <stk_mesh/base/GetEntities.hpp> // for count_entities
#include <stk_mesh/base/BulkData.hpp>    // for BulkData
#include <stk_mesh/base/MetaData.hpp>    // for MetaData
#include <stk_mesh/base/CreateFaces.hpp> // for handling faces
#include <stk_mesh/base/CreateEdges.hpp> // for handling faces
#include <stk_mesh/base/CoordinateSystems.hpp> // for Cartesian

#include "linalg/cl_XTK_Matrix.hpp"
#include "cl_Mesh_Parallel_Data_Stk.hpp"
#include "stk_mesh/base/Field.hpp"       // for Field
#include "stk_mesh/base/FieldBase.hpp"   // for field_data, etc
#include "Ioss_Region.h"                 // for Region, NodeBlockContainer
#include "Ioss_Field.h"                  // for Field

namespace mesh
{
template<typename Real, typename Integer, typename Real_Matrix, typename Integer_Matrix>
class Mesh_Data_Stk: public Mesh_Data<Real, Integer, Real_Matrix, Integer_Matrix>
{
public:

    Mesh_Data_Stk(MPI_Comm aCommunicator, int aParallelPoolSize,
                  stk::mesh::BulkData::AutomaticAuraOption aAuraOption,
                  bool aLoadedFromFile) :
            mLoadedFromFile(aLoadedFromFile),
            mStkMeshMetaData(std::make_shared<stk::mesh::MetaData>(3)),
            mStkMeshBulkData(std::make_shared<stk::mesh::BulkData>((*mStkMeshMetaData), aCommunicator, aAuraOption)),
            mStkMeshIoBroker(std::make_shared<stk::io::StkMeshIoBroker>(aCommunicator)),
            mMeshParallelData(aParallelPoolSize),
            mMeshExternalEntityData()
    {
    }

    // Required implementation for XTK


    moris::Mat_New<Integer, Integer_Matrix>
    get_entity_connected_to_entity_loc_inds(Integer aEntityIndex, enum EntityRank aInputEntityRank, enum EntityRank aOutputEntityRank) const
    {

        XTK_ASSERT(aInputEntityRank != aOutputEntityRank," Input and output entity rank cannot be the same (this is an invalid connectivity inside STK). Use get_element_to_element_loc_inds for element to element connectivity.");

        // Get Stk entity Id from local to global map
        Integer tId = this->get_glb_entity_id_from_entity_loc_index(aEntityIndex, aInputEntityRank);
        stk::mesh::EntityId tStkEntityId = { tId };

        // Call function that gets the connected entities
        stk::mesh::EntityRank tStkInputRank = this->get_stk_entity_rank(aInputEntityRank);
        stk::mesh::EntityRank tStkOutputRank = this->get_stk_entity_rank(aOutputEntityRank);
        stk::mesh::Entity tStkEntity = mStkMeshBulkData->get_entity(tStkInputRank, tStkEntityId);

        moris::Mat_New<Integer, Integer_Matrix> tEntitiesConnected = this->entities_connected_to_entity_stk_glb_id(&tStkEntity, tStkInputRank, tStkOutputRank);
        // Get the number of entities
        Integer tNumOutputEntities = tEntitiesConnected.n_cols();

        // Declare xtk::Mat that will contain the local entities
        Integer tFillValue = 0;
        moris::Mat_New<Integer, Integer_Matrix> tLocalIndices(1, tNumOutputEntities, tFillValue);

        // Fill local ids to xtk::Mat
        for (Integer i = 0; i < tNumOutputEntities; ++i)
        {

            tLocalIndices(0, i) = mStkMeshBulkData->local_id(mStkMeshBulkData->get_entity(tStkOutputRank, tEntitiesConnected(0, i)));
        }
        return tLocalIndices;
    }

    //
    moris::Mat_New<Integer, Integer_Matrix>
    get_element_connected_to_element_loc_inds(Integer aElementIndex) const
    {

        // First get faces connected to element
        // Get Stk entity Id from local to global map
        Integer tId = this->get_glb_entity_id_from_entity_loc_index(aElementIndex, EntityRank::ELEMENT);
        stk::mesh::EntityId tStkEntityId = { tId };

        // Call function that gets the connected entities
        stk::mesh::EntityRank tStkInputRank  = stk::topology::ELEMENT_RANK;
        stk::mesh::EntityRank tStkOutputRank = stk::topology::FACE_RANK;
        stk::mesh::Entity tStkEntity = mStkMeshBulkData->get_entity(tStkInputRank, tStkEntityId);

        moris::Mat_New<Integer, Integer_Matrix> tFacesInElem = this->entities_connected_to_entity_stk_glb_id(&tStkEntity, tStkInputRank, tStkOutputRank);


        XTK_ASSERT( ( tFacesInElem.n_rows() != 0 ) || ( tFacesInElem.n_cols() != 0 ),
                "No faces connected to element found. Maybe the CreateAllEdgesAndFaces flag is set to false. Check mesh struct." );

        // Then for each face get elements connected
        Integer tCounter  = 0;
        Integer tNumFaces = tFacesInElem.n_cols();

        moris::Mat_New<Integer, Integer_Matrix> tElemsConnectedToElem(2, tNumFaces);

        for ( Integer faceIt = 0; faceIt < tNumFaces; ++faceIt )
        {
            Integer tFaceId = tFacesInElem( 0, faceIt );
            stk::mesh::Entity aFaceEntity = mStkMeshBulkData->get_entity( stk::topology::FACE_RANK, tFaceId );
            moris::Mat_New<Integer, Integer_Matrix> tDummyConnectivity = this->entities_connected_to_entity_stk_glb_id( &aFaceEntity,stk::topology::FACE_RANK, stk::topology::ELEMENT_RANK );

            // Faces in mesh boundaries do not have more than one element
            if ( tDummyConnectivity.n_cols() > 0 )
            {
                if ( tDummyConnectivity( 0,0 ) != tId )
                {
                    tElemsConnectedToElem( 0,tCounter ) = get_loc_entity_index_from_entity_glb_id(tDummyConnectivity( 0,0 ),EntityRank::ELEMENT);
                    tElemsConnectedToElem( 1,tCounter ) = get_loc_entity_index_from_entity_glb_id(tFaceId,EntityRank::FACE);
                    tCounter++;
                }
            }

            if ( tDummyConnectivity.n_cols()  > 1 )
            {
                if ( tDummyConnectivity(0, 1 ) != tId )
                {
                    tElemsConnectedToElem( 0,tCounter ) = get_loc_entity_index_from_entity_glb_id(tDummyConnectivity( 0,1 ),EntityRank::ELEMENT);
                    tElemsConnectedToElem( 1,tCounter ) = get_loc_entity_index_from_entity_glb_id(tFaceId,EntityRank::FACE);
                    tCounter++;
                }
            }

            XTK_ASSERT( tDummyConnectivity.n_cols()  <= 2,
                    "For some reason face has more than 2 elements connected to it... Check get_elements_connected_to_element." );
        }

        // Resize to include only ids added above and get rid of initialized extra zeros
        tElemsConnectedToElem.resize( 2,tCounter );

        return tElemsConnectedToElem;
    }

    moris::Mat_New<Integer, Integer_Matrix>
    get_entity_connected_to_entity_glb_ids( Integer aEntityIndex,
                                            enum EntityRank aInputEntityRank,
                                            enum EntityRank aOutputEntityRank) const
        {

        // Get Stk entity Id from local to global map
        Integer tId = this->get_glb_entity_id_from_entity_loc_index(aEntityIndex, aInputEntityRank);
        stk::mesh::EntityId tStkEntityId = { tId };

        // Call function that gets the connected entities
        stk::mesh::EntityRank tStkInputRank = this->get_stk_entity_rank(aInputEntityRank);
        stk::mesh::EntityRank tStkOutputRank = this->get_stk_entity_rank(aOutputEntityRank);
        stk::mesh::Entity tStkEntity = mStkMeshBulkData->get_entity(tStkInputRank, tStkEntityId);

        moris::Mat_New<Integer, Integer_Matrix> tEntitiesConnected = this->entities_connected_to_entity_stk_glb_id(&tStkEntity, tStkInputRank, tStkOutputRank);

        return tEntitiesConnected;
        }

    Integer
    get_element_face_ordinal_loc_inds( Integer const & aElementIndex,
                                       Integer const & aFaceIndex) const
    {
        moris::Mat_New<Integer, Integer_Matrix> tElementFaces = get_entity_connected_to_entity_loc_inds(aElementIndex,EntityRank::ELEMENT,EntityRank::FACE);

        Integer tOrdinal = 50;
        for(Integer iOrd = 0; iOrd<tElementFaces.n_cols(); iOrd++)
        {
            if(tElementFaces(0,iOrd) == aFaceIndex)
            {
                tOrdinal = iOrd;
                return tOrdinal;
            }
        }

        return tOrdinal;

    }


    Integer
    get_num_entities(enum EntityRank aEntityRank) const
    {
        // Initialize
        stk::mesh::EntityRank requestedRank = this->get_stk_entity_rank(aEntityRank);

        // Get all entities from meta data
        stk::mesh::Selector allEntities = mStkMeshMetaData->universal_part();

        Integer tSTKEntities = stk::mesh::count_selected_entities(allEntities, mStkMeshBulkData->buckets(requestedRank));
        Integer tExternalEntities = mMeshExternalEntityData.get_num_entities_external_data(aEntityRank);

        return tSTKEntities + tExternalEntities;
    }

    moris::Mat_New<Real,Real_Matrix>
    get_selected_node_coordinates_loc_inds( moris::Mat_New<Integer, Integer_Matrix> const & aNodeIndices) const
    {
        // TODO: Add external entity check to see if xtk has the coordinate field or stk has it
        // Number of spatial dimensions
        Integer tSpatialDimension = 3;

        // Get number of nodes provided
        Integer tNumNodes = aNodeIndices.n_cols();

        // Define part where nodes information is stored
        stk::mesh::Selector tUniversalSelector = mStkMeshMetaData->universal_part();

        // Get the coordinate field from stk
        stk::mesh::FieldBase const * coord = mStkMeshMetaData->coordinate_field();

        // Initialize output matrix
        Real tFillValue = 0;
        moris::Mat_New<Real,Real_Matrix> tSelectedNodesCoords(tNumNodes, tSpatialDimension, tFillValue);
        // Loop over all nodes
        enum EntityRank tEntityRank = EntityRank::NODE;

        for (Integer n = 0; n < tNumNodes; ++n)
        {
            if (mMeshExternalEntityData.is_external_entity(aNodeIndices(0, n), tEntityRank))
            {
                moris::Mat_New<Real,Real_Matrix> const & tNodeCoords = mMeshExternalEntityData.get_selected_node_coordinates_loc_inds_external_data(aNodeIndices(0, n));

                for (Integer dim = 0; dim < tSpatialDimension; dim++)
                {
                    tSelectedNodesCoords(n, dim) = tNodeCoords(0, dim);
                }
            }

            else
            {
                // Get node id from provided index
                stk::mesh::EntityId tNodeId = { this->get_glb_entity_id_from_entity_loc_index(aNodeIndices(0, n), EntityRank::NODE) };

                // Declare node entity
                stk::mesh::Entity tNodeEntity = mStkMeshBulkData->get_entity(stk::topology::NODE_RANK, tNodeId);

                // Get coordinates of node n
                double *fieldValue = static_cast<double *>(stk::mesh::field_data(*coord, tNodeEntity));

                // Store coodinates into output matrix
                for (Integer dim = 0; dim < tSpatialDimension; dim++)
                {
                    tSelectedNodesCoords(n, dim) = fieldValue[dim];
                }
            }
        }

        return tSelectedNodesCoords;
    }

    moris::Mat_New<Real,Real_Matrix>
    get_all_node_coordinates_loc_inds() const
    {
        // Get counts inside of STK and in External Data (both happen in this function call)
        Integer tNumNodes = this->get_num_entities(EntityRank::NODE);

        moris::Mat_New<Real,Real_Matrix> tAllNodeCoordinates(tNumNodes,3,0);

        // Get node coordinates from STK
        // Number of spatial dimensions
        Integer k = 0;
        Integer tLocalId = 0;
        Integer tSpatialDim = 3;

        // Define part where nodes information is stored
        stk::mesh::Selector tUniversalSelector = mStkMeshMetaData->universal_part();

        // Get the coordinates field from stk
        stk::mesh::FieldBase const * tCoordinateField = mStkMeshMetaData->coordinate_field();

        // Get the entities corresponding to all nodes selected
        std::vector<stk::mesh::Entity> nodes;
        stk::mesh::get_entities( *mStkMeshBulkData , stk::topology::NODE_RANK , nodes );

        // Loop over all nodes
        for (Integer n=0; n<nodes.size(); ++n)
        {
            // Get node and corresponding Id
            stk::mesh::Entity node = nodes[n];

            // Get coordinates of node n
            double* coords = static_cast<double *>( stk::mesh::field_data( * tCoordinateField , node ) );

            tLocalId = mStkMeshBulkData->local_id(node); // Index return by the local index assigned to entity

            // Store coodinates into output matrix
            for (Integer dim = 0; dim < tSpatialDim ; ++dim )
            {
                tAllNodeCoordinates( tLocalId , dim ) = coords[dim];
            }
            k++;
        }

        // Get node coordinates from external entities
        mMeshExternalEntityData.get_all_node_coordinates_loc_inds_external_data(k,tAllNodeCoordinates);

        return tAllNodeCoordinates;
    }

    moris::Mat_New<Integer, Integer_Matrix> get_all_entity_indices(enum EntityRank aEntityRank) const
        {
            Integer tNumEntity = this->get_num_entities(aEntityRank);
            moris::Mat_New<Integer, Integer_Matrix> tEntityIndices(1,tNumEntity);

            for(Integer i = 0; i<tNumEntity; i++)
            {
                tEntityIndices(0,i) = i;
            }

            return tEntityIndices;
        }

    Integer get_glb_entity_id_from_entity_loc_index(Integer aEntityIndex, enum EntityRank aEntityRank) const
    {
        Integer tGlbId = 0;
        if (mMeshExternalEntityData.is_external_entity(aEntityIndex, aEntityRank))
        {
            tGlbId = mMeshExternalEntityData.get_glb_entity_id_from_entity_loc_index_external_data(aEntityIndex, aEntityRank);
        }
        else
        {
            tGlbId = mMeshParallelData.get_entity_glb_id_from_entity_index(aEntityIndex, aEntityRank);
        }
        return tGlbId;
    }

    Integer get_loc_entity_index_from_entity_glb_id(Integer aEntityId, enum EntityRank aEntityRank) const
    {
        return mStkMeshBulkData->local_id( mStkMeshBulkData->get_entity( get_stk_entity_rank(aEntityRank), aEntityId ) );
    }

    Integer get_first_available_index(enum EntityRank aEntityRank) const
    {
        return mMeshExternalEntityData.get_first_available_index_external_data(aEntityRank);
    }

    moris::Mat_New<Integer, Integer_Matrix> const & get_local_to_global_map(enum EntityRank aEntityRank) const
        {
        return mMeshParallelData.get_local_to_global_map_parallel_data(aEntityRank);
        }

    moris::Mat_New<Real,Real_Matrix> get_entity_field_value(moris::Mat_New<Integer, Integer_Matrix> const & aEntityIndex, std::string const & aFieldName, enum EntityRank aFieldEntityRank) const
    {
        XTK_ASSERT(aFieldEntityRank==EntityRank::NODE,"Only implemented for nodal scalar field");

        // Initialize Output
        Integer tNumEntities= aEntityIndex.n_cols();
        moris::Mat_New<Real,Real_Matrix> tFieldValues(1,tNumEntities);

        // Get field by name and entity rank
        stk::mesh::EntityRank tEntityRank = this->get_stk_entity_rank(aFieldEntityRank);
        stk::mesh::Field<Real> * tField = mStkMeshIoBroker->meta_data().get_field<stk::mesh::Field<Real>>(tEntityRank,aFieldName);

        // Loop over entities and access field value
        for (Integer i = 0; i < tNumEntities; i++ )
        {
            // Get global Id of current node and create "node entity" for stk mesh
            //stk::mesh::EntityId nodeGlobalId = node_i;
            if (mMeshExternalEntityData.is_external_entity(aEntityIndex(0, i), aFieldEntityRank))
            {
                 tFieldValues(0,i) = mMeshExternalEntityData.get_entity_field_value_external_data(aEntityIndex(0, i),aFieldEntityRank,aFieldName);
            }

            else
            {
                Integer tId = get_glb_entity_id_from_entity_loc_index(aEntityIndex(0,i),aFieldEntityRank);
                stk::mesh::Entity tEntity = mStkMeshBulkData->get_entity(tEntityRank, tId);

                // Store the coordinates of the current node
                Real* tFieldData = stk::mesh::field_data ( *tField, tEntity );
                tFieldValues(0,i) = tFieldData[0];
            }
        }

        return tFieldValues;
    }

    Real get_entity_field_value(Integer const & aEntityIndex, std::string const & aFieldName, enum EntityRank aFieldEntityRank)
    {

        // Initialize Output
        Real tFieldValues = 0;

        // Get field by name and entity rank
        stk::mesh::EntityRank tEntityRank = this->get_stk_entity_rank(aFieldEntityRank);
        stk::mesh::Field<Real> * tField = mStkMeshIoBroker->meta_data().get_field<stk::mesh::Field<Real>>(tEntityRank,aFieldName);

        if (mMeshExternalEntityData.is_external_entity(aEntityIndex, aFieldEntityRank))
        {
            tFieldValues = mMeshExternalEntityData.get_entity_field_value_external_data(aEntityIndex,aFieldEntityRank,aFieldName);
        }

        else
        {
            Integer tId = get_glb_entity_id_from_entity_loc_index(aEntityIndex,aFieldEntityRank);
            stk::mesh::Entity tEntity = mStkMeshBulkData->get_entity(tEntityRank, tId);

            // Store the coordinates of the current node
            Real* tFieldData = stk::mesh::field_data ( *tField, tEntity );
            tFieldValues = tFieldData[0];
        }

        return tFieldValues;
    }

    Integer get_entity_field_value_integer(Integer const & aEntityIndex, std::string const & aFieldName, enum EntityRank aFieldEntityRank)
        {
            XTK_ASSERT(aFieldEntityRank==EntityRank::ELEMENT,"Only implemented for ELEMENT scalar field");

            // Initialize Output
            Integer tFieldValues = 0;

            // Get field by name and entity rank
            stk::mesh::EntityRank tEntityRank = this->get_stk_entity_rank(aFieldEntityRank);
            stk::mesh::Field<Integer> * tField = mStkMeshIoBroker->meta_data().get_field<stk::mesh::Field<Integer>>(tEntityRank,aFieldName);

            if (mMeshExternalEntityData.is_external_entity(aEntityIndex, aFieldEntityRank))
            {
                tFieldValues = mMeshExternalEntityData.get_entity_field_value_external_data(aEntityIndex,aFieldEntityRank,aFieldName);
            }

            else
            {
                Integer tId = get_glb_entity_id_from_entity_loc_index(aEntityIndex,aFieldEntityRank);
                stk::mesh::Entity tEntity = mStkMeshBulkData->get_entity(tEntityRank, tId);

                // Store the coordinates of the current node
                Integer* tFieldData = stk::mesh::field_data ( *tField, tEntity );
                tFieldValues = tFieldData[0];
            }

            return tFieldValues;
        }

    void update_first_available_index(Integer aNewFirstAvailableIndex, enum EntityRank aEntityRank)
    {
        mMeshExternalEntityData.update_first_available_index_external_data(aNewFirstAvailableIndex, aEntityRank);
    }

    void batch_create_new_nodes(xtk::Cell<xtk::Pending_Node<Real, Integer, Real_Matrix, Integer_Matrix>> const & aPendingNodes)
    {
        mMeshExternalEntityData.batch_create_new_nodes_external_data(aPendingNodes,false);
        mMeshParallelData.update_mesh_parallel_data(aPendingNodes);

    }

    void batch_create_new_nodes_with_fields(xtk::Cell<xtk::Pending_Node<Real, Integer, Real_Matrix, Integer_Matrix>> const & aPendingNodes,
                                            xtk::Cell<std::string> const & aFieldNames)
    {
        mMeshExternalEntityData.batch_create_new_nodes_with_fields_external_data(aPendingNodes,aFieldNames);
        mMeshParallelData.update_mesh_parallel_data(aPendingNodes);
    }

    Integer get_entity_parallel_owner_rank(Integer aEntityIndex, enum EntityRank aEntityRank) const
    {

        XTK_ASSERT(aEntityRank != EntityRank::NODE, "Node ownership not finished for external entities");

        // Convert index to ID
        Integer tEntityId = mMeshParallelData.get_entity_glb_id_from_entity_index(aEntityIndex, aEntityRank);

        //Get entity Id
        stk::mesh::Entity tEntity = mStkMeshBulkData->get_entity(get_stk_entity_rank(aEntityRank), tEntityId);

        //processor rank that owns entity
        Integer tOwningProcessor = mStkMeshBulkData->parallel_owner_rank(tEntity);

        return tOwningProcessor;
    }

    Integer get_num_of_entities_shared_with_processor(Integer aProcessorRank, enum EntityRank aEntityRank, bool aSendFlag) const
    {
        return mMeshParallelData.get_num_of_entities_shared_with_processor_parallel_data(aProcessorRank, aEntityRank, aSendFlag);
    }

    void get_processors_whom_share_entity(Integer aEntityIndex,
                                          enum EntityRank aEntityRank,
                                          moris::Mat_New<Integer, Integer_Matrix> & aProcsWhomShareEntity) const
    {
        // Convert index to ID
        stk::mesh::EntityId tEntityId = { mMeshParallelData.get_entity_glb_id_from_entity_index(aEntityIndex, aEntityRank) };

        //Get entity
        stk::mesh::Entity tEntity = mStkMeshBulkData->get_entity(get_stk_entity_rank(aEntityRank), tEntityId);

        // Intialize shared procs
        std::vector<int> tSharedProcs;

        // get shared processor IDs
        mStkMeshBulkData->comm_procs(mStkMeshBulkData->entity_key(tEntity), tSharedProcs);

        if (tSharedProcs.size() == 0)
        {
            tSharedProcs.push_back(mStkMeshBulkData->parallel_owner_rank(tEntity));
        }

        // Initialize output
        aProcsWhomShareEntity.resize(1, tSharedProcs.size());

        // Cell to vector conversion
        for (Integer i = 0; i < tSharedProcs.size(); i++)
        {
            aProcsWhomShareEntity(0, i) = tSharedProcs[i];
        }
    }

    Integer allocate_entity_ids(Integer aNumIdstoAllocate, enum EntityRank aEntityRank) const
    {
        return mMeshExternalEntityData.allocate_entity_ids_external_entity_data(aNumIdstoAllocate, aEntityRank);
    }

    void add_mesh_field_data_loc_indices(std::string const &     aFieldName,
                                         enum EntityRank const & aFieldEntityRank,
                                         xtk::Cell<Real> const & aFieldData)
    {
        // Write Data to Field
        Integer tNumEntities = get_num_entities(aFieldEntityRank);

        // Get Field
        stk::mesh::EntityRank tEntityRank = this->get_stk_entity_rank(aFieldEntityRank);
        stk::mesh::Field<Real> * tField = mStkMeshMetaData->get_field<stk::mesh::Field<Real>>(tEntityRank,aFieldName);
        for (Integer i = 0; i < tNumEntities; i++ )
        {
            // Get global Id of current node and create "node entity" for stk mesh
            //stk::mesh::EntityId nodeGlobalId = node_i;
            Integer tId = get_glb_entity_id_from_entity_loc_index(i,aFieldEntityRank);
            stk::mesh::Entity tEntity = mStkMeshBulkData->get_entity(tEntityRank, tId);

            // Store the coordinates of the current node
            Real* tFieldData = stk::mesh::field_data ( *tField, tEntity );
            tFieldData[0] = aFieldData(i);
        }
    }


    void add_mesh_field_data_loc_indices_integer(std::string const &     aFieldName,
                                                 enum EntityRank const & aFieldEntityRank,
                                                 xtk::Cell<Integer> const &    aFieldData)
    {
        // Write Data to Field
        Integer tNumEntities = get_num_entities(aFieldEntityRank);

        // Get Field
        stk::mesh::EntityRank tEntityRank = this->get_stk_entity_rank(aFieldEntityRank);
        stk::mesh::Field<Integer> * tField = mStkMeshMetaData->get_field<stk::mesh::Field<Integer>>(tEntityRank,aFieldName);
        for (Integer i = 0; i < tNumEntities; i++ )
        {
            // Get global Id of current node and create "node entity" for stk mesh
            //stk::mesh::EntityId nodeGlobalId = node_i;
            Integer tId = get_glb_entity_id_from_entity_loc_index(i,aFieldEntityRank);
            stk::mesh::Entity tEntity = mStkMeshBulkData->get_entity(tEntityRank, tId);

            // Store the coordinates of the current node
            Integer* tFieldData = stk::mesh::field_data ( *tField, tEntity );
            tFieldData[0] = aFieldData(i);
        }
    }



    stk::mesh::Field<Real, stk::mesh::Cartesian3d> * get_coordinate_field()
    {
        return mStkMeshMetaData->get_field<stk::mesh::Field<Real, stk::mesh::Cartesian3d> >(stk::topology::NODE_RANK, "coordinates");
    }

    stk::mesh::FieldBase * get_field(enum EntityRank aFieldRank,
                                     std::string aFieldName)
    {
        return mStkMeshMetaData->get_field<stk::mesh::Field<Real> >(get_stk_entity_rank(aFieldRank), aFieldName);
    }

    void write_output_mesh(std::string            const & aMeshName,
                           xtk::Cell<std::string> const & aRealNodeFieldsToOutput       = {},
                           xtk::Cell<std::string> const & aIntNodeFieldsToOutput        = {},
                           xtk::Cell<std::string> const & aRealElementFieldsToOutput    = {},
                           xtk::Cell<std::string> const & aIntElementFieldsToOutput     = {},
                           xtk::Cell<std::string> const & aRealVectorNodeFieldsToOutput = {},
                           Integer aTime = 0) const
    {

        std::string tOutputMesh = aMeshName;
        std::cout<<"Output Path: "<< tOutputMesh<<std::endl;

        // If we loaded this mesh from a file then the io broker has been set up otherwise we set up the io broker now
        if(!mLoadedFromFile)
        {
            mStkMeshIoBroker->set_bulk_data(*mStkMeshBulkData);
        }


        Integer tFileHandle = mStkMeshIoBroker->create_output_mesh(tOutputMesh, stk::io::WRITE_RESULTS);
        mStkMeshIoBroker->write_output_mesh(tFileHandle);

        // Add node real fields to the output
        for(Integer i = 0; i <aRealNodeFieldsToOutput.size(); i++)
        {
            stk::mesh::Field<Real> * tField = mStkMeshIoBroker->meta_data().get_field<stk::mesh::Field<Real>>(stk::topology::NODE_RANK,aRealNodeFieldsToOutput(i));
            mStkMeshIoBroker->add_field(tFileHandle, *tField);
        }


        // Add node integer fields to the output
        for(Integer i = 0; i <aIntNodeFieldsToOutput.size(); i++)
        {
            stk::mesh::Field<Integer> * tField = mStkMeshIoBroker->meta_data().get_field<stk::mesh::Field<Integer>>(stk::topology::NODE_RANK,aIntNodeFieldsToOutput(i));
            mStkMeshIoBroker->add_field(tFileHandle, *tField);
        }


        // Add element real fields to the output
        for(Integer i = 0; i <aRealElementFieldsToOutput.size(); i++)
        {
            stk::mesh::Field<Real> * tField = mStkMeshIoBroker->meta_data().get_field<stk::mesh::Field<Real>>(stk::topology::ELEMENT_RANK,aRealElementFieldsToOutput(i));
            mStkMeshIoBroker->add_field(tFileHandle, *tField);
        }


        // Add element integer fields to the output
        for(Integer i = 0; i <aIntElementFieldsToOutput.size(); i++)
        {
            stk::mesh::Field<Integer> * tField = mStkMeshIoBroker->meta_data().get_field<stk::mesh::Field<Integer>>(stk::topology::ELEMENT_RANK,aIntElementFieldsToOutput(i));
            mStkMeshIoBroker->add_field(tFileHandle, *tField);
        }

        // Add node real cartesian vector fields to the output
        for(Integer i = 0; i <aRealVectorNodeFieldsToOutput.size(); i++)
        {
            stk::mesh::Field<Real,stk::mesh::Cartesian3d> * tField = mStkMeshIoBroker->meta_data().get_field<stk::mesh::Field<Real,stk::mesh::Cartesian3d>>(stk::topology::NODE_RANK,aRealVectorNodeFieldsToOutput(i));
            mStkMeshIoBroker->add_field(tFileHandle, *tField);
        }

        mStkMeshIoBroker->begin_output_step(tFileHandle, aTime);
        mStkMeshIoBroker->write_defined_output_fields(tFileHandle);
        mStkMeshIoBroker->end_output_step(tFileHandle);
    }



    Integer get_num_buckets(enum EntityRank aEntityRank) const
    {
        stk::mesh::EntityRank tEntityRank = this->get_stk_entity_rank(aEntityRank);
        stk::mesh::Selector tSelector(mStkMeshMetaData->locally_owned_part());
        stk::mesh::BucketVector const & tBuckets = mStkMeshBulkData->get_buckets(tEntityRank,tSelector);
        return tBuckets.size();
    }

    moris::Mat_New<Integer, Integer_Matrix> get_entities_in_bucket_loc_index(Integer aBucketOrdinal, enum EntityRank aEntityRank) const
    {
        /*
         * Get Bucket Vector
         */
        stk::mesh::EntityRank tEntityRank = this->get_stk_entity_rank(aEntityRank);
        stk::mesh::Selector tSelector(mStkMeshMetaData->locally_owned_part());
        stk::mesh::BucketVector const & tBuckets = mStkMeshBulkData->get_buckets(tEntityRank,tSelector);

        Integer tNumEntitiesInBucket = tBuckets[aBucketOrdinal]->size();
        moris::Mat_New<Integer, Integer_Matrix> tBucketEntityIndices(1,tNumEntitiesInBucket);

        /*
         * Iterate through specified bucket
         */

        Integer j = 0;
        for(stk::mesh::Bucket::iterator iEnt = tBuckets[aBucketOrdinal]->begin(); iEnt!=tBuckets[aBucketOrdinal]->end(); iEnt++)
        {
           tBucketEntityIndices(0,j) =  mStkMeshBulkData->local_id(*iEnt);
           j++;
        }

        return tBucketEntityIndices;
    }

    void get_entity_part_membership_ordinals(Integer const & aEntityIndex, enum EntityRank aEntityRank, xtk::Cell<Integer> & aPartOrdinals, bool const & aInduced = false) const
    {
        /*
         * Get global Id of current node and create "node entity" for stk mesh
         * stk::mesh::EntityId nodeGlobalId = node_i;
         */
        Integer tId = get_glb_entity_id_from_entity_loc_index(aEntityIndex,aEntityRank);
        stk::mesh::EntityRank tEntityRank = this->get_stk_entity_rank(aEntityRank);
        stk::mesh::Entity tEntity = mStkMeshBulkData->get_entity(tEntityRank, tId);

        /*
         * Get Entities bucket
         */
        stk::mesh::Bucket & tBucket = mStkMeshBulkData->bucket(tEntity);

        /*
         * Get part ordinals from the bucket
         * Note: this include internal parts
         */
        stk::mesh::OrdinalVector tPartOrdinals;
        tBucket.supersets(tPartOrdinals);

        /*
         * Copy to xtk::Cell
         */

        aPartOrdinals.clear();
        Integer tNumParts = tPartOrdinals.size();
        aPartOrdinals.reserve(tNumParts);

        if(aEntityRank !=EntityRank::FACE)
        {
            for(Integer i = 0; i<tNumParts; i++)
            {
                stk::mesh::Part & tPart = mStkMeshMetaData->get_part(tPartOrdinals[i]);
                if(!stk::mesh::impl::is_internal_part(tPart) && tEntityRank == tPart.primary_entity_rank())
                {

                    aPartOrdinals.push_back(tPartOrdinals[i]);
                }
            }
        }

        else
        {
            for(Integer i = 0; i<tNumParts -1 ; i++)
            {
                stk::mesh::Part & tPart = mStkMeshMetaData->get_part(tPartOrdinals[i]);

                if(!stk::mesh::impl::is_internal_part(tPart) && tEntityRank == tPart.primary_entity_rank())
                {

                    aPartOrdinals.push_back(tPartOrdinals[i]);
                }
            }
        }

    }


    void get_all_part_names(enum EntityRank const & aPrimaryRank,
                            xtk::Cell<std::string> & aPartNames) const
    {
        aPartNames.clear();
        stk::mesh::PartVector const tParts = mStkMeshMetaData->get_mesh_parts();
        aPartNames.reserve(tParts.size());

        stk::mesh::EntityRank tPrimaryRank = this->get_stk_entity_rank(aPrimaryRank);

        for(Integer iParts = 0; iParts<tParts.size(); iParts++)
        {
            if(tPrimaryRank == tParts[iParts]->primary_entity_rank())
            {
                aPartNames.push_back(tParts[iParts]->name());
            }
        }

    }

    void get_part_name_from_part_ordinals(xtk::Cell<Integer> const & aPartOrdinals, xtk::Cell<std::string> & aPartNames) const
    {

        aPartNames.clear();
        aPartNames.reserve(aPartOrdinals.size());

        for(Integer iParts = 0; iParts<aPartOrdinals.size(); iParts++)
        {
            stk::mesh::Part & tPart = mStkMeshMetaData->get_part(aPartOrdinals(iParts));
            aPartNames.push_back(tPart.name());
        }
    }


    /*
     * Returns the number of basis functions in the mesh
     */
    Integer
    get_num_basis_functions() const
    {
        // Initialize
        stk::mesh::EntityRank requestedRank = this->get_stk_entity_rank(EntityRank::NODE);

        // Get all entities from meta data
        stk::mesh::Selector allEntities = mStkMeshMetaData->universal_part();

        return stk::mesh::count_selected_entities(allEntities, mStkMeshBulkData->buckets(requestedRank));
    }


    /*
     *  Returns the elements in the support of a the basis function with index (aBasisIndex)
     */
    moris::Mat_New<Integer, Integer_Matrix>
    get_elements_in_basis_support(Integer aBasisIndex) const
    {
        return this->get_entity_connected_to_entity_loc_inds(aBasisIndex, EntityRank::NODE,EntityRank::ELEMENT);
    }


    // Used only in the mesh builder to access the member variables
    // These are not accessible via a ptr to the base class
    stk::io::StkMeshIoBroker & mesh_io_broker()
    {
        return (*mStkMeshIoBroker);
    }

    stk::mesh::MetaData & mesh_meta_data()
    {
        return (*mStkMeshMetaData);
    }

    stk::mesh::BulkData & mesh_bulk_data()
    {
        return (*mStkMeshBulkData);
    }

    void setup_parallel_mesh_data()
    {
        mMeshParallelData.setup_mesh_parallel_data(*mStkMeshBulkData, *mStkMeshMetaData);
    }

    void setup_mesh_external_entities_and_first_available_information()
    {
        Integer tNumNodes = get_num_entities(EntityRank::NODE);
        Integer tNumEdges = get_num_entities(EntityRank::EDGE);
        Integer tNumFaces = get_num_entities(EntityRank::FACE);
        Integer tNumElements = get_num_entities(EntityRank::ELEMENT);

        mMeshExternalEntityData.set_up_external_entity_data(tNumNodes, tNumEdges, tNumFaces, tNumElements, *mStkMeshBulkData);
    }

private:
    bool mLoadedFromFile;
    std::shared_ptr<stk::mesh::MetaData> mStkMeshMetaData;
    std::shared_ptr<stk::mesh::BulkData> mStkMeshBulkData;
    std::shared_ptr<stk::io::StkMeshIoBroker> mStkMeshIoBroker;

    // Xtk member variables
    mesh::Mesh_Parallel_Data_Stk<Real, Integer, Real_Matrix, Integer_Matrix> mMeshParallelData;
    mesh::Mesh_External_Entity_Data_Stk<Real, Integer, Real_Matrix, Integer_Matrix> mMeshExternalEntityData;

    // Private functions for implementation
private:
    // Functions needed to supply above functionality
    moris::Mat_New<Integer, Integer_Matrix>
    entities_connected_to_entity_stk_glb_id(stk::mesh::Entity* const aInputEntity,
                                            stk::mesh::EntityRank const aInputEntityRank,
                                            stk::mesh::EntityRank const aOutputEntityRank) const
    {
        // Declare the object where we are going to store the shared faces and handlers
        std::vector<stk::mesh::Entity> tDesiredEntitiesConnectedToInputEntities;

        // Check if the connectivity exists (i.e., was already generated and is stored in mesh data)
        if (mStkMeshBulkData->connectivity_map().valid(aInputEntityRank, aOutputEntityRank))
        {

            switch (aOutputEntityRank)
            {
            case stk::topology::NODE_RANK:

                // Fill entities connected
                if (mStkMeshBulkData->num_nodes(aInputEntity[0]) > 0)
                {
                    // Get pointers to the location of the connected nodes
                    stk::mesh::Entity const * tDesiredEntityStart = mStkMeshBulkData->begin_nodes(aInputEntity[0]);
                    stk::mesh::Entity const * tDesiredEntityEnd = mStkMeshBulkData->end_nodes(aInputEntity[0]);

                    // Store faces in output vector
                    tDesiredEntitiesConnectedToInputEntities.assign(tDesiredEntityStart, tDesiredEntityEnd);
                }
                break;

            case stk::topology::EDGE_RANK:

                // Fill entities connected
                if (mStkMeshBulkData->num_edges(aInputEntity[0]) > 0)
                {
                    // Get pointers to the location of the connected edges
                    stk::mesh::Entity const * tDesiredEntityStart = mStkMeshBulkData->begin_edges(aInputEntity[0]);
                    stk::mesh::Entity const * tDesiredEntityEnd = mStkMeshBulkData->end_edges(aInputEntity[0]);

                    // Store faces in output vector
                    tDesiredEntitiesConnectedToInputEntities.assign(tDesiredEntityStart, tDesiredEntityEnd);
                }
                break;

            case stk::topology::FACE_RANK:

                // Fill entities connected
                if (mStkMeshBulkData->num_faces(aInputEntity[0]) > 0)
                {
                    // Get pointers to the location of the connected faces
                    stk::mesh::Entity const * tDesiredEntityStart = mStkMeshBulkData->begin_faces(aInputEntity[0]);
                    stk::mesh::Entity const * tDesiredEntityEnd = mStkMeshBulkData->end_faces(aInputEntity[0]);

                    // Store faces in output vector
                    tDesiredEntitiesConnectedToInputEntities.assign(tDesiredEntityStart, tDesiredEntityEnd);
                }
                break;

            case stk::topology::ELEMENT_RANK:

                // Fill entities connected
                if (mStkMeshBulkData->num_elements(aInputEntity[0]) > 0)
                {
                    // Get pointers to the location of the connected elements
                    stk::mesh::Entity const * tDesiredEntityStart = mStkMeshBulkData->begin_elements(aInputEntity[0]);
                    stk::mesh::Entity const * tDesiredEntityEnd = mStkMeshBulkData->end_elements(aInputEntity[0]);

                    // Store faces in output vector
                    tDesiredEntitiesConnectedToInputEntities.assign(tDesiredEntityStart, tDesiredEntityEnd);
                }
                break;

            default:
                std::cerr << " wrong topology in entities_connected_to_entity_stk_glb_id ";
                break;
            }
        }
        else
        {
            std::cerr << " STK already has valid connectivity maps. Check if you are trying to access invalid connectivity (e.g., edge to edge)";
        }

        // Get number of connected entities
        Integer tNumEntities = tDesiredEntitiesConnectedToInputEntities.size();

        Integer tFillValue = 0;
        moris::Mat_New<Integer, Integer_Matrix> tEntitiesConnectedToGivenEntity(1, tNumEntities, tFillValue);

        // Loop over connected entities
        for (Integer i = 0; i < tNumEntities; ++i)
        {
            // Store the entity ID in output cell
            tEntitiesConnectedToGivenEntity(0, i) = (Integer) mStkMeshBulkData->identifier(tDesiredEntitiesConnectedToInputEntities[i]);
        }

        return tEntitiesConnectedToGivenEntity;

    }

    stk::mesh::EntityRank get_stk_entity_rank(enum EntityRank aXtkEntityRank) const
    {
        if (aXtkEntityRank == EntityRank::NODE)
        {
            return stk::topology::NODE_RANK;
        }
        else
            if (aXtkEntityRank == EntityRank::EDGE)
            {
                return stk::topology::EDGE_RANK;
            }

            else
                if (aXtkEntityRank == EntityRank::FACE)
                {
                    return stk::topology::FACE_RANK;
                }
                else
                    if (aXtkEntityRank == EntityRank::ELEMENT)
                    {
                        return stk::topology::ELEMENT_RANK;
                    }
                    else
                    {
                        return stk::topology::INVALID_RANK;
                    }
    }

};
}

#endif /* SRC_MESH_CL_MESH_DATA_STK_HPP_ */
