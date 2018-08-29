/*
 * cl_Mesh_Builder_Stk.hpp
 *
 *  Created on: Jun 20, 2017
 *      Author: ktdoble
 */

#ifndef SRC_MESH_CL_MESH_BUILDER_STK_HPP_
#define SRC_MESH_CL_MESH_BUILDER_STK_HPP_

#include "mesh/cl_Mesh_Builder.hpp"
#include "mesh/cl_Mesh_Entity.hpp"
#include "mesh/cl_Mesh_Enums.hpp"
#include "mesh/cl_Mesh_Data.hpp"
#include "mesh/cl_Mesh_Data_Stk.hpp"
#include "linalg/cl_XTK_Matrix.hpp"
#include "containers/cl_XTK_Cell.hpp"
#include "xtk/cl_XTK_Enums.hpp"


// TPL: STK Includes
#include <stk_io/StkMeshIoBroker.hpp>    // for StkMeshIoBroker
#include <stk_mesh/base/GetEntities.hpp> // for count_entities
#include <stk_mesh/base/BulkData.hpp>    // for BulkData
#include <stk_mesh/base/MetaData.hpp>    // for MetaData
#include <stk_mesh/base/CreateFaces.hpp> // for handling faces
#include <stk_mesh/base/CreateEdges.hpp> // for handling faces
#include <stk_mesh/base/FEMHelpers.hpp>
#include "stk_mesh/base/Field.hpp"       // for Field
#include "stk_mesh/base/FieldBase.hpp"   // for field_data, etc
#include "linalg/cl_XTK_Matrix.hpp"
#include "Ioss_Region.h"                 // for Region, NodeBlockContainer
#include "Ioss_Field.h"                  // for Field

// XTKL: Logging and Assertion Functions
#include "ios/cl_Logger.hpp"
#include "assert/fn_xtk_assert.hpp"
#include "tools/cl_MPI_Tools.hpp"

using namespace xtk;

namespace mesh
{

void calculate_shared_node_data(stk::mesh::BulkData &aBulkData,
                                stk::mesh::MetaData &aMetaData)
{
    int tNumProcs;
    int tMyRank;
    MPI_Comm_size(get_comm(), &tNumProcs);
    MPI_Comm_rank(get_comm(), &tMyRank);
    if(tNumProcs > 1)
    {
        std::vector<stk::mesh::Entity> aMyOwnedNodeEntities;
        aBulkData.get_entities(stk::topology::NODE_RANK, aMetaData.locally_owned_part(), aMyOwnedNodeEntities);
        int tMyNumNodes = aMyOwnedNodeEntities.size();
        std::vector<int> tMyGlobalNodeIds(tMyNumNodes, 0);
        for(int i = 0; i < tMyNumNodes; ++i)
        {
            tMyGlobalNodeIds[i] = aBulkData.identifier(aMyOwnedNodeEntities[i]);
        }
        std::vector<int> tLocalNumNodesOnProcs(tNumProcs, 0);
        std::vector<int> tGlobalNumNodesOnProcs(tNumProcs, 0);
        tLocalNumNodesOnProcs[tMyRank] = tMyNumNodes;
        MPI_Allreduce(tLocalNumNodesOnProcs.data(), tGlobalNumNodesOnProcs.data(), tNumProcs, MPI_INT, MPI_MAX, get_comm());
        int tMaxNumOnProc = 0;
        for(int i=0; i<tNumProcs; ++i)
        {
            if(tMaxNumOnProc < tGlobalNumNodesOnProcs[i])
                tMaxNumOnProc = tGlobalNumNodesOnProcs[i];
        }
        std::vector<int> tAllGlobalIds(tMaxNumOnProc*tNumProcs, 0);
        std::vector<int> tResizedGlobalIdsForThisProc(tMaxNumOnProc, 0);
        std::copy(tMyGlobalNodeIds.begin(), tMyGlobalNodeIds.end(), tResizedGlobalIdsForThisProc.begin());
        MPI_Allgather(tResizedGlobalIdsForThisProc.data(), tMaxNumOnProc, MPI_INT, tAllGlobalIds.data(), tMaxNumOnProc, MPI_INT, get_comm());
        // At this point every processor should have all of the global ids and info about which proc they are on.
        // See which nodes are on other procs and add sharing info for them.
        for(int i = 0; i < tMyNumNodes; ++i)
        {
            int tCurGlobalNodeId = aBulkData.identifier(aMyOwnedNodeEntities[i]);
            int tStartIndex = 0;
            for(int j=0; j<tNumProcs; ++j)
            {
                if(j != tMyRank)
                {
                    for(int k=0; k<tGlobalNumNodesOnProcs[j]; ++k)
                    {
                        if(tCurGlobalNodeId == tAllGlobalIds[tStartIndex+k])
                        {

                            aBulkData.add_node_sharing(aMyOwnedNodeEntities[i], j);
                            k=tGlobalNumNodesOnProcs[j];
                        }
                    }
                }
                tStartIndex += tMaxNumOnProc;
            }
        }
    }
}




template<typename Real, typename Integer, typename Real_Matrix, typename Integer_Matrix>
class Mesh_Builder_Stk : public Mesh_Builder<Real, Integer, Real_Matrix, Integer_Matrix>
{
public:
    // Build mesh from file
    Mesh_Builder_Stk()
    {
    }

    std::shared_ptr<mesh::Mesh_Data<Real, Integer,Real_Matrix, Integer_Matrix>> build_mesh_from_string(std::string const & aMeshFileName,
                                                                                                       Cell<std::string> const & aScalarFieldNames,
                                                                                                       bool aCreateFacesAndEdges)
    {

//        std::string tMeshFileName = prefix+aMeshFileName;
        std::string tMeshFileName = aMeshFileName;



        // Declare MPI communicator and get parallel size
        int tParallelPoolSize = 0;
        MPI_Comm tCommunicator = get_comm();
        MPI_Comm_size(tCommunicator, &tParallelPoolSize);

        // Specify that STK should use an automatic aura
        stk::mesh::BulkData::AutomaticAuraOption tAutoAuraOption = stk::mesh::BulkData::AutomaticAuraOption::NO_AUTO_AURA;

        // Initialize meta, bulk, mesh io broker via Mesh_Data_Stk
        std::shared_ptr<Mesh_Data_Stk<Real, Integer,Real_Matrix, Integer_Matrix>> tSTKMeshData =
                std::make_shared<Mesh_Data_Stk<Real, Integer, Real_Matrix, Integer_Matrix>>(tCommunicator, tParallelPoolSize, tAutoAuraOption,true);

        // Create mesh database using the IO broker
        // Create mesh database using the IO broker
        stk::io::StkMeshIoBroker & tMeshIo = tSTKMeshData->mesh_io_broker();
        tMeshIo.set_bulk_data(tSTKMeshData->mesh_bulk_data());
        tMeshIo.add_mesh_database(tMeshFileName, stk::io::READ_MESH);
        tMeshIo.create_input_mesh();


        tMeshIo.add_all_mesh_fields_as_input_fields();

        // Put fields on entire mesh
        for(Integer i = 0; i<aScalarFieldNames.size(); i++)
        {
            stk::mesh::Field<Real> & tNodeField = tMeshIo.meta_data().declare_field<stk::mesh::Field<Real> >(stk::topology::NODE_RANK, aScalarFieldNames(i), 1);

            stk::mesh::put_field_on_entire_mesh(tNodeField);

            stk::io::set_field_role(tNodeField, Ioss::Field::TRANSIENT);
        }

        // Include mesh fields and populate the database
        tSTKMeshData->mesh_io_broker().populate_bulk_data();

        // Create mesh edge entities
        stk::mesh::create_edges(tSTKMeshData->mesh_bulk_data());

        // Create mesh face entities
        stk::mesh::create_faces(tSTKMeshData->mesh_bulk_data(), true);

        // Tell the IO broker to read data if there is data
        Teuchos::RCP<Ioss::Region> tRegion = tMeshIo.get_input_io_region();
        int step_count = tRegion->get_property("state_count").get_int();

        if(step_count!=0)
        {
            tMeshIo.read_defined_input_fields(1);
        }

        // Setup the external entity implementation
        tSTKMeshData->setup_mesh_external_entities_and_first_available_information();

        // Tell the Mesh data to setup the information it needs
        tSTKMeshData->setup_parallel_mesh_data();

        // Dynamic cast to return a pointer to the base class
        std::shared_ptr<mesh::Mesh_Data<Real, Integer, Real_Matrix, Integer_Matrix>> tMeshData = std::dynamic_pointer_cast<
                mesh::Mesh_Data<Real, Integer, Real_Matrix, Integer_Matrix>>(tSTKMeshData);

        return tSTKMeshData;
    }




    std::shared_ptr<mesh::Mesh_Data<Real, Integer,Real_Matrix, Integer_Matrix>>
    build_mesh_from_data(Integer const                                                            aSpatialDimension,
                         Cell<Bucket<Integer,Integer_Matrix>> const &                             aElementBuckets,
                         Cell<Side_Set_Input<Integer, Integer_Matrix>> const &                    aSideSets,
                         Cell<Node_Set_Input<Real, Integer, Real_Matrix, Integer_Matrix>> const & aNodeSets,
                         moris::Mat_New<Real,Real_Matrix> const &                                            aNodeCoordinates,
                         moris::Mat_New<Integer, Integer_Matrix> const &                                     aLocaltoGlobalNodeMap,
                         Cell<std::string> const &                                                aElementPartNames,
                         Cell<enum EntityTopology> &                                              aElementPartTopologys,
                         Cell<std::string> const &                                                aSideSetNames,
                         Cell<std::string> const &                                                aNodeSetNames,
                         Cell<std::string> const &                                                aInterfaceNodeSetNames,
                         Cell<std::string> const &                                                aInterfaceSideSetNames,
                         Cell<std::string> const &                                                aRealScalarFields,
                         Cell<std::string> const &                                                aRealVectorFields,
                         Cell<std::string> const &                                                aElementRealFields,
                         Sensitivity<Real, Integer, Real_Matrix, Integer_Matrix> &                aSensitivityData,
                         bool const &                                                             aSetupDataForInternalUse = false) const
        {


        // State whether a local to global map is needed
            Integer tMapFlag = 1;
            if(aLocaltoGlobalNodeMap(0,0) == std::numeric_limits<Integer>::max())
            {
                tMapFlag = 0;
            }
            // Declare MPI communicator and get parallel size
            int tParallelPoolSize = 0;
            MPI_Comm tCommunicator = get_comm();
            MPI_Comm_size(tCommunicator, &tParallelPoolSize);

            // Specify that STK should use an automatic aura
            stk::mesh::BulkData::AutomaticAuraOption tAutoAuraOption = stk::mesh::BulkData::AutomaticAuraOption::NO_AUTO_AURA;

            // Declare and initialize Stk mesh
            std::shared_ptr<Mesh_Data_Stk<Real, Integer,Real_Matrix, Integer_Matrix>> tSTKMeshData = std::make_shared<Mesh_Data_Stk<Real, Integer, Real_Matrix, Integer_Matrix>>(tCommunicator, tParallelPoolSize, tAutoAuraOption,false);

            stk::mesh::MetaData & tMeshMeta = tSTKMeshData->mesh_meta_data();

            /*
             * Declare element parts
             * Add other parts later
             */

            Integer tNumParts = aElementPartNames.size();
            stk::topology::topology_t tSTKTopology = stk::topology::INVALID_TOPOLOGY;
            for(Integer iPart = 0;  iPart<aElementPartTopologys.size(); iPart++)
            {
                this->get_stk_topology(aElementPartTopologys(iPart),tSTKTopology);
                stk::mesh::Part & tPart = tMeshMeta.declare_part_with_topology(aElementPartNames(iPart),tSTKTopology);
                stk::io::put_io_part_attribute(tPart);
            }

            tNumParts = aSideSetNames.size();
            for(Integer iPart = 0;  iPart<tNumParts; iPart++)
            {
                stk::mesh::Part & tPart = tMeshMeta.declare_part(aSideSetNames(iPart),tMeshMeta.side_rank());
                stk::io::put_io_part_attribute(tPart);
            }


            Integer tNumNodeSets = aNodeSetNames.size();
            for(Integer iPart = 0;  iPart<tNumNodeSets; iPart++)
            {
                stk::mesh::Part & tPart = tMeshMeta.declare_part(aNodeSetNames(iPart),stk::topology::NODE_RANK);
                stk::io::put_io_part_attribute(tPart);
            }

            Integer tNumInterfaceNodeSets = aInterfaceNodeSetNames.size();
            for(Integer iPart = 0;  iPart<tNumInterfaceNodeSets; iPart++)
            {
                stk::mesh::Part & tPart = tMeshMeta.declare_part(aInterfaceNodeSetNames(iPart),stk::topology::NODE_RANK);
                stk::io::put_io_part_attribute(tPart);
            }

            /*
             * Field stuff
             */
            Real tVal = 0.0 ; // any value for initialization
            stk::mesh::Field<Real,stk::mesh::Cartesian3d>& coord_field = tMeshMeta.declare_field<stk::mesh::Field<double,stk::mesh::Cartesian3d> >(stk::topology::NODE_RANK, "coordinates");
            stk::mesh::put_field(coord_field, tMeshMeta.universal_part(), 3, &tVal);

            // Scalar  field declarations
            for(Integer iSField = 0; iSField<aRealScalarFields.size(); iSField++)
            {
                stk::mesh::Field<Real>& tNodeField = tMeshMeta.declare_field<stk::mesh::Field<Real> >(stk::topology::NODE_RANK, aRealScalarFields(iSField));
                stk::mesh::put_field(tNodeField, tMeshMeta.universal_part() , 1, &tVal);
            }

            for(Integer iVField = 0; iVField<aRealVectorFields.size(); iVField++)
            {
                stk::mesh::Field<Real,stk::mesh::Cartesian3d>& tNodeVectorField = tMeshMeta.declare_field<stk::mesh::Field<Real,stk::mesh::Cartesian3d> >(stk::topology::NODE_RANK, aRealVectorFields(iVField));
                stk::mesh::put_field(tNodeVectorField, tMeshMeta.universal_part() , 3, &tVal);
            }

            for(Integer iESField = 0; iESField<aElementRealFields.size(); iESField++)
            {
                stk::mesh::Field<Real>& tElementField = tMeshMeta.declare_field<stk::mesh::Field<Real> >(stk::topology::ELEMENT_RANK, aElementRealFields(iESField));
                stk::mesh::put_field(tElementField, tMeshMeta.universal_part() , &tVal);
            }


            Cell<stk::mesh::Field<Real,stk::mesh::Cartesian3d>*> tInterfaceFields;
            Cell<stk::mesh::Field<Integer>*> tADVIndexField;
//            stk::mesh::Field<Integer>* tNADVIndexField;
            // output sensitivity densely
            if( !aSensitivityData.output_sparesly()){


                for(Integer iG = 0; iG<aInterfaceNodeSetNames.size(); iG++)
                {
                    stk::mesh::Part* tInterfaceNodes = tMeshMeta.get_part(aInterfaceNodeSetNames(iG));
                    Cell<std::string> const & tADVFieldNames = aSensitivityData.get_all_field_names();
                    tInterfaceFields = Cell<stk::mesh::Field<Real,stk::mesh::Cartesian3d>*>(tADVFieldNames.size());

                    for(Integer i = 0; i<tADVFieldNames.size(); i++)
                    {
                        stk::mesh::Field<Real,stk::mesh::Cartesian3d>& tADVField = tMeshMeta.declare_field<stk::mesh::Field<Real,stk::mesh::Cartesian3d> >(stk::topology::NODE_RANK, tADVFieldNames(i));
                        stk::mesh::put_field(tADVField, *tInterfaceNodes , 3, &tVal);

                        tInterfaceFields(i) = &tADVField;
                    }
                }
            }

            // Sparse packing of sensitivity
            else
            {
                //TODO: MAX NUM OF ADVS ON A PROC RESULTS IN AN INCONSISTENT NUMBER OF PROCS DECLARING FIELDS.
                Integer tMaxNumADVs = aSensitivityData.get_max_num_advs();

                for(Integer iG = 0; iG<aInterfaceNodeSetNames.size(); iG++)
                {
                    stk::mesh::Part* tInterfaceNodes = tMeshMeta.get_part(aInterfaceNodeSetNames(iG));
                    Cell<std::string> const & tADVFieldNames    = aSensitivityData.get_all_field_names();
                    Cell<std::string> const & tADVIndFieldNames = aSensitivityData.get_adv_ind_field_name();
                    tInterfaceFields = Cell<stk::mesh::Field<Real,stk::mesh::Cartesian3d>*>(tMaxNumADVs);
                    tADVIndexField   = Cell<stk::mesh::Field<Integer>*>(tMaxNumADVs);

                    for(Integer i = 0; i<tMaxNumADVs; i++)
                    {
                        stk::mesh::Field<Real,stk::mesh::Cartesian3d>& tADVField = tMeshMeta.declare_field<stk::mesh::Field<Real,stk::mesh::Cartesian3d> >(stk::topology::NODE_RANK, tADVFieldNames(i));
                        stk::mesh::put_field(tADVField, *tInterfaceNodes , 3,&tVal);


                        stk::mesh::Field<Integer> & tIndexField = tMeshMeta.declare_field<stk::mesh::Field<Integer> >(stk::topology::NODE_RANK, tADVIndFieldNames(i));
                        stk::mesh::put_field(tIndexField, *tInterfaceNodes , 1);
                        tInterfaceFields(i) = &tADVField;
                        tADVIndexField(i)   = &tIndexField;
                    }

                    stk::mesh::Field<Integer> & tNAIndexField = tMeshMeta.declare_field<stk::mesh::Field<Integer> >(stk::topology::NODE_RANK, aSensitivityData.get_num_adv_ind_field_name());
                    stk::mesh::put_field(tNAIndexField, *tInterfaceNodes , 1);
//                    tNADVIndexField = &tNAIndexField;
                }

            }

            // Commit MetaData before populating the BulkData
            tMeshMeta.commit();

            // Create BulkData Object
            stk::mesh::BulkData & tMeshBulk = tSTKMeshData->mesh_bulk_data();


            // Start modification cycle
            tMeshBulk.modification_begin();

            // Declare set that will contain declared nodes for current processors (required for parallel meshes)
//            std::set<stk::mesh::EntityId> tNodesIDeclared;

            Integer tElementId = 1;
            Integer tNumBuckets = aElementBuckets.size();
            stk::mesh::EntityId tNodeGlobalId;
            for(Integer iBuck = 0; iBuck < tNumBuckets; iBuck++)
            {

                if(aElementBuckets(iBuck).bucket_has_entities())
                {
                    Integer tNumElems = aElementBuckets(iBuck).get_num_entities_in_bucket();
                    tNumParts = aElementBuckets(iBuck).get_num_parts();


                    Cell<std::string> const & tPartNames = aElementBuckets(iBuck).get_part_names();

                    stk::mesh::PartVector tBucketParts(tNumParts);

                    /*
                     * Get Parts
                     */


                    for(Integer iPart =0; iPart<tPartNames.size(); iPart++)
                    {
                        tBucketParts[iPart] = tMeshMeta.get_part(tPartNames(iPart));
                    }

                    Integer tNumNodesPerElem;
                    for (Integer elem_i = 0; elem_i < tNumElems; ++elem_i )
                    {


                        moris::Mat_New<Integer, Integer_Matrix> const & tElementToNodes = aElementBuckets(iBuck).get_entity(elem_i);

                        tElementId = aElementBuckets(iBuck).get_entity_id(elem_i);
                        tNumNodesPerElem = tElementToNodes.n_cols();

                        stk::mesh::Entity tElement = tMeshBulk.declare_entity(stk::topology::ELEM_RANK,tElementId,tBucketParts);

                        for (Integer node_i = 0; node_i < tNumNodesPerElem; ++node_i )
                        {
                            tNodeGlobalId = (stk::mesh::EntityId)tElementToNodes(0,node_i);


                            if(tNodeGlobalId != 0)
                            {
                                stk::mesh::Entity tNode = tMeshBulk.get_entity( stk::topology::NODE_RANK, tNodeGlobalId );

                                if(!tMeshBulk.is_valid( tNode ))
                                {
                                    tNode = tMeshBulk.declare_entity(stk::topology::NODE_RANK, tNodeGlobalId,tMeshMeta.universal_part());
                                }

                                tMeshBulk.declare_relation(tElement, tNode, node_i);
                            }
                        }
                    }
                }
            }


            /*
             * Add node sets
             */
            for(Integer iNS = 0; iNS<aNodeSets.size(); iNS++)
            {

                Cell<std::string> const & tNodeSetNames = aNodeSets(iNS).get_node_set_names();
                tNumParts = tNodeSetNames.size();
                stk::mesh::PartVector tNodeParts(tNumParts); // +1 for topology part

                for(Integer iPart =0; iPart<tNumParts; iPart++)
                {
                    tNodeParts[iPart] = tMeshMeta.get_part(tNodeSetNames(iPart));
                }

                Cell<Integer> const & tNodeIds = aNodeSets(iNS).get_node_ids();


                for(Integer iNode = 0; iNode <tNodeIds.size(); iNode++)
                {
                    stk::mesh::Entity tNode = tMeshBulk.get_entity( stk::topology::NODE_RANK, tNodeIds(iNode) );

                    if(tMeshBulk.is_valid(tNode))
                        tMeshBulk.change_entity_parts( tNode, tNodeParts);
                }
            }


            calculate_shared_node_data(tMeshBulk, tMeshMeta);

            //Close modification cycle
            tMeshBulk.modification_end();

            // Create mesh edges (not necessary for a basic mesh)
//            stk::mesh::create_edges( tMeshBulk );
//            // Create mesh faces (not necessary for a basic mesh)
//            stk::mesh::create_faces( tMeshBulk , true ); // boolean to specify if want to connect faces to edges


            // Add side sets
            tMeshBulk.modification_begin();

            Integer tNumSideBuckets = aSideSets.size();
            Integer tNumSidesInBucket = 0;
            Integer tSideOrdinal = 10;
            for(Integer iSet = 0; iSet<tNumSideBuckets; iSet++)
            {
                std::string const & tPartNames = aSideSets(iSet).get_side_set_name();
                tNumParts = aSideSets(iSet).get_num_side_sets();
                stk::mesh::PartVector tSideParts(tNumParts); // +1 for topology part

                /*
                 * Get Parts
                 */

                for(Integer iPart =0; iPart<1; iPart++)
                {
                    tSideParts[iPart] = tMeshMeta.get_part(tPartNames);
                }

                tNumSidesInBucket = aSideSets(iSet).get_num_of_sides();

                for(Integer iSide=0; iSide<tNumSidesInBucket;iSide++)
                {

                    tElementId = aSideSets(iSet).get_element_id(iSide);
                    stk::mesh::Entity tElement = tMeshBulk.get_entity( stk::topology::ELEMENT_RANK, tElementId  );

                    if(!tMeshBulk.is_valid(tElement))
                    {
                         tElement = tMeshBulk.get_entity( stk::topology::ELEMENT_RANK, tElementId + 1  );
                    }
                    if(tMeshBulk.is_valid(tElement))
                    {
                        tSideOrdinal = aSideSets(iSet).get_side_ordinal(iSide);
//                        moris::Mat_New<Integer, Integer_Matrix> const & tElementSideNodes = aSideSets(iSet).get_side_nodes(iSide);
//
//                        stk::mesh::EntityVector tSideNodesEntityVector;
//
//                        for(Integer i=0; i<tElementSideNodes.n_cols(); i++)
//                        {
//                            stk::mesh::Entity tNode = tMeshBulk.get_entity( stk::topology::NODE_RANK, tElementSideNodes(0,i)  );
//                            tSideNodesEntityVector.push_back(tNode);
//                        }



//                        stk::mesh::OrdinalAndPermutation tFaceOrd = stk::mesh::get_ordinal_and_permutation(tMeshBulk,tElement,stk::topology::FACE_RANK, tSideNodesEntityVector );
//
//                        tSTKTopology = tMeshBulk.bucket(tElement).topology();



//                        if(tFaceOrd.first != tSideOrdinal)
//                        {
////
//                            std::cout<<"Element ID: "<< tElementId <<" STK Ordinal: "<< tFaceOrd.first << " XTK Ordinal: "<< tSideOrdinal <<"Nodes: ";
//                            for(Integer i = 0; i<tElementSideNodes.n_cols(); i++)
//                            {
//                                std::cout<< " "<< tElementSideNodes(0,i);
//                            }
//                            std::cout<<std::endl;
//                            tMeshBulk.declare_element_side(tElement,tFaceOrd.first,tSideParts);
//                        }
//
//                        else
//                        {
                            tMeshBulk.declare_element_side(tElement,tSideOrdinal,tSideParts);
//                        }
                    }

                }
            }


            calculate_shared_node_data(tMeshBulk, tMeshMeta);

            tMeshBulk.modification_end();

            // Get the coordinates field from stk
            stk::mesh::FieldBase const * coord_field_i = tMeshMeta.coordinate_field();

            Integer tNumNodes = aNodeCoordinates.n_rows();

            // Loop over the number of nodes
            if(tMapFlag == 0)
            {
                for (Integer node_i = 0; node_i < tNumNodes; ++node_i )
                {
                    // Get global Id of current node and create "node entity" for stk mesh
                    //stk::mesh::EntityId nodeGlobalId = node_i;
                    stk::mesh::Entity node = tMeshBulk.get_entity( stk::topology::NODE_RANK, node_i+1  );

                    if(tMeshBulk.is_valid(node))
                    {
                        // Store the coordinates of the current node
                        Real* coord_data = static_cast <Real*> ( stk::mesh::field_data ( * coord_field_i, node ));

                        // Add coordinates information to the BulkData
                        for (Integer dim = 0; dim < aSpatialDimension; ++dim)
                        {
                            coord_data[dim] = (Real) aNodeCoordinates( node_i , dim );
                        }
                    }
                }
            }
            else // non-contiguous
            {
                for (Integer node_i = 0; node_i < tNumNodes; ++node_i )
                {
                    Integer tId = aLocaltoGlobalNodeMap(0,node_i);

                    if(tId!=0)
                    {
                    // Get global Id of current node and create "node entity" for stk mesh
                    stk::mesh::Entity node = tMeshBulk.get_entity( stk::topology::NODE_RANK, tId );

                    // Store the coordinates of the current node
                    if(tMeshBulk.is_valid(node))
                    {

                        Real* coord_data = static_cast <Real*> ( stk::mesh::field_data ( * coord_field_i, node ));

                        // Add coordinates information to the BulkData
                        for (Integer dim = 0; dim < aSpatialDimension; ++dim)
                        {
                            coord_data[dim] = aNodeCoordinates( node_i , dim );
                        }
                    }
                    }
                }
            }

//            tMeshIo.set_bulk_data(tSTKMeshData->mesh_bulk_data());


            // Add field data
            if(aNodeSets.size() != 0 )
            {
                //            std::cout<<"Adding Field"<<std::endl;

                Integer tInterfaceNodesIndex = aNodeSets.size()-1;
                Integer tNumInterfaceNodes = aNodeSets(tInterfaceNodesIndex).get_num_nodes_in_node_set();
                Integer tNumInterfaceNodeFields = aNodeSets(tInterfaceNodesIndex).get_num_real_fields();

                for(Integer iField = 0; iField<tNumInterfaceNodeFields; iField++)
                {

                    stk::mesh::Field<Real,stk::mesh::Cartesian3d> * tField = tMeshMeta.get_field<stk::mesh::Field<Real,stk::mesh::Cartesian3d>>(stk::topology::NODE_RANK,aRealVectorFields(iField));

                    for(Integer iNode = 0; iNode<tNumInterfaceNodes; iNode++)
                    {
                        moris::Mat_New<Real,Real_Matrix> const & tNodeData = aNodeSets(tInterfaceNodesIndex).get_real_field_data(iField,iNode);
                        Integer const & tNodeId = aNodeSets(tInterfaceNodesIndex).get_node_id(iNode);
                        // Get global Id of current node and create "node entity" for stk mesh
                        stk::mesh::Entity tEntity = tMeshBulk.get_entity(stk::topology::NODE_RANK, tNodeId);

                        // Store the coordinates of the current node
                        Real* tFieldData = stk::mesh::field_data ( *tField, tEntity );
                        for(Integer iS = 0; iS<aSpatialDimension; iS++)
                        {
                            tFieldData[iS] = tNodeData(0,iS);
                        }
                    }
                }
            }



            // Sensitivity Information
            std::unordered_map<Integer, moris::Mat_New<Integer, Integer_Matrix>> const & tDxDpMap = aSensitivityData.get_full_dxdp_map();

            // Iterate through map

            if(!aSensitivityData.output_sparesly())
            {
            for(auto it = tDxDpMap.begin(); it !=tDxDpMap.end(); ++it)
            {


                // Get the node entity
                Integer const & tNodeId = it->first;
                stk::mesh::Entity tEntity = tMeshBulk.get_entity(stk::topology::NODE_RANK, tNodeId);

                // get the adv data
                moris::Mat_New<Integer, Integer_Matrix> const & tADVIndices = (it->second);
                for(Integer j = 0; j<tADVIndices.n_cols(); j++)
                {
                    Integer const & tADVIndex = tADVIndices(0,j);
                    Integer const & tADVRow   = tADVIndices(1,j);
                    moris::Mat_New<Real,Real_Matrix> const &  tDxDpData = aSensitivityData.get_sensitivity_data(tADVIndex);

                    Real* tFieldData = stk::mesh::field_data ( *tInterfaceFields(tADVIndex), tEntity );

                    for(Integer iD = 0; iD<3; iD++)
                    {
                        tFieldData[iD] = tDxDpData(tADVRow,iD);
                    }

                }
            }
            }
            else
            {
                for(auto it = tDxDpMap.begin(); it !=tDxDpMap.end(); ++it)
                {

//                    // Get the node entity
//                    Integer const & tNodeId = it->first;
////                    stk::mesh::Field<Integer> * tField = tMeshIo.meta_data().get_field<stk::mesh::Field<Integer>>(stk::topology::NODE_RANK,"dxdp_indices");
//                    stk::mesh::Entity tEntity = tMeshBulk.get_entity(stk::topology::NODE_RANK, tNodeId);
//
//                    // get the adv data
//                    moris::Mat_New<Integer, Integer_Matrix> const & tADVIndices = *(it->second);
//                    Integer* tNADVIndices = stk::mesh::field_data(*tNADVIndexField, tEntity);
//                    tNADVIndices[0] = tADVIndices.n_cols();
//                    for(Integer j = 0; j<tADVIndices.n_cols(); j++)
//                    {
//                        Integer const & tADVIndex = tADVIndices(0,j);
//                        Integer const & tADVRow   = tADVIndices(1,j);
//                        moris::Mat_New<Real,Real_Matrix> const &  tDxDpData = aSensitivityData.get_sensitivity_data(j);
//
//                        Real*    tFieldData     = stk::mesh::field_data ( *tInterfaceFields(j), tEntity );
//                        Integer* tADVIndData = stk::mesh::field_data ( *tADVIndexField(j), tEntity );
//
//                        for(Integer iD = 0; iD<3; iD++)
//                        {
//                            tFieldData[iD] = tDxDpData(tADVRow,iD);
//                        }
//
//                        tADVIndData[0] = tADVIndex;
//                    }
                }
            }



            if(aSetupDataForInternalUse)
            {
                // Setup the external entity implementation
                tSTKMeshData->setup_mesh_external_entities_and_first_available_information();

                // Tell the Mesh data to setup the information it needs
                tSTKMeshData->setup_parallel_mesh_data();
            }

            // Dynamic cast to return a pointer to the base class
            std::shared_ptr<mesh::Mesh_Data<Real, Integer, Real_Matrix, Integer_Matrix>> tMeshData = std::dynamic_pointer_cast<
                    mesh::Mesh_Data<Real, Integer, Real_Matrix, Integer_Matrix>>(tSTKMeshData);

            return tSTKMeshData;
        }


    // Private Member Functions
private:

    void get_stk_topology(enum EntityRank const & aEntityRank,
                              Integer const & aNumNodes,
                              stk::topology::topology_t & aSTKTopology) const
    {
        if(aEntityRank == EntityRank::ELEMENT && aNumNodes == 8)
        {
            aSTKTopology = stk::topology::HEX_8;
        }

        else if(aEntityRank == EntityRank::ELEMENT && aNumNodes == 4)
        {
            aSTKTopology = stk::topology::TET_4;
        }

        else
        {
            XTK_ERROR<<"Topology not implemented";
        }
    }

    void get_stk_topology(enum EntityTopology const & aEntityTopology,
                          stk::topology::topology_t & aSTKTopology) const
    {
        if(aEntityTopology == EntityTopology::HEXA_8)
        {
            aSTKTopology = stk::topology::HEX_8;
        }

        else if(aEntityTopology == EntityTopology::TET_4)
        {
            aSTKTopology = stk::topology::TET_4;
        }

        else if(aEntityTopology == EntityTopology::TET_10)
        {
            aSTKTopology = stk::topology::TET_10;
        }
        else
        {
            XTK_ERROR<<"Topology not implemented";
        }
    }

};
}

#endif
