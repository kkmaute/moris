#ifndef MORIS_cl_XTK_Cut_Integration_Mesh_HPP_
#define MORIS_cl_XTK_Cut_Integration_Mesh_HPP_

#include "cl_Cell.hpp"
#include "cl_MTK_Cell.hpp"
#include "cl_MTK_Vertex.hpp"

#include "cl_MTK_Mesh.hpp"
#include "cl_XTK_Integration_Mesh_Generator.hpp"
#include "cl_Communication_Tools.hpp"
#include <stdio.h>

namespace xtk
{
    struct IG_Cell_Group
    {
        IG_Cell_Group(moris_index aNumCellsInGroup):
        mIgCellGroup(aNumCellsInGroup,nullptr)
        {

        }

        moris::Cell<moris::mtk::Cell*> mIgCellGroup;
    };

    struct IG_Vertex_Group
    {
        IG_Vertex_Group(moris_index aNumVerticesInGroup):
        mIgVertexGroup(0,nullptr)
        {
            mIgVertexGroup.reserve(aNumVerticesInGroup);
        }

        std::size_t
        size()
        {
            return mIgVertexGroup.size();
        }

        void
        reserve(std::size_t aReserveSize)
        {
            mIgVertexGroup.reserve(aReserveSize);
            mIgVertexLocalCoords.reserve(aReserveSize);
        }

        void
        add_vertex(
            moris::mtk::Vertex* aVertex,
            std::shared_ptr<Matrix<DDRMat>> aVertexLocalCoord)
        {
            moris_index tNewVertexOrdinal = (moris_index)mIgVertexGroup.size();

            MORIS_ASSERT(mIgVertexIndexToVertexOrdinal.find(aVertex->get_index()) == mIgVertexIndexToVertexOrdinal.end(),"Duplicate vertex in group");
            MORIS_ASSERT(mIgVertexGroup.size() == mIgVertexLocalCoords.size(),"Size issue");
            mIgVertexGroup.push_back(aVertex);
            mIgVertexLocalCoords.push_back(aVertexLocalCoord);
            mIgVertexIndexToVertexOrdinal[aVertex->get_index()] = tNewVertexOrdinal;
        }   

        void
        add_vertex_local_coord_pointers()
        {}

        moris::mtk::Vertex*
        get_vertex(moris_index aGroupVertexOrdinal)
        {
            return mIgVertexGroup(aGroupVertexOrdinal);
        }

        moris_index 
        get_vertex_group_ordinal(moris_index aVertex)
        {
            auto tIter = mIgVertexIndexToVertexOrdinal.find(aVertex);
            MORIS_ERROR(tIter != mIgVertexIndexToVertexOrdinal.end(),"Provided vertex not in group");
            return tIter->second;
        }

        std::shared_ptr<Matrix<DDRMat>>
        get_vertex_local_coords(moris_index aVertex)
        {
            return mIgVertexLocalCoords(this->get_vertex_group_ordinal(aVertex));
        }
        
        private:
        moris::Cell<moris::mtk::Vertex*>             mIgVertexGroup;
        std::unordered_map<moris_index,moris_index>  mIgVertexIndexToVertexOrdinal;
        moris::Cell<std::shared_ptr<Matrix<DDRMat>>> mIgVertexLocalCoords;

    };

    struct Edge_Based_Connectivity
    {
        moris::Cell<moris::Cell<moris::mtk::Vertex*>> mEdgeVertices;
        moris::Cell<moris::Cell<moris::mtk::Cell*>>   mEdgeToCell;
        moris::Cell<moris::Cell<moris::moris_index>>  mEdgeToCellEdgeOrdinal;
        moris::Cell<moris::Cell<moris_index>>         mCellToEdge;
    };
    
    struct Edge_Based_Ancestry
    {
        moris::Cell<moris::moris_index> mEdgeParentEntityIndex;
        moris::Cell<moris::moris_index> mEdgeParentEntityRank;
        moris::Cell<moris::moris_index> mEdgeParentEntityOrdinalWrtBackgroundCell;
    };

    

    class Child_Mesh_Experimental;

    class Cut_Integration_Mesh
    {
        friend class Integration_Mesh_Generator;
    protected:
        moris::uint mSpatialDimension;

        // integration cells
        moris_index mFirstControlledCellIndex;
        moris::Cell<moris::mtk::Cell*> mIntegrationCells;
        moris::Cell<std::shared_ptr<moris::mtk::Cell>> mControlledIgCells;

        // integration vertices
        moris_index mFirstControlledVertexIndex;
        moris::Cell<moris::mtk::Vertex*> mIntegrationVertices;
        moris::Cell<std::shared_ptr<moris::mtk::Vertex>> mControlledIgVerts;

        moris::Cell<moris::moris_index> mIgVertexParentEntityIndex;
        moris::Cell<moris::moris_index> mIgVertexParentEntityRank;  
        moris::Cell<moris::moris_index> mIgVertexConnectedCell;     // this is needed to deduce ancestry for an edge between two nodes.

        // vertex coords
        moris::Cell<std::shared_ptr<Matrix<DDRMat>>> mVertexCoordinates;

        // all data is stored in the current mesh. pointers are in the child mesh
        // as well as accessor function
        moris::Cell<std::shared_ptr<Child_Mesh_Experimental>> mChildMeshes;

        // communication map
        moris::Matrix<IdMat> mCommunicationMap;

        // group of all integration cells in a single parent cell
        // outer cell is the child mesh
        moris::Cell<std::shared_ptr<IG_Cell_Group>>               mIntegrationCellGroups;
        moris::Cell<std::shared_ptr<IG_Vertex_Group>>             mIntegrationVertexGroups;
        moris::Cell<moris::mtk::Cell*>                            mIntegrationCellGroupsParentCell;

        moris::Cell<moris::Cell<moris_index>> mIntegrationCellToCellGroupIndex;

        
        // bulk phase
        moris::Cell<moris::moris_index> mIntegrationCellBulkPhase;

        // subphase groupings
        // outer cell is the child mesh
        moris::Cell<moris::moris_index>             mSubphaseIds;
        moris::Cell<moris::moris_index>             mSubphaseBulkPhase;
        moris::Cell<moris::Cell<moris::mtk::Cell*>> mIntegrationCellSubphaseGroups;
        moris::Cell<moris::mtk::Cell*>              mIntegrationCellSubphaseGroupsParentCell;
        moris::Cell<moris::moris_index>             mParentCellSubphaseIndices;

        // vertex ancestry
        moris::Cell<moris::Cell<moris_index>> mVertexParentIndex;
        moris::Cell<moris::Cell<moris_index>> mVertexParentRank;

        // connectivity from vertex to child mesh
        // outer cell vertex
        // inner cell child meshes associated with the vertex
        moris::Cell<moris::Cell<moris_index>> mVertexToChildMeshIndex;

        // connectivity from vertex to child mesh
        // outer cell integration cell index
        // inner cell child meshes associated with the integration cell
        moris::Cell<moris::Cell<moris_index>> mCellToChildMeshIndex;

        // outer cell geomtry index
        // if vertex index is in the map then it is a member of the geometric interface
        moris::Cell<std::unordered_map<moris_index,moris_index>> mGeometryInterfaceVertexIndices;

        moris::moris_index mGlobalMaxVertexId;
        moris::moris_index mGlobalMaxCellId;

        std::unordered_map< moris_id, moris_index > mIntegrationCellIdToIndexMap;
        std::unordered_map< moris_id, moris_index > mIntegrationVertexIdToIndexMap;

        moris::Cell<moris::moris_index> mIntegrationCellIndexToId;
        moris::Cell<moris::moris_index> mIntegrationVertexIndexToId;

        moris::mtk::Mesh* mBackgroundMesh;

        std::shared_ptr<moris::mtk::Cell_Info> mChildCellInfo;

    
    public:
        Cut_Integration_Mesh(moris::mtk::Mesh* aBackgroundMesh);
        ~Cut_Integration_Mesh();

        moris::Cell<std::shared_ptr<Matrix<DDRMat>>>*
        get_all_vertex_coordinates_loc_inds()
        {
            return &mVertexCoordinates;
        }
        
        std::shared_ptr<Child_Mesh_Experimental>
        get_child_mesh(moris_index aChildMeshIndex)
        {
            MORIS_ERROR(mChildMeshes.size() > (uint)aChildMeshIndex,"Child mesh index out of bounds");
            return mChildMeshes(aChildMeshIndex);
        }

        mtk::Cell const &
        get_mtk_cell( moris_index aElementIndex ) const
        {
          return *mIntegrationCells(aElementIndex);
        }

        moris::mtk::Vertex*
        get_mtk_vertex_pointer(moris_index aVertexIndex)
        {
            return mIntegrationVertices(aVertexIndex);
        }

        std::shared_ptr<IG_Vertex_Group>
        get_vertex_group(moris_index aVertexGroupIndex)
        {
            return mIntegrationVertexGroups(aVertexGroupIndex);
        }

        void
        set_integration_cell(
            moris_index                       aCellIndex, 
            std::shared_ptr<moris::mtk::Cell> aNewCell,
            bool aNewCellBool )
        {
            if(!aNewCell)
            {
                MORIS_ERROR(aCellIndex >= mFirstControlledCellIndex,"Cannot set integration cell that I do not control."); 
            }

            // controlled index
            moris_index tIndexInControlledCells = aCellIndex - mFirstControlledCellIndex;

            mControlledIgCells(tIndexInControlledCells) = aNewCell;
            mIntegrationCells(aCellIndex) = mControlledIgCells(tIndexInControlledCells).get();
        }

        void
        add_cell_to_integration_mesh(
            moris_index aCellIndex, 
            moris_index aCellGroupIndex )
        {
            MORIS_ERROR(aCellGroupIndex   < (moris_index)mIntegrationCellGroups.size(),"Child Mesh Index out of bounds."); 
            MORIS_ERROR(aCellIndex < (moris_index)mIntegrationCells.size(),"Cell Index out of bounds."); 
            mIntegrationCellGroups(aCellGroupIndex)->mIgCellGroup.push_back(mIntegrationCells(aCellIndex));
            mIntegrationCellToCellGroupIndex(aCellIndex).push_back(aCellGroupIndex);
        }


        moris_id
        allocate_entity_ids( moris::size_t   aNumIdstoAllocate,
                             enum EntityRank aEntityRank)
        {
            MORIS_ERROR(aEntityRank == EntityRank::NODE || aEntityRank == EntityRank::ELEMENT,"Only node and element ids can be allocated with xtk.");
            int tProcRank = moris::par_rank();
            int tProcSize = moris::par_size();

            // size_t is defined as uint here because of aNumRequested
            //Initialize gathered information outputs (information which will be scattered across processors)
            moris::Cell<moris::moris_id> aGatheredInfo;
            moris::Cell<moris::moris_id> tFirstId(1);
            moris::Cell<moris::moris_id> tNumIdsRequested(1);

            tNumIdsRequested(0) = (moris::moris_id)aNumIdstoAllocate;

            moris::gather(tNumIdsRequested,aGatheredInfo);

            moris::Cell<moris::moris_id> tProcFirstID(tProcSize);


            if(tProcRank == 0)
            {
                // Loop over entities print the number of entities requested by each processor
                for (int iProc = 0; iProc < tProcSize; ++iProc)
                {
                    if(aEntityRank == EntityRank::NODE)
                    {
                        // Give each processor their desired amount of IDs
                        tProcFirstID(iProc) = mGlobalMaxVertexId;

                        // Increment the first available node ID
                        mGlobalMaxVertexId= mGlobalMaxVertexId+aGatheredInfo(iProc);
                    }

                    else if(aEntityRank == EntityRank::ELEMENT)
                    {
                        // Give each processor their desired amount of IDs
                        tProcFirstID(iProc) = mGlobalMaxCellId;

                        // Increment the first available node ID
                        mGlobalMaxCellId= mGlobalMaxCellId+aGatheredInfo(iProc);
                    }
                }
            }

            moris::scatter(tProcFirstID,tFirstId);

            return tFirstId(0);
        }

        moris::moris_index
        get_first_available_index(enum EntityRank aEntityRank) const
        {
            MORIS_ERROR(aEntityRank == EntityRank::NODE || aEntityRank == EntityRank::ELEMENT,"Only can handle this question for nodes and elements");

            if(aEntityRank == EntityRank::NODE)
            {
                return mIntegrationVertices.size();
            }

            if(aEntityRank == EntityRank::ELEMENT)
            {
                return mIntegrationCells.size();
            }

            return MORIS_INDEX_MAX;
        }

        std::shared_ptr<IG_Cell_Group>
        get_ig_cell_group(moris_index aGroupIndex)
        {
            return mIntegrationCellGroups(aGroupIndex);
        }

        moris::Cell<moris_index> const &
        get_ig_cell_group_memberships(moris_index aIgCellIndex)
        {
            return mIntegrationCellToCellGroupIndex(aIgCellIndex);
        }

        moris::mtk::Cell*
        get_ig_cell_group_parent_cell(moris_index aGroupIndex)
        {
            return mIntegrationCellGroupsParentCell(aGroupIndex);
        }

        uint
        get_num_entities( 
            enum EntityRank   aEntityRank, 
            const moris_index aIndex ) const
        {
            switch(aEntityRank)
            {
                case EntityRank::NODE:
                {
                    return mIntegrationVertices.size();
                    break;
                }
                case EntityRank::ELEMENT:
                {
                    return mIntegrationCells.size();
                    break;
                }
                default:
                {
                    MORIS_ERROR(0,"Only support get num entities for nodes and elements currently");
                }
                return 0;
            }
        }

    private:

        void
        setup_comm_map()
        {
            if(par_size() > 1)
            {
                std::unordered_map<moris_id,moris_id> tCommunicationMap;
                Cell<moris_index> tCellOfProcs;
        
                for(moris::uint i = 0; i < mBackgroundMesh->get_num_entities(EntityRank::NODE); i++)
                {
                    moris_index tOwner = mBackgroundMesh->get_entity_owner((moris_index)i,EntityRank::NODE);

                    if(tCommunicationMap.find(tOwner) == tCommunicationMap.end() && tOwner != par_rank())
                    {
                        tCellOfProcs.push_back(tOwner);
                        tCommunicationMap[tOwner] = 1;
                    }
                }
                mCommunicationMap.resize(1,tCellOfProcs.size());
                for(moris::uint i = 0; i < tCellOfProcs.size(); i++)
                {
                    mCommunicationMap(i) = tCellOfProcs(i);
                }

                Cell<Matrix<IndexMat>> tGatheredMats;
                moris_index tTag = 10009;
                if(tCellOfProcs.size() == 0)
                {
                    Matrix<IndexMat> tDummy(1,1,MORIS_INDEX_MAX);
                    all_gather_vector(tDummy,tGatheredMats,tTag,0,0);
                }
                else
                {
                    all_gather_vector(mCommunicationMap,tGatheredMats,tTag,0,0);
                }
                if(par_rank() == 0)
                {
                    Cell<Cell<uint>> tProcToProc(par_size());
                    for(moris::uint i = 0; i < tGatheredMats.size(); i++)
                    {
                        for(moris::uint j = 0; j < tGatheredMats(i).numel(); j ++)
                        {
                            if(tGatheredMats(i)(j) != MORIS_INDEX_MAX)
                            {
                                tProcToProc(tGatheredMats(i)(j)).push_back((moris_index)i);
                            }
                        }
                    }
                    Cell<Matrix<IndexMat>> tReturnMats(tProcToProc.size());
                    // convert to a matrix
                    for(moris::uint  i = 0; i < tProcToProc.size(); i++)
                    {
                        tReturnMats(i).resize(1,tProcToProc(i).size());

                        for(moris::uint j = 0; j < tProcToProc(i).size(); j++)
                        {
                            tReturnMats(i)(j) = tProcToProc(i)(j);
                        }
                        if(tProcToProc(i).size() == 0)
                        {
                            tReturnMats(i) = Matrix<IndexMat>(1,1,MORIS_INDEX_MAX);
                        }
                    }
                    // send them back
                    for(moris::uint i = 0; i < tReturnMats.size(); i++)
                    {
                        nonblocking_send(tReturnMats(i),tReturnMats(i).n_rows(),tReturnMats(i).n_cols(),i,tTag);
                    }
                }

                barrier();
                Matrix<IndexMat> tTempCommMap(1,1,0);
                receive(tTempCommMap,1,0,tTag);

                for(moris::uint i =0; i < tTempCommMap.numel(); i++)
                {
                    if(tTempCommMap(i) != MORIS_INDEX_MAX)
                    {
                    if(tCommunicationMap.find(tTempCommMap(i)) == tCommunicationMap.end() && tTempCommMap(i) != par_rank())
                    {
                        moris_index tIndex = mCommunicationMap.numel();
                        mCommunicationMap.resize(1,mCommunicationMap.numel()+1);

                        tCommunicationMap[tTempCommMap(i)] = 1;

                        mCommunicationMap(tIndex) = tTempCommMap(i);
                    }
                    }
                }

                barrier();
                // every proc tell the other procs that they communicate with them
            }
        }
    };
    
    Cut_Integration_Mesh::Cut_Integration_Mesh(moris::mtk::Mesh* aBackgroundMesh):mBackgroundMesh(aBackgroundMesh)
    {
        //
        mSpatialDimension = aBackgroundMesh->get_spatial_dim();

        // setup the integration cells
        moris::uint tNumBackgroundCells    = mBackgroundMesh->get_num_elems();
        moris::uint tNumBackgroundVertices = mBackgroundMesh->get_num_nodes();

        // cell setup
        mIntegrationCells.resize( tNumBackgroundCells, nullptr);
        mIntegrationCellIndexToId.resize( tNumBackgroundCells, MORIS_INDEX_MAX);
        mIntegrationCellToCellGroupIndex.resize(tNumBackgroundCells,0);
        for(moris::uint iCell = 0; iCell< tNumBackgroundCells; iCell++)
        {
            mIntegrationCells(iCell) = &mBackgroundMesh->get_mtk_cell((moris_index)iCell);
        }
        
        // vertex setup
        mIntegrationVertices.resize(tNumBackgroundVertices,nullptr);
        mIntegrationVertexIndexToId.resize(tNumBackgroundVertices, MORIS_INDEX_MAX);
        mVertexCoordinates.resize(tNumBackgroundVertices,nullptr);
        mIgVertexParentEntityIndex.resize(tNumBackgroundVertices);
        mIgVertexParentEntityRank.resize(tNumBackgroundVertices,0);

        for(moris::uint iV = 0; iV< tNumBackgroundVertices; iV++)
        {
            // get a vertex pointer into our data
            mIntegrationVertices(iV) = &mBackgroundMesh->get_mtk_vertex( (moris_index) iV );

            std::cout<<"iV = "<<iV << " | mIntegrationVertices(iV) index ="<<mIntegrationVertices(iV)->get_index()<<std::endl;

            // this vertex's parent is itself
            mIgVertexParentEntityIndex(iV) = iV;
            
            // check the vertex index lines up correctly
            MORIS_ERROR(mIntegrationVertices(iV)->get_index() == (moris_index) iV,"Vertex index mismatch");

            // add the vertex to the local to global vertex map
            mIntegrationVertexIndexToId(mIntegrationVertices(iV)->get_index()) = mIntegrationVertices(iV)->get_id();

            // store the coordinate
            mVertexCoordinates(mIntegrationVertices(iV)->get_index()) = std::make_shared<moris::Matrix<DDRMat>>(mIntegrationVertices(iV)->get_coords());

            // verify that we are not doubling up vertices in the id map
            MORIS_ERROR(mIntegrationVertexIdToIndexMap.find(mIntegrationVertices(iV)->get_id()) == mIntegrationVertexIdToIndexMap.end(), "Provided Vertex Id is already in the integration vertex map: Vertex Id =%uon process %u", mIntegrationVertices(iV)->get_id(), par_rank());

            // add the vertex to id to index map
            mIntegrationVertexIdToIndexMap[mIntegrationVertices(iV)->get_id()] = (moris_index) iV;
        }

        this->setup_comm_map();

        // max vertex and cell ids (to allocate new ones later)
        mGlobalMaxVertexId = mBackgroundMesh->get_max_entity_id(EntityRank::NODE)+1;
        mGlobalMaxCellId   = mBackgroundMesh->get_max_entity_id(EntityRank::ELEMENT)+1;
    
        mFirstControlledCellIndex = mIntegrationCells.size();
        mFirstControlledVertexIndex = mIntegrationVertices.size();


    }
    
    
    Cut_Integration_Mesh::~Cut_Integration_Mesh()
    {
    }
  

    
    class Child_Mesh_Experimental
    {
        friend class Cut_Integration_Mesh;
        friend class Integration_Mesh_Generator;
        public:
        moris::mtk::Cell*                mParentCell;
        moris::moris_index               mChildMeshIndex;
        std::shared_ptr<IG_Cell_Group>   mIgCells;
        std::shared_ptr<IG_Vertex_Group> mIgVerts;

        public:
        Child_Mesh_Experimental()
        {

        }


        moris::mtk::Cell*
        get_parent_cell()
        {
            return mParentCell;
        }

        moris::moris_index
        get_child_mesh_index()
        {
            return mChildMeshIndex;
        }
    };

}

#endif