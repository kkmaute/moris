#ifndef MORIS_cl_XTK_Cut_Integration_Mesh_HPP_
#define MORIS_cl_XTK_Cut_Integration_Mesh_HPP_

#include "cl_Cell.hpp"
#include "cl_MTK_Cell.hpp"
#include "cl_MTK_Vertex.hpp"
#include "cl_MTK_Writer_Exodus.hpp"

#include "cl_MTK_Mesh_Core.hpp"
#include "cl_Mesh_Enums.hpp"
#include "cl_MTK_Integration_Mesh.hpp"
#include "cl_MTK_Vertex_XTK_Impl.hpp"
#include "cl_Tracer.hpp"

#include "cl_Communication_Tools.hpp"
#include <stdio.h>
using namespace moris;

namespace xtk
{
    struct IG_Cell_Group
    {
        IG_Cell_Group(moris_index aNumCellsInGroup);

        IG_Cell_Group();

        moris::Cell<moris::mtk::Cell*> mIgCellGroup;
    };

    struct IG_Cell_Side_Group
    {
        IG_Cell_Side_Group(moris_index aEstimatedNumCells);

        moris::Cell<moris::mtk::Cell*> mIgCells;
        moris::Cell<moris_index>       mIgCellSideOrdinals;
    };

    struct IG_Cell_Double_Side_Group
    {
        IG_Cell_Double_Side_Group(moris_index aEstimatedNumCells);

        moris::Cell<moris::mtk::Cell*> mLeaderIgCells;
        moris::Cell<moris_index>       mLeaderIgCellSideOrdinals;

        moris::Cell<moris::mtk::Cell*> mFollowerIgCells;
        moris::Cell<moris_index>       mFollowerIgCellSideOrdinals;

        void
        print()
        {
            std::cout<<"Number of Leaders:   "<<mLeaderIgCells.size()<<std::endl;
            std::cout<<"Number of Followers: "<<mFollowerIgCells.size()<<std::endl;

            int tStrLen = std::string("Lead Cell Id   | ").size();

            std::cout<<"Lead Cell Id   | ";
            std::cout<<"Side Ord       | ";
            std::cout<<"Follow Cell Id | ";
            std::cout<<"Side Ord       "<<std::endl;
            // iterate through pairs
            for(moris::uint i = 0; i < mLeaderIgCells.size(); i++)
            {
                std::cout<<std::setw(tStrLen)<<mLeaderIgCells(i)->get_id();
                std::cout<<std::setw(tStrLen)<<mLeaderIgCellSideOrdinals(i);
                std::cout<<std::setw(tStrLen)<<mFollowerIgCells(i)->get_id();
                std::cout<<std::setw(tStrLen)<<mFollowerIgCellSideOrdinals(i)<<std::endl;
            }
        }
    };

    struct IG_Vertex_Group
    {
        IG_Vertex_Group(moris_index aNumVerticesInGroup);

        std::size_t
        size();

        void
        reserve(std::size_t aReserveSize);

        void
        add_vertex(
            moris::mtk::Vertex* aVertex,
            std::shared_ptr<Matrix<DDRMat>> aVertexLocalCoord);

        void
        add_vertex_local_coord_pointers();

        moris::mtk::Vertex*
        get_vertex(moris_index aGroupVertexOrdinal);

        moris_index 
        get_vertex_group_ordinal(moris_index aVertex);

        std::shared_ptr<Matrix<DDRMat>>
        get_vertex_local_coords(moris_index aVertex);

        void
        print();
        
        
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

    struct Facet_Based_Connectivity
    {
        moris::Cell<moris::Cell<moris::mtk::Vertex*>> mFacetVertices;
        moris::Cell<moris::Cell<moris::mtk::Cell*>>   mFacetToCell;
        moris::Cell<moris::Cell<moris::moris_index>>  mFacetToCellEdgeOrdinal;
        moris::Cell<moris::Cell<moris_index>>         mCellToFacet;
    };
    
    struct Facet_Based_Ancestry
    {
        moris::Cell<moris::moris_index> mFacetParentEntityIndex;
        moris::Cell<moris::moris_index> mFacetParentEntityRank;
        moris::Cell<moris::moris_index> mFacetParentEntityOrdinalWrtBackgroundCell;
    };

    struct Cell_Neighborhood_Connectivity
    {
        moris::Cell<std::shared_ptr<moris::Cell<moris::mtk::Cell*>>> mNeighborCells;
        moris::Cell<std::shared_ptr<moris::Cell<moris_index>>>       mMySideOrdinal;
        moris::Cell<std::shared_ptr<moris::Cell<moris_index>>>       mNeighborSideOrdinal;
    };

    struct Subphase_Neighborhood_Connectivity
    {
        moris::Cell<std::shared_ptr<moris::Cell<moris_index>>> mSubphaseToSubPhase;
        moris::Cell<std::shared_ptr<moris::Cell<moris_index>>> mSubphaseToSubPhaseMySideOrds;
        moris::Cell<std::shared_ptr<moris::Cell<moris_index>>> mSubphaseToSubPhaseNeighborSideOrds;
        moris::Cell<std::shared_ptr<moris::Cell<moris_index>>> mTransitionNeighborCellLocation;

        void
        print_subphase_neighborhood()
        {

            std::cout<<"Subphases"<<std::endl;
            for(moris::uint iC = 0; iC<mSubphaseToSubPhase.size(); iC++ )
            {
                std::cout<<std::setw(6)<<iC<<" | ";

                for(moris::uint iN = 0; iN< mSubphaseToSubPhase(iC)->size(); iN++)
                {
                    std::cout<<std::setw(6)<<(*mSubphaseToSubPhase(iC))(iN);
                }
                std::cout<<std::endl;
            }

            std::cout<<"Subphases My Side Ordinals"<<std::endl;
            for(moris::uint iC = 0; iC<mSubphaseToSubPhaseMySideOrds.size(); iC++ )
            {
                std::cout<<std::setw(6)<<iC<<" | ";

                for(moris::uint iN = 0; iN< mSubphaseToSubPhaseMySideOrds(iC)->size(); iN++)
                {
                    std::cout<<std::setw(6)<<(*mSubphaseToSubPhaseMySideOrds(iC))(iN);
                }
                std::cout<<std::endl;
            }

            std::cout<<"Subphases Neighbor Side Ordinals"<<std::endl;
            for(moris::uint iC = 0; iC<mSubphaseToSubPhaseNeighborSideOrds.size(); iC++ )
            {
                std::cout<<std::setw(6)<<iC<<" | ";

                for(moris::uint iN = 0; iN< mSubphaseToSubPhaseNeighborSideOrds(iC)->size(); iN++)
                {
                    std::cout<<std::setw(6)<<(*mSubphaseToSubPhaseNeighborSideOrds(iC))(iN);
                }
                std::cout<<std::endl;
            }


            std::cout<<"Transition Neighbor Locations"<<std::endl;
            for(moris::uint iC = 0; iC<mTransitionNeighborCellLocation.size(); iC++ )
            {
                std::cout<<std::setw(6)<<iC<<" | ";

                for(moris::uint iN = 0; iN< mTransitionNeighborCellLocation(iC)->size(); iN++)
                {
                    std::cout<<std::setw(12)<<(*mTransitionNeighborCellLocation(iC))(iN);
                }
                std::cout<<std::endl;
            }
        }
    };
    

    

    class Child_Mesh_Experimental;
    class Model;
    class Cell_XTK_No_CM;
    class Cut_Integration_Mesh : public moris::mtk::Mesh
    {
        friend class Integration_Mesh_Generator;
    protected:
        moris::uint mSpatialDimension;

        // integration cells
        moris_index mFirstControlledCellIndex;
        moris::Cell<moris::mtk::Cell*> mIntegrationCells;
        moris::Cell<std::shared_ptr<xtk::Cell_XTK_No_CM>> mControlledIgCells;

        // quantities related to integration cells
        moris::Cell<moris::Cell<moris_index>> mIntegrationCellToCellGroupIndex;
        moris::Cell<moris_index>              mIntegrationCellToSubphaseIndex;
        moris::Cell<moris::moris_index>       mIntegrationCellBulkPhase;

        // integration vertices
        moris_index mFirstControlledVertexIndex;
        moris::Cell<moris::mtk::Vertex*> mIntegrationVertices;
        moris::Cell<std::shared_ptr<moris::mtk::Vertex_XTK>> mControlledIgVerts;

        // vertex ancestry
        moris::Cell<moris::moris_index> mIgVertexParentEntityIndex;
        moris::Cell<moris::moris_index> mIgVertexParentEntityRank;  
        moris::Cell<moris::moris_index> mIgVertexConnectedCell;

        // vertex quantities
        moris::Cell<std::shared_ptr<Matrix<DDRMat>>> mVertexCoordinates;

        // all data is stored in the current mesh. pointers are in the child mesh
        // as well as accessor functions are provided there
        moris::Cell<std::shared_ptr<Child_Mesh_Experimental>> mChildMeshes;
        moris::Cell<std::shared_ptr<Child_Mesh_Experimental>> mOwnedChildMeshes;
        moris::Cell<std::shared_ptr<Child_Mesh_Experimental>> mNotOwnedChildMeshes;

        // communication map
        moris::Matrix<IdMat> mCommunicationMap;

        // group of all integration cells in a single parent cell
        moris::Cell<std::shared_ptr<IG_Cell_Group>>   mIntegrationCellGroups;
        moris::Cell<std::shared_ptr<IG_Vertex_Group>> mIntegrationVertexGroups;
        moris::Cell<moris::mtk::Cell*>                mIntegrationCellGroupsParentCell;
        moris::Cell<moris_index>                      mParentCellCellGroupIndex;
        
        moris::Cell<moris_index> mOwnedIntegrationCellGroupsInds;
        moris::Cell<moris_index> mNotOwnedIntegrationCellGroups;

        // subphase groupings
        moris::Cell<moris_index>                    mSubPhaseIds;
        moris::Cell<std::shared_ptr<IG_Cell_Group>> mSubPhaseCellGroups;
        moris::Cell<moris::moris_index>             mSubPhaseBulkPhase;
        moris::Cell<moris::mtk::Cell*>              mSubPhaseParentCell;
        moris::Cell<moris::Cell<moris_index>>       mParentCellToSubphase;
        moris::Cell<moris_index>                    mParentCellHasChildren;

        moris::Cell<moris_index> mOwnedSubphaseGroupsInds;
        moris::Cell<moris_index> mNotOwnedSubphaseGroupsInds;
        std::unordered_map<moris::moris_id,moris::moris_index> mGlobalToLocalSubphaseMap;

        std::shared_ptr<Subphase_Neighborhood_Connectivity> mSubphaseNeighborhood;

        // face connectivity
        std::shared_ptr<Facet_Based_Connectivity> mIgCellFaceConnectivity;

        // interface facets - indexed based on mIgCellFaceConnectivity facet indices
        moris::Cell<moris_index> mInterfaceFacets;
        
        // double side interface groups
        // outer cell - bulk phase 0
        // inner cell - bulk phase 1
        // IG_Cell_Double_Side_Group pairings between integration cells
        moris::Cell<moris::Cell<std::shared_ptr<IG_Cell_Double_Side_Group>>> mBptoBpDblSideInterfaces;

        // background facet to child facet connectivity
        moris::Cell<std::shared_ptr<moris::Cell<moris::moris_index>>> mBGFacetToChildFacet;

        // block set data
        std::unordered_map<std::string, moris_index> mBlockSetLabelToOrd;
        moris::Cell<std::string>                     mBlockSetNames;
        moris::Cell<std::shared_ptr<IG_Cell_Group>>  mBlockSetCellGroup;
        moris::Cell<enum CellTopology>               mBlockCellTopo;

        // Side Set Data
        std::unordered_map<std::string, moris_index>     mSideSideSetLabelToOrd;
        Cell<std::string>                                mSideSetLabels;
        moris::Cell<std::shared_ptr<IG_Cell_Side_Group>> mSideSetCellSides;


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
        Model*            mXTKModel;
        std::shared_ptr<moris::mtk::Cell_Info> mChildCellInfo;

    
    public:
        Cut_Integration_Mesh(moris::mtk::Mesh* aBackgroundMesh,
                             Model*            aXTKModel);

        ~Cut_Integration_Mesh();

        // Core Mesh Functions
        uint 
        get_spatial_dim() const
        {
            return mSpatialDimension;
        }

        MeshType 
        get_mesh_type() const
        {
            return MeshType::XTK;
        }

        uint
        get_num_entities( 
            enum EntityRank   aEntityRank, 
            const moris_index aIndex ) const;

        moris::uint get_num_sets() const
        {
            return 0;
        }

        Matrix< DDRMat >
        get_node_coordinate( moris_index aNodeIndex ) const 
        {
            return (*mVertexCoordinates(aNodeIndex));
        }

        uint 
        get_node_owner(moris_index aNodeIndex) const
        {
            return mIntegrationVertices(aNodeIndex)->get_owner();
        }

        uint 
        get_element_owner(moris_index aElementIndex) const 
        {
            return mIntegrationCells(aElementIndex)->get_owner();
        }

        Matrix< IdMat >
        get_communication_table() const
        {
            return mCommunicationMap;
        }

        Matrix<IndexMat> get_element_indices_in_block_set(uint aSetIndex)
        {
            // get cells in set
            std::shared_ptr<IG_Cell_Group> tCellsInBlock = mBlockSetCellGroup(aSetIndex);

            Matrix<IndexMat> tCellIndices(1,tCellsInBlock->mIgCellGroup.size());

            for(moris::uint iCell = 0; iCell < tCellsInBlock->mIgCellGroup.size(); iCell++)
            {
                tCellIndices(iCell) = tCellsInBlock->mIgCellGroup(iCell)->get_index();
            }

            return tCellIndices;
        }


        enum CellTopology
        get_blockset_topology(const std::string & aSetName)
        {
            // getindex
            moris_index tBlockIndex = this->get_block_set_index(aSetName);

            return mBlockCellTopo(tBlockIndex);
        }


        enum CellShape
        get_IG_blockset_shape(const std::string & aSetName)
        {
            return CellShape::INVALID;
        }


        enum CellShape
        get_IP_blockset_shape(const std::string & aSetName)
        {
            return CellShape::INVALID;
        }


        moris_id
        get_glb_entity_id_from_entity_loc_index(
                    moris_index        aEntityIndex,
                    enum EntityRank    aEntityRank,
                    const moris_index  aDiscretizationIndex = 0) const
        {
            MORIS_ERROR(aEntityRank == EntityRank::NODE || aEntityRank == EntityRank::ELEMENT,"Only supported for nodes and cells");
            if(aEntityRank == EntityRank::NODE)
            {
                return mIntegrationVertices(aEntityIndex)->get_index();
            }
            else if(aEntityRank == EntityRank::ELEMENT)
            {
                return mIntegrationCells(aEntityIndex)->get_index();
            }

            else
            {
                return 0;
            }
        }

        moris_index
        get_loc_entity_ind_from_entity_glb_id(
                moris_id        aEntityId,
                enum EntityRank aEntityRank) const
        {
            // warning element map is set up after integration mesh has been constructed
            MORIS_ERROR(aEntityRank == EntityRank::NODE || aEntityRank == EntityRank::ELEMENT,"Only a node map and element map is implemented in XTK");

            if(aEntityRank == EntityRank::NODE)
            {
                auto tIter = mIntegrationVertexIdToIndexMap.find(aEntityId);

                MORIS_ERROR(tIter!=mIntegrationVertexIdToIndexMap.end(),
                    "Provided Entity Id is not in the map, Has the map been initialized?: aEntityId =%u EntityRank = %u on process %u",
                    aEntityId, (uint)aEntityRank, par_rank());
                return tIter->second;
            }
            else if(aEntityRank == EntityRank::ELEMENT)
            {
                auto tIter = mIntegrationCellIdToIndexMap.find(aEntityId);

                MORIS_ERROR(tIter!=mIntegrationCellIdToIndexMap.end(),
                    "Provided Entity Id is not in the map, Has the map been initialized?: aEntityId =%u EntityRank = %u on process %u",
                    aEntityId, (uint)aEntityRank, par_rank());
                return tIter->second;
            }

            else
            {
                return 0;
            }

        }

        Matrix<IndexMat>
        get_entity_connected_to_entity_loc_inds(
                    moris_index        aEntityIndex,
                    enum EntityRank    aInputEntityRank,
                    enum EntityRank    aOutputEntityRank,
                    const moris_index  aDiscretizationIndex = 0) const
        {
            MORIS_ERROR(aInputEntityRank == EntityRank::ELEMENT && aOutputEntityRank == EntityRank::NODE,
                "Only support element to node connectivity");

            return this->get_mtk_cell(aEntityIndex).get_vertex_inds();
        }


        moris::Cell<std::string>
        get_set_names(enum EntityRank aSetEntityRank) const
        {
            switch(aSetEntityRank)
            {
                case EntityRank::NODE:
                {
                    return {};
                    break;
                }
                case EntityRank::EDGE:
                {
                    MORIS_ASSERT(this->get_facet_rank() == EntityRank::EDGE,"side sets are defined on edges in 2d");
                    return mSideSetLabels;
                    break;
                }
                case EntityRank::FACE:
                {
                    MORIS_ASSERT(this->get_facet_rank() == EntityRank::FACE,"side sets are defined on faces in 3d");
                    return mSideSetLabels;
                    break;
                }
                case EntityRank::ELEMENT:
                {
                    return mBlockSetNames;
                    break;
                }
                default:
                {
                    MORIS_ERROR(0,"Currently only supporting block, node and side sets in XTK enriched integration meshes");
                }
                return moris::Cell<std::string>(0);
                break;
            }
        }

        moris_index
        get_block_set_index(std::string aBlockSetLabel) const
        {
            auto tIter = mBlockSetLabelToOrd.find(aBlockSetLabel);

            MORIS_ERROR(tIter != mBlockSetLabelToOrd.end(),"block set set label not found");

            return tIter->second;
        }

        Matrix< IndexMat >
        get_block_entity_loc_inds( std::string     aSetName) const
        {
            // ge tindex
            moris_index tBlockIndex = this->get_block_set_index(aSetName);

            // get cells in set
            std::shared_ptr<IG_Cell_Group> tCellsInBlock = mBlockSetCellGroup(tBlockIndex);

            Matrix<IndexMat> tCellIndices(1,tCellsInBlock->mIgCellGroup.size());

            for(moris::uint iCell = 0; iCell < tCellsInBlock->mIgCellGroup.size(); iCell++)
            {
                tCellIndices(iCell) = tCellsInBlock->mIgCellGroup(iCell)->get_index();
            }

            return tCellIndices;
        }

        moris_index
        get_side_set_index(std::string aSideSetLabel) const
        {
            auto tIter = mSideSideSetLabelToOrd.find(aSideSetLabel);

            MORIS_ERROR(tIter != mSideSideSetLabelToOrd.end(),"side side set label not found");

            return tIter->second;
        }

        void
        get_sideset_elems_loc_inds_and_ords(
            const  std::string & aSetName,
            Matrix< IndexMat > & aElemIndices,
            Matrix< IndexMat > & aSidesetOrdinals ) const
        {
            // get the index
            moris_index tSideSetIndex = this->get_side_set_index(aSetName);

            std::shared_ptr<IG_Cell_Side_Group> tCellSidesInSet = mSideSetCellSides(tSideSetIndex);

            // iterate through side clusters and count number of sides in set
            moris::uint tNumSides = tCellSidesInSet->mIgCells.size();

            // size outputs
            aElemIndices.resize(1,tNumSides);
            aSidesetOrdinals.resize(1,tNumSides);

            for(moris::uint iSide =0; iSide < tNumSides; iSide++)
            {
                aElemIndices(iSide) = tCellSidesInSet->mIgCells(iSide)->get_index();
                aSidesetOrdinals(iSide) = tCellSidesInSet->mIgCellSideOrdinals(iSide);
            }

        }
        

        Matrix< IndexMat >
        get_set_entity_loc_inds(
                    enum EntityRank aSetEntityRank,
                    std::string     aSetName) const
        {
        switch(aSetEntityRank)
        {
            case EntityRank::NODE:
            {
                // // get the vertex set index
                // auto tSetIndex = mVertexSetLabelToOrd.find(aSetName);

                // moris::Cell<moris::mtk::Vertex*> tVerticesInSet = mVerticesInVertexSet(tSetIndex->second);
                // Matrix<IndexMat> tVerticesInSetMat(1,tVerticesInSet.size());
                // for(moris::uint i = 0; i < tVerticesInSet.size(); i++)
                // {
                //     tVerticesInSetMat(i) = tVerticesInSet(i)->get_index();
                // }

                return Matrix< IndexMat >(0,0);
                break;
            }
            case EntityRank::EDGE:
            {
                // MORIS_ASSERT(this->get_facet_rank() == EntityRank::EDGE,"side sets are defined on edges in 2d");
                return Matrix< IndexMat >(0,0);
                break;
            }
            case EntityRank::FACE:
            {
                // MORIS_ASSERT(this->get_facet_rank() == EntityRank::FACE,"side sets are defined on faces in 3d");
                return Matrix< IndexMat >(0,0);
                break;
            }
            case EntityRank::ELEMENT:
            {
                return this->get_block_entity_loc_inds(aSetName);
                return Matrix< IndexMat >(0,0);
                break;
            }
            default:
            {
                MORIS_ERROR(0,"Currently only supporting block, node and side sets in XTK enriched integration meshes");
                return Matrix< IndexMat >(0,0);
                break;
            }
        }
    }
        

        moris::Cell<std::shared_ptr<Matrix<DDRMat>>>*
        get_all_vertex_coordinates_loc_inds();

        moris::Cell<std::shared_ptr<IG_Cell_Group>> & 
        get_all_cell_groups();
        
        std::shared_ptr<Child_Mesh_Experimental>
        get_child_mesh(moris_index aChildMeshIndex);

        std::unordered_map<moris_id,moris_index>
        get_vertex_glb_id_to_loc_vertex_ind_map() const;


        mtk::Cell const &
        get_mtk_cell( moris_index aElementIndex ) const;

        mtk::Cell &
        get_mtk_cell( moris_index aElementIndex );

        mtk::Vertex & 
        get_mtk_vertex( moris_index aVertexIndex );

        mtk::Vertex const & 
        get_mtk_vertex( moris_index aVertexIndex ) const;


        moris::mtk::Vertex*
        get_mtk_vertex_pointer(moris_index aVertexIndex);

        std::shared_ptr<IG_Vertex_Group>
        get_vertex_group(moris_index aVertexGroupIndex);

        moris_index
        get_parent_cell_group_index( moris_index aParentCellIndex);
        
        void
        replace_controlled_ig_cell(
        moris_index aCellIndex,
        moris_id    aCellId,
        std::shared_ptr<moris::mtk::Cell_Info> aCellInfo,
        moris::Cell<moris::mtk::Vertex*> & aVertexPointers);

        

        void
        set_integration_cell(
            moris_index                          aCellIndex, 
            std::shared_ptr<xtk::Cell_XTK_No_CM> aNewCell );

        moris_index
        get_integration_cell_controlled_index(
            moris_index                       aCellIndex );

        void
        add_cell_to_integration_mesh(
            moris_index aCellIndex, 
            moris_index aCellGroupIndex );


        moris_id
        allocate_entity_ids( moris::size_t   aNumIdstoAllocate,
                             enum EntityRank aEntityRank);
        
        moris_id
        allocate_subphase_ids( moris::size_t   aNumIdstoAllocate );

        moris::moris_index
        get_first_available_index(enum EntityRank aEntityRank) const;

        moris::uint
        get_num_ig_cell_groups();

        std::shared_ptr<IG_Cell_Group>
        get_ig_cell_group(moris_index aGroupIndex);

        moris::Cell<moris_index> const &
        get_ig_cell_group_memberships(moris_index aIgCellIndex);

        moris::mtk::Cell*
        get_ig_cell_group_parent_cell(moris_index aGroupIndex);

        void
        set_child_mesh_subphase(
            moris_index aCMIndex,
            moris::Cell<moris_index> & aSubphasesGroups);

        moris::uint
        get_num_subphases();

        moris::Cell<std::shared_ptr<Child_Mesh_Experimental>> &
        get_owned_child_meshes();

        moris::Cell<moris_index> &
        get_owned_subphase_indices();

        moris::Cell<moris_index> &
        get_not_owned_subphase_indices();

        moris::mtk::Cell*
        get_subphase_parent_cell(moris_index aSubPhaseIndex);

        std::shared_ptr<IG_Cell_Group>
        get_subphase_ig_cells(moris_index aSubPhaseIndex);

        moris_index
        get_subphase_id(moris_index aSubPhaseIndex);

        moris_index
        get_subphase_index(moris_id aSubphaseId);

        moris_index
        get_subphase_bulk_phase(moris_index aSubPhaseIndex);

        moris::Cell<moris_index> const &
        get_parent_cell_subphases(moris_index aParentCellIndex);

        moris_index
        get_ig_cell_subphase_index(moris_index aIgCellIndex);

        bool
        parent_cell_has_children(moris_index aParentCellIndex);
        void
        finalize_cut_mesh_construction()
        {
            this->deduce_ig_cell_group_ownership();
            
            this->assign_controlled_ig_cell_ids();
        }

        void
        deduce_ig_cell_group_ownership();

        void
        assign_controlled_ig_cell_ids();

        void
        set_face_connectivity(std::shared_ptr<Facet_Based_Connectivity> aFaceConnectivity);

        std::shared_ptr<Facet_Based_Connectivity>
        get_face_connectivity();

        void
        set_interface_facets(moris::Cell<moris_index> & aInterfaces);

        moris::Cell<moris_index> const &
        get_interface_facets();

        void
        set_bulk_phase_to_bulk_phase_dbl_side_interface(moris::Cell<moris::Cell<std::shared_ptr<IG_Cell_Double_Side_Group>>> & aBptoBpDblSideInterfaces);

        moris::Cell<moris::Cell<std::shared_ptr<IG_Cell_Double_Side_Group>>> const &
        get_bulk_phase_to_bulk_phase_dbl_side_interface();

        void
        set_background_facet_to_child_facet_connectivity( moris::Cell<std::shared_ptr<moris::Cell<moris::moris_index>>> const & aBgtoChildFacet);


        moris::Cell<std::shared_ptr<moris::Cell<moris::moris_index>>> const &
        get_background_facet_to_child_facet_connectivity();

        void
        set_subphase_neighborhood(std::shared_ptr<Subphase_Neighborhood_Connectivity> aSubphaseNeighborhood);

        std::shared_ptr<Subphase_Neighborhood_Connectivity> 
        get_subphase_neighborhood();

        void
        setup_glob_to_loc_subphase_map();


        moris_index
        get_cell_bulk_phase( moris_index aCellIndex );

        // SET STUFF

        Cell<moris_index>
        register_side_set_names(moris::Cell<std::string> const & aSideSetNames)
        {
            uint tNumSetsToRegister = aSideSetNames.size();

            // block set ords
            Cell<moris_index> tSideSetOrds(tNumSetsToRegister);

            // iterate and add sets
            for(moris::uint i = 0; i < tNumSetsToRegister; i++)
            {
                tSideSetOrds(i) = mSideSetLabels.size();

                mSideSetLabels.push_back(aSideSetNames(i));
                mSideSetCellSides.push_back(nullptr);

                MORIS_ASSERT(mSideSideSetLabelToOrd.find(aSideSetNames(i)) ==  mSideSideSetLabelToOrd.end(),
                        "Duplicate side set in mesh");

                mSideSideSetLabelToOrd[aSideSetNames(i)] = tSideSetOrds(i) ;
            }
            return tSideSetOrds;
        }

        Cell<moris_index>
        register_block_set_names(moris::Cell<std::string> const & aBlockSetNames,
                                 enum CellTopology aCellTopo)
        {
            uint tNumSetsToRegister = aBlockSetNames.size();

            // block set ords
            Cell<moris_index> tBlockOrds(tNumSetsToRegister);

            // iterate and add sets
            for(moris::uint i = 0; i < tNumSetsToRegister; i++)
            {
                tBlockOrds(i) = mBlockSetNames.size();

                mBlockSetNames.push_back(aBlockSetNames(i));
                mBlockSetCellGroup.push_back(nullptr);
                mBlockCellTopo.push_back(aCellTopo);

                MORIS_ASSERT(mBlockSetLabelToOrd.find(aBlockSetNames(i)) ==  mBlockSetLabelToOrd.end(),
                        "Duplicate block set in mesh");

                mBlockSetLabelToOrd[aBlockSetNames(i)] = tBlockOrds(i) ;
            }

            return tBlockOrds;
        }

        void
        write_mesh(std::string aOutputPath,
                   std::string aOutputFile)
        {
            Tracer tTracer( "XTK", "Cut Integration Mesh", "Write mesh" );
            // get path to output XTK files to
            std::string tOutputPath = aOutputPath;
            std::string tOutputFile = aOutputFile;
            std::string tOutputBase = tOutputFile.substr(0,tOutputFile.find("."));
            std::string tOutputExt  = tOutputFile.substr(tOutputFile.find("."),tOutputFile.length());

            MORIS_ASSERT(tOutputExt == ".exo" || tOutputExt == ".e","Invalid file extension, needs to be .exo or .e");
            
            // Write mesh
            moris::mtk::Writer_Exodus writer( this );
 
            writer.write_mesh(
                "", tOutputPath + tOutputFile, 
                "", tOutputPath + "xtk_temp2.exo" );

            // Write the fields
            writer.set_time(0.0);
            writer.close_file();
        }




    private:

        void
        setup_comm_map();

        void
        prepare_child_element_identifier_requests(
            moris::Cell<moris::Cell<moris_id>>    & aNotOwnedChildMeshesToProcs,
            moris::Cell<moris::Matrix<IdMat>>     & aOwnedParentCellId,
            moris::Cell<moris::Matrix<IdMat>>     & aNumOwnedCellIdsOffsets,
            moris::Cell<uint>                     & aProcRanks,
            std::unordered_map<moris_id,moris_id> & aProcRankToDataIndex);
            
        void
        prepare_child_cell_id_answers(
            Cell<Matrix<IndexMat>> & aReceivedParentCellIds,
            Cell<Matrix<IndexMat>> & aReceivedParentCellNumChildren,
            Cell<Matrix<IndexMat>> & aChildCellIdOffset);
        void
        handle_received_child_cell_id_request_answers(
            Cell<Cell<moris_index>> const & aChildMeshesInInNotOwned,
            Cell<Matrix<IndexMat>>  const & aReceivedChildCellIdOffset);
    };
    

    
    class Child_Mesh_Experimental
    {
        friend class Cut_Integration_Mesh;
        friend class Integration_Mesh_Generator;
        public:
        moris::mtk::Cell*                mParentCell;
        moris::moris_index               mChildMeshIndex;
        std::shared_ptr<IG_Cell_Group>   mIgCells;
        std::shared_ptr<IG_Vertex_Group> mIgVerts;

        // subphases
        moris::Cell<std::shared_ptr<IG_Cell_Group>> mSubphaseCellGroups;

        public:
        Child_Mesh_Experimental()
        {

        }

        moris_index
        get_parent_element_index()
        {
            return this->get_parent_cell()->get_index();
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

        void
        set_subphase_groups( moris::Cell<std::shared_ptr<IG_Cell_Group>> & aSubphasesGroups)
        {
            mSubphaseCellGroups = aSubphasesGroups;
        }

        
    };

}

#endif