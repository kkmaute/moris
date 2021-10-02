#include "cl_XTK_Cut_Integration_Mesh.hpp"
#include "typedefs.hpp"
#include "cl_Cell.hpp"
#include "cl_MTK_Cell.hpp"
#include "cl_MTK_Vertex.hpp"

#include "cl_MTK_Mesh.hpp"
#include "cl_XTK_Model.hpp"
// #include "cl_XTK_Integration_Mesh_Generator.hpp"
#include "cl_Communication_Tools.hpp"
#include "cl_MPI_Tools.hpp"
#include "cl_XTK_Cell_No_CM.hpp"
using namespace moris;

namespace xtk
{
    IG_Cell_Group::IG_Cell_Group(moris_index aNumCellsInGroup):
    mIgCellGroup(aNumCellsInGroup,nullptr)
    {

    }

    IG_Cell_Group::IG_Cell_Group()
    {

    }

    IG_Cell_Side_Group::IG_Cell_Side_Group(moris_index aEstimatedNumCells)
    {
        mIgCells.reserve(aEstimatedNumCells);
        mIgCellSideOrdinals.reserve(aEstimatedNumCells);
    }

    IG_Cell_Double_Side_Group::IG_Cell_Double_Side_Group(moris_index aEstimatedNumCells)
    {
        mLeaderIgCells.reserve(aEstimatedNumCells);
        mLeaderIgCellSideOrdinals.reserve(aEstimatedNumCells);
        mFollowerIgCells.reserve(aEstimatedNumCells);
        mFollowerIgCellSideOrdinals.reserve(aEstimatedNumCells);
    }


    IG_Vertex_Group::IG_Vertex_Group(moris_index aNumVerticesInGroup):
    mIgVertexGroup(0,nullptr)
    {
        mIgVertexGroup.reserve(aNumVerticesInGroup);
    }

    std::size_t
    IG_Vertex_Group::size()
    {
        return mIgVertexGroup.size();
    }

    void
    IG_Vertex_Group::reserve(std::size_t aReserveSize)
    {
        mIgVertexGroup.reserve(aReserveSize);
        mIgVertexLocalCoords.reserve(aReserveSize);
    }

    void
    IG_Vertex_Group::add_vertex(
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
    IG_Vertex_Group::add_vertex_local_coord_pointers()
    {}

    moris::mtk::Vertex*
    IG_Vertex_Group::get_vertex(moris_index aGroupVertexOrdinal)
    {
        return mIgVertexGroup(aGroupVertexOrdinal);
    }

    moris_index 
    IG_Vertex_Group::get_vertex_group_ordinal(moris_index aVertex)
    {
        auto tIter = mIgVertexIndexToVertexOrdinal.find(aVertex);
        MORIS_ERROR(tIter != mIgVertexIndexToVertexOrdinal.end(),"Provided vertex not in group");
        return tIter->second;
    }

    std::shared_ptr<Matrix<DDRMat>>
    IG_Vertex_Group::get_vertex_local_coords(moris_index aVertex)
    {
        return mIgVertexLocalCoords(this->get_vertex_group_ordinal(aVertex));
    }

    void
    IG_Vertex_Group::print()
    {
        // iterate through vertices
        for(moris::uint iV = 0; iV < this->size(); iV++)
        {
            std::cout<<   "Vertex Id: "<<std::setw(8)<< mIgVertexGroup(iV)->get_id();
            std::cout<<" | Vertex Index: "<<std::setw(8)<<mIgVertexGroup(iV)->get_id();

            Matrix<DDRMat> tVertexCoords = mIgVertexGroup(iV)->get_coords();

            std::cout<<" | Coords:";
            for(moris::uint iSpatial = 0; iSpatial < tVertexCoords.numel(); iSpatial++ )
            {
                std::cout<<" "<<std::scientific<<tVertexCoords(iSpatial);
            }
            std::cout<<std::endl;
        }
    }

    // ----------------------------------------------------------------------------------
    Cut_Integration_Mesh::Cut_Integration_Mesh(
        moris::mtk::Mesh* aBackgroundMesh,
        xtk::Model* aXTKModel)
        :mBackgroundMesh(aBackgroundMesh),mXTKModel(aXTKModel)
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
        mParentCellCellGroupIndex.resize(tNumBackgroundCells,MORIS_INDEX_MAX);
        
        mIntegrationCellBulkPhase.resize(tNumBackgroundCells,MORIS_INDEX_MAX);

        for(moris::uint iCell = 0; iCell< tNumBackgroundCells; iCell++)
        {
            mIntegrationCells(iCell) = &mBackgroundMesh->get_mtk_cell((moris_index)iCell);

            // verify that we are not doubling up vertices in the id map
            MORIS_ERROR(mIntegrationCellIdToIndexMap.find(mIntegrationCells(iCell)->get_id()) == mIntegrationCellIdToIndexMap.end(), "Provided Cell Id is already in the integration vertex map: Vertex Id =%uon process %u", mIntegrationCells(iCell)->get_id(), par_rank());

            // add the vertex to id to index map
            mIntegrationCellIdToIndexMap[mIntegrationCells(iCell)->get_id()] = (moris_index) iCell;
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

    // ----------------------------------------------------------------------------------
    Cut_Integration_Mesh::~Cut_Integration_Mesh()
    {
    }
    // ----------------------------------------------------------------------------------
    moris::Cell<std::shared_ptr<Matrix<DDRMat>>>*
    Cut_Integration_Mesh::get_all_vertex_coordinates_loc_inds()
    {
        return &mVertexCoordinates;
    }

    moris::Cell<std::shared_ptr<IG_Cell_Group>> & 
    Cut_Integration_Mesh::get_all_cell_groups()
    {
        return mIntegrationCellGroups;
    }

    // ----------------------------------------------------------------------------------
    std::shared_ptr<Child_Mesh_Experimental>
    Cut_Integration_Mesh::get_child_mesh(moris_index aChildMeshIndex)
    {
        MORIS_ERROR(mChildMeshes.size() > (uint)aChildMeshIndex,"Child mesh index out of bounds");
        return mChildMeshes(aChildMeshIndex);
    }
    mtk::Vertex &
    Cut_Integration_Mesh::get_mtk_vertex( moris_index aVertexIndex )
    {
        MORIS_ASSERT(aVertexIndex < (moris_index)mIntegrationVertices.size(),
                "Vertex index out of bounds");

        return *mIntegrationVertices(aVertexIndex);
    }

    // ----------------------------------------------------------------------------

    mtk::Vertex const &
    Cut_Integration_Mesh::get_mtk_vertex( moris_index aVertexIndex ) const
    {
        MORIS_ASSERT(aVertexIndex < (moris_index)mIntegrationVertices.size(),
                "Vertex index out of bounds");

        return *mIntegrationVertices(aVertexIndex);
    }
    // ----------------------------------------------------------------------------------
    std::unordered_map<moris_id,moris_index>
    Cut_Integration_Mesh::get_vertex_glb_id_to_loc_vertex_ind_map() const
    {
        return mIntegrationVertexIdToIndexMap;
    }
    // ----------------------------------------------------------------------------------
    mtk::Cell const &
    Cut_Integration_Mesh::get_mtk_cell( moris_index aElementIndex ) const
    {
        return *mIntegrationCells(aElementIndex);
    }
    mtk::Cell &
    Cut_Integration_Mesh::get_mtk_cell( moris_index aElementIndex )
    {
        return *mIntegrationCells(aElementIndex);
    }
    // ----------------------------------------------------------------------------------
    moris::mtk::Vertex*
    Cut_Integration_Mesh::get_mtk_vertex_pointer(moris_index aVertexIndex)
    {
        return mIntegrationVertices(aVertexIndex);
    }
    // ----------------------------------------------------------------------------------
    std::shared_ptr<IG_Vertex_Group>
    Cut_Integration_Mesh::get_vertex_group(moris_index aVertexGroupIndex)
    {
        return mIntegrationVertexGroups(aVertexGroupIndex);
    }
    // ----------------------------------------------------------------------------------
    moris_index
    Cut_Integration_Mesh::get_parent_cell_group_index( moris_index aParentCellIndex)
    {
        return mParentCellCellGroupIndex(aParentCellIndex);
    }
    // ----------------------------------------------------------------------------------
    void
    Cut_Integration_Mesh::replace_controlled_ig_cell(
        moris_index aCellIndex,
        moris_id    aCellId,
        std::shared_ptr<moris::mtk::Cell_Info> aCellInfo,
        moris::Cell<moris::mtk::Vertex*> & aVertexPointers)
    {
        MORIS_ERROR(aCellIndex >= mFirstControlledCellIndex,"Cannot set integration cell that I do not control."); 

        // controlled index
        moris_index tIndexInControlledCells = aCellIndex - mFirstControlledCellIndex;

        mControlledIgCells(tIndexInControlledCells)->set_id(aCellId);
        mControlledIgCells(tIndexInControlledCells)->set_mtk_cell_info(aCellInfo);
        mControlledIgCells(tIndexInControlledCells)->set_vertex_pointers(aVertexPointers);

    }
    // ----------------------------------------------------------------------------------
    void
    Cut_Integration_Mesh::set_integration_cell(
            moris_index                          aCellIndex, 
            std::shared_ptr<xtk::Cell_XTK_No_CM> aNewCell )
    {
        MORIS_ERROR(aCellIndex >= mFirstControlledCellIndex,"Cannot set integration cell that I do not control."); 
        
        // controlled index
        moris_index tIndexInControlledCells = aCellIndex - mFirstControlledCellIndex;

        mControlledIgCells(tIndexInControlledCells) = aNewCell;
        mIntegrationCells(aCellIndex) = mControlledIgCells(tIndexInControlledCells).get();
    } 
       
    // ----------------------------------------------------------------------------------

    moris_index
    Cut_Integration_Mesh::get_integration_cell_controlled_index(
        moris_index aCellIndex )
    {
        MORIS_ERROR(aCellIndex >= mFirstControlledCellIndex,"Cell index provided I do not control"); 

        return aCellIndex - mFirstControlledCellIndex;
    }

    // ----------------------------------------------------------------------------------
    void
    Cut_Integration_Mesh::add_cell_to_integration_mesh(
            moris_index aCellIndex, 
            moris_index aCellGroupIndex )
    {
        MORIS_ERROR(aCellGroupIndex   < (moris_index)mIntegrationCellGroups.size(),"Child Mesh Index out of bounds."); 
        MORIS_ERROR(aCellIndex < (moris_index)mIntegrationCells.size(),"Cell Index out of bounds."); 
        mIntegrationCellGroups(aCellGroupIndex)->mIgCellGroup.push_back(mIntegrationCells(aCellIndex));
        mIntegrationCellToCellGroupIndex(aCellIndex).push_back(aCellGroupIndex);
    }
    // ----------------------------------------------------------------------------------
    moris_id
    Cut_Integration_Mesh::allocate_entity_ids( 
        moris::size_t   aNumIdstoAllocate,
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
    moris_id
    Cut_Integration_Mesh::allocate_subphase_ids( moris::size_t   aNumIdstoAllocate )
    {
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

        moris_index tFirstSubphaseId = mBackgroundMesh->get_max_entity_id(EntityRank::ELEMENT) + 1;


        if(tProcRank == 0)
        {
            // Loop over entities print the number of entities requested by each processor
            for (int iProc = 0; iProc < tProcSize; ++iProc)
            {

                // Give each processor their desired amount of IDs
                tProcFirstID(iProc) = tFirstSubphaseId;

                // Increment the first available node ID
                tFirstSubphaseId= tFirstSubphaseId+aGatheredInfo(iProc);

            }
        }

        moris::scatter(tProcFirstID,tFirstId);

        return tFirstId(0);
    }
    // ----------------------------------------------------------------------------------
    moris::moris_index
    Cut_Integration_Mesh::get_first_available_index(enum EntityRank aEntityRank) const
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
    moris::uint
    Cut_Integration_Mesh::get_num_ig_cell_groups()
    {
        return mIntegrationCellGroups.size();
    }
    // ----------------------------------------------------------------------------------
    std::shared_ptr<IG_Cell_Group>
    Cut_Integration_Mesh::get_ig_cell_group(moris_index aGroupIndex)
    {
        return mIntegrationCellGroups(aGroupIndex);
    }
    // ----------------------------------------------------------------------------------
    moris::Cell<moris_index> const &
    Cut_Integration_Mesh::get_ig_cell_group_memberships(moris_index aIgCellIndex)
    {
        return mIntegrationCellToCellGroupIndex(aIgCellIndex);
    }
    // ----------------------------------------------------------------------------------
    moris::mtk::Cell*
    Cut_Integration_Mesh::get_ig_cell_group_parent_cell(moris_index aGroupIndex)
    {
        return mIntegrationCellGroupsParentCell(aGroupIndex);
    }

    // ----------------------------------------------------------------------------------

    void
    Cut_Integration_Mesh::set_child_mesh_subphase(
            moris_index aCMIndex,
            moris::Cell<moris_index> & aSubphasesGroups)
    {
         moris::Cell<std::shared_ptr<IG_Cell_Group>> tIgCellSubphases(aSubphasesGroups.size());

         for(moris::uint i = 0; i < aSubphasesGroups.size(); i++)
         {
             tIgCellSubphases(i) = mSubPhaseCellGroups(aSubphasesGroups(i));
         }

        mChildMeshes(aCMIndex)->set_subphase_groups( tIgCellSubphases );
    }

    moris::uint
    Cut_Integration_Mesh::get_num_subphases()
    {
        return mSubPhaseCellGroups.size();
    }
    // ----------------------------------------------------------------------------------

    moris::Cell<std::shared_ptr<Child_Mesh_Experimental>> &
    Cut_Integration_Mesh::get_owned_child_meshes()
    {
        return mOwnedChildMeshes;
    }

    // ----------------------------------------------------------------------------------

    moris::Cell<moris_index> &
    Cut_Integration_Mesh::get_owned_subphase_indices()
    {
        return mOwnedSubphaseGroupsInds;
    }
   
    // ----------------------------------------------------------------------------------

    moris::Cell<moris_index> &
    Cut_Integration_Mesh::get_not_owned_subphase_indices()
    {
        return mNotOwnedSubphaseGroupsInds;
    }

    moris::mtk::Cell*
    Cut_Integration_Mesh::get_subphase_parent_cell(moris_index aSubPhaseIndex)
    {
        return mSubPhaseParentCell(aSubPhaseIndex);
    }

    std::shared_ptr<IG_Cell_Group>
    Cut_Integration_Mesh::get_subphase_ig_cells(moris_index aSubPhaseIndex)
    {
        return mSubPhaseCellGroups(aSubPhaseIndex);
    }

    moris_index
    Cut_Integration_Mesh::get_subphase_id(moris_index aSubPhaseIndex)
    {
        return mSubPhaseIds(aSubPhaseIndex);
    }

    moris_index
    Cut_Integration_Mesh::get_subphase_index(moris_id aSubphaseId)
    {
        auto tIter = mGlobalToLocalSubphaseMap.find(aSubphaseId);

        MORIS_ASSERT(tIter != mGlobalToLocalSubphaseMap.end(),"Subphase id not in map");

        return tIter->second;
    }

    moris_index
    Cut_Integration_Mesh::get_subphase_bulk_phase(moris_index aSubPhaseIndex)
    {
        return mSubPhaseBulkPhase(aSubPhaseIndex);
    }
    moris_index
    Cut_Integration_Mesh::get_ig_cell_subphase_index(moris_index aIgCellIndex)
    {
        return mIntegrationCellToSubphaseIndex(aIgCellIndex);
    }

    bool
    Cut_Integration_Mesh::parent_cell_has_children(moris_index aParentCellIndex)
    {
        return (bool)mParentCellHasChildren(aParentCellIndex);
    }
            
    moris::Cell<moris_index> const &
    Cut_Integration_Mesh::get_parent_cell_subphases(moris_index aParentCellIndex)
    {
        return mParentCellToSubphase(aParentCellIndex);
    }
    
    // ----------------------------------------------------------------------------------
    uint
    Cut_Integration_Mesh::get_num_entities( 
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

    std::shared_ptr<Facet_Based_Connectivity>
    Cut_Integration_Mesh::get_face_connectivity()
    {
        return mIgCellFaceConnectivity;
    }
    // ----------------------------------------------------------------------------------
    void
    Cut_Integration_Mesh::set_face_connectivity(std::shared_ptr<Facet_Based_Connectivity> aFaceConnectivity)
    {
        mIgCellFaceConnectivity = aFaceConnectivity;
    }

    void
    Cut_Integration_Mesh::set_interface_facets(moris::Cell<moris_index> & aInterfaces)
    {
        mInterfaceFacets = aInterfaces;
    }

    moris::Cell<moris_index> const &
    Cut_Integration_Mesh::get_interface_facets()
    {
        return mInterfaceFacets;
    }

    void
    Cut_Integration_Mesh::set_bulk_phase_to_bulk_phase_dbl_side_interface(moris::Cell<moris::Cell<std::shared_ptr<IG_Cell_Double_Side_Group>>> & aBptoBpDblSideInterfaces)
    {
        mBptoBpDblSideInterfaces = aBptoBpDblSideInterfaces;
    }

    moris::Cell<moris::Cell<std::shared_ptr<IG_Cell_Double_Side_Group>>> const &
    Cut_Integration_Mesh::get_bulk_phase_to_bulk_phase_dbl_side_interface()
    {
        return mBptoBpDblSideInterfaces;
    }


    void
    Cut_Integration_Mesh::set_background_facet_to_child_facet_connectivity( moris::Cell<std::shared_ptr<moris::Cell<moris::moris_index>>> const & aBgtoChildFacet)
    {
         mBGFacetToChildFacet = aBgtoChildFacet;
    }


    moris::Cell<std::shared_ptr<moris::Cell<moris::moris_index>>> const &
    Cut_Integration_Mesh::get_background_facet_to_child_facet_connectivity()
    {
        return mBGFacetToChildFacet;
    }

    void
    Cut_Integration_Mesh::set_subphase_neighborhood(std::shared_ptr<Subphase_Neighborhood_Connectivity> aSubphaseNeighborhood)
    {
        mSubphaseNeighborhood = aSubphaseNeighborhood;
    }

    std::shared_ptr<Subphase_Neighborhood_Connectivity> 
    Cut_Integration_Mesh::get_subphase_neighborhood()
    {
        return mSubphaseNeighborhood;
    }

    void
    Cut_Integration_Mesh::setup_glob_to_loc_subphase_map()
    {
        mGlobalToLocalSubphaseMap.clear();
        for(moris::uint i = 0; i < this->get_num_subphases(); i++)
        {
            moris_id tSubphaseId = this->get_subphase_id((moris_id)i);
            MORIS_ASSERT(mGlobalToLocalSubphaseMap.find(tSubphaseId) == mGlobalToLocalSubphaseMap.end(),"Subphase id already in map");
            mGlobalToLocalSubphaseMap[tSubphaseId] = i;
        }
    }

    moris_index
    Cut_Integration_Mesh::get_cell_bulk_phase( moris_index aCellIndex )
    {
        return mIntegrationCellBulkPhase(aCellIndex);
    }
    // ----------------------------------------------------------------------------------
    void
    Cut_Integration_Mesh::setup_comm_map()
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
        }
    }

    void
    Cut_Integration_Mesh::deduce_ig_cell_group_ownership()
    {
        mOwnedIntegrationCellGroupsInds.reserve(mIntegrationCellGroupsParentCell.size());
        mNotOwnedIntegrationCellGroups.reserve(mIntegrationCellGroupsParentCell.size());

        mOwnedChildMeshes.reserve(mChildMeshes.size());
        mNotOwnedChildMeshes.reserve(mChildMeshes.size());

        moris_index tParRank = moris::par_rank();
        // iterate through groups
        for(moris::uint i = 0; i < mIntegrationCellGroupsParentCell.size(); i++)
        {
            mIntegrationCellGroupsParentCell(i)->get_owner() == tParRank? 
                mOwnedIntegrationCellGroupsInds.push_back((moris_index) i): 
                mNotOwnedIntegrationCellGroups.push_back((moris_index) i);

            mIntegrationCellGroupsParentCell(i)->get_owner() == tParRank? 
                mOwnedChildMeshes.push_back(mChildMeshes(i)): 
                mNotOwnedChildMeshes.push_back(mChildMeshes(i));                
        }

        mOwnedIntegrationCellGroupsInds.shrink_to_fit();
        mNotOwnedIntegrationCellGroups.shrink_to_fit();
        mOwnedChildMeshes.shrink_to_fit();
        mNotOwnedChildMeshes.shrink_to_fit();
    }
    // ----------------------------------------------------------------------------------
    void
    Cut_Integration_Mesh::assign_controlled_ig_cell_ids()
    {
        
        // Set child element ids and indices
        moris::uint tNumControlledCellsInCutMesh = mControlledIgCells.size();

        // Allocate global element ids (these need to be give to the children meshes)
        moris_id    tElementIdOffset = this->allocate_entity_ids(tNumControlledCellsInCutMesh, moris::EntityRank::ELEMENT);

        // set child elements ids in the children meshes which I own and dont share
        for(moris::size_t i = 0; i<mOwnedIntegrationCellGroupsInds.size(); i++)
        {
            // tOwnedChildMeshes(i)->set_child_element_ids(tElementIdOffset);

            std::shared_ptr<IG_Cell_Group> tCellGroup = this->get_ig_cell_group(mOwnedIntegrationCellGroupsInds(i));

            // iterate through the child cell elements in the group
            for(moris::uint iCell = 0; iCell < tCellGroup->mIgCellGroup.size(); iCell++)
            {
                // get the controlled cell index
                moris_index tControlledIndex  = this->get_integration_cell_controlled_index(tCellGroup->mIgCellGroup(iCell)->get_index());

                mControlledIgCells(tControlledIndex)->set_id(tElementIdOffset++); 

                MORIS_ASSERT(mIntegrationCellIdToIndexMap.find(tCellGroup->mIgCellGroup(iCell)->get_id()) == mIntegrationCellIdToIndexMap.end(),"Id already in the map");

                mIntegrationCellIdToIndexMap[tCellGroup->mIgCellGroup(iCell)->get_id()] = tCellGroup->mIgCellGroup(iCell)->get_index();

            }   
        }

        // prepare outward requests
        Cell<Cell<moris_id>>        tNotOwnedChildMeshesToProcs;
        Cell<moris::Matrix<IdMat>>  tOwnedParentCellId;
        Cell<moris::Matrix<IdMat>>  tNumOwnedCellIdsOffsets;
        Cell<uint>                  tProcRanks;
        std::unordered_map<moris_id,moris_id>  tProcRankToDataIndex;

        this->prepare_child_element_identifier_requests(tNotOwnedChildMeshesToProcs, tOwnedParentCellId, tNumOwnedCellIdsOffsets, tProcRanks,tProcRankToDataIndex);

        // send requests
        moris::uint tMPITag = 141;
        mXTKModel->send_outward_requests(tMPITag, tProcRanks,tOwnedParentCellId);
        mXTKModel->send_outward_requests(tMPITag+1, tProcRanks, tNumOwnedCellIdsOffsets);

        barrier();

        // receive requests
        Cell<Matrix<IndexMat>> tReceivedParentCellIds;
        Cell<Matrix<IndexMat>> tReceivedParentCellNumChildren;
        Cell<uint> tProcsReceivedFrom1;
        Cell<uint> tProcsReceivedFrom2;
        mXTKModel->inward_receive_requests(tMPITag, 1, tReceivedParentCellIds, tProcsReceivedFrom1);
        mXTKModel->inward_receive_requests(tMPITag+1,1, tReceivedParentCellNumChildren, tProcsReceivedFrom2);

        MORIS_ASSERT(tProcsReceivedFrom1.size() == tProcsReceivedFrom2.size(),
                "Size mismatch between procs received from child cell ids and number of child cells");

        Cell<Matrix<IndexMat>> tChildIdOffsets;
        this->prepare_child_cell_id_answers(tReceivedParentCellIds,tReceivedParentCellNumChildren,tChildIdOffsets);

        // return information
        mXTKModel->return_request_answers(tMPITag+2, tChildIdOffsets, tProcsReceivedFrom1);

        // receive the information
        barrier();

        // receive the answers
        Cell<Matrix<IndexMat>> tReceivedChildIdOffsets;
        mXTKModel->inward_receive_request_answers(tMPITag+2,1,tProcRanks,tReceivedChildIdOffsets);

        // add child cell ids to not owned child meshes
        this->handle_received_child_cell_id_request_answers(tNotOwnedChildMeshesToProcs,tReceivedChildIdOffsets);

        barrier();
    }
    // ----------------------------------------------------------------------------------
    void
    Cut_Integration_Mesh::prepare_child_element_identifier_requests(
            Cell<Cell<moris_id>>       & aNotOwnedChildMeshesToProcs,
            Cell<moris::Matrix<IdMat>> & aOwnedParentCellId,
            Cell<moris::Matrix<IdMat>> & aNumOwnedCellIdsOffsets,
            Cell<uint>                 & aProcRanks,
            std::unordered_map<moris_id,moris_id> & aProcRankToDataIndex)
    {
        
        Cell<moris_id>            tCounts(0);
        moris_index tCurrentIndex = 0;

        // access the communication table
        Matrix<IdMat> tCommTable = this->get_communication_table();

        for(moris::size_t i = 0; i < tCommTable.numel(); i++)
        {
            aProcRankToDataIndex[tCommTable(i)] = tCurrentIndex;
            aProcRanks.push_back(tCommTable(i));
            aNotOwnedChildMeshesToProcs.push_back(Cell<moris_id>(0));
            tCounts.push_back(0);
            tCurrentIndex++;
        }

        for(moris::size_t i = 0; i < mNotOwnedIntegrationCellGroups.size(); i++)
        {
            moris_index tOwnerProc = mIntegrationCellGroupsParentCell(mNotOwnedIntegrationCellGroups(i))->get_owner();
            moris_index tProcDataIndex = aProcRankToDataIndex[tOwnerProc];
            aNotOwnedChildMeshesToProcs(tProcDataIndex).push_back(mNotOwnedIntegrationCellGroups(i));
        }

        aOwnedParentCellId.resize(aNotOwnedChildMeshesToProcs.size());
        aNumOwnedCellIdsOffsets.resize(aNotOwnedChildMeshesToProcs.size());

        // iterate through procs and child meshes shared with that processor
        for(moris::size_t i = 0; i < aNotOwnedChildMeshesToProcs.size(); i++)
        {
            // number of child meshes shared with this processor
            moris::uint tNumCM = aNotOwnedChildMeshesToProcs(i).size();

            // allocate matrix
            aOwnedParentCellId(i).resize(1,tNumCM);
            aNumOwnedCellIdsOffsets(i).resize(1,tNumCM);

            for(moris::uint j = 0; j < tNumCM; j++)
            {
                moris_index tIGCellGroupIndex = aNotOwnedChildMeshesToProcs(i)(j);
                aOwnedParentCellId(i)(j)      = mIntegrationCellGroupsParentCell(tIGCellGroupIndex)->get_id();
                aNumOwnedCellIdsOffsets(i)(j) = mIntegrationCellGroups(tIGCellGroupIndex)->mIgCellGroup.size();
            }
        }

        for(moris::size_t i = 0; i < tCommTable.numel(); i++)
        {
            if(aNotOwnedChildMeshesToProcs(i).size() == 0)
            {
                aOwnedParentCellId(i).resize(1,1);
                aNumOwnedCellIdsOffsets(i).resize(1,1);
                aOwnedParentCellId(i)(0,0) = MORIS_INDEX_MAX;
                aNumOwnedCellIdsOffsets(i)(0,0) = MORIS_INDEX_MAX;
            }
        }
    }
    // ----------------------------------------------------------------------------------
    void
    Cut_Integration_Mesh::prepare_child_cell_id_answers(
            Cell<Matrix<IndexMat>> & aReceivedParentCellIds,
            Cell<Matrix<IndexMat>> & aReceivedParentCellNumChildren,
            Cell<Matrix<IndexMat>> & aChildCellIdOffset)
    {
        MORIS_ASSERT(aReceivedParentCellIds.size() == aReceivedParentCellNumChildren.size(),
                "Mismatch in received parent cell ids and received parent cell number of children");

        // allocate answer size
        aChildCellIdOffset.resize(aReceivedParentCellIds.size());

        // iterate through received data
        for(moris::uint i = 0; i < aReceivedParentCellIds.size(); i++)
        {
            uint tNumReceivedReqs = aReceivedParentCellIds(i).n_cols();

            aChildCellIdOffset(i).resize(1,tNumReceivedReqs);

            if(aReceivedParentCellIds(i)(0) != MORIS_INDEX_MAX)
            {
                // iterate through received requests
                for(moris::uint j = 0; j < tNumReceivedReqs; j++)
                {
                    // parent cell information
                    moris_id tParentId           = aReceivedParentCellIds(i)(0,j);
                    moris_index tParentCellIndex = mBackgroundMesh->get_loc_entity_ind_from_entity_glb_id(tParentId,EntityRank::ELEMENT);

                    moris_index tCMIndex = mParentCellCellGroupIndex(tParentCellIndex);

                    // get child mesh
                    MORIS_ASSERT(tCMIndex != MORIS_INDEX_MAX, "Request is made for child element ids on a parent cell not intersected");

                    MORIS_ASSERT(par_rank() == mIntegrationCellGroupsParentCell(tCMIndex)->get_owner(),"I dont own this entity that had info requested.");
    
                    // place in return data
                    MORIS_ASSERT(mIntegrationCellGroups(tCMIndex)->mIgCellGroup.size() == (uint)aReceivedParentCellNumChildren(i)(j),
                            "Number of child cells in child meshes do not match number on other processor");


                    aChildCellIdOffset(i)(j) = mIntegrationCellGroups(tCMIndex)->mIgCellGroup(0)->get_id();
                }
            }
            else
            {
                aChildCellIdOffset(i)(0) =   MORIS_INDEX_MAX;
            }
        }
    }
    // ----------------------------------------------------------------------------------
    void
    Cut_Integration_Mesh::handle_received_child_cell_id_request_answers(
            Cell<Cell<moris_index>> const & aChildMeshesInInNotOwned,
            Cell<Matrix<IndexMat>>  const & aReceivedChildCellIdOffset)
    {
        // iterate through received data
        for(moris::uint i = 0; i < aChildMeshesInInNotOwned.size(); i++)
        {
            uint tNumReceivedReqs = aChildMeshesInInNotOwned(i).size();

            // iterate through received requests
            for(moris::uint j = 0; j < tNumReceivedReqs; j++)
            {
                moris_id tChildMeshInNotOwned = aChildMeshesInInNotOwned(i)(j);
                moris_id tChildCellFirstId = aReceivedChildCellIdOffset(i)(j);

                std::shared_ptr<IG_Cell_Group> tCellGroup = this->get_ig_cell_group(tChildMeshInNotOwned);

                for(moris::uint iCell = 0; iCell < tCellGroup->mIgCellGroup.size(); iCell++)
                {
                    // get the controlled cell index
                    moris_index tControlledIndex  = this->get_integration_cell_controlled_index(tCellGroup->mIgCellGroup(iCell)->get_index());

                    mControlledIgCells(tControlledIndex)->set_id(tChildCellFirstId++); 

                    MORIS_ASSERT(mIntegrationCellIdToIndexMap.find(tCellGroup->mIgCellGroup(iCell)->get_id()) == mIntegrationCellIdToIndexMap.end(),"Id already in the map");

                    mIntegrationCellIdToIndexMap[tCellGroup->mIgCellGroup(iCell)->get_id()] = tCellGroup->mIgCellGroup(iCell)->get_index();
                }   

            }
        }

    }
    // ----------------------------------------------------------------------------------
    // ----------------------------------------------------------------------------------
    // ----------------------------------------------------------------------------------
    // ----------------------------------------------------------------------------------



}
