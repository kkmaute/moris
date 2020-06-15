/*
 * cl_XTK_Enriched_Interpolation_Mesh.cpp
 *
 *  Created on: Jul 10, 2019
 *      Author: doble
 */



#include "cl_XTK_Enriched_Interpolation_Mesh.hpp"
#include "cl_Map.hpp"

#include "cl_XTK_Multigrid.hpp"

namespace xtk
{
    //------------------------------------------------------------------------------
    Enriched_Interpolation_Mesh::Enriched_Interpolation_Mesh(
            Model* aXTKModel)
    : mXTKModel(aXTKModel),
      mNumVerts(0),
      mNumVertsPerInterpCell(MORIS_UINT_MAX),
      mCellInfo(nullptr)
    {

    }

    Enriched_Interpolation_Mesh::~Enriched_Interpolation_Mesh()
    {
        if( mCellInfo!= nullptr ) { delete mCellInfo; }

        for(moris::uint i = 0; i < mInterpVertEnrichment.size(); i++)
        {
            for(auto it : mInterpVertEnrichment(i))
            {
                delete it;
            }

            mInterpVertEnrichment(i).clear();
        }

        mInterpVertEnrichment.clear();

        for(auto it : mEnrichedInterpCells)
        {
            delete it;
        }

        mEnrichedInterpCells.clear();
    }

    //------------------------------------------------------------------------------
    MeshType
    Enriched_Interpolation_Mesh::get_mesh_type() const
    {
        return MeshType::XTK;
    }
    //------------------------------------------------------------------------------
    moris::uint
    Enriched_Interpolation_Mesh::get_spatial_dim() const
    {
        return mXTKModel->get_spatial_dim();
    }
    //------------------------------------------------------------------------------
    uint
    Enriched_Interpolation_Mesh::get_num_entities(
            enum EntityRank aEntityRank,
            const moris_index aIndex ) const
    {
        switch(aEntityRank)
        {
            case(EntityRank::NODE):
                        {
                return mEnrichedInterpVerts.size();
                break;
                        }
            case(EntityRank::ELEMENT):
                        {
                return mEnrichedInterpCells.size();
                break;
                        }
            default:
                MORIS_ERROR(0,"Only support get num entities for nodes and elements currently");
                return 0;
        }
    }
    //------------------------------------------------------------------------------
    uint
    Enriched_Interpolation_Mesh::get_num_coeffs(const uint aBSplineMeshIndex) const
    {
        moris_index tLocalMeshIndex = this->get_local_mesh_index(aBSplineMeshIndex);

        return mEnrichCoeffLocToGlob(tLocalMeshIndex).numel();
    }

    //------------------------------------------------------------------------------
    Matrix<IndexMat>
    Enriched_Interpolation_Mesh::get_entity_connected_to_entity_loc_inds(
            moris_index       aEntityIndex,
            enum EntityRank   aInputEntityRank,
            enum EntityRank   aOutputEntityRank,
            const moris_index aIndex ) const
    {
        MORIS_ERROR(aInputEntityRank == EntityRank::ELEMENT && aOutputEntityRank == EntityRank::NODE,"Only support element to node connectivity");
        MORIS_ASSERT(aEntityIndex<(moris_index)mEnrichedInterpCells.size(),"Element index out of bounds");
        return mEnrichedInterpCells(aEntityIndex)->get_vertex_inds();
    }
    //------------------------------------------------------------------------------
    Matrix< IndexMat >
    Enriched_Interpolation_Mesh::get_elements_connected_to_element_and_face_ind_loc_inds(moris_index aElementIndex) const
    {
        MORIS_ERROR(0,"XTK ENRICHED MESH ERROR: get_elements_connected_to_element_and_face_ind_loc_inds no implemented");
        return Matrix<IndexMat>(0,0);
    }
    //------------------------------------------------------------------------------
    Cell<mtk::Vertex const *>
    Enriched_Interpolation_Mesh::get_all_vertices() const
    {
        moris::uint tNumNodes = this->get_num_entities(EntityRank::NODE);
        Cell<mtk::Vertex const *> tVertices(tNumNodes);
        for(moris::uint i = 0; i < tNumNodes; i++)
        {
            tVertices(i) = & mEnrichedInterpVerts(i);
        }
        return tVertices;
    }
    //------------------------------------------------------------------------------
    moris_id
    Enriched_Interpolation_Mesh::get_glb_entity_id_from_entity_loc_index(
            moris_index       aEntityIndex,
            enum EntityRank   aEntityRank,
            const moris_index aIndex ) const
    {
        moris::uint tMapIndex = (uint)aEntityRank;

        MORIS_ASSERT(aEntityRank == EntityRank::NODE || aEntityRank == EntityRank::ELEMENT,"XTK ENRICHED MESH ERROR: Only support glb to loc conversion for vertices and cells");
        MORIS_ASSERT(aEntityIndex<(moris_index)mLocalToGlobalMaps(tMapIndex).numel(),"Entityindex out of bounds");
        return mLocalToGlobalMaps(tMapIndex)(aEntityIndex);
    }
    //------------------------------------------------------------------------------
    moris_index
    Enriched_Interpolation_Mesh::get_loc_entity_ind_from_entity_glb_id(
            moris_id          aEntityId,
            enum EntityRank   aEntityRank,
            const moris_index aIndex ) const
    {
        moris::uint tMapIndex = (uint)aEntityRank;

        auto tIter = mGlobaltoLobalMaps(tMapIndex).find(aEntityId);

        MORIS_ASSERT(aEntityRank == EntityRank::NODE || aEntityRank == EntityRank::ELEMENT,"XTK ENRICHED MESH ERROR: Only support glb to loc conversion for vertices and cells");

        if(tIter ==  mGlobaltoLobalMaps(tMapIndex).end())
        {
            std::cout<<"Not Found  Entity Id = "<<aEntityId<<" | par_rank = "<<par_rank()<<std::endl;
        }

        MORIS_ASSERT(tIter!=  mGlobaltoLobalMaps(tMapIndex).end(),"Id does not appear in map");
        return tIter->second;
    }

    //------------------------------------------------------------------------------

    Matrix<IdMat>
    Enriched_Interpolation_Mesh::get_entity_connected_to_entity_glob_ids(
            moris_id          aEntityId,
            enum EntityRank   aInputEntityRank,
            enum EntityRank   aOutputEntityRank,
            const moris_index aIndex ) const
    {
        moris_index tEntityIndex = get_loc_entity_ind_from_entity_glb_id(aEntityId,aInputEntityRank);

        Matrix<IndexMat> tEntityToEntityLoc = this->get_entity_connected_to_entity_loc_inds(tEntityIndex,aInputEntityRank,aOutputEntityRank);

        return convert_indices_to_ids(tEntityToEntityLoc,aOutputEntityRank);
    }
    //------------------------------------------------------------------------------
    Matrix< DDRMat >
    Enriched_Interpolation_Mesh::get_node_coordinate( moris_index aNodeIndex ) const
    {
        mtk::Vertex const & tVertex = get_mtk_vertex(aNodeIndex);
        return tVertex.get_coords();
    }
    //------------------------------------------------------------------------------
    mtk::Vertex &
    Enriched_Interpolation_Mesh::get_mtk_vertex( moris_index aVertexIndex )
    {
        MORIS_ASSERT(aVertexIndex < (moris_index)mEnrichedInterpVerts.size(),"Vertex index out of bounds");
        return mEnrichedInterpVerts(aVertexIndex);
    }
    mtk::Vertex const &
    Enriched_Interpolation_Mesh::get_mtk_vertex( moris_index aVertexIndex ) const
    {
        MORIS_ASSERT(aVertexIndex < (moris_index)mEnrichedInterpVerts.size(),"Vertex index out of bounds");
        return mEnrichedInterpVerts(aVertexIndex);
    }
    //------------------------------------------------------------------------------
    mtk::Cell const &
    Enriched_Interpolation_Mesh::get_mtk_cell( moris_index aElementIndex ) const
    {
        MORIS_ASSERT(aElementIndex < (moris_index)mEnrichedInterpCells.size(),"Cell index out of bounds");
        return *mEnrichedInterpCells(aElementIndex);
    }
    //------------------------------------------------------------------------------
    mtk::Cell &
    Enriched_Interpolation_Mesh::get_writable_mtk_cell( moris_index aElementIndex )
    {
        MORIS_ASSERT(aElementIndex < (moris_index)mEnrichedInterpCells.size(),"Cell index out of bounds");
        return *mEnrichedInterpCells(aElementIndex);
    }
    //------------------------------------------------------------------------------
    Matrix< IdMat >
    Enriched_Interpolation_Mesh::get_communication_table() const
    {
        return mXTKModel->get_background_mesh().get_communication_table();
    }
    //------------------------------------------------------------------------------
    uint
    Enriched_Interpolation_Mesh::get_num_elements()
    {
        return mEnrichedInterpCells.size();
    }
    //------------------------------------------------------------------------------
    moris_id
    Enriched_Interpolation_Mesh::get_max_entity_id( enum EntityRank aEntityRank,const moris_index aIndex ) const
    {
        MORIS_ASSERT(aEntityRank == EntityRank::NODE || aEntityRank == EntityRank::ELEMENT,"Only Elements or Nodes have max id");

        moris::uint tNumEntities = this->get_num_entities(aEntityRank);

        moris_id tLocalMaxId = 0;

        for(moris::uint i = 0; i < tNumEntities; i++)
        {
            moris_id tId = this->get_glb_entity_id_from_entity_loc_index(i,aEntityRank);

            if(tId > tLocalMaxId)
            {
                tLocalMaxId = tId;
            }
        }

        moris_id tGlobalMaxId = 0;
        moris::max_all(tLocalMaxId,tGlobalMaxId);
        return tGlobalMaxId;
    }
    //------------------------------------------------------------------------------
    void
    Enriched_Interpolation_Mesh::get_adof_map( const uint aBSplineIndex, map< moris_id, moris_index > & aAdofMap ) const
    {
        aAdofMap.clear();

        moris_index tLocalMeshIndex = this->get_local_mesh_index(aBSplineIndex);

        for(moris::uint iB =0; iB <mEnrichCoeffLocToGlob(tLocalMeshIndex).numel(); iB++)
        {
            MORIS_ASSERT(!aAdofMap.key_exists(mEnrichCoeffLocToGlob(tLocalMeshIndex)(iB)),"Duplicate id in the basis map detected");

            aAdofMap[mEnrichCoeffLocToGlob(tLocalMeshIndex)(iB)] = (moris_index) iB;

        }


    }
    //------------------------------------------------------------------------------
    Matrix<IndexMat> const &
    Enriched_Interpolation_Mesh::get_enriched_coefficients_at_background_coefficient(moris_index const & aMeshIndex,
            moris_index aBackgroundCoeffIndex) const
    {
        moris_index tLocalMeshIndex = this->get_local_mesh_index(aMeshIndex);
        MORIS_ASSERT(aBackgroundCoeffIndex< (moris_index)mCoeffToEnrichCoeffs(tLocalMeshIndex).size(), "Background coefficient index out of bounds. Be sure this is not an enriched coefficient index passed in.");
        return mCoeffToEnrichCoeffs(tLocalMeshIndex)(aBackgroundCoeffIndex);
    }
    Cell<Matrix<IndexMat>> const &
    Enriched_Interpolation_Mesh::get_enriched_coefficients_to_background_coefficients(moris_index const & aMeshIndex) const
    {
        moris_index tLocalMeshIndex = this->get_local_mesh_index(aMeshIndex);
        return mCoeffToEnrichCoeffs(tLocalMeshIndex);
    }
    //------------------------------------------------------------------------------
    Matrix<IndexMat> const &
    Enriched_Interpolation_Mesh::get_enriched_coefficient_local_to_global_map(moris_index const & aMeshIndex) const
    {
        moris_index tLocalMeshIndex = this->get_local_mesh_index(aMeshIndex);
        return mEnrichCoeffLocToGlob(tLocalMeshIndex);
    }
    //------------------------------------------------------------------------------
    Matrix<IndexMat>
    Enriched_Interpolation_Mesh::get_background_coefficient_local_to_global_map() const
    {
        moris::mtk::Mesh & tBackgroundMeshData = mXTKModel->get_background_mesh().get_mesh_data();

        moris::uint tNumBackgroundCoeffs = tBackgroundMeshData.get_num_entities(mBasisRank);

        Matrix<IndexMat> tCoefficientLocalToGlobal(1,tNumBackgroundCoeffs);

        for(moris::uint i = 0; i < tNumBackgroundCoeffs; i++)
        {
            tCoefficientLocalToGlobal(i) = tBackgroundMeshData.get_glb_entity_id_from_entity_loc_index((moris_index) i, mBasisRank);
        }

        return tCoefficientLocalToGlobal;
    }
    //------------------------------------------------------------------------------
    uint
    Enriched_Interpolation_Mesh::get_num_background_coefficients() const
    {
        return mCoeffToEnrichCoeffs.size();
    }
    //------------------------------------------------------------------------------
    uint
    Enriched_Interpolation_Mesh::get_num_verts_per_interp_cell()
    {
        MORIS_ASSERT(mNumVertsPerInterpCell != MORIS_UINT_MAX,"Number of verts per interpolation cell not set");
        return mNumVertsPerInterpCell;
    }
    //------------------------------------------------------------------------------
    Interpolation_Vertex_Unzipped*
    Enriched_Interpolation_Mesh::get_unzipped_vertex_pointer(moris_index aVertexIndex)
    {
        MORIS_ASSERT(aVertexIndex < (moris_index)mEnrichedInterpVerts.size(),"Provided vertex index is out of bounds");
        return &mEnrichedInterpVerts(aVertexIndex);
    }

    //------------------------------------------------------------------------------
    moris_index
    Enriched_Interpolation_Mesh::add_vertex_enrichment(
            moris_index   const & aMeshIndex,
            mtk::Vertex *         aBaseInterpVertex,
            Vertex_Enrichment   & aVertexEnrichment,
            bool                & aNewVertex)
    {
        moris_index tLocalMeshIndex = this->get_local_mesh_index(aMeshIndex);

        // vertex index of the base interpolation vertex
        moris_index tBaseVertIndex = aBaseInterpVertex->get_index();

        // Number of enriched vertices related to the base vertex
        moris::uint tNumVertsEnrOnBaseVert = mBaseInterpVertToVertEnrichmentIndex(tLocalMeshIndex)(tBaseVertIndex).size();

        // not new until we make it to the end
        aNewVertex = false;

        // iterate through the enriched vertices related to the base vertex and see if any are equal
        for(moris::uint i = 0; i < tNumVertsEnrOnBaseVert; i++)
        {
            moris_index tVertEnrIndex = mBaseInterpVertToVertEnrichmentIndex(tLocalMeshIndex)(tBaseVertIndex)(i);

            if(aVertexEnrichment == *mInterpVertEnrichment(tLocalMeshIndex)(tVertEnrIndex))
            {
                return tVertEnrIndex;
            }
        }

        // if we make it through the loop without finding an enrichment vertex
        // make a new one
        aNewVertex = true;

        // index of the vertex enrichment
        moris_index tVertexEnrichmentIndex = mInterpVertEnrichment(tLocalMeshIndex).size();

        // add to member data
        mInterpVertEnrichment(tLocalMeshIndex).push_back(new Vertex_Enrichment(aVertexEnrichment));

        // add a dummy value to the parent vertex index of a vertex interpolation
        mVertexEnrichmentParentVertexIndex(tLocalMeshIndex).push_back(tVertexEnrichmentIndex);

        mBaseInterpVertToVertEnrichmentIndex(tLocalMeshIndex)(tBaseVertIndex).push_back(tVertexEnrichmentIndex);

        return tVertexEnrichmentIndex;
    }
    //------------------------------------------------------------------------------
    Vertex_Enrichment *
    Enriched_Interpolation_Mesh::get_vertex_enrichment(
            moris_index const & aMeshIndex,
            moris_index const & aVertexEnrichmentIndex)
    {
        moris_index tLocalMeshIndex = this->get_local_mesh_index(aMeshIndex);
        MORIS_ASSERT(aVertexEnrichmentIndex< (moris_index)mInterpVertEnrichment(tLocalMeshIndex).size(),"Provided vertex enrichment index out of bounds");
        return mInterpVertEnrichment(tLocalMeshIndex)(aVertexEnrichmentIndex);
    }
    //------------------------------------------------------------------------------
    moris_index
    Enriched_Interpolation_Mesh::get_vertex_related_to_vertex_enrichment(
            moris_index const & aMeshIndex,
            moris_index aVertexEnrichmentIndex) const
    {
        moris_index tLocalMeshIndex = this->get_local_mesh_index(aMeshIndex);
        MORIS_ASSERT(aVertexEnrichmentIndex< (moris_index)mVertexEnrichmentParentVertexIndex(tLocalMeshIndex).size(),"Provided vertex enrichment index out of bounds");
        return mVertexEnrichmentParentVertexIndex(tLocalMeshIndex)(aVertexEnrichmentIndex);
    }
    //------------------------------------------------------------------------------
    moris_index
    Enriched_Interpolation_Mesh::get_local_mesh_index(moris_index const & aMeshIndex) const
    {
        auto tIter = mMeshIndexToLocMeshIndex.find(aMeshIndex);
        MORIS_ASSERT(tIter != mMeshIndexToLocMeshIndex.end(),"Mesh index not in map");
        return tIter->second;
    }
    //------------------------------------------------------------------------------
    Cell<moris_index> const &
    Enriched_Interpolation_Mesh::get_not_owned_vertex_indices() const
    {
        return mNotOwnedVerts;
    }
    bool
    Enriched_Interpolation_Mesh::basis_exists_on_partition(
            moris_index const & aMeshIndex,
            moris_index const & aBasisId)
    {
        moris_index tMeshIndex = this->get_local_mesh_index(aMeshIndex);
        if(mGlobaltoLocalBasisMaps(tMeshIndex).find(aBasisId) == mGlobaltoLocalBasisMaps(tMeshIndex).end())
        {
            return false;
        }

        return true;
    }
    //------------------------------------------------------------------------------
    moris_index
    Enriched_Interpolation_Mesh::add_basis_function(
            moris_index const & aMeshIndex,
            moris_index const & aBasisIdToAdd,
            moris_index const & aBasisOwner)
    {
        MORIS_ASSERT(!this->basis_exists_on_partition(aMeshIndex,aBasisIdToAdd),"Basis that you are trying to add already exists in this mesh");

        moris_index tLocMesh = this->get_local_mesh_index(aMeshIndex);
        moris_index tNewIndex = mEnrichCoeffLocToGlob(tLocMesh).numel();

        // add a size of 1
        mEnrichCoeffLocToGlob(tLocMesh).resize(1,tNewIndex+1);
        mEnrichCoeffOwnership(tLocMesh).resize(1,tNewIndex+1);


        // add the local to glb map
        mEnrichCoeffLocToGlob(tLocMesh)(tNewIndex) = aBasisIdToAdd;
        mEnrichCoeffOwnership(tLocMesh)(tNewIndex) = aBasisOwner;

        // add to glb to local map
        mGlobaltoLocalBasisMaps(tLocMesh)[aBasisIdToAdd] = tNewIndex;

        return tNewIndex;
    }
    //------------------------------------------------------------------------------

    moris::Cell<Interpolation_Cell_Unzipped const*>
    Enriched_Interpolation_Mesh::get_enriched_cells_from_base_cells(moris::Cell<moris::mtk::Cell const*> const & aBaseCells) const
    {
        uint tNumBaseCells = aBaseCells.size();

        moris::Cell<Interpolation_Cell_Unzipped const*> tEnrichedCells;
        for(moris::uint  i = 0 ; i <tNumBaseCells; i++)
        {
            tEnrichedCells.append(this->get_enriched_cells_from_base_cell(aBaseCells(i)));
        }
        return tEnrichedCells;

    }
    //------------------------------------------------------------------------------
    Cell<Interpolation_Cell_Unzipped *> const &
    Enriched_Interpolation_Mesh::get_enriched_interpolation_cells() const
    {
        return mEnrichedInterpCells;
    }
    //------------------------------------------------------------------------------
    uint
    Enriched_Interpolation_Mesh::get_num_interpolation_types() const
    {
        return mMeshIndices.numel();
    }
    //------------------------------------------------------------------------------
    moris_index
    Enriched_Interpolation_Mesh::get_interpolation_index(moris_index const & aLocalInterpIndex) const
    {
        MORIS_ASSERT(aLocalInterpIndex < (moris_index) mMeshIndices.numel(),"Local interpolation index out of bounds");
        return mMeshIndices(aLocalInterpIndex);
    }
    //------------------------------------------------------------------------------
    moris_index
    Enriched_Interpolation_Mesh::get_basis_owner(moris_index aBasisIndex,
                                                 moris_index aMeshIndex)
    {
        moris_index tLocMesh = this->get_local_mesh_index(aMeshIndex);

        return mEnrichCoeffOwnership(tLocMesh)(aBasisIndex);
    }

    //------------------------------------------------------------------------------
    Cell<Interpolation_Cell_Unzipped *> &
    Enriched_Interpolation_Mesh::get_enriched_interpolation_cells()
    {
        return mEnrichedInterpCells;
    }
    //------------------------------------------------------------------------------
    void
    Enriched_Interpolation_Mesh::get_owned_and_not_owned_enriched_interpolation_cells(
            Cell<Interpolation_Cell_Unzipped *>       & aOwnedInterpCells,
            Cell<Cell<Interpolation_Cell_Unzipped *>> & aNotOwnedInterpCells,
            Cell<uint>                                & aProcRanks)
    {
        // get all interp cells
        Cell<Interpolation_Cell_Unzipped*> & tEnrInterpCells = this->get_enriched_interpolation_cells();

        // reserve space
        aOwnedInterpCells.resize(0);
        aOwnedInterpCells.reserve(tEnrInterpCells.size());
        aNotOwnedInterpCells.resize(0);
        // proc rank
        moris_index tParRank = par_rank();

        // counter
        moris::uint       tOwnerCount = 0;
        Cell<moris::uint> tCounts(0);

        // map
        std::unordered_map<moris_id,moris_id> tProcRankToDataIndex;

        // access the communication table
        Matrix<IdMat> tCommTable = this->get_communication_table();

        // resize proc ranks and setup map to comm table
        aProcRanks.resize(tCommTable.numel());
        for(moris::uint i = 0; i <tCommTable.numel(); i++)
        {
            tProcRankToDataIndex[tCommTable(i)] = i;
            aProcRanks(i) = (tCommTable(i));
            aNotOwnedInterpCells.push_back(Cell<Interpolation_Cell_Unzipped *>(0));
        }

        for(moris::uint i = 0; i <tEnrInterpCells.size(); i++)
        {
            moris_index tOwnerProc = tEnrInterpCells(i)->get_owner();

            if(tParRank == tOwnerProc)
            {
                aOwnedInterpCells.push_back( tEnrInterpCells(i) );
                tOwnerCount++;
            }
            else
            {
                moris_index tProcDataIndex = tProcRankToDataIndex[tOwnerProc];
                aNotOwnedInterpCells(tProcDataIndex).push_back( tEnrInterpCells(i) );
            }
        }
    }
    //------------------------------------------------------------------------------
    Interpolation_Vertex_Unzipped &
    Enriched_Interpolation_Mesh::get_xtk_interp_vertex(moris::uint aVertexIndex)
    {
        return mEnrichedInterpVerts(aVertexIndex);
    }
    //------------------------------------------------------------------------------

    void
    Enriched_Interpolation_Mesh::add_proc_to_comm_table(moris_index aProcRank)
    {
        mXTKModel->get_background_mesh().add_proc_to_comm_table(aProcRank);
    }

    //------------------------------------------------------------------------------
    moris_index
    Enriched_Interpolation_Mesh::get_enr_basis_index_from_enr_basis_id(
            moris_index const & aMeshIndex,
            moris_index const & aBasisId) const
    {
        moris_index tLocalMeshIndex = this->get_local_mesh_index(aMeshIndex);
        auto tIter = mGlobaltoLocalBasisMaps(tLocalMeshIndex).find(aBasisId);
        MORIS_ASSERT(tIter != mGlobaltoLocalBasisMaps(tLocalMeshIndex).end(),"Basis id not in map");
        return tIter->second;
    }
    //------------------------------------------------------------------------------
    moris::Cell<Interpolation_Cell_Unzipped const*>
    Enriched_Interpolation_Mesh::get_enriched_cells_from_base_cell(moris::mtk::Cell const * aBaseCells) const
    {
        moris_index tBaseIndex = aBaseCells->get_index();
        MORIS_ASSERT(tBaseIndex<(moris_index)mBaseCelltoEnrichedCell.size(),"Base Cell index is out of bounds. This index is related to the unenriched interpolation mesh. Make sure enriched cell is not passed into this function");

        uint tNumEnrichedCells = mBaseCelltoEnrichedCell(tBaseIndex).size();
        moris::Cell<Interpolation_Cell_Unzipped const*> tEnrichedCellPtrs(tNumEnrichedCells);

        for(moris::uint i = 0; i<tNumEnrichedCells; i++)
        {
            tEnrichedCellPtrs(i) = (mBaseCelltoEnrichedCell(tBaseIndex)(i));
        }

        return tEnrichedCellPtrs;
    }
    //------------------------------------------------------------------------------
    Interpolation_Vertex_Unzipped const &
    Enriched_Interpolation_Mesh::get_xtk_interp_vertex(moris::uint aVertexIndex) const
    {
        return mEnrichedInterpVerts(aVertexIndex);
    }
    //------------------------------------------------------------------------------
    Matrix<IdMat>
    Enriched_Interpolation_Mesh::convert_indices_to_ids(Matrix<IndexMat> const & aIndices,
            enum EntityRank          aEntityRank) const
    {
        moris::uint tNRow = aIndices.n_rows();
        moris::uint tNCol = aIndices.n_cols();
        Matrix<IdMat> tIds(tNRow,tNCol);
        for(moris::uint i = 0; i < tNRow; i++)
        {
            for(moris::uint j = 0; j<tNCol; j++)
            {
                tIds(i,j) = this->get_glb_entity_id_from_entity_loc_index(aIndices(i,j),aEntityRank);
            }
        }
        return tIds;
    }
    //------------------------------------------------------------------------------
    Matrix<IndexMat>
    Enriched_Interpolation_Mesh::convert_ids_to_indices(Matrix<IdMat> const & aIds,
            enum EntityRank       aEntityRank) const
    {
        moris::uint tNRow = aIds.n_rows();
        moris::uint tNCol = aIds.n_cols();
        Matrix<IdMat> tIndices(tNRow,tNCol);
        for(moris::uint i = 0; i < tNRow; i++)
        {
            for(moris::uint j = 0; j<tNCol; j++)
            {
                tIndices(i,j) = this->get_loc_entity_ind_from_entity_glb_id(aIds(i,j),aEntityRank);
            }
        }

        return tIndices;
    }

    //------------------------------------------------------------------------------
    void
    Enriched_Interpolation_Mesh::convert_enriched_basis_indices_to_ids(moris_index      const & aMeshIndex,
            Matrix<IndexMat> const & aEnrichedIndices,
            Matrix<IdMat>          & aEnrichedIds) const
    {
        moris_index tLocalMeshIndex = this->get_local_mesh_index(aMeshIndex);


        aEnrichedIds.resize(aEnrichedIndices.n_rows(),aEnrichedIndices.n_cols());

        for(moris::uint i = 0; i < aEnrichedIndices.n_rows(); i++)
        {
            for(moris::uint j = 0; j < aEnrichedIndices.n_cols(); j++)
            {
                aEnrichedIds(i,j) = mEnrichCoeffLocToGlob(tLocalMeshIndex)(aEnrichedIndices(i,j));
            }
        }
    }

    //------------------------------------------------------------------------------
    void
    Enriched_Interpolation_Mesh::convert_enriched_basis_indices_to_ids(moris_index            const & aMeshIndex,
            Cell<Matrix<IndexMat>> const & aEnrichedIndices,
            Cell<Matrix<IdMat>>          & aEnrichedIds) const
    {
        aEnrichedIds.resize(aEnrichedIndices.size());

        for(moris::uint i = 0; i < aEnrichedIndices.size(); i++)
        {
            this->convert_enriched_basis_indices_to_ids(aMeshIndex,aEnrichedIndices(i),aEnrichedIds(i));
        }
    }

    //------------------------------------------------------------------------------
    void
    Enriched_Interpolation_Mesh::print() const
    {
        std::cout<<"\n-------------------------------------"<<std::endl;
        this->print_vertex_maps();
        this->print_enriched_cell_maps();
        this->print_enriched_cells();
        this->print_basis_to_enriched_basis();
        std::cout<<"\n-------------------------------------"<<std::endl;
    }

    void Enriched_Interpolation_Mesh::print_enriched_cells() const
    {
        moris::uint tNumCells = this->get_num_entities(EntityRank::ELEMENT);
        std::cout<<"\nEnriched Interpolation Cells:"<<std::endl;
        for(moris::uint i =0; i <tNumCells; i++)
        {
            std::cout<<"   "<<*mEnrichedInterpCells(i);
        }
    }
    void Enriched_Interpolation_Mesh::print_vertex_maps() const
    {
        moris::uint tNumNodes = this->get_num_entities(EntityRank::NODE);
        std::cout<<"\nVertex Map:"<<std::endl;
        for(moris::uint i =0; i <tNumNodes; i++)
        {
            std::cout<<"    Vertex Index: "<<std::setw(9)<<i<< " | Vertex Id: "<<std::setw(9) <<mLocalToGlobalMaps(0)(i)<<std::endl;
        }
    }

    void Enriched_Interpolation_Mesh::print_enriched_cell_maps() const
    {
        moris::uint tNumCells = this->get_num_entities(EntityRank::ELEMENT);
        moris::uint tSpatialDim = mXTKModel->get_spatial_dim();

        std::cout<<"\nCell Map:"<<std::endl;
        for(moris::uint i =0; i <tNumCells; i++)
        {
            std::cout<<"    Cell Index: "<<std::setw(9)<<i<< " | Cell Id: "<<std::setw(9) <<mLocalToGlobalMaps(tSpatialDim)(i)<<std::endl;
        }
    }
    //------------------------------------------------------------------------------
    void Enriched_Interpolation_Mesh::print_basis_to_enriched_basis() const
    {
        for(moris::uint iM = 0; iM < mMeshIndices.numel(); iM++)
        {
            moris::uint tNumBasis = mCoeffToEnrichCoeffs(iM).size();

            std::cout<<"\nBackground Basis to Enriched Basis Indices For Mesh: " << mMeshIndices(iM) <<std::endl;

            for(moris::uint iB = 0; iB < tNumBasis; iB++)
            {
                std::cout<<"    Basis Index: "<< std::setw(9)<<iB<<" | Enriched Indices";

                for(moris::uint iEB = 0; iEB < mCoeffToEnrichCoeffs(iM)(iB).numel(); iEB++)
                {
                    std::cout<<std::setw(9)<<mCoeffToEnrichCoeffs(iM)(iB)(iEB);
                }
                std::cout<<std::endl;
            }
        }
    }
    //------------------------------------------------------------------------------
    void Enriched_Interpolation_Mesh::print_vertex_interpolation() const
    {
        moris::uint tNumVerts = this->get_num_entities(EntityRank::NODE);
        std::cout<<"\nVertex Interpolation:"<<std::endl;
        for(moris::moris_index i =0; i <(moris_index)tNumVerts; i++)
        {
            Interpolation_Vertex_Unzipped const & tVertex = this->get_xtk_interp_vertex(i);

            std::cout<<"\nVertex Id: "<<std::setw(9)<<tVertex.get_id()<<std::endl;
            std::cout<<*tVertex.get_xtk_interpolation(0)<<std::endl;
        }
    }
    void
    Enriched_Interpolation_Mesh::print_basis_information() const
    {
        std::cout<<"\nBasis Information on proc "<<par_rank()<<":"<<std::endl;
        for(moris::moris_index iM = 0 ; iM < (moris_index)mMeshIndices.numel(); iM++)
        {
            std::cout<<" Mesh Index: "<<mMeshIndices(iM)<<std::endl;
            for(moris::moris_index i =0; i <(moris_index)mEnrichCoeffLocToGlob(iM).numel(); i++)
            {
                moris_id tId = mEnrichCoeffLocToGlob(iM)(i);
                moris_index tIndex = this->get_enr_basis_index_from_enr_basis_id(mMeshIndices(iM),tId);
                std::cout<<"    Basis Id: "<<std::setw(9)<<tId<< " | Basis Index: "<<std::setw(9)<<tIndex<<std::endl;
            }
        }
    }

    //------------------------------------------------------------------------------
    void
    Enriched_Interpolation_Mesh::finalize_setup()
    {
        this->setup_local_to_global_maps();
        this->setup_not_owned_vertices();
        this->setup_basis_ownership();

        //    if(par_rank() == 1)
        //    {
        //        moris_id tVertexId = 5024;
        //
        //        moris_index tVertexInd= this->get_loc_entity_ind_from_entity_glb_id(tVertexId,EntityRank::NODE);
        //
        //        std::cout<<"tVertexInd = "<<tVertexInd<<std::endl;
        //
        //
        //        Interpolation_Vertex_Unzipped & tVertex = mEnrichedInterpVerts(tVertexInd);
        //
        //        mtk::Vertex const * tBGVertex = tVertex.get_base_vertex();
        //
        //        Matrix<DDRMat> tCoords = tBGVertex->get_coords();
        //
        //        std::cout<<"tBGVertex Id: "<<tBGVertex->get_id()<<" | Interpolation address: " <<tBGVertex->get_interpolation(0)<<" | Coords: "<<tCoords(0)<< ", "<<tCoords(1)<<std::endl;
        //
        //
        //        const mtk::Vertex_Interpolation * tInterp = tBGVertex->get_interpolation(0);
        //
        //        moris::print(*(tInterp->get_weights()),"tInterp->get_weights()");
        //        moris::print(tInterp->get_owners(),"tInterp->get_owners()");
        //        moris::print(tVertex.get_interpolation(0)->get_owners(),"tVertex.get_interpolation(0)->get_owners()");
        ////
        ////        Interpolation_Mesh & tIPMesh = mXTKModel->get_background_mesh().get_mesh_data();
        ////
        ////        mtk::Vertex const & tHMRVert = tIPMesh.get_mtk_vertex(tBGVertex-)
        //    }
    }
    //------------------------------------------------------------------------------
    void
    Enriched_Interpolation_Mesh::setup_local_to_global_maps()
    {
        // initialize local to global maps
        mLocalToGlobalMaps = Cell<Matrix< IdMat >>( 4 );
        mGlobaltoLobalMaps = Cell<std::unordered_map<moris_id,moris_index>>(4);

        this->setup_vertex_maps();
        this->setup_cell_maps();
        this->setup_basis_maps();
    }
    //------------------------------------------------------------------------------
    void
    Enriched_Interpolation_Mesh::setup_not_owned_vertices()
    {
        // my proc rank
        moris_index tMyProc = par_rank();

        for(moris::uint iV = 0; iV < this->get_num_entities(EntityRank::NODE); iV++)
        {
            moris::mtk::Vertex & tVert = this->get_mtk_vertex((moris_index)iV);

            if(tVert.get_owner() == tMyProc)
            {
                mNotOwnedVerts.push_back(tVert.get_index());
            }
        }
    }
    //------------------------------------------------------------------------------
    void
    Enriched_Interpolation_Mesh::setup_basis_ownership()
    {
        // size data
        mEnrichCoeffOwnership.resize(mMeshIndices.numel());

        // iterate through meshes
        for(moris::uint iM = 0; iM < mMeshIndices.numel(); iM++)
        {
            mEnrichCoeffOwnership(iM).resize(1,mEnrichCoeffLocToGlob(iM).numel());

            mEnrichCoeffOwnership(iM).fill(MORIS_INDEX_MAX);

            // iterate through basis functions
            for(moris::uint iB = 0; iB < mCoeffToEnrichCoeffs(iM).size(); iB++)
            {
                moris_index tOwner = mXTKModel->get_background_mesh().get_mesh_data().get_entity_owner((moris_index)iB,mBasisRank, mMeshIndices(iM));

                for(moris::uint iEB = 0; iEB < mCoeffToEnrichCoeffs(iM)(iB).numel(); iEB++)
                {
                    moris_index tEnrIndex = mCoeffToEnrichCoeffs(iM)(iB)(iEB);
                    mEnrichCoeffOwnership(iM)(tEnrIndex) = tOwner;
                }
            }
        }

    }
    //------------------------------------------------------------------------------
    void
    Enriched_Interpolation_Mesh::setup_vertex_maps()
    {
        moris::uint tNumNodes = this->get_num_entities(EntityRank::NODE);

        mLocalToGlobalMaps(0) = Matrix<IdMat>(tNumNodes,1);

        for(moris::uint i =0; i <tNumNodes; i++)
        {
            mLocalToGlobalMaps(0)(mEnrichedInterpVerts(i).get_index()) = mEnrichedInterpVerts(i).get_id();
            MORIS_ASSERT(mEnrichedInterpVerts(i).get_index() == (moris_index)i,"Index alignment issue in vertices");
            MORIS_ASSERT(mGlobaltoLobalMaps(0).find(mEnrichedInterpVerts(i).get_id()) == mGlobaltoLobalMaps(0).end(),"Duplicate id in the vertex map detected");

            mGlobaltoLobalMaps(0)[mEnrichedInterpVerts(i).get_id()] = mEnrichedInterpVerts(i).get_index();
        }
    }
    //------------------------------------------------------------------------------
    void
    Enriched_Interpolation_Mesh::setup_basis_maps()
    {

        mGlobaltoLocalBasisMaps.resize(mMeshIndices.numel());

        // iterate through meshes
        for(moris::uint iM = 0; iM < mEnrichCoeffLocToGlob.size(); iM++)
        {
            for(moris::uint iB =0; iB <mEnrichCoeffLocToGlob(iM).numel(); iB++)
            {
                MORIS_ASSERT(mGlobaltoLocalBasisMaps(iM).find(mEnrichCoeffLocToGlob(iM)(iB)) == mGlobaltoLocalBasisMaps(iM).end(),"Duplicate id in the basis map detected");

                mGlobaltoLocalBasisMaps(iM)[mEnrichCoeffLocToGlob(iM)(iB)] = (moris_index) iB;

                MORIS_ASSERT(this->get_enr_basis_index_from_enr_basis_id(mMeshIndices(iM),mEnrichCoeffLocToGlob(iM)(iB)) == (moris_index)iB, "Issue setting up the basis map");
            }
        }
    }
    //------------------------------------------------------------------------------
    void
    Enriched_Interpolation_Mesh::setup_cell_maps()
    {
        moris::uint tNumCells = this->get_num_entities(EntityRank::ELEMENT);

        mLocalToGlobalMaps(3) = Matrix<IdMat>(tNumCells,1);

        for(moris::uint i =0; i <tNumCells; i++)
        {
            mLocalToGlobalMaps(3)(mEnrichedInterpCells(i)->get_index()) = mEnrichedInterpCells(i)->get_id();

            MORIS_ASSERT(mEnrichedInterpCells(i)->get_index() == (moris_index)i,"Index alignment issue in cells");
            MORIS_ASSERT(mGlobaltoLobalMaps(3).find(mEnrichedInterpCells(i)->get_id()) == mGlobaltoLobalMaps(3).end(),"Duplicate id in the cell map detected");

            mGlobaltoLobalMaps(3)[mEnrichedInterpCells(i)->get_id()] = mEnrichedInterpCells(i)->get_index();
        }
    }

    void
    Enriched_Interpolation_Mesh::setup_mesh_index_map()
    {
        for(moris::uint i =0; i <mMeshIndices.numel(); i++)
        {
            MORIS_ASSERT(mMeshIndexToLocMeshIndex.find(mMeshIndices(i)) == mMeshIndexToLocMeshIndex.end(),"Duplicate id in the mesh index map detected");

            mMeshIndexToLocMeshIndex[mMeshIndices(i)] = i;
        }
    }

    //------------------------------------------------------------------------------
    moris_id
    Enriched_Interpolation_Mesh::allocate_entity_ids(
            moris::size_t   aNumReqs,
            enum EntityRank aEntityRank)
    {
        MORIS_ASSERT(aEntityRank == EntityRank::NODE || aEntityRank == EntityRank::ELEMENT,"Only Elements or Nodes have ids");

        moris_id tGlobalMax = this->get_max_entity_id(aEntityRank);

        int tProcRank = par_rank();
        int tProcSize = par_size();

        moris::Cell<moris::moris_id> aGatheredInfo;
        moris::Cell<moris::moris_id> tFirstId(1);
        moris::Cell<moris::moris_id> tNumIdsRequested(1);

        tNumIdsRequested(0) = (moris::moris_id)aNumReqs;

        xtk::gather(tNumIdsRequested,aGatheredInfo);

        moris::Cell<moris::moris_id> tProcFirstID(tProcSize);

        if(tProcRank == 0)
        {
            // Loop over entities print the number of entities requested by each processor
            for (int iProc = 0; iProc < tProcSize; ++iProc)
            {
                // Give each processor their desired amount of IDs
                tProcFirstID(iProc) = tGlobalMax;

                // Increment the first available node ID
                tGlobalMax = tGlobalMax+aGatheredInfo(iProc);
            }
        }

        xtk::scatter(tProcFirstID,tFirstId);

        return tFirstId(0);
    }

    //------------------------------------------------------------------------------

    //------------------------------------------------------------------------------
    // multigrid accessor functions
    //------------------------------------------------------------------------------

    uint Enriched_Interpolation_Mesh::get_num_interpolations()
    {
        return mXTKModel->get_multigrid_ptr()->get_num_interpolations();
    }

    //------------------------------------------------------------------------------

    uint Enriched_Interpolation_Mesh::get_max_level( const moris_index aInterpolationIndex )
    {
        return mXTKModel->get_multigrid_ptr()->get_max_level( aInterpolationIndex );
    }

    //------------------------------------------------------------------------------

    uint Enriched_Interpolation_Mesh::get_num_basis( const moris_index aInterpolationIndex )
    {
        return mXTKModel->get_multigrid_ptr()->get_num_basis( aInterpolationIndex );
    }

    //------------------------------------------------------------------------------

    uint Enriched_Interpolation_Mesh::get_basis_level( const moris_index aInterpolationIndex,
            const moris_index aBasisIndex )
    {
        return mXTKModel->get_multigrid_ptr()->get_basis_level( aInterpolationIndex,
                aBasisIndex );
    }

    //------------------------------------------------------------------------------

    uint Enriched_Interpolation_Mesh::get_num_coarse_basis_of_basis(
            const moris_index aInterpolationIndex,
            const moris_index aBasisIndex )
    {
        return mXTKModel->get_multigrid_ptr()->get_num_coarse_basis_of_basis( aInterpolationIndex,
                aBasisIndex );
    }

    //------------------------------------------------------------------------------

    uint Enriched_Interpolation_Mesh::get_coarse_basis_index_of_basis(
            const moris_index aInterpolationIndex,
            const moris_index aBasisIndex,
            const moris_index aCoarseParentIndex )
    {
        return mXTKModel->get_multigrid_ptr()->get_coarse_basis_index_of_basis( aInterpolationIndex,
                aBasisIndex,
                aCoarseParentIndex );
    }

    //------------------------------------------------------------------------------

    moris::Matrix< DDSMat > Enriched_Interpolation_Mesh::get_fine_basis_inds_of_basis(
            const moris_index aInterpolationIndex,
            const moris_index aBasisIndex )
    {
        return mXTKModel->get_multigrid_ptr()->get_fine_basis_inds_of_basis( aInterpolationIndex,
                aBasisIndex );
    }

    //------------------------------------------------------------------------------

    moris::Matrix< DDRMat > Enriched_Interpolation_Mesh::get_fine_basis_weights_of_basis(
            const moris_index aInterpolationIndex,
            const moris_index aBasisIndex )
    {
        return mXTKModel->get_multigrid_ptr()->get_fine_basis_weights_of_basis( aInterpolationIndex,
                aBasisIndex );
    }
    //------------------------------------------------------------------------------
    //------------------------------------------------------------------------------

}
