/*
 * cl_XTK_Enriched_Integration_Mesh.cpp
 *
 *  Created on: Jul 22, 2019
 *      Author: doble
 */

#include "cl_XTK_Enriched_Integration_Mesh.hpp"
#include "cl_XTK_Enriched_Interpolation_Mesh.hpp"
#include "cl_XTK_Model.hpp"
#include "cl_XTK_Cell_Cluster.hpp"
#include "cl_XTK_Side_Cluster.hpp"
#include "cl_MTK_Side_Set.hpp"
#include "cl_MTK_Double_Side_Set.hpp"
#include "cl_MTK_Cell_Info_Hex8.hpp"

#include <memory>

namespace xtk
{
//------------------------------------------------------------------------------
Enriched_Integration_Mesh::Enriched_Integration_Mesh(Model* aXTKModel,
                                                     moris::moris_index aInterpIndex):
                mModel(aXTKModel),
                mMeshIndexInModel(aInterpIndex),
                mCellClusters(0,nullptr),
                mFields(0),
                mFieldLabelToIndex(2)
{

    this->setup_cell_clusters();
    this->setup_blockset_with_cell_clusters();
    this->setup_side_set_clusters();
    this->setup_interface_side_sets();
    this->setup_double_side_set_clusters();
    this->collect_all_sets();
}
//------------------------------------------------------------------------------
Enriched_Integration_Mesh::~Enriched_Integration_Mesh()
{
    for(auto p : mListofBlocks)
    {
        delete p;
    }
    mListofBlocks.clear();

    for(auto p:mListofSideSets)
    {
        delete p;
    }
    mListofSideSets.clear();

    for(auto p:mListofDoubleSideSets)
    {
        delete p;
    }
    mListofDoubleSideSets.clear();


    for(auto p:mDoubleSideClusters)
    {
        delete p;
    }
    mDoubleSideClusters.clear();

    mDoubleSideSingleSideClusters.clear();

}
//------------------------------------------------------------------------------
MeshType
Enriched_Integration_Mesh::get_mesh_type() const
{
    return MeshType::XTK;
}
//------------------------------------------------------------------------------
moris::uint
Enriched_Integration_Mesh::get_spatial_dim() const
{
    return mModel->get_spatial_dim();
}
//------------------------------------------------------------------------------
uint
Enriched_Integration_Mesh::get_num_entities( enum EntityRank aEntityRank ) const
{
    switch(aEntityRank)
    {
        case(EntityRank::NODE):
            {
            return mModel->mBackgroundMesh.get_num_entities(EntityRank::NODE);
            break;
            }
        case(EntityRank::ELEMENT):
            {
            return mModel->get_num_elements_total();
            break;
            }
        default:
            MORIS_ERROR(0,"Only support get num entities for nodes and elements currently");
            return 0;
    }
}
//------------------------------------------------------------------------------
Matrix<IndexMat>
Enriched_Integration_Mesh::get_entity_connected_to_entity_loc_inds(moris_index  aEntityIndex,
                                                                   enum EntityRank aInputEntityRank,
                                                                   enum EntityRank aOutputEntityRank) const
{
    MORIS_ERROR(aInputEntityRank == EntityRank::ELEMENT && aOutputEntityRank == EntityRank::NODE,"Only support element to node connectivity");
    return this->get_mtk_cell(aEntityIndex).get_vertex_inds();
}
//------------------------------------------------------------------------------
Matrix< IndexMat >
Enriched_Integration_Mesh::get_elements_connected_to_element_and_face_ind_loc_inds(moris_index aElementIndex) const
{
    MORIS_ERROR(0,"XTK ENRICHED MESH ERROR: get_elements_connected_to_element_and_face_ind_loc_inds no implemented");
    return Matrix<IndexMat>(0,0);
}
//------------------------------------------------------------------------------
Cell<mtk::Vertex const *>
Enriched_Integration_Mesh::get_all_vertices() const
{
    moris::uint tNumNodes = this->get_num_entities(EntityRank::NODE);
    Cell<mtk::Vertex const *> tVertices(tNumNodes);

    for(moris::uint i = 0; i < tNumNodes; i++)
    {
        tVertices(i) = &mModel->mBackgroundMesh.get_mtk_vertex(i);
    }
    return tVertices;
}
//------------------------------------------------------------------------------
moris_id
Enriched_Integration_Mesh::get_glb_entity_id_from_entity_loc_index(moris_index     aEntityIndex,
                                                                   enum EntityRank aEntityRank) const
{
    return mModel->mBackgroundMesh.get_glb_entity_id_from_entity_loc_index(aEntityIndex,aEntityRank);
}
//------------------------------------------------------------------------------
moris_index
Enriched_Integration_Mesh::get_loc_entity_ind_from_entity_glb_id(moris_id        aEntityId,
                                                                 enum EntityRank aEntityRank) const
{
    return mModel->mBackgroundMesh.get_loc_entity_ind_from_entity_glb_id( aEntityId, aEntityRank );
}

//------------------------------------------------------------------------------

Matrix<IdMat>
Enriched_Integration_Mesh::get_entity_connected_to_entity_glob_ids( moris_id     aEntityId,
                                                                    enum EntityRank aInputEntityRank,
                                                                    enum EntityRank aOutputEntityRank) const
{
    moris_index tEntityIndex = get_loc_entity_ind_from_entity_glb_id(aEntityId,aInputEntityRank);

    Matrix<IndexMat> tEntityToEntityLoc = this->get_entity_connected_to_entity_loc_inds(tEntityIndex,aInputEntityRank,aOutputEntityRank);

    return convert_indices_to_ids(tEntityToEntityLoc,aOutputEntityRank);
}
//------------------------------------------------------------------------------
Matrix< DDRMat >
Enriched_Integration_Mesh::get_node_coordinate( moris_index aNodeIndex ) const
{
    mtk::Vertex const & tVertex = get_mtk_vertex(aNodeIndex);
    return tVertex.get_coords();
}
//------------------------------------------------------------------------------
mtk::Vertex &
Enriched_Integration_Mesh::get_mtk_vertex( moris_index aVertexIndex )
{
    return mModel->mBackgroundMesh.get_mtk_vertex(aVertexIndex);
}
mtk::Vertex const &
Enriched_Integration_Mesh::get_mtk_vertex( moris_index aVertexIndex ) const
{
    return mModel->mBackgroundMesh.get_mtk_vertex(aVertexIndex);
}
//------------------------------------------------------------------------------
mtk::Cell &
Enriched_Integration_Mesh::get_writable_mtk_cell( moris_index aElementIndex )
{
    return mModel->mBackgroundMesh.get_mtk_cell(aElementIndex);
}
//------------------------------------------------------------------------------

mtk::Cell &
Enriched_Integration_Mesh::get_mtk_cell( moris_index aElementIndex )
{
    return mModel->mBackgroundMesh.get_mtk_cell(aElementIndex);
}

mtk::Cell const &
Enriched_Integration_Mesh::get_mtk_cell( moris_index aElementIndex ) const
{
    return mModel->mBackgroundMesh.get_mtk_cell(aElementIndex);
}
//------------------------------------------------------------------------------
Matrix< IdMat >
Enriched_Integration_Mesh::get_communication_table() const
{
    return mModel->mBackgroundMesh.get_mesh_data().get_communication_table();
}

moris::Cell<std::string>
Enriched_Integration_Mesh::get_set_names(enum EntityRank aSetEntityRank) const
{
    switch(aSetEntityRank)
    {
        case(EntityRank::NODE):
            {
            return moris::Cell<std::string>(0);;
            break;
            }
        case(EntityRank::EDGE):
            {
            MORIS_ASSERT(this->get_facet_rank() == EntityRank::EDGE,"side sets are defined on edges in 2d");
            return mSideSetLabels;
            break;
            }
        case(EntityRank::FACE):
            {
            MORIS_ASSERT(this->get_facet_rank() == EntityRank::FACE,"side sets are defined on faces in 3d");
            return mSideSetLabels;
            break;
            }
        case(EntityRank::ELEMENT):
            {
            return mBlockSetNames;
            break;
            }
        default:
            MORIS_ERROR(0,"Currently only supporting block, node and side sets in XTK enriched integration meshes");

            return moris::Cell<std::string>(0);
            break;
    }
}
Matrix< IndexMat >
Enriched_Integration_Mesh::get_set_entity_loc_inds( enum EntityRank aSetEntityRank,
                                                    std::string     aSetName) const
{
    switch(aSetEntityRank)
    {
        case(EntityRank::NODE):
            {
            return Matrix< IndexMat >(0,0);
            break;
            }
        case(EntityRank::EDGE):
            {
            MORIS_ASSERT(this->get_facet_rank() == EntityRank::EDGE,"side sets are defined on edges in 2d");
            return Matrix< IndexMat >(0,0);
            break;
            }
        case(EntityRank::FACE):
            {
            MORIS_ASSERT(this->get_facet_rank() == EntityRank::FACE,"side sets are defined on faces in 3d");
            return Matrix< IndexMat >(0,0);
            break;
            }
        case(EntityRank::ELEMENT):
            {

            return this->get_block_entity_loc_inds(aSetName);

            break;
            }
        default:
            MORIS_ERROR(0,"Currently only supporting block, node and side sets in XTK enriched integration meshes");

            return Matrix< IndexMat >(0,0);
            break;
    }
}

void
Enriched_Integration_Mesh::get_sideset_elems_loc_inds_and_ords( const  std::string & aSetName,
                                                                Matrix< IndexMat > & aElemIndices,
                                                                Matrix< IndexMat > & aSidesetOrdinals ) const
{
    // get the index
    moris_index tSideSetIndex = this->get_side_set_index(aSetName);

    // get the cell clusters
    moris::Cell<mtk::Cluster const *> tSideClusters = this->get_side_set_cluster(tSideSetIndex);

    // iterate through side clusters and count number of sides in set
    moris::uint tNumSides = 0;
    for(auto iCluster:tSideClusters)
    {
        tNumSides = tNumSides + iCluster->get_num_primary_cells();
    }

    // size outputs
    aElemIndices.resize(1,tNumSides);
    aSidesetOrdinals.resize(1,tNumSides);

    // reset count
    tNumSides = 0;
    for(auto iCluster:tSideClusters)
    {
        Matrix<IndexMat> tSideOrdinals = iCluster->get_cell_side_ordinals();
        Matrix<IndexMat> tCellIndices  = iCluster->get_primary_cell_indices_in_cluster();

        aElemIndices({0,0},{tNumSides,tNumSides+tCellIndices.numel()-1}) = tCellIndices.matrix_data();
        aSidesetOrdinals({0,0},{tNumSides,tNumSides+tCellIndices.numel()-1}) = tSideOrdinals.matrix_data();

        tNumSides = tNumSides + tSideOrdinals.numel();
    }

}
//------------------------------------------------------------------------------
moris_id
Enriched_Integration_Mesh::get_max_entity_id( enum EntityRank aEntityRank ) const
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
Matrix< IndexMat >
Enriched_Integration_Mesh::get_block_entity_loc_inds( std::string     aSetName) const
{
    // ge tindex
    moris_index tBlockIndex = this->get_block_set_index(aSetName);

    // get clusters in block
    moris::Cell<mtk::Cluster const *> tCellClusters =  this->get_cell_clusters_in_set(tBlockIndex);

    // iterate through clusters and count all primary integration cells
    moris::uint tCount = 0;
    for(auto it: tCellClusters)
    {
        tCount = tCount + it->get_num_primary_cells();
    }

    // allocate output
    Matrix<IndexMat> tCellIndices(1,tCount);

    // reset count to use it for something else
    tCount = 0;

    // iterate through clusters and collect all primary integration cell indices
    for(auto it: tCellClusters)
    {
        // get primary cell clusters
        moris::Matrix<moris::IndexMat> tPrimaryCellIndices = it->get_primary_cell_indices_in_cluster();

        tCellIndices({0,0},{tCount,tCount+tPrimaryCellIndices.numel()-1}) = tPrimaryCellIndices.matrix_data();

        tCount = tCount + tPrimaryCellIndices.numel();
    }

    return tCellIndices;
}
//------------------------------------------------------------------------------
void
Enriched_Integration_Mesh::deactivate_empty_sets()
{
    this->deactivate_empty_block_sets();
    this->deactivate_empty_side_sets();
    this->collect_all_sets();
}

void
Enriched_Integration_Mesh::deactivate_empty_side_sets()
{
    // copy old data
    std::unordered_map<std::string, moris_index> tOldSetMap = mSideSideSetLabelToOrd;
    moris::Cell<std::string> tOldSetNames = mSideSetLabels;
    moris::Cell<moris::Cell<std::shared_ptr<xtk::Side_Cluster>>> tOldSetClusters = mSideSets;


    // clear member data
    mSideSideSetLabelToOrd.clear();
    mSideSetLabels.clear();
    mSideSets.clear();


    for(auto iB: mListofSideSets)
    {
        delete iB;
    }
    mListofSideSets.clear();

    // current index
    moris_index tSetIndex = 0;

    for(moris::uint i = 0; i < tOldSetClusters.size(); i++)
    {
        if ( tOldSetClusters(i).size() > 0)
        {
            mSideSetLabels.push_back(tOldSetNames(i));
            mSideSets.push_back(tOldSetClusters(i));

            MORIS_ASSERT(mSideSideSetLabelToOrd.find(tOldSetNames(i)) ==  mSideSideSetLabelToOrd.end(),"Duplicate block set in mesh");
            mSideSideSetLabelToOrd[tOldSetNames(i)] = tSetIndex ;
            this->commit_side_set(tSetIndex);
            tSetIndex++;

        }
    }

}


void
Enriched_Integration_Mesh::deactivate_empty_block_sets()
{
    std::unordered_map<std::string, moris_index>        tOldSetMap      = mBlockSetLabelToOrd;
    moris::Cell<std::string>                            tOldSetNames    = mBlockSetNames;
    moris::Cell<enum CellTopology>                      tOldSetTopo     = mBlockSetTopology;
    moris::Cell<moris::Cell<xtk::Cell_Cluster const *>> tOldSetClusters = mPrimaryBlockSetClusters;

    // clear member data
    mBlockSetLabelToOrd.clear();
    mBlockSetNames.clear();
    mBlockSetTopology.clear();
    mPrimaryBlockSetClusters.clear();

    for(auto iB: mListofBlocks)
    {
        delete iB;
    }
    mListofBlocks.clear();

    // current index
    moris_index tSetIndex = 0;

    for(moris::uint i = 0; i < tOldSetClusters.size(); i++)
    {


        if ( tOldSetClusters(i).size() > 0)
        {
            mBlockSetNames.push_back(tOldSetNames(i));
            mPrimaryBlockSetClusters.push_back(tOldSetClusters(i));
            mBlockSetTopology.push_back(tOldSetTopo(i));

            MORIS_ASSERT(mBlockSetLabelToOrd.find(tOldSetNames(i)) ==  mBlockSetLabelToOrd.end(),"Duplicate block set in mesh");
            mBlockSetLabelToOrd[tOldSetNames(i)] = tSetIndex ;
            this->commit_block_set(tSetIndex);
            tSetIndex++;
        }

    }

}
//------------------------------------------------------------------------------
enum CellTopology
Enriched_Integration_Mesh::get_blockset_topology(const  std::string & aSetName)
{
    moris_index tIndex = this->get_block_set_index(aSetName);
    return mBlockSetTopology(tIndex);
}

//------------------------------------------------------------------------------
Matrix<IdMat>
Enriched_Integration_Mesh::convert_indices_to_ids(Matrix<IndexMat> const & aIndices,
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
Enriched_Integration_Mesh::convert_ids_to_indices(Matrix<IdMat> const & aIds,
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
moris::Cell<moris::mtk::Cell const *>
Enriched_Integration_Mesh::get_mtk_cells_loc_inds(Matrix<IndexMat>  const &  aCellIndices)
{
    moris::uint tNumCells = aCellIndices.numel();
    moris::Cell<moris::mtk::Cell const *> tCells(tNumCells);
    for(moris::uint  i = 0; i < tNumCells; i++)
    {
        tCells(i) = &this->get_mtk_cell(aCellIndices(i));

    }

    return tCells;
}
//------------------------------------------------------------------------------
moris::Cell<moris::mtk::Vertex const *>
Enriched_Integration_Mesh::get_mtk_vertices_loc_inds(Matrix<IndexMat> const & aVertexIndices)
{
    moris::uint tNumVerts = aVertexIndices.numel();
    moris::Cell<moris::mtk::Vertex const *> tVertices(tNumVerts);
    for(moris::uint  i = 0; i < tNumVerts; i++)
    {
        tVertices(i) = &this->get_mtk_vertex((moris_index)aVertexIndices(i));
    }

    return tVertices;
}

//------------------------------------------------------------------------------
xtk::Cell_Cluster const &
Enriched_Integration_Mesh::get_xtk_cell_cluster(mtk::Cell const & aInterpCell) const
{
    return get_cell_cluster(aInterpCell.get_index());
}

//------------------------------------------------------------------------------
mtk::Cell_Cluster const &
Enriched_Integration_Mesh::get_cell_cluster(mtk::Cell const & aInterpCell) const
{
    return get_cell_cluster(aInterpCell.get_index());
}
//------------------------------------------------------------------------------
Cell_Cluster const &
Enriched_Integration_Mesh::get_cell_cluster(moris_index aInterpCellIndex) const
{
    MORIS_ASSERT(aInterpCellIndex<(moris_index)mCellClusters.size(),"Interpolation Cell index out of bounds");
    return *mCellClusters(aInterpCellIndex);
}
//------------------------------------------------------------------------------
moris::Cell<std::string>
Enriched_Integration_Mesh::get_block_set_names() const
{
    return mBlockSetNames;
}
std::string
Enriched_Integration_Mesh::get_block_set_label(moris_index aBlockSetOrdinal) const
{
    MORIS_ASSERT(aBlockSetOrdinal<(moris_index)mSideSetLabels.size(),"Block set ordinal out of bounds");
    return mBlockSetNames(aBlockSetOrdinal);
}
moris_index
Enriched_Integration_Mesh::get_block_set_index(std::string aBlockSetLabel) const
{
    auto tIter = mBlockSetLabelToOrd.find(aBlockSetLabel);

    MORIS_ERROR(tIter != mBlockSetLabelToOrd.end(),"block set set label not found");

    return tIter->second;
}

//------------------------------------------------------------------------------
moris::Cell<mtk::Cluster const *>
Enriched_Integration_Mesh::get_cell_clusters_in_set(moris_index aBlockSetOrdinal) const
{
    MORIS_ASSERT(aBlockSetOrdinal<(moris_index)mBlockSetNames.size(),"Requested block set ordinal out of bounds.");

    moris::Cell<xtk::Cell_Cluster const  *> const & tXTKClustersInSet = mPrimaryBlockSetClusters(aBlockSetOrdinal);

    moris::Cell<mtk::Cluster const *> tClusterInSet(tXTKClustersInSet.size());

    for(moris::uint i = 0; i <tXTKClustersInSet.size(); i++)
    {
        tClusterInSet(i) = tXTKClustersInSet(i);
    }

    return tClusterInSet;

}
//------------------------------------------------------------------------------
moris::Cell<mtk::Cluster const *>
Enriched_Integration_Mesh::get_side_set_cluster(moris_index aSideSetOrdinal) const
{
    MORIS_ASSERT(aSideSetOrdinal < (moris_index)mSideSets.size(), "Side set ordinal out of bounds");

    moris::uint tNumSideClustersInSet = mSideSets(aSideSetOrdinal).size();

    moris::Cell<mtk::Cluster const *> tSideClustersInSet(tNumSideClustersInSet);

    for(moris::uint i = 0; i <tNumSideClustersInSet; i++)
    {
        tSideClustersInSet(i) = mSideSets(aSideSetOrdinal)(i).get();
    }

    return tSideClustersInSet;
}
//------------------------------------------------------------------------------
uint
Enriched_Integration_Mesh::get_num_side_sets() const
{
    return mSideSets.size();
}
//------------------------------------------------------------------------------
std::string
Enriched_Integration_Mesh::get_side_set_label(moris_index aSideSetOrdinal) const
{
    MORIS_ASSERT(aSideSetOrdinal<(moris_index)mSideSetLabels.size(),"Side set ordinal out of bounds");
    return mSideSetLabels(aSideSetOrdinal);
}
//------------------------------------------------------------------------------
moris_index
Enriched_Integration_Mesh::get_side_set_index(std::string aSideSetLabel) const
{
    auto tIter = mSideSideSetLabelToOrd.find(aSideSetLabel);

    MORIS_ERROR(tIter != mSideSideSetLabelToOrd.end(),"side side set label not found");

    return tIter->second;
}
//------------------------------------------------------------------------------
uint
Enriched_Integration_Mesh::get_num_double_sided_sets() const
{
    return mDoubleSideSetLabels.size();
}
//------------------------------------------------------------------------------
std::string
Enriched_Integration_Mesh::get_double_sided_set_label(moris_index aSideSetOrdinal) const
{
    MORIS_ASSERT(aSideSetOrdinal<(moris_index)mDoubleSideSetLabels.size(),"Double side set ordinal out of bounds");
    return mDoubleSideSetLabels(aSideSetOrdinal);
}
//------------------------------------------------------------------------------
moris_index
Enriched_Integration_Mesh::get_double_sided_set_index(std::string aDoubleSideSetLabel) const
{
    auto tIter = mDoubleSideSetLabelToOrd.find(aDoubleSideSetLabel);

    MORIS_ERROR(tIter != mDoubleSideSetLabelToOrd.end(),"double side set label not found");

    return tIter->second;
}
//------------------------------------------------------------------------------
moris::Cell<mtk::Cluster const*>
Enriched_Integration_Mesh::get_double_side_set_cluster(moris_index aSideSetOrdinal) const
{
    MORIS_ASSERT(aSideSetOrdinal<(moris_index)mDoubleSideSetLabels.size(),"Double side set ordinal out of bounds");
    moris::Cell<mtk::Cluster const*> tDblSideClusters;
    tDblSideClusters.data() = std::vector<mtk::Cluster const *>(mDoubleSideSets(aSideSetOrdinal).cbegin(),
            mDoubleSideSets(aSideSetOrdinal).cend());

    return tDblSideClusters;
}
//------------------------------------------------------------------------------
uint
Enriched_Integration_Mesh::get_sidesets_num_faces( moris::Cell< moris_index > aSideSetIndex ) const
{
    moris::uint tNumSideSetFaces = 0;

    for( moris_index Ik=0; Ik < (moris_index)aSideSetIndex.size(); ++Ik )
    {
        MORIS_ASSERT(aSideSetIndex(Ik) < (moris_index)mSideSets.size(),"Side set index out of bounds");
        // add up the sideset number of faces
        tNumSideSetFaces = tNumSideSetFaces + mSideSets(aSideSetIndex(Ik)).size();
    }

    return tNumSideSetFaces;
}
void
Enriched_Integration_Mesh::print() const
{
    this->print_block_sets();
    this->print_side_sets();
    this->print_double_side_sets();
}
//------------------------------------------------------------------------------
void
Enriched_Integration_Mesh::print_cell_clusters(moris::uint aVerbosityLevel) const
{
    std::cout<<"\nCell Clusters:"<<std::endl;
    for(moris::uint i =0; i <mCellClusters.size(); i++)
    {
        xtk::Cell_Cluster* tCluster = mCellClusters(i).get();
        std::string tTrivialStr = "f";
        if(tCluster->is_trivial())
        {
            tTrivialStr = "t";
        }

        Interpolation_Cell_Unzipped const * tBaseInterpCell = tCluster->get_xtk_interpolation_cell();

        std::cout<<"    Cluster Index: "<<std::setw(9)<<i
        << " | Interp Cell Id: "     <<std::setw(9)<< tCluster->get_interpolation_cell_id()
        << " | Base Interp Cell Id: "<<std::setw(9)<< tBaseInterpCell->get_base_cell()->get_id()
        << " | Trivial: "                          << tTrivialStr
        << " | Num Primary: "        <<std::setw(9)<< tCluster->get_num_primary_cells()
        << " | Num Void: "                         << tCluster->get_num_void_cells()
        <<std::endl;

        if(aVerbosityLevel > 0)
        {
            moris::Cell<moris::mtk::Cell const *> const & tPrimaryCells = tCluster->get_primary_cells_in_cluster();
            std::cout<<"        Primary Integration Cell Ids: ";
            for(moris::uint  i = 0 ; i < tCluster->get_num_primary_cells(); i++)
            {
                std::cout<<std::setw(9)<<tPrimaryCells(i)->get_id();
            }
            std::cout<<"\n"<<std::endl;
        }
    }
}
//------------------------------------------------------------------------------
void
Enriched_Integration_Mesh::print_block_sets(moris::uint aVerbosityLevel) const
{
    std::cout<<"\nBlock Sets:"<<std::endl;
    std::cout<<"    Num Block Sets: "<<this->get_num_blocks();

    for(moris::uint iBS = 0; iBS < this->get_num_blocks(); iBS++)
    {
        std::cout<<"\n    Block Name: "     <<std::setw(20)<<mBlockSetNames(iBS)
                <<" | Block Set Ord: "    <<std::setw(9)<<iBS
                <<" | Num Cell Clusters: "<<std::setw(9)<<this->mPrimaryBlockSetClusters(iBS).size();

        if(aVerbosityLevel > 0)
        {
            moris::Cell<xtk::Cell_Cluster const *> tClusters = this->mPrimaryBlockSetClusters(iBS);
            std::cout<<"\n            Cluster in set\n";
            for(moris::uint i =0; i <tClusters.size(); i++)
            {
                xtk::Cell_Cluster const * tCluster = tClusters(i);
                std::string tTrivialStr = "f";
                if(tCluster->is_trivial())
                {
                    tTrivialStr = "t";
                }

                Interpolation_Cell_Unzipped const * tBaseInterpCell = tCluster->get_xtk_interpolation_cell();

                std::cout<<"            Cluster Index: "<<std::setw(9)<<i
                << " | Interp Cell Id: "     <<std::setw(9)<< tCluster->get_interpolation_cell_id()
                << " | Base Interp Cell Id: "<<std::setw(9)<< tBaseInterpCell->get_base_cell()->get_id()
                << " | Trivial: "                          << tTrivialStr
                << " | Num Primary: "        <<std::setw(9)<< tCluster->get_num_primary_cells()
                << " | Num Void: "                         << tCluster->get_num_void_cells()
                <<std::endl;
            }

        }

    }
}
//------------------------------------------------------------------------------
void
Enriched_Integration_Mesh::print_side_sets(moris::uint aVerbosityLevel) const
{
    std::cout<<"\nSide Sets:"<<std::endl;
    std::cout<<"    Num Side Sets: "<<this->get_num_side_sets()<<std::endl;

    for(moris::uint iSS = 0; iSS < this->get_num_side_sets(); iSS++)
    {
        std::cout<<"    Side Set Name: "<<std::setw(20)<<mSideSetLabels(iSS)<<" | Side Set Ord: "<<std::setw(9)<<iSS<<" | Num Cell Clusters: "<<std::setw(9)<<this->mSideSets(iSS).size()<<std::endl;
    }
}
//------------------------------------------------------------------------------
void
Enriched_Integration_Mesh::print_double_side_sets(moris::uint aVerbosityLevel) const
{
    std::cout<<"\nDouble Side Sets:"<<std::endl;
    std::cout<<"    Num Side Sets: "<<this->get_num_double_side_set()<<std::endl;

    for(moris::uint iSS = 0; iSS < this->get_num_double_side_set(); iSS++)
    {
        std::cout<<"    Dbl Side Set Name: "<<std::setw(20)<<mDoubleSideSetLabels(iSS)<<" | Dbl Side Set Ord: "<<std::setw(9)<<iSS<<" | Num Cell Clusters: "<<std::setw(9)<<this->mDoubleSideSets(iSS).size();

        if(aVerbosityLevel>0)
        {
            for(moris::uint  i = 0; i < mDoubleSideSets(iSS).size(); i++)
            {
             std::cout<<"\n      Master Interpolation Cell: "<<std::setw(9)<<mDoubleSideSets(iSS)(i)->get_interpolation_cell( mtk::Master_Slave::MASTER ).get_id();
             std::cout<<" | Slave Interpolation Cell: "<<std::setw(9)<<mDoubleSideSets(iSS)(i)->get_interpolation_cell( mtk::Master_Slave::SLAVE ).get_id();
            }
        }

        std::cout<<std::endl;
    }
}
void
Enriched_Integration_Mesh::print_double_side_clusters(moris::uint aVerbosityLevel) const
{
    std::cout<<"\nDouble Side Clusters:"<<std::endl;
    std::cout<<"    Num Double Side Clusters: "<<mDoubleSideClusters.size()<<std::endl;

    for(moris::uint i = 0; i < mDoubleSideClusters.size(); i++)
    {
        std::cout<<mDoubleSideClusters(i)<<std::endl;
    }
}
//------------------------------------------------------------------------------
moris_index
Enriched_Integration_Mesh::create_side_set_from_dbl_side_set(moris_index const & aDblSideSetIndex,
                                                             std::string const & aSideSetName)
{
    Cell<moris_index> tSideSetIndex = this->register_side_set_names({aSideSetName});

    moris::Cell<mtk::Double_Side_Cluster*> & tDblSideClusters = mDoubleSideSets(aDblSideSetIndex);

    moris::uint tCount = 0;

    for(moris::uint i = 0 ; i < tDblSideClusters.size(); i++)
    {
        // get the index
        moris_index tMasterIndex = mDoubleSideSetsMasterIndex(aDblSideSetIndex)(i);
        moris_index tSlaveIndex  = mDoubleSideSetsSlaveIndex(aDblSideSetIndex)(i);

        mSideSets(tSideSetIndex(0)).push_back( mDoubleSideSingleSideClusters(tMasterIndex) );
        mSideSets(tSideSetIndex(0)).push_back( mDoubleSideSingleSideClusters(tSlaveIndex) );
        tCount++;

    }

    this->commit_side_set(tSideSetIndex(0));

    this->collect_all_sets();

    return tSideSetIndex(0);
}
//------------------------------------------------------------------------------
moris_index
Enriched_Integration_Mesh::create_block_set_from_cells_of_side_set(moris_index const & aSideSetIndex,
                                                                   std::string const & aBlockSetName,
                                                                   enum CellTopology   aCellTopo)
{

    moris::Cell<std::shared_ptr<xtk::Side_Cluster>> & tSideClusters = mSideSets(aSideSetIndex);

    Cell<moris_index> tBlockSetIndex = this->register_block_set_names_with_cell_topo({aBlockSetName},aCellTopo);


    std::unordered_map<moris_index, moris_index> tIpCellInSet;
    for(moris::uint i = 0 ; i < tSideClusters.size(); i++)
    {
        // cast to xtk side cluster
        xtk::Side_Cluster* tSideCluster = tSideClusters(i).get();

        if(tIpCellInSet.find(tSideCluster->mIntegrationCells(0)->get_id()) == tIpCellInSet.end())
        {
            // create a new cell cluster
            std::shared_ptr<xtk::Cell_Cluster>  tCellCluster = std::make_shared< Cell_Cluster >();

            // get the ith enriched interpolation cell
            tCellCluster->mInterpolationCell = tSideCluster->mInterpolationCell;

            // mark as trivial
            tCellCluster->mTrivial = true;

            // add interp cell as integration cell
            tCellCluster->mPrimaryIntegrationCells.append(tSideCluster->mIntegrationCells);

            mCellClusters.push_back(tCellCluster);

            mPrimaryBlockSetClusters(tBlockSetIndex(0)).push_back( tCellCluster.get() );

            tIpCellInSet[tSideCluster->mIntegrationCells(0)->get_id()] = i;

        }
    }

    this->commit_block_set(tBlockSetIndex(0));

    this->collect_all_sets();

    return tBlockSetIndex(0);
}

//------------------------------------------------------------------------------
std::string
Enriched_Integration_Mesh::get_interface_side_set_name(moris_index aGeomIndex,
                                                       moris_index aBulkPhaseIndex0,
                                                       moris_index aBulkPhaseIndex1)
{
    MORIS_ASSERT(aGeomIndex< (moris_index)mModel->get_geom_engine().get_num_geometries(),"Geometry index out of bounds");
    MORIS_ASSERT(aBulkPhaseIndex0< (moris_index)mModel->get_geom_engine().get_num_bulk_phase(),"Bulk phase index 0 out of bounds");
    MORIS_ASSERT(aBulkPhaseIndex1< (moris_index)mModel->get_geom_engine().get_num_bulk_phase(),"Bulk phase index 1 out of bounds");

    return "iside_g_" + std::to_string(aGeomIndex) + "_b0_" + std::to_string(aBulkPhaseIndex0) + "_b1_" + std::to_string(aBulkPhaseIndex1);
}

std::string
Enriched_Integration_Mesh::get_dbl_interface_side_set_name(moris_index aBulkPhaseIndex0,
                                                           moris_index aBulkPhaseIndex1)
{
    MORIS_ASSERT(aBulkPhaseIndex0<aBulkPhaseIndex1,"Double side set names are defined from low index to high");
    MORIS_ASSERT(aBulkPhaseIndex0< (moris_index)mModel->get_geom_engine().get_num_bulk_phase(),"Bulk phase index 0 out of bounds");
    MORIS_ASSERT(aBulkPhaseIndex1< (moris_index)mModel->get_geom_engine().get_num_bulk_phase(),"Bulk phase index 1 out of bounds");

    return "dbl_iside_p0_" + std::to_string(aBulkPhaseIndex0) + "_p1_" + std::to_string(aBulkPhaseIndex1);
}
//------------------------------------------------------------------------------
moris::moris_index
Enriched_Integration_Mesh::create_field(std::string            aLabel,
                                        enum moris::EntityRank aEntityRank,
                                        moris::moris_index     aBulkPhaseIndex)
{
    MORIS_ASSERT(!field_exists(aLabel,aEntityRank),"Field already created");

    moris::moris_index tFieldIndex = mFields.size();
    mFieldLabelToIndex(this->get_entity_rank_field_index(aEntityRank))[aLabel] = tFieldIndex;
    mFields.push_back(Field(aLabel,aBulkPhaseIndex));

    return tFieldIndex;
}
//------------------------------------------------------------------------------
moris::moris_index
Enriched_Integration_Mesh::get_field_index(std::string              aLabel,
                                           enum moris::EntityRank   aEntityRank)
{
    MORIS_ASSERT(field_exists(aLabel,aEntityRank),"Field does not exist in mesh");

    moris_index tIndex = get_entity_rank_field_index(aEntityRank);
    auto tIter = mFieldLabelToIndex(tIndex).find(aLabel);
    return tIter->second;
}

//------------------------------------------------------------------------------
void
Enriched_Integration_Mesh::add_field_data(moris::moris_index       aFieldIndex,
                                          enum moris::EntityRank   aEntityRank,
                                          Matrix<DDRMat>  const  & aFieldData)
{
    mFields(aFieldIndex).mFieldData = aFieldData.copy();
}

//------------------------------------------------------------------------------
moris_id
Enriched_Integration_Mesh::allocate_entity_ids(moris::size_t  aNumReqs,
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
void
Enriched_Integration_Mesh::commit_double_side_set(moris_index const & aDoubleSideSetIndex)
{

    MORIS_ASSERT(mListofDoubleSideSets.size() == (uint)aDoubleSideSetIndex,"Committing double side set failed. aDoubleSideSetIndex needs to be equivalent to the size of the list of double side sets");
    mListofDoubleSideSets.resize( mListofDoubleSideSets.size()+1, nullptr );

    mListofDoubleSideSets( aDoubleSideSetIndex ) = new moris::mtk::Double_Side_Set(mDoubleSideSetLabels(aDoubleSideSetIndex), this->get_double_side_set_cluster( aDoubleSideSetIndex ), this->get_spatial_dim());
}
void
Enriched_Integration_Mesh::commit_side_set(moris_index const & aSideSetIndex)
{
    MORIS_ASSERT(mListofSideSets.size() == (uint)aSideSetIndex,"Committing side set failed. aSideSetIndex needs to be equivalent to the size of the list of double side sets");
    mListofSideSets.resize( mListofSideSets.size()+1, nullptr );

    mListofSideSets( aSideSetIndex ) = new moris::mtk::Side_Set(mSideSetLabels(aSideSetIndex), this->get_side_set_cluster( aSideSetIndex ), this->get_spatial_dim());
}
//------------------------------------------------------------------------------
void
Enriched_Integration_Mesh::commit_block_set(moris_index const & aBlockSetIndex)
{
    MORIS_ASSERT(mListofBlocks.size() == (uint)aBlockSetIndex,"Committing side set failed. aSideSetIndex needs to be equivalent to the size of the list of double side sets");
    mListofBlocks.resize( mListofBlocks.size()+1, nullptr );

    mListofBlocks( aBlockSetIndex ) = new moris::mtk::Block(mBlockSetNames(aBlockSetIndex), this->get_cell_clusters_in_set( aBlockSetIndex ), this->get_spatial_dim());
}
//------------------------------------------------------------------------------
void
Enriched_Integration_Mesh::setup_cell_clusters()
{
    Enriched_Interpolation_Mesh* tEnrInterpMesh  = mModel->mEnrichedInterpMesh(mMeshIndexInModel);
    Background_Mesh &            tBackgroundMesh = mModel->mBackgroundMesh;
    Cut_Mesh        &            tCutMesh        = mModel->mCutMesh;

    // Number of interpolation cells
    moris::uint tNumInterpCells = tEnrInterpMesh->get_num_entities(EntityRank::ELEMENT);

    // allocate cell cluster member data
    mCellClusters.resize(tNumInterpCells,nullptr);

    // Allocate subphase index to cluster index
    mSubphaseIndexToClusterIndex.resize(1,tNumInterpCells);

    // reference the enriched cells
    Cell<Interpolation_Cell_Unzipped*> const & tEnrichedInterpCells = tEnrInterpMesh->get_enriched_interpolation_cells();

    // iterate through interpolation cells to create cell clusters
    for(moris::uint i = 0; i <tNumInterpCells; i++)
    {
        // index
        moris_index tInterpCellIndex = tEnrichedInterpCells(i)->get_index();

        // create a new cell cluster
        mCellClusters(tInterpCellIndex) = std::make_shared< Cell_Cluster >();


        // get the ith enriched interpolation cell
        mCellClusters(tInterpCellIndex)->mInterpolationCell = tEnrichedInterpCells(i);


        // ad subphase index to cluster index to subphase index data
        mSubphaseIndexToClusterIndex(i) = mCellClusters(tInterpCellIndex)->get_xtk_interpolation_cell()->get_subphase_index();

        // base cell
        moris::mtk::Cell const * tBaseInterpCell = mCellClusters(tInterpCellIndex)->mInterpolationCell->get_base_cell();

        // ask background mesh if the base cell has children (the opposite answer to this question is the trivial flag)
        mCellClusters(tInterpCellIndex)->mTrivial = !tBackgroundMesh.entity_has_children(tBaseInterpCell->get_index(),EntityRank::ELEMENT);

        // if it has children get a pointer to the child mesh
        if(!mCellClusters(tInterpCellIndex)->mTrivial)
        {
            // subphase index
            moris_index tProcSubphaseIndex = mCellClusters(tInterpCellIndex)->mInterpolationCell->get_subphase_index();


            moris_index tChildMeshIndex = tBackgroundMesh.child_mesh_index(tBaseInterpCell->get_index(),EntityRank::ELEMENT);
            mCellClusters(tInterpCellIndex)->mChildMesh = & tCutMesh.get_child_mesh(tChildMeshIndex);

            moris_index tSubphaseIndex = mCellClusters(tInterpCellIndex)->mChildMesh->get_subphase_loc_index(tProcSubphaseIndex);

            // access the subphase information
            Cell<moris::Matrix< moris::IndexMat >> tSubPhaseGroups = mCellClusters(tInterpCellIndex)->mChildMesh->get_subphase_groups();

            // child cell proc inds
            Matrix<IndexMat> const & tChildCellInds = mCellClusters(tInterpCellIndex)->mChildMesh->get_element_inds();


            // convert to cell indices as child mesh returns them ordered by cm index
            for(moris::uint iSp = 0; iSp < tSubPhaseGroups.size(); iSp++ )
            {
                moris::uint tNumInSubPhase = tSubPhaseGroups(iSp).numel();
                for(moris::uint iC = 0; iC < tNumInSubPhase; iC++)
                {
                    tSubPhaseGroups(iSp)(iC) = tChildCellInds(tSubPhaseGroups(iSp)(iC) );
                }
            }

            // get cells in primary subphase
            mCellClusters(tInterpCellIndex)->mPrimaryIntegrationCells = this->get_mtk_cells_loc_inds(tSubPhaseGroups(tSubphaseIndex));

            // add other cells to void subphase
            for(moris::uint iSp = 0; iSp < tSubPhaseGroups.size(); iSp++ )
            {
                if(iSp != (uint) tSubphaseIndex)
                {
                    mCellClusters(tInterpCellIndex)->mVoidIntegrationCells.append(this->get_mtk_cells_loc_inds(tSubPhaseGroups(iSp)));
                }
            }

            mCellClusters(tInterpCellIndex)->mVerticesInCluster = this->get_mtk_vertices_loc_inds(mCellClusters(tInterpCellIndex)->mChildMesh->get_node_indices());

        }

        // trivial case, the enriched interpolation cell becomes the primary cell
        else
        {
            mCellClusters(tInterpCellIndex)->mPrimaryIntegrationCells.push_back(mCellClusters(tInterpCellIndex)->mInterpolationCell->get_base_cell());
        }
    }
}
//------------------------------------------------------------------------------
void
Enriched_Integration_Mesh::setup_blockset_with_cell_clusters()
{
    // get background mesh
    Background_Mesh & tBackgroundMesh = mModel->get_background_mesh();

    // enriched interpolation mesh
    Enriched_Interpolation_Mesh* tEnrInterpMesh = mModel->mEnrichedInterpMesh(mMeshIndexInModel);

    // get block sets (in background mesh data)
    Cell<std::string> tBlockSetsNames = tBackgroundMesh.get_mesh_data().get_set_names(EntityRank::ELEMENT);

    // for each block set construct
    for(moris::uint iBS = 0; iBS < tBlockSetsNames.size(); iBS++)
    {
        // split set into child and no child as we need to have the same type of integration cell in each set
        moris::Cell<std::string> tChildNoChildSetNames = this->split_set_name_by_child_no_child(tBlockSetsNames(iBS));

        // split child and no child sets by phases
        moris::Cell<std::string> tPhaseChildBlockSetNames = this->split_set_name_by_bulk_phase(tChildNoChildSetNames(0));
        moris::Cell<std::string> tPhaseNoChildBlockSetNames = this->split_set_name_by_bulk_phase(tChildNoChildSetNames(1));

        // topology enums
        enum CellTopology tChildTopo   = mModel->get_cut_mesh().get_child_element_topology();
        enum CellTopology tParentTopo  = mModel->get_background_mesh().get_parent_cell_topology();

        // add block set names to member data
        Cell<moris_index> tChildBlockSetOrds   = this->register_block_set_names_with_cell_topo(tPhaseChildBlockSetNames,tChildTopo);
        Cell<moris_index> tNoChildBlockSetOrds = this->register_block_set_names_with_cell_topo(tPhaseNoChildBlockSetNames,tParentTopo);

        // get the cells in this block
        moris::Cell<moris::mtk::Cell const*> tCellsInBlock = tBackgroundMesh.get_mesh_data().get_block_set_cells(tBlockSetsNames(iBS));

        // get the enriched interpolation cells in this block
        moris::Cell<xtk::Interpolation_Cell_Unzipped const * > tEnrichedCellsInBlock = tEnrInterpMesh->get_enriched_cells_from_base_cells(tCellsInBlock);

        // iterate through and add cluster associated with enriched cell to block set
        for(moris::uint iC = 0; iC < tEnrichedCellsInBlock.size(); iC++)
        {
            // get the bulkphase
            moris_index tBulkPhaseIndex = tEnrichedCellsInBlock(iC)->get_bulkphase_index();

            // get cluster associated with enriched cell
            xtk::Cell_Cluster const & tCluster = this->get_cell_cluster(tEnrichedCellsInBlock(iC)->get_index());

            // set ord
            moris_index tSetOrd = MORIS_INDEX_MAX;

            if(tCluster.is_trivial())
            {
                tSetOrd = tNoChildBlockSetOrds(tBulkPhaseIndex);
            }

            else
            {
                tSetOrd = tChildBlockSetOrds(tBulkPhaseIndex);
            }

            // add to member data
            mPrimaryBlockSetClusters(tSetOrd).push_back(&tCluster);
        }
    }

    // create mtk set information
    mListofBlocks.resize( mPrimaryBlockSetClusters.size(), nullptr );

    moris::Cell<std::string> tBSNames = this->get_block_set_names();

    for(moris::uint Ik = 0; Ik<mListofBlocks.size(); Ik++)
    {
        mListofBlocks( Ik ) = new moris::mtk::Block( tBSNames(Ik), this->get_cell_clusters_in_set( Ik ), this->get_spatial_dim() );
    }
}


//------------------------------------------------------------------------------
void
Enriched_Integration_Mesh::setup_side_set_clusters()
{
    // get data for easy access
    Enriched_Interpolation_Mesh* tEnrInterpMesh  = mModel->mEnrichedInterpMesh(mMeshIndexInModel);
    Background_Mesh &            tBackgroundMesh = mModel->mBackgroundMesh;
    Cut_Mesh        &            tCutMesh        = mModel->mCutMesh;

    // rank enum for facets
    enum EntityRank tFacetRank = mModel->mBackgroundMesh.get_mesh_data().get_facet_rank();

    // get side sets (in background mesh data)
    Cell<std::string> tSideSetNames = tBackgroundMesh.get_mesh_data().get_set_names(tFacetRank);

    tSideSetNames = mModel->check_for_and_remove_internal_seacas_side_sets(tSideSetNames);

    // my proc rank
    moris_index tParRank = par_rank();

    // for each side set construct
    for(moris::uint iSS = 0; iSS < tSideSetNames.size(); iSS++)
    {
        // split set into child and no child as we need to have the same type of integration cell in each set
        moris::Cell<std::string> tChildNoChildSetNames = this->split_set_name_by_child_no_child(tSideSetNames(iSS));

        // split child and no child sets by phases
        moris::Cell<std::string> tPhaseChildSideSetNames = this->split_set_name_by_bulk_phase(tChildNoChildSetNames(0));
        moris::Cell<std::string> tPhaseNoChildSideSetNames = this->split_set_name_by_bulk_phase(tChildNoChildSetNames(1));

        // add side set names to member data
        Cell<moris_index> tChildSideSetOrds = this->register_side_set_names(tPhaseChildSideSetNames);
        Cell<moris_index> tNoChildSideSetOrds = this->register_side_set_names(tPhaseNoChildSideSetNames);

        // get the cells in this side set and their side ordinals
        moris::Cell< mtk::Cell const * > tCellsInSideSet;
        Matrix< IndexMat >               tCellOrdsInSideSet;
        tBackgroundMesh.get_mesh_data().get_sideset_cells_and_ords(tSideSetNames(iSS),tCellsInSideSet,tCellOrdsInSideSet);

        // iterate through cells in side set
        for(moris::uint iC = 0; iC<tCellsInSideSet.size(); iC++)
        {
            mtk::Cell const * tBaseCell = tCellsInSideSet(iC);
            moris_index       tSideOrd  = tCellOrdsInSideSet(iC);

            if(tBaseCell->get_owner() == tParRank)
            {

                // ask background mesh if the base cell has children (the opposite answer to this question is the trivial flag)
                bool tTrivial = !tBackgroundMesh.entity_has_children(tBaseCell->get_index(),EntityRank::ELEMENT);

                // get the enriched interpolation cells associated with base cell
                moris::Cell<xtk::Interpolation_Cell_Unzipped const * > tEnrichedCellsOfBaseCell = tEnrInterpMesh->get_enriched_cells_from_base_cell(tBaseCell);

                // if there is a child mesh associated with this base cell
                if(!tTrivial)
                {
                    // get the face index associated with the side ordinal
                    moris_index tSideIndex = tBackgroundMesh.get_mesh_data().get_entity_connected_to_entity_loc_inds(tBaseCell->get_index(),EntityRank::ELEMENT,tFacetRank)(tSideOrd);

                    // get the child mesh
                    moris_index tChildMeshIndex   = tBackgroundMesh.child_mesh_index(tBaseCell->get_index(),EntityRank::ELEMENT);
                    Child_Mesh* tChildMesh = & tCutMesh.get_child_mesh(tChildMeshIndex);

                    MORIS_ASSERT(tEnrichedCellsOfBaseCell.size() == tChildMesh->get_num_subphase_bins(),"Number of enriched interpolation cells and subphase bins does not match");

                    // get child element indices and side ordinals on face
                    Matrix<IdMat>    tChildCellIdsOnFace;
                    Matrix<IndexMat> tChildCellsCMIndOnFace;
                    Matrix<IndexMat> tChildCellsOnFaceOrdinal;
                    tChildMesh->get_child_elements_connected_to_parent_facet(tSideIndex, tChildCellIdsOnFace, tChildCellsCMIndOnFace, tChildCellsOnFaceOrdinal);

                    // child cell indices
                    Matrix<IndexMat> const & tChildCellInds = tChildMesh->get_element_inds();

                    // get child cell pointers
                    moris::Cell<moris::mtk::Cell const *> tChildCells = this->get_mtk_cells_loc_inds(tChildCellInds);

                    // create a side cluster for each subphase in this child mesh
                    moris::Cell<std::shared_ptr<xtk::Side_Cluster>> tSideClustersForCM(tChildMesh->get_num_subphase_bins());
                    for(moris::uint  iSP = 0; iSP < tChildMesh->get_num_subphase_bins(); iSP++)
                    {
                        MORIS_ASSERT(tEnrichedCellsOfBaseCell(iSP)->get_subphase_index() == tChildMesh->get_subphase_indices()(iSP),"Enriched interpolation cell subphases associated with a base cell should be in ascending order.");

                        tSideClustersForCM(iSP) = std::make_shared< Side_Cluster >();
                        tSideClustersForCM(iSP)->mInterpolationCell = tEnrichedCellsOfBaseCell(iSP);
                        tSideClustersForCM(iSP)->mTrivial = false;
                        tSideClustersForCM(iSP)->mIntegrationCellSideOrdinals = Matrix<IndexMat>(1,tChildCellIdsOnFace.numel());
                        tSideClustersForCM(iSP)->mChildMesh = tChildMesh;
                        tSideClustersForCM(iSP)->mAssociatedCellCluster = &this->get_cell_cluster(*tEnrichedCellsOfBaseCell(iSP));
                    }


                    // iterate through child cells on face
                    for(moris::uint iF = 0; iF < tChildCellsCMIndOnFace.numel(); iF++)
                    {
                        moris_index tSubphaseGroup      = tChildMesh->get_element_subphase_index(tChildCellsCMIndOnFace(iF));
                        moris_index tCMCellIndex        = tChildCellsCMIndOnFace(iF);
                        moris_index tIndexInSideCluster = tSideClustersForCM(tSubphaseGroup)->mIntegrationCells.size();

                        // add information to side cluster
                        tSideClustersForCM(tSubphaseGroup)->mIntegrationCellSideOrdinals(tIndexInSideCluster) = tChildCellsOnFaceOrdinal(iF);
                        tSideClustersForCM(tSubphaseGroup)->mIntegrationCells.push_back(tChildCells(tCMCellIndex));
                    }

                    // iterate through, get rid of extra space in side ordinals and add to side set clusters data
                    for(moris::uint  iSP = 0; iSP < tChildMesh->get_num_subphase_bins(); iSP++)
                    {
                        // only add this side cluster to the set if it has at least one integration cell in it
                        if(tSideClustersForCM(iSP)->mIntegrationCells.size() > 0 )
                        {

                            tSideClustersForCM(iSP)->mIntegrationCellSideOrdinals.resize(1,tSideClustersForCM(iSP)->mIntegrationCells.size());

                            // bulk phase of this cluster
                            moris_index tBulkPhase = tSideClustersForCM(iSP)->mInterpolationCell->get_bulkphase_index();

                            // side set ordinal in mesh
                            moris_index tSideSetOrd = tChildSideSetOrds(tBulkPhase);

                            tSideClustersForCM(iSP)->finalize_setup();

                            // add
                            mSideSets(tSideSetOrd).push_back(tSideClustersForCM(iSP));
                        }

                    }
                }

                else
                {
                    MORIS_ASSERT(tEnrichedCellsOfBaseCell.size() == 1,"For the trivial case, a base cell should have 1 enriched interpolation cell associated with it");

                    // phase of cell
                    moris_index tBulkPhase = tEnrichedCellsOfBaseCell(0)->get_bulkphase_index();

                    // side set ordinal in mesh
                    moris_index tSideSetOrd = tNoChildSideSetOrds(tBulkPhase);

                    // create a new side cluster in the side set assoicate with this bulk phase
                    moris_index tIndex = mSideSets(tSideSetOrd).size();
                    mSideSets(tSideSetOrd).push_back(std::make_shared< Side_Cluster >());
                    std::shared_ptr<Side_Cluster> tSideCluster = mSideSets(tSideSetOrd)(tIndex);

                    // set trivial flag
                    tSideCluster->mTrivial = tTrivial;

                    // get the set enriched interpolation cell
                    tSideCluster->mInterpolationCell =tEnrichedCellsOfBaseCell(0);

                    // mark child mesh as nullptr
                    tSideCluster->mChildMesh = nullptr;

                    // integration cell is the same as the interpolation cell in this case
                    tSideCluster->mIntegrationCells = moris::Cell<moris::mtk::Cell const *>(1,tSideCluster->mInterpolationCell->get_base_cell());

                    // side ordinal
                    tSideCluster->mIntegrationCellSideOrdinals = Matrix<IndexMat>({{tSideOrd}});

                    // add vertices
                    tSideCluster->mVerticesInCluster.append(tSideCluster->mInterpolationCell->get_vertices_on_side_ordinal(tSideOrd));
                    tSideCluster->finalize_setup();

                    tSideCluster->mAssociatedCellCluster = &this->get_cell_cluster(*tEnrichedCellsOfBaseCell(0));
                }
            }
        }
    }

    mListofSideSets.resize( mSideSets.size(), nullptr );

    for(moris::uint Ik = 0; Ik<mListofSideSets.size(); Ik++)
    {
        mListofSideSets( Ik ) = new moris::mtk::Side_Set( mSideSetLabels(Ik),this->get_side_set_cluster( Ik ), this->get_spatial_dim() );
    }

}
//------------------------------------------------------------------------------
void
Enriched_Integration_Mesh::setup_double_side_set_clusters()
{
    this->setup_double_sided_interface_sides();
}
//------------------------------------------------------------------------------
void
Enriched_Integration_Mesh::setup_double_sided_interface_sides()
{
    this->declare_interface_double_side_sets();

    this->create_interface_double_side_sets_and_clusters();
}

void
Enriched_Integration_Mesh::declare_interface_double_side_sets()
{
    uint tNumBulkPhases = mModel->get_geom_engine().get_num_bulk_phase();

    Cell<std::string> tDoubleInterfaceSideNames;

    mBulkPhaseToDblSideIndex.resize(tNumBulkPhases,tNumBulkPhases);
    mBulkPhaseToDblSideIndex.fill(MORIS_INDEX_MAX);

    moris_index tCount = 0;

    for(moris::moris_index iP0 = 0; iP0 <(moris_index) tNumBulkPhases; iP0++)
    {
        for(moris::moris_index iP1 = iP0+1; iP1 < (moris_index)tNumBulkPhases; iP1++)
        {

            std::string tInterfaceSideSetName = this->get_dbl_interface_side_set_name(iP0,iP1);

            tDoubleInterfaceSideNames.push_back(tInterfaceSideSetName);

            mBulkPhaseToDblSideIndex(iP0,iP1) = tCount;
            mBulkPhaseToDblSideIndex(iP1,iP0) = tCount;
            tCount++;
        }
    }
    this->register_double_side_set_names(tDoubleInterfaceSideNames);
}

moris_index
Enriched_Integration_Mesh::get_dbl_side_set_index(moris_index aPhase0,
                                                  moris_index aPhase1)
{
    MORIS_ASSERT(aPhase0<aPhase1,"Double side sets are defined from low phase index to high");

    MORIS_ASSERT(mDoubleSideSetLabels(mBulkPhaseToDblSideIndex(aPhase0,aPhase1)) == this->get_dbl_interface_side_set_name(aPhase0,aPhase1),"Interface double side set not showing up in correct index");
    return mBulkPhaseToDblSideIndex(aPhase0,aPhase1);
}
//------------------------------------------------------------------------------
void
Enriched_Integration_Mesh::create_interface_double_side_sets_and_clusters()
{
    // access interpolation mesh
    Enriched_Interpolation_Mesh* tEnrInterpMesh  = mModel->mEnrichedInterpMesh(mMeshIndexInModel);

    // background mesh
    moris_index tMyProcRank = par_rank();
    moris::mtk::Interpolation_Mesh & tBGMesh = mModel->get_background_mesh().get_mesh_data();

    uint tNumChildMeshes = mModel->get_cut_mesh().get_num_child_meshes();
    // iterate through children meshes
    for(moris::uint iCM = 0; iCM < tNumChildMeshes; iCM++)
    {
        // get the child mesh
        Child_Mesh * tChildMesh = &mModel->get_cut_mesh().get_child_mesh((moris_index)iCM);

        if(tBGMesh.get_entity_owner(tChildMesh->get_parent_element_index(),EntityRank::ELEMENT) == tMyProcRank)
        {

            // tell the child mesh to construct its double side sets between subphases
            tChildMesh->construct_double_sides_between_subphases();

            // get the child mesh parent cell pointer
            moris::mtk::Cell const * tBaseCell = &this->get_mtk_cell(tChildMesh->get_parent_element_index());

            // get the enriched interpolation cells associated with base cell
            moris::Cell<xtk::Interpolation_Cell_Unzipped const * > tEnrichedCellsOfBaseCell = tEnrInterpMesh->get_enriched_cells_from_base_cell(tBaseCell);

            // child cell indices
            Matrix<IndexMat> const & tChildCellInds = tChildMesh->get_element_inds();

            // get child cell pointers
            moris::Cell<moris::mtk::Cell const *> tChildCells = this->get_mtk_cells_loc_inds(tChildCellInds);

            // create side clusters
            uint tNumDblSideClusters = tChildMesh->get_num_double_side_interfaces();

            // get the subphase bulkphase indices
            Cell<moris::moris_index> const & tSubphaseBulkPhases = tChildMesh->get_subphase_bin_bulk_phase();

            // iterate through clusters and construct side clusters, then add to appropriate side set
            for(moris::moris_index iDI = 0; iDI < (moris_index)tNumDblSideClusters; iDI++)
            {
                // access the dbl side cluster information
                moris::Cell<moris_index>              tDblSideClustSubphaseCMInds = tChildMesh->get_double_side_interface_subphase_indices(iDI);
                moris::Cell<moris::Cell<moris_index>> tDblSideClustCellCMInds     = tChildMesh->get_double_side_interface_cell_pairs(iDI);
                moris::Cell<moris::Cell<moris_index>> tDblSideClustCellFacetOrds  = tChildMesh->get_double_side_interface_cell_pairs_facet_ords(iDI);

                // create a new side cluster for each of the pairs
                std::shared_ptr<xtk::Side_Cluster> tLeftSideCluster  = std::make_shared< Side_Cluster >();
                std::shared_ptr<xtk::Side_Cluster> tRightSideCluster = std::make_shared< Side_Cluster >();

                // paired local child mesh subphase indices
                moris_index tSubphaseIndex0 = tDblSideClustSubphaseCMInds(0);
                moris_index tSubphaseIndex1 = tDblSideClustSubphaseCMInds(1);

                // add enriched interpolation cell
                tLeftSideCluster->mInterpolationCell = tEnrichedCellsOfBaseCell(tSubphaseIndex0);
                tRightSideCluster->mInterpolationCell = tEnrichedCellsOfBaseCell(tSubphaseIndex1);

                // flag both as non-trivial
                tLeftSideCluster->mTrivial = false;
                tRightSideCluster->mTrivial = false;

                // allocate space in integration cell side ordinals
                tLeftSideCluster->mIntegrationCellSideOrdinals = Matrix<IndexMat>(1,tDblSideClustCellCMInds.size());
                tRightSideCluster->mIntegrationCellSideOrdinals = Matrix<IndexMat>(1,tDblSideClustCellCMInds.size());

                // add child meshes to clusters
                tLeftSideCluster->mChildMesh = tChildMesh;
                tRightSideCluster->mChildMesh = tChildMesh;

                tLeftSideCluster->mAssociatedCellCluster = &this->get_cell_cluster(*tLeftSideCluster->mInterpolationCell);
                tRightSideCluster->mAssociatedCellCluster = &this->get_cell_cluster(*tRightSideCluster->mInterpolationCell );

                // add integration cells to cluster
                for(moris::uint iF = 0; iF < tDblSideClustCellCMInds.size(); iF++)
                {
                    // add side ordinals to cluster
                    tLeftSideCluster->mIntegrationCellSideOrdinals(iF) = tDblSideClustCellFacetOrds(iF)(0);
                    tRightSideCluster->mIntegrationCellSideOrdinals(iF) = tDblSideClustCellFacetOrds(iF)(1);

                    // add cells to cluster
                    tLeftSideCluster->mIntegrationCells.push_back(tChildCells(tDblSideClustCellCMInds(iF)(0)));
                    tRightSideCluster->mIntegrationCells.push_back(tChildCells(tDblSideClustCellCMInds(iF)(1)));
                }

                // index of double side set
                moris_index tDoubleSideSetIndex = this->get_dbl_side_set_index(tSubphaseBulkPhases(tSubphaseIndex0),tSubphaseBulkPhases(tSubphaseIndex1));

                // add to a place to store
                mDoubleSideSetsMasterIndex(tDoubleSideSetIndex).push_back(mDoubleSideSingleSideClusters.size());
                mDoubleSideSingleSideClusters.push_back(tLeftSideCluster);
                mDoubleSideSetsSlaveIndex(tDoubleSideSetIndex).push_back(mDoubleSideSingleSideClusters.size());
                mDoubleSideSingleSideClusters.push_back(tRightSideCluster);

                // create double side set
                mtk::Double_Side_Cluster* tDblSideCluster  = new mtk::Double_Side_Cluster(tLeftSideCluster.get(),tRightSideCluster.get(),tChildMesh->get_vertices());

                mDoubleSideClusters.push_back(tDblSideCluster);
                mDoubleSideSets(tDoubleSideSetIndex).push_back(tDblSideCluster);

            }
        }

        // remove from child mesh to not store twice
        tChildMesh->delete_double_sides_interface_sets();
    }

    mListofDoubleSideSets.resize( mDoubleSideSets.size(), nullptr );

    for(moris::uint Ik = 0; Ik<mListofDoubleSideSets.size(); Ik++)
    {
        mListofDoubleSideSets( Ik ) = new moris::mtk::Double_Side_Set(mDoubleSideSetLabels(Ik), this->get_double_side_set_cluster( Ik ), this->get_spatial_dim());
    }
}
//------------------------------------------------------------------------------
moris::Cell<std::string>
Enriched_Integration_Mesh::split_set_name_by_bulk_phase(std::string aBaseName)
{
    moris::uint tNumPhases = mModel->mGeometryEngine.get_num_bulk_phase();
    moris::Cell<std::string> tSetNames(tNumPhases);
    for(moris::uint  i = 0; i<tNumPhases; i++)
    {
        tSetNames(i) = aBaseName+"_p"+std::to_string(i);
    }

    return tSetNames;
}
//------------------------------------------------------------------------------

moris::Cell<std::string>
Enriched_Integration_Mesh::split_set_name_by_child_no_child(std::string aBaseName)
{
    moris::Cell<std::string> tSetNames(2);
    tSetNames(0) = aBaseName+"_c";
    tSetNames(1) = aBaseName+"_n";
    return tSetNames;
}
//------------------------------------------------------------------------------
Cell<moris_index>
Enriched_Integration_Mesh::register_block_set_names_with_cell_topo(moris::Cell<std::string> const & aBlockSetNames,
                                                                   enum CellTopology                aBlockTopology)
{
    uint tNumSetsToRegister = aBlockSetNames.size();

    // block set ords
    Cell<moris_index> tBlockSetOrds(tNumSetsToRegister);

    // iterate and add sets
    for(moris::uint i = 0; i < tNumSetsToRegister; i++)
    {
        tBlockSetOrds(i) = mBlockSetNames.size();

        mBlockSetNames.push_back(aBlockSetNames(i));
        mBlockSetTopology.push_back(aBlockTopology);
        MORIS_ASSERT(mBlockSetLabelToOrd.find(aBlockSetNames(i)) ==  mBlockSetLabelToOrd.end(),"Duplicate block set in mesh");
        mBlockSetLabelToOrd[aBlockSetNames(i)] = tBlockSetOrds(i) ;
    }

    mPrimaryBlockSetClusters.resize(mPrimaryBlockSetClusters.size() + tNumSetsToRegister);
    return tBlockSetOrds;
}
//------------------------------------------------------------------------------
Cell<moris_index>
Enriched_Integration_Mesh::register_side_set_names(moris::Cell<std::string> const & aSideSetNames)
{
    uint tNumSetsToRegister = aSideSetNames.size();

    // block set ords
    Cell<moris_index> tSideSetOrds(tNumSetsToRegister);

    // iterate and add sets
    for(moris::uint i = 0; i < tNumSetsToRegister; i++)
    {
        tSideSetOrds(i) = mSideSetLabels.size();

        mSideSetLabels.push_back(aSideSetNames(i));
        MORIS_ASSERT(mSideSideSetLabelToOrd.find(aSideSetNames(i)) ==  mSideSideSetLabelToOrd.end(),"Duplicate block set in mesh");
        mSideSideSetLabelToOrd[aSideSetNames(i)] = tSideSetOrds(i) ;
    }

    mSideSets.resize(mSideSets.size() + tNumSetsToRegister);

    return tSideSetOrds;
}
//------------------------------------------------------------------------------
Cell<moris_index>
Enriched_Integration_Mesh::register_double_side_set_names(moris::Cell<std::string> const & aDblSideSetNames)
{
    uint tNumSetsToRegister = aDblSideSetNames.size();

    // block set ords
    Cell<moris_index> tDblSideSetOrds(tNumSetsToRegister);

    // iterate and add sets
    for(moris::uint i = 0; i < tNumSetsToRegister; i++)
    {
        tDblSideSetOrds(i) = mDoubleSideSets.size()+i;

        mDoubleSideSetLabels.push_back(aDblSideSetNames(i));
        MORIS_ASSERT(mDoubleSideSetLabelToOrd.find(aDblSideSetNames(i)) ==  mDoubleSideSetLabelToOrd.end(),"Duplicate double side set in mesh");
        mDoubleSideSetLabelToOrd[aDblSideSetNames(i)] = tDblSideSetOrds(i) ;
    }

    mDoubleSideSets.resize(mDoubleSideSets.size() + tNumSetsToRegister);
    mDoubleSideSetsMasterIndex.resize(mDoubleSideSetsMasterIndex.size() + tNumSetsToRegister);
    mDoubleSideSetsSlaveIndex.resize(mDoubleSideSetsSlaveIndex.size() + tNumSetsToRegister);
    return tDblSideSetOrds;
}
//------------------------------------------------------------------------------
void
Enriched_Integration_Mesh::setup_interface_side_sets()
{
    this->declare_interface_side_sets();

    this->create_interface_side_sets_and_clusters();
}

void
Enriched_Integration_Mesh::declare_interface_side_sets()
{
    uint tNumGeometries = mModel->get_geom_engine().get_num_geometries();
    uint tNumBulkPhases = mModel->get_geom_engine().get_num_bulk_phase();

    Cell<std::string> tInterfaceSideNames;
    for(moris::moris_index iG = 0; iG < (moris_index)tNumGeometries; iG++)
    {
        for(moris::moris_index iP0 = 0; iP0 <(moris_index) tNumBulkPhases; iP0++)
        {
            for(moris::moris_index iP1 = 0; iP1 < (moris_index)tNumBulkPhases; iP1++)
            {
                if(iP1 != iP0)
                {
                    std::string tInterfaceSideSetName = get_interface_side_set_name(iG,iP0,iP1);

                    tInterfaceSideNames.push_back(tInterfaceSideSetName);
                }
            }
        }
    }

    register_side_set_names(tInterfaceSideNames);

}
//------------------------------------------------------------------------------
void
Enriched_Integration_Mesh::create_interface_side_sets_and_clusters()
{
    uint tNumGeometries  = mModel->get_geom_engine().get_num_geometries();
    uint tNumBulkPhases  = mModel->get_geom_engine().get_num_bulk_phase();
    uint tNumChildMeshes = mModel->get_cut_mesh().get_num_child_meshes();

    Enriched_Interpolation_Mesh* tEnrInterpMesh  = mModel->mEnrichedInterpMesh(mMeshIndexInModel);
    mtk::Interpolation_Mesh & tBGMesh = mModel->get_background_mesh().get_mesh_data();

    moris_index tMyProcRank = par_rank();

    for(moris::moris_index iG = 0; iG < (moris_index)tNumGeometries; iG++)
    {
        for(moris::moris_index iP0 = 0; iP0 <(moris_index) tNumBulkPhases; iP0++)
        {
            for(moris::moris_index iP1 = 0; iP1 < (moris_index)tNumBulkPhases; iP1++)
            {
                if(iP1 != iP0)
                {
                    std::string tInterfaceSideSetName = get_interface_side_set_name(iG,iP0,iP1);

                    moris_index tSideSetOrd = this->get_side_set_index(tInterfaceSideSetName);

                    // iterate through children meshes and package up interface side sets
                    // note only children meshes need to be considered as interfaces anywhere
                    // else do not make sense
                    for(moris::uint iCM = 0; iCM < tNumChildMeshes; iCM++)
                    {
                        // get child mesh
                        Child_Mesh * tChildMesh = &mModel->get_cut_mesh().get_child_mesh((moris_index)iCM);
                        if(tChildMesh->has_interface_along_geometry(iG) && tBGMesh.get_entity_owner(tChildMesh->get_parent_element_index(),EntityRank::ELEMENT) == tMyProcRank)
                        {
                            // package the interface sides
                            Matrix<IndexMat> tSideSetCellsAndOrds = tChildMesh->pack_interface_sides(iG,iP0,iP1,2);

                            // if there are cells along this boundary keep going
                            // there is chance that a geometry intersects but the boundary between different bulk phases
                            if(tSideSetCellsAndOrds.numel()>0)
                            {
                                // get the base cell
                                moris::mtk::Cell const * tBaseCell = &this->get_mtk_cell(tChildMesh->get_parent_element_index());

                                // get the enriched interpolation cells associated with base cell
                                moris::Cell<xtk::Interpolation_Cell_Unzipped const * > tEnrichedCellsOfBaseCell = tEnrInterpMesh->get_enriched_cells_from_base_cell(tBaseCell);

                                // create a side cluster for each subphase in this child mesh
                                moris::Cell<std::shared_ptr<xtk::Side_Cluster>> tSideClustersForCM(tChildMesh->get_num_subphase_bins());

                                // child cell indices
                                Matrix<IndexMat> const & tChildCellInds = tChildMesh->get_element_inds();

                                // get child cell pointers
                                moris::Cell<moris::mtk::Cell const *> tChildCells = this->get_mtk_cells_loc_inds(tChildCellInds);

                                for(moris::uint  iSP = 0; iSP < tChildMesh->get_num_subphase_bins(); iSP++)
                                {
                                    MORIS_ASSERT(tEnrichedCellsOfBaseCell(iSP)->get_subphase_index() == tChildMesh->get_subphase_indices()(iSP),"Enriched interpolation cell subphases associated with a base cell should be in ascending order.");

                                    tSideClustersForCM(iSP) = std::make_shared< Side_Cluster >();
                                    tSideClustersForCM(iSP)->mInterpolationCell = tEnrichedCellsOfBaseCell(iSP);
                                    tSideClustersForCM(iSP)->mTrivial = false;
                                    tSideClustersForCM(iSP)->mIntegrationCellSideOrdinals = Matrix<IndexMat>(1,tSideSetCellsAndOrds.n_rows());
                                    tSideClustersForCM(iSP)->mChildMesh = tChildMesh;
                                    tSideClustersForCM(iSP)->mAssociatedCellCluster = &this->get_cell_cluster(*tEnrichedCellsOfBaseCell(iSP));

                                }

                                // iterate through child cells on face
                                for(moris::uint iF = 0; iF < tSideSetCellsAndOrds.n_rows(); iF++)
                                {
                                    moris_index tSubphaseGroup      = tChildMesh->get_element_subphase_index(tSideSetCellsAndOrds(iF,0));
                                    moris_index tCMCellIndex        = tSideSetCellsAndOrds(iF,0);
                                    moris_index tIndexInSideCluster = tSideClustersForCM(tSubphaseGroup)->mIntegrationCells.size();

                                    // add information to side cluster
                                    tSideClustersForCM(tSubphaseGroup)->mIntegrationCellSideOrdinals(tIndexInSideCluster) = tSideSetCellsAndOrds(iF,1);
                                    tSideClustersForCM(tSubphaseGroup)->mIntegrationCells.push_back(tChildCells(tCMCellIndex));
                                }

                                // iterate through, get rid of extra space in side ordinals and add to side set clusters data
                                for(moris::uint  iSP = 0; iSP < tChildMesh->get_num_subphase_bins(); iSP++)
                                {
                                    // only add this side cluster to the set if it has at least one integration cell in it
                                    if(tSideClustersForCM(iSP)->mIntegrationCells.size() > 0 )
                                    {

                                        tSideClustersForCM(iSP)->mIntegrationCellSideOrdinals.resize(1,tSideClustersForCM(iSP)->mIntegrationCells.size());

                                        tSideClustersForCM(iSP)->finalize_setup();

                                        // add
                                        mSideSets(tSideSetOrd).push_back(tSideClustersForCM(iSP));
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }

    uint tCurrentSize = mListofSideSets.size();
    mListofSideSets.resize( mSideSets.size(), nullptr );

    for(moris::uint Ik = tCurrentSize; Ik<mListofSideSets.size(); Ik++)
    {
        mListofSideSets( Ik ) = new moris::mtk::Side_Set( mSideSetLabels(Ik), this->get_side_set_cluster( Ik ), this->get_spatial_dim() );
    }

}
//------------------------------------------------------------------------------
bool
Enriched_Integration_Mesh::field_exists(std::string              aLabel,
                                        enum moris::EntityRank   aEntityRank)
{
    moris::moris_index tIndex = this->get_entity_rank_field_index(aEntityRank);
    return mFieldLabelToIndex(tIndex).find(aLabel) != mFieldLabelToIndex(tIndex).end();

}
//------------------------------------------------------------------------------
moris_index
Enriched_Integration_Mesh::get_entity_rank_field_index(enum moris::EntityRank   aEntityRank)
{
    MORIS_ERROR(aEntityRank == EntityRank::NODE || aEntityRank == EntityRank::ELEMENT,"Only node and cell fields are supported");

    moris_index tIndex = MORIS_INDEX_MAX;
    if(aEntityRank == EntityRank::NODE)
    {
        tIndex = 0;
    }

    else if(aEntityRank == EntityRank::ELEMENT)
    {
        tIndex = 1;
    }

    return tIndex;

}
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------

}


