/*
 * cl_MTK_Intersection_Mesh.cpp
 *
 *  Created on: Aug 20, 2021
 *      Author: momo
 */

#include "cl_MTK_Intersection_Mesh.hpp"

namespace moris
{
    namespace mtk
    {

        Intersection_Mesh::Intersection_Mesh( moris::mtk::Integration_Mesh* aBackGroundMesh,
                moris::mtk::Intersection_Detect* aIntersectionDetect ):mBackGroundMesh(aBackGroundMesh),
                        mIntersectionDetect(aIntersectionDetect)
                        {}

        // ----------------------------------------------------------------------------

        moris::uint
        Intersection_Mesh::get_num_entities(enum EntityRank aEntityRank, const moris_index aIndex ) const
        {
            //number of vertices
            if(aEntityRank == EntityRank::NODE)
            {
                moris::size_t tNumBackgroundEntities = mBackGroundMesh->get_num_entities((moris::EntityRank)aEntityRank);
                moris::size_t tExternalEntities = mIntersectionDetect->mIntersectedMeshData.get_num_entities_external_data(aEntityRank);

                return tNumBackgroundEntities + tExternalEntities;
            }
            //number of cells
            else if(aEntityRank == EntityRank::ELEMENT)
            {
                moris::size_t tNumBackgroundEntities = mBackGroundMesh->get_num_entities((moris::EntityRank)aEntityRank);
                moris::size_t tExternalEntities = mIntersectionDetect->mIntersectedMeshData.get_num_entities_external_data(aEntityRank);

                return tNumBackgroundEntities + tExternalEntities;
            }
            else
            {
                MORIS_ERROR(0,"Only cells and verts supported.");
                return 0;
            }
        }

        // ----------------------------------------------------------------------------
        moris::Cell< std::string >
        Intersection_Mesh::get_set_names(enum EntityRank aSetEntityRank) const
        {
            //get the current list of sets
            moris::Cell< std::string > tSetNames = mBackGroundMesh->get_set_names(aSetEntityRank);

            //Distinguish between side sets and block sets
            if ( aSetEntityRank == EntityRank::FACE or aSetEntityRank == EntityRank::EDGE )
            {
                tSetNames.append(mIntersectionDetect->mSideSetLabels);
            }

            if( aSetEntityRank == EntityRank::ELEMENT)
            {
                tSetNames.append(mIntersectionDetect->mBlockSetLabels);
            }

            return tSetNames;
        }
        // ----------------------------------------------------------------------------

        bool
        Intersection_Mesh::is_intersection_side_set(const  std::string & aSetName) const
        {
            //iterate through the side set labels for intersection data
            auto tIter = mIntersectionDetect->mSideSideSetLabelToOrd.find(aSetName);

            //if not found
            if ( tIter == mIntersectionDetect->mSideSideSetLabelToOrd.end() )
            {
                return false;
            }

            //if found
            return true;
        }

        // ----------------------------------------------------------------------------

        bool
        Intersection_Mesh::is_intersection_block_set(const  std::string & aSetName) const
        {
            //iterate through the side set labels for intersection data
            auto tIter = mIntersectionDetect->mBlockSetLabelToOrd.find(aSetName);

            //if not found
            if ( tIter == mIntersectionDetect->mBlockSetLabelToOrd.end() )
            {
                return false;
            }

            //if found
            return true;
        }

        // ----------------------------------------------------------------------------

        void
        Intersection_Mesh::get_sideset_elems_loc_inds_and_ords(
                const  std::string & aSetName,
                Matrix< IndexMat > & aElemIndices,
                Matrix< IndexMat > & aSidesetOrdinals ) const
        {

            //Determine if it is a side set created in the intersection process
            if( this->is_intersection_side_set( aSetName ) )
            {
                //find the side set  and its index
                auto tIter = mIntersectionDetect->mSideSideSetLabelToOrd.find(aSetName);
                moris_index tSideSetIndex = tIter->second;

                // if the side set is NOT empty
                if( mIntersectionDetect->mMasterSideSets(tSideSetIndex).size() > 0)
                {

                    //access the clusters in the side set
                    moris::Cell<mtk::Cluster const *> tSideClusters = mIntersectionDetect->mMasterSideSets(tSideSetIndex);

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

                //if the side set is empty return an empty matrices
                else
                {
                    aElemIndices.resize(1,0);
                    aSidesetOrdinals.resize(1,0);
                }

            }

            //if the side set belongs to background mesh, let background mesh do the operation
            else
            {
                mBackGroundMesh->get_sideset_elems_loc_inds_and_ords( aSetName, aElemIndices, aSidesetOrdinals);
            }

        }

        // ----------------------------------------------------------------------------

        Matrix<IndexMat>
        Intersection_Mesh::get_element_indices_in_block_set(
                uint aSetIndex)
        {

            //Determine if it is a block set created in the intersection process
            if( this->is_intersection_block_set( this->get_set_names(EntityRank::ELEMENT)(aSetIndex) ) )
            {
                //If the block set is NOT empty
                if( mIntersectionDetect->mMasterSideCells.size() > 0)
                {
                    //Initialize the output
                    Matrix<IndexMat> tIndexMatrix(1, mIntersectionDetect->mMasterSideCells.size()) ;

                    //Iterate through the cells and obtain their index
                    for(uint i = 0 ; i < mIntersectionDetect->mMasterSideCells.size(); i++ )
                    {
                        tIndexMatrix(i) = mIntersectionDetect->mMasterSideCells(i)->get_index();
                    }

                    return tIndexMatrix;
                }

                //If the block set is empty return an empty matrix
                else
                {
                    return Matrix<IndexMat>(0,0);
                }

            }

            //If the block set belongs to background mesh, let background mesh do the operation
            else
            {
                return mBackGroundMesh->get_element_indices_in_block_set(aSetIndex) ;
            }

        }


        // ----------------------------------------------------------------------------

        Matrix<IndexMat>
        Intersection_Mesh::get_element_ids_in_block_set(
                uint aSetIndex) const
        {
            //Determine if it is a block set created in the intersection process
            if( this->is_intersection_block_set( this->get_set_names(EntityRank::ELEMENT)(aSetIndex) ) )
            {
                //If the block set is NOT empty
                if( mIntersectionDetect->mMasterSideCells.size() > 0)
                {
                    //Initialize the output
                    Matrix<IndexMat> tIndexMatrix(1, mIntersectionDetect->mMasterSideCells.size()) ;

                    //Iterate through the cells and obtain their index
                    for(uint i = 0 ; i < mIntersectionDetect->mMasterSideCells.size(); i++ )
                    {
                        tIndexMatrix(i) = mIntersectionDetect->mMasterSideCells(i)->get_id();
                    }

                    return tIndexMatrix;
                }

                //If the block set is empty return an empty matrix
                else
                {
                    std::cout<<"non"<<std::endl;
                    return Matrix<IndexMat>(0,0);
                }
            }

            //If the block set belongs to background mesh, let background mesh do the operation
            else
            {
                return mBackGroundMesh->get_element_ids_in_block_set(aSetIndex) ;
            }

        }

        // ----------------------------------------------------------------------------

        moris::Matrix< moris::DDRMat >
        Intersection_Mesh::get_node_coordinate(moris_index aNodeIndex) const
        {
            //Check if the node is added in the intersection process
            if( mIntersectionDetect->mIntersectedMeshData.is_external_entity(aNodeIndex, EntityRank::NODE) )
            {
                return mIntersectionDetect->mNewNodeCoords(mIntersectionDetect->mNodeIndexToCoordsMap[aNodeIndex]);
            }
            //If not an intersection node, refer it to background mesh
            else
            {
                return mBackGroundMesh->get_node_coordinate(aNodeIndex);
            }
        }

        // ----------------------------------------------------------------------------

        moris::moris_id
        Intersection_Mesh::get_glb_entity_id_from_entity_loc_index(moris_index   aEntityIndex,
                enum EntityRank aEntityRank,
                const moris_index aMeshIndex ) const
        {
            if( mIntersectionDetect->mIntersectedMeshData.is_external_entity( aEntityIndex, aEntityRank) )
            {
                return mIntersectionDetect->mEntityLocaltoGlobalMap((uint)aEntityRank)(aEntityIndex);
            }
            else
            {
                return mBackGroundMesh->get_glb_entity_id_from_entity_loc_index(aEntityIndex,aEntityRank );
            }
        }

        // ----------------------------------------------------------------------------

        enum CellTopology
        Intersection_Mesh::get_blockset_topology(const std::string & aSetName)
        {
            //All the cells created in intersection process are TET elements
            if( this->is_intersection_block_set( aSetName ) )
            {
                return CellTopology::TET4;
            }

            return mBackGroundMesh->get_blockset_topology(aSetName);
        }

        // ----------------------------------------------------------------------------

        Matrix< IndexMat >
        Intersection_Mesh::get_nodes_connected_to_element_loc_inds(moris_index aElementIndex) const
        {
            // If the element is an intersection element
            if( mIntersectionDetect->mIntersectedMeshData.is_external_entity(aElementIndex, EntityRank::ELEMENT) )
            {
                // Get the corresponding cell from list of the cells
                // based on the map from cell index to mtk::cell
                moris::mtk::Cell const * tCell = mIntersectionDetect->mMasterSideCells(mIntersectionDetect->mMasterCellIndextoCellMap[aElementIndex]);

                // Get vertices' indices attached to that cell
                return tCell->get_vertex_inds();

            }

            // If it is a background element
            else
            {
                return mBackGroundMesh->get_nodes_connected_to_element_loc_inds(aElementIndex);
            }

        }

    } /* end namespace mtk */

} /* end namespace moris */

