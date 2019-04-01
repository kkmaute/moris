/*
 * cl_XTK_Ghost_Stabilization.cpp
 *
 *  Created on: Mar 26, 2019
 *      Author: doble
 */

#include "cl_XTK_Ghost_Stabilization.hpp"

#include "fn_generate_shared_face_element_graph.hpp"

namespace xtk
{
    void
    Ghost_Stabilization::setup_ghost_stabilization_facets()
    {
        Ghost_Setup_Data tGhostSetupData;

        // identify candidate facets for ghost penalization
        identify_candidate_facets_for_ghost_penalization(tGhostSetupData);

        // iterate through bulk phases and construct ghost cells
        moris_index tNumBulkPhases = mXTKModel->get_geom_engine().get_num_bulk_phase();

        // resize member data
        mGhostCells.resize(tNumBulkPhases);

        for(moris_index  iBP = 0; iBP <tNumBulkPhases; iBP++)
        {
            // construct the ghost cells for this phase
            this->construct_ghost_cells(iBP,tGhostSetupData);
        }


    }


    void
    Ghost_Stabilization::identify_candidate_facets_for_ghost_penalization(Ghost_Setup_Data & aGhostSetupData)
    {
        Cut_Mesh &               tCutMesh    = mXTKModel->get_cut_mesh();
        Background_Mesh &        tBMMesh     = mXTKModel->get_background_mesh();
        mtk::Mesh const & tBMMeshData = tBMMesh.get_mesh_data();

        // iterate through children meshes
        for(uint iCM = 0 ; iCM <tCutMesh.get_num_child_meshes(); iCM++ )
        {
            // get parent element index of child mesh i
            moris_index tParentCellIndex = tCutMesh.get_parent_element_index(iCM);

            // get facets attached to this element
            Matrix<IndexMat> tElementFaces = tBMMeshData.get_entity_connected_to_entity_loc_inds(tParentCellIndex,EntityRank::ELEMENT,EntityRank::FACE);

            // add faces to candidate facets
            for(uint iF = 0; iF<tElementFaces.numel(); iF++)
            {
                aGhostSetupData.mCandidateFacets.push_back(tElementFaces(iF));
            }
        }

        // Make vector unique
        unique(aGhostSetupData.mCandidateFacets);
    }

    void
    Ghost_Stabilization::construct_ghost_cells( moris_index        aBulkPhase,
                                                Ghost_Setup_Data & aGhostSetupData)
    {
        // number of candidate facets
        moris_index tNumCandidateFacets = aGhostSetupData.mCandidateFacets.size();

        // iterate through candidate facets
        for(moris_index iF = 0; iF < tNumCandidateFacets; iF++)
        {
            // candidate facets
            moris_index tFacetIndex = aGhostSetupData.mCandidateFacets(iF);

            // figure out if we need to penalize this facet
            Cell<moris_index>      tCellsWithChildren;
            Cell<moris_index>      tCellsWithNoChildren;
            Cell<Matrix<IndexMat>> tChildMeshCellIndsOnFacet; /* local index to the specific child mesh*/
            Cell<Matrix<IndexMat>> tChildMeshCellOnFacetOrdinals;

            bool tPenalize = this->decide_whether_to_penalize_facet(aBulkPhase,tFacetIndex,tCellsWithChildren,tCellsWithNoChildren,tChildMeshCellIndsOnFacet,tChildMeshCellOnFacetOrdinals);

            // if we need to penalize this facet, then construct the ghost cell
            if(tPenalize)
            {
                this->construct_facet_ghost_cell(aBulkPhase,tFacetIndex,tCellsWithChildren,tCellsWithNoChildren,tChildMeshCellIndsOnFacet,tChildMeshCellOnFacetOrdinals);
            }



        }


    }

    bool
    Ghost_Stabilization::decide_whether_to_penalize_facet( moris_index aBulkPhase,
                                                           moris_index aFacetIndex,
                                                           Cell<moris_index>      & aCellsWithChildren,
                                                           Cell<moris_index>      & aCellsWithNoChildren,
                                                           Cell<Matrix<IndexMat>> & aChildMeshCellIndsOnFacet,
                                                           Cell<Matrix<IndexMat>> & aChildMeshCellOnFacetOrdinals)
    {
        // easy access member data
        Cut_Mesh &               tCutMesh    = mXTKModel->get_cut_mesh();
        Background_Mesh &        tBMMesh     = mXTKModel->get_background_mesh();
        mtk::Mesh const & tBMMeshData = tBMMesh.get_mesh_data();

        // Cells attached to facets
        Matrix<IndexMat> tFacetToCell = tBMMeshData.get_entity_connected_to_entity_loc_inds(aFacetIndex,EntityRank::FACE,EntityRank::ELEMENT);

        // If there are two cells attached to a face then we need to penalize the jump.
        // if there is only one then the jump operation doesn't really make sense.

        // flag indicating whether to penalize or not
        bool tPenalizeFacet = false;
        if(tFacetToCell.numel()>1)
        {
            MORIS_ASSERT(tFacetToCell.numel() == 2,"This section needs extra work to consider a facet shared by more than 2 elements although it has been written without this inherent assumption");

            aCellsWithChildren   = Cell<moris_index>();
            aCellsWithNoChildren = Cell<moris_index>();

            // iterate through the cells attached to the facet
            // to figure out which one has children and which do not
            for(uint iEl = 0; iEl<tFacetToCell.numel(); iEl++)
            {
                // if this cell has children then add it to the
                // cell tracking entities with children
                if(tBMMesh.entity_has_children(tFacetToCell(iEl),EntityRank::ELEMENT))
                {
                    aCellsWithChildren.push_back(tFacetToCell(iEl));
                }
                // if this cell does not have children add it to the other vector
                else
                {
                    aCellsWithNoChildren.push_back(tFacetToCell(iEl));
                }
            }


            tPenalizeFacet = true;

            // iterate through cells without children
            for(uint iWOChild = 0; iWOChild<aCellsWithNoChildren.size(); iWOChild++)
            {
                // get element bulk phase
                moris_index tElemBulkPhase = tBMMesh.get_element_phase_index(aCellsWithNoChildren(iWOChild));

                // if this element is not the current bulk phase we are constructing ghost for
                // then it does not need to be penalized
                if( tElemBulkPhase != aBulkPhase)
                {
                    tPenalizeFacet = false;
                }

            }

            aChildMeshCellIndsOnFacet     = Cell<Matrix<IndexMat>>(aCellsWithChildren.size());
            aChildMeshCellOnFacetOrdinals = Cell<Matrix<IndexMat>>(aCellsWithChildren.size());

            // If the parent element says to penalize the facet then there will
            // definitely be a ghost here. If it says not to penalize, then
            // there will not be a ghost, this is why we check the cells
            // with no children first.
            if(tPenalizeFacet)
            {
                // iterate through cells with children
                for(uint iWChild = 0; iWChild<aCellsWithChildren.size(); iWChild++)
                {
                    // get child mesh
                    moris_index tCMIndex    = tBMMesh.child_mesh_index(aCellsWithChildren(iWChild),EntityRank::ELEMENT);
                    Child_Mesh & tChildMesh = tCutMesh.get_child_mesh(tCMIndex);

                    // get elements attached to parent face
                    Matrix< IdMat > tChildElemIdsOnFace;
                    tChildMesh.get_child_elements_connected_to_parent_face(aFacetIndex,tChildElemIdsOnFace,aChildMeshCellIndsOnFacet(iWChild),aChildMeshCellOnFacetOrdinals(iWChild));


                    // reference the child mesh element indexs
                    Matrix<IndexMat> const &  tElementPhases = tChildMesh.get_element_phase_indices();

                    // iterate through child elements on the face and check if they have the bulk phase of interest here
                    for(uint iCEl = 0; iCEl<aChildMeshCellIndsOnFacet(iWChild).numel(); iCEl++)
                    {
                        moris_index tCMElemIndex = aChildMeshCellIndsOnFacet(iWChild)(iCEl);

                        if(tElementPhases(tCMElemIndex) == aBulkPhase )
                        {
                            tPenalizeFacet = true;
                            break;
                        }

                        else
                        {
                            tPenalizeFacet = false;
                        }
                    }
                }

            }
        }

        return tPenalizeFacet;
    }

    void
    Ghost_Stabilization::construct_facet_ghost_cell(moris_index              aBulkPhase,
                                                    moris_index              aFacetIndex,
                                                    Cell<moris_index>      & aCellsWithChildren,
                                                    Cell<moris_index>      & aCellsWithNoChildren,
                                                    Cell<Matrix<IndexMat>> & aChildMeshCellIndsOnFacet,
                                                    Cell<Matrix<IndexMat>> & aChildMeshCellOnFacetOrdinals)
    {
        //TODO: link the correct pairs
        MORIS_ASSERT(aCellsWithNoChildren.size() == 1 || aCellsWithNoChildren.size() == 0," Currently only tested for a single unintersected element in the configuration");
        MORIS_ASSERT(aCellsWithChildren.size() >= 1 ," There needs to be at least one cell that has a child mesh in it");

        // easy access
        Cut_Mesh &        tCutMesh    = mXTKModel->get_cut_mesh();
        Background_Mesh & tBMMesh     = mXTKModel->get_background_mesh();
        mtk::Mesh const & tBMMeshData = tBMMesh.get_mesh_data();

        if(aCellsWithChildren.size() == 2)
        {
            // child mesh indices
            moris_index tChildMeshIndex0 = tBMMesh.child_mesh_index(aCellsWithChildren(0),EntityRank::ELEMENT);
            moris_index tChildMeshIndex1 = tBMMesh.child_mesh_index(aCellsWithChildren(1),EntityRank::ELEMENT);


            Matrix< IndexMat > tBoundaryPairsOrdinals;
            Matrix< IndexMat > tBoundaryCellPairs = generate_shared_face_element_pairs(aFacetIndex,tChildMeshIndex0,tChildMeshIndex1,tCutMesh,tBoundaryPairsOrdinals);

            // child meshes
            Child_Mesh & tCM0 = tCutMesh.get_child_mesh(tChildMeshIndex0);
            Child_Mesh & tCM1 = tCutMesh.get_child_mesh(tChildMeshIndex1);

            // child mesh indices
            Matrix< IndexMat > const & tChildMeshElementInds0 = tCM0.get_element_inds();
            Matrix< IndexMat > const & tChildMeshElementInds1 = tCM1.get_element_inds();



            // for each pair construct a ghost
            for(uint i = 0 ; i <tBoundaryCellPairs.n_cols(); i++)
            {
                // get mtk cell pointer for first element
                moris_index tCellIndex0 = tChildMeshElementInds0(tBoundaryCellPairs(0,i));
                moris_index tCellIndex1 = tChildMeshElementInds1(tBoundaryCellPairs(1,i));

                // create ghost cell
                mGhostCells(aBulkPhase).push_back(Ghost_Cell( aFacetIndex,
                                                 &tBMMesh.get_mtk_cell(tCellIndex0),
                                                 tBoundaryPairsOrdinals(0,i),
                                                 &tBMMesh.get_mtk_cell(tCellIndex1),
                                                 tBoundaryPairsOrdinals(1,i)));

            }
        }

        // tets on hex interface (4-tris to 1 quad)
        else if(aCellsWithChildren.size() == 1 && aCellsWithNoChildren.size() == 1)
        {
            // in this case we assume all children elements on this boundary are the left and the unintersected is a
            moris_index tParentCellIndex = aCellsWithNoChildren(0);

            // Side ordinal of parent
            moris_index tParentOrdinal  = tBMMeshData.get_facet_ordinal_from_cell_and_facet_loc_inds(aFacetIndex,tParentCellIndex);

            // get the element inds
            moris_index tChildMeshIndex = tBMMesh.child_mesh_index(aCellsWithChildren(0),EntityRank::ELEMENT);
            Child_Mesh & tCM0 = tCutMesh.get_child_mesh(tChildMeshIndex);
            Matrix< IndexMat > const & tChildMeshElementInds = tCM0.get_element_inds();

            // iterate through children elements
            moris_index tNumChildrenOnFace = aChildMeshCellIndsOnFacet(0).numel();
            for(moris_index i = 0; i <tNumChildrenOnFace; i++)
            {
                // create ghost cell
                mGhostCells(aBulkPhase).push_back(Ghost_Cell( aFacetIndex,
                                                 &tBMMesh.get_mtk_cell(tParentCellIndex),
                                                 tParentOrdinal,
                                                 &tBMMesh.get_mtk_cell(tChildMeshElementInds(aChildMeshCellIndsOnFacet(0)(i))),
                                                 aChildMeshCellOnFacetOrdinals(0)(i)));
            }

        }
        else
        {
            MORIS_ASSERT(0,"Case not implemented in ghost");
        }

    }


}
