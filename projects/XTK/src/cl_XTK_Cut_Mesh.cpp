/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_XTK_Cut_Mesh.cpp
 *
 */

#include "cl_XTK_Cut_Mesh.hpp"
#include "cl_Matrix.hpp"
#include "cl_XTK_Model.hpp"
#include "cl_MPI_Tools.hpp"
#include "cl_TOL_Memory_Map.hpp"
namespace xtk
{
    Cut_Mesh:: Cut_Mesh() :
                    mNumberOfChildrenMesh(0),
                    mChildrenMeshes(0),
                    mConsistentCounts(false),
                    mNumEntities(4,0),
                    mChildElementTopo(mtk::CellTopology::TET4){}
    // ----------------------------------------------------------------------------------
    Cut_Mesh::Cut_Mesh(
            Model* aModel,
            moris::uint aModelDim) :
                    mModel(aModel),
                    mSpatialDim(aModelDim),
                    mNumberOfChildrenMesh(0),
                    mChildrenMeshes(0),
                    mConsistentCounts(false),
                    mNumEntities(4,0),
                    mChildElementTopo(mtk::CellTopology::TET4){}
    // ----------------------------------------------------------------------------------
    Cut_Mesh::Cut_Mesh(
            Model* aModel,
            moris::size_t aNumSimpleMesh,
            moris::size_t aModelDim) :
                    mModel(aModel),
                    mSpatialDim(aModelDim),
                    mNumberOfChildrenMesh(aNumSimpleMesh),
                    mChildrenMeshes(aNumSimpleMesh),
                    mConsistentCounts(false),
                    mNumEntities(4,0),
                    mChildElementTopo(mtk::CellTopology::TET4)
    {
        for(moris::size_t i = 0; i <aNumSimpleMesh; i++)
        {
            mChildrenMeshes(i) = new Child_Mesh();
        }
    }

    Cut_Mesh::~Cut_Mesh()
    {
        for(moris::size_t i = 0; i <mChildrenMeshes.size(); i++)
        {
            delete mChildrenMeshes(i);
        }
    }

    // ----------------------------------------------------------------------------------
    void
    Cut_Mesh::inititalize_new_child_meshes(
            moris::size_t aNumNewChildMesh,
            moris::size_t aModelDim)
    {
        mNumberOfChildrenMesh = aNumNewChildMesh + mNumberOfChildrenMesh;
        mChildrenMeshes.resize(mNumberOfChildrenMesh,nullptr);
    }
    // ----------------------------------------------------------------------------------
    void
    Cut_Mesh::generate_templated_mesh(
            moris::size_t     aChildMeshIndex,
            enum TemplateType aTemplate)
    {
        MORIS_ASSERT(aChildMeshIndex < mNumberOfChildrenMesh, "The requested mesh index is out of bounds. Consider allocating more space");
        mChildrenMeshes(aChildMeshIndex)->modify_child_mesh(aTemplate);

        // Meshes have changed so counts are wrong
        mConsistentCounts = false;
    }
    // ----------------------------------------------------------------------------------
    void
    Cut_Mesh::generate_templated_mesh(
            enum TemplateType aTemplate)
    {
        for (moris::size_t i = 0; i < mNumberOfChildrenMesh; i++)
        {
            mChildrenMeshes(i)->modify_child_mesh(aTemplate);
        }

        // Meshes have changed so counts are wrong
        mConsistentCounts = false;
    }
    // ----------------------------------------------------------------------------------
    void
    Cut_Mesh::generate_templated_mesh(
            Matrix< IndexMat > const & aChildMeshIndices,
            enum TemplateType          aTemplate)
    {
        for (moris::size_t i = 0; i < aChildMeshIndices.n_cols(); i++)
        {
            mChildrenMeshes(aChildMeshIndices(0,i))->modify_child_mesh(aTemplate);
        }

        // Meshes have changed so counts are wrong
        mConsistentCounts = false;
    }
    // ----------------------------------------------------------------------------------
    void
    Cut_Mesh::convert_cut_mesh_to_tet10s()
    {
        for(moris::size_t i = 0; i<get_num_child_meshes(); i++)
        {
            mChildrenMeshes(i)->convert_tet4_to_tet10_child();
        }

        mChildElementTopo = mtk::CellTopology::TET10;
    }
    // ----------------------------------------------------------------------------------
    void
    Cut_Mesh::set_child_mesh( moris_index aCMIndex,
                              Child_Mesh* aChildMesh)
    {
        MORIS_ASSERT(aCMIndex < (moris_index)mNumberOfChildrenMesh, "The requested mesh index is out of bounds.");
        MORIS_ASSERT(mChildrenMeshes(aCMIndex) == nullptr,"Overwriting existing child mesh not allowed");

        mChildrenMeshes(aCMIndex) = aChildMesh;
    }

    // ----------------------------------------------------------------------------------
    void
    Cut_Mesh::initialize_new_mesh_from_parent_element(
            moris::size_t              aChildMeshIndex,
            enum TemplateType          aTemplate,
            Matrix< IndexMat >       & aNodeIndices,
            Cell<Matrix< IndexMat >> & aParentEntities)
    {
        MORIS_ASSERT(aChildMeshIndex < mNumberOfChildrenMesh, "The requested mesh index is out of bounds.");

        // Construct a template and initialize this new mesh with the template

        // Set all node parent ranks to mtk::EntityRank::NODE which = 0;
        moris::Matrix< moris::DDSTMat > tElementNodeParentRanks(1,8);
        tElementNodeParentRanks.fill(0);

        // Set all edge parent ranks to mtk::EntityRank::EDGE which = 1;
        moris::Matrix< moris::DDSTMat > tParentEdgeRanks(1,aParentEntities(1).numel());

        // This is needed for HMR because edges in 2-d have rank 2
        if(mModel->get_background_mesh().get_mesh_type() == mtk::MeshType::HMR && mSpatialDim == 2)
        {
            tParentEdgeRanks.fill(2);
        }
        else
        {
            tParentEdgeRanks.fill(1);
        }

        // Set all node parent ranks to mtk::EntityRank::FACE which = 2;
        moris::Matrix< moris::DDSTMat > tParentFaceRanks(1,aParentEntities(2).numel());
        tParentFaceRanks.fill(2);

        // No interface sides
        moris::Matrix< moris::DDSTMat > tInterfaceSides(1,1);
        tInterfaceSides.fill(std::numeric_limits<moris::size_t>::max());

        // Note for this: child mesh node indices, parent node indices, and element to node connectivity are
        // the same thing. This is why aNodeIndices appears 3 times in this call.
        mChildrenMeshes(aChildMeshIndex) = new Child_Mesh(
                mSpatialDim,
                aParentEntities(3)(0,0),
                aNodeIndices,
                aNodeIndices,
                tElementNodeParentRanks,
                aNodeIndices,
                aParentEntities(1),
                tParentEdgeRanks,
                aParentEntities(2),
                tParentFaceRanks,
                tInterfaceSides);

        if(mModel->get_background_mesh().get_mesh_type() == mtk::MeshType::HMR)
        {
            mChildrenMeshes(aChildMeshIndex)->mark_as_hmr_child_mesh();
        }

        // set parent element parametric coordinate
        switch(aTemplate)
        {
            case TemplateType::HEX_8:
            {
                const moris::Matrix< moris::DDRMat > tParamCoords(
                    {{-1.0, -1.0, -1.0},
                    { 1.0, -1.0, -1.0},
                    { 1.0,  1.0, -1.0},
                    {-1.0,  1.0, -1.0},
                    {-1.0, -1.0,  1.0},
                    { 1.0, -1.0,  1.0},
                    { 1.0,  1.0,  1.0},
                    {-1.0,  1.0,  1.0}});

                // Add hex 8 parametric coordinates
                mChildrenMeshes(aChildMeshIndex)->allocate_parametric_coordinates(8,3);
                mChildrenMeshes(aChildMeshIndex)->add_node_parametric_coordinate(aNodeIndices,tParamCoords);

                break;
            }
            case TemplateType::TET_4:
            {
                const moris::Matrix< moris::DDRMat > tParamCoords(
                    {{ 1.0, 0.0, 0.0, 0.0},
                    { 0.0, 1.0, 0.0, 0.0},
                    { 0.0, 0.0, 1.0, 0.0},
                    { 0.0, 0.0, 0.0, 1.0}});

                // Add tetra parametric coordinates
                mChildrenMeshes(aChildMeshIndex)->allocate_parametric_coordinates(4,4);
                mChildrenMeshes(aChildMeshIndex)->add_node_parametric_coordinate(aNodeIndices,tParamCoords);
                break;
            }
            case TemplateType::QUAD_4:
            {
                const moris::Matrix< moris::DDRMat > tParamCoords = {
                        {-1.0, -1.0},
                        { 1.0, -1.0},
                        { 1.0,  1.0},
                        {-1.0,  1.0}};

                // Add quad 4 parametric coordinates
                mChildrenMeshes(aChildMeshIndex)->allocate_parametric_coordinates(4, 2);
                mChildrenMeshes(aChildMeshIndex)->add_node_parametric_coordinate(aNodeIndices,tParamCoords);

                break;
            }

            default:
            {
                MORIS_ERROR(0,"Parent element template type not recognized");
            }
        }
    }
    // ----------------------------------------------------------------------------------
    void
    Cut_Mesh::modify_templated_mesh(
            moris::size_t     aChildMeshIndex,
            enum TemplateType aTemplate)
    {
        MORIS_ASSERT(aChildMeshIndex < mNumberOfChildrenMesh, "The requested mesh index is out of bounds.");
        mChildrenMeshes(aChildMeshIndex)->modify_child_mesh(aTemplate);

        // Meshes have changed so counts are wrong
        mConsistentCounts = false;
    }
    // ----------------------------------------------------------------------------------
    void
    Cut_Mesh::add_entity_to_intersect_connectivity(
            moris::size_t aChildMeshIndex,
            moris::size_t aDPrime1Ind,
            moris::size_t aDPrime2Ind,
            moris::size_t aReturnType)
    {
        MORIS_ASSERT(aChildMeshIndex < mNumberOfChildrenMesh, "The requested mesh index is out of bounds.");
        mChildrenMeshes(aChildMeshIndex)->add_entity_to_intersect_connectivity(aDPrime1Ind, aDPrime2Ind, aReturnType);
    }
    // ----------------------------------------------------------------------------------
    void
    Cut_Mesh::add_interface_element(Interface_Element const & aInterfaceElement)
    {
        mInterfaceElements.push_back(aInterfaceElement);
    }
    // ----------------------------------------------------------------------------------
    moris::Cell<Interface_Element> &
    Cut_Mesh::get_interface_elements()
    {
        return mInterfaceElements;
    }
    // ----------------------------------------------------------------------------------
    moris::Matrix<moris::IdMat>
    Cut_Mesh::get_interface_element_ids()
    {
        moris::Matrix<moris::IdMat> tInterfaceElementIds(1,mInterfaceElements.size());

        moris::moris_index tMyProcRank = par_rank();

        moris::uint tCount = 0;

        for(moris::uint i = 0; i <mInterfaceElements.size(); i++)
        {
            if(mInterfaceElements(i).get_element_owner() == tMyProcRank)
            {
                tInterfaceElementIds(i) = mInterfaceElements(i).get_element_id();
                tCount++;
            }
        }

        tInterfaceElementIds.resize(1,tCount);

        return tInterfaceElementIds;
    }
    // ----------------------------------------------------------------------------------
    moris::Matrix<moris::IndexMat>
    Cut_Mesh::get_extracted_interface_elements_loc_inds()
    {
        MORIS_ASSERT(mChildElementTopo == mtk::CellTopology::TET4,"Interface unzipping only tested on child meshes with tet4s");

        // hardcoded for tet4s
        moris::uint tNumNodesPerElement = 6;

        moris::uint tNumInterfaceElements = mInterfaceElements.size();

        moris::Matrix<moris::IndexMat> tInterfaceElemToNode(tNumInterfaceElements,tNumNodesPerElement);

        moris::moris_index tMyProcRank = par_rank();

        moris::uint tCount = 0;

        for(moris::uint iInt =0; iInt < tNumInterfaceElements; iInt++)
        {
            if(mInterfaceElements(iInt).get_element_owner() == tMyProcRank)
            {
                tInterfaceElemToNode.get_row(iInt) = mInterfaceElements(iInt).extract_as_standard_element_loc_inds().get_row(0);
                tCount++;
            }
        }

        tInterfaceElemToNode.resize(tCount,tNumNodesPerElement);

        return tInterfaceElemToNode;
    }
    // ----------------------------------------------------------------------------------
    void
    Cut_Mesh::set_node_index(
            moris::size_t const & aChildMeshIndex,
            Matrix< IndexMat >  & aNodeInd)
    {
        MORIS_ASSERT(aChildMeshIndex < mNumberOfChildrenMesh, "The requested mesh index is out of bounds.");
        mChildrenMeshes(aChildMeshIndex)->add_node_indices(aNodeInd);

        // Meshes have changed so counts are wrong
        mConsistentCounts = false;
    }
    // ----------------------------------------------------------------------------------
    void
    Cut_Mesh::set_node_ids(
            moris::size_t const & aChildMeshIndex,
            moris::Matrix< moris::IdMat > & aNodeIds)
    {
        MORIS_ASSERT(aChildMeshIndex < mNumberOfChildrenMesh, "The requested mesh index is out of bounds.");
        mChildrenMeshes(aChildMeshIndex)->set_node_ids(aNodeIds);

        // Meshes have changed so counts are wrong and need to be updated before use
        mConsistentCounts = false;
    }
    // ----------------------------------------------------------------------------------
    void
    Cut_Mesh::add_node_ids(
            moris::size_t const & aChildMeshIndex,
            Matrix< IdMat >     & aNodeIds)
    {
        MORIS_ASSERT(aChildMeshIndex < mNumberOfChildrenMesh, "The requested mesh index is out of bounds.");
        mChildrenMeshes(aChildMeshIndex)->add_node_ids(aNodeIds);

        // Meshes have changed so counts are wrong and need to be updated before use
        mConsistentCounts = false;
    }
    // ----------------------------------------------------------------------------------
    void
    Cut_Mesh::set_child_element_ids(
            moris::size_t const &   aChildMeshIndex,
            moris::moris_id & aElementIdOffset )
    {
        mChildrenMeshes(aChildMeshIndex)->set_child_element_ids(aElementIdOffset);
    }
    // ----------------------------------------------------------------------------------
    void
    Cut_Mesh::set_child_element_inds(
            moris::size_t const & aChildMeshIndex,
            moris::moris_index &  aElementIndOffset )
    {
        mChildrenMeshes(aChildMeshIndex)->set_child_element_inds(aElementIndOffset);
    }
    // ----------------------------------------------------------------------------------
    moris::Matrix< moris::IdMat > const &
    Cut_Mesh::get_element_ids(
            moris::size_t const & aChildMeshIndex)
    {
        return  mChildrenMeshes(aChildMeshIndex)->get_element_ids();
    }
    // ----------------------------------------------------------------------------------
    moris::Matrix< moris::IdMat >
    Cut_Mesh::get_all_element_ids()
    {
        moris::size_t tNumElems = this->get_num_entities(mtk::EntityRank::ELEMENT);
        moris::Matrix< moris::IdMat > tElementIds(1,tNumElems);

        moris::size_t tCount = 0;
        for(moris::size_t iCM = 0; iCM<this->get_num_child_meshes(); iCM++)
        {
            moris::Matrix< moris::IdMat > const & tCMIds = this->get_element_ids(iCM);
            for(moris::size_t iE = 0; iE<tCMIds.numel(); iE++)
            {
                tElementIds(tCount) = tCMIds(iE);
                tCount++;
            }

        }
        return  tElementIds;
    }
    // ----------------------------------------------------------------------------------
    moris::Matrix< moris::IndexMat >
    Cut_Mesh::get_all_element_inds()
    {
        moris::size_t tNumElems = this->get_num_entities(mtk::EntityRank::ELEMENT);
        moris::Matrix< moris::IndexMat > tElementInds(1,tNumElems);

        moris::size_t tCount = 0;
        for(moris::size_t iCM = 0; iCM<this->get_num_child_meshes(); iCM++)
        {
            moris::Matrix< moris::IndexMat > const & tCMIds = this->get_element_inds(iCM);
            for(moris::size_t iE = 0; iE<tCMIds.numel(); iE++)
            {
                tElementInds(tCount) = tCMIds(iE);
                tCount++;
            }
        }
        return  tElementInds;
    }
    // ----------------------------------------------------------------------------------
    Cell<moris::Matrix< moris::IdMat >>
    Cut_Mesh::get_child_elements_by_phase(
            uint aNumPhases,
            moris::mtk::Mesh const & aBackgroundMeshData)
    {
        moris::size_t tNumElems = this->get_num_entities(mtk::EntityRank::ELEMENT);

        //Initialize output
        Cell<moris::Matrix<moris::IdMat>> tElementsByPhase(aNumPhases);
        moris::Matrix<moris::DDUMat> tPhaseCount(1,aNumPhases,0);
        for(uint i =0; i <aNumPhases; i++)
        {
            tElementsByPhase(i) = moris::Matrix<moris::IdMat>(1,tNumElems);
        }

        for(moris::size_t iCM = 0; iCM<this->get_num_child_meshes(); iCM++)
        {
            Child_Mesh const & tCM = get_child_mesh(iCM);

            moris::Matrix< moris::IdMat >    const & tCMIds     = tCM.get_element_ids();
            moris::Matrix< moris::IndexMat > const & tElemPhase = tCM.get_element_phase_indices();
            for(moris::size_t iE = 0; iE<tCMIds.numel(); iE++)
            {
                moris::moris_index tPhaseInd = tElemPhase(iE);
                uint tCount = tPhaseCount(tElemPhase(iE));
                tElementsByPhase(tPhaseInd)(tCount) = tCMIds(iE);
                tPhaseCount(tElemPhase(iE))++;
            }
        }

        // Resize
        for(uint i =0; i <aNumPhases; i++)
        {
            tElementsByPhase(i).resize(1,tPhaseCount(i));
        }

        return  tElementsByPhase;
    }
    // ----------------------------------------------------------------------------------
    moris::Matrix<moris::IndexMat> const &
    Cut_Mesh::get_element_inds(
            moris::size_t const & aChildMeshIndex) const
    {
        return  mChildrenMeshes(aChildMeshIndex)->get_element_inds();
    }
    // ----------------------------------------------------------------------------------
    moris::Matrix< moris::IndexMat > const &
    Cut_Mesh::get_node_indices(moris::size_t aChildMeshIndex)
    {
        MORIS_ASSERT(aChildMeshIndex < mNumberOfChildrenMesh, "The requested mesh index is out of bounds.");
        return mChildrenMeshes(aChildMeshIndex)->get_node_indices();
    }

    // ----------------------------------------------------------------------------------
    void
    Cut_Mesh::set_child_element_topology( mtk::CellTopology aChildCellTopo )
    {
        mChildElementTopo = aChildCellTopo;
    }
    // ----------------------------------------------------------------------------------
    mtk::CellTopology
    Cut_Mesh::get_child_element_topology()
    {
        return mChildElementTopo;
    }
    // ----------------------------------------------------------------------------------
    moris::size_t
    Cut_Mesh::get_num_child_meshes() const
    {
        return mNumberOfChildrenMesh;
    }
    // ----------------------------------------------------------------------------------
    moris::size_t
    Cut_Mesh::get_num_entities(
            moris::size_t aChildMeshIndex,
            mtk::EntityRank aEntityRank) const
    {
        MORIS_ASSERT(aChildMeshIndex < mNumberOfChildrenMesh, "The requested mesh index is out of bounds.");
        return mChildrenMeshes(aChildMeshIndex)->get_num_entities(aEntityRank);
    }
    // ----------------------------------------------------------------------------------
    moris::size_t
    Cut_Mesh::get_num_entities(
            mtk::EntityRank aEntityRank) const
    {
        // Make sure counts are up to date
        get_entity_counts();

        // Make sure the counts are up to date
        MORIS_ASSERT(mConsistentCounts, "Make sure to call get_entity_counts otherwise the mNumEntities variable is out dated and garbage is returned");

        return mNumEntities((moris::size_t) aEntityRank);
    }
    // ----------------------------------------------------------------------------------
    moris::moris_index
    Cut_Mesh::get_parent_element_index(
            moris::size_t aChildMeshIndex)
    {
        MORIS_ASSERT(aChildMeshIndex < mNumberOfChildrenMesh, "The requested mesh index is out of bounds.");
        return mChildrenMeshes(aChildMeshIndex)->get_parent_element_index();
    }
    // ----------------------------------------------------------------------------------

    void
    Cut_Mesh::get_child_elements_connected_to_parent_facet(
            moris::moris_index const         & aChildMeshIndex,
            moris::moris_index const         & aParentFaceIndex,
            moris::Matrix< moris::IdMat >    & aChildrenElementId,
            moris::Matrix< moris::IndexMat > & aChildrenElementCMInd,
            moris::Matrix< moris::IndexMat > & aFaceOrdinal) const
    {
        // estimate maximum number of elements on face
         const uint tMaxElemOnFace = 100;

        // matrices used throughout routine
        aChildrenElementId.set_size(1,tMaxElemOnFace);
        aChildrenElementCMInd.set_size(1,tMaxElemOnFace);
        aFaceOrdinal.set_size(1,tMaxElemOnFace);

        // define variable for actual number of child elements on face
        uint tNumberOfChildElemsOnFace;

        mChildrenMeshes(aChildMeshIndex)->get_child_elements_connected_to_parent_facet(
                aParentFaceIndex,
                tNumberOfChildElemsOnFace,
                aChildrenElementId,
                aChildrenElementCMInd,
                aFaceOrdinal);

        // resize matrices matrices used throughout routine
        aChildrenElementId.resize(0,tNumberOfChildElemsOnFace);
        aChildrenElementCMInd.resize(0,tNumberOfChildElemsOnFace);
        aFaceOrdinal.resize(0,tNumberOfChildElemsOnFace);
    }

    // ----------------------------------------------------------------------------------
    void
    Cut_Mesh::pack_cut_mesh_by_phase(
            moris::size_t const & aMeshIndex,
            moris::size_t const & aNumPhases,
            Cell<moris::Matrix< moris::IdMat >> & aElementIds,
            Cell<moris::Matrix< moris::IdMat >> & aElementCMInds) const
    {

        mChildrenMeshes(aMeshIndex)->pack_child_mesh_by_phase(aNumPhases,aElementCMInds,aElementIds);
    }
    // ----------------------------------------------------------------------------------
    moris::Matrix< moris::IdMat >
    Cut_Mesh::pack_interface_sides(
            moris_index aGeometryIndex,
            moris_index aPhaseIndex0,
            moris_index aPhaseIndex1,
            moris_index aIndexFlag) const
    {
        uint tNumElem = this->get_num_entities(mtk::EntityRank::ELEMENT);

        moris::Matrix< moris::IdMat > tFullElementIdAndSideOrdinals(tNumElem,2);

        uint tCount = 0;
        uint tNumElemsFromCM = 0;
        for(uint i = 0; i <this->get_num_child_meshes(); i++)
        {
            moris::Matrix< moris::IdMat > tSingleCMElementIdAndSideOrds = mChildrenMeshes(i)->pack_interface_sides(aGeometryIndex,aPhaseIndex0,aPhaseIndex1,aIndexFlag);

            tNumElemsFromCM = tSingleCMElementIdAndSideOrds.n_rows();
            if(tNumElemsFromCM!=0)
            {
                tFullElementIdAndSideOrdinals({tCount,tCount+tNumElemsFromCM-1},{0,1}) = tSingleCMElementIdAndSideOrds({0,tNumElemsFromCM-1},{0,1});
            }

            tCount = tCount + tNumElemsFromCM;
        }

        tFullElementIdAndSideOrdinals.resize(tCount,2);
        return tFullElementIdAndSideOrdinals;
    }

    // ----------------------------------------------------------------------------------
    moris::Matrix< moris::IndexMat >
    Cut_Mesh::pack_interface_sides_loc_inds() const
    {
        uint tNumElem = this->get_num_entities(mtk::EntityRank::ELEMENT);

        moris::Matrix< moris::IdMat > tFullElementIdAndSideOrdinals(tNumElem,2);

        uint tCount = 0;
        uint tNumElemsFromCM = 0;
        for(uint i = 0; i <this->get_num_child_meshes(); i++)
        {
            moris::Matrix< moris::IdMat > tSingleCMElementIdAndSideOrds = mChildrenMeshes(i)->pack_interface_sides_loc_inds();

            tNumElemsFromCM = tSingleCMElementIdAndSideOrds.n_rows();
            if(tNumElemsFromCM!=0)
            {
                tFullElementIdAndSideOrdinals({tCount,tCount+tNumElemsFromCM-1},{0,1}) = tSingleCMElementIdAndSideOrds({0,tNumElemsFromCM-1},{0,1});
            }

            tCount = tCount + tNumElemsFromCM;
        }

        tFullElementIdAndSideOrdinals.resize(tCount,2);
        return tFullElementIdAndSideOrdinals;
    }
    // ----------------------------------------------------------------------------------
    moris::Matrix<moris::IdMat>
    Cut_Mesh::get_full_element_to_node_glob_ids()
    {
        mtk::CellTopology tChildElementTopo = this->get_child_element_topology();
        moris::size_t     tNumElements = this->get_num_entities(mtk::EntityRank::ELEMENT);
        moris::size_t tNumNodesPerElem = 0;
        if(tChildElementTopo == mtk::CellTopology::TET4)
        {
            tNumNodesPerElem = 4;
        }
        else
        {
            MORIS_ERROR(0,"Not implemented");
        }

        moris::Matrix<moris::IdMat> tElementToNodeIds(tNumElements,tNumNodesPerElem);
        moris::size_t tCount = 0;
        for(moris::size_t i = 0; i<this->get_num_child_meshes(); i++)
        {
            moris::Matrix<moris::IdMat> tElementToNodeIdsCM = mChildrenMeshes(i)->get_element_to_node_global();
            for(moris::size_t j = 0; j<tElementToNodeIdsCM.n_rows(); j++)
            {
                tElementToNodeIds.set_row(tCount, tElementToNodeIdsCM.get_row(j));
                tCount++;
            }

        }

        return tElementToNodeIds;
    }
    // ----------------------------------------------------------------------------------
    moris::Matrix<moris::IdMat>
    Cut_Mesh::get_full_element_to_node_loc_inds()
    {
        moris::Matrix<moris::IdMat> tElementToNodeInds;
        if(mChildrenMeshes.size() > 0)
        {
            moris::mtk::Cell_Info const * tCellInfo = mChildrenMeshes(0)->get_cell_info();

            moris::size_t  tNumElements = this->get_num_entities(mtk::EntityRank::ELEMENT);
            moris::uint    tNumNodesPerElem = tCellInfo->get_num_verts();

            tElementToNodeInds.resize(tNumElements,tNumNodesPerElem);
            moris::size_t tCount = 0;
            for(moris::size_t i = 0; i<this->get_num_child_meshes(); i++)
            {
                moris::Matrix<moris::IdMat> tElementToNodeIndsCM = mChildrenMeshes(i)->get_element_to_node();
                for(moris::size_t j = 0; j<tElementToNodeIndsCM.n_rows(); j++)
                {
                    tElementToNodeInds.set_row(tCount, tElementToNodeIndsCM.get_row(j));
                    tCount++;
                }

            }
        }
        return tElementToNodeInds;
    }
    // ----------------------------------------------------------------------------------
    moris::Cell<moris::Matrix<moris::IdMat>>
    Cut_Mesh::get_full_element_to_node_by_phase_glob_ids(moris::uint aNumPhases,
            moris::mtk::Mesh & aBackgroundMeshData)
    {
        mtk::CellTopology tChildElementTopo = this->get_child_element_topology();
        moris::size_t     tNumElements = this->get_num_entities(mtk::EntityRank::ELEMENT);

        moris::size_t tNumNodesPerElem = 0;
        if(tChildElementTopo == mtk::CellTopology::TET4)
        {
            tNumNodesPerElem = 4;
        }
        else if(tChildElementTopo == mtk::CellTopology::TRI3)
        {
            tNumNodesPerElem = 3;
        }
        else
        {
            MORIS_ERROR(0,"Not implemented");
        }

        moris::Cell<moris::Matrix<moris::IdMat>> tElementToNodeIdsByPhase(aNumPhases , moris::Matrix<moris::IdMat>(tNumElements,tNumNodesPerElem));
        moris::Cell<moris::size_t> tCount(aNumPhases,0);

        for(moris::size_t i = 0; i<this->get_num_child_meshes(); i++)
        {
            moris::Matrix< moris::IdMat >  tElementToNodeIdsCM = mChildrenMeshes(i)->get_element_to_node_global();
            moris::Matrix<moris::IndexMat> tElementPhase       = mChildrenMeshes(i)->get_element_phase_indices();
            for(moris::size_t j = 0; j<tElementToNodeIdsCM.n_rows(); j++)
            {
                tElementToNodeIdsByPhase(tElementPhase(j)).set_row(tCount(tElementPhase(j)), tElementToNodeIdsCM.get_row(j));
                tCount(tElementPhase(j))++;
            }
        }

        // size out extra space
        for(moris::uint i = 0; i <aNumPhases; i++)
        {
            tElementToNodeIdsByPhase(i).resize(tCount(i),tNumNodesPerElem);
        }

        return tElementToNodeIdsByPhase;
    }
    // ----------------------------------------------------------------------------------

    Child_Mesh const &
    Cut_Mesh::get_child_mesh(
            moris::size_t const & aChildMeshIndex) const
    {
        return *mChildrenMeshes(aChildMeshIndex);
    }
    // ----------------------------------------------------------------------------------
    Child_Mesh &
    Cut_Mesh::get_child_mesh(
            moris::size_t const & aChildMeshIndex)
    {
        return *mChildrenMeshes(aChildMeshIndex);
    }
    // ----------------------------------------------------------------------------------
    void
    Cut_Mesh::set_num_subphases(
            moris::uint aNumSubPhases)
    {
        mNumSubPhases = aNumSubPhases;
    }
    // ----------------------------------------------------------------------------------
    moris::uint
    Cut_Mesh::get_num_subphases()
    {
        return mNumSubPhases;
    }
    // ----------------------------------------------------------------------------------
    void
    Cut_Mesh::setup_subphase_to_child_mesh_connectivity()
    {
        mSubPhaseIndexToChildMesh              = Matrix<IndexMat>(1,this->get_num_subphases(),MORIS_INDEX_MAX);
        mSubPhaseIndexToChildMeshSubphaseIndex = Matrix<IndexMat>(1,this->get_num_subphases(),MORIS_INDEX_MAX);

        // iterate through child meshes
        for(moris::uint iCM = 0; iCM < this->get_num_child_meshes(); iCM++)
        {
            Child_Mesh const & tChildMesh = this->get_child_mesh(iCM);

            moris::Cell<moris::moris_index> const & tCMSubphaseIndices = tChildMesh.get_subphase_indices();

            // iterate through subphase on child meshes
            for(moris::uint iSP = 0; iSP < tCMSubphaseIndices.size(); iSP++)
            {
                MORIS_ASSERT( mSubPhaseIndexToChildMesh(tCMSubphaseIndices(iSP)) == MORIS_INDEX_MAX, "Subphase to child mesh index already set");

                mSubPhaseIndexToChildMesh(tCMSubphaseIndices(iSP)) = iCM;
                mSubPhaseIndexToChildMeshSubphaseIndex(tCMSubphaseIndices(iSP)) = iSP;
            }
        }
    }
    // ----------------------------------------------------------------------------------
    uint
    Cut_Mesh::get_bulk_phase_index(
            moris_index aSubPhaseIndex)
    {
        moris_index tCMIndex = mSubPhaseIndexToChildMesh(aSubPhaseIndex);
        moris_index tCMSubphaseIndex = mSubPhaseIndexToChildMeshSubphaseIndex(aSubPhaseIndex);

        Cell<moris::moris_index> const & tSubphaseBulk =  this->get_child_mesh(tCMIndex).get_subphase_bin_bulk_phase();

        return tSubphaseBulk(tCMSubphaseIndex);
    }
    // ----------------------------------------------------------------------------------
    void
    Cut_Mesh::populate_subphase_vector(
            moris::Matrix<moris::IndexMat> & aSubphase)
    {
        moris::uint tNumChildMeshes = this->get_num_child_meshes();

        for(moris::uint i = 0 ; i <tNumChildMeshes; i++)
        {
            Child_Mesh & tChildMesh = this->get_child_mesh(i);
            moris::Matrix<moris::IndexMat> const & tElementInds = tChildMesh.get_element_inds();
            Cell<moris_index> const & tSubphaseInds = tChildMesh.get_subphase_indices();

            for(moris::uint iCE = 0; iCE < tElementInds.numel(); iCE++)
            {
                aSubphase(tElementInds(iCE)) = tSubphaseInds(tChildMesh.get_element_subphase_index(iCE));
            }
        }
    }

    // ----------------------------------------------------------------------------------

    moris::Memory_Map
    Cut_Mesh::get_memory_usage()
    {
        // memory map
        moris::Memory_Map tMemoryMap;

        // memory map of child mesh
        moris::Memory_Map tMemoryMapCMs;

        // number of child meshes
        moris::uint tNumChildMeshes = this->get_num_child_meshes();

        //sum up all child mesh data
        for (moris::uint i = 0; i < tNumChildMeshes; i++)
        {
            tMemoryMapCMs = tMemoryMapCMs + mChildrenMeshes(i)->get_memory_usage();
        }

        tMemoryMap.mMemoryMapData["Child Meshes"] = tMemoryMapCMs.sum();
        tMemoryMap.mMemoryMapData["mSpatialDim"]  = sizeof(mSpatialDim);
        tMemoryMap.mMemoryMapData["mOwnedChildrenMeshes"] = mOwnedChildrenMeshes.capacity();
        tMemoryMap.mMemoryMapData["mNotOwnedChildrenMeshes"] = mNotOwnedChildrenMeshes.capacity();
        tMemoryMap.mMemoryMapData["mNotOwnedOwningProc"] = mNotOwnedOwningProc.capacity();
        tMemoryMap.mMemoryMapData["mNumberOfChildrenMesh"] = sizeof(mNumberOfChildrenMesh);
        tMemoryMap.mMemoryMapData["mConsistentCounts"] = sizeof(mConsistentCounts);
        tMemoryMap.mMemoryMapData["mNumEntities"] = mNumEntities.capacity();
        tMemoryMap.mMemoryMapData["mChildElementTopo"] = sizeof(mChildElementTopo);
        tMemoryMap.mMemoryMapData["mNumSubPhases"] = sizeof(mNumSubPhases);
        tMemoryMap.mMemoryMapData["mSubPhaseIndexToChildMesh"] = mSubPhaseIndexToChildMesh.capacity();
        tMemoryMap.mMemoryMapData["mSubPhaseIndexToChildMeshSubphaseIndex"] = mSubPhaseIndexToChildMeshSubphaseIndex.capacity();

        return tMemoryMap;
    }

    // ----------------------------------------------------------------------------------
    Matrix<IndexMat> const &
    Cut_Mesh::get_subphase_to_child_mesh_connectivity()
    {
        return mSubPhaseIndexToChildMesh;
    }
    // ----------------------------------------------------------------------------------
    moris_id
    Cut_Mesh::get_subphase_id(
            moris_index aSubphaseIndex)
    {
        moris_index tChildMeshIndex   = mSubPhaseIndexToChildMesh(aSubphaseIndex);
        moris_index tChildMeshSPIndex = mSubPhaseIndexToChildMeshSubphaseIndex(aSubphaseIndex);
        return mChildrenMeshes(tChildMeshIndex)->get_subphase_ids()(tChildMeshSPIndex);
    }
    // ----------------------------------------------------------------------------------
    void
    Cut_Mesh::add_child_mesh_groups(
            Cell<Child_Mesh*>   & tOwnedChildrenMeshes,
            Cell<Child_Mesh*>   & tNotOwnedChildrenMeshes,
            Cell<moris_id>      & tNotOwnedOwningProc)
    {
        mOwnedChildrenMeshes    = tOwnedChildrenMeshes;
        mNotOwnedChildrenMeshes = tNotOwnedChildrenMeshes;
        mNotOwnedOwningProc     = tNotOwnedOwningProc;
    }
    // ----------------------------------------------------------------------------------
    Cell<Child_Mesh*> &
    Cut_Mesh::get_owned_child_meshes()
    {
        return mOwnedChildrenMeshes;
    }
    // ----------------------------------------------------------------------------------
    Cell<Child_Mesh*> &
    Cut_Mesh::get_not_owned_child_meshes()
    {
        return mNotOwnedChildrenMeshes;
    }
    // ----------------------------------------------------------------------------------
    Cell<moris_id> &
    Cut_Mesh::get_not_owned_child_owners()
    {
        return mNotOwnedOwningProc;
    }

    void
    Cut_Mesh::remove_all_child_meshes_but_selected(Cell<moris::uint> const & aMeshesToKeep,
                                                   Cell<moris::uint> const & aMeshesToDelete,
                                                   Cell<moris_index> & aCellsToRemoveFromMesh)
    {
        uint tNumCells = 0;
        // count the cells to delete
        for(moris::uint iCM = 0; iCM < aMeshesToDelete.size(); iCM++)
        {
            tNumCells = tNumCells+mChildrenMeshes(aMeshesToDelete(iCM))->get_num_entities(mtk::EntityRank::ELEMENT);
        }

        aCellsToRemoveFromMesh.reserve(tNumCells);

        Cell<Child_Mesh*> tCMToKeep(aMeshesToKeep.size());

        // collect a temporary list then populate member data
        for(moris::uint iCM = 0; iCM < aMeshesToKeep.size(); iCM++)
        {
            tCMToKeep(iCM) = mChildrenMeshes(aMeshesToKeep(iCM));
        }

        // delete the others
        for(moris::uint iCM = 0; iCM < aMeshesToDelete.size(); iCM++)
        {

            Matrix<IndexMat> const & tCellInds = mChildrenMeshes(aMeshesToDelete(iCM))->get_element_inds();
            for(moris::uint iCell = 0; iCell<mChildrenMeshes(aMeshesToDelete(iCM))->get_num_entities(mtk::EntityRank::ELEMENT);iCell++)
            {
                aCellsToRemoveFromMesh.push_back(tCellInds(iCell));
            }

            delete mChildrenMeshes(aMeshesToDelete(iCM));
        }

        mChildrenMeshes.resize(aMeshesToKeep.size());
        for(moris::uint iCM = 0; iCM < tCMToKeep.size(); iCM++)
        {
            mChildrenMeshes(iCM) = tCMToKeep(iCM);
        }

        mNumberOfChildrenMesh = mChildrenMeshes.size();

        // recount the mesh
        mConsistentCounts = false;
        this->get_entity_counts();

        this->setup_subphase_to_child_mesh_connectivity();

    }

    void
    Cut_Mesh::reindex_cells(Cell<moris_index> & aOldIndexToNewCellIndex)
    {
        for(moris::uint i = 0; i < this->get_num_child_meshes(); i++)
        {
            mChildrenMeshes(i)->reindex_cells(aOldIndexToNewCellIndex);
        }
    }

    // ----------------------------------------------------------------------------------
    void
    Cut_Mesh::get_entity_counts() const
    {
        // If the counts aren't up to date then recount them.
        if(!mConsistentCounts)
        {
            moris::Cell<moris::size_t> tRanks = {0,1,2,3};

            if(mSpatialDim == 2)
            {
                tRanks = {0,1,3};
            }

            mtk::EntityRank tRank = mtk::EntityRank::NODE;
            for(auto iR:tRanks)
            {
                moris::size_t tCount = 0;

                // Cast r to enum
                tRank = static_cast<mtk::EntityRank>(iR);
                for(moris::size_t m = 0; m<mNumberOfChildrenMesh; m++)
                {
                    tCount += mChildrenMeshes(m)->get_num_entities(tRank);
                }
                mNumEntities(iR) = tCount;
            }
        }

        // This should be the only spot where this is changed to true
        mConsistentCounts = true;

    }
    // ----------------------------------------------------------------------------------
}

