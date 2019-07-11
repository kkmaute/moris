/*
 * cl_XTK_Cut_Mesh.cpp
 *
 *  Created on: Jul 4, 2019
 *      Author: doble
 */

#include "cl_XTK_Cut_Mesh.hpp"


namespace xtk
{
// ----------------------------------------------------------------------------------
Cut_Mesh::Cut_Mesh(moris::uint aModelDim) :
                   mNumberOfChildrenMesh(0),
                   mChildrenMeshes(1, Child_Mesh()),
                   mNumEntities(4,0),
                   mChildElementTopo(CellTopology::TET4){}
// ----------------------------------------------------------------------------------
Cut_Mesh::Cut_Mesh(moris::size_t aNumSimpleMesh,
                   moris::size_t aModelDim) :
                mNumberOfChildrenMesh(aNumSimpleMesh),
                mChildrenMeshes(aNumSimpleMesh),
                mChildElementTopo(CellTopology::TET4)
{
    for(moris::size_t i = 0; i <aNumSimpleMesh; i++)
    {
        mChildrenMeshes(i) = Child_Mesh();
    }
}
// ----------------------------------------------------------------------------------
void
Cut_Mesh::inititalize_new_child_meshes(moris::size_t aNumNewChildMesh, moris::size_t aModelDim)
{

    mNumberOfChildrenMesh = aNumNewChildMesh + mNumberOfChildrenMesh;
    mChildrenMeshes.resize(mNumberOfChildrenMesh, Child_Mesh());
}
// ----------------------------------------------------------------------------------
void
Cut_Mesh::generate_templated_mesh(moris::size_t     aChildMeshIndex,
                                  enum TemplateType aTemplate)
{
    MORIS_ASSERT(aChildMeshIndex < mNumberOfChildrenMesh, "The requested mesh index is out of bounds. Consider allocating more space");
    mChildrenMeshes(aChildMeshIndex).modify_child_mesh(aTemplate);


    // Meshes have changed so counts are wrong
    mConsistentCounts = false;
}
// ----------------------------------------------------------------------------------
void
Cut_Mesh::generate_templated_mesh(enum TemplateType aTemplate)
{
    for (moris::size_t i = 0; i < mNumberOfChildrenMesh; i++)
    {
        mChildrenMeshes(i).modify_child_mesh(aTemplate);
    }

    // Meshes have changed so counts are wrong
    mConsistentCounts = false;
}
// ----------------------------------------------------------------------------------
void
Cut_Mesh::generate_templated_mesh(Matrix< IndexMat > const & aChildMeshIndices,
                                  enum TemplateType          aTemplate)
{
    for (moris::size_t i = 0; i < aChildMeshIndices.n_cols(); i++)
    {
        mChildrenMeshes(aChildMeshIndices(0,i)).modify_child_mesh(aTemplate);
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
        mChildrenMeshes(i).convert_tet4_to_tet10_child();
    }

    mChildElementTopo = CellTopology::TET10;
}
// ----------------------------------------------------------------------------------
void
Cut_Mesh::initialize_new_mesh_from_parent_element(moris::size_t              aChildMeshIndex,
                                                  enum TemplateType          aTemplate,
                                                  Matrix< IndexMat >       & aNodeIndices,
                                                  Cell<Matrix< IndexMat >> & aParentEntities)
 {
     MORIS_ASSERT(aChildMeshIndex < mNumberOfChildrenMesh, "The requested mesh index is out of bounds.");

     // Construct a template and initialize this new mesh with the template

     // Set all node parent ranks to EntityRank::NODE which = 0;
     moris::Matrix< moris::DDSTMat > tElementNodeParentRanks(1,8);
     tElementNodeParentRanks.fill(0);

     // Set all edge parent ranks to EntityRank::EDGE which = 1;
     moris::Matrix< moris::DDSTMat > tParentEdgeRanks(1,aParentEntities(1).numel());
     tParentEdgeRanks.fill(1);

     // Set all node parent ranks to EntityRank::FACE which = 2;
     moris::Matrix< moris::DDSTMat > tParentFaceRanks(1,aParentEntities(2).numel());
     tParentFaceRanks.fill(2);

     // No interface sides
     moris::Matrix< moris::DDSTMat > tInterfaceSides(1,1);
     tInterfaceSides.fill(std::numeric_limits<moris::size_t>::max());

     // Note for this: child mesh node indices, parent node indices, and element to node connectivity are
     // the same thing. This is why aNodeIndices appears 3 times in this call.
     mChildrenMeshes(aChildMeshIndex) = Child_Mesh(aParentEntities(3)(0,0),
                                                        aNodeIndices,
                                                        aNodeIndices,
                                                        tElementNodeParentRanks,
                                                        aNodeIndices,
                                                        aParentEntities(1),
                                                        tParentEdgeRanks,
                                                        aParentEntities(2),
                                                        tParentFaceRanks,
                                                        tInterfaceSides);

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
             mChildrenMeshes(aChildMeshIndex).allocate_parametric_coordinates(8,3);
             mChildrenMeshes(aChildMeshIndex).add_node_parametric_coordinate(aNodeIndices,tParamCoords);


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
             mChildrenMeshes(aChildMeshIndex).allocate_parametric_coordinates(4,4);
             mChildrenMeshes(aChildMeshIndex).add_node_parametric_coordinate(aNodeIndices,tParamCoords);
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
Cut_Mesh::modify_templated_mesh(moris::size_t     aChildMeshIndex,
                                enum TemplateType aTemplate)
{
    MORIS_ASSERT(aChildMeshIndex < mNumberOfChildrenMesh, "The requested mesh index is out of bounds.");
    mChildrenMeshes(aChildMeshIndex).modify_child_mesh(aTemplate);

    // Meshes have changed so counts are wrong
    mConsistentCounts = false;
}
// ----------------------------------------------------------------------------------
void
Cut_Mesh::add_entity_to_intersect_connectivity(moris::size_t aChildMeshIndex,
                                               moris::size_t aDPrime1Ind,
                                               moris::size_t aDPrime2Ind,
                                               moris::size_t aReturnType)
{
    MORIS_ASSERT(aChildMeshIndex < mNumberOfChildrenMesh, "The requested mesh index is out of bounds.");
    mChildrenMeshes(aChildMeshIndex).add_entity_to_intersect_connectivity(aDPrime1Ind, aDPrime2Ind, aReturnType);
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
    MORIS_ASSERT(mChildElementTopo == CellTopology::TET4,"Interface unzipping only tested on child meshes with tet4s");

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
Cut_Mesh::set_node_index(moris::size_t const & aChildMeshIndex,
               Matrix< IndexMat >  & aNodeInd)
{
    MORIS_ASSERT(aChildMeshIndex < mNumberOfChildrenMesh, "The requested mesh index is out of bounds.");
    mChildrenMeshes(aChildMeshIndex).add_node_indices(aNodeInd);

    // Meshes have changed so counts are wrong
    mConsistentCounts = false;
}
// ----------------------------------------------------------------------------------
void
Cut_Mesh::set_node_ids(moris::size_t const & aChildMeshIndex,
             moris::Matrix< moris::IdMat > & aNodeIds)
{
    MORIS_ASSERT(aChildMeshIndex < mNumberOfChildrenMesh, "The requested mesh index is out of bounds.");
    mChildrenMeshes(aChildMeshIndex).set_node_ids(aNodeIds);

    // Meshes have changed so counts are wrong and need to be updated before use
    mConsistentCounts = false;
}
// ----------------------------------------------------------------------------------
void
Cut_Mesh::add_node_ids(moris::size_t const & aChildMeshIndex,
             Matrix< IdMat >     & aNodeIds)
{
    MORIS_ASSERT(aChildMeshIndex < mNumberOfChildrenMesh, "The requested mesh index is out of bounds.");
    mChildrenMeshes(aChildMeshIndex).add_node_ids(aNodeIds);

    // Meshes have changed so counts are wrong and need to be updated before use
    mConsistentCounts = false;
}
// ----------------------------------------------------------------------------------
void
Cut_Mesh::set_child_element_ids(moris::size_t const &   aChildMeshIndex,
                      moris::moris_id & aElementIdOffset )
{
    mChildrenMeshes(aChildMeshIndex).set_child_element_ids(aElementIdOffset);
}
// ----------------------------------------------------------------------------------
void
Cut_Mesh::set_child_element_inds(moris::size_t const & aChildMeshIndex,
                       moris::moris_index &  aElementIndOffset )
{
    mChildrenMeshes(aChildMeshIndex).set_child_element_inds(aElementIndOffset);
}
// ----------------------------------------------------------------------------------
moris::Matrix< moris::IdMat > const &
Cut_Mesh::get_element_ids(moris::size_t const & aChildMeshIndex)
{
    return  mChildrenMeshes(aChildMeshIndex).get_element_ids();
}
// ----------------------------------------------------------------------------------
moris::Matrix< moris::IdMat >
Cut_Mesh::get_all_element_ids()
{
   moris::size_t tNumElems = this->get_num_entities(EntityRank::ELEMENT);
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
    moris::size_t tNumElems = this->get_num_entities(EntityRank::ELEMENT);
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
Cut_Mesh::get_child_elements_by_phase(uint aNumPhases,
                            moris::mtk::Mesh const & aBackgroundMeshData)
{
    moris::size_t tNumElems = this->get_num_entities(EntityRank::ELEMENT);

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
Cut_Mesh::get_element_inds(moris::size_t const & aChildMeshIndex) const
{
    return  mChildrenMeshes(aChildMeshIndex).get_element_inds();
}
// ----------------------------------------------------------------------------------
moris::Matrix< moris::IndexMat > const &
Cut_Mesh::get_node_indices(moris::size_t aChildMeshIndex)
{
    MORIS_ASSERT(aChildMeshIndex < mNumberOfChildrenMesh, "The requested mesh index is out of bounds.");
    return mChildrenMeshes(aChildMeshIndex).get_node_indices();
}

// ----------------------------------------------------------------------------------
enum CellTopology
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
Cut_Mesh::get_num_entities(moris::size_t aChildMeshIndex, enum EntityRank aEntityRank) const
{
    MORIS_ASSERT(aChildMeshIndex < mNumberOfChildrenMesh, "The requested mesh index is out of bounds.");
    return mChildrenMeshes(aChildMeshIndex).get_num_entities(aEntityRank);
}
// ----------------------------------------------------------------------------------
moris::size_t
Cut_Mesh::get_num_entities(enum EntityRank aEntityRank) const
{
    // Make sure counts are up to date
    get_entity_counts();

    // Make sure the counts are up to date
    MORIS_ASSERT(mConsistentCounts, "Make sure to call get_entity_counts otherwise the mNumEntities variable is out dated and garbage is returned");

    return mNumEntities((moris::size_t) aEntityRank);
}
// ----------------------------------------------------------------------------------
moris::moris_index
Cut_Mesh::get_parent_element_index(moris::size_t aChildMeshIndex)
{
    MORIS_ASSERT(aChildMeshIndex < mNumberOfChildrenMesh, "The requested mesh index is out of bounds.");
    return mChildrenMeshes(aChildMeshIndex).get_parent_element_index();
}
// ----------------------------------------------------------------------------------
void
Cut_Mesh::get_child_elements_connected_to_parent_face(moris::moris_index const & aChildMeshIndex,
                                                 moris::moris_index const & aParentFaceIndex,
                                                 moris::Matrix< moris::IdMat > & aChildrenElementId,
                                                 moris::Matrix< moris::IndexMat > & aChildrenElementCMInd,
                                                 moris::Matrix< moris::IndexMat > & aFaceOrdinal) const
{
    mChildrenMeshes(aChildMeshIndex).get_child_elements_connected_to_parent_face(aParentFaceIndex,aChildrenElementId,aChildrenElementCMInd,aFaceOrdinal);
}

// ----------------------------------------------------------------------------------
void
Cut_Mesh::pack_cut_mesh_by_phase(moris::size_t const & aMeshIndex,
                            moris::size_t const & aNumPhases,
                            Cell<moris::Matrix< moris::IdMat >> & aElementIds,
                            Cell<moris::Matrix< moris::IdMat >> & aElementCMInds) const
{

    mChildrenMeshes(aMeshIndex).pack_child_mesh_by_phase(aNumPhases,aElementCMInds,aElementIds);
}
// ----------------------------------------------------------------------------------
moris::Matrix< moris::IdMat >
Cut_Mesh::pack_interface_sides(bool aIndexFlag,
                               moris_index aPhaseIndex) const
{
    uint tNumElem = this->get_num_entities(EntityRank::ELEMENT);

    moris::Matrix< moris::IdMat > tFullElementIdAndSideOrdinals(tNumElem,2);

    uint tCount = 0;
    uint tNumElemsFromCM = 0;
    for(uint i = 0; i <this->get_num_child_meshes(); i++)
    {
       moris::Matrix< moris::IdMat > tSingleCMElementIdAndSideOrds = mChildrenMeshes(i).pack_interface_sides(aIndexFlag,aPhaseIndex);

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
    uint tNumElem = this->get_num_entities(EntityRank::ELEMENT);

    moris::Matrix< moris::IdMat > tFullElementIdAndSideOrdinals(tNumElem,2);

    uint tCount = 0;
    uint tNumElemsFromCM = 0;
    for(uint i = 0; i <this->get_num_child_meshes(); i++)
    {
       moris::Matrix< moris::IdMat > tSingleCMElementIdAndSideOrds = mChildrenMeshes(i).pack_interface_sides_loc_inds();

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
    enum CellTopology tChildElementTopo = this->get_child_element_topology();
    moris::size_t     tNumElements = this->get_num_entities(EntityRank::ELEMENT);
    moris::size_t tNumNodesPerElem = 0;
    if(tChildElementTopo == CellTopology::TET4)
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
       moris::Matrix<moris::IdMat> tElementToNodeIdsCM = mChildrenMeshes(i).get_element_to_node_global();
       for(moris::size_t j = 0; j<tElementToNodeIdsCM.n_rows(); j++)
       {
           tElementToNodeIds.set_row(tCount, tElementToNodeIdsCM.get_row(j));
           tCount++;
       }

    }

    return tElementToNodeIds;
}
// ----------------------------------------------------------------------------------
moris::Cell<moris::Matrix<moris::IdMat>>
Cut_Mesh::get_full_element_to_node_by_phase_glob_ids(moris::uint aNumPhases,
                                           moris::mtk::Mesh & aBackgroundMeshData)
{
    enum CellTopology tChildElementTopo = this->get_child_element_topology();
    moris::size_t     tNumElements = this->get_num_entities(EntityRank::ELEMENT);

    moris::size_t tNumNodesPerElem = 0;
    if(tChildElementTopo == CellTopology::TET4)
    {
        tNumNodesPerElem = 4;
    }
    else
    {
        MORIS_ERROR(0,"Not implemented");
    }

    moris::Cell<moris::Matrix<moris::IdMat>> tElementToNodeIdsByPhase(aNumPhases , moris::Matrix<moris::IdMat>(tNumElements,tNumNodesPerElem));
    moris::Cell<moris::size_t> tCount(aNumPhases,0);

    for(moris::size_t i = 0; i<this->get_num_child_meshes(); i++)
    {
        moris::Matrix< moris::IdMat >  tElementToNodeIdsCM = mChildrenMeshes(i).get_element_to_node_global();
        moris::Matrix<moris::IndexMat> tElementPhase       = mChildrenMeshes(i).get_element_phase_indices();
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
Cut_Mesh::get_child_mesh(moris::size_t const & aChildMeshIndex) const
    {
        return mChildrenMeshes(aChildMeshIndex);
    }
// ----------------------------------------------------------------------------------
Child_Mesh &
Cut_Mesh::get_child_mesh(moris::size_t const & aChildMeshIndex)
{
    return mChildrenMeshes(aChildMeshIndex);
}
// ----------------------------------------------------------------------------------
void
Cut_Mesh::add_child_mesh_groups(    Cell<Child_Mesh*>   & tOwnedNotSharedChildrenMeshes,
                          Cell<Child_Mesh*>   & tOwnedSharedChildrenMeshes,
                          Cell<Child_Mesh*>   & tNotOwnedSharedChildrenMeshes,
                          Cell<Matrix<IdMat>> & tOwnedSharedOtherProcs,
                          Cell<moris_id>      & tNotOwnedSharedOwningProc)
{
    mOwnedNotSharedChildrenMeshes = tOwnedNotSharedChildrenMeshes;
    mOwnedSharedChildrenMeshes    = tOwnedSharedChildrenMeshes;
    mNotOwnedSharedChildrenMeshes = tNotOwnedSharedChildrenMeshes;
    mOwnedSharedOtherProcs        = tOwnedSharedOtherProcs;
    mNotOwnedSharedOwningProc     = tNotOwnedSharedOwningProc;
}
// ----------------------------------------------------------------------------------
Cell<Child_Mesh*> &
Cut_Mesh::get_owned_not_shared_child_meshes()
{
    return mOwnedNotSharedChildrenMeshes;
}
// ----------------------------------------------------------------------------------
Cell<Child_Mesh*> &
Cut_Mesh::get_owned_shared_child_meshes()
{
    return mOwnedSharedChildrenMeshes;
}
// ----------------------------------------------------------------------------------
Cell<Matrix<IdMat>> &
Cut_Mesh::get_owned_shared_child_meshes_sharing_procs()
{
    return mOwnedSharedOtherProcs;
}
// ----------------------------------------------------------------------------------
Cell<Child_Mesh*> &
Cut_Mesh::get_not_owned_shared_child_meshes()
{
    return mNotOwnedSharedChildrenMeshes;
}
// ----------------------------------------------------------------------------------
Cell<moris_id> &
Cut_Mesh::get_not_owned_shared_child_owners()
{
    return mNotOwnedSharedOwningProc;
}
// ----------------------------------------------------------------------------------
void
Cut_Mesh::get_entity_counts() const
{
    // If the counts aren't up to date then recount them.
    if(!mConsistentCounts)
    {
        enum EntityRank tRank = EntityRank::NODE;
        for(moris::size_t r = 0; r<(moris::size_t)EntityRank::ELEMENT+1; r++)
        {
            moris::size_t tCount = 0;
            // Cast r to enum
            tRank = static_cast<EntityRank>(r);
            for(moris::size_t m = 0; m<mNumberOfChildrenMesh; m++)
            {
                tCount += mChildrenMeshes(m).get_num_entities(tRank);
            }
            mNumEntities(r) = tCount;
        }

        // This should be the only spot where this is changed to true
        mConsistentCounts = true;
    }
}
// ----------------------------------------------------------------------------------
}


