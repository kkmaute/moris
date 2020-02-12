/*
 * UT_XTK_Model_Parallel.cpp
 *
 *  Created on: Nov 4, 2019
 *      Author: doble
 */
#include "catch.hpp"
#include "cl_XTK_Model.hpp"
#include "cl_XTK_Enums.hpp"
#include "cl_XTK_Cut_Mesh.hpp"
#include "cl_XTK_Enrichment.hpp"
#include "cl_XTK_Enriched_Interpolation_Mesh.hpp"
#include "cl_MGE_Geometry_Engine.hpp"
#include "cl_Sphere.hpp"
#include "cl_MPI_Tools.hpp"
#include "fn_trans.hpp"
#include "fn_all_true.hpp"
#include "fn_sort.hpp"
#include "op_equal_equal.hpp"


#include "cl_GEN_Sphere.hpp"

namespace xtk
{

template<typename T>
void
convert_vector_to_sparse(Cell<Matrix<T>> const & aCellOfMatrixToConvert,
                              Matrix<T> &        aData,
                              Matrix<T> &        aOffset)
{
    // figure out how big the data is
    moris::uint tSize = 0;
    for(moris::uint i = 0; i<aCellOfMatrixToConvert.size(); i++)
    {
        tSize += aCellOfMatrixToConvert(i).numel();

        MORIS_ASSERT(isvector(aCellOfMatrixToConvert(i)),"only works on vectors");
    }

    // size outputs
    aData.resize(1,tSize);
    aOffset.resize(1,aCellOfMatrixToConvert.size()+1);
    aOffset(0) = 0;

    moris::uint tStart = 0;
    moris::uint tEnd   = 0;
    moris::uint tCount = 1;
    for(moris::uint i = 0; i<aCellOfMatrixToConvert.size(); i++)
    {
        // figure out end
        tEnd = tStart + aCellOfMatrixToConvert(i).numel() - 1;

        // set data and offsets
        aData({0,0},{tStart,tEnd}) = aCellOfMatrixToConvert(i).matrix_data();
        aOffset(tCount) = tEnd;

        // increment
        tStart = tEnd+1;
        tCount++;
    }

}

template<typename T>
void
convert_sparse_to_vector(Matrix<T> const & aData,
                         Matrix<T> const & aOffset,
                         Cell<Matrix<T>> & aCellOfMatrixToConvert)
{
    moris::uint tSize = aOffset.numel() - 1;
    aCellOfMatrixToConvert.resize(tSize);

    moris::uint tCount = 0;
    moris::uint tStart = 0;
    moris::uint tEnd   = 0;
    for(moris::uint i = 1; i < tSize+1; i++)
    {
        tEnd = aOffset(i);

        aCellOfMatrixToConvert(i-1).resize(1,tEnd-tStart);
        aCellOfMatrixToConvert(i-1).matrix_data() = aData({0,0},{tStart,tEnd});

        tStart = tEnd+1;
    }

}

/*!
 * Gathers all coordinates on proc 0 and verifies that all nodes with the same id have the same coordinate
 * This function is for test purposes and is not very efficient
 */
inline
bool
verify_parallel_node_coordinates(Model & aModel)
{
    bool tValidCoords = true;

    moris::Matrix<DDRMat>   tMyNodeCoords = aModel.get_background_mesh().get_all_node_coordinates_loc_inds();
    moris::Matrix<IndexMat> tMyNodeMap    = aModel.get_background_mesh().get_local_to_global_map(EntityRank::NODE);
    if(par_rank() > 0)
    {
        nonblocking_send(tMyNodeCoords,tMyNodeCoords.n_rows(),tMyNodeCoords.n_cols(),0,11);
        nonblocking_send(tMyNodeMap,tMyNodeMap.n_rows(),tMyNodeMap.n_cols(),0,12);
    }

    barrier();

    if(par_rank() == 0)
    {
        moris::Cell<moris::Matrix<DDRMat>>   tNodeCoords(par_size(),moris::Matrix<DDRMat>(1,3));
        moris::Cell<moris::Matrix<IndexMat>> tNodeMaps(par_size(),moris::Matrix<IndexMat>(1,1));

        tNodeCoords(0) = tMyNodeCoords;
        tNodeMaps(0)   = tMyNodeMap;

        // total number of nodes with duplicate
        moris::uint tNumNodesTotalDup = tMyNodeMap.numel();

        // go and receive all processors node coords and maps
        for(moris::uint i= 1; i <(uint)par_size(); i++)
        {
            receive_col_known(tNodeCoords(i), 3, i,11);
            receive(tNodeMaps(i), 1, i,12);

            tNumNodesTotalDup = tNumNodesTotalDup + tNodeMaps(i).numel();
        }
        // concatenate matrices into 1 which is easier to work with
        moris::Matrix<moris::DDRMat> tAllNodeCoords(tNumNodesTotalDup,3);
        moris::Matrix<moris::IndexMat> tAllNodeMaps(1,tNumNodesTotalDup);
        moris_index tStart = 0;
        moris_index tEnd = 0;

        for(moris::uint i= 0; i <(uint)par_size(); i++)
        {
            tEnd = tStart + tNodeCoords(i).n_rows() - 1;

            tAllNodeCoords({tStart,tEnd},{0,2}) = tNodeCoords(i).matrix_data();
            tAllNodeMaps({0,0},{tStart,tEnd}) = tNodeMaps(i).matrix_data();
            tStart = tEnd+1;
        }
        // check coordinates
        std::unordered_map<moris_id,moris_index> tNodeIdToIndMap;
        for( moris::uint i = 0 ; i < tAllNodeMaps.numel(); i++)
        {
            auto tIter = tNodeIdToIndMap.find(tAllNodeMaps(i));
            if(tIter == tNodeIdToIndMap.end())
            {
                tNodeIdToIndMap[tAllNodeMaps(i)] = i;
            }

            else
            {
                moris_index tIndex = tIter->second;
                if(!all_true(tAllNodeCoords.get_row(i) == tAllNodeCoords.get_row(tIndex)))
                {
                    tValidCoords = false;
                }
            }
        }
    }

    barrier();

    return tValidCoords;

}

inline
void
collect_child_elements_and_maps(Model & aModel,
                                Matrix<IndexMat> & aAllChildElements,
                                Matrix<IndexMat> & aAllChildElementIdMap)
{
    // data structures
    Cut_Mesh & tCutMesh = aModel.get_cut_mesh();
    moris::ge::GEN_Geometry_Engine & tGeometryEngine = aModel.get_geom_engine();
    Background_Mesh & tBackgroundMesh = aModel.get_background_mesh();

    // Children element nodes connected to elements
    moris::Cell<Matrix<IdMat>>  tElementToNodeChildrenByPhase = tCutMesh.get_full_element_to_node_by_phase_glob_ids(tGeometryEngine.get_num_bulk_phase(), tBackgroundMesh.get_mesh_data());

    // Child element ids
    moris::Cell<Matrix<IdMat>>  tChildElementsByPhase = tCutMesh.get_child_elements_by_phase(tGeometryEngine.get_num_bulk_phase(),tBackgroundMesh.get_mesh_data());

    moris::uint tNumNodesPerElem = tElementToNodeChildrenByPhase(0).n_cols();
    moris::uint tNumElements   = 0;
    for(moris::uint i =0; i < tChildElementsByPhase.size(); i++)
    {
        tNumElements = tNumElements + tChildElementsByPhase(i).numel();
    }

    Matrix<IndexMat> tMyElementToNode(tNumElements,tNumNodesPerElem);
    Matrix<IndexMat> tMyElementIds(1,tNumElements);
    moris_index tStart = 0;
    moris_index tEnd = 0;
    for(moris::uint i =0; i < tChildElementsByPhase.size(); i++)
    {

        tEnd = tStart + tElementToNodeChildrenByPhase(i).n_rows() - 1;
        tMyElementToNode({tStart,tEnd},{0,tNumNodesPerElem-1}) = tElementToNodeChildrenByPhase(i).matrix_data();
        tMyElementIds({0,0},{tStart,tEnd}) = tChildElementsByPhase(i).matrix_data();
        tStart = tEnd+1;
    }

    if(par_rank() > 0)
    {
        nonblocking_send(tMyElementToNode,tMyElementToNode.n_rows(),tMyElementToNode.n_cols(),0,11);
        nonblocking_send(tMyElementIds,tMyElementIds.n_rows(),tMyElementIds.n_cols(),0,12);
    }

    barrier();

    if(par_rank() == 0)
    {
        moris::Cell<moris::Matrix<IndexMat>> tElementToNodeCells(par_size(),moris::Matrix<IndexMat>(1,1));
        moris::Cell<moris::Matrix<IndexMat>> tElementMapsCells(par_size(),moris::Matrix<IndexMat>(1,1));

        tElementToNodeCells(0) = tMyElementToNode;
        tElementMapsCells(0)   = tMyElementIds;

        // total number of nodes with duplicate
        moris::uint tNumElementsTotal = tMyElementIds.numel();

        // go and receive all processors node coords and maps
        for(moris::uint i= 1; i <(uint)par_size(); i++)
        {
            receive_col_known(tElementToNodeCells(i), tNumNodesPerElem, i,11);
            receive(tElementMapsCells(i), 1, i,12);

            tNumElementsTotal = tNumElementsTotal + tElementMapsCells(i).numel();
        }
        // concatenate matrices into 1 which is easier to work with
        aAllChildElements.resize(tNumElementsTotal,tNumNodesPerElem);
        aAllChildElementIdMap.resize(1,tNumElementsTotal);
        tStart = 0;
        tEnd = 0;

        for(moris::uint i= 0; i <(uint)par_size(); i++)
        {
            tEnd = tStart + tElementToNodeCells(i).n_rows() - 1;

            aAllChildElements({tStart,tEnd},{0,tNumNodesPerElem-1}) = tElementToNodeCells(i).matrix_data();
            aAllChildElementIdMap({0,0},{tStart,tEnd}) = tElementMapsCells(i).matrix_data();
            tStart = tEnd+1;
        }
    }
}

inline
bool
verify_parallel_elements(Model & aModel)
{
    bool tValidElements = true;
    Matrix<IndexMat> tAllChildElements;
    Matrix<IndexMat> tAllChildElementIdMap;
    collect_child_elements_and_maps(aModel, tAllChildElements, tAllChildElementIdMap);

    // verify connectivity
    if(par_rank() == 0)
    {
        // check coordinates
        std::unordered_map<moris_id,moris_index> tElementIdToIndMap;
        for( moris::uint i = 0 ; i < tAllChildElementIdMap.numel(); i++)
        {
            auto tIter = tElementIdToIndMap.find(tAllChildElementIdMap(i));
            if(tIter == tElementIdToIndMap.end())
            {
                tElementIdToIndMap[tAllChildElementIdMap(i)] = i;
            }
            else
            {
                moris_index tIndex = tIter->second;
                if(!all_true(tAllChildElements.get_row(i) == tAllChildElements.get_row(tIndex)))
                {
                    tValidElements = false;
                }
            }
        }
    }
    barrier();

    return tValidElements;
}

inline
bool
verify_parallel_subphase_ids(Model & aModel)
{
//    moris_index tId = 4;
//    moris_index tIndex = aModel.get_background_mesh().get_mesh_data().get_loc_entity_ind_from_entity_glb_id(tId,EntityRank::ELEMENT);
//
//    moris_index tChildMeshIndex = aModel.get_background_mesh().child_mesh_index(tIndex,EntityRank::ELEMENT);
//
//    Child_Mesh & tCM = aModel.get_cut_mesh().get_child_mesh(tChildMeshIndex);
//
//    moris::Matrix<moris::IndexMat> const & tSubphaseIds = tCM.get_subphase_ids();
//
//    if(par_rank() == 0)
//    {
//        moris::print(tSubphaseIds,"tSubphaseIds");
//    }
//
//    barrier();
//
//    if(par_rank() == 1)
//    {
//        moris::print(tSubphaseIds,"tSubphaseIds");
//    }

    return true;
}


inline
void
collect_enrichment_data(Matrix<IndexMat>       const & aEnrichedBasisLocToGlobal,
                        Cell<Matrix<IndexMat>> const & aSubphaseIdsInSupportOfEnrichedBasis,
                        Matrix<IndexMat>             & aAllEnrichedBasisLocToGlobal,
                        Cell<Matrix<IndexMat>>       & aAllSubphaseIdsInSupportOfEnrichedBasis)
{
    // concatenate into a vector
    Matrix<IndexMat> tSubphaseIdsInEnrichedBasisSuppData;
    Matrix<IndexMat> tSubphaseIdsInEnrichedBasisSuppOffsets;
    convert_vector_to_sparse(aSubphaseIdsInSupportOfEnrichedBasis,tSubphaseIdsInEnrichedBasisSuppData,tSubphaseIdsInEnrichedBasisSuppOffsets);

    // send data to proc 0
    if(par_rank() > 0)
    {
        // send map
        nonblocking_send(aEnrichedBasisLocToGlobal,aEnrichedBasisLocToGlobal.n_rows(),aEnrichedBasisLocToGlobal.n_cols(),0,10);

        // send sparse form of subphase ids
        nonblocking_send(tSubphaseIdsInEnrichedBasisSuppData,tSubphaseIdsInEnrichedBasisSuppData.n_rows(),tSubphaseIdsInEnrichedBasisSuppData.n_cols(),0,11);
        nonblocking_send(tSubphaseIdsInEnrichedBasisSuppOffsets,tSubphaseIdsInEnrichedBasisSuppOffsets.n_rows(),tSubphaseIdsInEnrichedBasisSuppOffsets.n_cols(),0,12);
    }

    barrier();

    if(par_rank() == 0)
    {
        moris::Cell<moris::Matrix<IndexMat>> tAllEnrichedBasisLocToGlobal(par_size(),moris::Matrix<IndexMat>(1,1));
        moris::Cell<moris::Matrix<IndexMat>> tAllSubphaseIdsInEnrichedBasisSuppData(par_size(),moris::Matrix<IndexMat>(1,1));
        moris::Cell<moris::Matrix<IndexMat>> tAllSubphaseIdsInEnrichedBasisSuppOffsets(par_size(),moris::Matrix<IndexMat>(1,1));

        // number of enriched basis function
        moris::uint tNumEnrichedBasisFuncs = aEnrichedBasisLocToGlobal.numel();

        // go and receive all processors node coords and maps
        for(moris::uint i= 1; i <(uint)par_size(); i++)
        {
            receive(tAllEnrichedBasisLocToGlobal(i), 1, i,10);
            receive(tAllSubphaseIdsInEnrichedBasisSuppData(i), 1, i,11);
            receive(tAllSubphaseIdsInEnrichedBasisSuppOffsets(i), 1, i,12);

            tNumEnrichedBasisFuncs = tNumEnrichedBasisFuncs + tAllEnrichedBasisLocToGlobal(i).numel();
        }

        // size return data
        aAllEnrichedBasisLocToGlobal.resize(1,tNumEnrichedBasisFuncs);

        // add my data to the all list
        aAllSubphaseIdsInSupportOfEnrichedBasis.append(aSubphaseIdsInSupportOfEnrichedBasis);
        aAllEnrichedBasisLocToGlobal({0,0},{0,aEnrichedBasisLocToGlobal.numel()-1}) = aEnrichedBasisLocToGlobal.matrix_data();

        moris::uint tStart = aEnrichedBasisLocToGlobal.numel();
        moris::uint tEnd = 0;

        // sorted matrix
        Matrix<IndexMat> tSorted;

        for(moris::uint i= 0; i <(uint)par_size(); i++)
        {
            tEnd = tStart + tAllEnrichedBasisLocToGlobal(i).numel() - 1;

            if(tStart != tEnd)
            {
                aAllEnrichedBasisLocToGlobal({0,0},{tStart,tEnd}) = tAllEnrichedBasisLocToGlobal(i).matrix_data();
                Cell<Matrix<IndexMat>> tProcsSubphaseIdsInSupport(0);
                convert_sparse_to_vector(tAllSubphaseIdsInEnrichedBasisSuppData(i), tAllSubphaseIdsInEnrichedBasisSuppOffsets(i), tProcsSubphaseIdsInSupport);

                // iterate through and sort
                for(auto it:tProcsSubphaseIdsInSupport)
                {
                    moris::sort(it, tSorted, "ascend", 1 );

                    aAllSubphaseIdsInSupportOfEnrichedBasis.push_back(tSorted);
                }

                tStart = tEnd+1;
            }
        }
    }
}

inline
bool
verify_subphases_in_basis_support(Matrix<IndexMat>       const & aAllEnrichedBasisLocToGlobal,
                                  Cell<Matrix<IndexMat>> const & aAllSubphaseIdsInSupportOfEnrichedBasis)
{
    // INDEX FOR THIS FUNCTION = ID - 1

    bool tValid = true;

    if(par_rank() == 0)
    {
        moris_id tMaxId = aAllEnrichedBasisLocToGlobal.max();

        // assign owner to cell which has the most interpolation cells in support (indicating it has the entire support)
        moris::Cell<moris::uint> tOwnerIndices(tMaxId+1,MORIS_UINT_MAX);
        moris::Cell<moris::Cell<moris::moris_index>> tBasisIdToIndicesInLocToGlobalMap(tMaxId+1,0);

        // iterate through loc to global basis map
        for( moris::uint i = 0 ; i < aAllEnrichedBasisLocToGlobal.numel(); i++)
        {
            moris_id tId = aAllEnrichedBasisLocToGlobal(i);
            tBasisIdToIndicesInLocToGlobalMap(tId-1).push_back(i);
        }


        //asssign owners
        for(moris::uint i = 0; i <tBasisIdToIndicesInLocToGlobalMap.size(); i++)
        {
            moris::uint tMaxSupportNumInSupport = 0;
            for(moris::uint j = 0; j < tBasisIdToIndicesInLocToGlobalMap(i).size(); j++)
            {
                moris_index tBasisIndexInCell = tBasisIdToIndicesInLocToGlobalMap(i)(j);
                if(aAllSubphaseIdsInSupportOfEnrichedBasis(tBasisIndexInCell).numel() > tMaxSupportNumInSupport)
                {
                    moris_id tId = aAllEnrichedBasisLocToGlobal(tBasisIndexInCell);

                    tOwnerIndices(tId - 1) = tBasisIndexInCell;

                    tMaxSupportNumInSupport = aAllSubphaseIdsInSupportOfEnrichedBasis(tBasisIndexInCell).numel();
                }
            }
        }

        for( moris::uint i = 0 ; i < aAllEnrichedBasisLocToGlobal.numel(); i++)
        {
            moris_index tBasisId = aAllEnrichedBasisLocToGlobal(i);
            moris_index tOwnerIndex = tOwnerIndices(tBasisId-1);

            if(tOwnerIndex != (moris_index)i)
            {
                CHECK(aAllSubphaseIdsInSupportOfEnrichedBasis(tOwnerIndex).numel() >= aAllSubphaseIdsInSupportOfEnrichedBasis(i).numel());

                // iterate through my subphase ids in this basis support
                for(moris::uint  iMySp = 0; iMySp < aAllSubphaseIdsInSupportOfEnrichedBasis(i).numel(); iMySp++)
                {
                    moris_index tMySubphaseId = aAllSubphaseIdsInSupportOfEnrichedBasis(i)(iMySp);

                    bool tFoundIt = false;
                    // iterate through owner and find the desi
                    for(moris::uint iOwnSp = 0; iOwnSp < aAllSubphaseIdsInSupportOfEnrichedBasis(tOwnerIndex).numel(); iOwnSp++)
                    {
                        if(aAllSubphaseIdsInSupportOfEnrichedBasis(tOwnerIndex)(iOwnSp) == tMySubphaseId)
                        {
                            tFoundIt = true;
                        }
                    }

                    // verify we found it
                    CHECK(tFoundIt);
                }

            }
        }
    }
    barrier();

    return tValid;
}

inline
bool
collect_parallel_basis_functions_and_support(Model & aModel,
                                             moris::Matrix<moris::IndexMat> & aBasisToEnrichedBasisIds,
                                             moris::Matrix<moris::IndexMat> & aBasisToEnrichedBasisOffsets,
                                             moris::Matrix<moris::IndexMat> & aBasisIds,
                                             moris::Matrix<moris::IndexMat> & aCellsInEnrichedBasis,
                                             moris::Matrix<moris::IndexMat> & aCellsInEnrichedBasisOffsets,
                                             moris::Matrix<moris::IndexMat> & aEnrichedBasisIds)
{
    // data structures
    Cut_Mesh & tCutMesh = aModel.get_cut_mesh();
//    Background_Mesh & tBackgroundMesh = aModel.get_background_mesh();
    Enriched_Interpolation_Mesh const & tEnrInterpMesh = aModel.get_enriched_interp_mesh();
    Enrichment const & tEnrichment = aModel.get_basis_enrichment();

    // get the enriched coefficient local to global ids
    moris::Matrix<moris::IdMat> const & tLocToGlobalEnrCoeffs = tEnrInterpMesh.get_enriched_coefficient_local_to_global_map();

    // get the background coefficient local to global ids
    moris::Matrix<moris::IdMat> tLocToGlobalBackgroundCoeffs = tEnrInterpMesh.get_background_coefficient_local_to_global_map();

    // get the background coefficient to enriched coefficient data
    Cell<moris::Matrix<moris::IdMat>> const & tBaseCoeffToEnrichedCoeffInds = tEnrInterpMesh.get_enriched_coefficients_to_background_coefficients();

    // get the subphases ids in support of enrichde basis function
    Cell<moris::Matrix< moris::IndexMat >> tSubphaseIdsInEnrichedBasisSupp = tEnrichment.get_subphases_glb_id_in_enriched_basis();


    Cell<Matrix<IndexMat>> tAllSubphaseIdsInSupportOfEnrichedBasis(0);
    Matrix<IndexMat>       tAllEnrichedBasisLocToGlobal(0,0);
    collect_enrichment_data(tLocToGlobalEnrCoeffs, tSubphaseIdsInEnrichedBasisSupp, tAllEnrichedBasisLocToGlobal, tAllSubphaseIdsInSupportOfEnrichedBasis);

    bool tFlag = verify_subphases_in_basis_support( tAllEnrichedBasisLocToGlobal, tAllSubphaseIdsInSupportOfEnrichedBasis );

    return tFlag;
}

inline
bool
verify_parallel_basis_functions(Model & aModel)
{
    moris::Matrix<moris::IndexMat> aBasisToEnrichedBasisIds;
    moris::Matrix<moris::IndexMat> aBasisToEnrichedBasisOffsets;
    moris::Matrix<moris::IndexMat> aBasisIds;
    moris::Matrix<moris::IndexMat> aCellsInEnrichedBasis;
    moris::Matrix<moris::IndexMat> aCellsInEnrichedBasisOffsets;
    moris::Matrix<moris::IndexMat> aEnrichedBasisIds;

    bool tFlag = collect_parallel_basis_functions_and_support(aModel,aBasisToEnrichedBasisIds,aBasisToEnrichedBasisOffsets,aBasisIds,aCellsInEnrichedBasis,aCellsInEnrichedBasisOffsets,aEnrichedBasisIds);


    return tFlag;
}

TEST_CASE("Regular Subdivision Method Parallel","[REG_SUB_PARALLEL]")
{
    //     Geometry Engine Setup -----------------------
    //     Using a Levelset Sphere as the Geometry
    real tRadius = 0.65;
    real tXCenter = 1.0;
    real tYCenter = 1.0;
    real tZCenter = 2.0;
    moris::ge::Sphere  tLevelsetSphere(tRadius, tXCenter, tYCenter, tZCenter);
    moris::ge::GEN_Phase_Table tPhaseTable (1,  Phase_Table_Structure::EXP_BASE_2);
    moris::ge::GEN_Geometry_Engine tGeometryEngine(tLevelsetSphere,tPhaseTable);

    // Create Mesh ---------------------------------
    std::string tMeshFileName = "generated:1x2x4";
    moris::mtk::Interpolation_Mesh* tMeshData = moris::mtk::create_interpolation_mesh( MeshType::STK, tMeshFileName);

    // Setup XTK Model -----------------------------
    size_t tModelDimension = 3;
    Model tXTKModel(tModelDimension,tMeshData,tGeometryEngine);

    //Specify your decomposition methods and start cutting
    Cell<enum Subdivision_Method> tDecompositionMethods = {Subdivision_Method::NC_REGULAR_SUBDIVISION_HEX8};
    tXTKModel.decompose(tDecompositionMethods);

    Cut_Mesh const & tCutMesh = tXTKModel.get_cut_mesh();

    // check nodes
    CHECK(verify_parallel_node_coordinates(tXTKModel));

    // check cells
    CHECK(verify_parallel_elements(tXTKModel));

    delete tMeshData;

}
TEST_CASE("Regular Subdivision and Node Hierarchy Method Parallel","[CONF_PARALLEL]")
{
    //     Geometry Engine Setup -----------------------
    //     Using a Levelset Sphere as the Geometry
    real tRadius = 0.25;
    real tXCenter = 1.0;
    real tYCenter = 1.0;
    real tZCenter = 2.0;
    moris::ge::Sphere tLevelsetSphere(tRadius, tXCenter, tYCenter, tZCenter);
    moris::ge::GEN_Phase_Table tPhaseTable (1,  Phase_Table_Structure::EXP_BASE_2);
    moris::ge::GEN_Geometry_Engine tGeometryEngine(tLevelsetSphere,tPhaseTable);

    // Create Mesh ---------------------------------
    std::string tMeshFileName = "generated:1x1x4";
    moris::mtk::Interpolation_Mesh* tMeshData = moris::mtk::create_interpolation_mesh( MeshType::STK, tMeshFileName);

    // Setup XTK Model -----------------------------
    size_t tModelDimension = 3;
    Model tXTKModel(tModelDimension,tMeshData,tGeometryEngine);
    tXTKModel.mVerbose = false;

    //Specify your decomposition methods and start cutting
    Cell<enum Subdivision_Method> tDecompositionMethods = {Subdivision_Method::NC_REGULAR_SUBDIVISION_HEX8,Subdivision_Method::C_HIERARCHY_TET4};
    tXTKModel.decompose(tDecompositionMethods);

    Cut_Mesh const & tCutMesh = tXTKModel.get_cut_mesh();

    CHECK(verify_parallel_node_coordinates(tXTKModel));
    CHECK(verify_parallel_elements(tXTKModel));
    CHECK(verify_parallel_subphase_ids(tXTKModel));

    // perform basis enrichment
    tXTKModel.perform_basis_enrichment(EntityRank::NODE);

    CHECK(verify_parallel_basis_functions(tXTKModel));

    // setup output mesh options with cell enrichment fields
    Output_Options tOutputOptions;
//    tOutputOptions.mRealElementExternalFieldNames = tEnrichmentFieldNames;
//    tOutputOptions.mAddParallelFields = true;

    moris::mtk::Mesh* tCutMeshData = tXTKModel.get_output_mesh(tOutputOptions);
//    tEnrichment.write_cell_enrichment_to_fields(tEnrichmentFieldNames,tCutMeshData);
    std::string tMeshOutputFile ="./xtk_exo/conformalparallel.e";
    tCutMeshData->create_output_mesh(tMeshOutputFile);


    delete tMeshData;
    delete tCutMeshData;
}
}

