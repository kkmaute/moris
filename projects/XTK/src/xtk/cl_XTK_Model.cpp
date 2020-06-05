/*
 * cl_XTK_Model.cpp
 *
 *  Created on: Feb 18, 2019
 *      Author: doble
 */

#include "cl_XTK_Model.hpp"
#include "cl_XTK_Background_Mesh.hpp"
#include "cl_MTK_Integration_Mesh.hpp"
#include "cl_MTK_Interpolation_Mesh.hpp"
#include "cl_XTK_Enriched_Interpolation_Mesh.hpp"
#include "cl_XTK_Enriched_Integration_Mesh.hpp"
#include "cl_XTK_Ghost_Stabilization.hpp"
#include "cl_XTK_Multigrid.hpp"
#include "fn_all_true.hpp"
#include "fn_unique.hpp"
#include "op_equal_equal.hpp"
#include "fn_sort.hpp"
#include "fn_iscol.hpp"
#include "fn_trans.hpp"
#include "fn_equal_to.hpp"
#include "fn_generate_element_to_element.hpp"
#include "fn_create_faces_from_element_to_node.hpp"
#include "fn_create_edges_from_element_to_node.hpp"
#include "HDF5_Tools.hpp"
#include "cl_MTK_Visualization_STK.hpp"
#include "cl_MTK_Cell_Info_Factory.hpp"
#include "cl_MTK_Cell_Info.hpp"
#include "cl_MTK_Writer_Exodus.hpp"
#include "fn_Parsing_Tools.hpp"
//#include "cl_XTK_Enrichment.hpp"
using namespace moris;

namespace xtk
{
// ----------------------------------------------------------------------------------
// Constructor/Deconstructor Source code
// ----------------------------------------------------------------------------------
Model::~Model()
{
    if(mEnrichment         != nullptr ) { delete mEnrichment;         }
    if(mGhostStabilization != nullptr ) { delete mGhostStabilization; }


    for(auto tIt:mEnrichedInterpMesh)
    {
        if(tIt !=nullptr)
        {
            delete tIt;
        }
    }
    mEnrichedInterpMesh.clear();

    for(auto tIt:mEnrichedIntegMesh)
    {
        if(tIt !=nullptr)
        {
            delete tIt;
        }
    }
    mEnrichedIntegMesh.clear();

}

/*
 * using the general geometry engine
 */
Model::Model(uint aModelDimension,
             moris::mtk::Interpolation_Mesh* aMeshData,
             moris::ge::Geometry_Engine* aGeometryEngine,
             bool aLinkGeometryOnConstruction ) :
                                          mSameMesh(false),
                                          mModelDimension(aModelDimension),
                                          mBackgroundMesh(aMeshData,aGeometryEngine),
                                          mCutMesh(this,mModelDimension),
                                          mGeometryEngine(aGeometryEngine),
                                          mEnrichment(nullptr),
                                          mGhostStabilization(nullptr),
                                          mEnrichedInterpMesh(0,nullptr),
                                          mEnrichedIntegMesh(0,nullptr),
                                          mConvertedToTet10s(false)
{
    // flag this as a non-parameter list based run
    mParameterList.insert("has_parameter_list", false);

    if(aLinkGeometryOnConstruction == true)
    {
        link_background_mesh_to_geometry_objects();
    }

    mBackgroundMesh.initialize_interface_node_flags(mBackgroundMesh.get_num_entities(EntityRank::NODE),mGeometryEngine->get_num_geometries());
}

Model::Model( moris::ParameterList const & aParameterList):
        mParameterList(aParameterList)
{
    // flag this as a paramter list based run
    mParameterList.insert("has_parameter_list", true);
}

void
Model::set_geometry_engine(moris::ge::Geometry_Engine* aGeometryEngine)
{
    mGeometryEngine = aGeometryEngine;
}

void
Model::set_mtk_background_mesh(moris::mtk::Interpolation_Mesh* aMesh)
{
    this->initialize( aMesh );

    mInitializeCalled = true;
}

void
Model::set_input_performer( std::shared_ptr< mtk::Mesh_Manager > aMTKPerformer )
{
	mMTKInputPerformer = aMTKPerformer;
}

void
Model::set_output_performer( std::shared_ptr< mtk::Mesh_Manager > aMTKPerformer )
{
	mMTKOutputPerformer = aMTKPerformer;
}

void
Model::initialize( moris::mtk::Interpolation_Mesh* aMesh )
{
    mSameMesh           = false;
    mModelDimension     = aMesh->get_spatial_dim();
    mCutMesh            = Cut_Mesh(this,mModelDimension),
    mEnrichment         = nullptr;
    mGhostStabilization = nullptr;
    mEnrichedInterpMesh = Cell<Enriched_Interpolation_Mesh*>(0, nullptr);
    mEnrichedIntegMesh  = Cell<Enriched_Integration_Mesh*>(0, nullptr);
    mConvertedToTet10s  = false;

    mBackgroundMesh = Background_Mesh(aMesh,mGeometryEngine);
    link_background_mesh_to_geometry_objects();
    mBackgroundMesh.initialize_interface_node_flags(mBackgroundMesh.get_num_entities(EntityRank::NODE),mGeometryEngine->get_num_geometries());
}

void
Model::perform()
{
    if( !mInitializeCalled )
    {
        MORIS_ERROR( mMTKInputPerformer != nullptr ,"xtk::Model::perform(), mMTKInputPerformer not set!");

        //FIXME hardcodes to mesh pair index 0
        moris::mtk::Interpolation_Mesh* tMesh = mMTKInputPerformer->get_interpolation_mesh( 0 );

        this->initialize( tMesh );
    }

    MORIS_ASSERT(this->has_parameter_list(),"Perform can only be called on a parameter list based XTK");
    MORIS_ERROR(this->valid_parameters(),"Invalid parameters detected in XTK.");

    mVerbose = mParameterList.get<bool>("verbose");

    if(mParameterList.get<bool>("decompose"))
    {
        Cell<enum Subdivision_Method> tSubdivisionMethods = this->get_subdivision_methods();
        this->decompose(tSubdivisionMethods);
    }

    if(mParameterList.get<bool>("enrich"))
    {
        enum EntityRank tBasisRank = get_entity_rank_from_str(mParameterList.get<std::string>("basis_rank"));

        Matrix<IndexMat> tMeshIndexCell;
        moris::string_to_mat(mParameterList.get< std::string >( "enrich_mesh_indices" ), tMeshIndexCell);

        this->perform_basis_enrichment(tBasisRank,tMeshIndexCell);

        // if high to low double side sets need to be created
        if(mParameterList.get<bool>("high_to_low_dbl_side_sets"))
        {
            for(moris::uint i = 0; i < mGeometryEngine->get_num_bulk_phase(); i++)
            {
                for(moris::uint j = 0; j < mGeometryEngine->get_num_bulk_phase(); j++)
                {
                    if(i > j)
                    {
                        mEnrichedIntegMesh(0)->create_dbl_sided_interface_set( i, j );
                    }
                }
            }
        }

    }

    if(mParameterList.get<bool>("ghost_stab"))
    {
        this->construct_face_oriented_ghost_penalization_cells();
    }

    if( mParameterList.get<bool>("multigrid") )
    {
        this->construct_multigrid();
    }

    // get meshes
    xtk::Enriched_Interpolation_Mesh & tEnrInterpMesh = this->get_enriched_interp_mesh();
    xtk::Enriched_Integration_Mesh   & tEnrIntegMesh  = this->get_enriched_integ_mesh();

    // place the pair in mesh manager
    mMTKOutputPerformer->register_mesh_pair( &tEnrInterpMesh, &tEnrIntegMesh );

    if( mParameterList.get<bool>("print_enriched_ig_mesh") )
    {
        tEnrIntegMesh.print();
    }

    if( mParameterList.get<bool>("exodus_output_XTK_ig_mesh") )
    {
//        tEnrIntegMesh.deactivate_empty_sets();

        // Write mesh
        moris::mtk::Writer_Exodus writer( &tEnrIntegMesh );
        writer.write_mesh("", "./xtk_temp.exo");

        // Write the fields
        writer.set_time(0.0);
        writer.close_file();
    }
}

bool
Model::has_parameter_list()
{
    return mParameterList.get<bool>("has_parameter_list");
}

bool
Model::valid_parameters()
{
    bool tDecompose = mParameterList.get<bool>("decompose");
    bool tEnrich    = mParameterList.get<bool>("enrich");
    bool tGhost     = mParameterList.get<bool>("ghost_stab");
    bool tMultigrid     = mParameterList.get<bool>("multigrid");

    if(tEnrich == true)
    {
        MORIS_ERROR(tDecompose, "To perform basis enrichment, decomposition is also required.");
    }

    if(tGhost == true)
    {
        MORIS_ERROR(tDecompose && tEnrich, "To perform ghost stabilization, decomposition and enrichment are also required.");
    }

    if(tMultigrid == true)
    {
        MORIS_ERROR(tDecompose && tEnrich, "To perform multigrid, decomposition and enrichment are also required.");
    }

    return true;
}

Cell<enum Subdivision_Method>
Model::get_subdivision_methods()
{
    MORIS_ASSERT(this->has_parameter_list(),"Perform can only be called on a parameter list based XTK");


    moris::uint       tSpatialDimension = this->get_spatial_dim();
    enum CellTopology tBGCellTopo       = mBackgroundMesh.get_parent_cell_topology();
    std::string       tDecompStr        = mParameterList.get<std::string>("decomposition_type");

    // determine if we are going conformal or not
    bool tConformal = true;
    if(tDecompStr.compare("conformal") == 0 )
    {
        tConformal = true;
    }
    else if(tDecompStr.compare("nonconformal") == 0 )
    {
        tConformal = true;
    }
    else
    {
        MORIS_ERROR(0,"Invalid decomposition_type provided. Recognized Options: Conformal and Nonconformal");
    }


    if(tSpatialDimension == 2 )
    {
        if(tBGCellTopo == CellTopology::QUAD4  && tConformal)
        {
            return {Subdivision_Method::NC_REGULAR_SUBDIVISION_QUAD4, Subdivision_Method::C_TRI3};
        }
        else if(tBGCellTopo == CellTopology::QUAD4  && !tConformal)
        {
            return {Subdivision_Method::NC_REGULAR_SUBDIVISION_QUAD4};
        }
    }
    else if ( tSpatialDimension == 3 )
    {
        if(tBGCellTopo == CellTopology::HEX8  && tConformal)
        {
            return {Subdivision_Method::NC_REGULAR_SUBDIVISION_HEX8, Subdivision_Method::C_HIERARCHY_TET4};
        }
        else if(tBGCellTopo == CellTopology::HEX8  && !tConformal)
        {
            return {Subdivision_Method::NC_REGULAR_SUBDIVISION_HEX8};
        }
        else if(tBGCellTopo == CellTopology::TET4  && tConformal)
        {
            return {Subdivision_Method::C_HIERARCHY_TET4};
        }
    }
    else
    {
        MORIS_ASSERT(0,"Invalid spatial dimension");
    }

    MORIS_ERROR(0,"Failed determining subdivision methods");

    return Cell<enum Subdivision_Method>(0);

}


// ----------------------------------------------------------------------------------
// Decomposition Source code
// ----------------------------------------------------------------------------------
void
Model::decompose(Cell<enum Subdivision_Method> aMethods)
{
    // Start clock
    std::clock_t tTotalTime = std::clock();

    // Assert that there has been a link between geometry model and background mesh
    MORIS_ERROR(mLinkedBackground, "Geometry model and background mesh have not been linked via call to link_background_mesh_to_geometry_objects");

    // Process for a decomposition
    uint tNumDecompositions = aMethods.size();
    uint tNumGeometries = mGeometryEngine->get_num_geometries();

    print_decompsition_preamble(aMethods);

    // Tell the subdivision to assign node Ids if it is the only subdivision method (critical for outputting)
    // This is usually only going to happen in test cases
    // Note: the Conformal subdivision methods dependent on node ids for subdivision routine, the node Ids are set regardless of the below boolean

    bool tNonConformingMeshFlag = false;
    bool tSetPhase = true;
    if(aMethods.size() == 1)
    {
        tNonConformingMeshFlag = true;
    }

    // Loop over each geometry and have an active child mesh indices list for each
    for(moris::size_t iGeom = 0; iGeom<tNumGeometries; iGeom++)
    {
        bool tFirstSubdivisionFlag = true;
        moris::Matrix< moris::IndexMat > tActiveChildMeshIndices(1,1,0);

        for (moris::size_t iDecomp = 0; iDecomp < tNumDecompositions; iDecomp++)
        {
            // start timing on this decomposition
            std::clock_t start = std::clock();

            // Perform subdivision
            this->decompose_internal(aMethods(iDecomp), iGeom, tActiveChildMeshIndices, tFirstSubdivisionFlag, tNonConformingMeshFlag);

            // Change the first subdivision flag as false
            tFirstSubdivisionFlag = false;

            // print timing
            if(moris::par_rank() == 0 && mVerbose)
            {
                std::cout<<"XTK: Decomposition "<<get_enum_str(aMethods(iDecomp))<<" for geometry "<<iGeom<< " completed in " <<(std::clock() - start) / (double)(CLOCKS_PER_SEC)<<" s."<<std::endl;
                std::cout<<"XTK: Decomposition "<<get_enum_str(aMethods(iDecomp))<<" for geometry "<<iGeom<< " had "<<  tActiveChildMeshIndices.numel()<<" intersected background elements."<<std::endl;
            }
        }
        // If it's not the last geometry tell the geometry engine we're moving on
        if(iGeom!= tNumGeometries-1)
        {
            mGeometryEngine->advance_geometry_index();
        }
    }

    // Tell the xtk mesh to set all necessary information to finalize decomposition allowing
    // i.e set element ids, indices for children elements
    this->finalize_decomp_in_xtk_mesh(tSetPhase);

    if(moris::par_rank() == 0 && mVerbose)
    {
        std::cout<<"XTK: Decomposition completed in " <<(std::clock() - tTotalTime) / (double)(CLOCKS_PER_SEC)<<" s."<<std::endl;
    }
}

void
Model::decompose_internal(enum Subdivision_Method    const & aSubdivisionMethod,
                          moris::uint                        aGeomIndex,
                          moris::Matrix< moris::IndexMat > & aActiveChildMeshIndices,
                          bool const &                       aFirstSubdivision,
                          bool const &                       aSetIds)
{
    switch (aSubdivisionMethod)
    {
        case (Subdivision_Method::NC_REGULAR_SUBDIVISION_HEX8):
                {
            //            MORIS_ASSERT(tXTKMeshData.get_entity_connected_to_entity_loc_inds(0, moris::EntityRank::ELEMENT, moris::EntityRank::NODE).numel() == 8, "NC_REGULAR_SUBDIVISION_HEX8 is for HEX8 meshes only.");
            MORIS_ASSERT(aFirstSubdivision,"NC_REGULAR_SUBDIVISION_HEX8 needs to be the first subdivision routine for each geometry");
            MORIS_ASSERT(mModelDimension == 3,"NC_REGULAR_SUBDIVISION_HEX8 needs to be done on a 3D mesh");

            // Runs the first cut routine to get the new active child mesh indices and indicate which are new and need to be regularly subdivided and which ones dont
            moris::Matrix< moris::IndexMat > tNewPairBool;
            run_first_cut_routine(TemplateType::HEX_8, aGeomIndex, 8,  aActiveChildMeshIndices,tNewPairBool);

            // set the child cell topology as tet 3s
            mCutMesh.set_child_element_topology(CellTopology::TET4);

            // initialize a struct of all the data we are keeping track of in this decomposition
            // intended to reduce the clutter of function inputs etc
            Decomposition_Data tDecompData;
            tDecompData.mSubdivisionMethod = Subdivision_Method::NC_REGULAR_SUBDIVISION_HEX8;

            // number of intersected elements
            moris::uint tIntersectedCount = aActiveChildMeshIndices.n_cols();

            // make node requests for each intersected element
            this->decompose_internal_reg_sub_hex8_make_requests(aActiveChildMeshIndices,tNewPairBool,tDecompData);

            // specify a dummy secondary id (not really needed for this type of decomposition)
            tDecompData.tSecondaryIdentifiers = Cell<moris_index>(tDecompData.tNewNodeParentIndex.size(), MORIS_INDEX_MAX);

            moris_index tMessageTag = 60000; /*arbitrary tag for regular subdivision*/
            assign_node_requests_identifiers(tDecompData,tMessageTag);

            // Allocate interface flag space in XTK mesh even though these are not interface nodes
            mBackgroundMesh.allocate_space_in_interface_node_flags(tDecompData.tNewNodeIndex.size(),mGeometryEngine->get_num_geometries());

            // add nodes to the background mesh
            mBackgroundMesh.batch_create_new_nodes(tDecompData.tNewNodeId,tDecompData.tNewNodeIndex, tDecompData.tNewNodeOwner,tDecompData.tNewNodeCoordinate);

            // add nodes to child mesh
            this->decompose_internal_set_new_nodes_in_child_mesh_reg_sub(aActiveChildMeshIndices,tNewPairBool, 3, tDecompData);

            // associate new nodes with geometry objects
            create_new_node_association_with_geometry(tDecompData);

            for(moris::size_t i = 0; i< tIntersectedCount; i++)
            {
                if(tNewPairBool(0,i) == 0)
                {
                    mCutMesh.generate_templated_mesh(aActiveChildMeshIndices(i),TemplateType::REGULAR_SUBDIVISION_HEX8);
                }
            }

            break;
                }
        case (Subdivision_Method::NC_REGULAR_SUBDIVISION_QUAD4):
            {
            MORIS_ASSERT(aFirstSubdivision, "NC_REGULAR_SUBDIVISION_QUAD4 needs to be the first subdivision routine for each geometry.");
            MORIS_ASSERT(mModelDimension == 2, "NC_REGULAR_SUBDIVISION_QUAD4 needs to be done on a 2D mesh.");

            // Runs the first cut routine to get the new active child mesh indices and indicate which are new and need to be regularly subdivided and which ones don't
            moris::Matrix< moris::IndexMat > tNewPairBool;
            run_first_cut_routine(TemplateType::QUAD_4, aGeomIndex, 4, aActiveChildMeshIndices, tNewPairBool);

            // mark child cells as tri 3s
            mCutMesh.set_child_element_topology(CellTopology::TRI3);

            // initialize a struct of all the data we are keeping track of in this decomposition
            // intended to reduce the clutter of function inputs etc
            Decomposition_Data tDecompData;
            tDecompData.mSubdivisionMethod = Subdivision_Method::NC_REGULAR_SUBDIVISION_QUAD4;

            // number of intersected elements
            moris::uint tIntersectedCount = aActiveChildMeshIndices.n_cols();

            // make node requests for each intersected element
            this->decompose_internal_reg_sub_quad4_make_requests(aActiveChildMeshIndices, tNewPairBool, tDecompData);

            // specify a dummy secondary id (not really needed for this type of decomposition)
            tDecompData.tSecondaryIdentifiers = Cell<moris_index>(tDecompData.tNewNodeParentIndex.size(), MORIS_INDEX_MAX);

            moris_index tMessageTag = 60000; /*arbitrary tag for regular subdivision*/
            assign_node_requests_identifiers(tDecompData,tMessageTag);

            // Allocate interface flag space in XTK mesh even though these are not interface nodes
            mBackgroundMesh.allocate_space_in_interface_node_flags(tDecompData.tNewNodeIndex.size(),mGeometryEngine->get_num_geometries());

            // add nodes to the background mesh
            mBackgroundMesh.batch_create_new_nodes(tDecompData.tNewNodeId,tDecompData.tNewNodeIndex, tDecompData.tNewNodeOwner,tDecompData.tNewNodeCoordinate);

            // crate nodes in child mesh
            this->decompose_internal_set_new_nodes_in_child_mesh_reg_sub(aActiveChildMeshIndices,tNewPairBool,2,tDecompData);

            // associate new nodes with geometry objects
            create_new_node_association_with_geometry(tDecompData);


            for(moris::size_t i = 0; i< tIntersectedCount; i++)
            {
                if(tNewPairBool(0,i) == 0)
                {
                    mCutMesh.generate_templated_mesh(aActiveChildMeshIndices(i),TemplateType::REGULAR_SUBDIVISION_QUAD4);
                }
            }
            break;

            }
        case (Subdivision_Method::C_HIERARCHY_TET4):
                {

            // If it the first subdivision we need to find the intersected before placing the conformal nodes
            // Intersected elements are flagged via the Geometry_Engine
            if(aFirstSubdivision)
            {
                moris::Matrix< moris::IndexMat > tNewPairBool;
                run_first_cut_routine(TemplateType::TET_4, aGeomIndex, 4, aActiveChildMeshIndices,tNewPairBool);

                for(moris::size_t i = 0; i<aActiveChildMeshIndices.n_cols(); i++)
                {
                    Child_Mesh & tChildMesh = mCutMesh.get_child_mesh(aActiveChildMeshIndices(0,i));
                    tChildMesh.generate_connectivities(true,true,true);
                }

                // set the child cell topology as tet 4s
                mCutMesh.set_child_element_topology(CellTopology::TET4);

            }

            // For hex background meshes we have a three dimension parametric coordinate
            moris::size_t tDimParamCoord = 3;

            // For tet background meshes we have a 4-d parametric coordinate
            if(aFirstSubdivision)
            {
                tDimParamCoord = 4;
            }

            // initialize a struct of all the data we are keeping track of in this decomposition
            // intended to reduce the clutter of function inputs etc
            Decomposition_Data tDecompData;
            tDecompData.mSubdivisionMethod = Subdivision_Method::C_HIERARCHY_TET4;
            tDecompData.mConformalDecomp = true;
            tDecompData.mHasSecondaryIdentifier = true;
            tDecompData.mFirstSubdivision = aFirstSubdivision;

            // Initialize
            moris::size_t tEdgeInd = MORIS_INDEX_MAX;

            // Initialize topologies used in this method (all local coordinates are with respect to an edge)
            Edge_Topology tEdgeTopology;

            // initialize a couple of commonly used matrices in this method
            moris::Matrix< moris::DDRMat > tLocalCoordRelativeToEdge(1,1, 0); // ALong an edge
            moris::Matrix< moris::DDRMat > tGlobalCoord(1,3, 0); // ALong an edge
            moris::Matrix< moris::DDRMat > tEdgeNodeParamCoordinates(2,tDimParamCoord); // parametric coordinate of end nodes wrt parent element

            // Check type specified as conformal (could change this to enum)
            moris::size_t tCheckType = 1;
            moris::Matrix< moris::DDRMat > tNodeCoords = mBackgroundMesh.get_all_node_coordinates_loc_inds();

            // get the underlying background mesh data
            moris::mtk::Mesh & tMeshData = mBackgroundMesh.get_mesh_data();

            // resize child mesh to new node information
            tDecompData.tCMNewNodeLoc.resize(aActiveChildMeshIndices.n_cols());
            tDecompData.tCMNewNodeParamCoord.resize(aActiveChildMeshIndices.n_cols());

            for (moris::size_t j = 0; j < aActiveChildMeshIndices.n_cols(); j++)
            {

                // Initialize geometry objects
                Cell<moris::ge::GEN_Geometry_Object> tGeoObjects;

                // Get the child mesh that is active
                Child_Mesh & tChildMesh = mCutMesh.get_child_mesh(aActiveChildMeshIndices(0,j));

                // edge to node connectivity from child mesh
                moris::Matrix< moris::IndexMat > const & tEdgeToNode = tChildMesh.get_edge_to_node();

                // Ask geometry engine which edges are intersected (Simple mesh local indexed edges)
                mGeometryEngine->is_intersected(tNodeCoords, tEdgeToNode, tCheckType, tGeoObjects);

                // Initialize node index pointers based on number of intersected edges and parametric coordinates
                uint tNumNewNodes = 0;
                Cell<moris::moris_index*>      tNodeInds(tGeoObjects.size());
                moris::Matrix< moris::DDRMat > tParametricCoordsRelativeToParentElem(tGeoObjects.size(),tDimParamCoord);

                // get reference to child mesh edge parent information
                moris::Matrix< moris::IndexMat > const & tEdgeParentIndices = tChildMesh.get_edge_parent_inds();
                moris::Matrix< moris::DDSTMat >  const & tEdgeParentRanks   = tChildMesh.get_edge_parent_ranks();

                for (moris::size_t k = 0; k < tGeoObjects.size(); k++)
                {
                    if(!tGeoObjects(k).has_parent_nodes_on_interface())
                    {
                        // Local index to XTK Mesh
                        tEdgeInd = tGeoObjects(k).get_parent_entity_index();

                        // get a local coordinate along the intersected edge [-1,1]
                        tLocalCoordRelativeToEdge(0,0) = tGeoObjects(k).get_interface_lcl_coord();

                        // get the interpolated global coordinate
                        tGlobalCoord = tGeoObjects(k).get_interface_glb_coord();

                        // Add edge to the entity intersection connectivity
                        mCutMesh.add_entity_to_intersect_connectivity(aActiveChildMeshIndices(0,j), tNumNewNodes, tEdgeInd, 0);

                        // Edge nodes
                        moris::Matrix<moris::IndexMat> tEdgeNodes = tEdgeToNode.get_row(tEdgeInd);

                        // Compute new node parametric coordinate with respect to the current parent element
                        tEdgeNodeParamCoordinates.set_row(0, tChildMesh.get_parametric_coordinates(tEdgeNodes(0)));
                        tEdgeNodeParamCoordinates.set_row(1, tChildMesh.get_parametric_coordinates(tEdgeNodes(1)));
                        moris::Matrix< moris::DDRMat > tParametricCoordsRelativeToParentElem = Interpolation::linear_interpolation_location(tEdgeNodeParamCoordinates,tLocalCoordRelativeToEdge);

                        // Parent edge information
                        moris::size_t      tParentRank  = tEdgeParentRanks(0, tEdgeInd);
                        moris::moris_index tParentIndex = tEdgeParentIndices(0, tEdgeInd);

                        // get the owning processor for an entity
                        moris::moris_index tOwningProc = tMeshData.get_entity_owner(tParentIndex, (enum EntityRank)tParentRank);

                        // Convert to global id using mesh
                        tEdgeNodes(0, 0) = mBackgroundMesh.get_glb_entity_id_from_entity_loc_index(tEdgeNodes(0, 0), EntityRank::NODE);
                        tEdgeNodes(0, 1) = mBackgroundMesh.get_glb_entity_id_from_entity_loc_index(tEdgeNodes(0, 1), EntityRank::NODE);

                        // Order the nodes in ascending order
                        if(tEdgeNodes(0, 1) < tEdgeNodes(0, 0))
                        {
                            moris::size_t tSwap = tEdgeNodes(0, 0);
                            tEdgeNodes(0, 0) = tEdgeNodes(0, 1);
                            tEdgeNodes(0, 1) = tSwap;
                        }

                        // Intersected edge is an existing  edge
                        // Make request in edge requests
                        // This does not require a supplemental identifier
                        // TODO: ADD OVERFLOW CHECK IN CANTOR PAIRING!!!!!!
                        moris::moris_index tSecondaryId = xtk::cantor_pairing(tEdgeNodes(0, 0),tEdgeNodes(0, 1));
                        moris_index tNewNodeIndexInSubdivision = MORIS_INDEX_MAX;
                        bool tRequestExist = tDecompData.request_exists(tParentIndex,tSecondaryId,(enum EntityRank)tParentRank,tNewNodeIndexInSubdivision);

                        // location for this face in the map
                        if(!tRequestExist)
                        {
                            tNewNodeIndexInSubdivision = tDecompData.register_new_request(tParentIndex,
                                                                                          tSecondaryId,
                                                                                          tOwningProc,
                                                                                          (enum EntityRank)tParentRank,
                                                                                          tGlobalCoord,
                                                                                          new Edge_Topology(tEdgeToNode.get_row(tEdgeInd)), /*Note: this is deleted in the decomp data deconstructor*/
                                                                                          tLocalCoordRelativeToEdge.get_row(0));
                        }

                        // add to pending node pointers for child mesh
                        tDecompData.tCMNewNodeLoc(j).push_back(tNewNodeIndexInSubdivision);

                        // add parametric coordinate to decomp data
                        tDecompData.tCMNewNodeParamCoord(j).push_back(tParametricCoordsRelativeToParentElem);

                        // Creating a new node add 1 to count
                        tNumNewNodes++;
                    }

                    else if(tGeoObjects(k).all_parent_nodes_on_interface())
                    {
                        moris::moris_index tParentIndex = tGeoObjects(k).get_parent_entity_index();

                        // Tell the child mesh this edge is actually on the interface already
                        tChildMesh.mark_edge_as_on_interface(tParentIndex);

                        // Tell the xtk mesh that these edge nodes are interface nodes
                        mBackgroundMesh.mark_node_as_interface_node(tEdgeToNode(tParentIndex,0),mGeometryEngine->get_active_geometry_index());
                        mBackgroundMesh.mark_node_as_interface_node(tEdgeToNode(tParentIndex,1),mGeometryEngine->get_active_geometry_index());
                    }
                } // geometry object

                tChildMesh.mark_interface_faces_from_interface_coincident_faces();
            } // XTK Mesh loop

            moris_index tMessageTag = 60001; /*arbitrary tag for regular subdivision*/
            assign_node_requests_identifiers(tDecompData,tMessageTag);

            // Allocate interface flag space in XTK mesh even though these are not interface nodes
            mBackgroundMesh.allocate_space_in_interface_node_flags(tDecompData.tNewNodeIndex.size(),mGeometryEngine->get_num_geometries());

            // add nodes to the background mesh
            mBackgroundMesh.batch_create_new_nodes(tDecompData.tNewNodeId,tDecompData.tNewNodeIndex,tDecompData.tNewNodeOwner,tDecompData.tNewNodeCoordinate);

            // add nodes to child mesh
            this->decompose_internal_set_new_nodes_in_child_mesh_nh(aActiveChildMeshIndices,tDecompData);

            // associate new nodes with geometry objects
            create_new_node_association_with_geometry(tDecompData);

            // mark nodes as interface nodes
            moris_index tGeomIndex = mGeometryEngine->get_active_geometry_index();
            for(moris::uint i = 0; i <tDecompData.tNewNodeId.size(); i++)
            {
                mBackgroundMesh.mark_node_as_interface_node(tDecompData.tNewNodeIndex(i),tGeomIndex);

                // determine if this vertex is on other interfaces
                for(moris::uint j = 0; j < mGeometryEngine->get_num_geometries(); j++)
                {
                    moris::real const & tPhaseVal = mGeometryEngine->get_entity_phase_val(tDecompData.tNewNodeIndex(i),(moris_index)j);
                    if(moris::equal_to(0.0,tPhaseVal))
                    {
                        mBackgroundMesh.mark_node_as_interface_node(tDecompData.tNewNodeIndex(i),j);
                    }
                }
            }


            // Set Node Ids and tell the child mesh to update
            for (moris::size_t j = 0; j < aActiveChildMeshIndices.n_cols(); j++)
            {
                moris::Matrix< moris::IndexMat > const & tNodeIndices = mCutMesh.get_node_indices(aActiveChildMeshIndices(0,j));
                moris::Matrix< moris::IdMat > tNodeIds = mBackgroundMesh.get_glb_entity_id_from_entity_loc_index_range(tNodeIndices, EntityRank::NODE);

                mCutMesh.set_node_ids(aActiveChildMeshIndices(0,j), tNodeIds);
                mCutMesh.modify_templated_mesh(aActiveChildMeshIndices(0,j), TemplateType::HIERARCHY_TET4);
            }

            break;
                }
        case (Subdivision_Method::C_TRI3):
    {
            // If it the first subdivision we need to find the intersected before placing the conformal nodes
            // Intersected elements are flagged via the Geometry_Engine
            if(aFirstSubdivision)
            {

                moris::Matrix< moris::IndexMat > tNewPairBool;
                run_first_cut_routine(TemplateType::TRI_3, aGeomIndex, 3, aActiveChildMeshIndices,tNewPairBool);

                for(moris::size_t i = 0; i<aActiveChildMeshIndices.n_cols(); i++)
                {
                    Child_Mesh & tChildMesh = mCutMesh.get_child_mesh(aActiveChildMeshIndices(0,i));
                    tChildMesh.generate_connectivities(false,true,true);
                }
                // set the child cell topology as tet 4s
                mCutMesh.set_child_element_topology(CellTopology::TRI3);

            }

            // For quad background meshes we have a 2 dimension parametric coordinate
            moris::size_t tDimParamCoord = 2;

            // For tri background meshes we have a 3-d parametric coordinate
            if(aFirstSubdivision)
            {
                tDimParamCoord = 3;
            }

            // initialize a struct of all the data we are keeping track of in this decomposition
            // intended to reduce the clutter of function inputs etc
            Decomposition_Data tDecompData;
            tDecompData.mSubdivisionMethod = Subdivision_Method::C_TRI3;
            tDecompData.mConformalDecomp = true;
            tDecompData.mHasSecondaryIdentifier = true;
            tDecompData.mFirstSubdivision = aFirstSubdivision;

            // Initialize
            moris::size_t tEdgeInd = MORIS_INDEX_MAX;

            // Initialize topologies used in this method (all local coordinates are with respect to an edge)
            Edge_Topology tEdgeTopology;

            // initialize a couple of commonly used matrices in this method
            moris::Matrix< moris::DDRMat > tLocalCoordRelativeToEdge(1,1, 0); // ALong an edge
            moris::Matrix< moris::DDRMat > tGlobalCoord(1,2, 0); // ALong an edge
            moris::Matrix< moris::DDRMat > tEdgeNodeParamCoordinates(2,tDimParamCoord); // parametric coordinate of end nodes wrt parent element

            // Check type specified as conformal (could change this to enum)
            moris::size_t tCheckType = 1;
            moris::Matrix< moris::DDRMat > tNodeCoords = mBackgroundMesh.get_all_node_coordinates_loc_inds();

            // get the underlying background mesh data
            moris::mtk::Mesh & tMeshData = mBackgroundMesh.get_mesh_data();

            // resize child mesh to new node information
            tDecompData.tCMNewNodeLoc.resize(aActiveChildMeshIndices.n_cols());
            tDecompData.tCMNewNodeParamCoord.resize(aActiveChildMeshIndices.n_cols());

            for (moris::size_t j = 0; j < aActiveChildMeshIndices.n_cols(); j++)
            {

                // Initialize geometry objects
                Cell<moris::ge::GEN_Geometry_Object> tGeoObjects;

                // Get the child mesh that is active
                Child_Mesh & tChildMesh = mCutMesh.get_child_mesh(aActiveChildMeshIndices(0,j));

                // edge to node connectivity from child mesh
                moris::Matrix< moris::IndexMat > const & tEdgeToNode = tChildMesh.get_edge_to_node();

                // Ask geometry engine which edges are intersected (Simple mesh local indexed edges)
                mGeometryEngine->is_intersected(tNodeCoords, tEdgeToNode, tCheckType, tGeoObjects);

                // Initialize node index pointers based on number of intersected edges and parametric coordinates
                uint tNumNewNodes = 0;
                Cell<moris::moris_index*>      tNodeInds(tGeoObjects.size());
                moris::Matrix< moris::DDRMat > tParametricCoordsRelativeToParentElem(tGeoObjects.size(),tDimParamCoord);

                // get reference to child mesh edge parent information
                moris::Matrix< moris::IndexMat > const & tEdgeParentIndices = tChildMesh.get_edge_parent_inds();
                moris::Matrix< moris::DDSTMat >  const & tEdgeParentRanks   = tChildMesh.get_edge_parent_ranks();

                for (moris::size_t k = 0; k < tGeoObjects.size(); k++)
                {
                    if(!tGeoObjects(k).has_parent_nodes_on_interface())
                    {
                        // Local index to XTK Mesh
                        tEdgeInd = tGeoObjects(k).get_parent_entity_index();

                        // get a local coordinate along the intersected edge [-1,1]
                        tLocalCoordRelativeToEdge(0,0) = tGeoObjects(k).get_interface_lcl_coord();

                        // get the interpolated global coordinate
                        tGlobalCoord = tGeoObjects(k).get_interface_glb_coord();

                        // Add edge to the entity intersection connectivity
                        mCutMesh.add_entity_to_intersect_connectivity(aActiveChildMeshIndices(0,j), tNumNewNodes, tEdgeInd, 0);

                        // Edge nodes
                        moris::Matrix<moris::IndexMat> tEdgeNodes = tEdgeToNode.get_row(tEdgeInd);

                        // Compute new node parametric coordinate with respect to the current parent element
                        tEdgeNodeParamCoordinates.set_row(0, tChildMesh.get_parametric_coordinates(tEdgeNodes(0)));
                        tEdgeNodeParamCoordinates.set_row(1, tChildMesh.get_parametric_coordinates(tEdgeNodes(1)));
                        moris::Matrix< moris::DDRMat > tParametricCoordsRelativeToParentElem = Interpolation::linear_interpolation_location(tEdgeNodeParamCoordinates,tLocalCoordRelativeToEdge);

                        // Parent edge information
                        moris::size_t      tParentRank  = tEdgeParentRanks(0, tEdgeInd);
                        moris::moris_index tParentIndex = tEdgeParentIndices(0, tEdgeInd);

                        // get the owning processor for an entity
                        moris::moris_index tOwningProc = tMeshData.get_entity_owner(tParentIndex, (enum EntityRank)tParentRank);


                        // Convert to global id using mesh
                        tEdgeNodes(0, 0) = mBackgroundMesh.get_glb_entity_id_from_entity_loc_index(tEdgeNodes(0, 0), EntityRank::NODE);
                        tEdgeNodes(0, 1) = mBackgroundMesh.get_glb_entity_id_from_entity_loc_index(tEdgeNodes(0, 1), EntityRank::NODE);

                        // Order the nodes in ascending order
                        if(tEdgeNodes(0, 1) < tEdgeNodes(0, 0))
                        {
                            moris::size_t tSwap = tEdgeNodes(0, 0);
                            tEdgeNodes(0, 0) = tEdgeNodes(0, 1);
                            tEdgeNodes(0, 1) = tSwap;
                        }

                        // Intersected edge is an existing  edge
                        // Make request in edge requests
                        // This does not require a supplemental identifier
                        // TODO: ADD OVERFLOW CHECK IN CANTOR PAIRING!!!!!!
                        moris::moris_index tSecondaryId = xtk::cantor_pairing(tEdgeNodes(0, 0),tEdgeNodes(0, 1));
                        moris_index tNewNodeIndexInSubdivision = MORIS_INDEX_MAX;
                        bool tRequestExist = tDecompData.request_exists(tParentIndex,tSecondaryId,(enum EntityRank)tParentRank,tNewNodeIndexInSubdivision);

                        // location for this face in the map
                        if(!tRequestExist)
                        {
                            tNewNodeIndexInSubdivision = tDecompData.register_new_request(tParentIndex,
                                                                                          tSecondaryId,
                                                                                          tOwningProc,
                                                                                          (enum EntityRank)tParentRank,
                                                                                          tGlobalCoord,
                                                                                          new Edge_Topology(tEdgeToNode.get_row(tEdgeInd)), /*Note: this is deleted in the decomp data deconstructor*/
                                                                                          tLocalCoordRelativeToEdge.get_row(0));
                        }

                        // add to pending node pointers for child mesh
                        tDecompData.tCMNewNodeLoc(j).push_back(tNewNodeIndexInSubdivision);

                        // add parametric coordinate to decomp data
                        tDecompData.tCMNewNodeParamCoord(j).push_back(tParametricCoordsRelativeToParentElem);

                        // Creating a new node add 1 to count
                        tNumNewNodes++;
                    }

                    else if(tGeoObjects(k).all_parent_nodes_on_interface())
                    {
                        MORIS_ERROR(0,"All parents on interface not handled in 2D yet");
                    }
                } // geometry object

            } // XTK Mesh loop

            moris_index tMessageTag = 60001; /*arbitrary tag for regular subdivision*/
            assign_node_requests_identifiers(tDecompData,tMessageTag);

            // Allocate interface flag space in XTK mesh even though these are not interface nodes
            mBackgroundMesh.allocate_space_in_interface_node_flags(tDecompData.tNewNodeIndex.size(),mGeometryEngine->get_num_geometries());

            // add nodes to the background mesh
            mBackgroundMesh.batch_create_new_nodes(tDecompData.tNewNodeId,tDecompData.tNewNodeIndex,tDecompData.tNewNodeOwner,tDecompData.tNewNodeCoordinate);

            // add nodes to child mesh
            this->decompose_internal_set_new_nodes_in_child_mesh_nh(aActiveChildMeshIndices,tDecompData);

            // associate new nodes with geometry objects
            create_new_node_association_with_geometry(tDecompData);

            // mark nodes as interface nodes
            moris_index tGeomIndex = mGeometryEngine->get_active_geometry_index();
            for(moris::uint i = 0; i <tDecompData.tNewNodeId.size(); i++)
            {
                // this node is always on the geometry interface of current so mark this
                mBackgroundMesh.mark_node_as_interface_node(tDecompData.tNewNodeIndex(i),tGeomIndex);

                // determine if this vertex is on other interfaces
                for(moris::uint j = 0; j < mGeometryEngine->get_num_geometries(); j++)
                {
                    moris::real const & tPhaseVal = mGeometryEngine->get_entity_phase_val(tDecompData.tNewNodeIndex(i),(moris_index)j);
                    if(moris::equal_to(0.0,tPhaseVal))
                    {
                        mBackgroundMesh.mark_node_as_interface_node(tDecompData.tNewNodeIndex(i),j);
                    }
                }
            }

            // Set Node Ids and tell the child mesh to update
            for (moris::size_t j = 0; j < aActiveChildMeshIndices.n_cols(); j++)
            {
                moris::Matrix< moris::IndexMat > const & tNodeIndices = mCutMesh.get_node_indices(aActiveChildMeshIndices(0,j));
                moris::Matrix< moris::IdMat > tNodeIds = mBackgroundMesh.get_glb_entity_id_from_entity_loc_index_range(tNodeIndices, EntityRank::NODE);

                mCutMesh.set_node_ids(aActiveChildMeshIndices(0,j), tNodeIds);
                mCutMesh.modify_templated_mesh(aActiveChildMeshIndices(0,j), TemplateType::CONFORMAL_TRI3);
            }



            break;
    }
        default:
        {
            moris::size_t breaker = 0;
            MORIS_ERROR(breaker != 0, "formulate_node_request should not enter the default case, check to see if your aCheckType is undefined.");
        }
    }
}


void
Model::decompose_internal_reg_sub_hex8_make_requests(moris::Matrix< moris::IndexMat > & aActiveChildMeshIndices,
                                                     moris::Matrix< moris::IndexMat > & tNewPairBool,
                                                     Decomposition_Data & tDecompData)
{
    // mesh data accessor
    moris::mtk::Mesh & tXTKMeshData = mBackgroundMesh.get_mesh_data();

    // number of intersected elements
    moris::uint tIntersectedCount = aActiveChildMeshIndices.n_cols();

    // allocate child mesh to new node location
    tDecompData.tCMNewNodeLoc.resize(tIntersectedCount,7);
    tDecompData.tCMNewNodeParamCoord.resize(tIntersectedCount);

    // parametric coordinates relative to hex where we put the nodes
    // Parametric coordinates for this subdivision routine
    const moris::Matrix< moris::DDRMat > tParamCoordsRelativeToElem(
            {{ 0.0, -1.0,  0.0},
        { 1.0,  0.0,  0.0},
        { 0.0,  1.0,  0.0},
        {-1.0,  0.0,  0.0},
        { 0.0,  0.0, -1.0},
        { 0.0,  0.0,  1.0},
        { 0.0,  0.0,  0.0}});

    const moris::Matrix< moris::DDRMat > tParamCoordsRelativeToFace(
            {{ 0.0,  0.0},
        { 0.0,  0.0},
        { 0.0,  0.0},
        { 0.0,  0.0},
        { 0.0,  0.0},
        { 0.0,  0.0}});


    // get the underlying background mesh data
    moris::mtk::Mesh & tMeshData = mBackgroundMesh.get_mesh_data();

    // setup Child mesh to new node location
    for (moris::size_t i = 0; i < tIntersectedCount; i++)
    {
        if(tNewPairBool(0,i) == 0)
        {

            // Get element index
            moris::moris_index tElemInd = mCutMesh.get_parent_element_index(aActiveChildMeshIndices(0,i));

            // Get local index of faces connected to element using local element index
            moris::Matrix<moris::IndexMat> tFaceIndices = tXTKMeshData.get_entity_connected_to_entity_loc_inds(tElemInd, moris::EntityRank::ELEMENT, moris::EntityRank::FACE);

            // Loop over faces (6 in a hex 8) and set a node request.
            // Request will return a pointer to where the created node index will be placed
            for (moris::size_t fi = 0; fi < 6; fi++)
            {

                moris_index tRequestLoc = MORIS_INDEX_MAX;
                bool tRequestExists = tDecompData.request_exists(tFaceIndices(fi),EntityRank::FACE,tRequestLoc);

                // if we haven't created a node on this face then create one
                if(!tRequestExists)
                {
                    // node indices attached to face fi
                    moris::Matrix<moris::IndexMat> tFaceNodes = tXTKMeshData.get_entity_connected_to_entity_loc_inds(tFaceIndices(fi), moris::EntityRank::FACE, moris::EntityRank::NODE);

                    // face owner
                    moris::moris_index tOwningProc = tMeshData.get_entity_owner(tFaceIndices(fi), EntityRank::FACE);


                    // coordinates of nodes attached to the nodes of this face
                    moris::Matrix<moris::DDRMat> tCoordinates = mBackgroundMesh.get_selected_node_coordinates_loc_inds(tFaceNodes);

                    // bilinearly interpolate to the center of this face fi
                    moris::Matrix<moris::DDRMat> tNewNodeCoordinates;
                    xtk::Interpolation::bilinear_interpolation(tCoordinates, tParamCoordsRelativeToFace.get_row(fi),tNewNodeCoordinates);

                    // location for this face in the map
                    moris_index tNewNodeIndexInSubdivision = tDecompData.register_new_request(tFaceIndices(fi),
                                                                                              tOwningProc,
                                                                                              EntityRank::FACE,
                                                                                              tNewNodeCoordinates,
                                                                                              new Quad_4_Topology(tFaceNodes), /*Note: this is deleted in the decomp data deconstructor*/
                                                                                              tParamCoordsRelativeToFace.get_row(fi));
                    // add to pending node pointers for child mesh
                    tDecompData.tCMNewNodeLoc(i)(fi) = tNewNodeIndexInSubdivision;

                    // add parametric coordinate to decomp data
                    tDecompData.tCMNewNodeParamCoord(i).push_back(tParamCoordsRelativeToElem.get_row(fi));
                }

                // if debug check the coordinate will be the same
                else
                {
                    tDecompData.tCMNewNodeLoc(i)(fi) = tRequestLoc;
                    tDecompData.tCMNewNodeParamCoord(i).push_back(tParamCoordsRelativeToElem.get_row(fi));
#ifdef DEBUG
                    moris::uint tNewNodeIndexInSubdivision = tRequestLoc;

                    // node indices attached to face fi
                    moris::Matrix<moris::IndexMat> tFaceNodes = tXTKMeshData.get_entity_connected_to_entity_loc_inds(tFaceIndices(fi), moris::EntityRank::FACE, moris::EntityRank::NODE);

                    // coordinates of nodes attached to the nodes of this face
                    moris::Matrix<moris::DDRMat> tCoordinates = mBackgroundMesh.get_selected_node_coordinates_loc_inds(tFaceNodes);

                    // bilinearly interpolate to the center of this face fi
                    moris::Matrix<moris::DDRMat> tNewNodeCoordinates;
                    xtk::Interpolation::bilinear_interpolation(tCoordinates, tParamCoordsRelativeToFace.get_row(fi),tNewNodeCoordinates);

                    // other coordinate
                    moris::Matrix<moris::DDRMat> tExistingNodeCoordinate = tDecompData.tNewNodeCoordinate(tNewNodeIndexInSubdivision);

                    MORIS_ASSERT(all_true(tNewNodeCoordinates == tExistingNodeCoordinate) ,"Node coordinates created on same face do not match");
#endif

                }
            }

            // Place node at center of element
            // get the nodes attached to the element
            moris::Matrix<moris::IndexMat>tElementNodes = tXTKMeshData.get_entity_connected_to_entity_loc_inds(tElemInd, moris::EntityRank::ELEMENT, moris::EntityRank::NODE);

            // coordinates of nodes attached to element
            moris::Matrix<moris::DDRMat> tCoordinates = mBackgroundMesh.get_selected_node_coordinates_loc_inds(tElementNodes);

            // trilinearly interpolate to the center of the element
            moris::Matrix<moris::DDRMat> tNewNodeCoordinates;
            xtk::Interpolation::trilinear_interpolation(tCoordinates, tParamCoordsRelativeToElem.get_row(6), tNewNodeCoordinates);

            // add the new node at center of element to the map
            // location for this face in the map
            moris_index tNewNodeIndexInSubdivision = MORIS_INDEX_MAX;

            // owner of element
            moris::moris_index tOwningProc = tMeshData.get_entity_owner(tElemInd, EntityRank::ELEMENT);


            MORIS_ASSERT(!tDecompData.request_exists(tElemInd,EntityRank::ELEMENT,tNewNodeIndexInSubdivision),"All element requests should be unique, therefore tNewRequest is expected to be true here");


            tNewNodeIndexInSubdivision = tDecompData.register_new_request(tElemInd,
                                                                          tOwningProc,
                                                                          EntityRank::ELEMENT,
                                                                          tNewNodeCoordinates,
                                                                          new Hexahedron_8_Topology(tElementNodes), /*Note: this is deleted in the decomp data deconstructor*/
                                                                          tParamCoordsRelativeToElem.get_row(6));

            // add child mesh new node location and parametric coordinate relative to element
            tDecompData.tCMNewNodeLoc(i)(6) = tNewNodeIndexInSubdivision;
            tDecompData.tCMNewNodeParamCoord(i).push_back(tParamCoordsRelativeToElem.get_row(6));

        }

    }

}

void
Model::decompose_internal_reg_sub_quad4_make_requests(moris::Matrix< moris::IndexMat > & aActiveChildMeshIndices,
                                                      moris::Matrix< moris::IndexMat > & tNewPairBool,
                                                      Decomposition_Data                & tDecompData)
{
    // Get access to mesh data
    moris::mtk::Mesh & tXTKMeshData = mBackgroundMesh.get_mesh_data();

    // Get number of intersected elements
    moris::uint tIntersectedCount = aActiveChildMeshIndices.n_cols();

    // Allocate child mesh to new node location
    tDecompData.tCMNewNodeLoc.resize(tIntersectedCount,1);
    tDecompData.tCMNewNodeParamCoord.resize(tIntersectedCount);

    // Parametric coordinates relative to quad where we put the nodes
    const moris::Matrix< moris::DDRMat > tParamCoordsRelativeToElem(
            {{ 0.0,  0.0}});

    // get the underlying background mesh data
    moris::mtk::Mesh & tMeshData = mBackgroundMesh.get_mesh_data();

    // setup Child mesh to new node location
    for (moris::size_t i = 0; i < tIntersectedCount; i++)
    {
        if(tNewPairBool(0,i) == 0)
        {

            // Get element index
            moris::moris_index tElemInd = mCutMesh.get_parent_element_index(aActiveChildMeshIndices(0,i));

            // Place node at center of element
            // get the nodes attached to the element
            moris::Matrix<moris::IndexMat>tElementNodes = tXTKMeshData.get_entity_connected_to_entity_loc_inds(tElemInd, moris::EntityRank::ELEMENT, moris::EntityRank::NODE);

            // coordinates of nodes attached to element
            moris::Matrix<moris::DDRMat> tCoordinates = mBackgroundMesh.get_selected_node_coordinates_loc_inds(tElementNodes);

            // trilinearly interpolate to the center of the element
            moris::Matrix<moris::DDRMat> tNewNodeCoordinates;
            xtk::Interpolation::bilinear_interpolation(tCoordinates, tParamCoordsRelativeToElem.get_row(0), tNewNodeCoordinates);

            // add the new node at center of element to the map
            // location for this face in the map
            moris_index tNewNodeIndexInSubdivision = MORIS_INDEX_MAX;

            // owner of element
            moris::moris_index tOwningProc = tMeshData.get_entity_owner(tElemInd, EntityRank::ELEMENT);


            MORIS_ASSERT(!tDecompData.request_exists(tElemInd,EntityRank::ELEMENT,tNewNodeIndexInSubdivision),"All element requests should be unique, therefore tNewRequest is expected to be true here");


            tNewNodeIndexInSubdivision = tDecompData.register_new_request(tElemInd,
                                                                          tOwningProc,
                                                                          EntityRank::ELEMENT,
                                                                          tNewNodeCoordinates,
                                                                          new Quad_4_Topology(tElementNodes), /*Note: this is deleted in the decomp data deconstructor*/
                                                                          tParamCoordsRelativeToElem.get_row(0));

            // add child mesh new node location and parametric coordinate relative to element
            tDecompData.tCMNewNodeLoc(i)(0) = tNewNodeIndexInSubdivision;
            tDecompData.tCMNewNodeParamCoord(i).push_back(tParamCoordsRelativeToElem.get_row(0));

        }

    }
}


void
Model::decompose_internal_set_new_nodes_in_child_mesh_reg_sub(moris::Matrix< moris::IndexMat > & aActiveChildMeshIndices,
                                                              moris::Matrix< moris::IndexMat > & tNewPairBool,
                                                              moris::real                        tNumParamCoords,
                                                              Decomposition_Data &               tDecompData)
{
    // number of intersected elements
    moris::uint tIntersectedCount = aActiveChildMeshIndices.n_cols();

    // iterate through active child mesh indices
    for(moris::uint i = 0 ; i <tIntersectedCount; i++)
    {
        // only regularly subdivide if it hasnt already been regularly subdivided
        if(tNewPairBool(0,i) == 0)
        {

            // number of new nodes for child mesh i
            moris::uint tNumNewNodesForCM = tDecompData.tCMNewNodeLoc(i).size();

            // matrix of new node indices
            moris::Matrix<IndexMat> tCMNewNodeInds(1,tNumNewNodesForCM);

            // matrix of new node ids
            moris::Matrix<IdMat> tCMNewNodeIds(1,tNumNewNodesForCM);

            // iterate through new nodes for child mesh i to collect index and id
            for(moris::uint iN =0; iN< tNumNewNodesForCM; iN++)
            {
                // location relative to the decomposition data
                moris::moris_index tNodeIndexInRequestVect = tDecompData.tCMNewNodeLoc(i)(iN);

                // retreive node index and id
                tCMNewNodeInds(iN) = tDecompData.tNewNodeIndex(tNodeIndexInRequestVect);
                tCMNewNodeIds(iN)  = tDecompData.tNewNodeId(tNodeIndexInRequestVect);
            }

            // retrieve child mesh
            Child_Mesh & tChildMesh = mCutMesh.get_child_mesh(aActiveChildMeshIndices(i));

            // add node indices, ids, and vertices to child mesh
            tChildMesh.add_node_indices(tCMNewNodeInds);
            tChildMesh.add_node_ids(tCMNewNodeIds);

            // allocate space for parametric coordinates
            tChildMesh.allocate_parametric_coordinates(tNumNewNodesForCM,tNumParamCoords);

            // iterate through nods and add parametric coordinate
            for(moris::uint iN =0; iN< tNumNewNodesForCM; iN++)
            {
                tChildMesh.add_node_parametric_coordinate( tCMNewNodeInds(iN),tDecompData.tCMNewNodeParamCoord(i)(iN));
            }
        }
    }

}



void
Model::decompose_internal_set_new_nodes_in_child_mesh_nh(moris::Matrix< moris::IndexMat > & aActiveChildMeshIndices,
                                                         Decomposition_Data &               tDecompData)
{
    // number of intersected elements
    moris::uint tIntersectedCount = aActiveChildMeshIndices.n_cols();

    // iterate through active child mesh indices
    for(moris::uint i = 0 ; i <tIntersectedCount; i++)
    {
        // number of new nodes for child mesh i
        moris::uint tNumNewNodesForCM = tDecompData.tCMNewNodeLoc(i).size();

        // matrix of new node indices
        moris::Matrix<IndexMat> tCMNewNodeInds(1,tNumNewNodesForCM);

        // matrix of new node ids
        moris::Matrix<IdMat> tCMNewNodeIds(1,tNumNewNodesForCM);

        // iterate through new nodes for child mesh i to collect index and id
        for(moris::uint iN =0; iN< tNumNewNodesForCM; iN++)
        {
            // location relative to the decomposition data
            moris::moris_index tNodeIndexInRequestVect = tDecompData.tCMNewNodeLoc(i)(iN);

            // retreive node index and id
            tCMNewNodeInds(iN) = tDecompData.tNewNodeIndex(tNodeIndexInRequestVect);
            tCMNewNodeIds(iN)  = tDecompData.tNewNodeId(tNodeIndexInRequestVect);
        }

        // retrieve child mesh
        Child_Mesh & tChildMesh = mCutMesh.get_child_mesh(aActiveChildMeshIndices(i));


        // add node indices, ids, and vertices to child mesh
        tChildMesh.add_node_indices(tCMNewNodeInds);
        tChildMesh.add_node_ids(tCMNewNodeIds);

        // allocate space for parametric coordinate

        // For hex background meshes we have a three dimension parametric coordinate
        enum CellTopology tBackgroundTopo = mBackgroundMesh.get_parent_cell_topology();

        moris::size_t tDimParamCoord = 0;

        if(tBackgroundTopo == CellTopology::HEX8)
        {
            tDimParamCoord =3;
        }
        else if(tBackgroundTopo == CellTopology::TET4)
        {
            tDimParamCoord =4;
        }
        else if(tBackgroundTopo == CellTopology::QUAD4)
        {
            tDimParamCoord = 2;
        }
        else if(tBackgroundTopo == CellTopology::TRI3)
        {
            tDimParamCoord = 3;
        }
        else
        {
            MORIS_ERROR(0,"Invalid background cell topo");
        }


        tChildMesh.allocate_parametric_coordinates(tNumNewNodesForCM,tDimParamCoord);

        // iterate through nods and add parametric coordinate
        for(moris::uint iN =0; iN< tNumNewNodesForCM; iN++)
        {
            tChildMesh.add_node_parametric_coordinate( tCMNewNodeInds(iN),tDecompData.tCMNewNodeParamCoord(i)(iN));
        }
    }
}

void
Model::create_new_node_association_with_geometry(Decomposition_Data & tDecompData)
{
    // create geometry objects for each node
    mGeometryEngine->create_new_node_geometry_objects(tDecompData.tNewNodeIndex,
                                                     tDecompData.mConformalDecomp,
                                                     tDecompData.tNewNodeParentTopology,
                                                     tDecompData.tParamCoordRelativeToParent,
                                                     tDecompData.tNewNodeCoordinate);
}


void
Model::assign_node_requests_identifiers(Decomposition_Data & aDecompData,
                                        moris::moris_index   aMPITag)
{
    barrier();
    // asserts
    MORIS_ASSERT(aDecompData.tNewNodeId.size() == aDecompData.tNewNodeIndex.size(),      "Dimension mismatch in assign_node_requests_identifiers");
    MORIS_ASSERT(aDecompData.tNewNodeId.size() == aDecompData.tNewNodeParentRank.size(), "Dimension mismatch in assign_node_requests_identifiers");
    MORIS_ASSERT(aDecompData.tNewNodeId.size() == aDecompData.tNewNodeParentIndex.size(),"Dimension mismatch in assign_node_requests_identifiers");

    // owned requests and shared requests sorted by owning proc
    Cell<uint> tOwnedRequest;
    Cell<Cell<uint>> tNotOwnedRequests;
    Cell<uint> tProcRanks;
    std::unordered_map<moris_id,moris_id> tProcRankToDataIndex;
    this->sort_new_node_requests_by_owned_and_not_owned(aDecompData,tOwnedRequest,tNotOwnedRequests,tProcRanks,tProcRankToDataIndex);

    // allocate ids for nodes I own
    moris::moris_id tNodeId  = mBackgroundMesh.allocate_entity_ids(aDecompData.tNewNodeId.size(), EntityRank::NODE);

    // get first available index
    moris::moris_id tNodeInd = mBackgroundMesh.get_first_available_index(EntityRank::NODE);

    // Assign owned request identifiers
    this->assign_owned_request_identifiers(aDecompData, tOwnedRequest, tNodeInd, tNodeId);

    // prepare node information request data
    Cell<Matrix<IndexMat>> tOutwardRequests;
    this->setup_outward_requests(aDecompData, tNotOwnedRequests, tProcRanks, tProcRankToDataIndex, tOutwardRequests);


    // send requests to owning processor
    this->send_outward_requests(aMPITag,tProcRanks,tOutwardRequests);

    // hold on to make sure everyone has sent all their information
    barrier();



    // receive the requests
    Cell<Matrix<IndexMat>> tReceivedRequests;
    Cell<uint> tProcsReceivedFrom;
    this->inward_receive_requests(aMPITag, 3, tReceivedRequests, tProcsReceivedFrom);

    // Prepare request answers
    Cell<Matrix<IndexMat>> tRequestAnwers;
    this->prepare_request_answers(aDecompData,tReceivedRequests,tRequestAnwers);

    // send the answers back
    this->return_request_answers(aMPITag+1, tRequestAnwers, tProcsReceivedFrom);

    barrier();

    // receive the answers
    Cell<Matrix<IndexMat>> tReceivedRequestsAnswers;
    this->inward_receive_request_answers(aMPITag+1,1,tProcRanks,tReceivedRequestsAnswers);

    // handle received information
    this->handle_received_request_answers(aDecompData,tOutwardRequests,tReceivedRequestsAnswers,tNodeInd,tNodeId);

    // return index to update
    mBackgroundMesh.update_first_available_index(tNodeInd,EntityRank::NODE);

    MORIS_ASSERT(this->verify_successful_node_assignment(aDecompData),"Unsuccesssful node assignment detected.");

    barrier();
}


void
Model::sort_new_node_requests_by_owned_and_not_owned(Decomposition_Data                    & tDecompData,
                                                     Cell<uint>                            & aOwnedRequests,
                                                     Cell<Cell<uint>>                      & aNotOwnedRequests,
                                                     Cell<uint>                            & aProcRanks,
                                                     std::unordered_map<moris_id,moris_id> & aProcRankToIndexInData)
{
    // access the mesh data behind the background mesh
    moris::mtk::Mesh & tMeshData = mBackgroundMesh.get_mesh_data();

    // access the communication
    Matrix<IdMat> tCommTable = mBackgroundMesh.get_communication_table();


    // number of new nodes
    moris::uint tNumNewNodes = tDecompData.tNewNodeParentIndex.size();

    // Par rank
    moris::moris_index tParRank = par_rank();

    // index
    moris::uint tCurrentIndex = 0;

    // resize proc ranks and setup map to comm table
    aProcRanks.resize(tCommTable.numel());
    for(moris::uint i = 0; i <tCommTable.numel(); i++)
    {
        aProcRankToIndexInData[tCommTable(i)] = i;
        aProcRanks(i) = (tCommTable(i));
        aNotOwnedRequests.push_back(Cell<uint>(0));
    }

    // iterate through each node request and figure out the owner
    for(moris::uint i = 0; i <tNumNewNodes; i++)
    {
        // Parent Rank
        enum EntityRank    tParentRank  = tDecompData.tNewNodeParentRank(i);
        moris::moris_index tParentIndex = tDecompData.tNewNodeParentIndex(i);

        // get the owner processor
        moris::moris_index tOwnerProc = tMeshData.get_entity_owner(tParentIndex,tParentRank);

        // If i own the request keep track of the index
        if(tOwnerProc == tParRank)
        {
            aOwnedRequests.push_back(i);
        }
        else
        {
            moris_index tIndex = aProcRankToIndexInData[tOwnerProc];

            aNotOwnedRequests(tIndex).push_back(i);

        }
    }
}

void
Model::assign_owned_request_identifiers(Decomposition_Data & aDecompData,
                                        Cell<uint> const &   aOwnedRequest,
                                        moris::moris_id &    aNodeInd,
                                        moris::moris_id &    aNodeId)
{
    for(moris::uint i = 0; i < aOwnedRequest.size(); i++)
    {
        moris_index tRequestIndex = aOwnedRequest(i);

        // set the new node index
        aDecompData.tNewNodeIndex(tRequestIndex) = aNodeInd;
        aNodeInd++;

        // set the new node id
        aDecompData.tNewNodeId(tRequestIndex) = aNodeId;
        aNodeId++;

        // increment number of new nodes with set ids (for assertion purposes)
        aDecompData.mNumNewNodesWithIds++;
    }
}


void
Model::setup_outward_requests(Decomposition_Data              const & aDecompData,
                              Cell<Cell<uint>>                const & aNotOwnedRequests,
                              Cell<uint>                      const & aProcRanks,
                              std::unordered_map<moris_id,moris_id> & aProcRankToIndexInData,
                              Cell<Matrix<IndexMat>>                & aOutwardRequests)
{
    // size data
    aOutwardRequests.resize(aProcRanks.size());

    // iterate through the processors we need information from and package the matrix
    for(moris::uint i = 0; i < aProcRanks.size(); i++)
    {
        uint tProcRank = aProcRanks(i);

        MORIS_ASSERT(aProcRankToIndexInData.find(tProcRank) != aProcRankToIndexInData.end(),"Proc rank not in map");
        uint tIndexInData = aProcRankToIndexInData[tProcRank];

        uint tNumRequests = aNotOwnedRequests(tIndexInData).size();

        // size the sending matrix
        // column - request
        //   r0 - parent entity id
        //   r1 - parent entity rank
        //   r2 - Secondary id
        if(tNumRequests > 0)
        {
            aOutwardRequests(i) = moris::Matrix<IndexMat>(3,tNumRequests);
        }

        else
        {
            aOutwardRequests(i) = moris::Matrix<IndexMat>(3,1,MORIS_INDEX_MAX);
        }

        // populate matrix to send;
        for(moris::uint j = 0; j < tNumRequests; j++)
        {
            moris_index     tRequestIndex = aNotOwnedRequests(tIndexInData)(j);
            moris_index     tParentIndex  = aDecompData.tNewNodeParentIndex(tRequestIndex);
            moris_index     tSecondaryId  = aDecompData.tSecondaryIdentifiers(tRequestIndex);
            enum EntityRank tParentRank   = aDecompData.tNewNodeParentRank(tRequestIndex);

            aOutwardRequests(i)(0,j) = mBackgroundMesh.get_glb_entity_id_from_entity_loc_index(tParentIndex,tParentRank);
            aOutwardRequests(i)(1,j) = (moris_index)tParentRank;
            aOutwardRequests(i)(2,j) = tSecondaryId;
        }
    }
}

bool
Model::verify_successful_node_assignment(Decomposition_Data & aDecompData)
{
    for(moris::uint i = 0; i < aDecompData.tNewNodeId.size(); i++)
    {
        if(aDecompData.tNewNodeId(i) == MORIS_INDEX_MAX)
        {
            return false;
        }
    }

    return true;

}

void
Model::send_outward_requests(moris_index            const & aMPITag,
                             Cell<uint>             const & aProcRanks,
                             Cell<Matrix<IndexMat>> & aOutwardRequests)
{
    // Cell of requests
    Cell<MPI_Request> tRequests(aProcRanks.size());
    MPI_Status tStatus;

    // iterate through owned requests and send
    for(moris::uint i = 0; i < aProcRanks.size(); i++)
    {
        tRequests(i) = nonblocking_send(aOutwardRequests(i),aOutwardRequests(i).n_rows(),aOutwardRequests(i).n_cols(),aProcRanks(i),aMPITag);
    }

}

void
Model::inward_receive_requests(
        moris_index            const & aMPITag,
        moris::uint                    aNumRows,
        Cell<Matrix<IndexMat>> &       aReceivedData,
        Cell<uint>             &       aProcRanksReceivedFrom)
{

    // ensure the sizes are correct.
    aReceivedData.resize(0);
    aProcRanksReceivedFrom.resize(0);

    // access the mesh data behind the background mesh
    moris::mtk::Mesh & tMeshData = mBackgroundMesh.get_mesh_data();

    // access the communication table
    Matrix<IdMat> tCommTable = mBackgroundMesh.get_communication_table();

    moris::moris_index tParRank = par_rank();
    moris::uint tCount = 0;
    MPI_Status tStatus;
    for(moris::uint i = 0; i<tCommTable.numel(); i++)
    {
        aReceivedData.push_back(Matrix<IndexMat>(1,1));
        aProcRanksReceivedFrom.push_back(tCommTable(i));
        receive(aReceivedData(tCount),aNumRows, tCommTable(i),aMPITag);
        tCount++;
    }
}

void
Model::inward_receive_request_answers(moris_index            const & aMPITag,
                                      moris::uint            const & aNumRows,
                                      Cell<uint>             const & aProcRanks,
                                      Cell<Matrix<IndexMat>> &       aReceivedRequestAnswers)
{

    MPI_Status tStatus;

    for(moris::uint i = 0; i<aProcRanks.size(); i++)
    {

        bool tFlag = sent_message_exists(aProcRanks(i),aMPITag,tStatus);
        while(tFlag == false)
        {
            tFlag = sent_message_exists(aProcRanks(i),aMPITag,tStatus);
        }

        aReceivedRequestAnswers.push_back(Matrix<IndexMat>(1,1));
        receive(aReceivedRequestAnswers(i),aNumRows, aProcRanks(i),aMPITag);

    }
}

void
Model::handle_received_request_answers(Decomposition_Data & aDecompData,
                                        Cell<Matrix<IndexMat>> const & aRequests,
                                        Cell<Matrix<IndexMat>> const & aRequestAnswers,
                                        moris::moris_id &    aNodeInd,
                                        moris::moris_id &    aNodeId)
{
    // iterate through received data
    for(moris::uint i = 0; i < aRequests.size(); i++)
    {
        uint tNumReceivedReqs = aRequests(i).n_cols();


        // avoid the dummy message
        if(aRequests(i)(0,0) != MORIS_INDEX_MAX)
        {

            // iterate through received requests
            for(moris::uint j = 0; j < tNumReceivedReqs; j++)
            {
                moris_id        tParentId      = aRequests(i)(0,j);
                enum EntityRank tParentRank    = (enum EntityRank) aRequests(i)(1,j);
                moris_id        tSecondaryId   = aRequests(i)(2,j);
                moris_index     tParentInd     = mBackgroundMesh.get_mesh_data().get_loc_entity_ind_from_entity_glb_id(tParentId,tParentRank);
                bool            tRequestExists = false;
                moris_index     tRequestIndex  = MORIS_INDEX_MAX;

                if(aDecompData.mHasSecondaryIdentifier)
                {
                    tRequestExists = aDecompData.request_exists(tParentInd,tSecondaryId,(EntityRank)tParentRank,tRequestIndex);
                }

                else
                {
                    tRequestExists = aDecompData.request_exists(tParentInd,(EntityRank)tParentRank,tRequestIndex);
                }

                if(tRequestExists && aRequestAnswers(i)(j))
                {

                    moris_id tNodeId =aRequestAnswers(i)(j);

                    // meaning the owning processor expected this and gave an answer
                    if(tNodeId < MORIS_ID_MAX && aDecompData.tNewNodeId(tRequestIndex) == MORIS_INDEX_MAX)
                    {
                        aDecompData.tNewNodeId(tRequestIndex) = tNodeId;
                        aDecompData.tNewNodeIndex(tRequestIndex) = aNodeInd;
                        aNodeInd++;

                        // set the new node id
                        aDecompData.tNewNodeId(tRequestIndex) = tNodeId;

                        aDecompData.mNumNewNodesWithIds++;

                    }

                    // This is a hanging
                    else
                    {
                        //                    MORIS_ASSERT(aDecompData.mSubdivisionMethod == Subdivision_Method::NC_REGULAR_SUBDIVISION_HEX8 && (int)tParentRank == 2," The only known case where hanging nodes are allowed to show up is on a face during NC_REGULAR_SUBDIVISION_HEX8");
                        //                    aDecompData.tNewNodeOwner(tRequestIndex) = par_rank();

                        //                    aDecompData.tNewNodeId(tRequestIndex) = tNodeId;
                        //                    aDecompData.tNewNodeIndex(tRequestIndex) = aNodeInd;
                        //                    aNodeInd++;
                        //
                        //                    // set the new node id
                        //                    aDecompData.tNewNodeId(tRequestIndex) = tNodeId;
                        //
                        //                    aDecompData.mNumNewNodesWithIds++;

                    }
                }
                else
                {
                    MORIS_ASSERT(0,"Request does not exist.");
                }
            }
        }
    }
}

void
Model::send_outward_requests_reals(
        moris_index const    & aMPITag,
        Cell<uint>  const    & aProcRanks,
        Cell<Matrix<DDRMat>> & aOutwardRequests)
{
    // iterate through owned requests and send
    for(moris::uint i = 0; i < aProcRanks.size(); i++)
    {
        nonblocking_send(aOutwardRequests(i),aOutwardRequests(i).n_rows(),aOutwardRequests(i).n_cols(),aProcRanks(i),aMPITag);
    }
}

void
Model::inward_receive_requests_reals(
        moris_index const &    aMPITag,
        moris::uint            aNumRows,
        Cell<Matrix<DDRMat>> & aReceivedData,
        Cell<uint>           & aProcRanksReceivedFrom)
{
    {
        moris::moris_index tParRank = par_rank();
        moris::uint tCount = 0;
        MPI_Status tStatus;
        for(moris::uint i = 0; i<(moris::uint)par_size(); i++)
        {
            if((moris_index)i != tParRank)
            {
                // if there is a sent message from a processor go receive it
                if(sent_message_exists(i,aMPITag,tStatus))
                {
                    aReceivedData.push_back(Matrix<DDRMat>(1,1));
                    aProcRanksReceivedFrom.push_back(i);
                    receive(aReceivedData(tCount),aNumRows, i,aMPITag);
                    tCount++;
                }
            }
        }
    }
}

void
Model::return_request_answers_reals(
        moris_index const & aMPITag,
        Cell<Matrix<DDRMat>> const & aRequestAnswers,
        Cell<uint>              const & aProcRanks)
{
    // iterate through owned requests and send
    for(moris::uint i = 0; i < aProcRanks.size(); i++)
    {
        nonblocking_send(aRequestAnswers(i),aRequestAnswers(i).n_rows(),aRequestAnswers(i).n_cols(),aProcRanks(i),aMPITag);
    }
}

void
Model::inward_receive_request_answers_reals(moris_index            const & aMPITag,
        moris::uint            const & aNumRows,
        Cell<uint>             const & aProcRanks,
        Cell<Matrix<DDRMat>> &       aReceivedData)
{
    for(moris::uint i = 0; i<aProcRanks.size(); i++)
    {
        aReceivedData.push_back(Matrix<DDRMat>(1,1));
        receive(aReceivedData(i),aNumRows, aProcRanks(i),aMPITag);
    }
}


void
Model::prepare_request_answers( Decomposition_Data & aDecompData,
                                Cell<Matrix<IndexMat>> const & aReceiveData,
                                Cell<Matrix<IndexMat>>       & aRequestAnswers)
{
    // allocate answer size
    aRequestAnswers.resize(aReceiveData.size());

    // iterate through received data
    for(moris::uint i = 0; i < aReceiveData.size(); i++)
    {
        uint tNumReceivedReqs = aReceiveData(i).n_cols();

        aRequestAnswers(i).resize(1,tNumReceivedReqs);

        aRequestAnswers(i)(0) = MORIS_INDEX_MAX;

        // avoid the dummy message
        if(aReceiveData(i)(0,0) != MORIS_INDEX_MAX)
        {

            // iterate through received requests
            for(moris::uint j = 0; j < tNumReceivedReqs; j++)
            {
                moris_id        tParentId      = aReceiveData(i)(0,j);
                enum EntityRank tParentRank    = (enum EntityRank) aReceiveData(i)(1,j);
                moris_id        tSecondaryId   = aReceiveData(i)(2,j);
                moris_index     tParentInd     = mBackgroundMesh.get_mesh_data().get_loc_entity_ind_from_entity_glb_id(tParentId,tParentRank);
                bool            tRequestExists = false;
                moris_index     tRequestIndex  = MORIS_INDEX_MAX;

                if(aDecompData.mHasSecondaryIdentifier)
                {
                    tRequestExists = aDecompData.request_exists(tParentInd,tSecondaryId,(EntityRank)tParentRank,tRequestIndex);
                }

                else
                {
                    tRequestExists = aDecompData.request_exists(tParentInd,(EntityRank)tParentRank,tRequestIndex);
                }

                if(tRequestExists)
                {
                    moris_id tNodeId =aDecompData.tNewNodeId(tRequestIndex);
                    aRequestAnswers(i)(j) = tNodeId;

                    if(tNodeId == MORIS_ID_MAX)
                    {
                        std::cout<<"tParentId = "<<tParentId<<" | Rank "<<(uint)tParentRank<<std::endl;
                        //                    MORIS_ERROR(0,"Max node");
                    }
                }
                else
                {
                    aRequestAnswers(i)(j) = MORIS_ID_MAX;
                    MORIS_ASSERT(0,"Request does not exist. Need to handle hanging node in reg sub of hex8");
                }
            }
        }


    }

}

void
Model::return_request_answers(  moris_index const & aMPITag,
                                Cell<Matrix<IndexMat>> const & aRequestAnswers,
                                Cell<uint> const & aProcRanks)
{
    // access the mesh data behind the background mesh
    moris::mtk::Mesh & tMeshData = mBackgroundMesh.get_mesh_data();

    // access the communication table
    Matrix<IdMat> tCommTable = mBackgroundMesh.get_communication_table();

    // iterate through owned requests and send
    for(moris::uint i = 0; i < tCommTable.numel(); i++)
    {
        nonblocking_send(aRequestAnswers(i),aRequestAnswers(i).n_rows(),aRequestAnswers(i).n_cols(),tCommTable(i),aMPITag);
    }
}

void
Model::link_background_mesh_to_geometry_objects()
{
    // initialize geometry objects
    moris::Matrix< moris::DDRMat > tNodeCoords = mBackgroundMesh.get_all_node_coordinates_loc_inds();

    mGeometryEngine->initialize_geometry_objects_for_background_mesh_nodes(tNodeCoords.n_rows());
    mGeometryEngine->initialize_geometry_object_phase_values(tNodeCoords);
    mLinkedBackground = true;
}

void
Model::finalize_decomp_in_xtk_mesh(bool aSetPhase)
{
    // Change XTK model decomposition state flag
    mDecomposed = true;

    // Sort the children meshes into groups
    this->sort_children_meshes_into_groups();

    // give each child cell its id (parallel consistent) and index (not parallel consistent)
    this->assign_child_element_identifiers();

    // add child element to local to global map
    this->add_child_elements_to_local_to_global_map();

    // Associate nodes created during decomposition to their child meshes
    this->associate_nodes_created_during_decomp_to_child_meshes();

    // creates mtk cells for all child elements (parent elements are assumed to have mtk cells in the mtk mesh)
    this->create_child_element_mtk_cells();

    // add vertices to child meshes
    this->add_vertices_to_child_meshes();

    // set the glb to loc map for all cells
    this->setup_cell_glb_to_local_map();

    // Compute the child element phase using the geometry engine
    // a case where the phase may not be set is when we only do a
    // non-conformal decomposition
    this->set_element_phases();

    // identify local subphases in child mesh
    this->identify_local_subphase_clusters_in_child_meshes();

    // assign subphase ids
    this->assign_subphase_glob_ids();

    // setup global to local subphase map
    this->setup_glob_to_loc_subphase_map();

}

void
Model::assign_child_element_identifiers()
{
    // Set child element ids and indices
    moris::size_t tNumElementsInCutMesh = mCutMesh.get_num_entities(EntityRank::ELEMENT);

    // Allocate global element ids (these need to be give to the children meshes)
    moris_id    tElementIdOffset = mBackgroundMesh.allocate_entity_ids(tNumElementsInCutMesh, moris::EntityRank::ELEMENT);
    moris_index tElementIndOffset = mBackgroundMesh.get_first_available_index(EntityRank::ELEMENT);

    // set child elements ids in the children meshes which I own and dont share
    Cell<Child_Mesh*> const & tOwnedChildMeshes = mCutMesh.get_owned_child_meshes();
    for(moris::size_t i = 0; i<tOwnedChildMeshes.size(); i++)
    {
        tOwnedChildMeshes(i)->set_child_element_ids(tElementIdOffset);
        tOwnedChildMeshes(i)->set_child_element_inds(tElementIndOffset);
    }

    // set the indices of the child elements in not owned but not the ids
    // set child elements ids in the children meshes which I own and dont share
    Cell<Child_Mesh*> const & tNotOwnedChildMeshes     = mCutMesh.get_not_owned_child_meshes();
    for(moris::size_t i = 0; i<tNotOwnedChildMeshes.size(); i++)
    {
        tNotOwnedChildMeshes(i)->set_child_element_inds(tElementIndOffset);
    }

    // prepare outward requests
    Cell<Cell<moris_id>>        tNotOwnedChildMeshesToProcs;
    Cell<moris::Matrix<IdMat>>  tOwnedParentCellId;
    Cell<moris::Matrix<IdMat>>  tNumOwnedCellIdsOffsets;
    Cell<uint>           tProcRanks;
    std::unordered_map<moris_id,moris_id>  tProcRankToDataIndex;
    this->prepare_child_element_identifier_requests(tNotOwnedChildMeshesToProcs, tOwnedParentCellId, tNumOwnedCellIdsOffsets, tProcRanks,tProcRankToDataIndex);

    // send requests
    moris::uint tMPITag = 141;
    this->send_outward_requests(tMPITag, tProcRanks,tOwnedParentCellId);
    this->send_outward_requests(tMPITag+1, tProcRanks, tNumOwnedCellIdsOffsets);

    barrier();

    // receive requests
    Cell<Matrix<IndexMat>> tReceivedParentCellIds;
    Cell<Matrix<IndexMat>> tReceivedParentCellNumChildren;
    Cell<uint> tProcsReceivedFrom1;
    Cell<uint> tProcsReceivedFrom2;
    this->inward_receive_requests(tMPITag, 1, tReceivedParentCellIds, tProcsReceivedFrom1);
    this->inward_receive_requests(tMPITag+1,1, tReceivedParentCellNumChildren, tProcsReceivedFrom2);

    MORIS_ASSERT(tProcsReceivedFrom1.size() == tProcsReceivedFrom2.size(),"Size mismatch between procs received from child cell ids and number of child cells");
    Cell<Matrix<IndexMat>> tChildIdOffsets;
    this->prepare_child_cell_id_answers(tReceivedParentCellIds,tReceivedParentCellNumChildren,tChildIdOffsets);

    // return information
    this->return_request_answers(tMPITag+2, tChildIdOffsets, tProcsReceivedFrom1);

    // receive the information
    barrier();

    // receive the answers
    Cell<Matrix<IndexMat>> tReceivedChildIdOffsets;
    this->inward_receive_request_answers(tMPITag+2,1,tProcRanks,tReceivedChildIdOffsets);

    // add child cell ids to not owned child meshes
    this->handle_received_child_cell_id_request_answers(tNotOwnedChildMeshesToProcs,tReceivedChildIdOffsets,tElementIndOffset);

    // tell the background mesh about the new first available index
    mBackgroundMesh.update_first_available_index(tElementIndOffset,EntityRank::ELEMENT);

    barrier();

}

void
Model::prepare_child_element_identifier_requests(Cell<Cell<moris_id>>       & aNotOwnedChildMeshesToProcs,
                                                 Cell<moris::Matrix<IdMat>> & aOwnedParentCellId,
                                                 Cell<moris::Matrix<IdMat>> & aNumOwnedCellIdsOffsets,
                                                 Cell<uint>                 & aProcRanks,
                                                 std::unordered_map<moris_id,moris_id> & aProcRankToDataIndex)
{
    // ask owning processor about child element ids
    Cell<Child_Mesh*> const & tNotOwnedChildMeshes     = mCutMesh.get_not_owned_child_meshes();
    Cell<moris_id>    const & tNotOwnedChildMeshOwners = mCutMesh.get_not_owned_child_owners();
    Cell<moris_id>            tCounts(0);
    moris_index tCurrentIndex = 0;

    //set up the porcs
    // access the mesh data behind the background mesh
    moris::mtk::Mesh & tMeshData = mBackgroundMesh.get_mesh_data();

    // access the communication table
    Matrix<IdMat> tCommTable = mBackgroundMesh.get_communication_table();

    for(moris::size_t i = 0; i < tCommTable.numel(); i++)
    {
            aProcRankToDataIndex[tCommTable(i)] = tCurrentIndex;
            aProcRanks.push_back(tCommTable(i));
            aNotOwnedChildMeshesToProcs.push_back(Cell<moris_id>(0));
            tCounts.push_back(0);
            tCurrentIndex++;
    }


    // sort child meshes by owner
    for(moris::size_t i = 0; i < tNotOwnedChildMeshes.size(); i++)
    {
        moris_index tOwnerProc = tNotOwnedChildMeshOwners(i);
        moris_index tProcDataIndex = aProcRankToDataIndex[tOwnerProc];
        aNotOwnedChildMeshesToProcs(tProcDataIndex).push_back(i);

        moris_id tParentCellId = mBackgroundMesh.get_mesh_data().get_glb_entity_id_from_entity_loc_index(tNotOwnedChildMeshes(i)->get_parent_element_index(),EntityRank::ELEMENT);
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
            Child_Mesh* tCM = tNotOwnedChildMeshes(aNotOwnedChildMeshesToProcs(i)(j));

            aOwnedParentCellId(i)(j)      = mBackgroundMesh.get_mesh_data().get_glb_entity_id_from_entity_loc_index(tCM->get_parent_element_index(),EntityRank::ELEMENT);
            aNumOwnedCellIdsOffsets(i)(j) = tCM->get_num_entities(EntityRank::ELEMENT);
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

void
Model::prepare_child_cell_id_answers(Cell<Matrix<IndexMat>> & aReceivedParentCellIds,
                                     Cell<Matrix<IndexMat>> & aReceivedParentCellNumChildren,
                                     Cell<Matrix<IndexMat>> & aChildCellIdOffset)
{
    MORIS_ASSERT(aReceivedParentCellIds.size() == aReceivedParentCellNumChildren.size(),"Mismatch in received parent cell ids and received parent cell number of children");

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
                moris_index tParentCellIndex = mBackgroundMesh.get_mesh_data().get_loc_entity_ind_from_entity_glb_id(tParentId,EntityRank::ELEMENT);

                // get child mesh
                MORIS_ASSERT(mBackgroundMesh.entity_has_children(tParentCellIndex,EntityRank::ELEMENT),"Request is made for child element ids on a parent cell not intersected");
                moris_index tCMIndex = mBackgroundMesh.child_mesh_index(tParentCellIndex,EntityRank::ELEMENT);
                Child_Mesh & tCM = mCutMesh.get_child_mesh(tCMIndex);

                MORIS_ASSERT(par_rank() == mBackgroundMesh.get_mesh_data().get_entity_owner(tParentCellIndex,EntityRank::ELEMENT),"I dont own this entity that had info requestsed.");

                // place in return data
                MORIS_ASSERT(tCM.get_num_entities(EntityRank::ELEMENT) == (uint)aReceivedParentCellNumChildren(i)(j),"Number of child cells in child mesh do not match number on other processor");

                // since hmr ownership is not correct
                if(tCM.get_element_ids().numel()>0)
                {
                    aChildCellIdOffset(i)(j) = tCM.get_element_ids()(0);
                }
            }
        }
        else
        {
            aChildCellIdOffset(i)(0) =   MORIS_INDEX_MAX;
        }

    }
}

void
Model::handle_received_child_cell_id_request_answers(
                                        Cell<Cell<moris_index>> const & aChildMeshesInInNotOwned,
                                        Cell<Matrix<IndexMat>> const & aReceivedChildCellIdOffset,
                                        moris::moris_id              & aCellInd)
{
    Cell<Child_Mesh*> const & tNotOwnedChildMeshes     = mCutMesh.get_not_owned_child_meshes();

    // iterate through received data
    for(moris::uint i = 0; i < aChildMeshesInInNotOwned.size(); i++)
    {
        uint tNumReceivedReqs = aChildMeshesInInNotOwned(i).size();

        // iterate through received requests
        for(moris::uint j = 0; j < tNumReceivedReqs; j++)
        {
            moris_id tChildMeshInNotOwned = aChildMeshesInInNotOwned(i)(j);
            Child_Mesh* tCM = tNotOwnedChildMeshes(tChildMeshInNotOwned);
            moris_id tChildCellFirstId = aReceivedChildCellIdOffset(i)(j);

            tCM->set_child_element_ids(tChildCellFirstId);
        }
    }
}

void
Model::add_child_elements_to_local_to_global_map()
{
    // get all child element ids and indexes
    Matrix<IndexMat> tChildElementInds = mCutMesh.get_all_element_inds();
    Matrix<IndexMat> tChildElementIds  = mCutMesh.get_all_element_ids();

    mBackgroundMesh.add_cells_to_global_to_local_map(tChildElementInds,tChildElementIds);


}

void
Model::sort_children_meshes_into_groups()
{
    // background mesh data
    moris::mtk::Mesh & tMeshData = mBackgroundMesh.get_mesh_data();

    // my proc rank
    moris_index tProcRank = par_rank();

    // number of children meshes
    uint tNumChildrenMeshes = mCutMesh.get_num_child_meshes();

    // allocate data
    Cell<Child_Mesh*>   tOwnedChildrenMeshes(tNumChildrenMeshes);
    Cell<Child_Mesh*>   tNotOwnedChildrenMeshes(tNumChildrenMeshes);
    Cell<moris_id>      tNotOwnedOwningProc(tNumChildrenMeshes);

    // keep track of the number in each group
    uint tOwnedCount    = 0;
    uint tNotOwnedCount = 0;

    for(moris::size_t i = 0; i<mCutMesh.get_num_child_meshes(); i++)
    {
        Child_Mesh & tChildMesh = mCutMesh.get_child_mesh(i);

        moris_index tParentCellInd = tChildMesh.get_parent_element_index();

        // get owner of parent cell
        moris_index tOwnerProc = tMeshData.get_entity_owner(tParentCellInd,EntityRank::ELEMENT);

        // if this processor does not own the element add it to the not owned shared list
        if(tOwnerProc != tProcRank)
        {
            tNotOwnedChildrenMeshes(tNotOwnedCount) = &tChildMesh;
            tNotOwnedOwningProc(tNotOwnedCount) = tOwnerProc;
            tNotOwnedCount++;
        }

        else
        {
            tOwnedChildrenMeshes(tOwnedCount) = &tChildMesh;
            tOwnedCount++;
        }
    }

    // size out extra space
    tOwnedChildrenMeshes.resize(tOwnedCount);
    tNotOwnedChildrenMeshes.resize(tNotOwnedCount);
    tNotOwnedOwningProc.resize(tNotOwnedCount);

    // add to cut mesh
    mCutMesh.add_child_mesh_groups( tOwnedChildrenMeshes, tNotOwnedChildrenMeshes, tNotOwnedOwningProc);
}

void
Model::associate_nodes_created_during_decomp_to_child_meshes()
{
    // Initialize the data in the XTK mesh
    mBackgroundMesh.allocate_external_node_to_child_mesh_associations();

    // Number of children meshes
    size_t tNumCM = mCutMesh.get_num_child_meshes();
    for(size_t i = 0 ; i < tNumCM; i++)
    {
        // Get reference to the child mesh
        Child_Mesh const & tChildMesh = mCutMesh.get_child_mesh(i);

        // Get reference to the nods in the child mesh node indices
        moris::Matrix<moris::IndexMat> const & tNodeIndices = tChildMesh.get_node_indices();

        // Associate these node indices with their child mesh index
        mBackgroundMesh.associate_external_nodes_to_child_mesh(i,tNodeIndices);
    }
}

void
Model::set_element_phases()
{
    // Set element phase indices
    mBackgroundMesh.initialize_element_phase_indices(this->get_num_elements_total());

    moris::size_t tNumElem = mBackgroundMesh.get_mesh_data().get_num_entities(EntityRank::ELEMENT);

    for(moris::size_t i = 0; i<tNumElem; i++)
    {
        if(mBackgroundMesh.entity_has_children(i,EntityRank::ELEMENT))
        {
            moris::size_t tChildMeshIndex = mBackgroundMesh.child_mesh_index(i,EntityRank::ELEMENT);

            Child_Mesh & tChildMesh = mCutMesh.get_child_mesh(tChildMeshIndex);

            moris::Matrix< moris::IndexMat > tElemToNode = tChildMesh.get_element_to_node();

            moris::Matrix< moris::IndexMat > const & tElemInds  = tChildMesh.get_element_inds();

            tChildMesh.initialize_element_phase_mat();

            moris::size_t tNumElem = tChildMesh.get_num_entities(EntityRank::ELEMENT);

            for( moris::size_t j = 0; j<tNumElem; j++)
            {
                moris::size_t tElemPhaseIndex = determine_element_phase_index(j,tElemToNode);
                mBackgroundMesh.set_element_phase_index(tElemInds(0,j),tElemPhaseIndex);
                tChildMesh.set_element_phase_index(j,tElemPhaseIndex);
            }
        }

        else
        {
            moris::Matrix< moris::IndexMat > tElementNodes = mBackgroundMesh.get_mesh_data().get_entity_connected_to_entity_loc_inds(i,moris::EntityRank::ELEMENT, moris::EntityRank::NODE);

            if(iscol(tElementNodes))
            {
                tElementNodes = trans(tElementNodes);
            }

            moris::size_t tElemPhaseIndex = determine_element_phase_index(0,tElementNodes);

            mBackgroundMesh.set_element_phase_index(i,tElemPhaseIndex);
        }


    }
}

void Model::set_downward_inheritance()
{
    moris::size_t tNumChildMesh = mCutMesh.get_num_child_meshes();
    Cell<std::pair<moris::moris_index,moris::moris_index>> tXTKElementToCutMeshPairs(tNumChildMesh);

    for(moris::size_t iMesh = 0; iMesh<tNumChildMesh; iMesh++)
    {
        tXTKElementToCutMeshPairs(iMesh) = std::pair<moris::moris_index,moris::moris_index> (mCutMesh.get_parent_element_index(iMesh),iMesh);
    }

    mBackgroundMesh.register_new_downward_inheritance(tXTKElementToCutMeshPairs);
}


void  Model::run_first_cut_routine(enum TemplateType const &          aTemplateType,
                                   moris::uint                        aGeomIndex,
                                   moris::size_t const &              aNumNodesPerElement,
                                   moris::Matrix< moris::IndexMat > & aActiveChildMeshIndices,
                                   moris::Matrix< moris::IndexMat > & aNewPairBool)
{
    // Note this method is independent of node ids for this reason Background_Mesh is not given the node Ids during this subdivision
    moris::mtk::Mesh & tXTKMeshData = mBackgroundMesh.get_mesh_data();

    // Package up node to element connectivity
    moris::moris_index tParentElementIndex = MORIS_INDEX_MAX;
    moris::size_t tNumElements = mBackgroundMesh.get_num_entities(EntityRank::ELEMENT);

    moris::Matrix< moris::IndexMat > tElementToNodeConnInd (tNumElements, aNumNodesPerElement);
    moris::Matrix< moris::IndexMat > tElementToNodeVector (1, aNumNodesPerElement);

    for (moris::size_t i = 0; i < tNumElements; i++)
    {
        tElementToNodeVector = tXTKMeshData.get_entity_connected_to_entity_loc_inds(i, moris::EntityRank::ELEMENT, moris::EntityRank::NODE);

        for(moris::uint j = 0; j < aNumNodesPerElement; j++)
        {
            tElementToNodeConnInd(i,j) = tElementToNodeVector(j);
        }

    }

    // Get the Node Coordinates
    moris::Matrix< moris::DDRMat > tAllNodeCoords = mBackgroundMesh.get_all_node_coordinates_loc_inds();

    // Intersected elements are flagged via the Geometry_Engine
    Cell<moris::ge::GEN_Geometry_Object> tGeoObjects;
    mGeometryEngine->is_intersected(tAllNodeCoords, tElementToNodeConnInd, 0,tGeoObjects);

    // Count number intersected
    moris::size_t tIntersectedCount = tGeoObjects.size();

    // Loop over and determine how many new meshes that need to be registered (Avoids dynamic allocation in the child mesh)
    // Also register active mesh pairs
    Cell<std::pair<moris::moris_index,moris::moris_index>> tNewChildElementPair;
    aNewPairBool = moris::Matrix< moris::IndexMat >(1,tIntersectedCount,0);
    tNewChildElementPair.reserve(tIntersectedCount);

    moris::size_t tNumNewChildMeshes = 0;
    moris::moris_index tNewIndex = 0;
    aActiveChildMeshIndices.resize(1,tIntersectedCount);
    for (moris::size_t j = 0; j < tIntersectedCount; j++)
    {
        tParentElementIndex = tGeoObjects(j).get_parent_entity_index();
        if(!mBackgroundMesh.entity_has_children(tParentElementIndex,EntityRank::ELEMENT))
        {
            tNewIndex = tNumNewChildMeshes+mCutMesh.get_num_child_meshes();
            tNewChildElementPair.push_back( std::pair<moris::moris_index,moris::moris_index>(tParentElementIndex, tNewIndex));
            aActiveChildMeshIndices(0,j) = tNewIndex;
            tNumNewChildMeshes++;
        }

        else
        {
            aActiveChildMeshIndices(0,j) = mBackgroundMesh.child_mesh_index(tParentElementIndex,EntityRank::ELEMENT);
            aNewPairBool(0,j) = 1;
        }
    }


    // Add the downward pair to the mesh for all the newly created element pairs
    mBackgroundMesh.register_new_downward_inheritance(tNewChildElementPair);

    // Allocate space for more simple meshes in XTK mesh
    mCutMesh.inititalize_new_child_meshes(tNumNewChildMeshes, mModelDimension);


    moris::Matrix< moris::IndexMat > tPlaceHolder(1, 1);
    for (moris::size_t j = 0; j < tIntersectedCount; j++)
    {
        if(aNewPairBool(0,j) == 0)
        {
            tParentElementIndex = tGeoObjects(j).get_parent_entity_index();

            // Get information to provide ancestry
            // This could be replaced with a proper topology implementation that knows faces, edges based on parent element nodes
            Matrix< IndexMat > tNodetoElemConnVec = tElementToNodeConnInd.get_row(tParentElementIndex);
            Matrix< IndexMat > tEdgetoElemConnInd = tXTKMeshData.get_entity_connected_to_entity_loc_inds(tParentElementIndex, moris::EntityRank::ELEMENT, moris::EntityRank::EDGE);
            Matrix< IndexMat > tFacetoElemConnInd = tXTKMeshData.get_entity_connected_to_entity_loc_inds(tParentElementIndex, moris::EntityRank::ELEMENT, moris::EntityRank::FACE);
            Matrix< IndexMat > tElementMat        = {{tParentElementIndex}};


            // Set parent element, nodes, and entity ancestry
            moris::Matrix< moris::IndexMat > tElemToNodeIndices(tNodetoElemConnVec);

            Cell<moris::Matrix< moris::IndexMat >> tAncestorInformation = {tPlaceHolder, tEdgetoElemConnInd, tFacetoElemConnInd, tElementMat};
            mCutMesh.initialize_new_mesh_from_parent_element(aActiveChildMeshIndices(0,j), aTemplateType, tNodetoElemConnVec, tAncestorInformation);

            // add node ids
            mBackgroundMesh.convert_loc_entity_ind_to_glb_entity_ids(EntityRank::NODE,tNodetoElemConnVec);

            // get child mesh
            moris::moris_index tCMIndex = mBackgroundMesh.child_mesh_index(tParentElementIndex,EntityRank::ELEMENT);

            // get child mesh
            Child_Mesh & tChildMesh = mCutMesh.get_child_mesh(tCMIndex);
            tChildMesh.add_new_geometry_interface(aGeomIndex);

            // add node ids
            tChildMesh.add_node_ids(tNodetoElemConnVec);

        }
        else
        {
            // get parent element index
            tParentElementIndex = tGeoObjects(j).get_parent_entity_index();

            // get child mesh
            moris::moris_index tCMIndex = mBackgroundMesh.child_mesh_index(tParentElementIndex,EntityRank::ELEMENT);

            // get child mesh
            Child_Mesh & tChildMesh = mCutMesh.get_child_mesh(tCMIndex);
            tChildMesh.add_new_geometry_interface(aGeomIndex);
        }
    }
}

void
Model::create_child_element_mtk_cells()
{

    moris::uint tNumChildMeshes = mCutMesh.get_num_child_meshes();

    moris::mtk::Mesh const & tMeshData = mBackgroundMesh.get_mesh_data();

    for(moris::uint i=0; i<tNumChildMeshes; i++)
    {
        Child_Mesh & tChildMesh = mCutMesh.get_child_mesh(i);

        moris::moris_index tOwnerProc = tMeshData.get_entity_owner(tChildMesh.get_parent_element_index(),EntityRank::ELEMENT);

        moris::Matrix< moris::IdMat >    const & tElementIds  = tChildMesh.get_element_ids();
        moris::Matrix< moris::IndexMat > const & tElementInds = tChildMesh.get_element_inds();
        // Iterate over elements
        for(moris::uint j = 0; j<tElementIds.numel(); j++)
        {
            mBackgroundMesh.add_child_element_to_mtk_cells(tElementInds(j),tElementIds(j),tOwnerProc,j, &tChildMesh);
        }
    }
}

void
Model::add_vertices_to_child_meshes()
{
  moris::uint tNumChildMeshes = mCutMesh.get_num_child_meshes();

  for(moris::uint i=0; i<tNumChildMeshes; i++)
  {
      Child_Mesh & tChildMesh = mCutMesh.get_child_mesh(i);
      Cell<moris::mtk::Vertex const *> tVertices = mBackgroundMesh.get_mtk_vertices(tChildMesh.get_node_indices());
      tChildMesh.add_vertices(tVertices);

  }

}

void
Model::setup_cell_glb_to_local_map()
{
    for(moris::uint i = 0; i < this->get_num_elements_total(); i++)
    {
        moris_id tId = mBackgroundMesh.get_glb_entity_id_from_entity_loc_index((moris_index)i,EntityRank::ELEMENT);
        MORIS_ASSERT(mCellGlbToLocalMap.find(tId) == mCellGlbToLocalMap.end(),"Id already in map");
        mCellGlbToLocalMap[tId] = (moris_index) i;
    }
}


void
Model::identify_local_subphase_clusters_in_child_meshes()
{

    // get the number of children meshes
    moris::size_t tNumChildMeshes =  mCutMesh.get_num_child_meshes();

    // first subphase index (this is incremented by child meshes on call to set_elemental_subphase)
    // (proc local subphase index, never global)
    moris::moris_index tSubPhaseIndex = mBackgroundMesh.get_mesh_data().get_num_entities(EntityRank::ELEMENT);

    // iterate over children meshes and perform local flood-fill
    for(moris::size_t i = 0; i<tNumChildMeshes; i++)
    {
        // Get child mesh index
        Child_Mesh & tChildMesh = mCutMesh.get_child_mesh(i);

        // Perform local flood-fill on child mesh to identify subphase
        moris::Matrix< moris::IndexMat > tLocalFloodFill = local_child_mesh_flood_fill(tChildMesh);

        // Set the local floodfill data as the elemental subphase values in the child mesh
        // The child mesh then sorts the elements into bins
        tChildMesh.set_elemental_subphase(tSubPhaseIndex,tLocalFloodFill);
    }

    // tell the cut mesh how many subphases there are
    mCutMesh.set_num_subphases(tSubPhaseIndex);

    // tell the cut mesh to setup subphase to child mesh connectivity
    mCutMesh.setup_subphase_to_child_mesh_connectivity();
}

void
Model::assign_subphase_glob_ids()
{
    // Get the number of subphases
    moris_id tNumSubphases = (moris_id)mCutMesh.get_num_subphases();

    // Allocate global element ids starting at the maximum id in the background mesh (these need to be give to the children meshes)
    moris::moris_id tSubphaseIdOffset = mBackgroundMesh.allocate_entity_ids(tNumSubphases, EntityRank::ELEMENT);

    // set subphase ids in the children meshes which I own
    Cell<Child_Mesh*> const & tOwnedChildMeshes = mCutMesh.get_owned_child_meshes();
    for(moris::size_t i = 0; i<tOwnedChildMeshes.size(); i++)
    {
        moris_id tCellId = mBackgroundMesh.get_mesh_data().get_glb_entity_id_from_entity_loc_index(tOwnedChildMeshes(i)->get_parent_element_index(),EntityRank::ELEMENT);

        // iterate through subphase ids
        tOwnedChildMeshes(i)->set_subphase_id(0,tCellId);
        for(moris::uint j = 1; j <tOwnedChildMeshes(i)->get_num_subphase_bins(); j++)
        {
            tOwnedChildMeshes(i)->set_subphase_id(j,tSubphaseIdOffset);
            tSubphaseIdOffset++;
        }
    }

    // prepare outward requests
    Cell<Cell<moris_id>>        tNotOwnedSubphasesToProcs;
    Cell<Cell<moris_id>>        tCMSubphaseIndices;
    Cell<moris::Matrix<IdMat>>  tParentCellIds;
    Cell<moris::Matrix<IdMat>>  tChildCellIds;
    Cell<uint>                  tProcRanks;
    std::unordered_map<moris_id,moris_id>  tProcRankToDataIndex;
    this->prepare_subphase_identifier_requests(tNotOwnedSubphasesToProcs, tCMSubphaseIndices, tParentCellIds, tChildCellIds, tProcRanks,tProcRankToDataIndex);

    // send requests
    moris::uint tMPITag = 221;
    this->send_outward_requests(tMPITag, tProcRanks,tParentCellIds);
    this->send_outward_requests(tMPITag+1, tProcRanks, tChildCellIds);

    barrier();

    // receive requests
    Cell<Matrix<IndexMat>> tReceivedParentCellIds;
    Cell<Matrix<IndexMat>> tFirstChildCellIds;
    Cell<uint> tProcsReceivedFrom1;
    Cell<uint> tProcsReceivedFrom2;
    this->inward_receive_requests(tMPITag, 1, tReceivedParentCellIds, tProcsReceivedFrom1);
    this->inward_receive_requests(tMPITag+1,1, tFirstChildCellIds, tProcsReceivedFrom2);
    MORIS_ASSERT(tProcsReceivedFrom1.size() == tProcsReceivedFrom2.size(),"Size mismatch between procs received from child cell ids and number of child cells");

    // prepare answers
    Cell<Matrix<IndexMat>> tSubphaseIds;
    this->prepare_subphase_id_answers(tReceivedParentCellIds,tFirstChildCellIds,tSubphaseIds);

    // return information
    this->return_request_answers(tMPITag+2, tSubphaseIds, tProcsReceivedFrom1);

    barrier();

    // receive the answers
    Cell<Matrix<IndexMat>> tReceivedSubphaseIds;
    this->inward_receive_request_answers(tMPITag+2,1,tProcRanks,tReceivedSubphaseIds);

    // add child cell ids to not owned child meshes
    this->handle_received_subphase_id_request_answers(tNotOwnedSubphasesToProcs,tCMSubphaseIndices,tReceivedSubphaseIds);

    barrier();

}

void
Model::prepare_subphase_identifier_requests(Cell<Cell<moris_id>>       & aNotOwnedSubphasesToProcs,
                                            Cell<Cell<moris_id>>       & aSubphaseCMIndices,
                                            Cell<moris::Matrix<IdMat>> & aParentCellIds,
                                            Cell<moris::Matrix<IdMat>> & aChildCellIds,
                                            Cell<uint>                 & aProcRanks,
                                            std::unordered_map<moris_id,moris_id> & aProcRankToDataIndex)
{
    // access the mesh data behind the background mesh
    moris::mtk::Mesh & tMeshData = mBackgroundMesh.get_mesh_data();

    // access the communication table
    Matrix<IdMat> tCommTable = mBackgroundMesh.get_communication_table();

    // resize proc ranks and setup map to comm table
    aProcRanks.resize(tCommTable.numel());
    for(moris::uint i = 0; i <tCommTable.numel(); i++)
    {
        aProcRankToDataIndex[tCommTable(i)] = i;
        aProcRanks(i) = (tCommTable(i));
        aNotOwnedSubphasesToProcs.push_back(Cell<moris_id>(0));
        aSubphaseCMIndices.push_back(Cell<moris_id>(0));
    }


    // ask owning processor about child element ids
    Cell<Child_Mesh*> const & tNotOwnedChildMeshes     = mCutMesh.get_not_owned_child_meshes();
    Cell<moris_id>    const & tNotOwnedChildMeshOwners = mCutMesh.get_not_owned_child_owners();
    moris_index tCurrentIndex = 0;

    // sort child meshes by owner
    for(moris::size_t i = 0; i < tNotOwnedChildMeshes.size(); i++)
    {
        moris_index tOwnerProc = tNotOwnedChildMeshOwners(i);

        moris_index tProcDataIndex = aProcRankToDataIndex[tOwnerProc];

        // iterate through subphases
        Child_Mesh* tCM = tNotOwnedChildMeshes(i);

        for(moris::uint iS = 0; iS < tCM->get_num_subphase_bins();iS++)
        {
            aNotOwnedSubphasesToProcs(tProcDataIndex).push_back(i);
            aSubphaseCMIndices(tProcDataIndex).push_back(iS);
        }

    }

    aParentCellIds.resize(aNotOwnedSubphasesToProcs.size());
    aChildCellIds.resize(aNotOwnedSubphasesToProcs.size());

    // iterate through procs and child meshes shared with that processor
    for(moris::size_t i = 0; i < aNotOwnedSubphasesToProcs.size(); i++)
    {
        // number of child meshes shared with this processor
        moris::uint tNumCM = aNotOwnedSubphasesToProcs(i).size();

        // allocate matrix
        aParentCellIds(i).resize(1,tNumCM);
        aChildCellIds(i).resize(1,tNumCM);

        if(tNumCM == 0)
        {
            aParentCellIds(i).resize(1,1);
            aChildCellIds(i).resize(1,1);
            aParentCellIds(i)(0) = MORIS_INDEX_MAX;
            aChildCellIds(i)(0) = MORIS_INDEX_MAX;
        }

        for(moris::uint j = 0; j < tNumCM; j++)
        {
            Child_Mesh* tCM = tNotOwnedChildMeshes(aNotOwnedSubphasesToProcs(i)(j));

            // subphase groups
            Cell<moris::Matrix< moris::IndexMat >> const & tSubphaseClusters = tCM->get_subphase_groups();

            // element ids
            Matrix<IdMat> const & tCellIds  = tCM->get_element_ids();

            // subphase index
            moris_index tSPIndex = aSubphaseCMIndices(i)(j);

            // cell index
            moris_index tCMCellInd = tSubphaseClusters(tSPIndex)(0);

            aParentCellIds(i)(j) = mBackgroundMesh.get_glb_entity_id_from_entity_loc_index(tCM->get_parent_element_index(),EntityRank::ELEMENT);
            aChildCellIds(i)(j)  = tCellIds(tCMCellInd);
        }
    }
}

void
Model::prepare_subphase_id_answers(Cell<Matrix<IndexMat>> & aReceivedParentCellIds,
                                   Cell<Matrix<IndexMat>> & aFirstChildCellIds,
                                   Cell<Matrix<IndexMat>> & aSubphaseIds)
{
    MORIS_ASSERT(aReceivedParentCellIds.size() == aFirstChildCellIds.size(),"Mismatch in received parent cell ids and received parent cell number of children");

    // allocate answer size
    aSubphaseIds.resize(aReceivedParentCellIds.size());

    // iterate through received data
    for(moris::uint i = 0; i < aReceivedParentCellIds.size(); i++)
    {
        uint tNumReceivedReqs = aReceivedParentCellIds(i).n_cols();

        aSubphaseIds(i).resize(1,tNumReceivedReqs);

        if(aReceivedParentCellIds(i)(0,0) != MORIS_INDEX_MAX)
        {
            // iterate through received requests
            for(moris::uint j = 0; j < tNumReceivedReqs; j++)
            {
                // parent cell information
                moris_id tParentId           = aReceivedParentCellIds(i)(0,j);
                moris_index tParentCellIndex = mBackgroundMesh.get_mesh_data().get_loc_entity_ind_from_entity_glb_id(tParentId,EntityRank::ELEMENT);

                // Child cell information
                moris_id tChildCellId = aFirstChildCellIds(i)(j);

                // get child mesh
                MORIS_ASSERT(mBackgroundMesh.entity_has_children(tParentCellIndex,EntityRank::ELEMENT),"Request is made for child element ids on a parent cell not intersected");
                moris_index tCMIndex = mBackgroundMesh.child_mesh_index(tParentCellIndex,EntityRank::ELEMENT);
                Child_Mesh & tCM = mCutMesh.get_child_mesh(tCMIndex);


                // figure out which subphase this child cell belongs to in the given child mesh
                Matrix<IdMat> const & tChildMeshCellIds = tCM.get_element_ids();
                moris_id tSubphaseId = MORIS_ID_MAX;
                for(moris::uint iCM = 0; iCM < tChildMeshCellIds.numel(); iCM++)
                {
                    if(tChildMeshCellIds(iCM) == tChildCellId)
                    {
                        tSubphaseId = tCM.get_element_subphase_id(iCM);
                    }
                }

                MORIS_ERROR(tSubphaseId!= MORIS_ID_MAX,"Child cell id not found in child mesh");

                // place in return data
                aSubphaseIds(i)(j) = tSubphaseId;
            }
        }
    }
}

void
Model::handle_received_subphase_id_request_answers(
                                        Cell<Cell<moris_index>>    const & aChildMeshesInNotOwned,
                                        Cell<Cell<moris_index>>    const & aCMSubphaseIndices,
                                        Cell<Matrix<IndexMat>>     const & aReceivedSubphaseIds)
{
    Cell<Child_Mesh*> const & tNotOwnedChildMeshes     = mCutMesh.get_not_owned_child_meshes();

    // iterate through received data
    for(moris::uint i = 0; i < aChildMeshesInNotOwned.size(); i++)
    {
        uint tNumReceivedReqs = aChildMeshesInNotOwned(i).size();

        // iterate through received requests
        for(moris::uint j = 0; j < tNumReceivedReqs; j++)
        {
            moris_id tChildMeshInNotOwned = aChildMeshesInNotOwned(i)(j);
            Child_Mesh* tCM = tNotOwnedChildMeshes(tChildMeshInNotOwned);
            moris_id tSubphaseIndex = aCMSubphaseIndices(i)(j);
            moris_id tSubphaseId = aReceivedSubphaseIds(i)(j);

            tCM->set_subphase_id(tSubphaseIndex,tSubphaseId);
        }
    }
}


void
Model::setup_glob_to_loc_subphase_map()
{

    for(moris::uint i = 0; i < mCutMesh.get_num_subphases(); i++)
    {
        moris_id tSubphaseId = this->get_subphase_id((moris_id)i);
        MORIS_ASSERT(mGlobalToLocalSubphaseMap.find(tSubphaseId) == mGlobalToLocalSubphaseMap.end(),"Subphase id already in map");
        mGlobalToLocalSubphaseMap[tSubphaseId] = i;
    }
}

// ----------------------------------------------------------------------------------
// Sensitivity Source code
// ----------------------------------------------------------------------------------
void
Model::compute_sensitivity()
{
    // Start the clock
    std::clock_t start = std::clock();

    // verify the state of the xtk model
    MORIS_ERROR(mDecomposed,"Prior to computing sensitivity, the decomposition process must be called");
    MORIS_ERROR(!mConvertedToTet10s,"Prior to computing sensitivity, the convert tet4 to tet10 process was called");
    MORIS_ERROR(!mSensitivity,"Calling compute interface sensitivity twice is not supported");

    // Compute interface sensitivity
    compute_interface_sensitivity_internal();

    // Change the sensitivity computation state flag
    mSensitivity = true;

    if(moris::par_rank() == 0 && mVerbose)
    {
        std::cout<<"XTK: Sensitivity computation completed in " <<(std::clock() - start) / (double)(CLOCKS_PER_SEC)<<" s."<<std::endl;
    }
}

// ----------------------------------------------------------------------------------
// Unzipping Child Mesh Source code
// ----------------------------------------------------------------------------------
void
Model::unzip_child_mesh()
{
    // start the clock
    std::clock_t start = std::clock();

    MORIS_ERROR(mDecomposed,"Prior to unzip_child_mesh, the decomposition process must be called");

    // unzip the interface
    unzip_child_mesh_internal();

    mUnzipped = true;
    if(moris::par_rank() == 0  && mVerbose)
    {
        std::cout<<"XTK: Child mesh unzipping completed in "<< (std::clock() - start) / (double)(CLOCKS_PER_SEC)<<" s."<<std::endl;
    }
}

void
Model::unzip_child_mesh_internal()
{
    // get the number of children meshes
    moris::size_t tNumChildMeshes =  mCutMesh.get_num_child_meshes();

    for(moris::size_t i = 0; i<tNumChildMeshes; i++)
    {
        // Get child mesh index
        //        Child_Mesh & tChildMesh = mCutMesh.get_child_mesh(i);
    }
}

// ----------------------------------------------------------------------------------
// Unzipping Interface Source code
// ----------------------------------------------------------------------------------
void
Model::unzip_interface()
{
    // start the clock
    std::clock_t start = std::clock();

    // unzip the interface
    unzip_interface_internal();

    mUnzipped = true;
    if(moris::par_rank() == 0  && mVerbose)
    {
        std::cout<<"XTK: Interface unzipping completed in "<< (std::clock() - start) / (double)(CLOCKS_PER_SEC)<<" s."<<std::endl;
    }
}

void
Model::unzip_interface_internal()
{
    // Get the number of geometries (we need to unzip each interface)
    uint tNumGeoms = mGeometryEngine->get_num_geometries();

    // Get the interface nodes (wrt all geometries)
    Cell<moris::Matrix<moris::IndexMat>> tAllInterfaceNodeIds = mBackgroundMesh.get_interface_nodes_glob_ids();
    Cell<moris::Matrix<moris::IndexMat>> tAllInterfaceNodeInds = mBackgroundMesh.get_interface_nodes_loc_inds();

    // Keep count of which interface node index in matrices were in
    for(uint iG = 0; iG<tNumGeoms; iG++)
    {
        // Interface node wrt geometry iG
        moris::Matrix<moris::IndexMat> const & tInterfaceNodeIds = tAllInterfaceNodeIds(iG);
        moris::Matrix<moris::IndexMat> const & tInterfaceNodeInds = tAllInterfaceNodeInds(iG);

        // Number of interface nodes wrt geometry iG (and check the sizes of the interface information)
        uint tNumInterfaceNodes = tInterfaceNodeIds.numel();
        MORIS_ASSERT(tInterfaceNodeIds.numel() == tInterfaceNodeInds.numel(), "Interface Ids and Indices dimension mismatch");

        // Assign node ids and indices ( row - 0 Node ids, row 1 - New node ids)
        moris::Matrix<moris::IdMat> tNewUnzippedNodeIds((size_t)tNumInterfaceNodes);
        moris::Matrix<moris::IndexMat> tNewUnzippedNodeInds((size_t)tNumInterfaceNodes);
        this->unzip_interface_internal_assign_node_identifiers(tNumInterfaceNodes,tNewUnzippedNodeInds,tNewUnzippedNodeIds);

        // Add new nodes to the mesh (as a copy of the existing node)
        mBackgroundMesh.batch_create_new_nodes_as_copy_of_other_nodes(tInterfaceNodeInds,tNewUnzippedNodeIds,tNewUnzippedNodeInds);

        // Allocate space in background mesh interface node flags
        mBackgroundMesh.allocate_space_in_interface_node_flags(tNumInterfaceNodes, mGeometryEngine->get_num_geometries());

        // Mark the newly created nodes as interface nodes
        mBackgroundMesh.mark_nodes_as_interface_node_loc_inds(tNewUnzippedNodeInds,iG);

        // Link the new nodes to the geometry object of the node they were copied to
        mGeometryEngine->link_new_nodes_to_existing_geometry_objects(tInterfaceNodeInds,tNewUnzippedNodeInds);

        // unzip_child_mesh_index
        this->unzip_interface_internal_modify_child_mesh(iG,tInterfaceNodeInds,tNewUnzippedNodeInds,tNewUnzippedNodeIds);

    }

}
// ----------------------------------------------------------------------------------
void
Model::unzip_interface_internal_assign_node_identifiers(moris::uint aNumNodes,
                                                        moris::Matrix<moris::IdMat> & aUnzippedNodeIndices,
                                                        moris::Matrix<moris::IdMat> & aUnzippedNodeIds)
{
    // Verify sizes
    MORIS_ASSERT(aUnzippedNodeIndices.numel() == aNumNodes, "Size mismatch between aNumNodes and aUnzippedNodeIndices. Please pre-allocate these matrices ");
    MORIS_ASSERT(aUnzippedNodeIds.numel() == aNumNodes, "Size mismatch between aNumNodes and aUnzippedNodeIds.  Please pre-allocate these matrices ");

    // Ask the mesh for new node ids
    moris::moris_index tNodeIndexOffset = mBackgroundMesh.get_first_available_index(EntityRank::NODE);
    moris::moris_id    tNodeIdOffset    = mBackgroundMesh.allocate_entity_ids(aNumNodes,EntityRank::NODE);

    // Iterate new nodes and assign new node ids
    for( uint iN = 0; iN<aNumNodes; iN++)
    {

        // TODO: ADD PARALLEL OWNERSHIP STUFF HERE TO ASSIGN CONSISTENT NODE IDS
        // Give node global ids
        aUnzippedNodeIds(iN) = tNodeIdOffset;

        // increase the id offset
        tNodeIdOffset++;

        // Give nodes processor indices
        aUnzippedNodeIndices(iN) = tNodeIndexOffset;

        // increase the node index offset
        tNodeIndexOffset++;
    }

    // update the first available node index in our background mesh
    mBackgroundMesh.update_first_available_index(tNodeIndexOffset,EntityRank::NODE);

}
// ----------------------------------------------------------------------------------
void
Model::unzip_interface_internal_modify_child_mesh(moris::uint                         aGeometryIndex,
                                                  moris::Matrix<moris::IdMat> const & aInterfaceNodeIndices,
                                                  moris::Matrix<moris::IdMat> const & aUnzippedNodeIndices,
                                                  moris::Matrix<moris::IdMat> const & aUnzippedNodeIds)
{

    // from interface node indices, figure out which interface nodes live in which interface
    moris::Cell<moris::Cell< moris::moris_index >> tChildMeshInterfaceNodes = unzip_interface_internal_collect_child_mesh_to_interface_node(aInterfaceNodeIndices,aUnzippedNodeIndices,aUnzippedNodeIds);

    // Flag indicating there is an interface without an element pair
    bool tNoPairFlag = false;

    // Iterate through child meshes and add new unzipped nodes
    uint tCMIndex = 0;
    for(auto iCM = tChildMeshInterfaceNodes.begin(); iCM != tChildMeshInterfaceNodes.end(); ++iCM)
    {
        // number of interface nodes in this child mesh
        uint tNumCMInterfaceNodes = iCM->size();

        // Allocate matrices of interface node indices
        moris::Matrix< moris::IndexMat > tCMInterfaceNodeIndices(1,tNumCMInterfaceNodes);
        moris::Matrix< moris::IndexMat > tCMUnzippedInterfaceNodeIndices(1,tNumCMInterfaceNodes);
        moris::Matrix< moris::IdMat >    tCMUnzippedInterfaceNodeIds(1,tNumCMInterfaceNodes);

        // Collect information on the interface nodes on this child mesh
        for(moris::uint iN = 0; iN<tNumCMInterfaceNodes; iN++)
        {
            // node index local to the numbering scheme in interface nodes
            moris::moris_index tInterfaceLocInd = (*iCM)(iN);

            tCMInterfaceNodeIndices(iN)         = aInterfaceNodeIndices(tInterfaceLocInd);
            tCMUnzippedInterfaceNodeIndices(iN) = aUnzippedNodeIndices(tInterfaceLocInd);
            tCMUnzippedInterfaceNodeIds(iN)     = aUnzippedNodeIds(tInterfaceLocInd);
        }

        // Tell the child mesh to unzip it's interface
        Child_Mesh & tChildMesh = mCutMesh.get_child_mesh(tCMIndex);

        // initialize the unzipping (which basically allocated node to element connectivity
        // in the child mesh because it is not typically needed)
        tChildMesh.initialize_unzipping();

        // Ask the child mesh to construct interface element pairs
        moris::Matrix<moris::IndexMat> tInterfaceElementPairsCMIndex;
        moris::Matrix<moris::IndexMat> tInterfaceSideOrdinals;

        tChildMesh.unzip_child_mesh_interface_get_interface_element_pairs(aGeometryIndex,tNoPairFlag,tInterfaceElementPairsCMIndex,tInterfaceSideOrdinals);

        // Convert the pairs to processor local indices because we need to be able to access the element phase index
        moris::Matrix<moris::IndexMat> tInterfaceElementPairs = tChildMesh.convert_to_proc_local_elem_inds(tInterfaceElementPairsCMIndex);

        // TODO: Add method to resolve cross child mesh element pairs for when the interface coincides with a parent face
        // NOTE: By using the sign of the geometry value, it really shouldnt take a whole lot of work to accomodated
        MORIS_ERROR(!tNoPairFlag," in unzip_interface_internal_modify_child_mesh, interface detected on a child mesh boundary. Currently, no method is implemented to resolve this");

        // Take the child mesh pairs and determine who gets which id
        // This output is either a 0 or 1, meaning the first or second element of the pair gets the unzipped nodes
        moris::Matrix< moris::IndexMat > tElementWhichKeepsOriginalNodes =
                this->unzip_interface_internal_assign_which_element_uses_unzipped_nodes(aGeometryIndex,tInterfaceElementPairs);

        // Get the elements on the boundary
        tChildMesh.unzip_child_mesh_interface(aGeometryIndex,
                                              tInterfaceElementPairsCMIndex,
                                              tElementWhichKeepsOriginalNodes,
                                              tCMInterfaceNodeIndices,
                                              tCMUnzippedInterfaceNodeIndices,
                                              tCMUnzippedInterfaceNodeIds);

        // Construct interface elements
        unzip_interface_construct_interface_elements(aGeometryIndex,
                                                     tInterfaceElementPairs,
                                                     tInterfaceSideOrdinals);


        tChildMesh.finalize_unzipping();

        tCMIndex++;
    }

    unzip_interface_assign_element_identifiers();


}
// ----------------------------------------------------------------------------------
moris::Matrix< moris::IndexMat >
Model::unzip_interface_internal_assign_which_element_uses_unzipped_nodes( moris::moris_index aGeometryIndex,
                                                                          moris::Matrix< moris::IndexMat > const & aInterfaceElementPairs )
{

    // specify which geometry sign gets to keep the nodes
    moris::moris_index tValWhichUsesUnzipped = 0;

    // The rule used here is whichever phase gets a 1 from the phase table with respect to the current geometry index
    // gets to keep the original node. The other element changes it's nodes to the unzipped indices.
    moris::Matrix< moris::IndexMat > tElementWhichKeepsUsesUnzippedNodes(aInterfaceElementPairs.n_cols());

    // number of pairs
    moris::uint tNumPairs = aInterfaceElementPairs.n_cols();

    // allocate
    moris::moris_index tElement0           = MORIS_INDEX_MAX;
    moris::moris_index tElement1           = MORIS_INDEX_MAX;
    moris::moris_index tElement0PhaseIndex = MORIS_INDEX_MAX;
    moris::moris_index tElement1PhaseIndex = MORIS_INDEX_MAX;

    // iterate through pairs
    for(moris::uint iP = 0; iP<tNumPairs; iP++)
    {
        // set this back to moris index max so we dont run into issues if there is actually only one element in the pair
        moris::moris_index tElement0GeomSign   = MORIS_INDEX_MAX; // sign of phase of element wrt to aGeometryIndex (0 for negative, 1 for positive)
        moris::moris_index tElement1GeomSign   = MORIS_INDEX_MAX;// sign of phase of element wrt to aGeometryIndex (0 for negative, 1 for positive)

        // Element indices
        tElement0 = aInterfaceElementPairs(0,iP);
        tElement1 = aInterfaceElementPairs(1,iP);

        // Figure out which element in the pair gets to keep the original
        if(tElement0 != MORIS_INDEX_MAX)
        {
            tElement0PhaseIndex = mBackgroundMesh.get_element_phase_index(tElement0);
            tElement0GeomSign = mGeometryEngine->get_phase_sign_of_given_phase_and_geometry(tElement0PhaseIndex,aGeometryIndex);

            if(tElement0GeomSign == tValWhichUsesUnzipped)
            {
                tElementWhichKeepsUsesUnzippedNodes(iP) = 0;
            }
        }

        if(tElement1 != MORIS_INDEX_MAX)
        {
            tElement1PhaseIndex = mBackgroundMesh.get_element_phase_index(tElement1);
            tElement1GeomSign = mGeometryEngine->get_phase_sign_of_given_phase_and_geometry(tElement1PhaseIndex,aGeometryIndex);
            if(tElement1GeomSign == tValWhichUsesUnzipped)
            {
                tElementWhichKeepsUsesUnzippedNodes(iP) = 1;
            }
        }

        // Make sure both don't end up with the same sign
        MORIS_ASSERT(tElement1GeomSign != tElement0GeomSign,"Both elements in an interface pair returned the same phase sign");

        MORIS_ASSERT(tElement0GeomSign != MORIS_INDEX_MAX,"tElement0GeomSign no pair");
        MORIS_ASSERT(tElement1GeomSign != MORIS_INDEX_MAX,"tElement1GeomSign no pair");

    }

    return tElementWhichKeepsUsesUnzippedNodes;

}


// ----------------------------------------------------------------------------------
moris::Cell<moris::Cell< moris::moris_index >>
Model::unzip_interface_internal_collect_child_mesh_to_interface_node(moris::Matrix<moris::IdMat> const & aInterfaceNodeIndices,
                                                                     moris::Matrix<moris::IdMat> const & aUnzippedNodeIndices,
                                                                     moris::Matrix<moris::IdMat> const & aUnzippedNodeIds)
{
    // Allocate cell to keep track of the node indices of each child mesh
    moris::uint tNumChildMeshes = mCutMesh.get_num_child_meshes();

    moris::Cell<moris::Cell< moris::moris_index >> tChildMeshInterfaceNodes(tNumChildMeshes);

    // Iterate through interface node indices
    for(moris::uint iN = 0; iN < aInterfaceNodeIndices.numel(); iN++)
    {
        // Get the child meshes which have the interface node
        moris::Matrix< moris::IndexMat > tNodeChildMeshIndices = mBackgroundMesh.get_node_child_mesh_assocation(aInterfaceNodeIndices(iN));

        // Iterate through child meshes and mark interface node indices in this child mesh
        for(moris::uint iCM = 0; iCM<tNodeChildMeshIndices.numel(); iCM++)
        {
            moris::moris_index tCMIndex = tNodeChildMeshIndices(iCM);
            tChildMeshInterfaceNodes(tCMIndex).push_back(iN);
        }
    }

    return tChildMeshInterfaceNodes;
}
// ----------------------------------------------------------------------------------
void
Model::unzip_interface_construct_interface_elements(moris::uint aGeometryIndex,
                                                    moris::Matrix< moris::IndexMat > const & aElementPairs,
                                                    moris::Matrix< moris::IndexMat > const & aSideOrdinalPairs)
{
    moris::uint tNumPairs = aElementPairs.n_cols();

    moris::Cell<const moris::mtk::Cell*> tPairCells(2);

    for(moris::uint i = 0; i <tNumPairs; i++)
    {
        tPairCells(0) = mBackgroundMesh.get_child_element_mtk_cell_ptr(aElementPairs(0,i));
        tPairCells(1) = mBackgroundMesh.get_child_element_mtk_cell_ptr(aElementPairs(1,i));

        // give ownership information to left element
        moris_index tOwnerOfElem0 = tPairCells(0)->get_owner();

        // construct an interface element
        Interface_Element tInterfaceElement;
        tInterfaceElement.set_element_pair_and_side_ordinal(tPairCells,aSideOrdinalPairs.get_column(i));
        tInterfaceElement.set_element_owner(tOwnerOfElem0);

        // Add interface element to cut mesh
        mCutMesh.add_interface_element(tInterfaceElement);
    }
}
// ----------------------------------------------------------------------------------
void
Model::unzip_interface_assign_element_identifiers()
{
    moris::Cell<Interface_Element> & tInterfaceElements = mCutMesh.get_interface_elements();
    moris::uint tNumInterfaceElements = tInterfaceElements.size();

    // Allocate ids
    moris::moris_index tIdOffset    = mBackgroundMesh.allocate_entity_ids(tNumInterfaceElements,EntityRank::ELEMENT);
    moris::moris_index tIndexOffset = mBackgroundMesh.get_first_available_index(EntityRank::ELEMENT);

    for(moris::uint i = 0; i <tNumInterfaceElements; i++)
    {
        tInterfaceElements(i).set_element_identifiers(tIndexOffset,tIdOffset);
        tIndexOffset++;
        tIdOffset++;
    }

    mBackgroundMesh.update_first_available_index(tIndexOffset,EntityRank::ELEMENT);

}



// ----------------------------------------------------------------------------------
// Enrichment Source code
// ----------------------------------------------------------------------------------
void
Model::perform_basis_enrichment(enum EntityRank  const & aBasisRank,
                                moris_index      const & aMeshIndex)
{
    // Start the clock
    std::clock_t start = std::clock();

    MORIS_ERROR(mDecomposed,"Prior to computing basis enrichment, the decomposition process must be called");

    // allocate some new enriched interpolation and integration meshes
    mEnrichedIntegMesh.resize(aMeshIndex+1,nullptr);
    mEnrichedInterpMesh.resize(aMeshIndex+1,nullptr);

    this->perform_basis_enrichment_internal(aBasisRank,{{aMeshIndex}});

    // Change the enrichment flag
    mEnriched = true;

    if(moris::par_rank() == 0 && mVerbose)
    {
        std::cout<<"XTK: Basis enrichment computation completed in " <<(std::clock() - start) / (double)(CLOCKS_PER_SEC)<<" s."<<std::endl;
        std::cout<<"XTK: Basis enrichment performed on mesh index: "<< aMeshIndex<<std::endl;
    }
}
// ----------------------------------------------------------------------------------
void
Model::perform_basis_enrichment(enum EntityRank  const & aBasisRank,
                                Matrix<IndexMat> const & aMeshIndex)
{
    // Start the clock
    std::clock_t start = std::clock();

    MORIS_ERROR(mDecomposed,"Prior to computing basis enrichment, the decomposition process must be called");

    // allocate some new enriched interpolation and integration meshes
    mEnrichedIntegMesh.resize(aMeshIndex.numel()+1,nullptr);
    mEnrichedInterpMesh.resize(aMeshIndex.numel()+1,nullptr);

    this->perform_basis_enrichment_internal(aBasisRank,aMeshIndex);

    // Change the enrichment flag
    mEnriched = true;

    if(moris::par_rank() == 0 && mVerbose)
    {
        std::cout<<"XTK: Basis enrichment computation completed in " <<(std::clock() - start) / (double)(CLOCKS_PER_SEC)<<" s."<<std::endl;
        std::cout<<"XTK: Basis enrichment performed on meshes:";
        for(moris::uint i = 0; i < aMeshIndex.numel(); i++)
        {
            std::cout<<std::setw(6)<<aMeshIndex(i);
        }
        std::cout<<std::endl;
    }
}
// ----------------------------------------------------------------------------------

Enrichment const &
Model::get_basis_enrichment()
{
    MORIS_ASSERT(mEnriched,"Cannot get basis enrichment from an XTK model which has not called perform_basis_enrichment ");
    return *mEnrichment;
}
Enriched_Interpolation_Mesh &
Model::get_enriched_interp_mesh(moris::moris_index aIndex)
{
    MORIS_ASSERT(mEnriched,"Cannot get enriched interpolation mesh from an XTK model which has not called perform_basis_enrichment ");
    return *(mEnrichedInterpMesh(aIndex));
}
Enriched_Integration_Mesh &
Model::get_enriched_integ_mesh(moris::moris_index aIndex)
{
    MORIS_ASSERT(mEnriched,"Cannot get enriched integration mesh from an XTK model which has not called perform_basis_enrichment ");
    return *(mEnrichedIntegMesh(aIndex));
}


void
Model::perform_basis_enrichment_internal(enum EntityRank  const & aBasisRank,
                                         Matrix<IndexMat> const & aMeshIndex)
{
    // initialize enrichment (ptr because of circular dependency)
    mEnrichment = new Enrichment(Enrichment_Method::USE_INTERPOLATION_CELL_BASIS,
    							 aBasisRank,
    							aMeshIndex,
    							mGeometryEngine->get_num_phases(),
    							this,
    							&mCutMesh,
    							&mBackgroundMesh);

    // Set verbose flag to match XTK.
    mEnrichment->mVerbose = mVerbose;

    // perform the enrichment
    mEnrichment->perform_enrichment();

}


void
Model::construct_face_oriented_ghost_penalization_cells()
{
    MORIS_ERROR(mDecomposed,"Mesh needs to be decomposed prior to calling ghost penalization");
    MORIS_ERROR(!mGhost,"Ghost penalization has already been called");

    std::clock_t start = std::clock();

    mGhostStabilization = new Ghost_Stabilization(this);

    mGhostStabilization->setup_ghost_stabilization();

    mGhost = true;

    if(moris::par_rank() == 0  && mVerbose)
    {
        std::cout<<"XTK: Ghost stabilization setup completed in "<< (std::clock() - start) / (double)(CLOCKS_PER_SEC)<<" s."<<std::endl;
    }
}

Ghost_Stabilization &
Model::get_ghost_stabilization(moris::moris_index  aIndex)
{
    MORIS_ERROR(mGhost,"Ghost has not been constructed on this model.");
    return *mGhostStabilization;
}

void Model::construct_multigrid()
{
    mMultigrid = std::make_shared< xtk::Multigrid >( this );

    mMultigrid->build_enriched_coeff_to_background_coeff_map();

    mMultigrid->create_fine_to_coarse_relationship();

    mMultigrid->create_coarse_to_fine_relationship();

    mMultigrid->create_coarse_to_fine_weights();

    std::string tName = "Enriched_bspline_1.exo";
    mMultigrid->build_basis_exodus_information(tName);
}

// ----------------------------------------------------------------------------------
// Tet 10 conversion Source code
// ----------------------------------------------------------------------------------
void
Model::convert_mesh_tet4_to_tet10()
{
    MORIS_ASSERT(0,"not currently working");
    mConvertedToTet10s = true;

    // start timing on this decomposition
    std::clock_t start = std::clock();
    //        convert_mesh_tet4_to_tet10_internal();

    if(moris::par_rank() == 0)
    {
        std::cout<<"Tet4 to Tet10 conversion completed in "<< (std::clock() - start) / (double)(CLOCKS_PER_SEC)<<" s."<<std::endl;
    }
}

// ----------------------------------------------------------------------------------
// Export mesh Source code
// ----------------------------------------------------------------------------------
void
Model::extract_surface_mesh_to_obj(std::string                      aOutputFile,
                                   size_t                           aPhaseIndex,
                                   moris::Cell<std::string> const & aBoundingSideSets)
{
    // start timing on this decomposition
    std::clock_t start = std::clock();

    // create the output mtk mesh
    extract_surface_mesh_to_obj_internal(aOutputFile,aPhaseIndex,aBoundingSideSets);

    if(moris::par_rank() == 0  && mVerbose)
    {
        std::cout<<"XTK: Extract surface to obj completed in "<< (std::clock() - start) / (double)(CLOCKS_PER_SEC)<<" s."<<std::endl;
        std::cout<<"XTK: OBJ File: "<<aOutputFile<<std::endl;
    }
}
//------------------------------------------------------------------------------

void
Model::construct_neighborhood()
{

    mElementToElement.resize(this->get_num_elements_total());

    // add uncut neighborhood to connectivity
    // keep track of boundaries with an uncut mesh next to a cut mesh
    moris::Cell<moris::Cell<moris_index>> tCutToUncutFace;
    construct_uncut_neighborhood(tCutToUncutFace);

    // since there are hanging nodes in HMR we need to do a little extra
    // to get the neighborhood correct
    if(mBackgroundMesh.get_mesh_data().get_mesh_type() == MeshType::HMR)
    {
        construct_cut_mesh_simple_neighborhood();

        construct_cut_mesh_to_uncut_mesh_neighborhood(tCutToUncutFace);

        construct_complex_neighborhood();
    }
    else
    {
        // create the simple relationship neighborhood between child meshes
        construct_cut_mesh_simple_neighborhood();

        // create the link between uncut background cells and their neighboring children cells
        construct_cut_mesh_to_uncut_mesh_neighborhood(tCutToUncutFace);
    }

    //    print_neighborhood();
}


void
Model::construct_subphase_neighborhood()
{
    // get the interpolation mesh
    moris::mtk::Interpolation_Mesh & tInterpMesh = mBackgroundMesh.get_mesh_data();

    // allocate subphase to subphase connectivity
    mSubphaseToSubPhase = moris::Cell<moris::Cell<moris::moris_index>>(mCutMesh.get_num_subphases());
    mSubphaseToSubPhaseMySideOrds = moris::Cell<moris::Cell<moris::moris_index>>(mCutMesh.get_num_subphases());
    mSubphaseToSubPhaseNeighborSideOrds = moris::Cell<moris::Cell<moris::moris_index>>(mCutMesh.get_num_subphases());
    mTransitionNeighborCellLocation = moris::Cell<moris::Cell<moris::moris_index>>(mCutMesh.get_num_subphases());

    // non unique temporary data
    moris::Cell<moris::Cell<moris::moris_index>> tNonUniqueSubphaseToSubphase(mCutMesh.get_num_subphases());
    moris::Cell<moris::Cell<moris::moris_index>> tNonUniqueSubphaseToSubPhaseMySideOrds(mCutMesh.get_num_subphases());
    moris::Cell<moris::Cell<moris::moris_index>> tNonUniqueSubphaseToSubPhaseNeighborSideOrds(mCutMesh.get_num_subphases());
    moris::Cell<moris::Cell<moris::moris_index>> tNonUniqueTransitionLocation(mCutMesh.get_num_subphases());

    //iterate through facets
    moris::uint tNumFacets = tInterpMesh.get_num_entities(EntityRank::ELEMENT);
    for(moris::moris_index iC = 0; iC < (moris::moris_index)tNumFacets; iC++)
    {
        // current cell
        mtk::Cell const * tCurrentCell = & tInterpMesh.get_mtk_cell(iC);

        // get the cells attached to the facet
        Matrix<IndexMat> tCellToCellSideIndex = tInterpMesh.get_elements_connected_to_element_and_face_ind_loc_inds(iC);
        Matrix<IndexMat> tCellToCellSideOrd  = tInterpMesh.get_elements_connected_to_element_and_face_ord_loc_inds(iC);

        // get the neighboring cells
        Cell<mtk::Cell const *> tCells(tCellToCellSideOrd.numel());
        tInterpMesh.get_mtk_cells(tCellToCellSideOrd.get_row(0),tCells);

        // iterate through neighbor
        for(moris::uint iN = 0; iN < tCellToCellSideOrd.n_cols(); iN++)
        {
            // facet ordinal shared for current neighbors
            moris_index tMyOrdinal  = tCellToCellSideOrd(1,iN);
            moris_index tNeighborOrdinal = tCellToCellSideOrd(2,iN);
            moris_index tTransitionCellLocation = tCellToCellSideOrd(3,iN);

            // neighbor cell
            mtk::Cell const * tOtherCell =  & tInterpMesh.get_mtk_cell(tCellToCellSideIndex(0,iN));

            // get the subphase indices attached to the facet which is connected to the current cell
            Cell<moris::moris_index> tMyCellSubphaseIndices(0);
            Cell<moris::moris_index> tMyCellSubphaseBulkIndices(0);
            this->collect_subphases_attached_to_facet_on_cell( tCurrentCell->get_index(), tMyOrdinal, tMyCellSubphaseIndices, tMyCellSubphaseBulkIndices);

            // get the subphase indices attached to the facet which is connected to the other cell in the neighborhood
            Cell<moris::moris_index> tNeighborSubphaseIndices(0);
            Cell<moris::moris_index> tNeighborSubphaseBulkIndices(0);
            this->collect_subphases_attached_to_facet_on_cell( tOtherCell->get_index(), tNeighborOrdinal, tNeighborSubphaseIndices, tNeighborSubphaseBulkIndices);

            // iterate over subphases and add to neighborhood
            for(moris::uint i = 0; i < tMyCellSubphaseIndices.size(); i++)
            {
                moris_index tMyBulkIndex = tMyCellSubphaseBulkIndices(i);
                moris_index tMySubphaseIndex = tMyCellSubphaseIndices(i);

                for(moris::uint j = 0; j < tNeighborSubphaseIndices.size(); j++)
                {
                    moris_index tNeighborBulkIndex = tNeighborSubphaseBulkIndices(j);
                    moris_index tNeighborSubphaseIndex = tNeighborSubphaseIndices(j);

                    if(tMyBulkIndex == tNeighborBulkIndex)
                    {
                        mSubphaseToSubPhase(tMySubphaseIndex).push_back(tNeighborSubphaseIndex);
                        mSubphaseToSubPhaseMySideOrds(tMySubphaseIndex).push_back(tMyOrdinal);
                        mSubphaseToSubPhaseNeighborSideOrds(tMySubphaseIndex).push_back(tNeighborOrdinal);
                        mTransitionNeighborCellLocation(tMySubphaseIndex).push_back(tTransitionCellLocation);

//                        tNonUniqueSubphaseToSubphase(tNeighborSubphaseIndex).push_back(tMySubphaseIndex);
//                        tNonUniqueSubphaseToSubPhaseMySideOrds(tNeighborSubphaseIndex).push_back(tNeighborOrdinal);
//                        tNonUniqueSubphaseToSubPhaseNeighborSideOrds(tNeighborSubphaseIndex).push_back(tMyOrdinal);
//                        tNonUniqueTransitionLocation(tNeighborSubphaseIndex).push_back(MORIS_INDEX_MAX);
                    }
                }
            }

        }
    }

    //    // make unique (We do this this way because HMR only has neighbors from larger cells to lower cells so
    //    // a criteria to check cell ids does not work)
    //    for(moris::uint i = 0; i < mSubphaseToSubPhase.size(); i++)
    //    {
    //        Cell<moris::moris_index> tUniqueIndices = unique_index( tNonUniqueSubphaseToSubphase(i) );
    //
    //
    //        std::cout<<" i = "<< i<<std::endl;
    //
    //        // resize member data
    //        mSubphaseToSubPhase(i).resize(tUniqueIndices.size());
    //        mSubphaseToSubPhaseMySideOrds(i).resize(tUniqueIndices.size());
    //        mSubphaseToSubPhaseNeighborSideOrds(i).resize(tUniqueIndices.size());
    //        mTransitionNeighborCellLocation(i).resize(tUniqueIndices.size());
    //
    //        moris::print(tNonUniqueTransitionLocation(i),"tNonUniqueTransitionLocation");
    //        moris::print(tNonUniqueSubphaseToSubphase(i),"tNonUniqueSubphaseToSubphase");
    //        moris::print(tUniqueIndices,"tUniqueIndices");
    //        // commit unique entries to member data
    //        for(moris::uint j = 0; j < tUniqueIndices.size(); j ++)
    //        {
    //            mSubphaseToSubPhase(i)(j) = tNonUniqueSubphaseToSubphase(i)(tUniqueIndices(j));
    //            mSubphaseToSubPhaseMySideOrds(i)(j) = tNonUniqueSubphaseToSubPhaseMySideOrds(i)(tUniqueIndices(j));
    //            mSubphaseToSubPhaseNeighborSideOrds(i)(j) = tNonUniqueSubphaseToSubPhaseNeighborSideOrds(i)(tUniqueIndices(j));
    //            mTransitionNeighborCellLocation(i)(j) = tNonUniqueTransitionLocation(i)(tUniqueIndices(j));
    //        }
    //
    //        moris::print(mTransitionNeighborCellLocation(i),"mTransitionNeighborCellLocation(i)");
    //
    //   }
}

void
Model::collect_subphases_attached_to_facet_on_cell(moris::moris_index aCellIndex,
                                                   moris::moris_index aFacetOrdinal,
                                                   Cell<moris::moris_index> & aCellSubphaseIndices,
                                                   Cell<moris::moris_index> & aCellSubphaseBulkIndices)
{
    // set pointer to child mesh if it has children
    if(mBackgroundMesh.entity_has_children(aCellIndex, EntityRank::ELEMENT))
    {
        // get the child mesh ptr
        moris::moris_index tCMIndex = mBackgroundMesh.child_mesh_index(aCellIndex,EntityRank::ELEMENT);
        Child_Mesh const * tCMCell = & mCutMesh.get_child_mesh(tCMIndex);

        Matrix<IndexMat> tCellFacets = mBackgroundMesh.get_mesh_data().get_entity_connected_to_entity_loc_inds(aCellIndex,EntityRank::ELEMENT,mBackgroundMesh.get_mesh_data().get_facet_rank());

        moris_index tFacetIndex = tCellFacets(aFacetOrdinal);

        // get subphases attached to facet
        Cell<moris::moris_index> tCellCMSubphaseIndices;
        tCMCell->get_subphases_attached_to_facet(tFacetIndex, tCellCMSubphaseIndices);

        // reference to subphase indices and bulkphases of subphases
        Cell<moris::moris_index> const & tCMSubphaseBulkIndices = tCMCell->get_subphase_bin_bulk_phase();
        Cell<moris::moris_index> const & tCMSubphaseIndices     = tCMCell->get_subphase_indices();

        // put the information in processor local indices with bulk phase
        aCellSubphaseIndices.resize(tCellCMSubphaseIndices.size());
        aCellSubphaseBulkIndices.resize(tCellCMSubphaseIndices.size());

        for(moris::uint i = 0; i < tCellCMSubphaseIndices.size(); i++)
        {
            aCellSubphaseBulkIndices(i) = tCMSubphaseBulkIndices(tCellCMSubphaseIndices(i));
            aCellSubphaseIndices(i) = tCMSubphaseIndices(tCellCMSubphaseIndices(i));
        }

    }
    else
    {
        // get the cell subphase indices
        aCellSubphaseBulkIndices = {{mBackgroundMesh.get_element_phase_index(aCellIndex)}};
        aCellSubphaseIndices     = {{aCellIndex}};
    }

}

bool
Model::subphase_is_in_child_mesh(moris_index aSubphaseIndex)
{
    if(aSubphaseIndex >= (moris_index)mBackgroundMesh.get_mesh_data().get_num_entities(EntityRank::ELEMENT))
    {
        return true;
    }

    else if (mBackgroundMesh.entity_has_children(aSubphaseIndex,EntityRank::ELEMENT))
    {
        return true;
    }

    return false;
}

void
Model::construct_cut_mesh_simple_neighborhood()
{
    // Collect element to node of child mesh and create a full element to element graph for child mesh
    Matrix<IndexMat> tCMElementToNode = mCutMesh.get_full_element_to_node_loc_inds();

    // generate connectivities (face then  element conn) we do not keep faces around though
    // face connectivity
    CellTopology tCellTopo = mCutMesh.get_child_element_topology();
    Matrix<IndexMat> tElementToFace;
    Matrix<IndexMat> tFaceToNode;
    Matrix<IndexMat> tNodeToFace;
    Matrix<IndexMat> tFaceToElement;

    moris::mtk::Cell_Info_Factory tFactory;
    moris::mtk::Cell_Info* tCellInfo = tFactory.create_cell_info(tCellTopo);
    xtk::create_faces_from_element_to_node(tCellInfo, mBackgroundMesh.get_num_entities(EntityRank::NODE), tCMElementToNode, tElementToFace, tFaceToNode, tNodeToFace, tFaceToElement);

    // element connectivity
    moris::size_t tNumFacePerElem = tCellInfo->get_num_facets();

    // generate the element to element connectivity. this only captures the easy neighborhood relationships
    // where we consider elements neighbors if they share a full face
    Matrix<IndexMat> tElementToElementMat = generate_element_to_element(tFaceToElement, mCutMesh.get_num_entities(EntityRank::ELEMENT), tNumFacePerElem, MORIS_INDEX_MAX);

    moris::Matrix< moris::IndexMat > tCMElementInds = mCutMesh.get_all_element_inds();

    moris::Cell<moris::mtk::Cell*> tCMCellPtrs(tCMElementInds.numel());
    for(moris::uint iC = 0; iC<tCMCellPtrs.size(); iC++ )
    {
        tCMCellPtrs(iC) = &mBackgroundMesh.get_mtk_cell(tCMElementInds(iC));
    }

    // add to member data
    for(moris::uint iC = 0; iC<tElementToElementMat.n_rows(); iC++ )
    {
        for(moris::uint iN = 0; iN< tElementToElementMat.n_cols(); iN++)
        {
            if(tElementToElementMat(iC,iN) == MORIS_INDEX_MAX)
            {
                continue;
            }
            mElementToElement(tCMElementInds(iC)).push_back(tCMCellPtrs(tElementToElementMat(iC,iN)));
        }
    }

    delete tCellInfo;
}

void
Model::construct_complex_neighborhood()
{

}


void
Model::construct_cut_mesh_to_uncut_mesh_neighborhood(moris::Cell<moris::Cell<moris_index>> & aCutToUncutFace)
{
    // matrices used throughout routine
    moris::Matrix< moris::IdMat >    tChildElemsIdsOnFace;
    moris::Matrix< moris::IndexMat > tChildElemsCMIndOnFace;
    moris::Matrix< moris::IndexMat > tChildElemOnFaceOrdinal;

    // iterate through the cut to uncut relationships
    for(moris::uint iR = 0; iR < aCutToUncutFace.size(); iR++)
    {
        // get data for readability from aCutToUncutFace
        moris_index tCellUnCutInd = aCutToUncutFace(iR)(0);
        moris_index tCellCutInd   = aCutToUncutFace(iR)(1);
        moris_index tFaceIndex    = aCutToUncutFace(iR)(2);

        // uncut cell
        moris::mtk::Cell* tCellUnCut = &mBackgroundMesh.get_mtk_cell(tCellUnCutInd);

        Child_Mesh &  tChildMesh = mCutMesh.get_child_mesh(mBackgroundMesh.child_mesh_index(tCellCutInd,EntityRank::ELEMENT));
        Matrix<IndexMat> const & tChildCellInds = tChildMesh.get_element_inds();
        tChildMesh.get_child_elements_connected_to_parent_facet(tFaceIndex, tChildElemsIdsOnFace, tChildElemsCMIndOnFace, tChildElemOnFaceOrdinal);

        // get child cells and add cut cell to neighborhood, cut cells to neighborhood of uncut
        for(moris::uint  i = 0; i < tChildElemsCMIndOnFace.numel(); i++ )
        {
            mElementToElement(tChildCellInds(tChildElemsCMIndOnFace(i))).push_back(tCellUnCut);
            mElementToElement(tCellUnCutInd).push_back(&mBackgroundMesh.get_mtk_cell(tChildCellInds(tChildElemsCMIndOnFace(i))));
        }
    }
}

void
Model::construct_uncut_neighborhood(moris::Cell<moris::Cell<moris_index>> & aCutToUncutFace)
{
    moris::mtk::Interpolation_Mesh & tInterpMesh = mBackgroundMesh.get_mesh_data();
    moris::uint tNumCells                        = tInterpMesh.get_num_entities(EntityRank::ELEMENT);

    // iterate through background cells
    for(moris::uint iC = 0; iC < tNumCells; iC++)
    {
        // if i have no children
        if(!mBackgroundMesh.entity_has_children((moris_index)iC,EntityRank::ELEMENT))
        {
            // get the neighbors
            Matrix<IndexMat> tElementNeighors = tInterpMesh.get_elements_connected_to_element_and_face_ind_loc_inds(iC);

            // iterate through neighbors
            for(moris::uint iN = 0; iN<tElementNeighors.n_cols(); iN++ )
            {
                // if the neighbor has no children
                if(!mBackgroundMesh.entity_has_children(tElementNeighors(0,iN),EntityRank::ELEMENT))
                {
                    mElementToElement(iC).push_back(&mBackgroundMesh.get_mtk_cell(tElementNeighors(0,iN)));
                }

                // mark as a cut to uncut boundary
                else
                {
                    aCutToUncutFace.push_back( {(moris_index)iC , tElementNeighors(0,iN), tElementNeighors(1,iN)} );
                }
            }
        }
    }
}

void
Model::print_neighborhood()
{
    for(moris::uint iE = 0; iE < mElementToElement.size(); iE++)
    {
        moris::mtk::Cell & tCell = mBackgroundMesh.get_mtk_cell(iE);

        std::cout<<"Element Id: "<<std::setw(8)<<tCell.get_id()<<" | ";

        for(moris::uint iN = 0; iN < mElementToElement(iE).size(); iN++)
        {
            std::cout<<std::setw(8)<<mElementToElement(iE)(iN)->get_id();
        }
        std::cout<<std::endl;
    }
}


void
Model::print_cells()
{
    for(moris::uint  i = 0; i < this->get_num_elements_total(); i++ )
    {
        mtk::Cell & tCell = mBackgroundMesh.get_mtk_cell((moris_index)i);

        moris::Cell<moris::mtk::Vertex *> tVertices = tCell.get_vertex_pointers();

        std::cout<<"Cell Id: "<<std::setw(8)<<tCell.get_id();
        std::cout<<" | Cell Index: "<<std::setw(8)<<tCell.get_index();
        std::cout<<" | Phase: "<<std::setw(8)<<mBackgroundMesh.get_element_phase_index(i) << " | Verts: ";
        for(moris::uint j = 0; j < tVertices.size(); j++)
        {
            std::cout<<std::setw(8)<<tVertices(j)->get_id();
        }
        std::cout<<std::endl;
    }
}

void
Model::print_vertex_geometry()
{
        for(moris::uint  i = 0; i < mBackgroundMesh.get_num_entities(EntityRank::NODE); i++)
    {
        moris::mtk::Vertex & tVertex = mBackgroundMesh.get_mtk_vertex(i);

        std::cout<<"Vertex Id: "<<std::setw(8)<<tVertex.get_id()<<" | ";
        for(moris::uint j = 0; j < mGeometryEngine->get_num_geometries(); j++)
        {
            std::cout<<std::setw(12)<<mGeometryEngine->get_entity_phase_val(tVertex.get_index(),j)<<" , ";
        }
        std::cout<<std::endl;
    }
}


void
Model::print_interface_vertices()
{
    mBackgroundMesh.print_interface_node_flags();
}
void
Model::print_subphase_neighborhood()
{

    std::cout<<"Subphases"<<std::endl;
    for(moris::uint iC = 0; iC<mSubphaseToSubPhase.size(); iC++ )
    {
        std::cout<<std::setw(6)<<iC<<" | ";


        for(moris::uint iN = 0; iN< mSubphaseToSubPhase(iC).size(); iN++)
        {
            std::cout<<std::setw(6)<<mSubphaseToSubPhase(iC)(iN);
        }
        std::cout<<std::endl;
    }

    std::cout<<"Subphases My Side Ordinals"<<std::endl;
    for(moris::uint iC = 0; iC<mSubphaseToSubPhaseMySideOrds.size(); iC++ )
    {
        std::cout<<std::setw(6)<<iC<<" | ";


        for(moris::uint iN = 0; iN< mSubphaseToSubPhaseMySideOrds(iC).size(); iN++)
        {
            std::cout<<std::setw(6)<<mSubphaseToSubPhaseMySideOrds(iC)(iN);
        }
        std::cout<<std::endl;
    }


    std::cout<<"Subphases Neighbor Side Ordinals"<<std::endl;
    for(moris::uint iC = 0; iC<mSubphaseToSubPhaseNeighborSideOrds.size(); iC++ )
    {
        std::cout<<std::setw(6)<<iC<<" | ";


        for(moris::uint iN = 0; iN< mSubphaseToSubPhaseNeighborSideOrds(iC).size(); iN++)
        {
            std::cout<<std::setw(6)<<mSubphaseToSubPhaseNeighborSideOrds(iC)(iN);
        }
        std::cout<<std::endl;
    }

    if(mBackgroundMesh.get_mesh_data().get_mesh_type() == MeshType::HMR)
    {
        std::cout<<"Transition Neighbor Locations"<<std::endl;
        for(moris::uint iC = 0; iC<mTransitionNeighborCellLocation.size(); iC++ )
        {
            std::cout<<std::setw(6)<<iC<<" | ";


            for(moris::uint iN = 0; iN< mTransitionNeighborCellLocation(iC).size(); iN++)
            {
                std::cout<<std::setw(12)<<mTransitionNeighborCellLocation(iC)(iN);
            }
            std::cout<<std::endl;
        }
    }
}

//------------------------------------------------------------------------------

void
Model::extract_surface_mesh_to_obj_internal(std::string                      aOutputFile,
                                            size_t                           aPhaseIndex,
                                            moris::Cell<std::string> const & aBoundingSideSets)
{
    MORIS_ERROR(aBoundingSideSets.size() == 6," There needs to be 6 side sets which skin the mesh to extract the surface");

    // allocate side set data
    moris::Cell<moris::Matrix<IndexMat>> tElementIndAndSideOrds(2*6);
    moris::Cell<moris::mtk::MtkSideSetInfo> tBackgroundSideSets(2*6);

    // background mesh data
    moris::mtk::Mesh const & tMeshData = mBackgroundMesh.get_mesh_data();

    // setup outputting options
    Output_Options tOutputOptions;

    // Specify there are 2 possible phases
    size_t tNumPhases = mGeometryEngine->get_num_bulk_phase();

    // Say I only want to output phase 1
    Cell<size_t> tPhasesToOutput = {aPhaseIndex};
    tOutputOptions.change_phases_to_output(tNumPhases,tPhasesToOutput);

    // Initialize Sets information structure
    moris::mtk::MtkSetsInfo tMtkMeshSets;

    moris::uint tNumChildCells   = 0;
    moris::uint tNumNoChildCells = 0;
    for(moris::uint i = 0; i <aBoundingSideSets.size(); i++)
    {
        moris::uint tNoChildInd = 2*i;
        moris::uint tChildInd = 2*i+1;

        this->propogate_background_side_set(aBoundingSideSets(i),tNoChildInd,tChildInd,tElementIndAndSideOrds,tBackgroundSideSets,tOutputOptions,true);

        tNumNoChildCells = tNumNoChildCells + tElementIndAndSideOrds(tNoChildInd).n_rows();
        tNumChildCells = tNumChildCells + tElementIndAndSideOrds(tChildInd).n_rows();
    }

    // get the interface information
    Cell<moris::Matrix<moris::IndexMat>> tInterfaceNodes = mBackgroundMesh.get_interface_nodes_glob_ids();

    // interface sides
    moris::Matrix<moris::IdMat> tInterfaceElemIndandSideOrd = mCutMesh.pack_interface_sides(0,aPhaseIndex,1);


    // tri 3s
    moris::Matrix<moris::IdMat> tTri3ElemToNode(tNumChildCells + tInterfaceElemIndandSideOrd.n_rows() + tNumNoChildCells*4,3);

    // keep track of nodes that are in skinned mesh
    moris::uint tNodeCount = 0;
    moris::Matrix<moris::IdMat> tNodesOnBoundary((tNumNoChildCells*5+tNumChildCells*3 + tInterfaceNodes(aPhaseIndex).numel()) ,1);
    tNodesOnBoundary.fill(MORIS_ID_MAX);

    // collect tri nodes
    moris_index tTriCount = 0;
    for(moris::uint  i = 0; i<6; i++)
    {
        moris::uint tChildInd = 2*i+1;


        // iterate through cells in this side set
        for(moris::uint  j = 0; j <tElementIndAndSideOrds(tChildInd).n_rows(); j++)
        {
            // get the mtk cell
            moris::mtk::Cell const & tCell = mBackgroundMesh.get_mtk_cell(tElementIndAndSideOrds(tChildInd)(j,0));

            moris::Cell<moris::mtk::Vertex const *> tVerticesOnSide = tCell.get_vertices_on_side_ordinal(tElementIndAndSideOrds(tChildInd)(j,1));

            Matrix<IdMat> tVertexIds( 1, tVerticesOnSide.size());
            for(moris::uint k =0; k < 3; k++)
            {
                tVertexIds(k) = tVerticesOnSide(k)->get_id();
                tNodesOnBoundary(tNodeCount) = tVertexIds(k);
                tNodeCount++;
            }

            tTri3ElemToNode.get_row(tTriCount) = tVertexIds.get_row(0);
            tTriCount++;
        }
    }

    // Collect the nodes on the surface and construct a quad 4 / tri 3 mesh
    // quad 4s
    uint tQuadCount = 0;
    moris::Matrix<moris::IdMat> tQuad4ElemToNode(tNumNoChildCells,4);

    for(moris::uint i = 0; i < 6; i++)
    {
        moris::uint tNoChildInd = 2*i;

        for(moris::uint  j = 0; j <tElementIndAndSideOrds(tNoChildInd).n_rows(); j++)
        {
            // get the mtk cell
            moris::mtk::Cell const & tCell = mBackgroundMesh.get_mtk_cell(tElementIndAndSideOrds(tNoChildInd)(j,0));

            moris::Cell<moris::mtk::Vertex const *> tVerticesOnSide = tCell.get_vertices_on_side_ordinal(tElementIndAndSideOrds(tNoChildInd)(j,1));

            Matrix<IdMat> tVertexIds( 1, tVerticesOnSide.size());

            for(moris::uint k =0; k < 4; k++)
            {
                tVertexIds(k) = tVerticesOnSide(k)->get_id();
                tNodesOnBoundary(tNodeCount) = tVertexIds(k);
                tNodeCount++;
            }
            tQuad4ElemToNode.get_row(tQuadCount) = tVertexIds.get_row(0);
            tQuadCount++;
        }
    }

    // triangulate quad 4s
    //
    // allocate ids for triangulation
    moris_id tNodeId = mBackgroundMesh.allocate_entity_ids(tQuad4ElemToNode.n_rows(),EntityRank::NODE);

    // template for splitting the quad 4 into tris
    moris::Matrix<moris::IdMat> tTriangulatedQuadMap = {{0,1,4},{1,2,4},{2,3,4},{3,0,4}};

    moris::Matrix<moris::IdMat> tNodeIdsForTemplate(1,5);

    moris::Matrix<moris::IdMat>  tNewVertexId(tQuad4ElemToNode.n_rows(),1);
    moris::Matrix<moris::DDRMat> tNewVertexCoords(tQuad4ElemToNode.n_rows(),3);
    moris::Matrix<moris::DDRMat> tParamCoordCenterQuad( {{ 0.0,  0.0}} );

    for(moris::uint i = 0; i<tQuad4ElemToNode.n_rows(); i++)
    {
        // assign ids needed for template
        tNodeIdsForTemplate({0,0},{0,3})              = tQuad4ElemToNode.get_row(i);
        tNodeIdsForTemplate(4)                        = tNodeId;
        tNewVertexId(i)                               = tNodeId;
        tNodeId++;

        // turn quad into triangles
        moris::Matrix<moris::IdMat> tTriangulatedQuad = reindex_matrix(tTriangulatedQuadMap,0,tNodeIdsForTemplate);

        // compute vertex coordinate
        moris::Matrix<moris::DDRMat> tNodeCoordsOnFace(4,3);
        for(moris::uint j = 0; j < 4; j++)
        {
            moris_index tNodeIndex = tMeshData.get_loc_entity_ind_from_entity_glb_id(tNodeIdsForTemplate(j),EntityRank::NODE);
            tNodeCoordsOnFace.get_row(j) = tMeshData.get_node_coordinate(tNodeIndex).get_row(0);
        }

        // bilinearly interpolate to the center of this face fi
        moris::Matrix<moris::DDRMat> tNewNodeCoordinates;
        xtk::Interpolation::bilinear_interpolation(tNodeCoordsOnFace, tParamCoordCenterQuad, tNewNodeCoordinates);
        tNewVertexCoords.get_row(i) = tNewNodeCoordinates.get_row(0);

        // add quad to tri 3s
        for(moris::uint  j = 0; j < 4; j++)
        {
            tTri3ElemToNode.get_row(tTriCount) = tTriangulatedQuad.get_row(j);
            tTriCount++;
        }
    }

    // add interface nodes to nodes on boundary
    tNodesOnBoundary({tNodeCount, tNodeCount + tInterfaceNodes(aPhaseIndex).numel()-1},{0,0}) = moris::trans(tInterfaceNodes(aPhaseIndex));
    tNodeCount = tNodeCount + tInterfaceNodes(aPhaseIndex).numel();

    // add interface facets

    // iterate through cells in interface
    for(moris::uint  i = 0; i <tInterfaceElemIndandSideOrd.n_rows(); i++)
    {
        // get the mtk cell

        moris::mtk::Cell const & tCell = mBackgroundMesh.get_mtk_cell(tInterfaceElemIndandSideOrd(i,0));
        moris::Cell<moris::mtk::Vertex const *> tVerticesOnSide = tCell.get_vertices_on_side_ordinal(tInterfaceElemIndandSideOrd(i,1));

        Matrix<IdMat> tVertexIds( 1, tVerticesOnSide.size());
        for(moris::uint k =0; k < 3; k++)
        {
            tVertexIds(k) = tVerticesOnSide(k)->get_id();
        }

        tTri3ElemToNode.get_row(tTriCount) = tVertexIds.get_row(0);
        tTriCount++;
    }

    // assign element ids
    //fixme: Make these ids unique across processors
    moris::Matrix<moris::IdMat> tChildFaceElementIds(tTriCount,1);
    for(moris::uint  i = 0; i <tChildFaceElementIds.numel(); i++)
    {
        tChildFaceElementIds(i) = (moris_id)i+1;
    }


    // Interface elements
    moris::mtk::MtkBlockSetInfo tTri3Block;
    tTri3Block.mCellIdsInSet = &tChildFaceElementIds;
    tTri3Block.mBlockSetName = "tri_surf";
    tTri3Block.mBlockSetTopo = CellTopology::TRI3;
    tMtkMeshSets.add_block_set(&tTri3Block);

    // Convert to vertex indices
    tNodesOnBoundary.resize(tNodeCount,1);
    moris::Matrix<moris::IdMat> tNodeMap;

    moris::unique( tNodesOnBoundary,tNodeMap);
    moris::Matrix<moris::IndexMat> tNodeIndices(tNodeMap.numel(),1);

    for(moris::uint i = 0; i <tNodeMap.numel(); i++)
    {
        tNodeIndices(i) = mBackgroundMesh.get_loc_entity_ind_from_entity_glb_id(tNodeMap(i),EntityRank::NODE);
    }

    // Get the node coordinates
    moris::Matrix<moris::DDRMat> tNodeCoordinates = mBackgroundMesh.get_selected_node_coordinates_loc_inds(tNodeIndices);

    // add nodes created during quad 4
    uint tNumNodeNoMidside = tNodeCoordinates.n_rows();
    uint tNumMidsideNodes = tNewVertexId.numel();
    uint tTotalNumNodes  =  tNumNodeNoMidside + tNumMidsideNodes;
    tNodeCoordinates.resize(tNumNodeNoMidside + tNumMidsideNodes,3);
    tNodeCoordinates({tNumNodeNoMidside,tTotalNumNodes-1},{0,2}) = tNewVertexCoords.matrix_data();

    tNodeMap.resize(tTotalNumNodes, 1);
    tNodeMap({tNumNodeNoMidside, tTotalNumNodes-1},{0,0}) = tNewVertexId.matrix_data();

    // make vertices consecutive for obj output
    std::unordered_map<moris_index, moris_index> tObjIndex;
    for(moris::uint i = 0; i < tNodeMap.numel(); i++)
    {
        tObjIndex[tNodeMap(i)] = i+1;
    }

    // renumber vertex
    for(moris::uint i  = 0; i <tTri3ElemToNode.n_rows(); i++)
    {
        for(moris::uint  j = 0; j<tTri3ElemToNode.n_cols(); j++)
        {
            auto tIter = tObjIndex.find(tTri3ElemToNode(i,j));
            tTri3ElemToNode(i,j) = tIter->second;
        }
    }

    // write to an obj file
    std::string tParObjPath =  moris::make_path_parallel( aOutputFile );

    std::ofstream tOFS (tParObjPath, std::ofstream::out);

    tOFS << std::setprecision (15);

    tOFS << "# XTK OBJ EXTRACTION"<<std::endl;

    // add vertices
    for(moris::uint i = 0 ; i < tNodeCoordinates.n_rows(); i++)
    {
        tOFS << std::scientific << "v "<< tNodeCoordinates(i,0) <<" "<< tNodeCoordinates(i,1) << " "<< tNodeCoordinates(i,2) <<std::endl;
    }

    // add facets
    for(moris::uint i = 0; i < tTri3ElemToNode.n_rows(); i++)
    {
        tOFS << "f "<< tTri3ElemToNode(i,0) <<" "<< tTri3ElemToNode(i,1)<< " "<< tTri3ElemToNode(i,2) <<std::endl;
    }

    tOFS.close();
}

moris::mtk::Integration_Mesh*
Model::get_output_mesh(Output_Options const & aOutputOptions)

{
    // start timing on this output
    std::clock_t start = std::clock();

    // create the output mtk mesh
    moris::mtk::Integration_Mesh* tOutputMesh = construct_output_mesh(aOutputOptions);

    if(moris::par_rank() == 0  && mVerbose)
    {
        std::cout<<"XTK: Mesh output completed in "<< (std::clock() - start) / (double)(CLOCKS_PER_SEC)<<" s."<<std::endl;
    }
    return tOutputMesh;
}

//------------------------------------------------------------------------------

moris::uint
Model::get_spatial_dim() const
{
    return mModelDimension;
}

//------------------------------------------------------------------------------

moris::uint
Model::get_num_elements_total()
{
    MORIS_ASSERT(mDecomposed,"Prior to using get_num_elements, the decomposition process must be finished");

    return mBackgroundMesh.get_num_entities(EntityRank::ELEMENT);
}

moris::uint
Model::get_num_elements_unzipped()
{
    return this->get_num_elements_total() - mCutMesh.get_num_child_meshes();
}
//------------------------------------------------------------------------------
moris_index
Model::get_cell_xtk_index(moris_id aCellId)
{
    auto tIter = mCellGlbToLocalMap.find(aCellId);
    MORIS_ASSERT(tIter != mCellGlbToLocalMap.end(),"Id not in map");
    return tIter->second;
}
//------------------------------------------------------------------------------
moris::Matrix<moris::IndexMat>
Model::get_element_to_subphase()
{
    // child mesh subphases
    moris::uint tNumElem = this->get_num_elements_total();
    Matrix<IndexMat> tSubphase(1,tNumElem);
    mCutMesh.populate_subphase_vector(tSubphase);

    // populate the non intersected background cells
    for(moris::uint i = 0; i < mBackgroundMesh.get_mesh_data().get_num_entities(EntityRank::ELEMENT); i++)
    {
        tSubphase(i) = i;
    }

    return tSubphase;
}

//------------------------------------------------------------------------------

moris_id
Model::get_subphase_id(moris_id aSubphaseIndex)
{
    if(this->subphase_is_in_child_mesh(aSubphaseIndex))
    {
        return mCutMesh.get_subphase_id(aSubphaseIndex);
    }
    else
    {
        return mBackgroundMesh.get_glb_entity_id_from_entity_loc_index(aSubphaseIndex,EntityRank::ELEMENT);
    }
}

moris_index
Model::get_subphase_index(moris_id aSubphaseId)
{

    auto tIter = mGlobalToLocalSubphaseMap.find(aSubphaseId);
    MORIS_ASSERT(tIter != mGlobalToLocalSubphaseMap.end(),"Subphase id not in map");
    return tIter->second;
}

//------------------------------------------------------------------------------

moris::mtk::Integration_Mesh*
Model::construct_output_mesh( Output_Options const & aOutputOptions )
{

    // start timing on this decomposition
    std::clock_t start = std::clock();

    // Get mesh information ready for outputting
    // Package element to Node Connectivity
    moris::uint tSpatialDim = mBackgroundMesh.get_mesh_data().get_spatial_dim();

    // Children element nodes connected to elements
    moris::Cell<moris::Matrix<moris::IdMat>>  tElementToNodeChildrenByPhase = mCutMesh.get_full_element_to_node_by_phase_glob_ids(mGeometryEngine->get_num_bulk_phase(),mBackgroundMesh.get_mesh_data());

    // Child element ids
    moris::Cell<moris::Matrix<moris::IdMat>>  tChildElementsByPhase = mCutMesh.get_child_elements_by_phase(mGeometryEngine->get_num_bulk_phase(),mBackgroundMesh.get_mesh_data());

    // Parent elements without children
    Cell<moris::Matrix<moris::IdMat>>  tElementNoChildrenIdsByPhase = mBackgroundMesh.get_all_non_intersected_elements_by_phase(mGeometryEngine->get_num_bulk_phase());

    // Connectivity of parent elements without children
    Cell<moris::Matrix<moris::IdMat>>  tElementToNodeNoChildrenByPhase = mBackgroundMesh.get_non_intersected_element_to_node_by_phase(mGeometryEngine->get_num_bulk_phase());

    // Node map  of nodes in the phase we are interested in
    moris::Matrix<moris::IndexMat> tOutputtedNodeInds;
    moris::Matrix<moris::IdMat>  tLocalToGlobalNodeMap = this->get_node_map_restricted_to_output_phases(aOutputOptions,tOutputtedNodeInds);

    // All node coordinates
    moris::Matrix<moris::DDRMat> tNodeCoordinates = mBackgroundMesh.get_selected_node_coordinates_loc_inds(tOutputtedNodeInds);


    // Number of bulk phases
    uint tNumBulkPhases = mGeometryEngine->get_num_phases();

    // Get non-interescted parent elements by phase
    Cell<moris::Matrix<moris::IdMat>> tNoChildElementsByPhase = mBackgroundMesh.get_all_non_intersected_elements_by_phase(tNumBulkPhases);

    // combination of the elements by phase (if specified)
    Cell<moris::Matrix<moris::IdMat>> tCombinedElementsByPhase(tNoChildElementsByPhase.size());
    if(!aOutputOptions.mSeparateInterfaceBlock)
    {
        MORIS_ASSERT(mBackgroundMesh.get_parent_cell_topology() == CellTopology::TET4," Combining the interface block w/ non-interface block only valid on tet background mesh");

        tCombinedElementsByPhase = combine_interface_and_non_interface_blocks(tChildElementsByPhase,tNoChildElementsByPhase);

    }

    // Interface nodes
    Cell<moris::Matrix<moris::IndexMat>> tInterfaceNodes = mBackgroundMesh.get_interface_nodes_glob_ids();

    // Assemble geometry data as field for mesh output
    moris::Cell< moris::Matrix < moris::DDRMat > > tGeometryFieldData = assemble_geometry_data_as_mesh_field(tOutputtedNodeInds);

    // Give the geometry data a name
    moris::Cell<std::string> tGeometryFieldNames = assign_geometry_data_names();

    // Get rank of the geometry data field
    moris::Cell < enum moris::EntityRank > tFieldRanks =  assign_geometry_data_field_ranks();

    // number of phases being output
    moris::uint tNumPhasesOutput = get_num_phases_to_output(aOutputOptions);

    // Set up field data structure for MTK input
    moris::mtk::MtkFieldsInfo tFieldsInfo;
    moris::Cell<moris::mtk::Scalar_Field_Info<DDRMat>> tGeometryFields(tGeometryFieldData.size());

    for(uint i = 0; i <tGeometryFieldData.size(); i++)
    {
        tGeometryFields(i).set_field_name(tGeometryFieldNames(i));
        tGeometryFields(i).set_field_entity_rank(moris::EntityRank::NODE);
        tGeometryFields(i).add_field_data(&tLocalToGlobalNodeMap, &tGeometryFieldData(i));
        tFieldsInfo.mRealScalarFields.push_back(&tGeometryFields(i));
    }

    // External Fields - real cell fields
    uint tNumExtRealCellScalarFields = aOutputOptions.mRealElementExternalFieldNames.size();
    moris::Cell<moris::mtk::Scalar_Field_Info<DDRMat>> tExternalRealCellScalarFields(tNumExtRealCellScalarFields);
    for(uint i = 0; i<tNumExtRealCellScalarFields; i++)
    {
        tExternalRealCellScalarFields(i).set_field_name(aOutputOptions.mRealElementExternalFieldNames(i));
        tExternalRealCellScalarFields(i).set_field_entity_rank(moris::EntityRank::ELEMENT);
        add_field_for_mesh_input(&tExternalRealCellScalarFields(i),tFieldsInfo);
    }

    // External Fields - real vertex fields
    uint tNumExtRealVertexScalarFields = aOutputOptions.mRealNodeExternalFieldNames.size();
    moris::Cell<moris::mtk::Scalar_Field_Info<DDRMat>> tExternalRealVertexScalarFields(tNumExtRealVertexScalarFields);
    for(uint i = 0; i<tNumExtRealVertexScalarFields; i++)
    {
        tExternalRealVertexScalarFields(i).set_field_name(aOutputOptions.mRealNodeExternalFieldNames(i));
        tExternalRealVertexScalarFields(i).set_field_entity_rank(moris::EntityRank::NODE);
        add_field_for_mesh_input(&tExternalRealVertexScalarFields(i),tFieldsInfo);
    }


    // sensitivity fields
    moris::Cell<moris::Matrix<DDRMat>> adxdpData;
    moris::Cell<std::string>           adxdpNames;
    moris::Cell<moris::Matrix<DDRMat>> aDesVars;
    moris::Cell<std::string>           aDesVarsName;
    moris::Matrix<moris::DDRMat>       aNumDesVars;
    std::string                        aNumDesVarsName;
    moris::Cell<moris::mtk::Scalar_Field_Info<DDRMat>> tdxdpDataFields;
    moris::Cell<moris::mtk::Scalar_Field_Info<DDRMat>> tDesVarFields;
    moris::Cell<moris::mtk::Scalar_Field_Info<DDRMat>> tNumDesVarsField;
    if(aOutputOptions.mPackageDxDpSparsely && mSensitivity)
    {
        this->extract_interface_sensitivity_sparse(tOutputtedNodeInds,adxdpData,adxdpNames,aDesVars,aDesVarsName,aNumDesVars,aNumDesVarsName);

        tdxdpDataFields.resize(adxdpData.size());
        tDesVarFields.resize(aDesVars.size());
        tDesVarFields.resize(1);

        // place into a field
        for(moris::uint  i = 0; i <tdxdpDataFields.size(); i++)
        {
            tdxdpDataFields(i).set_field_name(adxdpNames(i));
            tdxdpDataFields(i).set_field_entity_rank(moris::EntityRank::NODE);
            tdxdpDataFields(i).add_field_data( &tLocalToGlobalNodeMap, &adxdpData(i));
            add_field_for_mesh_input(&tdxdpDataFields(i),tFieldsInfo);
        }

        for(moris::uint  i = 0; i <tDesVarFields.size(); i++)
        {
            tDesVarFields(i).set_field_name(aDesVarsName(i));
            tDesVarFields(i).set_field_entity_rank(moris::EntityRank::NODE);
            tDesVarFields(i).add_field_data( &tLocalToGlobalNodeMap, &aDesVars(i));
            add_field_for_mesh_input(&tDesVarFields(i),tFieldsInfo);
        }

        tNumDesVarsField(0).set_field_name(aNumDesVarsName);
        tNumDesVarsField(0).set_field_entity_rank(moris::EntityRank::NODE);
        tNumDesVarsField(0).add_field_data( &tLocalToGlobalNodeMap, &aNumDesVars);
        add_field_for_mesh_input(&tNumDesVarsField(0),tFieldsInfo);

    }

    if(aOutputOptions.mPackageDxDpDensely && mSensitivity)
    {
        this->extract_interface_sensitivity_dense(tOutputtedNodeInds,adxdpData,adxdpNames);

        tdxdpDataFields.resize(adxdpData.size());

        // place into a field
        for(moris::uint  i = 0; i <tdxdpDataFields.size(); i++)
        {
            tdxdpDataFields(i).set_field_name(adxdpNames(i));
            tdxdpDataFields(i).set_field_entity_rank(moris::EntityRank::NODE);
            tdxdpDataFields(i).add_field_data( &tLocalToGlobalNodeMap, &adxdpData(i));
            add_field_for_mesh_input(&tdxdpDataFields(i),tFieldsInfo);
        }
    }


    //TODO: implement node owner (currently set to owned by this proc)
    //    moris::Matrix<moris::IdMat> tNodeOwner(1,tOutputtedNodeInds.numel(),moris::par_rank());

    moris::Matrix<moris::IdMat> tNodeOwner = mBackgroundMesh.get_vertices_owner(tOutputtedNodeInds);

    // Set up mesh sets
    // Initialize Sets information structure
    moris::mtk::MtkSetsInfo tMtkMeshSets;

    //
    moris::uint tNumBlocksPerPhase = 2;
    if(!aOutputOptions.mSeparateInterfaceBlock)
    {
        MORIS_ASSERT(mBackgroundMesh.get_parent_cell_topology() == CellTopology::TET4," Combining the interface block w/ non-interface block only valid on tet background mesh");
        tNumBlocksPerPhase = 1;
    }

    // Setup block sets
    Cell<moris::mtk::MtkBlockSetInfo> tBlockSets(tNumPhasesOutput*tNumBlocksPerPhase);
    uint tCount= 0;

    for(uint i = 0; i <tNumBulkPhases; i++)
    {
        if(aOutputOptions.output_phase(i) && aOutputOptions.mSeparateInterfaceBlock)
        {
            // Children of material phase i
            tBlockSets(tCount).mCellIdsInSet = &tChildElementsByPhase(i);
            tBlockSets(tCount).mBlockSetName = "child_"+std::to_string(i);
            tBlockSets(tCount).mBlockSetTopo = mCutMesh.get_child_element_topology();

            tMtkMeshSets.add_block_set(&tBlockSets(tCount));
            tCount++;

            // Children of material phase i
            tBlockSets(tCount).mCellIdsInSet = &tNoChildElementsByPhase(i);
            tBlockSets(tCount).mBlockSetName = "parent_"+std::to_string(i);
            tBlockSets(tCount).mBlockSetTopo = mBackgroundMesh.get_parent_cell_topology();

            tMtkMeshSets.add_block_set(&tBlockSets(tCount));
            tCount++;
        }

        else if(aOutputOptions.output_phase(i) && !aOutputOptions.mSeparateInterfaceBlock)
        {
            // Children of material phase i
            tBlockSets(tCount).mCellIdsInSet = &tCombinedElementsByPhase(i);
            tBlockSets(tCount).mBlockSetName = "phase_"+std::to_string(i);
            tBlockSets(tCount).mBlockSetTopo = mBackgroundMesh.get_parent_cell_topology();

            tMtkMeshSets.add_block_set(&tBlockSets(tCount));
            tCount++;
        }
    }

    // Interface elements
    moris::Matrix<moris::IndexMat> tInterfaceElements(0,0);
    moris::Matrix<moris::IndexMat> tInterfaceElementIds(0,0);
    moris::mtk::MtkBlockSetInfo tUnzippedInterfaceBlockSet;
    if(mUnzipped && aOutputOptions.mHaveInterface)
    {
        // get the interface elements local node index element connectivity
        tInterfaceElements = mCutMesh.get_extracted_interface_elements_loc_inds();

        mBackgroundMesh.convert_loc_entity_ind_to_glb_entity_ids(EntityRank::NODE,tInterfaceElements);

        tInterfaceElementIds = mCutMesh.get_interface_element_ids();

        tUnzippedInterfaceBlockSet.mCellIdsInSet = &tInterfaceElementIds;
        tUnzippedInterfaceBlockSet.mBlockSetName = "interface";
        tUnzippedInterfaceBlockSet.mBlockSetTopo = CellTopology::PRISM6;
        tMtkMeshSets.add_block_set(&tUnzippedInterfaceBlockSet);
    }

    // propogate background mesh node sets
    moris::Cell<moris::Matrix<IndexMat>> tBackgroundNodeSetData;
    moris::Cell<moris::mtk::MtkNodeSetInfo> tBackgroundNodeSets;
    if(aOutputOptions.mAddNodeSets)
    {
        tBackgroundNodeSets = propogate_background_node_sets(tBackgroundNodeSetData,aOutputOptions);

        for(moris::uint i = 0; i<tBackgroundNodeSets.size(); i++)
        {
            tMtkMeshSets.add_node_set(&tBackgroundNodeSets(i));
        }
    }


    moris::Cell<moris::mtk::MtkNodeSetInfo> tInterfaceNodeSets(tInterfaceNodes.size());
    if(aOutputOptions.mHaveInterface)
    {


        for(uint i = 0; i<tInterfaceNodes.size(); i++)
        {
            tInterfaceNodeSets(i).mNodeIds     = &tInterfaceNodes(i);
            tInterfaceNodeSets(i).mNodeSetName = "inodes_" +std::to_string(i) ;
            tMtkMeshSets.add_node_set(&tInterfaceNodeSets(i));
        }
    }

    // Get the packaged interface side sets from the cut mesh
    Cell<moris::Matrix<moris::IdMat>> tInterfaceElemIdandSideOrd;
    Cell<std::string> tInterfaceNames;
    Cell<moris::mtk::MtkSideSetInfo> tInterfaceSideSets;
    if(aOutputOptions.mHaveInterface)
    {
        this->setup_interface_single_side_sets(aOutputOptions,tInterfaceElemIdandSideOrd,tInterfaceNames);

        tInterfaceSideSets.resize(tInterfaceElemIdandSideOrd.size());
        for(moris::uint i = 0; i < tInterfaceElemIdandSideOrd.size(); i++)
        {
            tInterfaceSideSets(i).mElemIdsAndSideOrds = &tInterfaceElemIdandSideOrd(i);
            tInterfaceSideSets(i).mSideSetName        = tInterfaceNames(i);
            tMtkMeshSets.add_side_set(&tInterfaceSideSets(i));
        }
    }

    // propogate side sets from background mesh
    moris::Cell<moris::Matrix<IndexMat>> tSideSetData;
    moris::Cell<moris::mtk::MtkSideSetInfo> tBackgroundSideSets;
    if(aOutputOptions.mAddSideSets)
    {
        // collect information about background side set
        tBackgroundSideSets = this->propogate_background_side_sets(tSideSetData,aOutputOptions);

        // add to mesh input structure
        for(moris::uint i = 0; i<tBackgroundSideSets.size(); i++)
        {
            tMtkMeshSets.add_side_set(&tBackgroundSideSets(i));
        }
    }

    // Mesh data input structure (with multi element type mesh)
    moris::mtk::MtkMeshData tMeshDataInput;

    moris::uint tNumElemTypes = tNumPhasesOutput*2;
    if(mUnzipped)
    {
        tNumElemTypes = tNumElemTypes + 1;
    }

    tMeshDataInput.ElemConn             = moris::Cell<moris::Matrix<IdMat>*>(tNumElemTypes);
    tMeshDataInput.LocaltoGlobalElemMap = moris::Cell<moris::Matrix<IdMat>*>(tNumElemTypes);
    tMeshDataInput.CellTopology         = moris::Cell<enum CellTopology>(tNumElemTypes,CellTopology::INVALID);

    tMeshDataInput.Verbose                 = mVerbose;
    tMeshDataInput.SpatialDim              = &tSpatialDim;

    tCount = 0;
    for(moris::uint  i = 0 ; i <mGeometryEngine->get_num_bulk_phase(); i++)
    {
        if(aOutputOptions.output_phase(i))
        {
            tMeshDataInput.ElemConn(tCount)             = &tElementToNodeChildrenByPhase(i);
            tMeshDataInput.LocaltoGlobalElemMap(tCount) = &tChildElementsByPhase(i);
            tMeshDataInput.CellTopology(tCount)         = mCutMesh.get_child_element_topology();
            tCount++;
            tMeshDataInput.ElemConn(tCount)             = &tElementToNodeNoChildrenByPhase(i);
            tMeshDataInput.LocaltoGlobalElemMap(tCount) = &tElementNoChildrenIdsByPhase(i);
            tMeshDataInput.CellTopology(tCount)         = mBackgroundMesh.get_parent_cell_topology();
            tCount++;
        }
    }

    tMeshDataInput.NodeCoords              = &tNodeCoordinates;
    tMeshDataInput.FieldsInfo              = &tFieldsInfo;
    tMeshDataInput.LocaltoGlobalNodeMap    = &tLocalToGlobalNodeMap;
    tMeshDataInput.SetsInfo                = &tMtkMeshSets;
    tMeshDataInput.MarkNoBlockForIO        = false;
    tMeshDataInput.CreateAllEdgesAndFaces  = true;
    tMeshDataInput.AutoAuraOptionInSTK     = false;

    //Add clustering information
    moris::mtk::Cell_Cluster_Input tCellClusterInput;                // Cell clusters
    moris::Cell<Matrix<IdMat>>     tClusterCellIds;                  // Cell cluster Ids
    moris::mtk::Side_Cluster_Input tInterfaceSideClusterInput;       // Side clusters
    moris::Cell<Matrix<IdMat>>     tInterfaceCellIdsandSideOrds;     // side cluster ids and side ordinals
    moris::Cell<Matrix<DDRMat>>    tInterfaceSideClusterParamCoords; // side cluster vertex parametric coordinates
//
    if(aOutputOptions.mAddClusters)
    {
        // cell clustering
        this->setup_cell_clusters_for_output(tCellClusterInput,aOutputOptions,tClusterCellIds);

        tMeshDataInput.CellClusterInput = &tCellClusterInput;

        MORIS_ASSERT(mGeometryEngine->get_num_geometries() == 1,"This has not been setup for multi geometry problems.");

        //fixme: support changes to interface side
        this->setup_interface_side_cluster(tInterfaceNames(0), tInterfaceSideClusterInput, aOutputOptions, tInterfaceCellIdsandSideOrds, tInterfaceSideClusterParamCoords);

        tMeshDataInput.SideClusterInput = &tInterfaceSideClusterInput;
    }

    // Interface elements
    if(mUnzipped)
    {
        tMeshDataInput.ElemConn(tNumElemTypes-1)             = &tInterfaceElements;
        tMeshDataInput.LocaltoGlobalElemMap(tNumElemTypes-1) = &tInterfaceElementIds;
        tMeshDataInput.CellTopology(tNumElemTypes-1)         = CellTopology::PRISM6;
    }

    // add parallel information
    moris::mtk::Visualization_STK tVizTool;
    if(aOutputOptions.mAddParallelFields && par_size()>1)
    {
        moris::mtk::MtkFieldsInfo* tParFields = tVizTool.setup_parallel_cell_fields_for_declaration();
        tFieldsInfo.combine_fields_info(*tParFields);
    }


    if(moris::par_rank() == 0 && mVerbose)
    {
        std::cout<<"XTK: Mesh data setup completed in " <<(std::clock() - start) / (double)(CLOCKS_PER_SEC)<<" s."<<std::endl;
    }


    start = std::clock();

    // cast background mesh to an interpolation mesh and pass in
    moris::mtk::Integration_Mesh* tMeshData = nullptr;
    if(aOutputOptions.mAddClusters)
    {
        moris::mtk::Interpolation_Mesh* tInterpMesh = dynamic_cast<moris::mtk::Interpolation_Mesh*>(&mBackgroundMesh.get_mesh_data());

        tMeshData = moris::mtk::create_integration_mesh( MeshType::STK, tMeshDataInput, tInterpMesh );
    }
    else
    {
        tMeshData = moris::mtk::create_integration_mesh( MeshType::STK, tMeshDataInput);
    }

    if(aOutputOptions.mAddParallelFields && par_size()>1)
    {
        tVizTool.populate_parallel_cell_fields_on_mesh(tMeshData);
    }

    if(moris::par_rank() == 0 && mVerbose)
    {
        std::cout<<"XTK: Write to MTK completed in " <<(std::clock() - start) / (double)(CLOCKS_PER_SEC)<<" s."<<std::endl;
        tMeshDataInput.print_summary();
    }

    return tMeshData;

}

//------------------------------------------------------------------------------

moris::Matrix<moris::IndexMat>
Model::get_node_map_restricted_to_output_phases(Output_Options const & aOutputOptions,
                                                moris::Matrix<moris::IndexMat> & aOutputtedNodeInds)
{
    moris::Matrix<moris::IndexMat> tNodeMap = mBackgroundMesh.get_local_to_global_map(EntityRank::NODE);

    moris_index tMyProcRank = par_rank();

    aOutputtedNodeInds.resize(tNodeMap.n_rows(),tNodeMap.n_cols());
    moris::uint tNumNodes = tNodeMap.numel();
    // if we are returning all phases there is no reason to restrict the map
    if(aOutputOptions.output_all_phases())
    {

        for(moris::uint i = 0; i <tNumNodes; i++)
        {
            aOutputtedNodeInds(i) = mBackgroundMesh.get_loc_entity_ind_from_entity_glb_id(tNodeMap(i),EntityRank::NODE);
        }

        return tNodeMap;
    }

    else
    {
        moris::Matrix<moris::IndexMat> tRestrictedNodeMap(tNodeMap.n_rows(),tNodeMap.n_cols());

        moris::uint tCount = 0;
        for(moris::uint i = 0; i <tNumNodes; i++)
        {
            if(output_node(i,aOutputOptions) && mBackgroundMesh.get_vertex_owner(i) == tMyProcRank)
            {

                moris::size_t tPhaseIndex = 0;
                moris_index   tVertexIndex = mBackgroundMesh.get_loc_entity_ind_from_entity_glb_id(tNodeMap(i),EntityRank::NODE);
                mGeometryEngine->get_phase_index(tVertexIndex,tPhaseIndex);
                aOutputtedNodeInds(tCount) = tVertexIndex;
                tRestrictedNodeMap(tCount) = tNodeMap(i);

                tCount++;
            }
        }


        tRestrictedNodeMap.resize(1,tCount);
        aOutputtedNodeInds.resize(1,tCount);

        return tRestrictedNodeMap;
    }

}

//------------------------------------------------------------------------------
void
Model::setup_interface_single_side_sets(Output_Options const & aOutputOptions,
                                        Cell<moris::Matrix<moris::IdMat>> & aCellIdsAndSideOrds,
                                        Cell<std::string> &                 aInterfaceSetNames)
{
    std::string tSetNameBase = "iside_";

    for(moris_index iG = 0; iG< (moris_index)mGeometryEngine->get_num_geometries(); iG++)
    {
        for(moris_index iP0 = 0; iP0 < (moris_index)mGeometryEngine->get_num_bulk_phase(); iP0++)
        {
            for(moris_index iP1 = iP0+1; iP1 < (moris_index)mGeometryEngine->get_num_bulk_phase(); iP1++)
            {
                if(aOutputOptions.output_phase(iP0) && aOutputOptions.output_phase(iP1))
                {
                    std::string tSetName = tSetNameBase+"g_"+ std::to_string(iG) + "_p0_"+std::to_string(iP0)+"_p1_"+std::to_string(iP1);

                    aInterfaceSetNames.push_back(tSetName);
                    aCellIdsAndSideOrds.push_back(mCutMesh.pack_interface_sides(iG,iP0,iP1));
                }
            }
        }
    }
}

//------------------------------------------------------------------------------

moris::Cell<moris::mtk::MtkNodeSetInfo>
Model::propogate_background_node_sets(moris::Cell<moris::Matrix<IndexMat>> & aNodeSetData,
                                      Output_Options const & aOutputOptions)
{
    // access background mesh data
    moris::mtk::Mesh & tMeshData = mBackgroundMesh.get_mesh_data();

    // get all node set names in background mesh
    moris::Cell<std::string> tSetNames = tMeshData.get_set_names(EntityRank::NODE);

    // allocate output
    aNodeSetData = moris::Cell<moris::Matrix<IndexMat>>(tSetNames.size());
    moris::Cell<moris::mtk::MtkNodeSetInfo> tNodeSetInfo(tSetNames.size());
    for(moris::uint i = 0; i <tSetNames.size(); i++)
    {
        moris::uint tCount = 0;
        moris::Cell<moris::mtk::Vertex const *> tNodesInSetInds = tMeshData.get_vertices_in_vertex_set_no_aura(tSetNames(i));

        aNodeSetData(i) = moris::Matrix<moris::IndexMat>(tNodesInSetInds.size(),1);

        for(moris::uint iNode =0; iNode<tNodesInSetInds.size(); iNode++)
        {
            moris_index tVertexInd = tNodesInSetInds(iNode)->get_index();

            if(this->output_node(tVertexInd,aOutputOptions))
            {
                aNodeSetData(i)(tCount) = tNodesInSetInds(iNode)->get_index();
                tCount++;
            }
        }


        aNodeSetData(i).resize(tCount,1);

        mBackgroundMesh.convert_loc_entity_ind_to_glb_entity_ids(EntityRank::NODE,aNodeSetData(i));

        tNodeSetInfo(i).mNodeIds = &aNodeSetData(i);
        tNodeSetInfo(i).mNodeSetName = tSetNames(i);

    }
    return tNodeSetInfo;

}

//------------------------------------------------------------------------------

moris::Cell<moris::mtk::MtkSideSetInfo>
Model::propogate_background_side_sets(moris::Cell<moris::Matrix<IndexMat>> & aSideSetElemIdsAndOrds,
                                      Output_Options const & aOutputOptions)
{
    // access background mesh data
    moris::mtk::Mesh & tMeshData = mBackgroundMesh.get_mesh_data();

    // get all side set names in background mesh
    moris::Cell<std::string> tSetNames = tMeshData.get_set_names(tMeshData.get_facet_rank());

    // remove internal side sets which show up with a generated string
    tSetNames = check_for_and_remove_internal_seacas_side_sets(tSetNames);

    // allocate output side sets
    moris::Cell<moris::mtk::MtkSideSetInfo> tSideSets(2*tSetNames.size());
    aSideSetElemIdsAndOrds = moris::Cell<moris::Matrix<IndexMat>>(2*tSetNames.size());


    for(moris::uint i = 0; i <tSetNames.size(); i++)
    {
        moris::uint tNoChildInd = 2*i;
        moris::uint tChildInd = 2*i+1;

        this->propogate_background_side_set(tSetNames(i),tNoChildInd,tChildInd,aSideSetElemIdsAndOrds,tSideSets,aOutputOptions,false);
    }


    return tSideSets;
}

void
Model::propogate_background_side_set( std::string             const &             aSideSetName,
                                      moris::moris_index                          aNoChildIndex,
                                      moris::moris_index                          aChildIndex,
                                      moris::Cell<moris::Matrix<IndexMat>>      & aElementIdsAndSideOrd,
                                      moris::Cell<moris::mtk::MtkSideSetInfo>   & aSideSetData,
                                      Output_Options          const             & aOutputOptions,
                                      bool                                        aOutputIndices)
{

    // access background mesh data
    moris::mtk::Mesh & tMeshData = mBackgroundMesh.get_mesh_data();

    // declare matrices used through
    moris::Matrix< moris::IndexMat > tElementsAttachedToFace(1,1);
    moris::Matrix< moris::IdMat >    tChildElemsIdsOnFace;
    moris::Matrix< moris::IndexMat > tChildElemsCMIndOnFace;
    moris::Matrix< moris::IndexMat > tChildElemOnFaceOrdinal;

    moris::uint tElementIndex          = 0;
    moris::uint tPhaseIndex            = 0;
    moris::uint tFaceOrdinal           = 0;
    moris::moris_index tChildMeshIndex = 0;
    moris::moris_id    tElementId      = 0;
    moris::moris_index tMyProcRank     = par_rank();
    bool        tHasChildren = false;


    // get cells and sides in side set
    moris::Cell< mtk::Cell const * > tCellsInSideSet;
    moris::Matrix< moris::IndexMat > tSideSetOrdinals;
    tMeshData.get_sideset_cells_and_ords(aSideSetName, tCellsInSideSet, tSideSetOrdinals );

    // side set data non-intersected
    aElementIdsAndSideOrd(aNoChildIndex)   = Matrix<IndexMat>(tCellsInSideSet.size()*2,2);

    // intersected data
    //TODO: FIGURE OUT MAXIMUM VALUE
    aElementIdsAndSideOrd(aChildIndex) = Matrix<IndexMat>(tCellsInSideSet.size()*2*10,2);

    // keep count
    moris::Cell<moris::uint> tCount(2,0);

    // iterate through sides in set i
    for(moris::uint iSide= 0; iSide<tCellsInSideSet.size(); iSide++)
    {
        tElementIndex = tCellsInSideSet(iSide)->get_index();

        // sides attached to cell
        moris::Matrix<moris::IdMat> tElementFaces = tMeshData.get_entity_connected_to_entity_loc_inds(tElementIndex,EntityRank::ELEMENT,EntityRank::FACE);

        moris_index tSideIndex = tElementFaces(tSideSetOrdinals(iSide));

        if(tMeshData.get_entity_owner(tElementIndex, EntityRank::ELEMENT) == tMyProcRank)
        {
            tHasChildren = mBackgroundMesh.entity_has_children(tElementIndex,EntityRank::ELEMENT);
            // get the faces from the child mesh
            if(tHasChildren)
            {
                tChildMeshIndex = mBackgroundMesh.child_mesh_index(tElementIndex,EntityRank::ELEMENT);

                Child_Mesh const & tChildMesh = mCutMesh.get_child_mesh(tChildMeshIndex);

                tChildMesh.get_child_elements_connected_to_parent_facet(tSideIndex,
                                                                       tChildElemsIdsOnFace,
                                                                       tChildElemsCMIndOnFace,
                                                                       tChildElemOnFaceOrdinal);

                moris::Matrix< moris::IndexMat > const & tChildElementPhaseIndices = tChildMesh.get_element_phase_indices();
                moris::Matrix< moris::IndexMat > const & tChildElementIndices = tChildMesh.get_element_inds();
                moris::Matrix< moris::IndexMat > const & tElementIds = tChildMesh.get_element_ids();

                for(moris::moris_index iCElem  = 0; iCElem < (moris::moris_index)tChildElemsCMIndOnFace.numel(); iCElem++)
                {
                    tPhaseIndex = tChildElementPhaseIndices(0,tChildElemsCMIndOnFace(0,iCElem));
                    if(aOutputOptions.output_phase(tPhaseIndex))
                    {
                        // Child Element Id
                        if(!aOutputIndices)
                        {
                            tElementId = tElementIds(tChildElemsCMIndOnFace(iCElem));
                            tFaceOrdinal   = tChildElemOnFaceOrdinal(iCElem);
                            aElementIdsAndSideOrd(aChildIndex)(tCount(1),0) = tElementId;
                            aElementIdsAndSideOrd(aChildIndex)(tCount(1),1) = tFaceOrdinal;
                            tCount(1)++;
                        }
                        else
                        {
                            tFaceOrdinal   = tChildElemOnFaceOrdinal(iCElem);
                            aElementIdsAndSideOrd(aChildIndex)(tCount(1),0) = tChildElementIndices(tChildElemsCMIndOnFace(iCElem));
                            aElementIdsAndSideOrd(aChildIndex)(tCount(1),1) = tFaceOrdinal;
                            tCount(1)++;
                        }
                    }
                }
            }

            else
            {

                tFaceOrdinal = tSideSetOrdinals(iSide);
                tElementId   = tCellsInSideSet(iSide)->get_id();

                if(aOutputOptions.output_phase(mBackgroundMesh.get_element_phase_index(tElementIndex)))
                {
                    if(!aOutputIndices)
                    {
                        aElementIdsAndSideOrd(aNoChildIndex)(tCount(0),0) = tElementId;
                        aElementIdsAndSideOrd(aNoChildIndex)(tCount(0),1) = tFaceOrdinal;
                        tCount(0)++;
                    }
                    else
                    {
                        aElementIdsAndSideOrd(aNoChildIndex)(tCount(0),0) = tElementIndex;
                        aElementIdsAndSideOrd(aNoChildIndex)(tCount(0),1) = tFaceOrdinal;
                        tCount(0)++;
                    }
                }

            }
        }

    }

    // resize data
    aElementIdsAndSideOrd(aChildIndex).resize(tCount(1),2);
    aElementIdsAndSideOrd(aNoChildIndex).resize(tCount(0),2);

    // Add data to side set info
    // no child
    aSideSetData(aNoChildIndex).mElemIdsAndSideOrds = &aElementIdsAndSideOrd(aNoChildIndex);
    aSideSetData(aNoChildIndex).mSideSetName        = aSideSetName;
    aSideSetData(aNoChildIndex).mSideTopology       = CellTopology::QUAD4;
    aSideSetData(aChildIndex).mElemIdsAndSideOrds   = &aElementIdsAndSideOrd(aChildIndex);
    aSideSetData(aChildIndex).mSideSetName          = aSideSetName + "_i";
    aSideSetData(aChildIndex).mSideTopology         = CellTopology::TRI3;


}

moris::Cell<std::string>
Model::check_for_and_remove_internal_seacas_side_sets(moris::Cell<std::string> & aSideSetNames)
{
    for(std::vector<std::string>::iterator iSet = aSideSetNames.begin(); iSet != aSideSetNames.end(); ++iSet)
    {
        if(iSet->compare("surface_1_quad4") == 0)
        {
            aSideSetNames.data().erase(iSet--);
        }

        else if(iSet->compare("surface_2_quad4") == 0)
        {
            aSideSetNames.data().erase(iSet--);
        }
        else if(iSet->compare("surface_3_quad4") == 0)
        {
            aSideSetNames.data().erase(iSet--);
        }
        else if(iSet->compare("surface_4_quad4") == 0)
        {
            aSideSetNames.data().erase(iSet--);
        }
        else if(iSet->compare("surface_5_quad4") == 0)
        {
            aSideSetNames.data().erase(iSet--);
        }
        else if(iSet->compare("surface_6_quad4") == 0)
        {
            aSideSetNames.data().erase(iSet--);
        }
        else if(iSet->compare("surface_hex8_quad_1") == 0)
        {
            aSideSetNames.data().erase(iSet--);
        }
        else if(iSet->compare("surface_hex8_quad_2") == 0)
        {
            aSideSetNames.data().erase(iSet--);
        }
        else if(iSet->compare("surface_hex8_quad4_1") == 0)
        {
            aSideSetNames.data().erase(iSet--);
        }
        else if(iSet->compare("surface_hex8_quad4_2") == 0)
        {
            aSideSetNames.data().erase(iSet--);
        }
        else if(iSet->compare("surface_hex8_quad4_3") == 0)
        {
            aSideSetNames.data().erase(iSet--);
        }
        else if(iSet->compare("surface_hex8_quad4_4") == 0)
        {
            aSideSetNames.data().erase(iSet--);
        }
    }
    return aSideSetNames;

}

//------------------------------------------------------------------------------

Cell<moris::Matrix<moris::IdMat>>
Model::combine_interface_and_non_interface_blocks(Cell<moris::Matrix<moris::IdMat>> & aChildElementsByPhase,
                                                  Cell<moris::Matrix<moris::IdMat>> & aNoChildElementsByPhase)
{
    moris::uint tNumPhase = aChildElementsByPhase.size();

    Cell<moris::Matrix<moris::IdMat>> tCombinedElementsByPhase(tNumPhase);

    for(moris::uint i =0; i<tNumPhase; i++)
    {
        moris::uint tNumChildElems   = aChildElementsByPhase(i).numel();
        moris::uint tNumNoChildElems = aNoChildElementsByPhase(i).numel();

        tCombinedElementsByPhase(i) = moris::Matrix<moris::IdMat>(1,tNumChildElems + tNumNoChildElems);

        tCombinedElementsByPhase(i)({0,0},{0,tNumChildElems-1}) = aChildElementsByPhase(i).get_row(0);
        tCombinedElementsByPhase(i)({0,0},{tNumChildElems,tNumChildElems + tNumNoChildElems -1}) = aNoChildElementsByPhase(i).get_row(0);
    }

    return tCombinedElementsByPhase;

}

//------------------------------------------------------------------------------

uint
Model::get_num_phases_to_output(Output_Options const & aOutputOptions)
{
    uint tNumPhasesOutput = 0;
    if(aOutputOptions.output_all_phases())
    {
        tNumPhasesOutput = mGeometryEngine->get_num_bulk_phase();
    }
    else
    {
        tNumPhasesOutput = aOutputOptions.num_phases_to_output();
    }

    return tNumPhasesOutput;
}

//------------------------------------------------------------------------------

void
Model::setup_cell_clusters_for_output(moris::mtk::Cell_Cluster_Input & aCellClusterInput,
                                      Output_Options const & aOutputOptions,
                                      moris::Cell<Matrix<IdMat>> & aCellIds)
{
    // iterate through child meshes and construct cells
    uint tNumChildMeshes = mCutMesh.get_num_child_meshes();

    for(moris::uint i = 0; i < tNumChildMeshes; i ++)
    {
        // Get child mesh
        Child_Mesh const & tChildMesh = mCutMesh.get_child_mesh(i);

        // pack the element ids into phase grouping
        Cell<moris::Matrix< moris::IdMat >> tElementIds;
        Cell<moris::Matrix< moris::IdMat >> tCMElementInds;
        tChildMesh.pack_child_mesh_by_phase(mGeometryEngine->get_num_bulk_phase(), tElementIds, tCMElementInds);

        // add them to cell to keep in scope
        aCellIds.push_back(tElementIds(0));
        aCellIds.push_back(tElementIds(1));
    }

    for(moris::uint i = 0; i < tNumChildMeshes; i ++)
    {
        // Get child mesh
        Child_Mesh const & tChildMesh = mCutMesh.get_child_mesh(i);

        // primary index
        moris_index tPrimaryCellIndex = 2*i;
        moris_index tVoidCellIndex    = 2*i+1;

        // parent index
        moris_index tParentCellIndex = tChildMesh.get_parent_element_index();

        // access the parent element from the background mesh
        moris::mtk::Cell* tInterpCell = &mBackgroundMesh.get_mesh_data().get_mtk_cell(tParentCellIndex);

        // add to cluster
        aCellClusterInput.add_cluster_data(tInterpCell,&aCellIds(tPrimaryCellIndex),&aCellIds(tVoidCellIndex),&tChildMesh.get_node_ids(),&tChildMesh.get_parametric_coordinates());

    }
}

void
Model::setup_interface_side_cluster(std::string                      aInterfaceSideLabelBase,
                                    moris::mtk::Side_Cluster_Input & aSideClusterInput,
                                    Output_Options           const & aOutputOptions,
                                    moris::Cell<Matrix<IdMat>>     & aCellIdsandSideOrds,
                                    moris::Cell<Matrix<DDRMat>>    & aParametricCoordinates)
{
    moris::uint tNumPhases = mGeometryEngine->get_num_bulk_phase();

    moris::uint tNumChildMeshes = mCutMesh.get_num_child_meshes();

    for(moris::uint  iP = 0; iP<tNumPhases; iP++)
    {
        // if we are outputting this phase
        //        if(aOutputOptions.output_phase((size_t)iP))
        if(iP == 0)
        {
            // add side set to output
            std::string tSetName = aInterfaceSideLabelBase;

            //iterate through children meshes
            for(moris::uint iC = 0; iC < tNumChildMeshes; iC ++)
            {
                // Get child mesh
                Child_Mesh const & tChildMesh = mCutMesh.get_child_mesh(iC);

                // package this child element by bulk phase
                moris::Matrix< moris::IdMat > tInterfaceElementIdsAndSideOrd = tChildMesh.pack_interface_sides( 0, 0, 1 );

                // add to data which will stay in scope
                aCellIdsandSideOrds.push_back(tInterfaceElementIdsAndSideOrd);

            }
        }
    }

    uint tCount = 0;
    for(moris::uint  iP = 0; iP<tNumPhases; iP++)
    {
        // if we are outputting this phase
        //        if(aOutputOptions.output_phase((size_t)iP))
        if(iP == 0)
        {
            // add side set to output
            std::string tSetName = aInterfaceSideLabelBase;
            moris_index tSideSetOrd = aSideClusterInput.add_side_set_label(tSetName);

            //iterate through children meshes
            for(moris::uint iC = 0; iC < tNumChildMeshes; iC ++)
            {
                // Get child mesh
                Child_Mesh const & tChildMesh = mCutMesh.get_child_mesh(iC);

                // parent cell index
                moris_index tParentCellIndex = tChildMesh.get_parent_element_index();

                // access the parent element from the background mesh
                moris::mtk::Cell* tInterpCell = &mBackgroundMesh.get_mesh_data().get_mtk_cell(tParentCellIndex);

                // add to cluster input data
                //fixme: Add only vertex indices on the interface to cluster. Adding all.
                aSideClusterInput.add_cluster_data(false,tSideSetOrd,tInterpCell,&aCellIdsandSideOrds(tCount),&tChildMesh.get_node_ids(),&tChildMesh.get_parametric_coordinates());

                tCount++;
            }
        }
    }
}
//------------------------------------------------------------------------------


bool
Model::output_node(moris::moris_index aNodeIndex,
                   Output_Options const & aOutputOptions)
{
    bool tIsInterface = mBackgroundMesh.is_interface_node(aNodeIndex,0);
    moris::size_t tPhaseIndex = 0;
    mGeometryEngine->get_phase_index(aNodeIndex,tPhaseIndex);


    if(aOutputOptions.output_phase(tPhaseIndex) && !tIsInterface)
    {
        return true;
    }
    else if(tIsInterface)
    {
        return true;
    }

    return false;
}

//------------------------------------------------------------------------------
moris::size_t
Model::determine_element_phase_index(moris::size_t aRowIndex,
                                     moris::Matrix< moris::IndexMat > const & aElementToNodeIndex)
{
    moris::size_t tNumGeom = mGeometryEngine->get_num_geometries();
    moris::size_t tNumNodesPerElem = aElementToNodeIndex.n_cols();
    moris::Matrix< moris::IndexMat > tNodalPhaseVals(1,tNumGeom,MORIS_INDEX_MAX);


    for(moris::size_t i = 0; i<tNumGeom; i++)
    {
        bool tFoundNonInterfaceNode = false;
        for( moris::size_t j = 0; j<tNumNodesPerElem; j++)
        {
            if(!mBackgroundMesh.is_interface_node(aElementToNodeIndex(aRowIndex,j),i))
            {
                tNodalPhaseVals(0,i) = mGeometryEngine->get_node_phase_index_wrt_a_geometry(aElementToNodeIndex(aRowIndex,j),i);
                tFoundNonInterfaceNode = true;
                break;
            }
        }

        if(!tFoundNonInterfaceNode)
        {
            std::cout<<"Did not find a non-interface node for this element"<<std::endl;
            tNodalPhaseVals(0,i) = 1001;
        }
    }


    moris::moris_index tElemPhaseVal = mGeometryEngine->get_elem_phase_index(tNodalPhaseVals);

    return tElemPhaseVal;
}
void
Model::print_decompsition_preamble(Cell<enum Subdivision_Method> aMethods)
{
    // Only process with rank 0 prints the preamble


    if(moris::par_rank() == 0 && mVerbose)
    {
        std::cout<<"--------------------------------------------------------"<<std::endl;
        std::cout<<"XTK: Specified Decomposition Routines: ";

        for(moris::size_t i = 0 ; i<aMethods.size(); i++)
        {
            std::cout<<"["<<get_enum_str(aMethods(i))<<  "] ";
        }

        std::cout<<std::endl;
    }
}


//------------------------------------------------------------------------------
void
Model::compute_interface_sensitivity_internal()
{
    // Number of geometries in the geometry engine (we need to compute sensitivity wrt each)
    uint tNumGeoms = mGeometryEngine->get_num_geometries();

    // Node coordinates
    moris::Matrix< moris::DDRMat > tNodeCoords = mBackgroundMesh.get_all_node_coordinates_loc_inds();

    for(uint iGeo = 0; iGeo <tNumGeoms; iGeo++)
    {
        // Get interface nodes
        moris::Matrix< moris::IndexMat > tInterfaceNodes = mBackgroundMesh.get_interface_nodes_loc_inds(iGeo);

        // Compute interface sensitivity
        mGeometryEngine->compute_interface_sensitivity(tInterfaceNodes,tNodeCoords,iGeo);
    }
}


void
Model::extract_interface_sensitivity_sparse(moris::Matrix<moris::IndexMat> const & aNodeIndsToOutput,
                                            moris::Cell<moris::Matrix<DDRMat>>   & adxdpData,
                                            moris::Cell<std::string>             & adxdpNames,
                                            moris::Cell<moris::Matrix<DDRMat>>   & aDesVars,
                                            moris::Cell<std::string>             & aDesVarsName,
                                            moris::Matrix<moris::DDRMat>         & aNumDesVars,
                                            std::string                          & aNumDesVarsName) const
{
    // names of sparsely packaged fields
    moris::uint tNumFields = 6;
    adxdpNames = moris::Cell<std::string>({{"dx0dp0"},
        {"dx1dp0"},
        {"dx2dp0"},
        {"dx0dp1"},
        {"dx1dp1"},
        {"dx2dp1"}});



    moris::uint tNumNodes = aNodeIndsToOutput.numel();
    adxdpData = moris::Cell<moris::Matrix<moris::DDRMat>>(tNumFields,moris::Matrix<moris::DDRMat>(tNumNodes,1,0.0));

    //TODO: hardcoded to 2
    tNumFields = 2;
    aDesVarsName = moris::Cell<std::string>({{"DesVar0"},{"DesVar1"}});
    aDesVars = moris::Cell<moris::Matrix<moris::DDRMat>>(tNumFields,moris::Matrix<moris::DDRMat>(tNumNodes,1));

    aNumDesVarsName = "NumDesVar";
    aNumDesVars = moris::Matrix<moris::DDRMat>(1,tNumNodes,0);

    for(moris::uint iNode = 0; iNode<tNumNodes; iNode++)
    {
        moris::moris_index tNodeIndex = aNodeIndsToOutput(iNode);

        if(mBackgroundMesh.is_interface_node(tNodeIndex,0))
        {
            moris::ge::GEN_Geometry_Object const & tNodeGeoObj = mGeometryEngine->get_geometry_object(tNodeIndex);

            moris::Matrix< moris::DDRMat > const & tdxdp = tNodeGeoObj.get_sensitivity_dx_dp();

            MORIS_ASSERT(tdxdp.n_rows() == 2,"Invalid dxdp size for sparse packing, This function only works on tet meshes with discrete fields at the moment");
            MORIS_ASSERT(tdxdp.n_cols() == 3,"Invalid dxdp size for sparse packing, This function only works on tet meshes with discrete fields at the moment");

            adxdpData(0)(iNode) = tdxdp(0,0);
            adxdpData(1)(iNode) = tdxdp(0,1);
            adxdpData(2)(iNode) = tdxdp(0,2);
            adxdpData(3)(iNode) = tdxdp(1,0);
            adxdpData(4)(iNode) = tdxdp(1,1);
            adxdpData(5)(iNode) = tdxdp(1,2);

            moris::Matrix< moris::IndexMat > const & tDesVarInds = tNodeGeoObj.get_node_adv_indices();

            aDesVars(0)(iNode) = (moris::real)mBackgroundMesh.get_mesh_data().get_glb_entity_id_from_entity_loc_index(tDesVarInds(0),EntityRank::NODE);
            aDesVars(1)(iNode) = (moris::real)mBackgroundMesh.get_mesh_data().get_glb_entity_id_from_entity_loc_index(tDesVarInds(1),EntityRank::NODE);

            aNumDesVars(iNode) = tDesVarInds.numel();
        }

    }
}

void
Model::extract_interface_sensitivity_dense(moris::Matrix<moris::IndexMat> const & aNodeIndsToOutput,
                                            moris::Cell<moris::Matrix<DDRMat>>   & adxdpData,
                                            moris::Cell<std::string>             & adxdpNames) const
{
    // names of sparsely packaged fields
    moris::uint tNumDVs = mGeometryEngine->get_num_design_variables();

    // spatial dimension (used often here)
    moris::uint tSpatialDim = this->get_spatial_dim();

    for(moris::uint i = 0; i< tNumDVs; i++)
    {
        for(moris::uint j = 0; j < tSpatialDim; j++)
        {
            adxdpNames.push_back("dx" + std::to_string(j)+"dp"+std::to_string(i));
        }
    }

    moris::uint tNumNodes = aNodeIndsToOutput.numel();
    adxdpData = moris::Cell<moris::Matrix<moris::DDRMat>>(tNumDVs*tSpatialDim,moris::Matrix<moris::DDRMat>(tNumNodes,1,0.0));

    for(moris::uint iNode = 0; iNode<tNumNodes; iNode++)
    {
        moris::moris_index tNodeIndex = aNodeIndsToOutput(iNode);

        if(mBackgroundMesh.is_interface_node(tNodeIndex,0))
        {
            moris::ge::GEN_Geometry_Object const & tNodeGeoObj = mGeometryEngine->get_geometry_object(tNodeIndex);

            moris::Matrix< moris::DDRMat > const & tdxdp = tNodeGeoObj.get_sensitivity_dx_dp();

            MORIS_ASSERT(tdxdp.n_rows() == tNumDVs,"Invalid dxdp size for dense packing");
            MORIS_ASSERT(tdxdp.n_cols() == tSpatialDim,"Invalid dxdp size for dense packing");

            moris::uint tCount = 0;
            for(moris::uint i = 0; i< tNumDVs; i++)
            {
                for(moris::uint j = 0; j < tSpatialDim; j++)
                {
                    adxdpData(tCount)(iNode) = tdxdp(i,j);
                    tCount++;
                }
            }
        }
    }
}


//------------------------------------------------------------------------------

moris::Cell< moris::Matrix < moris::DDRMat > >
Model::assemble_geometry_data_as_mesh_field(moris::Matrix<moris::IndexMat> const & aNodeIndsToOutput)
{
    uint tNumGeometries = mGeometryEngine->get_num_geometries();
    uint tNumNodes      = aNodeIndsToOutput.numel();


    // Allocate output data
    moris::Cell< moris::Matrix < moris::DDRMat > > tGeometryData(tNumGeometries, moris::Matrix<moris::DDRMat>(tNumNodes,1));

    //Iterate through geometries
    for(uint iG = 0; iG <tNumGeometries; iG++)
    {
        // Iterate through nodes
        for(uint iN = 0; iN<tNumNodes; iN++)
        {
            tGeometryData(iG)(iN) = mGeometryEngine->get_entity_phase_val(aNodeIndsToOutput(iN),iG);
        }
    }

    return tGeometryData;
}

//------------------------------------------------------------------------------

moris::Cell<std::string>
Model::assign_geometry_data_names()
{
    uint tNumGeometries = mGeometryEngine->get_num_geometries();

    // base string of geometry data
    std::string tBaseName = "gd_";

    // Allocate output
    moris::Cell<std::string> tGeometryFieldName(tNumGeometries);

    //Iterate through geometries
    for(uint iG = 0; iG <tNumGeometries; iG++)
    {
        tGeometryFieldName(iG) = tBaseName+std::to_string(iG);
    }

    return tGeometryFieldName;
}

//------------------------------------------------------------------------------

moris::Cell < enum moris::EntityRank >
Model::assign_geometry_data_field_ranks()
{
    uint tNumGeometries = mGeometryEngine->get_num_geometries();

    // base string of geometry data
    std::string tBaseName = "gd_";

    // Allocate output
    // Note: for now this is a nodal field always
    moris::Cell<enum moris::EntityRank> tGeometryFieldRank(tNumGeometries,moris::EntityRank::NODE);

    return tGeometryFieldRank;
}
//------------------------------------------------------------------------------



}
