/*
 * cl_XTK_Model.hpp
 *
 *  Created on: Jul 2, 2017
 *      Author: ktdoble
 */

#ifndef SRC_XTK_CL_XTK_MODEL_HPP_
#define SRC_XTK_CL_XTK_MODEL_HPP_

// Standard Include
#include <limits>
#include <unordered_map>
#include <mpi.h>
#include <ctime>

// XTKL: Mesh Includes
#include "cl_MTK_Mesh_Core.hpp"
#include "cl_MTK_Interpolation_Mesh.hpp"

// XTKL: Container includes
#include "cl_Cell.hpp"

// XTKL: Tools includes
#include "cl_XTK_Enums.hpp"
#include "cl_XTK_Cut_Integration_Mesh.hpp"
#include "cl_XTK_Cut_Mesh.hpp"
#include "cl_XTK_Enrichment.hpp"
#include "cl_XTK_Ghost_Stabilization.hpp"
#include "cl_XTK_Background_Mesh.hpp"
#include "cl_XTK_Decomposition_Data.hpp"
#include "cl_Interpolaton.hpp"

#include "cl_XTK_Output_Options.hpp"

// Linalg Includes
#include "cl_Matrix.hpp"
#include "typedefs.hpp"
#include "fn_print.hpp"

#include "cl_Communication_Tools.hpp"

// Topology
//TODO: MOVE THESE WITH CUTTING METHODS SOMEWHERE ELSE
#include "cl_XTK_Topology.hpp"
#include "cl_XTK_Edge_Topology.hpp"
#include "cl_XTK_Quad_4_Topology.hpp"
#include "cl_XTK_Hexahedron_8_Topology.hpp"

// general geometry engine class
#include "cl_GEN_Geometry_Engine.hpp"

#include "cl_Param_List.hpp"

#include "cl_MTK_Intersection_Detect.hpp"
#include "cl_MTK_Intersection_Detect_2D.hpp"

namespace xtk
{
class Integration_Mesh_Generator;
class Enrichment;
class Enrichment_Parameters;
class Ghost_Stabilization;
class Enriched_Interpolation_Mesh;
class Enriched_Integration_Mesh;
class Ghost_Stabilization;
class Multigrid;
}// namespace xtk


using namespace moris;

namespace xtk
{
class Model
{
  public:
    // Public member functions/data
    bool mVerbose = false;
    moris::uint mVerboseLevel = 0;
    // friend classes
    friend class Integration_Mesh_Generator;
    friend class Cut_Integration_Mesh;
    friend class Enrichment;
    friend class Enriched_Interpolation_Mesh;
    friend class Enriched_Integration_Mesh;
    friend class Cut_Integration_Mesh;
    friend class Ghost_Stabilization;
    friend class Multigrid;

    //--------------------------------------------------------------------------------
    // Initialization
    //--------------------------------------------------------------------------------
    Model(){};

    //--------------------------------------------------------------------------------
    /**
   * @brief Constructor for standalone use of XTK
   *
   */
    Model( uint                         aModelDimension,
        moris::mtk::Interpolation_Mesh* aMeshData,
        moris::ge::Geometry_Engine*     aGeometryEngine,
        bool                            aLinkGeometryOnConstruction = true );

    //--------------------------------------------------------------------------------

    /**
   * @brief Parameter list based XTK initialization
   */
    Model( moris::ParameterList const& aParameterList );

    //--------------------------------------------------------------------------------

    ~Model();


    //--------------------------------------------------------------------------------
    /*
   * @brief Set the pointer to the geometry engine
   */
    void
    set_geometry_engine( moris::ge::Geometry_Engine* aGeometryEngine );

    //--------------------------------------------------------------------------------
    /**
   * @brief Set the pointer to the background interpolation mesh
   */
    void
    set_mtk_background_mesh( moris::mtk::Interpolation_Mesh* aMesh );

    //--------------------------------------------------------------------------------

    /**
   * @brief Set the pointer to the output performer
   */
    void
    set_output_performer( std::shared_ptr< mtk::Mesh_Manager > aMTKPerformer );

    //--------------------------------------------------------------------------------

    /**
   * @brief Set the pointer to the input performer
   */
    void
    set_input_performer( std::shared_ptr< mtk::Mesh_Manager > aMTKPerformer );


    void
    set_cut_ig_mesh( std::shared_ptr< Cut_Integration_Mesh > aCutIgMesh );

    //--------------------------------------------------------------------------------

    /**
   * @brief Initialize data using the interpolation mesh
   */
    void
    initialize( moris::mtk::Interpolation_Mesh* aMesh );


    //--------------------------------------------------------------------------------
    // Operations
    //--------------------------------------------------------------------------------
    /**
   * @brief Workflow based perform call. Uses parameter list to decide operations to
   * perform
   */
    bool
    perform();

    //--------------------------------------------------------------------------------
    /**
   * Decomposes a mesh to conform to a geometric interface
   * @param[in] aMethods Subdivision Methods
   */
    bool decompose( Cell< enum Subdivision_Method > aMethods );

    //--------------------------------------------------------------------------------
    /**
   * Remove child meshes that have an intersection but all cells are in the same bulk phase
   */
    void
    set_cleanup_cut_mesh_flag( bool aCleanupMesh );


    //--------------------------------------------------------------------------------
    /**
   * @brief get the memory usage of XTK
   */
    moris::Memory_Map
    get_memory_usage();

    // ----------------------------------------------------------------------------------

    /**
   * Perform generalized heaviside  enrichment.
   */
    void
    perform_basis_enrichment(
        enum EntityRank const& aBasisRank,
        moris_index const&     aMeshIndex = 0 );

    // ----------------------------------------------------------------------------------

    /**
   * @brief Perform basis enrichment
   * @param[in] aBasisRank Entity rank of basis functions
   * @param[in] aMeshIndex Mesh index
   */
    void
    perform_basis_enrichment(
        enum EntityRank const&    aBasisRank,
        Matrix< IndexMat > const& aMeshIndex );

    // ----------------------------------------------------------------------------------

    void perform_hanging_node_identification();

    // ----------------------------------------------------------------------------------
    /**
   * @brief Probes and prints the information about a background cell
   */
    void
    probe_bg_cell( Matrix< IndexMat > const& tBGCellIds );

    // ----------------------------------------------------------------------------------

    Cut_Integration_Mesh*
    get_cut_integration_mesh();
    /**
   * @return Basis enrichment
   */
    Enrichment const&
    get_basis_enrichment();

    // ----------------------------------------------------------------------------------
    /**
   * @param[in] aIndex Interpolation Index
   * @return Enriched interpolation mesh
   */
    Enriched_Interpolation_Mesh&
    get_enriched_interp_mesh( moris::moris_index aIndex = 0 );

    // ----------------------------------------------------------------------------------
    /**
   * @param[in] aIndex Interpolation Index
   * @return Enriched integration mesh
   */
    Enriched_Integration_Mesh&
    get_enriched_integ_mesh( moris::moris_index aIndex = 0 );


    // ----------------------------------------------------------------------------------
    /**
   * @brief Constructs the face oriented ghost penalization
   */
    void
    construct_face_oriented_ghost_penalization_cells();

    // ----------------------------------------------------------------------------------
    /**
   * @return Ghost stabilization
   */
    Ghost_Stabilization&
    get_ghost_stabilization( moris::moris_index aIndex = 0 );

    // ----------------------------------------------------------------------------------

    /**
   * @brief Construct multigrid information
   */
    void construct_multigrid();

    //--------------------------------------------------------------------------------
    // Member data access functions
    //--------------------------------------------------------------------------------
    /**
   * @return Cut mesh (collection of child meshes)
   */
    Cut_Mesh&
    get_cut_mesh();

    // ----------------------------------------------------------------------------------
    /**
   * @return Const cut mesh (collection of child meshes)
   */
    Cut_Mesh const&
    get_cut_mesh() const;

    // ----------------------------------------------------------------------------------
    /**
   * @return Background mesh
   */
    moris::mtk::Interpolation_Mesh&
    get_background_mesh();

    // ----------------------------------------------------------------------------------
    /**
   * @return Constant background mesh
   */
    moris::mtk::Interpolation_Mesh const&
    get_background_mesh() const;

    // ----------------------------------------------------------------------------------
    /**
   * @return Geometry engine pointer
   */
    moris::ge::Geometry_Engine*
    get_geom_engine();

    // ----------------------------------------------------------------------------------
    // Outputting functions
    //-----------------------------------------------------------------------------------
    /**
   * @brief Outputs the Mesh to a mesh data which can then be written to exodus files as desired.
   * @param[in] aOutputOptions Mesh output options
   */
    moris::mtk::Integration_Mesh*
    get_output_mesh( Output_Options const& aOutputOptions = Output_Options() );


    /**
   * get spatial dimension of model
   */
    moris::uint
    get_spatial_dim() const;

    //-----------------------------------------------------------------------------------

    /**
   * @returns element to element connectivity
   */
    moris::Cell< moris::Cell< moris::mtk::Cell* > > const&
    get_element_to_element()
    {
        return mElementToElement;
    };

    //-----------------------------------------------------------------------------------

    /**
   * @returns element to element my side ordinals
   */
    moris::Cell< moris::Cell< moris_index > > const&
    get_element_to_element_my_side_ords()
    {
        return mElementToElementSideOrds;
    };

    //-----------------------------------------------------------------------------------

    /**
   * @returns element to element neighbor side ordinals
   */
    moris::Cell< moris::Cell< moris_index > > const&
    get_element_to_element_neighbor_side_ords()
    {
        return mElementToElementNeighborSideOrds;
    };

    //-----------------------------------------------------------------------------------

    moris_index
    get_cell_xtk_index( moris_id aCellId );

    //-----------------------------------------------------------------------------------

    moris::Cell< moris::Cell< moris_index > > const&
    get_subphase_to_subphase();

    //-----------------------------------------------------------------------------------

    moris::Cell< moris::Cell< moris_index > > const&
    get_subphase_to_subphase_my_side_ords();

    //-----------------------------------------------------------------------------------

    moris::Cell< moris::Cell< moris_index > > const&
    get_subphase_to_subphase_neighbor_side_ords();

    //-----------------------------------------------------------------------------------
    moris::Matrix< moris::IndexMat >
    get_num_subphase_neighbors();

    //-----------------------------------------------------------------------------------

    moris::Cell< moris::Cell< moris_index > > const&
    get_subphase_to_subphase_transition_loc();

    //-----------------------------------------------------------------------------------

    bool
    subphase_is_in_child_mesh( moris_index aSubphaseIndex );

    //-----------------------------------------------------------------------------------

    moris::Matrix< moris::IndexMat >
    get_element_to_subphase();

    //-----------------------------------------------------------------------------------

    moris_id
    get_subphase_id( moris_id aSubphaseIndex );

    //-----------------------------------------------------------------------------------

    moris_index
    get_subphase_index( moris_id aSubphaseId );

    //-----------------------------------------------------------------------------------

    std::shared_ptr< Multigrid > get_multigrid_ptr();

    //--------------------------------------------------------------------------------
    // Printing Functions
    //--------------------------------------------------------------------------------
    void
    print_cells();
    //------------------------------------------------------------------------------

    moris::ParameterList&
    get_parameter_list();

    void
    setup_diagnostics(
        bool               aDiagnostics,
        std::string const& aDiagnosticPath,
        std::string const& aDiagnosticLabel );

    std::string
    get_diagnostic_file_name( std::string const& aLabel ) const;

      /**
     * @brief Tells the integration mesh to either decompose all background cells or not
     * 
     * @return true - decomposition mesh will triangulate/tessalate all cells
     * @return false 
     */
    bool
    triangulate_all();

    enum CellTopology
    get_parent_cell_topology() const
    {
        mtk::Cell&                     tCell     = mBackgroundMesh->get_mtk_cell( 0 );
        enum moris::mtk::Geometry_Type tGeomType = tCell.get_geometry_type();

        enum CellTopology tTopo = CellTopology::INVALID;
        if ( tGeomType == moris::mtk::Geometry_Type::HEX )
        {
            tTopo = CellTopology::HEX8;
        }
        else if ( tGeomType == moris::mtk::Geometry_Type::TET )
        {
            tTopo = CellTopology::TET4;
        }
        else if ( tGeomType == moris::mtk::Geometry_Type::QUAD )
        {
            tTopo = CellTopology::QUAD4;
        }
        else
        {
            MORIS_ERROR( 0, "Provided parent cell topology not implemented." );
        }


        return tTopo;
    }

    Matrix< IdMat >
    get_communication_table() const
    {
        return mCutIntegrationMesh->get_communication_table();
    }


  protected:
    moris::ParameterList mParameterList;

    uint                                    mModelDimension;
    moris::mtk::Interpolation_Mesh*         mBackgroundMesh;
    Cut_Mesh                                mCutMesh;
    moris::ge::Geometry_Engine*             mGeometryEngine;
    std::shared_ptr< Cut_Integration_Mesh > mCutIntegrationMesh;
    Enrichment*                             mEnrichment;
    Ghost_Stabilization*                    mGhostStabilization;
    Cell< Enriched_Interpolation_Mesh* >    mEnrichedInterpMesh;
    Cell< Enriched_Integration_Mesh* >      mEnrichedIntegMesh;
    std::shared_ptr< xtk::Multigrid >       mMultigrid;

  private:
    // XTK Model State Flags
    bool mDecomposed        = false;// Model has been decomposed
    bool mConvertedToTet10s = false;// Model has been converted from tet4's to tet10s
    bool mMeshDataFinalized = false;
    bool mEnriched          = false;// Model has been enriched
    bool mUnzipped          = false;// Model has been unzipped
    bool mGhost             = false;// Model has setup ghost stabilization

    // Flag to cleanup mesh at end of decomposition
    bool mTriangulateAll = false;// Triangulate all background cells
    bool mCleanupMesh    = false;// Cleanup the mesh

    // diagnostic information
    bool        mDiagnostics    = false;
    std::string mDiagnosticPath = "";
    std::string mDiagnosticId   = "";

    // communication table
    moris::Matrix<IdMat> mCommunicationMap;


    // cell map
    std::map< moris_id, moris_index > mCellGlbToLocalMap;


    // element to element neighborhood
    moris::Cell< moris::Cell< moris::mtk::Cell* > > mElementToElement;
    moris::Cell< moris::Cell< moris_index > >       mElementToElementSideOrds;
    moris::Cell< moris::Cell< moris_index > >       mElementToElementNeighborSideOrds;
    moris::Cell< moris::Cell< moris_index > >       mSubphaseToSubPhase;
    moris::Cell< moris::Cell< moris_index > >       mSubphaseToSubPhaseMySideOrds;
    moris::Cell< moris::Cell< moris_index > >       mSubphaseToSubPhaseNeighborSideOrds;

    // in the case of a hierarchically refined mesh, there are transitions with hanging nodes
    // this data flags the transition from a large facet to a smaller facet. (this is trivial
    // for non hmr type meshes)
    moris::Cell< moris::Cell< moris_index > > mTransitionNeighborCellLocation;

    // local to global subphase map
    std::unordered_map< moris::moris_id, moris::moris_index > mGlobalToLocalSubphaseMap;

    std::shared_ptr< mtk::Mesh_Manager > mMTKInputPerformer  = nullptr;
    std::shared_ptr< mtk::Mesh_Manager > mMTKOutputPerformer = nullptr;

    bool mInitializeCalled = false;

    // Intersection algorithm data members
    // They contain added clusters,cells, vertices
    moris::mtk::Intersection_Detect*    mIntersectionDetect   = nullptr;
    moris::mtk::Intersection_Detect_2D* mIntersectionDetect2D = nullptr;

    // Private Functions
  private:
    //------------------------------------------------------------------------------

    // Functions that take parameter inputs and setup XTK inputs
    bool
    has_parameter_list();


    //------------------------------------------------------------------------------

    /**
   * Verifys provided parameter list
   */
    bool
    valid_parameters();


    //------------------------------------------------------------------------------

    /*
   * Using the parameter list, figure out the cell of subdivision methods
   */
    Cell< enum Subdivision_Method >
    get_subdivision_methods();


    //------------------------------------------------------------------------------

    void
    create_new_node_association_with_geometry( Decomposition_Data& tDecompData );

    //------------------------------------------------------------------------------

    /**
    * Verifies that the nodes in a decomposition have all been assigned node ids
    */
    bool
    verify_successful_node_assignment( Decomposition_Data& aDecompData );

    //------------------------------------------------------------------------------

  protected:
    void
    send_outward_requests(
        moris_index const&          aMPITag,
        Cell< uint > const&         aProcRanks,
        Cell< Matrix< IndexMat > >& aOutwardRequests );

    //------------------------------------------------------------------------------

    void
    inward_receive_requests(
        moris_index const&          aMPITag,
        moris::uint                 aNumRows,
        Cell< Matrix< IndexMat > >& aReceivedData,
        Cell< uint >&               aProcRanksReceivedFrom );

    //------------------------------------------------------------------------------

    void
    prepare_request_answers(
        Decomposition_Data&               aDecompData,
        Cell< Matrix< IndexMat > > const& aReceiveData,
        Cell< Matrix< IndexMat > >&       aReceivedRequestAnswers );

    //------------------------------------------------------------------------------

    void
    return_request_answers(
        moris_index const&                aMPITag,
        Cell< Matrix< IndexMat > > const& aRequestAnswers,
        Cell< uint > const&               aProcRanks );

    //------------------------------------------------------------------------------

    void
    inward_receive_request_answers(
        moris_index const&          aMPITag,
        moris::uint const&          aNumRows,
        Cell< uint > const&         aProcRanks,
        Cell< Matrix< IndexMat > >& aReceivedData );

    //------------------------------------------------------------------------------

    void
    handle_received_request_answers(
        Decomposition_Data&               aDecompData,
        Cell< Matrix< IndexMat > > const& aRequests,
        Cell< Matrix< IndexMat > > const& aRequestAnswers,
        moris::moris_id&                  aNodeId );

    //------------------------------------------------------------------------------
    // moris real versions of above
    void
    send_outward_requests_reals(
        moris_index const&        aMPITag,
        Cell< uint > const&       aProcRanks,
        Cell< Matrix< DDRMat > >& aOutwardRequests );

    //------------------------------------------------------------------------------

    void
    inward_receive_requests_reals(
        moris_index const&        aMPITag,
        moris::uint               aNumRows,
        Cell< Matrix< DDRMat > >& aReceivedData,
        Cell< uint >&             aProcRanksReceivedFrom );

    //------------------------------------------------------------------------------

    void
    return_request_answers_reals(
        moris_index const&              aMPITag,
        Cell< Matrix< DDRMat > > const& aRequestAnswers,
        Cell< uint > const&             aProcRanks );
    //------------------------------------------------------------------------------

    /*
               * For nodes that are created during the decomposition process, tell
               * the XTK mesh about where they live in child meshes.
               */
    void
    associate_nodes_created_during_decomp_to_child_meshes();

    //------------------------------------------------------------------------------

    /*
               * Set element phase index
               */
    void
    set_element_phases();

    //------------------------------------------------------------------------------

    /*
               * Tells the XTK mesh about where it's children live in the cut mesh
               */
    void
    set_downward_inheritance();

    //------------------------------------------------------------------------------
    /*
               * This algorithm sets up the active child mesh indices and registers new pairs in the downward inheritance
               */

    void
    run_first_cut_routine(
        moris::uint                       aGeomIndex,
        moris::Matrix< moris::IndexMat >& aActiveChildMeshIndices,
        moris::Matrix< moris::IndexMat >& aNewPairBool );

    //------------------------------------------------------------------------------
    /*
               * Returns a flag that all intersected cells are on the same level if false the decomposition will return false
               */
    bool
    all_child_meshes_on_same_level();

    //------------------------------------------------------------------------------

    /*!
               * Constructs the output mesh using provided Output_Options
               */
    moris::mtk::Integration_Mesh*
    construct_output_mesh( Output_Options const& aOutputOptions );

    //------------------------------------------------------------------------------

    /*!
               * Setup interface single sided side sets
               */
    void
    setup_interface_single_side_sets(
        Output_Options const&                  aOutputOptions,
        Cell< moris::Matrix< moris::IdMat > >& aCellIdsAndSideOrds,
        Cell< std::string >&                   aInterfaceSetNames );

    //------------------------------------------------------------------------------

    //------------------------------------------------------------------------------

    moris::Cell< std::string >
    assign_geometry_data_names();

    //------------------------------------------------------------------------------


    moris::Cell< enum moris::EntityRank >
    assign_geometry_data_field_ranks();

    //------------------------------------------------------------------------------

    //------------------------------------------------------------------------------

    /*!
               * Sets up background side sets for mesh output. Propogates the side set from
               * the background mesh to the output mesh
               */
    moris::Cell< moris::mtk::MtkSideSetInfo >
    propogate_background_side_sets(
        moris::Cell< moris::Matrix< IndexMat > >& aSideSetData,
        Output_Options const&                     aOutputOptions );

    //------------------------------------------------------------------------------

    /*!
               * Add a single side set from the background mesh
               */
    void
    propogate_background_side_set(
        std::string const&                         aSideSetName,
        moris::moris_index                         aNoChildIndex,
        moris::moris_index                         aChildIndex,
        moris::Cell< moris::Matrix< IndexMat > >&  aElementIdsAndSideOrd,
        moris::Cell< moris::mtk::MtkSideSetInfo >& aSideSetData,
        Output_Options const&                      aOutputOptions,
        bool                                       aOutputIndices );

    //------------------------------------------------------------------------------

    /*!
               * This function checks for the side sets which appear in a mesh that comes
               * from a SEACAS generated string. For example:
               * "generated:1x1x1|sideset:xXyYzZ"
               */
    moris::Cell< std::string >
    check_for_and_remove_internal_seacas_side_sets( moris::Cell< std::string >& aSideSetNames );

    //------------------------------------------------------------------------------

    /*!
               * Combine interface and non-interface blocks
               */
    Cell< moris::Matrix< moris::IdMat > >
    combine_interface_and_non_interface_blocks(
        Cell< moris::Matrix< moris::IdMat > >& tChildElementsByPhase,
        Cell< moris::Matrix< moris::IdMat > >& tNoChildElementsByPhase );

    //------------------------------------------------------------------------------

    /*
               * Returns the number of phases to output given the output options
               */
    uint
    get_num_phases_to_output( Output_Options const& aOutputOptions );

    //------------------------------------------------------------------------------

    /*!
               * Setup clustering data
               */
    void
    setup_cell_clusters_for_output(
        moris::mtk::Cell_Cluster_Input& aCellClusterInput,
        Output_Options const&           aOutputOptions,
        moris::Cell< Matrix< IdMat > >& tCellIds );

    //------------------------------------------------------------------------------

    void
    setup_interface_side_cluster(
        std::string                      aInterfaceSideLabelBase,
        moris::mtk::Side_Cluster_Input&  aCellClusterInput,
        Output_Options const&            aOutputOptions,
        moris::Cell< Matrix< IdMat > >&  tCellIdsandSideOrds,
        moris::Cell< Matrix< DDRMat > >& tParametricCoordinates );

    //------------------------------------------------------------------------------

    bool
    output_node(
        moris::moris_index    aNodeIndex,
        Output_Options const& aOutputOptions );

    //------------------------------------------------------------------------------

    /*
               * Prints the method of decomposition, type of background mesh,
               */
    void
    print_decompsition_preamble( Cell< enum Subdivision_Method > aMethods );

    //------------------------------------------------------------------------------

    moris::size_t
    determine_element_phase_index(
        moris::size_t                           aRowIndex,
        moris::Matrix< moris::IndexMat > const& aElementToNodeIndex );

    //------------------------------------------------------------------------------

    void
    collect_subphases_attached_to_facet_on_cell(
        moris::moris_index          aCellIndex,
        moris::moris_index          aFacetOrdinal,
        Cell< moris::moris_index >& aCellSubphaseIndices,
        Cell< moris::moris_index >& aCellSubphaseBulkIndices,
        Cell< moris::moris_index >& aRepresentativeCellInd,
        Cell< moris::moris_index >& aRepresentativeCellSideOrdinal );


    //------------------------------------------------------------------------------

    void
    inward_receive_request_answers_reals(
        moris_index const&        aMPITag,
        moris::uint const&        aNumRows,
        Cell< uint > const&       aProcRanks,
        Cell< Matrix< DDRMat > >& aReceivedData );

    //------------------------------------------------------------------------------


  private:
    //------------------------------------------------------------------------------
    // Enrichment computation functions
    //------------------------------------------------------------------------------

    void perform_basis_enrichment_internal( enum EntityRank const& aBasisRank, Matrix< IndexMat > const& aMeshIndex );
};
}// namespace xtk

#endif /* SRC_XTK_CL_XTK_MODEL_HPP_ */
