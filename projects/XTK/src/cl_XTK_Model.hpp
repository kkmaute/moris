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
}


using namespace moris;

namespace xtk
{
    class Model
    {
        public:
            // Public member functions/data
            bool mVerbose = false;
            bool mDiagnostics = true;

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
            /*!
             * @brief Constructor for standalone use of XTK
             *
             */
            Model(uint aModelDimension,
                    moris::mtk::Interpolation_Mesh* aMeshData,
                    moris::ge::Geometry_Engine* aGeometryEngine,
                    bool aLinkGeometryOnConstruction = true);

            //--------------------------------------------------------------------------------

            /*!
             * @brief Parameter list based XTK initialization
             */
            Model(moris::ParameterList const & aParameterList);

            //--------------------------------------------------------------------------------

            ~Model();

            //--------------------------------------------------------------------------------
            /*!
             * @brief Set the pointer to the geometry engine
             */
            void
            set_geometry_engine(moris::ge::Geometry_Engine* aGeometryEngine);

            //--------------------------------------------------------------------------------
            /*!
             * @brief Set the pointer to the background interpolation mesh
             */
            void
            set_mtk_background_mesh(moris::mtk::Interpolation_Mesh* aMesh);

            //--------------------------------------------------------------------------------

            /*!
             * @brief Set the pointer to the output performer
             */
            void
            set_output_performer( std::shared_ptr< mtk::Mesh_Manager > aMTKPerformer );

            //--------------------------------------------------------------------------------

            /*!
             * @brief Set the pointer to the input performer
             */
            void
            set_input_performer( std::shared_ptr< mtk::Mesh_Manager > aMTKPerformer );

            //--------------------------------------------------------------------------------

            /*!
             * @brief Initialize data using the interpolation mesh
             */
            void
            initialize( moris::mtk::Interpolation_Mesh* aMesh );


            //--------------------------------------------------------------------------------
            // Operations
            //--------------------------------------------------------------------------------
            /*!
             * @brief Workflow based perform call. Uses parameter list to decide operations to
             * perform
             */
            bool
            perform();

            //--------------------------------------------------------------------------------
            /*!
             * Decomposes a mesh to conform to a geometric interface
             * @param[in] aMethods Subdivision Methods
             */
            bool decompose(Cell<enum Subdivision_Method> aMethods);

            //--------------------------------------------------------------------------------
            /*!
            * Remove child meshes that have an intersection but all cells are in the same bulk phase
            */
            void
            set_cleanup_cut_mesh_flag( bool aCleanupMesh );


            //--------------------------------------------------------------------------------
            /*!
            * @brief get the memory usage of XTK
            */
            moris::Memory_Map
            get_memory_usage();

            //--------------------------------------------------------------------------------
            /*!
             * Unzipps the interface (must be called after decompose but prior
             * to enrichment)
             */
            void
            unzip_interface();

            // ----------------------------------------------------------------------------------
 
            /*!
             * Perform generalized heaviside  enrichment.
             */
            void
            perform_basis_enrichment(
                    enum EntityRank  const & aBasisRank,
                    moris_index      const & aMeshIndex = 0);

            // ----------------------------------------------------------------------------------
 
            /*!
             * @brief Perform basis enrichment
             * @param[in] aBasisRank Entity rank of basis functions
             * @param[in] aMeshIndex Mesh index
             */
            void
            perform_basis_enrichment(
                    enum EntityRank  const & aBasisRank,
                    Matrix<IndexMat> const & aMeshIndex);

            // ----------------------------------------------------------------------------------

            void perform_hanging_node_identification();

            // ----------------------------------------------------------------------------------
            /*!
             * @brief Probes and prints the information about a background cell
             */
            void
            probe_bg_cell(Matrix<IndexMat> const & tBGCellIds);

            // ----------------------------------------------------------------------------------

            Cut_Integration_Mesh*
            get_cut_integration_mesh();
            /*!
             * @return Basis enrichment
             */
            Enrichment const &
            get_basis_enrichment();

            // ----------------------------------------------------------------------------------
            /*!
             * @param[in] aIndex Interpolation Index
             * @return Enriched interpolation mesh
             */
            Enriched_Interpolation_Mesh &
            get_enriched_interp_mesh(moris::moris_index aIndex = 0);

            // ----------------------------------------------------------------------------------
            /*!
             * @param[in] aIndex Interpolation Index
             * @return Enriched integration mesh
             */
            Enriched_Integration_Mesh &
            get_enriched_integ_mesh(moris::moris_index  aIndex = 0);


            // ----------------------------------------------------------------------------------
            /*!
             * @brief Constructs the face oriented ghost penalization
             */
            void
            construct_face_oriented_ghost_penalization_cells();

            // ----------------------------------------------------------------------------------
            /*!
             * @return Ghost stabilization
             */
            Ghost_Stabilization &
            get_ghost_stabilization(moris::moris_index  aIndex = 0);

            // ----------------------------------------------------------------------------------
            //FIXME: REMOVE
            /*!
             * Convert Tet4 elements to Tet10 Elements
             */
            void
            convert_mesh_tet4_to_tet10();

            // ----------------------------------------------------------------------------------

            /*!
             * @brief Construct multigrid information
             */
            void construct_multigrid();

            //--------------------------------------------------------------------------------
            // Member data access functions
            //--------------------------------------------------------------------------------
            /*!
             * @return Cut mesh (collection of child meshes)
             */
            Cut_Mesh &
            get_cut_mesh();

            // ----------------------------------------------------------------------------------
            /*!
             * @return Const cut mesh (collection of child meshes)
             */
            Cut_Mesh const &
            get_cut_mesh() const;

            // ----------------------------------------------------------------------------------
            /*!
             * @return Background mesh
             */
            Background_Mesh &
            get_background_mesh();

            // ----------------------------------------------------------------------------------
            /*!
             * @return Constant background mesh
             */
            Background_Mesh const &
            get_background_mesh() const;

            // ----------------------------------------------------------------------------------
            /*!
             * @return Geometry engine pointer
             */
            moris::ge::Geometry_Engine*
            get_geom_engine();

            // ----------------------------------------------------------------------------------
            // Outputting functions
            //-----------------------------------------------------------------------------------
            /*!
             * @brief Outputs the Mesh to a mesh data which can then be written to exodus files as desired.
             * @param[in] aOutputOptions Mesh output options
             */
            moris::mtk::Integration_Mesh*
            get_output_mesh(Output_Options const & aOutputOptions = Output_Options());

            //--------------------------------------------------------------------------------
            // Cell Neighborhood creation and access
            //--------------------------------------------------------------------------------

            //-----------------------------------------------------------------------------------


            //-----------------------------------------------------------------------------------
            /*!
             * @briefHanging node considerations for cell neighborhood
             */
            void
            construct_complex_neighborhood();

            //-----------------------------------------------------------------------------------
            // Data access functions
            //--------------------------------------------------------------------------------

            /*!
             * get spatial dimension of model
             */
            moris::uint
            get_spatial_dim() const;

            //-----------------------------------------------------------------------------------

            /*!
             * @returns element to element connectivity
             */
            moris::Cell<moris::Cell<moris::mtk::Cell*>> const &
            get_element_to_element(){ return mElementToElement; };

             /*!
             * @returns element to element my side ordinals
             */
            moris::Cell<moris::Cell<moris_index>> const &
            get_element_to_element_my_side_ords(){ return mElementToElementSideOrds; };

            /*!
             * @returns element to element neighbor side ordinals
             */
            moris::Cell<moris::Cell<moris_index>> const &
            get_element_to_element_neighbor_side_ords(){ return mElementToElementNeighborSideOrds; };

            //-----------------------------------------------------------------------------------

            moris_index
            get_cell_xtk_index(moris_id aCellId);

            //-----------------------------------------------------------------------------------

            moris::Cell<moris::Cell<moris_index>>  const &
            get_subphase_to_subphase();

            //-----------------------------------------------------------------------------------

            moris::Cell<moris::Cell<moris_index>>  const &
            get_subphase_to_subphase_my_side_ords();

            //-----------------------------------------------------------------------------------

            moris::Cell<moris::Cell<moris_index>>  const &
            get_subphase_to_subphase_neighbor_side_ords();

            //-----------------------------------------------------------------------------------

            moris::Cell<moris::Cell<moris_index>>  const &
            get_subphase_to_subphase_transition_loc();

            //-----------------------------------------------------------------------------------

            bool
            subphase_is_in_child_mesh(moris_index aSubphaseIndex);

            //-----------------------------------------------------------------------------------

            moris::Matrix<moris::IndexMat>
            get_element_to_subphase();

            //-----------------------------------------------------------------------------------

            moris_id
            get_subphase_id(moris_id aSubphaseIndex);

            //-----------------------------------------------------------------------------------

            moris_index
            get_subphase_index(moris_id aSubphaseId);

            //-----------------------------------------------------------------------------------

            std::shared_ptr< Multigrid > get_multigrid_ptr();

            //--------------------------------------------------------------------------------
            // Printing Functions
            //--------------------------------------------------------------------------------
            void
            print_cells();

            //------------------------------------------------------------------------------

        protected:
            moris::ParameterList               mParameterList;

            uint                               mModelDimension;
            Background_Mesh                    mBackgroundMesh;
            Cut_Mesh                           mCutMesh;

            moris::ge::Geometry_Engine*        mGeometryEngine;

            std::shared_ptr<Cut_Integration_Mesh> mCutIntegrationMesh;
            Enrichment*                           mEnrichment;
            Ghost_Stabilization*                  mGhostStabilization;
            Cell<Enriched_Interpolation_Mesh*>    mEnrichedInterpMesh;
            Cell<Enriched_Integration_Mesh*>      mEnrichedIntegMesh;

            std::shared_ptr< xtk::Multigrid >  mMultigrid;

        private:

            // XTK Model State Flags
            bool mDecomposed        = false; // Model has been decomposed
            bool mConvertedToTet10s = false; // Model has been converted from tet4's to tet10s
            bool mMeshDataFinalized = false;
            bool mEnriched          = false; // Model has been enriched
            bool mUnzipped          = false; // Model has been unzipped
            bool mGhost             = false; // Model has setup ghost stabilization

            // Flag to cleanup mesh at end of decomposition
            bool mTriangulateAll = false; // Triangulate all background cells
            bool mCleanupMesh    = false; // Cleanup the mesh 

            // cell map
            std::map< moris_id, moris_index> mCellGlbToLocalMap;


            // element to element neighborhood
            moris::Cell<moris::Cell<moris::mtk::Cell*>> mElementToElement;
            moris::Cell<moris::Cell<moris_index>>       mElementToElementSideOrds;
            moris::Cell<moris::Cell<moris_index>>       mElementToElementNeighborSideOrds;
            moris::Cell<moris::Cell<moris_index>>       mSubphaseToSubPhase;
            moris::Cell<moris::Cell<moris_index>>       mSubphaseToSubPhaseMySideOrds;
            moris::Cell<moris::Cell<moris_index>>       mSubphaseToSubPhaseNeighborSideOrds;

            // in the case of a hierarchically refined mesh, there are transitions with hanging nodes
            // this data flags the transition from a large facet to a smaller facet. (this is trivial
            // for non hmr type meshes)
            moris::Cell<moris::Cell<moris_index>> mTransitionNeighborCellLocation;

            // local to global subphase map
            std::unordered_map<moris::moris_id,moris::moris_index> mGlobalToLocalSubphaseMap;

            std::shared_ptr< mtk::Mesh_Manager > mMTKInputPerformer  = nullptr;
            std::shared_ptr< mtk::Mesh_Manager > mMTKOutputPerformer = nullptr;

            bool mInitializeCalled = false;

            // Intersection algorithm data members
            // They contain added clusters,cells, vertices
            moris::mtk::Intersection_Detect* mIntersectionDetect=nullptr;
            moris::mtk::Intersection_Detect_2D* mIntersectionDetect2D=nullptr;


            // Private Functions
        private:
            //------------------------------------------------------------------------------

            // Functions that take parameter inputs and setup XTK inputs
            bool
            has_parameter_list();


            //------------------------------------------------------------------------------

            /*!
             * Verifys provided parameter list
             */
            bool
            valid_parameters();

            //------------------------------------------------------------------------------

            /*
             * Using the parameter list, figure out the cell of subdivision methods
             */
            Cell<enum Subdivision_Method>
            get_subdivision_methods();


            //------------------------------------------------------------------------------

            void
            create_new_node_association_with_geometry(Decomposition_Data & tDecompData);

             //------------------------------------------------------------------------------

             /*!
              * Verifies that the nodes in a decomposition have all been assigned node ids
              */
             bool
             verify_successful_node_assignment(Decomposition_Data & aDecompData);

             //------------------------------------------------------------------------------

        protected:
             void
             send_outward_requests(
                     moris_index            const & aMPITag,
                     Cell<uint>             const & aProcRanks,
                     Cell<Matrix<IndexMat>> & aOutwardRequests);

             //------------------------------------------------------------------------------

             void
             inward_receive_requests(
                     moris_index            const & aMPITag,
                     moris::uint                    aNumRows,
                     Cell<Matrix<IndexMat>> &       aReceivedData,
                     Cell<uint>             &       aProcRanksReceivedFrom);

             //------------------------------------------------------------------------------

             void
             prepare_request_answers(
                     Decomposition_Data           & aDecompData,
                     Cell<Matrix<IndexMat>> const & aReceiveData,
                     Cell<Matrix<IndexMat>>       & aReceivedRequestAnswers);

             //------------------------------------------------------------------------------

             void
             return_request_answers(
                     moris_index            const & aMPITag,
                     Cell<Matrix<IndexMat>> const & aRequestAnswers,
                     Cell<uint>             const & aProcRanks);

             //------------------------------------------------------------------------------

             void
             inward_receive_request_answers(
                     moris_index            const & aMPITag,
                     moris::uint            const & aNumRows,
                     Cell<uint>             const & aProcRanks,
                     Cell<Matrix<IndexMat>> &       aReceivedData);

             //------------------------------------------------------------------------------

             void
             handle_received_request_answers(
                     Decomposition_Data           & aDecompData,
                     Cell<Matrix<IndexMat>> const & aRequests,
                     Cell<Matrix<IndexMat>> const & aRequestAnswers,
                     moris::moris_id              & aNodeId);

             //------------------------------------------------------------------------------
             // moris real versions of above
             void
             send_outward_requests_reals(
                     moris_index const    & aMPITag,
                     Cell<uint>  const    & aProcRanks,
                     Cell<Matrix<DDRMat>> & aOutwardRequests);

             //------------------------------------------------------------------------------

             void
             inward_receive_requests_reals(
                     moris_index const &    aMPITag,
                     moris::uint            aNumRows,
                     Cell<Matrix<DDRMat>> & aReceivedData,
                     Cell<uint>           & aProcRanksReceivedFrom);

             //------------------------------------------------------------------------------

             void
             return_request_answers_reals(
                     moris_index const & aMPITag,
                     Cell<Matrix<DDRMat>> const & aRequestAnswers,
                     Cell<uint>              const & aProcRanks);

             //------------------------------------------------------------------------------

             void
             inward_receive_request_answers_reals(
                     moris_index          const & aMPITag,
                     moris::uint          const & aNumRows,
                     Cell<uint>           const & aProcRanks,
                     Cell<Matrix<DDRMat>>       & aReceivedData);

            //------------------------------------------------------------------------------

            moris::Cell<std::string>
            check_for_and_remove_internal_seacas_side_sets(moris::Cell<std::string> & aSideSetNames);


        private:

              //------------------------------------------------------------------------------
              // Enrichment computation functions
              //------------------------------------------------------------------------------

              void
              perform_basis_enrichment_internal(
                      enum EntityRank   const & aBasisRank,
                      Matrix<IndexMat> const & aMeshIndex);



    };
}

#endif /* SRC_XTK_CL_XTK_MODEL_HPP_ */
