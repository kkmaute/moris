
/*
 * cl_XTK_Enrichment.hpp
 *
 *  Created on: Feb 23, 2018
 *      Author: ktdoble
 */

#ifndef XTK_SRC_XTK_CL_XTK_ENRICHMENT_HPP_
#define XTK_SRC_XTK_CL_XTK_ENRICHMENT_HPP_

// XTKL: Linalg Includes
#include "cl_Matrix.hpp"
#include "cl_XTK_Matrix_Base_Utilities.hpp"

// Std includes
#include <limits>

// XTKL: XTK Includes
#include "cl_Cell.hpp"
#include "cl_XTK_Model.hpp"
#include "cl_XTK_Child_Mesh.hpp"
#include "cl_XTK_Cut_Mesh.hpp"
#include "cl_XTK_Subphase_Group.hpp"
#include "fn_mesh_flood_fill.hpp"
#include "fn_prune_element_to_element.hpp"
#include "fn_generate_element_to_element.hpp"
#include "fn_local_child_mesh_flood_fill.hpp"
#include "fn_mesh_flood_fill.hpp"
#include "fn_Pairing.hpp"
#include "fn_equal_to.hpp"

// Mesh includes
#include "cl_MTK_Mesh.hpp"
#include "cl_MTK_Cell.hpp"
#include "cl_Mesh_Enums.hpp"
#include "cl_XTK_Background_Mesh.hpp"
#include "cl_Mesh_Enums.hpp"

#include "fn_unique.hpp"

#include "cl_XTK_Vertex_Enrichment.hpp"
#include "cl_MTK_Vertex_Interpolation.hpp"

#include "cl_TOL_Memory_Map.hpp"
/*
 * This class provides all the functions to perform the enrichment strategy on a child mesh
 */
namespace xtk
{
    // ----------------------------------------------------------------------------------

    enum class Enrichment_Method
    {
            USE_INTEGRATION_CELL_BASIS,// this one computes directly the vertex interpolation (uses the basis of tetrahedral cell)
            USE_INTERPOLATION_CELL_BASIS,// This one constructs an interpolation cell which interpolates into each subphase. (Uses basis of interpolation cell)
            INVALID
    };

    // ----------------------------------------------------------------------------------

    class Enrichment_Parameters
    {
        public:
            Enrichment_Parameters(){};

            enum moris::EntityRank mBasisToEnrich = EntityRank::NODE; /*For lagrange mesh this is node, for HMR this may be bsplines*/
    };

    // ----------------------------------------------------------------------------------

    class Enrichment_Data
    {
        public:
            Enrichment_Data( moris_index aNumSubphases ) :
                mSubphaseBGBasisIndices( aNumSubphases ),
                mSubphaseBGBasisEnrLev( aNumSubphases ){};

            // Enrichment Data ordered by basis function indices
            // For each basis function, the element indices and elemental enrichment levels //?
            Cell< moris::Matrix< moris::IndexMat > > mElementEnrichmentLevel;
            Cell< moris::Matrix< moris::IndexMat > > mElementIndsInBasis;

            // For each enriched basis function, the subphase indices in support
            Cell< moris::Matrix< moris::IndexMat > > mSubphaseIndsInEnrichedBasis;

            // Basis enrichment level indices
            moris::Cell< moris::Matrix< moris::IndexMat > > mBasisEnrichmentIndices; //???
            moris::Matrix< moris::IndexMat >                mEnrichedBasisIndexToId;

            // Unintersected Parent Cell, BackBasis interpolating in them and corresponding enrichment level
            // outer cell corresponds to interp cell index
            // inner cell corresponds to basis/enr-lev in interp cell
            moris::Cell< moris::Cell< moris_index > > mSubphaseBGBasisIndices;
            moris::Cell< moris::Cell< moris_index > > mSubphaseBGBasisEnrLev;

            // total number of basis enrichment levels (all basis functions)
            moris::uint mNumEnrichmentLevels;

            // Vertex interpolations for this enrichment ordered by background vertex index
            Cell< mtk::Vertex_Interpolation* > mBGVertexInterpolations;

            // Basis bulk measures
    };

    // ----------------------------------------------------------------------------------

    class Model;
    class Cut_Integration_Mesh;
    class Integration_Mesh_Generator;

    // ----------------------------------------------------------------------------------

    class Enrichment
    {
        public:
            Enrichment(){};

            Enrichment(
                    enum Enrichment_Method const& aMethod,
                    enum EntityRank const&        aBasisRank,
                    Matrix< IndexMat > const&     aInterpIndex,
                    moris::moris_index const&     aNumBulkPhases,
                    xtk::Model*                   aXTKModelPtr,
                    moris::mtk::Mesh*             aBackgroundMeshPtr,
                    bool                          aSortBasisEnrichmentLevels);

            // ----------------------------------------------------------------------------------

            ~Enrichment();
            // ----------------------------------------------------------------------------------

            bool mVerbose = false; /* Verbose output of enrichment information flag */

            // ----------------------------------------------------------------------------------

            /**
             * Perform basis function enrichment
             */
            void
            perform_enrichment();

            // ----------------------------------------------------------------------------------
            // Accessing enrichment data
            // ----------------------------------------------------------------------------------

            /**
             * Returns the subphases indices in enriched basis functions support
             * Note: Used in XTK Multigrid
             */
            Cell< moris::Matrix< moris::IndexMat > > const&
            get_subphases_loc_inds_in_enriched_basis( moris_index const& aEnrichmentDataIndex = 0 ) const;

            // ----------------------------------------------------------------------------------
            // Field writing for debugging
            // ----------------------------------------------------------------------------------
            /*
             * Provided an MTK mesh, writes the cell enrichment data onto the mesh. (fields declared from get_cell_enrichment_field_names )
             */
            void
            write_cell_enrichment_to_fields(
                    Cell< std::string >& aEnrichmentFieldStrs,
                    mtk::Mesh*           aMeshWithEnrFields ) const;

            // ----------------------------------------------------------------------------------

            /**
             * Returns a vector of cell fields names to declare in STK mesh if you want to visualize the cell level
             * enrichment fields. cells within each subphase of a given basis function. One field per basis function, one field per child mesh
             */
            Cell< std::string >
            get_cell_enrichment_field_names() const;

            // ----------------------------------------------------------------------------------

            /**
             * Returns the element inds in a basis support constructed in call to perform_enrichment. These are indexed by basis function index.
             * Note: Only used for debug outputting in function write_cell_enrichment_to_fields.
             */
            Cell< moris::Matrix< moris::IdMat > > const&
            get_element_inds_in_basis_support( moris_index const& aEnrichmentDataIndex = 0 ) const;

            // ----------------------------------------------------------------------------------

            /**
             * Returns the element enrichment levels in a basis support constructed in call to perform_enrichment. These are indexed by basis function index.
             * Correspond to the element inds found at the same index in mElementIndsInBasis.
             * Note: Only used for debug outputting in function write_cell_enrichment_to_fields.
             */
            Cell< moris::Matrix< moris::IndexMat > > const&
            get_element_enrichment_levels_in_basis_support( moris_index const& aEnrichmentDataIndex = 0 ) const;

            // ----------------------------------------------------------------------------------

            /**
             * @brief get the memory usage of enrichment
             */
            moris::Memory_Map
            get_memory_usage();

            // ----------------------------------------------------------------------------------

            void
            write_diagnostics();

            void
            print_enriched_basis_to_subphase_id(
                    const moris_index & aMeshIndex,
                    std::string aFileName);

            typedef std::unordered_map< moris::moris_index, moris::moris_index > IndexMap;
            friend class Multigrid;

        private:

            // enrichment method
            enum Enrichment_Method mEnrichmentMethod;

            // basis rank
            enum EntityRank mBasisRank; /*Entity rank of the basis functions*/

            // index of interpolation
            Matrix< IndexMat > mMeshIndices; /* Mesh indices to perform enrichment on*/

            // number of bulk-phases possible in model
            moris::size_t mNumBulkPhases;

            // Pointers to necessary classes
            Model*                      mXTKModelPtr       = nullptr;
            Cut_Integration_Mesh*       mCutIgMesh         = nullptr;
            Integration_Mesh_Generator* mIgMeshTools       = nullptr;
            moris::mtk::Mesh*           mBackgroundMeshPtr = nullptr;

            // enrichment strategy data (outer cell - mesh index, inner cell - necessary data for enrichment of mesh index)
            Cell< Enrichment_Data > mEnrichmentData;

            // flag whether to sort basis enrichment levels
            bool mSortBasisEnrichmentLevels;

            // B-spline - Lagrange Mesh information
            moris::Cell< Bspline_Mesh_Info * > mBsplineMeshInfos;

            // ----------------------------------------------------------------------------------

            /*
             * Performs enrichment on elements in support of full basis cluster. This enrichment includes all children elements of parents in
             * the basis cluster and parent elements with no children
             */
            void
            perform_basis_cluster_enrichment();

            // ----------------------------------------------------------------------------------

            /*
             * Constructs the subphase neighborhood in the XTK model
             */
            void
            construct_neighborhoods();

            // ----------------------------------------------------------------------------------

            /*
             * Verifies vertex interpolation is present and communicates non-existing vertex interpolation.
             * The vertex interpolation on the aura needs to be communicated if we are using HMR
             */
            void
            setup_background_vertex_interpolations();

            // ----------------------------------------------------------------------------------

            Matrix< IndexMat >
            get_subphase_clusters_in_support( moris::Matrix< moris::IndexMat > const& aElementsInSupport );

            // ----------------------------------------------------------------------------------

            Matrix< IndexMat >
            get_subphase_clusters_in_bspline_cell( moris::Cell< mtk::Cell* > const& aCellsInGroup );

            // ----------------------------------------------------------------------------------

            void
            construct_subphase_in_support_map(
                    moris::Matrix< moris::IndexMat > const& aSubphaseClusterIndicesInSupport,
                    IndexMap&                               aSubPhaseIndexToSupportIndex );

            // ----------------------------------------------------------------------------------
            void
            construct_subphase_in_bspline_cell_map(
                    moris::Matrix< moris::IndexMat > const& aSubphaseIndicesInBsplineCell,
                    IndexMap&                               aSubphaseIndexToBsplineCellIndex );

            // ----------------------------------------------------------------------------------

            /**
             * @brief Construct full subphase neighbor graph in basis support and the corresponding shared faces
             * 
             * @param aSubphasesInSupport 
             * @param aSubPhaseIndexToSupportIndex 
             * @param aPrunedSubPhaseToSubphase 
             */
            void
            generate_pruned_subphase_graph_in_basis_support(
                    moris::Matrix< moris::IndexMat > const& aSubphasesInSupport,
                    IndexMap&                               aSubPhaseIndexToSupportIndex,
                    moris::Matrix< moris::IndexMat >&       aPrunedSubPhaseToSubphase );

            // ----------------------------------------------------------------------------------

            void
            generate_pruned_subphase_graph_in_bspline_cell(
                    moris::Matrix< moris::IndexMat > const& aSubphasesInBsplineCell,
                    IndexMap&                               aSubphaseIndicesToBspline,
                    moris::Matrix< moris::IndexMat >&       aPrunedBsplineSubphaseToSubphase );

            // ----------------------------------------------------------------------------------

            void
            assign_subphase_bin_enrichment_levels_in_basis_support(
                    moris::Matrix< moris::IndexMat > const& aSubphasesInSupport,
                    IndexMap&                               aSubPhaseIndexToSupportIndex,
                    moris::Matrix< moris::IndexMat > const& aPrunedSubPhaseToSubphase,
                    moris::Matrix< moris::IndexMat >&       aSubPhaseBinEnrichmentVals,
                    moris_index&                            aMaxEnrichmentLevel );

            // ----------------------------------------------------------------------------------

            void
            create_subphase_groups(
                    moris::Cell< mtk::Cell* >& aIpCellsInGroup,
                    moris::size_t              aBsplineMeshIndex );

            // ----------------------------------------------------------------------------------

            /**
             * @brief splits up the bin of subphases into lists of subphase indices with the same color
             * 
             * @param aSubphaseBin 
             * @param aSubphaseIndicesInBsplineCell 
             * @return moris::Cell< moris::Cell< moris_index > > 
             */
            moris::Cell< moris::Cell< moris_index > > 
            split_flood_fill_bin( 
                    moris::Matrix< moris::IndexMat > aSubphaseBin, 
                    moris::Matrix< moris::IndexMat > aSubphaseIndicesInBsplineCell,
                    uint                             aNumSPGs );

            // ----------------------------------------------------------------------------------

            moris::Cell< bool >
            collect_subphase_group_ligament_side_ordinals( 
                    moris::Cell< moris_index >& aSPsInGroup,
                    IndexMap&                   aSubphaseIndicesToBspline );

            // ----------------------------------------------------------------------------------

            /*
             * @brief Sort enrichment levels by distance of centroid  of subphases with same
             *        enrichment level from origin
             */

            void
            sort_enrichment_levels_in_basis_support(
                    moris::Matrix<moris::IndexMat> const& aSubphasesInSupport,
                    moris::Matrix<moris::IndexMat>&       aSubPhaseBinEnrichmentVals,
                    moris_index const                     aMaxEnrichmentLevel);

            // ----------------------------------------------------------------------------------

            void
            unzip_subphase_bin_enrichment_into_element_enrichment(
                    moris_index const&                      aEnrichmentDataIndex,
                    moris_index const&                      aBasisIndex,
                    moris::Matrix< moris::IndexMat > const& aParentElementsInSupport,
                    moris::Matrix< moris::IndexMat > const& aSubphasesInSupport,
                    IndexMap&                               aSubPhaseIndexToSupportIndex,
                    moris::Matrix< moris::IndexMat > const& aPrunedSubPhaseToSubphase,
                    moris::Matrix< moris::IndexMat >&       aSubPhaseBinEnrichmentVals );

            // ----------------------------------------------------------------------------------

            void
            construct_enriched_basis_to_subphase_connectivity(
                    moris_index const&                                     aEnrichmentDataIndex,
                    moris::Cell< moris::Matrix< moris::IndexMat > > const& aSubPhaseBinEnrichment,
                    moris::Cell< moris::Matrix< moris::IndexMat > > const& aSubphaseClusterIndicesInSupport,
                    moris::Cell< moris_index > const&                      aMaxEnrichmentLevel );

            // ----------------------------------------------------------------------------------

            /*!
             * Assign the enrichment level local identifiers
             */
            void
            assign_enriched_coefficients_identifiers(
                    moris_index const&                aEnrichmentDataIndex,
                    moris::Cell< moris_index > const& aMaxEnrichmentLevel );

            // ----------------------------------------------------------------------------------

            void
            communicate_basis_information_with_owner(
                    moris_index const&                        aEnrichmentDataIndex,
                    Cell< Cell< moris_index > > const&        aBasisIdToBasisOwner,
                    Cell< Cell< moris_index > > const&        aMaxSubphaseIdToBasisOwner,
                    Cell< moris_index > const&                aProcRanks,
                    std::unordered_map< moris_id, moris_id >& aProcRankToIndexInData,
                    Cell< moris::Matrix< moris::IndexMat > >& aEnrichedBasisId );

            // ----------------------------------------------------------------------------------

            void
            set_received_enriched_basis_ids(
                    moris_index const&                              aEnrichmentDataIndex,
                    Cell< moris::Matrix< moris::IndexMat > > const& aReceivedEnrichedIds,
                    Cell< Cell< moris_index > > const&              aBasisIndexToBasisOwner,
                    Cell< Cell< moris_index > > const&              aSubphaseIdInSupport );

            // ----------------------------------------------------------------------------------

            moris::size_t
            count_elements_in_support( moris::Matrix< moris::IndexMat > const& aParentElementsInSupport );

            // ----------------------------------------------------------------------------------
            bool
            subphase_is_in_support(
                    moris_index const& aEnrichmentDataIndex,
                    moris_index        aSubphaseIndex,
                    moris_index        aEnrichedBasisIndex );
            // ----------------------------------------------------------------------------------
            void
            print_basis_support_debug(
                    moris_index                             aBasisIndex,
                    moris::Matrix< moris::IndexMat > const& aParentElementsInSupport,
                    moris::Matrix< moris::IndexMat > const& aSubphasesInSupport,
                    IndexMap&                               aSubPhaseIndexToSupportIndex,
                    moris::Matrix< moris::IndexMat > const& aPrunedSubPhaseToSubphase,
                    moris::Matrix< moris::IndexMat >&       aSubPhaseBinEnrichmentVals );

            // ----------------------------------------------------------------------------------
            // Setup enriched interpolation mesh
            // ----------------------------------------------------------------------------------
            void
            construct_enriched_interpolation_mesh();

            // ----------------------------------------------------------------------------------

            void
            construct_enriched_integration_mesh();

            // ----------------------------------------------------------------------------------

            void
            allocate_interpolation_cells();

            // ----------------------------------------------------------------------------------

            void
            construct_enriched_interpolation_vertices_and_cells();

            // ----------------------------------------------------------------------------------

            void
            construct_enriched_vertex_interpolation(
                    moris_index const&                        aEnrichmentDataIndex,
                    mtk::Vertex_Interpolation*                aBaseVertexInterp,
                    Cell< moris_index > const&                aSubPhaseBasisEnrLev,
                    std::unordered_map< moris_id, moris_id >& aMapBasisIndexToLocInSubPhase,
                    Vertex_Enrichment&                        aVertexEnrichment );
            // ----------------------------------------------------------------------------------

            std::unordered_map< moris_id, moris_id >
            construct_subphase_basis_to_basis_map( Cell< moris_id > const& aSubPhaseBasisIndex );

            // ----------------------------------------------------------------------------------
            // functions that collects the appropriate vertex interpolation for the enrichment strategy only
            moris::Cell< mtk::Vertex_Interpolation* >
            get_vertex_interpolations(
                    moris::mtk::Cell& aParentCell,
                    const uint        aMeshIndex ) const;

            // ----------------------------------------------------------------------------------
            /*!
             *   @param[in] aMeshIndex Mesh index
             *   @param[in] aNumIdsToAllocate Number of basis ids needed
             *   @return Basis index offset
             */
            moris_index
            allocate_basis_ids( moris_index const& aMeshIndex,
                    moris_index const&                 aNumIdsToAllocate );
            // ----------------------------------------------------------------------------------

            /*!
             *   @param[in] aMeshIndex Mesh index
             *   @return Maxmimum basis id on this proc
             */
            moris_index
            get_max_basis_id( moris_index const& aMeshIndex );
    };
}// namespace xtk
#endif /* XTK_SRC_XTK_CL_XTK_ENRICHMENT_HPP_ */
