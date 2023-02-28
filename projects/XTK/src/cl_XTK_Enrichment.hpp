/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_XTK_Enrichment.hpp
 *
 */

#ifndef XTK_SRC_XTK_CL_XTK_ENRICHMENT_HPP_
#define XTK_SRC_XTK_CL_XTK_ENRICHMENT_HPP_

// XTKL: Linalg Includes
#include "cl_Matrix.hpp"
#include "cl_XTK_Matrix_Base_Utilities.hpp"

// Std includes
#include <limits>
#include <unordered_set>

// XTKL: XTK Includes
#include "cl_Cell.hpp"
#include "cl_XTK_Model.hpp"
#include "cl_XTK_Child_Mesh.hpp"
#include "cl_XTK_Cut_Mesh.hpp"
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
#include "cl_XTK_Subphase_Group.hpp"

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
            Enrichment_Data( 
                    moris_index aNumSubphases ) :
                mSubphaseBGBasisIndices( aNumSubphases ),
                mSubphaseBGBasisEnrLev( aNumSubphases ){};

            // FIXME: this needs to go once everything is moved over to SPG bases Enrichment
            void
            reinitialize_for_SPGs( moris_index aNumSPGs )
            {
                // resize maps to number of SPGs
                mSubphaseGroupBGBasisIndices.resize( aNumSPGs );
                mSubphaseGroupBGBasisEnrLev.resize ( aNumSPGs );

                // set flag whether SP or SPG based enrichment is used
                mUseSPGs = true;
            } 

            // Enrichment Data ordered by basis function indices
            // For each basis function, the element indices and elemental enrichment levels //?
            Cell< moris::Matrix< moris::IndexMat > > mElementIndsInBasis; // input: non-enriched BF index || output: list of IG cells within the basis support
            Cell< moris::Matrix< moris::IndexMat > > mElementEnrichmentLevel; // input: non-enriched BF index || output: which enrichment level of current basis is active on the IG cell

            // For each enriched basis function, the subphase indices in support
            // input: enriched BF index || output: list of subphase indices in which enriched BF is active
            Cell< moris::Matrix< moris::IndexMat > > mSubphaseIndsInEnrichedBasis;
            
            // FIXME: for SPG based enrichment, the above needs to go eventually
            Cell< moris::Matrix< moris::IndexMat > > mSubphaseGroupIndsInEnrichedBasis; // input: enriched BF index || output: list of SPG indices in which enriched BF is active
            Matrix< IdMat > mBulkPhaseInEnrichedBasis; // input: enriched BF index || output: list of SPG indices in which enriched BF is active

            // Basis enrichment level indices
            moris::Cell< moris::Matrix< moris::IndexMat > > mBasisEnrichmentIndices; // input1: non-enriched BF index, input2: enrichment level || output: enriched BF index
            moris::Matrix< moris::IndexMat >                mEnrichedBasisIndexToId; // input: enriched BF index || output: global ID of that enriched BF 

            // Unintersected Parent Cell, BackBasis interpolating in them and corresponding enrichment level
            // outer cell corresponds to interp cell index
            // inner cell corresponds to basis/enr-lev in interp cell
            moris::Cell< moris::Cell< moris_index > > mSubphaseBGBasisIndices; // input: subphase index || output: list of non-enriched BF indices active on given subphase
            moris::Cell< moris::Cell< moris_index > > mSubphaseBGBasisEnrLev; // input: subphase index || output: list of enrichment levels for the respective BFs active on the given subphase
            
            // FIXME: for SPG based enrichment, the above needs to go eventually
            moris::Cell< moris::Cell< moris_index > > mSubphaseGroupBGBasisIndices; // input: SPG index || output: list of non-enriched BF indices active on given subphase group
            moris::Cell< moris::Cell< moris_index > > mSubphaseGroupBGBasisEnrLev; // input: SPG index || output: list of enrichment levels for the respective BFs active on the given subphase group

            // total number of basis enrichment levels (i.e. number of enriched Basis functions)
            moris::uint mNumEnrichmentLevels;

            // Vertex interpolations for this enrichment ordered by background vertex index
            Cell< mtk::Vertex_Interpolation* > mBGVertexInterpolations;

            // flag whether SPG or SP based Enrichment is used
            bool mUseSPGs = false;
    };

    // ----------------------------------------------------------------------------------

    class Model;
    class Cut_Integration_Mesh;
    class Integration_Mesh_Generator;

    // ----------------------------------------------------------------------------------

    class Enrichment
    {
        public: 
            typedef std::unordered_map< moris::moris_index, moris::moris_index > IndexMap;
            friend class Multigrid;

        private:

            // enrichment method
            enum Enrichment_Method mEnrichmentMethod;

            // basis rank
            enum EntityRank mBasisRank; /*Entity rank of the basis functions*/

            // index of interpolation
            Matrix< IndexMat > mMeshIndices; /* Mesh indices to perform enrichment on*/
            IndexMap mMeshIndexToListIndexMap; /* map that contains the position of a given mesh index in the above list */
            bool mMeshIndexMapIsSet = false; /* flag that marks whether the above index map has been computed or not */

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

            // quick access to the Cut integration mesh's Bspline Mesh Infos (the Bspline to Lagrange mesh relationships)
            moris::Cell< Bspline_Mesh_Info* > mBsplineMeshInfos;

            // Maps tracking how an IP cell gets unzipped
            moris_index mNumEnrIpCells; 
            moris::Cell< uint > mNumUnzippingsOnIpCell; // input: IP-cell index || output: number of enr. IP-cells and clusters to be created
            moris::Cell< moris::Cell<  moris::Cell< moris_index > > > mMaterialSpgsUnzippedOnIpCell; // input: B-spline mesh index, IP-cell index || output: list of SPG indices wrt. which clusters containing material need to be created
            moris::Cell< moris::Cell<  moris::Cell< moris_index > > > mVoidSpgsUnzippedOnIpCell; // input: B-spline mesh index, IP-cell index || output: list of SPG indices wrt. which void clusters need to be created

            // indices of enriched IP cells as a function of the base IP cell and the local SPG index
            moris::Cell< moris::Cell< moris_index > > mEnrIpCellIndices; // input: base IP cell index, index of unzipping || output: index of enr. IP cell
            moris::Cell< moris_index > mUipcUnzippingIndices; // input: UIPC index || output: unzipping index

            // map allowing correct UIPC to be grabbed based on base IP cell and SPG index (for Ghost)
            // input: B-spline mesh index, base IP cell index, SPG index local to corresponding B-spline element || output: Enr. IP cell index
            Cell< Cell< Cell< moris_index > > > mBaseIpCellAndSpgToUnzipping;

            // map allowing to enr. IP cell / non-void cluster associated with a given subphase
            // input: subphase index || output: enriched interpolation cell index
            Cell< moris_index > mSubphaseIndexToEnrIpCellIndex; 

            // maps associating UIPCs / Cell clusters with SPGs
            // input: B-spline mesh list index, SPG index || output: List of enr. IP cells / cell clusters on SPG
            Cell< Cell< Cell< moris_index > > > mSpgToUipcIndex;

            // input: B-spline mesh list index, Cluster/UIPC index || output: SPG index this Cluster/UIPC is in
            Cell< Cell< moris_index > > mUipcToSpgIndex;

            // input: B-spline mesh index, base IP cell index, index of unzipping || output: SPG index the UIPC belongs to
            Cell< Cell< Cell< moris_index > > > mUnzippingToSpgIndex;

            // ----------------------------------------------------------------------------------

        public:

            // ----------------------------------------------------------------------------------
            
            Enrichment(){};

            Enrichment(
                    enum Enrichment_Method const& aMethod,
                    enum EntityRank const&        aBasisRank,
                    Matrix< IndexMat > const&     aInterpIndex,
                    moris::moris_index const&     aNumBulkPhases,
                    xtk::Model*                   aXTKModelPtr,
                    moris::mtk::Mesh*             aBackgroundMeshPtr,
                    bool                          aSortBasisEnrichmentLevels,
                    bool                          aUseSpgBasedEnrichment = false );

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

            void
            perform_enrichment_new();

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

            // ----------------------------------------------------------------------------------

            void
            print_enriched_basis_to_subphase_id(
                    const moris_index & aMeshIndex,
                    std::string aFileName);

            // ----------------------------------------------------------------------------------
            
            Matrix< IndexMat >
            get_mesh_indices() const 
            {
                return mMeshIndices;
            }

            // ----------------------------------------------------------------------------------
            
            moris::Cell< moris::Cell< moris_index > > const&
            get_unzipped_IP_cell_to_base_IP_cell_and_SPG( const moris_index aBsplineMeshListIndex ) const 
            {
                return mBaseIpCellAndSpgToUnzipping( aBsplineMeshListIndex );
            }

            // ----------------------------------------------------------------------------------

            /**
             * @brief Get the bspline mesh info for a given B-spline mesh
             * 
             * @param aMeshIndexInList position of the B-spline mesh in the list of discretization meshes
             * @return Bspline_Mesh_Info const* 
             */
            Bspline_Mesh_Info const* 
            get_bspline_mesh_info_for_list_index( moris_index aMeshIndexInList ) const
            {
                return mBsplineMeshInfos( aMeshIndexInList );
            }

            // ----------------------------------------------------------------------------------

            uint
            get_num_unzippings_of_base_ip_cell( moris_index aBaseIpCellIndex ) const
            {
                return mEnrIpCellIndices( aBaseIpCellIndex ).size();
            }

            // ----------------------------------------------------------------------------------
            
            moris::Cell< moris_index > const&
            get_subphase_to_UIPC_map() const 
            {
                return mSubphaseIndexToEnrIpCellIndex;
            }

            // ----------------------------------------------------------------------------------

            uint
            get_num_enr_ip_cells() const
            {
                return (uint) mNumEnrIpCells;
            }

            // ----------------------------------------------------------------------------------

            moris::Cell< moris_index > const&
            get_enr_ip_cell_indices_on_base_ip_cell( moris_index aBaseIpCellIndex ) const
            {
                return mEnrIpCellIndices( aBaseIpCellIndex );
            }

            // ----------------------------------------------------------------------------------

            moris::Cell< moris::Cell< moris::Cell< moris_index > > > const&
            get_SPG_to_UIPC_map() const
            {
                return mSpgToUipcIndex;
            }

            // ----------------------------------------------------------------------------------

            moris::Cell< moris_index > const&
            get_UIPC_indices_on_SPG( 
                    const moris_index aBspMeshListIndex, 
                    const moris_index aSpgIndex ) const
            {
                return mSpgToUipcIndex( aBspMeshListIndex )( aSpgIndex );
            }

            // ----------------------------------------------------------------------------------
            
            moris::Cell< moris::Cell< moris_index > > const&
            get_UIPC_to_SPG_map() const
            {
                return mUipcToSpgIndex;
            }

            // ----------------------------------------------------------------------------------

            moris_index
            get_SPG_on_UIPC( 
                    const moris_index aBspMeshListIndex, 
                    const moris_index aEnrIpCellIndex ) const
            {
                return mUipcToSpgIndex( aBspMeshListIndex )( aEnrIpCellIndex );
            }

            // ----------------------------------------------------------------------------------

        private:

            // ----------------------------------------------------------------------------------

            /**
             * @brief Get the position of a given mesh index in the list of mesh indices mMeshIndices
             * 
             * @param aMeshIndex mesh index
             * @return moris_index position of the mesh index in the list of MeshIndices
             */
            moris_index
            get_list_index_for_mesh_index( moris_index aMeshIndex );

            // ----------------------------------------------------------------------------------

            /**
             * @brief Performs enrichment on elements in support of full basis cluster using the 
             * previously construced Subphase Groups. This enrichment includes all children 
             * elements of parents in the basis cluster, parent elements with no children  
             */
            void
            perform_basis_cluster_enrichment_new();

            // ----------------------------------------------------------------------------------

            /**
             * @brief Performs enrichment on elements in support of full basis cluster. This 
             * enrichment includes all children elements of parents in the basis cluster and 
             * parent elements with no children
             */
            void
            perform_basis_cluster_enrichment();

            // ----------------------------------------------------------------------------------

            /**
             * @brief Constructs the subphase neighborhood in the XTK model
             */
            void
            construct_neighborhoods();

            // ----------------------------------------------------------------------------------

            /**
             * @brief Verifies vertex interpolation is present and communicates non-existing 
             * vertex interpolation. The vertex interpolation on the aura needs to be 
             * communicated if we are using HMR
             */
            void
            setup_background_vertex_interpolations();

            // ----------------------------------------------------------------------------------

            Matrix< IndexMat >
            get_subphase_clusters_in_support( moris::Matrix< moris::IndexMat > const& aElementsInSupport );

            // ----------------------------------------------------------------------------------

            void
            get_subphase_groups_in_support( 
                moris_index                             aMeshIndexPosition, 
                moris::Matrix< moris::IndexMat > const& aLagElementsInSupport,
                moris::Matrix< moris::IndexMat >&       aSubphaseGroupIndicesInSupport,
                IndexMap&                               aSubphaseGroupIndexToSupportIndex );

            // ----------------------------------------------------------------------------------

            void
            construct_subphase_in_support_map(
                    moris::Matrix< moris::IndexMat > const& aSubphaseClusterIndicesInSupport,
                    IndexMap&                               aSubPhaseIndexToSupportIndex );

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

            /**
             * @brief Construct subphase group neighbor graph in basis support and the corresponding shared faces
             * 
             * @param aSubphaseGroupIndicesInSupport 
             * @param aSubphaseGroupIndexToSupportIndex 
             * @param aPrunedSubphaseGroupToSubphase 
             */
            void
            generate_pruned_subphase_group_graph_in_basis_support(
                    const moris_index                       aMeshIndex,
                    moris::Matrix< moris::IndexMat > const& aSubphaseGroupIndicesInSupport,
                    IndexMap&                               aSubphaseGroupIndexToSupportIndex,
                    moris::Matrix< moris::IndexMat >&       aPrunedSubphaseGroupToSubphase );

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
            assign_subphase_group_bin_enrichment_levels_in_basis_support(
                    moris::Matrix< moris::IndexMat > const & aSpgsInSupport,
                    moris::Matrix< moris::IndexMat > const & aPrunedSpgToSpg,
                    moris::Matrix< moris::IndexMat >&        aSpgBinEnrichmentVals,
                    moris_index&                             aMaxEnrichmentLevel );

            // ----------------------------------------------------------------------------------

            /**
             * @brief Sort enrichment levels by distance of centroid  of subphases with same
             *        enrichment level from origin
             */

            void
            sort_enrichment_levels_in_basis_support(
                    moris::Matrix<moris::IndexMat> const& aSubphasesInSupport,
                    moris::Matrix<moris::IndexMat>&       aSubPhaseBinEnrichmentVals,
                    moris_index const                     aMaxEnrichmentLevel);

            // ----------------------------------------------------------------------------------

            /**
             * @brief go through all subphases in the support of a given basis funtion and save 
             * the current basis function's index to them, and which enrichment level of this 
             * given basis function is active on them. This data is stored in the Enrichment_Data
             * 
             * @param aEnrichmentDataIndex 
             * @param aBasisIndex 
             * @param aParentElementsInSupport 
             * @param aSubphasesInSupport 
             * @param aSubPhaseIndexToSupportIndex 
             * @param aPrunedSubPhaseToSubphase 
             * @param aSubPhaseBinEnrichmentVals 
             */
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

            /**
             * @brief go through all SPGs in the support of a given basis funtion and save the 
             * current BFs's index to them, and which enrichment level of this given BF is active 
             * on them. This data is stored in the Enrichment_Data
             * 
             * @param aEnrichmentDataIndex 
             * @param aBasisIndex 
             * @param aParentElementsInSupport 
             * @param aSpgsInSupport 
             * @param aSpgBinEnrichmentVals 
             */
            void
            unzip_subphase_group_bin_enrichment_into_element_enrichment(
                    moris_index const &                      aEnrichmentDataIndex,
                    moris_index const &                      aBasisIndex,
                    moris::Matrix< moris::IndexMat > const & aParentElementsInSupport,
                    moris::Matrix< moris::IndexMat > const & aSpgsInSupport,
                    moris::Matrix< moris::IndexMat >&        aSpgBinEnrichmentVals );

            // ----------------------------------------------------------------------------------

            /**
             * @brief counts number of enriched basis function indices and finds which subphases
             * a given enriched BF index is active on. The data is stored in the Enrichment_Data
             * 
             * @param aEnrichmentDataIndex 
             * @param aSubPhaseBinEnrichment 
             * @param aSubphaseClusterIndicesInSupport 
             * @param aMaxEnrichmentLevel 
             */
            void
            construct_enriched_basis_to_subphase_connectivity(
                    moris_index const&                                     aEnrichmentDataIndex,
                    moris::Cell< moris::Matrix< moris::IndexMat > > const& aSubPhaseBinEnrichment,
                    moris::Cell< moris::Matrix< moris::IndexMat > > const& aSubphaseClusterIndicesInSupport,
                    moris::Cell< moris_index > const&                      aMaxEnrichmentLevel );

            // ----------------------------------------------------------------------------------
            /**
             * @brief counts number of enriched basis function indices and finds which SPGs a
             * given enriched BF index is active on. The data is stored in the Enrichment_Data 
             * 
             * @param aEnrichmentDataIndex 
             * @param aSpgBinEnrichment 
             * @param aSpgIndicesInSupport 
             * @param aMaxEnrichmentLevel 
             */
            void
            construct_enriched_basis_to_subphase_group_connectivity(
                    moris_index const &                                     aEnrichmentDataIndex,
                    moris::Cell< moris::Matrix< moris::IndexMat > > const & aSpgBinEnrichment,
                    moris::Cell< moris::Matrix< moris::IndexMat > > const & aSpgIndicesInSupport,
                    moris::Cell< moris_index > const &                      aMaxEnrichmentLevel );

            // ----------------------------------------------------------------------------------

            /**
             * @brief Assign the enrichment level local identifiers (legacy enrichment procedure)
             */
            void
            assign_enriched_coefficients_identifiers(
                    moris_index const&                aEnrichmentDataIndex,
                    moris::Cell< moris_index > const& aMaxEnrichmentLevel );

            /**
             * @brief Assign the enrichment level local identifiers (in SPG based enrichment)
             */
            void
            assign_enriched_coefficients_identifiers_new(
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

            /**
             * @brief 
             * 
             * @param aEnrichmentDataIndex B-spline mesh index
             * @param aBasisIdToBasisOwner request lists: outer cell: proc index in comm table to request from, inner cell: list of non-enriched Basis IDs to request
             * @param aSubphaseGroupIdInSupport inner cell: list of SPG IDs for unique identification purposes
             * @param aProcRanks comm table: proc IDs associated with outer cell indices in above lists
             * @param aProcRankToIndexInData map corresponding to comm table
             * @param aEnrichedBasisId output: list of Basis IDs as requested above
             */
            void
            communicate_basis_information_with_owner_new(
                    moris_index const &                       aEnrichmentDataIndex,         
                    Cell< Cell< moris_index > > const &       aBasisIdToBasisOwner,         
                    Cell< Cell< moris_index > > const &       aSubphaseGroupIdInSupport,    
                    Cell< moris_index > const &               aProcRanks,                   
                    std::unordered_map< moris_id, moris_id >& aProcRankToIndexInData,       
                    Cell< moris::Matrix< moris::IndexMat > >& aEnrichedBasisId );

            // ----------------------------------------------------------------------------------

            void
            set_received_enriched_basis_ids(
                    moris_index const&                              aEnrichmentDataIndex,
                    Cell< moris::Matrix< moris::IndexMat > > const& aReceivedEnrichedIds,
                    Cell< Cell< moris_index > > const&              aBasisIndexToBasisOwner );

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

            bool
            subphase_group_is_in_support(
                    moris_index const& aBsplineMeshListIndex,
                    moris_index        aSubphaseGroupIndex,
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

            void
            construct_enriched_interpolation_mesh_new();

            // ----------------------------------------------------------------------------------

            void
            construct_enriched_integration_mesh();

            void
            construct_enriched_integration_mesh( const Matrix< IndexMat > aBsplineMeshIndices );
            
            // ----------------------------------------------------------------------------------

            void
            allocate_interpolation_cells();

            void
            allocate_interpolation_cells_based_on_SPGs();

            void
            allocate_interpolation_cells_based_on_SPGs_new();

            // ----------------------------------------------------------------------------------

            /**
             * @brief figure out if the basis associated with an existing SPG can be used for an enriched IP cell wrt. to a certain B-spline mesh
             * 
             * @param aIpCellIndex index of the base IP cell / Lagrange element
             * @param aUnzippingOnIpCell index of the unzipping (i.e. position of the enr. IP cell (UIPC) in the stack of UIPCs on this base IP cell)
             * @param aBsplineMeshIndex list index of the B-spline mesh whose basis is considered
             * @return moris_index index of the SPG whose basis will be used, return -1 if there is no SPG from which the basis can be used
             */
            
            moris_index
            get_SPG_for_basis_extension( 
                    const moris_index aIpCellIndex, 
                    const moris_index aUnzippingOnIpCell, 
                    const moris_index aBsplineMeshIndex ) const;

            // ----------------------------------------------------------------------------------

            /**
             * @brief Construct an averaged T-matrix from all material enr. IP cells for a given 
             * material sub-domain on the coarsest B-spline element containing a certain Lagrange 
             * element wrt. to a given B-spline mesh index.
             * 
             * @param aIpCellIndex index of the base IP cell / Lagrange element
             * @param aUnzippingOnIpCell index of the unzipping (i.e. position of the enr. IP cell (UIPC) in the stack of UIPCs on this base IP cell)
             * @param aBsplineMeshIndex list index of the B-spline mesh whose basis is considered
             * @param aAverageEnrichedTmatrix average T-matrix 
             */
            void
            construct_averaged_T_matrix_for_extension_on_enriched_IP_cell(
                    const moris_index  aIpCellIndex,
                    const moris_index  aUnzippingIndex,
                    const moris_index  aBsplineMeshIndex,
                    Vertex_Enrichment& aAverageEnrichedTmatrix ) const;

            // ----------------------------------------------------------------------------------
            
            /**
             * @brief computes a weighted average from multiple T-matrices
             * 
             * @param aAverageTmatrices list of T-matrices to average
             * @param aRelativeIpCellVolumes corresponding weights
             * @param aAverageEnrichedTmatrix output the averaged T-matrix
             */
            void
            average_T_matrices(
                    Cell< Vertex_Enrichment* > const & aAverageTmatrices,
                    Cell< real > const &               aWeights,
                    Vertex_Enrichment&                 aAverageEnrichedTmatrix ) const;

            // ----------------------------------------------------------------------------------
            
            /**
             * @brief fill a enriched T-matrix with NANs to mark it as a dummy object
             * 
             * @param aDummyEnrichedTmatrix T-matrix to fill with dummy information
             */
            void
            fill_T_matrix_dummy( Vertex_Enrichment& aDummyEnrichedTmatrix ) const;

            // ----------------------------------------------------------------------------------

            /**
             * @brief find and collect the Subphases in the coarsest B-spline mesh that are associated 
             * with the same material sub-domain as the enriched IP cell specified by the inputs
             * 
             * @param aIpCellIndex index of the base IP cell / Lagrange element
             * @param aUnzippingOnIpCell index of the unzipping (i.e. position of the enr. IP cell (UIPC) in the stack of UIPCs on this base IP cell)
             * @return Cell< moris_index > const& list of subphases
             */
            Cell< moris_index > const&
            collect_SPs_in_coarsest_element_in_same_material_subdomain_as_unzipped_cell(
                    const moris_index    aIpCellIndex,
                    const moris_index    aUnzippingIndex ) const;

            // ----------------------------------------------------------------------------------

            void
            get_averaged_vertex_enrichment_for_subphase(
                    const moris_index  aSubphaseIndex,
                    const moris_index  aBsplineMeshIndex,
                    Vertex_Enrichment* aAveragedTmatrix ) const;

            // ----------------------------------------------------------------------------------

            void
            construct_enriched_interpolation_vertices_and_cells();

            void
            construct_enriched_interpolation_vertices_and_cells_based_on_SPGs();

            void
            construct_enriched_interpolation_vertices_and_cells_based_on_SPGs_new();

            // ----------------------------------------------------------------------------------

            void
            communicate_unzipped_ip_cells();

            // ----------------------------------------------------------------------------------

            /**
             * @brief retrieve the unzipping for every enriched IP cell and store it in the member variable list
             */
            void
            construct_UIPC_to_unzipping_index();

            // ----------------------------------------------------------------------------------

            /**
             * @brief compute the maximum number of times an IP cell may need to be unzipped 
             * across all discretization mesh indices (DMIs). Recall that depending on the 
             * configuration of the SPGs void clusters may need to be created for the Ghost. The 
             * SPG layout depends on the DMI. Therefore the number of unzippings may also depend 
             * on the DMI. This function returns the maximum for a given IP cell.
             * 
             * @param aIpCellIndex index of the IP cell
             * @return uint maximum number of unzippings for that IP cell
             */
            uint
            maximum_number_of_unzippings_for_IP_cell( moris_index aIpCellIndex );

            // ----------------------------------------------------------------------------------

            /**
             * @brief finds SPs within an IP element that are also within the same SPG
             * 
             * @param aMeshIndex index of the B-spline mesh to be treated
             */
            void
            establish_IP_SPG_SP_relationship( const moris_index aMeshIndex );

            // ----------------------------------------------------------------------------------

            /**
             * @brief find out the enriched BF indices and IDs interpolating into the current vertex 
             * and store it along side the weights and other remaining interpolation information in 
             * the Vertex_Enrichment object
             * 
             * @param aEnrichmentDataIndex 
             * @param aBaseVertexInterp 
             * @param aSubPhaseBasisEnrLev 
             * @param aMapBasisIndexToLocInSubPhase 
             * @param aVertexEnrichment output: Vertex_Enrichment
             */
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

            /**
             * @brief function that collects the appropriate vertex interpolation for the enrichment strategy only
             * 
             * @param aParentCell 
             * @param aMeshIndex 
             * @return moris::Cell< mtk::Vertex_Interpolation* > 
             */
            moris::Cell< mtk::Vertex_Interpolation* >
            get_vertex_interpolations(
                    moris::mtk::Cell& aParentCell,
                    const uint        aMeshIndex ) const;

            // ----------------------------------------------------------------------------------

            /**
             * @brief ???
             * 
             *   @param[in] aMeshIndex Mesh index
             *   @param[in] aNumIdsToAllocate Number of basis ids needed
             *   @return Basis index offset
             */
            moris_index
            allocate_basis_ids( moris_index const& aMeshIndex,
                    moris_index const&                 aNumIdsToAllocate );

            // ----------------------------------------------------------------------------------

            /**
             * @brief ???
             *   
             * @param[in] aMeshIndex Mesh index
             * @return Maxmimum basis id on this proc
             */
            moris_index
            get_max_basis_id( moris_index const& aMeshIndex );

            // ----------------------------------------------------------------------------------            
    };
}// namespace xtk
#endif /* XTK_SRC_XTK_CL_XTK_ENRICHMENT_HPP_ */
