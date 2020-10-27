// MRS
#include "fn_Parsing_Tools.hpp"

// GEN
#include "cl_GEN_Geometry_Engine.hpp"
#include "fn_GEN_create_geometries.hpp"
#include "cl_GEN_Level_Set.hpp"
#include "fn_GEN_create_properties.hpp"
#include "cl_GEN_Interpolation.hpp"
#include "cl_GEN_Child_Node.hpp"
#include "cl_GEN_Intersection_Node.hpp"

// MTK
#include "cl_MTK_Integration_Mesh.hpp"
#include "cl_MTK_Interpolation_Mesh.hpp"
#include "cl_MTK_Writer_Exodus.hpp"

// XTK FIXME
#include "cl_XTK_Topology.hpp"

// SOL FIXME
#include "cl_SOL_Matrix_Vector_Factory.hpp"
#include "cl_SOL_Dist_Map.hpp"

namespace moris
{
    namespace ge
    {

        //--------------------------------------------------------------------------------------------------------------
        // PUBLIC
        //--------------------------------------------------------------------------------------------------------------

        Geometry_Engine::Geometry_Engine(
                Cell<Cell<ParameterList>> aParameterLists,
                std::shared_ptr<Library_IO> aLibrary):

                // Level set options
                mIsocontourThreshold(aParameterLists(0)(0).get<real>("isocontour_threshold")),
                mErrorFactor(aParameterLists(0)(0).get<real>("isocontour_error_factor")),

                // ADVs/IQIs
                mADVs(aParameterLists(0)(0).get<sint>("advs_size")
                        ? Matrix<DDRMat>((uint)aParameterLists(0)(0).get<sint>("advs_size"), 1, aParameterLists(0)(0).get<real>("initial_advs_fill"))
                        : string_to_mat<DDRMat>(aParameterLists(0)(0).get<std::string>("initial_advs"))),
                mLowerBounds(aParameterLists(0)(0).get<sint>("advs_size")
                        ? Matrix<DDRMat>((uint)aParameterLists(0)(0).get<sint>("advs_size"), 1, aParameterLists(0)(0).get<real>("lower_bounds_fill"))
                        : string_to_mat<DDRMat>(aParameterLists(0)(0).get<std::string>("lower_bounds"))),
                mUpperBounds(aParameterLists(0)(0).get<sint>("advs_size")
                        ? Matrix<DDRMat>((uint)aParameterLists(0)(0).get<sint>("advs_size"), 1, aParameterLists(0)(0).get<real>("upper_bounds_fill"))
                        : string_to_mat<DDRMat>(aParameterLists(0)(0).get<std::string>("upper_bounds"))),
                mRequestedIQIs(string_to_cell<std::string>(aParameterLists(0)(0).get<std::string>("IQI_types"))),

                // Library
                mLibrary(aLibrary),

                // Geometries
                mGeometries(create_geometries(aParameterLists(1), mADVs, mLibrary)),
                mGeometryParameterLists(aParameterLists(1)),
                mGeometryFieldFile(aParameterLists(0)(0).get<std::string>("geometry_field_file")),
                mOutputMeshFile(aParameterLists(0)(0).get<std::string>("output_mesh_file")),

                // Properties
                mProperties(create_properties(aParameterLists(2), mADVs, mGeometries, mLibrary)),
                mPropertyParameterLists(aParameterLists(2)),
                
                // phase table
                mPhaseTable(mGeometries.size())
        {
            // Set requested PDVs
            Cell<std::string> tRequestedPdvNames = string_to_cell<std::string>(aParameterLists(0)(0).get<std::string>("PDV_types"));
            Cell<PDV_Type> tRequestedPdvTypes(tRequestedPdvNames.size());
            map<std::string, PDV_Type> tPdvTypeMap = get_pdv_type_map();
            for (uint tPdvTypeIndex = 0; tPdvTypeIndex < tRequestedPdvTypes.size(); tPdvTypeIndex++)
            {
                tRequestedPdvTypes(tPdvTypeIndex) = tPdvTypeMap[tRequestedPdvNames(tPdvTypeIndex)];
            }
            mPdvHostManager.set_requested_interpolation_pdv_types(tRequestedPdvTypes);

            // Initialize PDV type list
            this->initialize_pdv_type_list();
            
            // Set map if its a non-standard phase table (i.e. the map is not 1-1 between index and bulk phase).
            if (aParameterLists(0)(0).get<std::string>("phase_table").length() > 0)
            {
                Matrix<IndexMat> tPhaseTable = string_to_mat<IndexMat>(aParameterLists(0)(0).get<std::string>("phase_table"));
                mPhaseTable.set_index_to_bulk_phase_map(tPhaseTable);
            }

            // print the phase table if requested (since GEN doesn't have a perform operator this is going here)
            if (aParameterLists(0)(0).get<bool>("print_phase_table") and par_rank() == 0)
            {
                mPhaseTable.print();
            }
        }

        //--------------------------------------------------------------------------------------------------------------

        Geometry_Engine::Geometry_Engine(
                Cell< std::shared_ptr<Geometry> > aGeometry,
                Phase_Table                       aPhaseTable,
                mtk::Interpolation_Mesh*          aMesh,
                Matrix<DDRMat>                    aADVs,
                real                              aIsocontourThreshold,
                real                              aErrorFactor)
                : mIsocontourThreshold(aIsocontourThreshold),
                  mErrorFactor(aErrorFactor),
                  mADVs(aADVs),
                  mGeometries(aGeometry),
                  mPhaseTable(aPhaseTable)
        {
            this->compute_level_set_data(aMesh);
        }

        //--------------------------------------------------------------------------------------------------------------

        Geometry_Engine::~Geometry_Engine()
        {
            delete mOwnedADVs;
        }

        //--------------------------------------------------------------------------------------------------------------

        void Geometry_Engine::set_advs(Matrix<DDRMat> aNewADVs)
        {
            // Set new ADVs
            mOwnedADVs->vec_put_scalar(0);
            mOwnedADVs->replace_global_values(mFullADVIds, aNewADVs);
            mOwnedADVs->vector_global_asembly();

            // Reset info related to the mesh
            mPdvHostManager.reset();
            mActiveGeometryIndex = 0;
            for (uint tGeometryIndex = 0; tGeometryIndex < mGeometries.size(); tGeometryIndex++)
            {
                mGeometries(tGeometryIndex)->import_advs(mOwnedADVs);
                mGeometries(tGeometryIndex)->reset_child_nodes();
            }
        }

        //--------------------------------------------------------------------------------------------------------------

        Matrix<DDRMat>& Geometry_Engine::get_advs()
        {
            // Create full ADVs
            sol::Matrix_Vector_Factory tDistributedFactory;
            std::shared_ptr<sol::Dist_Map> tFullMap = tDistributedFactory.create_map(mFullADVIds);
            moris::sol::Dist_Vector* tFullVector = tDistributedFactory.create_vector(tFullMap);

            // Import ADVs
            tFullVector->import_local_to_global(*mOwnedADVs);

            // Extract copy
            tFullVector->extract_copy(mADVs);

            // Delete full ADVs/map
            delete tFullVector;

            return mADVs;
        }

        //--------------------------------------------------------------------------------------------------------------

        Matrix<DDRMat>& Geometry_Engine::get_lower_bounds()
        {
            print(mLowerBounds, "lower");
            return mLowerBounds;
        }

        //--------------------------------------------------------------------------------------------------------------

        Matrix<DDRMat>& Geometry_Engine::get_upper_bounds()
        {
            print(mUpperBounds, "upper");
            return mUpperBounds;
        }

        //--------------------------------------------------------------------------------------------------------------

        void Geometry_Engine::communicate_requested_IQIs()
        {
            mPdvHostManager.set_requested_IQIs(mRequestedIQIs);
        }

        //--------------------------------------------------------------------------------------------------------------

        void Geometry_Engine::communicate_requested_IQIs(Cell<std::string> aIQINames)
        {
            mPdvHostManager.set_requested_IQIs(aIQINames);
        }

        //--------------------------------------------------------------------------------------------------------------

        Matrix<DDRMat> Geometry_Engine::get_dcriteria_dadv()
        {
            return mPdvHostManager.compute_diqi_dadv(mFullADVIds);
        }

        //--------------------------------------------------------------------------------------------------------------

        MSI::Design_Variable_Interface* Geometry_Engine::get_design_variable_interface()
        {
            return &mPdvHostManager;
        }

        //--------------------------------------------------------------------------------------------------------------

        real Geometry_Engine::get_geometry_field_value(
                uint                   aNodeIndex,
                const Matrix<DDRMat> & aCoordinates,
                uint                   aGeometryIndex)
        {
            return mGeometries(aGeometryIndex)->get_field_value(aNodeIndex, aCoordinates);
        }

        //--------------------------------------------------------------------------------------------------------------

        bool Geometry_Engine::is_intersected(const Matrix<IndexMat>& aNodeIndices, const Matrix<DDRMat>& aNodeCoordinates)
        {
            // Check input
            MORIS_ASSERT(aNodeIndices.length() == aNodeCoordinates.n_rows(),
                    "Geometry engine must be provided the same number of node indices as node coordinates for "
                    "determining if an element is intersected or not.");
            MORIS_ASSERT(aNodeIndices.length() > 0,
                    "Geometry engine must be provided at least 1 node to determine if an element is intersected or not.");

            // Initialize by evaluating the first node
            real tMin = mGeometries(mActiveGeometryIndex)->get_field_value(aNodeIndices(0), aNodeCoordinates.get_row(0));
            real tMax = tMin;

            // Evaluate the rest of the nodes
            for (uint tNodeCount = 0; tNodeCount < aNodeIndices.length(); tNodeCount++)
            {
                real tEval = mGeometries(mActiveGeometryIndex)->get_field_value(
                        aNodeIndices(tNodeCount),
                        aNodeCoordinates.get_row(tNodeCount));
                tMin = std::min(tMin, tEval);
                tMax = std::max(tMax, tEval);
            }

            // Return result
            return (tMax >= mIsocontourThreshold and tMin <= mIsocontourThreshold);
        }

        //--------------------------------------------------------------------------------------------------------------

        bool Geometry_Engine::queue_intersection(
                uint aFirstNodeIndex,
                uint aSecondNodeIndex,
                const Matrix<DDRMat>& aFirstNodeCoordinates,
                const Matrix<DDRMat>& aSecondNodeCoordinates)
        {
            // Determine if edge is intersected
            bool tEdgeIsIntersected = mGeometries(mActiveGeometryIndex)->get_field_value(aFirstNodeIndex, aFirstNodeCoordinates)
                    * mGeometries(mActiveGeometryIndex)->get_field_value(aSecondNodeIndex, aSecondNodeCoordinates) <= 0;

            // If edge is intersected, queue intersection node
            if (tEdgeIsIntersected)
            {
                mQueuedIntersectionNode = std::make_shared<Intersection_Node>(
                        aFirstNodeIndex,
                        aSecondNodeIndex,
                        aFirstNodeCoordinates,
                        aSecondNodeCoordinates,
                        mGeometries(mActiveGeometryIndex),
                        mIsocontourThreshold);
            }

            return tEdgeIsIntersected;
        }

        //--------------------------------------------------------------------------------------------------------------

        bool Geometry_Engine::queued_intersection_first_parent_on_interface()
        {
            return mQueuedIntersectionNode->first_parent_on_interface();
        }

        //--------------------------------------------------------------------------------------------------------------

        bool Geometry_Engine::queued_intersection_second_parent_on_interface()
        {
            return mQueuedIntersectionNode->second_parent_on_interface();
        }

        //--------------------------------------------------------------------------------------------------------------

        real Geometry_Engine::get_queued_intersection_local_coordinate()
        {
            return mQueuedIntersectionNode->get_local_coordinates()(0);
        }

        //--------------------------------------------------------------------------------------------------------------

        Matrix<DDRMat> Geometry_Engine::get_queued_intersection_global_coordinates()
        {
            return mQueuedIntersectionNode->get_global_coordinates();
        }

        //--------------------------------------------------------------------------------------------------------------

        void Geometry_Engine::admit_queued_intersection(uint aNodeIndex)
        {
            // Assign as PDV host
            if (mGeometries(mActiveGeometryIndex)->depends_on_advs())
            {
                mPdvHostManager.set_intersection_node(aNodeIndex, mQueuedIntersectionNode);
            }
            else
            {
                mPdvHostManager.set_intersection_node(aNodeIndex, nullptr);
            }

            // Assign as child node
            for (uint tGeometryIndex = 0; tGeometryIndex < mGeometries.size(); tGeometryIndex++)
            {
                mGeometries(tGeometryIndex)->add_child_node(aNodeIndex, mQueuedIntersectionNode);
            }
        }

        //--------------------------------------------------------------------------------------------------------------

        void Geometry_Engine::update_queued_intersection(
                const moris_index & aNodeIndex,
                const moris_index & aNodeId,
                const moris_index & aNodeOwner )
        {
            mPdvHostManager.update_intersection_node(aNodeIndex, aNodeId, aNodeOwner);
        }

        //--------------------------------------------------------------------------------------------------------------

        void Geometry_Engine::create_new_child_nodes(
                const Cell<moris_index>&    aNewNodeIndices,
                const Cell<xtk::Topology*>& aParentTopo,
                const Cell<Matrix<DDRMat>>& aParamCoordRelativeToParent,
                const Matrix<DDRMat>&       aGlobalNodeCoord )
        {
            // Loop over nodes
            for (uint tNode = 0; tNode < aNewNodeIndices.size(); tNode++)
            {
                // Create child node
                Matrix<DDUMat> tParentNodeIndices(aParentTopo(tNode)->get_node_indices().length(), 1);
                Cell<Matrix<DDRMat>> tParentNodeCoordinates(tParentNodeIndices.length());
                for (uint tParentNode = 0; tParentNode < tParentNodeIndices.length(); tParentNode++)
                {
                    tParentNodeIndices(tParentNode) = aParentTopo(tNode)->get_node_indices()(tParentNode);
                    tParentNodeCoordinates(tParentNode) = aGlobalNodeCoord.get_row(tParentNodeIndices(tParentNode));
                }
                std::shared_ptr<Child_Node> tChildNode = std::make_shared<Child_Node>(
                        tParentNodeIndices, tParentNodeCoordinates, aParentTopo(tNode)->get_basis_function(), aParamCoordRelativeToParent(tNode));

                // Assign to geometries
                for (uint tGeometryIndex = 0; tGeometryIndex < mGeometries.size(); tGeometryIndex++)
                {
                    mGeometries(tGeometryIndex)->add_child_node(aNewNodeIndices(tNode), tChildNode);
                }
            }

            // Set max node index
            if (aNewNodeIndices.size() > 0)
            {
                mPdvHostManager.set_num_background_nodes(aNewNodeIndices(aNewNodeIndices.size() - 1) + 1);
            }
        }

        //--------------------------------------------------------------------------------------------------------------

        size_t Geometry_Engine::get_num_phases()
        {
            return mPhaseTable.get_num_phases();
        }

        //--------------------------------------------------------------------------------------------------------------

        moris_index Geometry_Engine::get_phase_sign_of_given_phase_and_geometry(
                moris_index aPhaseIndex,
                moris_index aGeometryIndex )
        {
            return mPhaseTable.get_phase_sign_of_given_phase_and_geometry( aPhaseIndex,aGeometryIndex );
        }

        //--------------------------------------------------------------------------------------------------------------

        size_t Geometry_Engine::get_phase_index(
                moris_index            aNodeIndex,
                const Matrix<DDRMat> & aCoordinates)
        {
            // 0 for neg 1 for pos
            real tNodePhaseValue = 0;
            Matrix< IndexMat > tPhaseOnOff(1, this->get_num_geometries());

            for (uint tGeometryIndex = 0; tGeometryIndex < mGeometries.size(); tGeometryIndex++)
            {
                tNodePhaseValue = this->get_geometry_field_value(aNodeIndex, aCoordinates, tGeometryIndex);

                // Negative
                if (tNodePhaseValue < mIsocontourThreshold)
                {
                    tPhaseOnOff(0, tGeometryIndex) = 0;
                }
                else
                {
                    tPhaseOnOff(0, tGeometryIndex) = 1;
                }
            }

            return mPhaseTable.get_phase_index(tPhaseOnOff);
        }

        //--------------------------------------------------------------------------------------------------------------

        moris_index Geometry_Engine::get_elem_phase_index(Matrix< IndexMat > const & aElemOnOff)
        {
            return mPhaseTable.get_phase_index(aElemOnOff);
        }

        //--------------------------------------------------------------------------------------------------------------

        size_t Geometry_Engine::get_node_phase_index_wrt_a_geometry(
                uint                   aNodeIndex,
                const Matrix<DDRMat> & aCoordinates,
                uint                   aGeometryIndex)
        {
            real tNodePhaseValue = this->get_geometry_field_value(
                    aNodeIndex,
                    aCoordinates,
                    aGeometryIndex);

            size_t tPhaseOnOff = 1;

            if (tNodePhaseValue < mIsocontourThreshold)
            {
                tPhaseOnOff = 0;
            }

            return tPhaseOnOff;
        }

        //--------------------------------------------------------------------------------------------------------------

        size_t Geometry_Engine::get_num_geometries()
        {
            return mGeometries.size();
        }

        //--------------------------------------------------------------------------------------------------------------

        size_t Geometry_Engine::get_num_bulk_phase()
        {
            return mPhaseTable.get_num_phases();
        }

        //--------------------------------------------------------------------------------------------------------------

        size_t Geometry_Engine::get_active_geometry_index()
        {
            return mActiveGeometryIndex;
        }

        //--------------------------------------------------------------------------------------------------------------

        void Geometry_Engine::advance_geometry_index()
        {
            MORIS_ASSERT(mActiveGeometryIndex < mGeometries.size(),
                    "Trying to advance past the number of geometries in the geometry engine");
            mActiveGeometryIndex += 1;
        }

        //--------------------------------------------------------------------------------------------------------------

        uint Geometry_Engine::get_num_refinement_fields()
        {
            return mGeometries.size();
        }

        //--------------------------------------------------------------------------------------------------------------

        bool Geometry_Engine::refinement_needed(
                uint aFieldIndex,
                uint aRefinementIndex)
        {
            return ((sint)aRefinementIndex < mGeometries(aFieldIndex)->get_num_refinements());
        }

        //--------------------------------------------------------------------------------------------------------------

        real Geometry_Engine::get_field_value(
                uint                  aFieldIndex,
                uint                  aNodeIndex,
                const Matrix<DDRMat>& aCoordinates)
        {
            // TODO can return property field too
            return mGeometries(aFieldIndex)->get_field_value(aNodeIndex, aCoordinates);
        }

        //--------------------------------------------------------------------------------------------------------------

        sint Geometry_Engine::get_refinement_function_index(
                uint aFieldIndex,
                uint aRefinementIndex)
        {
            return mGeometries(aFieldIndex)->get_refinement_function_index();
        }

        //--------------------------------------------------------------------------------------------------------------

        void Geometry_Engine::create_pdvs(std::shared_ptr<mtk::Mesh_Manager> aMeshManager)
        {
            // Get integration mesh
            mtk::Integration_Mesh* tIntegrationMesh = aMeshManager->get_integration_mesh(0);

            // Initialize PDV type groups and mesh set info
            Cell<Cell<Cell<PDV_Type>>> tPdvTypes(tIntegrationMesh->get_num_sets());
            Cell<PDV_Type> tPdvTypeGroup(1);

            Cell<std::string> tMeshSetNames(0);
            Matrix<DDUMat> tMeshSetIndices(0, 0);

            // PDV type map
            map< std::string, PDV_Type > tPdvTypeMap = get_pdv_type_map();

            // Loop over properties to create PDVs
            for (uint tPropertyIndex = 0; tPropertyIndex < mPropertyParameterLists.size(); tPropertyIndex++)
            {
                // PDV type and mesh set names/indices from parameter list
                tPdvTypeGroup(0) = tPdvTypeMap[mPropertyParameterLists(tPropertyIndex).get<std::string>("pdv_type")];

                string_to_cell(mPropertyParameterLists(tPropertyIndex).get<std::string>("pdv_mesh_set_names"), tMeshSetNames);
                string_to_mat(mPropertyParameterLists(tPropertyIndex).get<std::string>("pdv_mesh_set_indices"), tMeshSetIndices);

                // Convert mesh set names to indices
                uint tNumSetIndices = tMeshSetIndices.length();
                tMeshSetIndices.resize(tNumSetIndices + tMeshSetNames.size(), 1);

                for (uint tIndex = tNumSetIndices; tIndex < tMeshSetIndices.length(); tIndex++)
                {
                    tMeshSetIndices(tIndex) = tIntegrationMesh->get_set_index_by_name(tMeshSetNames(tIndex - tNumSetIndices));
                }

                // Assign PDV types
                for (uint tIndex = 0; tIndex < tMeshSetIndices.length(); tIndex++)
                {
                    tPdvTypes(tMeshSetIndices(tIndex)).push_back(tPdvTypeGroup);
                }
            }

            Matrix< IdMat > tCommTable = aMeshManager->get_interpolation_mesh(0)->get_communication_table();
            std::unordered_map<moris_id,moris_index> tIPVertexGlobaToLocalMap =
                    aMeshManager->get_interpolation_mesh(0)->get_vertex_glb_id_to_loc_vertex_ind_map();
            std::unordered_map<moris_id,moris_index> tIGVertexGlobaToLocalMap =
                    tIntegrationMesh->get_vertex_glb_id_to_loc_vertex_ind_map();

            mPdvHostManager.set_communication_table( tCommTable );
            mPdvHostManager.set_vertex_global_to_local_maps( tIPVertexGlobaToLocalMap, tIGVertexGlobaToLocalMap);

            // Create PDV hosts
            this->create_interpolation_pdv_hosts(
                    aMeshManager->get_interpolation_mesh(0),
                    tIntegrationMesh,
                    tPdvTypes);

            // Set integration PDV types
            if (mShapeSensitivities)
            {
                this->set_integration_pdv_types(tIntegrationMesh);
            }

            // Loop over properties to assign PDVs
            for (uint tPropertyIndex = 0; tPropertyIndex < mPropertyParameterLists.size(); tPropertyIndex++)
            {
                // Assign PDVs
                if (mPropertyParameterLists(tPropertyIndex).get<std::string>("pdv_mesh_type") == "interpolation")
                {
                    this->assign_property_to_pdv_hosts(mProperties(tPropertyIndex), tPdvTypeGroup(0), tIntegrationMesh, tMeshSetIndices);
                }
                else
                {
                    MORIS_ERROR(false, "Assignment of PDVs is only supported with an interpolation mesh right now.");
                }
            }

            // Create PDV IDs
            mPdvHostManager.create_pdv_ids();

            // Get rid of parameter lists
            mPropertyParameterLists.resize(0);
        }

        //--------------------------------------------------------------------------------------------------------------

        void Geometry_Engine::initialize_pdv_type_list()
        {
            // Reserve of temporary pdv type list
            moris::Cell< enum PDV_Type > tTemporaryPdvTypeList;
            tTemporaryPdvTypeList.reserve( static_cast< int >( PDV_Type::UNDEFINED ) + 1 );
            Matrix< DDUMat > tListToCheckIfEnumExist( (static_cast< int >(PDV_Type::UNDEFINED) + 1), 1, 0 );

            // PDV type map
            map< std::string, PDV_Type > tPdvTypeMap = get_pdv_type_map();

            // Loop over properties to build parallel consitent pdv list
             for (uint tPropertyIndex = 0; tPropertyIndex < mPropertyParameterLists.size(); tPropertyIndex++)
             {
                 // PDV type and mesh set names/indices from parameter list
                 enum PDV_Type tPdvType = tPdvTypeMap[mPropertyParameterLists(tPropertyIndex).get<std::string>("pdv_type")];

                 if ( tListToCheckIfEnumExist(static_cast< int >(tPdvType) , 0) == 0)
                 {
                     // Set 1 at position of the enum value
                     tListToCheckIfEnumExist( static_cast< int >(tPdvType) ,0 ) = 1;

                     tTemporaryPdvTypeList.push_back( tPdvType );
                 }
             }

            // Shrink pdv type list to fit
             tTemporaryPdvTypeList.shrink_to_fit();

            // Communicate dof types so that all processors have the same unique list
             mPdvHostManager.communicate_dof_types( tTemporaryPdvTypeList );

            // Create a map
             mPdvHostManager.create_dv_type_map();
        }

        //--------------------------------------------------------------------------------------------------------------

        void Geometry_Engine::compute_level_set_data(mtk::Interpolation_Mesh* aMesh)
        {
            // Register spatial dimension
            mSpatialDim = aMesh->get_spatial_dim();

            // Set number of background nodes
            mPdvHostManager.set_num_background_nodes(aMesh->get_num_nodes());

            // Build distributed ADVs
            if (mOwnedADVs == nullptr)
            {
                //------------------------------------//
                // Determine owned and shared ADV IDs //
                //------------------------------------//

                // Loop over all geometries to get number of new ADVs
                uint tNumNewOwnedADVs = 0;
                for (uint tGeometryIndex = 0; tGeometryIndex < mGeometries.size(); tGeometryIndex++)
                {
                    // Determine if level set will be created
                    if (mGeometries(tGeometryIndex)->conversion_to_bsplines())
                    {
                        // Get number of coefficients
                        sint tBSplineMeshIndex = mGeometries(tGeometryIndex)->get_bspline_mesh_index();
                        uint tNumCoefficients = aMesh->get_num_coeffs(tBSplineMeshIndex);

                        // Loop over B-spline coefficients
                        for (uint tBSplineIndex = 0; tBSplineIndex < tNumCoefficients; tBSplineIndex++)
                        {
                            // If this processor owns this coefficient
                            if (par_rank() == aMesh->get_entity_owner(tBSplineIndex, EntityRank::BSPLINE, tBSplineMeshIndex))
                            {
                                // New ADV needs to be created on this processor
                                tNumNewOwnedADVs++;
                            }
                        }
                    }
                }

                // Set number of owned ADVs
                uint tNumOwnedADVs = mADVs.length() + tNumNewOwnedADVs;
                mLowerBounds.resize(tNumOwnedADVs, 1);
                mUpperBounds.resize(tNumOwnedADVs, 1);

                // Resize ADV IDs and set primitive IDs
                Matrix<DDSMat> tOwnedADVIds(mADVs.length(), 1);
                for (uint tADVIndex = 0; tADVIndex < mADVs.length(); tADVIndex++)
                {
                    tOwnedADVIds(tADVIndex) = tADVIndex;
                }
                Matrix<DDSMat> tPrimitiveIds = tOwnedADVIds;
                tOwnedADVIds.resize(tNumOwnedADVs, 1);

                // Cell of shared ADV IDs
                Cell<Matrix<DDSMat>> tSharedADVIds(mGeometries.size());

                // Loop over all geometry parameter lists to set B-spline ADV bounds and IDs
                uint tOwnedADVIndex = mADVs.length();
                uint tIdOffset = tOwnedADVIndex; // FIXME this needs to be updated for multiple level set fields
                for (uint tGeometryIndex = 0; tGeometryIndex < mGeometries.size(); tGeometryIndex++)
                {
                    // Determine if level set will be created
                    if (mGeometries(tGeometryIndex)->conversion_to_bsplines())
                    {
                        // Get bounds
                        real tBSplineLowerBound = mGeometries(tGeometryIndex)->get_bspline_lower_bound();
                        real tBSplineUpperBound = mGeometries(tGeometryIndex)->get_bspline_upper_bound();

                        // Get number of coefficients
                        sint tBSplineMeshIndex = mGeometries(tGeometryIndex)->get_bspline_mesh_index();
                        uint tNumCoefficients = aMesh->get_num_coeffs(tBSplineMeshIndex);

                        // Resize shared ADV IDs
                        tSharedADVIds(tGeometryIndex).resize(tNumCoefficients, 1);

                        // Loop over B-spline coefficients
                        for (uint tBSplineIndex = 0; tBSplineIndex < tNumCoefficients; tBSplineIndex++)
                        {
                            // Get new ADV ID and set as shared ID
                            sint tNewADVId = tIdOffset + aMesh->get_glb_entity_id_from_entity_loc_index(
                                    tBSplineIndex,
                                    EntityRank::BSPLINE,
                                    tBSplineMeshIndex);
                            tSharedADVIds(tGeometryIndex)(tBSplineIndex) = tNewADVId;

                            // If this processor owns this coefficient set to owned list and set bounds
                            if (par_rank() == aMesh->get_entity_owner(tBSplineIndex, EntityRank::BSPLINE, tBSplineMeshIndex))
                            {
                                // Bounds
                                mLowerBounds(tOwnedADVIndex) = tBSplineLowerBound;
                                mUpperBounds(tOwnedADVIndex) = tBSplineUpperBound;

                                // New ADV ID
                                tOwnedADVIds(tOwnedADVIndex++) = tNewADVId;
                            }
                        }
                    }
                }

                // Set owned ADV IDs
                mPdvHostManager.set_owned_adv_ids(tOwnedADVIds);

                //----------------------------------------//
                // Create owned ADV vector                //
                //----------------------------------------//

                // Create factory for distributed ADV vector
                sol::Matrix_Vector_Factory tDistributedFactory;

                // Create map for distributed vector
                std::shared_ptr<sol::Dist_Map> tOwnedADVMap = tDistributedFactory.create_map(tOwnedADVIds);

                // Create vector
                mOwnedADVs = tDistributedFactory.create_vector(tOwnedADVMap);

                // Assign primitive IDs
                mOwnedADVs->replace_global_values(tPrimitiveIds, mADVs);

                // Global assembly
                mOwnedADVs->vector_global_asembly();

                // Build geometries and properties from parameter lists using distributed vector
                // TODO augmented copy constructor for fields
                if (mGeometryParameterLists.size() > 0)
                {
                    // Build geometries and properties
                    mGeometries = create_geometries(mGeometryParameterLists, mOwnedADVs, mLibrary);
                    mProperties = create_properties(mPropertyParameterLists, mOwnedADVs, mGeometries, mLibrary);
                    mGeometryParameterLists.resize(0);
                }

                //----------------------------------------//
                // Convert geometries to level sets       //
                //----------------------------------------//

                // Determine if conversion to level sets are needed and if shape sensitivities are needed
                for (uint tGeometryIndex = 0; tGeometryIndex < mGeometries.size(); tGeometryIndex++)
                {
                    // Shape sensitivities logic in case of no level sets
                    mShapeSensitivities = (mShapeSensitivities or mGeometries(tGeometryIndex)->depends_on_advs());

                    // Convert to level set if needed
                    if (mGeometries(tGeometryIndex)->conversion_to_bsplines())
                    {
                        // Always have shape sensitivities if level set
                        mShapeSensitivities = true;

                        // Create level set FIXME for multiple level set fields
                        mGeometries(tGeometryIndex) = std::make_shared<Level_Set>(
                                mOwnedADVs,
                                tOwnedADVIds,
                                tSharedADVIds(tGeometryIndex),
                                tPrimitiveIds.length(),
                                aMesh,
                                mGeometries(tGeometryIndex));
                    }
                }

                //----------------------------------------//
                // Communicate all ADV IDs to processor 0 //
                //----------------------------------------//

                // Sending mats
                Cell<Matrix<DDSMat>> tSendingIDs(0);
                Cell<Matrix<DDRMat>> tSendingLowerBounds(0);
                Cell<Matrix<DDRMat>> tSendingUpperBounds(0);

                // Receiving mats
                Cell<Matrix<DDSMat>> tReceivingIDs(0);
                Cell<Matrix<DDRMat>> tReceivingLowerBounds(0);
                Cell<Matrix<DDRMat>> tReceivingUpperBounds(0);

                // Set up communication list for communicating ADV IDs
                Matrix<IdMat> tCommunicationList(1, 1, 0);
                if (par_rank() == 0)
                {
                    // Resize communication list and sending mats
                    tCommunicationList.resize(par_size() - 1, 1);
                    tSendingIDs.resize(par_size() - 1);
                    tSendingLowerBounds.resize(par_size() - 1);
                    tSendingUpperBounds.resize(par_size() - 1);

                    // Assign communication list
                    for (uint tProcessorIndex = 1; tProcessorIndex < (uint)par_size(); tProcessorIndex++)
                    {
                        tCommunicationList(tProcessorIndex - 1) = tProcessorIndex;
                    }
                }
                else
                {
                    tSendingIDs = {tOwnedADVIds};
                    tSendingLowerBounds = {mLowerBounds};
                    tSendingUpperBounds = {mUpperBounds};
                }

                // Communicate mats
                communicate_mats(tCommunicationList, tSendingIDs, tReceivingIDs);
                communicate_mats(tCommunicationList, tSendingLowerBounds, tReceivingLowerBounds);
                communicate_mats(tCommunicationList, tSendingUpperBounds, tReceivingUpperBounds);

                // Assemble full ADVs/bounds
                if (par_rank() == 0)
                {
                    // Start full IDs with owned IDs on processor 0
                    mFullADVIds = tOwnedADVIds;

                    // Assemble additional IDs/bounds from other processors
                    for (uint tProcessorIndex = 1; tProcessorIndex < (uint)par_size(); tProcessorIndex++)
                    {
                        // Get number of received ADVs
                        uint tFullADVsFilled = mFullADVIds.length();
                        uint tNumReceivedADVs = tReceivingIDs(tProcessorIndex - 1).length();

                        // Resize full ADV IDs and bounds
                        mFullADVIds.resize(tFullADVsFilled + tNumReceivedADVs, 1);
                        mLowerBounds.resize(tFullADVsFilled + tNumReceivedADVs, 1);
                        mUpperBounds.resize(tFullADVsFilled + tNumReceivedADVs, 1);

                        // Assign received ADV IDs
                        for (uint tADVIndex = 0; tADVIndex < tNumReceivedADVs; tADVIndex++)
                        {
                            mFullADVIds(tFullADVsFilled + tADVIndex) = tReceivingIDs(tProcessorIndex - 1)(tADVIndex);
                            mLowerBounds(tFullADVsFilled + tADVIndex) =
                                    tReceivingLowerBounds(tProcessorIndex - 1)(tADVIndex);
                            mUpperBounds(tFullADVsFilled + tADVIndex) =
                                    tReceivingUpperBounds(tProcessorIndex - 1)(tADVIndex);
                        }
                    }
                }
                else
                {
                    mLowerBounds.set_size(0, 0);
                    mUpperBounds.set_size(0, 0);
                }
            }
        }

        //--------------------------------------------------------------------------------------------------------------

        void Geometry_Engine::output_fields(mtk::Mesh* aMesh)
        {
            this->output_fields_on_mesh(aMesh, mOutputMeshFile);
            this->write_geometry_fields(aMesh, mGeometryFieldFile);
        }

        //--------------------------------------------------------------------------------------------------------------

        void Geometry_Engine::output_fields_on_mesh(mtk::Mesh* aMesh, std::string aExodusFileName)
        {
            if (aExodusFileName != "")
            {
                // Write mesh
                mtk::Writer_Exodus tWriter(aMesh);
                tWriter.write_mesh("", aExodusFileName, "", "gen_temp.exo");

                // Setup fields
                Cell<std::string> tFieldNames(mGeometries.size());
                for (uint tGeometryIndex = 0; tGeometryIndex < mGeometries.size(); tGeometryIndex++)
                {
                    tFieldNames(tGeometryIndex) = "Geometry " + std::to_string(tGeometryIndex);
                }
                tWriter.set_nodal_fields(tFieldNames);

                // Get all node coordinates
                Cell<Matrix<DDRMat>> tNodeCoordinates(aMesh->get_num_nodes());
                for (uint tNodeIndex = 0; tNodeIndex < aMesh->get_num_nodes(); tNodeIndex++)
                {
                    tNodeCoordinates(tNodeIndex) = aMesh->get_node_coordinate(tNodeIndex);
                }

                // Loop over geometries
                for (uint tGeometryIndex = 0; tGeometryIndex < mGeometries.size(); tGeometryIndex++)
                {
                    // Create field vector
                    Matrix<DDRMat> tFieldData(aMesh->get_num_nodes(), 1);

                    // Assign field to vector
                    for (uint tNodeIndex = 0; tNodeIndex < aMesh->get_num_nodes(); tNodeIndex++)
                    {
                        tFieldData(tNodeIndex) = mGeometries(tGeometryIndex)->get_field_value(
                                tNodeIndex,
                                tNodeCoordinates(tNodeIndex));
                    }

                    // Create field on mesh
                    tWriter.write_nodal_field("Geometry " + std::to_string(tGeometryIndex), tFieldData);
                }

                // Finalize
                tWriter.close_file(true);
            }
        }

        //--------------------------------------------------------------------------------------------------------------

        void Geometry_Engine::write_geometry_fields(mtk::Mesh* aMesh, std::string aBaseFileName)
        {
            if (aBaseFileName != "")
            {
                // Get all node coordinates
                Cell<Matrix<DDRMat>> tNodeCoordinates(aMesh->get_num_nodes());
                for (uint tNodeIndex = 0; tNodeIndex < aMesh->get_num_nodes(); tNodeIndex++)
                {
                    tNodeCoordinates(tNodeIndex) = aMesh->get_node_coordinate(tNodeIndex);
                }

                // Loop over geometries
                for (uint tGeometryIndex = 0; tGeometryIndex < mGeometries.size(); tGeometryIndex++)
                {
                    // Create file
                    std::ofstream tOutFile(aBaseFileName + "_" + std::to_string(tGeometryIndex) + ".txt");

                    // Write to file
                    for (uint tNodeIndex = 0; tNodeIndex < aMesh->get_num_nodes(); tNodeIndex++)
                    {
                        // Coordinates
                        for (uint tDimension = 0; tDimension < mSpatialDim; tDimension++)
                        {
                            tOutFile << tNodeCoordinates(tNodeIndex)(tDimension) << ", ";
                        }

                        // Fill unused dimensions with zeros
                        for (uint tDimension = mSpatialDim; tDimension < 3; tDimension++)
                        {
                            tOutFile << 0.0 << ", ";
                        }

                        // Level-set field
                        tOutFile << mGeometries(tGeometryIndex)->get_field_value(
                                tNodeIndex,
                                tNodeCoordinates(tNodeIndex)) << std::endl;
                    }

                    // Close file
                    tOutFile.close();
                }
            }
        }

        //--------------------------------------------------------------------------------------------------------------
        // PRIVATE
        //--------------------------------------------------------------------------------------------------------------

        void Geometry_Engine::create_interpolation_pdv_hosts(
                mtk::Interpolation_Mesh     * aInterpolationMesh,
                mtk::Integration_Mesh       * aIntegrationMesh,
                Cell<Cell<Cell<PDV_Type>>>    aPdvTypes)
        {
            // Get information from integration mesh
            //uint tNumSets = tInterpolationMesh->get_num_sets(); FIXME
            uint tNumSets  = aPdvTypes.size();
            uint tNumNodes = aInterpolationMesh->get_num_nodes();

            Cell<Matrix<DDSMat>> tNodeIndicesPerSet(tNumSets);
            Cell<Matrix<DDSMat>> tNodeIdsPerSet(tNumSets);
            Cell<Matrix<DDSMat>> tNodeOwnersPerSet(tNumSets);
            Cell<Matrix<DDRMat>> tNodeCoordinates(tNumNodes);

            // Loop through sets
            for (uint tMeshSetIndex = 0; tMeshSetIndex < tNumSets; tMeshSetIndex++)
            {
                uint tCurrentNode = 0;
                mtk::Set* tSet = aIntegrationMesh->get_set_by_index(tMeshSetIndex);

                // Clusters per set
                for (uint tClusterIndex = 0; tClusterIndex < tSet->get_num_clusters_on_set(); tClusterIndex++)
                {
                    const mtk::Cluster* tCluster = tSet->get_clusters_by_index(tClusterIndex);

                    // Indices on cluster // FIXME this is really bad and slow. especially when building the pdvs
                    Matrix<IndexMat> tNodeIndicesInCluster = tCluster->get_interpolation_cell().get_vertex_inds();
                    Matrix<IndexMat> tNodeIdsInCluster     = tCluster->get_interpolation_cell().get_vertex_ids();
                    Matrix<IndexMat> tNodeOwnersInCluster = tCluster->get_interpolation_cell().get_vertex_owners();

                    // FIXME don't undersand this resize. it's really slow
                    tNodeIndicesPerSet(tMeshSetIndex).resize(tNodeIndicesPerSet(tMeshSetIndex).length() + tNodeIndicesInCluster.length(), 1);
                    tNodeIdsPerSet(tMeshSetIndex).resize(tNodeIdsPerSet(tMeshSetIndex).length() + tNodeIdsInCluster.length(), 1);
                    tNodeOwnersPerSet(tMeshSetIndex).resize(tNodeOwnersPerSet(tMeshSetIndex).length() + tNodeOwnersInCluster.length(), 1);

                    // FIXME we have nodes up to 8 tims in this list in 3d
                    for (uint tNodeInCluster = 0; tNodeInCluster < tNodeIndicesInCluster.length(); tNodeInCluster++)
                    {
                        tNodeIndicesPerSet(tMeshSetIndex)(tCurrentNode)   = tNodeIndicesInCluster(tNodeInCluster);
                        tNodeIdsPerSet(tMeshSetIndex)(tCurrentNode)   = tNodeIdsInCluster(tNodeInCluster);
                        tNodeOwnersPerSet(tMeshSetIndex)(tCurrentNode++) = tNodeOwnersInCluster(tNodeInCluster);
                    }
                }
            }

            // Get node coordinates
            for (uint tNodeIndex = 0; tNodeIndex < tNumNodes; tNodeIndex++)
            {
                tNodeCoordinates(tNodeIndex) = aInterpolationMesh->get_node_coordinate(tNodeIndex);
            }

            // Create hosts
            mPdvHostManager.create_interpolation_pdv_hosts(
                    tNodeIndicesPerSet,
                    tNodeIdsPerSet,
                    tNodeOwnersPerSet,
                    tNodeCoordinates,
                    aPdvTypes);
        }

        //--------------------------------------------------------------------------------------------------------------

        void Geometry_Engine::set_integration_pdv_types(mtk::Integration_Mesh* aIntegrationMesh)
        {
            // Get information from integration mesh
            uint tNumSets = aIntegrationMesh->get_num_sets();

            // Cell of IG PDV_Type types
            Cell<PDV_Type> tCoordinatePdvs(mSpatialDim);

            switch(mSpatialDim)
            {
                case 2:
                {
                    tCoordinatePdvs(0) = PDV_Type::X_COORDINATE;
                    tCoordinatePdvs(1) = PDV_Type::Y_COORDINATE;
                    break;
                }
                case 3:
                {
                    tCoordinatePdvs(0) = PDV_Type::X_COORDINATE;
                    tCoordinatePdvs(1) = PDV_Type::Y_COORDINATE;
                    tCoordinatePdvs(2) = PDV_Type::Z_COORDINATE;
                    break;
                }
                default:
                {
                    MORIS_ERROR( false, "Geometry Engine only works for 2D and 3D models." );
                }
            }

            // Loop through sets
            Cell<Cell<Cell<PDV_Type>>> tPdvTypes(tNumSets);
            for (uint tMeshSetIndex = 0; tMeshSetIndex < tNumSets; tMeshSetIndex++)
            {
                // PDV_Type types per set
                tPdvTypes(tMeshSetIndex).resize(1);
                tPdvTypes(tMeshSetIndex)(0) = tCoordinatePdvs;
            }

            // Set PDV types
            mPdvHostManager.set_integration_pdv_types(tPdvTypes);
            mPdvHostManager.set_requested_integration_pdv_types(tCoordinatePdvs);

        }

        //--------------------------------------------------------------------------------------------------------------

        void Geometry_Engine::assign_property_to_pdv_hosts(
                std::shared_ptr<Property> aPropertyPointer,
                PDV_Type                  aPdvType,
                mtk::Integration_Mesh*    aIntegrationMesh,
                Matrix<DDUMat>            aSetIndices)
        {
            for (uint tSet = 0; tSet < aSetIndices.length(); tSet++)
            {
                // get the mesh set from index
                mtk::Set* tSetPointer = aIntegrationMesh->get_set_by_index( aSetIndices(tSet) );

                // get the list of cluster on mesh set
                Cell< mtk::Cluster const * > tClusterPointers = tSetPointer->get_clusters_on_set();

                // get number of clusters on mesh set
                uint tNumClusters = tClusterPointers.size();

                // loop over the clusters on mesh set
                for(uint iClust=0; iClust<tNumClusters; iClust++)
                {
                    // get the IP cell from cluster
                    mtk::Cell const & tIPCell = tClusterPointers(iClust)->get_interpolation_cell();

                    // get the vertices from IP cell
                    Cell< mtk::Vertex * > tVertices = tIPCell.get_vertex_pointers();

                    // get the number of vertices on IP cell
                    uint tNumVerts = tVertices.size();

                    // loop over vertices on IP cell
                    for(uint iVert = 0; iVert < tNumVerts; iVert++)
                    {
                        // get the vertex index
                        moris_index tVertIndex = tVertices(iVert)->get_index();

                        // ask pdv host manager to assign to vertex a pdv type and a property
                        mPdvHostManager.create_interpolation_pdv( uint(tVertIndex), aPdvType, aPropertyPointer );
                    }
                }
            }
        }

        //--------------------------------------------------------------------------------------------------------------

    }
}
