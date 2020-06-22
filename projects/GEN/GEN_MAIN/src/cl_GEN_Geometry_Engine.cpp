// GEN
#include "cl_GEN_Geometry_Engine.hpp"
#include "fn_GEN_Matrix_Base_Utilities.hpp"
#include "fn_GEN_create_geometry.hpp"
#include "fn_GEN_create_properties.hpp"
#include "cl_GEN_Interpolation.hpp"
#include "cl_GEN_Child_Node.hpp"

// LINALG
#include "cl_Matrix.hpp"
#include "fn_trans.hpp"
#include "op_equal_equal.hpp"
#include "op_times.hpp"
#include "linalg_typedefs.hpp"

// HMR
#include "cl_HMR_Mesh.hpp"
#include "cl_HMR.hpp"

// PRM
#include "cl_PRM_HMR_Parameters.hpp"

// MRS/IOS
#include "fn_Parsing_Tools.hpp"

namespace moris
{
    namespace ge
    {
        //--------------------------------------------------------------------------------------------------------------

        Geometry_Engine::Geometry_Engine(Cell<Cell<ParameterList>> aParameterLists,
                                         std::shared_ptr<moris::Library_IO> aLibrary):
                // User options
                mThresholdValue(aParameterLists(0)(0).get<real>("threshold_value")),
                mPerturbationValue(aParameterLists(0)(0).get<real>("perturbation_value")),
                mNumRefinements(aParameterLists(0)(0).get<int>("HMR_refinements")),
                mUserDefinedFunc(nullptr),

                // ADVs/IQIs
                mADVs((uint)aParameterLists(0)(0).get<int>("advs_size"), 1, aParameterLists(0)(0).get<real>("initial_advs_fill")),
                mLowerBounds((uint)aParameterLists(0)(0).get<int>("advs_size"), 1, aParameterLists(0)(0).get<real>("lower_bounds_fill")),
                mUpperBounds((uint)aParameterLists(0)(0).get<int>("advs_size"), 1, aParameterLists(0)(0).get<real>("upper_bounds_fill")),
                mRequestedIQIs(string_to_cell<std::string>(aParameterLists(0)(0).get<std::string>("IQI_types"))),

                // Phase table
                mPhaseTable(string_to_mat<IndexMat>(aParameterLists(0)(0).get<std::string>("phase_table")).numel()
                      ? Phase_Table(string_to_mat<IndexMat>(aParameterLists(0)(0).get<std::string>("phase_table")), aParameterLists(0)(0).get<std::string>("phase_table_structure"))
                      : Phase_Table(aParameterLists(1).size(), aParameterLists(0)(0).get<std::string>("phase_table_structure")))
        {
            // Explicit ADVs and bounds
            Matrix<DDRMat> tInitialADVs = string_to_mat<DDRMat>(aParameterLists(0)(0).get<std::string>("initial_advs"));
            Matrix<DDRMat> tLowerBounds = string_to_mat<DDRMat>(aParameterLists(0)(0).get<std::string>("lower_bounds"));
            Matrix<DDRMat> tUpperBounds = string_to_mat<DDRMat>(aParameterLists(0)(0).get<std::string>("upper_bounds"));

            // Resize if needed
            if (tInitialADVs.length() > mADVs.length())
            {
                mADVs.resize(tInitialADVs.length(), 1);
                mLowerBounds.resize(tInitialADVs.length(), 1);
                mUpperBounds.resize(tInitialADVs.length(), 1);
            }

            // Copy over values
            for (uint tADVIndex = 0; tADVIndex < tInitialADVs.length(); tADVIndex++)
            {
                mADVs(tADVIndex) = tInitialADVs(tADVIndex);
            }
            for (uint tADVIndex = 0; tADVIndex < tLowerBounds.length(); tADVIndex++)
            {
                mLowerBounds(tADVIndex) = tLowerBounds(tADVIndex);
            }
            for (uint tADVIndex = 0; tADVIndex < tUpperBounds.length(); tADVIndex++)
            {
                mUpperBounds(tADVIndex) = tUpperBounds(tADVIndex);
            }

            // Build geometry (just analytic for right now)
            if (aParameterLists(1).size() > 0)
            {
                mGeometry.resize(aParameterLists(1).size());
                for (uint tGeometryIndex = 0; tGeometryIndex < aParameterLists(1).size(); tGeometryIndex++)
                {
                    mGeometry(tGeometryIndex) = create_geometry(aParameterLists(1)(tGeometryIndex), mADVs, aLibrary);
                }
            }

            // Create function pointer for user defined refinement function
            std::string tUserDefinedFunc = aParameterLists(0)(0).get<std::string>("user_defined_refinement_function");

            if ( tUserDefinedFunc.size() > 1 )
            {
                mUserDefinedFunc = aLibrary->load_user_defined_refinement_functions( tUserDefinedFunc );
            }

            // Create properties
            mProperties = create_properties(aParameterLists(2), mADVs, aLibrary);
            mPropertyParameterLists = aParameterLists(2);

            // Set requested PDVs
            Cell<std::string> tRequestedPdvNames = string_to_cell<std::string>(aParameterLists(0)(0).get<std::string>("PDV_types"));
            Cell<PDV_Type> tRequestedPdvTypes(tRequestedPdvNames.size());
            moris::map<std::string, PDV_Type> tPdvTypeMap = get_pdv_type_map();
            for (uint tPdvTypeIndex = 0; tPdvTypeIndex < tRequestedPdvTypes.size(); tPdvTypeIndex++)
            {
                tRequestedPdvTypes(tPdvTypeIndex) = tPdvTypeMap[tRequestedPdvNames(tPdvTypeIndex)];
            }
            mPdvHostManager.set_ip_requested_dv_types(tRequestedPdvTypes);
        }

        //--------------------------------------------------------------------------------------------------------------

        Geometry_Engine::Geometry_Engine(Cell<std::shared_ptr<Geometry>>    aGeometry,
                                         Phase_Table                        aPhaseTable,
                                         uint                               aSpatialDim,
                                         real                               aThresholdValue,
                                         real                               aPerturbationValue)
            : mSpatialDim(aSpatialDim),
              mThresholdValue(aThresholdValue),
              mPerturbationValue(aPerturbationValue),
              mActiveGeometryIndex(0),
              mGeometry(aGeometry),
              mPhaseTable(aPhaseTable)
        {
        }

        //--------------------------------------------------------------------------------------------------------------

        Geometry_Engine::~Geometry_Engine()
        {
        }

        //--------------------------------------------------------------------------------------------------------------

        void Geometry_Engine::set_advs(Matrix<DDRMat> aNewADVs)
        {
            mADVs = aNewADVs;
        }

        //--------------------------------------------------------------------------------------------------------------

        Matrix<DDRMat>& Geometry_Engine::get_advs()
        {
            return mADVs;
        }

        //--------------------------------------------------------------------------------------------------------------

        Matrix<DDRMat>& Geometry_Engine::get_lower_bounds()
        {
            return mLowerBounds;
        }

        //--------------------------------------------------------------------------------------------------------------

        Matrix<DDRMat>& Geometry_Engine::get_upper_bounds()
        {
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
            return mPdvHostManager.compute_diqi_dadv();
        }

        //--------------------------------------------------------------------------------------------------------------

        MSI::Design_Variable_Interface* Geometry_Engine::get_design_variable_interface()
        {
            return &mPdvHostManager;
        }

        //--------------------------------------------------------------------------------------------------------------

        real Geometry_Engine::get_geometry_field_value(      uint            aNodeIndex,
                                                       const Matrix<DDRMat>& aCoordinates,
                                                             uint            aGeometryIndex)
        {
            return mGeometry(aGeometryIndex)->evaluate_field_value(aNodeIndex, aCoordinates);
        }

        //--------------------------------------------------------------------------------------------------------------

        void Geometry_Engine::create_new_node_geometry_objects( const Cell<moris_index>&    aNewNodeIndices,
                                                                bool                        aStoreParentTopo,
                                                                const Cell<xtk::Topology*>& aParentTopo,
                                                                const Cell<Matrix<DDRMat>>& aParamCoordRelativeToParent,
                                                                const Matrix<DDRMat>&       aGlobalNodeCoord )
        {
            // Allocate space
            moris::size_t tNumNewNodes = aNewNodeIndices.size();
            Cell<GEN_Geometry_Object> tGeometryObjects(tNumNewNodes);

            moris::Matrix< moris::IndexMat > tNodeIndices(1,tNumNewNodes);
            for (moris::size_t i = 0; i < tNumNewNodes; i++)
            {
                tGeometryObjects(i).mGeometryIndex = mActiveGeometryIndex;
                tNodeIndices(0, i) = aNewNodeIndices(i);
                if (aStoreParentTopo)
                {
                    tGeometryObjects(i).set_parent_entity_topology(aParentTopo(i)->copy());
                }
            }

            for (uint tNode = 0; tNode < tNumNewNodes; tNode++ )
            {
                // Create child node
                Matrix<DDUMat> tParentNodeIndices(aParentTopo(tNode)->get_node_indices().length(), 1);
                Cell<Matrix<DDRMat>> tParentNodeCoordinates(tParentNodeIndices.length());
                for (uint tParentNode = 0; tParentNode < tParentNodeIndices.length(); tParentNode++)
                {
                    tParentNodeIndices(tParentNode) = aParentTopo(tNode)->get_node_indices()(tParentNode);
                    tParentNodeCoordinates(tParentNode) = aGlobalNodeCoord.get_row(tParentNodeIndices(tParentNode));
                }
                Child_Node tChildNode(tParentNodeIndices, tParentNodeCoordinates, aParentTopo(tNode)->get_basis_function(), aParamCoordRelativeToParent(tNode));

                // Assign to geometries
                for (uint tGeometryIndex = 0; tGeometryIndex < this->get_num_geometries(); tGeometryIndex++)
                {
                    mGeometry(tGeometryIndex)->add_child_node(aNewNodeIndices(tNode), tChildNode);
                }
            }

            if (tNumNewNodes != 0)
            {
                mGeometryObjectManager.store_geometry_objects(tNodeIndices, tGeometryObjects);
            }
        }

        //--------------------------------------------------------------------------------------------------------------

        void Geometry_Engine::link_new_nodes_to_existing_geometry_objects( Matrix< IndexMat > const & aNodesIndicesWithGeomObj,
                                                                               Matrix< IndexMat > const & aNodesIndicesToLink )
        {
            // Assert lengths match
            MORIS_ASSERT(aNodesIndicesWithGeomObj.numel() == aNodesIndicesToLink.numel(),
            "Length of nodes with geometry objects does not match length of list with node indices to link  ");

            // Number of nodes to link
            uint tNumNodes = aNodesIndicesWithGeomObj.numel();

            // Iterate through nodes and create the link
            for(uint i = 0; i <tNumNodes; i++)
            {
                mGeometryObjectManager.link_to_node_to_another_nodes_geometry_object(aNodesIndicesWithGeomObj(i),aNodesIndicesToLink(i));
            }

        }

        //--------------------------------------------------------------------------------------------------------------

        void Geometry_Engine::is_intersected( const Matrix<DDRMat>&      aNodeCoords,
                                              const Matrix<IndexMat>&    aNodetoEntityConn,
                                              moris::size_t              aCheckType,
                                              Cell<GEN_Geometry_Object>& aGeometryObjects )
        {
            //Get information for loops
            moris::size_t tNumEntities = aNodetoEntityConn.n_rows(); // Number of entities provided to the geometry engine

            //Initialize
            moris::size_t tIntersectedCount = 0;    // Intersected element counter
            aGeometryObjects.clear();
            aGeometryObjects.resize(tNumEntities,GEN_Geometry_Object());

            //Loop over elements and determine if the element has an intersection
            for(moris::moris_index i = 0; i < (moris::moris_index)tNumEntities; i++)
            {

                //Populate the intersection flag of this element with a bool
                moris::Matrix< moris::IndexMat > tRow = aNodetoEntityConn.get_row(i);
                moris::Matrix< moris::IndexMat > tNodeADVIndices;
                bool tIsIntersected = compute_intersection_info( i,tRow, aNodeCoords, aCheckType, tNodeADVIndices, aGeometryObjects(tIntersectedCount) );

                if(tIsIntersected)
                {
                    tIntersectedCount++;
                }
            }

            // resize
            aGeometryObjects.resize( tIntersectedCount, GEN_Geometry_Object() );
        }

        //--------------------------------------------------------------------------------------------------------------

        void Geometry_Engine::set_interface_nodes( Matrix< IndexMat > const & aInterfaceNodeIndices)
        {
            mInterfaceNodeIndices = aInterfaceNodeIndices;
        }

        //--------------------------------------------------------------------------------------------------------------

        void Geometry_Engine::get_intersection_location( real                    aIsocontourThreshold,
                                                         real                    aPerturbationThreshold,
                                                         const Matrix<DDRMat>&   aGlobalNodeCoordinates,
                                                         const Matrix<DDRMat>&   aEntityNodeVars,
                                                         const Matrix<IndexMat>& aEntityNodeIndices,
                                                         Matrix<DDRMat>&         aIntersectionLocalCoordinates,
                                                         Matrix<DDRMat>&         aIntersectionGlobalCoordinates,
                                                         bool                    aCheckLocalCoordinate,
                                                         bool                    aComputeGlobalCoordinate )
        {

            // compute the local coordinate where the intersection occurs
            Interpolation::linear_interpolation_value(aEntityNodeVars, aIsocontourThreshold, aIntersectionLocalCoordinates);

            // Perturb away from node if necessary
            if(aCheckLocalCoordinate)
            {
                if(aIntersectionLocalCoordinates(0, 0) >= 1-aPerturbationThreshold)
                {
                    aIntersectionLocalCoordinates(0, 0) = aIntersectionLocalCoordinates(0, 0) - aPerturbationThreshold;
                }

                if(aIntersectionLocalCoordinates(0, 0) <= -1+aPerturbationThreshold)
                {
                    aIntersectionLocalCoordinates(0, 0) = aIntersectionLocalCoordinates(0, 0) + aPerturbationThreshold;
                }
            }

            // Compute the global coordinate only if you plan to use it
            if(aComputeGlobalCoordinate)
            {
                // Place only the entity coordinates in a matrix
                moris::Matrix< moris::DDRMat > tEntityCoordinates(2,mSpatialDim);
                replace_row(aEntityNodeIndices(0,0), aGlobalNodeCoordinates,0,tEntityCoordinates);
                replace_row(aEntityNodeIndices(0,1), aGlobalNodeCoordinates,1,tEntityCoordinates);

                // compute the global coordinate
                Interpolation::linear_interpolation_location(tEntityCoordinates,aIntersectionLocalCoordinates,aIntersectionGlobalCoordinates);
            }
        }

        //--------------------------------------------------------------------------------------------------------------

        moris::size_t Geometry_Engine::get_num_phases()
        {
            return mPhaseTable.get_num_phases();
        }

        //--------------------------------------------------------------------------------------------------------------

        moris::moris_index
        Geometry_Engine::get_phase_sign_of_given_phase_and_geometry( moris::moris_index aPhaseIndex,
                                                                         moris::moris_index aGeometryIndex )
        {
            return mPhaseTable.get_phase_sign_of_given_phase_and_geometry( aPhaseIndex,aGeometryIndex );
        }

        //--------------------------------------------------------------------------------------------------------------

        size_t Geometry_Engine::get_phase_index(moris_index aNodeIndex, const Matrix<DDRMat>& aCoordinates)
        {
            // 0 for neg 1 for pos
            moris::real tNodePhaseValue = 0;
            moris::Matrix< moris::IndexMat > tPhaseOnOff(1, this->get_num_geometries());

            for (uint tGeometryIndex = 0; tGeometryIndex < mGeometry.size(); tGeometryIndex++)
            {
                tNodePhaseValue = this->get_geometry_field_value(aNodeIndex, aCoordinates, tGeometryIndex);

                // Negative
                if (tNodePhaseValue < mThresholdValue)
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

        moris_index Geometry_Engine::get_elem_phase_index(moris::Matrix< moris::IndexMat > const & aElemOnOff)
        {
            return mPhaseTable.get_phase_index(aElemOnOff);
        }

        //--------------------------------------------------------------------------------------------------------------

        size_t Geometry_Engine::get_node_phase_index_wrt_a_geometry(      uint            aNodeIndex,
                                                                    const Matrix<DDRMat>& aCoordinates,
                                                                          uint            aGeometryIndex)
        {
            real tNodePhaseValue = this->get_geometry_field_value(aNodeIndex, aCoordinates, aGeometryIndex);

            moris::size_t tPhaseOnOff = 1;
            if (tNodePhaseValue < mThresholdValue)
            {
                tPhaseOnOff = 0;
            }

            return tPhaseOnOff;
        }

        //--------------------------------------------------------------------------------------------------------------

        moris::size_t Geometry_Engine::get_num_geometries()
        {
            return mGeometry.size();
        }

        //--------------------------------------------------------------------------------------------------------------

        moris::size_t Geometry_Engine::get_num_bulk_phase()
        {
            return mPhaseTable.get_num_phases();
        }

        //--------------------------------------------------------------------------------------------------------------

        moris::size_t Geometry_Engine::get_active_geometry_index()
        {
            return mActiveGeometryIndex;
        }

        //--------------------------------------------------------------------------------------------------------------

        void Geometry_Engine::advance_geometry_index()
        {
            MORIS_ASSERT(mActiveGeometryIndex < mGeometry.size(),
                    "Trying to advance past the number of geometries in the geometry engine");
            mActiveGeometryIndex += 1;
        }

        //--------------------------------------------------------------------------------------------------------------

        void Geometry_Engine::register_mesh(std::shared_ptr<mtk::Mesh_Manager> aMeshManager)
        {
            mMeshManager = aMeshManager;
            mSpatialDim = mMeshManager->get_interpolation_mesh(0)->get_spatial_dim();
        }

        //--------------------------------------------------------------------------------------------------------------

        void Geometry_Engine::perform_refinement(std::shared_ptr<hmr::HMR> aHMRPerformer, uint aNumRefinements)
        {
            // Create mesh
            std::shared_ptr<hmr::Mesh> tMesh = aHMRPerformer->create_mesh(0);

            // Determine number of refinements
            if (aNumRefinements == 0)
            {
                aNumRefinements = mNumRefinements;
            }

            // FIXME
            ParameterList tParameters = prm::create_hmr_parameter_list();

            // Loop over set number of refinement levels
            for (uint tRefinement = 0; tRefinement < aNumRefinements; tRefinement++)
            {
                // Loop over geometries to get field values and put on queue
                for (uint tGeometryIndex = 0; tGeometryIndex < mGeometry.size(); tGeometryIndex++)
                {
                    Matrix<DDRMat> tFieldValues(tMesh->get_num_nodes(), 1);
                    for (uint tNodeIndex = 0; tNodeIndex < tMesh->get_num_nodes(); tNodeIndex++)
                    {
                        tFieldValues(tNodeIndex) = this->get_geometry_field_value(tNodeIndex, tMesh->get_node_coordinate(tNodeIndex), tGeometryIndex);
                    }

                    // Call either user defined refinement function or default function
                    if (mUserDefinedFunc != nullptr )
                    {
                        Cell<Matrix<DDRMat>> tCellFieldValues({tFieldValues});
                        aHMRPerformer->user_defined_flagging( 0, tCellFieldValues, tParameters, 0 );
                    }
                    else
                    {
                        aHMRPerformer->based_on_field_put_elements_on_queue(tFieldValues, 0);
                    }
                }

                // Perform refinement
                aHMRPerformer->perform_refinement_based_on_working_pattern( 0, false );
            }
        }

        //--------------------------------------------------------------------------------------------------------------

        bool Geometry_Engine::compute_intersection_info( moris::moris_index               const & aEntityIndex,
                                                         moris::Matrix< moris::IndexMat > const & aEntityNodeInds,
                                                         moris::Matrix< moris::DDRMat >   const & aNodeCoords,
                                                         moris::size_t const &                    aCheckType,
                                                         moris::Matrix< moris::IndexMat > &       aNodeADVIndices,
                                                         GEN_Geometry_Object & aGeometryObject )
        {
            //Initialize
            bool tIsIntersected = false;

            moris::real tMax = 0;
            moris::real tMin = 0;
            moris::uint tMaxLocRow = 0;
            moris::uint tMaxLocCol = 0;
            moris::uint tMinLocRow = 0;
            moris::uint tMinLocCol = 0;

            moris::size_t tNodeInd  = 0;
            moris::size_t tNumNodes = aEntityNodeInds.numel();
            moris::Matrix< moris::DDRMat > tEntityNodeVars(tNumNodes, 1);

            // Loop through nodes and get levelset values from precomputed values in aNodeVars or in the levelset mesh
            for (moris::size_t n = 0; n < tNumNodes; n++)
            {
                tNodeInd = aEntityNodeInds(n);

                tEntityNodeVars(n) = this->get_geometry_field_value(tNodeInd, aNodeCoords.get_row(tNodeInd), mActiveGeometryIndex); //FIXME Wrong
            }

            //get the max and minimum levelset value for the entity
            tMax = tEntityNodeVars.max(tMaxLocRow,tMaxLocCol);
            tMin = tEntityNodeVars.min(tMinLocRow,tMinLocCol);

            //    If there is a sign change in element node variables return true, else return false

            //TODO: intersection flag should not be a moris::real (needs to be a bool) split this function
            moris::Matrix< moris::DDRMat > tIntersection(1, 2, 0.0);// Initialize as false

            moris::real tErrorFactor = 1;
            // If the max is also the threshold value, figure out which node is on the interface

            if( moris::ge::approximate(tMin, mThresholdValue, tErrorFactor) && moris::ge::approximate(tMax, mThresholdValue,tErrorFactor))
            {
                aGeometryObject.set_parent_entity_index(aEntityIndex);
                aGeometryObject.mark_all_nodes_as_on_interface();
                tIsIntersected = true;
            }

            else if(moris::ge::approximate(tMax,mThresholdValue, tErrorFactor))
            {
                aGeometryObject.set_parent_entity_index(aEntityIndex);
                aGeometryObject.mark_node_as_on_interface(tMaxLocRow);
                tIsIntersected = true;
            }

            // If the min is also the threshold value, figure out which node is on the interface
            else if(moris::ge::approximate(tMin,mThresholdValue, tErrorFactor))
            {
                aGeometryObject.set_parent_entity_index(aEntityIndex);
                aGeometryObject.mark_node_as_on_interface(tMinLocRow);
                tIsIntersected = true;
            }

            else if((tMax > mThresholdValue) &&
               (tMin < mThresholdValue))
            {
                aGeometryObject.set_parent_entity_index(aEntityIndex);
                aGeometryObject.mark_nodes_as_not_on_interface();
                tIsIntersected = true;
                if(aCheckType == 1)
                {
                    moris::Matrix< moris::DDRMat > tIntersectLocalCoordinate(1,1);
                    moris::Matrix< moris::DDRMat > tIntersectGlobalCoordinate(1,mSpatialDim);

                    get_intersection_location(mThresholdValue,
                                              mPerturbationValue,
                                              aNodeCoords,
                                              tEntityNodeVars,
                                              aEntityNodeInds,
                                              tIntersectLocalCoordinate,
                                              tIntersectGlobalCoordinate,
                                              true,
                                              true);

                    aGeometryObject.set_interface_loc_coord(tIntersectLocalCoordinate(0));
                    aGeometryObject.set_interface_glb_coord(tIntersectGlobalCoordinate);
               }
            }

            return tIsIntersected;

        }

        //--------------------------------------------------------------------------------------------------------------

        void Geometry_Engine::assign_ip_hosts_by_set_name(std::string                 aSetName,
                                                          std::shared_ptr< Property > aPropertyPointer,
                                                          PDV_Type                    aPdvType,
                                                          moris_index                 aWhichMesh)
        {
            // get the mesh set from name
            moris::mtk::Set* tSetPointer = mMeshManager->get_integration_mesh( aWhichMesh )->get_set_by_name( aSetName );

            // get the list of cluster on mesh set
            moris::Cell< mtk::Cluster const * > tClusterPointers = tSetPointer->get_clusters_on_set();

            // get number of clusters on mesh set
            uint tNumClusters = tClusterPointers.size();

            // loop over the clusters on mesh set
            for(uint iClust=0; iClust<tNumClusters; iClust++)
            {
                // get the IP cell from cluster
                moris::mtk::Cell const & tIPCell = tClusterPointers(iClust)->get_interpolation_cell();

                // get the vertices from IP cell
                moris::Cell< moris::mtk::Vertex * > tVertices = tIPCell.get_vertex_pointers();

                // get the number of vertices on IP cell
                uint tNumVerts = tVertices.size();

                // loop over vertices on IP cell
                for(uint iVert = 0; iVert < tNumVerts; iVert++)
                {
                    // get the vertex index
                    moris_index tVertIndex = tVertices(iVert)->get_index();

                    // ask pdv host manager to assign to vertex a pdv type and a property
                    mPdvHostManager.create_ip_pdv( uint(tVertIndex), aPdvType, aPropertyPointer);
                }
            }
        }

        //--------------------------------------------------------------------------------------------------------------

        void Geometry_Engine::assign_ip_hosts_by_set_index( moris_index                 aSetIndex,
                                                            std::shared_ptr< Property > aPropertyPointer,
                                                            PDV_Type                    aPdvType,
                                                            moris_index                 aWhichMesh)
        {
            // get the mesh set from index
            moris::mtk::Set* tSetPointer = mMeshManager->get_integration_mesh( aWhichMesh )->get_set_by_index( aSetIndex );

            // get the list of cluster on mesh set
            moris::Cell< mtk::Cluster const * > tClusterPointers = tSetPointer->get_clusters_on_set();

            // get number of clusters on mesh set
            uint tNumClusters = tClusterPointers.size();

            // loop over the clusters on mesh set
            for(uint iClust=0; iClust<tNumClusters; iClust++)
            {
                // get the IP cell from cluster
                moris::mtk::Cell const & tIPCell = tClusterPointers(iClust)->get_interpolation_cell();

                // get the vertices from IP cell
                moris::Cell< moris::mtk::Vertex * > tVertices = tIPCell.get_vertex_pointers();

                // get the number of vertices on IP cell
                uint tNumVerts = tVertices.size();

                // loop over vertices on IP cell
                for(uint iVert = 0; iVert < tNumVerts; iVert++)
                {
                    // get the vertex index
                    moris_index tVertIndex = tVertices(iVert)->get_index();

                    // ask pdv host manager to assign to vertex a pdv type and a property
                    mPdvHostManager.create_ip_pdv( uint(tVertIndex), aPdvType, aPropertyPointer );
                }
            }
        }

        //--------------------------------------------------------------------------------------------------------------

        void Geometry_Engine::create_ip_pdv_hosts(Cell<Cell<Cell<PDV_Type>>> aPdvTypes, moris_index aMeshIndex)
        {
            // Get information from integration mesh
            mtk::Integration_Mesh* tIntegrationMesh = mMeshManager->get_integration_mesh(aMeshIndex);
            //uint tNumSets = tInterpolationMesh->get_num_sets(); FIXME
            uint tNumSets = aPdvTypes.size();
            uint tNumNodes = mMeshManager->get_interpolation_mesh(aMeshIndex)->get_num_nodes();
            Cell<Matrix<DDSMat>> tNodeIndicesPerSet(tNumSets);
            Cell<Matrix<F31RMat>> tNodeCoordinates(tNumNodes);

            // Loop through sets
            Cell<Cell<Cell<PDV_Type>>> tPdvTypes(tNumSets);
            mtk::Set* tSet;
            const mtk::Cluster* tCluster;
            uint tCurrentNode;
            Matrix<IndexMat> tNodeIndicesInCluster;
            for (uint tMeshSetIndex = 0; tMeshSetIndex < tNumSets; tMeshSetIndex++)
            {
                tCurrentNode = 0;
                tSet = tIntegrationMesh->get_set_by_index(tMeshSetIndex);

                // Clusters per set
                for (uint tClusterIndex = 0; tClusterIndex < tSet->get_num_clusters_on_set(); tClusterIndex++)
                {
                    tCluster = tSet->get_clusters_by_index(tClusterIndex);

                    // Indices on cluster
                    tNodeIndicesInCluster = tCluster->get_interpolation_cell().get_vertex_inds();
                    tNodeIndicesPerSet(tMeshSetIndex).resize(tNodeIndicesPerSet(tMeshSetIndex).length() + tNodeIndicesInCluster.length(), 1);

                    for (uint tNodeInCluster = 0; tNodeInCluster < tNodeIndicesInCluster.length(); tNodeInCluster++)
                    {
                        tNodeIndicesPerSet(tMeshSetIndex)(tCurrentNode++) = tNodeIndicesInCluster(tNodeInCluster);
                    }
                }
            }

            // Get node coordinates
            for (uint tNodeIndex = 0; tNodeIndex < tNumNodes; tNodeIndex++)
            {
                tNodeCoordinates(tNodeIndex) = mMeshManager->get_interpolation_mesh(aMeshIndex)->get_node_coordinate(tNodeIndex);
            }

            // Create hosts
            mPdvHostManager.create_ip_pdv_hosts(tNodeIndicesPerSet, tNodeCoordinates, aPdvTypes);
        }

        //--------------------------------------------------------------------------------------------------------------

        void Geometry_Engine::create_ig_pdv_hosts(moris_index aMeshIndex)
        {
            // Get information from integration mesh
            mtk::Integration_Mesh* tIntegrationMesh = mMeshManager->get_integration_mesh(aMeshIndex);
            uint tNumSets = tIntegrationMesh->get_num_sets();
            uint tNumNodes = tIntegrationMesh->get_num_nodes();
            Cell<Matrix<DDSMat>> tNodeIndicesPerSet(tNumSets);
            Cell<Matrix<F31RMat>> tNodeCoordinates(tNumNodes);

            // Cell of IG PDV_Type types
            Cell<PDV_Type> tCoordinatePdvs(mSpatialDim);

            switch(mSpatialDim)
            {
                case(2):
                {
                    tCoordinatePdvs(0) = PDV_Type::X_COORDINATE;
                    tCoordinatePdvs(1) = PDV_Type::Y_COORDINATE;
                    break;
                }
                case(3):
                {
                    tCoordinatePdvs(0) = PDV_Type::X_COORDINATE;
                    tCoordinatePdvs(1) = PDV_Type::Y_COORDINATE;
                    tCoordinatePdvs(2) = PDV_Type::Z_COORDINATE;
                    break;
                }
                default:
                {
                    MORIS_ERROR( false, "Geometry_Engine::initialize_integ_pdv_host_list() - Geometry Engine only works for 2D and 3D models." );
                }
            }

            // Loop through sets
            Cell<Cell<Cell<PDV_Type>>> tPdvTypes(tNumSets);
            for (uint tMeshSetIndex = 0; tMeshSetIndex < tNumSets; tMeshSetIndex++)
            {
                // Node indices per set
                tNodeIndicesPerSet(tMeshSetIndex) = tIntegrationMesh->get_set_by_index(tMeshSetIndex)->get_vertices_inds_on_block(false);

                // PDV_Type types per set
                tPdvTypes(tMeshSetIndex).resize(1);
                tPdvTypes(tMeshSetIndex)(0) = tCoordinatePdvs;
            }

            // Get node coordinates
            for (uint tNodeIndex = 0; tNodeIndex < tNumNodes; tNodeIndex++)
            {
                tNodeCoordinates(tNodeIndex) = mMeshManager->get_integration_mesh(aMeshIndex)->get_node_coordinate(tNodeIndex);
            }

            // Get intersection dependencies
            Cell<Intersection_Info> tIntersectionInfo(mInterfaceNodeIndices.length());
            for (uint tNodeIndex = 0; tNodeIndex < mInterfaceNodeIndices.length(); tNodeIndex++)
            {
                // Get geometry object TODO will be removed in future
                GEN_Geometry_Object* tGeometryObject = mGeometryObjectManager.get_geometry_object( mInterfaceNodeIndices(tNodeIndex));
                xtk::Topology const & tParentEdge = tGeometryObject->get_parent_entity_topology();

                // Set info
                tIntersectionInfo(tNodeIndex).mGeometry = this->mGeometry(tGeometryObject->mGeometryIndex);
                tIntersectionInfo(tNodeIndex).mNodeIndex = mInterfaceNodeIndices(tNodeIndex);
                tIntersectionInfo(tNodeIndex).mParentNodeIndices = tParentEdge.get_node_indices();;
            }

            // Create hosts
            mPdvHostManager.create_ig_pdv_hosts(tNodeIndicesPerSet, tNodeCoordinates, tPdvTypes, tIntersectionInfo);
        }

        //--------------------------------------------------------------------------------------------------------------

        void Geometry_Engine::assign_pdv_hosts()
        {
            // Initialize
            mtk::Integration_Mesh* tIntegrationMesh = mMeshManager->get_integration_mesh(0);
            Cell<Cell<Cell<PDV_Type>>> tPdvTypes(tIntegrationMesh->get_num_sets());
            Cell<PDV_Type> tPdvTypeGroup(1);
            Cell<std::string> tMeshSetNames(0);
            Matrix<DDUMat> tMeshSetIndices(0, 0);

            // PDV type map
            moris::map< std::string, PDV_Type > tPdvTypeMap = get_pdv_type_map();

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
            this->create_ip_pdv_hosts(tPdvTypes);

            // Loop over properties to assign PDVs
            for (uint tPropertyIndex = 0; tPropertyIndex < mPropertyParameterLists.size(); tPropertyIndex++)
            {
                // PDV type and mesh set names/indices from parameter list
                tPdvTypeGroup(0) = tPdvTypeMap[mPropertyParameterLists(tPropertyIndex).get<std::string>("pdv_type")];
                string_to_cell(mPropertyParameterLists(tPropertyIndex).get<std::string>("pdv_mesh_set_names"), tMeshSetNames);
                string_to_mat(mPropertyParameterLists(tPropertyIndex).get<std::string>("pdv_mesh_set_indices"), tMeshSetIndices);

                // Assign PDVs
                if (mPropertyParameterLists(tPropertyIndex).get<std::string>("pdv_mesh_type") == "interpolation")
                {
                    // Set names
                    for (uint tNameIndex = 0; tNameIndex < tMeshSetNames.size(); tNameIndex++)
                    {
                        this->assign_ip_hosts_by_set_name(tMeshSetNames(tNameIndex), mProperties(tPropertyIndex), tPdvTypeGroup(0));
                    }

                    // Set indices
                    for (uint tIndex = 0; tIndex < tMeshSetIndices.length(); tIndex++)
                    {
                        this->assign_ip_hosts_by_set_index(tMeshSetIndices(tIndex), mProperties(tPropertyIndex), tPdvTypeGroup(0));
                    }
                }
                else
                {
                    MORIS_ERROR(false, "Assignment of PDVs is only supported with an interpolation mesh right now.");
                }
            }
        }

        //--------------------------------------------------------------------------------------------------------------

    }   // end ge namespace
}   // end moris namespace
