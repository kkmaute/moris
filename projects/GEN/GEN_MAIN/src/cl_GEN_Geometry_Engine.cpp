// GEN
#include "cl_GEN_Geometry_Engine.hpp"
#include "fn_GEN_Matrix_Base_Utilities.hpp"
#include "fn_GEN_create_geometry.hpp"
#include "fn_GEN_create_properties.hpp"

// LINALG
#include "cl_Matrix.hpp"
#include "fn_trans.hpp"
#include "op_equal_equal.hpp"
#include "op_times.hpp"
#include "linalg_typedefs.hpp"

// HMR
#include "cl_HMR_Mesh.hpp"
#include "cl_HMR.hpp"

// MRS/IOS
#include "fn_Parsing_Tools.hpp"

namespace moris
{
    namespace ge
    {
        //--------------------------------------------------------------------------------------------------------------
        
        Geometry_Engine::Geometry_Engine(moris::Cell<moris::Cell<ParameterList>> aParameterLists, std::shared_ptr<moris::Library_IO> aLibrary) :
        // User options
        mSpatialDim(aParameterLists(0)(0).get<int>("spatial_dimensions")),
        mThresholdValue(aParameterLists(0)(0).get<real>("threshold_value")),
        mPerturbationValue(aParameterLists(0)(0).get<real>("perturbation_value")),
        mNumRefinements(aParameterLists(0)(0).get<int>("HMR_refinements")),

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
                mADVs(tADVIndex) = tLowerBounds(tADVIndex);
            }
            for (uint tADVIndex = 0; tADVIndex < tUpperBounds.length(); tADVIndex++)
            {
                mADVs(tADVIndex) = tUpperBounds(tADVIndex);
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

        Matrix<DDRMat> Geometry_Engine::get_dcriteria_dadv() // TODO
        {
            return mPdvHostManager.compute_diqi_dadv();
        }

        //--------------------------------------------------------------------------------------------------------------

        MSI::Design_Variable_Interface* Geometry_Engine::get_design_variable_interface()
        {
            return &mPdvHostManager;
        }

        //--------------------------------------------------------------------------------------------------------------
        
        void Geometry_Engine::initialize_geometry_objects_for_background_mesh_nodes(moris::size_t const & aNumNodes)
        {
            // Allocate space
            mNodePhaseVals = moris::Matrix< moris::DDRMat >(aNumNodes, this->get_num_geometries(), 0);
        
            // Allocate geometry object
            Cell< GEN_Geometry_Object > tGeometryObjects( aNumNodes );
        
            // Associate each geometry object with a row in phase val matrix (note phase val computed later)
            moris::Matrix< moris::IndexMat > tNodeIndex( 1, aNumNodes );
            for(moris::size_t i = 0; i < aNumNodes; i++)
            {
                tGeometryObjects(i).set_phase_val_row(i);
                tNodeIndex(0, i) = i;
            }
        
            // Place these in the geometry object manager
            mGeometryObjectManager.store_geometry_objects(tNodeIndex, tGeometryObjects);
        }

        //--------------------------------------------------------------------------------------------------------------
        
        void Geometry_Engine::initialize_geometry_object_phase_values( moris::Matrix< moris::DDRMat > const & aNodeCoords )
        {
            // Allocate space
            size_t tNumNodes = aNodeCoords.n_rows();
            uint tNumGeometries = mGeometry.size();
        
            // Loop through each geometry and then each node and compute the level set field value
            for (moris::size_t tGeometryIndex = 0; tGeometryIndex < tNumGeometries; tGeometryIndex++) // Analytic
            {
                for (moris::size_t tNodeIndex = 0; tNodeIndex < tNumNodes; tNodeIndex++)
                {
                    mNodePhaseVals(tNodeIndex, tGeometryIndex) =
                            mGeometry(tGeometryIndex)->evaluate_field_value(tNodeIndex, aNodeCoords.get_row(tNodeIndex));
                }
            }
        }

        //--------------------------------------------------------------------------------------------------------------
        
        void Geometry_Engine::associate_new_nodes_with_geometry_object( Cell<Pending_Node> & aNewNodes,
                                                                            bool                 aInterfaceNodes )
        {
            // Allocate space
            moris::size_t tNumNewNodes = aNewNodes.size();
            moris::size_t tNumCurrNodes = mNodePhaseVals.n_rows();
        
            // add space to the node phase value table
            mNodePhaseVals.resize(tNumNewNodes+tNumCurrNodes, this->get_num_geometries());
        
            Cell<GEN_Geometry_Object> tGeometryObjects(tNumNewNodes);
        
            moris::Matrix< moris::IndexMat > tNodeIndex(1,tNumNewNodes);
            for(moris::size_t i = 0; i<tNumNewNodes; i++)
            {
                tGeometryObjects(i).set_phase_val_row(i+tNumCurrNodes);
                tNodeIndex(0,i) = aNewNodes(i).get_node_index();
                if(aInterfaceNodes)
                {
                    tGeometryObjects(i).set_parent_entity_topology(aNewNodes(i).get_parent_topology_ptr());
                }
            }
        
            if(tNumNewNodes !=0)
            {
                mGeometryObjectManager.store_geometry_objects(tNodeIndex,tGeometryObjects);
            }
        
            // Compute and store level set value of this node for each new node
            for(moris::size_t j = 0; j < this->get_num_geometries(); j++)
            {
        
                for(moris::size_t i = 0; i < tNumNewNodes; i++ )
                {
                    // Ask the pending node about its parent
                    // This information is needed to know what to interpolate based on
                    moris::Matrix< moris::DDRMat > const & tLocalCoordinate = aNewNodes(i).get_local_coordinate_relative_to_parent();
                    moris::Matrix< moris::DDRMat >  tLevelSetValues(1,1);
                    xtk::Topology const & tParentTopology = aNewNodes(i).get_parent_topology();
        
                    // Interpolate all level set values to node
                    this->interpolate_level_set_value_to_child_node_location(tParentTopology, j,tLocalCoordinate,tLevelSetValues);
                    mNodePhaseVals(i+tNumCurrNodes,j) = tLevelSetValues(0,0);
                }
        
            }
        }
        
        //--------------------------------------------------------------------------------------------------------------
        
        void Geometry_Engine::create_new_node_geometry_objects( Cell< moris_index >  const & aNewNodeIndices,
                                                                    bool                         aStoreParentTopo,
                                                                    Cell<xtk::Topology*> const & aParentTopo,
                                                                    Cell<Matrix<DDRMat>> const & aParamCoordRelativeToParent,
                                                                    Cell<Matrix<DDRMat>> const & aGlobalNodeCoord )
        {
            // Allocate space
            moris::size_t tNumNewNodes = aNewNodeIndices.size();
            moris::size_t tNumCurrNodes = mNodePhaseVals.n_rows();
        
            // add space to the node phase value table
            mNodePhaseVals.resize(tNumNewNodes+tNumCurrNodes, this->get_num_geometries());
        
            Cell<GEN_Geometry_Object> tGeometryObjects(tNumNewNodes);
        
            moris::Matrix< moris::IndexMat > tNodeIndex(1,tNumNewNodes);
            for(moris::size_t i = 0; i<tNumNewNodes; i++)
            {
                tGeometryObjects(i).set_phase_val_row(i+tNumCurrNodes);
                tNodeIndex(0,i) = aNewNodeIndices(i);
                if(aStoreParentTopo)
                {
                    tGeometryObjects(i).set_parent_entity_topology(aParentTopo(i)->copy());
                    mIntegNodeIndices.push_back( aNewNodeIndices(i) );
                }
            }
        
            if(tNumNewNodes !=0)
            {
                mGeometryObjectManager.store_geometry_objects(tNodeIndex,tGeometryObjects);
            }
        
            for(moris::size_t j = 0; j < this->get_num_geometries(); j++)
            {
        
                for(moris::size_t i = 0; i<tNumNewNodes; i++ )
                {
                    // Ask the pending node about its parent
                    // This information is needed to know what to interpolate based on
                    moris::Matrix< moris::DDRMat >  tLevelSetValues(1,1);
                    this->interpolate_level_set_value_to_child_node_location(*aParentTopo(i), j, aParamCoordRelativeToParent(i),tLevelSetValues);
        
                    mNodePhaseVals(i+tNumCurrNodes,j) = tLevelSetValues(0,0);
                }
        
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
        
        void Geometry_Engine::is_intersected( moris::Matrix< moris::DDRMat > const &   aNodeCoords,
                                                  moris::Matrix< moris::IndexMat > const & aNodetoEntityConn,
                                                  moris::size_t                            aCheckType,
                                                  Cell< GEN_Geometry_Object > &            aGeometryObjects )
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
        
        void Geometry_Engine::compute_interface_sensitivity( Matrix< IndexMat > const & aInterfaceNodeIndices,
                                                                 Matrix< DDRMat >   const & aNodeCoords,
                                                                 moris_index                aGeomIndex,
                                                                 bool               const   aGlbCoord )
        {
            // Figure out how many entities to compute sensitivity for
            uint tNumEntities = aInterfaceNodeIndices.numel();
        
            // iterate through node indices and compute sensitivity for each
            for( uint iEnt = 0; iEnt<tNumEntities; iEnt++ )
            {
                // get the node index
                moris::moris_index tNodeIndex = aInterfaceNodeIndices( iEnt );
        
                // Get the node geometry object
                GEN_Geometry_Object & tGeoObj = this->get_geometry_object( tNodeIndex );
        
                // Get the parent topology that this node was created on
                xtk::Topology const & tParentEdge = tGeoObj.get_parent_entity_topology( );
        
                MORIS_ASSERT(tParentEdge.get_topology_type() == xtk::Topology_Type::EDGE,"Only supporting interface sensitivity computation on an edge");
        
                // Get the node indices from the topology
                Matrix< IndexMat > const & tParentEntityNodes = tParentEdge.get_node_indices( );
        
                // Initialize sensitivity
                Matrix< DDRMat > tDxDp(1,1,0.0);
        
                // Get the node vars of the parent edge nodes
                Matrix< DDRMat > tEntityNodeVars( tParentEntityNodes.numel(), 1 );
                for(uint i = 0; i < tParentEntityNodes.numel(); i++)
                {
                    tEntityNodeVars(i) = this->get_entity_phase_val( tParentEntityNodes(i), aGeomIndex );
                }
        
                // Recompute local intersection (This could be stored instead)
                Matrix< DDRMat > tIntersectLocalCoordinate( 1, 1, 0.0 );
                Matrix< DDRMat > tIntersectGlobalCoordinate( 1, 1, 0.0 );
                get_intersection_location(mThresholdValue,
                                          mPerturbationValue,
                                          aNodeCoords,
                                          tEntityNodeVars,
                                          tParentEntityNodes,
                                          tIntersectLocalCoordinate,
                                          tIntersectGlobalCoordinate,
                                          true,
                                          aGlbCoord);
        
                // FIXME: Parent edge nodes need to not be the ADVs
                Matrix< IndexMat > tADVIndices;
        
                compute_dx_dp_for_an_intersection( tParentEntityNodes, aNodeCoords, tIntersectLocalCoordinate, tEntityNodeVars, tDxDp, tADVIndices );
        
                tGeoObj.set_sensitivity_dx_dp( tDxDp );
                tGeoObj.set_node_adv_indices( tParentEntityNodes );
              }
        }

        //--------------------------------------------------------------------------------------------------------------
        
        void Geometry_Engine::get_intersection_location( moris::real const &                      aIsocontourThreshold,
                                                             moris::real const &                      aPerturbationThreshold,
                                                             moris::Matrix< moris::DDRMat > const &   aGlobalNodeCoordinates,
                                                             moris::Matrix< moris::DDRMat > const &   aEntityNodeVars,
                                                             moris::Matrix< moris::IndexMat > const & aEntityNodeIndices,
                                                             moris::Matrix< moris::DDRMat > &         aIntersectionLocalCoordinates,
                                                             moris::Matrix< moris::DDRMat > &         aIntersectionGlobalCoordinates,
                                                             bool                                     aCheckLocalCoordinate,
                                                             bool                                     aComputeGlobalCoordinate )
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
        
        void Geometry_Engine::compute_dx_dp_finite_difference( moris::real                      const & aPerturbationVal,
                                                                   moris::Matrix< moris::DDRMat >   const & aGlobalNodeCoordinates,
                                                                   moris::Matrix< moris::DDRMat >   const & aEntityNodeCoordinates,
                                                                   moris::Matrix< moris::DDRMat >   const & aIntersectionLclCoordinate,
                                                                   moris::Matrix< moris::IndexMat > const & aEntityNodeIndices,
                                                                   moris::Matrix< moris::DDRMat >         & aEntityNodeVars,
                                                                   moris::Matrix< moris::DDRMat >         & aDxDp )
        {
        
            moris::size_t tNumNodeVars = aEntityNodeVars.n_rows();
            MORIS_ASSERT(tNumNodeVars == 2, "Currently compute_dx_dp_finite_difference has only been tested on edges");
            aDxDp.resize(2, 3);
        
            moris::real tPerturbationLen = 2 * aPerturbationVal;
            moris::real tScale   = 1/tPerturbationLen;
            Cell<moris::real>  tPerturbationSign = {1, -1};
        
            moris::Matrix< moris::DDRMat >       tDxDp(1, 3);
            moris::Matrix< moris::DDRMat >       tPerturbedLocalCoordinate(1, 1);
            Cell<moris::Matrix< moris::DDRMat >> tPerturbedGlobCoordinates = {moris::Matrix< moris::DDRMat >(1, 3),
                                                                     moris::Matrix< moris::DDRMat >(1, 3)};
            // Loop over all the nodes and perturb up and down
            for(moris::size_t i = 0; i < tNumNodeVars; i++)
            {
                // Perturb up and down
                for(moris::size_t j = 0; j < 2; j++)
                {
        
                    moris::real tPerturb = tPerturbationSign(j) * aPerturbationVal;
                    // Perturb
                    aEntityNodeVars(i, 0) = aEntityNodeVars(i, 0) + tPerturb;
        
                    // Locate perturbed interface
                    get_intersection_location(mThresholdValue, aPerturbationVal, aGlobalNodeCoordinates, aEntityNodeVars, aEntityNodeIndices, tPerturbedLocalCoordinate, tPerturbedGlobCoordinates(j),false, true);
        
                    // Reverse perturb
                    aEntityNodeVars(i, 0) = aEntityNodeVars(i, 0) - tPerturb;
        
                }
        
                tDxDp.matrix_data() = tScale * (tPerturbedGlobCoordinates(1).matrix_data() - tPerturbedGlobCoordinates(0).matrix_data());
        
                replace_row(0, tDxDp, i, aDxDp);
        
        
            }
        }

        //--------------------------------------------------------------------------------------------------------------
        
        void Geometry_Engine::compute_dx_dp_for_an_intersection( moris::Matrix< moris::IndexMat > const & aEntityNodeIndices,
                                                                     moris::Matrix< moris::DDRMat >   const & aGlobalNodeCoordinates,
                                                                     moris::Matrix< moris::DDRMat >   const & aIntersectionLclCoordinate,
                                                                     moris::Matrix< moris::DDRMat >         & aEntityNodeVars,
                                                                     moris::Matrix< moris::DDRMat >         & aDxDp,
                                                                     moris::Matrix< moris::IndexMat >       & aADVIndices )
        {
            moris::size_t tNumNodes = aEntityNodeIndices.n_cols();
        
            MORIS_ASSERT(tNumNodes == 2, "Currently, compute_dx_dp_for_an_intersection is only supported on edges");
        
            // Initialize
            moris::Cell< moris::Matrix< moris::DDRMat > > tDPhiiDp(2);// = { moris::Matrix< moris::DDRMat >(0,0), moris::Matrix< moris::DDRMat >(0,0) };
            uint tDim = aGlobalNodeCoordinates.n_cols();
            moris::Matrix< moris::DDRMat > tEntityNodeCoordinates( tNumNodes,tDim );
        
            // Assemble the entity local coordinates
            replace_row(aEntityNodeIndices(0,0), aGlobalNodeCoordinates,0,tEntityNodeCoordinates);
            replace_row(aEntityNodeIndices(0,1), aGlobalNodeCoordinates,1,tEntityNodeCoordinates);
        
            // Get information from geometry
            if (mGeometry(mActiveGeometryIndex)->sensitivities_available())
            {
                // Analytic
                aADVIndices = get_node_adv_indices_analytic();
                for(moris::size_t i = 0; i < tNumNodes; i++)
                {
                    mGeometry(mActiveGeometryIndex)->evaluate_sensitivity(aEntityNodeIndices(0, i),
                                                                          aGlobalNodeCoordinates.get_row(aEntityNodeIndices(0, i)),
                                                                          tDPhiiDp(i));
                }
                compute_dx_dp_with_linear_basis( tDPhiiDp(0), tDPhiiDp(1), tEntityNodeCoordinates, aEntityNodeVars, aDxDp );
            }
            else
            {
                // Discrete
                aADVIndices = aEntityNodeIndices;
                compute_dx_dp_finite_difference( mPerturbationValue, aGlobalNodeCoordinates, tEntityNodeCoordinates,
                        aIntersectionLclCoordinate, aEntityNodeIndices, aEntityNodeVars, aDxDp );
            }
        }

        //--------------------------------------------------------------------------------------------------------------
        
        GEN_Geometry_Object &
        Geometry_Engine::get_geometry_object(moris::size_t const & aNodeIndex)
        {
           return mGeometryObjectManager.get_geometry_object_from_manager(aNodeIndex);
        }

        //--------------------------------------------------------------------------------------------------------------
        
        GEN_Geometry_Object const &
        Geometry_Engine::get_geometry_object(moris::size_t const & aNodeIndex) const
        {
           return mGeometryObjectManager.get_geometry_object_from_manager(aNodeIndex);
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
        
        moris::real
        Geometry_Engine::get_entity_phase_val( moris::size_t const & aNodeIndex,
                                                   moris::size_t const & aGeomIndex )
        {
            GEN_Geometry_Object & tNodesGeoObj = get_geometry_object(aNodeIndex);
            moris::size_t tNodeRowIndex = tNodesGeoObj.get_phase_val_row();
        
            MORIS_ASSERT(tNodeRowIndex<mNodePhaseVals.n_rows(),"Entity row index out of bounds in the nodal phase val matrix");
        
            return mNodePhaseVals(tNodeRowIndex,aGeomIndex);
        }

        //--------------------------------------------------------------------------------------------------------------
        
        moris::Matrix< moris::DDRMat > const &
        Geometry_Engine::get_node_dx_dp( moris::size_t const & aNodeIndex ) const
        {
            GEN_Geometry_Object const & tNodesGeoObj = get_geometry_object( aNodeIndex );
            return tNodesGeoObj.get_sensitivity_dx_dp();
        }

        //--------------------------------------------------------------------------------------------------------------
        
        moris::Matrix< moris::IndexMat > const &
        Geometry_Engine::get_node_adv_indices( moris::size_t const & aNodeIndex ) const
        {
            GEN_Geometry_Object const & tNodesGeoObj = get_geometry_object( aNodeIndex );
            return tNodesGeoObj.get_node_adv_indices();
        }

        //--------------------------------------------------------------------------------------------------------------
        
        void Geometry_Engine::get_phase_index( moris::Matrix< moris::DDSTMat > const & aNodeIndex,
                                                   moris::Matrix< moris::DDSTMat > & aNodePhaseIndex )
        {
            // 0 for neg 1 for pos
            moris::real tNodePhaseValue = 0;
            moris::Matrix< moris::IndexMat > tPhaseOnOff(1, this->get_num_geometries());
        
            for(moris::size_t i = 0; i<aNodeIndex.n_cols(); i++)
            {
                GEN_Geometry_Object & tNodesGeoObj = get_geometry_object(aNodeIndex(0,i));
                moris::size_t tNodeRowIndex = tNodesGeoObj.get_phase_val_row();
        
                for(moris::size_t iG = 0; iG < this->get_num_geometries(); iG++)
                {
                    tNodePhaseValue =  mNodePhaseVals(tNodeRowIndex, iG);
        
                    // Negative
                    if(tNodePhaseValue<mThresholdValue)
                    {
                        tPhaseOnOff(0, iG) = 0;
                    }
        
                    else
                    {
                        tPhaseOnOff(0, iG) = 1;
                    }
                }
                aNodePhaseIndex(i, 0) = mPhaseTable.get_phase_index(tPhaseOnOff);
            }
        }

        //--------------------------------------------------------------------------------------------------------------
        
        void Geometry_Engine::get_phase_index( moris::moris_index const & aNodeIndex,
                                                   moris::size_t & aNodePhaseIndex )
        {
            // 0 for neg 1 for pos
            moris::real tNodePhaseValue = 0;
            moris::Matrix< moris::IndexMat > tPhaseOnOff(1, this->get_num_geometries());
        
            GEN_Geometry_Object & tNodesGeoObj = get_geometry_object(aNodeIndex);
            moris::size_t tNodeRowIndex = tNodesGeoObj.get_phase_val_row();
        
            for(moris::size_t iG = 0; iG < this->get_num_geometries(); iG++)
            {
                tNodePhaseValue =  mNodePhaseVals(tNodeRowIndex,iG);
        
                // Negative
                if(tNodePhaseValue<mThresholdValue)
                {
                    tPhaseOnOff(0,iG) = 0;
                }
        
                else
                {
                    tPhaseOnOff(0,iG) = 1;
                }
            }
            aNodePhaseIndex = mPhaseTable.get_phase_index(tPhaseOnOff);
        }

        //--------------------------------------------------------------------------------------------------------------

        moris::moris_index
        Geometry_Engine::get_elem_phase_index(moris::Matrix< moris::IndexMat > const & aElemOnOff)
        {
            return mPhaseTable.get_phase_index(aElemOnOff);
        }

        //--------------------------------------------------------------------------------------------------------------
        
        moris::size_t
        Geometry_Engine::get_node_phase_index_wrt_a_geometry( moris::size_t aNodeIndex,
                                                                  moris::size_t aGeometryIndex )
        {
            GEN_Geometry_Object & tNodesGeoObj = get_geometry_object(aNodeIndex);
            moris::size_t tNodeRowIndex = tNodesGeoObj.get_phase_val_row();
        
            moris::real tNodePhaseVal = mNodePhaseVals(tNodeRowIndex,aGeometryIndex);
        
            moris::size_t tPhaseOnOff = 10000;
        
            if(tNodePhaseVal < mThresholdValue)
            {
                tPhaseOnOff = 0;
            }
            else
            {
                tPhaseOnOff = 1;
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
        
        moris::Matrix< moris::IndexMat >
        Geometry_Engine::get_node_adv_indices_analytic()
        {
            moris::size_t tNumDVS = get_num_design_variables();
            moris::Matrix< moris::IndexMat > tMatrix(1,tNumDVS);
            for(moris::size_t i = 0; i<tNumDVS; i++)
            {
                tMatrix(0,i) = (moris::moris_index)i;
            }
            return tMatrix;
        }

        //--------------------------------------------------------------------------------------------------------------

        moris::size_t Geometry_Engine::get_num_design_variables()
        {
            return mADVs.length();
        }

        //--------------------------------------------------------------------------------------------------------------
        
        Geometry_Object_Manager* Geometry_Engine::get_all_geom_obj()
        {
            return &mGeometryObjectManager;
        }
        
        //--------------------------------------------------------------------------------------------------------------
        
        void Geometry_Engine::register_mesh( mtk::Mesh_Manager* aMesh )
        {
            mMeshManager = aMesh;
        
            mSpatialDim = mMeshManager->get_interpolation_mesh(0 )->get_spatial_dim();	// assuming there is only one pair in the manager
        }

        //--------------------------------------------------------------------------------------------------------------
        
        moris_index Geometry_Engine::register_mesh( std::shared_ptr< moris::hmr::Mesh > aMesh )   //FIXME: this needs to be deleted and the GE should only be able to register a mesh pair
        {
            mMesh_HMR.push_back( aMesh );
        
            mSpatialDim = mMesh_HMR( mMesh_HMR.size()-1 )->get_spatial_dim();
        
            return mMesh_HMR.size()-1;
        }

        //--------------------------------------------------------------------------------------------------------------

        moris_index Geometry_Engine::register_mesh( std::shared_ptr< hmr::HMR > aHMR )   //FIXME: same as above
        {
            mHMRPerformer.push_back(aHMR);
            mMesh_HMR.push_back( aHMR->create_mesh(0) );

            mSpatialDim = mMesh_HMR( mMesh_HMR.size()-1 )->get_spatial_dim();

            return mMesh_HMR.size()-1;
        }

        //--------------------------------------------------------------------------------------------------------------

        void Geometry_Engine::perform( )
        {
            this->perform_refinement();
        }

        //--------------------------------------------------------------------------------------------------------------
        
        void Geometry_Engine::perform_refinement()
        {
            for( uint Ik = 0; Ik < mNumRefinements; ++Ik )
            {
                moris::Cell< moris::Matrix< DDRMat > > tValues;
        
                this->get_field_values_for_all_geometries( tValues );
        
                for( uint Ij = 0; Ij < tValues.size(); ++Ij )
                {
                    mHMRPerformer( 0 )->based_on_field_put_elements_on_queue( tValues( Ij ), 0 );
                }

                mHMRPerformer( 0 )->perform_refinement_based_on_working_pattern( 0, false );
            }
        }

        //--------------------------------------------------------------------------------------------------------------

        //Matrix< DDRMat > Geometry_Engine::get_cylinder_vals( moris_index aWhichMesh,
        //                                                         GEN_CylinderWithEndCaps* aFiber,
        //                                                         uint aNumberOfFibers )  //FIXME this is currently only setup to work with an HMR member mesh
        //{
        //    uint tNumOfIPNodes = mMesh_HMR( aWhichMesh )->get_num_nodes();
        //
        //    Matrix< DDRMat > tLSVals(tNumOfIPNodes, 1, 1.0); // FIXME: 10.0 needs to be replaced with problem dependent value
        //
        //    for( uint k=0; k<aNumberOfFibers; ++k )
        //    {
        //        uint tNumCylinders = aFiber->get_number_of_cylinders(k);
        //
        //        for( uint l=0; l<tNumCylinders; ++l )
        //        {
        //            Matrix<DDRMat> tMidPoint;
        //            Matrix<DDRMat> tLength;
        //            aFiber->midPoint_and_BB_dims( k, l, tMidPoint, tLength, 1.0 ); // FIXME: 1.0 needs to be the size of the coarsest element ( this is the bounding box buffer value )
        //
        //            Matrix<IndexMat> tNodeIndices;
        //            mMesh_HMR( aWhichMesh )->get_nodes_indices_in_bounding_box( tMidPoint,
        //                                                                        { { tLength(0) },{ tLength(1) } ,{ tLength(2) }},
        //                                                                        tNodeIndices );
        //            for( uint i=0; i<tNodeIndices.numel(); ++i )
        //            {
        //                Matrix<DDRMat> tVertexCoords = mMesh_HMR( aWhichMesh )->get_mtk_vertex( tNodeIndices( i ) ).get_coords();
        //                tLSVals( tNodeIndices( i ) ) = std::min( tLSVals( tNodeIndices( i ) ),
        //                                                         aFiber->create_cylinder( tVertexCoords, k, l ) );
        //            }
        //
        //        }
        //    }
        //
        //      return tLSVals;
        //}


        size_t Geometry_Engine::analytic_geometry_index(size_t aGlobalGeometryIndex)
        {
            return aGlobalGeometryIndex;
        }

        //--------------------------------------------------------------------------------------------------------------

        size_t Geometry_Engine::discrete_geometry_index(size_t aGlobalGeometryIndex)
        {
            return aGlobalGeometryIndex - mGeometry.size();
        }

        //--------------------------------------------------------------------------------------------------------------

        bool
        Geometry_Engine::compute_intersection_info( moris::moris_index               const & aEntityIndex,
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
            moris::Matrix< moris::DDRMat > tInterpLocationCoords(1,1);
        
            // Loop through nodes and get levelset values from precomputed values in aNodeVars or in the levelset mesh
            for(moris::size_t n = 0; n < tNumNodes; n++)
            {   //Get node id n
                tNodeInd = aEntityNodeInds(n);
        
                GEN_Geometry_Object & tGeoObj = get_geometry_object(tNodeInd);
                moris::size_t tPhaseValRowIndex = tGeoObj.get_phase_val_row();
                tEntityNodeVars(n) = mNodePhaseVals(tPhaseValRowIndex, mActiveGeometryIndex);
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
        
        //        MORIS_ASSERT(tMax != mThresholdValue && tMin != mThresholdValue, "Threshold levelset value at all nodes! There is no handling of this inside XTK currently.");
        
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
                    if(mComputeDxDp)
                    {
                        moris::Matrix< moris::DDRMat > tDxDp(1,1,100.0);
                        compute_dx_dp_for_an_intersection(aEntityNodeInds,aNodeCoords,tIntersectLocalCoordinate,tEntityNodeVars, tDxDp, aNodeADVIndices);
                        aGeometryObject.set_sensitivity_dx_dp(tDxDp);
                        aGeometryObject.set_node_adv_indices(aNodeADVIndices);
                    }
               }
            }
        
            return tIsIntersected;
        
        }

        //--------------------------------------------------------------------------------------------------------------

        void
        Geometry_Engine::interpolate_level_set_value_to_child_node_location( xtk::Topology const & aParentTopology,
                                                                                 moris::size_t const &                  aGeometryIndex,
                                                                                 moris::Matrix< moris::DDRMat > const & aNodeLocalCoordinate,
                                                                                 moris::Matrix< moris::DDRMat >       & aLevelSetValues  )
        {
                     // Get node indices attached to parent (These are indices relative to another mesh and may need to be mapped)
            moris::Matrix< moris::IndexMat > const & tNodesAttachedToParent = aParentTopology.get_node_indices();
        
            // Get number of nodes attached to parent
            moris::size_t tNumNodesAttachedToParent = tNodesAttachedToParent.numel();
            moris::Matrix< moris::DDRMat > tNodesLevelSetValues(1, tNumNodesAttachedToParent);

            for(moris::size_t i = 0; i < tNumNodesAttachedToParent; i++)
            {
                GEN_Geometry_Object & tGeoObj = get_geometry_object(tNodesAttachedToParent(i));
                moris::size_t tPhaseRow = tGeoObj.get_phase_val_row();

                tNodesLevelSetValues(0,i) = mNodePhaseVals(tPhaseRow,aGeometryIndex);
            }
        
            // Ask the topology how to interpolate
            moris::Matrix< moris::DDRMat > tBasisValues(1,1);
            xtk::Basis_Function const & tParentBasisFunctions = aParentTopology.get_basis_function();
        
            // Evaluate basis function
            tParentBasisFunctions.evaluate_basis_function(aNodeLocalCoordinate,tBasisValues);
        
            // Compute \phi = Ni.\phi_i
            aLevelSetValues = tBasisValues*moris::trans(tNodesLevelSetValues);
        
         }

        //--------------------------------------------------------------------------------------------------------------
        
        void Geometry_Engine::get_field_values_for_all_geometries( moris::Cell< Matrix< DDRMat > > & aAllFieldVals,
                                                  const moris_index                 aWhichMesh )
        {
            //TODO: implement for a mesh manager (rather than just an HMR mesh)
            uint tNumVertices = mMesh_HMR( aWhichMesh )->get_num_nodes();
        
            aAllFieldVals.resize(this->get_num_geometries());

            // Evaluate field values
            for ( uint Ik = 0; Ik < mGeometry.size(); Ik++ )
            {
                aAllFieldVals( Ik ).set_size( tNumVertices, 1, - MORIS_REAL_MAX );
                for( uint iVert = 0; iVert <tNumVertices; iVert++)
                {
                    Matrix< DDRMat > tCoord = mMesh_HMR( aWhichMesh )->get_mtk_vertex( iVert ).get_coords();
                    aAllFieldVals( Ik )( iVert ) = mGeometry(Ik)->evaluate_field_value(iVert, tCoord);
        
                    // FIXME will not work in parallel. Ind are not consistent because of aura
                }
            }
        }

        //--------------------------------------------------------------------------------------------------------------

        void Geometry_Engine::assign_ip_hosts_by_set_name(std::string                     aSetName,
                                                              std::shared_ptr< Property > aPropertyPointer,
                                                              PDV_Type                     aPdvType,
                                                              moris_index                     aWhichMesh)
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

        void Geometry_Engine::assign_ip_hosts_by_set_index( moris_index                     aSetIndex,
                                                                std::shared_ptr< Property > aPropertyPointer,
                                                                PDV_Type                     aPdvType,
                                                                moris_index                     aWhichMesh)
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

                    // Indices on cluster FIXME get rid of copying over indices
                    tNodeIndicesInCluster = tCluster->get_interpolation_cell().get_vertex_inds();
                    tNodeIndicesPerSet(tMeshSetIndex).resize(tNodeIndicesPerSet(tMeshSetIndex).length() + tNodeIndicesInCluster.length(), 1);
                    
                    for (uint tNodeInCluster = 0; tNodeInCluster < tNodeIndicesInCluster.length(); tNodeInCluster++)
                    {
                        tNodeIndicesPerSet(tMeshSetIndex)(tCurrentNode++) = tNodeIndicesInCluster(tNodeInCluster);
                    }
                }
                //tNodeIndicesPerSet(tMeshSetIndex) = tInterpolationMesh->get_set_by_index(tMeshSetIndex)->get_vertices_inds_on_block(false);
            }

            // Create hosts
            mPdvHostManager.create_ip_pdv_hosts(tNumNodes, tNodeIndicesPerSet, aPdvTypes);
        }

        //--------------------------------------------------------------------------------------------------------------

        void Geometry_Engine::create_ig_pdv_hosts(moris_index aMeshIndex)
        {
            // Get information from integration mesh
            mtk::Integration_Mesh* tIntegrationMesh = mMeshManager->get_integration_mesh(aMeshIndex);
            uint tNumSets = tIntegrationMesh->get_num_sets();
            uint tNumNodes = tIntegrationMesh->get_num_nodes();
            Cell<Matrix<DDSMat>> tNodeIndicesPerSet(tNumSets);

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
                    MORIS_ERROR( false,"Geometry_Engine::initialize_integ_pdv_host_list() - Geometry Engine only works for 2D and 3D models." );
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

            // Create hosts
            mPdvHostManager.create_ig_pdv_hosts(tNumNodes, tNodeIndicesPerSet, tPdvTypes);

            // Assign PDVs
            Matrix<F31RMat> tCoordinates;
            for(uint tNodeIndex = 0; tNodeIndex < tNumNodes; tNodeIndex++)
            {
                // Get coordinates
                tCoordinates = tIntegrationMesh->get_node_coordinate(tNodeIndex);

                // Create PDVs
                mPdvHostManager.create_ig_pdv(tNodeIndex, PDV_Type::X_COORDINATE, tCoordinates(0));
                mPdvHostManager.create_ig_pdv(tNodeIndex, PDV_Type::Y_COORDINATE, tCoordinates(1));
                if (mSpatialDim == 3)
                {
                    mPdvHostManager.create_ig_pdv(tNodeIndex, PDV_Type::Z_COORDINATE, tCoordinates(2));
                }
            }
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


