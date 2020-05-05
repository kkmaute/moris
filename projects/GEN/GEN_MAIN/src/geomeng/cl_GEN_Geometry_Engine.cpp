// GEN
#include "cl_GEN_Geometry_Engine.hpp"
#include "fn_GEN_create_geometry.hpp"

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

        // ADVs
        mADVs(string_to_mat<DDRMat>(aParameterLists(0)(0).get<std::string>("initial_advs"))),
        mLowerBounds(string_to_mat<DDRMat>(aParameterLists(0)(0).get<std::string>("lower_bounds"))),
        mUpperBounds(string_to_mat<DDRMat>(aParameterLists(0)(0).get<std::string>("upper_bounds"))),

        // Phase table
        mPhaseTable(string_to_mat<IndexMat>(aParameterLists(0)(0).get<std::string>("phase_table")).numel()
              ? Phase_Table(string_to_mat<IndexMat>(aParameterLists(0)(0).get<std::string>("phase_table")), aParameterLists(0)(0).get<std::string>("phase_table_structure"))
              : Phase_Table(aParameterLists(1).size(), aParameterLists(0)(0).get<std::string>("phase_table_structure")))

        {
            // Build geometry (just analytic for right now)
            if (aParameterLists(1).size() > 0)
            {
                mGeometryAnalytic.resize(aParameterLists(1).size());
                for (uint tGeometryIndex = 0; tGeometryIndex < aParameterLists(1).size(); tGeometryIndex++)
                {
                    mGeometryAnalytic(tGeometryIndex) = create_geometry(aParameterLists(1)(tGeometryIndex), mADVs, aLibrary);
                }
            }
        }

        //--------------------------------------------------------------------------------------------------------------

        Geometry_Engine::Geometry_Engine(Cell<std::shared_ptr<Geometry_Analytic>>   aGeometry,
                                         Phase_Table                                aPhaseTable,
                                         uint                                       aSpatialDim,
                                         real                                       aThresholdValue,
                                         real                                       aPerturbationValue)
            : mSpatialDim(aSpatialDim),
              mThresholdValue(aThresholdValue),
              mPerturbationValue(aPerturbationValue),
              mActiveGeometryIndex(0),
              mGeometryAnalytic(aGeometry),
              mPhaseTable(aPhaseTable)
        {
        }

        //--------------------------------------------------------------------------------------------------------------

        Geometry_Engine::Geometry_Engine(Cell<std::shared_ptr<Geometry_Discrete>>   aGeometry,
                                                 Phase_Table                                aPhaseTable,
                                                 uint                                       aSpatialDim,
                                                 real                                       aThresholdValue,
                                                 real                                       aPerturbationValue)
            : mSpatialDim(aSpatialDim),
              mThresholdValue(aThresholdValue),
              mPerturbationValue(aPerturbationValue),
              mActiveGeometryIndex(0),
              mGeometryDiscrete(aGeometry),
              mPhaseTable(aPhaseTable)
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
            uint tNumAnalyticGeometries = mGeometryAnalytic.size();
        
            // Loop through each geometry and then each node and compute the level set field value
            // add value to phase value matrix
            for (moris::size_t tGeometryIndex = 0; tGeometryIndex < tNumAnalyticGeometries; tGeometryIndex++) // Analytic
            {
                for (moris::size_t tNodeIndex = 0; tNodeIndex < tNumNodes; tNodeIndex++)
                {
                    mNodePhaseVals(tNodeIndex, tGeometryIndex) =
                            mGeometryAnalytic(this->analytic_geometry_index(tGeometryIndex))->evaluate_field_value(aNodeCoords.get_row(tNodeIndex));
                }
            }

            for (moris::size_t tGeometryIndex = tNumAnalyticGeometries; tGeometryIndex < this->get_num_geometries(); tGeometryIndex++) // Discrete
            {
                for (moris::size_t tNodeIndex = 0; tNodeIndex < tNumNodes; tNodeIndex++ )
                {
                    mNodePhaseVals(tNodeIndex, tGeometryIndex) =
                            mGeometryDiscrete(this->discrete_geometry_index(tGeometryIndex))->evaluate_field_value(tNodeIndex);
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
        
            // Get information from a analytic geometry
            if (mActiveGeometryIndex < mGeometryAnalytic.size())
            {
                aADVIndices = get_node_adv_indices_analytic();
                for(moris::size_t i = 0; i < tNumNodes; i++)
                {
                    tDPhiiDp(i) = mGeometryAnalytic(this->analytic_geometry_index(mActiveGeometryIndex))
                            ->evaluate_sensitivity(aGlobalNodeCoordinates.get_row(aEntityNodeIndices(0, i)));
                }

                compute_dx_dp_with_linear_basis( tDPhiiDp(0), tDPhiiDp(1), tEntityNodeCoordinates, aEntityNodeVars, aDxDp );
            }

            // Get information from a discrete geometry
            else
            {
                aADVIndices = aEntityNodeIndices; // FIXME I don't undestand what noah originally intended
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
            return mGeometryAnalytic.size() + mGeometryDiscrete.size();
        }
        
        //--------------------------------------------------------------------------------------------------------------
        
        bool Geometry_Engine::is_geometry_analytic()
        {
            return (mActiveGeometryIndex < mGeometryAnalytic.size());
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
            MORIS_ASSERT(mActiveGeometryIndex < this->get_num_geometries(), "Trying to advance past the number of geometries in the geometry engine");
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
        
        Pdv_Host_Manager* Geometry_Engine::get_pdv_host_manager(  )
        {
            return &mPdvHostManager;
        }

        //--------------------------------------------------------------------------------------------------------------
        
        Geometry_Object_Manager* Geometry_Engine::get_all_geom_obj()
        {
            return &mGeometryObjectManager;
        }
        
        //--------------------------------------------------------------------------------------------------------------
        
        void Geometry_Engine::register_mesh( mtk::Mesh_Manager* aMesh )
        {
            mMesh = aMesh;
        
            mSpatialDim = mMesh->get_interpolation_mesh( 0 )->get_spatial_dim();	// assuming there is only one pair in the manager
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
        
        //--------------------------------------------------------------------------------------------------------------
        
        
        
        
        //--------------------------------------------------------------------------------------------------------------
        // private functions
        //--------------------------------------------------------------------------------------------------------------

        size_t Geometry_Engine::analytic_geometry_index(size_t aGlobalGeometryIndex)
        {
            return aGlobalGeometryIndex;
        }

        //--------------------------------------------------------------------------------------------------------------

        size_t Geometry_Engine::discrete_geometry_index(size_t aGlobalGeometryIndex)
        {
            return aGlobalGeometryIndex - mGeometryAnalytic.size();
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

            // Analytic
            for ( uint Ik = 0; Ik < mGeometryAnalytic.size(); Ik++ )
            {
                aAllFieldVals( Ik ).set_size( tNumVertices, 1, - MORIS_REAL_MAX );
                for( uint iVert = 0; iVert <tNumVertices; iVert++)
                {
                    Matrix< DDRMat > tCoord = mMesh_HMR( aWhichMesh )->get_mtk_vertex( iVert ).get_coords();
                    aAllFieldVals( Ik )( iVert ) = mGeometryAnalytic(this->analytic_geometry_index(Ik))->evaluate_field_value(tCoord);
        
                    // FIXME will not work in parallel. Ind are not consistent because of aura
                }
            }

            // Discrete
            for ( uint Ik = mGeometryAnalytic.size(); Ik < this->get_num_geometries(); Ik++ )
            {
                aAllFieldVals( Ik ).set_size( tNumVertices, 1, - MORIS_REAL_MAX );
                for( uint iVert = 0; iVert <tNumVertices; iVert++)
                {
                    aAllFieldVals( Ik )( iVert ) = mGeometryDiscrete(this->discrete_geometry_index(Ik))->evaluate_field_value(iVert);
                }
            }
        }

        //--------------------------------------------------------------------------------------------------------------

        void Geometry_Engine::set_pdv_types(Cell<Cell<Cell<GEN_DV>>> aPdvTypes)
        {
            // set the set dv type flag to true
            mTypesSet = true;

            // set the dv type list for the pdv host manager
            mPdvHostManager.set_ip_pdv_types(aPdvTypes);
        }

        //--------------------------------------------------------------------------------------------------------------

        void Geometry_Engine::initialize_interp_pdv_host_list(moris_index aWhichMesh)
        {
            // check if the dv type list was set
            MORIS_ASSERT( mTypesSet, "Geometry_Engine::initialize_interp_pdv_host_list() - set_pdv_types() must be called before this function." );

            // get number of vertices on the IP mesh
            uint tTotalNumVertices = mMesh->get_interpolation_mesh(aWhichMesh)->get_num_nodes();

            // ask pdv host manager to init host
            mPdvHostManager.initialize_ip_hosts( tTotalNumVertices );
        }

        //--------------------------------------------------------------------------------------------------------------

        void Geometry_Engine::assign_ip_hosts_by_set_name( std::string                  aSetName,
                                                               std::shared_ptr< GEN_Field > aFieldPointer,
                                                               enum GEN_DV                  aPdvType,
                                                               moris_index                  aWhichMesh)
        {
            // get the mesh set from name
            moris::mtk::Set* tSetPointer = mMesh->get_integration_mesh( aWhichMesh )->get_set_by_name( aSetName );

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
                for(uint iVert=0; iVert<tNumVerts; iVert++)
                {
                    // get the vertex index
                    moris_index tVertIndex = tVertices(iVert)->get_index();

                    // ask pdv host manager to assign to vertex a pdv type and a property
                    mPdvHostManager.assign_field_to_pdv_type_by_vertex_index( aFieldPointer, aPdvType, tVertIndex );
                }
            }

            // ask pdv host manager to update local to global dv type map
            mPdvHostManager.update_ip_local_to_global_dv_type_map();

            // mark this DV type as unchanging TODO
            // mPdvHostManager.mark_ip_pdv_as_unchanging( aPdvType );
        }

        //--------------------------------------------------------------------------------------------------------------

        void Geometry_Engine::assign_ip_hosts_by_set_name(std::string                     aSetName,
                                                              std::shared_ptr< GEN_Property > aPropertyPointer,
                                                              enum GEN_DV                     aPdvType,
                                                              moris_index                     aWhichMesh)
        {
            // get the mesh set from name
            moris::mtk::Set* tSetPointer = mMesh->get_integration_mesh( aWhichMesh )->get_set_by_name( aSetName );

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
                for(uint iVert=0; iVert<tNumVerts; iVert++)
                {
                    // get the vertex index
                    moris_index tVertIndex = tVertices(iVert)->get_index();

                    // ask pdv host manager to assign to vertex a pdv type and a property
                    mPdvHostManager.assign_property_to_pdv_type_by_vertex_index( aPropertyPointer, aPdvType, tVertIndex );
                }
            }

            // ask pdv host manager to update local to global dv type map
            mPdvHostManager.update_ip_local_to_global_dv_type_map();
        }

        //--------------------------------------------------------------------------------------------------------------

        void Geometry_Engine::assign_ip_hosts_by_set_index(moris_index                  aSetIndex,
                                                               std::shared_ptr< GEN_Field > aFieldPointer,
                                                               enum GEN_DV                  aPdvType,
                                                               moris_index                  aWhichMesh)
        {
            // get the mesh set from index
            moris::mtk::Set* tSetPointer = mMesh->get_integration_mesh( aWhichMesh )->get_set_by_index( aSetIndex );

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
                for(uint iVert=0; iVert<tNumVerts; iVert++)
                {
                    // get the vertex index
                    moris_index tVertIndex = tVertices(iVert)->get_index();

                    // ask pdv host manager to assign to vertex a pdv type and a property
                    mPdvHostManager.assign_field_to_pdv_type_by_vertex_index( aFieldPointer, aPdvType, tVertIndex );
                }
            }

            // ask pdv host manager to update local to global dv type map
            mPdvHostManager.update_ip_local_to_global_dv_type_map();

            // mark this DV type as unchanging TODO
            // mPdvHostManager.mark_ip_pdv_as_unchanging( aPdvType );
        }

        //--------------------------------------------------------------------------------------------------------------

        void Geometry_Engine::assign_ip_hosts_by_set_index( moris_index                     aSetIndex,
                                                                std::shared_ptr< GEN_Property > aPropertyPointer,
                                                                enum GEN_DV                     aPdvType,
                                                                moris_index                     aWhichMesh)
        {
            // get the mesh set from index
            moris::mtk::Set* tSetPointer = mMesh->get_integration_mesh( aWhichMesh )->get_set_by_index( aSetIndex );

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
                for(uint iVert=0; iVert<tNumVerts; iVert++)
                {
                    // get the vertex index
                    moris_index tVertIndex = tVertices(iVert)->get_index();

                    // ask pdv host manager to assign to vertex a pdv type and a property
                    mPdvHostManager.assign_property_to_pdv_type_by_vertex_index( aPropertyPointer, aPdvType, tVertIndex );
                }
            }

            // ask pdv host manager to update local to global dv type map
            mPdvHostManager.update_ip_local_to_global_dv_type_map();
        }

        //--------------------------------------------------------------------------------------------------------------

        void Geometry_Engine::initialize_integ_pdv_host_list(moris_index aWhichMesh) //FIXME
        {
            Cell<Cell<Cell<GEN_DV>>> tDimDvList(1);
            tDimDvList(0).resize(1);
            tDimDvList(0)(0).resize(mSpatialDim);

            switch(mSpatialDim)
            {
                case(2):
                {
                    tDimDvList(0)(0)(0) = GEN_DV::XCOORD;
                    tDimDvList(0)(0)(1) = GEN_DV::YCOORD;
                    break;
                }
                case(3):
                {
                    tDimDvList(0)(0)(0) = GEN_DV::XCOORD;
                    tDimDvList(0)(0)(1) = GEN_DV::YCOORD;
                    tDimDvList(0)(0)(2) = GEN_DV::ZCOORD;
                    break;
                }
                default:
                {
                    MORIS_ERROR( false,"Geometry_Engine::initialize_integ_pdv_host_list() - Geometry Engine only works for 2D and 3D models." );
                }
            }

            mPdvHostManager.set_ig_pdv_types( tDimDvList );

            // get number of vertices on the IG mesh
            uint tTotalNumVertices = mMesh->get_integration_mesh(aWhichMesh)->get_num_nodes();

            // ask pdv host manager to initialize hosts
            mPdvHostManager.initialize_ig_hosts( tTotalNumVertices );

            // check to make sure there are IG nodes to put hosts on
            MORIS_ASSERT( mIntegNodeIndices.size() != 0, "Geometry_Engine::initialize_integ_pdv_host_list() - no integration node indices stored, has the XTK model performed decomposition already?" );

            uint tNumIndices = mIntegNodeIndices.size();

            mPdvHostManager.update_ig_local_to_global_dv_type_map();

            for( uint iInd=0; iInd<tNumIndices; iInd++ )
            {
                mPdvHostManager.create_ig_pdv_host( mSpatialDim, mIntegNodeIndices(iInd) );

                //mPdvHostManager.get_ig_pdv_host( mIntegNodeIndices(iInd) )->update_pdv_list( mSpatialDim ); FIXME

                Matrix< DDRMat > tTempCoords = mMesh->get_integration_mesh( aWhichMesh )->get_node_coordinate( mIntegNodeIndices(iInd) );

                for(uint iDim=0; iDim<mSpatialDim; iDim++)
                {
                    //mPdvHostManager.get_ig_pdv_host( mIntegNodeIndices(iInd) )->create_pdv( tTempCoords(iDim), tDimDvList(iDim), mPdvHostManager.get_ig_global_map() ); FIXME
                }
            }
        }

        //--------------------------------------------------------------------------------------------------------------

        void Geometry_Engine::mark_ig_pdv_as_inactive(moris_index aNodeIndex, GEN_DV aPdvType)
        {
            mPdvHostManager.mark_ig_pdv_as_inactive(aNodeIndex, aPdvType);
        }

        //--------------------------------------------------------------------------------------------------------------


    }   // end ge namespace
}   // end moris namespace


