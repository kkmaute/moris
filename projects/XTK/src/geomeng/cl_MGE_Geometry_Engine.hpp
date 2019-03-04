/*
 * cl_Geometry_Engine.hpp
 *
 *  Created on: Jun 21, 2017
 *      Author: ktdoble
 */

#ifndef SRC_GEOMENG_CL_MGE_GEOMETRY_ENGINE_HPP_
#define SRC_GEOMENG_CL_MGE_GEOMETRY_ENGINE_HPP_

// Standard library includes
#include <memory> // for shared_ptr
#include <math.h>

// Linear algebra includes
#include "cl_Matrix.hpp"
#include "fn_trans.hpp"
#include "op_times.hpp"
#include "linalg_typedefs.hpp"

// base class for geometry
#include "cl_Geometry.hpp"

// MTK includes
#include "cl_Mesh_Enums.hpp"

// Geometry Engine Includes
#include "cl_MGE_Enums.hpp"
#include "cl_MGE_Geometry_Object.hpp"
#include "cl_MGE_Geometry_Object_Manager.hpp"

// XTKL: General Includes

#include "cl_Logger.hpp"

// XTKL: Container includes
#include "cl_Cell.hpp"

// XTKL: Tools includes
#include "cl_Interpolaton.hpp"

//XTKL: Geometry Engine Includes
#include "cl_XTK_Phase_Table.hpp"
#include "cl_XTK_Pending_Node.hpp"

//XTKL: Topology
#include "cl_XTK_Topology.hpp"
#include "cl_XTK_Basis_Function.hpp"


using namespace moris;

namespace xtk
{


/*
 *
 * $\frac{\partial{\phi_A}}{\partial{p}}$ (change in phi with respect to a design variable
 * See for more detailed description of this function:
 *
 */
inline
void compute_dx_dp_with_linear_basis(moris::Matrix< moris::DDRMat >  & aDPhiADp,
                                     moris::Matrix< moris::DDRMat >  & aDPhiBDp,
                                     moris::Matrix< moris::DDRMat >  & aEdgeCoordinates,
                                     moris::Matrix< moris::DDRMat >  & aEdgeNodePhi,
                                     moris::Matrix< moris::DDRMat >  & aDxDp)
{

    MORIS_ASSERT(aDPhiADp.n_rows() != 0,"dPhi/dp not implemented in geometry would cause a seg fault here");
    MORIS_ASSERT(aDPhiBDp.n_rows() != 0,"dPhi/dp not implemented in geometry would cause a seg fault here");
    moris::real const & tPhiA = aEdgeNodePhi(0,0);
    moris::real const & tPhiB = aEdgeNodePhi(1,0);

    // Initialize
    moris::Matrix< moris::DDRMat > tXa = aEdgeCoordinates.get_row(0);

    moris::Matrix< moris::DDRMat > tXb = aEdgeCoordinates.get_row(1);

    // Compute $\frac{\partial x_{\Gamma}}{\partial \phi}$
    moris::DDRMat tDxgammaDphiA = -(tPhiB)/std::pow((tPhiA-tPhiB),2)*(tXb.matrix_data()-tXa.matrix_data());
    moris::DDRMat tDxgammaDphiB =  (tPhiA)/std::pow((tPhiA-tPhiB),2)*(tXb.matrix_data()-tXa.matrix_data());

    moris::Matrix< moris::DDRMat > tDxgDphiAMat(tDxgammaDphiA);
    moris::Matrix< moris::DDRMat > tDxgDphiBMat(tDxgammaDphiB);

    // Compute dx/dp
    moris::DDRMat tDxDp = aDPhiADp * moris::trans(tDxgDphiAMat) +  aDPhiBDp * moris::trans(tDxgDphiBMat);
    aDxDp = moris::Matrix< moris::DDRMat >(tDxDp);

}


class Geometry_Engine
{

public:

    // Single geometry constructor
    Geometry_Engine(Geometry  & aGeometry,
                    Phase_Table const & aPhaseTable) :
        mThresholdValue(0),
        mComputeDxDp(false),
        mActiveGeometryIndex(0),
        mGeometryObjects(),
        mPhaseTable(aPhaseTable)
    {
        mPerturbationValue = 0.0000005;
        mGeometry.push_back(&aGeometry);
    }

    // geometry vector constructor
    Geometry_Engine(Cell<Geometry*> const  & aGeometry,
                    Phase_Table const & aPhaseTable) :
        mThresholdValue(0),
        mComputeDxDp(false),
        mActiveGeometryIndex(0),
        mGeometry(aGeometry),
        mPhaseTable(aPhaseTable)
    {
        mPerturbationValue = 0.0000005;
    }


    Geometry_Engine()
    {
    }

    // Options which the user can change (all are given defaults)
    moris::real mThresholdValue;
    moris::real mPerturbationValue;
    bool        mComputeDxDp; // Should be turned off if a sensitivity has not been implemented


    /*
     * Initiali allocation of geometry objects, this creates a geometry object for each node coordinate.
     * In this case, aNodeCoords needs to be ordered by proc indices.
     */
    void
    initialize_geometry_objects_for_background_mesh_nodes(moris::size_t const & aNumNodes)
    {
        // Allocate space
        mNodePhaseVals = moris::Matrix< moris::DDRMat >(aNumNodes,mGeometry.size(),0);

        // Allocate geometry object
        Cell<Geometry_Object> tGeometryObjects(aNumNodes);

        // Associate each geometry object with a row in phase val matrix (note phase val computed later)
        moris::Matrix< moris::IndexMat > tNodeIndex(1,aNumNodes);
        for(moris::size_t i = 0; i<aNumNodes; i++)
        {
            tGeometryObjects(i).set_phase_val_row(i);
            tNodeIndex(0,i) = i;
        }

        // Place these in the geometry object manager
        mGeometryObjects.store_geometry_objects(tNodeIndex,tGeometryObjects);

    }

    void
    initialize_geometry_object_phase_values(moris::Matrix< moris::DDRMat > const & aNodeCoords)
    {
        // Allocate space
        moris::size_t tNumNodes = aNodeCoords.n_rows();

        // Loop through each geometry and then each node and compute the level set field value
        // add value to phase value matrix

        for(moris::size_t j = 0; j<get_num_geometries(); j++)
        {
            bool tAnalyticFlag = mGeometry(j)->is_analytic();

            // Analytic
            if(tAnalyticFlag)
            {
                for(moris::size_t i = 0; i<tNumNodes; i++ )
                {
                    mNodePhaseVals(i,j) = mGeometry(j)->evaluate_field_value_with_coordinate(i,aNodeCoords);
                }
            }

            // Discrete
            else
            {
                for(moris::size_t i = 0; i<tNumNodes; i++ )
                {
                    mNodePhaseVals(i,j) = mGeometry(j)->access_field_value_with_entity_index(i, moris::EntityRank::NODE);
                }
            }
        }
    }


    /*
     * Creates a geometry object association for pending nodes
     * These nodes have node indices and parent information
     */
    void
    associate_new_nodes_with_geometry_object(Cell<Pending_Node> & aNewNodes,
                                             bool aInterfaceNodes)
    {

        // Allocate space
        moris::size_t tNumNewNodes = aNewNodes.size();
        moris::size_t tNumCurrNodes = mNodePhaseVals.n_rows();

        // add space to the node phase value table
        mNodePhaseVals.resize(tNumNewNodes+tNumCurrNodes,mGeometry.size());

        Cell<Geometry_Object> tGeometryObjects(tNumNewNodes);

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
            mGeometryObjects.store_geometry_objects(tNodeIndex,tGeometryObjects);
        }

        // Compute and store level set value of this node for each new node
        for(moris::size_t j = 0; j<get_num_geometries(); j++)
        {

            for(moris::size_t i = 0; i<tNumNewNodes; i++ )
            {
                // Ask the pending node about its parent
                // This information is needed to know what to interpolate based on
                moris::Matrix< moris::DDRMat > const & tLocalCoordinate = aNewNodes(i).get_local_coordinate_relative_to_parent();
                moris::Matrix< moris::DDRMat >  tLevelSetValues(1,1);
                Topology const & tParentTopology = aNewNodes(i).get_parent_topology();

                // Interpolate all level set values to node
                this->interpolate_level_set_value_to_child_node_location(tParentTopology, j,tLocalCoordinate,tLevelSetValues);
                mNodePhaseVals(i+tNumCurrNodes,j) = tLevelSetValues(0,0);
            }

        }
    }

    /**
     * Links new nodes with an existing geometry object. This is used for unzipped interfaces
     * where more than one node is at the same location
     * @param[in] aNodesIndicesWithGeomObj - Node indices which already have a geometry object
     * @param[in] aNodesIndicesToLink - Node indices to link to the corresponding nodes in aNodesIndicesWithGeomObj
     */

    void
    link_new_nodes_to_existing_geometry_objects( Matrix< IndexMat > const & aNodesIndicesWithGeomObj,
                                                 Matrix< IndexMat > const & aNodesIndicesToLink)
    {
        // Assert lengths match
        MORIS_ASSERT(aNodesIndicesWithGeomObj.numel() == aNodesIndicesToLink.numel(),
        "Length of nodes with geometry objects does not match length of list with node indices to link  ");

        // Number of nodes to link
        uint tNumNodes = aNodesIndicesWithGeomObj.numel();

        // Iterate through nodes and create the link
        for(uint i = 0; i <tNumNodes; i++)
        {
            mGeometryObjects.link_to_node_to_another_nodes_geometry_object(aNodesIndicesWithGeomObj(i),aNodesIndicesToLink(i));
        }


    }


    /** is_intersected checks to see if an entity provided to it intersects a geometry field. Intersects in this context
     * means a geometry crosses a certain threshold (typically 0). For levelset fields, this can be thought of as a phase change
     *
     * @param[in] aNodeCoords       - Node coordinate
     * @param[in] aNodeToEntityConn - Connectivity between nodes and parent entity
     * @param[in] aCheckType        - Specifies what type of intersection check is to be performed
     *                                   0 - No information on interface required
     *                                   1 - information on interface required
     */
    void is_intersected(moris::Matrix< moris::DDRMat > const &   aNodeCoords,
                        moris::Matrix< moris::IndexMat > const & aNodetoEntityConn,
                        moris::size_t                            aCheckType,
                        Cell<Geometry_Object> &                  aGeometryObjects)
    {
        //Get information for loops
        moris::size_t tNumEntities = aNodetoEntityConn.n_rows(); // Number of entities provided to the geometry engine

        //Initialize
        moris::size_t tIntersectedCount = 0;    // Intersected element counter
        aGeometryObjects.clear();
        aGeometryObjects.resize(tNumEntities,Geometry_Object());

        //Loop over elements and determine if the element has an intersection
        for(moris::moris_index i = 0; i < (moris::moris_index)tNumEntities; i++)
        {

            //Populate the intersection flag of this element with a bool
            moris::Matrix< moris::IndexMat > tRow = aNodetoEntityConn.get_row(i);
            moris::Matrix< moris::IndexMat > tNodeADVIndices;
            bool tIsIntersected = compute_intersection_info(i,tRow, aNodeCoords, aCheckType,tNodeADVIndices,aGeometryObjects(tIntersectedCount));

            if(tIsIntersected)
            {
                tIntersectedCount++;
            }
        }

        // resize
        aGeometryObjects.resize(tIntersectedCount, Geometry_Object());
    }


    /*!
     * Computes the interface sensitivity of the provided node indices. After this call,
     * the sensitivity information of these interface nodes can be accessed through the interface
     * nodes respective geometry object.
     * @param[in] aInterfaceNodeIndices - Interface Node Indices (should be interface nodes wrt geometry index provided)
     * @param[in] aNodeCoords -  Node coordinates with location corresponding to indices of aIntefaceNodeIndices.
     * @param[in] aGeomIndex - Geometry Index
     */
    void
    compute_interface_sensitivity( Matrix< IndexMat > const & aInterfaceNodeIndices,
                                   Matrix< DDRMat >   const & aNodeCoords,
                                   moris_index                aGeomIndex)
    {
        // Figure out how many entities to compute sensitivity for
        uint tNumEntities = aInterfaceNodeIndices.numel();

        // iterate through node indices and compute sensitivity for each
        for(uint iEnt = 0; iEnt<tNumEntities; iEnt++)
        {
            // get the node index
            moris::moris_index tNodeIndex = aInterfaceNodeIndices(iEnt);

            // Get the node geometry object
            Geometry_Object & tGeoObj = this->get_geometry_object(tNodeIndex);

            // Get the parent topology that this node was created on
            Topology const & tParentEdge = tGeoObj.get_parent_entity_topology();

            MORIS_ASSERT(tParentEdge.get_topology_type() == Topology_Type::EDGE,"Only supporting interface sensitivity computation on an edge");

            // Get the node indices from the topology
            Matrix< IndexMat > const & tParentEntityNodes = tParentEdge.get_node_indices();

            // Initialize sensitivity
            Matrix< DDRMat > tDxDp(1,1,0.0);

            // Get the node vars of the parent edge nodes
            Matrix< DDRMat > tEntityNodeVars(tParentEntityNodes.numel(),1);
            for(uint i = 0; i < tParentEntityNodes.numel(); i++)
            {
                tEntityNodeVars(i) = this->get_entity_phase_val(tParentEntityNodes(i),aGeomIndex);
            }


            // Recompute local intersection (This could be stored instead)
            Matrix< DDRMat > tIntersectLocalCoordinate(1,1,0.0);
            Matrix< DDRMat > tIntersectGlobalCoordinate(1,1,0.0);
            get_intersection_location(mThresholdValue,
                                      mPerturbationValue,
                                      aNodeCoords,
                                      tEntityNodeVars,
                                      tParentEntityNodes,
                                      tIntersectLocalCoordinate,
                                      tIntersectGlobalCoordinate,
                                      true,
                                      false);

            // FIXME: Parent edge nodes need to not be the ADVs
            Matrix< IndexMat > tADVIndices;

            compute_dx_dp_for_an_intersection(tParentEntityNodes,aNodeCoords,tIntersectLocalCoordinate,tEntityNodeVars, tDxDp, tADVIndices);

            tGeoObj.set_sensitivity_dx_dp(tDxDp);
            tGeoObj.set_node_adv_indices(tParentEntityNodes);
        }

    }

    /*
     * Computes the intersection of an isocountour with an entity and returning the local coordinate relative to the parent
     * and the global coordinate if needed
     */
    void
    get_intersection_location(moris::real const &                      aIsocontourThreshold,
                              moris::real const &                      aPerturbationThreshold,
                              moris::Matrix< moris::DDRMat > const &   aGlobalNodeCoordinates,
                              moris::Matrix< moris::DDRMat > const &   aEntityNodeVars,
                              moris::Matrix< moris::IndexMat > const & aEntityNodeIndices,
                              moris::Matrix< moris::DDRMat > &         aIntersectionLocalCoordinates,
                              moris::Matrix< moris::DDRMat > &         aIntersectionGlobalCoordinates,
                              bool                                     aCheckLocalCoordinate = true,
                              bool                                     aComputeGlobalCoordinate = false)
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
            moris::Matrix< moris::DDRMat > tEntityCoordinates(2,3);
            replace_row(aEntityNodeIndices(0,0), aGlobalNodeCoordinates,0,tEntityCoordinates);
            replace_row(aEntityNodeIndices(0,1), aGlobalNodeCoordinates,1,tEntityCoordinates);

            // compute the global coordinate
            Interpolation::linear_interpolation_location(tEntityCoordinates,aIntersectionLocalCoordinates,aIntersectionGlobalCoordinates);
        }
    }


    void
    compute_dx_dp_finite_difference(moris::real                      const & aPerturbationVal,
                                    moris::Matrix< moris::DDRMat >   const & aGlobalNodeCoordinates,
                                    moris::Matrix< moris::DDRMat >   const & aEntityNodeCoordinates,
                                    moris::Matrix< moris::DDRMat >   const & aIntersectionLclCoordinate,
                                    moris::Matrix< moris::IndexMat > const & aEntityNodeIndices,
                                    moris::Matrix< moris::DDRMat >         & aEntityNodeVars,
                                    moris::Matrix< moris::DDRMat >         & aDxDp)
    {

        moris::size_t tNumNodeVars = aEntityNodeVars.n_rows();
        MORIS_ASSERT(tNumNodeVars == 2,"Currently compute_dx_dp_finite_difference has only been tested on edges");
        aDxDp.resize(2,3);

        moris::real tPerturbationLen = 2*aPerturbationVal;
        moris::real tScale   = 1/tPerturbationLen;
        Cell<moris::real>  tPerturbationSign = {1,-1};

        moris::Matrix< moris::DDRMat >       tDxDp(1,3);
        moris::Matrix< moris::DDRMat >       tPerturbedLocalCoordinate(1,1);
        Cell<moris::Matrix< moris::DDRMat >> tPerturbedGlobCoordinates = {moris::Matrix< moris::DDRMat >(1,3),
                                                                 moris::Matrix< moris::DDRMat >(1,3)};
        // Loop over all the nodes and perturb up and down
        for(moris::size_t i = 0; i<tNumNodeVars; i++)
        {
            // Perturb up and down
            for(moris::size_t j = 0; j<2; j++)
            {

                moris::real tPerturb = tPerturbationSign(j)*aPerturbationVal;
                // Perturb
                aEntityNodeVars(i,0) = aEntityNodeVars(i,0) + tPerturb;

                // Locate perturbed interface
                get_intersection_location(mThresholdValue, aPerturbationVal, aGlobalNodeCoordinates, aEntityNodeVars, aEntityNodeIndices, tPerturbedLocalCoordinate, tPerturbedGlobCoordinates(j),false, true);

                // Reverse perturb
                aEntityNodeVars(i,0) = aEntityNodeVars(i,0) - tPerturb;

            }

            tDxDp.matrix_data() = tScale*(tPerturbedGlobCoordinates(1).matrix_data() - tPerturbedGlobCoordinates(0).matrix_data());

            replace_row(0,tDxDp,i,aDxDp);


        }
    }


    void
    compute_dx_dp_for_an_intersection(moris::Matrix< moris::IndexMat > const & aEntityNodeIndices,
                                      moris::Matrix< moris::DDRMat >   const & aGlobalNodeCoordinates,
                                      moris::Matrix< moris::DDRMat >   const & aIntersectionLclCoordinate,
                                      moris::Matrix< moris::DDRMat >         & aEntityNodeVars,
                                      moris::Matrix< moris::DDRMat >         & aDxDp,
                                      moris::Matrix< moris::IndexMat >       & aADVIndices)
    {
        moris::real tPerturbationLength = 0.005;
        moris::size_t tNumNodes = aEntityNodeIndices.n_cols();

        MORIS_ASSERT(tNumNodes == 2,"Currently, this function is only supported on edges");

        // Initialize
        Cell<moris::Matrix< moris::DDRMat >> tDPhiiDp = {moris::Matrix< moris::DDRMat >(0,0),
                                                moris::Matrix< moris::DDRMat >(0,0)};

        moris::Matrix< moris::DDRMat > tEntityNodeCoordinates(2,3);

        // Assemble the entity local coordinates
        replace_row(aEntityNodeIndices(0,0), aGlobalNodeCoordinates,0,tEntityNodeCoordinates);
        replace_row(aEntityNodeIndices(0,1), aGlobalNodeCoordinates,1,tEntityNodeCoordinates);

        // Get information from a analytic geometry
        if(ActiveGeometry().is_analytic())
        {
            aADVIndices = get_node_adv_indices_analytic();
            for(moris::size_t i = 0; i < tNumNodes; i++)
            {
                tDPhiiDp(i) = ActiveGeometry().evaluate_sensitivity_dphi_dp_with_coordinate(aEntityNodeIndices(0, i),aGlobalNodeCoordinates);
            }

            compute_dx_dp_with_linear_basis(tDPhiiDp(0), tDPhiiDp(1), tEntityNodeCoordinates, aEntityNodeVars, aDxDp);
        }

        // Get information from a discrete geometry
        else
        {
            aADVIndices = get_node_adv_indices_discrete(aEntityNodeIndices);
            compute_dx_dp_finite_difference(tPerturbationLength, aGlobalNodeCoordinates, tEntityNodeCoordinates, aIntersectionLclCoordinate,aEntityNodeIndices,aEntityNodeVars, aDxDp);
        }
    }

    /**
     * Returns a reference to the geometry object at the provided index
     */
    Geometry_Object &
    get_geometry_object(moris::size_t const & aNodeIndex)
    {
       return mGeometryObjects.get_geometry_object_from_manager(aNodeIndex);
    }

    /**
     * Returns a reference to the geometry object at the provided index
     */
    Geometry_Object const &
    get_geometry_object(moris::size_t const & aNodeIndex) const
    {
       return mGeometryObjects.get_geometry_object_from_manager(aNodeIndex);
    }

    /*
     * Get the total number of phases in the phase table
     */
    moris::size_t get_num_phases()
    {
        return mPhaseTable.get_num_phases();
    }

    /*
     * Get the 0 or 1 value associated with a given phase and geometry index
     */
    moris::moris_index
    get_phase_sign_of_given_phase_and_geometry(moris::moris_index aPhaseIndex,
                                               moris::moris_index aGeometryIndex)
    {
        return mPhaseTable.get_phase_sign_of_given_phase_and_geometry(aPhaseIndex,aGeometryIndex);
    }

    /*
     * Get phase value for a given node and geometry index
     */
    moris::real
    get_entity_phase_val(moris::size_t const & aNodeIndex,
                         moris::size_t const & aGeomIndex)
    {
        Geometry_Object & tNodesGeoObj = get_geometry_object(aNodeIndex);
        moris::size_t tNodeRowIndex = tNodesGeoObj.get_phase_val_row();

        return mNodePhaseVals(tNodeRowIndex,aGeomIndex);
    }

    /*
     * Get dxdp for a node
     */
    moris::Matrix< moris::DDRMat > const &
    get_node_dx_dp(moris::size_t const & aNodeIndex) const
    {
        Geometry_Object const & tNodesGeoObj = get_geometry_object(aNodeIndex);
        return tNodesGeoObj.get_sensitivity_dx_dp();
    }

    moris::Matrix< moris::IndexMat > const &
    get_node_adv_indices(moris::size_t const & aNodeIndex) const
    {
        Geometry_Object const & tNodesGeoObj = get_geometry_object(aNodeIndex);
        return tNodesGeoObj.get_node_adv_indices();
    }



    /*
     * For a given, node index return the phase index relative to each geometry (i.e. inside/outside indicator)
     */
    void get_phase_index(moris::Matrix< moris::DDSTMat > const & aNodeIndex,
                         moris::Matrix< moris::DDSTMat > & aNodePhaseIndex)
    {
        // 0 for neg 1 for pos
        moris::real tNodePhaseValue = 0;
        moris::Matrix< moris::IndexMat > tPhaseOnOff(1,mGeometry.size());

        for(moris::size_t i = 0; i<aNodeIndex.n_cols(); i++)
        {
            Geometry_Object & tNodesGeoObj = get_geometry_object(aNodeIndex(0,i));
            moris::size_t tNodeRowIndex = tNodesGeoObj.get_phase_val_row();

            for(moris::size_t iG = 0; iG<mGeometry.size(); iG++)
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
            aNodePhaseIndex(i,0) = mPhaseTable.get_phase_index(tPhaseOnOff);
        }
    }

    /*
      * For a given, node index return the phase index relative to each geometry (i.e. inside/outside indicator)
      */
     void get_phase_index(moris::moris_index const & aNodeIndex,
                          moris::size_t & aNodePhaseIndex)
     {
         // 0 for neg 1 for pos
         moris::real tNodePhaseValue = 0;
         moris::Matrix< moris::IndexMat > tPhaseOnOff(1,mGeometry.size());

         Geometry_Object & tNodesGeoObj = get_geometry_object(aNodeIndex);
         moris::size_t tNodeRowIndex = tNodesGeoObj.get_phase_val_row();

         for(moris::size_t iG = 0; iG<mGeometry.size(); iG++)
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

    /*
     * Provided the inside and out phase values for an entity, return the phase index
     */
    moris::moris_index
    get_elem_phase_index(moris::Matrix< moris::IndexMat > const & aElemOnOff)
    {
        return mPhaseTable.get_phase_index(aElemOnOff);
    }

    /*
     * Returns whether a node is inside or outside wrt to a given geometry index
     */
    moris::size_t
    get_node_phase_index_wrt_a_geometry(moris::size_t aNodeIndex,
                                        moris::size_t aGeometryIndex)
    {
        Geometry_Object & tNodesGeoObj = get_geometry_object(aNodeIndex);
        moris::size_t tNodeRowIndex = tNodesGeoObj.get_phase_val_row();

        moris::real tNodePhaseVal = mNodePhaseVals(tNodeRowIndex,aGeometryIndex);

        moris::size_t tPhaseOnOff = 10000;

        if(tNodePhaseVal<mThresholdValue)
        {
            tPhaseOnOff = 0;
        }
        else
        {
            tPhaseOnOff = 1;
        }

        return tPhaseOnOff;
    }

    /*
     * Returns whether the active geometry is analytic
     */
    bool is_geometry_analytic()
    {
        return ActiveGeometry().is_analytic();
    }

    /*
     * Returns the number of geometries
     */
    moris::size_t get_num_geometries()
    {
        return mGeometry.size();
    }

    /*
     * Returns the number of geometries
     */
    moris::size_t get_num_bulk_phase()
    {
        return mPhaseTable.get_num_phases();
    }

    /*
     * Returns the active geometry index
     */
    moris::size_t get_active_geometry_index()
    {
        return mActiveGeometryIndex;
    }

    /*
     * Advance the active geometry index
     */
    void advance_geometry_index()
    {
        MORIS_ASSERT(mActiveGeometryIndex<get_num_geometries(),"Trying to advance past the number of geometries in the geometry engine");
        mActiveGeometryIndex += 1;
    }

    moris::Matrix< moris::IndexMat >
    get_node_adv_indices_analytic()
    {
        moris::size_t tNumDVS = get_num_design_vars_analytic();
        moris::Matrix< moris::IndexMat > tMatrix(1,tNumDVS);
        for(moris::size_t i = 0; i<tNumDVS; i++)
        {
            tMatrix(0,i) = (moris::moris_index)i;
        }
        return tMatrix;
    }


    /*
     * Returns the ADV indices of the provided nodes
     */
    moris::Matrix< moris::IndexMat >
    get_node_adv_indices_discrete(moris::Matrix< moris::IndexMat > const & aEntityNodes)
    {
        //FIXME: use this line
//        moris::Matrix< moris::IndexMat > tNodeADVIndices = ActiveGeometry().get_node_adv_indices(aEntityNodes);

        return aEntityNodes.copy();
    }

    moris::size_t
    get_num_design_vars_analytic()
    {
        moris::size_t tNumRows = 0;
        moris::size_t tNumCols = 0;
        moris::size_t tNumDVs = 0;
        for(moris::size_t i = 0; i<mGeometry.size(); i++)
        {
            mGeometry(i)->get_dphi_dp_size(tNumRows,tNumCols);
            tNumDVs += tNumCols;
        }

        return tNumDVs;
    }


private:
    moris::size_t mActiveGeometryIndex;
    Cell<Geometry*> mGeometry;

    // Contains all the geometry objects
    Geometry_Object_Manager mGeometryObjects;

    // Phase Table
    Phase_Table mPhaseTable;

    // Node Entity Phase Vals
    // Only analytic phase values are stored here to prevent duplicate storage of discrete geometries
    moris::Matrix< moris::DDRMat > mNodePhaseVals;



private:
    Geometry &
    ActiveGeometry() const
    {
        return (*mGeometry(mActiveGeometryIndex));
    }




    /**
     * compute_intersection_info, calculates the relevant intersection information placed in the geometry object
     * @param[in]  aEntityNodeInds - node to entity connectivity
     * @param[in]  aNodeVars      - node level set values
     * @param[in]  aCheckType     - if a entity local location is necessary 1, else 0.
     * @param[out] Returns an intersection flag and local coordinates if aCheckType 1 in cell 1 and node sensitivity information in cell 2 if intersection point located
     **/
    bool
    compute_intersection_info(moris::moris_index               const & aEntityIndex,
                              moris::Matrix< moris::IndexMat > const & aEntityNodeInds,
                              moris::Matrix< moris::DDRMat >   const & aNodeCoords,
                              moris::size_t const &                    aCheckType,
                              moris::Matrix< moris::IndexMat > &       aNodeADVIndices,
                              Geometry_Object & aGeometryObject)
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
        moris::size_t tNumNodes = aEntityNodeInds.n_cols();
        moris::Matrix< moris::DDRMat > tEntityNodeVars(tNumNodes, 1);
        moris::Matrix< moris::DDRMat > tInterpLocationCoords(1,1);


        // Loop through nodes and get levelset values from precomputed values in aNodeVars or in the levelset mesh
        for(moris::size_t n = 0; n < tNumNodes; n++)
        {   //Get node id n
            tNodeInd = aEntityNodeInds(0, n);

            Geometry_Object & tGeoObj = get_geometry_object(tNodeInd);
            moris::size_t tPhaseValRowIndex = tGeoObj.get_phase_val_row();
            tEntityNodeVars(n, 0) = mNodePhaseVals(tPhaseValRowIndex, mActiveGeometryIndex);
        }

        //get the max and minimum levelset value for the entity
        tMax = tEntityNodeVars.max(tMaxLocRow,tMaxLocCol);
        tMin = tEntityNodeVars.min(tMinLocRow,tMinLocCol);

        //    If there is a sign change in element node variables return true, else return false

        //TODO: intersection flag should not be a moris::real (needs to be a bool) split this function
        moris::Matrix< moris::DDRMat > tIntersection(1, 2, 0.0);// Initialize as false

        moris::real tErrorFactor = 1;
        // If the max is also the threshold value, figure out which node is on the interface
        // If both are the interface, the location may not be correct (sincent
        if( approximate(tMin, mThresholdValue,tErrorFactor) && approximate(tMax, mThresholdValue,tErrorFactor))
        {
            aGeometryObject.set_parent_entity_index(aEntityIndex);
            aGeometryObject.mark_all_nodes_as_on_interface();
            tIsIntersected = true;
        }

        else if(approximate(tMax,mThresholdValue,tErrorFactor))
        {
            aGeometryObject.set_parent_entity_index(aEntityIndex);
            aGeometryObject.mark_node_as_on_interface(tMaxLocRow);
            tIsIntersected = true;
        }

        // If the min is also the threshold value, figure out which node is on the interface
        else if(approximate(tMin,mThresholdValue,tErrorFactor))
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
                moris::Matrix< moris::DDRMat > tIntersectGlobalCoordinate(1,3);
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

    void
    interpolate_level_set_value_to_child_node_location(Topology const & aParentTopology,
                                                       moris::size_t const &                                              aGeometryIndex,
                                                       moris::Matrix< moris::DDRMat > const &                                aNodeLocalCoordinate,
                                                       moris::Matrix< moris::DDRMat > & aLevelSetValues)
     {

         // Get node indices attached to parent (These are indices relative to another mesh and may need to be mapped)
         moris::Matrix< moris::IndexMat > const & tNodesAttachedToParent = aParentTopology.get_node_indices();

         // Get number of nodes attached to parent
         moris::size_t tNumNodesAttachedToParent = tNodesAttachedToParent.numel();
         moris::Matrix< moris::DDRMat > tNodesLevelSetValues(1, tNumNodesAttachedToParent);

         for(moris::size_t i = 0; i < tNumNodesAttachedToParent; i++)
         {
             Geometry_Object & tGeoObj = get_geometry_object(tNodesAttachedToParent(i));
             moris::size_t tPhaseRow = tGeoObj.get_phase_val_row();

             tNodesLevelSetValues(0,i) = mNodePhaseVals(tPhaseRow,aGeometryIndex);
         }

         // Ask the topology how to interpolate
         moris::Matrix< moris::DDRMat > tBasisValues(1,1);
         Basis_Function const & tParentBasisFunctions = aParentTopology.get_basis_function();

         // Evaluate basis function
         tParentBasisFunctions.evaluate_basis_function(aNodeLocalCoordinate,tBasisValues);

         // Compute \phi = Ni.\phi_i
         aLevelSetValues = tBasisValues*moris::trans(tNodesLevelSetValues);

     }


};
}

#endif /* SRC_GEOMENG_CL_MGE_GEOMETRY_ENGINE_HPP_ */
