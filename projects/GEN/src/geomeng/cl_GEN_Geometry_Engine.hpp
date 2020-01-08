/*
 * cl_GEN_Geometry_Engine.hpp
 *
 *  Created on: Jun 21, 2017
 *      Author: ktdoble
 */

#ifndef PROJECTS_GEN_SRC_NEW_GEOMENG_CL_GEN_GEOMETRY_ENGINE_HPP_
#define PROJECTS_GEN_SRC_NEW_GEOMENG_CL_GEN_GEOMETRY_ENGINE_HPP_

// Standard library includes
#include <memory> // for shared_ptr
#include <math.h>

#include "cl_Cell.hpp"
#include "cl_Logger.hpp"

// Linear algebra includes
#include "cl_Matrix.hpp"
#include "fn_trans.hpp"
#include "op_equal_equal.hpp"
#include "op_times.hpp"
#include "linalg_typedefs.hpp"

// GE
#include "../additional/cl_GEN_Basis_Function.hpp"
#include "../additional/cl_GEN_Interpolaton.hpp"
#include "../additional/cl_GEN_Pending_Node.hpp"
#include "../additional/cl_GEN_Phase_Table.hpp"
#include "../additional/fn_GEN_approximate.hpp"

#include "../field/cl_GEN_Field.hpp"
#include "../field/cl_GEN_Field_User_Defined.hpp"

#include "cl_GEN_Geometry_Object.hpp"
#include "cl_GEN_Geometry_Object_Manager.hpp"
#include "cl_GEN_Pdv_Host.hpp"
#include "cl_GEN_Pdv_Host_Manager.hpp"

#include "../geometry/cl_GEN_Geometry.hpp"
#include "../geometry/cl_GEN_Cylinder_With_End_Caps.hpp"
#include "../geometry/cl_GEN_Geometry.hpp"

#include "../property/cl_GEN_Property.hpp"

// MTK
#include "cl_MTK_Mesh_Manager.hpp"
#include "cl_Mesh_Enums.hpp"

// HMR
#include "../projects/HMR/src/cl_HMR_Mesh.hpp"

namespace moris
{
namespace ge
{


/*
 *
 * $\frac{\partial{\phi_A}}{\partial{p}}$ (change in phi with respect to a design variable
 * See for more detailed description of this function:
 *
 */
inline
void compute_dx_dp_with_linear_basis( moris::Matrix< moris::DDRMat >  & aDPhiADp,
                                      moris::Matrix< moris::DDRMat >  & aDPhiBDp,
                                      moris::Matrix< moris::DDRMat >  & aEdgeCoordinates,
                                      moris::Matrix< moris::DDRMat >  & aEdgeNodePhi,
                                      moris::Matrix< moris::DDRMat >  & aDxDp )
{

  MORIS_ASSERT(aDPhiADp.n_rows() != 0,"dPhi/dp not implemented in geometry would cause a seg fault here");
  MORIS_ASSERT(aDPhiBDp.n_rows() != 0,"dPhi/dp not implemented in geometry would cause a seg fault here");
  moris::real const & tPhiA = aEdgeNodePhi(0,0);
  moris::real const & tPhiB = aEdgeNodePhi(1,0);

  // Initialize
  moris::Matrix< moris::DDRMat > tXa = aEdgeCoordinates.get_row(0);

  moris::Matrix< moris::DDRMat > tXb = aEdgeCoordinates.get_row(1);

  // ------- Compute $\frac{\partial x_{\Gamma}}{\partial \phi}$ -------
  moris::DDRMat tDxgammaDphiA = -(tPhiB)/std::pow((tPhiA-tPhiB),2)*(tXb.matrix_data()-tXa.matrix_data());
  moris::DDRMat tDxgammaDphiB =  (tPhiA)/std::pow((tPhiA-tPhiB),2)*(tXb.matrix_data()-tXa.matrix_data());

  moris::Matrix< moris::DDRMat > tDxgDphiAMat(tDxgammaDphiA);
  moris::Matrix< moris::DDRMat > tDxgDphiBMat(tDxgammaDphiB);
  // ------------------------------ end --------------------------------

  // Compute dx/dp
  moris::DDRMat tDxDp = aDPhiADp * moris::trans(tDxgDphiAMat) +  aDPhiBDp * moris::trans(tDxgDphiBMat);
  aDxDp = moris::Matrix< moris::DDRMat >(tDxDp);

}
//------------------------------------------------------------------------------

class GEN_Geometry_Engine
{
public:

    // Single geometry constructor
    GEN_Geometry_Engine( moris::ge::GEN_Geometry  & aGeometry,
                         moris::ge::GEN_Phase_Table const & aPhaseTable,
                         moris::uint aSpatialDim = 3 );

    // geometry vector constructor
    GEN_Geometry_Engine( Cell< GEN_Geometry* > const  & aGeometry,
                         moris::ge::GEN_Phase_Table const & aPhaseTable,
                         moris::uint aSpatialDim = 3 );

    GEN_Geometry_Engine(  )
    {
    }

    // Options which the user can change (all are given defaults)
    moris::real mThresholdValue;
    moris::real mPerturbationValue;
    bool        mComputeDxDp; // Should be turned off if a sensitivity has not been implemented
    moris::uint mSpatialDim;


    //------------------------------------------------------------------------------
    /*
     * @brief set the pdv type and associated property lists
     *
     * The user needs to be sure that the entry in aPdvList is associated with the same entry in aPropertyList.
     *
     * e.g.
     * aPdvList = { GEN_PDV::DENSITY, GEN_PDV::MODULUS }
     * aPropertyList = { property 1, property 2 }
     *
     * property 1 must correspond to DENSITY and property 2 must correspond to MODULUS
     */

    void set_pdv_property_list( const Cell< enum GEN_PDV >  aPdvList,
                                const Cell< GEN_Property* > aPropertyList );
    //------------------------------------------------------------------------------
    /*
     * @brief Initial allocation of geometry objects, this creates a geometry object for each node coordinate.
     * In this case, aNodeCoords needs to be ordered by proc indices.
     */
    void initialize_geometry_objects_for_background_mesh_nodes(moris::size_t const & aNumNodes);
    //------------------------------------------------------------------------------
    /*
     * @brief initial allocation of pdv hosts, creates a pdv host for each node coordinate
     */
    void initialize_pdv_hosts_for_background_mesh_nodes(moris::size_t const & aNumNodes);
    //------------------------------------------------------------------------------
    void initialize_geometry_object_phase_values(moris::Matrix< moris::DDRMat > const & aNodeCoords);
    //------------------------------------------------------------------------------
    /*
     * @brief Creates a geometry object association for pending nodes
     * These nodes have node indices and parent information
     */
    void associate_new_nodes_with_geometry_object( Cell<Pending_Node> & aNewNodes,
                                                   bool                 aInterfaceNodes );
    //------------------------------------------------------------------------------
    void associate_new_nodes_with_pdv_hosts( Cell< Pending_Node > & aNewNodes,
                                             bool                   aInterfaceNodes );

    //------------------------------------------------------------------------------
    void create_new_node_geometry_objects(Cell< moris_index >  const & aNewNodeIndices,
                                          bool                         aStoreParentTopo,
                                          Cell<xtk::Topology*> const & aParentTopo,
                                          Cell<Matrix<DDRMat>> const & aParamCoordRelativeToParent,
                                          Cell<Matrix<DDRMat>> const & aGlobalNodeCoord);
    //------------------------------------------------------------------------------
    void create_new_node_pdv_hosts( Cell< moris_index >  const & aNewNodeIndices,
                                    bool                         aStoreParentTopo,
                                    Cell<xtk::Topology*> const & aParentTopo,
                                    Cell<Matrix<DDRMat>> const & aGlobalNodeCoord );
    //------------------------------------------------------------------------------
    /**
     * @brief Links new nodes with an existing geometry object. This is used for unzipped interfaces
     * where more than one node is at the same location
     * @param[in] aNodesIndicesWithGeomObj - Node indices which already have a geometry object
     * @param[in] aNodesIndicesToLink - Node indices to link to the corresponding nodes in aNodesIndicesWithGeomObj
     */

    void link_new_nodes_to_existing_geometry_objects( Matrix< IndexMat > const & aNodesIndicesWithGeomObj,
                                                      Matrix< IndexMat > const & aNodesIndicesToLink );
    //------------------------------------------------------------------------------
    /**
     * @brief is_intersected checks to see if an entity provided to it intersects a geometry field. Intersects in this context
     * means a geometry crosses a certain threshold (typically 0). For levelset fields, this can be thought of as a phase change
     *
     * @param[in] aNodeCoords       - Node coordinate
     * @param[in] aNodeToEntityConn - Connectivity between nodes and parent entity
     * @param[in] aCheckType        - Specifies what type of intersection check is to be performed
     *                                   0 - No information on interface required
     *                                   1 - information on interface required
     */
    void is_intersected( moris::Matrix< moris::DDRMat > const &   aNodeCoords,
                         moris::Matrix< moris::IndexMat > const & aNodetoEntityConn,
                         moris::size_t                            aCheckType,
                         Cell<GEN_Geometry_Object> &              aGeometryObjects );
    //------------------------------------------------------------------------------
    /*!
     * @brief Computes the interface sensitivity of the provided node indices. After this call,
     * the sensitivity information of these interface nodes can be accessed through the interface
     * nodes respective geometry object.
     * @param[in] aInterfaceNodeIndices - Interface Node Indices (should be interface nodes wrt geometry index provided)
     * @param[in] aNodeCoords -  Node coordinates with location corresponding to indices of aIntefaceNodeIndices.
     * @param[in] aGeomIndex - Geometry Index
     * @param[in] aGlbCoord  - bool to calculate the global coordinate of the intersection point
     */
    void compute_interface_sensitivity( Matrix< IndexMat > const & aInterfaceNodeIndices,
                                        Matrix< DDRMat >   const & aNodeCoords,
                                        moris_index                aGeomIndex,
                                        bool               const   aGlbCoord = false );
    //------------------------------------------------------------------------------
    /*
     * @brief Computes the intersection of an isocountour with an entity and returning the local coordinate relative to the parent
     * and the global coordinate if needed
     */
    void get_intersection_location( moris::real const &                      aIsocontourThreshold,
                                    moris::real const &                      aPerturbationThreshold,
                                    moris::Matrix< moris::DDRMat > const &   aGlobalNodeCoordinates,
                                    moris::Matrix< moris::DDRMat > const &   aEntityNodeVars,
                                    moris::Matrix< moris::IndexMat > const & aEntityNodeIndices,
                                    moris::Matrix< moris::DDRMat > &         aIntersectionLocalCoordinates,
                                    moris::Matrix< moris::DDRMat > &         aIntersectionGlobalCoordinates,
                                    bool                                     aCheckLocalCoordinate = true,
                                    bool                                     aComputeGlobalCoordinate = false );
    //------------------------------------------------------------------------------
    void compute_dx_dp_finite_difference( moris::real                      const & aPerturbationVal,
                                          moris::Matrix< moris::DDRMat >   const & aGlobalNodeCoordinates,
                                          moris::Matrix< moris::DDRMat >   const & aEntityNodeCoordinates,
                                          moris::Matrix< moris::DDRMat >   const & aIntersectionLclCoordinate,
                                          moris::Matrix< moris::IndexMat > const & aEntityNodeIndices,
                                          moris::Matrix< moris::DDRMat >         & aEntityNodeVars,
                                          moris::Matrix< moris::DDRMat >         & aDxDp );
    //------------------------------------------------------------------------------
    void compute_dx_dp_for_an_intersection( moris::Matrix< moris::IndexMat > const & aEntityNodeIndices,
                                            moris::Matrix< moris::DDRMat >   const & aGlobalNodeCoordinates,
                                            moris::Matrix< moris::DDRMat >   const & aIntersectionLclCoordinate,
                                            moris::Matrix< moris::DDRMat >         & aEntityNodeVars,
                                            moris::Matrix< moris::DDRMat >         & aDxDp,
                                            moris::Matrix< moris::IndexMat >       & aADVIndices );
    //------------------------------------------------------------------------------
    /**
     * @brief Returns a reference to the geometry object at the provided index
     */
    GEN_Geometry_Object &
    get_geometry_object(moris::size_t const & aNodeIndex);
    //------------------------------------------------------------------------------
    /**
     * @brief Returns a reference to the geometry object at the provided index
     */
    GEN_Geometry_Object const &
    get_geometry_object(moris::size_t const & aNodeIndex) const;
    //------------------------------------------------------------------------------
    /*
     * @brief Get the total number of phases in the phase table
     */
    moris::size_t get_num_phases();
    //------------------------------------------------------------------------------
    /*
     * @brief Get the 0 or 1 value associated with a given phase and geometry index
     */
    moris::moris_index
    get_phase_sign_of_given_phase_and_geometry( moris::moris_index aPhaseIndex,
                                                moris::moris_index aGeometryIndex );
    //------------------------------------------------------------------------------
    /*
     * @brief Get phase value for a given node and geometry index
     */
    moris::real
    get_entity_phase_val( moris::size_t const & aNodeIndex,
                          moris::size_t const & aGeomIndex );
    //------------------------------------------------------------------------------
    /*
     * @brief Get dxdp for a node
     */
    moris::Matrix< moris::DDRMat > const &
    get_node_dx_dp(moris::size_t const & aNodeIndex) const;
    //------------------------------------------------------------------------------
    /*
     * @brief get adv indices for a node
     */
    moris::Matrix< moris::IndexMat > const &
    get_node_adv_indices( moris::size_t const & aNodeIndex ) const;
    //------------------------------------------------------------------------------
    /*
     * @brief For a given node index, return the phase index relative to each geometry (i.e. inside/outside indicator)
     */
    void get_phase_index( moris::Matrix< moris::DDSTMat > const & aNodeIndex,
                          moris::Matrix< moris::DDSTMat > & aNodePhaseIndex );
    //------------------------------------------------------------------------------
    /*
      * @brief For a given node index, return the phase index relative to each geometry (i.e. inside/outside indicator)
      */
     void get_phase_index( moris::moris_index const & aNodeIndex,
                           moris::size_t & aNodePhaseIndex );
     //------------------------------------------------------------------------------
    /*
     * @brief Provided the inside and out phase values for an entity, return the phase index
     */
    moris::moris_index
    get_elem_phase_index(moris::Matrix< moris::IndexMat > const & aElemOnOff)
    {
        return mPhaseTable.get_phase_index(aElemOnOff);
    }
    //------------------------------------------------------------------------------
    /*
     * @brief Returns whether a node is inside or outside wrt to a given geometry index
     */
    moris::size_t
    get_node_phase_index_wrt_a_geometry(moris::size_t aNodeIndex,
                                        moris::size_t aGeometryIndex);
    //------------------------------------------------------------------------------
    /*
     * @brief Returns whether the active geometry is analytic
     */
    bool is_geometry_analytic();
    //------------------------------------------------------------------------------
    /*
     * @brief Returns the number of geometries
     */
    moris::size_t get_num_geometries();
    //------------------------------------------------------------------------------
    /*
     * @brief Returns the number of phases
     */
    moris::size_t get_num_bulk_phase();
    //------------------------------------------------------------------------------
    /*
     * @brief Returns the active geometry index
     */
    moris::size_t get_active_geometry_index();
    //------------------------------------------------------------------------------
    /*
     * @brief Advance the active geometry index
     */
    void advance_geometry_index();
    //------------------------------------------------------------------------------
    moris::Matrix< moris::IndexMat > get_node_adv_indices_analytic();

    //------------------------------------------------------------------------------
    moris::uint get_num_design_variables() const;
    //------------------------------------------------------------------------------
    /*
     * @brief Returns the ADV indices of the provided nodes
     */
    moris::Matrix< moris::IndexMat >
    get_node_adv_indices_discrete(moris::Matrix< moris::IndexMat > const & aEntityNodes);
    //------------------------------------------------------------------------------
    moris::size_t
    get_num_design_vars_analytic();
    //------------------------------------------------------------------------------
    Pdv_Host_Manager* get_pdv_hosts();

    Geometry_Object_Manager* get_all_geom_obj();
    //------------------------------------------------------------------------------
    /*
     * @brief register a mesh to be used for later computation(s)
     */
    moris_index register_mesh( mtk::Mesh_Manager* aMesh );

    moris_index register_mesh( std::shared_ptr< moris::hmr::Mesh > aMesh ); //fixme: this needs to be deleted and the GE should only be able to register an mtk mesh pair
    //------------------------------------------------------------------------------
    /*
     * @brief register a field or cell of fields for later computation
     */
    moris_index register_field( GEN_Field* aField );

    void set_field_cell( moris::Cell< GEN_Field* > aFieldCell );
    //------------------------------------------------------------------------------
    /*
     * @brief calculates the field values at all nodes
     */
    //fixme set the field up in a similar way as the geometries (phase value table, field object manager, etc.)
    void calc_field_vals_at_nodes( moris_index        aMeshIndex,
                                   moris_index        aFieldIndex,
                                   Matrix< DDRMat > & aNodeVals,
                                   moris_index        aMeshIndexInManager = 0 );
    //------------------------------------------------------------------------------
    /*
     * @brief function specific to fiber problem
     */
    Matrix< DDRMat > get_cylinder_vals( moris_index aWhichMesh,
                                        GEN_CylinderWithEndCaps* aFiber,
                                        uint aNumberOfFibers ); //FIXME this is currently only setup to work with an HMR member mesh

    //------------------------------------------------------------------------------

    //------------------------------------------------------------------------------


private:
    GEN_Geometry & ActiveGeometry() const;
    //------------------------------------------------------------------------------
    /**
     * @brief compute_intersection_info, calculates the relevant intersection information placed in the geometry object
     * @param[in]  aEntityNodeInds - node to entity connectivity
     * @param[in]  aNodeVars       - node level set values
     * @param[in]  aCheckType      - if a entity local location is necessary 1, else 0.
     * @param[out] Returns an intersection flag and local coordinates if aCheckType 1 in cell 1 and node sensitivity information in cell 2 if intersection point located
     **/
    bool compute_intersection_info( moris::moris_index               const & aEntityIndex,
                                    moris::Matrix< moris::IndexMat > const & aEntityNodeInds,
                                    moris::Matrix< moris::DDRMat >   const & aNodeCoords,
                                    moris::size_t                    const & aCheckType,
                                    moris::Matrix< moris::IndexMat >       & aNodeADVIndices,
                                    GEN_Geometry_Object                    & aGeometryObject );
    //------------------------------------------------------------------------------
    void interpolate_level_set_value_to_child_node_location( xtk::Topology                  const & aParentTopology,
                                                             moris::size_t                  const & aGeometryIndex,
                                                             moris::Matrix< moris::DDRMat > const & aNodeLocalCoordinate,
                                                             moris::Matrix< moris::DDRMat >       & aLevelSetValues );
    //------------------------------------------------------------------------------
private:    // member data
    moris::size_t mActiveGeometryIndex;
    Cell< GEN_Geometry* > mGeometry;

    // Contains all the geometry objects
    Geometry_Object_Manager mGeometryObjects;

    // Contains all the pdv hosts
    Pdv_Host_Manager mPdvHosts;

    // Phase Table
    moris::ge::GEN_Phase_Table mPhaseTable;

    // Node Entity Phase Vals
    // Only analytic phase values are stored here to prevent duplicate storage of discrete geometries
    moris::Matrix< moris::DDRMat > mNodePhaseVals;
    //------------------------------------------------------------------------------
    moris::Cell< mtk::Mesh_Manager* > mMesh;

    moris::Cell< std::shared_ptr< moris::hmr::Mesh > > mMesh_HMR; //FIXME delete this one
    //------------------------------------------------------------------------------

    moris::Cell< GEN_Field* > mFields;
    //------------------------------------------------------------------------------
    Cell< enum GEN_PDV >  mPdvList{{ GEN_PDV::END_ENUM }};
    Cell< GEN_Property* > mPropertyList{{ nullptr }};
    //------------------------------------------------------------------------------

};
}
}

#endif /* PROJECTS_GEN_SRC_NEW_GEOMENG_CL_GEN_GEOMETRY_ENGINE_HPP_ */
