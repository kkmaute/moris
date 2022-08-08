#include <iostream>
// FEM/INT/src
#include "cl_FEM_Element_Double_Sideset.hpp"
#include "cl_FEM_Set.hpp"
#include "cl_FEM_Model.hpp"
#include "cl_FEM_Field_Interpolator_Manager.hpp"
#include "fn_FEM_Side_Coordinate_Map.hpp"
// FEM/MSI/src
#include "cl_MSI_Equation_Model.hpp"

namespace moris
{
    namespace fem
    {

        //------------------------------------------------------------------------------

        Element_Double_Sideset::Element_Double_Sideset(
                mtk::Cell const *  aMasterIGCell,
                mtk::Cell const *  aSlaveIGCell,
                Set*               aSet,
                Cluster*           aCluster,
                moris::moris_index aCellIndexInCluster )
                : Element(
                        aMasterIGCell,
                        aSlaveIGCell,
                        aSet,
                        aCluster,
                        aCellIndexInCluster )
        {
        }

        //------------------------------------------------------------------------------

        Element_Double_Sideset::~Element_Double_Sideset() {}

        //------------------------------------------------------------------------------

        void
        Element_Double_Sideset::init_ig_geometry_interpolator(
                uint aMasterSideOrdinal,
                uint aSlaveSideOrdinal )
        {
            // get master IG geometry interpolator
            Geometry_Interpolator* tMasterIGGI =
                    mSet->get_field_interpolator_manager( mtk::Master_Slave::MASTER )->    //
                    get_IG_geometry_interpolator();

            // get slave IG geometry interpolator
            Geometry_Interpolator* tSlaveIGGI =
                    mSet->get_field_interpolator_manager( mtk::Master_Slave::SLAVE )->    //
                    get_IG_geometry_interpolator();

            // get master physical space and time coordinates for IG element
            Matrix< DDRMat > tMasterIGPhysSpaceCoords =
                    mMasterCell->get_cell_physical_coords_on_side_ordinal( aMasterSideOrdinal );
            Matrix< DDRMat > tMasterIGPhysTimeCoords =
                    mCluster->mInterpolationElement->get_time();

            // get slave physical space and time coordinates for IG element
            Matrix< DDRMat > tSlaveIGPhysSpaceCoords =
                    mSlaveCell->get_cell_physical_coords_on_side_ordinal( aSlaveSideOrdinal );
            Matrix< DDRMat > tSlaveIGPhysTimeCoords =
                    mCluster->mInterpolationElement->get_time();

            // get master parametric space and time coordinates for IG element
            Matrix< DDRMat > tMasterIGParamSpaceCoords =
                    mCluster->get_cell_local_coords_on_side_wrt_interp_cell(
                            mCellIndexInCluster,
                            aMasterSideOrdinal,
                            mtk::Master_Slave::MASTER );
            // FIXME not true if time is not linear
            Matrix< DDRMat > tMasterIGParamTimeCoords = { { -1.0 }, { 1.0 } };

            // get slave parametric space and time coordinates for IG element
            Matrix< DDRMat > tSlaveIGParamSpaceCoords =
                    mCluster->get_cell_local_coords_on_side_wrt_interp_cell(
                            mCellIndexInCluster,
                            aSlaveSideOrdinal,
                            mtk::Master_Slave::SLAVE );
            // FIXME not true if time is not linear
            Matrix< DDRMat > tSlaveIGParamTimeCoords = { { -1.0 }, { 1.0 } };

            // set physical space and time coefficients for master IG element GI
            tMasterIGGI->set_space_coeff( tMasterIGPhysSpaceCoords );
            tMasterIGGI->set_time_coeff( tMasterIGPhysTimeCoords );

            // set physical space and time coefficients for slave IG element GI
            tSlaveIGGI->set_space_coeff( tSlaveIGPhysSpaceCoords );
            tSlaveIGGI->set_time_coeff( tSlaveIGPhysTimeCoords );

            // set parametric space and time coefficients for master IG element GI
            tMasterIGGI->set_space_param_coeff( tMasterIGParamSpaceCoords );
            tMasterIGGI->set_time_param_coeff( tMasterIGParamTimeCoords );

            // set parametric space and time coefficients for slave IG element GI
            tSlaveIGGI->set_space_param_coeff( tSlaveIGParamSpaceCoords );
            tSlaveIGGI->set_time_param_coeff( tSlaveIGParamTimeCoords );
        }

        //------------------------------------------------------------------------------

        void
        Element_Double_Sideset::init_ig_geometry_interpolator(
                uint              aMasterSideOrdinal,
                uint              aSlaveSideOrdinal,
                Matrix< DDSMat >& aGeoLocalAssembly )
        {
            // get master IG geometry interpolator
            Geometry_Interpolator* tMasterIGGI =
                    mSet->get_field_interpolator_manager( mtk::Master_Slave::MASTER )->    //
                    get_IG_geometry_interpolator();

            // get slave IG geometry interpolator
            Geometry_Interpolator* tSlaveIGGI =
                    mSet->get_field_interpolator_manager( mtk::Master_Slave::SLAVE )->    //
                    get_IG_geometry_interpolator();

            // get master physical space and time coordinates for IG element
            Matrix< DDRMat > tMasterIGPhysSpaceCoords =
                    mMasterCell->get_cell_physical_coords_on_side_ordinal( aMasterSideOrdinal );

            Matrix< DDRMat > tMasterIGPhysTimeCoords =
                    mCluster->mInterpolationElement->get_time();

            // get slave physical space and time coordinates for IG element
            Matrix< DDRMat > tSlaveIGPhysSpaceCoords =
                    mSlaveCell->get_cell_physical_coords_on_side_ordinal( aSlaveSideOrdinal );

            Matrix< DDRMat > tSlaveIGPhysTimeCoords =
                    mCluster->mInterpolationElement->get_time();

            // get master parametric space and time coordinates for IG element
            Matrix< DDRMat > tMasterIGParamSpaceCoords =
                    mCluster->get_cell_local_coords_on_side_wrt_interp_cell(
                            mCellIndexInCluster,
                            aMasterSideOrdinal,
                            mtk::Master_Slave::MASTER );

            // FIXME not true if time is not linear
            Matrix< DDRMat > tMasterIGParamTimeCoords = { { -1.0 }, { 1.0 } };

            // get slave parametric space and time coordinates for IG element
            Matrix< DDRMat > tSlaveIGParamSpaceCoords =
                    mCluster->get_cell_local_coords_on_side_wrt_interp_cell(
                            mCellIndexInCluster,
                            aSlaveSideOrdinal,
                            mtk::Master_Slave::SLAVE );

            // FIXME not true if time is not linear
            Matrix< DDRMat > tSlaveIGParamTimeCoords = { { -1.0 }, { 1.0 } };

            // get the local cluster assembly indices
            if ( mSet->get_geo_pdv_assembly_flag() )
            {
                // get the vertices indices
                Matrix< IndexMat > tMasterVertexIndices =
                        mMasterCell->get_vertices_ind_on_side_ordinal( aMasterSideOrdinal );
                
                Matrix< IndexMat > tSlaveVertexIndices =
                        mSlaveCell->get_vertices_ind_on_side_ordinal( aSlaveSideOrdinal );

                // get the requested geo pdv types
                moris::Cell< enum PDV_Type > tGeoPdvType;
                mSet->get_ig_unique_dv_types_for_set( tGeoPdvType );

                // get local assembly indices
                mSet->get_equation_model()->get_integration_xyz_pdv_assembly_indices(
                        tMasterVertexIndices,
                        tGeoPdvType,
                        aGeoLocalAssembly );
            }

            // set physical space and time coefficients for master IG element GI
            tMasterIGGI->set_space_coeff( tMasterIGPhysSpaceCoords );
            tMasterIGGI->set_time_coeff( tMasterIGPhysTimeCoords );

            // set physical space and time coefficients for slave IG element GI
            tSlaveIGGI->set_space_coeff( tSlaveIGPhysSpaceCoords );
            tSlaveIGGI->set_time_coeff( tSlaveIGPhysTimeCoords );

            // set parametric space and time coefficients for master IG element GI
            tMasterIGGI->set_space_param_coeff( tMasterIGParamSpaceCoords );
            tMasterIGGI->set_time_param_coeff( tMasterIGParamTimeCoords );

            // set parametric space and time coefficients for slave IG element GI
            tSlaveIGGI->set_space_param_coeff( tSlaveIGParamSpaceCoords );
            tSlaveIGGI->set_time_param_coeff( tSlaveIGParamTimeCoords );
        }

        //------------------------------------------------------------------------------

        void
        Element_Double_Sideset::compute_residual()
        {
            // get number of IWGs
            uint tNumIWGs = mSet->get_number_of_requested_IWGs();

            // check for active IWGs
            if ( tNumIWGs == 0 )
            {
                return;
            }

            // get treated side ordinal on the master and on the slave
            uint tMasterSideOrd = mCluster->mMasterListOfSideOrdinals( mCellIndexInCluster );
            uint tSlaveSideOrd  = mCluster->mSlaveListOfSideOrdinals( mCellIndexInCluster );

            // set the master/slave ig geometry interpolator physical/parametric space and time coefficients
            this->init_ig_geometry_interpolator( tMasterSideOrd, tSlaveSideOrd );

            // get first corresponding node from master to slave
            moris::mtk::Vertex const * tSlaveNode =
                    mCluster->get_left_vertex_pair(
                            mMasterCell->get_vertices_on_side_ordinal( tMasterSideOrd )( 0 ) );

            moris_index tSlaveNodeOrdOnSide =
                    mCluster->get_right_vertex_ordinal_on_facet( mCellIndexInCluster, tSlaveNode );

            // loop over the integration points
            uint tNumIntegPoints = mSet->get_number_of_integration_points();

            for ( uint iGP = 0; iGP < tNumIntegPoints; iGP++ )
            {
                // get local integration point for the master integration cell
                const Matrix< DDRMat >& tMasterLocalIntegPoint = mSet->get_integration_points().get_column( iGP );

                // get copy of local integration point for the slave integration cell
                Matrix< DDRMat > tSlaveLocalIntegPoint =
                        side_coordinate_map(
                                mSet->get_IG_geometry_type(),
                                tSlaveNodeOrdOnSide,
                                tMasterLocalIntegPoint );

                // set evaluation point for master and slave interpolators
                mSet->get_field_interpolator_manager( mtk::Master_Slave::MASTER )->    //
                        set_space_time_from_local_IG_point( tMasterLocalIntegPoint );

                mSet->get_field_interpolator_manager( mtk::Master_Slave::SLAVE )->    //
                        set_space_time_from_local_IG_point( tSlaveLocalIntegPoint );

                // compute detJ of integration domain
                real tDetJ = mSet->get_field_interpolator_manager()->get_IG_geometry_interpolator()->det_J();

                // skip if detJ smaller than threshold
                if ( tDetJ < Geometry_Interpolator::sDetJInvJacLowerLimit )
                {
                    continue;
                }

                // compute integration point weight
                real tWStar = mSet->get_integration_weights()( iGP ) * tDetJ;

                // get the normal from mesh
                Matrix< DDRMat > tNormal = mCluster->get_side_normal( mMasterCell, tMasterSideOrd );

                // loop over the IWGs
                for ( uint iIWG = 0; iIWG < tNumIWGs; iIWG++ )
                {
                    // get requested IWG
                    const std::shared_ptr< IWG >& tReqIWG =
                            mSet->get_requested_IWGs()( iIWG );

                    // reset IWG
                    tReqIWG->reset_eval_flags();

                    // set the normal for the IWG
                    tReqIWG->set_normal( tNormal );

                    // compute residual at integration point
                    tReqIWG->compute_residual( tWStar );

                    // compute Jacobian at integration point
                    // compute off-diagonal Jacobian for staggered solve
                    ( this->*m_compute_jacobian )( tReqIWG, tWStar );
                }
            }
        }

        //------------------------------------------------------------------------------

        void
        Element_Double_Sideset::compute_jacobian()
        {
            // get number of IWGs
            uint tNumIWGs = mSet->get_number_of_requested_IWGs();

            // check for active IWGs
            if ( tNumIWGs == 0 )
            {
                return;
            }

            // get treated side ordinal on the master and on the slave
            uint tMasterSideOrd = mCluster->mMasterListOfSideOrdinals( mCellIndexInCluster );
            uint tSlaveSideOrd  = mCluster->mSlaveListOfSideOrdinals( mCellIndexInCluster );

            // set the master/slave ig geometry interpolator physical/parametric space and time coefficients
            this->init_ig_geometry_interpolator( tMasterSideOrd, tSlaveSideOrd );

            // get first corresponding node from master to slave
            moris::mtk::Vertex const * tSlaveNode =
                    mCluster->get_left_vertex_pair(
                            mMasterCell->get_vertices_on_side_ordinal( tMasterSideOrd )( 0 ) );

            moris_index tSlaveNodeOrdOnSide =
                    mCluster->get_right_vertex_ordinal_on_facet( mCellIndexInCluster, tSlaveNode );

            // loop over the integration points
            uint tNumIntegPoints = mSet->get_number_of_integration_points();

            for ( uint iGP = 0; iGP < tNumIntegPoints; iGP++ )
            {
                // get local integration point for the master integration cell
                const Matrix< DDRMat >& tMasterLocalIntegPoint = mSet->get_integration_points().get_column( iGP );

                // get copy of local integration point for the slave integration cell
                Matrix< DDRMat > tSlaveLocalIntegPoint =
                        side_coordinate_map(
                                mSet->get_IG_geometry_type(),
                                tSlaveNodeOrdOnSide,
                                tMasterLocalIntegPoint );

                // set evaluation point for master and slave interpolators
                mSet->get_field_interpolator_manager( mtk::Master_Slave::MASTER )->    //
                        set_space_time_from_local_IG_point( tMasterLocalIntegPoint );

                mSet->get_field_interpolator_manager( mtk::Master_Slave::SLAVE )->    //
                        set_space_time_from_local_IG_point( tSlaveLocalIntegPoint );

                // compute detJ of integration domain
                real tDetJ = mSet->get_field_interpolator_manager()->get_IG_geometry_interpolator()->det_J();

                // skip if detJ smaller than threshold
                if ( tDetJ < Geometry_Interpolator::sDetJInvJacLowerLimit )
                {
                    continue;
                }

                // compute integration point weight
                real tWStar = mSet->get_integration_weights()( iGP ) * tDetJ;

                // get the normal from mesh and set if for the IWG
                Matrix< DDRMat > tNormal = mCluster->get_side_normal( mMasterCell, tMasterSideOrd );

                // loop over the IWGs
                for ( uint iIWG = 0; iIWG < tNumIWGs; iIWG++ )
                {
                    // get requested IWG
                    const std::shared_ptr< IWG >& tReqIWG =
                            mSet->get_requested_IWGs()( iIWG );

                    // reset IWG
                    tReqIWG->reset_eval_flags();

                    // set the normal for the IWG
                    tReqIWG->set_normal( tNormal );

                    // compute residual at integration point
                    ( this->*m_compute_jacobian )( tReqIWG, tWStar );
                }
            }
        }

        //------------------------------------------------------------------------------

        void
        Element_Double_Sideset::compute_jacobian_and_residual()
        {
            // get number of IWGs
            uint tNumIWGs = mSet->get_number_of_requested_IWGs();

            // check for active IWGs
            if ( tNumIWGs == 0 )
            {
                return;
            }

            // get treated side ordinal on the master and on the slave
            uint tMasterSideOrd = mCluster->mMasterListOfSideOrdinals( mCellIndexInCluster );
            uint tSlaveSideOrd  = mCluster->mSlaveListOfSideOrdinals( mCellIndexInCluster );

            // set the master/slave ig geometry interpolator physical/parametric space and time coefficients
            this->init_ig_geometry_interpolator( tMasterSideOrd, tSlaveSideOrd );

            // get first corresponding node from master to slave
            moris::mtk::Vertex const * tSlaveNode =
                    mCluster->get_left_vertex_pair( mMasterCell->get_vertices_on_side_ordinal( tMasterSideOrd )( 0 ) );
            moris_index tSlaveNodeOrdOnSide =
                    mCluster->get_right_vertex_ordinal_on_facet( mCellIndexInCluster, tSlaveNode );

            // loop over the integration points
            uint tNumIntegPoints = mSet->get_number_of_integration_points();

            for ( uint iGP = 0; iGP < tNumIntegPoints; iGP++ )
            {
                // get local integration point for the master integration cell
                const Matrix< DDRMat >& tMasterLocalIntegPoint = mSet->get_integration_points().get_column( iGP );

                // get copy of local integration point for the slave integration cell
                Matrix< DDRMat > tSlaveLocalIntegPoint =
                        side_coordinate_map(
                                mSet->get_IG_geometry_type(),
                                tSlaveNodeOrdOnSide,
                                tMasterLocalIntegPoint );

                // set evaluation point for master and slave interpolators
                mSet->get_field_interpolator_manager( mtk::Master_Slave::MASTER )->    //
                        set_space_time_from_local_IG_point( tMasterLocalIntegPoint );

                mSet->get_field_interpolator_manager( mtk::Master_Slave::SLAVE )->    //
                        set_space_time_from_local_IG_point( tSlaveLocalIntegPoint );

                // compute detJ of integration domain
                real tDetJ = mSet->get_field_interpolator_manager()->get_IG_geometry_interpolator()->det_J();

                // skip if detJ smaller than threshold
                if ( tDetJ < Geometry_Interpolator::sDetJInvJacLowerLimit )
                {
                    continue;
                }

                // compute integration point weight
                real tWStar = mSet->get_integration_weights()( iGP ) * tDetJ;

                // get the normal from mesh
                Matrix< DDRMat > tNormal = mCluster->get_side_normal( mMasterCell, tMasterSideOrd );

                // loop over the IWGs
                for ( uint iIWG = 0; iIWG < tNumIWGs; iIWG++ )
                {
                    // get requested IWG
                    const std::shared_ptr< IWG >& tReqIWG =
                            mSet->get_requested_IWGs()( iIWG );

                    // reset IWG
                    tReqIWG->reset_eval_flags();

                    // set the normal for the IWG
                    tReqIWG->set_normal( tNormal );

                    if ( mSet->mEquationModel->get_is_forward_analysis() )
                    {
                        // compute residual at integration point
                        tReqIWG->compute_residual( tWStar );
                    }

                    // compute jacobian at integration point
                    ( this->*m_compute_jacobian )( tReqIWG, tWStar );
                }

                mSet->mFemModel->mDoubleSidedSideSetsGaussPoints++;
            }
        }

        //------------------------------------------------------------------------------

        void
        Element_Double_Sideset::compute_dRdp()
        {
            // get number of IWGs
            uint tNumIWGs = mSet->get_number_of_requested_IWGs();

            // check for active IWGs
            if ( tNumIWGs == 0 )
            {
                return;
            }

            // get treated side ordinal on the master and on the slave
            uint tMasterSideOrd = mCluster->mMasterListOfSideOrdinals( mCellIndexInCluster );
            uint tSlaveSideOrd  = mCluster->mSlaveListOfSideOrdinals( mCellIndexInCluster );

            // set the master/slave ig geometry interpolator physical/parametric space and time coefficients
            Matrix< DDSMat > tGeoLocalAssembly;
            this->init_ig_geometry_interpolator(
                    tMasterSideOrd,
                    tSlaveSideOrd,
                    tGeoLocalAssembly );

            // get first corresponding node from master to slave
            moris::mtk::Vertex const * tSlaveNode =
                    mCluster->get_left_vertex_pair(
                            mMasterCell->get_vertices_on_side_ordinal( tMasterSideOrd )( 0 ) );

            moris_index tSlaveNodeOrdOnSide =
                    mCluster->get_right_vertex_ordinal_on_facet( mCellIndexInCluster, tSlaveNode );

            // get the vertices indices
            moris::Cell< Matrix< IndexMat > > tVertexIndices( 2 );
            tVertexIndices( 0 ) =
                    mMasterCell->get_vertices_ind_on_side_ordinal( tMasterSideOrd );
            tVertexIndices( 1 ) =
                    mSlaveCell->get_vertices_ind_on_side_ordinal( tSlaveSideOrd );

            // loop over integration points
            uint tNumIntegPoints = mSet->get_number_of_integration_points();

            for ( uint iGP = 0; iGP < tNumIntegPoints; iGP++ )
            {
                // get local integration point for the master integration cell
                const Matrix< DDRMat > tMasterLocalIntegPoint = mSet->get_integration_points().get_column( iGP );

                // get local integration point for the slave integration cell
                const Matrix< DDRMat > tSlaveLocalIntegPoint =
                        side_coordinate_map(
                                mSet->get_IG_geometry_type(),
                                tSlaveNodeOrdOnSide,
                                tMasterLocalIntegPoint );

                // set evaluation point for master and slave interpolators
                mSet->get_field_interpolator_manager( mtk::Master_Slave::MASTER )->    //
                        set_space_time_from_local_IG_point( tMasterLocalIntegPoint );

                mSet->get_field_interpolator_manager( mtk::Master_Slave::SLAVE )->    //
                        set_space_time_from_local_IG_point( tSlaveLocalIntegPoint );

                // compute detJ of integration domain
                real tDetJ = mSet->get_field_interpolator_manager()->get_IG_geometry_interpolator()->det_J();

                // skip if detJ smaller than threshold
                if ( tDetJ < Geometry_Interpolator::sDetJInvJacLowerLimit )
                {
                    continue;
                }

                // compute integration point weight
                real tWStar = mSet->get_integration_weights()( iGP ) * tDetJ;

                // get the normal from mesh and set if for the IWG
                Matrix< DDRMat > tNormal = mCluster->get_side_normal( mMasterCell, tMasterSideOrd );

                // loop over the IWGs
                for ( uint iIWG = 0; iIWG < tNumIWGs; iIWG++ )
                {
                    // get requested IWG
                    const std::shared_ptr< IWG >& tReqIWG =
                            mSet->get_requested_IWGs()( iIWG );

                    // reset IWG
                    tReqIWG->reset_eval_flags();

                    // set the normal for the IWG
                    tReqIWG->set_normal( tNormal );

                    // compute dRdpMat at evaluation point
                    ( this->*m_compute_dRdp )( tReqIWG, tWStar, tGeoLocalAssembly, tVertexIndices );
                }
            }
        }

        //------------------------------------------------------------------------------

        void
        Element_Double_Sideset::compute_dRdp_and_dQIdp()
        {
            // get number of IWGs
            uint tNumIWGs = mSet->get_number_of_requested_IWGs();

            // check for active IWGs
            if ( tNumIWGs == 0 )
            {
                return;
            }

            // get treated side ordinal on the master and on the slave
            uint tMasterSideOrd = mCluster->mMasterListOfSideOrdinals( mCellIndexInCluster );
            uint tSlaveSideOrd  = mCluster->mSlaveListOfSideOrdinals( mCellIndexInCluster );

            // set the master/slave ig geometry interpolator physical/parametric space and time coefficients
            Matrix< DDSMat > tGeoLocalAssembly;
            this->init_ig_geometry_interpolator(
                    tMasterSideOrd,
                    tSlaveSideOrd,
                    tGeoLocalAssembly );

            // get first corresponding node from master to slave
            moris::mtk::Vertex const * tSlaveNode =
                    mCluster->get_left_vertex_pair(
                            mMasterCell->get_vertices_on_side_ordinal( tMasterSideOrd )( 0 ) );

            moris_index tSlaveNodeOrdOnSide =
                    mCluster->get_right_vertex_ordinal_on_facet( mCellIndexInCluster, tSlaveNode );

            // get the vertices indices
            moris::Cell< Matrix< IndexMat > > tVertexIndices( 2 );
            tVertexIndices( 0 ) =
                    mMasterCell->get_vertices_ind_on_side_ordinal( tMasterSideOrd );
            tVertexIndices( 1 ) =
                    mSlaveCell->get_vertices_ind_on_side_ordinal( tSlaveSideOrd );

            // loop over integration points
            uint tNumIntegPoints = mSet->get_number_of_integration_points();

            for ( uint iGP = 0; iGP < tNumIntegPoints; iGP++ )
            {
                // get local integration point for the master integration cell
                const Matrix< DDRMat > tMasterLocalIntegPoint = mSet->get_integration_points().get_column( iGP );

                // get local integration point for the slave integration cell
                const Matrix< DDRMat > tSlaveLocalIntegPoint =
                        side_coordinate_map(
                                mSet->get_IG_geometry_type(),
                                tSlaveNodeOrdOnSide,
                                tMasterLocalIntegPoint );

                // set evaluation point for master and slave interpolators
                mSet->get_field_interpolator_manager( mtk::Master_Slave::MASTER )->    //
                        set_space_time_from_local_IG_point( tMasterLocalIntegPoint );

                mSet->get_field_interpolator_manager( mtk::Master_Slave::SLAVE )->    //
                        set_space_time_from_local_IG_point( tSlaveLocalIntegPoint );

                // compute detJ of integration domain
                const real tDetJ = mSet->get_field_interpolator_manager()->get_IG_geometry_interpolator()->det_J();

                // skip if detJ smaller than threshold
                if ( tDetJ < Geometry_Interpolator::sDetJInvJacLowerLimit )
                {
                    continue;
                }

                // compute integration point weight
                const real tWStar = mSet->get_integration_weights()( iGP ) * tDetJ;

                // get the normal from mesh and set if for the IWG
                const Matrix< DDRMat > tNormal = mCluster->get_side_normal( mMasterCell, tMasterSideOrd );

                // loop over the IWGs
                for ( uint iIWG = 0; iIWG < tNumIWGs; iIWG++ )
                {
                    // get requested IWG
                    const std::shared_ptr< IWG >& tReqIWG =
                            mSet->get_requested_IWGs()( iIWG );

                    // reset IWG
                    tReqIWG->reset_eval_flags();

                    // set the normal for the IWG
                    tReqIWG->set_normal( tNormal );

                    // compute dRdpMat at evaluation point
                    ( this->*m_compute_dRdp )( tReqIWG, tWStar, tGeoLocalAssembly, tVertexIndices );
                }

                // FIXME add part over IQIs
            }
        }

        //------------------------------------------------------------------------------

        real
        Element_Double_Sideset::compute_volume( mtk::Master_Slave aIsMaster )
        {
            // get treated side ordinal on the master and on the slave
            uint tMasterSideOrd = mCluster->mMasterListOfSideOrdinals( mCellIndexInCluster );
            uint tSlaveSideOrd  = mCluster->mSlaveListOfSideOrdinals( mCellIndexInCluster );

            // set the master/slave ig geometry interpolator physical/parametric space and time coefficients
            this->init_ig_geometry_interpolator( tMasterSideOrd, tSlaveSideOrd );

            // get number of integration points
            uint tNumOfIntegPoints = mSet->get_number_of_integration_points();

            // initialize volume
            real tVolume = 0;

            // get geometry interpolator
            Geometry_Interpolator* tIGGI =
                    mSet->get_field_interpolator_manager( aIsMaster )->get_IG_geometry_interpolator();

            // loop over integration points
            for ( uint iGP = 0; iGP < tNumOfIntegPoints; iGP++ )
            {
                // set integration point for geometry interpolator
                tIGGI->set_space_time( mSet->get_integration_points().get_column( iGP ) );

                // compute and add integration point contribution to volume
                tVolume += tIGGI->det_J() * mSet->get_integration_weights()( iGP );
            }

            // get time step
            real tTimeStep = tIGGI->get_time_step();

            // return the volume value
            return tVolume / tTimeStep;
        }

        //------------------------------------------------------------------------------

        void
        Element_Double_Sideset::compute_QI()
        {
            // get number of IQIs
            uint tNumIQIs = mSet->get_number_of_requested_IQIs();

            // check for active IQIs
            if ( tNumIQIs == 0 )
            {
                return;
            }

            // get treated side ordinal on the master and on the slave
            uint tMasterSideOrd = mCluster->mMasterListOfSideOrdinals( mCellIndexInCluster );
            uint tSlaveSideOrd  = mCluster->mSlaveListOfSideOrdinals( mCellIndexInCluster );

            // set the master/slave ig geometry interpolator physical/parametric space and time coefficients
            this->init_ig_geometry_interpolator( tMasterSideOrd, tSlaveSideOrd );

            // get first corresponding node from master to slave
            moris::mtk::Vertex const * tSlaveNode =
                    mCluster->get_left_vertex_pair( mMasterCell->get_vertices_on_side_ordinal( tMasterSideOrd )( 0 ) );
            moris_index tSlaveNodeOrdOnSide =
                    mCluster->get_right_vertex_ordinal_on_facet( mCellIndexInCluster, tSlaveNode );

            // loop over the integration points
            uint tNumIntegPoints = mSet->get_number_of_integration_points();

            for ( uint iGP = 0; iGP < tNumIntegPoints; iGP++ )
            {

                // get local integration point for the master integration cell
                const Matrix< DDRMat >& tMasterLocalIntegPoint = mSet->get_integration_points().get_column( iGP );

                // get copy of local integration point for the slave integration cell
                const Matrix< DDRMat > tSlaveLocalIntegPoint = side_coordinate_map(
                        mSet->get_IG_geometry_type(),
                        tSlaveNodeOrdOnSide,
                        tMasterLocalIntegPoint );

                // set evaluation point for master and slave interpolators
                mSet->get_field_interpolator_manager( mtk::Master_Slave::MASTER )->    //
                        set_space_time_from_local_IG_point( tMasterLocalIntegPoint );

                mSet->get_field_interpolator_manager( mtk::Master_Slave::SLAVE )->    //
                        set_space_time_from_local_IG_point( tSlaveLocalIntegPoint );

                // compute detJ of integration domain
                const real tDetJ = mSet->get_field_interpolator_manager()->get_IG_geometry_interpolator()->det_J();

                // skip if detJ smaller than threshold
                if ( tDetJ < Geometry_Interpolator::sDetJInvJacLowerLimit )
                {
                    continue;
                }

                // compute integration point weight
                const real tWStar = mSet->get_integration_weights()( iGP ) * tDetJ;

                // get the normal from mesh
                const Matrix< DDRMat > tNormal = mCluster->get_side_normal( mMasterCell, tMasterSideOrd );

                // loop over the IQIs
                for ( uint iIQI = 0; iIQI < tNumIQIs; iIQI++ )
                {
                    // get requested IQI
                    const std::shared_ptr< IQI >& tReqIQI =
                            mSet->get_requested_IQIs()( iIQI );

                    // reset IQI
                    tReqIQI->reset_eval_flags();

                    // set the normal for the IQI
                    tReqIQI->set_normal( tNormal );

                    // compute QI at evaluation point
                    tReqIQI->compute_QI( tWStar );
                }
            }
        }

        //------------------------------------------------------------------------------

    } /* namespace fem */
} /* namespace moris */
