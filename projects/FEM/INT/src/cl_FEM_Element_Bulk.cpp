#include <iostream>
//FEM/INT/src
#include "cl_FEM_Element_Bulk.hpp"
#include "cl_FEM_Field_Interpolator_Manager.hpp"
#include "cl_FEM_Set.hpp"
//FEM/MSI/src
#include "cl_MSI_Design_Variable_Interface.hpp"
#include "cl_MSI_Equation_Model.hpp"

namespace moris
{
    namespace fem
    {
        //----------------------------------------------------------------------

        Element_Bulk::Element_Bulk(
                mtk::Cell const    * aCell,
                Set                * aSet,
                Cluster            * aCluster,
                moris::moris_index   aCellIndexInCluster )
        : Element( aCell, aSet, aCluster, aCellIndexInCluster )
        {}

        //----------------------------------------------------------------------

        Element_Bulk::~Element_Bulk(){}

        //----------------------------------------------------------------------

        void Element_Bulk::init_ig_geometry_interpolator(
                Matrix< DDSMat > & aGeoLocalAssembly )
        {
            // get geometry interpolator for IG element
            Geometry_Interpolator * tIGGI =
                    mSet->get_field_interpolator_manager()->get_IG_geometry_interpolator();

            // get master physical space and time coordinates for IG element
            Matrix< DDRMat > tIGPhysSpaceCoords =
                    mMasterCell->get_vertex_coords();
            Matrix< DDRMat > tIGPhysTimeCoords =
                    mCluster->mInterpolationElement->get_time();

            // get master parametric space and time coordinates for IG element
            Matrix< DDRMat > tIGParamSpaceCoords =
                    mCluster->get_primary_cell_local_coords_on_side_wrt_interp_cell( mCellIndexInCluster );
            // FIXME not true if time is not linear
            Matrix< DDRMat > tIGParamTimeCoords = { { -1.0 }, { 1.0 } };

            // get the local cluster assembly indices
            if( mSet->get_geo_pdv_assembly_flag() )
            {
                // get the vertices indices for IG element
                Matrix< IndexMat > tVertexIndices = mMasterCell->get_vertex_inds();

                // get the requested geo pdv types
                moris::Cell < enum PDV_Type > tGeoPdvType;
                mSet->get_ig_unique_dv_types_for_set( tGeoPdvType );

                // get local assembly indices
                mSet->get_equation_model()->get_integration_xyz_pdv_assembly_indices(
                        tVertexIndices,
                        tGeoPdvType,
                        aGeoLocalAssembly );
            }

            // set physical space and time coefficients for IG element GI
            tIGGI->set_space_coeff( tIGPhysSpaceCoords );
            tIGGI->set_time_coeff(  tIGPhysTimeCoords );

            // set parametric space and time coefficients for IG element GI
            tIGGI->set_space_param_coeff( tIGParamSpaceCoords );
            tIGGI->set_time_param_coeff(  tIGParamTimeCoords );
        }

        //----------------------------------------------------------------------

        void Element_Bulk::init_ig_geometry_interpolator()
        {
            // get geometry interpolator for IG element
            Geometry_Interpolator * tIGGI =
                    mSet->get_field_interpolator_manager()->get_IG_geometry_interpolator();

            // get master physical space and time coordinates for IG element
            Matrix< DDRMat > tIGPhysSpaceCoords =
                    mMasterCell->get_vertex_coords();
            Matrix< DDRMat > tIGPhysTimeCoords =
                    mCluster->mInterpolationElement->get_time();

            // get master parametric space and time coordinates for IG element
            Matrix< DDRMat > tIGParamSpaceCoords =
                    mCluster->get_primary_cell_local_coords_on_side_wrt_interp_cell( mCellIndexInCluster );
            // FIXME not true if time is not linear
            Matrix< DDRMat > tIGParamTimeCoords = { { -1.0 }, { 1.0 } };

            // set physical space and time coefficients for IG element GI
            tIGGI->set_space_coeff( tIGPhysSpaceCoords );
            tIGGI->set_time_coeff(  tIGPhysTimeCoords );

            // set parametric space and time coefficients for IG element GI
            tIGGI->set_space_param_coeff( tIGParamSpaceCoords );
            tIGGI->set_time_param_coeff(  tIGParamTimeCoords );
        }

        //------------------------------------------------------------------------------

        void Element_Bulk::compute_residual()
        {
            // set physical and parametric space and time coefficients for IG element
            this->init_ig_geometry_interpolator();

            // get number of IWGs
            uint tNumIWGs = mSet->get_number_of_requested_IWGs();

            // loop over integration points
            uint tNumIntegPoints = mSet->get_number_of_integration_points();
            for( uint iGP = 0; iGP < tNumIntegPoints; iGP++ )
            {
                // get the current integration point in the IG param space
                const Matrix< DDRMat > & tLocalIntegPoint =
                        mSet->get_integration_points().get_column( iGP );

                // set evaluation point for interpolators (FIs and GIs)
                mSet->get_field_interpolator_manager()->set_space_time_from_local_IG_point( tLocalIntegPoint );

                // compute detJ of integration domain
                real tDetJ = mSet->get_field_interpolator_manager()->get_IG_geometry_interpolator()->det_J();

                // skip if detJ smaller than threshold
                if ( tDetJ < Geometry_Interpolator::sDetJInvJacLowerLimit )
                {
//                    MORIS_LOG_INFO("Element_Bulk::compute_residual: Skip quadrature point evaluation due to small detJ (%e)\n.",tDetJ);
                    continue;
                }

                // compute integration point weight
                real tWStar = mSet->get_integration_weights()( iGP ) * tDetJ;

                // loop over the IWGs
                for( uint iIWG = 0; iIWG < tNumIWGs; iIWG++ )
                {
                    // get requested IWG
                    const std::shared_ptr< IWG > & tReqIWG = mSet->get_requested_IWGs()( iIWG );

                    // reset IWG
                    tReqIWG->reset_eval_flags();

                    // FIXME: enforced nodal weak bcs
                    tReqIWG->set_nodal_weak_bcs(
                            mCluster->mInterpolationElement->get_weak_bcs() );

                    // compute residual at evaluation point
                    tReqIWG->compute_residual( tWStar );

                    // compute off-diagonal Jacobian for staggered solve
                    tReqIWG->compute_jacobian( tWStar );
                }
            }
        }

        //------------------------------------------------------------------------------

        void Element_Bulk::compute_jacobian()
        {
            // set physical and parametric space and time coefficients for IG element
            this->init_ig_geometry_interpolator();

            // get number of IWGs
            uint tNumIWGs = mSet->get_number_of_requested_IWGs();

            // loop over integration points
            uint tNumIntegPoints = mSet->get_number_of_integration_points();

            for( uint iGP = 0; iGP < tNumIntegPoints; iGP++ )
            {
                // get the ith integration point in the IG param space
                const Matrix< DDRMat > & tLocalIntegPoint =
                        mSet->get_integration_points().get_column( iGP );

                // set evaluation point for interpolators (FIs and GIs)
                mSet->get_field_interpolator_manager()->set_space_time_from_local_IG_point( tLocalIntegPoint );

                // compute detJ of integration domain
                real tDetJ = mSet->get_field_interpolator_manager()->get_IG_geometry_interpolator()->det_J();

                // skip if detJ smaller than threshold
                if ( tDetJ < Geometry_Interpolator::sDetJInvJacLowerLimit )
                {
//                    MORIS_LOG_INFO("Element_Bulk::compute_jacobian: Skip quadrature point evaluation due to small detJ (%e)\n.",tDetJ);
                    continue;
                }

                // compute integration point weight
                real tWStar = mSet->get_integration_weights()( iGP ) * tDetJ;

                // loop over the IWGs
                for( uint iIWG = 0; iIWG < tNumIWGs; iIWG++ )
                {
                    // get requested IWG
                    const std::shared_ptr< IWG > & tReqIWG = mSet->get_requested_IWGs()( iIWG );

                    // reset IWG
                    tReqIWG->reset_eval_flags();

                    // FIXME set nodal weak BCs
                    tReqIWG->set_nodal_weak_bcs(
                            mCluster->mInterpolationElement->get_weak_bcs() );

                    // compute Jacobian at evaluation point
                    tReqIWG->compute_jacobian( tWStar );
                }
            }
        }

        //------------------------------------------------------------------------------

        void Element_Bulk::compute_jacobian_and_residual()
        {
            // set physical and parametric space and time coefficients for IG element
            this->init_ig_geometry_interpolator();

            // get number of IWGs
            uint tNumIWGs = mSet->get_number_of_requested_IWGs();

            // get number of IQIs
            uint tNumIQIs = mSet->get_number_of_requested_IQIs();

            // loop over integration points
            uint tNumIntegPoints = mSet->get_number_of_integration_points();

            for( uint iGP = 0; iGP < tNumIntegPoints; iGP++ )
            {
                // get the ith integration point in the IG param space
                const Matrix< DDRMat > & tLocalIntegPoint =
                        mSet->get_integration_points().get_column( iGP );

                // set evaluation point for interpolators (FIs and GIs)
                mSet->get_field_interpolator_manager()->set_space_time_from_local_IG_point( tLocalIntegPoint );

                // compute detJ of integration domain
                real tDetJ = mSet->get_field_interpolator_manager()->get_IG_geometry_interpolator()->det_J();

                // skip if detJ smaller than threshold
                if ( tDetJ < Geometry_Interpolator::sDetJInvJacLowerLimit )
                {
//                    MORIS_LOG_INFO("Element_Bulk::compute_jacobian_and_residual: Skip quadrature point evaluation due to small detJ (%e)\n.",tDetJ);
                    continue;
                }

                // compute integration point weight
                real tWStar = mSet->get_integration_weights()( iGP ) * tDetJ;

                // loop over the IWGs
                for( uint iIWG = 0; iIWG < tNumIWGs; iIWG++ )
                {
                    // get requested IWG
                    const std::shared_ptr< IWG > & tReqIWG = mSet->get_requested_IWGs()( iIWG );

                    // reset IWG
                    tReqIWG->reset_eval_flags();

                    // FIXME set nodal weak BCs
                    tReqIWG->set_nodal_weak_bcs(
                            mCluster->mInterpolationElement->get_weak_bcs() );

                    if( mSet->mEquationModel->get_is_forward_analysis() )
                    {
                        // compute residual at evaluation point
                        tReqIWG->compute_residual( tWStar );
                    }

                    // compute Jacobian at evaluation point
                    tReqIWG->compute_jacobian( tWStar );
                }

                if( ( !mSet->mEquationModel->get_is_forward_analysis() ) && ( tNumIQIs > 0 ) )
                {
                    // loop over the IQIs
                    for( uint iIQI = 0; iIQI < tNumIQIs; iIQI++ )
                    {
                        // reset IWG
                        mSet->get_requested_IQIs()( iIQI )->reset_eval_flags();

                        // compute dQIdu at evaluation point
                        mSet->get_requested_IQIs()( iIQI )->add_dQIdu_on_set( tWStar );
                    }
                }
            }
        }

        //------------------------------------------------------------------------------

        void Element_Bulk::compute_dRdp()
        {
            // set physical and parametric space and time coefficients for IG element
            Matrix< DDSMat > tGeoLocalAssembly;
            this->init_ig_geometry_interpolator( tGeoLocalAssembly );

            // get number of IWGs
            uint tNumIWGs = mSet->get_number_of_requested_IWGs();

            // loop over integration points
            uint tNumIntegPoints = mSet->get_number_of_integration_points();
            for( uint iGP = 0; iGP < tNumIntegPoints; iGP++ )
            {
                // get the ith integration point in the IG param space
                const Matrix< DDRMat > & tLocalIntegPoint =
                        mSet->get_integration_points().get_column( iGP );

                // set evaluation point for interpolators (FIs and GIs)
                mSet->get_field_interpolator_manager()->set_space_time_from_local_IG_point( tLocalIntegPoint );

                // compute detJ of integration domain
                real tDetJ = mSet->get_field_interpolator_manager()->get_IG_geometry_interpolator()->det_J();

                // skip if detJ smaller than threshold
                if ( tDetJ < Geometry_Interpolator::sDetJInvJacLowerLimit )
                {
//                    MORIS_LOG_INFO("Element_Bulk::compute_dRdp: Skip quadrature point evaluation due to small detJ (%e)\n.",tDetJ);
                    continue;
                }

                // compute integration point weight
                real tWStar = mSet->get_integration_weights()( iGP ) * tDetJ;

                // loop over the IWGs
                for( uint iIWG = 0; iIWG < tNumIWGs; iIWG++ )
                {
                    // get requested IWG
                    const std::shared_ptr< IWG > & tReqIWG = mSet->get_requested_IWGs()( iIWG );

                    // reset IWG
                    tReqIWG->reset_eval_flags();

                    // FIXME set nodal weak BCs
                    tReqIWG->set_nodal_weak_bcs(
                            mCluster->mInterpolationElement->get_weak_bcs() );

                    // compute dRdpMat at evaluation point
                    tReqIWG->compute_dRdp( tWStar );
                }
            }
        }

        //------------------------------------------------------------------------------

        void Element_Bulk::compute_dRdp_FD()
        {
            // get finite difference scheme type
            fem::FDScheme_Type tFDScheme =
                    mSet->get_finite_difference_scheme_for_sensitivity_analysis();

            // get the finite difference perturbation size
            real tFDPerturbation = mSet->get_finite_difference_perturbation_size();

            // get the vertices indices
            Matrix< IndexMat > tVertexIndices = mMasterCell->get_vertex_inds();

            // set physical and parametric space and time coefficients for IG element
            Matrix< DDSMat > tGeoLocalAssembly;
            this->init_ig_geometry_interpolator( tGeoLocalAssembly );

            // get number of IWGs
            uint tNumIWGs = mSet->get_number_of_requested_IWGs();

            // loop over integration points
            uint tNumIntegPoints = mSet->get_number_of_integration_points();
            for( uint iGP = 0; iGP < tNumIntegPoints; iGP++ )
            {
                // get the ith integration point in the IG param space
                const Matrix< DDRMat > & tLocalIntegPoint =
                        mSet->get_integration_points().get_column( iGP );

                // set evaluation point for interpolators (FIs and GIs)
                mSet->get_field_interpolator_manager()->
                        set_space_time_from_local_IG_point( tLocalIntegPoint );

                // compute detJ of integration domain
                real tDetJ = mSet->get_field_interpolator_manager()->get_IG_geometry_interpolator()->det_J();

                // skip if detJ smaller than threshold
                if ( tDetJ < Geometry_Interpolator::sDetJInvJacLowerLimit )
                {
//                    MORIS_LOG_INFO("Element_Bulk::compute_dRdp_FD: Skip quadrature point evaluation due to small detJ (%e)\n.",tDetJ);
                    continue;
                }

                // compute integration point weight
                real tWStar = mSet->get_integration_weights()( iGP ) * tDetJ;

                // loop over the IWGs
                for( uint iIWG = 0; iIWG < tNumIWGs; iIWG++ )
                {
                    // get requested IWG
                    const std::shared_ptr< IWG > & tReqIWG = mSet->get_requested_IWGs()( iIWG );

                    // reset IWG
                    tReqIWG->reset_eval_flags();

                    // FIXME set nodal weak BCs
                    tReqIWG->set_nodal_weak_bcs(
                            mCluster->mInterpolationElement->get_weak_bcs() );

                    // compute dRdpMat at evaluation point
                    tReqIWG->compute_dRdp_FD_material(
                            tWStar,
                            tFDPerturbation,
                            tFDScheme );

                    // compute dRdpGeo at evaluation point
                    if( mSet->get_geo_pdv_assembly_flag() )
                    {
                        tReqIWG->compute_dRdp_FD_geometry(
                                tWStar,
                                tFDPerturbation,
                                tGeoLocalAssembly,
                                tFDScheme );
                    }
                }
            }
        }

        //------------------------------------------------------------------------------

        void Element_Bulk::compute_QI()
        {
            // set physical and parametric space and time coefficients for IG element
            this->init_ig_geometry_interpolator();

            // get number of IQIs
            uint tNumIQIs = mSet->get_number_of_requested_IQIs();

            // loop over integration points
            uint tNumIntegPoints = mSet->get_number_of_integration_points();
            for( uint iGP = 0; iGP < tNumIntegPoints; iGP++ )
            {
                // get the ith integration point in the IG param space
                const Matrix< DDRMat > & tLocalIntegPoint =
                        mSet->get_integration_points().get_column( iGP );

                // set evaluation point for interpolators (FIs and GIs)
                mSet->get_field_interpolator_manager()->set_space_time_from_local_IG_point( tLocalIntegPoint );

                // compute detJ of integration domain
                real tDetJ = mSet->get_field_interpolator_manager()->get_IG_geometry_interpolator()->det_J();

                // skip if detJ smaller than threshold
                if ( tDetJ < Geometry_Interpolator::sDetJInvJacLowerLimit )
                {
//                    MORIS_LOG_INFO("Element_Bulk::compute_QI: Skip quadrature point evaluation due to small detJ (%e)\n.",tDetJ);
                    continue;
                }

                // compute integration point weight
                real tWStar = mSet->get_integration_weights()( iGP ) * tDetJ;

                // loop over the IQIs
                for( uint iIQI = 0; iIQI < tNumIQIs; iIQI++ )
                {
                    // reset IQI
                    mSet->get_requested_IQIs()( iIQI )->reset_eval_flags();

                    // compute QI at evaluation point
                    mSet->get_requested_IQIs()( iIQI )->add_QI_on_set( tWStar );
                }
            }
        }

        //------------------------------------------------------------------------------

        void Element_Bulk::compute_dQIdp_explicit()
        {
            // set physical and parametric space and time coefficients for IG element
            Matrix< DDSMat > tGeoLocalAssembly;
            this->init_ig_geometry_interpolator( tGeoLocalAssembly );

            // get number of IWGs
            uint tNumIQIs = mSet->get_number_of_requested_IQIs();

            // loop over integration points
            uint tNumIntegPoints = mSet->get_number_of_integration_points();
            for( uint iGP = 0; iGP < tNumIntegPoints; iGP++ )
            {
                // get the ith integration point in the IG param space
                const Matrix< DDRMat > & tLocalIntegPoint =
                        mSet->get_integration_points().get_column( iGP );

                // set evaluation point for interpolators (FIs and GIs)
                mSet->get_field_interpolator_manager()->set_space_time_from_local_IG_point( tLocalIntegPoint );

                // compute detJ of integration domain
                real tDetJ = mSet->get_field_interpolator_manager()->get_IG_geometry_interpolator()->det_J();

                // skip if detJ smaller than threshold
                if ( tDetJ < Geometry_Interpolator::sDetJInvJacLowerLimit )
                {
//                    MORIS_LOG_INFO("Element_Bulk::compute_dQIdp_explicit: Skip quadrature point evaluation due to small detJ (%e)\n.",tDetJ);
                    continue;
                }

                // compute integration point weight
                real tWStar = mSet->get_integration_weights()( iGP ) * tDetJ;

                // loop over the IQIs
                for( uint iIQI = 0; iIQI < tNumIQIs; iIQI++ )
                {
                    // reset IWG
                    mSet->get_requested_IQIs()( iIQI )->reset_eval_flags();

                    // compute dQIdpMat at evaluation point
                    mSet->get_requested_IQIs()( iIQI )->compute_dQIdp( tWStar );
                }
            }
        }

        //------------------------------------------------------------------------------

        void Element_Bulk::compute_dQIdp_explicit_FD()
        {
            // get finite difference scheme type
            fem::FDScheme_Type tFDScheme =
                    mSet->get_finite_difference_scheme_for_sensitivity_analysis();

            // get the finite difference perturbation size
            real tFDPerturbation = mSet->get_finite_difference_perturbation_size();

            // get the vertices indices
            Matrix< IndexMat > tVertexIndices = mMasterCell->get_vertex_inds();

            // set physical and parametric space and time coefficients for IG element
            Matrix< DDSMat > tGeoLocalAssembly;
            this->init_ig_geometry_interpolator( tGeoLocalAssembly );

            // get number of IWGs
            uint tNumIQIs = mSet->get_number_of_requested_IQIs();

            // loop over integration points
            uint tNumIntegPoints = mSet->get_number_of_integration_points();
            for( uint iGP = 0; iGP < tNumIntegPoints; iGP++ )
            {
                // get the ith integration point in the IG param space
                const Matrix< DDRMat > & tLocalIntegPoint =
                        mSet->get_integration_points().get_column( iGP );

                // set evaluation point for interpolators (FIs and GIs)
                mSet->get_field_interpolator_manager()->set_space_time_from_local_IG_point( tLocalIntegPoint );

                // compute detJ of integration domain
                real tDetJ = mSet->get_field_interpolator_manager()->get_IG_geometry_interpolator()->det_J();

                // skip if detJ smaller than threshold
                if ( tDetJ < Geometry_Interpolator::sDetJInvJacLowerLimit )
                {
//                    MORIS_LOG_INFO("Element_Bulk::compute_dQIdp_explicit_FD: Skip quadrature point evaluation due to small detJ (%e)\n.",tDetJ);
                    continue;
                }

                // compute integration point weight
                real tWStar = mSet->get_integration_weights()( iGP ) * tDetJ;

                // loop over the IQIs
                for( uint iIQI = 0; iIQI < tNumIQIs; iIQI++ )
                {
                    // reset IWG
                    mSet->get_requested_IQIs()( iIQI )->reset_eval_flags();

                    // compute dQIdpMat at evaluation point
                    mSet->get_requested_IQIs()( iIQI )->compute_dQIdp_FD_material(
                            tWStar,
                            tFDPerturbation,
                            tFDScheme );

                    // compute dQIdpGeo at evaluation point
                    if( mSet->get_geo_pdv_assembly_flag() )
                    {
                        mSet->get_requested_IQIs()( iIQI )->compute_dQIdp_FD_geometry(
                                tWStar,
                                tFDPerturbation,
                                tGeoLocalAssembly,
                                tFDScheme );
                    }
                }
            }
        }

        //------------------------------------------------------------------------------

        void Element_Bulk::compute_dQIdu()
        {
            // set physical and parametric space and time coefficients for IG element
            this->init_ig_geometry_interpolator();

            // get number of IQIs
            uint tNumIQIs = mSet->get_number_of_requested_IQIs();

            // loop over integration points
            uint tNumIntegPoints = mSet->get_number_of_integration_points();
            for( uint iGP = 0; iGP < tNumIntegPoints; iGP++ )
            {
                // get the ith integration point in the IG param space
                const Matrix< DDRMat > & tLocalIntegPoint =
                        mSet->get_integration_points().get_column( iGP );

                // set evaluation point for interpolators (FIs and GIs)
                mSet->get_field_interpolator_manager()->set_space_time_from_local_IG_point( tLocalIntegPoint );

                // compute detJ of integration domain
                real tDetJ = mSet->get_field_interpolator_manager()->get_IG_geometry_interpolator()->det_J();

                // skip if detJ smaller than threshold
                if ( tDetJ < Geometry_Interpolator::sDetJInvJacLowerLimit )
                {
//                    MORIS_LOG_INFO("Element_Bulk::compute_dQIdu: Skip quadrature point evaluation due to small detJ (%e)\n.",tDetJ);
                    continue;
                }

                // compute integration point weight
                real tWStar = mSet->get_integration_weights()( iGP ) * tDetJ;

                // loop over the IQIs
                for( uint iIQI = 0; iIQI < tNumIQIs; iIQI++ )
                {
//                    // check if IQI has dof dependencies
//                    moris::Cell< moris::Cell< MSI::Dof_Type > > aDofTypeList =
//                            mSet->get_requested_IQIs()( iIQI )->get_global_dof_type_list();

                    // if there are dof dependencies
                    //if( aDofTypeList.size() > 0 )
                    //{
                        // reset IWG
                        mSet->get_requested_IQIs()( iIQI )->reset_eval_flags();

                        // compute dQIdu at evaluation point
                        mSet->get_requested_IQIs()( iIQI )->add_dQIdu_on_set( tWStar );
                    //}
                }
            }
        }

        //------------------------------------------------------------------------------

        void Element_Bulk::compute_quantity_of_interest_global(
                const uint                         aMeshIndex,
                const moris::Cell< std::string > & aQINames )
        {
            // set physical and parametric space and time coefficients for IG element
            this->init_ig_geometry_interpolator();

            // loop over integration points
            uint tNumIntegPoints = mSet->get_number_of_integration_points();
            for( uint iGP = 0; iGP < tNumIntegPoints; iGP++ )
            {
                // get the ith integration point in the IG param space
                const Matrix< DDRMat > & tLocalIntegPoint =
                        mSet->get_integration_points().get_column( iGP );

                // set evaluation point for interpolators (FIs and GIs)
                mSet->get_field_interpolator_manager()->set_space_time_from_local_IG_point( tLocalIntegPoint );

                // compute detJ of integration domain
                real tDetJ = mSet->get_field_interpolator_manager()->get_IG_geometry_interpolator()->det_J();

                // skip if detJ smaller than threshold
                if ( tDetJ < Geometry_Interpolator::sDetJInvJacLowerLimit )
                {
//                    MORIS_LOG_INFO("Element_Bulk::compute_quantity_of_interest_global: Skip quadrature point evaluation due to small detJ (%e)\n.",tDetJ);
                    continue;
                }

                // compute integration point weight
                real tWStar = mSet->get_integration_weights()( iGP ) * tDetJ;

                // loop over IQI
                uint tNumIQIs = aQINames.size();
                for( uint iIQI = 0; iIQI < tNumIQIs; iIQI++ )
                {
                    if( mSet->mIQINameToIndexMap.key_exists( aQINames( iIQI ) ) )
                    {
                        // get the set local index
                        moris_index tIQISetLocalIndex =
                                mSet->mIQINameToIndexMap.find( aQINames( iIQI ) );

                        // reset the requested IQI
                        mSet->mIQIs( tIQISetLocalIndex )->reset_eval_flags();

                        // compute quantity of interest at evaluation point
                        Matrix< DDRMat > tQIValue;
                        mSet->mIQIs( tIQISetLocalIndex )->compute_QI( tQIValue );

                        // assemble the global QI value on the set
                        ( *( mSet->mSetGlobalValues ) )( iIQI ) += tQIValue( 0 ) * tWStar;
                    }
                }
            }
        }

        //------------------------------------------------------------------------------

        void Element_Bulk::compute_quantity_of_interest_elemental(
                const uint                         aMeshIndex,
                const moris::Cell< std::string > & aQINames )
        {
            // set physical and parametric space and time coefficients for IG element
            this->init_ig_geometry_interpolator();

            // loop over integration points
            uint tNumIntegPoints = mSet->get_number_of_integration_points();
            for( uint iGP = 0; iGP < tNumIntegPoints; iGP++ )
            {
                // get the ith integration point in the IG param space
                const Matrix< DDRMat > & tLocalIntegPoint =
                        mSet->get_integration_points().get_column( iGP );

                // set evaluation point for interpolators (FIs and GIs)
                mSet->get_field_interpolator_manager()->set_space_time_from_local_IG_point( tLocalIntegPoint );

                // compute detJ of integration domain
                real tDetJ = mSet->get_field_interpolator_manager()->get_IG_geometry_interpolator()->det_J();

                // skip if detJ smaller than threshold
                if ( tDetJ < Geometry_Interpolator::sDetJInvJacLowerLimit )
                {
//                    MORIS_LOG_INFO("Element_Bulk::compute_quantity_of_interest_elemental: Skip quadrature point evaluation due to small detJ (%e)\n.",tDetJ);
                    continue;
                }

                // compute integration point weight
                real tWStar = mSet->get_integration_weights()( iGP ) * tDetJ;

                // loop over IQI
                for( uint iIQI = 0; iIQI < aQINames.size(); iIQI++ )
                {
                    // if IQI defined
                    if( mSet->mIQINameToIndexMap.key_exists( aQINames( iIQI ) ) )
                    {
                        // get the set local index
                        moris_index tIQISetLocalIndex =
                                mSet->mIQINameToIndexMap.find( aQINames( iIQI ) );

                        // reset the requested IQI
                        mSet->mIQIs( tIQISetLocalIndex )->reset_eval_flags();

                        // compute quantity of interest at evaluation point
                        Matrix< DDRMat > tQIValue;
                        mSet->mIQIs( tIQISetLocalIndex )->compute_QI( tQIValue );

                        // assemble the nodal QI value on the set
                        ( *mSet->mSetElementalValues )( mSet->mCellAssemblyMap( aMeshIndex )( mMasterCell->get_index() ), iIQI ) +=
                                tQIValue( 0 ) * tWStar / tNumIntegPoints;
                    }
                }
            }
        }

        //------------------------------------------------------------------------------

        real Element_Bulk::compute_volume( mtk::Master_Slave aIsMaster )
        {
            // set physical and parametric space and time coefficients for IG element
            this->init_ig_geometry_interpolator();

            //get number of integration points
            uint tNumOfIntegPoints = mSet->get_number_of_integration_points();

            // init volume
            real tVolume = 0;

            // get geometry interpolator
            Geometry_Interpolator * tIGGI =
                    mSet->get_field_interpolator_manager()->get_IG_geometry_interpolator();

            // loop over integration points
            for( uint iGP = 0; iGP < tNumOfIntegPoints; iGP++ )
            {
                // set integration point for geometry interpolator
                tIGGI->set_space_time( mSet->get_integration_points().get_column( iGP ) );

                // compute and add integration point contribution to volume
                tVolume += tIGGI->det_J() * mSet->get_integration_weights()( iGP );
            }

            // return the volume value
            return tVolume;
        }

        //------------------------------------------------------------------------------

    } /* namespace fem */
} /* namespace moris */
