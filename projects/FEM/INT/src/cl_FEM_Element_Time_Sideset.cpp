#include <iostream>
//FEM/INT/src
#include "cl_FEM_Element_Time_Sideset.hpp"
#include "cl_FEM_Field_Interpolator_Manager.hpp"
#include "cl_FEM_Set.hpp"
//FEM/MSI/src
#include "cl_MSI_Equation_Model.hpp"
#include "cl_MSI_Design_Variable_Interface.hpp"

namespace moris
{
    namespace fem
    {

        //------------------------------------------------------------------------------

        Element_Time_Sideset::Element_Time_Sideset(
                mtk::Cell const  * aCell,
                Set              * aElementBlock,
                Cluster          * aCluster,
                moris::moris_index aCellIndexInCluster )
        : Element( aCell, aElementBlock, aCluster, aCellIndexInCluster )
        {}

        //------------------------------------------------------------------------------

        Element_Time_Sideset::~Element_Time_Sideset(){}

        //------------------------------------------------------------------------------

        void Element_Time_Sideset::init_ig_geometry_interpolator(
                moris::Cell< Matrix< DDSMat > > & aIsActiveDv )
        {
            // get geometry interpolator
            Geometry_Interpolator * tIGGI =
                    mSet->get_field_interpolator_manager()->get_IG_geometry_interpolator();
            Geometry_Interpolator * tPreviousIGGI =
                    mSet->get_field_interpolator_manager_previous_time()->get_IG_geometry_interpolator();

            // get physical space and current and previous time coordinates for IG element
            Matrix< DDRMat > tIGPhysSpaceCoords = mMasterCell->get_vertex_coords();
            Matrix< DDRMat > tIGPhysTimeCoords( 1, 1,
                    mCluster->mInterpolationElement->get_time()( 0 ) );
            Matrix< DDRMat > tIGPhysPreviousTimeCoords( 1, 1,
                    mCluster->mInterpolationElement->get_previous_time()( 1 ) );

            // get master parametric space and current and previous time coordinates for IG element
            Matrix< DDRMat > tIGParamSpaceCoords =
                    mCluster->get_primary_cell_local_coords_on_side_wrt_interp_cell( mCellIndexInCluster );
            Matrix< DDRMat > tIGParamTimeCoords( 1, 1, -1.0 );
            Matrix< DDRMat > tIGParamPreviousTimeCoords( 1, 1, 1.0 );

            // get the requested geo pdv types
            moris::Cell < enum PDV_Type > tGeoPdvType;
            mSet->get_ig_unique_dv_types_for_set( tGeoPdvType );

            // determine if there are IG pdvs
            if ( tGeoPdvType.size() )
            {
                // get space dimension
                uint tSpaceDim = tIGPhysSpaceCoords.n_cols();

                // reshape the XYZ values into a cell of vectors
                moris::Cell< Matrix< DDRMat > > tPdvValueList( tSpaceDim );
                for( uint iSpaceDim = 0; iSpaceDim < tSpaceDim; iSpaceDim++ )
                {
                    tPdvValueList( iSpaceDim ) = tIGPhysSpaceCoords.get_column( iSpaceDim );
                }

                // get the vertices indices for IG element
                Matrix< IndexMat > tVertexIndices = mMasterCell->get_vertex_inds();

                // get the pdv values from the MSI/GEN interface
                mSet->get_equation_model()->get_design_variable_interface()->get_ig_pdv_value(
                        tVertexIndices,
                        tGeoPdvType,
                        tPdvValueList,
                        aIsActiveDv );

                // reshape the cell of vectors tPdvValueList into a matrix tIGPhysSpaceCoords
                tIGPhysSpaceCoords.set_size( 0, 0 );
                mSet->get_equation_model()->get_design_variable_interface()->reshape_pdv_values(
                        tPdvValueList,
                        tIGPhysSpaceCoords );
            }

            // set physical space and current time coefficients for IG element GI
            tIGGI->set_space_coeff( tIGPhysSpaceCoords );
            tIGGI->set_time_coeff(  tIGPhysTimeCoords );

            // set physical space and previous time coefficients for IG element GI
            tPreviousIGGI->set_space_coeff( tIGPhysSpaceCoords );
            tPreviousIGGI->set_time_coeff(  tIGPhysPreviousTimeCoords );

            // set parametric space and current time coefficients for IG element GI
            tIGGI->set_space_param_coeff( tIGParamSpaceCoords );
            tIGGI->set_time_param_coeff(  tIGParamTimeCoords );

            // set parametric space and previous time coefficients for IG element GI
            tPreviousIGGI->set_space_param_coeff( tIGParamSpaceCoords );
            tPreviousIGGI->set_time_param_coeff(  tIGParamPreviousTimeCoords );
        }

        //------------------------------------------------------------------------------

        void Element_Time_Sideset::compute_residual()
        {
            // set physical and parametric space and time coefficients for IG element
            moris::Cell< Matrix< DDSMat > > tIsActiveDv;
            this->init_ig_geometry_interpolator( tIsActiveDv );

            // get number of IWGs
            uint tNumIWGs = mSet->get_number_of_requested_IWGs();

            // loop over integration points
            uint tNumIntegPoints = mSet->get_number_of_integration_points();
            for( uint iGP = 0; iGP < tNumIntegPoints; iGP++ )
            {
                // get integration point location in the reference surface
                Matrix< DDRMat > tLocalIntegPoint = mSet->get_integration_points().get_column( iGP );

                // set evaluation point for interpolators (FIs and GIs)
                mSet->get_field_interpolator_manager()->
                        set_space_time_from_local_IG_point( tLocalIntegPoint );
                mSet->get_field_interpolator_manager_previous_time()->
                        set_space_time_from_local_IG_point( tLocalIntegPoint );

                // compute integration point weight
                real tWStar = mSet->get_integration_weights()( iGP ) *
                        mSet->get_field_interpolator_manager()->get_IG_geometry_interpolator()->det_J();

                // loop over the IWGs
                for( uint iIWG = 0; iIWG < tNumIWGs; iIWG++ )
                {
                    // reset IWG
                    mSet->get_requested_IWGs()( iIWG )->reset_eval_flags();

                    // compute jacobian at evaluation point
                    mSet->get_requested_IWGs()( iIWG )->compute_residual( tWStar );

                    // compute jacobian at evaluation point
                    // compute off-diagonal jacobian for staggered solve
                    mSet->get_requested_IWGs()( iIWG )->compute_jacobian( tWStar );
                }
            }
        }

        //------------------------------------------------------------------------------

        void Element_Time_Sideset::compute_jacobian()
        {
            // set physical and parametric space and time coefficients for IG element
            moris::Cell< Matrix< DDSMat > > tIsActiveDv;
            this->init_ig_geometry_interpolator( tIsActiveDv );

            // get number of IWGs
            uint tNumIWGs = mSet->get_number_of_requested_IWGs();

            // loop over integration points
            uint tNumIntegPoints = mSet->get_number_of_integration_points();
            for( uint iGP = 0; iGP < tNumIntegPoints; iGP++ )
            {
                // get local integration point location
                Matrix< DDRMat > tLocalIntegPoint = mSet->get_integration_points().get_column( iGP );

                // set evaluation point for interpolators (FIs and GIs)
                mSet->get_field_interpolator_manager()->
                        set_space_time_from_local_IG_point( tLocalIntegPoint );
                mSet->get_field_interpolator_manager_previous_time()->
                        set_space_time_from_local_IG_point( tLocalIntegPoint );

                // compute integration point weight
                real tWStar = mSet->get_integration_weights()( iGP ) *
                        mSet->get_field_interpolator_manager()->get_IG_geometry_interpolator()->det_J();

                // loop over the IWGs
                for( uint iIWG = 0; iIWG < tNumIWGs; iIWG++ )
                {
                    // reset IWG
                    mSet->get_requested_IWGs()( iIWG )->reset_eval_flags();

                    if( mSet->mEquationModel->get_is_adjoint_off_diagonal_time_contribution() )
                    {
                        mSet->get_requested_IWGs()( iIWG )->compute_jacobian_previous( tWStar );
                    }
                    else
                    {
                        // compute jacobian at evaluation point
                        mSet->get_requested_IWGs()( iIWG )->compute_jacobian( tWStar );
                    }
                }
            }
        }

        //------------------------------------------------------------------------------

        void Element_Time_Sideset::compute_jacobian_and_residual()
        {
            // set physical and parametric space and time coefficients for IG element
            moris::Cell< Matrix< DDSMat > > tIsActiveDv;
            this->init_ig_geometry_interpolator( tIsActiveDv );

            // get number of IWGs
            uint tNumIWGs = mSet->get_number_of_requested_IWGs();

            // loop over integration points
            uint tNumIntegPoints = mSet->get_number_of_integration_points();
            for( uint iGP = 0; iGP < tNumIntegPoints; iGP++ )
            {
                // get local integration point location
                Matrix< DDRMat > tLocalIntegPoint = mSet->get_integration_points().get_column( iGP );

                // set evaluation point for interpolators (FIs and GIs)
                mSet->get_field_interpolator_manager()->
                        set_space_time_from_local_IG_point( tLocalIntegPoint );
                mSet->get_field_interpolator_manager_previous_time()->
                        set_space_time_from_local_IG_point( tLocalIntegPoint );

                // compute integration point weight
                real tWStar = mSet->get_integration_weights()( iGP ) *
                        mSet->get_field_interpolator_manager()->get_IG_geometry_interpolator()->det_J();

                // loop over the IWGs
                for( uint iIWG = 0; iIWG < tNumIWGs; iIWG++ )
                {
                    // reset IWG
                    mSet->get_requested_IWGs()( iIWG )->reset_eval_flags();

                    // compute jacobian at evaluation point
                    mSet->get_requested_IWGs()( iIWG )->compute_residual( tWStar );

                    // compute jacobian at evaluation point
                    mSet->get_requested_IWGs()( iIWG )->compute_jacobian( tWStar );
                }
            }
        }

        //------------------------------------------------------------------------------

        void Element_Time_Sideset::compute_dRdp()
        {
            // get the vertices indices
            Matrix< IndexMat > tVertexIndices = mMasterCell->get_vertex_inds();

            // set physical and parametric space and time coefficients for IG element
            moris::Cell< Matrix< DDSMat > > tIsActiveDv;
            this->init_ig_geometry_interpolator( tIsActiveDv );

            // get number of IWGs
            uint tNumIWGs = mSet->get_number_of_requested_IWGs();

            // loop over integration points
            uint tNumIntegPoints = mSet->get_number_of_integration_points();
            for( uint iGP = 0; iGP < tNumIntegPoints; iGP++ )
            {
                // get integration point location in the reference surface
                Matrix< DDRMat > tLocalIntegPoint = mSet->get_integration_points().get_column( iGP );

                // set evaluation point for interpolators (FIs and GIs)
                mSet->get_field_interpolator_manager()->
                        set_space_time_from_local_IG_point( tLocalIntegPoint );
                mSet->get_field_interpolator_manager_previous_time()->
                        set_space_time_from_local_IG_point( tLocalIntegPoint );

                // compute integration point weight
                real tWStar = mSet->get_integration_weights()( iGP ) *
                        mSet->get_field_interpolator_manager()->get_IG_geometry_interpolator()->det_J();

                // loop over the IWGs
                for( uint iIWG = 0; iIWG < tNumIWGs; iIWG++ )
                {
                    // reset IWG
                    mSet->get_requested_IWGs()( iIWG )->reset_eval_flags();

                    // set a perturbation size
                    real tPerturbation = 1E-6;

                    // compute dRdpMat at evaluation point
                    mSet->get_requested_IWGs()( iIWG )->compute_dRdp_FD_material(
                            tWStar,
                            tPerturbation );

                    // compute dRdpGeo at evaluation point
                    if( tIsActiveDv.size() != 0 )
                    {
                        mSet->get_requested_IWGs()( iIWG )->compute_dRdp_FD_geometry(
                                tWStar,
                                tPerturbation,
                                tIsActiveDv,
                                tVertexIndices );
                    }
                }
            }
        }

        //------------------------------------------------------------------------------

        void Element_Time_Sideset::compute_QI()
        {
            // set physical and parametric space and time coefficients for IG element
            moris::Cell< Matrix< DDSMat > > tIsActiveDv;
            this->init_ig_geometry_interpolator( tIsActiveDv );

            // get number of IQIs
            uint tNumIQIs = mSet->get_number_of_requested_IQIs();

            // FIXME create a cell of IQI values
            moris::Cell< Matrix< DDRMat > > tQIValues( tNumIQIs );

            // loop over integration points
            uint tNumIntegPoints = mSet->get_number_of_integration_points();
            for( uint iGP = 0; iGP < tNumIntegPoints; iGP++ )
            {
                // get the ith integration point in the IG param space
                Matrix< DDRMat > tLocalIntegPoint = mSet->get_integration_points().get_column( iGP );

                // set evaluation point for interpolators (FIs and GIs)
                mSet->get_field_interpolator_manager()->
                        set_space_time_from_local_IG_point( tLocalIntegPoint );
                mSet->get_field_interpolator_manager_previous_time()->
                        set_space_time_from_local_IG_point( tLocalIntegPoint );

                // compute integration point weight
                real tWStar = mSet->get_integration_weights()( iGP ) *
                        mSet->get_field_interpolator_manager()->get_IG_geometry_interpolator()->det_J();

                // loop over the IQIs
                for( uint iIQI = 0; iIQI < tNumIQIs; iIQI++ )
                {
                    // reset IQI
                    mSet->get_requested_IQIs()( iIQI )->reset_eval_flags();

                    // compute QI at evaluation point
                    mSet->get_requested_IQIs()( iIQI )->compute_QI( tWStar );
                }
            }
        }

        //------------------------------------------------------------------------------

        void Element_Time_Sideset::compute_dQIdu()
        {
            // set physical and parametric space and time coefficients for IG element
            moris::Cell< Matrix< DDSMat > > tIsActiveDv;
            this->init_ig_geometry_interpolator( tIsActiveDv );

            // get number of IQIs
            uint tNumIQIs = mSet->get_number_of_requested_IQIs();

            // loop over integration points
            uint tNumIntegPoints = mSet->get_number_of_integration_points();
            for( uint iGP = 0; iGP < tNumIntegPoints; iGP++ )
            {
                // get the ith integration point in the IG param space
                Matrix< DDRMat > tLocalIntegPoint = mSet->get_integration_points().get_column( iGP );

                // set evaluation point for interpolators (FIs and GIs)
                mSet->get_field_interpolator_manager()->
                        set_space_time_from_local_IG_point( tLocalIntegPoint );
                mSet->get_field_interpolator_manager_previous_time()->
                        set_space_time_from_local_IG_point( tLocalIntegPoint );

                // compute integration point weight
                real tWStar = mSet->get_integration_weights()( iGP ) *
                        mSet->get_field_interpolator_manager()->get_IG_geometry_interpolator()->det_J();

                // loop over the IQIs
                for( uint iIQI = 0; iIQI < tNumIQIs; iIQI++ )
                {
                    // reset IWG
                    mSet->get_requested_IQIs()( iIQI )->reset_eval_flags();

                    // compute dQIdu at evaluation point
                    mSet->get_requested_IQIs()( iIQI )->compute_dQIdu( tWStar );
                }
            }
        }

        //------------------------------------------------------------------------------

        void Element_Time_Sideset::compute_dQIdp_explicit()
        {
            // get the vertices indices
            Matrix< IndexMat > tVertexIndices = mMasterCell->get_vertex_inds();

            // set physical and parametric space and time coefficients for IG element
            moris::Cell< Matrix< DDSMat > > tIsActiveDv;
            this->init_ig_geometry_interpolator( tIsActiveDv );

            // get number of IWGs
            uint tNumIQIs = mSet->get_number_of_requested_IQIs();

            // loop over integration points
            uint tNumIntegPoints = mSet->get_number_of_integration_points();
            for( uint iGP = 0; iGP < tNumIntegPoints; iGP++ )
            {
                // get the ith integration point in the IG param space
                Matrix< DDRMat > tLocalIntegPoint = mSet->get_integration_points().get_column( iGP );

                // set evaluation point for interpolators (FIs and GIs)
                mSet->get_field_interpolator_manager()->
                        set_space_time_from_local_IG_point( tLocalIntegPoint );
                mSet->get_field_interpolator_manager_previous_time()->
                        set_space_time_from_local_IG_point( tLocalIntegPoint );

                // compute integration point weight
                real tWStar = mSet->get_integration_weights()( iGP ) *
                        mSet->get_field_interpolator_manager()->get_IG_geometry_interpolator()->det_J();

                // loop over the IQIs
                for( uint iIQI = 0; iIQI < tNumIQIs; iIQI++ )
                {
                    // reset IWG
                    mSet->get_requested_IQIs()( iIQI )->reset_eval_flags();

                    // relative perturbation size
                    real tPerturbation = 1E-6;

                    // compute dQIdpMat at evaluation point
                    mSet->get_requested_IQIs()( iIQI )->compute_dQIdp_FD_material(
                            tWStar,
                            tPerturbation );

                    // compute dQIdpGeo at evaluation point
                    if( tIsActiveDv.size() != 0 )
                    {
                        mSet->get_requested_IQIs()( iIQI )->compute_dQIdp_FD_geometry(
                                tWStar,
                                tPerturbation,
                                tIsActiveDv,
                                tVertexIndices );
                    }
                }
            }
        }


        //------------------------------------------------------------------------------

    } /* namespace fem */
} /* namespace moris */
