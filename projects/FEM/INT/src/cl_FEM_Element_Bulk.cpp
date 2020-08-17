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
        //------------------------------------------------------------------------------

        Element_Bulk::Element_Bulk(
                mtk::Cell const    * aCell,
                Set                * aSet,
                Cluster            * aCluster,
                moris::moris_index   aCellIndexInCluster )
        : Element( aCell, aSet, aCluster, aCellIndexInCluster )
        {}

        //------------------------------------------------------------------------------

        Element_Bulk::~Element_Bulk(){}

        //----------------------------------------------------------------------

        void Element_Bulk::init_ig_geometry_interpolator(
                moris::Cell< Matrix< DDSMat > > & aIsActiveDv )
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

            // get the requested geo pdv types
            moris::Cell < enum PDV_Type > tGeoPdvType;
            mSet->get_ig_unique_dv_types_for_set( tGeoPdvType );

            // Determine if there are IG PDVs
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
            moris::Cell< Matrix< DDSMat > > tIsActiveDv;
            this->init_ig_geometry_interpolator( tIsActiveDv );

            // get number of IWGs
            uint tNumIWGs = mSet->get_number_of_requested_IWGs();

            // loop over integration points
            uint tNumIntegPoints = mSet->get_number_of_integration_points();
            for( uint iGP = 0; iGP < tNumIntegPoints; iGP++ )
            {
                // get the current integration point in the IG param space
                Matrix< DDRMat > tLocalIntegPoint = mSet->get_integration_points().get_column( iGP );

                // set evaluation point for interpolators (FIs and GIs)
                mSet->get_field_interpolator_manager()->set_space_time_from_local_IG_point( tLocalIntegPoint );

                // compute integration point weight
                real tWStar = mSet->get_integration_weights()( iGP ) *
                        mSet->get_field_interpolator_manager()->get_IG_geometry_interpolator()->det_J();

                // loop over the IWGs
                for( uint iIWG = 0; iIWG < tNumIWGs; iIWG++ )
                {
                    // reset IWG
                    mSet->get_requested_IWGs()( iIWG )->reset_eval_flags();

                    // FIXME: enforced nodal weak bcs
                    mSet->get_requested_IWGs()( iIWG )->set_nodal_weak_bcs(
                            mCluster->mInterpolationElement->get_weak_bcs() );

                    // compute residual at evaluation point
                    mSet->get_requested_IWGs()( iIWG )->compute_residual( tWStar );

                    // compute off-diagonal jacobian for staggered solve
                    mSet->get_requested_IWGs()( iIWG )->compute_jacobian( tWStar );
                }
            }
        }

        //------------------------------------------------------------------------------

        void Element_Bulk::compute_jacobian()
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
                // get the ith integration point in the IG param space
                Matrix< DDRMat > tLocalIntegPoint = mSet->get_integration_points().get_column( iGP );

                // set evaluation point for interpolators (FIs and GIs)
                mSet->get_field_interpolator_manager()->set_space_time_from_local_IG_point( tLocalIntegPoint );

                // compute integration point weight
                real tWStar = mSet->get_integration_weights()( iGP ) *
                        mSet->get_field_interpolator_manager()->get_IG_geometry_interpolator()->det_J();

                // loop over the IWGs
                for( uint iIWG = 0; iIWG < tNumIWGs; iIWG++ )
                {
                    // reset IWG
                    mSet->get_requested_IWGs()( iIWG )->reset_eval_flags();

                    // FIXME set nodal weak BCs
                    mSet->get_requested_IWGs()( iIWG )->set_nodal_weak_bcs(
                            mCluster->mInterpolationElement->get_weak_bcs() );

                    // compute jacobian at evaluation point
                    mSet->get_requested_IWGs()( iIWG )->compute_jacobian( tWStar );
                }
            }
        }

        //------------------------------------------------------------------------------

        void Element_Bulk::compute_jacobian_and_residual()
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
                // get the ith integration point in the IG param space
                Matrix< DDRMat > tLocalIntegPoint = mSet->get_integration_points().get_column( iGP );

                // set evaluation point for interpolators (FIs and GIs)
                mSet->get_field_interpolator_manager()->set_space_time_from_local_IG_point( tLocalIntegPoint );

                // compute integration point weight
                real tWStar = mSet->get_integration_weights()( iGP ) *
                        mSet->get_field_interpolator_manager()->get_IG_geometry_interpolator()->det_J();

                // loop over the IWGs
                for( uint iIWG = 0; iIWG < tNumIWGs; iIWG++ )
                {
                    // reset IWG
                    mSet->get_requested_IWGs()( iIWG )->reset_eval_flags();

                    // FIXME set nodal weak BCs
                    mSet->get_requested_IWGs()( iIWG )->set_nodal_weak_bcs(
                            mCluster->mInterpolationElement->get_weak_bcs() );

                    // compute residual at evaluation point
                    mSet->get_requested_IWGs()( iIWG )->compute_residual( tWStar );

                    // compute jacobian at evaluation point
                    mSet->get_requested_IWGs()( iIWG )->compute_jacobian( tWStar );
                }
            }
        }

        //------------------------------------------------------------------------------

        void Element_Bulk::compute_dRdp()
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
                // get the ith integration point in the IG param space
                Matrix< DDRMat > tLocalIntegPoint = mSet->get_integration_points().get_column( iGP );

                // set evaluation point for interpolators (FIs and GIs)
                mSet->get_field_interpolator_manager()->set_space_time_from_local_IG_point( tLocalIntegPoint );

                // compute integration point weight
                real tWStar = mSet->get_integration_weights()( iGP ) *
                        mSet->get_field_interpolator_manager()->get_IG_geometry_interpolator()->det_J();

                // loop over the IWGs
                for( uint iIWG = 0; iIWG < tNumIWGs; iIWG++ )
                {
                    // reset IWG
                    mSet->get_requested_IWGs()( iIWG )->reset_eval_flags();

                    // FIXME set nodal weak BCs
                    mSet->get_requested_IWGs()( iIWG )->set_nodal_weak_bcs(
                            mCluster->mInterpolationElement->get_weak_bcs() );

                    // compute dRdpMat at evaluation point
                    mSet->get_requested_IWGs()( iIWG )->compute_dRdp( tWStar );
                }
            }
        }

        //------------------------------------------------------------------------------

        void Element_Bulk::compute_dRdp_FD()
        {
            // get finite difference scheme type
            fem::FDScheme_Type tFDScheme =
                    mSet->get_finite_difference_scheme_for_sensitivity_analysis();

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
                // get the ith integration point in the IG param space
                Matrix< DDRMat > tLocalIntegPoint = mSet->get_integration_points().get_column( iGP );

                // set evaluation point for interpolators (FIs and GIs)
                mSet->get_field_interpolator_manager()->set_space_time_from_local_IG_point( tLocalIntegPoint );

                // compute integration point weight
                real tWStar = mSet->get_integration_weights()( iGP ) *
                        mSet->get_field_interpolator_manager()->get_IG_geometry_interpolator()->det_J();

                // loop over the IWGs
                for( uint iIWG = 0; iIWG < tNumIWGs; iIWG++ )
                {
                    // reset IWG
                    mSet->get_requested_IWGs()( iIWG )->reset_eval_flags();

                    // FIXME set nodal weak BCs
                    mSet->get_requested_IWGs()( iIWG )->set_nodal_weak_bcs(
                            mCluster->mInterpolationElement->get_weak_bcs() );

                    // set a perturbation size
                    real tPerturbation = 1E-6;

                    // compute dRdpMat at evaluation point
                    mSet->get_requested_IWGs()( iIWG )->compute_dRdp_FD_material(
                            tWStar,
                            tPerturbation,
                            tFDScheme );

                    // compute dRdpGeo at evaluation point
                    if( tIsActiveDv.size() != 0 )
                    {
                        mSet->get_requested_IWGs()( iIWG )->compute_dRdp_FD_geometry(
                                tWStar,
                                tPerturbation,
                                tIsActiveDv,
                                tVertexIndices,
                                tFDScheme );
                    }
                }
            }
        }

        //------------------------------------------------------------------------------

        void Element_Bulk::compute_QI()
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
                mSet->get_field_interpolator_manager()->set_space_time_from_local_IG_point( tLocalIntegPoint );

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

        void Element_Bulk::compute_dQIdp_explicit()
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
                mSet->get_field_interpolator_manager()->set_space_time_from_local_IG_point( tLocalIntegPoint );

                // compute integration point weight
                real tWStar = mSet->get_integration_weights()( iGP ) *
                        mSet->get_field_interpolator_manager()->get_IG_geometry_interpolator()->det_J();

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
                mSet->get_field_interpolator_manager()->set_space_time_from_local_IG_point( tLocalIntegPoint );

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
                            tPerturbation,
                            tFDScheme );

                    // compute dQIdpGeo at evaluation point
                    if( tIsActiveDv.size() != 0 )
                    {
                        mSet->get_requested_IQIs()( iIQI )->compute_dQIdp_FD_geometry(
                                tWStar,
                                tPerturbation,
                                tIsActiveDv,
                                tVertexIndices,
                                tFDScheme );
                    }
                }
            }
        }

        //------------------------------------------------------------------------------

        void Element_Bulk::compute_dQIdu()
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
                mSet->get_field_interpolator_manager()->set_space_time_from_local_IG_point( tLocalIntegPoint );

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

        void Element_Bulk::compute_quantity_of_interest_global(
                const uint             aMeshIndex,
                enum  vis::Output_Type aOutputType )
        {
            // set physical and parametric space and time coefficients for IG element
            moris::Cell< Matrix< DDSMat > > tIsActiveDv;
            this->init_ig_geometry_interpolator( tIsActiveDv );

            // loop over integration points
            uint tNumIntegPoints = mSet->get_number_of_integration_points();
            for( uint iGP = 0; iGP < tNumIntegPoints; iGP++ )
            {
                // get the ith integration point in the IG param space
                Matrix< DDRMat > tLocalIntegPoint = mSet->get_integration_points().get_column( iGP );

                // set evaluation point for interpolators (FIs and GIs)
                mSet->get_field_interpolator_manager()->set_space_time_from_local_IG_point( tLocalIntegPoint );

                // compute integration point weight
                real tWStar = mSet->get_integration_weights()( iGP ) *
                        mSet->get_field_interpolator_manager()->get_IG_geometry_interpolator()->det_J();

                // reset the requested IQI
                mSet->get_IQI_for_vis( aOutputType )->reset_eval_flags();

                // compute quantity of interest at evaluation point
                Matrix< DDRMat > tQIValue;
                mSet->get_IQI_for_vis( aOutputType )->compute_QI( tQIValue );

                // FIXME assemble on the set here or inside the compute QI?
                *( mSet->mSetGlobalValues ) += tQIValue( 0 ) * tWStar;
            }
        }

        //------------------------------------------------------------------------------

        void Element_Bulk::compute_quantity_of_interest_nodal(
                const uint aMeshIndex,
                enum vis::Output_Type aOutputType )
        {
            // set physical and parametric space and time coefficients for IG element
            moris::Cell< Matrix< DDSMat > > tIsActiveDv;
            this->init_ig_geometry_interpolator( tIsActiveDv );

            // get the vertices
            moris::Cell< mtk::Vertex * > tVertices = mMasterCell->get_vertex_pointers();

            // loop over the vertices
            uint tNumNodes = tVertices.size();
            for( uint iVertex = 0; iVertex < tNumNodes; iVertex++ )
            {
                // get the ith vertex coordinates in the IP param space
                Matrix< DDRMat > tGlobalIntegPoint =
                        mCluster->get_primary_cell_local_coords_on_side_wrt_interp_cell( mCellIndexInCluster ).get_row( iVertex );

                tGlobalIntegPoint.resize( 1, tGlobalIntegPoint.numel() + 1 );
                tGlobalIntegPoint( tGlobalIntegPoint.numel() - 1 ) =
                        mCluster->mInterpolationElement->get_time()( 0 ); //FIXME plot at beginning of slab
                tGlobalIntegPoint = trans( tGlobalIntegPoint );

                // set vertex coordinates for field interpolator
                mSet->get_field_interpolator_manager()->set_space_time( tGlobalIntegPoint );

                // reset the requested IQI
                mSet->get_IQI_for_vis( aOutputType )->reset_eval_flags();

                // compute quantity of interest at evaluation point
                Matrix< DDRMat > tQIValue;
                mSet->get_IQI_for_vis( aOutputType )->compute_QI( tQIValue );

                // FIXME assemble on the set here or inside the compute QI?
                // FIXME add up on shared node and divide or overwrite
                (*mSet->mSetNodalValues)( tVertices( iVertex )->get_index(), 0 ) += tQIValue( 0 );

                mSet->mSetNodalCounter( tVertices( iVertex )->get_index(), 0 ) += 1;
            }
        }

        //------------------------------------------------------------------------------

        void Element_Bulk::compute_quantity_of_interest_elemental(
                const uint aMeshIndex,
                enum vis::Output_Type aOutputType )
        {
            // set physical and parametric space and time coefficients for IG element
            moris::Cell< Matrix< DDSMat > > tIsActiveDv;
            this->init_ig_geometry_interpolator( tIsActiveDv );

            // loop over integration points
            uint tNumIntegPoints = mSet->get_number_of_integration_points();
            for( uint iGP = 0; iGP < tNumIntegPoints; iGP++ )
            {
                // get the ith integration point in the IG param space
                Matrix< DDRMat > tLocalIntegPoint = mSet->get_integration_points().get_column( iGP );

                // set evaluation point for interpolators (FIs and GIs)
                mSet->get_field_interpolator_manager()->set_space_time_from_local_IG_point( tLocalIntegPoint );

                // compute integration point weight
                real tWStar = mSet->get_integration_weights()( iGP ) *
                        mSet->get_field_interpolator_manager()->get_IG_geometry_interpolator()->det_J();

                // reset the requested IQI
                mSet->get_IQI_for_vis( aOutputType )->reset_eval_flags();

                // compute quantity of interest at evaluation point
                Matrix< DDRMat > tQIValue;
                mSet->get_IQI_for_vis( aOutputType )->compute_QI( tQIValue );

                // FIXME assemble on the set here or inside the compute QI?
                ( *mSet->mSetElementalValues )( mSet->mCellAssemblyMap( aMeshIndex )( mMasterCell->get_index() ), 0 ) += 
                        tQIValue( 0 ) * tWStar / tNumIntegPoints;
            }
        }

        //------------------------------------------------------------------------------

        real Element_Bulk::compute_volume( mtk::Master_Slave aIsMaster )
        {
            // set physical and parametric space and time coefficients for IG element
            moris::Cell< Matrix< DDSMat > > tIsActiveDv;
            this->init_ig_geometry_interpolator( tIsActiveDv );

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
