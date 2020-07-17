#include <iostream>
#include "cl_FEM_Element_Sideset.hpp" //FEM/INT/src
#include "cl_FEM_Field_Interpolator_Manager.hpp" //FEM/INT/src
#include "cl_FEM_Set.hpp"   //FEM/INT/src
#include "cl_MSI_Equation_Model.hpp"

namespace moris
{
    namespace fem
    {

        //------------------------------------------------------------------------------
        Element_Sideset::Element_Sideset(
                mtk::Cell const    * aCell,
                Set                * aSet,
                Cluster            * aCluster,
                moris::moris_index   aCellIndexInCluster)
        : Element( aCell, aSet, aCluster, aCellIndexInCluster )
        {}

        //------------------------------------------------------------------------------
        Element_Sideset::~Element_Sideset(){}

        //------------------------------------------------------------------------------
        void Element_Sideset::init_ig_geometry_interpolator( uint aSideOrdinal )
        {
            // get IG geometry interpolator
            Geometry_Interpolator * tIGGI = mSet->get_field_interpolator_manager()->get_IG_geometry_interpolator();

            // set the geometry interpolator physical space and time coefficients for integration cell
            tIGGI->set_space_coeff( mMasterCell->get_cell_physical_coords_on_side_ordinal( aSideOrdinal ) );
            tIGGI->set_time_coeff( mCluster->mInterpolationElement->get_time() );

            // set the geometry interpolator param space and time coefficients for integration cell
            tIGGI->set_space_param_coeff( mCluster->get_cell_local_coords_on_side_wrt_interp_cell( mCellIndexInCluster, aSideOrdinal ) );
            tIGGI->set_time_param_coeff( {{-1.0}, {1.0}} ); //fixme
        }

        //------------------------------------------------------------------------------
        void Element_Sideset::init_ig_geometry_interpolator_with_pdv(
                uint                              aSideOrdinal,
                moris::Cell< Matrix< DDSMat > > & aIsActiveDv )
        {
            // get the vertices indices
            Matrix< IndexMat > tVertexIndices = mMasterCell->get_vertices_ind_on_side_ordinal( aSideOrdinal );

            // get the geometry XYZ values
            Matrix< DDRMat > tXYZValues =
                    mMasterCell->get_cell_physical_coords_on_side_ordinal( aSideOrdinal );

            // get the requested geo pdv types
            moris::Cell < enum PDV_Type > tGeoPdvType;
            mSet->get_ig_unique_dv_types_for_set( tGeoPdvType );

            if (tGeoPdvType.size())
            {
                // get space dimension
                uint tSpaceDim = tXYZValues.n_cols();

                // reshape the XYZ values into a cell of vectors
                moris::Cell< Matrix< DDRMat > > tPdvValueList( tSpaceDim );
                for( uint iSpaceDim = 0; iSpaceDim < tSpaceDim; iSpaceDim++ )
                {
                    tPdvValueList( iSpaceDim ) = tXYZValues.get_column( iSpaceDim );
                }

                // get the pdv values from the MSI/GEN interface
                mSet->get_equation_model()->get_design_variable_interface()->get_ig_pdv_value(
                        tVertexIndices,
                        tGeoPdvType,
                        tPdvValueList,
                        aIsActiveDv );

                // reshape the cell of vectors tPdvValueList into a matrix tPdvValues
                Matrix< DDRMat > tPdvValues;
                mSet->get_equation_model()->get_design_variable_interface()->reshape_pdv_values(
                        tPdvValueList,
                        tPdvValues );

                // get IG geometry interpolator
                Geometry_Interpolator * tIGGI = mSet->get_field_interpolator_manager()->get_IG_geometry_interpolator();

                // set the geometry interpolator physical space and time coefficients for integration cell
                tIGGI->set_space_coeff( tPdvValues );
                tIGGI->set_time_coeff ( mCluster->mInterpolationElement->get_time() );

                // set the geometry interpolator param space and time coefficients for integration cell
                tIGGI->set_space_param_coeff( mCluster->get_cell_local_coords_on_side_wrt_interp_cell( mCellIndexInCluster, aSideOrdinal ) );
                tIGGI->set_time_param_coeff( {{-1.0}, {1.0}} ); // FIXME
            }
            else
            {
                this->init_ig_geometry_interpolator( aSideOrdinal );
            }
        }

        //------------------------------------------------------------------------------
        void Element_Sideset::compute_residual()
        {
            // get treated side ordinal
            uint tSideOrd = mCluster->mMasterListOfSideOrdinals( mCellIndexInCluster );

            // set the ig geometry interpolator physical/param space and time coefficients
            this->init_ig_geometry_interpolator( tSideOrd );

            // get number of IWGs
            uint tNumIWGs = mSet->get_number_of_requested_IWGs();

            // loop over integration points
            uint tNumIntegPoints = mSet->get_number_of_integration_points();

            for( uint iGP = 0; iGP < tNumIntegPoints; iGP++ )
            {
                // get integration point location in the reference surface
                Matrix< DDRMat > tLocalIntegPoint = mSet->get_integration_points().get_column( iGP );

                // set evaluation point for interpolators (FIs and GIs)
                mSet->get_field_interpolator_manager()->set_space_time_from_local_IG_point( tLocalIntegPoint );

                // compute the integration point weight
                real tWStar = mSet->get_integration_weights()( iGP ) *
                        mSet->get_field_interpolator_manager()->get_IG_geometry_interpolator()->det_J();

                // get the normal from mesh
                Matrix< DDRMat > tNormal = mCluster->get_side_normal( mMasterCell, tSideOrd );

                // loop over the IWGs
                for( uint iIWG = 0; iIWG < tNumIWGs; iIWG++ )
                {
                    // reset IWG
                    mSet->get_requested_IWGs()( iIWG )->reset_eval_flags();

                    // FIXME
                    mSet->get_requested_IWGs()( iIWG )->set_nodal_weak_bcs( mCluster->mInterpolationElement->get_weak_bcs() );

                    // set the normal for the IWG
                    mSet->get_requested_IWGs()( iIWG )->set_normal( tNormal );

                    // compute residual at integration point
                    mSet->get_requested_IWGs()( iIWG )->compute_residual( tWStar );

                    // compute jacobian at evaluation point
                    // compute off-diagonal jacobian for staggered solve
                    mSet->get_requested_IWGs()( iIWG )->compute_jacobian( tWStar );
                }
            }
        }

        //------------------------------------------------------------------------------
        void Element_Sideset::compute_jacobian()
        {
            // get treated side ordinal
            uint tSideOrd = mCluster->mMasterListOfSideOrdinals( mCellIndexInCluster );

            // set the ig geometry interpolator physical/param space and time coefficients
            this->init_ig_geometry_interpolator( tSideOrd );

            // get number of IWGs
            uint tNumIWGs = mSet->get_number_of_requested_IWGs();

            // loop over integration points
            uint tNumIntegPoints = mSet->get_number_of_integration_points();

            for( uint iGP = 0; iGP < tNumIntegPoints; iGP++ )
            {
                // get integration point location in the reference surface
                Matrix< DDRMat > tLocalIntegPoint = mSet->get_integration_points().get_column( iGP );

                // set evaluation point for interpolators (FIs and GIs)
                mSet->get_field_interpolator_manager()->set_space_time_from_local_IG_point( tLocalIntegPoint );

                // compute integration point weight
                real tWStar = mSet->get_integration_weights()( iGP ) *
                        mSet->get_field_interpolator_manager()->get_IG_geometry_interpolator()->det_J();

                // get the normal from mesh
                Matrix< DDRMat > tNormal = mCluster->get_side_normal( mMasterCell, tSideOrd );

                // loop over the IWGs
                for( uint iIWG = 0; iIWG < tNumIWGs; iIWG++ )
                {
                    // reset IWG
                    mSet->get_requested_IWGs()( iIWG )->reset_eval_flags();

                    // FIXME set BCs
                    mSet->get_requested_IWGs()( iIWG )->set_nodal_weak_bcs( mCluster->mInterpolationElement->get_weak_bcs() );

                    // set the normal for the IWG
                    mSet->get_requested_IWGs()( iIWG )->set_normal( tNormal );

                    // compute jacobian at evaluation point
                    mSet->get_requested_IWGs()( iIWG )->compute_jacobian( tWStar );
                }
            }
        }

        //------------------------------------------------------------------------------
        void Element_Sideset::compute_jacobian_and_residual()
        {
            // get treated side ordinal
            uint tSideOrd = mCluster->mMasterListOfSideOrdinals( mCellIndexInCluster );

            // set the ig geometry interpolator physical/param space and time coefficients
            this->init_ig_geometry_interpolator( tSideOrd );

            // get number of IWGs
            uint tNumIWGs = mSet->get_number_of_requested_IWGs();

            // loop over integration points
            uint tNumIntegPoints = mSet->get_number_of_integration_points();
            for( uint iGP = 0; iGP < tNumIntegPoints; iGP++ )
            {
                // get integration point location in the reference surface
                Matrix< DDRMat > tLocalIntegPoint = mSet->get_integration_points().get_column( iGP );

                // set evaluation point for interpolators (FIs and GIs)
                mSet->get_field_interpolator_manager()->set_space_time_from_local_IG_point( tLocalIntegPoint );

                // compute integration point weight
                real tWStar = mSet->get_integration_weights()( iGP ) *
                        mSet->get_field_interpolator_manager()->get_IG_geometry_interpolator()->det_J();

                // get the normal from mesh
                Matrix< DDRMat > tNormal = mCluster->get_side_normal( mMasterCell, tSideOrd );

                // loop over the IWGs
                for( uint iIWG = 0; iIWG < tNumIWGs; iIWG++ )
                {
                    // reset IWG
                    mSet->get_requested_IWGs()( iIWG )->reset_eval_flags();

                    // FIXME set BCs
                    mSet->get_requested_IWGs()( iIWG )->set_nodal_weak_bcs( mCluster->mInterpolationElement->get_weak_bcs() );

                    // set the normal for the IWG
                    mSet->get_requested_IWGs()( iIWG )->set_normal( tNormal );

                    // compute residual at integration point
                    mSet->get_requested_IWGs()( iIWG )->compute_residual( tWStar );

                    // compute jacobian at evaluation point
                    mSet->get_requested_IWGs()( iIWG )->compute_jacobian( tWStar );
                }
            }
        }

        //------------------------------------------------------------------------------
        void Element_Sideset::compute_dRdp()
        {
            // get treated side ordinal
            uint tSideOrd = mCluster->mMasterListOfSideOrdinals( mCellIndexInCluster );

            // get the vertices indices
            Matrix< IndexMat > tVertexIndices = mMasterCell->get_vertices_ind_on_side_ordinal( tSideOrd );

            // set the IG geometry interpolator physical/param space and time coefficients
            moris::Cell< Matrix< DDSMat > > tIsActiveDv;
            this->init_ig_geometry_interpolator_with_pdv(
                    tSideOrd,
                    tIsActiveDv );

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
                    mSet->get_requested_IWGs()( iIWG )->set_nodal_weak_bcs( mCluster->mInterpolationElement->get_weak_bcs() );

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
        void Element_Sideset::compute_QI()
        {
            // get treated side ordinal
            uint tSideOrd = mCluster->mMasterListOfSideOrdinals( mCellIndexInCluster );

            // set the ig geometry interpolator physical/param space and time coefficients
            this->init_ig_geometry_interpolator( tSideOrd );

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
        void Element_Sideset::compute_dQIdp_explicit()
        {
            // get treated side ordinal
            uint tSideOrd = mCluster->mMasterListOfSideOrdinals( mCellIndexInCluster );

            // get the vertices indices
            Matrix< IndexMat > tVertexIndices = mMasterCell->get_vertices_ind_on_side_ordinal( tSideOrd );

            // set the IG geometry interpolator physical/param space and time coefficients
            moris::Cell< Matrix< DDSMat > > tIsActiveDv;
            this->init_ig_geometry_interpolator_with_pdv(
                    tSideOrd,
                    tIsActiveDv );

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
        void Element_Sideset::compute_quantity_of_interest_global(
                const uint             aMeshIndex,
                enum  vis::Output_Type aOutputType )
        {
            // get treated side ordinal
            uint tSideOrd = mCluster->mMasterListOfSideOrdinals( mCellIndexInCluster );

            // set the ig geometry interpolator physical/param space and time coefficients
            this->init_ig_geometry_interpolator( tSideOrd );

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
        void Element_Sideset::compute_quantity_of_interest_nodal(
                const uint             aMeshIndex,
                enum  vis::Output_Type aOutputType )
        {
            // get treated side ordinal
            uint tSideOrd = mCluster->mMasterListOfSideOrdinals( mCellIndexInCluster );

            // set the ig geometry interpolator physical/param space and time coefficients
            this->init_ig_geometry_interpolator( tSideOrd );

            // get the vertices
            moris::Cell< mtk::Vertex const * > tVertices = mMasterCell->get_vertices_on_side_ordinal( tSideOrd );

            // loop over the vertices
            uint tNumNodes = tVertices.size();
            for( uint iVertex = 0; iVertex < tNumNodes; iVertex++ )
            {
                // get the ith vertex coordinates in the IP param space
                Matrix< DDRMat > tGlobalIntegPoint = mCluster->get_primary_cell_local_coords_on_side_wrt_interp_cell( mCellIndexInCluster ).get_row( iVertex );
                tGlobalIntegPoint.resize( 1, tGlobalIntegPoint.numel() + 1 );
                tGlobalIntegPoint( tGlobalIntegPoint.numel() - 1 ) = mCluster->mInterpolationElement->get_time()( 0 );
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
        void Element_Sideset::compute_quantity_of_interest_elemental(
                const uint             aMeshIndex,
                enum  vis::Output_Type aOutputType )
        {
            // get treated side ordinal
            uint tSideOrd = mCluster->mMasterListOfSideOrdinals( mCellIndexInCluster );

            // set the ig geometry interpolator physical/param space and time coefficients
            this->init_ig_geometry_interpolator( tSideOrd );

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
        real Element_Sideset::compute_volume( mtk::Master_Slave aIsMaster )
        {
            // get treated side ordinal
            uint tSideOrd = mCluster->mMasterListOfSideOrdinals( mCellIndexInCluster );

            // set the ig geometry interpolator physical/param space and time coefficients
            this->init_ig_geometry_interpolator( tSideOrd );

            //get number of integration points
            uint tNumOfIntegPoints = mSet->get_number_of_integration_points();

            // init volume
            real tVolume = 0;

            // get IG geometry interpolator
            Geometry_Interpolator * tIGGI =
                    mSet->get_field_interpolator_manager( aIsMaster )->get_IG_geometry_interpolator();

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
