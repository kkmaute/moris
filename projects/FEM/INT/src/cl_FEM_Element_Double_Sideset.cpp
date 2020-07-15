#include <iostream>

#include "cl_FEM_Element_Double_Sideset.hpp" //FEM/INT/src
#include "cl_FEM_Set.hpp"                    //FEM/INT/src
#include "cl_FEM_Field_Interpolator_Manager.hpp" //FEM/INT/src
#include "fn_FEM_Rotation_Matrix.hpp"        //FEM/INT/src
#include "cl_MSI_Equation_Model.hpp"

namespace moris
{
    namespace fem
    {

        //------------------------------------------------------------------------------

        Element_Double_Sideset::Element_Double_Sideset(
                mtk::Cell const  * aMasterIGCell,
                mtk::Cell const  * aSlaveIGCell,
                Set              * aSet,
                Cluster          * aCluster,
                moris::moris_index aCellIndexInCluster )
        : Element(
                aMasterIGCell,
                aSlaveIGCell,
                aSet,
                aCluster,
                aCellIndexInCluster )
        {}

        //------------------------------------------------------------------------------

        Element_Double_Sideset::~Element_Double_Sideset(){}

        //------------------------------------------------------------------------------

        void Element_Double_Sideset::init_ig_geometry_interpolator(
                uint aMasterSideOrdinal,
                uint aSlaveSideOrdinal )
        {
            // get master IG geometry interpolator
            Geometry_Interpolator * tMasterIGGI =
                    mSet->get_field_interpolator_manager( mtk::Master_Slave::MASTER )->
                    get_IG_geometry_interpolator();

            // get slave IG geometry interpolator
            Geometry_Interpolator * tSlaveIGGI =
                    mSet->get_field_interpolator_manager( mtk::Master_Slave::SLAVE )->
                    get_IG_geometry_interpolator();

            // set the geometry interpolator physical space and time coefficients for master integration cell
            tMasterIGGI->set_space_coeff( mMasterCell->get_cell_physical_coords_on_side_ordinal( aMasterSideOrdinal ) );
            tMasterIGGI->set_time_coeff(  mCluster->mInterpolationElement->get_time() );

            // set the geometry interpolator physical space and time coefficients for slave integration cell
            tSlaveIGGI->set_space_coeff( mSlaveCell->get_cell_physical_coords_on_side_ordinal( aSlaveSideOrdinal ) );
            tSlaveIGGI->set_time_coeff( mCluster->mInterpolationElement->get_time() );

            // set the geometry interpolator param space and time coefficients for master integration cell
            tMasterIGGI->set_space_param_coeff( mCluster->get_cell_local_coords_on_side_wrt_interp_cell(
                    mCellIndexInCluster,
                    aMasterSideOrdinal,
                    mtk::Master_Slave::MASTER ) );

            tMasterIGGI->set_time_param_coeff( {{-1.0}, {1.0}} ); //fixme

            // set the geometry interpolator param space and time coefficients for slave integration cell
            tSlaveIGGI->set_space_param_coeff( mCluster->get_cell_local_coords_on_side_wrt_interp_cell(
                    mCellIndexInCluster,
                    aSlaveSideOrdinal,
                    mtk::Master_Slave::SLAVE ) );

            tSlaveIGGI->set_time_param_coeff( {{-1.0}, {1.0}} ); //fixme
        }

        //------------------------------------------------------------------------------
        void Element_Double_Sideset::init_ig_geometry_interpolator_with_pdv(
                uint                              aMasterSideOrdinal,
                uint                              aSlaveSideOrdinal,
                moris::Cell< Matrix< DDSMat > > & aMasterIsActiveDv,
                moris::Cell< Matrix< DDSMat > > & aSlaveIsActiveDv )
        {
            // get the vertices indices
            Matrix< IndexMat > tMasterVertexIndices = mMasterCell->get_vertex_inds();
            Matrix< IndexMat > tSlaveVertexIndices  = mSlaveCell->get_vertex_inds();

            // get the geometry XYZ values
            Matrix< DDRMat > tMasterXYZValues =
                    mMasterCell->get_cell_physical_coords_on_side_ordinal( aMasterSideOrdinal );
            Matrix< DDRMat > tSlaveXYZValues =
                    mSlaveCell->get_cell_physical_coords_on_side_ordinal( aSlaveSideOrdinal );

            // get the requested geo pdv types
            moris::Cell < enum PDV_Type > tGeoPdvType;
            mSet->get_ig_unique_dv_types_for_set( tGeoPdvType );

            // get space dimension
            uint tSpaceDim = tMasterXYZValues.n_cols();

            // reshape the XYZ values into a cell of vectors
            moris::Cell< Matrix< DDRMat > > tMasterPdvValueList( tSpaceDim );
            moris::Cell< Matrix< DDRMat > > tSlavePdvValueList( tSpaceDim );
            for( uint iSpaceDim = 0; iSpaceDim < tSpaceDim; iSpaceDim++ )
            {
                tMasterPdvValueList( iSpaceDim ) = tMasterXYZValues.get_column( iSpaceDim );
                tSlavePdvValueList( iSpaceDim )  = tSlaveXYZValues.get_column( iSpaceDim );
            }

            // get the pdv values from the MSI/GEN interface
            mSet->get_equation_model()->get_design_variable_interface()->get_ig_pdv_value(
                    tMasterVertexIndices,
                    tGeoPdvType,
                    tMasterPdvValueList,
                    aMasterIsActiveDv );

            mSet->get_equation_model()->get_design_variable_interface()->get_ig_pdv_value(
                    tSlaveVertexIndices,
                    tGeoPdvType,
                    tSlavePdvValueList,
                    aSlaveIsActiveDv );

            // reshape the cell of vectors tPdvValueList into a matrix tPdvValues
            Matrix< DDRMat > tMasterPdvValues;
            Matrix< DDRMat > tSlavePdvValues;

            mSet->get_equation_model()->get_design_variable_interface()->reshape_pdv_values(
                    tMasterPdvValueList,
                    tMasterPdvValues );

            mSet->get_equation_model()->get_design_variable_interface()->reshape_pdv_values(
                    tSlavePdvValueList,
                    tSlavePdvValues );

            // get master IG geometry interpolator
            Geometry_Interpolator * tMasterIGGI =
                    mSet->get_field_interpolator_manager( mtk::Master_Slave::MASTER )->
                    get_IG_geometry_interpolator();

            // get slave IG geometry interpolator
            Geometry_Interpolator * tSlaveIGGI =
                    mSet->get_field_interpolator_manager( mtk::Master_Slave::SLAVE )->
                    get_IG_geometry_interpolator();

            // set the geometry interpolator physical space and time coefficients for master integration cell
            tMasterIGGI->set_space_coeff( tMasterPdvValues );
            tMasterIGGI->set_time_coeff( mCluster->mInterpolationElement->get_time() );

            // set the geometry interpolator physical space and time coefficients for slave integration cell
            tSlaveIGGI->set_space_coeff( tSlavePdvValues );
            tSlaveIGGI->set_time_coeff( mCluster->mInterpolationElement->get_time() );

            // set the geometry interpolator param space and time coefficients for master integration cell
            tMasterIGGI->set_space_param_coeff( mCluster->get_cell_local_coords_on_side_wrt_interp_cell(
                    mCellIndexInCluster,
                    aMasterSideOrdinal,
                    mtk::Master_Slave::MASTER ) );

            tMasterIGGI->set_time_param_coeff( {{-1.0}, {1.0}} ); //fixme

            // set the geometry interpolator param space and time coefficients for slave integration cell
            tSlaveIGGI->set_space_param_coeff( mCluster->get_cell_local_coords_on_side_wrt_interp_cell(
                    mCellIndexInCluster,
                    aSlaveSideOrdinal,
                    mtk::Master_Slave::SLAVE ) );

            tSlaveIGGI->set_time_param_coeff( {{-1.0}, {1.0}} ); //fixme
        }

        //------------------------------------------------------------------------------

        void Element_Double_Sideset::compute_residual()
        {
            // get treated side ordinal on the master and on the slave
            uint tMasterSideOrd = mCluster->mMasterListOfSideOrdinals( mCellIndexInCluster );
            uint tSlaveSideOrd  = mCluster->mSlaveListOfSideOrdinals( mCellIndexInCluster );

            // set the master/slave ig geometry interpolator physical/parm space and time coefficients
            this->init_ig_geometry_interpolator( tMasterSideOrd, tSlaveSideOrd );

            // get first corresponding node from master to slave
            moris::mtk::Vertex const * tSlaveNode =
                    mCluster->get_left_vertex_pair(
                            mMasterCell->get_vertices_on_side_ordinal( tMasterSideOrd )( 0 ) );

            moris_index tSlaveNodeOrdOnSide =
                    mCluster->get_right_vertex_ordinal_on_facet( mCellIndexInCluster, tSlaveNode );

            // get rotation matrix from left to right
            Matrix< DDRMat> tR;
            rotation_matrix( mSet->get_IG_geometry_type(), tSlaveNodeOrdOnSide, tR );

            // get number of IWGs
            uint tNumIWGs = mSet->get_number_of_requested_IWGs();

            // loop over the integration points
            uint tNumIntegPoints = mSet->get_number_of_integration_points();

            for( uint iGP = 0; iGP < tNumIntegPoints; iGP++ )
            {
                // get local integration point for the master integration cell
                Matrix< DDRMat > tMasterLocalIntegPoint = mSet->get_integration_points().get_column( iGP );

                // get local integration point for the slave integration cell
                Matrix< DDRMat > tSlaveLocalIntegPoint = tMasterLocalIntegPoint;
                tSlaveLocalIntegPoint({0,tSlaveLocalIntegPoint.numel()-2},{0,0}) =
                        tR * tMasterLocalIntegPoint({0,tSlaveLocalIntegPoint.numel()-2},{0,0}); //fixme better way?

                // set evaluation point for master and slave interpolators
                mSet->get_field_interpolator_manager( mtk::Master_Slave::MASTER )->
                        set_space_time_from_local_IG_point( tMasterLocalIntegPoint );
                mSet->get_field_interpolator_manager( mtk::Master_Slave::SLAVE )->
                        set_space_time_from_local_IG_point( tSlaveLocalIntegPoint );

                // compute the integration point weight
                real tWStar = mSet->get_integration_weights()( iGP ) *
                        mSet->get_field_interpolator_manager()->get_IG_geometry_interpolator()->det_J();

                // get the normal from mesh
                Matrix< DDRMat > tNormal = mCluster->get_side_normal( mMasterCell, tMasterSideOrd );

                // loop over the IWGs
                for( uint iIWG = 0; iIWG < tNumIWGs; iIWG++ )
                {
                    // reset IWG
                    mSet->get_requested_IWGs()( iIWG )->reset_eval_flags();

                    // set the normal for the IWG
                    mSet->get_requested_IWGs()( iIWG )->set_normal( tNormal );

                    // compute residual at integration point
                    mSet->get_requested_IWGs()( iIWG )->compute_residual( tWStar );

                    // compute jacobian at integration point
                    // compute off-diagonal jacobian for staggered solve
                    mSet->get_requested_IWGs()( iIWG )->compute_jacobian( tWStar );
                }
            }
        }

        //------------------------------------------------------------------------------
        void Element_Double_Sideset::compute_jacobian()
        {
            // get treated side ordinal on the master and on the slave
            uint tMasterSideOrd = mCluster->mMasterListOfSideOrdinals( mCellIndexInCluster );
            uint tSlaveSideOrd  = mCluster->mSlaveListOfSideOrdinals( mCellIndexInCluster );

            // set the master/slave ig geometry interpolator physical/parm space and time coefficients
            this->init_ig_geometry_interpolator( tMasterSideOrd, tSlaveSideOrd );

            // get first corresponding node from master to slave
            moris::mtk::Vertex const * tSlaveNode =
                    mCluster->get_left_vertex_pair(
                            mMasterCell->get_vertices_on_side_ordinal( tMasterSideOrd )( 0 ) );

            moris_index tSlaveNodeOrdOnSide =
                    mCluster->get_right_vertex_ordinal_on_facet(mCellIndexInCluster,tSlaveNode);

            // get rotation matrix from left to right
            Matrix< DDRMat> tR;
            rotation_matrix( mSet->get_IG_geometry_type(), tSlaveNodeOrdOnSide, tR );

            // get number of IWGs
            uint tNumIWGs = mSet->get_number_of_requested_IWGs();

            // loop over the integration points
            uint tNumIntegPoints = mSet->get_number_of_integration_points();

            for( uint iGP = 0; iGP < tNumIntegPoints; iGP++ )
            {
                // get local integration point for the master integration cell
                Matrix< DDRMat > tMasterLocalIntegPoint = mSet->get_integration_points().get_column( iGP );

                // get local integration point for the slave integration cell
                Matrix< DDRMat > tSlaveLocalIntegPoint = tMasterLocalIntegPoint;
                tSlaveLocalIntegPoint({0,tMasterLocalIntegPoint.numel()-2},{0,0}) =
                        tR * tMasterLocalIntegPoint({0,tSlaveLocalIntegPoint.numel()-2},{0,0}); //fixme better way?

                // set evaluation point for master and slave interpolators
                mSet->get_field_interpolator_manager( mtk::Master_Slave::MASTER )->
                        set_space_time_from_local_IG_point( tMasterLocalIntegPoint );
                mSet->get_field_interpolator_manager( mtk::Master_Slave::SLAVE )->
                        set_space_time_from_local_IG_point( tSlaveLocalIntegPoint );

                // compute the integration point weight
                real tWStar = mSet->get_integration_weights()( iGP ) *
                        mSet->get_field_interpolator_manager()->get_IG_geometry_interpolator()->det_J();

                // get the normal from mesh and set if for the IWG
                Matrix< DDRMat > tNormal = mCluster->get_side_normal( mMasterCell, tMasterSideOrd );

                // loop over the IWGs
                for( uint iIWG = 0; iIWG < tNumIWGs; iIWG++ )
                {
                    // reset IWG
                    mSet->get_requested_IWGs()( iIWG )->reset_eval_flags();

                    // set the normal for the IWG
                    mSet->get_requested_IWGs()( iIWG )->set_normal( tNormal );

                    // compute residual at integration point
                    mSet->get_requested_IWGs()( iIWG )->compute_jacobian( tWStar );
                }
            }
        }

        //------------------------------------------------------------------------------
        void Element_Double_Sideset::compute_jacobian_and_residual()
        {
            // get treated side ordinal on the master and on the slave
            uint tMasterSideOrd = mCluster->mMasterListOfSideOrdinals( mCellIndexInCluster );
            uint tSlaveSideOrd  = mCluster->mSlaveListOfSideOrdinals( mCellIndexInCluster );

            // set the master/slave ig geometry interpolator physical/parm space and time coefficients
            this->init_ig_geometry_interpolator( tMasterSideOrd, tSlaveSideOrd );

            // get first corresponding node from master to slave
            moris::mtk::Vertex const * tSlaveNode =
                    mCluster->get_left_vertex_pair( mMasterCell->get_vertices_on_side_ordinal( tMasterSideOrd )( 0 ) );
            moris_index tSlaveNodeOrdOnSide =
                    mCluster->get_right_vertex_ordinal_on_facet( mCellIndexInCluster, tSlaveNode );

            // get rotation matrix from left to right
            Matrix< DDRMat> tR;
            rotation_matrix( mSet->get_IG_geometry_type(), tSlaveNodeOrdOnSide, tR );

            // get number of IWGs
            uint tNumIWGs = mSet->get_number_of_requested_IWGs();

            // loop over the integration points
            uint tNumIntegPoints = mSet->get_number_of_integration_points();
            for( uint iGP = 0; iGP < tNumIntegPoints; iGP++ )
            {
                // get local integration point for the master integration cell
                Matrix< DDRMat > tMasterLocalIntegPoint = mSet->get_integration_points().get_column( iGP );

                // get local integration point for the slave integration cell
                Matrix< DDRMat > tSlaveLocalIntegPoint = tMasterLocalIntegPoint;
                tSlaveLocalIntegPoint({0,tSlaveLocalIntegPoint.numel()-2},{0,0}) =
                        tR * tMasterLocalIntegPoint({0,tSlaveLocalIntegPoint.numel()-2},{0,0}); //fixme better way?

                // set evaluation point for master and slave interpolators
                mSet->get_field_interpolator_manager( mtk::Master_Slave::MASTER )->
                        set_space_time_from_local_IG_point( tMasterLocalIntegPoint );
                mSet->get_field_interpolator_manager( mtk::Master_Slave::SLAVE )->
                        set_space_time_from_local_IG_point( tSlaveLocalIntegPoint );

                // compute the integration point weight
                real tWStar = mSet->get_integration_weights()( iGP ) *
                        mSet->get_field_interpolator_manager()->get_IG_geometry_interpolator()->det_J();

                // get the normal from mesh
                Matrix< DDRMat > tNormal = mCluster->get_side_normal( mMasterCell, tMasterSideOrd );

                // loop over the IWGs
                for( uint iIWG = 0; iIWG < tNumIWGs; iIWG++ )
                {
                    // reset IWG
                    mSet->get_requested_IWGs()( iIWG )->reset_eval_flags();

                    // set the normal for the IWG
                    mSet->get_requested_IWGs()( iIWG )->set_normal( tNormal );

                    // compute residual at integration point
                    mSet->get_requested_IWGs()( iIWG )->compute_residual( tWStar );

                    // compute jacobian at integration point
                    mSet->get_requested_IWGs()( iIWG )->compute_jacobian( tWStar );
                }
            }
        }

        //------------------------------------------------------------------------------

        void Element_Double_Sideset::compute_dRdp()
        {
            // get treated side ordinal on the master and on the slave
            uint tMasterSideOrd = mCluster->mMasterListOfSideOrdinals( mCellIndexInCluster );
            uint tSlaveSideOrd  = mCluster->mSlaveListOfSideOrdinals( mCellIndexInCluster );

            // set the master/slave ig geometry interpolator physical/parm space and time coefficients
            moris::Cell< Matrix< DDSMat > > tMasterIsActiveDv;
            moris::Cell< Matrix< DDSMat > > tSlaveIsActiveDv;

            this->init_ig_geometry_interpolator_with_pdv(
                    tMasterSideOrd,
                    tSlaveSideOrd,
                    tMasterIsActiveDv,
                    tSlaveIsActiveDv );

            // get first corresponding node from master to slave
            moris::mtk::Vertex const * tSlaveNode =
                    mCluster->get_left_vertex_pair(
                            mMasterCell->get_vertices_on_side_ordinal( tMasterSideOrd )( 0 ) );

            moris_index tSlaveNodeOrdOnSide =
                    mCluster->get_right_vertex_ordinal_on_facet(mCellIndexInCluster,tSlaveNode );

            // get rotation matrix from left to right
            Matrix< DDRMat> tR;
            rotation_matrix( mSet->get_IG_geometry_type(), tSlaveNodeOrdOnSide, tR );

            // get the vertices indices
            Matrix< IndexMat > tMasterVertexIndices =
                    mMasterCell->get_vertices_ind_on_side_ordinal( tMasterSideOrd );
            Matrix< IndexMat > tSlaveVertexIndices  =
                    mSlaveCell->get_vertices_ind_on_side_ordinal( tSlaveSideOrd );

            // get number of IWGs
            uint tNumIWGs = mSet->get_number_of_requested_IWGs();

            // loop over integration points
            uint tNumIntegPoints = mSet->get_number_of_integration_points();

            for( uint iGP = 0; iGP < tNumIntegPoints; iGP++ )
            {
                // get local integration point for the master integration cell
                Matrix< DDRMat > tMasterLocalIntegPoint = mSet->get_integration_points().get_column( iGP );

                // get local integration point for the slave integration cell
                Matrix< DDRMat > tSlaveLocalIntegPoint = tMasterLocalIntegPoint;
                tSlaveLocalIntegPoint({0,tMasterLocalIntegPoint.numel()-2},{0,0}) =
                        tR * tMasterLocalIntegPoint({0,tSlaveLocalIntegPoint.numel()-2},{0,0}); //fixme better way?

                // set evaluation point for master and slave interpolators
                mSet->get_field_interpolator_manager( mtk::Master_Slave::MASTER )->
                        set_space_time_from_local_IG_point( tMasterLocalIntegPoint );

                mSet->get_field_interpolator_manager( mtk::Master_Slave::SLAVE )->
                        set_space_time_from_local_IG_point( tSlaveLocalIntegPoint );

                // compute the integration point weight
                real tWStar = mSet->get_integration_weights()( iGP ) *
                        mSet->get_field_interpolator_manager()->get_IG_geometry_interpolator()->det_J();

                // get the normal from mesh and set if for the IWG
                Matrix< DDRMat > tNormal = mCluster->get_side_normal( mMasterCell, tMasterSideOrd );

                // loop over the IWGs
                for( uint iIWG = 0; iIWG < tNumIWGs; iIWG++ )
                {
                    // reset IWG
                    mSet->get_requested_IWGs()( iIWG )->reset_eval_flags();

                    // set a perturbation size
                    real tPerturbation = 1E-6;

                    // compute dRdpMat at evaluation point
                    mSet->get_requested_IWGs()( iIWG )->compute_dRdp_FD_material_double(
                            tWStar,
                            tPerturbation );

                    // compute dRdpGeo at evaluation point
                    mSet->get_requested_IWGs()( iIWG )->compute_dRdp_FD_geometry_double(
                            tWStar,
                            tPerturbation,
                            tMasterIsActiveDv,
                            tMasterVertexIndices,
                            tSlaveIsActiveDv,
                            tSlaveVertexIndices );
                }
            }
        }

        //------------------------------------------------------------------------------

        real Element_Double_Sideset::compute_volume( mtk::Master_Slave aIsMaster )
        {
            // get treated side ordinal on the master and on the slave
            uint tMasterSideOrd = mCluster->mMasterListOfSideOrdinals( mCellIndexInCluster );
            uint tSlaveSideOrd  = mCluster->mSlaveListOfSideOrdinals( mCellIndexInCluster );

            // set the master/slave ig geometry interpolator physical/parm space and time coefficients
            this->init_ig_geometry_interpolator( tMasterSideOrd, tSlaveSideOrd );

            //get number of integration points
            uint tNumOfIntegPoints = mSet->get_number_of_integration_points();

            // initialize volume
            real tVolume = 0;

            // get geometry interpolator
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
