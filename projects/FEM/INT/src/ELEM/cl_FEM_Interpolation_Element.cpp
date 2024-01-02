/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_FEM_Interpolation_Element.cpp
 *
 */

#include <iostream>
// FEM/INT/src
#include "cl_FEM_Element.hpp"
#include "cl_FEM_Interpolation_Element.hpp"
#include "cl_FEM_Field_Interpolator_Manager.hpp"
#include "cl_MSI_Design_Variable_Interface.hpp"
#include "cl_FEM_Cluster.hpp"
#include "cl_FEM_Set.hpp"
#include "cl_FEM_Model.hpp"
#include "cl_FEM_Field.hpp"
// SOL/src
#include "cl_SOL_Dist_Vector.hpp"
// LINALG/src
#include "fn_isfinite.hpp"

namespace moris
{
    namespace fem
    {
        //------------------------------------------------------------------------------

        Interpolation_Element::Interpolation_Element(
                const Element_Type               aElementType,
                const Vector< const mtk::Cell* >&  aInterpolationCell,
                const Vector< Node_Base* >& aNodes,
                Set*                             aSet )
                : MSI::Equation_Object( aSet )
                , mSet( aSet )
                , mElementType( aElementType )
        {
            // fill the leader interpolation cell
            mLeaderInterpolationCell = aInterpolationCell( 0 );

            // get vertices from cell
            Vector< mtk::Vertex* > tVertices =
                    mLeaderInterpolationCell->get_vertex_pointers();

            // get number of vertices from cell
            uint tNumOfVertices = tVertices.size();

            // assign node object
            mNodeObj.resize( 1 );
            mNodeObj( 0 ).resize( tNumOfVertices, nullptr );

            // fill leader node objects
            for ( uint iVertex = 0; iVertex < tNumOfVertices; iVertex++ )
            {
                mNodeObj( 0 )( iVertex ) = aNodes( tVertices( iVertex )->get_index() );
            }

            // if double sided sideset
            if ( mElementType == fem::Element_Type::DOUBLE_SIDESET )
            {
                // fill the follower interpolation cell
                mFollowerInterpolationCell = aInterpolationCell( 1 );

                // get vertices from cell
                Vector< mtk::Vertex* > tFollowerVertices =
                        mFollowerInterpolationCell->get_vertex_pointers();

                // get number of vertices from cell
                uint tNumOfFollowerVertices = tFollowerVertices.size();

                // assign node object
                mNodeObj.resize( 2 );
                mNodeObj( 1 ).resize( tNumOfFollowerVertices, nullptr );

                // fill follower node objects
                for ( uint iVertex = 0; iVertex < tNumOfFollowerVertices; iVertex++ )
                {
                    mNodeObj( 1 )( iVertex ) = aNodes( tFollowerVertices( iVertex )->get_index() );
                }
            }
        }

        //------------------------------------------------------------------------------

        void
        Interpolation_Element::set_cluster(
                std::shared_ptr< fem::Cluster > aCluster,
                const uint                      aFemMeshIndex )
        {
            // if mesh index is 0 (i.e., forward analysis mesh, IG mesh)
            if ( aFemMeshIndex == 0 )
            {
                // fem cluster with index 0 should be set only once and shall not be changed
                MORIS_ASSERT( !( mFemCluster.size() >= 1 ),
                        "Interpolation_Element::set_cluster() - first fem cluster is already set" );
            }

            // get max size for fem cluster list
            sint tSize = std::max( (sint)mFemCluster.size(), (sint)aFemMeshIndex + 1 );

            // resize fem cluster list
            mFemCluster.resize( tSize );

            // add the fem cluster to the list
            mFemCluster( aFemMeshIndex ) = aCluster;
        }

        //------------------------------------------------------------------------------

        const std::shared_ptr< fem::Cluster >&
        Interpolation_Element::get_cluster( const uint aIndex )
        {
            MORIS_ERROR( aIndex < mFemCluster.size(),
                    "Interpolation_Element::get_cluster - index out of bounds.\n" );

            return mFemCluster( aIndex );
        }

        //------------------------------------------------------------------------------

        void
        Interpolation_Element::set_field_interpolators_coefficients()
        {
            // dof field interpolators------------------------------------------

            // get leader dof type list from set
            Vector< Vector< MSI::Dof_Type > >& tLeaderDofTypeList =
                    mSet->get_dof_type_list( mtk::Leader_Follower::LEADER );

            // get number of leader dof types
            uint tLeaderNumDofTypes = tLeaderDofTypeList.size();

            // loop on the leader dof types
            for ( uint iDOF = 0; iDOF < tLeaderNumDofTypes; iDOF++ )
            {
                // get the ith dof type group
                Vector< MSI::Dof_Type >& tDofTypeGroup = tLeaderDofTypeList( iDOF );

                // get the pdof values for the ith dof type group
                Vector< Vector< Matrix< DDRMat > > > tCoeff_Original;
                this->get_my_pdof_values( mSet->mPdofValues, tDofTypeGroup, tCoeff_Original );

                // reshape tCoeffs into the order the cluster expects them
                Matrix< DDRMat > tCoeff;
                this->reshape_pdof_values( tCoeff_Original( 0 ), tCoeff );

                // set field interpolator coefficients
                mSet->get_field_interpolator_manager()->set_coeff_for_type( tDofTypeGroup( 0 ), tCoeff );

                // if previous solution
                if ( mSet->get_time_continuity() )
                {
                    // get the pdof values for the ith dof type group
                    Vector< Vector< Matrix< DDRMat > > > tCoeff_Original;
                    this->get_my_pdof_values( mSet->mPreviousPdofValues, tDofTypeGroup, tCoeff_Original );

                    // reshape tCoeffs into the order the cluster expects them
                    Matrix< DDRMat > tCoeff;
                    this->reshape_pdof_values( tCoeff_Original( 0 ), tCoeff );

                    // set field interpolator coefficients
                    mSet->get_field_interpolator_manager_previous_time()->set_coeff_for_type( tDofTypeGroup( 0 ), tCoeff );
                }

                // if eigen vectors
                if ( mSet->mEigenVectorPdofValues.size() > 0 )
                {
                    // get the pdof values for the ith dof type group
                    Vector< Vector< Matrix< DDRMat > > > tCoeff_Original;
                    this->get_my_pdof_values( mSet->mEigenVectorPdofValues, tDofTypeGroup, tCoeff_Original );

                    // check for consistency of number of eigen vectors
                    MORIS_ASSERT( tCoeff_Original.size() == mSet->mNumEigenVectors,
                            "Interpolation_Element::set_field_interpolators_coefficients - inconsistent number of eigen vectors" );

                    // loop over all eigen vectors
                    for ( uint iv = 0; iv < tCoeff_Original.size(); ++iv )
                    {
                        // reshape tCoeffs into the order the cluster expects them
                        Matrix< DDRMat > tCoeff;
                        this->reshape_pdof_values( tCoeff_Original( iv ), tCoeff );

                        // set field interpolator coefficients
                        mSet->get_field_interpolator_manager_eigen_vectors()->set_coeff_for_type( tDofTypeGroup( 0 ), tCoeff, iv );
                    }
                }
            }

            // get follower dof type list from set
            Vector< Vector< MSI::Dof_Type > >& tFollowerDofTypeList =
                    mSet->get_dof_type_list( mtk::Leader_Follower::FOLLOWER );

            // get number of follower dof types
            uint tFollowerNumDofTypes = tFollowerDofTypeList.size();

            // loop on the follower dof types
            for ( uint iDOF = 0; iDOF < tFollowerNumDofTypes; iDOF++ )
            {
                // get the ith dof type group
                Vector< MSI::Dof_Type >& tDofTypeGroup = tFollowerDofTypeList( iDOF );

                // get the pdof values for the ith dof type group
                Vector< Vector< Matrix< DDRMat > > > tCoeff_Original;
                this->get_my_pdof_values( mSet->mPdofValues, tDofTypeGroup, tCoeff_Original, mtk::Leader_Follower::FOLLOWER );

                // reshape tCoeffs into the order the cluster expects them
                Matrix< DDRMat > tCoeff;
                this->reshape_pdof_values( tCoeff_Original( 0 ), tCoeff );

                // set the field coefficients
                mSet->get_field_interpolator_manager( mtk::Leader_Follower::FOLLOWER )->set_coeff_for_type( tDofTypeGroup( 0 ), tCoeff );
            }

            // dv field interpolators------------------------------------------

            // get leader dv type list from set
            const Vector< Vector< PDV_Type > >& tLeaderDvTypeList =
                    mSet->get_dv_type_list( mtk::Leader_Follower::LEADER );

            // get number of leader dv types
            uint tLeaderNumDvTypes = tLeaderDvTypeList.size();

            // loop on the leader dv types
            for ( uint iDv = 0; iDv < tLeaderNumDvTypes; iDv++ )
            {
                // get the dv type group
                const Vector< PDV_Type >& tDvTypeGroup = tLeaderDvTypeList( iDv );

                // get the pdv values for the ith dv type group
                // FIXME: the underlying use of the base cell needs to be hidden within PDV
                Vector< Matrix< DDRMat > > tCoeff_Original;
                mSet->get_equation_model()->get_design_variable_interface()->get_ip_pdv_value(
                        mLeaderInterpolationCell->get_base_cell()->get_vertex_inds(),
                        tDvTypeGroup,
                        tCoeff_Original );

                // reshape tCoeffs into the order the FI expects them
                Matrix< DDRMat > tCoeff;
                mSet->get_equation_model()->get_design_variable_interface()->reshape_pdv_values( tCoeff_Original, tCoeff );

                // set field interpolator coefficients
                mSet->get_field_interpolator_manager()->set_coeff_for_type( tDvTypeGroup( 0 ), tCoeff );
            }

            // get follower dv type list from set
            const Vector< Vector< PDV_Type > >& tFollowerDvTypeList =
                    mSet->get_dv_type_list( mtk::Leader_Follower::FOLLOWER );

            // get number of follower dv types
            uint tFollowerNumDvTypes = tFollowerDvTypeList.size();

            // loop on the follower dv types
            for ( uint iDv = 0; iDv < tFollowerNumDvTypes; iDv++ )
            {
                // get the dv type group
                const Vector< PDV_Type >& tDvTypeGroup = tFollowerDvTypeList( iDv );

                // get the pdv values for the ith dv type group
                Vector< Matrix< DDRMat > > tCoeff_Original;
                mSet->get_equation_model()->get_design_variable_interface()->get_ip_pdv_value(
                        mFollowerInterpolationCell->get_base_cell()->get_vertex_inds(),
                        tDvTypeGroup,
                        tCoeff_Original );

                // reshape tCoeffs into the order the FI expects them
                Matrix< DDRMat > tCoeff;
                mSet->get_equation_model()->get_design_variable_interface()->reshape_pdv_values( tCoeff_Original, tCoeff );

                // set the field coefficients
                mSet->get_field_interpolator_manager( mtk::Leader_Follower::FOLLOWER )->set_coeff_for_type( tDvTypeGroup( 0 ), tCoeff );
            }

            // field field interpolators------------------------------------------

            // get leader field type list from set
            const Vector< Vector< mtk::Field_Type > >& tLeaderFieldTypeList =
                    mSet->get_field_type_list( mtk::Leader_Follower::LEADER );

            // get number of leader field types
            uint tLeaderNumFieldTypes = tLeaderFieldTypeList.size();

            // loop on the leader field types
            for ( uint iFi = 0; iFi < tLeaderNumFieldTypes; iFi++ )
            {
                // get the field type group
                const Vector< mtk::Field_Type >& tFieldTypeGroup = tLeaderFieldTypeList( iFi );

                Matrix< IndexMat > tIPCellIndices = mLeaderInterpolationCell->get_vertex_inds();

                Matrix< DDRMat > tCoeff;
                mSet->get_fem_model()->get_field( tFieldTypeGroup( 0 ) )->get_values( tIPCellIndices, tCoeff, tFieldTypeGroup );

                // FIXME implement reshape for vector fields

                // set field interpolator coefficients
                mSet->get_field_interpolator_manager()->set_coeff_for_type( tFieldTypeGroup( 0 ), tCoeff );
            }

            // get follower field type list from set
            const Vector< Vector< mtk::Field_Type > >& tFollowerFieldTypeList =
                    mSet->get_field_type_list( mtk::Leader_Follower::FOLLOWER );

            // get number of follower field types
            uint tFollowerNumFieldTypes = tFollowerFieldTypeList.size();

            // loop on the follower field types
            for ( uint iFi = 0; iFi < tFollowerNumFieldTypes; iFi++ )
            {
                // get the field type group
                const Vector< mtk::Field_Type >& tFieldTypeGroup = tFollowerFieldTypeList( iFi );

                Matrix< IndexMat > tIPCellIndices = mFollowerInterpolationCell->get_vertex_inds();

                Matrix< DDRMat > tCoeff;
                mSet->get_fem_model()->get_field( tFieldTypeGroup( 0 ) )->get_values( tIPCellIndices, tCoeff, tFieldTypeGroup );

                // FIXME implement reshape for vector fields

                // set field interpolator coefficients
                mSet->get_field_interpolator_manager( mtk::Leader_Follower::FOLLOWER )->set_coeff_for_type( tFieldTypeGroup( 0 ), tCoeff );
            }

            // geometry interpolators------------------------------------------
            // set the IP geometry interpolator physical space and time coefficients for the leader
            mSet->get_field_interpolator_manager( mtk::Leader_Follower::LEADER )->get_IP_geometry_interpolator()->set_space_coeff( mLeaderInterpolationCell->get_vertex_coords() );
            mSet->get_field_interpolator_manager( mtk::Leader_Follower::LEADER )->get_IP_geometry_interpolator()->set_param_coeff();
            mSet->get_field_interpolator_manager( mtk::Leader_Follower::LEADER )->get_IP_geometry_interpolator()->set_time_coeff( this->get_time() );

            // if double sideset
            if ( mElementType == fem::Element_Type::DOUBLE_SIDESET )
            {
                // set the IP geometry interpolator physical space and time coefficients for the follower
                mSet->get_field_interpolator_manager( mtk::Leader_Follower::FOLLOWER )->get_IP_geometry_interpolator()->set_space_coeff( mFollowerInterpolationCell->get_vertex_coords() );
                mSet->get_field_interpolator_manager( mtk::Leader_Follower::FOLLOWER )->get_IP_geometry_interpolator()->set_time_coeff( this->get_time() );
            }

            // if time sideset
            if ( mElementType == fem::Element_Type::TIME_SIDESET )
            {
                // set the IP geometry interpolator physical space and time coefficients for the previous
                mSet->get_field_interpolator_manager_previous_time( mtk::Leader_Follower::LEADER )->get_IP_geometry_interpolator()->set_space_coeff( mLeaderInterpolationCell->get_vertex_coords() );
                mSet->get_field_interpolator_manager_previous_time( mtk::Leader_Follower::LEADER )->get_IP_geometry_interpolator()->set_time_coeff( this->get_previous_time() );
            }

            // if eigen vectors
            if ( mSet->mEigenVectorPdofValues.size() > 0 )
            {
                // set the IP geometry interpolator physical space and time coefficients for eigen vectors
                mSet->get_field_interpolator_manager_eigen_vectors( mtk::Leader_Follower::LEADER )->get_IP_geometry_interpolator()->set_space_coeff( mLeaderInterpolationCell->get_vertex_coords() );
                mSet->get_field_interpolator_manager_eigen_vectors( mtk::Leader_Follower::LEADER )->get_IP_geometry_interpolator()->set_time_coeff( this->get_time() );
            }
        }

        //------------------------------------------------------------------------------

        void
        Interpolation_Element::fill_mat_pdv_assembly_vector()
        {
            // get the design variable interface
            MSI::Design_Variable_Interface* tDVInterface =
                    mSet->mEquationModel->get_design_variable_interface();

            // get the list of requested dv types by the opt solver
            Vector< Vector< enum PDV_Type > > tRequestedDvTypes;
            mSet->get_ip_dv_types_for_set( tRequestedDvTypes );

            // reset material pdv assembly vector
            mSet->get_mat_pdv_assembly_vector().fill( -1 );

            // init pdv counter
            uint tCounter = 0;

            // get leader vertices from cell
            Matrix< IndexMat > tLeaderVerticesInds =
                    mLeaderInterpolationCell->get_base_cell()->get_vertex_inds();

            // loop over the dv types
            for ( uint Ik = 0; Ik < tRequestedDvTypes.size(); Ik++ )
            {
                // get dv ids for this type and node indices
                Vector< Matrix< IdMat > > tPdvIds;

                // get the pdv ids for requested vertices and pdv type
                tDVInterface->get_ip_dv_ids_for_type_and_ind(
                        tLeaderVerticesInds,
                        tRequestedDvTypes( Ik ),
                        tPdvIds );

                // get number of coefficients
                uint tNumCoeff = tPdvIds( 0 ).numel();

                // check that coefficients are larger or equal 0
                MORIS_ASSERT( tPdvIds( 0 ).min() > -1,
                        "Interpolation_Element::fill_mat_pdv_assembly_vector - %s",
                        "Coefficient vector includes negative numbers.\n" );

                // fill the assembly vector with pdv ids
                mSet->get_mat_pdv_assembly_vector()(
                        { tCounter, tCounter + tNumCoeff - 1 } ) = tPdvIds( 0 ).matrix_data();

                // update pdv counter
                tCounter += tNumCoeff;
            }

            // double sided
            if ( mElementType == fem::Element_Type::DOUBLE_SIDESET )
            {
                // get follower vertices from cell
                Matrix< IndexMat > tFollowerVerticesInds =
                        mFollowerInterpolationCell->get_base_cell()->get_vertex_inds();

                // get the list of requested dv types by the opt solver for the follower side
                mSet->get_ip_dv_types_for_set( tRequestedDvTypes, mtk::Leader_Follower::FOLLOWER );

                // loop over the dv types
                for ( uint Ik = 0; Ik < tRequestedDvTypes.size(); Ik++ )
                {
                    // get dv ids for this type and node indices
                    Vector< Matrix< IdMat > > tPdvIds;

                    // get the pdv ids for requested vertices and pdv type
                    tDVInterface->get_ip_dv_ids_for_type_and_ind(
                            tFollowerVerticesInds,
                            tRequestedDvTypes( Ik ),
                            tPdvIds );

                    // get number of coefficients
                    uint tNumCoeff = tPdvIds( 0 ).numel();

                    // check that coefficients are larger or equal 0
                    MORIS_ASSERT( tPdvIds( 0 ).min() > -1,
                            "Interpolation_Element::fill_mat_pdv_assembly_vector - %s",
                            "Coefficient vector includes negative numbers.\n" );

                    // fill the assembly vector with pdv ids
                    mSet->get_mat_pdv_assembly_vector()(
                            { tCounter, tCounter + tNumCoeff - 1 } ) = tPdvIds( 0 ).matrix_data();

                    // update pdv counter
                    tCounter += tNumCoeff;
                }
            }
        }

        //------------------------------------------------------------------------------

        void
        Interpolation_Element::compute_jacobian()
        {
            // compute pdof values
            // FIXME do this only once
            this->compute_my_pdof_values();

            // if time continuity set
            if ( mSet->get_time_continuity() )
            {
                // compute pdof values for previous time step
                // FIXME do this only once
                this->compute_previous_pdof_values();
            }

            // initialize the residual
            mSet->initialize_mResidual();

            // initialize the jacobian
            mSet->initialize_mJacobian();

            // set the field interpolators coefficients
            this->set_field_interpolators_coefficients();

            // FIXME should not be like this
            mSet->set_IWG_field_interpolator_managers();

            // set cluster for stabilization parameter
            mSet->set_IWG_cluster_for_stabilization_parameters( mFemCluster( 0 ).get() );

            // ask cluster to compute jacobian
            mFemCluster( 0 )->compute_jacobian();
        }

        //------------------------------------------------------------------------------

        void
        Interpolation_Element::compute_residual()
        {
            // Fixme do this only once
            this->compute_my_pdof_values();

            // if time continuity set
            if ( mSet->get_time_continuity() )
            {
                // compute pdof values for previous time step
                // FIXME do this only once
                this->compute_previous_pdof_values();
            }

            // initialize the residual
            mSet->initialize_mResidual();

            // initialize the jacobian
            mSet->initialize_mJacobian();

            // set the field interpolators coefficients
            this->set_field_interpolators_coefficients();

            // FIXME should not be like this
            mSet->set_IWG_field_interpolator_managers();

            // set cluster for stabilization parameter
            mSet->set_IWG_cluster_for_stabilization_parameters( mFemCluster( 0 ).get() );

            if ( mSet->mEquationModel->get_is_forward_analysis() )
            {
                // FIXME should not be like this
                mSet->set_IWG_field_interpolator_managers();

                // set cluster for stabilization parameter
                mSet->set_IWG_cluster_for_stabilization_parameters( mFemCluster( 0 ).get() );

                // ask cluster to compute residual
                mFemCluster( 0 )->compute_residual();
            }
            else if ( ( !mSet->mEquationModel->get_is_forward_analysis() ) && ( mSet->get_number_of_requested_IQIs() > 0 ) )
            {
                // FIXME should not be like this
                mSet->set_IQI_field_interpolator_managers();

                // set cluster for stabilization parameter
                mSet->set_IQI_cluster_for_stabilization_parameters( mFemCluster( 0 ).get() );

                // ask cluster to compute jacobian
                mFemCluster( 0 )->compute_dQIdu();
            }
        }

        //------------------------------------------------------------------------------

        void
        Interpolation_Element::compute_jacobian_and_residual()
        {
            // Fixme do this only once
            this->compute_my_pdof_values();

            // if time continuity set
            if ( mSet->get_time_continuity() )
            {
                // compute pdof values for previous time step
                // FIXME do this only once
                this->compute_previous_pdof_values();
            }

            // initialize the Jacobian
            mSet->initialize_mJacobian();

            // initialize the residual
            mSet->initialize_mResidual();

            // set the field interpolators coefficients
            this->set_field_interpolators_coefficients();

            // FIXME should not be like this
            mSet->set_IWG_field_interpolator_managers();

            // set cluster for stabilization parameter
            mSet->set_IWG_cluster_for_stabilization_parameters( mFemCluster( 0 ).get() );

            if ( ( !mSet->mEquationModel->get_is_forward_analysis() ) && ( mSet->get_number_of_requested_IQIs() > 0 ) )
            {
                // FIXME should not be like this
                mSet->set_IQI_field_interpolator_managers();

                // set cluster for stabilization parameter
                mSet->set_IQI_cluster_for_stabilization_parameters( mFemCluster( 0 ).get() );
            }

            // ask cluster to compute Jacobian and residual
            mFemCluster( 0 )->compute_jacobian_and_residual();
        }

        //------------------------------------------------------------------------------

        void
        Interpolation_Element::compute_dRdp()
        {
            // compute pdof values
            // FIXME do this only once
            this->compute_my_pdof_values();

            // if time continuity set
            if ( mSet->get_time_continuity() )
            {
                // compute pdof values for previous time step
                // FIXME do this only once
                this->compute_previous_pdof_values();
            }

            // init geo pdv assembly map
            mSet->create_geo_pdv_assembly_map( mFemCluster( 0 ) );

            // init dRdp
            mSet->initialize_mdRdpMat();
            mSet->initialize_mdRdpGeo( mFemCluster( 0 ) );

            // as long as dRdp is computed with FD we need this
            mSet->initialize_mResidual();

            // set the field interpolators coefficients
            this->set_field_interpolators_coefficients();

            // FIXME should not be like this
            mSet->set_IWG_field_interpolator_managers();

            // set cluster for stabilization parameter
            mSet->set_IWG_cluster_for_stabilization_parameters( mFemCluster( 0 ).get() );

            // ask cluster to compute jacobian
            mFemCluster( 0 )->compute_dRdp();
        }

        //------------------------------------------------------------------------------

        void
        Interpolation_Element::compute_dQIdp_explicit()
        {
            // initialize IP and IG pdv assembly maps
            this->fill_mat_pdv_assembly_vector();
            mSet->create_geo_pdv_assembly_map( mFemCluster( 0 ) );

            // get the assembly vector
            const Matrix< DDSMat >& tLocalToGlobalIdsIPPdv =
                    mEquationSet->get_mat_pdv_assembly_vector();

            // get the assembly vector
            const Matrix< DDSMat >& tLocalToGlobalIdsIGPdv =
                    mEquationSet->get_geo_pdv_assembly_vector();

            // if there is no pdv defined, return
            if ( tLocalToGlobalIdsIPPdv.numel() == 0 && tLocalToGlobalIdsIGPdv.numel() == 0 )
            {
                return;
            }

            // compute pdof values
            // FIXME do this only once
            this->compute_my_pdof_values();

            // if time continuity set
            if ( mSet->get_time_continuity() )
            {
                // compute pdof values for previous time step
                // FIXME do this only once
                this->compute_previous_pdof_values();
            }

            // set the field interpolators coefficients
            this->set_field_interpolators_coefficients();

            // FIXME should not be like this
            mSet->set_IQI_field_interpolator_managers();

            // set cluster for stabilization parameter
            mSet->set_IQI_cluster_for_stabilization_parameters( mFemCluster( 0 ).get() );

            // initialize dQIdp
            mSet->initialize_mdQIdpMat();
            mSet->initialize_mdQIdpGeo( mFemCluster( 0 ) );

            // ask cluster to compute jacobian
            mFemCluster( 0 )->compute_dQIdp_explicit();

            // Assembly for the IP pdv
            //----------------------------------------------------------------------------------------
            // if assembly vector is not empty
            if ( tLocalToGlobalIdsIPPdv.numel() != 0 )
            {
                // loop over the IP pdv
                uint tNumIQIs = mSet->mdQIdp( 0 ).size();
                for ( uint Ik = 0; Ik < tNumIQIs; Ik++ )
                {
                    // assemble explicit dQIdpMat into multivector
                    mEquationSet->get_equation_model()->get_explicit_dQidp()->sum_into_global_values(
                            tLocalToGlobalIdsIPPdv,
                            mSet->mdQIdp( 0 )( Ik ),
                            Ik );
                }
            }

            // Assembly for the IG pdv
            //----------------------------------------------------------------------------------------
            // if assembly vector is not empty
            if ( tLocalToGlobalIdsIGPdv.numel() != 0 )
            {
                // loop over the IG pdv
                uint tNumIQIs = mSet->mdQIdp( 1 ).size();
                for ( uint Ik = 0; Ik < tNumIQIs; Ik++ )
                {
                    // assemble explicit dQIdpGeo into multivector
                    mEquationSet->get_equation_model()->get_explicit_dQidp()->sum_into_global_values(
                            tLocalToGlobalIdsIGPdv,
                            mSet->mdQIdp( 1 )( Ik ),
                            Ik );
                }
            }
        }

        //------------------------------------------------------------------------------

        void
        Interpolation_Element::compute_dQIdp_explicit_implicit()
        {
            // fill IP pdv assembly vector
            this->fill_mat_pdv_assembly_vector();

            // make sure this is done on the FEM clusters and not on the VIS clusters
            MORIS_ASSERT( !mFemCluster( 0 )->is_VIS_cluster(),
                    "FEM::Set::compute_dQIdp_explicit_implicit() - "
                    "Trying to compute sensitivities on a VIS cluster. This shouldn't happen." );

            // init IG pdv assembly map
            mSet->create_geo_pdv_assembly_map( mFemCluster( 0 ) );

            // get the IP pdv assembly vector
            const Matrix< DDSMat >& tLocalToGlobalIdsIPPdv =
                    mEquationSet->get_mat_pdv_assembly_vector();

            // get the IG pdv assembly vector
            const Matrix< DDSMat >& tLocalToGlobalIdsIGPdv =
                    mEquationSet->get_geo_pdv_assembly_vector();

            // if there is no pdv defined, return
            if ( tLocalToGlobalIdsIPPdv.numel() == 0 && tLocalToGlobalIdsIGPdv.numel() == 0 )
            {
                return;
            }

            // init dRdp
            mSet->initialize_mdRdpMat();
            mSet->initialize_mdRdpGeo( mFemCluster( 0 ) );

            // initialize dQIdp
            mSet->initialize_mdQIdpMat();
            mSet->initialize_mdQIdpGeo( mFemCluster( 0 ) );

            // as long as dRdp is computed with FD,
            // we need to init the residual storage
            mSet->initialize_mResidual();

            // compute pdof values
            // FIXME do this only once
            this->compute_my_pdof_values();

            // if time continuity set
            if ( mSet->get_time_continuity() )
            {
                // compute pdof values for previous time step
                // FIXME do this only once
                this->compute_previous_pdof_values();
            }

            // set the field interpolators coefficients
            this->set_field_interpolators_coefficients();

            // FIXME should not be like this
            mSet->set_IWG_field_interpolator_managers();

            // set cluster for stabilization parameter
            mSet->set_IWG_cluster_for_stabilization_parameters( mFemCluster( 0 ).get() );

            // FIXME should not be like this
            mSet->set_IQI_field_interpolator_managers();

            // set cluster for stabilization parameter
            mSet->set_IQI_cluster_for_stabilization_parameters( mFemCluster( 0 ).get() );

            // ask cluster to compute jacobian
            mFemCluster( 0 )->compute_dRdp_and_dQIdp();

            // get reference to computed dRdp
            Vector< Matrix< DDRMat > >& tdRdp = mEquationSet->get_drdp();

            // extract adjoint values for this equation object
            this->compute_my_adjoint_values();

            // get number of RHS
            uint tNumRHS = mSet->mAdjointPdofValues.size();

            // get number of pdof values
            uint tNumPdofValues = tdRdp( 0 ).n_rows();

            // reorder adjoint values following the requested dof types order
            Matrix< DDRMat > tAdjointPdofValuesReordered;

            // get leader dof type list from set
            const Vector< Vector< MSI::Dof_Type > >& tLeaderDofTypeGroup =
                    mSet->get_dof_type_list( mtk::Leader_Follower::LEADER );

            // get number of leader dof types
            uint tNumLeaderDofTypes = tLeaderDofTypeGroup.size();

            // get follower dof type list from set
            const Vector< Vector< MSI::Dof_Type > >& tFollowerDofTypeGroup =
                    mSet->get_dof_type_list( mtk::Leader_Follower::FOLLOWER );

            // get number of follower dof types
            uint tNumFollowerDofTypes = tFollowerDofTypeGroup.size();

            // loop over the RHS
            for ( uint Ik = 0; Ik < tNumRHS; Ik++ )
            {
                // set size for reordered adjoint values
                tAdjointPdofValuesReordered.set_size( tNumPdofValues, 1, 0.0 );

                // loop over the leader dof types
                for ( uint Ia = 0; Ia < tNumLeaderDofTypes; Ia++ )
                {
                    // get the adjoint values for the ith dof type group
                    Vector< Vector< Matrix< DDRMat > > > tLeaderAdjointOriginal;
                    this->get_my_pdof_values(
                            mSet->mAdjointPdofValues,
                            tLeaderDofTypeGroup( Ia ),
                            tLeaderAdjointOriginal,
                            mtk::Leader_Follower::LEADER );

                    // reshape adjoint values
                    Matrix< DDRMat > tLeaderAdjointCoeff;
                    this->reshape_pdof_values_vector( tLeaderAdjointOriginal( Ik ), tLeaderAdjointCoeff );

                    // get indices for begin and end
                    uint tDofIndex   = mSet->get_dof_index_for_type( tLeaderDofTypeGroup( Ia )( 0 ), mtk::Leader_Follower::LEADER );
                    uint tStartIndex = mSet->get_res_dof_assembly_map()( tDofIndex )( 0, 0 );
                    uint tStopIndex  = mSet->get_res_dof_assembly_map()( tDofIndex )( 0, 1 );

                    // fill reordered adjoint pdof values
                    tAdjointPdofValuesReordered( { tStartIndex, tStopIndex } ) =
                            tLeaderAdjointCoeff.matrix_data();
                }

                // loop over the follower dof types
                for ( uint Ia = 0; Ia < tNumFollowerDofTypes; Ia++ )
                {
                    // get the adjoint values for the ith dof type group
                    Vector< Vector< Matrix< DDRMat > > > tFollowerAdjointOriginal;
                    this->get_my_pdof_values(
                            mSet->mAdjointPdofValues,
                            tFollowerDofTypeGroup( Ia ),
                            tFollowerAdjointOriginal,
                            mtk::Leader_Follower::FOLLOWER );

                    // reshape adjoint values
                    Matrix< DDRMat > tFollowerAdjointCoeff;
                    this->reshape_pdof_values_vector( tFollowerAdjointOriginal( Ik ), tFollowerAdjointCoeff );

                    // get indices for begin and end
                    uint tDofIndex   = mSet->get_dof_index_for_type( tFollowerDofTypeGroup( Ia )( 0 ), mtk::Leader_Follower::FOLLOWER );
                    uint tStartIndex = mSet->get_res_dof_assembly_map()( tDofIndex )( 0, 0 );
                    uint tStopIndex  = mSet->get_res_dof_assembly_map()( tDofIndex )( 0, 1 );

                    // fill reordered adjoint pdof values
                    tAdjointPdofValuesReordered( { tStartIndex, tStopIndex } ) =
                            tFollowerAdjointCoeff.matrix_data();
                }

                // Assembly for the IP pdv
                //----------------------------------------------------------------------------------------
                // if the assembly vector is not empty
                if ( tLocalToGlobalIdsIPPdv.numel() != 0 )
                {
                    // assemble explicit dQIdpMat into multivector
                    mEquationSet->get_equation_model()->get_explicit_dQidp()->sum_into_global_values(
                            tLocalToGlobalIdsIPPdv,
                            mSet->mdQIdp( 0 )( Ik ),
                            Ik );

                    // post multiplication of adjoint values time dRdp
                    Matrix< DDRMat > tLocalIPdQiDp =
                            -1.0 * trans( tAdjointPdofValuesReordered ) * tdRdp( 0 );

                    // assemble implicit dQidp into multivector
                    mEquationSet->get_equation_model()->get_implicit_dQidp()->sum_into_global_values(
                            tLocalToGlobalIdsIPPdv,
                            tLocalIPdQiDp,
                            Ik );
                }

                // Assembly for the IG pdv
                //----------------------------------------------------------------------------------------
                // if assembly vector is not empty
                if ( tLocalToGlobalIdsIGPdv.numel() != 0 )
                {
                    // assemble explicit dQIdpGeo into multivector
                    mEquationSet->get_equation_model()->get_explicit_dQidp()->sum_into_global_values(
                            tLocalToGlobalIdsIGPdv,
                            mSet->mdQIdp( 1 )( Ik ),
                            Ik );

                    // post multiplication of adjoint values time dRdp
                    Matrix< DDRMat > tLocalIGdQiDp =
                            -1.0 * trans( tAdjointPdofValuesReordered ) * tdRdp( 1 );

                    // assemble implicit dQidp into multivector
                    mEquationSet->get_equation_model()->get_implicit_dQidp()->sum_into_global_values(
                            tLocalToGlobalIdsIGPdv,
                            tLocalIGdQiDp,
                            Ik );
                }
            }
        }

        //-------------------------------------------------------------------------------------------------

        void
        Interpolation_Element::compute_dQIdp_implicit()
        {
            // fill IP pdv assembly vector
            this->fill_mat_pdv_assembly_vector();

            // init IG pdv assembly map
            mSet->create_geo_pdv_assembly_map( mFemCluster( 0 ) );

            // get the IP pdv assembly vector
            const Matrix< DDSMat >& tLocalToGlobalIdsIPPdv =
                    mEquationSet->get_mat_pdv_assembly_vector();

            // get the IG pdv assembly vector
            const Matrix< DDSMat >& tLocalToGlobalIdsIGPdv =
                    mEquationSet->get_geo_pdv_assembly_vector();

            // if there is no pdv defined, return
            if ( tLocalToGlobalIdsIPPdv.numel() == 0 && tLocalToGlobalIdsIGPdv.numel() == 0 )
            {
                return;
            }

            // init dRdp
            mSet->initialize_mdRdpMat();
            mSet->initialize_mdRdpGeo( mFemCluster( 0 ) );

            // as long as dRdp is computed with FD,
            // we need to init the residual storage
            mSet->initialize_mResidual();

            // compute pdof values
            // FIXME do this only once
            this->compute_my_pdof_values();

            // if time continuity set
            if ( mSet->get_time_continuity() )
            {
                // compute pdof values for previous time step
                // FIXME do this only once
                this->compute_previous_pdof_values();
            }

            // set the field interpolators coefficients
            this->set_field_interpolators_coefficients();

            // FIXME should not be like this
            mSet->set_IWG_field_interpolator_managers();

            // set cluster for stabilization parameter
            mSet->set_IWG_cluster_for_stabilization_parameters( mFemCluster( 0 ).get() );

            // ask cluster to compute dRdp
            mFemCluster( 0 )->compute_dRdp();

            // get reference to computed dRdp
            Vector< Matrix< DDRMat > >& tdRdp = mEquationSet->get_drdp();

            // extract adjoint values for this equation object
            this->compute_my_adjoint_values();

            // get number of  RHS
            uint tNumRHS = mSet->mAdjointPdofValues.size();

            // get number of pdof values
            uint tNumPdofValues = tdRdp( 0 ).n_rows();

            // reorder adjoint values following the requested dof types order
            Matrix< DDRMat > tAdjointPdofValuesReordered;

            // get leader dof type list from set
            const Vector< Vector< MSI::Dof_Type > >& tLeaderDofTypeGroup =
                    mSet->get_dof_type_list( mtk::Leader_Follower::LEADER );

            // get number of leader dof types
            uint tNumLeaderDofTypes = tLeaderDofTypeGroup.size();

            // get follower dof type list from set
            const Vector< Vector< MSI::Dof_Type > >& tFollowerDofTypeGroup =
                    mSet->get_dof_type_list( mtk::Leader_Follower::FOLLOWER );

            // get number of follower dof types
            uint tNumFollowerDofTypes = tFollowerDofTypeGroup.size();

            // loop over the RHS
            for ( uint Ik = 0; Ik < tNumRHS; Ik++ )
            {
                // set size for reordered adjoint values
                tAdjointPdofValuesReordered.set_size( tNumPdofValues, 1, 0.0 );

                // loop over the leader dof types
                for ( uint Ia = 0; Ia < tNumLeaderDofTypes; Ia++ )
                {
                    // get the adjoint values for the ith dof type group
                    Vector< Vector< Matrix< DDRMat > > > tLeaderAdjointOriginal;
                    this->get_my_pdof_values(
                            mSet->mAdjointPdofValues,
                            tLeaderDofTypeGroup( Ia ),
                            tLeaderAdjointOriginal,
                            mtk::Leader_Follower::LEADER );

                    // reshape adjoint values
                    Matrix< DDRMat > tLeaderAdjointCoeff;
                    this->reshape_pdof_values_vector( tLeaderAdjointOriginal( Ik ), tLeaderAdjointCoeff );

                    // get indices for begin and end
                    uint tDofIndex   = mSet->get_dof_index_for_type( tLeaderDofTypeGroup( Ia )( 0 ), mtk::Leader_Follower::LEADER );
                    uint tStartIndex = mSet->get_res_dof_assembly_map()( tDofIndex )( 0, 0 );
                    uint tStopIndex  = mSet->get_res_dof_assembly_map()( tDofIndex )( 0, 1 );

                    // fill reordered adjoint pdof values
                    tAdjointPdofValuesReordered( { tStartIndex, tStopIndex } ) =
                            tLeaderAdjointCoeff.matrix_data();
                }

                // loop over the follower dof types
                for ( uint Ia = 0; Ia < tNumFollowerDofTypes; Ia++ )
                {
                    // get the adjoint values for the ith dof type group
                    Vector< Vector< Matrix< DDRMat > > > tFollowerAdjointOriginal;
                    this->get_my_pdof_values(
                            mSet->mAdjointPdofValues,
                            tFollowerDofTypeGroup( Ia ),
                            tFollowerAdjointOriginal,
                            mtk::Leader_Follower::FOLLOWER );

                    // reshape adjoint values
                    Matrix< DDRMat > tFollowerAdjointCoeff;
                    this->reshape_pdof_values_vector( tFollowerAdjointOriginal( Ik ), tFollowerAdjointCoeff );

                    // get indices for begin and end
                    uint tDofIndex   = mSet->get_dof_index_for_type( tFollowerDofTypeGroup( Ia )( 0 ), mtk::Leader_Follower::FOLLOWER );
                    uint tStartIndex = mSet->get_res_dof_assembly_map()( tDofIndex )( 0, 0 );
                    uint tStopIndex  = mSet->get_res_dof_assembly_map()( tDofIndex )( 0, 1 );

                    // fill reordered adjoint pdof values
                    tAdjointPdofValuesReordered( { tStartIndex, tStopIndex } ) =
                            tFollowerAdjointCoeff.matrix_data();
                }

                // Assembly for the IP pdv
                //----------------------------------------------------------------------------------------
                // if the assembly vector is not empty
                if ( tLocalToGlobalIdsIPPdv.numel() != 0 )
                {
                    // post multiplication of adjoint values time dRdp
                    Matrix< DDRMat > tLocalIPdQiDp =
                            -1.0 * trans( tAdjointPdofValuesReordered ) * tdRdp( 0 );

                    // assemble implicit dQidp into multivector
                    mEquationSet->get_equation_model()->get_implicit_dQidp()->sum_into_global_values(
                            tLocalToGlobalIdsIPPdv,
                            tLocalIPdQiDp,
                            Ik );
                }

                // Assembly for the IG pdv
                //----------------------------------------------------------------------------------------
                // if assembly vector is not empty
                if ( tLocalToGlobalIdsIGPdv.numel() != 0 )
                {
                    // post multiplication of adjoint values time dRdp
                    Matrix< DDRMat > tLocalIGdQiDp =
                            -1.0 * trans( tAdjointPdofValuesReordered ) * tdRdp( 1 );

                    // assemble implicit dQidp into multivector
                    mEquationSet->get_equation_model()->get_implicit_dQidp()->sum_into_global_values(
                            tLocalToGlobalIdsIGPdv,
                            tLocalIGdQiDp,
                            Ik );
                }
            }
        }

        //------------------------------------------------------------------------------

        void
        Interpolation_Element::compute_dQIdu()
        {
            // compute pdof values
            // FIXME do this only once
            this->compute_my_pdof_values();

            // if time continuity set
            if ( mSet->get_time_continuity() )
            {
                // compute pdof values for previous time step
                // FIXME do this only once
                this->compute_previous_pdof_values();
            }

            // set the field interpolators coefficients
            this->set_field_interpolators_coefficients();

            // FIXME should not be like this
            mSet->set_IQI_field_interpolator_managers();

            // set cluster for stabilization parameter
            mSet->set_IQI_cluster_for_stabilization_parameters( mFemCluster( 0 ).get() );

            // ask cluster to compute jacobian
            mFemCluster( 0 )->compute_dQIdu();
        }

        //------------------------------------------------------------------------------

        void
        Interpolation_Element::compute_QI()
        {
            // compute pdof values
            // FIXME do this only once
            this->compute_my_pdof_values();

            // if time continuity set
            if ( mSet->get_time_continuity() )
            {
                // compute pdof values for previous time step
                // FIXME do this only once
                this->compute_previous_pdof_values();
            }

            // if eigen vectors
            if ( mSet->mNumEigenVectors )
            {
                // compute pdof values for previous time step
                // FIXME do this only once
                this->compute_my_eigen_vector_values();
            }

            // initialize IQI
            mSet->initialize_mQI();

            // set the field interpolators coefficients
            this->set_field_interpolators_coefficients();

            // FIXME should not be like this
            mSet->set_IQI_field_interpolator_managers();

            // set cluster for stabilization parameter
            mSet->set_IQI_cluster_for_stabilization_parameters( mFemCluster( 0 ).get() );

            // ask cluster to compute quantity of interest
            mFemCluster( 0 )->compute_QI();
        }

        //------------------------------------------------------------------------------

        void
        Interpolation_Element::compute_quantity_of_interest(
                const uint           aFemMeshIndex,
                enum vis::Field_Type aFieldType )
        {
            // compute elemental coefficients and set FIs for evaluating IQIs
            this->setup_IQI_computation();

            // if nodal field
            if ( aFieldType == vis::Field_Type::NODAL )
            {
                // get the type of the current cluster
                enum fem::Element_Type tElementType = mFemCluster( aFemMeshIndex )->get_element_type();

                // compute the nodal values of the requested QIs
                // double sided clusters require additional treatment of the follower side and therefore need special consideration
                if ( tElementType == fem::Element_Type::DOUBLE_SIDESET )
                {
                    this->compute_nodal_QIs_double_sided( aFemMeshIndex );
                }
                else    // on BULK or SIDESET elements use the default
                {
                    this->compute_nodal_QIs_standard( aFemMeshIndex );
                }
            }

            // if elemental or global field
            else
            {
                // ask cluster to compute quantity of interest
                mFemCluster( aFemMeshIndex )->compute_quantity_of_interest( aFemMeshIndex, aFieldType );
            }
        }

        //------------------------------------------------------------------------------

        void
        Interpolation_Element::setup_IQI_computation()
        {
            // compute pdof values
            this->compute_my_pdof_values();    // FIXME: do this only once

            // if time continuity set
            if ( mSet->get_time_continuity() )
            {
                // compute pdof values for previous time step
                this->compute_previous_pdof_values();    // FIXME: do this only once
            }

            // if eigen vectors
            if ( mSet->mNumEigenVectors )
            {
                // compute pdof values for previous time step
                this->compute_my_eigen_vector_values();    // FIXME: do this only once
            }

            // set the field interpolators coefficients
            this->set_field_interpolators_coefficients();

            // give IQI access to the field interpolator
            mSet->set_IQI_field_interpolator_managers();    // FIXME: should not be like this

            // set cluster for stabilization parameter
            mSet->set_IQI_cluster_for_stabilization_parameters( mFemCluster( 0 ).get() );
        }

        //------------------------------------------------------------------------------

        void
        Interpolation_Element::compute_nodal_QIs_standard(
                const uint           aFemMeshIndex,
                mtk::Leader_Follower aLeaderOrFollowerSide )
        {
            // get the vertex indices on the mesh cluster
            Matrix< IndexMat > tVertexIndices;
            mFemCluster( aFemMeshIndex )->get_vertex_indices_in_cluster_for_visualization( tVertexIndices, aLeaderOrFollowerSide );

            // get the vertices local coordinates on the mesh cluster
            Matrix< moris::DDRMat > tVertexLocalCoords =
                    mFemCluster( aFemMeshIndex )->get_vertices_local_coordinates_wrt_interp_cell( aLeaderOrFollowerSide );

            // get number of vertices on the treated mesh cluster
            uint tNumVertices = tVertexLocalCoords.n_rows();

            // loop over the vertices on the treated mesh cluster
            for ( uint iVertex = 0; iVertex < tNumVertices; iVertex++ )
            {
                // get the i-th vertex's coordinates in the IP parametric space
                Matrix< DDRMat > tNodalPointLocalCoords = tVertexLocalCoords.get_row( iVertex );
                tNodalPointLocalCoords.resize( 1, tNodalPointLocalCoords.numel() + 1 );
                tNodalPointLocalCoords( tNodalPointLocalCoords.numel() - 1 ) = -1.0;
                tNodalPointLocalCoords                                       = trans( tNodalPointLocalCoords );

                // set vertex coordinates for field interpolator
                mSet->get_field_interpolator_manager()->set_space_time( tNodalPointLocalCoords );

                // set vertex coordinates for field interpolator of eigen vectors
                if ( mSet->mNumEigenVectors )
                {
                    mSet->get_field_interpolator_manager_eigen_vectors()->set_space_time( tNodalPointLocalCoords );
                }

                // get the current vertex's coordinates
                moris_index tVertexIndex = tVertexIndices( iVertex );

                // get number of active local IQIs
                uint tNumLocalIQIs = mSet->get_number_of_requested_nodal_IQIs_for_visualization();

                // loop over IQI
                for ( uint iIQI = 0; iIQI < tNumLocalIQIs; iIQI++ )
                {
                    // get requested IQI
                    const std::shared_ptr< IQI >& tReqIQI =
                            mSet->get_requested_nodal_IQIs_for_visualization()( iIQI );

                    // get IQI global index
                    moris_index tGlobalIqiIndex =
                            mSet->get_requested_nodal_IQIs_global_indices_for_visualization()( iIQI );

                    // reset the requested IQI
                    tReqIQI->reset_eval_flags();

                    // compute quantity of interest at evaluation point
                    Matrix< DDRMat > tQINodal( 1, 1, 0.0 );
                    tReqIQI->compute_QI( tQINodal );

                    // assemble the nodal QI value on the set
                    ( *mSet->mSetNodalValues )( tVertexIndex, tGlobalIqiIndex ) = tQINodal( 0 );
                }

            }    // end for: vertices on cluster

        }        // end function: Interpolation_Element::compute_nodal_QIs()

        //------------------------------------------------------------------------------

        void
        Interpolation_Element::compute_nodal_QIs_double_sided( const uint aFemMeshIndex )
        {
            // make sure this function is only called on double sided side clusters
            MORIS_ASSERT(
                    mFemCluster( aFemMeshIndex )->get_element_type() == fem::Element_Type::DOUBLE_SIDESET,
                    "Interpolation_Element::compute_nodal_QIs_double_sided() - "
                    "This function can only be called on double sided side clusters." );

            // get the VIS vertex indices on the mesh cluster
            Matrix< IndexMat > tLeaderVisVertexIndices;
            Matrix< IndexMat > tFollowerVisVertexIndices;
            mFemCluster( aFemMeshIndex )
                    ->get_vertex_indices_in_cluster_for_visualization( tLeaderVisVertexIndices, mtk::Leader_Follower::LEADER );
            mFemCluster( aFemMeshIndex )
                    ->get_vertex_indices_in_cluster_for_visualization( tFollowerVisVertexIndices, mtk::Leader_Follower::FOLLOWER );

            // get the vertices' local coordinates on the respective mesh clusters
            Matrix< moris::DDRMat > tLeaderVertexLocalCoords =
                    mFemCluster( aFemMeshIndex )->get_vertices_local_coordinates_wrt_interp_cell( mtk::Leader_Follower::LEADER );
            Matrix< moris::DDRMat > tFollowerVertexLocalCoords =
                    mFemCluster( aFemMeshIndex )->get_vertices_local_coordinates_wrt_interp_cell( mtk::Leader_Follower::FOLLOWER );

            // get number of vertices on the treated mesh cluster
            uint tNumVertices = tLeaderVertexLocalCoords.n_rows();

            // deduce the number of spatial dimensions
            uint tNumSpatialDims = tLeaderVertexLocalCoords.n_cols();

            // make sure the number of vertices on both clusters is the same (since the assumption is made that they correspond to each other)
            MORIS_ASSERT(
                    tNumVertices == tFollowerVertexLocalCoords.n_rows(),
                    "Interpolation_Element::compute_nodal_QIs_double_sided() - "
                    "The number of vertices on the Leader and Follower VIS clusters don't match. "
                    "Here the assumption is made that the two hold the same (interface only) vertices." );

            // loop over the interfaces vertices on the treated dbl side cluster cluster
            for ( uint iVertex = 0; iVertex < tNumVertices; iVertex++ )
            {
                // get the vertex's coordinates with respect to both the leader and follower cells
                Matrix< DDRMat > tNodalPointLeaderLocalCoords   = tLeaderVertexLocalCoords.get_row( iVertex );
                Matrix< DDRMat > tNodalPointFollowerLocalCoords = tFollowerVertexLocalCoords.get_row( iVertex );

                // add the time local coordinate and set it to -1 (beginning of time slab)
                tNodalPointLeaderLocalCoords.resize( 1, tNumSpatialDims + 1 ); // add time dimension
                tNodalPointFollowerLocalCoords.resize( 1, tNumSpatialDims + 1 );
                tNodalPointLeaderLocalCoords( tNumSpatialDims ) = -1.0; // set time location to -1
                tNodalPointFollowerLocalCoords( tNumSpatialDims ) = -1.0;

                // field interpolator takes the transpose
                tNodalPointLeaderLocalCoords = trans( tNodalPointLeaderLocalCoords );
                tNodalPointFollowerLocalCoords = trans( tNodalPointFollowerLocalCoords );

                // set vertex coordinates for field interpolators on the leader and follower sides
                mSet->get_field_interpolator_manager( mtk::Leader_Follower::LEADER  )->set_space_time( tNodalPointLeaderLocalCoords );
                mSet->get_field_interpolator_manager( mtk::Leader_Follower::FOLLOWER )->set_space_time( tNodalPointFollowerLocalCoords );

                // set vertex coordinates for field interpolator of eigen vectors
                if ( mSet->mNumEigenVectors )
                {
                    mSet->get_field_interpolator_manager_eigen_vectors( mtk::Leader_Follower::LEADER  )
                            ->set_space_time( tNodalPointLeaderLocalCoords );
                    mSet->get_field_interpolator_manager_eigen_vectors( mtk::Leader_Follower::FOLLOWER )
                            ->set_space_time( tNodalPointFollowerLocalCoords );
                }

                // get the current vertex's coordinates
                moris_index tLeaderVertexIndex = tLeaderVisVertexIndices( iVertex );

                // get number of active local IQIs
                uint tNumLocalIQIs = mSet->get_number_of_requested_nodal_IQIs_for_visualization();

                // loop over IQI
                for ( uint iIQI = 0; iIQI < tNumLocalIQIs; iIQI++ )
                {
                    // get requested IQI
                    const std::shared_ptr< IQI >& tReqIQI =
                            mSet->get_requested_nodal_IQIs_for_visualization()( iIQI );

                    // get IQI global index
                    moris_index tGlobalIqiIndex =
                            mSet->get_requested_nodal_IQIs_global_indices_for_visualization()( iIQI );

                    // reset the requested IQI
                    tReqIQI->reset_eval_flags();

                    // compute quantity of interest at evaluation point
                    Matrix< DDRMat > tQINodal( 1, 1, 0.0 );
                    tReqIQI->compute_QI( tQINodal );

                    // assemble the nodal QI value on the set
                    // NOTE: output of the dbl. sided elements by the VIS mesh is defined on the leader side; hence, the nodal value is outputted to the leader vertex
                    ( *mSet->mSetNodalValues )( tLeaderVertexIndex, tGlobalIqiIndex ) = tQINodal( 0 );
                }

            }    // end for: vertices on cluster
        }

        //------------------------------------------------------------------------------

        void
        Interpolation_Element::populate_fields(
                Vector< std::shared_ptr< fem::Field > >& aFields,
                Vector< std::string > const &            tFieldIQINames )
        {
            // compute pdof values
            // FIXME do this only once
            this->compute_my_pdof_values();

            // if time continuity set
            if ( mSet->get_time_continuity() )
            {
                // compute pdof values for previous time step
                // FIXME do this only once
                this->compute_previous_pdof_values();
            }

            // set the field interpolators coefficients
            this->set_field_interpolators_coefficients();

            // FIXME should not be like this
            mSet->set_IQI_field_interpolator_managers();

            // set cluster for stabilization parameter
            mSet->set_IQI_cluster_for_stabilization_parameters( mFemCluster( 0 ).get() );

            const Vector< std::shared_ptr< IQI > >& tIQI =
                    mSet->get_requested_field_IQIs();

            // get number of active local IQIs
            uint tNumIQIs = tIQI.size();

            for ( uint iIQI = 0; iIQI < tNumIQIs; iIQI++ )
            {
                moris_index tGlobalIndex =
                        mSet->get_requested_field_IQIs_global_indices()( iIQI );

                // if nodal field
                if ( aFields( tGlobalIndex )->get_field_entity_type() == mtk::Field_Entity_Type::NODAL )
                {
                    // // get the leader vertices indices on the mesh cluster
                    // Matrix< IndexMat > tVertexIndices = mFemCluster( 0 )->
                    //         get_mesh_cluster()->
                    //         get_interpolation_cell().
                    //         get_vertex_inds();

                    Matrix< IndexMat > tVertexIndices =
                            mLeaderInterpolationCell->get_vertex_inds();

                    // get the leader vertices local coordinates on the interpolation element
                    Geometry_Interpolator* tIPGI =
                            mSet->get_field_interpolator_manager()->get_IP_geometry_interpolator();

                    const Matrix< DDRMat > tIGLocalCoords = tIPGI->get_space_param_coeff();

                    // get number of vertices on the treated mesh cluster
                    uint tNumNodes = tIGLocalCoords.n_rows();

                    // loop over the vertices on the treated mesh cluster
                    for ( uint iVertex = 0; iVertex < tNumNodes; iVertex++ )
                    {
                        // get the ith vertex coordinates in the IP param space
                        Matrix< DDRMat > tIntegPoint = tIGLocalCoords.get_row( iVertex );
                        tIntegPoint.resize( 1, tIntegPoint.numel() + 1 );
                        tIntegPoint( tIntegPoint.numel() - 1 ) = -1.0;
                        tIntegPoint                            = trans( tIntegPoint );

                        // set vertex coordinates for field interpolator
                        mSet->get_field_interpolator_manager()->set_space_time( tIntegPoint );

                        // reset the requested IQI
                        tIQI( iIQI )->reset_eval_flags();

                        // compute quantity of interest at evaluation point
                        Matrix< DDRMat > tQINodal;
                        tIQI( iIQI )->compute_QI( tQINodal );

                        aFields( tGlobalIndex )->set_field_value( tVertexIndices( iVertex ), tQINodal );
                    }
                }
                else
                {
                    // get mtk interpolation element index
                    moris_index tIndex = mLeaderInterpolationCell->get_index();

                    Matrix< DDRMat > tValues( 1, 1, 0.0 );

                    MORIS_ASSERT( aFields( tGlobalIndex )->get_value( tIndex, 0 ) == MORIS_REAL_MIN, "elemental field values previously set." );

                    real tSpaceTimeVolume = 0.0;

                    // ask cluster to compute quantity of interest
                    mFemCluster( 0 )->compute_quantity_of_interest(
                            tValues,
                            aFields( tGlobalIndex )->get_field_entity_type(),
                            iIQI,
                            tSpaceTimeVolume );

                    tValues( 0 ) = tValues( 0 ) / tSpaceTimeVolume;

                    aFields( tGlobalIndex )->set_field_value( tIndex, tValues );
                }

            }    // end for: each requested IQIs for field output

        }        // end function: Interpolation_Element::populate_fields()

        //------------------------------------------------------------------------------
    } /* namespace fem */
} /* namespace moris */
