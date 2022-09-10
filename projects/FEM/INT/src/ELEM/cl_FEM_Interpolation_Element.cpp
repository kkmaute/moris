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
        const Cell< const mtk::Cell* >&  aInterpolationCell,
        const moris::Cell< Node_Base* >& aNodes,
        Set*                             aSet ) :
        MSI::Equation_Object( aSet ),
        mSet( aSet ),
        mElementType( aElementType )
    {
        // fill the master interpolation cell
        mMasterInterpolationCell = aInterpolationCell( 0 );

        // get vertices from cell
        moris::Cell< mtk::Vertex* > tVertices =
            mMasterInterpolationCell->get_vertex_pointers();

        // get number of vertices from cell
        uint tNumOfVertices = tVertices.size();

        // assign node object
        mNodeObj.resize( 1 );
        mNodeObj( 0 ).resize( tNumOfVertices, nullptr );

        // fill master node objects
        for ( uint iVertex = 0; iVertex < tNumOfVertices; iVertex++ )
        {
            mNodeObj( 0 )( iVertex ) = aNodes( tVertices( iVertex )->get_index() );
        }

        // if double sided sideset
        if ( mElementType == fem::Element_Type::DOUBLE_SIDESET )
        {
            // fill the slave interpolation cell
            mSlaveInterpolationCell = aInterpolationCell( 1 );

            // get vertices from cell
            moris::Cell< mtk::Vertex* > tSlaveVertices =
                mSlaveInterpolationCell->get_vertex_pointers();

            // get number of vertices from cell
            uint tNumOfSlaveVertices = tSlaveVertices.size();

            // assign node object
            mNodeObj.resize( 2 );
            mNodeObj( 1 ).resize( tNumOfSlaveVertices, nullptr );

            // fill slave node objects
            for ( uint iVertex = 0; iVertex < tNumOfSlaveVertices; iVertex++ )
            {
                mNodeObj( 1 )( iVertex ) = aNodes( tSlaveVertices( iVertex )->get_index() );
            }
        }
    }

    //------------------------------------------------------------------------------

    void
    Interpolation_Element::set_cluster(
        std::shared_ptr< fem::Cluster > aCluster,
        const uint                      aMeshIndex )
    {
        // if mesh index is 0 (i.e., forward analysis mesh, IG mesh)
        if ( aMeshIndex == 0 )
        {
            // fem cluster with index 0 should be set only once and shall not be changed
            MORIS_ASSERT( !( mFemCluster.size() >= 1 ),
                "Interpolation_Element::set_cluster() - first fem cluster is already set" );
        }

        // get max size for fem cluster list
        sint tSize = std::max( (sint)mFemCluster.size(), (sint)aMeshIndex + 1 );

        // resize fem cluster list
        mFemCluster.resize( tSize );

        // add the fem cluster to the list
        mFemCluster( aMeshIndex ) = aCluster;
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

        // get master dof type list from set
        Cell< Cell< MSI::Dof_Type > >& tMasterDofTypeList =
            mSet->get_dof_type_list( mtk::Master_Slave::MASTER );

        // get number of master dof types
        uint tMasterNumDofTypes = tMasterDofTypeList.size();

        // loop on the master dof types
        for ( uint iDOF = 0; iDOF < tMasterNumDofTypes; iDOF++ )
        {
            // get the ith dof type group
            moris::Cell< MSI::Dof_Type >& tDofTypeGroup = tMasterDofTypeList( iDOF );

            // get the pdof values for the ith dof type group
            Cell< Cell< Matrix< DDRMat > > > tCoeff_Original;
            this->get_my_pdof_values( mSet->mPdofValues, tDofTypeGroup, tCoeff_Original );

            // reshape tCoeffs into the order the cluster expects them
            Matrix< DDRMat > tCoeff;
            this->reshape_pdof_values( tCoeff_Original( 0 ), tCoeff );

            // set field interpolator coefficients
            mSet->get_field_interpolator_manager()->set_coeff_for_type( tDofTypeGroup( 0 ), tCoeff );

            if ( mSet->get_time_continuity() )
            {
                // get the pdof values for the ith dof type group
                Cell< Cell< Matrix< DDRMat > > > tCoeff_Original;
                this->get_my_pdof_values( mSet->mPreviousPdofValues, tDofTypeGroup, tCoeff_Original );

                // reshape tCoeffs into the order the cluster expects them
                Matrix< DDRMat > tCoeff;
                this->reshape_pdof_values( tCoeff_Original( 0 ), tCoeff );

                // set field interpolator coefficients
                mSet->get_field_interpolator_manager_previous_time()->set_coeff_for_type( tDofTypeGroup( 0 ), tCoeff );
            }
        }

        // get slave dof type list from set
        Cell< Cell< MSI::Dof_Type > >& tSlaveDofTypeList =
            mSet->get_dof_type_list( mtk::Master_Slave::SLAVE );

        // get number of slave dof types
        uint tSlaveNumDofTypes = tSlaveDofTypeList.size();

        // loop on the slave dof types
        for ( uint iDOF = 0; iDOF < tSlaveNumDofTypes; iDOF++ )
        {
            // get the ith dof type group
            moris::Cell< MSI::Dof_Type >& tDofTypeGroup = tSlaveDofTypeList( iDOF );

            // get the pdof values for the ith dof type group
            Cell< Cell< Matrix< DDRMat > > > tCoeff_Original;
            this->get_my_pdof_values( mSet->mPdofValues, tDofTypeGroup, tCoeff_Original, mtk::Master_Slave::SLAVE );

            // reshape tCoeffs into the order the cluster expects them
            Matrix< DDRMat > tCoeff;
            this->reshape_pdof_values( tCoeff_Original( 0 ), tCoeff );

            // set the field coefficients
            mSet->get_field_interpolator_manager( mtk::Master_Slave::SLAVE )->set_coeff_for_type( tDofTypeGroup( 0 ), tCoeff );
        }

        // dv field interpolators------------------------------------------

        // get master dv type list from set
        const Cell< Cell< PDV_Type > >& tMasterDvTypeList =
            mSet->get_dv_type_list( mtk::Master_Slave::MASTER );

        // get number of master dv types
        uint tMasterNumDvTypes = tMasterDvTypeList.size();

        // loop on the master dv types
        for ( uint iDv = 0; iDv < tMasterNumDvTypes; iDv++ )
        {
            // get the dv type group
            const moris::Cell< PDV_Type >& tDvTypeGroup = tMasterDvTypeList( iDv );

            // get the pdv values for the ith dv type group
            // FIXME: the underlying use of the base cell needs to be hidden within PDV
            Cell< Matrix< DDRMat > > tCoeff_Original;
            mSet->get_equation_model()->get_design_variable_interface()->get_ip_pdv_value(
                mMasterInterpolationCell->get_base_cell()->get_vertex_inds(),
                tDvTypeGroup,
                tCoeff_Original );

            // reshape tCoeffs into the order the FI expects them
            Matrix< DDRMat > tCoeff;
            mSet->get_equation_model()->get_design_variable_interface()->reshape_pdv_values( tCoeff_Original, tCoeff );

            // set field interpolator coefficients
            mSet->get_field_interpolator_manager()->set_coeff_for_type( tDvTypeGroup( 0 ), tCoeff );
        }

        // get slave dv type list from set
        const Cell< Cell< PDV_Type > >& tSlaveDvTypeList =
            mSet->get_dv_type_list( mtk::Master_Slave::SLAVE );

        // get number of slave dv types
        uint tSlaveNumDvTypes = tSlaveDvTypeList.size();

        // loop on the slave dv types
        for ( uint iDv = 0; iDv < tSlaveNumDvTypes; iDv++ )
        {
            // get the dv type group
            const moris::Cell< PDV_Type >& tDvTypeGroup = tSlaveDvTypeList( iDv );

            // get the pdv values for the ith dv type group
            Cell< Matrix< DDRMat > > tCoeff_Original;
            mSet->get_equation_model()->get_design_variable_interface()->get_ip_pdv_value(
                mSlaveInterpolationCell->get_base_cell()->get_vertex_inds(),
                tDvTypeGroup,
                tCoeff_Original );

            // reshape tCoeffs into the order the FI expects them
            Matrix< DDRMat > tCoeff;
            mSet->get_equation_model()->get_design_variable_interface()->reshape_pdv_values( tCoeff_Original, tCoeff );

            // set the field coefficients
            mSet->get_field_interpolator_manager( mtk::Master_Slave::SLAVE )->set_coeff_for_type( tDvTypeGroup( 0 ), tCoeff );
        }

        // field field interpolators------------------------------------------

        // get master field type list from set
        const Cell< Cell< mtk::Field_Type > >& tMasterFieldTypeList =
            mSet->get_field_type_list( mtk::Master_Slave::MASTER );

        // get number of master field types
        uint tMasterNumFieldTypes = tMasterFieldTypeList.size();

        // loop on the master field types
        for ( uint iFi = 0; iFi < tMasterNumFieldTypes; iFi++ )
        {
            // get the field type group
            const moris::Cell< mtk::Field_Type >& tFieldTypeGroup = tMasterFieldTypeList( iFi );

            Matrix< IndexMat > tIPCellIndices = mMasterInterpolationCell->get_vertex_inds();

            Matrix< DDRMat > tCoeff;
            mSet->get_fem_model()->get_field( tFieldTypeGroup( 0 ) )->get_values( tIPCellIndices, tCoeff, tFieldTypeGroup );

            // FIXME implement reshape for vector fields

            // set field interpolator coefficients
            mSet->get_field_interpolator_manager()->set_coeff_for_type( tFieldTypeGroup( 0 ), tCoeff );
        }

        // get slave field type list from set
        const Cell< Cell< mtk::Field_Type > >& tSlaveFieldTypeList =
            mSet->get_field_type_list( mtk::Master_Slave::SLAVE );

        // get number of slave field types
        uint tSlaveNumFieldTypes = tSlaveFieldTypeList.size();

        // loop on the slave field types
        for ( uint iFi = 0; iFi < tSlaveNumFieldTypes; iFi++ )
        {
            // get the field type group
            const moris::Cell< mtk::Field_Type >& tFieldTypeGroup = tSlaveFieldTypeList( iFi );

            Matrix< IndexMat > tIPCellIndices = mSlaveInterpolationCell->get_vertex_inds();

            Matrix< DDRMat > tCoeff;
            mSet->get_fem_model()->get_field( tFieldTypeGroup( 0 ) )->get_values( tIPCellIndices, tCoeff, tFieldTypeGroup );

            // FIXME implement reshape for vector fields

            // set field interpolator coefficients
            mSet->get_field_interpolator_manager( mtk::Master_Slave::SLAVE )->set_coeff_for_type( tFieldTypeGroup( 0 ), tCoeff );
        }

        // geometry interpolators------------------------------------------
        // set the IP geometry interpolator physical space and time coefficients for the master
        mSet->get_field_interpolator_manager( mtk::Master_Slave::MASTER )->get_IP_geometry_interpolator()->set_space_coeff( mMasterInterpolationCell->get_vertex_coords() );
        mSet->get_field_interpolator_manager( mtk::Master_Slave::MASTER )->get_IP_geometry_interpolator()->set_param_coeff();
        mSet->get_field_interpolator_manager( mtk::Master_Slave::MASTER )->get_IP_geometry_interpolator()->set_time_coeff( this->get_time() );

        // if double sideset
        if ( mElementType == fem::Element_Type::DOUBLE_SIDESET )
        {
            // set the IP geometry interpolator physical space and time coefficients for the slave
            mSet->get_field_interpolator_manager( mtk::Master_Slave::SLAVE )->get_IP_geometry_interpolator()->set_space_coeff( mSlaveInterpolationCell->get_vertex_coords() );
            mSet->get_field_interpolator_manager( mtk::Master_Slave::SLAVE )->get_IP_geometry_interpolator()->set_time_coeff( this->get_time() );
        }

        // if time sideset
        if ( mElementType == fem::Element_Type::TIME_SIDESET )
        {
            // set the IP geometry interpolator physical space and time coefficients for the previous
            mSet->get_field_interpolator_manager_previous_time( mtk::Master_Slave::MASTER )->get_IP_geometry_interpolator()->set_space_coeff( mMasterInterpolationCell->get_vertex_coords() );
            mSet->get_field_interpolator_manager_previous_time( mtk::Master_Slave::MASTER )->get_IP_geometry_interpolator()->set_time_coeff( this->get_previous_time() );
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
        moris::Cell< moris::Cell< enum PDV_Type > > tRequestedDvTypes;
        mSet->get_ip_dv_types_for_set( tRequestedDvTypes );

        // reset material pdv assembly vector
        mSet->get_mat_pdv_assembly_vector().fill( -1 );

        // init pdv counter
        uint tCounter = 0;

        // get master vertices from cell
        Matrix< IndexMat > tMasterVerticesInds =
            mMasterInterpolationCell->get_base_cell()->get_vertex_inds();

        // loop over the dv types
        for ( uint Ik = 0; Ik < tRequestedDvTypes.size(); Ik++ )
        {
            // get dv ids for this type and node indices
            moris::Cell< moris::Matrix< IdMat > > tPdvIds;

            // get the pdv ids for requested vertices and pdv type
            tDVInterface->get_ip_dv_ids_for_type_and_ind(
                tMasterVerticesInds,
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
            // get slave vertices from cell
            Matrix< IndexMat > tSlaveVerticesInds =
                mSlaveInterpolationCell->get_base_cell()->get_vertex_inds();

            // get the list of requested dv types by the opt solver for the slave side
            mSet->get_ip_dv_types_for_set( tRequestedDvTypes, mtk::Master_Slave::SLAVE );

            // loop over the dv types
            for ( uint Ik = 0; Ik < tRequestedDvTypes.size(); Ik++ )
            {
                // get dv ids for this type and node indices
                moris::Cell< moris::Matrix< IdMat > > tPdvIds;

                // get the pdv ids for requested vertices and pdv type
                tDVInterface->get_ip_dv_ids_for_type_and_ind(
                    tSlaveVerticesInds,
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
        moris::Cell< Matrix< DDRMat > >& tdRdp = mEquationSet->get_drdp();

        // extract adjoint values for this equation object
        this->compute_my_adjoint_values();

        // get number of RHS
        uint tNumRHS = mSet->mAdjointPdofValues.size();

        // get number of pdof values
        uint tNumPdofValues = tdRdp( 0 ).n_rows();

        // reorder adjoint values following the requested dof types order
        Matrix< DDRMat > tAdjointPdofValuesReordered;

        // get master dof type list from set
        const Cell< Cell< MSI::Dof_Type > >& tMasterDofTypeGroup =
            mSet->get_dof_type_list( mtk::Master_Slave::MASTER );

        // get number of master dof types
        uint tNumMasterDofTypes = tMasterDofTypeGroup.size();

        // get slave dof type list from set
        const Cell< Cell< MSI::Dof_Type > >& tSlaveDofTypeGroup =
            mSet->get_dof_type_list( mtk::Master_Slave::SLAVE );

        // get number of slave dof types
        uint tNumSlaveDofTypes = tSlaveDofTypeGroup.size();

        // loop over the RHS
        for ( uint Ik = 0; Ik < tNumRHS; Ik++ )
        {
            // set size for reordered adjoint values
            tAdjointPdofValuesReordered.set_size( tNumPdofValues, 1, 0.0 );

            // loop over the master dof types
            for ( uint Ia = 0; Ia < tNumMasterDofTypes; Ia++ )
            {
                // get the adjoint values for the ith dof type group
                Cell< Cell< Matrix< DDRMat > > > tMasterAdjointOriginal;
                this->get_my_pdof_values(
                    mSet->mAdjointPdofValues,
                    tMasterDofTypeGroup( Ia ),
                    tMasterAdjointOriginal,
                    mtk::Master_Slave::MASTER );

                // reshape adjoint values
                Matrix< DDRMat > tMasterAdjointCoeff;
                this->reshape_pdof_values_vector( tMasterAdjointOriginal( Ik ), tMasterAdjointCoeff );

                // get indices for begin and end
                uint tDofIndex   = mSet->get_dof_index_for_type( tMasterDofTypeGroup( Ia )( 0 ), mtk::Master_Slave::MASTER );
                uint tStartIndex = mSet->get_res_dof_assembly_map()( tDofIndex )( 0, 0 );
                uint tStopIndex  = mSet->get_res_dof_assembly_map()( tDofIndex )( 0, 1 );

                // fill reordered adjoint pdof values
                tAdjointPdofValuesReordered( { tStartIndex, tStopIndex } ) =
                    tMasterAdjointCoeff.matrix_data();
            }

            // loop over the slave dof types
            for ( uint Ia = 0; Ia < tNumSlaveDofTypes; Ia++ )
            {
                // get the adjoint values for the ith dof type group
                Cell< Cell< Matrix< DDRMat > > > tSlaveAdjointOriginal;
                this->get_my_pdof_values(
                    mSet->mAdjointPdofValues,
                    tSlaveDofTypeGroup( Ia ),
                    tSlaveAdjointOriginal,
                    mtk::Master_Slave::SLAVE );

                // reshape adjoint values
                Matrix< DDRMat > tSlaveAdjointCoeff;
                this->reshape_pdof_values_vector( tSlaveAdjointOriginal( Ik ), tSlaveAdjointCoeff );

                // get indices for begin and end
                uint tDofIndex   = mSet->get_dof_index_for_type( tSlaveDofTypeGroup( Ia )( 0 ), mtk::Master_Slave::SLAVE );
                uint tStartIndex = mSet->get_res_dof_assembly_map()( tDofIndex )( 0, 0 );
                uint tStopIndex  = mSet->get_res_dof_assembly_map()( tDofIndex )( 0, 1 );

                // fill reordered adjoint pdof values
                tAdjointPdofValuesReordered( { tStartIndex, tStopIndex } ) =
                    tSlaveAdjointCoeff.matrix_data();
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
                moris::Matrix< DDRMat > tLocalIPdQiDp =
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
                moris::Matrix< DDRMat > tLocalIGdQiDp =
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
        moris::Cell< Matrix< DDRMat > >& tdRdp = mEquationSet->get_drdp();

        // extract adjoint values for this equation object
        this->compute_my_adjoint_values();

        // get number of  RHS
        uint tNumRHS = mSet->mAdjointPdofValues.size();

        // get number of pdof values
        uint tNumPdofValues = tdRdp( 0 ).n_rows();

        // reorder adjoint values following the requested dof types order
        Matrix< DDRMat > tAdjointPdofValuesReordered;

        // get master dof type list from set
        const Cell< Cell< MSI::Dof_Type > >& tMasterDofTypeGroup =
            mSet->get_dof_type_list( mtk::Master_Slave::MASTER );

        // get number of master dof types
        uint tNumMasterDofTypes = tMasterDofTypeGroup.size();

        // get slave dof type list from set
        const Cell< Cell< MSI::Dof_Type > >& tSlaveDofTypeGroup =
            mSet->get_dof_type_list( mtk::Master_Slave::SLAVE );

        // get number of slave dof types
        uint tNumSlaveDofTypes = tSlaveDofTypeGroup.size();

        // loop over the RHS
        for ( uint Ik = 0; Ik < tNumRHS; Ik++ )
        {
            // set size for reordered adjoint values
            tAdjointPdofValuesReordered.set_size( tNumPdofValues, 1, 0.0 );

            // loop over the master dof types
            for ( uint Ia = 0; Ia < tNumMasterDofTypes; Ia++ )
            {
                // get the adjoint values for the ith dof type group
                Cell< Cell< Matrix< DDRMat > > > tMasterAdjointOriginal;
                this->get_my_pdof_values(
                    mSet->mAdjointPdofValues,
                    tMasterDofTypeGroup( Ia ),
                    tMasterAdjointOriginal,
                    mtk::Master_Slave::MASTER );

                // reshape adjoint values
                Matrix< DDRMat > tMasterAdjointCoeff;
                this->reshape_pdof_values_vector( tMasterAdjointOriginal( Ik ), tMasterAdjointCoeff );

                // get indices for begin and end
                uint tDofIndex   = mSet->get_dof_index_for_type( tMasterDofTypeGroup( Ia )( 0 ), mtk::Master_Slave::MASTER );
                uint tStartIndex = mSet->get_res_dof_assembly_map()( tDofIndex )( 0, 0 );
                uint tStopIndex  = mSet->get_res_dof_assembly_map()( tDofIndex )( 0, 1 );

                // fill reordered adjoint pdof values
                tAdjointPdofValuesReordered( { tStartIndex, tStopIndex } ) =
                    tMasterAdjointCoeff.matrix_data();
            }

            // loop over the slave dof types
            for ( uint Ia = 0; Ia < tNumSlaveDofTypes; Ia++ )
            {
                // get the adjoint values for the ith dof type group
                Cell< Cell< Matrix< DDRMat > > > tSlaveAdjointOriginal;
                this->get_my_pdof_values(
                    mSet->mAdjointPdofValues,
                    tSlaveDofTypeGroup( Ia ),
                    tSlaveAdjointOriginal,
                    mtk::Master_Slave::SLAVE );

                // reshape adjoint values
                Matrix< DDRMat > tSlaveAdjointCoeff;
                this->reshape_pdof_values_vector( tSlaveAdjointOriginal( Ik ), tSlaveAdjointCoeff );

                // get indices for begin and end
                uint tDofIndex   = mSet->get_dof_index_for_type( tSlaveDofTypeGroup( Ia )( 0 ), mtk::Master_Slave::SLAVE );
                uint tStartIndex = mSet->get_res_dof_assembly_map()( tDofIndex )( 0, 0 );
                uint tStopIndex  = mSet->get_res_dof_assembly_map()( tDofIndex )( 0, 1 );

                // fill reordered adjoint pdof values
                tAdjointPdofValuesReordered( { tStartIndex, tStopIndex } ) =
                    tSlaveAdjointCoeff.matrix_data();
            }

            // Assembly for the IP pdv
            //----------------------------------------------------------------------------------------
            // if the assembly vector is not empty
            if ( tLocalToGlobalIdsIPPdv.numel() != 0 )
            {
                // post multiplication of adjoint values time dRdp
                moris::Matrix< DDRMat > tLocalIPdQiDp =
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
                moris::Matrix< DDRMat > tLocalIGdQiDp =
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
        const uint           aMeshIndex,
        enum vis::Field_Type aFieldType )
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

        // if nodal field
        if ( aFieldType == vis::Field_Type::NODAL )
        {
            // get the master vertices indices on the mesh cluster
            moris::Matrix< moris::IndexMat > tVertexIndices;
            mFemCluster( aMeshIndex )->get_vertex_indices_in_cluster_for_visualization( tVertexIndices );

            // get the master vertices local coordinates on the mesh cluster
            moris::Matrix< moris::DDRMat > tVertexLocalCoords =
                mFemCluster( aMeshIndex )->get_vertices_local_coordinates_wrt_interp_cell();

            // get number of vertices on the treated mesh cluster
            uint tNumNodes = tVertexLocalCoords.n_rows();

            // loop over the vertices on the treated mesh cluster
            for ( uint iVertex = 0; iVertex < tNumNodes; iVertex++ )
            {
                // get the ith vertex coordinates in the IP param space
                Matrix< DDRMat > tGlobalIntegPoint = tVertexLocalCoords.get_row( iVertex );
                tGlobalIntegPoint.resize( 1, tGlobalIntegPoint.numel() + 1 );
                tGlobalIntegPoint( tGlobalIntegPoint.numel() - 1 ) = -1.0;
                tGlobalIntegPoint                                  = trans( tGlobalIntegPoint );

                // set vertex coordinates for field interpolator
                mSet->get_field_interpolator_manager()->set_space_time( tGlobalIntegPoint );

                // get number of active local IQIs
                uint tNumLocalIQIs = mSet->get_number_of_requested_nodal_IQIs_for_visualization();

                // loop over IQI
                for ( uint iIQI = 0; iIQI < tNumLocalIQIs; iIQI++ )
                {
                    // get requested IQI
                    const std::shared_ptr< IQI >& tReqIQI =
                        mSet->get_requested_nodal_IQIs_for_visualization()( iIQI );

                    // get IQI global index
                    moris_index tGlobalIndex =
                        mSet->get_requested_nodal_IQIs_global_indices_for_visualization()( iIQI );

                    // reset the requested IQI
                    tReqIQI->reset_eval_flags();

                    // compute quantity of interest at evaluation point
                    Matrix< DDRMat > tQINodal( 1, 1, 0.0 );
                    tReqIQI->compute_QI( tQINodal );

                    // assemble the nodal QI value on the set
                    ( *mSet->mSetNodalValues )( tVertexIndices( iVertex ), tGlobalIndex ) = tQINodal( 0 );
                }
            }
        }
        else
        {
            // ask cluster to compute quantity of interest
            mFemCluster( aMeshIndex )->compute_quantity_of_interest( aMeshIndex, aFieldType );
        }
    }

    //------------------------------------------------------------------------------

    void
    Interpolation_Element::populate_fields(
        moris::Cell< std::shared_ptr< fem::Field > >& aFields,
        moris::Cell< std::string > const&             tFieldIQINames )
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

        const moris::Cell< std::shared_ptr< IQI > >& tIQI =
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
                //                // get the master vertices indices on the mesh cluster
                //                moris::Matrix< moris::IndexMat > tVertexIndices = mFemCluster( 0 )->
                //                        get_mesh_cluster()->
                //                        get_interpolation_cell().
                //                        get_vertex_inds();

                Matrix< IndexMat > tVertexIndices =
                    mMasterInterpolationCell->get_vertex_inds();

                // get the master vertices local coordinates on the interpolation element
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
                moris_index tIndex = mMasterInterpolationCell->get_index();

                moris::Matrix< DDRMat > tValues( 1, 1, 0.0 );

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
        }
    }

    //------------------------------------------------------------------------------
} /* namespace fem */
} /* namespace moris */

