/*
 * cl_FEM_Set.cpp
 *
 *  Created on: Apr 11, 2019
 *      Author: schmidt/noel
 */
#include <iostream>

#include "cl_MSI_Model_Solver_Interface.hpp" //FEM/MSI/src
#include "cl_FEM_Set.hpp"                    //FEM/INT/src
#include "cl_FEM_Set_User_Info.hpp"                    //FEM/INT/src
#include "cl_FEM_Element_Factory.hpp"        //FEM/INT/src
#include "cl_FEM_Integrator.hpp"             //FEM/INT/src

#include "cl_MTK_Set.hpp"             //FEM/INT/src
#include "fn_equal_to.hpp"

namespace moris
{
    namespace fem
    {
//------------------------------------------------------------------------------
    Set::Set( moris::mtk::Set           * aMeshSet,
              fem::Set_User_Info        & aSetInfo,
              moris::Cell< Node_Base* > & aIPNodes )
    : mMeshSet( aMeshSet ),
      mNodes( aIPNodes ),
      mIWGs( aSetInfo.get_IWGs() ),
      mElementType( aSetInfo.get_set_type() )
    {
        // get mesh clusters on set
        mMeshClusterList = mMeshSet->get_clusters_on_set();

        // get number of mesh clusters on set
        uint tNumMeshClusters = mMeshClusterList.size();

        // set size for the equation objects list
        mEquationObjList.resize( tNumMeshClusters, nullptr);

        // create a fem cluster factory
        fem::Element_Factory tClusterFactory;

        // loop over mesh clusters on set
        for( luint iCluster = 0; iCluster < tNumMeshClusters; iCluster++ )
        {
            // create a fem cluster
            mEquationObjList( iCluster ) = tClusterFactory.create_cluster( mElementType,
                                                                           mMeshClusterList( iCluster ),
                                                                           mNodes,
                                                                           this );
        }

        // get spatial dimension
        mSpaceDim = mMeshSet->get_spatial_dim();

        // bool true is master IG cell are trivial
        mIsTrivialMaster = mMeshSet->is_trivial( mtk::Master_Slave::MASTER );

        // get interpolation geometry type
        mIPGeometryType = mMeshSet->get_interpolation_cell_geometry_type();

        // get integration geometry type
        mIGGeometryType = mMeshSet->get_integration_cell_geometry_type();

        // get space interpolation order for IP cells
        // FIXME if different for different fields
        mIPSpaceInterpolationOrder = mMeshSet->get_interpolation_cell_interpolation_order();

        // get space interpolation order for IG cells
        // FIXME if different for different fields
        mIGSpaceInterpolationOrder = mMeshSet->get_integration_cell_interpolation_order();

        // create geometry interpolation rule for IP cells
        Interpolation_Rule tIPGeometryInterpolationRule( mIPGeometryType,
                                                         Interpolation_Type::LAGRANGE,
                                                         mIPSpaceInterpolationOrder,
                                                         Interpolation_Type::LAGRANGE,
                                                         mtk::Interpolation_Order::LINEAR ); // FIXME not linear?

        // create geometry interpolation rule for IG cells
        Interpolation_Rule tIGGeometryInterpolationRule( mIGGeometryType,
                                                         Interpolation_Type::LAGRANGE,
                                                         mIGSpaceInterpolationOrder,
                                                         Interpolation_Type::LAGRANGE,
                                                         mtk::Interpolation_Order::LINEAR ); // FIXME not linear?

        // switch on set type
        switch ( mElementType )
        {
            // if block-set
            case ( fem::Element_Type::BULK ):
            {
                // create a geometry interpolator for IP cells
                mMasterIPGeometryInterpolator = new Geometry_Interpolator( tIPGeometryInterpolationRule, false );

                // create a geometry interpolator for IG cells
                mMasterIGGeometryInterpolator = new Geometry_Interpolator( tIGGeometryInterpolationRule, false );

                break;
            }

            // if side-set
            case( fem::Element_Type::SIDESET ):
            {
                // create a geometry interpolator for IP cells
                mMasterIPGeometryInterpolator = new Geometry_Interpolator( tIPGeometryInterpolationRule, true );

                // create a geometry interpolator for IG cells
                mMasterIGGeometryInterpolator = new Geometry_Interpolator( tIGGeometryInterpolationRule, true );

                break;
            }

            // if double side-set
            case( fem::Element_Type::DOUBLE_SIDESET ):
            {
                // bool true is slave IG cell are trivial
                mIsTrivialSlave = mMeshSet->is_trivial( mtk::Master_Slave::SLAVE );

                // create a geometry interpolator for master and slave IP cells
                mMasterIPGeometryInterpolator = new Geometry_Interpolator( tIPGeometryInterpolationRule, true );
                mSlaveIPGeometryInterpolator  = new Geometry_Interpolator( tIPGeometryInterpolationRule, true );

                // create a geometry interpolator for master and slave IG cells
                mMasterIGGeometryInterpolator = new Geometry_Interpolator( tIGGeometryInterpolationRule, true );
                mSlaveIGGeometryInterpolator  = new Geometry_Interpolator( tIGGeometryInterpolationRule, true );

                break;
            }

            // if none of the above
            default:
            {
                MORIS_ERROR(false, "Set::Set - unknown element type");
                break;
            }
        }

        // create a unique dof type list for solver
        this->create_unique_dof_type_list();

        // create a dof type list
        this->create_dof_type_list();

        // create a dof type map
        this->create_dof_type_map();

        // create an interpolation rule
        Integration_Rule tIntegrationRule = Integration_Rule( mIGGeometryType,
                                                              Integration_Type::GAUSS,
                                                              this->get_auto_integration_order( mIGGeometryType ),
                                                              Integration_Type::GAUSS,
                                                              Integration_Order::BAR_1 ); // fixme time order

        // create an integrator
        Integrator tIntegrator( tIntegrationRule );

        // get integration points
        tIntegrator.get_points( mIntegPoints );

        // get integration weights
        tIntegrator.get_weights( mIntegWeights );
    }

//------------------------------------------------------------------------------
    Set::~Set()
    {
        // delete the equation object pointers
        for( MSI::Equation_Object* tEquationObj : mEquationObjList )
        {
            delete tEquationObj;
        }
        mEquationObjList.clear();

        // delete the master interpolation geometry interpolator pointer
        if ( mMasterIPGeometryInterpolator != nullptr )
        {
            delete mMasterIPGeometryInterpolator;
        }

        // delete the slave interpolation geometry interpolator pointer
        if ( mSlaveIPGeometryInterpolator != nullptr )
        {
            delete mSlaveIPGeometryInterpolator;
        }

        // delete the master integration geometry interpolator pointer
        if ( mMasterIGGeometryInterpolator != nullptr )
        {
            delete mMasterIGGeometryInterpolator;
        }

        // delete the slave integration geometry interpolator pointer
        if ( mSlaveIGGeometryInterpolator != nullptr )
        {
            delete mSlaveIGGeometryInterpolator;
        }

        // delete the field interpolator pointers
        this->delete_pointers();

    }

//------------------------------------------------------------------------------

    void Set::delete_pointers()
    {
        // delete the master field interpolator pointers
        for( Field_Interpolator* tMasterFieldInterpolator : mMasterFI )
        {
            delete tMasterFieldInterpolator;
        }
        mMasterFI.clear();

        // delete the slave field interpolator pointers
        for( Field_Interpolator* tSlaveFieldInterpolator : mSlaveFI )
        {
            delete tSlaveFieldInterpolator;
        }
        mSlaveFI.clear();
    }

//------------------------------------------------------------------------------
    void Set::finalize( MSI::Model_Solver_Interface * aModelSolverInterface )
    {
        // delete the field interpolator pointers
        this->delete_pointers();

        // create the field interpolators
        this->create_field_interpolators( aModelSolverInterface );

        // create dof assembly map
        this->create_dof_assembly_map();  //FIXME delete this function

        this->create_num_dof_map();

        // set field interpolators for the IWGs
        this->set_IWG_field_interpolators();

        // create the IWG info for the set
        this->create_IWG_dof_assembly_map();

        this->create_IWG_dof_assembly_map_2();

//        // set field interpolators for the properties
//        this->set_properties_field_interpolators();
//
//        // set field interpolators for the contitutive models
//        this->set_CM_field_interpolators();

    }

//------------------------------------------------------------------------------
    void Set::create_unique_dof_type_list()
    {
        // set the size of the dof type list for the set
        uint tCounter = 0;

        // loop over the IWGs
        for ( std::shared_ptr< IWG > tIWG : mIWGs )
        {
            // loop over the IWG master dof type
            for ( uint iDOF = 0; iDOF< tIWG->get_global_dof_type_list( mtk::Master_Slave::MASTER ).size(); iDOF++ )
            {
                // add the number of dof type
                tCounter += tIWG->get_global_dof_type_list( mtk::Master_Slave::MASTER )( iDOF ).size();
            }

            // loop over the IWG slave dof type
            for ( uint iDOF = 0; iDOF< tIWG->get_global_dof_type_list( mtk::Master_Slave::SLAVE ).size(); iDOF++ )
            {
                // add the number of dof type
                tCounter += tIWG->get_global_dof_type_list( mtk::Master_Slave::SLAVE )( iDOF ).size();
            }
        }

        // set max size for the dof type list
        mEqnObjDofTypeList.reserve( tCounter );

        // loop over the IWGs
        for ( std::shared_ptr< IWG > tIWG : mIWGs )
        {
            // loop over the IWG master dof type
            for ( uint iDOF = 0; iDOF < tIWG->get_global_dof_type_list( mtk::Master_Slave::MASTER ).size(); iDOF++ )
            {
                // put the master dof type in the dof type list
                mEqnObjDofTypeList.append( tIWG->get_global_dof_type_list( mtk::Master_Slave::MASTER )( iDOF ) );
            }

            // loop over the IWG slave dof type
            for ( uint iDOF = 0; iDOF < tIWG->get_global_dof_type_list( mtk::Master_Slave::SLAVE ).size(); iDOF++ )
            {
                // put the slave dof type in the dof type list
                mEqnObjDofTypeList.append( tIWG->get_global_dof_type_list( mtk::Master_Slave::SLAVE )( iDOF ) );
            }
        }

        // make the dof type list unique
        std::sort( ( mEqnObjDofTypeList.data() ).data(),
                   ( mEqnObjDofTypeList.data() ).data() + mEqnObjDofTypeList.size());
        auto last = std::unique( ( mEqnObjDofTypeList.data() ).data(),
                                 ( mEqnObjDofTypeList.data() ).data() + mEqnObjDofTypeList.size() );
        auto pos  = std::distance( ( mEqnObjDofTypeList.data() ).data(), last );
        mEqnObjDofTypeList.resize( pos );
    }

//------------------------------------------------------------------------------
    void Set::create_dof_type_list()
    {
        // set counter for master and slave dof types
        uint tMasterDofCounter = 0;
        uint tSlaveDofCounter = 0;

        // loop over the IWGs
        for ( std::shared_ptr< IWG > tIWG : mIWGs )
        {
           // add the number of master dof types
           tMasterDofCounter += tIWG->get_global_dof_type_list( mtk::Master_Slave::MASTER ).size();

           // add the number of master dof types
           tSlaveDofCounter += tIWG->get_global_dof_type_list( mtk::Master_Slave::SLAVE ).size();
        }
        // set size for master and slave dof type list
        mMasterDofTypes.resize( tMasterDofCounter );
        mSlaveDofTypes.resize( tSlaveDofCounter );

        // create a list to check if dof type is already in the list
        moris::Cell< sint > tMasterCheckList( tMasterDofCounter, -1 );
        moris::Cell< sint > tSlaveCheckList( tSlaveDofCounter, -1 );

        // reset master and slave dof type counter
        tMasterDofCounter = 0;
        tSlaveDofCounter  = 0;

        // loop over the IWGs
        for ( std::shared_ptr< IWG > tIWG : mIWGs )
        {
            // loop over the IWG active dof type
            for ( uint iDOF = 0; iDOF < tIWG->get_global_dof_type_list().size(); iDOF++ )
            {
                // check enum is not already in the list
                bool tMasterCheck = false;

                // loop over the check list
                for( uint i = 0; i < tMasterDofCounter; i++ )
                {
                    // check dof type is not already in list
                    tMasterCheck = tMasterCheck || equal_to( tMasterCheckList( i ), static_cast< uint >( tIWG->get_global_dof_type_list()( iDOF )( 0 ) ) );
                }

                // if dof type not in the list yet
                if ( !tMasterCheck )
                {
                    // add dof type to the check list
                    tMasterCheckList( tMasterDofCounter ) = static_cast< uint >( tIWG->get_global_dof_type_list()( iDOF )( 0 ) );

                    // add dof type to the dof list
                    mMasterDofTypes( tMasterDofCounter ) = tIWG->get_global_dof_type_list()( iDOF );

                    // update dof type counter
                    tMasterDofCounter++;
                }
            }

            // loop over the IWG active dof type
            for ( uint iDOF = 0; iDOF < tIWG->get_global_dof_type_list( mtk::Master_Slave::SLAVE ).size(); iDOF++ )
            {
                // check enum is not already in the list
                bool tSlaveCheck = false;

                // loop over the check list
                for( uint i = 0; i < tSlaveDofCounter; i++ )
                {
                    // check dof type is not already in list
                    tSlaveCheck = tSlaveCheck || equal_to( tSlaveCheckList( i ), static_cast< uint >( tIWG->get_global_dof_type_list( mtk::Master_Slave::SLAVE )( iDOF )( 0 ) ) );
                }

                // if dof type not in the list yet
                if ( !tSlaveCheck )
                {
                    // add dof type to the check list
                    tSlaveCheckList( tSlaveDofCounter ) = static_cast< uint >( tIWG->get_global_dof_type_list( mtk::Master_Slave::SLAVE )( iDOF )( 0 ) );

                    // add dof type to the dof list
                    mSlaveDofTypes( tSlaveDofCounter ) = tIWG->get_global_dof_type_list(mtk::Master_Slave::SLAVE )( iDOF );

                    // update dof type counter
                    tSlaveDofCounter++;
                }
            }
        }

        // resize the matser and slave dof type lists
        mMasterDofTypes.resize( tMasterDofCounter );
        mSlaveDofTypes.resize( tSlaveDofCounter );
    }

//------------------------------------------------------------------------------
    void Set::create_dof_type_map()
    {
        //MASTER------------------------------------------------------------------------------
        // get number of master dof types
        uint tMasterNumDofs = this->get_number_of_field_interpolators();

        // get maximal dof type enum
        sint tMaxEnum = -1;

        // loop over the IWGs
        for ( uint iDOF = 0; iDOF < tMasterNumDofs; iDOF++ )
        {
            // get the highest dof type enum
            tMaxEnum = std::max( tMaxEnum, static_cast< int >( mMasterDofTypes( iDOF )( 0 ) ) );
        }
        // +1 since start at 0
        tMaxEnum++;

        // set size of dof type map    // FIXME replace with map
        mMasterDofTypeMap.set_size( tMaxEnum, 1, -1 );

        // loop over dof types
        for ( uint iDOF = 0; iDOF < tMasterNumDofs; iDOF++ )
        {
            mMasterDofTypeMap( static_cast< int >( mMasterDofTypes( iDOF )( 0 ) ), 0 ) = iDOF;
        }

        //SLAVE------------------------------------------------------------------------------
        // get number of slave dof types
        uint tSlaveNumDofs =  this->get_number_of_field_interpolators( mtk::Master_Slave::SLAVE );

        // get maximal dof type enum
        tMaxEnum = -1;

        // loop over the IWGs
        for ( uint iDOF = 0; iDOF < tSlaveNumDofs; iDOF++ )
        {
            // get the highest dof type enum
            tMaxEnum = std::max( tMaxEnum, static_cast< int >( mSlaveDofTypes( iDOF )( 0 ) ) );
        }
        // +1 since start at 0
        tMaxEnum++;

        // set size of dof type map
        mSlaveDofTypeMap.set_size( tMaxEnum, 1, -1 );

        // loop over dof types
        for ( uint iDOF = 0; iDOF < tSlaveNumDofs; iDOF++ )
        {
            mSlaveDofTypeMap( static_cast< int >( mSlaveDofTypes( iDOF )( 0 ) ), 0 ) = iDOF;
        }
    }

//------------------------------------------------------------------------------
    void Set::create_dof_assembly_map( )
    {
        // get master number of field interpolators
        uint tMasterNumFI = this->get_number_of_field_interpolators();
        uint tSlaveNumFI  = this->get_number_of_field_interpolators( mtk::Master_Slave::SLAVE );

        // set size of assembly map
        mDofAssemblyMap.set_size( tMasterNumFI + tSlaveNumFI, 2, -1 );

        //MASTER------------------------------------------------------------------------------
        // init dof counter
        uint tDofCounter = 0;

        // loop on the dof type groups
        for( uint iFI = 0; iFI < tMasterNumFI; iFI++ )
        {
            // fill the assembly map with start dof counter
            mDofAssemblyMap( iFI, 0 ) = tDofCounter;

            // update dof counter
            tDofCounter += mMasterFI( iFI )->get_number_of_space_time_coefficients() - 1;

            // fill the assembly map with stop dof counter
            mDofAssemblyMap( iFI, 1 ) = tDofCounter;

            // update dof counter
            tDofCounter++;
        }

        //SLAVE------------------------------------------------------------------------------
        // loop on the dof type groups
        for( uint iFI = 0; iFI < tSlaveNumFI; iFI++ )
        {
            // fill the assembly map with start dof counter
            mDofAssemblyMap( tMasterNumFI + iFI, 0 ) = tDofCounter;

            // update dof counter
            tDofCounter += mSlaveFI( iFI )->get_number_of_space_time_coefficients() - 1;

            // fill the assembly map with stop dof counter
            mDofAssemblyMap( tMasterNumFI + iFI, 1 ) = tDofCounter;

            // update dof counter
            tDofCounter++;
        }

        // set mTotalDof
        mTotalDof = tDofCounter;
    }

//------------------------------------------------------------------------------

    void Set::create_num_dof_map( )
    {
        // get master number of field interpolators
        uint tMasterNumFI = this->get_number_of_field_interpolators();
        uint tSlaveNumFI  = this->get_number_of_field_interpolators( mtk::Master_Slave::SLAVE );

        uint tMaxNumFI = std::max( tMasterNumFI, tSlaveNumFI );
        // set size of assembly map
        mNumDofMap.set_size( tMaxNumFI, 2, -1 );

        uint mTotalDof = 0;

        // loop on the dof type groups
        for( uint iFI = 0; iFI < tMasterNumFI; iFI++ )
        {
            // get number of dofs
            uint tNumDofs = mMasterFI( iFI )->get_number_of_space_time_coefficients();

            mTotalDof += tNumDofs;

            // fill the assembly map with start dof counter
            mNumDofMap( iFI, 0 ) = tNumDofs;
        }

        // loop on the dof type groups
        for( uint iFI = 0; iFI < tSlaveNumFI; iFI++ )
        {
            // get number of dofs
            uint tNumDofs = mSlaveFI( iFI )->get_number_of_space_time_coefficients();

            mTotalDof += tNumDofs;

            // fill the assembly map with start dof counter
            mNumDofMap( iFI, 1 ) = tNumDofs;
        }
    }

//------------------------------------------------------------------------------

        void Set::set_jacobian_list_size()
        {
            // get master number of field interpolators
            uint tMasterNumFI = this->get_number_of_field_interpolators();
            uint tSlaveNumFI  = this->get_number_of_field_interpolators( mtk::Master_Slave::SLAVE );

            uint tNumFI = tMasterNumFI + tSlaveNumFI;

            mJacobians.resize( tNumFI );

            for( uint Ik = 0; Ik < tNumFI; Ik++ )
            {
                mJacobians( Ik ).resize( tNumFI );
            }

            // loop on the dof type groups
            for( uint Ik = 0; Ik < tMasterNumFI; Ik++ )
            {
                uint tNumRow = mMasterFI( Ik )->get_number_of_space_time_coefficients();

                for( uint Ii = 0; Ii < tMasterNumFI; Ii++ )
                {
                    uint tNumCols = mMasterFI( Ii )->get_number_of_space_time_coefficients();

                    mJacobians( Ik )( Ii ).set_size( tNumRow, tNumCols, 0.0 );
                }

                for( uint Ia = 0; Ia < tSlaveNumFI; Ia++ )
                {
                    uint tNumCols = mSlaveFI( Ia )->get_number_of_space_time_coefficients();

                    mJacobians( Ik )( Ia + tMasterNumFI ).set_size( tNumRow, tNumCols, 0.0 );
                }
            }

            // loop on the dof type groups
            for( uint Ik = 0; Ik < tSlaveNumFI; Ik++ )
            {
                uint tNumRow = mSlaveFI( Ik )->get_number_of_space_time_coefficients();

                for( uint Ii = 0; Ii < tMasterNumFI; Ii++ )
                {
                    uint tNumCols = mMasterFI( Ii )->get_number_of_space_time_coefficients();

                    mJacobians( Ik + tMasterNumFI + 1 )( Ii ).set_size( tNumRow, tNumCols, 0.0 );
                }

                for( uint Ia = 0; Ia < tSlaveNumFI; Ia++ )
                {
                    uint tNumCols = mSlaveFI( Ia )->get_number_of_space_time_coefficients();

                    mJacobians( Ik + tMasterNumFI + 1 )( Ia + tMasterNumFI ).set_size( tNumRow, tNumCols, 0.0 );
                }
            }
        }

//------------------------------------------------------------------------------

            void Set::set_residual_list_size()
            {
                // get master number of field interpolators
                uint tMasterNumFI = this->get_number_of_field_interpolators();
                uint tSlaveNumFI  = this->get_number_of_field_interpolators( mtk::Master_Slave::SLAVE );

                uint tNumFI = tMasterNumFI + tSlaveNumFI;

                mResiduals.resize( tNumFI );

                // loop on the dof type groups
                for( uint Ik = 0; Ik < tMasterNumFI; Ik++ )
                {
                    uint tNumRow = mMasterFI( Ik )->get_number_of_space_time_coefficients();

                    mResiduals( Ik ).set_size( tNumRow, 1, 0.0 );
                }

                // loop on the dof type groups
                for( uint Ik = 0; Ik < tSlaveNumFI; Ik++ )
                {
                    uint tNumRow = mSlaveFI( Ik )->get_number_of_space_time_coefficients();

                    mResiduals( Ik + tMasterNumFI + 1 ).set_size( tNumRow, 1, 0.0 );
                }
            }

//-----------------------------------------------------------------------------
    void Set::create_field_interpolators( MSI::Model_Solver_Interface * aModelSolverInterface )
    {
        //MASTER------------------------------------------------------------------------
        // get number of master dof types
        uint tMasterNumDOF = this->get_number_of_field_interpolators();

        // set the size of the cell of field interpolators
        mMasterFI.resize( tMasterNumDOF, nullptr );

        // loop over the master dof type groups
        for( uint iDOF = 0; iDOF < tMasterNumDOF; iDOF++ )
        {
            // get the number of time level for the dof type group
            uint tNumTimeNodes = aModelSolverInterface->get_time_levels_for_type( mMasterDofTypes( iDOF )( 0 ) );

            // create the field interpolation rule for the ith dof type group
            Interpolation_Rule tFieldInterpolationRule( mIPGeometryType,
                                                        Interpolation_Type::LAGRANGE,
                                                        mIPSpaceInterpolationOrder,
                                                        this->get_auto_time_interpolation_type( tNumTimeNodes ), // fixme
                                                        // If interpolation type CONSTANT, iInterpolation order is not used
                                                        this->get_auto_interpolation_order( tNumTimeNodes, mtk::Geometry_Type::LINE ) ); //fixme

            // create a field interpolator for the dof type group
            mMasterFI( iDOF ) = new Field_Interpolator( mMasterDofTypes( iDOF ).size(),
                                                        tFieldInterpolationRule,
                                                        mMasterIPGeometryInterpolator,
                                                        mMasterDofTypes( iDOF ) );
        }

        //SLAVE-------------------------------------------------------------------------
        // get number of slave dof types
        uint tSlaveNumDOF = this->get_number_of_field_interpolators( mtk::Master_Slave::SLAVE );

        // set the size of the cell of field interpolators
        mSlaveFI.resize( tSlaveNumDOF, nullptr );

        // loop over the master dof type groups
        for( uint iDOF = 0; iDOF < tSlaveNumDOF; iDOF++ )
        {
            // get the number of time level for the dof type group
            uint tNumTimeNodes = aModelSolverInterface->get_time_levels_for_type( mSlaveDofTypes( iDOF )( 0 ) );

            // create the field interpolation rule for the ith dof type group
            Interpolation_Rule tFieldInterpolationRule( mIPGeometryType,
                                                        Interpolation_Type::LAGRANGE,
                                                        mIPSpaceInterpolationOrder,
                                                        this->get_auto_time_interpolation_type( tNumTimeNodes ), // fixme
                                                        // If interpolation type CONSTANT, iInterpolation order is not used
                                                        this->get_auto_interpolation_order( tNumTimeNodes, mtk::Geometry_Type::LINE ) ); //fixme

            // create a field interpolator for the dof type group
            mSlaveFI( iDOF ) = new Field_Interpolator( mSlaveDofTypes( iDOF ).size(),
                                                        tFieldInterpolationRule,
                                                        mSlaveIPGeometryInterpolator,
                                                        mSlaveDofTypes( iDOF ) );
        }
    }

//------------------------------------------------------------------------------

//    void Set::create_properties( const moris::Cell< moris::Cell< fem::Property_User_Defined_Info > > & aPropertyUserDefinedInfo )
//    {
//        // build a master input property map--------------------------------------------
//
//        //get the max enum for the properties to create a property map
//        sint tMaxEnum = 0;
//
//        // number of input properties
//        uint tMasterInputProps = aPropertyUserDefinedInfo( 0 ).size();
//
//        // loop over the property types
//        for( uint iProp = 0; iProp < tMasterInputProps; iProp++ )
//        {
//            // check if enum is greater
//            tMaxEnum = std::max( tMaxEnum, static_cast< int >( aPropertyUserDefinedInfo( 0 )( iProp ).get_property_type() ) );
//        }
//        // +1 because 0 based
//        tMaxEnum++;
//
//        // set the size of the property map
//        Matrix< DDSMat > tMasterPropTypeMap( tMaxEnum, 1, -1 );
//
//        // loop over the property types
//        for( uint iProp = 0; iProp < tMasterInputProps; iProp++ )
//        {
//            // fill the property map
//            tMasterPropTypeMap( static_cast< int >( aPropertyUserDefinedInfo( 0 )( iProp ).get_property_type() ), 0 ) = iProp;
//        }
//
//        //MASTER------------------------------------------------------------------------
//        // get number of master property types
//        uint tMasterNumProp = this->get_number_of_properties();
//
//        // set size for the cell of properties
//        mMasterProperties.resize( tMasterNumProp, nullptr );
//
//        // create a property for each property type
//        for( uint iProp = 0; iProp < tMasterNumProp; iProp++ )
//        {
//            // get the index of the treated property in the info
//            sint propIndex = tMasterPropTypeMap( static_cast< int >( mMasterPropTypes( iProp ) ) );
//
//            MORIS_ASSERT( propIndex != -1, "Set::create_properties(), tMasterPropTypeMap returned -1" );
//
//            // create an property for the treated property type
//            mMasterProperties( iProp ) = new Property( mMasterPropTypes( iProp ),
//                                                       aPropertyUserDefinedInfo( 0 )( propIndex ).get_property_dof_type_list(),
//                                                       aPropertyUserDefinedInfo( 0 )( propIndex ).get_property_dv_type_list(),
//                                                       aPropertyUserDefinedInfo( 0 )( propIndex ).get_property_param_list(),
//                                                       aPropertyUserDefinedInfo( 0 )( propIndex ).get_property_valFunc(),
//                                                       aPropertyUserDefinedInfo( 0 )( propIndex ).get_property_dof_derFunc_list(),
//                                                       aPropertyUserDefinedInfo( 0 )( propIndex ).get_property_dv_derFunc_list(),
//                                                       mMasterIPGeometryInterpolator );
//        }
//
//        // build a slave input property map--------------------------------------------
//
//        //get the max enum for the properties to create a property map
//        sint tSlaveMaxEnum = 0;
//
//        // number of input properties
//        uint tSlaveInputProps = 0;
//        if ( aPropertyUserDefinedInfo.size() > 1 )
//        {
//            tSlaveInputProps = aPropertyUserDefinedInfo( 1 ).size();
//        }
//
//        // loop over the property types
//        for( uint iProp = 0; iProp < tSlaveInputProps; iProp++ )
//        {
//            // check if enum is greater
//            tSlaveMaxEnum = std::max( tSlaveMaxEnum, static_cast< int >( aPropertyUserDefinedInfo( 1 )( iProp ).get_property_type() ) );
//        }
//        // +1 because 0 based
//        tSlaveMaxEnum++;
//
//        // set the size of the property map
//        Matrix< DDSMat > tSlavePropTypeMap( tSlaveMaxEnum, 1, -1 );
//
//        // loop over the property types
//        for( uint iProp = 0; iProp < tSlaveInputProps; iProp++ )
//        {
//            // fill the property map
//            tSlavePropTypeMap( static_cast< int >( aPropertyUserDefinedInfo( 1 )( iProp ).get_property_type() ), 0 ) = iProp;
//        }
//
//        //SLAVE-------------------------------------------------------------------------
//        // get number of slave property types
//        uint tSlaveNumProp = this->get_number_of_properties( mtk::Master_Slave::SLAVE );
//
//        // set size for the cell of properties
//        mSlaveProperties.resize( tSlaveNumProp, nullptr );
//
//        // create a property for each property type
//        for( uint iProp = 0; iProp < tSlaveNumProp; iProp++ )
//        {
//            // get the index of the treated property in the info
//            uint propIndex = tSlavePropTypeMap( static_cast< int >( mSlavePropTypes( iProp ) ) );
//
//            // create an property for the treated property type
//            mSlaveProperties( iProp ) = new Property( mSlavePropTypes( iProp ),
//                                                      aPropertyUserDefinedInfo( 1 )( propIndex ).get_property_dof_type_list(),
//                                                      aPropertyUserDefinedInfo( 1 )( propIndex ).get_property_dv_type_list(),
//                                                      aPropertyUserDefinedInfo( 1 )( propIndex ).get_property_param_list(),
//                                                      aPropertyUserDefinedInfo( 1 )( propIndex ).get_property_valFunc(),
//                                                      aPropertyUserDefinedInfo( 1 )( propIndex ).get_property_dof_derFunc_list(),
//                                                      aPropertyUserDefinedInfo( 1 )( propIndex ).get_property_dv_derFunc_list(),
//                                                      mSlaveIPGeometryInterpolator );
//        }
//    }
//
////------------------------------------------------------------------------------
//    void Set::create_constitutive_models( const moris::Cell< moris::Cell< fem::Constitutive_User_Defined_Info > > & aConstitutiveUserDefinedInfo )
//    {
//        // build a master input constitutive map--------------------------------------------------
//        //get the max enum for the constitutive type to create a constitutive map
//        sint tMasterMaxEnum = 0;
//
//        // number of input constitutive
//        uint tMasterInputCMs = aConstitutiveUserDefinedInfo( 0 ).size();
//
//        // loop over the constitutive types
//        for( uint iCM = 0; iCM < tMasterInputCMs; iCM++ )
//        {
//            // check if enum is greater
//            tMasterMaxEnum = std::max( tMasterMaxEnum, static_cast< int >( aConstitutiveUserDefinedInfo( 0 )( iCM ).get_constitutive_type() ) );
//        }
//        // +1 because 0 based
//        tMasterMaxEnum++;
//
//        // set the size of the constitutive map
//        Matrix< DDSMat > tMasterConstitutiveTypeMap( tMasterMaxEnum, 1, -1 );
//
//        // loop over the constitutive types
//        for( uint iCM = 0; iCM < tMasterInputCMs; iCM++ )
//        {
//            // fill the constitutive map
//            tMasterConstitutiveTypeMap( static_cast< int >( aConstitutiveUserDefinedInfo( 0 )( iCM ).get_constitutive_type() ), 0 ) = iCM;
//        }
//
//        // create a constitutive model factory
//        fem::CM_Factory tCMFactory;
//
//        //MASTER------------------------------------------------------------------------
//        // get number of master constitutive types
//        uint tMasterNumCM = this->get_number_of_constitutive_models();
//
//        // set size for the cell of constitutive models
//        mMasterCM.resize( tMasterNumCM, nullptr );
//
//        // create a property for each constitutive type
//        for( uint iCM = 0; iCM < tMasterNumCM; iCM++ )
//        {
//            // get the index of the treated constitutive model in the info
//            uint tCMIndex = tMasterConstitutiveTypeMap( static_cast< int >( mMasterConstitutiveTypes( iCM ) ) );
//
//            // create an property for the treated property type
//            mMasterCM( iCM ) =  tCMFactory.create_CM( mMasterConstitutiveTypes( iCM ) );
//
//            // set space dim
//            mMasterCM( iCM )->set_space_dim( mSpaceDim );
//
//            // set dof types
//            mMasterCM( iCM )->set_dof_type_list( aConstitutiveUserDefinedInfo( 0 )( tCMIndex ).get_constitutive_dof_type_list() );
//
//            // set dv types
//            mMasterCM( iCM )->set_dv_type_list( aConstitutiveUserDefinedInfo( 0 )( tCMIndex ).get_constitutive_dv_type_list() );
//
//            // set property type
//            mMasterCM( iCM )->set_property_type_list( aConstitutiveUserDefinedInfo( 0 )( tCMIndex ).get_constitutive_property_type_list() );
//        }
//
//        // build a master input constitutive map--------------------------------------------------
//        //get the max enum for the constitutive type to create a constitutive map
//        sint tSlaveMaxEnum = 0;
//
//        // number of input constitutive
//        uint tSlaveInputCMs = 0;
//        if( aConstitutiveUserDefinedInfo.size() > 1 )
//        {
//            tSlaveInputCMs = aConstitutiveUserDefinedInfo( 1 ).size();
//        }
//
//        // loop over the constitutive types
//        for( uint iCM = 0; iCM < tSlaveInputCMs; iCM++ )
//        {
//            // check if enum is greater
//            tSlaveMaxEnum = std::max( tSlaveMaxEnum, static_cast< int >( aConstitutiveUserDefinedInfo( 1 )( iCM ).get_constitutive_type() ) );
//        }
//        // +1 because 0 based
//        tSlaveMaxEnum++;
//
//        // set the size of the constitutive map
//        Matrix< DDSMat > tSlaveConstitutiveTypeMap( tSlaveMaxEnum, 1, -1 );
//
//        // loop over the constitutive types
//        for( uint iCM = 0; iCM < tSlaveInputCMs; iCM++ )
//        {
//            // fill the constitutive map
//            tSlaveConstitutiveTypeMap( static_cast< int >( aConstitutiveUserDefinedInfo( 1 )( iCM ).get_constitutive_type() ), 0 ) = iCM;
//        }
//
//        //SLAVE-------------------------------------------------------------------------
//        // get number of slave constitutive types
//        uint tSlaveNumCM = this->get_number_of_constitutive_models( mtk::Master_Slave::SLAVE );
//
//        // set size for the cell of constitutive models
//        mSlaveCM.resize( tSlaveNumCM, nullptr );
//
//        // create a constitutive model for each constitutive type
//        for( uint iCM = 0; iCM < tSlaveNumCM; iCM++ )
//        {
//            // get the index of the treated property in the info
//            uint tCMIndex = tSlaveConstitutiveTypeMap( static_cast< int >( mSlaveConstitutiveTypes( iCM ) ) );
//
//            // create an property for the treated property type
//            mSlaveCM( iCM ) = tCMFactory.create_CM( mSlaveConstitutiveTypes( iCM ) );
//
//            // set space dim
//            mSlaveCM( iCM )->set_space_dim( mSpaceDim );
//
//            // set dof types
//            mSlaveCM( iCM )->set_dof_type_list( aConstitutiveUserDefinedInfo( 1 )( tCMIndex ).get_constitutive_dof_type_list() );
//
//            // set dv types
//            mSlaveCM( iCM )->set_dv_type_list( aConstitutiveUserDefinedInfo( 1 )( tCMIndex ).get_constitutive_dv_type_list() );
//
//            // set property type
//            mSlaveCM( iCM )->set_property_type_list( aConstitutiveUserDefinedInfo( 1 )( tCMIndex ).get_constitutive_property_type_list() );
//        }
//    }
//
////------------------------------------------------------------------------------
//    void Set::set_IWG_properties()
//    {
//         // loop over the IWGs
//         for ( IWG * tIWG : mIWGs )
//         {
//             //MASTER------------------------------------------------------------------------
//             // get the number of master property type
//             uint tMasterIWGNumProp = tIWG->get_global_property_type_list().size();
//
//             // set size for cell of master properties
//             moris::Cell< fem::Property* > tIWGProp( tMasterIWGNumProp, nullptr );
//
//             // loop over the property type
//             for( uint iProp = 0; iProp < tMasterIWGNumProp; iProp++ )
//             {
//                 // get the property index
//                 uint propIndex = mMasterPropTypeMap( static_cast< int >( tIWG->get_global_property_type_list()( iProp ) ) );
//
//                 // collect the properties that are active for the IWG
//                 tIWGProp( iProp ) = mMasterProperties( propIndex );
//             }
//
//             // set the IWG properties
//             tIWG->set_properties( tIWGProp );
//
//             //SLAVE-------------------------------------------------------------------------
//             // get the number of slave property type
//             uint tSlaveIWGNumProp = tIWG->get_global_property_type_list( mtk::Master_Slave::SLAVE ).size();
//
//             // reset size for cell of slave properties
//             tIWGProp.resize( tSlaveIWGNumProp, nullptr );
//
//             // loop over the property type
//             for( uint iProp = 0; iProp < tSlaveIWGNumProp; iProp++ )
//             {
//                 // get the property index
//                 uint propIndex = mSlavePropTypeMap( static_cast< int >( tIWG->get_global_property_type_list( mtk::Master_Slave::SLAVE )( iProp ) ) );
//
//                 // collect the properties that are active for the IWG
//                 tIWGProp( iProp ) = mSlaveProperties( propIndex );
//             }
//
//             // set the IWG properties
//             tIWG->set_properties( tIWGProp, mtk::Master_Slave::SLAVE );
//        }
//    }

//------------------------------------------------------------------------------

    void Set::set_IWG_field_interpolators()
    {
        // loop over the IWGs
        for ( std::shared_ptr< IWG > tIWG : mIWGs )
        {
            //MASTER------------------------------------------------------------------------
            // get number of master dof type for IWG
            uint tMasterIWGNumDof = tIWG->get_global_dof_type_list().size();

            // set size for cell of field interpolators
            Cell< Field_Interpolator* > tIWGFI( tMasterIWGNumDof, nullptr );

            // select associated active interpolators
            for( uint iDOF = 0; iDOF < tMasterIWGNumDof; iDOF++ )
            {
                // find the index of active dof type in the list of element dof type
                uint tIWGDofIndex = mMasterDofTypeMap( static_cast< int >( tIWG->get_global_dof_type_list()( iDOF )( 0 ) ), 0 );

                // select the corresponding interpolator
                tIWGFI( iDOF ) = mMasterFI( tIWGDofIndex );
            }

            // set IWG field interpolators
            tIWG->set_dof_field_interpolators( tIWGFI );

            //SLAVE------------------------------------------------------------------------
            // get number of slave dof type for IWG
            uint tSlaveIWGNumDof = tIWG->get_global_dof_type_list( mtk::Master_Slave::SLAVE ).size();

            // reset size for cell of field interpolators
            tIWGFI.resize( tSlaveIWGNumDof, nullptr );

            // select associated active interpolators
            for( uint iDOF = 0; iDOF < tSlaveIWGNumDof; iDOF++ )
            {
                // find the index of active dof type in the list of element dof type
                uint tIWGDofIndex = mSlaveDofTypeMap( static_cast< int >( tIWG->get_global_dof_type_list( mtk::Master_Slave::SLAVE )( iDOF )( 0 ) ), 0 );

                // select the corresponding interpolator
                tIWGFI( iDOF ) = mSlaveFI( tIWGDofIndex );
            }

            // set IWG field interpolators
            tIWG->set_dof_field_interpolators( tIWGFI, mtk::Master_Slave::SLAVE );
        }
    }

//------------------------------------------------------------------------------
    void Set::set_IWG_geometry_interpolators()
    {
        // loop over the IWGs
        for ( std::shared_ptr< IWG > tIWG : mIWGs )
        {
            //MASTER------------------------------------------------------------------------
            // set IWG geometry interpolators
            tIWG->set_geometry_interpolator( mMasterIPGeometryInterpolator );

            //SLAVE------------------------------------------------------------------------
            // set IWG field interpolators
            if( mElementType == fem::Element_Type::DOUBLE_SIDESET )
            {
                tIWG->set_geometry_interpolator( mSlaveIPGeometryInterpolator, mtk::Master_Slave::SLAVE );
            }
        }
    }

//------------------------------------------------------------------------------
    void Set::create_IWG_dof_assembly_map()
    {
        // get number of IWGs
        uint tNumOfIWG = this->get_number_of_IWGs();

        // set info size
        mIWGResDofAssemblyMap.resize( tNumOfIWG );
        mIWGJacDofAssemblyMap.resize( tNumOfIWG );

        // loop over the IWGs
        for ( uint iIWG = 0; iIWG < tNumOfIWG; iIWG++ )
        {
            // get IWG
            std::shared_ptr< IWG > tIWG = mIWGs( iIWG );

            // get the number of dof type
            uint tMasterIWGNumDof = tIWG->get_global_dof_type_list().size();
            uint tSlaveIWGNumDof  = tIWG->get_global_dof_type_list( mtk::Master_Slave::SLAVE ).size();
            uint tTotalIWGNumDof  = tMasterIWGNumDof + tSlaveIWGNumDof;

            // set size
            mIWGResDofAssemblyMap( iIWG ).set_size( 1, 2 );
            mIWGJacDofAssemblyMap( iIWG ).set_size( tTotalIWGNumDof, 2 );

            // select associated active interpolators
            Matrix< DDUMat > tStartJDof( tTotalIWGNumDof, 1 );
            Matrix< DDUMat > tStopJDof ( tTotalIWGNumDof, 1 );

            // loop over the master and slave dof types
            for( uint iDOF = 0; iDOF < tTotalIWGNumDof; iDOF++ )
            {
                // init index of the dof type in the list of dof types
                uint tIWGDofIndex;

                // if master dof type
                if( iDOF < tMasterIWGNumDof )
                {
                    // find the index of dof type in the list of master dof type
                    tIWGDofIndex = mMasterDofTypeMap( static_cast< int >( tIWG->get_global_dof_type_list()( iDOF )( 0 ) ), 0 );
                }
                // if slave dof type
                else
                {
                    // find the index of dof type in the list of slave dof type
                    tIWGDofIndex = this->get_number_of_field_interpolators()
                                 + mSlaveDofTypeMap( static_cast< int >( tIWG->get_global_dof_type_list( mtk::Master_Slave::SLAVE )( iDOF - tMasterIWGNumDof )( 0 ) ), 0 );
                }
                // get the assembly indices for the dof type
                tStartJDof( iDOF ) = mDofAssemblyMap( tIWGDofIndex, 0 );
                tStopJDof( iDOF )  = mDofAssemblyMap( tIWGDofIndex, 1 );
            }

            // find the index of residual dof type in the list of master dof types
            uint tMasterResDofIndex = mMasterDofTypeMap( static_cast< int >( tIWG->get_residual_dof_type()( 0 ) ), 0 );

            // fill the start and stop indices for the master residual dofs
            mIWGResDofAssemblyMap( iIWG )( 0, 0 ) = mDofAssemblyMap( tMasterResDofIndex, 0 );
            mIWGResDofAssemblyMap( iIWG )( 0, 1 ) = mDofAssemblyMap( tMasterResDofIndex, 1 );

            // fill the map
            mIWGJacDofAssemblyMap( iIWG ).get_column( 0 ) = tStartJDof.matrix_data();
            mIWGJacDofAssemblyMap( iIWG ).get_column( 1 ) = tStopJDof.matrix_data();

            // if there is a slave
            if( tSlaveIWGNumDof > 0 )
            {
                // resize the residual dof assembly map to account for the slave
                mIWGResDofAssemblyMap( iIWG ).resize( 2, 2 );

                // find the index of residual dof type in the list of slave dof types
                uint tSlaveResDofIndex = this->get_number_of_field_interpolators()
                                       + mSlaveDofTypeMap( static_cast< int >( tIWG->get_residual_dof_type()( 0 ) ), 0 );

                // fill the start and stop indices for the slave residual dofs
                mIWGResDofAssemblyMap( iIWG )( 1, 0 ) = mDofAssemblyMap( tSlaveResDofIndex, 0 );
                mIWGResDofAssemblyMap( iIWG )( 1, 1 ) = mDofAssemblyMap( tSlaveResDofIndex, 1 );
            }
        }
    }

//------------------------------------------------------------------------------

    void Set::create_IWG_dof_assembly_map_2()
    {
        // get number of IWGs
        uint tNumOfIWG = this->get_number_of_IWGs();

        // set info size
        mIWGResDofAssemblyMap_2.resize( tNumOfIWG );
        mIWGJacDofAssemblyMap_2.resize( tNumOfIWG );

        // loop over the IWGs
        for ( uint iIWG = 0; iIWG < tNumOfIWG; iIWG++ )
        {
            // get IWG
            std::shared_ptr< IWG > tIWG = mIWGs( iIWG );

            // get the number of dof type
            uint tMasterIWGNumDof = tIWG->get_global_dof_type_list().size();
            uint tSlaveIWGNumDof  = tIWG->get_global_dof_type_list( mtk::Master_Slave::SLAVE ).size();
            uint tTotalIWGNumDof  = tMasterIWGNumDof + tSlaveIWGNumDof;

            // set size
            mIWGResDofAssemblyMap_2( iIWG ).set_size( 1, 1 );
            mIWGJacDofAssemblyMap_2( iIWG ).set_size( tTotalIWGNumDof, 1 );

            // select associated active interpolators
            Matrix< DDUMat > tPosDof( tTotalIWGNumDof, 1 );
//            Matrix< DDUMat > tStopJDof ( tTotalIWGNumDof, 1 );

            // loop over the master and slave dof types
            for( uint iDOF = 0; iDOF < tTotalIWGNumDof; iDOF++ )
            {
                // init index of the dof type in the list of dof types
                uint tIWGDofIndex;

                // if master dof type
                if( iDOF < tMasterIWGNumDof )
                {
                    // find the index of dof type in the list of master dof type
                    tIWGDofIndex = mMasterDofTypeMap( static_cast< int >( tIWG->get_global_dof_type_list()( iDOF )( 0 ) ), 0 );
                }
                // if slave dof type
//                else
//                {
//                    // find the index of dof type in the list of slave dof type
//                    tIWGDofIndex = this->get_number_of_field_interpolators()
//                                 + mSlaveDofTypeMap( static_cast< int >( tIWG->get_global_dof_type_list( mtk::Master_Slave::SLAVE )( iDOF - tMasterIWGNumDof )( 0 ) ), 0 );
//                }
                // get the assembly indices for the dof type
                tPosDof( iDOF ) = tIWGDofIndex;
//                tStopJDof( iDOF )  = mDofAssemblyMap( tIWGDofIndex, 1 );
            }

            // find the index of residual dof type in the list of master dof types
//            uint tMasterResDofIndex = mMasterDofTypeMap( static_cast< int >( tIWG->get_residual_dof_type()( 0 ) ), 0 );

            // fill the start and stop indices for the master residual dofs
//            mIWGResDofAssemblyMap_2( iIWG )( 0, 0 ) = mDofAssemblyMap( tMasterResDofIndex, 0 );
//            mIWGResDofAssemblyMap_2( iIWG )( 0, 1 ) = mDofAssemblyMap( tMasterResDofIndex, 1 );

            // fill the map
            mIWGJacDofAssemblyMap_2( iIWG ).get_column( 0 ) = tPosDof.matrix_data();
//            mIWGJacDofAssemblyMap_2( iIWG ).get_column( 1 ) = tStopJDof.matrix_data();

//            // if there is a slave
//            if( tSlaveIWGNumDof > 0 )
//            {
//                // resize the residual dof assembly map to account for the slave
//                mIWGResDofAssemblyMap( iIWG ).resize( 2, 2 );
//
//                // find the index of residual dof type in the list of slave dof types
//                uint tSlaveResDofIndex = this->get_number_of_field_interpolators()
//                                       + mSlaveDofTypeMap( static_cast< int >( tIWG->get_residual_dof_type()( 0 ) ), 0 );
//
//                // fill the start and stop indices for the slave residual dofs
//                mIWGResDofAssemblyMap( iIWG )( 1, 0 ) = mDofAssemblyMap( tSlaveResDofIndex, 0 );
//                mIWGResDofAssemblyMap( iIWG )( 1, 1 ) = mDofAssemblyMap( tSlaveResDofIndex, 1 );
//            }
        }
    }

//------------------------------------------------------------------------------
//    void Set::set_properties_field_interpolators()
//    {
//        //MASTER-----------------------------------------------------------------------------
//        // loop over the properties
//        for ( Property* tProperty: mMasterProperties )
//        {
//            // get number of master dof type for property
//            uint tMasterNumDof = tProperty->get_dof_type_list().size();
//
//            // set size for cell of field interpolator for master properties
//            moris::Cell< Field_Interpolator* > tPropMasterFI( tMasterNumDof, nullptr );
//
//            // select associated active interpolators
//            for( uint iDOF = 0; iDOF < tMasterNumDof; iDOF++ )
//            {
//                // find the index of the dof type in the list of set dof type
//                uint tDofIndex = mMasterDofTypeMap( static_cast< int >( tProperty->get_dof_type_list()( iDOF )( 0 ) ) );
//
//                // add the corresponding field interpolator to the property list
//                tPropMasterFI( iDOF ) = mMasterFI( tDofIndex );
//            }
//                // set field interpolators for the property
//                tProperty->set_dof_field_interpolators( tPropMasterFI );
//        }
//
//        //SLAVE------------------------------------------------------------------------------
//
//        // loop over the properties
//        for ( Property* tProperty: mSlaveProperties )
//        {
//            // get number of master dof type for property
//            uint tSlaveNumDof = tProperty->get_dof_type_list().size();
//
//            // set size for cell of field interpolator for slave properties
//            moris::Cell< Field_Interpolator* > tPropSlaveFI( tSlaveNumDof, nullptr );
//
//            // select associated active interpolators
//            for( uint iDOF = 0; iDOF < tSlaveNumDof; iDOF++ )
//            {
//                // find the index of the dof type in the list of set dof type
//                uint tDofIndex = mSlaveDofTypeMap( static_cast< int >( tProperty->get_dof_type_list()( iDOF )( 0 ) ) );
//
//                // add the corresponding field interpolator to the property list
//                tPropSlaveFI( iDOF ) = mSlaveFI( tDofIndex );
//            }
//                // set field interpolators for the property
//                tProperty->set_dof_field_interpolators( tPropSlaveFI );
//        }
//    }

//------------------------------------------------------------------------------

    Field_Interpolator* Set::get_dof_type_field_interpolators ( enum MSI::Dof_Type aDofType,
                                                                mtk::Master_Slave  aIsMaster )
    {
        switch( aIsMaster )
        {
            case ( mtk::Master_Slave::MASTER ):
            {
                // get the index of the dof type in the master dof type list
                uint tDofIndex = mMasterDofTypeMap( static_cast< int >( aDofType ) );

                // return corresponding master field interpolator
                return mMasterFI( tDofIndex );
            }
            case ( mtk::Master_Slave::SLAVE ):
            {
                // get the index of the dof type in the slave dof type list
                uint tDofIndex = mSlaveDofTypeMap( static_cast< int >( aDofType ) );

                // return corresponding slave field interpolator
                return mSlaveFI( tDofIndex );
            }
            default:
            {
                MORIS_ERROR( false, "Set::get_dof_type_field_interpolators - can only be MASTER or SLAVE.");
                return nullptr;
            }
        }
    }

//------------------------------------------------------------------------------
    void Set::initialize_mJacobian()
    {
        if ( !mJacobianExist )
        {
            uint tTotalDof = this->get_total_number_of_dofs();
            mJacobian.set_size( tTotalDof, tTotalDof, 0.0 );

            this->set_jacobian_list_size();

            mJacobianExist = true;
        }
        else
        {
            MORIS_ASSERT( mJacobian.numel() > 0, "Set::initialize_mJacobian() - Jacobian not properly initialized.");

            mJacobian.fill( 0.0 );

            for( uint Ik = 0; Ik < mJacobians.size(); Ik++ )
            {
                for( uint Ii = 0; Ii < mJacobians( Ik ).size(); Ii++ )
                {
                     mJacobians( Ik )( Ii ).fill( 0.0 );
                }
            }
        }
    }

//------------------------------------------------------------------------------
    void Set::initialize_mResidual()
    {
        if ( !mResidualExist )
        {
            mResidual.set_size( this->get_total_number_of_dofs(), 1, 0.0 );

            this->set_residual_list_size();

            mResidualExist = true;
        }
        else
        {
            MORIS_ASSERT( mResidual.numel() > 0, "Set::initialize_mResidual() - Residual not properly initialized.");

            mResidual.fill( 0.0 );

            for( uint Ik = 0; Ik < mResiduals.size(); Ik++ )
            {
                mResiduals( Ik ).fill( 0.0 );
            }
        }
    }

//------------------------------------------------------------------------------
    mtk::Interpolation_Order Set::get_auto_interpolation_order( const moris::uint        aNumVertices,
                                                                const mtk::Geometry_Type aGeometryType )
    {
        switch( aGeometryType )
        {
            case( mtk::Geometry_Type::LINE ) :
                switch( aNumVertices )
                {
                   case( 1 ) :
                       return mtk::Interpolation_Order::UNDEFINED;
                       break;
                   case( 2 ) :
                       return mtk::Interpolation_Order::LINEAR;
                       break;

                   case( 3 ) :
                       return mtk::Interpolation_Order::QUADRATIC;
                       break;

                   default :
                       MORIS_ERROR( false, " Element::get_auto_interpolation_order - not defined for LINE and number of vertices. ");
                       return mtk::Interpolation_Order::UNDEFINED;
                       break;
                }

            case( mtk::Geometry_Type::QUAD ) :
                switch( aNumVertices )
                {
                    case( 4 ) :
                        return mtk::Interpolation_Order::LINEAR;
                        break;

                    case( 8 ) :
                        return mtk::Interpolation_Order::SERENDIPITY;
                        break;

                    case( 9 ) :
                        return mtk::Interpolation_Order::QUADRATIC;
                        break;

                    case( 16 ) :
                        return mtk::Interpolation_Order::CUBIC;
                        break;

                    default :
                        MORIS_ERROR( false, " Element::get_auto_interpolation_order - not defined for QUAD and number of vertices. ");
                        return mtk::Interpolation_Order::UNDEFINED;
                        break;
                }

            case( mtk::Geometry_Type::HEX ) :
                switch( aNumVertices )
                {
                    case( 8 ) :
                        return mtk::Interpolation_Order::LINEAR;
                        break;

                    case( 20 ) :
                        return mtk::Interpolation_Order::SERENDIPITY;
                        break;

                    case( 27 ) :
                        return mtk::Interpolation_Order::QUADRATIC;
                        break;

                    case( 64 ) :
                        return mtk::Interpolation_Order::CUBIC;
                        break;

                    default :
                        MORIS_ERROR( false, " Element::get_auto_interpolation_order - not defined for HEX and number of vertices. ");
                        return mtk::Interpolation_Order::UNDEFINED;
                        break;
                }

                case( mtk::Geometry_Type::TET ) :
                switch( aNumVertices )
                {
                    case( 4 ) :
                        return mtk::Interpolation_Order::LINEAR;
                        break;

                    case( 10 ) :
                        return mtk::Interpolation_Order::QUADRATIC;
                        break;

                    case( 20 ) :
                        return mtk::Interpolation_Order::CUBIC;
                        break;

                    default :
                        MORIS_ERROR( false, " Element::get_auto_interpolation_order - not defined for TET and number of vertices. ");
                        return mtk::Interpolation_Order::UNDEFINED;
                        break;
                }

            default :
                MORIS_ERROR( false, " Element::get_auto_interpolation_order - not defined for this geometry type. ");
                return mtk::Interpolation_Order::UNDEFINED;
                break;
        }
    }

//------------------------------------------------------------------------------
    fem::Interpolation_Type Set::get_auto_time_interpolation_type( const moris::uint aNumVertices )
    {
        switch( aNumVertices )
        {
          case( 1 ) :
              return Interpolation_Type::CONSTANT;
              break;
          case( 2 ) :
          case( 3 ) :
          case( 4 ) :
              return Interpolation_Type::LAGRANGE;
              break;
          default :
              MORIS_ERROR( false, " Element::get_auto_time_interpolation_type - not defined this number of time vertices. ");
              return Interpolation_Type::UNDEFINED;
              break;
        }
    }

//------------------------------------------------------------------------------
    fem::Integration_Order Set::get_auto_integration_order( const mtk::Geometry_Type aGeometryType )
    {
        switch( aGeometryType )
        {
            case( mtk::Geometry_Type::LINE ) :
                return fem::Integration_Order::BAR_3;
                break;

            case( mtk::Geometry_Type::QUAD ) :
                 return fem::Integration_Order::QUAD_2x2;         //FIXME
                 break;

            case( mtk::Geometry_Type::HEX ) :
                return fem::Integration_Order::HEX_3x3x3;
                break;

            case( mtk::Geometry_Type::TRI ) :
                return fem::Integration_Order::TRI_6;
                break;

            case( mtk::Geometry_Type::TET ) :
                return fem::Integration_Order::TET_5;
                break;

            default :
                MORIS_ERROR( false, " Element::get_auto_integration_order - not defined for this geometry type. ");
                return Integration_Order::UNDEFINED;
                break;
        }
    }

//------------------------------------------------------------------------------

    } /* namespace fem */
} /* namespace moris */
