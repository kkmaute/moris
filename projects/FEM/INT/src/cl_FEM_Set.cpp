/*
 * cl_FEM_Set.cpp
 *
 *  Created on: Apr 11, 2019
 *      Author: schmidt/noel
 */
#include <iostream>

#include "cl_MSI_Model_Solver_Interface.hpp"     //FEM/MSI/src
#include "cl_MSI_Solver_Interface.hpp"           //FEM/MSI/src
#include "cl_FEM_Model.hpp"                      //FEM/INT/src
#include "cl_FEM_Set.hpp"                        //FEM/INT/src
#include "cl_FEM_Set_User_Info.hpp"              //FEM/INT/src
#include "cl_FEM_Element_Factory.hpp"            //FEM/INT/src
#include "cl_FEM_Integrator.hpp"                 //FEM/INT/src
#include "cl_FEM_Field_Interpolator_Manager.hpp" //FEM/INT/src
#include "cl_FEM_Interpolation_Element.hpp"      //FEM/INT/src
#include "cl_FEM_Cluster.hpp"                    //FEM/INT/src
#include "cl_MTK_Set.hpp"                        //FEM/INT/src
#include "fn_equal_to.hpp"

namespace moris
{
    namespace fem
    {
        //------------------------------------------------------------------------------

        Set::Set(
                fem::FEM_Model                  * aFemModel,
                moris::mtk::Set                 * aMeshSet,
                const fem::Set_User_Info        & aSetInfo,
                const moris::Cell< Node_Base* > & aIPNodes )
        : mFemModel( aFemModel ),
          mMeshSet( aMeshSet ),
          mNodes( aIPNodes ),
          mIWGs( aSetInfo.get_IWGs() ),
          mIQIs( aSetInfo.get_IQIs() ),
          mTimeContinuity( aSetInfo.get_time_continuity() ),
          mIsAnalyticalSA( aSetInfo.get_is_analytical_sensitivity_analysis() ),
          mFDSchemeForSA( aSetInfo.get_finite_difference_scheme_for_sensitivity_analysis() ),
          mFDPerturbation( aSetInfo.get_finite_difference_perturbation_size() )
        {
            // get the set type (BULK, SIDESET, DOUBLE_SIDESET, TIME_SIDESET)
            this->determine_set_type();

            // loop over the IWGs on the set
            for(  std::shared_ptr< IWG > tIWG : mIWGs )
            {
                // set the fem set pointer to the IWG
                tIWG->set_set_pointer( this );
            }

            // loop over the IQIs on the set
            for(  std::shared_ptr< IQI > tIQI : mIQIs )
            {
                // set the fem set pointer to the IQI
                tIQI->set_set_pointer( this );
            }

            // get mesh clusters on set
            mMeshClusterList = mMeshSet->get_clusters_on_set();

            // get number of mesh clusters on set
            uint tNumMeshClusters = mMeshClusterList.size();

            // set size for the equation objects list
            mEquationObjList.resize( tNumMeshClusters, nullptr );

            // create a fem cluster factory
            fem::Element_Factory tClusterFactory;

            // loop over mesh clusters on set
            for( luint iCluster = 0; iCluster < tNumMeshClusters; iCluster++ )
            {
                // init list of pointers to IP mesh cell
                moris::Cell< const mtk::Cell * > tInterpolationCell;

                // switch on set type
                switch ( mElementType )
                {
                    // if bulk or sideset
                    case fem::Element_Type::BULK:
                    case fem::Element_Type::SIDESET:
                    case fem::Element_Type::TIME_SIDESET:
                    case fem::Element_Type::TIME_BOUNDARY:
                    {
                        tInterpolationCell.resize( 1, &mMeshClusterList( iCluster )->get_interpolation_cell() );
                        break;
                    }
                    // if double sideset
                    case fem::Element_Type::DOUBLE_SIDESET:
                    {                                         
                        tInterpolationCell.resize( 2 );
                        tInterpolationCell( 0 ) = &mMeshClusterList( iCluster )->get_interpolation_cell( mtk::Master_Slave::MASTER );
                        tInterpolationCell( 1 ) = &mMeshClusterList( iCluster )->get_interpolation_cell( mtk::Master_Slave::SLAVE );
                        break;
                    }                                               
                    // if none of the above
                    default:
                    {
                        MORIS_ERROR(false, "Set::Set - unknown element type");
                    }
                }

                // create an interpolation element
                mEquationObjList( iCluster ) = new fem::Interpolation_Element( mElementType, tInterpolationCell, mNodes, this );

                // create a fem cluster
                std::shared_ptr< fem::Cluster > tCluster = std::make_shared< fem::Cluster >(
                        mElementType,
                        mMeshClusterList( iCluster ),
                        this,
                        mEquationObjList( iCluster ) );

                // set the cluster to the interpolation element
                reinterpret_cast< fem::Interpolation_Element* >( mEquationObjList( iCluster ) )->set_cluster( tCluster, 0 );
            }

            // geometry and interpolation info----------------------------------------------
            //------------------------------------------------------------------------------
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

            // dof and dv dependencies info-------------------------------------------------
            //------------------------------------------------------------------------------
            // create a unique dof and dv type lists for solver
            this->create_unique_dof_and_dv_type_lists();

            // create a unique dof and dv type maps
            this->create_unique_dof_and_dv_type_maps();

            // create a dof and dv type lists
            this->create_dof_and_dv_type_lists();

            // create a dof and dv type maps
            this->create_dof_and_dv_type_maps();

            // integration info-------------------------------------------------------------
            //------------------------------------------------------------------------------

            // set default time integration order
            mtk::Geometry_Type tTimeGeometryType         = mtk::Geometry_Type::LINE;
            fem::Integration_Order tTimeIntegrationOrder = fem::Integration_Order::BAR_2;

            // if a time set
            if( mTimeContinuity || mTimeBoundary )
            {
                tTimeGeometryType     = mtk::Geometry_Type::POINT;
                tTimeIntegrationOrder = fem::Integration_Order::POINT;
            }

            // create an integration rule
            Integration_Rule tIntegrationRule = Integration_Rule( mIGGeometryType,
                    Integration_Type::GAUSS,
                    this->get_auto_integration_order(
                            mElementType,
                            mIGGeometryType,
                            mIPSpaceInterpolationOrder ),
                            tTimeGeometryType,
                            Integration_Type::GAUSS,
                            tTimeIntegrationOrder );

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

            // delete the field interpolator pointers
            this->delete_pointers();
        }

        //------------------------------------------------------------------------------

        void Set::initialize_set( const bool aIsResidual )
        {
            if ( !mIsEmptySet )    //FIXME this flag is a hack. find better solution
            {
                mIsResidual = aIsResidual;

                this->create_residual_dof_assembly_map();

                this->create_dof_assembly_map( aIsResidual );

                this->create_mat_pdv_assembly_map();

                this->create_requested_IWG_list();

                this->create_requested_IQI_list();

                this->create_requested_IQI_type_map();

                this->build_requested_IWG_dof_type_list( aIsResidual );

                this->build_requested_IQI_dof_type_list();

                // set fem set pointer to IWGs FIXME still needed done in constructor?
                for(  std::shared_ptr< IWG > tIWG : mRequestedIWGs )
                {
                    tIWG->set_set_pointer( this );
                }

                // set fem set pointer to IQIs FIXME still needed done in constructor?
                for(  std::shared_ptr< IQI > tIQI : mRequestedIQIs )
                {
                    tIQI->set_set_pointer( this );
                }

                //            if( !aIsForward )
                //            {
                //                this->create_requested_dv_assembly_map();
                //            }
            }
        }

        //------------------------------------------------------------------------------

        void Set::finalize( MSI::Model_Solver_Interface * aModelSolverInterface )
        {
            if ( !mIsEmptySet )    //FIXME this flag is a hack. find better solution
            {
                // delete the field interpolator pointers
                this->delete_pointers();

                // create the field interpolators
                this->create_field_interpolator_managers( aModelSolverInterface );

                // set field interpolator managers for the IWGs
                this->set_IWG_field_interpolator_managers();

                // set field interpolator managers for the IQIs
                this->set_IQI_field_interpolator_managers();
            }
        }

        //------------------------------------------------------------------------------

        void Set::free_memory()
        {
            for(  std::shared_ptr< IWG > tIWG : mIWGs )
            {
                tIWG->free_memory();
            }
        }

        //------------------------------------------------------------------------------

        void Set::delete_pointers()
        {
            if( mMasterFIManager != nullptr )
            {
                delete mMasterFIManager;
                mMasterFIManager = nullptr;
            }
            if( mSlaveFIManager != nullptr )
            {
                delete mSlaveFIManager;
                mSlaveFIManager = nullptr;
            }
            if( mMasterPreviousFIManager != nullptr )
            {
                delete mMasterPreviousFIManager;
                mMasterPreviousFIManager = nullptr;
            }
        }

        //------------------------------------------------------------------------------

        void Set::create_unique_dof_and_dv_type_lists()
        {
            // init dof and dv type counter
            uint tMasterDofCounter = 0;
            uint tMasterDvCounter  = 0;
            uint tSlaveDofCounter = 0;
            uint tSlaveDvCounter  = 0;

            // loop over the IWGs
            for ( std::shared_ptr< IWG > tIWG : mIWGs )
            {
                // get an IWG non unique dof and dv types
                moris::Cell< moris::Cell< MSI::Dof_Type > > tActiveDofType;
                moris::Cell< moris::Cell< PDV_Type > >      tActiveDvType;
                tIWG->get_non_unique_dof_and_dv_types( tActiveDofType, tActiveDvType );

                // update dof and dv type counters
                tMasterDofCounter += tActiveDofType( 0 ).size();
                tMasterDvCounter  += tActiveDvType( 0 ).size();
                tSlaveDofCounter += tActiveDofType( 1 ).size();
                tSlaveDvCounter  += tActiveDvType( 1 ).size();
            }

            // loop over the IQIs
            for ( std::shared_ptr< IQI > tIQI : mIQIs )
            {
                // get an IWG non unique dof and dv types
                moris::Cell< moris::Cell< MSI::Dof_Type > > tActiveDofType;
                moris::Cell< moris::Cell< PDV_Type > >      tActiveDvType;
                tIQI->get_non_unique_dof_and_dv_types( tActiveDofType, tActiveDvType );

                // update dof and dv type counter
                tMasterDofCounter += tActiveDofType( 0 ).size();
                tMasterDvCounter  += tActiveDvType( 0 ).size();
                tSlaveDofCounter += tActiveDofType( 1 ).size();
                tSlaveDvCounter  += tActiveDvType( 1 ).size();
            }

            mUniqueDofTypeListMasterSlave.resize( 2 );
            mUniqueDvTypeListMasterSlave.resize( 2 );

            mUniqueDofTypeListMasterSlave( 0 ).reserve( tMasterDofCounter );
            mUniqueDofTypeListMasterSlave( 1 ).reserve( tSlaveDofCounter );
            mUniqueDvTypeListMasterSlave ( 0 ).reserve( tMasterDvCounter );
            mUniqueDvTypeListMasterSlave ( 1 ).reserve( tSlaveDvCounter );

            // set max size for the unique dof and dv type lists
            mUniqueDofTypeList.reserve( tMasterDofCounter + tSlaveDofCounter );
            mUniqueDvTypeList.reserve( tMasterDvCounter + tSlaveDvCounter );

            // loop over the IWGs
            for ( std::shared_ptr< IWG > tIWG : mIWGs )
            {
                // get non unique dof and dv types
                moris::Cell< moris::Cell< MSI::Dof_Type > > tActiveDofType;
                moris::Cell< moris::Cell< PDV_Type > >      tActiveDvType;
                tIWG->get_non_unique_dof_and_dv_types( tActiveDofType, tActiveDvType );

                // populate the corresponding unique dof and dv type lists
                mUniqueDofTypeListMasterSlave( 0 ).append( tActiveDofType( 0 ) );
                mUniqueDofTypeListMasterSlave( 1 ).append( tActiveDofType( 1 ) );
                mUniqueDvTypeListMasterSlave( 0 ).append( tActiveDvType( 0 ) );
                mUniqueDvTypeListMasterSlave( 1 ).append( tActiveDvType( 1 ) );

                mUniqueDofTypeList.append( tActiveDofType( 0 ) );
                mUniqueDofTypeList.append( tActiveDofType( 1 ) );
                mUniqueDvTypeList.append( tActiveDvType( 0 ) );
                mUniqueDvTypeList.append( tActiveDvType( 1 ) );
            }

            // loop over the IQIs
            for ( std::shared_ptr< IQI > tIQI : mIQIs )
            {
                // get non unique dof and dv types
                moris::Cell< moris::Cell< MSI::Dof_Type > > tActiveDofType;
                moris::Cell< moris::Cell< PDV_Type > >      tActiveDvType;
                tIQI->get_non_unique_dof_and_dv_types( tActiveDofType, tActiveDvType );

                // populate the corresponding unique dof and dv type lists
                mUniqueDofTypeListMasterSlave( 0 ).append( tActiveDofType( 0 ) );
                mUniqueDofTypeListMasterSlave( 1 ).append( tActiveDofType( 1 ) );
                mUniqueDvTypeListMasterSlave( 0 ).append( tActiveDvType( 0 ) );
                mUniqueDvTypeListMasterSlave( 1 ).append( tActiveDvType( 1 ) );

                mUniqueDofTypeList.append( tActiveDofType( 0 ) );
                mUniqueDofTypeList.append( tActiveDofType( 1 ) );
                mUniqueDvTypeList.append( tActiveDvType( 0 ) );
                mUniqueDvTypeList.append( tActiveDvType( 1 ) );
            }

            {
                // make the dof type list unique
                std::sort( ( mUniqueDofTypeListMasterSlave( 0 ).data() ).data(),
                        ( mUniqueDofTypeListMasterSlave( 0 ).data() ).data() + mUniqueDofTypeListMasterSlave( 0 ).size());
                auto last = std::unique( ( mUniqueDofTypeListMasterSlave( 0 ).data() ).data(),
                        ( mUniqueDofTypeListMasterSlave( 0 ).data() ).data() + mUniqueDofTypeListMasterSlave( 0 ).size() );
                auto pos  = std::distance( ( mUniqueDofTypeListMasterSlave( 0 ).data() ).data(), last );
                mUniqueDofTypeListMasterSlave( 0 ).resize( pos );
            }

            {
                // make the dof type list unique
                std::sort( ( mUniqueDofTypeListMasterSlave( 1 ).data() ).data(),
                        ( mUniqueDofTypeListMasterSlave( 1 ).data() ).data() + mUniqueDofTypeListMasterSlave( 1 ).size());
                auto last = std::unique( ( mUniqueDofTypeListMasterSlave( 1 ).data() ).data(),
                        ( mUniqueDofTypeListMasterSlave( 1 ).data() ).data() + mUniqueDofTypeListMasterSlave( 1 ).size() );
                auto pos  = std::distance( ( mUniqueDofTypeListMasterSlave( 1 ).data() ).data(), last );
                mUniqueDofTypeListMasterSlave( 1 ).resize( pos );
            }

            {
                // make the dof type list unique
                std::sort( ( mUniqueDvTypeListMasterSlave( 0 ).data() ).data(),
                        ( mUniqueDvTypeListMasterSlave( 0 ).data() ).data() + mUniqueDvTypeListMasterSlave( 0 ).size());
                auto last = std::unique( ( mUniqueDvTypeListMasterSlave( 0 ).data() ).data(),
                        ( mUniqueDvTypeListMasterSlave( 0 ).data() ).data() + mUniqueDvTypeListMasterSlave( 0 ).size() );
                auto pos  = std::distance( ( mUniqueDvTypeListMasterSlave( 0 ).data() ).data(), last );
                mUniqueDvTypeListMasterSlave( 0 ).resize( pos );
            }

            {
                // make the dof type list unique
                std::sort( ( mUniqueDvTypeListMasterSlave( 1 ).data() ).data(),
                        ( mUniqueDvTypeListMasterSlave( 1 ).data() ).data() + mUniqueDvTypeListMasterSlave( 1 ).size());
                auto last = std::unique( ( mUniqueDvTypeListMasterSlave( 1 ).data() ).data(),
                        ( mUniqueDvTypeListMasterSlave( 1 ).data() ).data() + mUniqueDvTypeListMasterSlave( 1 ).size() );
                auto pos  = std::distance( ( mUniqueDvTypeListMasterSlave( 1 ).data() ).data(), last );
                mUniqueDvTypeListMasterSlave( 1 ).resize( pos );
            }

            {
                // make the dof type list unique
                std::sort( ( mUniqueDofTypeList.data() ).data(),
                        ( mUniqueDofTypeList.data() ).data() + mUniqueDofTypeList.size());
                auto last = std::unique( ( mUniqueDofTypeList.data() ).data(),
                        ( mUniqueDofTypeList.data() ).data() + mUniqueDofTypeList.size() );
                auto pos  = std::distance( ( mUniqueDofTypeList.data() ).data(), last );
                mUniqueDofTypeList.resize( pos );
            }

            {
                // make the dv type list unique
                std::sort( ( mUniqueDvTypeList.data() ).data(),
                        ( mUniqueDvTypeList.data() ).data() + mUniqueDvTypeList.size());
                auto last = std::unique( ( mUniqueDvTypeList.data() ).data(),
                        ( mUniqueDvTypeList.data() ).data() + mUniqueDvTypeList.size() );
                auto pos  = std::distance( ( mUniqueDvTypeList.data() ).data(), last );
                mUniqueDvTypeList.resize( pos );
            }
        }

        //------------------------------------------------------------------------------

        void Set::create_dof_and_dv_type_lists()
        {
            // get number of dof and dv types
            uint tNumDofTypes = this->get_num_unique_dof_types();
            uint tNumDvTypes  = this->get_num_unique_dv_types();

            // set size for the global dof type list
            mMasterDofTypes.reserve( tNumDofTypes );
            mSlaveDofTypes .reserve( tNumDofTypes );
            mMasterDvTypes.reserve( tNumDvTypes );
            mSlaveDvTypes .reserve( tNumDvTypes );

            // create a list to check if dof type is already in the list
            Matrix< DDSMat > tMasterCheckList( tNumDofTypes, 1, -1 );
            Matrix< DDSMat > tSlaveCheckList ( tNumDofTypes, 1, -1 );
            Matrix< DDSMat > tMasterDvCheckList( tNumDvTypes, 1, -1 );
            Matrix< DDSMat > tSlaveDvCheckList ( tNumDvTypes, 1, -1 );

            // loop over the IWGs
            for ( std::shared_ptr< IWG > tIWG : mIWGs )
            {
                // get master dof and dv types for the IWG
                moris::Cell< moris::Cell< MSI::Dof_Type > > tDofTypeMaster =
                        tIWG->get_global_dof_type_list();
                moris::Cell< moris::Cell< PDV_Type > >  tDvTypeMaster =
                        tIWG->get_global_dv_type_list();

                // loop over the IWG active master dof type
                for ( uint iDOF = 0; iDOF < tDofTypeMaster.size(); iDOF++ )
                {
                    // get set index for the treated master dof type
                    sint tDofTypeindex = this->get_index_from_unique_dof_type_map( tDofTypeMaster( iDOF )( 0 ) );

                    // if dof enum not in the list
                    if ( tMasterCheckList( tDofTypeindex ) != 1 )
                    {
                        // put the dof type in the checklist
                        tMasterCheckList( tDofTypeindex ) = 1;

                        // put the dof type in the global type list
                        mMasterDofTypes.push_back( tDofTypeMaster( iDOF ) );
                    }
                }

                // loop over the IWG active master dv type
                for ( uint iDv = 0; iDv < tDvTypeMaster.size(); iDv++ )
                {
                    // get set index for the treated master dof type
                    sint tDvTypeindex = this->get_index_from_unique_dv_type_map( tDvTypeMaster( iDv )( 0 ) );

                    // if dv enum not in the list
                    if ( tMasterDvCheckList( tDvTypeindex ) != 1 )
                    {
                        // put the dof type in the checklist
                        tMasterDvCheckList( tDvTypeindex ) = 1;

                        // put the dof type in the global type list
                        mMasterDvTypes.push_back( tDvTypeMaster( iDv ) );
                    }
                }

                // get slave dof and dv types for the IWG
                moris::Cell< moris::Cell< MSI::Dof_Type > > tDofTypeSlave =
                        tIWG->get_global_dof_type_list( mtk::Master_Slave::SLAVE );

                moris::Cell< moris::Cell< PDV_Type > >  tDvTypeSlave =
                        tIWG->get_global_dv_type_list( mtk::Master_Slave::SLAVE );

                // loop over the IWG active slave dof type
                for ( uint iDOF = 0; iDOF < tDofTypeSlave.size(); iDOF++ )
                {
                    // get set index for the treated slave dof type
                    sint tDofTypeindex = this->get_index_from_unique_dof_type_map( tDofTypeSlave( iDOF )( 0 ) );

                    // if dof enum not in the list
                    if ( tSlaveCheckList( tDofTypeindex ) != 1 )
                    {
                        // put the dof type in the checklist
                        tSlaveCheckList( tDofTypeindex ) = 1;

                        // put the dof type in the global type list
                        mSlaveDofTypes.push_back( tDofTypeSlave( iDOF ) );
                    }
                }

                // loop over the IWG active slave dv type
                for ( uint iDv = 0; iDv < tDvTypeSlave.size(); iDv++ )
                {
                    // get set index for the treated slave dv type
                    sint tDvTypeindex = this->get_index_from_unique_dv_type_map( tDvTypeSlave( iDv )( 0 ) );

                    // if dv enum not in the list
                    if ( tSlaveDvCheckList( tDvTypeindex ) != 1 )
                    {
                        // put the dv type in the checklist
                        tSlaveDvCheckList( tDvTypeindex ) = 1;

                        // put the dv type in the global type list
                        mSlaveDvTypes.push_back( tDvTypeSlave( iDv ) );
                    }
                }
            }

            // loop over the IQIs
            for ( std::shared_ptr< IQI > tIQI : mIQIs )
            {
                // get master dof and dv types for the IWG
                moris::Cell< moris::Cell< MSI::Dof_Type > > tDofTypeMaster =
                        tIQI->get_global_dof_type_list();
                moris::Cell< moris::Cell< PDV_Type > >  tDvTypeMaster =
                        tIQI->get_global_dv_type_list();

                // loop over the IQI active master dof type
                for ( uint iDOF = 0; iDOF < tDofTypeMaster.size(); iDOF++ )
                {
                    // get set index for the treated master dof type
                    sint tDofTypeindex = this->get_index_from_unique_dof_type_map( tDofTypeMaster( iDOF )( 0 ) );

                    // if dof enum not in the list
                    if ( tMasterCheckList( tDofTypeindex ) != 1 )
                    {
                        // put the dof type in the checklist
                        tMasterCheckList( tDofTypeindex ) = 1;

                        // put the dof type in the global type list
                        mMasterDofTypes.push_back( tDofTypeMaster( iDOF ) );
                    }
                }

                // loop over the IQI active master dv type
                for ( uint iDv = 0; iDv < tDvTypeMaster.size(); iDv++ )
                {
                    // get set index for the treated master dv type
                    sint tDvTypeindex = this->get_index_from_unique_dv_type_map( tDvTypeMaster( iDv )( 0 ) );

                    // if dv enum not in the list
                    if ( tMasterDvCheckList( tDvTypeindex ) != 1 )
                    {
                        // put the dv type in the checklist
                        tMasterDvCheckList( tDvTypeindex ) = 1;

                        // put the dv type in the global type list
                        mMasterDvTypes.push_back( tDvTypeMaster( iDv ) );
                    }
                }

                // get slave dof and dv types for the IWG
                moris::Cell< moris::Cell< MSI::Dof_Type > > tDofTypeSlave =
                        tIQI->get_global_dof_type_list( mtk::Master_Slave::SLAVE );
                moris::Cell< moris::Cell< PDV_Type > >  tDvTypeSlave =
                        tIQI->get_global_dv_type_list( mtk::Master_Slave::SLAVE );

                // loop over the IWG active slave dof type
                for ( uint iDOF = 0; iDOF < tDofTypeSlave.size(); iDOF++ )
                {
                    // get set index for the treated slave dof type
                    sint tDofTypeindex = this->get_index_from_unique_dof_type_map( tDofTypeSlave( iDOF )( 0 ) );

                    // if dof enum not in the list
                    if ( tSlaveCheckList( tDofTypeindex ) != 1 )
                    {
                        // put the dof type in the checklist
                        tSlaveCheckList( tDofTypeindex ) = 1;

                        // put the dof type in the global type list
                        mSlaveDofTypes.push_back( tDofTypeSlave( iDOF ) );
                    }
                }

                // loop over the IWG active slave dv type
                for ( uint iDv = 0; iDv < tDvTypeSlave.size(); iDv++ )
                {
                    // get set index for the treated slave dv type
                    sint tDvTypeindex = this->get_index_from_unique_dv_type_map( tDvTypeSlave( iDv )( 0 ) );

                    // if dv enum not in the list
                    if ( tSlaveDvCheckList( tDvTypeindex ) != 1 )
                    {
                        // put the dv type in the checklist
                        tSlaveDvCheckList( tDvTypeindex ) = 1;

                        // put the dv type in the global type list
                        mSlaveDvTypes.push_back( tDvTypeSlave( iDv ) );
                    }
                }
            }

            // shrink list to fit to number of unique dof and dv types
            mMasterDofTypes.shrink_to_fit();
            mSlaveDofTypes .shrink_to_fit();
            mMasterDvTypes.shrink_to_fit();
            mSlaveDvTypes .shrink_to_fit();
        }

        //------------------------------------------------------------------------------
        void Set::create_unique_dof_and_dv_type_maps()
        {
            // dof types
            //------------------------------------------------------------------------------
            // Create temporary dof type list
            const moris::Cell< enum MSI::Dof_Type > & tDofType = get_unique_dof_type_list();

            //Get number of unique adofs of this equation object
            moris::uint tNumUniqueDofTypes = tDofType.size();

            // Get maximal dof type enum number
            moris::sint tMaxDofTypeEnumNumber = 0;

            // Loop over all pdof types to get the highest enum index
            for ( moris::uint Ii = 0; Ii < tNumUniqueDofTypes; Ii++ )
            {
                tMaxDofTypeEnumNumber = std::max( tMaxDofTypeEnumNumber, static_cast< int >( tDofType( Ii ) ) );
            }

            // +1 because c++ is 0 based
            tMaxDofTypeEnumNumber++;

            // Set size of mapping matrix
            mUniqueDofTypeMap.set_size( tMaxDofTypeEnumNumber, 1, -1 );

            // Loop over all pdof types to create the mapping matrix
            for ( moris::uint Ii = 0; Ii < tNumUniqueDofTypes; Ii++ )
            {
                mUniqueDofTypeMap( static_cast< int >( tDofType( Ii ) ), 0 ) = Ii;
            }

            // dv types
            //------------------------------------------------------------------------------
            // Create temporary dv type list
            moris::Cell< enum PDV_Type > tDvType = get_unique_dv_type_list();

            //Get number of unique dvs of this equation object
            moris::uint tNumUniqueDvTypes = tDvType.size();

            // Get maximal dv type enum number
            moris::sint tMaxDvTypeEnumNumber = 0;

            // Loop over all dv types to get the highest enum index
            for ( moris::uint Ii = 0; Ii < tNumUniqueDvTypes; Ii++ )
            {
                tMaxDvTypeEnumNumber =
                        std::max( tMaxDvTypeEnumNumber, static_cast< int >( tDvType( Ii ) ) );
            }

            // +1 because c++ is 0 based
            tMaxDvTypeEnumNumber++;

            // Set size of mapping matrix
            mUniqueDvTypeMap.set_size( tMaxDvTypeEnumNumber, 1, -1 );

            // Loop over all dv types to create the mapping matrix
            for ( moris::uint Ii = 0; Ii < tNumUniqueDvTypes; Ii++ )
            {
                mUniqueDvTypeMap( static_cast< int >( tDvType( Ii ) ), 0 ) = Ii;
            }
        }

        //------------------------------------------------------------------------------
        void Set::create_dof_and_dv_type_maps()
        {
            // get number of master dof types
            uint tMasterNumDofs = this->get_dof_type_list().size();

            // get maximal dof type enum
            sint tMaxEnum = -1;

            // loop over the IWGs
            for ( uint iDOF = 0; iDOF < tMasterNumDofs; iDOF++ )
            {
                for ( uint Ik = 0; Ik < mMasterDofTypes( iDOF ).size(); Ik++ )
                {
                    // get the highest dof type enum
                    tMaxEnum = std::max( tMaxEnum, static_cast< int >( mMasterDofTypes( iDOF )( Ik ) ) );
                }
            }

            // get number of slave dof types
            uint tSlaveNumDofs = this->get_dof_type_list( mtk::Master_Slave::SLAVE ).size();

            // loop over the IWGs
            for ( uint iDOF = 0; iDOF < tSlaveNumDofs; iDOF++ )
            {
                for ( uint Ik = 0; Ik < mSlaveDofTypes( iDOF ).size(); Ik++ )
                {
                    // get the highest dof type enum
                    tMaxEnum = std::max( tMaxEnum, static_cast< int >( mSlaveDofTypes( iDOF )( Ik ) ) );
                }
            }
            // +1 since start at 0
            tMaxEnum++;

            MORIS_ASSERT( tMaxEnum != -1, "Set::create_dof_and_dv_type_maps(), no information to build dof type map" );

            tMaxEnum = static_cast< int >( MSI::Dof_Type::END_ENUM );

            // set size of dof type map    // FIXME replace with map
            mMasterDofTypeMap.set_size( tMaxEnum, 1, -1 );

            // loop over dof types
            for ( uint iDOF = 0; iDOF < tMasterNumDofs; iDOF++ )
            {
                mMasterDofTypeMap( static_cast< int >( mMasterDofTypes( iDOF )( 0 ) ), 0 ) = iDOF;
            }

            // set size of dof type map
            mSlaveDofTypeMap.set_size( tMaxEnum, 1, -1 );

            // loop over dof types
            for ( uint iDOF = 0; iDOF < tSlaveNumDofs; iDOF++ )
            {
                mSlaveDofTypeMap( static_cast< int >( mSlaveDofTypes( iDOF )( 0 ) ), 0 ) = iDOF;
            }

            // get number of master dv types
            uint tMasterNumDvs  = this->get_dv_type_list().size();

            // get maximal dv type enum
            tMaxEnum = -1;

            // loop over the dv types
            for ( uint iDv = 0; iDv < tMasterNumDvs; iDv++ )
            {
                for ( uint Ik = 0; Ik < mMasterDvTypes( iDv ).size(); Ik++ )
                {
                    // get the highest dof type enum
                    tMaxEnum = std::max( tMaxEnum, static_cast< int >( mMasterDvTypes( iDv )( Ik ) ) );
                }
            }

            // get number of slave dv types
            uint tSlaveNumDvs =  this->get_dv_type_list( mtk::Master_Slave::SLAVE ).size();

            // loop over the IWGs
            for ( uint iDOF = 0; iDOF < tSlaveNumDvs; iDOF++ )
            {
                for ( uint Ik = 0; Ik < mSlaveDvTypes( iDOF ).size(); Ik++ )
                {
                    // get the highest dof type enum
                    tMaxEnum = std::max( tMaxEnum, static_cast< int >( mSlaveDvTypes( iDOF )( Ik ) ) );
                }
            }
            // +1 since start at 0
            tMaxEnum++;

            MORIS_ASSERT( tMaxEnum != -1, "Set::create_dv_type_map(), no information to build dv type map" );

            // set size of dv type map    // FIXME replace with map
            mMasterDvTypeMap.set_size( tMaxEnum, 1, -1 );

            // loop over dv types
            for ( uint iDv = 0; iDv < tMasterNumDvs; iDv++ )
            {
                mMasterDvTypeMap( static_cast< int >( mMasterDvTypes( iDv )( 0 ) ), 0 ) = iDv;
            }

            // set size of dv type map
            mSlaveDvTypeMap.set_size( tMaxEnum, 1, -1 );

            // loop over dv types
            for ( uint iDv = 0; iDv < tSlaveNumDvs; iDv++ )
            {
                mSlaveDvTypeMap( static_cast< int >( mSlaveDvTypes( iDv )( 0 ) ), 0 ) = iDv;
            }
        }

        //-----------------------------------------------------------------------------

        void Set::create_field_interpolator_managers(
                MSI::Model_Solver_Interface * aModelSolverInterface )
        {
            // create the master field interpolator manager
            mMasterFIManager = new Field_Interpolator_Manager(
                    mMasterDofTypes,
                    mMasterDvTypes,
                    this );

            // create the geometry interpolators on the master FI manager
            mMasterFIManager->create_geometry_interpolators();

            // create the field interpolators on the master FI manager
            mMasterFIManager->create_field_interpolators( aModelSolverInterface );

            // create the slave field interpolator manager
            mSlaveFIManager = new Field_Interpolator_Manager(
                    mSlaveDofTypes,
                    mSlaveDvTypes,
                    this,
                    mtk::Master_Slave::SLAVE );

            // create the geometry interpolators on the slave FI manager
            mSlaveFIManager->create_geometry_interpolators();

            // create the field interpolators on the slave FI manager
            mSlaveFIManager->create_field_interpolators( aModelSolverInterface );

            // if time sideset
            if( mElementType == fem::Element_Type::TIME_SIDESET )
            {
                // create the master field interpolator manager
                mMasterPreviousFIManager = new Field_Interpolator_Manager(
                        mMasterDofTypes,
                        mMasterDvTypes,
                        this );

                // create the geometry interpolators on the master FI manager
                mMasterPreviousFIManager->create_geometry_interpolators();

                // create the field interpolators on the master FI manager
                mMasterPreviousFIManager->create_field_interpolators( aModelSolverInterface );
            }
        }

        //------------------------------------------------------------------------------

        Field_Interpolator_Manager * Set::get_field_interpolator_manager(
                mtk::Master_Slave aIsMaster )
        {
            switch ( aIsMaster )
            {
                case mtk::Master_Slave::MASTER :
                    return mMasterFIManager;

                case mtk::Master_Slave::SLAVE :
                    return mSlaveFIManager;

                default :
                    MORIS_ERROR( false, "Set::get_field_interpolator_manager - can only be master or slave.");
                    return mMasterFIManager;
            }
        }

        //------------------------------------------------------------------------------

        Field_Interpolator_Manager * Set::get_field_interpolator_manager_previous_time(
                mtk::Master_Slave aIsMaster )
        {
            switch ( aIsMaster )
            {
                case mtk::Master_Slave::MASTER :
                    return mMasterPreviousFIManager;

                default :
                    MORIS_ERROR( false, "Set::get_field_interpolator_manager - can only be master.");
                    return mMasterPreviousFIManager;
            }
        }

        //------------------------------------------------------------------------------

        void Set::set_IWG_field_interpolator_managers()
        {
            // loop over the IWGs
            for ( std::shared_ptr< IWG > tIWG : mIWGs )
            {
                // set the master FI manager
                tIWG->set_field_interpolator_manager( mMasterFIManager );

                // if double sideset, set slave
                if( mElementType == fem::Element_Type::DOUBLE_SIDESET )
                {
                    // set IWG slave field interpolator manager
                    tIWG->set_field_interpolator_manager(
                            mSlaveFIManager,
                            mtk::Master_Slave::SLAVE );
                }

                // if time sideset, set previous
                if( mElementType == fem::Element_Type::TIME_SIDESET )
                {
                    // set IWG master field interpolator manager for previous time step
                    tIWG->set_field_interpolator_manager_previous_time(
                            mMasterPreviousFIManager,
                            mtk::Master_Slave::MASTER );
                }
            }
        }

        //------------------------------------------------------------------------------

        void Set::set_IWG_cluster_for_stabilization_parameters( fem::Cluster * aCluster )
        {
            // loop over the IWGs
            for( auto tIWG : mIWGs)
            {
                // get the SP from the IWG
                moris::Cell< std::shared_ptr< Stabilization_Parameter > > & tSPs =
                        tIWG->get_stabilization_parameters();

                // loop over the SP
                for( auto tSP : tSPs )
                {
                    // check if SP is null
                    if( tSP != nullptr )
                    {
                        // set the fem cluster
                        tSP->set_cluster( aCluster );
                    }
                }
            }
        }

        //------------------------------------------------------------------------------

        void Set::set_IQI_cluster_for_stabilization_parameters( fem::Cluster * aCluster )
        {
            // loop over the IQIs
            for( auto tIQI : mIQIs)
            {
                // get the SP from the IQI
                moris::Cell< std::shared_ptr< Stabilization_Parameter > > & tSPs =
                        tIQI->get_stabilization_parameters();

                // loop over the SPs
                for( auto tSP : tSPs )
                {
                    // check if SP is null
                    if( tSP != nullptr )
                    {
                        // set the fem cluster
                        tSP->set_cluster( aCluster );
                    }
                }
            }
        }

        //------------------------------------------------------------------------------

        void Set::set_IQI_field_interpolator_managers()
        {
            // loop over the IQIs
            for ( std::shared_ptr< IQI > tIQI : mIQIs )
            {
                // set IQI master FI manager
                tIQI->set_field_interpolator_manager( mMasterFIManager );

                // if double sideset, set slave
                if( mElementType == fem::Element_Type::DOUBLE_SIDESET )
                {
                    // set IQI slave FI manager
                    tIQI->set_field_interpolator_manager(
                            mSlaveFIManager,
                            mtk::Master_Slave::SLAVE );
                }
            }
        }

        //------------------------------------------------------------------------------

        void Set::create_residual_dof_assembly_map()
        {
            // get the list of requested dof types by the solver
            moris::Cell < enum MSI::Dof_Type > tRequestedDofTypes =
                    this->get_model_solver_interface()->
                    get_solver_interface()->
                    get_requested_dof_types();

            // init the max index for dof types
            sint tMaxDofIndex = -1;

            // loop over the requested dof types
            for( uint Ik = 0; Ik < tRequestedDofTypes.size(); Ik++ )
            {
                // get the set index for the requested master dof type
                sint tDofIndex = this->get_dof_index_for_type(
                        tRequestedDofTypes( Ik ),
                        mtk::Master_Slave::MASTER );

                // if the index was set (and is different from -1)
                if( tDofIndex != -1 )
                {
                    // update the max index for dof type
                    tMaxDofIndex = std::max( tMaxDofIndex, tDofIndex );
                }

                // get the set index for the requested slave dof type
                tDofIndex = this->get_dof_index_for_type(
                        tRequestedDofTypes( Ik ),
                        mtk::Master_Slave::SLAVE );

                // if the index was set (and is different -1)
                if( tDofIndex != -1 )
                {
                    // update the max index for dof type
                    tMaxDofIndex = std::max( tMaxDofIndex, tDofIndex );
                }
            }
            // add +1 to the max index for dof type (since 0 based)
            tMaxDofIndex++;

            // set size for the residual assembly map
            mResDofAssemblyMap.resize( tMaxDofIndex );

            // init the residual assembly map
            for( uint Ik = 0; Ik < mResDofAssemblyMap.size(); Ik++ )
            {
                mResDofAssemblyMap( Ik ).set_size( 1, 2, -1 );
            }

            // init dof coefficients counter
            uint tCounter = 0;

            // loop over the requested dof types
            for( uint Ik = 0; Ik < tRequestedDofTypes.size(); Ik++ )
            {
                // get the set index for the requested master dof type
                sint tDofIndex = this->get_dof_index_for_type(
                        tRequestedDofTypes( Ik ),
                        mtk::Master_Slave::MASTER );

                // if the index was set (and is different from -1)
                if( tDofIndex != -1 )
                {
                    // get the number of coefficients re;ated to the master dof type
                    uint tNumCoeff = mMasterFIManager->
                            get_field_interpolators_for_type( tRequestedDofTypes( Ik ) )->
                            get_number_of_space_time_coefficients();

                    // fill the residual assembly map with starting and ending indices for the master dof type
                    mResDofAssemblyMap( tDofIndex )( 0, 0 ) = tCounter;
                    mResDofAssemblyMap( tDofIndex )( 0, 1 ) = tCounter + tNumCoeff - 1;

                    // update the dof coefficient counter
                    tCounter += tNumCoeff;
                }
            }

            // loop over the requested dof types
            for( uint Ik = 0; Ik < tRequestedDofTypes.size(); Ik++ )
            {
                //get the set index for the requested slave dof type
                sint tDofIndex = this->get_dof_index_for_type( tRequestedDofTypes( Ik ),
                        mtk::Master_Slave::SLAVE );

                // if the dof type was set (its set index is different from -1)
                if( tDofIndex != -1 )
                {
                    // get the number of coefficients for the slave dof type
                    uint tNumCoeff = mSlaveFIManager->
                            get_field_interpolators_for_type( tRequestedDofTypes( Ik ) )->
                            get_number_of_space_time_coefficients();

                    // fill the residual assembly map with starting and ending indices for the slave dof type
                    mResDofAssemblyMap( tDofIndex )( 0, 0 ) = tCounter;
                    mResDofAssemblyMap( tDofIndex )( 0, 1 ) = tCounter + tNumCoeff - 1;

                    // update the dof coefficient counter
                    tCounter += tNumCoeff;
                }
            }
        }

        //------------------------------------------------------------------------------

        void Set::create_dof_assembly_map( const bool aIsResidual )
        {
            if( aIsResidual )
            {
                this->create_staggered_jacobian_dof_assembly_map();
            }
            else
            {
                this->create_jacobian_dof_assembly_map();
            }
        }

        //------------------------------------------------------------------------------

        void Set::create_jacobian_dof_assembly_map()
        {
            // get list of requested dof types (by the solver)
            moris::Cell < enum MSI::Dof_Type >  tRequestedDofTypes =
                    this->get_model_solver_interface()->
                    get_solver_interface()->
                    get_requested_dof_types();

            sint tMaxDofIndex = -1;

            // master
            for( uint Ik = 0; Ik < tRequestedDofTypes.size(); Ik++ )
            {
                sint tDofIndex = this->get_dof_index_for_type(
                        tRequestedDofTypes( Ik ),
                        mtk::Master_Slave::MASTER );

                if( tDofIndex != -1 )
                {
                    tMaxDofIndex = std::max( tMaxDofIndex, tDofIndex );
                }

                tDofIndex = this->get_dof_index_for_type(
                        tRequestedDofTypes( Ik ),
                        mtk::Master_Slave::SLAVE );

                if( tDofIndex != -1 )
                {
                    tMaxDofIndex = std::max( tMaxDofIndex, tDofIndex );
                }
            }

            tMaxDofIndex++;

            mJacDofAssemblyMap.resize( tMaxDofIndex );

            for( uint Ik = 0; Ik < mResDofAssemblyMap.size(); Ik++ )
            {
                mJacDofAssemblyMap( Ik ).set_size( tMaxDofIndex, 2, -1 );
            }

            for( uint Ik = 0; Ik < tRequestedDofTypes.size(); Ik++ )
            {
                sint tDofIndex = this->get_dof_index_for_type(
                        tRequestedDofTypes( Ik ),
                        mtk::Master_Slave::MASTER );

                if( tDofIndex != -1 )
                {
                    uint tCounter_2 = 0;

                    for( uint Ii = 0; Ii < tRequestedDofTypes.size(); Ii++ )
                    {
                        sint tDofIndex_2 = this->get_dof_index_for_type(
                                tRequestedDofTypes( Ii ),
                                mtk::Master_Slave::MASTER );

                        if( tDofIndex_2 != -1 )
                        {
                            uint tNumCoeff_2 = mMasterFIManager->
                                    get_field_interpolators_for_type( tRequestedDofTypes( Ii ) )->
                                    get_number_of_space_time_coefficients();

                            mJacDofAssemblyMap( tDofIndex )( tDofIndex_2, 0 ) = tCounter_2;
                            mJacDofAssemblyMap( tDofIndex )( tDofIndex_2, 1 ) = tCounter_2 + tNumCoeff_2 - 1;

                            tCounter_2 += tNumCoeff_2;
                        }
                    }

                    // slave
                    for( uint Ii = 0; Ii < tRequestedDofTypes.size(); Ii++ )
                    {
                        sint tDofIndex_2 = this->get_dof_index_for_type(
                                tRequestedDofTypes( Ii ),
                                mtk::Master_Slave::SLAVE );

                        if( tDofIndex_2 != -1 )
                        {
                            uint tNumCoeff_2 = mSlaveFIManager->
                                    get_field_interpolators_for_type( tRequestedDofTypes( Ii ) )->
                                    get_number_of_space_time_coefficients();

                            mJacDofAssemblyMap( tDofIndex )( tDofIndex_2, 0 ) = tCounter_2;
                            mJacDofAssemblyMap( tDofIndex )( tDofIndex_2, 1 ) = tCounter_2 + tNumCoeff_2 - 1;

                            tCounter_2 += tNumCoeff_2;
                        }
                    }
                }
            }

            for( uint Ik = 0; Ik < tRequestedDofTypes.size(); Ik++ )
            {
                sint tDofIndex = this->get_dof_index_for_type(
                        tRequestedDofTypes( Ik ),
                        mtk::Master_Slave::SLAVE  );

                if( tDofIndex != -1 )
                {
                    uint tCounter_2 = 0;

                    for( uint Ii = 0; Ii < tRequestedDofTypes.size(); Ii++ )
                    {
                        sint tDofIndex_2 = this->get_dof_index_for_type(
                                tRequestedDofTypes( Ii ),
                                mtk::Master_Slave::MASTER  );

                        if( tDofIndex_2 != -1 )
                        {
                            uint tNumCoeff_2 = mMasterFIManager->
                                    get_field_interpolators_for_type( tRequestedDofTypes( Ii ) )->
                                    get_number_of_space_time_coefficients();

                            mJacDofAssemblyMap( tDofIndex )( tDofIndex_2, 0 ) = tCounter_2;
                            mJacDofAssemblyMap( tDofIndex )( tDofIndex_2, 1 ) = tCounter_2 + tNumCoeff_2 - 1;

                            tCounter_2 += tNumCoeff_2;
                        }
                    }

                    for( uint Ii = 0; Ii < tRequestedDofTypes.size(); Ii++ )
                    {
                        sint tDofIndex_2 = this->get_dof_index_for_type(
                                tRequestedDofTypes( Ii ),
                                mtk::Master_Slave::SLAVE  );

                        if( tDofIndex_2 != -1 )
                        {
                            uint tNumCoeff_2 = mSlaveFIManager->
                                    get_field_interpolators_for_type( tRequestedDofTypes( Ii ) )->
                                    get_number_of_space_time_coefficients();

                            mJacDofAssemblyMap( tDofIndex )( tDofIndex_2, 0 ) = tCounter_2;
                            mJacDofAssemblyMap( tDofIndex )( tDofIndex_2, 1 ) = tCounter_2 + tNumCoeff_2 - 1;

                            tCounter_2 += tNumCoeff_2;
                        }
                    }
                }
            }
        }

        //------------------------------------------------------------------------------

        void Set::create_staggered_jacobian_dof_assembly_map()
        {
            // get list of requested dof types
            moris::Cell < enum MSI::Dof_Type > tRequestedDofTypes =
                    this->get_model_solver_interface()->
                    get_solver_interface()->
                    get_requested_dof_types();

            moris::Cell< moris::Cell < enum MSI::Dof_Type > > tSecondaryDofTypes =
                    this->get_model_solver_interface()->
                    get_solver_interface()->
                    get_secondary_dof_types();

            sint tMaxDofIndex = -1;

            for( uint Ik = 0; Ik < tRequestedDofTypes.size(); Ik++ )
            {
                sint tDofIndex = this->get_dof_index_for_type(
                        tRequestedDofTypes( Ik ),
                        mtk::Master_Slave::MASTER );

                if( tDofIndex != -1 )
                {
                    tMaxDofIndex = std::max( tMaxDofIndex, tDofIndex );
                }

                tDofIndex = this->get_dof_index_for_type(
                        tRequestedDofTypes( Ik ),
                        mtk::Master_Slave::SLAVE );

                if( tDofIndex != -1 )
                {
                    tMaxDofIndex = std::max( tMaxDofIndex, tDofIndex );
                }
            }
            tMaxDofIndex++;

            sint tMaxDofIndexSec = -1;

            for( uint Ik = 0; Ik < tSecondaryDofTypes.size(); Ik++ )
            {
                sint tDofIndex = this->get_dof_index_for_type(
                        tSecondaryDofTypes( Ik )( 0 ),
                        mtk::Master_Slave::MASTER );

                if( tDofIndex != -1 )
                {
                    tMaxDofIndexSec = std::max( tMaxDofIndexSec, tDofIndex );
                }

                tDofIndex = this->get_dof_index_for_type(
                        tSecondaryDofTypes( Ik )( 0 ),
                        mtk::Master_Slave::SLAVE );

                if( tDofIndex != -1 )
                {
                    tMaxDofIndexSec = std::max( tMaxDofIndexSec, tDofIndex );
                }
            }
            tMaxDofIndexSec++;

            mJacDofAssemblyMap.resize( tMaxDofIndex );

            for( uint Ik = 0; Ik < mResDofAssemblyMap.size(); Ik++ )
            {
                mJacDofAssemblyMap( Ik ).set_size( tMaxDofIndexSec, 2, -1 );
            }

            for( uint Ik = 0; Ik < tRequestedDofTypes.size(); Ik++ )
            {
                sint tDofIndex = this->get_dof_index_for_type(
                        tRequestedDofTypes( Ik ),
                        mtk::Master_Slave::MASTER );

                if( tDofIndex != -1 )
                {
                    uint tCounter_2 = 0;

                    for( uint Ii = 0; Ii < tSecondaryDofTypes.size(); Ii++ )
                    {
                        sint tDofIndex_2 = this->get_dof_index_for_type(
                                tSecondaryDofTypes( Ii )( 0 ),
                                mtk::Master_Slave::MASTER );

                        if( tDofIndex_2 != -1 )
                        {
                            uint tNumCoeff_2 = mMasterFIManager->
                                    get_field_interpolators_for_type( tSecondaryDofTypes( Ii )( 0 ) )->
                                    get_number_of_space_time_coefficients();

                            mJacDofAssemblyMap( tDofIndex )( tDofIndex_2, 0 ) = tCounter_2;
                            mJacDofAssemblyMap( tDofIndex )( tDofIndex_2, 1 ) = tCounter_2 + tNumCoeff_2 - 1;

                            tCounter_2 += tNumCoeff_2;
                        }
                    }

                    for( uint Ii = 0; Ii < tSecondaryDofTypes.size(); Ii++ )
                    {
                        sint tDofIndex_2 = this->get_dof_index_for_type(
                                tSecondaryDofTypes( Ii )( 0 ),
                                mtk::Master_Slave::SLAVE );

                        if( tDofIndex_2 != -1 )
                        {
                            uint tNumCoeff_2 = mSlaveFIManager->
                                    get_field_interpolators_for_type(  tSecondaryDofTypes( Ii )( 0 ) )->
                                    get_number_of_space_time_coefficients();

                            mJacDofAssemblyMap( tDofIndex )( tDofIndex_2, 0 ) = tCounter_2;
                            mJacDofAssemblyMap( tDofIndex )( tDofIndex_2, 1 ) = tCounter_2 + tNumCoeff_2 - 1;

                            tCounter_2 += tNumCoeff_2;
                        }
                    }
                }
            }

            for( uint Ik = 0; Ik < tRequestedDofTypes.size(); Ik++ )
            {
                sint tDofIndex = this->get_dof_index_for_type(
                        tRequestedDofTypes( Ik ),
                        mtk::Master_Slave::SLAVE  );

                if( tDofIndex != -1 )
                {
                    uint tCounter_2 = 0;

                    for( uint Ii = 0; Ii < tSecondaryDofTypes.size(); Ii++ )
                    {
                        sint tDofIndex_2 = this->get_dof_index_for_type(
                                tSecondaryDofTypes( Ii )( 0 ),
                                mtk::Master_Slave::MASTER  );

                        if( tDofIndex_2 != -1 )
                        {
                            uint tNumCoeff_2 = mMasterFIManager->
                                    get_field_interpolators_for_type( tSecondaryDofTypes( Ii )( 0 ) )
                                    ->get_number_of_space_time_coefficients();

                            mJacDofAssemblyMap( tDofIndex )( tDofIndex_2, 0 ) = tCounter_2;
                            mJacDofAssemblyMap( tDofIndex )( tDofIndex_2, 1 ) = tCounter_2 + tNumCoeff_2 - 1;

                            tCounter_2 += tNumCoeff_2;
                        }
                    }

                    for( uint Ii = 0; Ii < tSecondaryDofTypes.size(); Ii++ )
                    {
                        sint tDofIndex_2 = this->get_dof_index_for_type(
                                tSecondaryDofTypes( Ii )( 0 ),
                                mtk::Master_Slave::SLAVE  );

                        if( tDofIndex_2 != -1 )
                        {
                            uint tNumCoeff_2 = mSlaveFIManager->
                                    get_field_interpolators_for_type(  tSecondaryDofTypes( Ii )( 0 ) )->
                                    get_number_of_space_time_coefficients();

                            mJacDofAssemblyMap( tDofIndex )( tDofIndex_2, 0 ) = tCounter_2;
                            mJacDofAssemblyMap( tDofIndex )( tDofIndex_2, 1 ) = tCounter_2 + tNumCoeff_2 - 1;

                            tCounter_2 += tNumCoeff_2;
                        }
                    }
                }
            }
        }

        //------------------------------------------------------------------------------

        void Set::create_mat_pdv_assembly_map()
        {
            // get the list of requested dv types by the opt solver
            moris::Cell< moris::Cell< enum PDV_Type > > tRequestedDvTypes;
            this->get_ip_dv_types_for_set( tRequestedDvTypes );

            // init the max index for dv types
            sint tMaxDvIndex = -1;

            // loop over the dv types
            for( uint Ik = 0; Ik < tRequestedDvTypes.size(); Ik++ )
            {
                // get the set index for the requested master dof type
                sint tDvIndex = this->get_dv_index_for_type(
                        tRequestedDvTypes( Ik )( 0 ),
                        mtk::Master_Slave::MASTER );

                // if the index was set (and is different from -1)
                if( tDvIndex != -1 )
                {
                    // update the max index for dv type
                    tMaxDvIndex = std::max( tMaxDvIndex, tDvIndex );
                }

                // get the set index for the requested slave slave type
                tDvIndex = this->get_dv_index_for_type(
                        tRequestedDvTypes( Ik )( 0 ),
                        mtk::Master_Slave::SLAVE );

                // if the index was set (and is different -1)
                if( tDvIndex != -1 )
                {
                    // update the max index for dv type
                    tMaxDvIndex = std::max( tMaxDvIndex, tDvIndex );
                }
            }
            // add +1 to the max index for dv type (since 0 based)
            tMaxDvIndex++;

            // set size for the dv assembly map
            mPdvMatAssemblyMap.resize( tMaxDvIndex );

            // init the dv assembly map
            for( uint Ik = 0; Ik < mPdvMatAssemblyMap.size(); Ik++ )
            {
                mPdvMatAssemblyMap( Ik ).set_size( 1, 2, -1 );
            }

            // init dv coefficients counter
            uint tCounter = 0;

            // loop over the dv types
            for( uint Ik = 0; Ik < tRequestedDvTypes.size(); Ik++ )
            {
                // get the set index for the requested master dv type
                sint tDvIndex = this->get_dv_index_for_type(
                        tRequestedDvTypes( Ik )( 0 ),
                        mtk::Master_Slave::MASTER );

                // if the index was set (and is different from -1)
                if( tDvIndex != -1 )
                {
                    // get the FI related to the master pdv type
                    Field_Interpolator * tFI = mMasterFIManager->
                            get_field_interpolators_for_type( tRequestedDvTypes( Ik )( 0 ) );

                    // get the number of coefficients related to the master dv type
                    uint tNumCoeff = tFI->get_number_of_space_time_coefficients();

                    // fill the dv assembly map with starting and ending indices
                    // for the master dv type
                    mPdvMatAssemblyMap( tDvIndex )( 0, 0 ) = tCounter;
                    mPdvMatAssemblyMap( tDvIndex )( 0, 1 ) = tCounter + tNumCoeff - 1;

                    // update the dv coefficient counter
                    tCounter += tNumCoeff;
                }
            }

            // loop over the slave dv types
            for( uint Ik = 0; Ik < tRequestedDvTypes.size(); Ik++ )
            {
                //get the set index for the slave dv type
                sint tDvIndex = this->get_dv_index_for_type(
                        tRequestedDvTypes( Ik )( 0 ),
                        mtk::Master_Slave::SLAVE );

                // if the dv type was set (its set index is different from -1)
                if( tDvIndex != -1 )
                {
                    // get the FI related to the master pdv type
                    Field_Interpolator * tFI = mSlaveFIManager->
                            get_field_interpolators_for_type( tRequestedDvTypes( Ik )( 0 ) );

                    // get the number of coefficients for the slave dv type
                    uint tNumCoeff = tFI->get_number_of_space_time_coefficients();

                    // fill the residual assembly map with starting and ending indices for the slave dof type
                    mPdvMatAssemblyMap( tDvIndex )( 0, 0 ) = tCounter;
                    mPdvMatAssemblyMap( tDvIndex )( 0, 1 ) = tCounter + tNumCoeff - 1;

                    // update the dof coefficient counter
                    tCounter += tNumCoeff;
                }
            }

            // set size for assembly vector
            mPdvMatAssemblyVector.set_size( tCounter, 1, -1 );

        }

        //--------------------------------------------------------------------------

        void Set::create_geo_pdv_assembly_map( std::shared_ptr< fem::Cluster > aFemCluster )
        {
            // get the design variable interface
            MSI::Design_Variable_Interface * tDVInterface =
                    mEquationModel->get_design_variable_interface();

            // get the geo dv types requested by the opt
            moris::Cell < enum PDV_Type > tRequestedDvTypes;
            tDVInterface->get_ig_unique_dv_types_for_set(
                    mMeshSet->get_set_index(),
                    tRequestedDvTypes );

            // get node indices on cluster
            moris::Matrix< moris::IndexMat > tNodeIndicesOnCluster;
            aFemCluster->get_vertex_indices_in_cluster_for_sensitivity( tNodeIndicesOnCluster );

            // clean up assembly map
            mPdvGeoAssemblyMap.clear();

            // clean up assembly vector
            Matrix< DDSMat > tPdvGeoAssemblyVectorTemp( tNodeIndicesOnCluster.numel() * tRequestedDvTypes.size(), 1, -1 );

            // init active geo pdv counter
            uint tActiveGeoPdvCounter = 0;

            // loop over the requested pdv types
            for( uint iGeoPdv = 0; iGeoPdv < tRequestedDvTypes.size(); iGeoPdv++ )
            {
                // get treated geo pdv type
                PDV_Type tGeoPdvType = tRequestedDvTypes( iGeoPdv );

                // build geo pdv type cell
                moris::Cell< PDV_Type > tGeoPdvMat = { tGeoPdvType };

                // loop over the ig nodes on cluster
                uint tNumIGNodes = tNodeIndicesOnCluster.numel();
                for( uint iIGNode = 0; iIGNode < tNumIGNodes; iIGNode++ )
                {
                    // get treated node index
                    moris_index tNodeIndex = tNodeIndicesOnCluster( iIGNode );

                    // build node index matrix
                    Matrix< IndexMat > tNodeIndexMat = {{ tNodeIndex }};

                    // check if active
                    moris::Cell< moris::Matrix< DDRMat > > tPdvValue( 1 );
                    tPdvValue( 0 ).set_size( 1, 1 );
                    moris::Cell< moris::Matrix< DDSMat > > tIsActivePdv;
                    tDVInterface->get_ig_pdv_value(
                            tNodeIndexMat,
                            tGeoPdvMat,
                            tPdvValue,
                            tIsActivePdv );

                    // create key pair
                    std::pair< moris_index, PDV_Type > tKeyPair = std::make_pair( tNodeIndex, tGeoPdvType );

                    // if active and not set in the map
                    if( tIsActivePdv( 0 )( 0 ) && ( mPdvGeoAssemblyMap.find( tKeyPair ) == mPdvGeoAssemblyMap.end() ) )
                    {
                        // fill the map
                        mPdvGeoAssemblyMap[ tKeyPair ] = tActiveGeoPdvCounter;

                        // get the id associated to the pdv
                        moris::Cell< moris::Matrix< DDSMat > > tPdvId;
                        tDVInterface->get_ig_dv_ids_for_type_and_ind(
                                tNodeIndexMat,
                                tGeoPdvMat,
                                tPdvId );

                        // fill the global assembly vector
                        tPdvGeoAssemblyVectorTemp( tActiveGeoPdvCounter ) = tPdvId( 0 )( 0 );

                        // update active geo pdv counter
                        tActiveGeoPdvCounter++;
                    }
                }
            }

            mPdvGeoAssemblyVector.set_size( 0, 0 );
            if( tActiveGeoPdvCounter > 0 )
            {
                // fill assembly vector
                mPdvGeoAssemblyVector.set_size( tActiveGeoPdvCounter, 1, -1 );
                mPdvGeoAssemblyVector = tPdvGeoAssemblyVectorTemp( { 0, tActiveGeoPdvCounter - 1 }, { 0, 0 } );
            }
        }

        //------------------------------------------------------------------------------

        void Set::create_requested_IWG_list()
        {
            // get list of requested dof types from solver
            moris::Cell < enum MSI::Dof_Type > tRequestedDofTypes =
                    this->get_model_solver_interface()->get_solver_interface()->get_requested_dof_types();

            // clear requested IWG list
            mRequestedIWGs.clear();

            // reserve max size for requested IWG list
            mRequestedIWGs.reserve( mIWGs.size() );

            if( mEquationModel->get_is_forward_analysis() )
            {
                // loop over the requested dof types
                for( MSI::Dof_Type tDofType : tRequestedDofTypes )
                {
                    // loop over the IWG in set IWG list
                    for( uint iIWG = 0; iIWG < mIWGs.size(); iIWG++ )
                    {
                        // if the IWG residual dof type is requested
                        if( mIWGs( iIWG )->get_residual_dof_type()( 0 ) == tDofType )
                        {
                            // add the IWg to the requested IWG list
                            mRequestedIWGs.push_back( mIWGs( iIWG ) );
                        }
                    }
                }
            }
            else
            {
                // loop over the requested dof types
                for( MSI::Dof_Type tDofType : tRequestedDofTypes )
                {
                    // loop over the IWG in set IWG list
                    for( uint iIWG = 0; iIWG < mIWGs.size(); iIWG++ )
                    {
                        // if the IWG residual dof type is requested
                        if( mIWGs( iIWG )->get_residual_dof_type()( 0 ) == tDofType )
                        {
                            if( mEquationModel->get_is_adjoint_off_diagonal_time_contribution() )
                            {
                                if( mIWGs( iIWG )->get_IWG_type() == moris::fem::IWG_Type::TIME_CONTINUITY_DOF )
                                {
                                    // add the IWg to the requested IWG list
                                    mRequestedIWGs.push_back( mIWGs( iIWG ) );
                                }
                            }
                            else
                            {
                                // add the IWg to the requested IWG list
                                mRequestedIWGs.push_back( mIWGs( iIWG ) );
                            }
                        }
                    }
                }
            }
            // reduce the size of requested IWG list to fit
            mRequestedIWGs.shrink_to_fit();
        }

        //------------------------------------------------------------------------------

        void Set::create_requested_IQI_list()
        {
            // create IQI map
            this->create_IQI_map();

            // clear requested IQI list
            mRequestedIQIs.clear();

            // Get names of potential requested IQIs
            const moris::Cell< std::string > & tRequestedIQINames =
                    mEquationModel->get_requested_IQI_names();

            // get number of potential requested IQIs
            uint tNumREquestedIQINames = tRequestedIQINames.size();

            // reserve memory
            mRequestedIQIs.reserve( tNumREquestedIQINames );

            // loop over requested IQI names
            for( uint Ik = 0; Ik < tNumREquestedIQINames; Ik++ )
            {
                // check if this set has the requested IQI
                if( mIQINameToIndexMap.key_exists( tRequestedIQINames( Ik ) ) )
                {
                    // get the set local index
                    moris_index tIQISetLocalIndex =
                            mIQINameToIndexMap.find( tRequestedIQINames( Ik ) );

                    // put IQI in requested IQI list
                    mRequestedIQIs.push_back( mIQIs( tIQISetLocalIndex ) );
                }
            }

            // reduce memory to used space
            mRequestedIQIs.shrink_to_fit();
        }

        //------------------------------------------------------------------------------

        void Set::create_IQI_map()
        {
            // erase the content of the map
            mIQINameToIndexMap.clear();

            uint tCounter = 0;

            // loop over all IQIs and build a name to index map
            for( auto tIQI : mIQIs )
            {
                std::string tIQIName = tIQI->get_name();

                mIQINameToIndexMap[ tIQIName ] = tCounter++;
            }
        }

        //------------------------------------------------------------------------------

        void Set::build_requested_IWG_dof_type_list( const bool aIsResidual )
        {
            for( auto tIWG : mRequestedIWGs )
            {
                tIWG->build_requested_dof_type_list( aIsResidual );
            }
        }

        //------------------------------------------------------------------------------

        void Set::build_requested_IQI_dof_type_list()
        {
            for( auto tIQI : mRequestedIQIs )
            {
                tIQI->build_requested_dof_type_list();
            }
        }

        //------------------------------------------------------------------------------

        void Set::initialize_mJacobian()
        {
            // if residual not initialized before
            if ( !mJacobianExist )
            {
                // get the dof types requested by the solver
                moris::Cell< enum MSI::Dof_Type > tRequestedDofTypes =
                        this->get_requested_dof_types();

                // init dof coefficient counter
                uint tNumCols = 0;

                // loop over the requested dof types
                for ( uint Ik = 0; Ik < tRequestedDofTypes.size(); Ik++ )
                {
                    // get the set index for the master dof type
                    sint tDofIndex = this->get_dof_index_for_type(
                            tRequestedDofTypes( Ik ),
                            mtk::Master_Slave::MASTER);

                    // if this master dof is active
                    if ( tDofIndex != -1 )
                    {
                        // update number of dof coefficients
                        tNumCols += mMasterFIManager->
                                get_field_interpolators_for_type( tRequestedDofTypes( Ik ) )->
                                get_number_of_space_time_coefficients();
                    }

                    // get the set index for the slave dof type
                    tDofIndex = this->get_dof_index_for_type(
                            tRequestedDofTypes(Ik),
                            mtk::Master_Slave::SLAVE );

                    // if this slave dof is active
                    if ( tDofIndex != -1 )
                    {
                        // update number of dof coefficients
                        tNumCols += mSlaveFIManager->
                                get_field_interpolators_for_type( tRequestedDofTypes( Ik ) )->
                                get_number_of_space_time_coefficients();
                    }
                }

                // if for residual evaluation
                if( !mIsResidual )
                {
                    // set size for the jacobian matrix
                    mJacobian.set_size( tNumCols, tNumCols, 0.0 );
                }
                // if for Jacobian evaluation
                else
                {
                    // get the secondary dof types from the solver
                    moris::Cell< moris::Cell< enum MSI::Dof_Type > > tSecDofTypes =
                            this->get_secondary_dof_types();

                    // init dof coefficient counter for rows
                    uint tNumRows = 0;

                    // loop over the groups of secondary dof types
                    for ( auto tSecDofTypesI : tSecDofTypes )
                    {
                        // loop over the secondary dof types in group
                        for (uint Ik = 0; Ik < tSecDofTypesI.size(); Ik++)
                        {
                            // get the set index for the master dof type
                            sint tDofIndex = this->get_dof_index_for_type(
                                    tSecDofTypesI( Ik ),
                                    mtk::Master_Slave::MASTER );

                            // if this master dof is active
                            if ( tDofIndex != -1 )
                            {
                                // update number of dof coefficients
                                tNumRows += mMasterFIManager->
                                        get_field_interpolators_for_type( tSecDofTypesI( Ik ) )->
                                        get_number_of_space_time_coefficients();
                            }

                            // get the set index for the slave dof type
                            tDofIndex = this->get_dof_index_for_type(
                                    tSecDofTypesI(Ik),
                                    mtk::Master_Slave::SLAVE );

                            // if this slave dof is active
                            if ( tDofIndex != -1 )
                            {
                                // update number of dof coefficients
                                tNumRows += mSlaveFIManager->
                                        get_field_interpolators_for_type( tSecDofTypesI( Ik ) )->
                                        get_number_of_space_time_coefficients();
                            }
                        }
                    }

                    // set size for the jacobian matrix
                    mJacobian.set_size( tNumCols, tNumRows, 0.0 );
                }
                // set the jacobian initialization flag to true
                mJacobianExist = true;
            }
            else
            {
                // fill the jacobian matrix with zeros
                mJacobian.fill( 0.0 );
            }
        }

        //------------------------------------------------------------------------------

        void Set::initialize_mResidual()
        {
            // if residual not initialized before
            if ( !mResidualExist )
            {
                if( !mIsResidual )
                {
                    // get the dof types requested by the solver
                    moris::Cell < enum MSI::Dof_Type >tRequestedDofTypes =
                            this->get_requested_dof_types();

                    // init dof coefficient counter
                    uint tNumCoeff = 0;

                    // loop over the requested dof types
                    for( uint Ik = 0; Ik < tRequestedDofTypes.size(); Ik++ )
                    {
                        // get the set index for the master dof type
                        sint tDofIndex = this->get_dof_index_for_type(
                                tRequestedDofTypes( Ik ),
                                mtk::Master_Slave::MASTER );

                        // if this master dof is active
                        if( tDofIndex != -1 )
                        {
                            // update number of dof coefficients
                            tNumCoeff += mMasterFIManager->
                                    get_field_interpolators_for_type( tRequestedDofTypes( Ik ) )->
                                    get_number_of_space_time_coefficients();
                        }

                        // get the set index for the slave dof type
                        tDofIndex = this->get_dof_index_for_type(
                                tRequestedDofTypes( Ik ),
                                mtk::Master_Slave::SLAVE  );

                        // if this slave dof is active
                        if( tDofIndex != -1 )
                        {
                            // update number of dof coefficients
                            tNumCoeff += mSlaveFIManager->
                                    get_field_interpolators_for_type( tRequestedDofTypes( Ik ) )->
                                    get_number_of_space_time_coefficients();
                        }
                    }

                    // get the number of rhs
                    uint tNumRHS = mEquationModel->get_num_rhs();

                    // set size for the list of dQIdu vectors
                    mResidual.resize( tNumRHS );

                    // loop over the dQIdu vectors
                    for( auto & tRes : mResidual )
                    {
                        // set size for the dQIdu vector
                        tRes.set_size( tNumCoeff, 1, 0.0 );
                    }
                }
                else
                {
                    // get the dof types requested by the solver FIXME
                    moris::Cell < enum MSI::Dof_Type >tRequestedDofTypes =
                            this->get_secondary_dof_types()(0);

                    // init dof coefficient counter
                    uint tNumCoeff = 0;

                    // loop over the requested dof types
                    for( uint Ik = 0; Ik < tRequestedDofTypes.size(); Ik++ )
                    {
                        // get the set index for the master dof type
                        sint tDofIndex = this->get_dof_index_for_type(
                                tRequestedDofTypes( Ik ),
                                mtk::Master_Slave::MASTER );

                        // if this master dof is active
                        if( tDofIndex != -1 )
                        {
                            // update number of dof coefficients
                            tNumCoeff += mMasterFIManager->
                                    get_field_interpolators_for_type( tRequestedDofTypes( Ik ) )->
                                    get_number_of_space_time_coefficients();
                        }

                        // get the set index for the slave dof type
                        tDofIndex = this->get_dof_index_for_type(
                                tRequestedDofTypes( Ik ),
                                mtk::Master_Slave::SLAVE  );

                        // if this slave dof is active
                        if( tDofIndex != -1 )
                        {
                            // update number of dof coefficients
                            tNumCoeff += mSlaveFIManager->
                                    get_field_interpolators_for_type( tRequestedDofTypes( Ik ) )->
                                    get_number_of_space_time_coefficients();
                        }
                    }

                    // get the number of rhs
                    uint tNumRHS = mEquationModel->get_num_rhs();

                    // set size for the list of dQIdu vectors
                    mResidual.resize( tNumRHS );

                    // loop over the dQIdu vectors
                    for( auto & tRes : mResidual )
                    {
                        // set size for the dQIdu vector
                        tRes.set_size( tNumCoeff, 1, 0.0 );
                    }
                }

                // set the residual initialization flag to true
                mResidualExist = true;
            }
            // if residual initialized before
            else
            {
                // loop over the residual vectors
                for( auto & tRes : mResidual )
                {
                    // fill the residual vector with zeros
                    tRes.fill( 0.0 );
                }
            }
        }

        //------------------------------------------------------------------------------

        void Set::initialize_mQI()
        {
            // if list of QI values not initialized before
            if ( !mQIExist )
            {
                uint tNumQI = mEquationModel->get_requested_IQI_names().size();

                // set size for the list of QI values
                mQI.resize( tNumQI );

                for( auto & tQI : mQI )
                {
                    // set size for the QI value
                    // FIXME assumed scalar
                    tQI.set_size( 1, 1, 0.0 );
                }

                // set the QI initialization flag to true
                mQIExist = true;
            }
            // if list of QI values initialized before
            else
            {
                // loop over the QI values
                for( auto & tQI : mQI )
                {
                    // fill the QI value vector with zero
                    tQI.fill( 0.0 );
                }
            }
        }

        //------------------------------------------------------------------------------

        void Set::initialize_mdQIdpMat()
        {
            // if dRdpMap not initialized before
            if ( !mdQIdpMatExist )
            {
                // set size for dQIdp
                mdQIdp.resize( 2 );

                // get the number of requested IQIs on the model
                uint tNumRequestedIQIs = this->get_equation_model()->get_requested_IQI_names().size();

                // set size for dQIdp
                mdQIdp( 0 ).resize( tNumRequestedIQIs );

                // get the requested pdv types
                moris::Cell< moris::Cell < enum PDV_Type > > tRequestedDvTypes;
                this->get_ip_dv_types_for_set( tRequestedDvTypes );

                // get the number of requested pdv types
                uint tNumRequestedPdvTypes = tRequestedDvTypes.size();

                // init pdv coefficient counter
                uint tNumPdvCoefficients = 0;

                // loop over the requested dv types
                for( uint Ik = 0; Ik < tNumRequestedPdvTypes; Ik++ )
                {
                    // get the set index for the master dof type
                    sint tDvIndex = this->get_dv_index_for_type(
                            tRequestedDvTypes( Ik )( 0 ),
                            mtk::Master_Slave::MASTER );

                    // if this master dv is active
                    if( tDvIndex != -1 )
                    {
                        // update number of dof coefficients
                        tNumPdvCoefficients += mMasterFIManager->
                                get_field_interpolators_for_type( tRequestedDvTypes( Ik )( 0 ) )->
                                get_number_of_space_time_coefficients();
                    }

                    // get the set index for the slave dv type
                    tDvIndex = this->get_dv_index_for_type(
                            tRequestedDvTypes( Ik )( 0 ),
                            mtk::Master_Slave::SLAVE  );

                    // if this slave dv is active
                    if( tDvIndex != -1 )
                    {
                        // update number of dof coefficients
                        tNumPdvCoefficients += mSlaveFIManager->
                                get_field_interpolators_for_type( tRequestedDvTypes( Ik )( 0 ) )->
                                get_number_of_space_time_coefficients();
                    }
                }
                // loop over requested IQIs
                for ( auto & tdQIdp : mdQIdp( 0 ) )
                {
                    // set size
                    tdQIdp.set_size( 1, tNumPdvCoefficients, 0.0 );
                }

                // set exist flag to true
                mdQIdpMatExist = true;
            }
            else
            {
                // loop over requested IQIs
                for ( auto & tdQIdp : mdQIdp( 0 ) )
                {
                    // fill the dQIdp vector with zero
                    tdQIdp.fill( 0.0 );
                }
            }
        }

        //------------------------------------------------------------------------------

        void Set::initialize_mdQIdpGeo( std::shared_ptr< fem::Cluster > aFemCluster )
        {
            // get the number of requested IQIs
            uint tNumRequestedIQIs = this->get_equation_model()->get_requested_IQI_names().size();

            // set size for dQIdp
            mdQIdp( 1 ).resize( tNumRequestedIQIs );

            // loop over requested IQIs
            for ( auto & tdQIdp : mdQIdp( 1 ) )
            {
                // fill the dQIdp vector with zero
                tdQIdp.set_size( 1, mPdvGeoAssemblyVector.numel(), 0.0 );
            }
        }

        //----------------------------------------------------------------------

        void Set::initialize_mdRdpMat()
        {
            // if dRdpMap not initialized before
            if ( !mdRdpMatExist )
            {
                // set size for dRdp
                mdRdp.resize( 2 );

                // get the dof types requested by the solver
                moris::Cell < enum MSI::Dof_Type >tRequestedDofTypes =
                        this->get_requested_dof_types();

                // init dof coefficient counter
                uint tNumRows = 0;

                // loop over the requested dof types
                for( uint Ik = 0; Ik < tRequestedDofTypes.size(); Ik++ )
                {
                    // get the set index for the master dof type
                    sint tDofIndex = this->get_dof_index_for_type(
                            tRequestedDofTypes( Ik ),
                            mtk::Master_Slave::MASTER );

                    // if this master dof is active
                    if( tDofIndex != -1 )
                    {
                        // update number of dof coefficients
                        tNumRows += mMasterFIManager->
                                get_field_interpolators_for_type( tRequestedDofTypes( Ik ) )->
                                get_number_of_space_time_coefficients();
                    }

                    // get the set index for the slave dof type
                    tDofIndex = this->get_dof_index_for_type(
                            tRequestedDofTypes( Ik ),
                            mtk::Master_Slave::SLAVE  );

                    // if this slave dof is active
                    if( tDofIndex != -1 )
                    {
                        // update number of dof coefficients
                        tNumRows += mSlaveFIManager->
                                get_field_interpolators_for_type( tRequestedDofTypes( Ik ) )->
                                get_number_of_space_time_coefficients();
                    }
                }

                // get the dv types requested by the opt
                moris::Cell< moris::Cell < enum PDV_Type > > tRequestedDvTypes;
                this->get_ip_dv_types_for_set( tRequestedDvTypes );

                // init dv coefficient counter
                uint tNumCols = 0;

                // loop over the requested dv types
                for( uint Ik = 0; Ik < tRequestedDvTypes.size(); Ik++ )
                {
                    // get the set index for the master dof type
                    sint tDvIndex = this->get_dv_index_for_type(
                            tRequestedDvTypes( Ik )( 0 ),
                            mtk::Master_Slave::MASTER );

                    // if this master dv is active
                    if( tDvIndex != -1 )
                    {
                        // update number of dof coefficients
                        tNumCols += mMasterFIManager->
                                get_field_interpolators_for_type( tRequestedDvTypes( Ik )( 0 ) )->
                                get_number_of_space_time_coefficients();
                    }

                    // get the set index for the slave dv type
                    tDvIndex = this->get_dv_index_for_type(
                            tRequestedDvTypes( Ik )( 0 ),
                            mtk::Master_Slave::SLAVE  );

                    // if this slave dv is active
                    if( tDvIndex != -1 )
                    {
                        // update number of dof coefficients
                        tNumCols += mSlaveFIManager->
                                get_field_interpolators_for_type( tRequestedDvTypes( Ik )( 0 ) )->
                                get_number_of_space_time_coefficients();
                    }
                }

                mdRdp( 0 ).set_size( tNumRows, tNumCols, 0.0 );

                // set the dRdpMat initialization flag to true
                mdRdpMatExist = true;
            }
            // if dRdpMat initialized before
            else
            {
                mdRdp( 0 ).fill( 0.0 );
            }
        }

        //----------------------------------------------------------------------

        void Set::initialize_mdRdpGeo( std::shared_ptr< fem::Cluster > aFemCluster )
        {
            // get the dof types requested by the solver
            moris::Cell < enum MSI::Dof_Type >tRequestedDofTypes =
                    this->get_requested_dof_types();

            // init dof coefficient counter
            uint tNumRows = 0;

            // loop over the requested dof types
            for( uint Ik = 0; Ik < tRequestedDofTypes.size(); Ik++ )
            {
                // get the set index for the master dof type
                sint tDofIndex = this->get_dof_index_for_type(
                        tRequestedDofTypes( Ik ),
                        mtk::Master_Slave::MASTER );

                // if this master dof is active
                if( tDofIndex != -1 )
                {
                    // update number of dof coefficients
                    tNumRows += mMasterFIManager->
                            get_field_interpolators_for_type( tRequestedDofTypes( Ik ) )->
                            get_number_of_space_time_coefficients();
                }

                // get the set index for the slave dof type
                tDofIndex = this->get_dof_index_for_type(
                        tRequestedDofTypes( Ik ),
                        mtk::Master_Slave::SLAVE );

                // if this slave dof is active
                if( tDofIndex != -1 )
                {
                    // update number of dof coefficients
                    tNumRows += mSlaveFIManager->
                            get_field_interpolators_for_type( tRequestedDofTypes( Ik ) )->
                            get_number_of_space_time_coefficients();
                }
            }

            // set size for dRdpgeo
            mdRdp( 1 ).set_size( tNumRows, mPdvGeoAssemblyVector.numel(), 0.0 );
        }

        //------------------------------------------------------------------------------

        mtk::Interpolation_Order Set::get_auto_interpolation_order(
                const moris::uint        aNumVertices,
                const mtk::Geometry_Type aGeometryType )
        {
            switch( aGeometryType )
            {
                case mtk::Geometry_Type::LINE:
                {
                    switch( aNumVertices )
                    {
                        case 1:
                            return mtk::Interpolation_Order::UNDEFINED;

                        case 2:
                            return mtk::Interpolation_Order::LINEAR;

                        case 3:
                            return mtk::Interpolation_Order::QUADRATIC;

                        default:
                            MORIS_ERROR( false, " Set::get_auto_interpolation_order - not defined for LINE and number of vertices. ");
                            return mtk::Interpolation_Order::UNDEFINED;
                    }
                }
                case mtk::Geometry_Type::QUAD:
                {
                    switch( aNumVertices )
                    {
                        case 4:
                            return mtk::Interpolation_Order::LINEAR;

                        case 8:
                            return mtk::Interpolation_Order::SERENDIPITY;

                        case 9:
                            return mtk::Interpolation_Order::QUADRATIC;

                        case 16:
                            return mtk::Interpolation_Order::CUBIC;

                        default:
                            MORIS_ERROR( false, " Set::get_auto_interpolation_order - not defined for QUAD and number of vertices. ");
                            return mtk::Interpolation_Order::UNDEFINED;
                    }
                }
                case mtk::Geometry_Type::HEX:
                {
                    switch( aNumVertices )
                    {
                        case 8:
                            return mtk::Interpolation_Order::LINEAR;

                        case 20:
                            return mtk::Interpolation_Order::SERENDIPITY;

                        case 27:
                            return mtk::Interpolation_Order::QUADRATIC;

                        case 64:
                            return mtk::Interpolation_Order::CUBIC;

                        default:
                            MORIS_ERROR( false, " Set::get_auto_interpolation_order - not defined for HEX and number of vertices. ");
                            return mtk::Interpolation_Order::UNDEFINED;
                    }
                }
                case mtk::Geometry_Type::TET:
                {
                    switch( aNumVertices )
                    {
                        case 4:
                            return mtk::Interpolation_Order::LINEAR;

                        case 10:
                            return mtk::Interpolation_Order::QUADRATIC;

                        case 20:
                            return mtk::Interpolation_Order::CUBIC;

                        default:
                            MORIS_ERROR( false, " Set::get_auto_interpolation_order - not defined for TET and number of vertices. ");
                            return mtk::Interpolation_Order::UNDEFINED;
                    }
                }
                default:
                {
                    MORIS_ERROR( false, " Set::get_auto_interpolation_order - not defined for this geometry type. ");
                    return mtk::Interpolation_Order::UNDEFINED;
                }
            }
        }

        //------------------------------------------------------------------------------

        fem::Interpolation_Type Set::get_auto_time_interpolation_type(
                const moris::uint aNumVertices )
        {
            switch( aNumVertices )
            {
                case 1:
                    return Interpolation_Type::CONSTANT;

                case 2:
                case 3:
                case 4:
                    return Interpolation_Type::LAGRANGE;

                default:
                    MORIS_ERROR( false, " Element::get_auto_time_interpolation_type - not defined this number of time vertices. ");
                    return Interpolation_Type::UNDEFINED;
            }
        }

        //------------------------------------------------------------------------------

        fem::Integration_Order Set::get_auto_integration_order(
                const fem::Element_Type        aSetType,
                const mtk::Geometry_Type       aGeometryType,
                const mtk::Interpolation_Order aInterpolationOrder )
        {
            switch( aSetType )
            {
                case fem::Element_Type::BULK :
                case fem::Element_Type::TIME_SIDESET :
                case fem::Element_Type::TIME_BOUNDARY :
                {
                    switch( aGeometryType )
                    {
                        case mtk::Geometry_Type::LINE:
                        {
                            switch( aInterpolationOrder )
                            {
                                case mtk::Interpolation_Order::LINEAR:
                                    return fem::Integration_Order::BAR_3;

                                case mtk::Interpolation_Order::QUADRATIC:
                                    return fem::Integration_Order::BAR_4;

                                case mtk::Interpolation_Order::CUBIC:
                                    return fem::Integration_Order::BAR_5;

                                default:
                                    MORIS_ERROR( false, "Set::get_auto_integration_order - Unknown or unsupported interpolation order.");
                                    return fem::Integration_Order::UNDEFINED;
                            }
                            break;
                        }

                        case mtk::Geometry_Type::QUAD :
                        {
                            switch( aInterpolationOrder )
                            {
                                case mtk::Interpolation_Order::LINEAR:
                                    return fem::Integration_Order::QUAD_2x2;

                                case mtk::Interpolation_Order::SERENDIPITY:
                                    return fem::Integration_Order::QUAD_3x3;

                                case mtk::Interpolation_Order::QUADRATIC:
                                    return fem::Integration_Order::QUAD_3x3;

                                case mtk::Interpolation_Order::CUBIC:
                                    return fem::Integration_Order::QUAD_5x5;

                                default:
                                    MORIS_ERROR( false, "Set::get_auto_integration_order - Unknown or unsupported interpolation order.");
                                    return fem::Integration_Order::UNDEFINED;
                            }
                            break;
                        }

                        case mtk::Geometry_Type::HEX:
                        {
                            switch( aInterpolationOrder )
                            {
                                case mtk::Interpolation_Order::LINEAR:
                                    return fem::Integration_Order::HEX_2x2x2;

                                case mtk::Interpolation_Order::SERENDIPITY:
                                    return fem::Integration_Order::HEX_4x4x4;

                                case mtk::Interpolation_Order::QUADRATIC:
                                    return fem::Integration_Order::HEX_4x4x4;

                                case mtk::Interpolation_Order::CUBIC:
                                    return fem::Integration_Order::HEX_5x5x5;

                                default:
                                    MORIS_ERROR( false, "Set::get_auto_integration_order - Unknown or unsupported interpolation order.");
                                    return fem::Integration_Order::UNDEFINED;
                            }
                            break;
                        }

                        case mtk::Geometry_Type::TRI:
                        {
                            switch( aInterpolationOrder )
                            {
                                case mtk::Interpolation_Order::LINEAR:
                                    return fem::Integration_Order::TRI_7;

                                case mtk::Interpolation_Order::QUADRATIC:
                                    return fem::Integration_Order::TRI_12;

                                case mtk::Interpolation_Order::CUBIC:
                                    return fem::Integration_Order::TRI_25;

                                default:
                                    MORIS_ERROR( false, "Set::get_auto_integration_order - Unknown or unsupported interpolation order.");
                                    return fem::Integration_Order::UNDEFINED;
                            }
                            break;
                        }

                        case mtk::Geometry_Type::TET:
                        {
                            switch( aInterpolationOrder )
                            {
                                case mtk::Interpolation_Order::LINEAR:
                                    return fem::Integration_Order::TET_11;

                                case mtk::Interpolation_Order::QUADRATIC:
                                    return fem::Integration_Order::TET_35;

                                case mtk::Interpolation_Order::CUBIC:
                                    return fem::Integration_Order::TET_56;

                                default:
                                    MORIS_ERROR( false, "Set::get_auto_integration_order - Unknown or unsupported interpolation order.");
                                    return fem::Integration_Order::UNDEFINED;
                            }
                            break;
                        }

                        default :
                            MORIS_ERROR( false, " Set::get_auto_integration_order - Unknown or unsupported geometry type. ");
                            return Integration_Order::UNDEFINED;
                    }
                    break;
                }
                case fem::Element_Type::SIDESET :
                case fem::Element_Type::DOUBLE_SIDESET :
                {
                    switch( aGeometryType )
                    {
                        case mtk::Geometry_Type::LINE:
                        {
                            switch( aInterpolationOrder )
                            {
                                case mtk::Interpolation_Order::LINEAR:
                                    return fem::Integration_Order::BAR_3;

                                case mtk::Interpolation_Order::QUADRATIC:
                                    return fem::Integration_Order::BAR_4;

                                case mtk::Interpolation_Order::CUBIC:
                                    return fem::Integration_Order::BAR_5;

                                default:
                                    MORIS_ERROR( false, "Set::get_auto_integration_order - Unknown or unsupported interpolation order.");
                                    return fem::Integration_Order::UNDEFINED;
                            }
                            break;
                        }

                        case mtk::Geometry_Type::QUAD :
                        {
                            switch( aInterpolationOrder )
                            {
                                case mtk::Interpolation_Order::LINEAR:
                                    return fem::Integration_Order::QUAD_2x2;

                                case mtk::Interpolation_Order::SERENDIPITY:
                                    return fem::Integration_Order::QUAD_4x4;

                                case mtk::Interpolation_Order::QUADRATIC:
                                    return fem::Integration_Order::QUAD_4x4;

                                case mtk::Interpolation_Order::CUBIC:
                                    return fem::Integration_Order::QUAD_5x5;

                                default:
                                    MORIS_ERROR( false, "Set::get_auto_integration_order - Unknown or unsupported interpolation order.");
                                    return fem::Integration_Order::UNDEFINED;
                            }
                            break;
                        }

                        case mtk::Geometry_Type::TRI:
                        {
                            switch( aInterpolationOrder )
                            {
                                case mtk::Interpolation_Order::LINEAR:
                                    return fem::Integration_Order::TRI_7;

                                case mtk::Interpolation_Order::QUADRATIC:
                                    return fem::Integration_Order::TRI_12;

                                case mtk::Interpolation_Order::CUBIC:
                                    return fem::Integration_Order::TRI_25;

                                default:
                                    MORIS_ERROR( false, "Set::get_auto_integration_order - Unknown or unsupported interpolation order.");
                                    return fem::Integration_Order::UNDEFINED;
                            }
                            break;
                        }

                        default :
                            MORIS_ERROR( false, " Set::get_auto_integration_order - Unknown or unsupported geometry type. ");
                            return Integration_Order::UNDEFINED;
                    }
                    break;
                }

                default:
                    MORIS_ERROR( false, "Set::get_auto_integration_order - unknown set type.");
                    return fem::Integration_Order::UNDEFINED;
            }
        }

        //------------------------------------------------------------------------------

        void Set::set_visualization_set(
                const uint         aMeshIndex,
                moris::mtk::Set  * aVisMeshSet,
                const bool         aOnlyPrimayCells )
        {
            uint tNumClustersOnSets = aVisMeshSet->get_num_clusters_on_set();

            // set vis clusters to clusters
            for( uint Ik = 0; Ik < tNumClustersOnSets; Ik++ )
            {
                // create a fem cluster
                std::shared_ptr< fem::Cluster > tCluster =
                        std::make_shared< fem::Cluster >(
                                mElementType,
                                aVisMeshSet->get_clusters_by_index( Ik ),
                                this,
                                mEquationObjList( Ik ) );

                reinterpret_cast< fem::Interpolation_Element* >( mEquationObjList( Ik ) )->set_cluster(
                        tCluster,
                        aMeshIndex );
            }

            // build set element map
            uint tNumCells = aVisMeshSet->get_num_cells_on_set( aOnlyPrimayCells );

            moris::Matrix< DDSMat > tCellIndex = aVisMeshSet->get_cell_inds_on_block( aOnlyPrimayCells );

            sint tSize = std::max( ( sint )mCellAssemblyMap.size(), ( sint )aMeshIndex + 1 );

            mCellAssemblyMap.resize( tSize );
            mMtkIgCellOnSet .resize( tSize );

            mMtkIgCellOnSet( aMeshIndex )= tNumCells;

            if(tNumCells>0)
            {
                sint tMaxIndex = tCellIndex.max();
                //             sint tMinIndex = tCellIndex.min();

                mCellAssemblyMap( aMeshIndex ).set_size( tMaxIndex + 1, 1, -1 );

                for( uint Ik = 0; Ik < tNumCells; Ik++ )
                {
                    mCellAssemblyMap( aMeshIndex )( tCellIndex( Ik ) ) = Ik;
                }
            }
        }

        //------------------------------------------------------------------------------

        void Set::compute_quantity_of_interest(
                const uint              aMeshIndex,
                Matrix< DDRMat >      * aElementFieldValues,
                Matrix< DDRMat >      * aNodalFieldValues,
                moris::real           * aGlobalScalar,
                const std::string     & aQIName,
                enum vis::Field_Type    aFieldType )
        {
            mSetElementalValues = aElementFieldValues;
            mSetNodalValues     = aNodalFieldValues;
            mSetGlobalValues    = aGlobalScalar;

            mSetNodalCounter.set_size( (*mSetNodalValues).numel(), 1, 0 );

            mSetElementalValues->set_size( mMtkIgCellOnSet( aMeshIndex ), 1, 0.0 );

            // check if this set has the requested IQI
            if( mIQINameToIndexMap.key_exists( aQIName ) )
            {
                for( uint Ik = 0; Ik < mEquationObjList.size(); Ik++ )
                {
                    mEquationObjList( Ik )->compute_quantity_of_interest(
                            aMeshIndex,
                            aQIName,
                            aFieldType );
                }

                //FIXME I do not like this at all. someone change it
                for( uint Ik = 0; Ik < mSetNodalValues->numel(); Ik++ )
                {
                    if( mSetNodalCounter(Ik) != 0)
                    {
                        (*mSetNodalValues)(Ik) = (*mSetNodalValues)(Ik)/mSetNodalCounter(Ik);
                    }
                }
            }
        }

        //------------------------------------------------------------------------------

        void Set::determine_set_type()
        {
            enum moris::SetType tMtkSetType = mMeshSet->get_set_type();

            switch( tMtkSetType )
            {
                case moris::SetType::BULK:
                {
                    mElementType = fem::Element_Type::BULK;

                    // if time continuity
                    if ( mTimeContinuity )
                    {
                        mElementType = fem::Element_Type::TIME_SIDESET;
                    }

                    // if time boundary
                    if( mTimeBoundary )
                    {
                        mElementType = fem::Element_Type::TIME_BOUNDARY;
                    }
                    break;
                }

                case moris::SetType::SIDESET:
                    mElementType = fem::Element_Type::SIDESET;
                    break;

                case moris::SetType::DOUBLE_SIDED_SIDESET:
                    mElementType = fem::Element_Type::DOUBLE_SIDESET;
                    break;

                default :
                    MORIS_ERROR( false, "Set::determine_set_type() - not defined for this set type. ");
            }
        }

        //------------------------------------------------------------------------------

        void Set::get_ig_unique_dv_types_for_set(
                moris::Cell < enum PDV_Type > & aGeoPdvType )
        {
            // get design variable interface
            MSI::Design_Variable_Interface * tPdvInterface =
                    mEquationModel->get_design_variable_interface();

            // if the pdv interface is set
            if( tPdvInterface )
            {
                // get ig unique pdv types for set
                tPdvInterface->get_ig_unique_dv_types_for_set(
                        mMeshSet->get_set_index(),
                        aGeoPdvType );
            }
        }

        //------------------------------------------------------------------------------

        void Set::get_ip_dv_types_for_set(
                moris::Cell< moris::Cell< enum PDV_Type > > & aMatPdvType )
        {
            //
            aMatPdvType = mMasterDvTypes;

            // FIXME
            //            // get design variable interface
            //            MSI::Design_Variable_Interface * tPdvInterface =
            //                    mEquationModel->get_design_variable_interface();
            //
            //            // get ip pdv types for set
            //            tPdvInterface->get_ip_dv_types_for_set(
            //                    mMeshSet->get_set_index(),
            //                    aMatPdvType );
        }

        //------------------------------------------------------------------------------

    } /* namespace fem */
} /* namespace moris */
