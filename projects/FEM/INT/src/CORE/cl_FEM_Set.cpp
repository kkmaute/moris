/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_FEM_Set.cpp
 *
 */
#include <iostream>

#include "cl_MSI_Model_Solver_Interface.hpp"        //FEM/MSI/src
#include "cl_MSI_Solver_Interface.hpp"              //FEM/MSI/src
#include "cl_FEM_Model.hpp"                         //FEM/INT/src
#include "cl_FEM_Set.hpp"                           //FEM/INT/src
#include "cl_FEM_Set_User_Info.hpp"                 //FEM/INT/src
#include "cl_FEM_Element_Factory.hpp"               //FEM/INT/src
#include "cl_MTK_Integrator.hpp"                    //MTK/src
#include "cl_FEM_Field_Interpolator_Manager.hpp"    //FEM/INT/src
#include "cl_FEM_Interpolation_Element.hpp"         //FEM/INT/src
#include "cl_FEM_Cluster.hpp"                       //FEM/INT/src
#include "cl_MTK_Set.hpp"                           //FEM/INT/src
#include "cl_MTK_Cell_Info.hpp"
#include "fn_equal_to.hpp"

namespace moris
{
    namespace fem
    {
        //------------------------------------------------------------------------------

        Set::Set(
                fem::FEM_Model*                  aFemModel,
                moris::mtk::Set*                 aMeshSet,
                const fem::Set_User_Info&        aSetInfo,
                const moris::Cell< Node_Base* >& aIPNodes )
                : mFemModel( aFemModel )
                , mMeshSet( aMeshSet )
                , mIWGs( aSetInfo.get_IWGs() )
                , mIQIs( aSetInfo.get_IQIs() )
                , mTimeContinuity( aSetInfo.get_time_continuity() )
                , mIsAnalyticalFA( aSetInfo.get_is_analytical_forward_analysis() )
                , mFDSchemeForFA( aSetInfo.get_finite_difference_scheme_for_forward_analysis() )
                , mFDPerturbationFA( aSetInfo.get_finite_difference_perturbation_size_for_forward_analysis() )
                , mIsAnalyticalSA( aSetInfo.get_is_analytical_sensitivity_analysis() )
                , mFDSchemeForSA( aSetInfo.get_finite_difference_scheme_for_sensitivity_analysis() )
                , mFDPerturbation( aSetInfo.get_finite_difference_perturbation_size() )
                , mPerturbationStrategy( aSetInfo.get_perturbation_strategy() )
        {
            // get the set type (BULK, SIDESET, DOUBLE_SIDESET, TIME_SIDESET)
            this->determine_set_type();

            // loop over the IWGs on the set
            for ( const std::shared_ptr< IWG >& tIWG : mIWGs )
            {
                // set the fem set pointer to the IWG
                tIWG->set_set_pointer( this );
            }

            // loop over the IQIs on the set
            for ( const std::shared_ptr< IQI >& tIQI : mIQIs )
            {
                // set the fem set pointer to the IQI
                tIQI->set_set_pointer( this );
            }

            // get mesh clusters on set
            moris::Cell< mtk::Cluster const * > tMeshClusterList = mMeshSet->get_clusters_on_set();

            // get number of mesh clusters on set
            uint tNumMeshClusters = tMeshClusterList.size();

            // set size for the equation objects list
            mEquationObjList.resize( tNumMeshClusters, nullptr );

            // get cluster measures used on set
            this->build_cluster_measure_tuples_and_map();

            // create a fem cluster factory
            fem::Element_Factory tClusterFactory;

            // loop over mesh clusters on set
            for ( luint iCluster = 0; iCluster < tNumMeshClusters; iCluster++ )
            {
                // init list of pointers to IP mesh cell
                moris::Cell< const mtk::Cell* > tInterpolationCell;

                // switch on set type
                switch ( mElementType )
                {
                    // if bulk or sideset
                    case fem::Element_Type::BULK:
                    case fem::Element_Type::SIDESET:
                    case fem::Element_Type::TIME_SIDESET:
                    case fem::Element_Type::TIME_BOUNDARY:
                    {
                        tInterpolationCell.resize( 1, &tMeshClusterList( iCluster )->get_interpolation_cell() );
                        break;
                    }
                    // if double sideset
                    case fem::Element_Type::DOUBLE_SIDESET:
                    {
                        tInterpolationCell.resize( 2 );
                        tInterpolationCell( 0 ) = &tMeshClusterList( iCluster )->get_interpolation_cell( mtk::Leader_Follower::LEADER );
                        tInterpolationCell( 1 ) = &tMeshClusterList( iCluster )->get_interpolation_cell( mtk::Leader_Follower::FOLLOWER );
                        break;
                    }
                    // if none of the above
                    default:
                    {
                        MORIS_ERROR( false, "Set::Set - unknown element type" );
                    }
                }

                // create an interpolation element
                mEquationObjList( iCluster ) = new fem::Interpolation_Element(
                        mElementType,
                        tInterpolationCell,
                        aIPNodes,
                        this );

                // create a fem cluster
                std::shared_ptr< fem::Cluster > tCluster = std::make_shared< fem::Cluster >(
                        mElementType,
                        tMeshClusterList( iCluster ),
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
            this->create_unique_dof_dv_and_field_type_maps();

            // create a dof and dv type lists
            this->create_dof_and_dv_type_lists();

            // create a dof and dv type maps
            this->create_dof_and_dv_type_maps();

            // create IQI map
            this->create_IQI_map();
        }

        //------------------------------------------------------------------------------

        Set::~Set()
        {
            // delete the equation object pointers
            for ( MSI::Equation_Object* tEquationObj : mEquationObjList )
            {
                delete tEquationObj;
            }
            mEquationObjList.clear();

            // delete the field interpolator pointers
            this->delete_pointers();
        }

        //------------------------------------------------------------------------------

        void
        Set::initialize_set(
                const bool                 aIsStaggered,
                const Time_Continuity_Flag aTimeContinuityOnlyFlag )
        {
            if ( !mIsEmptySet )    // FIXME this flag is a hack. find better solution
            {
                // std::cout << "Initialize FEM Set with Name: " << mMeshSet->get_set_name() << std::endl;

                mIsStaggered = aIsStaggered;

                this->create_residual_dof_assembly_map();

                this->create_dof_assembly_map( aIsStaggered );

                this->create_mat_pdv_assembly_map();

                this->create_requested_IWG_list( aTimeContinuityOnlyFlag );

                this->create_requested_IQI_list();

                this->create_requested_IQI_type_map();

                this->build_requested_IWG_dof_type_list( aIsStaggered );

                this->build_requested_IQI_dof_type_list();

                // set fem set pointer to IWGs FIXME still needed done in constructor?
                for ( const std::shared_ptr< IWG >& tIWG : mRequestedIWGs )
                {
                    tIWG->set_set_pointer( this );
                }

                // set fem set pointer to IQIs FIXME still needed done in constructor?
                for ( const std::shared_ptr< IQI >& tIQI : mRequestedIQIs )
                {
                    tIQI->set_set_pointer( this );
                }
            }
        }

        //------------------------------------------------------------------------------

        void
        Set::finalize( MSI::Model_Solver_Interface* aModelSolverInterface )
        {
            if ( !mIsEmptySet )    // FIXME this flag is a hack. find better solution
            {
                // delete the field interpolator pointers
                this->delete_pointers();

                // create integration information
                this->create_integrator( aModelSolverInterface );

                // create the field interpolators
                this->create_field_interpolator_managers( aModelSolverInterface );

                // set field interpolator managers for the IWGs
                this->set_IWG_field_interpolator_managers();

                // set field interpolator managers for the IQIs
                this->set_IQI_field_interpolator_managers();
            }
        }

        //------------------------------------------------------------------------------

        void
        Set::free_memory()
        {
            for ( const std::shared_ptr< IWG >& tIWG : mIWGs )
            {
                tIWG->free_memory();
            }
        }

        //------------------------------------------------------------------------------

        std::string
        Set::get_set_name()
        {
            return mMeshSet->get_set_name();
        }
        //------------------------------------------------------------------------------

        void
        Set::delete_pointers()
        {
            if ( mLeaderFIManager != nullptr )
            {
                delete mLeaderFIManager;
                mLeaderFIManager = nullptr;
            }
            if ( mFollowerFIManager != nullptr )
            {
                delete mFollowerFIManager;
                mFollowerFIManager = nullptr;
            }
            if ( mLeaderPreviousFIManager != nullptr )
            {
                delete mLeaderPreviousFIManager;
                mLeaderPreviousFIManager = nullptr;
            }
            if ( mLeaderEigenFIManager != nullptr )
            {
                delete mLeaderEigenFIManager;
                mLeaderEigenFIManager = nullptr;
            }
        }

        //------------------------------------------------------------------------------

        void
        Set::create_integrator( MSI::Model_Solver_Interface* aModelSolverInterface )
        {
            // get time levels from model solver interface
            const Matrix< DDUMat >& tTimeLevels = aModelSolverInterface->get_dof_manager()->get_time_levels();

            uint tMaxTimeLevels = tTimeLevels.max();

            // initialize time geometry type
            mtk::Geometry_Type tTimeGeometryType = mtk::Geometry_Type::UNDEFINED;

            // init time integration order
            mtk::Integration_Order tTimeIntegrationOrder = mtk::Integration_Order::UNDEFINED;

            // switch on maximum time level
            switch ( tMaxTimeLevels )
            {
                case 1:
                {
                    tTimeGeometryType     = mtk::Geometry_Type::LINE;
                    tTimeIntegrationOrder = mtk::Integration_Order::BAR_1;
                    break;
                }
                case 2:
                {
                    tTimeGeometryType     = mtk::Geometry_Type::LINE;
                    tTimeIntegrationOrder = mtk::Integration_Order::BAR_2;
                    break;
                }
                default:
                {
                    MORIS_ERROR( false, "Set::create_integrator - only 1 or 2 time levels handled so far." );
                }
            }

            // if a time sideset or boundary
            if ( mTimeContinuity || mTimeBoundary )
            {
                tTimeGeometryType     = mtk::Geometry_Type::POINT;
                tTimeIntegrationOrder = mtk::Integration_Order::POINT;
            }

            // create an integration rule
            mtk::Integration_Rule tIntegrationRule(
                    mIGGeometryType,
                    mtk::Integration_Type::GAUSS,
                    this->get_auto_integration_order(
                            mElementType,
                            mIGGeometryType,
                            mIPSpaceInterpolationOrder ),
                    tTimeGeometryType,
                    mtk::Integration_Type::GAUSS,
                    tTimeIntegrationOrder );

            // create an integrator
            mtk::Integrator tIntegrator( tIntegrationRule );

            // get integration points
            tIntegrator.get_points( mIntegPoints );

            // get integration weights
            tIntegrator.get_weights( mIntegWeights );
        }

        //------------------------------------------------------------------------------

        void
        Set::create_unique_dof_and_dv_type_lists()
        {
            // init dof and dv type counter
            uint tLeaderDofCounter     = 0;
            uint tLeaderDvCounter      = 0;
            uint tLeaderFieldCounter   = 0;
            uint tFollowerDofCounter   = 0;
            uint tFollowerDvCounter    = 0;
            uint tFollowerFieldCounter = 0;

            // loop over the IWGs
            for ( const std::shared_ptr< IWG >& tIWG : mIWGs )
            {
                // get an IWG non unique dof and dv types
                moris::Cell< moris::Cell< MSI::Dof_Type > >   tActiveDofType;
                moris::Cell< moris::Cell< PDV_Type > >        tActiveDvType;
                moris::Cell< moris::Cell< mtk::Field_Type > > tActiveFieldType;

                tIWG->get_non_unique_dof_dv_and_field_types( tActiveDofType, tActiveDvType, tActiveFieldType );

                // update dof and dv type counters
                tLeaderDofCounter += tActiveDofType( 0 ).size();
                tLeaderDvCounter += tActiveDvType( 0 ).size();
                tLeaderFieldCounter += tActiveFieldType( 0 ).size();
                tFollowerDofCounter += tActiveDofType( 1 ).size();
                tFollowerDvCounter += tActiveDvType( 1 ).size();
                tFollowerFieldCounter += tActiveFieldType( 1 ).size();
            }

            // loop over the IQIs
            for ( const std::shared_ptr< IQI >& tIQI : mIQIs )
            {
                // get an IWG non unique dof and dv types
                moris::Cell< moris::Cell< MSI::Dof_Type > >   tActiveDofType;
                moris::Cell< moris::Cell< PDV_Type > >        tActiveDvType;
                moris::Cell< moris::Cell< mtk::Field_Type > > tActiveFieldType;

                tIQI->get_non_unique_dof_dv_and_field_types( tActiveDofType, tActiveDvType, tActiveFieldType );

                // update dof and dv type counter
                tLeaderDofCounter += tActiveDofType( 0 ).size();
                tLeaderDvCounter += tActiveDvType( 0 ).size();
                tLeaderFieldCounter += tActiveFieldType( 0 ).size();
                tFollowerDofCounter += tActiveDofType( 1 ).size();
                tFollowerDvCounter += tActiveDvType( 1 ).size();
                tFollowerFieldCounter += tActiveFieldType( 1 ).size();
            }

            mUniqueDofTypeListLeaderFollower.resize( 2 );
            mUniqueDvTypeListLeaderFollower.resize( 2 );
            mUniqueFieldTypeListLeaderFollower.resize( 2 );

            mUniqueDofTypeListLeaderFollower( 0 ).reserve( tLeaderDofCounter );
            mUniqueDofTypeListLeaderFollower( 1 ).reserve( tFollowerDofCounter );
            mUniqueDvTypeListLeaderFollower( 0 ).reserve( tLeaderDvCounter );
            mUniqueDvTypeListLeaderFollower( 1 ).reserve( tFollowerDvCounter );
            mUniqueFieldTypeListLeaderFollower( 0 ).reserve( tLeaderFieldCounter );
            mUniqueFieldTypeListLeaderFollower( 1 ).reserve( tFollowerFieldCounter );

            // set max size for the unique dof and dv type lists
            mUniqueDofTypeList.reserve( tLeaderDofCounter + tFollowerDofCounter );
            mUniqueDvTypeList.reserve( tLeaderDvCounter + tFollowerDvCounter );
            mUniqueFieldTypeList.reserve( tLeaderFieldCounter + tFollowerFieldCounter );

            // loop over the IWGs
            for ( const std::shared_ptr< IWG >& tIWG : mIWGs )
            {
                // get non unique dof and dv types
                moris::Cell< moris::Cell< MSI::Dof_Type > >   tActiveDofType;
                moris::Cell< moris::Cell< PDV_Type > >        tActiveDvType;
                moris::Cell< moris::Cell< mtk::Field_Type > > tActiveFieldType;

                tIWG->get_non_unique_dof_dv_and_field_types( tActiveDofType, tActiveDvType, tActiveFieldType );

                // populate the corresponding unique dof and dv type lists
                mUniqueDofTypeListLeaderFollower( 0 ).append( tActiveDofType( 0 ) );
                mUniqueDofTypeListLeaderFollower( 1 ).append( tActiveDofType( 1 ) );
                mUniqueDvTypeListLeaderFollower( 0 ).append( tActiveDvType( 0 ) );
                mUniqueDvTypeListLeaderFollower( 1 ).append( tActiveDvType( 1 ) );
                mUniqueFieldTypeListLeaderFollower( 0 ).append( tActiveFieldType( 0 ) );
                mUniqueFieldTypeListLeaderFollower( 1 ).append( tActiveFieldType( 1 ) );

                mUniqueDofTypeList.append( tActiveDofType( 0 ) );
                mUniqueDofTypeList.append( tActiveDofType( 1 ) );
                mUniqueDvTypeList.append( tActiveDvType( 0 ) );
                mUniqueDvTypeList.append( tActiveDvType( 1 ) );
                mUniqueFieldTypeList.append( tActiveFieldType( 0 ) );
                mUniqueFieldTypeList.append( tActiveFieldType( 1 ) );
            }

            // loop over the IQIs
            for ( const std::shared_ptr< IQI >& tIQI : mIQIs )
            {
                // get non unique dof and dv types
                moris::Cell< moris::Cell< MSI::Dof_Type > >   tActiveDofType;
                moris::Cell< moris::Cell< PDV_Type > >        tActiveDvType;
                moris::Cell< moris::Cell< mtk::Field_Type > > tActiveFieldType;

                tIQI->get_non_unique_dof_dv_and_field_types( tActiveDofType, tActiveDvType, tActiveFieldType );

                // populate the corresponding unique dof and dv type lists
                mUniqueDofTypeListLeaderFollower( 0 ).append( tActiveDofType( 0 ) );
                mUniqueDofTypeListLeaderFollower( 1 ).append( tActiveDofType( 1 ) );
                mUniqueDvTypeListLeaderFollower( 0 ).append( tActiveDvType( 0 ) );
                mUniqueDvTypeListLeaderFollower( 1 ).append( tActiveDvType( 1 ) );
                mUniqueFieldTypeListLeaderFollower( 0 ).append( tActiveFieldType( 0 ) );
                mUniqueFieldTypeListLeaderFollower( 1 ).append( tActiveFieldType( 1 ) );

                mUniqueDofTypeList.append( tActiveDofType( 0 ) );
                mUniqueDofTypeList.append( tActiveDofType( 1 ) );
                mUniqueDvTypeList.append( tActiveDvType( 0 ) );
                mUniqueDvTypeList.append( tActiveDvType( 1 ) );
                mUniqueFieldTypeList.append( tActiveFieldType( 0 ) );
                mUniqueFieldTypeList.append( tActiveFieldType( 1 ) );
            }

            {
                // make the dof type list unique
                std::sort( ( mUniqueDofTypeListLeaderFollower( 0 ).data() ).data(),
                        ( mUniqueDofTypeListLeaderFollower( 0 ).data() ).data() + mUniqueDofTypeListLeaderFollower( 0 ).size() );
                auto last = std::unique( ( mUniqueDofTypeListLeaderFollower( 0 ).data() ).data(),
                        ( mUniqueDofTypeListLeaderFollower( 0 ).data() ).data() + mUniqueDofTypeListLeaderFollower( 0 ).size() );
                auto pos  = std::distance( ( mUniqueDofTypeListLeaderFollower( 0 ).data() ).data(), last );
                mUniqueDofTypeListLeaderFollower( 0 ).resize( pos );
            }

            {
                // make the dof type list unique
                std::sort( ( mUniqueDofTypeListLeaderFollower( 1 ).data() ).data(),
                        ( mUniqueDofTypeListLeaderFollower( 1 ).data() ).data() + mUniqueDofTypeListLeaderFollower( 1 ).size() );
                auto last = std::unique( ( mUniqueDofTypeListLeaderFollower( 1 ).data() ).data(),
                        ( mUniqueDofTypeListLeaderFollower( 1 ).data() ).data() + mUniqueDofTypeListLeaderFollower( 1 ).size() );
                auto pos  = std::distance( ( mUniqueDofTypeListLeaderFollower( 1 ).data() ).data(), last );
                mUniqueDofTypeListLeaderFollower( 1 ).resize( pos );
            }

            {
                // make the dv type list unique
                std::sort( ( mUniqueDvTypeListLeaderFollower( 0 ).data() ).data(),
                        ( mUniqueDvTypeListLeaderFollower( 0 ).data() ).data() + mUniqueDvTypeListLeaderFollower( 0 ).size() );
                auto last = std::unique( ( mUniqueDvTypeListLeaderFollower( 0 ).data() ).data(),
                        ( mUniqueDvTypeListLeaderFollower( 0 ).data() ).data() + mUniqueDvTypeListLeaderFollower( 0 ).size() );
                auto pos  = std::distance( ( mUniqueDvTypeListLeaderFollower( 0 ).data() ).data(), last );
                mUniqueDvTypeListLeaderFollower( 0 ).resize( pos );
            }

            {
                // make the dv type list unique
                std::sort( ( mUniqueDvTypeListLeaderFollower( 1 ).data() ).data(),
                        ( mUniqueDvTypeListLeaderFollower( 1 ).data() ).data() + mUniqueDvTypeListLeaderFollower( 1 ).size() );
                auto last = std::unique( ( mUniqueDvTypeListLeaderFollower( 1 ).data() ).data(),
                        ( mUniqueDvTypeListLeaderFollower( 1 ).data() ).data() + mUniqueDvTypeListLeaderFollower( 1 ).size() );
                auto pos  = std::distance( ( mUniqueDvTypeListLeaderFollower( 1 ).data() ).data(), last );
                mUniqueDvTypeListLeaderFollower( 1 ).resize( pos );
            }

            {
                // make the field type list unique
                std::sort( ( mUniqueFieldTypeListLeaderFollower( 0 ).data() ).data(),
                        ( mUniqueFieldTypeListLeaderFollower( 0 ).data() ).data() + mUniqueFieldTypeListLeaderFollower( 0 ).size() );
                auto last = std::unique( ( mUniqueFieldTypeListLeaderFollower( 0 ).data() ).data(),
                        ( mUniqueFieldTypeListLeaderFollower( 0 ).data() ).data() + mUniqueFieldTypeListLeaderFollower( 0 ).size() );
                auto pos  = std::distance( ( mUniqueFieldTypeListLeaderFollower( 0 ).data() ).data(), last );
                mUniqueFieldTypeListLeaderFollower( 0 ).resize( pos );
            }

            {
                // make the field type list unique
                std::sort( ( mUniqueFieldTypeListLeaderFollower( 1 ).data() ).data(),
                        ( mUniqueFieldTypeListLeaderFollower( 1 ).data() ).data() + mUniqueFieldTypeListLeaderFollower( 1 ).size() );
                auto last = std::unique( ( mUniqueFieldTypeListLeaderFollower( 1 ).data() ).data(),
                        ( mUniqueFieldTypeListLeaderFollower( 1 ).data() ).data() + mUniqueFieldTypeListLeaderFollower( 1 ).size() );
                auto pos  = std::distance( ( mUniqueFieldTypeListLeaderFollower( 1 ).data() ).data(), last );
                mUniqueFieldTypeListLeaderFollower( 1 ).resize( pos );
            }

            {
                // make the dof type list unique
                std::sort( ( mUniqueDofTypeList.data() ).data(),
                        ( mUniqueDofTypeList.data() ).data() + mUniqueDofTypeList.size() );
                auto last = std::unique( ( mUniqueDofTypeList.data() ).data(),
                        ( mUniqueDofTypeList.data() ).data() + mUniqueDofTypeList.size() );
                auto pos  = std::distance( ( mUniqueDofTypeList.data() ).data(), last );
                mUniqueDofTypeList.resize( pos );
            }

            {
                // make the dv type list unique
                std::sort( ( mUniqueDvTypeList.data() ).data(),
                        ( mUniqueDvTypeList.data() ).data() + mUniqueDvTypeList.size() );
                auto last = std::unique( ( mUniqueDvTypeList.data() ).data(),
                        ( mUniqueDvTypeList.data() ).data() + mUniqueDvTypeList.size() );
                auto pos  = std::distance( ( mUniqueDvTypeList.data() ).data(), last );
                mUniqueDvTypeList.resize( pos );
            }

            {
                // make the field type list unique
                std::sort( ( mUniqueFieldTypeList.data() ).data(),
                        ( mUniqueFieldTypeList.data() ).data() + mUniqueFieldTypeList.size() );
                auto last = std::unique( ( mUniqueFieldTypeList.data() ).data(),
                        ( mUniqueFieldTypeList.data() ).data() + mUniqueFieldTypeList.size() );
                auto pos  = std::distance( ( mUniqueFieldTypeList.data() ).data(), last );
                mUniqueFieldTypeList.resize( pos );
            }
        }

        //------------------------------------------------------------------------------

        void
        Set::create_dof_and_dv_type_lists()
        {
            // get number of dof and dv types
            uint tNumDofTypes   = this->get_num_unique_dof_types();
            uint tNumDvTypes    = this->get_num_unique_dv_types();
            uint tNumFieldTypes = this->get_num_unique_field_types();

            // set size for the global dof type list
            mLeaderDofTypes.reserve( tNumDofTypes );
            mFollowerDofTypes.reserve( tNumDofTypes );
            mLeaderDvTypes.reserve( tNumDvTypes );
            mFollowerDvTypes.reserve( tNumDvTypes );
            mLeaderFieldTypes.reserve( tNumFieldTypes );
            mFollowerFieldTypes.reserve( tNumFieldTypes );

            // create a list to check if dof type is already in the list
            Matrix< DDSMat > tLeaderCheckList( tNumDofTypes, 1, -1 );
            Matrix< DDSMat > tFollowerCheckList( tNumDofTypes, 1, -1 );
            Matrix< DDSMat > tLeaderDvCheckList( tNumDvTypes, 1, -1 );
            Matrix< DDSMat > tFollowerDvCheckList( tNumDvTypes, 1, -1 );
            Matrix< DDSMat > tLeaderFieldCheckList( tNumFieldTypes, 1, -1 );
            Matrix< DDSMat > tFollowerFieldCheckList( tNumFieldTypes, 1, -1 );

            // loop over the IWGs
            for ( const std::shared_ptr< IWG >& tIWG : mIWGs )
            {
                // get leader dof and dv types for the IWG
                const moris::Cell< moris::Cell< MSI::Dof_Type > >& tDofTypeLeader =
                        tIWG->get_global_dof_type_list();

                const moris::Cell< moris::Cell< PDV_Type > >& tDvTypeLeader =
                        tIWG->get_global_dv_type_list();

                const moris::Cell< moris::Cell< mtk::Field_Type > >& tFieldTypeLeader =
                        tIWG->get_global_field_type_list();

                // loop over the IWG active leader dof type
                for ( uint iDOF = 0; iDOF < tDofTypeLeader.size(); iDOF++ )
                {
                    // get set index for the treated leader dof type
                    sint tDofTypeindex = this->get_index_from_unique_dof_type_map( tDofTypeLeader( iDOF )( 0 ) );

                    // if dof enum not in the list
                    if ( tLeaderCheckList( tDofTypeindex ) != 1 )
                    {
                        // put the dof type in the checklist
                        tLeaderCheckList( tDofTypeindex ) = 1;

                        // put the dof type in the global type list
                        mLeaderDofTypes.push_back( tDofTypeLeader( iDOF ) );
                    }
                }

                // loop over the IWG active leader dv type
                for ( uint iDv = 0; iDv < tDvTypeLeader.size(); iDv++ )
                {
                    // get set index for the treated leader dof type
                    sint tDvTypeindex = this->get_index_from_unique_dv_type_map( tDvTypeLeader( iDv )( 0 ) );

                    // if dv enum not in the list
                    if ( tLeaderDvCheckList( tDvTypeindex ) != 1 )
                    {
                        // put the dof type in the checklist
                        tLeaderDvCheckList( tDvTypeindex ) = 1;

                        // put the dof type in the global type list
                        mLeaderDvTypes.push_back( tDvTypeLeader( iDv ) );
                    }
                }

                // loop over the IWG active leader field type
                for ( uint iFi = 0; iFi < tFieldTypeLeader.size(); iFi++ )
                {
                    // get set index for the treated leader dof type
                    sint tFieldTypeindex = this->get_index_from_unique_field_type_map( tFieldTypeLeader( iFi )( 0 ) );

                    // if dv enum not in the list
                    if ( tLeaderFieldCheckList( tFieldTypeindex ) != 1 )
                    {
                        // put the dof type in the checklist
                        tLeaderFieldCheckList( tFieldTypeindex ) = 1;

                        // put the dof type in the global type list
                        mLeaderFieldTypes.push_back( tFieldTypeLeader( iFi ) );
                    }
                }

                // get follower dof and dv types for the IWG
                const moris::Cell< moris::Cell< MSI::Dof_Type > >& tDofTypeFollower =
                        tIWG->get_global_dof_type_list( mtk::Leader_Follower::FOLLOWER );

                const moris::Cell< moris::Cell< PDV_Type > >& tDvTypeFollower =
                        tIWG->get_global_dv_type_list( mtk::Leader_Follower::FOLLOWER );

                const moris::Cell< moris::Cell< mtk::Field_Type > >& tFieldTypeFollower =
                        tIWG->get_global_field_type_list( mtk::Leader_Follower::FOLLOWER );

                // loop over the IWG active follower dof type
                for ( uint iDOF = 0; iDOF < tDofTypeFollower.size(); iDOF++ )
                {
                    // get set index for the treated follower dof type
                    sint tDofTypeindex = this->get_index_from_unique_dof_type_map( tDofTypeFollower( iDOF )( 0 ) );

                    // if dof enum not in the list
                    if ( tFollowerCheckList( tDofTypeindex ) != 1 )
                    {
                        // put the dof type in the check list
                        tFollowerCheckList( tDofTypeindex ) = 1;

                        // put the dof type in the global type list
                        mFollowerDofTypes.push_back( tDofTypeFollower( iDOF ) );
                    }
                }

                // loop over the IWG active follower dv type
                for ( uint iDv = 0; iDv < tDvTypeFollower.size(); iDv++ )
                {
                    // get set index for the treated follower dv type
                    sint tDvTypeindex = this->get_index_from_unique_dv_type_map( tDvTypeFollower( iDv )( 0 ) );

                    // if dv enum not in the list
                    if ( tFollowerDvCheckList( tDvTypeindex ) != 1 )
                    {
                        // put the dv type in the checklist
                        tFollowerDvCheckList( tDvTypeindex ) = 1;

                        // put the dv type in the global type list
                        mFollowerDvTypes.push_back( tDvTypeFollower( iDv ) );
                    }
                }

                // loop over the IWG active leader field type
                for ( uint iFi = 0; iFi < tFieldTypeFollower.size(); iFi++ )
                {
                    // get set index for the treated leader dof type
                    sint tFieldTypeindex = this->get_index_from_unique_field_type_map( tFieldTypeFollower( iFi )( 0 ) );

                    // if dv enum not in the list
                    if ( tFollowerFieldCheckList( tFieldTypeindex ) != 1 )
                    {
                        // put the dof type in the checklist
                        tFollowerFieldCheckList( tFieldTypeindex ) = 1;

                        // put the dof type in the global type list
                        mFollowerFieldTypes.push_back( tFieldTypeFollower( iFi ) );
                    }
                }
            }

            // loop over the IQIs
            for ( const std::shared_ptr< IQI >& tIQI : mIQIs )
            {
                // get leader dof and dv types for the IWG
                const moris::Cell< moris::Cell< MSI::Dof_Type > >& tDofTypeLeader =
                        tIQI->get_global_dof_type_list();

                const moris::Cell< moris::Cell< PDV_Type > >& tDvTypeLeader =
                        tIQI->get_global_dv_type_list();

                const moris::Cell< moris::Cell< mtk::Field_Type > >& tFieldTypeLeader =
                        tIQI->get_global_field_type_list();

                // loop over the IQI active leader dof type
                for ( uint iDOF = 0; iDOF < tDofTypeLeader.size(); iDOF++ )
                {
                    // get set index for the treated leader dof type
                    sint tDofTypeindex = this->get_index_from_unique_dof_type_map( tDofTypeLeader( iDOF )( 0 ) );

                    // if dof enum not in the list
                    if ( tLeaderCheckList( tDofTypeindex ) != 1 )
                    {
                        // put the dof type in the checklist
                        tLeaderCheckList( tDofTypeindex ) = 1;

                        // put the dof type in the global type list
                        mLeaderDofTypes.push_back( tDofTypeLeader( iDOF ) );
                    }
                }

                // loop over the IQI active leader dv type
                for ( uint iDv = 0; iDv < tDvTypeLeader.size(); iDv++ )
                {
                    // get set index for the treated leader dv type
                    sint tDvTypeindex = this->get_index_from_unique_dv_type_map( tDvTypeLeader( iDv )( 0 ) );

                    // if dv enum not in the list
                    if ( tLeaderDvCheckList( tDvTypeindex ) != 1 )
                    {
                        // put the dv type in the checklist
                        tLeaderDvCheckList( tDvTypeindex ) = 1;

                        // put the dv type in the global type list
                        mLeaderDvTypes.push_back( tDvTypeLeader( iDv ) );
                    }
                }

                // loop over the IWG active leader field type
                for ( uint iFi = 0; iFi < tFieldTypeLeader.size(); iFi++ )
                {
                    // get set index for the treated leader dof type
                    sint tFieldTypeindex = this->get_index_from_unique_field_type_map( tFieldTypeLeader( iFi )( 0 ) );

                    // if dv enum not in the list
                    if ( tLeaderFieldCheckList( tFieldTypeindex ) != 1 )
                    {
                        // put the dof type in the checklist
                        tLeaderFieldCheckList( tFieldTypeindex ) = 1;

                        // put the dof type in the global type list
                        mLeaderFieldTypes.push_back( tFieldTypeLeader( iFi ) );
                    }
                }

                // get follower dof and dv types for the IWG
                const moris::Cell< moris::Cell< MSI::Dof_Type > >& tDofTypeFollower =
                        tIQI->get_global_dof_type_list( mtk::Leader_Follower::FOLLOWER );
                const moris::Cell< moris::Cell< PDV_Type > >& tDvTypeFollower =
                        tIQI->get_global_dv_type_list( mtk::Leader_Follower::FOLLOWER );
                const moris::Cell< moris::Cell< mtk::Field_Type > >& tFieldTypeFollower =
                        tIQI->get_global_field_type_list( mtk::Leader_Follower::FOLLOWER );

                // loop over the IWG active follower dof type
                for ( uint iDOF = 0; iDOF < tDofTypeFollower.size(); iDOF++ )
                {
                    // get set index for the treated follower dof type
                    sint tDofTypeindex = this->get_index_from_unique_dof_type_map( tDofTypeFollower( iDOF )( 0 ) );

                    // if dof enum not in the list
                    if ( tFollowerCheckList( tDofTypeindex ) != 1 )
                    {
                        // put the dof type in the checklist
                        tFollowerCheckList( tDofTypeindex ) = 1;

                        // put the dof type in the global type list
                        mFollowerDofTypes.push_back( tDofTypeFollower( iDOF ) );
                    }
                }

                // loop over the IWG active follower dv type
                for ( uint iDv = 0; iDv < tDvTypeFollower.size(); iDv++ )
                {
                    // get set index for the treated follower dv type
                    sint tDvTypeindex = this->get_index_from_unique_dv_type_map( tDvTypeFollower( iDv )( 0 ) );

                    // if dv enum not in the list
                    if ( tFollowerDvCheckList( tDvTypeindex ) != 1 )
                    {
                        // put the dv type in the checklist
                        tFollowerDvCheckList( tDvTypeindex ) = 1;

                        // put the dv type in the global type list
                        mFollowerDvTypes.push_back( tDvTypeFollower( iDv ) );
                    }
                }

                // loop over the IWG active leader field type
                for ( uint iFi = 0; iFi < tFieldTypeFollower.size(); iFi++ )
                {
                    // get set index for the treated leader dof type
                    sint tFieldTypeindex = this->get_index_from_unique_field_type_map( tFieldTypeFollower( iFi )( 0 ) );

                    // if dv enum not in the list
                    if ( tFollowerFieldCheckList( tFieldTypeindex ) != 1 )
                    {
                        // put the dof type in the checklist
                        tFollowerFieldCheckList( tFieldTypeindex ) = 1;

                        // put the dof type in the global type list
                        mFollowerFieldTypes.push_back( tFieldTypeFollower( iFi ) );
                    }
                }
            }

            // shrink list to fit to number of unique dof and dv types
            mLeaderDofTypes.shrink_to_fit();
            mFollowerDofTypes.shrink_to_fit();
            mLeaderDvTypes.shrink_to_fit();
            mFollowerDvTypes.shrink_to_fit();
            mLeaderFieldTypes.shrink_to_fit();
            mFollowerFieldTypes.shrink_to_fit();
        }

        //------------------------------------------------------------------------------

        void
        Set::create_unique_dof_dv_and_field_type_maps()
        {
            // dof types
            //------------------------------------------------------------------------------
            // Create temporary dof type list
            const moris::Cell< enum MSI::Dof_Type >& tDofType = get_unique_dof_type_list();

            // Get number of unique adofs of this equation object
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
            const moris::Cell< enum PDV_Type >& tDvType = get_unique_dv_type_list();

            // Get number of unique dvs of this equation object
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

            // field types
            //------------------------------------------------------------------------------
            // Create temporary field type list
            const moris::Cell< mtk::Field_Type >& tFieldType = get_unique_field_type_list();

            // Get number of unique dvs of this equation object
            moris::uint tNumUniqueFieldTypes = tFieldType.size();

            // Get maximal dv type enum number
            moris::sint tMaxFieldTypeEnumNumber = 0;

            // Loop over all dv types to get the highest enum index
            for ( moris::uint Ii = 0; Ii < tNumUniqueFieldTypes; Ii++ )
            {
                tMaxFieldTypeEnumNumber =
                        std::max( tMaxFieldTypeEnumNumber, static_cast< int >( tFieldType( Ii ) ) );
            }

            // +1 because c++ is 0 based
            tMaxFieldTypeEnumNumber++;

            // Set size of mapping matrix
            mUniqueFieldTypeMap.set_size( tMaxFieldTypeEnumNumber, 1, -1 );

            // Loop over all dv types to create the mapping matrix
            for ( moris::uint Ii = 0; Ii < tNumUniqueFieldTypes; Ii++ )
            {
                mUniqueFieldTypeMap( static_cast< int >( tFieldType( Ii ) ), 0 ) = Ii;
            }
        }

        //------------------------------------------------------------------------------

        void
        Set::create_dof_and_dv_type_maps()
        {
            // get number of leader dof types
            uint tLeaderNumDofs = this->get_dof_type_list().size();

            // get maximal dof type enum
            sint tMaxEnum = -1;

            // loop over the IWGs
            for ( uint iDOF = 0; iDOF < tLeaderNumDofs; iDOF++ )
            {
                for ( uint Ik = 0; Ik < mLeaderDofTypes( iDOF ).size(); Ik++ )
                {
                    // get the highest dof type enum
                    tMaxEnum = std::max( tMaxEnum, static_cast< int >( mLeaderDofTypes( iDOF )( Ik ) ) );
                }
            }

            // get number of follower dof types
            uint tFollowerNumDofs = this->get_dof_type_list( mtk::Leader_Follower::FOLLOWER ).size();

            // loop over the IWGs
            for ( uint iDOF = 0; iDOF < tFollowerNumDofs; iDOF++ )
            {
                for ( uint Ik = 0; Ik < mFollowerDofTypes( iDOF ).size(); Ik++ )
                {
                    // get the highest dof type enum
                    tMaxEnum = std::max( tMaxEnum, static_cast< int >( mFollowerDofTypes( iDOF )( Ik ) ) );
                }
            }
            // +1 since start at 0
            tMaxEnum++;

            MORIS_ASSERT( tMaxEnum != -1, "Set::create_dof_and_dv_type_maps(), no information to build dof type map" );

            tMaxEnum = static_cast< int >( MSI::Dof_Type::END_ENUM );

            // set size of dof type map    // FIXME replace with map
            mLeaderDofTypeMap.set_size( tMaxEnum, 1, -1 );

            // loop over dof types
            for ( uint iDOF = 0; iDOF < tLeaderNumDofs; iDOF++ )
            {
                mLeaderDofTypeMap( static_cast< int >( mLeaderDofTypes( iDOF )( 0 ) ), 0 ) = iDOF;
            }

            // set size of dof type map
            mFollowerDofTypeMap.set_size( tMaxEnum, 1, -1 );

            // loop over dof types
            for ( uint iDOF = 0; iDOF < tFollowerNumDofs; iDOF++ )
            {
                mFollowerDofTypeMap( static_cast< int >( mFollowerDofTypes( iDOF )( 0 ) ), 0 ) = iDOF;
            }

            // dv type maps -------------------------------------------------------
            // get number of leader dv types
            uint tLeaderNumDvs = this->get_dv_type_list().size();

            // get maximal dv type enum
            tMaxEnum = -1;

            // loop over the dv types
            for ( uint iDv = 0; iDv < tLeaderNumDvs; iDv++ )
            {
                for ( uint Ik = 0; Ik < mLeaderDvTypes( iDv ).size(); Ik++ )
                {
                    // get the highest dof type enum
                    tMaxEnum = std::max( tMaxEnum, static_cast< int >( mLeaderDvTypes( iDv )( Ik ) ) );
                }
            }

            // get number of follower dv types
            uint tFollowerNumDvs = this->get_dv_type_list( mtk::Leader_Follower::FOLLOWER ).size();

            // loop over the IWGs
            for ( uint iDOF = 0; iDOF < tFollowerNumDvs; iDOF++ )
            {
                for ( uint Ik = 0; Ik < mFollowerDvTypes( iDOF ).size(); Ik++ )
                {
                    // get the highest dof type enum
                    tMaxEnum = std::max( tMaxEnum, static_cast< int >( mFollowerDvTypes( iDOF )( Ik ) ) );
                }
            }
            // +1 since start at 0
            tMaxEnum++;

            MORIS_ASSERT( tMaxEnum != -1, "Set::create_dv_type_map(), no information to build dv type map" );

            // set size of dv type map    // FIXME replace with map
            mLeaderDvTypeMap.set_size( tMaxEnum, 1, -1 );

            // loop over dv types
            for ( uint iDv = 0; iDv < tLeaderNumDvs; iDv++ )
            {
                mLeaderDvTypeMap( static_cast< int >( mLeaderDvTypes( iDv )( 0 ) ), 0 ) = iDv;
            }

            // set size of dv type map
            mFollowerDvTypeMap.set_size( tMaxEnum, 1, -1 );

            // loop over dv types
            for ( uint iDv = 0; iDv < tFollowerNumDvs; iDv++ )
            {
                mFollowerDvTypeMap( static_cast< int >( mFollowerDvTypes( iDv )( 0 ) ), 0 ) = iDv;
            }

            // field type maps -------------------------------------------------------
            // get number of leader field types
            uint tLeaderNumFields = this->get_field_type_list().size();

            // get maximal field type enum
            tMaxEnum = -1;

            // loop over the field types
            for ( uint iFi = 0; iFi < tLeaderNumFields; iFi++ )
            {
                for ( uint Ik = 0; Ik < mLeaderFieldTypes( iFi ).size(); Ik++ )
                {
                    // get the highest field type enum
                    tMaxEnum = std::max( tMaxEnum, static_cast< int >( mLeaderFieldTypes( iFi )( Ik ) ) );
                }
            }

            // get number of follower dv types
            uint tFollowerNumFields = this->get_field_type_list( mtk::Leader_Follower::FOLLOWER ).size();

            // loop over the IWGs
            for ( uint iFi = 0; iFi < tFollowerNumFields; iFi++ )
            {
                for ( uint Ik = 0; Ik < mFollowerFieldTypes( iFi ).size(); Ik++ )
                {
                    // get the highest dof type enum
                    tMaxEnum = std::max( tMaxEnum, static_cast< int >( mFollowerFieldTypes( iFi )( Ik ) ) );
                }
            }
            // +1 since start at 0
            tMaxEnum++;

            MORIS_ASSERT( tMaxEnum != -1, "Set::create_field_type_map(), no information to build field type map" );

            // set size of field type map    // FIXME replace with map
            mLeaderFieldTypeMap.set_size( tMaxEnum, 1, -1 );

            // loop over field types
            for ( uint iFi = 0; iFi < tLeaderNumFields; iFi++ )
            {
                mLeaderFieldTypeMap( static_cast< int >( mLeaderFieldTypes( iFi )( 0 ) ), 0 ) = iFi;
            }

            // set size of field type map
            mFollowerFieldTypeMap.set_size( tMaxEnum, 1, -1 );

            // loop over field types
            for ( uint iFi = 0; iFi < tFollowerNumFields; iFi++ )
            {
                mFollowerFieldTypeMap( static_cast< int >( mFollowerFieldTypes( iFi )( 0 ) ), 0 ) = iFi;
            }
        }

        //-----------------------------------------------------------------------------

        void
        Set::create_field_interpolator_managers(
                MSI::Model_Solver_Interface* aModelSolverInterface )
        {
            // create the leader field interpolator manager
            mLeaderFIManager = new Field_Interpolator_Manager(
                    mLeaderDofTypes,
                    mLeaderDvTypes,
                    mLeaderFieldTypes,
                    this );

            // assign cell shape to the field interpolator manager
            mLeaderFIManager->set_IG_cell_shape( mMeshSet->get_IG_cell_shape() );
            mLeaderFIManager->set_IP_cell_shape( mMeshSet->get_IP_cell_shape() );

            // create the geometry interpolators on the leader FI manager
            mLeaderFIManager->create_geometry_interpolators();

            // create the field interpolators on the leader FI manager
            mLeaderFIManager->create_field_interpolators( aModelSolverInterface );

            // create the follower field interpolator manager
            mFollowerFIManager = new Field_Interpolator_Manager(
                    mFollowerDofTypes,
                    mFollowerDvTypes,
                    mFollowerFieldTypes,
                    this,
                    mtk::Leader_Follower::FOLLOWER );

            // assign cell shape to the field interpolator manager
            mFollowerFIManager->set_IG_cell_shape( mMeshSet->get_IG_cell_shape() );
            mFollowerFIManager->set_IP_cell_shape( mMeshSet->get_IP_cell_shape() );

            // create the geometry interpolators on the follower FI manager
            mFollowerFIManager->create_geometry_interpolators();

            // create the field interpolators on the follower FI manager
            mFollowerFIManager->create_field_interpolators( aModelSolverInterface );

            // if time sideset
            if ( mElementType == fem::Element_Type::TIME_SIDESET )
            {
                // create the leader field interpolator manager
                mLeaderPreviousFIManager = new Field_Interpolator_Manager(
                        mLeaderDofTypes,
                        mLeaderDvTypes,
                        mLeaderFieldTypes,
                        this );

                // assign cell shape to the field interpolator manager
                mLeaderPreviousFIManager->set_IG_cell_shape( mMeshSet->get_IG_cell_shape() );
                mLeaderPreviousFIManager->set_IP_cell_shape( mMeshSet->get_IP_cell_shape() );

                // create the geometry interpolators on the leader FI manager
                mLeaderPreviousFIManager->create_geometry_interpolators();

                // create the field interpolators on the leader FI manager
                mLeaderPreviousFIManager->create_field_interpolators( aModelSolverInterface );
            }

            // if eigen vectors
            mNumEigenVectors = aModelSolverInterface->get_num_eigen_vectors();

            if ( mNumEigenVectors > 0 )
            {
                // create the leader field interpolator manager
                mLeaderEigenFIManager = new Field_Interpolator_Manager(
                        mLeaderDofTypes,
                        mLeaderDvTypes,
                        mLeaderFieldTypes,
                        this );

                // assign cell shape to the field interpolator manager
                mLeaderEigenFIManager->set_IG_cell_shape( mMeshSet->get_IG_cell_shape() );
                mLeaderEigenFIManager->set_IP_cell_shape( mMeshSet->get_IP_cell_shape() );

                // create the geometry interpolators on the leader FI manager
                mLeaderEigenFIManager->create_geometry_interpolators();

                // create the field interpolators on the leader FI manager
                mLeaderEigenFIManager->create_field_interpolators( aModelSolverInterface, mNumEigenVectors );
            }
        }

        //------------------------------------------------------------------------------

        Field_Interpolator_Manager*
        Set::get_field_interpolator_manager(
                mtk::Leader_Follower aIsLeader )
        {
            switch ( aIsLeader )
            {
                case mtk::Leader_Follower::LEADER:
                    return mLeaderFIManager;

                case mtk::Leader_Follower::FOLLOWER:
                    return mFollowerFIManager;

                default:
                    MORIS_ERROR( false, "Set::get_field_interpolator_manager - can only be leader or follower." );
                    return mLeaderFIManager;
            }
        }

        //------------------------------------------------------------------------------

        Field_Interpolator_Manager*
        Set::get_field_interpolator_manager_previous_time(
                mtk::Leader_Follower aIsLeader )
        {
            switch ( aIsLeader )
            {
                case mtk::Leader_Follower::LEADER:
                    return mLeaderPreviousFIManager;

                default:
                    MORIS_ERROR( false, "Set::get_field_interpolator_manager - can only be leader." );
                    return mLeaderPreviousFIManager;
            }
        }

        //------------------------------------------------------------------------------

        Field_Interpolator_Manager*
        Set::get_field_interpolator_manager_eigen_vectors(
                mtk::Leader_Follower aIsLeader )
        {
            switch ( aIsLeader )
            {
                case mtk::Leader_Follower::LEADER:
                    return mLeaderEigenFIManager;

                default:
                    MORIS_ERROR( false, "Set::get_field_interpolator_manager - can only be leader." );
                    return mLeaderPreviousFIManager;
            }
        }

        //------------------------------------------------------------------------------

        void
        Set::set_IWG_field_interpolator_managers()
        {
            // loop over the IWGs
            for ( const std::shared_ptr< IWG >& tIWG : mIWGs )
            {
                // set the leader FI manager
                tIWG->set_field_interpolator_manager( mLeaderFIManager );

                // if double sideset, set follower
                if ( mElementType == fem::Element_Type::DOUBLE_SIDESET )
                {
                    // set IWG follower field interpolator manager
                    tIWG->set_field_interpolator_manager(
                            mFollowerFIManager,
                            mtk::Leader_Follower::FOLLOWER );
                }

                // if time sideset, set previous
                if ( mElementType == fem::Element_Type::TIME_SIDESET )
                {
                    // set IWG leader field interpolator manager for previous time step
                    tIWG->set_field_interpolator_manager_previous_time(
                            mLeaderPreviousFIManager,
                            mtk::Leader_Follower::LEADER );
                }
            }
        }

        //------------------------------------------------------------------------------

        void
        Set::set_IWG_cluster_for_stabilization_parameters( fem::Cluster* aCluster )
        {
            // loop over the IWGs
            for ( auto tIWG : mIWGs )
            {
                // set the fem cluster to IWG
                tIWG->set_cluster_pointer( aCluster );

                // get the SP from the IWG
                moris::Cell< std::shared_ptr< Stabilization_Parameter > >& tSPs =
                        tIWG->get_stabilization_parameters();

                // loop over the SP
                for ( const std::shared_ptr< Stabilization_Parameter >& tSP : tSPs )
                {
                    // check if SP is null
                    if ( tSP != nullptr )
                    {
                        // set the fem cluster
                        tSP->set_cluster( aCluster );
                    }
                }
            }
        }

        //------------------------------------------------------------------------------

        void
        Set::set_IQI_cluster_for_stabilization_parameters( fem::Cluster* aCluster )
        {
            // loop over the IQIs
            for ( auto tIQI : mIQIs )
            {
                // set the fem cluster to IQI
                tIQI->set_cluster_pointer( aCluster );

                // get the SP from the IQI
                moris::Cell< std::shared_ptr< Stabilization_Parameter > >& tSPs =
                        tIQI->get_stabilization_parameters();

                // loop over the SPs
                for ( const std::shared_ptr< Stabilization_Parameter >& tSP : tSPs )
                {
                    // check if SP is null
                    if ( tSP != nullptr )
                    {
                        // set the fem cluster
                        tSP->set_cluster( aCluster );
                    }
                }
            }
        }

        //------------------------------------------------------------------------------

        void
        Set::build_cluster_measure_tuples_and_map()
        {
            // init cluster measure counter
            uint tCMEACounter = 0;

            // loop over the IWGs
            for ( auto tIWG : mIWGs )
            {
                // get the SP from the IWG
                moris::Cell< std::shared_ptr< Stabilization_Parameter > >& tSPs =
                        tIWG->get_stabilization_parameters();

                // loop over the SP
                for ( const std::shared_ptr< Stabilization_Parameter >& tSP : tSPs )
                {
                    // check if SP is null
                    if ( tSP != nullptr )
                    {
                        // get list of cluster measure tuple from SP
                        moris::Cell< std::tuple<
                                fem::Measure_Type,
                                mtk::Primary_Void,
                                mtk::Leader_Follower > >
                                tClusterMEASPTuples =
                                        tSP->get_cluster_measure_tuple_list();

                        // add number of cluster measure tuple to counter
                        tCMEACounter += tClusterMEASPTuples.size();
                    }
                }
            }

            // loop over the IQIs
            for ( auto tIQI : mIQIs )
            {
                // get the SP from the IQI
                moris::Cell< std::shared_ptr< Stabilization_Parameter > >& tSPs =
                        tIQI->get_stabilization_parameters();

                // loop over the SPs
                for ( const std::shared_ptr< Stabilization_Parameter >& tSP : tSPs )
                {
                    // check if SP is null
                    if ( tSP != nullptr )
                    {
                        // get list of cluster measure tuple from SP
                        moris::Cell< std::tuple<
                                fem::Measure_Type,
                                mtk::Primary_Void,
                                mtk::Leader_Follower > >
                                tClusterMEASPTuples =
                                        tSP->get_cluster_measure_tuple_list();

                        // add number of cluster measure tuple to counter
                        tCMEACounter += tClusterMEASPTuples.size();
                    }
                }
            }
            // init size for cell of cluster measure tuples
            mClusterMEATuples.resize( tCMEACounter );

            // reset CMECounter
            tCMEACounter = 0;

            // loop over the IWGs
            for ( auto tIWG : mIWGs )
            {
                // get the SP from the IWG
                moris::Cell< std::shared_ptr< Stabilization_Parameter > >& tSPs =
                        tIWG->get_stabilization_parameters();

                // loop over the SP
                for ( const std::shared_ptr< Stabilization_Parameter >& tSP : tSPs )
                {
                    // check if SP is null
                    if ( tSP != nullptr )
                    {
                        // get list of cluster measure tuple from SP
                        moris::Cell< std::tuple<
                                fem::Measure_Type,
                                mtk::Primary_Void,
                                mtk::Leader_Follower > >
                                tClusterMEASPTuples =
                                        tSP->get_cluster_measure_tuple_list();

                        // loop over the cluster measure tuples from SP
                        for ( uint iCMEA = 0; iCMEA < tClusterMEASPTuples.size(); iCMEA++ )
                        {
                            // check if the cluster measure tuple already in map
                            if ( mClusterMEAMap.find( tClusterMEASPTuples( iCMEA ) ) == mClusterMEAMap.end() )
                            {
                                // add cluster measure tuple in map
                                mClusterMEAMap[ tClusterMEASPTuples( iCMEA ) ] = tCMEACounter;

                                // add cluster measure tuple to list
                                mClusterMEATuples( tCMEACounter ) = tClusterMEASPTuples( iCMEA );

                                // update counter
                                tCMEACounter++;
                            }
                        }
                    }
                }
            }

            // loop over the IQIs
            for ( auto tIQI : mIQIs )
            {
                // get the SP from the IQI
                moris::Cell< std::shared_ptr< Stabilization_Parameter > >& tSPs =
                        tIQI->get_stabilization_parameters();

                // loop over the SPs
                for ( const std::shared_ptr< Stabilization_Parameter >& tSP : tSPs )
                {
                    // check if SP is null
                    if ( tSP != nullptr )
                    {
                        // get list of cluster measure tuple from SP
                        moris::Cell< std::tuple<
                                fem::Measure_Type,
                                mtk::Primary_Void,
                                mtk::Leader_Follower > >
                                tClusterMEASPTuples =
                                        tSP->get_cluster_measure_tuple_list();

                        // loop over the cluster measure tuples from SP
                        for ( uint iCMEA = 0; iCMEA < tClusterMEASPTuples.size(); iCMEA++ )
                        {
                            // check if the cluster measure tuple already in map
                            if ( mClusterMEAMap.find( tClusterMEASPTuples( iCMEA ) ) == mClusterMEAMap.end() )
                            {
                                // add cluster measure tuple in map
                                mClusterMEAMap[ tClusterMEASPTuples( iCMEA ) ] = tCMEACounter;

                                // add cluster measure tuple to list
                                mClusterMEATuples( tCMEACounter ) = tClusterMEASPTuples( iCMEA );

                                // update counter
                                tCMEACounter++;
                            }
                        }
                    }
                }
            }
            // shrink to fit as list of tuples shorter than what was reserved
            mClusterMEATuples.resize( tCMEACounter );

            // set build bool to true
            mBuildClusterMEA = true;
        }

        //------------------------------------------------------------------------------

        moris::Cell< std::tuple<
                fem::Measure_Type,
                mtk::Primary_Void,
                mtk::Leader_Follower > >&
        Set::get_cluster_measure_tuples()
        {
            // check that tuples were built
            if ( mBuildClusterMEA == false )
            {
                // if not build list
                this->build_cluster_measure_tuples_and_map();
            }

            return mClusterMEATuples;
        }

        //------------------------------------------------------------------------------

        std::map< std::tuple<
                          fem::Measure_Type,
                          mtk::Primary_Void,
                          mtk::Leader_Follower >,
                uint >&
        Set::get_cluster_measure_map()
        {
            // check that map was built
            if ( mBuildClusterMEA == false )
            {
                // if not build map
                this->build_cluster_measure_tuples_and_map();
            }

            return mClusterMEAMap;
        }

        //------------------------------------------------------------------------------

        void
        Set::set_IQI_field_interpolator_managers()
        {
            // loop over the IQIs
            for ( const std::shared_ptr< IQI >& tIQI : mIQIs )
            {
                // set IQI leader FI manager
                tIQI->set_field_interpolator_manager( mLeaderFIManager );

                // if double sideset, set follower
                if ( mElementType == fem::Element_Type::DOUBLE_SIDESET )
                {
                    // set IQI follower FI manager
                    tIQI->set_field_interpolator_manager(
                            mFollowerFIManager,
                            mtk::Leader_Follower::FOLLOWER );
                }

                // set IQI leader FI manager for eigen vectors
                if ( mLeaderEigenFIManager )
                {
                    tIQI->set_field_interpolator_manager_eigen_vector( mLeaderEigenFIManager );
                }
            }
        }

        //------------------------------------------------------------------------------

        void
        Set::create_residual_dof_assembly_map()
        {
            // get the list of requested dof types by the solver
            const moris::Cell< enum MSI::Dof_Type >& tRequestedDofTypes =
                    this->get_requested_dof_types();

            // init the max index for dof types
            sint tMaxDofIndex = -1;

            // loop over the requested dof types
            for ( uint Ik = 0; Ik < tRequestedDofTypes.size(); Ik++ )
            {
                // get the set index for the requested leader dof type
                sint tDofIndex = this->get_dof_index_for_type(
                        tRequestedDofTypes( Ik ),
                        mtk::Leader_Follower::LEADER );

                // if the index was set (and is different from -1)
                if ( tDofIndex != -1 )
                {
                    // update the max index for dof type
                    tMaxDofIndex = std::max( tMaxDofIndex, tDofIndex );
                }

                // get the set index for the requested follower dof type
                tDofIndex = this->get_dof_index_for_type(
                        tRequestedDofTypes( Ik ),
                        mtk::Leader_Follower::FOLLOWER );

                // if the index was set (and is different -1)
                if ( tDofIndex != -1 )
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
            for ( uint Ik = 0; Ik < mResDofAssemblyMap.size(); Ik++ )
            {
                mResDofAssemblyMap( Ik ).set_size( 1, 2, -1 );
            }

            // initialize dof coefficients counter
            uint tCounter = 0;

            // loop over the requested dof types
            for ( uint Ik = 0; Ik < tRequestedDofTypes.size(); Ik++ )
            {
                // get the set index for the requested leader dof type
                sint tDofIndex = this->get_dof_index_for_type(
                        tRequestedDofTypes( Ik ),
                        mtk::Leader_Follower::LEADER );

                // if the index was set (and is different from -1)
                if ( tDofIndex != -1 )
                {
                    // get the number of coefficients related to the leader dof type
                    uint tNumCoeff = mLeaderFIManager->get_field_interpolators_for_type( tRequestedDofTypes( Ik ) )->get_number_of_space_time_coefficients();

                    // fill the residual assembly map with starting and ending indices for the leader dof type
                    mResDofAssemblyMap( tDofIndex )( 0 ) = tCounter;
                    mResDofAssemblyMap( tDofIndex )( 1 ) = tCounter + tNumCoeff - 1;

                    // update the dof coefficient counter
                    tCounter += tNumCoeff;
                }
            }

            // loop over the requested dof types
            for ( uint Ik = 0; Ik < tRequestedDofTypes.size(); Ik++ )
            {
                // get the set index for the requested follower dof type
                sint tDofIndex = this->get_dof_index_for_type( tRequestedDofTypes( Ik ),
                        mtk::Leader_Follower::FOLLOWER );

                // if the dof type was set (its set index is different from -1)
                if ( tDofIndex != -1 )
                {
                    // get the number of coefficients for the follower dof type
                    uint tNumCoeff = mFollowerFIManager->get_field_interpolators_for_type( tRequestedDofTypes( Ik ) )->get_number_of_space_time_coefficients();

                    // fill the residual assembly map with starting and ending indices for the follower dof type
                    mResDofAssemblyMap( tDofIndex )( 0 ) = tCounter;
                    mResDofAssemblyMap( tDofIndex )( 1 ) = tCounter + tNumCoeff - 1;

                    // update the dof coefficient counter
                    tCounter += tNumCoeff;
                }
            }
        }

        //------------------------------------------------------------------------------

        void
        Set::create_dof_assembly_map( const bool aIsStaggered )
        {
            if ( aIsStaggered )
            {
                this->create_staggered_jacobian_dof_assembly_map();
            }
            else
            {
                this->create_jacobian_dof_assembly_map();
            }
        }

        //------------------------------------------------------------------------------

        void
        Set::create_jacobian_dof_assembly_map()
        {
            // get list of requested dof types (by the solver)
            const moris::Cell< enum MSI::Dof_Type >& tRequestedDofTypes =
                    this->get_requested_dof_types();

            sint tMaxDofIndex = -1;

            // leader
            for ( uint Ik = 0; Ik < tRequestedDofTypes.size(); Ik++ )
            {
                sint tDofIndex = this->get_dof_index_for_type(
                        tRequestedDofTypes( Ik ),
                        mtk::Leader_Follower::LEADER );

                if ( tDofIndex != -1 )
                {
                    tMaxDofIndex = std::max( tMaxDofIndex, tDofIndex );
                }

                tDofIndex = this->get_dof_index_for_type(
                        tRequestedDofTypes( Ik ),
                        mtk::Leader_Follower::FOLLOWER );

                if ( tDofIndex != -1 )
                {
                    tMaxDofIndex = std::max( tMaxDofIndex, tDofIndex );
                }
            }

            tMaxDofIndex++;

            mJacDofAssemblyMap.resize( tMaxDofIndex );

            for ( uint Ik = 0; Ik < mResDofAssemblyMap.size(); Ik++ )
            {
                mJacDofAssemblyMap( Ik ).set_size( tMaxDofIndex, 2, -1 );
            }

            for ( uint Ik = 0; Ik < tRequestedDofTypes.size(); Ik++ )
            {
                sint tDofIndex = this->get_dof_index_for_type(
                        tRequestedDofTypes( Ik ),
                        mtk::Leader_Follower::LEADER );

                if ( tDofIndex != -1 )
                {
                    uint tCounter_2 = 0;

                    for ( uint Ii = 0; Ii < tRequestedDofTypes.size(); Ii++ )
                    {
                        sint tDofIndex_2 = this->get_dof_index_for_type(
                                tRequestedDofTypes( Ii ),
                                mtk::Leader_Follower::LEADER );

                        if ( tDofIndex_2 != -1 )
                        {
                            uint tNumCoeff_2 = mLeaderFIManager->get_field_interpolators_for_type( tRequestedDofTypes( Ii ) )->get_number_of_space_time_coefficients();

                            mJacDofAssemblyMap( tDofIndex )( tDofIndex_2, 0 ) = tCounter_2;
                            mJacDofAssemblyMap( tDofIndex )( tDofIndex_2, 1 ) = tCounter_2 + tNumCoeff_2 - 1;

                            tCounter_2 += tNumCoeff_2;
                        }
                    }

                    // follower
                    for ( uint Ii = 0; Ii < tRequestedDofTypes.size(); Ii++ )
                    {
                        sint tDofIndex_2 = this->get_dof_index_for_type(
                                tRequestedDofTypes( Ii ),
                                mtk::Leader_Follower::FOLLOWER );

                        if ( tDofIndex_2 != -1 )
                        {
                            uint tNumCoeff_2 = mFollowerFIManager->get_field_interpolators_for_type( tRequestedDofTypes( Ii ) )->get_number_of_space_time_coefficients();

                            mJacDofAssemblyMap( tDofIndex )( tDofIndex_2, 0 ) = tCounter_2;
                            mJacDofAssemblyMap( tDofIndex )( tDofIndex_2, 1 ) = tCounter_2 + tNumCoeff_2 - 1;

                            tCounter_2 += tNumCoeff_2;
                        }
                    }
                }
            }

            for ( uint Ik = 0; Ik < tRequestedDofTypes.size(); Ik++ )
            {
                sint tDofIndex = this->get_dof_index_for_type(
                        tRequestedDofTypes( Ik ),
                        mtk::Leader_Follower::FOLLOWER );

                if ( tDofIndex != -1 )
                {
                    uint tCounter_2 = 0;

                    for ( uint Ii = 0; Ii < tRequestedDofTypes.size(); Ii++ )
                    {
                        sint tDofIndex_2 = this->get_dof_index_for_type(
                                tRequestedDofTypes( Ii ),
                                mtk::Leader_Follower::LEADER );

                        if ( tDofIndex_2 != -1 )
                        {
                            uint tNumCoeff_2 = mLeaderFIManager->get_field_interpolators_for_type( tRequestedDofTypes( Ii ) )->get_number_of_space_time_coefficients();

                            mJacDofAssemblyMap( tDofIndex )( tDofIndex_2, 0 ) = tCounter_2;
                            mJacDofAssemblyMap( tDofIndex )( tDofIndex_2, 1 ) = tCounter_2 + tNumCoeff_2 - 1;

                            tCounter_2 += tNumCoeff_2;
                        }
                    }

                    for ( uint Ii = 0; Ii < tRequestedDofTypes.size(); Ii++ )
                    {
                        sint tDofIndex_2 = this->get_dof_index_for_type(
                                tRequestedDofTypes( Ii ),
                                mtk::Leader_Follower::FOLLOWER );

                        if ( tDofIndex_2 != -1 )
                        {
                            uint tNumCoeff_2 = mFollowerFIManager->get_field_interpolators_for_type( tRequestedDofTypes( Ii ) )->get_number_of_space_time_coefficients();

                            mJacDofAssemblyMap( tDofIndex )( tDofIndex_2, 0 ) = tCounter_2;
                            mJacDofAssemblyMap( tDofIndex )( tDofIndex_2, 1 ) = tCounter_2 + tNumCoeff_2 - 1;

                            tCounter_2 += tNumCoeff_2;
                        }
                    }
                }
            }
        }

        //------------------------------------------------------------------------------

        void
        Set::create_staggered_jacobian_dof_assembly_map()
        {
            // get list of requested dof types
            const moris::Cell< enum MSI::Dof_Type >& tRequestedDofTypes =
                    this->get_requested_dof_types();

            const moris::Cell< enum MSI::Dof_Type >& tSecondaryDofTypes =
                    this->get_secondary_dof_types();

            sint tMaxDofIndex = -1;

            for ( uint Ik = 0; Ik < tRequestedDofTypes.size(); Ik++ )
            {
                sint tDofIndex = this->get_dof_index_for_type(
                        tRequestedDofTypes( Ik ),
                        mtk::Leader_Follower::LEADER );

                if ( tDofIndex != -1 )
                {
                    tMaxDofIndex = std::max( tMaxDofIndex, tDofIndex );
                }

                tDofIndex = this->get_dof_index_for_type(
                        tRequestedDofTypes( Ik ),
                        mtk::Leader_Follower::FOLLOWER );

                if ( tDofIndex != -1 )
                {
                    tMaxDofIndex = std::max( tMaxDofIndex, tDofIndex );
                }
            }
            tMaxDofIndex++;

            sint tMaxDofIndexSec = -1;

            for ( uint Ik = 0; Ik < tSecondaryDofTypes.size(); Ik++ )
            {
                sint tDofIndex = this->get_dof_index_for_type(
                        tSecondaryDofTypes( Ik ),
                        mtk::Leader_Follower::LEADER );

                if ( tDofIndex != -1 )
                {
                    tMaxDofIndexSec = std::max( tMaxDofIndexSec, tDofIndex );
                }

                tDofIndex = this->get_dof_index_for_type(
                        tSecondaryDofTypes( Ik ),
                        mtk::Leader_Follower::FOLLOWER );

                if ( tDofIndex != -1 )
                {
                    tMaxDofIndexSec = std::max( tMaxDofIndexSec, tDofIndex );
                }
            }
            tMaxDofIndexSec++;

            mJacDofAssemblyMap.resize( tMaxDofIndex );

            for ( uint Ik = 0; Ik < mResDofAssemblyMap.size(); Ik++ )
            {
                mJacDofAssemblyMap( Ik ).set_size( tMaxDofIndexSec, 2, -1 );
            }

            for ( uint Ik = 0; Ik < tRequestedDofTypes.size(); Ik++ )
            {
                sint tDofIndex = this->get_dof_index_for_type(
                        tRequestedDofTypes( Ik ),
                        mtk::Leader_Follower::LEADER );

                if ( tDofIndex != -1 )
                {
                    uint tCounter_2 = 0;

                    for ( uint Ii = 0; Ii < tSecondaryDofTypes.size(); Ii++ )
                    {
                        sint tDofIndex_2 = this->get_dof_index_for_type(
                                tSecondaryDofTypes( Ii ),
                                mtk::Leader_Follower::LEADER );

                        if ( tDofIndex_2 != -1 )
                        {
                            uint tNumCoeff_2 = mLeaderFIManager->get_field_interpolators_for_type( tSecondaryDofTypes( Ii ) )->get_number_of_space_time_coefficients();

                            mJacDofAssemblyMap( tDofIndex )( tDofIndex_2, 0 ) = tCounter_2;
                            mJacDofAssemblyMap( tDofIndex )( tDofIndex_2, 1 ) = tCounter_2 + tNumCoeff_2 - 1;

                            tCounter_2 += tNumCoeff_2;
                        }
                    }

                    for ( uint Ii = 0; Ii < tSecondaryDofTypes.size(); Ii++ )
                    {
                        sint tDofIndex_2 = this->get_dof_index_for_type(
                                tSecondaryDofTypes( Ii ),
                                mtk::Leader_Follower::FOLLOWER );

                        if ( tDofIndex_2 != -1 )
                        {
                            uint tNumCoeff_2 = mFollowerFIManager->get_field_interpolators_for_type( tSecondaryDofTypes( Ii ) )->get_number_of_space_time_coefficients();

                            mJacDofAssemblyMap( tDofIndex )( tDofIndex_2, 0 ) = tCounter_2;
                            mJacDofAssemblyMap( tDofIndex )( tDofIndex_2, 1 ) = tCounter_2 + tNumCoeff_2 - 1;

                            tCounter_2 += tNumCoeff_2;
                        }
                    }
                }
            }

            for ( uint Ik = 0; Ik < tRequestedDofTypes.size(); Ik++ )
            {
                sint tDofIndex = this->get_dof_index_for_type(
                        tRequestedDofTypes( Ik ),
                        mtk::Leader_Follower::FOLLOWER );

                if ( tDofIndex != -1 )
                {
                    uint tCounter_2 = 0;

                    for ( uint Ii = 0; Ii < tSecondaryDofTypes.size(); Ii++ )
                    {
                        sint tDofIndex_2 = this->get_dof_index_for_type(
                                tSecondaryDofTypes( Ii ),
                                mtk::Leader_Follower::LEADER );

                        if ( tDofIndex_2 != -1 )
                        {
                            uint tNumCoeff_2 = mLeaderFIManager->get_field_interpolators_for_type( tSecondaryDofTypes( Ii ) )
                                                       ->get_number_of_space_time_coefficients();

                            mJacDofAssemblyMap( tDofIndex )( tDofIndex_2, 0 ) = tCounter_2;
                            mJacDofAssemblyMap( tDofIndex )( tDofIndex_2, 1 ) = tCounter_2 + tNumCoeff_2 - 1;

                            tCounter_2 += tNumCoeff_2;
                        }
                    }

                    for ( uint Ii = 0; Ii < tSecondaryDofTypes.size(); Ii++ )
                    {
                        sint tDofIndex_2 = this->get_dof_index_for_type(
                                tSecondaryDofTypes( Ii ),
                                mtk::Leader_Follower::FOLLOWER );

                        if ( tDofIndex_2 != -1 )
                        {
                            uint tNumCoeff_2 = mFollowerFIManager->get_field_interpolators_for_type( tSecondaryDofTypes( Ii ) )->get_number_of_space_time_coefficients();

                            mJacDofAssemblyMap( tDofIndex )( tDofIndex_2, 0 ) = tCounter_2;
                            mJacDofAssemblyMap( tDofIndex )( tDofIndex_2, 1 ) = tCounter_2 + tNumCoeff_2 - 1;

                            tCounter_2 += tNumCoeff_2;
                        }
                    }
                }
            }
        }

        //------------------------------------------------------------------------------

        void
        Set::create_mat_pdv_assembly_map()
        {
            // get the list of requested dv types by the opt solver
            moris::Cell< enum PDV_Type > tRequestedDvTypes = mUniqueDvTypeList;

            // init the max index for dv types
            moris_index tMaxDvIndex = -1;

            // loop over the dv types
            for ( uint iPdvType = 0; iPdvType < tRequestedDvTypes.size(); iPdvType++ )
            {
                // get the current PDV type
                PDV_Type tPdvType = tRequestedDvTypes( iPdvType );

                // get the set index for the requested leader dof type
                moris_index tDvIndex = this->get_dv_index_for_type(
                        tPdvType,
                        mtk::Leader_Follower::LEADER );

                // if the index was set (and is different from -1)
                if ( tDvIndex != -1 )
                {
                    // update the max index for dv type
                    tMaxDvIndex = std::max( tMaxDvIndex, tDvIndex );
                }

                // get the set index for the requested follower follower type
                tDvIndex = this->get_dv_index_for_type(
                        tPdvType,
                        mtk::Leader_Follower::FOLLOWER );

                // if the index was set (and is different -1)
                if ( tDvIndex != -1 )
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
            for ( uint Ik = 0; Ik < mPdvMatAssemblyMap.size(); Ik++ )
            {
                mPdvMatAssemblyMap( Ik ).set_size( 1, 2, -1 );
            }

            // init dv coefficients counter
            uint tCounter = 0;

            // loop over the dv types
            for ( uint iPdvType = 0; iPdvType < tRequestedDvTypes.size(); iPdvType++ )
            {
                // get the current PDV type
                PDV_Type tPdvType = tRequestedDvTypes( iPdvType );

                // get the set index for the requested leader dv type
                moris_index tDvIndex = this->get_dv_index_for_type(
                        tPdvType,
                        mtk::Leader_Follower::LEADER );

                // if the index was set (and is different from -1)
                if ( tDvIndex != -1 )
                {
                    // get the FI related to the leader pdv type
                    Field_Interpolator* tFI = mLeaderFIManager->get_field_interpolators_for_type( tPdvType );

                    // get the number of coefficients related to the leader dv type
                    uint tNumCoeff = tFI->get_number_of_space_time_coefficients();

                    // fill the dv assembly map with starting and ending indices
                    // for the leader dv type
                    mPdvMatAssemblyMap( tDvIndex )( 0, 0 ) = tCounter;
                    mPdvMatAssemblyMap( tDvIndex )( 0, 1 ) = tCounter + tNumCoeff - 1;

                    // update the dv coefficient counter
                    tCounter += tNumCoeff;
                }
            }

            // loop over the follower dv types
            for ( uint iPdvType = 0; iPdvType < tRequestedDvTypes.size(); iPdvType++ )
            {
                // get the current PDV type
                PDV_Type tPdvType = tRequestedDvTypes( iPdvType );

                // get the set index for the follower dv type
                moris_index tDvIndex = this->get_dv_index_for_type(
                        tPdvType,
                        mtk::Leader_Follower::FOLLOWER );

                // if the dv type was set (its set index is different from -1)
                if ( tDvIndex != -1 )
                {
                    // get the FI related to the leader pdv type
                    Field_Interpolator* tFI = mFollowerFIManager->get_field_interpolators_for_type( tPdvType );

                    // get the number of coefficients for the follower dv type
                    uint tNumCoeff = tFI->get_number_of_space_time_coefficients();

                    // fill the residual assembly map with starting and ending indices for the follower dof type
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

        void
        Set::create_geo_pdv_assembly_map(
                std::shared_ptr< fem::Cluster > aFemCluster )
        {
            // make sure this is done on the FEM clusters and not on the VIS clusters
            MORIS_ASSERT( !aFemCluster->is_VIS_cluster(),
                    "FEM::Set::create_geo_pdv_assembly_map() - "
                    "Trying to set PDV assembly map on a VIS cluster. This shouldn't happen." );

            // get the design variable interface
            MSI::Design_Variable_Interface* tDVInterface =
                    mEquationModel->get_design_variable_interface();

            // get the geo dv types requested by the opt
            moris::Cell< enum PDV_Type > tRequestedDvTypes;
            moris_index                  tMeshSetIndex = mMeshSet->get_set_index();
            tDVInterface->get_ig_unique_dv_types_for_set(
                    tMeshSetIndex,
                    tRequestedDvTypes );

            // get node indices on cluster
            moris::Matrix< moris::IndexMat > tNodeIndicesOnCluster;
            aFemCluster->get_vertex_indices_in_cluster_for_sensitivity( tNodeIndicesOnCluster );

            // get access to the FEM model
            MSI::Equation_Model* tFemModel = this->get_equation_model();

            // clean up assembly vector
            mPdvGeoAssemblyVector.set_size( 0, 0 );

            // get the pdv active flags and ids from the FEM IG nodes
            Matrix< DDSMat > tIsActivePdv;
            Matrix< DDSMat > tPdvIds;
            tFemModel->get_integration_xyz_pdv_active_flags_and_ids(
                    tNodeIndicesOnCluster,
                    tRequestedDvTypes,
                    tIsActivePdv,
                    tPdvIds );

            // init active geo pdv counter
            uint tActiveGeoPdvCounter = 0;

            // set flag for active pdv
            mPdvGeoAssemblyFlag = ( sum( tIsActivePdv ) > 0 );

            // if there are no active pdvs, skip the rest
            if ( !mPdvGeoAssemblyFlag )
            {
                return;
            }

            // clean up assembly map
            mPdvGeoAssemblyMap.clear();

            // get number of pdv types
            uint tNumPdvTypes = tRequestedDvTypes.size();

            // get number of nodes on cluster
            uint tNumIGNodes = tNodeIndicesOnCluster.numel();

            // reset assembly indices on nodes
            tFemModel->reset_integration_xyz_pdv_assembly_indices( tNodeIndicesOnCluster );

            // clean up assembly vector
            Matrix< DDSMat > tPdvGeoAssemblyVectorTemp( tNumIGNodes * tNumPdvTypes, 1, -1 );

            // loop over the requested pdv types
            for ( uint iGeoPdv = 0; iGeoPdv < tNumPdvTypes; iGeoPdv++ )
            {
                // get treated geo pdv type
                PDV_Type tGeoPdvType = tRequestedDvTypes( iGeoPdv );

                // get treated geo pdv type index
                // moris_index tGeoPdvIndex = static_cast< uint >( tGeoPdvType );

                // loop over the ig nodes on cluster
                for ( uint iIGNode = 0; iIGNode < tNumIGNodes; iIGNode++ )
                {
                    // get treated node index
                    moris_index tNodeIndex = tNodeIndicesOnCluster( iIGNode );

                    // create key pair
                    std::pair< moris_index, PDV_Type > tKeyPair = std::make_pair( tNodeIndex, tGeoPdvType );

                    // if active and not set in the map
                    bool tPdvIsActive   = tIsActivePdv( iIGNode, iGeoPdv );
                    bool tPdvIsNotInMap = ( mPdvGeoAssemblyMap.find( tKeyPair ) == mPdvGeoAssemblyMap.end() );
                    if ( tPdvIsActive && tPdvIsNotInMap )
                    {
                        // fill the map
                        mPdvGeoAssemblyMap[ tKeyPair ] = tActiveGeoPdvCounter;

                        // set node local index on cluster
                        tFemModel->set_integration_xyz_pdv_assembly_index(
                                tNodeIndex,
                                tGeoPdvType,
                                tActiveGeoPdvCounter );

                        // fill the global assembly vector
                        moris_id tPdvId                                   = tPdvIds( iIGNode, iGeoPdv );
                        tPdvGeoAssemblyVectorTemp( tActiveGeoPdvCounter ) = tPdvId;

                        // update active geo pdv counter
                        tActiveGeoPdvCounter++;
                    }

                }    // end for: each IG vertex

            }        // end for: each PDV type

            if ( tActiveGeoPdvCounter > 0 )
            {
                // fill assembly vector
                mPdvGeoAssemblyVector.set_size( tActiveGeoPdvCounter, 1, -1 );
                mPdvGeoAssemblyVector = tPdvGeoAssemblyVectorTemp( { 0, tActiveGeoPdvCounter - 1 }, { 0, 0 } );
            }

        }    // end function: fem::Set::create_geo_pdv_assembly_map()

        //------------------------------------------------------------------------------

        void
        Set::create_requested_IWG_list( const Time_Continuity_Flag aTimeContinuityOnlyFlag )
        {
            // get list of requested dof types from solver
            moris::Cell< enum MSI::Dof_Type > tRequestedDofTypes =
                    this->get_requested_dof_types();

            // clear requested IWG list
            mRequestedIWGs.clear();

            // reserve max size for requested IWG list
            mRequestedIWGs.reserve( mIWGs.size() );

            if ( mEquationModel->get_is_forward_analysis() )
            {
                // create identifier list, marking which IWGs are active
                Matrix< DDBMat > tActiveIWGs( mIWGs.size(), 1, false );

                // loop over the requested dof types
                for ( MSI::Dof_Type tDofType : tRequestedDofTypes )
                {
                    // loop over the IWG in set IWG list
                    for ( uint iIWG = 0; iIWG < mIWGs.size(); iIWG++ )
                    {
                        // residual dof types for current IWG
                        const moris::Cell< moris::Cell< MSI::Dof_Type > >& tResDofType =
                                mIWGs( iIWG )->get_residual_dof_type();

                        // number of residual dof types
                        const uint tNumResDofTypes = tResDofType.size();

                        // loop over all residual dof types
                        for ( uint iType = 0; iType < tNumResDofTypes; ++iType )
                        {
                            // if the IWG residual dof type is requested
                            if ( ( tResDofType( iType )( 0 ) == tDofType ) and ( !tActiveIWGs( iIWG ) ) )
                            {
                                // check whether only time continuity IWG should be considered
                                if ( aTimeContinuityOnlyFlag == Time_Continuity_Flag::TIME_CONTINUITY_ONLY )
                                {
                                    if ( mIWGs( iIWG )->get_IWG_type() == IWG_Type::TIME_CONTINUITY_DOF )
                                    {
                                        mRequestedIWGs.push_back( mIWGs( iIWG ) );
                                    }
                                }
                                else if ( aTimeContinuityOnlyFlag == Time_Continuity_Flag::NO_TIME_CONTINUITY )
                                {
                                    if ( mIWGs( iIWG )->get_IWG_type() != IWG_Type::TIME_CONTINUITY_DOF )
                                    {
                                        mRequestedIWGs.push_back( mIWGs( iIWG ) );
                                    }
                                }
                                else if ( aTimeContinuityOnlyFlag == Time_Continuity_Flag::DEFAULT )
                                {
                                    // add the IWG to the requested IWG list
                                    mRequestedIWGs.push_back( mIWGs( iIWG ) );
                                }

                                // mark IWG as active
                                tActiveIWGs( iIWG ) = true;
                            }
                        }
                    }
                }
            }
            else
            {
                // create identifier list, marking which IWGs are active
                Matrix< DDBMat > tActiveIWGs( mIWGs.size(), 1, false );

                // loop over the requested dof types
                for ( MSI::Dof_Type tDofType : tRequestedDofTypes )
                {
                    // loop over the IWG in set IWG list
                    for ( uint iIWG = 0; iIWG < mIWGs.size(); iIWG++ )
                    {
                        // residual dof types for current IWG
                        const moris::Cell< moris::Cell< MSI::Dof_Type > >& tResDofType =
                                mIWGs( iIWG )->get_residual_dof_type();

                        // number of residual dof types
                        const uint tNumResDofTypes = tResDofType.size();

                        // loop over all residual dof types
                        for ( uint iType = 0; iType < tNumResDofTypes; ++iType )
                        {
                            // if the IWG residual dof type is requested
                            if ( ( tResDofType( iType )( 0 ) == tDofType ) and ( !tActiveIWGs( iIWG ) ) )
                            {
                                if ( mEquationModel->get_is_adjoint_off_diagonal_time_contribution() )
                                {
                                    if ( mIWGs( iIWG )->get_IWG_type() == moris::fem::IWG_Type::TIME_CONTINUITY_DOF )
                                    {
                                        // add the IWg to the requested IWG list
                                        mRequestedIWGs.push_back( mIWGs( iIWG ) );
                                    }
                                }
                                else
                                {
                                    // add the IWG to the requested IWG list
                                    mRequestedIWGs.push_back( mIWGs( iIWG ) );
                                }

                                // mark IWG as active
                                tActiveIWGs( iIWG ) = true;
                            }
                        }
                    }
                }
            }

            // reduce the size of requested IWG list to fit
            mRequestedIWGs.shrink_to_fit();
        }

        //------------------------------------------------------------------------------

        void
        Set::create_requested_IQI_list()
        {
            // clear requested IQI list
            mRequestedIQIs.clear();

            // Get names of potential requested IQIs
            const moris::Cell< std::string >& tRequestedIQINames =
                    mEquationModel->get_requested_IQI_names();

            // get number of potential requested IQIs
            uint tNumREquestedIQINames = tRequestedIQINames.size();

            // reserve memory
            mRequestedIQIs.reserve( tNumREquestedIQINames );

            // loop over requested IQI names
            for ( uint Ik = 0; Ik < tNumREquestedIQINames; Ik++ )
            {
                // check if this set has the requested IQI
                if ( mIQINameToIndexMap.key_exists( tRequestedIQINames( Ik ) ) )
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

        void
        Set::create_IQI_map()
        {
            // erase the content of the map
            mIQINameToIndexMap.clear();

            uint tCounter = 0;

            // loop over all IQIs and build a name to index map
            for ( const std::shared_ptr< IQI >& tIQI : mIQIs )
            {
                std::string tIQIName = tIQI->get_name();

                mIQINameToIndexMap[ tIQIName ] = tCounter++;
            }
        }

        //------------------------------------------------------------------------------

        void
        Set::build_requested_IWG_dof_type_list( const bool aIsStaggered )
        {
            for ( const std::shared_ptr< IWG >& tIWG : mRequestedIWGs )
            {
                tIWG->build_requested_dof_type_list( aIsStaggered );
            }
        }

        //------------------------------------------------------------------------------

        void
        Set::build_requested_IQI_dof_type_list()
        {
            for ( const std::shared_ptr< IQI >& tIQI : mRequestedIQIs )
            {
                tIQI->build_requested_dof_type_list();
            }
        }

        //------------------------------------------------------------------------------

        void
        Set::initialize_mJacobian()
        {
            // if residual not initialized before
            if ( !mJacobianExist )
            {
                // get the dof types requested by the solver
                const moris::Cell< enum MSI::Dof_Type >& tRequestedDofTypes =
                        this->get_requested_dof_types();

                // init dof coefficient counter
                uint tNumCols = 0;

                // loop over the requested dof types
                for ( uint Ik = 0; Ik < tRequestedDofTypes.size(); Ik++ )
                {
                    // get the set index for the leader dof type
                    sint tDofIndex = this->get_dof_index_for_type(
                            tRequestedDofTypes( Ik ),
                            mtk::Leader_Follower::LEADER );

                    // if this leader dof is active
                    if ( tDofIndex != -1 )
                    {
                        // update number of dof coefficients
                        tNumCols += mLeaderFIManager->get_field_interpolators_for_type( tRequestedDofTypes( Ik ) )->get_number_of_space_time_coefficients();
                    }

                    // get the set index for the follower dof type
                    tDofIndex = this->get_dof_index_for_type(
                            tRequestedDofTypes( Ik ),
                            mtk::Leader_Follower::FOLLOWER );

                    // if this follower dof is active
                    if ( tDofIndex != -1 )
                    {
                        // update number of dof coefficients
                        tNumCols += mFollowerFIManager->get_field_interpolators_for_type( tRequestedDofTypes( Ik ) )->get_number_of_space_time_coefficients();
                    }
                }

                // if for residual evaluation
                if ( !mIsStaggered )
                {
                    // set size for the jacobian matrix
                    mJacobian.set_size( tNumCols, tNumCols, 0.0 );
                }
                // if for Jacobian evaluation
                else
                {
                    // get the secondary dof types from the solver
                    moris::Cell< enum MSI::Dof_Type > tSecDofTypes =
                            this->get_secondary_dof_types();

                    // init dof coefficient counter for rows
                    uint tNumRows = 0;

                    // loop over the groups of secondary dof types
                    for ( auto tSecDofTypesI : tSecDofTypes )
                    {
                        // get the set index for the leader dof type
                        sint tDofIndex = this->get_dof_index_for_type(
                                tSecDofTypesI,
                                mtk::Leader_Follower::LEADER );

                        // if this leader dof is active
                        if ( tDofIndex != -1 )
                        {
                            // update number of dof coefficients
                            tNumRows += mLeaderFIManager->get_field_interpolators_for_type( tSecDofTypesI )->get_number_of_space_time_coefficients();
                        }

                        // get the set index for the follower dof type
                        tDofIndex = this->get_dof_index_for_type(
                                tSecDofTypesI,
                                mtk::Leader_Follower::FOLLOWER );

                        // if this follower dof is active
                        if ( tDofIndex != -1 )
                        {
                            // update number of dof coefficients
                            tNumRows += mFollowerFIManager->get_field_interpolators_for_type( tSecDofTypesI )->get_number_of_space_time_coefficients();
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

        void
        Set::initialize_mResidual()
        {
            // if residual not initialized before
            if ( !mResidualExist )
            {
                moris::Cell< enum MSI::Dof_Type > tRequestedDofTypes;

                if ( !mIsStaggered )
                {
                    // get the dof types requested by the solver
                    tRequestedDofTypes = this->get_requested_dof_types();
                }
                else
                {
                    // get the dof types. sec and requested dof types are flipped
                    tRequestedDofTypes = this->get_secondary_dof_types();
                }

                // init dof coefficient counter
                uint tNumCoeff = 0;

                // loop over the requested dof types
                for ( uint Ik = 0; Ik < tRequestedDofTypes.size(); Ik++ )
                {
                    // get the set index for the leader dof type
                    sint tDofIndex = this->get_dof_index_for_type(
                            tRequestedDofTypes( Ik ),
                            mtk::Leader_Follower::LEADER );

                    // if this leader dof is active
                    if ( tDofIndex != -1 )
                    {
                        // update number of dof coefficients
                        tNumCoeff += mLeaderFIManager->get_field_interpolators_for_type( tRequestedDofTypes( Ik ) )->get_number_of_space_time_coefficients();
                    }

                    // get the set index for the follower dof type
                    tDofIndex = this->get_dof_index_for_type(
                            tRequestedDofTypes( Ik ),
                            mtk::Leader_Follower::FOLLOWER );

                    // if this follower dof is active
                    if ( tDofIndex != -1 )
                    {
                        // update number of dof coefficients
                        tNumCoeff += mFollowerFIManager->get_field_interpolators_for_type( tRequestedDofTypes( Ik ) )->get_number_of_space_time_coefficients();
                    }


                    // get the number of rhs
                    uint tNumRHS = mEquationModel->get_num_rhs();

                    // set size for the list of dQIdu vectors
                    mResidual.resize( tNumRHS );

                    // loop over the dQIdu vectors
                    for ( auto& tRes : mResidual )
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
                for ( auto& tRes : mResidual )
                {
                    // fill the residual vector with zeros
                    tRes.fill( 0.0 );
                }
            }
        }

        //------------------------------------------------------------------------------

        void
        Set::initialize_mQI()
        {
            // if list of QI values not initialized before
            if ( !mQIExist )
            {
                uint tNumQI = mEquationModel->get_requested_IQI_names().size();

                // set size for the list of QI values
                mQI.resize( tNumQI );

                // counter to iterate through the IQIs
                uint iIQICounter = 0;

                // loop over mQIs to set the size
                for ( auto& tQI : mQI )
                {
                    // get number of rows and columns for the IQI from the equation model
                    uint tRowNum = mEquationModel->get_IQI_values()( iIQICounter ).n_rows();
                    uint tColNum = mEquationModel->get_IQI_values()( iIQICounter ).n_cols();

                    // set the size of the IQI
                    tQI.set_size( tRowNum, tColNum, 0.0 );

                    // increase the counter
                    iIQICounter++;
                }

                // set the QI initialization flag to true
                mQIExist = true;
            }
            // if list of QI values initialized before
            else
            {
                // loop over the QI values
                for ( auto& tQI : mQI )
                {
                    // fill the QI value vector with zero
                    tQI.fill( 0.0 );
                }
            }
        }

        //------------------------------------------------------------------------------

        void
        Set::initialize_mdQIdpMat()
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
                moris::Cell< moris::Cell< enum PDV_Type > > tRequestedDvTypes;
                this->get_ip_dv_types_for_set( tRequestedDvTypes );

                // init pdv coefficient counter
                uint tNumPdvCoefficients = 0;

                // loop over the requested dv types
                for ( uint Ik = 0; Ik < tRequestedDvTypes.size(); Ik++ )
                {
                    // get the set index for the leader dof type
                    sint tDvIndex = this->get_dv_index_for_type(
                            tRequestedDvTypes( Ik )( 0 ),
                            mtk::Leader_Follower::LEADER );

                    // if this leader dv is active
                    if ( tDvIndex != -1 )
                    {
                        // update number of dof coefficients
                        tNumPdvCoefficients += mLeaderFIManager->get_field_interpolators_for_type( tRequestedDvTypes( Ik )( 0 ) )->get_number_of_space_time_coefficients();
                    }
                }
                // get the requested pdv types for the follower side
                this->get_ip_dv_types_for_set( tRequestedDvTypes, mtk::Leader_Follower::FOLLOWER );

                // loop over the requested dv types
                for ( uint Ik = 0; Ik < tRequestedDvTypes.size(); Ik++ )
                {
                    // get the set index for the follower dv type
                    sint tDvIndex = this->get_dv_index_for_type(
                            tRequestedDvTypes( Ik )( 0 ),
                            mtk::Leader_Follower::FOLLOWER );

                    // if this follower dv is active
                    if ( tDvIndex != -1 )
                    {
                        // update number of dof coefficients
                        tNumPdvCoefficients += mFollowerFIManager->get_field_interpolators_for_type( tRequestedDvTypes( Ik )( 0 ) )->get_number_of_space_time_coefficients();
                    }
                }
                // loop over requested IQIs
                for ( auto& tdQIdp : mdQIdp( 0 ) )
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
                for ( auto& tdQIdp : mdQIdp( 0 ) )
                {
                    // fill the dQIdp vector with zero
                    tdQIdp.fill( 0.0 );
                }
            }
        }

        //------------------------------------------------------------------------------

        void
        Set::initialize_mdQIdpGeo( std::shared_ptr< fem::Cluster > aFemCluster )
        {
            // get the number of requested IQIs
            uint tNumRequestedIQIs = this->get_equation_model()->get_requested_IQI_names().size();

            // set size for dQIdp
            mdQIdp( 1 ).resize( tNumRequestedIQIs );

            // loop over requested IQIs
            for ( auto& tdQIdp : mdQIdp( 1 ) )
            {
                // fill the dQIdp vector with zero
                tdQIdp.set_size( 1, mPdvGeoAssemblyVector.numel(), 0.0 );
            }
        }

        //----------------------------------------------------------------------

        void
        Set::initialize_mdRdpMat()
        {
            // if dRdpMap not initialized before
            if ( !mdRdpMatExist )
            {
                // set size for dRdp
                mdRdp.resize( 2 );

                // get the dof types requested by the solver
                const moris::Cell< enum MSI::Dof_Type >& tRequestedDofTypes =
                        this->get_requested_dof_types();

                // init dof coefficient counter
                uint tNumRows = 0;

                // loop over the requested dof types
                for ( uint Ik = 0; Ik < tRequestedDofTypes.size(); Ik++ )
                {
                    // get the set index for the leader dof type
                    sint tDofIndex = this->get_dof_index_for_type(
                            tRequestedDofTypes( Ik ),
                            mtk::Leader_Follower::LEADER );

                    // if this leader dof is active
                    if ( tDofIndex != -1 )
                    {
                        // update number of dof coefficients
                        tNumRows += mLeaderFIManager->get_field_interpolators_for_type( tRequestedDofTypes( Ik ) )->get_number_of_space_time_coefficients();
                    }

                    // get the set index for the follower dof type
                    tDofIndex = this->get_dof_index_for_type(
                            tRequestedDofTypes( Ik ),
                            mtk::Leader_Follower::FOLLOWER );

                    // if this follower dof is active
                    if ( tDofIndex != -1 )
                    {
                        // update number of dof coefficients
                        tNumRows += mFollowerFIManager->get_field_interpolators_for_type( tRequestedDofTypes( Ik ) )->get_number_of_space_time_coefficients();
                    }
                }

                // get the dv types requested by the opt
                moris::Cell< moris::Cell< enum PDV_Type > > tRequestedDvTypes;
                this->get_ip_dv_types_for_set( tRequestedDvTypes );

                // init dv coefficient counter
                uint tNumCols = 0;

                // loop over the requested dv types
                for ( uint Ik = 0; Ik < tRequestedDvTypes.size(); Ik++ )
                {
                    // get the set index for the leader dof type
                    sint tDvIndex = this->get_dv_index_for_type(
                            tRequestedDvTypes( Ik )( 0 ),
                            mtk::Leader_Follower::LEADER );

                    // if this leader dv is active
                    if ( tDvIndex != -1 )
                    {
                        // update number of dof coefficients
                        tNumCols += mLeaderFIManager->get_field_interpolators_for_type( tRequestedDvTypes( Ik )( 0 ) )->get_number_of_space_time_coefficients();
                    }
                }

                // get the follower side dv types
                this->get_ip_dv_types_for_set( tRequestedDvTypes, mtk::Leader_Follower::FOLLOWER );

                // loop over the requested dv types
                for ( uint Ik = 0; Ik < tRequestedDvTypes.size(); Ik++ )
                {
                    // get the set index for the leader dof type
                    sint tDvIndex = this->get_dv_index_for_type(
                            tRequestedDvTypes( Ik )( 0 ),
                            mtk::Leader_Follower::FOLLOWER );

                    // if this follower dv is active
                    if ( tDvIndex != -1 )
                    {
                        // update number of dof coefficients
                        tNumCols += mFollowerFIManager->get_field_interpolators_for_type( tRequestedDvTypes( Ik )( 0 ) )->get_number_of_space_time_coefficients();
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

        void
        Set::initialize_mdRdpGeo( std::shared_ptr< fem::Cluster > aFemCluster )
        {
            // get the dof types requested by the solver
            const moris::Cell< enum MSI::Dof_Type >& tRequestedDofTypes =
                    this->get_requested_dof_types();

            // init dof coefficient counter
            uint tNumRows = 0;

            // loop over the requested dof types
            for ( uint Ik = 0; Ik < tRequestedDofTypes.size(); Ik++ )
            {
                // get the set index for the leader dof type
                sint tDofIndex = this->get_dof_index_for_type(
                        tRequestedDofTypes( Ik ),
                        mtk::Leader_Follower::LEADER );

                // if this leader dof is active
                if ( tDofIndex != -1 )
                {
                    // update number of dof coefficients
                    tNumRows += mLeaderFIManager->get_field_interpolators_for_type( tRequestedDofTypes( Ik ) )->get_number_of_space_time_coefficients();
                }

                // get the set index for the follower dof type
                tDofIndex = this->get_dof_index_for_type(
                        tRequestedDofTypes( Ik ),
                        mtk::Leader_Follower::FOLLOWER );

                // if this follower dof is active
                if ( tDofIndex != -1 )
                {
                    // update number of dof coefficients
                    tNumRows += mFollowerFIManager->get_field_interpolators_for_type( tRequestedDofTypes( Ik ) )->get_number_of_space_time_coefficients();
                }
            }

            // set size for dRdpgeo
            mdRdp( 1 ).set_size( tNumRows, mPdvGeoAssemblyVector.numel(), 0.0 );
        }

        //------------------------------------------------------------------------------

        mtk::Interpolation_Order
        Set::get_auto_interpolation_order(
                const moris::uint        aNumVertices,
                const mtk::Geometry_Type aGeometryType )
        {
            switch ( aGeometryType )
            {
                case mtk::Geometry_Type::LINE:
                {
                    switch ( aNumVertices )
                    {
                        case 1:
                            return mtk::Interpolation_Order::UNDEFINED;

                        case 2:
                            return mtk::Interpolation_Order::LINEAR;

                        case 3:
                            return mtk::Interpolation_Order::QUADRATIC;

                        default:
                            MORIS_ERROR( false, " Set::get_auto_interpolation_order - not defined for LINE and number of vertices. " );
                            return mtk::Interpolation_Order::UNDEFINED;
                    }
                }
                case mtk::Geometry_Type::QUAD:
                {
                    switch ( aNumVertices )
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
                            MORIS_ERROR( false, " Set::get_auto_interpolation_order - not defined for QUAD and number of vertices. " );
                            return mtk::Interpolation_Order::UNDEFINED;
                    }
                }
                case mtk::Geometry_Type::HEX:
                {
                    switch ( aNumVertices )
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
                            MORIS_ERROR( false, " Set::get_auto_interpolation_order - not defined for HEX and number of vertices. " );
                            return mtk::Interpolation_Order::UNDEFINED;
                    }
                }
                case mtk::Geometry_Type::TET:
                {
                    switch ( aNumVertices )
                    {
                        case 4:
                            return mtk::Interpolation_Order::LINEAR;

                        case 10:
                            return mtk::Interpolation_Order::QUADRATIC;

                        case 20:
                            return mtk::Interpolation_Order::CUBIC;

                        default:
                            MORIS_ERROR( false, " Set::get_auto_interpolation_order - not defined for TET and number of vertices. " );
                            return mtk::Interpolation_Order::UNDEFINED;
                    }
                }
                default:
                {
                    MORIS_ERROR( false, " Set::get_auto_interpolation_order - not defined for this geometry type. " );
                    return mtk::Interpolation_Order::UNDEFINED;
                }
            }
        }

        //------------------------------------------------------------------------------

        mtk::Interpolation_Type
        Set::get_auto_time_interpolation_type(
                const moris::uint aNumVertices )
        {
            switch ( aNumVertices )
            {
                case 1:
                    return mtk::Interpolation_Type::CONSTANT;

                case 2:
                case 3:
                case 4:
                    return mtk::Interpolation_Type::LAGRANGE;

                default:
                    MORIS_ERROR( false, " Element::get_auto_time_interpolation_type - not defined this number of time vertices. " );
                    return mtk::Interpolation_Type::UNDEFINED;
            }
        }

        //------------------------------------------------------------------------------

        mtk::Integration_Order
        Set::get_auto_integration_order(
                const fem::Element_Type        aSetType,
                const mtk::Geometry_Type       aGeometryType,
                const mtk::Interpolation_Order aInterpolationOrder )
        {
            switch ( aSetType )
            {
                case fem::Element_Type::BULK:
                case fem::Element_Type::TIME_SIDESET:
                case fem::Element_Type::TIME_BOUNDARY:
                {
                    switch ( aGeometryType )
                    {
                        case mtk::Geometry_Type::LINE:
                        {
                            switch ( aInterpolationOrder )
                            {
                                case mtk::Interpolation_Order::LINEAR:
                                    return mtk::Integration_Order::BAR_3;

                                case mtk::Interpolation_Order::QUADRATIC:
                                    return mtk::Integration_Order::BAR_4;

                                case mtk::Interpolation_Order::CUBIC:
                                    return mtk::Integration_Order::BAR_5;

                                default:
                                    MORIS_ERROR( false, "Set::get_auto_integration_order - Unknown or unsupported interpolation order." );
                                    return mtk::Integration_Order::UNDEFINED;
                            }
                            break;
                        }

                        case mtk::Geometry_Type::QUAD:
                        {
                            switch ( aInterpolationOrder )
                            {
                                case mtk::Interpolation_Order::LINEAR:
                                    return mtk::Integration_Order::QUAD_2x2;

                                case mtk::Interpolation_Order::SERENDIPITY:
                                    return mtk::Integration_Order::QUAD_3x3;

                                case mtk::Interpolation_Order::QUADRATIC:
                                    return mtk::Integration_Order::QUAD_3x3;

                                case mtk::Interpolation_Order::CUBIC:
                                    return mtk::Integration_Order::QUAD_5x5;

                                default:
                                    MORIS_ERROR( false, "Set::get_auto_integration_order - Unknown or unsupported interpolation order." );
                                    return mtk::Integration_Order::UNDEFINED;
                            }
                            break;
                        }

                        case mtk::Geometry_Type::HEX:
                        {
                            switch ( aInterpolationOrder )
                            {
                                case mtk::Interpolation_Order::LINEAR:
                                    return mtk::Integration_Order::HEX_2x2x2;

                                case mtk::Interpolation_Order::SERENDIPITY:
                                    return mtk::Integration_Order::HEX_3x3x3;

                                case mtk::Interpolation_Order::QUADRATIC:
                                    return mtk::Integration_Order::HEX_3x3x3;

                                case mtk::Interpolation_Order::CUBIC:
                                    return mtk::Integration_Order::HEX_4x4x4;

                                default:
                                    MORIS_ERROR( false, "Set::get_auto_integration_order - Unknown or unsupported interpolation order." );
                                    return mtk::Integration_Order::UNDEFINED;
                            }
                            break;
                        }

                        case mtk::Geometry_Type::TRI:
                        {
                            switch ( aInterpolationOrder )
                            {
                                case mtk::Interpolation_Order::LINEAR:
                                    return mtk::Integration_Order::TRI_7;

                                case mtk::Interpolation_Order::QUADRATIC:
                                    return mtk::Integration_Order::TRI_12;

                                case mtk::Interpolation_Order::CUBIC:
                                    return mtk::Integration_Order::TRI_25;

                                default:
                                    MORIS_ERROR( false, "Set::get_auto_integration_order - Unknown or unsupported interpolation order." );
                                    return mtk::Integration_Order::UNDEFINED;
                            }
                            break;
                        }

                        case mtk::Geometry_Type::TET:
                        {
                            switch ( aInterpolationOrder )
                            {
                                case mtk::Interpolation_Order::LINEAR:
                                    return mtk::Integration_Order::TET_11;

                                case mtk::Interpolation_Order::QUADRATIC:
                                    return mtk::Integration_Order::TET_35;

                                case mtk::Interpolation_Order::CUBIC:
                                    return mtk::Integration_Order::TET_56;

                                default:
                                    MORIS_ERROR( false, "Set::get_auto_integration_order - Unknown or unsupported interpolation order." );
                                    return mtk::Integration_Order::UNDEFINED;
                            }
                            break;
                        }

                        default:
                            MORIS_ERROR( false, " Set::get_auto_integration_order - Unknown or unsupported geometry type. " );
                            return mtk::Integration_Order::UNDEFINED;
                    }
                    break;
                }
                case fem::Element_Type::SIDESET:
                case fem::Element_Type::DOUBLE_SIDESET:
                {
                    switch ( aGeometryType )
                    {
                        case mtk::Geometry_Type::LINE:
                        {
                            switch ( aInterpolationOrder )
                            {
                                case mtk::Interpolation_Order::LINEAR:
                                    return mtk::Integration_Order::BAR_3;

                                case mtk::Interpolation_Order::SERENDIPITY:
                                    return mtk::Integration_Order::BAR_4;

                                case mtk::Interpolation_Order::QUADRATIC:
                                    return mtk::Integration_Order::BAR_4;

                                case mtk::Interpolation_Order::CUBIC:
                                    return mtk::Integration_Order::BAR_5;

                                default:
                                    MORIS_ERROR( false, "Set::get_auto_integration_order - Unknown or unsupported interpolation order." );
                                    return mtk::Integration_Order::UNDEFINED;
                            }
                            break;
                        }

                        case mtk::Geometry_Type::QUAD:
                        {
                            switch ( aInterpolationOrder )
                            {
                                case mtk::Interpolation_Order::LINEAR:
                                    return mtk::Integration_Order::QUAD_2x2;

                                case mtk::Interpolation_Order::SERENDIPITY:
                                    return mtk::Integration_Order::QUAD_4x4;

                                case mtk::Interpolation_Order::QUADRATIC:
                                    return mtk::Integration_Order::QUAD_4x4;

                                case mtk::Interpolation_Order::CUBIC:
                                    return mtk::Integration_Order::QUAD_5x5;

                                default:
                                    MORIS_ERROR( false, "Set::get_auto_integration_order - Unknown or unsupported interpolation order." );
                                    return mtk::Integration_Order::UNDEFINED;
                            }
                            break;
                        }

                        case mtk::Geometry_Type::TRI:
                        {
                            switch ( aInterpolationOrder )
                            {
                                case mtk::Interpolation_Order::LINEAR:
                                    return mtk::Integration_Order::TRI_7;

                                case mtk::Interpolation_Order::QUADRATIC:
                                    return mtk::Integration_Order::TRI_12;

                                case mtk::Interpolation_Order::CUBIC:
                                    return mtk::Integration_Order::TRI_25;

                                default:
                                    MORIS_ERROR( false, "Set::get_auto_integration_order - Unknown or unsupported interpolation order." );
                                    return mtk::Integration_Order::UNDEFINED;
                            }
                            break;
                        }

                        default:
                            MORIS_ERROR( false, " Set::get_auto_integration_order - Unknown or unsupported geometry type. " );
                            return mtk::Integration_Order::UNDEFINED;
                    }
                    break;
                }

                default:
                    MORIS_ERROR( false, "Set::get_auto_integration_order - unknown set type." );
                    return mtk::Integration_Order::UNDEFINED;
            }
        }

        //------------------------------------------------------------------------------

        void
        Set::set_visualization_set(
                const uint       aVisMeshIndex,
                moris::mtk::Set* aVisMeshSet,
                const bool       aOnlyPrimaryCells )
        {
            uint tNumClustersOnSets = aVisMeshSet->get_num_clusters_on_set();

            // the FEM mesh index corresponding to the current VIS mesh index is always +1
            moris_index tFemMeshIndex = aVisMeshIndex + 1;

            // set vis clusters to fem clusters
            for ( uint iCluster = 0; iCluster < tNumClustersOnSets; iCluster++ )
            {
                // get the VIS cluster
                const mtk::Cluster* tVisCluster = aVisMeshSet->get_clusters_by_index( iCluster );

                // get access to the corresponding FEM element
                MSI::Equation_Object* tFemElement = mEquationObjList( iCluster );

                // create a fem cluster
                std::shared_ptr< fem::Cluster > tCluster =
                        std::make_shared< fem::Cluster >(
                                mElementType,
                                tVisCluster,
                                this,
                                tFemElement,
                                true );    // mark this cluster to be a visualization cluster

                // reference the FEM cluster created from the VIS cluster on the FEM element for use in output
                reinterpret_cast< fem::Interpolation_Element* >( tFemElement )->set_cluster( tCluster, tFemMeshIndex );
            }

            // get the set type
            mtk::SetType tSetType = aVisMeshSet->get_set_type();

            // the below steps are only needed for block sets
            if ( tSetType == mtk::SetType::BULK )
            {
                this->construct_cell_assembly_map_for_VIS_set( aVisMeshIndex, aVisMeshSet, aOnlyPrimaryCells );
            }
            else if ( tSetType == mtk::SetType::SIDESET || tSetType == mtk::SetType::DOUBLE_SIDED_SIDESET )
            {
                this->construct_facet_assembly_map_for_VIS_set( aVisMeshIndex, aVisMeshSet );
            }
            else
            {
                MORIS_ERROR( false, "FEM::Set::set_visualization_set() - Unknown cluster type." );
            }


        }    // end function: Set::set_visualization_set()

        //------------------------------------------------------------------------------

        void
        Set::construct_cell_assembly_map_for_VIS_set(
                const uint       aVisMeshIndex,
                moris::mtk::Set* aVisMeshSet,
                const bool       aOnlyPrimaryCells )
        {
            // make sure this function is only called on block sets
            MORIS_ASSERT(
                    aVisMeshSet->get_set_type() == mtk::SetType::BULK,
                    "fem::Set::construct_cell_assembly_map_for_VIS_set() - Function can only be called for BULK sets." );

            // get the cells that are on the VIS cluster the fem cluster is supposed to output to
            Matrix< DDSMat > tCellIndices = aVisMeshSet->get_cell_inds_on_block( aOnlyPrimaryCells );
            uint             tNumCells    = aVisMeshSet->get_num_cells_on_set( aOnlyPrimaryCells );

            // resize arrays for current VIS mesh information, if array is too small
            if ( mCellAssemblyMap.size() < aVisMeshIndex + 1 )
            {
                mCellAssemblyMap.resize( aVisMeshIndex + 1 );
                mNumIgCellsOnSet.resize( aVisMeshIndex + 1 );
            }

            // store how many cells are in the current mesh set
            mNumIgCellsOnSet( aVisMeshIndex ) = tNumCells;

            // skip the rest for empty sets
            if ( tNumCells == 0 )
            {
                return;
            }

            // initialize the map relating mtk cell index in the mesh to the
            moris_index tMaxIndex = tCellIndices.max();
            mCellAssemblyMap( aVisMeshIndex ).set_size( tMaxIndex + 1, 1, -1 );

            // relate the mtk cell index to the position of the cell in the list of cells on a given block
            for ( uint iCell = 0; iCell < tNumCells; iCell++ )
            {
                moris_index tCellIndex                          = tCellIndices( iCell );
                mCellAssemblyMap( aVisMeshIndex )( tCellIndex ) = iCell;
            }

        }    // end function: fem::Set::construct_cell_assembly_map_for_VIS_set()

        //------------------------------------------------------------------------------

        void
        Set::construct_facet_assembly_map_for_VIS_set(
                const uint       aVisMeshIndex,
                moris::mtk::Set* aVisMeshSet )
        {
            // get the set type
            mtk::SetType tSetType = aVisMeshSet->get_set_type();

            // make sure this function is only called on (dbl) side sets
            MORIS_ASSERT(
                    tSetType == mtk::SetType::SIDESET || tSetType == mtk::SetType::DOUBLE_SIDED_SIDESET,
                    "fem::Set::construct_facet_assembly_map_for_VIS_set() - "
                    "Function can only be called for SIDESETs or DOUBLE_SIDED_SIDESETs." );

            // -------------
            // get initialization data for facet assembly map

            // get the clusters on the current set
            Cell< mtk::Cluster const * > const & tClustersOnSet = aVisMeshSet->get_clusters_on_set();

            // get the number of facets on the set's element type in the loop below
            uint tNumFacetsOnElementType = 0;

            // find the maximum (leader) cell index on the current set and count number of facets
            moris_index tMaxLeaderCellIndex = 0;
            uint        tNumFacetsOnSet     = 0;
            for ( uint iClusterOnSet = 0; iClusterOnSet < tClustersOnSet.size(); iClusterOnSet++ )
            {
                // get pointer for the current cluster
                mtk::Cluster const * tCluster = tClustersOnSet( iClusterOnSet );

                // get the Leader side
                mtk::Cluster const * tLeaderSideCluster;
                if ( tSetType == mtk::SetType::DOUBLE_SIDED_SIDESET )
                {
                    tLeaderSideCluster = &tCluster->get_leader_side_cluster();
                }
                else if ( tSetType == mtk::SetType::SIDESET )
                {
                    tLeaderSideCluster = tCluster;
                }
                else
                {
                    MORIS_ERROR( false,
                            "fem::Set::construct_facet_assembly_map_for_VIS_set() - "
                            "Function can only be called for SIDESETs or DOUBLE_SIDED_SIDESETs." );
                    tLeaderSideCluster = tCluster;
                }

                // get the info for a representative cell
                if ( iClusterOnSet == 0 )
                {
                    tNumFacetsOnElementType = tLeaderSideCluster->get_primary_cells_in_cluster()( 0 )->get_cell_info()->get_num_facets();
                }

                // get the IG cell indices on the leader cluster
                Matrix< IndexMat > tLeaderCellIndices = tLeaderSideCluster->get_primary_cell_indices_in_cluster();

                // get the maximum out of all the IG cell indices and count up number of facets
                for ( uint iCellOnCluster = 0; iCellOnCluster < tLeaderCellIndices.numel(); iCellOnCluster++ )
                {
                    moris_index tLeaderCellIndex = tLeaderCellIndices( iCellOnCluster );
                    tMaxLeaderCellIndex          = std::max( tLeaderCellIndex, tMaxLeaderCellIndex );
                    tNumFacetsOnSet++;
                }
            }

            // -------------
            // initialize facet assembly map

            // resize arrays for current VIS mesh information, if array is too small
            if ( mFacetAssemblyMap.size() < aVisMeshIndex + 1 )
            {
                mFacetAssemblyMap.resize( aVisMeshIndex + 1 );
                mNumFacetsOnSet.resize( aVisMeshIndex + 1 );
            }

            // store how many facets are in the current mesh set
            mNumFacetsOnSet( aVisMeshIndex ) = tNumFacetsOnSet;

            // skip the rest for empty sets
            if ( tClustersOnSet.size() == 0 )
            {
                return;
            }

            // set size of the map and initialize with default value
            mFacetAssemblyMap( aVisMeshIndex ).set_size( tMaxLeaderCellIndex + 1, tNumFacetsOnElementType, -1 );

            // -------------
            // construct facet assembly map

            // initialize index counter for facets
            moris_index tFacetIndexOnSet = 0;

            // go over clusters and facets within to construct map
            for ( uint iClusterOnSet = 0; iClusterOnSet < tClustersOnSet.size(); iClusterOnSet++ )
            {
                // get pointer for the current cluster
                mtk::Cluster const * tCluster = tClustersOnSet( iClusterOnSet );

                // get the Leader side
                mtk::Cluster const * tLeaderSideCluster;
                if ( tSetType == mtk::SetType::DOUBLE_SIDED_SIDESET )
                {
                    tLeaderSideCluster = &tCluster->get_leader_side_cluster();
                }
                else if ( tSetType == mtk::SetType::SIDESET )
                {
                    tLeaderSideCluster = tCluster;
                }
                else
                {
                    MORIS_ERROR( false,
                            "fem::Set::construct_facet_assembly_map_for_VIS_set() - "
                            "Function can only be called for SIDESETs or DOUBLE_SIDED_SIDESETs." );
                    tLeaderSideCluster = tCluster;
                }

                // get the IG cell indices on the leader cluster
                Matrix< IndexMat > tLeaderCellIndices = tLeaderSideCluster->get_primary_cell_indices_in_cluster();

                // get the corresponding side ordinals
                Matrix< IndexMat > tLeaderSideOrdinals = tLeaderSideCluster->get_cell_side_ordinals();

                // associate each facet in the cluster with its position in the list of facets in the set
                for ( uint iCellOnCluster = 0; iCellOnCluster < tLeaderCellIndices.numel(); iCellOnCluster++ )
                {
                    // get the info for the current facet
                    moris_index tLeaderCellIndex   = tLeaderCellIndices( iCellOnCluster );
                    moris_index tLeaderSideOrdinal = tLeaderSideOrdinals( iCellOnCluster );

                    // fill map
                    mFacetAssemblyMap( aVisMeshIndex )( tLeaderCellIndex, tLeaderSideOrdinal ) = tFacetIndexOnSet;

                    // increment facet index
                    tFacetIndexOnSet++;
                }

            }    // end for: clusters on (dbl) side set

        }        // end function: fem::Set::construct_facet_assembly_map_for_VIS_set()

        //------------------------------------------------------------------------------

        void
        Set::compute_quantity_of_interest_nodal(
                const uint                        aVisMeshIndex,
                Matrix< DDRMat >*                 aNodalFieldValues,
                const moris::Cell< std::string >& aQINames )
        {
            // FEM mesh index is VIS mesh index +1
            moris_index tFemMeshIndex = aVisMeshIndex + 1;

            // set the nodal set values to the ones provided
            mSetNodalValues = aNodalFieldValues;

            this->gather_requested_IQIs( aQINames, mRequestedNodalIQIs, mRequestedNodalIQIsGlobalIndices );

            // loop over equation objects
            uint tNumElements = mEquationObjList.size();
            for ( uint iElement = 0; iElement < tNumElements; iElement++ )
            {
                // compute quantity of interest
                mEquationObjList( iElement )->compute_quantity_of_interest( tFemMeshIndex, vis::Field_Type::NODAL );
            }
        }

        //------------------------------------------------------------------------------

        const moris::Cell< std::shared_ptr< fem::IQI > >&
        Set::get_requested_nodal_IQIs_for_visualization()
        {
            return mRequestedNodalIQIs;
        }

        //------------------------------------------------------------------------------

        uint
        Set::get_number_of_requested_nodal_IQIs_for_visualization()
        {
            return mRequestedNodalIQIs.size();
        }

        //------------------------------------------------------------------------------

        const moris::Cell< moris_index >&
        Set::get_requested_nodal_IQIs_global_indices_for_visualization()
        {
            return mRequestedNodalIQIsGlobalIndices;
        }

        //------------------------------------------------------------------------------

        void
        Set::compute_quantity_of_interest_global(
                const uint                        aVisMeshIndex,
                Matrix< DDRMat >*                 aGlobalFieldValues,
                const moris::Cell< std::string >& aQINames )
        {
            // FEM mesh index is VIS mesh index +1
            moris_index tFemMeshIndex = aVisMeshIndex + 1;

            // set the global set values to the ones provided
            mSetGlobalValues = aGlobalFieldValues;

            this->gather_requested_IQIs( aQINames, mRequestedGlobalIQIs, mRequestedGlobalIQIsGlobalIndices );

            // loop over equation objects
            uint tNumEqObjs = mEquationObjList.size();
            for ( uint iElement = 0; iElement < tNumEqObjs; iElement++ )
            {
                // compute quantity of interest
                mEquationObjList( iElement )->compute_quantity_of_interest( tFemMeshIndex, vis::Field_Type::GLOBAL );
            }
        }

        //------------------------------------------------------------------------------

        const moris::Cell< std::shared_ptr< fem::IQI > >&
        Set::get_requested_global_IQIs_for_visualization()
        {
            return mRequestedGlobalIQIs;
        }

        //------------------------------------------------------------------------------

        uint
        Set::get_number_of_requested_global_IQIs_for_visualization()
        {
            return mRequestedGlobalIQIs.size();
        }

        //------------------------------------------------------------------------------

        const moris::Cell< moris_index >&
        Set::get_requested_global_IQIs_global_indices_for_visualization()
        {
            return mRequestedGlobalIQIsGlobalIndices;
        }

        //------------------------------------------------------------------------------

        void
        Set::compute_quantity_of_interest_elemental(
                const uint                        aVisMeshIndex,
                Matrix< DDRMat >*                 aElementalFieldValues,
                const moris::Cell< std::string >& aQINames,
                const bool                        aOutputAverageValue )
        {
            // FEM mesh index is VIS mesh index +1
            moris_index tFemMeshIndex = aVisMeshIndex + 1;

            // set the elemental set values to the ones provided
            mSetElementalValues = aElementalFieldValues;

            // get the number of elements or facets on the set depending on the set type
            uint tNumElementsOnSet;
            if ( this->get_element_type() == fem::Element_Type::BULK )
            {
                tNumElementsOnSet = mNumIgCellsOnSet( aVisMeshIndex );
            }
            else    // SIDESET or DOUBLE_SIDESET
            {
                tNumElementsOnSet = mNumFacetsOnSet( aVisMeshIndex );
            }

            // initialize output vector
            // mSetElementalValues->set_size( tNumElementsOnSet, aQINames.size(), std::numeric_limits< real >::quiet_NaN() );
            mSetElementalValues->set_size( tNumElementsOnSet, aQINames.size(), 0.0 );

            this->gather_requested_IQIs( aQINames, mRequestedElementalIQIs, mRequestedElementalIQIsGlobalIndices );    // get the used IQIs on the current set

            // compute averages if requested by user
            vis::Field_Type tOutputType = vis::Field_Type::ELEMENTAL_INT;
            if ( aOutputAverageValue )
            {
                tOutputType = vis::Field_Type::ELEMENTAL_AVG;
            }

            // loop over equation objects
            uint tNumEqObjs = mEquationObjList.size();
            for ( uint iElement = 0; iElement < tNumEqObjs; iElement++ )
            {
                // compute quantity of interest
                mEquationObjList( iElement )->compute_quantity_of_interest( tFemMeshIndex, tOutputType );
            }
        }

        //------------------------------------------------------------------------------

        const moris::Cell< std::shared_ptr< fem::IQI > >&
        Set::get_requested_elemental_IQIs_for_visualization()
        {
            return mRequestedElementalIQIs;
        }

        //------------------------------------------------------------------------------

        uint
        Set::get_number_of_requested_elemental_IQIs_for_visualization()
        {
            return mRequestedElementalIQIs.size();
        }

        //------------------------------------------------------------------------------

        const moris::Cell< moris_index >&
        Set::get_requested_elemental_IQIs_global_indices_for_visualization()
        {
            return mRequestedElementalIQIsGlobalIndices;
        }

        //------------------------------------------------------------------------------

        const moris::Cell< std::shared_ptr< fem::IQI > >&
        Set::get_requested_field_IQIs()
        {
            return mRequestedFieldIQIs;
        }

        //------------------------------------------------------------------------------

        uint
        Set::get_number_of_requested_field_IQIs()
        {
            return mRequestedFieldIQIs.size();
        }

        //------------------------------------------------------------------------------

        const moris::Cell< moris_index >&
        Set::get_requested_field_IQIs_global_indices()
        {
            return mRequestedFieldIQIsGlobalIndices;
        }

        //------------------------------------------------------------------------------

        void
        Set::gather_requested_IQIs(
                moris::Cell< std::string > const &     aNames,
                moris::Cell< std::shared_ptr< IQI > >& aListOfRequestedIQIs,
                moris::Cell< moris_index >&            aListOfIQIGlobalIndices )
        {
            // get number of potential requested IQIs
            uint tNumRequestedElementalIQINames = aNames.size();

            // clear requested IQI list
            aListOfRequestedIQIs.clear();
            aListOfIQIGlobalIndices.clear();

            // reserve memory
            aListOfRequestedIQIs.reserve( tNumRequestedElementalIQINames );
            aListOfIQIGlobalIndices.reserve( tNumRequestedElementalIQINames );

            // loop over requested IQI names
            for ( uint Ik = 0; Ik < tNumRequestedElementalIQINames; Ik++ )
            {
                // check if this set has the requested IQI
                if ( mIQINameToIndexMap.key_exists( aNames( Ik ) ) )
                {
                    // get the set local index
                    moris_index tIQISetLocalIndex =
                            mIQINameToIndexMap.find( aNames( Ik ) );

                    // put global model index in list
                    aListOfIQIGlobalIndices.push_back( Ik );

                    // put IQI in requested IQI list
                    aListOfRequestedIQIs.push_back( mIQIs( tIQISetLocalIndex ) );
                }
            }

            // reduce memory to used space
            aListOfRequestedIQIs.shrink_to_fit();
            aListOfIQIGlobalIndices.shrink_to_fit();
        }

        //------------------------------------------------------------------------------

        void
        Set::determine_set_type()
        {
            mtk::SetType tMtkSetType = mMeshSet->get_set_type();

            switch ( tMtkSetType )
            {
                case mtk::SetType::BULK:
                {
                    mElementType = fem::Element_Type::BULK;

                    // if time continuity
                    if ( mTimeContinuity )
                    {
                        mElementType = fem::Element_Type::TIME_SIDESET;
                    }

                    // if time boundary
                    if ( mTimeBoundary )
                    {
                        mElementType = fem::Element_Type::TIME_BOUNDARY;
                    }
                    break;
                }

                case mtk::SetType::SIDESET:
                    mElementType = fem::Element_Type::SIDESET;
                    break;

                case mtk::SetType::DOUBLE_SIDED_SIDESET:
                    mElementType = fem::Element_Type::DOUBLE_SIDESET;
                    break;

                default:
                    MORIS_ERROR( false, "Set::determine_set_type() - not defined for this set type. " );
            }
        }

        //------------------------------------------------------------------------------

        void
        Set::set_set_type( fem::Element_Type aElementType )
        {
            mElementType = aElementType;
        }

        //------------------------------------------------------------------------------

        void
        Set::get_ig_unique_dv_types_for_set(
                moris::Cell< enum PDV_Type >& aGeoPdvType )
        {
            // get design variable interface
            MSI::Design_Variable_Interface* tPdvInterface =
                    mEquationModel->get_design_variable_interface();

            // if the pdv interface is set
            if ( tPdvInterface )
            {
                // get ig unique pdv types for set
                tPdvInterface->get_ig_unique_dv_types_for_set(
                        mMeshSet->get_set_index(),
                        aGeoPdvType );
            }
        }

        //------------------------------------------------------------------------------

        void
        Set::get_ip_dv_types_for_set(
                moris::Cell< moris::Cell< enum PDV_Type > >& aMatPdvType,
                mtk::Leader_Follower                         aIsLeader )
        {
            // choose based on the leader, follower type
            // the output here is gather from fem
            switch ( aIsLeader )
            {
                case mtk::Leader_Follower::LEADER:
                {
                    aMatPdvType = mLeaderDvTypes;
                    break;
                }

                case mtk::Leader_Follower::FOLLOWER:
                {
                    aMatPdvType = mFollowerDvTypes;
                    break;
                }

                default:
                    MORIS_ERROR( false, "Set::get_ip_dv_types_for_set - can only be leader or follower." );
            }

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

        void
        Set::populate_fields(
                moris::Cell< std::shared_ptr< fem::Field > >& aFieldToPopulate,
                moris::Cell< std::string > const &            aFieldIQINames )
        {
            this->gather_requested_IQIs(
                    aFieldIQINames,
                    mRequestedFieldIQIs,
                    mRequestedFieldIQIsGlobalIndices );

            for ( uint Ik = 0; Ik < mEquationObjList.size(); Ik++ )
            {
                mEquationObjList( Ik )->populate_fields( aFieldToPopulate, aFieldIQINames );
            }
        }

        //------------------------------------------------------------------------------

    } /* namespace fem */
} /* namespace moris */
