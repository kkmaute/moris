/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_FEM_Field_Interpolator_Manager.cpp
 *
 */

#include "cl_FEM_Field_Interpolator_Manager.hpp"    //FEM/INT/src
#include "cl_FEM_Set.hpp"                           //FEM/INT/src

namespace moris
{
    namespace fem
    {

        //------------------------------------------------------------------------------

        Field_Interpolator_Manager::Field_Interpolator_Manager(
                const moris::Cell< moris::Cell< enum MSI::Dof_Type > >& aDofTypes,
                MSI::Equation_Set*                                      aEquationSet,
                mtk::Leader_Follower                                       aIsLeader )
                : mDofTypes( aDofTypes )
                , mEquationSet( aEquationSet )
                , mIsLeader( aIsLeader )
        {
            // set the dof type map
            mDofTypeMap = mEquationSet->get_dof_type_map( aIsLeader );

            // maximum number of dof field interpolators
            mMaxNumDofFI = mEquationSet->get_num_unique_dof_types();

            // FIXME default
            mMaxNumDvFI = 0;
        }

        Field_Interpolator_Manager::Field_Interpolator_Manager(
                const moris::Cell< moris::Cell< enum MSI::Dof_Type > >&   aDofTypes,
                const moris::Cell< moris::Cell< enum PDV_Type > >&        aDvTypes,
                const moris::Cell< moris::Cell< enum mtk::Field_Type > >& aFieldTypes,
                MSI::Equation_Set*                                        aEquationSet,
                mtk::Leader_Follower                                         aIsLeader )
                : mDofTypes( aDofTypes )
                , mEquationSet( aEquationSet )
                , mIsLeader( aIsLeader )
                , mDvTypes( aDvTypes )
                , mFieldTypes( aFieldTypes )
        {
            // set the dof type map
            mDofTypeMap = mEquationSet->get_dof_type_map( aIsLeader );

            // maximum number of dof field interpolators
            mMaxNumDofFI = mEquationSet->get_num_unique_dof_types();

            // FIXME maximum number of dv field interpolators
            mMaxNumDvFI = 3;    // FIXME FIXME FIXME

            mDvTypeMap = mEquationSet->get_dv_type_map( aIsLeader );

            mFieldTypeMap = mEquationSet->get_field_type_map( aIsLeader );

            mMaxNumFieldFI = mEquationSet->get_num_unique_field_types();
        }

        Field_Interpolator_Manager::Field_Interpolator_Manager(
                const moris::Cell< moris::Cell< enum MSI::Dof_Type > >& aDofTypes,
                MSI::Equation_Set*                                      aEquationSet,
                MSI::Model_Solver_Interface*                            aModelSolverInterface,
                mtk::Leader_Follower                                       aIsLeader )
                : mDofTypes( aDofTypes )
                , mEquationSet( aEquationSet )
                , mIsLeader( aIsLeader )
        {
            // set the dof type map
            mDofTypeMap = mEquationSet->get_dof_type_map( aIsLeader );

            // maximum number of dof field interpolators
            mMaxNumDofFI = mEquationSet->get_num_unique_dof_types();

            // FIXME default
            mMaxNumDvFI = 0;
        }

        //------------------------------------------------------------------------------

        Field_Interpolator_Manager::~Field_Interpolator_Manager()
        {
            // delete pointers on the FI manager
            this->delete_pointers();
        }

        //------------------------------------------------------------------------------

        void
        Field_Interpolator_Manager::delete_pointers()
        {
            // delete the dof field interpolator pointers
            for ( Field_Interpolator* tFI : mFI )
            {
                delete tFI;
            }
            mFI.clear();

            // delete the dv field interpolator pointers
            for ( Field_Interpolator* tFI : mDvFI )
            {
                delete tFI;
            }
            mDvFI.clear();

            // delete the field field interpolator pointers
            for ( Field_Interpolator* tFI : mFieldFI )
            {
                delete tFI;
            }
            mFieldFI.clear();

            // delete the IP geometry interpolator pointer
            if ( mIPGeometryInterpolator != nullptr && mGeometryInterpolatorOwned )
            {
                delete mIPGeometryInterpolator;
                mIPGeometryInterpolator = nullptr;
            }

            // delete the IG geometry interpolator pointer
            if ( mIGGeometryInterpolator != nullptr && mGeometryInterpolatorOwned )
            {
                delete mIGGeometryInterpolator;
                mIGGeometryInterpolator = nullptr;
            }
        }

        //------------------------------------------------------------------------------

        void
        Field_Interpolator_Manager::create_field_interpolators(
                MSI::Model_Solver_Interface* aModelSolverInterface,
                uint                         aNumSolutionSets )
        {
            // dof field interpolators------------------------------------------

            // store number of solution sets
            mNumSolutionSets = aNumSolutionSets;

            // set the size of the cell of field interpolators
            mFI.resize( mMaxNumDofFI * mNumSolutionSets, nullptr );

            // loop over the dof type groups
            for ( uint iDof = 0; iDof < mDofTypes.size(); iDof++ )
            {
                // get the number of time level for the dof type group
                uint tNumTimeNodes = aModelSolverInterface->get_time_levels_for_type( mDofTypes( iDof )( 0 ) );

                // get the set index for the dof type group
                uint tDofIndex = mEquationSet->get_dof_index_for_type_1( mDofTypes( iDof )( 0 ), mIsLeader );

                // create the field interpolation rule for the dof type group
                mtk::Interpolation_Rule tFieldInterpolationRule(
                        reinterpret_cast< Set* >( mEquationSet )->mIPGeometryType,
                        mtk::Interpolation_Type::LAGRANGE,
                        reinterpret_cast< Set* >( mEquationSet )->mIPSpaceInterpolationOrder,
                        reinterpret_cast< Set* >( mEquationSet )->get_auto_time_interpolation_type( tNumTimeNodes ),    // fixme
                        // If interpolation type CONSTANT, iInterpolation order is not used
                        reinterpret_cast< Set* >( mEquationSet )->get_auto_interpolation_order( tNumTimeNodes, mtk::Geometry_Type::LINE ) );    // fixme

                // get the discretization mesh index for the DoF type
                moris_index tMeshIndexForDof = aModelSolverInterface->get_adof_index_for_type( mDofTypes( iDof )( 0 ) );

                // loop over all solution sets
                for ( uint is = 0; is < mNumSolutionSets; ++is )
                {
                    // compute field interpolator index
                    uint tFiIndex = tDofIndex + is * mMaxNumDofFI;

                    // check if the field interpolator was created previously
                    MORIS_ASSERT( mFI( tFiIndex ) == nullptr,
                            "Field_Interpolator_Manager::create_field_interpolators - Field interpolator was created previously" );

                    // create a field interpolator for the dof type group
                    mFI( tFiIndex ) = new Field_Interpolator(
                            mDofTypes( iDof ).size(),
                            tFieldInterpolationRule,
                            mIPGeometryInterpolator,
                            mDofTypes( iDof ) );

                    // discretization mesh index it in the field interpolator
                    mFI( tFiIndex )->set_discretization_mesh_index( tMeshIndexForDof );
                }
            }

            // dv field interpolators------------------------------------------

            // set the size of the cell of field interpolators
            mDvFI.resize( mMaxNumDvFI, nullptr );

            // loop over the dv type groups
            for ( uint iDv = 0; iDv < mDvTypes.size(); iDv++ )
            {
                // get the number of time level for the dv type group
                // FIXME where do we get this info
                uint tNumTimeNodes = 1;

                // get the set index for the dv type group
                uint tDvIndex = mEquationSet->get_dv_index_for_type_1( mDvTypes( iDv )( 0 ), mIsLeader );

                // create the field interpolation rule for the dv type group
                mtk::Interpolation_Rule tFieldInterpolationRule(
                        reinterpret_cast< Set* >( mEquationSet )->mIPGeometryType,
                        mtk::Interpolation_Type::LAGRANGE,
                        reinterpret_cast< Set* >( mEquationSet )->mIPSpaceInterpolationOrder,
                        reinterpret_cast< Set* >( mEquationSet )->get_auto_time_interpolation_type( tNumTimeNodes ),    // fixme
                        // If interpolation type CONSTANT, iInterpolation order is not used
                        reinterpret_cast< Set* >( mEquationSet )->get_auto_interpolation_order( tNumTimeNodes, mtk::Geometry_Type::LINE ) );    // fixme

                // check if the field interpolator was created previously
                MORIS_ASSERT( mDvFI( tDvIndex ) == nullptr,
                        "Field_Interpolator_Manager::create_field_interpolators - Field interpolator was created previously." );

                // create a field interpolator for the dof type group
                mDvFI( tDvIndex ) = new Field_Interpolator(
                        mDvTypes( iDv ).size(),
                        tFieldInterpolationRule,
                        mIPGeometryInterpolator,
                        mDvTypes( iDv ) );
            }

            // field field interpolators------------------------------------------

            // set the size of the cell of field interpolators
            mFieldFI.resize( mMaxNumFieldFI, nullptr );

            // loop over the field type groups
            for ( uint iFi = 0; iFi < mFieldTypes.size(); iFi++ )
            {
                // get the number of time level for the dv type group
                // FIXME where do we get this info
                uint tNumTimeNodes = 1;

                // get the set index for the dv type group
                uint tFieldIndex = mEquationSet->get_field_index_for_type_1( mFieldTypes( iFi )( 0 ), mIsLeader );

                // create the field interpolation rule for the dv type group
                mtk::Interpolation_Rule tFieldInterpolationRule(
                        reinterpret_cast< Set* >( mEquationSet )->mIPGeometryType,
                        mtk::Interpolation_Type::LAGRANGE,
                        reinterpret_cast< Set* >( mEquationSet )->mIPSpaceInterpolationOrder,
                        reinterpret_cast< Set* >( mEquationSet )->get_auto_time_interpolation_type( tNumTimeNodes ),    // fixme
                        // If interpolation type CONSTANT, iInterpolation order is not used
                        reinterpret_cast< Set* >( mEquationSet )->get_auto_interpolation_order( tNumTimeNodes, mtk::Geometry_Type::LINE ) );    // fixme

                // check if the field interpolator was created previously
                MORIS_ASSERT( mFieldFI( tFieldIndex ) == nullptr,
                        "Field_Interpolator_Manager::create_field_interpolators - Field interpolator was created previously." );

                // create a field interpolator for the dof type group
                mFieldFI( tFieldIndex ) = new Field_Interpolator(
                        mFieldTypes( iFi ).size(),
                        tFieldInterpolationRule,
                        mIPGeometryInterpolator,
                        mFieldTypes( iFi ) );
            }
        }

        //------------------------------------------------------------------------------

        void
        Field_Interpolator_Manager::create_geometry_interpolators()
        {
            // get element type for set
            fem::Element_Type tElementType = reinterpret_cast< Set* >( mEquationSet )->mElementType;

            // bool true if time sideset
            bool tIsTimeSide = ( tElementType == fem::Element_Type::TIME_SIDESET );

            // bool true if sideset or double sideset
            bool tIsSide = ( tElementType != fem::Element_Type::BULK ) &&    //
                           ( tElementType != fem::Element_Type::TIME_SIDESET );

            // create geometry interpolation rule for IP elements
            mtk::Interpolation_Rule tIPGeometryInterpolationRule(
                    reinterpret_cast< Set* >( mEquationSet )->mIPGeometryType,
                    mtk::Interpolation_Type::LAGRANGE,
                    reinterpret_cast< Set* >( mEquationSet )->mIPSpaceInterpolationOrder,
                    mtk::Geometry_Type::LINE,
                    mtk::Interpolation_Type::LAGRANGE,
                    mtk::Interpolation_Order::LINEAR );    // FIXME not linear?

            // FIXME default should be given by the MSI
            mtk::Geometry_Type       tIGTimeGeometryType = mtk::Geometry_Type::LINE;
            mtk::Interpolation_Type  tIGTimeInterpType   = mtk::Interpolation_Type::LAGRANGE;
            mtk::Interpolation_Order tIGTimeInterpOrder  = mtk::Interpolation_Order::LINEAR;

            // if time sideset
            if ( tIsTimeSide )
            {
                tIGTimeGeometryType = mtk::Geometry_Type::POINT;
                tIGTimeInterpType   = mtk::Interpolation_Type::CONSTANT;
                tIGTimeInterpOrder  = mtk::Interpolation_Order::CONSTANT;
            }

            // create geometry interpolation rule for IG elements
            mtk::Interpolation_Rule tIGGeometryInterpolationRule(
                    reinterpret_cast< Set* >( mEquationSet )->mIGGeometryType,
                    mtk::Interpolation_Type::LAGRANGE,
                    reinterpret_cast< Set* >( mEquationSet )->mIGSpaceInterpolationOrder,
                    tIGTimeGeometryType,
                    tIGTimeInterpType,
                    tIGTimeInterpOrder );

            // check that geometry interpolators are not already initialized
            MORIS_ASSERT( mIPGeometryInterpolator == nullptr && mIGGeometryInterpolator == nullptr,
                    "Field_Interpolator_Manager::create_geometry_interpolators - geometry interpolators already initialized" );

            // create a geometry interpolator for IP cells
            mIPGeometryInterpolator = new Geometry_Interpolator(
                    tIPGeometryInterpolationRule,
                    mIPCellShape,
                    false,
                    false );

            // create a geometry interpolator for IG cells
            // IG interpolation requires knowledge about  the IP Element in
            // order to have an appropriate mapping for integration points,
            // tIPGeometryInterpolationRule
            mIGGeometryInterpolator = new Geometry_Interpolator(
                    tIGGeometryInterpolationRule,
                    tIPGeometryInterpolationRule,
                    mIGCellShape,
                    tIsSide,
                    tIsTimeSide );

            // set flag that Field_Interpolator_Manager owns pointers to GeometryInterpolators
            mGeometryInterpolatorOwned = true;
        }

        //------------------------------------------------------------------------------

        Field_Interpolator*
        Field_Interpolator_Manager::get_field_interpolators_for_type(
                const enum MSI::Dof_Type aDofType,
                const uint               aSolutionSetIndex )
        {
            // check of the equation set pointer was set for the FI manager
            MORIS_ASSERT( mEquationSet != nullptr,
                    "Field_Interpolator_Manager::get_field_interpolators_for_type - Equation Set pointer not set" );

            // get the set index for the requested dof type
            sint tDofIndex = mEquationSet->get_dof_index_for_type_1( aDofType, mIsLeader );

            // if the index was set for the equation set
            if ( tDofIndex != -1 )
            {
                // compute field interpolator index
                uint tFiIndex = tDofIndex + aSolutionSetIndex * mMaxNumDofFI;

                // check if the FI exists for the FI manager
                MORIS_ASSERT( mFI.size() > tFiIndex,
                        "Field_Interpolator_Manager::get_field_interpolators_for_type - field interpolator does not exist" );

                // return the FI
                return mFI( tFiIndex );
            }
            else
            {
                return nullptr;
            }
        }

        //------------------------------------------------------------------------------

        Field_Interpolator*
        Field_Interpolator_Manager::get_field_interpolators_for_type(
                enum PDV_Type aDvType )
        {
            // get the set index for the requested dv type
            sint tDvIndex = mEquationSet->get_dv_index_for_type_1( aDvType, mIsLeader );

            // if the index was set for the equation set
            if ( tDvIndex != -1 )
            {
                // check if the FI exists for the FI manager
                MORIS_ASSERT( (sint)mDvFI.size() > tDvIndex,
                        "Field_Interpolator_Manager::get_field_interpolators_for_type - field interpolator does not exist" );

                // return the FI
                return mDvFI( tDvIndex );
            }
            else
            {
                return nullptr;
            }
        }

        //------------------------------------------------------------------------------

        Field_Interpolator*
        Field_Interpolator_Manager::get_field_interpolators_for_type(
                enum mtk::Field_Type aFieldType )
        {
            // get the set index for the requested dv type
            sint tFieldIndex = mEquationSet->get_field_index_for_type_1( aFieldType, mIsLeader );

            // if the index was set for the equation set
            if ( tFieldIndex != -1 )
            {
                // check if the FI exists for the FI manager
                MORIS_ASSERT( (sint)mFieldFI.size() > tFieldIndex,
                        "Field_Interpolator_Manager::get_field_interpolators_for_type - field interpolator does not exist" );

                // return the FI
                return mFieldFI( tFieldIndex );
            }
            else
            {
                return nullptr;
            }
        }

        //------------------------------------------------------------------------------

        void
        Field_Interpolator_Manager::set_space_time(
                const Matrix< DDRMat >& aParamPoint )
        {
            // loop over the dof field interpolators
            for ( uint iDofFI = 0; iDofFI < mDofTypes.size(); iDofFI++ )
            {
                // get the set index for the dof type
                sint tDofIndex = mDofTypeMap( static_cast< uint >( mDofTypes( iDofFI )( 0 ) ) );

                // loop over all solution sets
                for ( uint is = 0; is < mNumSolutionSets; ++is )
                {
                    // compute field interpolator index
                    uint tFiIndex = tDofIndex + is * mMaxNumDofFI;

                    // check if the FI exists for the FI manager
                    MORIS_ASSERT( mFI.size() > tFiIndex,
                            "Field_Interpolator_Manager::get_field_interpolators_for_type - field interpolator does not exist" );

                    // set the evaluation point
                    mFI( tFiIndex )->set_space_time( aParamPoint );
                }
            }

            // loop over the dv field interpolators
            for ( uint iDvFI = 0; iDvFI < mDvTypes.size(); iDvFI++ )
            {
                // get the set index for the dv type
                sint tDvIndex = mDvTypeMap( static_cast< uint >( mDvTypes( iDvFI )( 0 ) ) );

                // set the evaluation point
                mDvFI( tDvIndex )->set_space_time( aParamPoint );
            }

            // loop over the field field interpolators
            for ( uint iFieldFI = 0; iFieldFI < mFieldTypes.size(); iFieldFI++ )
            {
                // get the set index for the field type
                sint tFieldIndex = mFieldTypeMap( static_cast< uint >( mFieldTypes( iFieldFI )( 0 ) ) );

                // set the evaluation point
                mFieldFI( tFieldIndex )->set_space_time( aParamPoint );
            }

            // IP geometry interpolator
            mIPGeometryInterpolator->set_space_time( aParamPoint );
        }

        //------------------------------------------------------------------------------

        void
        Field_Interpolator_Manager::set_space_time_from_local_IG_point(
                const Matrix< DDRMat >& aLocalParamPoint )
        {
            // set evaluation point in the IG param space for IG geometry interpolator
            mIGGeometryInterpolator->set_space_time( aLocalParamPoint );

            // bring evaluation point in the IP param space
            const Matrix< DDRMat >& tGlobalParamPoint =
                    mIGGeometryInterpolator->map_integration_point();

            // set evaluation point for interpolators (FIs and IP GI)
            this->set_space_time( tGlobalParamPoint );
        }

        //------------------------------------------------------------------------------

        void
        Field_Interpolator_Manager::set_coeff_for_type(
                const enum MSI::Dof_Type aDofType,
                const Matrix< DDRMat >&  aCoeff,
                const uint               aSolutionSetIndex )
        {
            // get field interpolator for dof type and set coefficients
            this->get_field_interpolators_for_type( aDofType, aSolutionSetIndex )->set_coeff( aCoeff );
        }

        //------------------------------------------------------------------------------

        void
        Field_Interpolator_Manager::set_coeff_for_type(
                enum PDV_Type           aDvType,
                const Matrix< DDRMat >& aCoeff )
        {
            // get field interpolator for dof type and set coefficients
            this->get_field_interpolators_for_type( aDvType )->set_coeff( aCoeff );
        }

        //------------------------------------------------------------------------------

        void
        Field_Interpolator_Manager::set_coeff_for_type(
                enum mtk::Field_Type    aFieldType,
                const Matrix< DDRMat >& aCoeff )
        {
            // get field interpolator for dof type and set coefficients
            this->get_field_interpolators_for_type( aFieldType )->set_coeff( aCoeff );
        }

        //------------------------------------------------------------------------------

        void
        Field_Interpolator_Manager::set_IG_cell_shape( enum CellShape aCellShape )
        {
            mIGCellShape = aCellShape;
        }

        //------------------------------------------------------------------------------

        void
        Field_Interpolator_Manager::set_IP_cell_shape( enum CellShape aCellShape )
        {
            mIPCellShape = aCellShape;
        }

        //------------------------------------------------------------------------------

    } /* namespace fem */
} /* namespace moris */
