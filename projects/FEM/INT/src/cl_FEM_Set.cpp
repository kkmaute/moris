/*
 * cl_FEM_Set.cpp
 *
 *  Created on: Apr 11, 2019
 *      Author: schmidt
 */
#include <iostream>

#include "cl_MSI_Model_Solver_Interface.hpp" //FEM/MSI/src
#include "cl_FEM_Set.hpp"                    //FEM/INT/src
#include "cl_FEM_Element_Factory.hpp"        //FEM/INT/src
#include "cl_FEM_Integrator.hpp"             //FEM/INT/src

#include "cl_MTK_Set.hpp"             //FEM/INT/src

namespace moris
{
    namespace fem
    {
//------------------------------------------------------------------------------
        Set::Set( moris::mtk::Set                          * aSet,
                  enum fem::Element_Type                     aElementType,
                  moris::Cell< IWG* >                      & aIWGs,
                  moris::Cell< Node_Base* >                & aIPNodes) : mMeshSet(aSet),
                                                                         mNodes(aIPNodes),
                                                                         mIWGs( aIWGs ),
                                                                         mElementType( aElementType )
        {
            mMeshClusterList = aSet->get_clusters_on_set();
            // create a unique dof type list
            this->create_unique_dof_type_lists();

            // create a dof type list for field interpolators
            this->create_dof_type_lists();

            // create a unique property type list for field interpolators
            this->create_unique_property_type_list();

            // init the equation object list
            mEquationObjList.resize( mMeshClusterList.size(), nullptr);

            // create a factory to create fem cluster
            fem::Element_Factory tClusterFactory;

            for( luint k = 0; k < mMeshClusterList.size(); ++k )
            {
                // create a cluster
                mEquationObjList( k ) = tClusterFactory.create_cluster( mElementType,
                                                                        mMeshClusterList( k ),
                                                                        mNodes,
                                                                        this );
            }
        }

//------------------------------------------------------------------------------

    Set::~Set()
    {
        this->delete_pointers();

        for( auto tEquationObj : mEquationObjList )
        {
            delete tEquationObj;
        }
        mEquationObjList.clear();
    }

//------------------------------------------------------------------------------

    void Set::delete_pointers()
    {
        // delete the interpolation geometry interpolator pointer
        if ( mMasterIPGeometryInterpolator != nullptr )
        {
            delete mMasterIPGeometryInterpolator;
        }

        if ( mSlaveIPGeometryInterpolator != nullptr )
        {
            delete mSlaveIPGeometryInterpolator;
        }

        // delete the integration geometry interpolator pointer
        if ( mMasterIGGeometryInterpolator != nullptr )
        {
            delete mMasterIGGeometryInterpolator;
        }

        if ( mSlaveIGGeometryInterpolator != nullptr )
        {
            delete mSlaveIGGeometryInterpolator;
        }

        // delete the list of field interpolator pointers
        for( auto tMasterFieldInterpolators : mMasterFieldInterpolators )
        {
            delete tMasterFieldInterpolators;        }

        mMasterFieldInterpolators.clear();

        for( auto tSlaveFieldInterpolators : mSlaveFieldInterpolators )
        {
            delete tSlaveFieldInterpolators;
        }
        mSlaveFieldInterpolators.clear();
    }

//------------------------------------------------------------------------------

    void Set::finalize( MSI::Model_Solver_Interface * aModelSolverInterface )
    {
        this->delete_pointers();

        mIsTrivialMaster = mMeshSet->is_trivial( mtk::Master_Slave::MASTER );

        mIPGeometryType = mMeshSet->get_interpolation_cell_geometry_type();

        mIGGeometryType = mMeshSet->get_integration_cell_geometry_type();

        mIPSpaceInterpolationOrder = mMeshSet->get_interpolation_cell_interpolation_order();

        mIGSpaceInterpolationOrder = mMeshSet->get_integration_cell_interpolation_order();

        // time interpolation order for IP cells fixme not linear
        mIPTimeInterpolationOrder = mtk::Interpolation_Order::LINEAR;

        // time interpolation order for IG cells fixme not linear
        mIGTimeInterpolationOrder = mtk::Interpolation_Order::LINEAR;

        // geometry interpolation rule for interpolation cells
        Interpolation_Rule tIPGeometryInterpolationRule( mIPGeometryType,
                                                         Interpolation_Type::LAGRANGE,
                                                         mIPSpaceInterpolationOrder,
                                                         Interpolation_Type::LAGRANGE,
                                                         mIPTimeInterpolationOrder );

        // geometry interpolation rule for integration cells
        Interpolation_Rule tIGGeometryInterpolationRule( mIGGeometryType,
                                                         Interpolation_Type::LAGRANGE,
                                                         mIGSpaceInterpolationOrder,
                                                         Interpolation_Type::LAGRANGE,
                                                         mIGTimeInterpolationOrder );

        // if block-set
        switch ( mElementType )
        {
            case ( fem::Element_Type::BULK ):
            {
                // create an interpolation geometry intepolator
                mMasterIPGeometryInterpolator = new Geometry_Interpolator( tIPGeometryInterpolationRule, false );

                // create an integration geometry intepolator
                mMasterIGGeometryInterpolator = new Geometry_Interpolator( tIGGeometryInterpolationRule, false );

                // create the element field interpolators
                this->create_field_interpolators( aModelSolverInterface );

                // create the element dof assembly map
                this->create_dof_assembly_map();

                break;
            }

            // if side-set
            case( fem::Element_Type::SIDESET ):
            {
                // create an interpolation geometry intepolator
                mMasterIPGeometryInterpolator = new Geometry_Interpolator( tIPGeometryInterpolationRule, true );

                // create an integration geometry intepolator
                mMasterIGGeometryInterpolator = new Geometry_Interpolator( tIGGeometryInterpolationRule, true );

                // create the element field interpolators
                this->create_field_interpolators( aModelSolverInterface );

                // create the element dof assembly map
                this->create_dof_assembly_map();

                break;
            }

            // if double side-set
            case( fem::Element_Type::DOUBLE_SIDESET ):
            {
                mIsTrivialSlave = mMeshSet->is_trivial( mtk::Master_Slave::SLAVE );

                // create an interpolation geometry intepolator
                mMasterIPGeometryInterpolator = new Geometry_Interpolator( tIPGeometryInterpolationRule, true );
                mSlaveIPGeometryInterpolator  = new Geometry_Interpolator( tIPGeometryInterpolationRule, true );

                // create an integration geometry intepolator
                mMasterIGGeometryInterpolator = new Geometry_Interpolator( tIGGeometryInterpolationRule, true );
                mSlaveIGGeometryInterpolator  = new Geometry_Interpolator( tIGGeometryInterpolationRule, true );

                // create the element field interpolators
                this->create_field_interpolators_double( aModelSolverInterface );

                // create the element dof assembly map
                this->create_dof_assembly_map_double();

                break;
            }
            default:
            {
                MORIS_ERROR(false, "Set::finalize - unknown element type");
                break;
            }
        }

        // create an interpolation rule for the side
        Integration_Rule tIntegrationRule = Integration_Rule( mIGGeometryType,
                                                              Integration_Type::GAUSS,
                                                              this->get_auto_integration_order( mIGGeometryType ),
                                                              Integration_Type::GAUSS,
                                                              Integration_Order::BAR_1 ); // fixme time order

        // create an integrator
        Integrator tIntegrator( tIntegrationRule );

        //get number of integration points
        mNumOfIntegPoints = tIntegrator.get_number_of_points();

        // get integration points
        mIntegPoints = tIntegrator.get_points();

        // get integration weights
        mIntegWeights = tIntegrator.get_weights();
    }

//------------------------------------------------------------------------------

    void Set::create_unique_dof_type_lists()
    {
        // set the size of the element active dof type list
        uint tCounter = 0;
        for ( IWG * tIWG : mIWGs )
        {
            tCounter = tCounter + tIWG->get_residual_dof_type().size();
        }
        mEqnObjDofTypeList.reserve( tCounter );

        // loop over the IWGs
        tCounter = 0;
        for ( IWG * tIWG : mIWGs )
        {
            mEqnObjDofTypeList.append( tIWG->get_residual_dof_type() );
        }

        auto last = std::unique( ( mEqnObjDofTypeList.data() ).data(),
                                 ( mEqnObjDofTypeList.data() ).data() + mEqnObjDofTypeList.size() );
        auto pos  = std::distance( ( mEqnObjDofTypeList.data() ).data(), last );
        mEqnObjDofTypeList.resize( pos );
    }

//------------------------------------------------------------------------------

    void Set::create_dof_type_lists()
    {
        // get the number of IWGs
        uint tNumOfIWGs = this->get_num_IWG();

        // create a list of the groups of dof types provided by the IWGs----------------
        // FIXME works as long as the dof type are always grouped in the same way
        moris::Cell< MSI::Dof_Type > tInterpDofTypeListBuild( tNumOfIWGs );

        // loop over the IWGs
        for ( uint i = 0; i < tNumOfIWGs; i++ )
        {
            // get the first dof type of each group
            tInterpDofTypeListBuild( i ) = mIWGs( i )->get_residual_dof_type()( 0 );
        }

        // get a unique list of the first dof type of each group
        Cell< moris::moris_index > tUniqueDofTypeGroupsIndices = unique_index( tInterpDofTypeListBuild );

        // get the number of unique dof type groups
        uint tNumOfUniqueDofTypeGroupsIndices = tUniqueDofTypeGroupsIndices.size();

        // set the size of the list of unique dof type groups
        mInterpDofTypeList.resize( tNumOfUniqueDofTypeGroupsIndices );

        // loop over the list of unique dof type groups
        for ( uint i = 0; i < tNumOfUniqueDofTypeGroupsIndices; i++ )
        {
            // get the unique residual dof type groups
            mInterpDofTypeList( i ) = mIWGs( tUniqueDofTypeGroupsIndices( i ) )->get_residual_dof_type();
        }

        // create a map of the element active dof type list------------------------
        // set number of unique pdof type of the element
        mNumOfInterp = tNumOfUniqueDofTypeGroupsIndices;

        // get maximal dof type enum number
        sint tMaxDofTypeEnumNumber = 0;

        // loop over all pdof types to get the highest enum index
        for ( uint i = 0; i < mNumOfInterp; i++ )
        {
            tMaxDofTypeEnumNumber = std::max( tMaxDofTypeEnumNumber, static_cast< int >( mInterpDofTypeList( i )( 0 ) ) );
        }

        // +1 because c++ is 0 based
        tMaxDofTypeEnumNumber = tMaxDofTypeEnumNumber + 1;

        // set size of mapping matrix
        mInterpDofTypeMap.set_size( tMaxDofTypeEnumNumber, 1, -1 );

        // loop over all dof types to create the mapping matrix
        for ( uint i = 0; i < mNumOfInterp; i++ )
        {
            mInterpDofTypeMap( static_cast< int >( mInterpDofTypeList( i )( 0 ) ), 0 ) = i;
        }
    }

//------------------------------------------------------------------------------

    void Set::create_unique_property_type_list()
       {
           // set the size of the active mp type list
           uint tCounter = 0;
           for ( IWG * tIWG : mIWGs )
           {
               tCounter = tCounter + tIWG->get_active_property_types().size();
           }
           mPropertyTypeList.reserve( tCounter );

           // loop over the IWGs
           tCounter = 0;
           for ( IWG * tIWG : mIWGs )
           {
               mPropertyTypeList.append( tIWG->get_active_property_types() );
           }

           auto last = std::unique( ( mPropertyTypeList.data() ).data(), ( mPropertyTypeList.data() ).data() + mPropertyTypeList.size() );
           auto pos  = std::distance( ( mPropertyTypeList.data() ).data(), last );
           mPropertyTypeList.resize( pos );
       }

//------------------------------------------------------------------------------

    void Set::create_dof_assembly_map( )
    {
        // set size of assembly mapping matrix
        mInterpDofAssemblyMap.set_size( mNumOfInterp, 2, -1 );

        // init dof counter
        uint tDofCounter = 0;

        // loop on the dof type groups and create a field interpolator for each
        for( uint i = 0; i < mNumOfInterp; i++ )
        {
            // fill the assembly map with starting dof counter
            mInterpDofAssemblyMap( i, 0 ) = tDofCounter;

            // update dof counter
            tDofCounter = tDofCounter + mMasterFieldInterpolators( i )->get_number_of_space_time_coefficients()-1;

            // fill the assembly map with starting dof counter
            mInterpDofAssemblyMap( i, 1 ) = tDofCounter;

            // update dof counter
            tDofCounter = tDofCounter + 1;
        }

        // set mTotalDof
        mTotalDof = tDofCounter;
    }

    void Set::create_dof_assembly_map_double( )
    {
        // set size of assembly mapping matrix
        mInterpDofAssemblyMap.set_size( mNumOfInterp, 2, -1 );

        // init dof counter
        uint tDofCounter = 0;

        // loop on the dof type groups and create a field interpolator for each
        for( uint i = 0; i < mNumOfInterp; i++ )
        {
            // fill the assembly map with starting dof counter
            mInterpDofAssemblyMap( i, 0 ) = tDofCounter;

            // update dof counter
            //FIXME works if left and right field interpolators are the same
            tDofCounter = tDofCounter + 2 * mMasterFieldInterpolators( i )->get_number_of_space_time_coefficients()-1;

            // fill the assembly map with starting dof counter
            mInterpDofAssemblyMap( i, 1 ) = tDofCounter;

            // update dof counter
            tDofCounter = tDofCounter + 1;
        }

        // set mTotalDof (two times what we have on the left side)
        mTotalDof = tDofCounter;
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

    void Set::create_field_interpolators( MSI::Model_Solver_Interface * aModelSolverInterface )
    {
        // cell of field interpolators
        mMasterFieldInterpolators.resize( mInterpDofTypeList.size() + mPropertyTypeList.size(), nullptr );

        // create a field interpolator for each dof type
        for( uint i = 0; i < mInterpDofTypeList.size(); i++ )
        {
            // get the ith dof type group
            Cell< MSI::Dof_Type > tDofTypeGroup = mInterpDofTypeList( i );

            moris::uint tNumTimeNodes = aModelSolverInterface->get_time_levels_for_type( tDofTypeGroup( 0 ) );

            // create the field interpolation rule for the ith dof type group
            Interpolation_Rule tFieldInterpolationRule( mIPGeometryType,
                                                        Interpolation_Type::LAGRANGE,
                                                        mIPSpaceInterpolationOrder,
                                                        this->get_auto_time_interpolation_type( tNumTimeNodes ), // fixme
                                                        // If interpolation type CONSTANT, iInterpolation order is not used
                                                        this->get_auto_interpolation_order( tNumTimeNodes, mtk::Geometry_Type::LINE ) ); //fixme

            // get number of field interpolated by the ith field interpolator
            uint tNumOfFields = tDofTypeGroup.size();

            // create an interpolator for the ith dof type group
            mMasterFieldInterpolators( i ) = new Field_Interpolator( tNumOfFields,
                                                                     tFieldInterpolationRule,
                                                                     mMasterIPGeometryInterpolator );
        }

        // create a field interpolator for each property type
        for( uint i = 0; i < mPropertyTypeList.size(); i++ )
        {
            // get the ith property type
            // fem::Property_Type tPropertyType = mPropertyTypeList( i );

            //moris::uint tNumTimeNodes = aModelSolverInterface->get_time_levels_for_type( tPropertyType );

            // create the field interpolation rule for the ith dof type group
            Interpolation_Rule tFieldInterpolationRule( mIPGeometryType,
                                                        Interpolation_Type::LAGRANGE,
                                                        mIPSpaceInterpolationOrder,
                                                        this->get_auto_time_interpolation_type( 2 ), // fixme
                                                        // If interpolation type CONSTANT, iInterpolation order is not used
                                                        this->get_auto_interpolation_order( 2, mtk::Geometry_Type::LINE ) ); //fixme

            // get number of field interpolated by the ith field interpolator
            uint tNumOfFields = 1; // fixme

            // create an interpolator for the ith dof type group
            mMasterFieldInterpolators( mInterpDofTypeList.size() + i ) = new Field_Interpolator( tNumOfFields,
                                                                                                 tFieldInterpolationRule,
                                                                                                 mMasterIPGeometryInterpolator );
        }
    }

    void Set::create_field_interpolators_double( MSI::Model_Solver_Interface * aModelSolverInterface )
    {
        // init the lists of field interpolators for left and right IP cells
        mMasterFieldInterpolators.resize( mNumOfInterp, nullptr );
        mSlaveFieldInterpolators.resize( mNumOfInterp, nullptr );

        // loop on the dof type groups and create a field interpolator for each
        for( uint i = 0; i < mNumOfInterp; i++ )
        {
            // get the ith dof type group
            Cell< MSI::Dof_Type > tDofTypeGroup = mInterpDofTypeList( i );

            moris::uint tNumTimeNodes = aModelSolverInterface->get_time_levels_for_type( tDofTypeGroup( 0 ) );

            // create the field interpolation rule for the ith dof type group
            Interpolation_Rule tFieldInterpolationRule( mIPGeometryType,
                                                        Interpolation_Type::LAGRANGE,
                                                        mIPSpaceInterpolationOrder,
                                                        this->get_auto_time_interpolation_type( tNumTimeNodes ), // fixme
                                                        // If interpolation type CONSTANT, iInterpolation order is not used
                                                        this->get_auto_interpolation_order( tNumTimeNodes, mtk::Geometry_Type::LINE ) ); //fixme

            // get number of field interpolated by the ith field interpolator
            uint tNumOfFields = tDofTypeGroup.size();

            // create an interpolator for the ith dof type group
            mMasterFieldInterpolators( i ) = new Field_Interpolator( tNumOfFields,
                                                                     tFieldInterpolationRule,
                                                                     mMasterIPGeometryInterpolator );
            mSlaveFieldInterpolators( i )  = new Field_Interpolator( tNumOfFields,
                                                                     tFieldInterpolationRule,
                                                                     mSlaveIPGeometryInterpolator );
        }
    }

//------------------------------------------------------------------------------

    moris::Cell< Field_Interpolator* > Set::get_IWG_field_interpolators ( IWG*                               & aIWG,
                                                                          moris::Cell< Field_Interpolator* > & aFieldInterpolators )
    {
        // ask the IWG for its active dof types
        moris::Cell< moris::Cell< MSI::Dof_Type > > tIWGActiveDof = aIWG->get_active_dof_types();

        // number of active dof type for the IWG
        uint tNumOfIWGActiveDof = tIWGActiveDof.size();

        // select associated active interpolators
        moris::Cell< Field_Interpolator* > tIWGFieldInterpolators( tNumOfIWGActiveDof, nullptr );
        for( uint i = 0; i < tNumOfIWGActiveDof; i++ )
        {
            // find the index of active dof type in the list of element dof type
            uint tIWGDofIndex = mInterpDofTypeMap( static_cast< int >( tIWGActiveDof( i )( 0 ) ) );

            // select the corresponding interpolator
            tIWGFieldInterpolators( i ) = aFieldInterpolators( tIWGDofIndex );
        }
        return tIWGFieldInterpolators;
    }

//------------------------------------------------------------------------------

    Field_Interpolator*
    Set::get_dof_type_field_interpolators ( enum MSI::Dof_Type aDofType )
    {
        uint tIndex = mInterpDofTypeMap( static_cast< int >( aDofType ) );

        return mMasterFieldInterpolators(tIndex);
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
                 return fem::Integration_Order::QUAD_3x3;
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

    void Set::initialize_mJacobian()
    {
        if ( !mJacobianExist )
        {
            uint tTotalDof = this->get_total_number_of_dofs();
            mJacobian.set_size( tTotalDof, tTotalDof, 0.0 );

            mJacobianExist = true;
        }
        else
        {
            MORIS_ASSERT( mJacobian.numel() > 0, "Set::initialize_mJacobian() - Jacobian not properly initialized.");

            mJacobian.fill( 0.0 );
        }
    }

//------------------------------------------------------------------------------

    void Set::initialize_mResidual()
    {
        if ( !mResidualExist )
        {
            mResidual.set_size( this->get_total_number_of_dofs(), 1, 0.0 );

            mResidualExist = true;
        }
        else
        {
            MORIS_ASSERT( mResidual.numel() > 0, "Set::initialize_mResidual() - Residual not properly initialized.");

            mResidual.fill( 0.0 );
        }
    }

//------------------------------------------------------------------------------


    } /* namespace fem */
} /* namespace moris */
