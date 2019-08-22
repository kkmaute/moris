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

#include "fn_equal_to.hpp"

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
            delete tMasterFieldInterpolators;
        }
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

                // create the element IWG info
                this->create_IWG_set_info();

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

                // create the element IWG info
                this->create_IWG_set_info();

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

                // create the element IWG info
                this->create_IWG_set_info_double();

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
        tIntegrator.get_points( mIntegPoints );

        // get integration weights
        tIntegrator.get_weights( mIntegWeights );
    }

//------------------------------------------------------------------------------

//    void Set::create_unique_dof_type_lists()
//    {
//        // set the size of the element active dof type list
//        uint tCounter = 0;
//        for ( IWG * tIWG : mIWGs )
//        {
//            tCounter = tCounter + tIWG->get_residual_dof_type().size();
//        }
//        mEqnObjDofTypeList.reserve( tCounter );
//
//        // loop over the IWGs
//        tCounter = 0;
//        for ( IWG * tIWG : mIWGs )
//        {
//            mEqnObjDofTypeList.append( tIWG->get_residual_dof_type() );
//        }
//
//        auto last = std::unique( ( mEqnObjDofTypeList.data() ).data(),
//                                 ( mEqnObjDofTypeList.data() ).data() + mEqnObjDofTypeList.size() );
//        auto pos  = std::distance( ( mEqnObjDofTypeList.data() ).data(), last );
//        mEqnObjDofTypeList.resize( pos );
//    }

    void Set::create_unique_dof_type_lists()
    {
        // set the size of the element active dof type list
        uint tCounter = 0;
        for ( IWG * tIWG : mIWGs )
        {
            // get active dof type
            Cell< Cell< MSI::Dof_Type > > tActiveDofType = tIWG->get_active_dof_types();

            for ( uint iDOF = 0; iDOF< tActiveDofType.size(); iDOF++ )
            {
                tCounter += tActiveDofType( iDOF ).size();
            }
        }
        mEqnObjDofTypeList.reserve( tCounter );

        // loop over the IWGs
        tCounter = 0;
        for ( IWG * tIWG : mIWGs )
        {
            // get active dof type
            Cell< Cell< MSI::Dof_Type > > tActiveDofType = tIWG->get_active_dof_types();

            for ( uint iDOF = 0; iDOF< tActiveDofType.size(); iDOF++ )
            {
                mEqnObjDofTypeList.append( tActiveDofType( iDOF ) );
            }
        }

        auto last = std::unique( ( mEqnObjDofTypeList.data() ).data(),
                                 ( mEqnObjDofTypeList.data() ).data() + mEqnObjDofTypeList.size() );
        auto pos  = std::distance( ( mEqnObjDofTypeList.data() ).data(), last );
        mEqnObjDofTypeList.resize( pos );
    }

//------------------------------------------------------------------------------

//    void Set::create_dof_type_lists()
//    {
//        // get the number of IWGs
//        uint tNumOfIWGs = this->get_num_IWG();
//
//        // create a list of the groups of dof types provided by the IWGs----------------
//        // FIXME works as long as the dof type are always grouped in the same way
//        moris::Cell< MSI::Dof_Type > tInterpDofTypeListBuild( tNumOfIWGs );
//
//        // loop over the IWGs
//        for ( uint i = 0; i < tNumOfIWGs; i++ )
//        {
//            // get the first dof type of each group
//            tInterpDofTypeListBuild( i ) = mIWGs( i )->get_residual_dof_type()( 0 );
//        }
//
//        // get a unique list of the first dof type of each group
//        Cell< moris::moris_index > tUniqueDofTypeGroupsIndices = unique_index( tInterpDofTypeListBuild );
//
//        // get the number of unique dof type groups, i.e. the number of interpolators
//        mNumOfInterp = tUniqueDofTypeGroupsIndices.size();
//
//        // set the size of the list of unique dof type groups
//        mInterpDofTypeList.resize( mNumOfInterp );
//
//        // get maximal dof type enum
//        sint tMaxEnum = 0;
//
//        // loop over the list of unique dof type groups
//        for ( uint i = 0; i < mNumOfInterp; i++ )
//        {
//            // get the unique residual dof type groups
//            mInterpDofTypeList( i ) = mIWGs( tUniqueDofTypeGroupsIndices( i ) )->get_residual_dof_type();
//
//            // get the highest dof type enum
//            tMaxEnum = std::max( tMaxEnum, static_cast< int >( mInterpDofTypeList( i )( 0 ) ) );
//        }
//
//        // +1 because c++ is 0 based
//        tMaxEnum++;
//
//        // create a map of the set active dof type list------------------------
//        // set size of mapping matrix
//        mInterpDofTypeMap.set_size( tMaxEnum, 1, -1 );
//
//        // loop over all dof types to create the mapping matrix
//        for ( uint i = 0; i < mNumOfInterp; i++ )
//        {
//            mInterpDofTypeMap( static_cast< int >( mInterpDofTypeList( i )( 0 ) ), 0 ) = i;
//        }
//    }

    void Set::create_dof_type_lists()
    {
        // set the size of the element active dof type list
        uint tCounterMax = 0;
        for ( IWG * tIWG : mIWGs )
        {
            tCounterMax += tIWG->get_active_dof_types().size();
        }
        mInterpDofTypeList.resize( tCounterMax );
        moris::Cell< sint > tCheckList( tCounterMax, -1 );

        // get maximal dof type enum
        sint tMaxEnum = 0;

        uint tCounter = 0;
        for ( IWG * tIWG : mIWGs )
        {
            // get active dof type
            Cell< Cell< MSI::Dof_Type > > tActiveDofType = tIWG->get_active_dof_types();

            for ( uint iDOF = 0; iDOF < tActiveDofType.size(); iDOF++ )
            {
                // check enum is not already in the list
                bool tCheck = false;
                for( uint i = 0; i < tCounter; i++ )
                {
                    tCheck = tCheck || equal_to( tCheckList( i ), static_cast< uint >( tActiveDofType( iDOF )( 0 ) ) );
                }

                // if dof enum not in the list
                if ( !tCheck )
                {
                    tCheckList( tCounter ) = static_cast< uint >( tActiveDofType( iDOF )( 0 ) );
                    mInterpDofTypeList( tCounter ) = tActiveDofType( iDOF );
                    tCounter++;

                    // get the highest dof type enum
                    tMaxEnum = std::max( tMaxEnum, static_cast< int >( tActiveDofType( iDOF )( 0 ) ) );
                }
            }
        }

        // get the number of unique dof type groups, i.e. the number of interpolators
        mNumOfInterp = tCounter;
        mInterpDofTypeList.resize( mNumOfInterp );
        tMaxEnum++;

        // create a map of the set active dof type list------------------------
        // set size of mapping matrix
        mInterpDofTypeMap.set_size( tMaxEnum, 1, -1 );

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
            tDofCounter += mMasterFieldInterpolators( i )->get_number_of_space_time_coefficients() - 1;

            // fill the assembly map with starting dof counter
            mInterpDofAssemblyMap( i, 1 ) = tDofCounter;

            // update dof counter
            tDofCounter++;
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
        mMasterFieldInterpolators.resize( mInterpDofTypeList.size(), nullptr );

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

    void Set::create_properties()
    {
        // cell of properties
        mMasterProperties.resize( mPropertyTypeList.size(), nullptr );

        // create a field interpolator for each dof type
        for( uint i = 0; i < mPropertyTypeList.size(); i++ )
        {
            // fixme create/get property active dof type list
            Cell< Cell< MSI::Dof_Type > >  tActiveDofTypes;

            // fixme create/get property val function
            std::function< Matrix< DDRMat > ( Cell< Matrix< DDRMat > >    & aCoeff,
                                              Cell< Field_Interpolator* > & aFieldInterpolator ) > tValFunction;

            // fixme create/get property derivatives
            Cell< std::function< Matrix< DDRMat > ( Cell< Matrix< DDRMat > >    & aCoeff,
                                                    Cell< Field_Interpolator* > & aFieldInterpolator ) > > tDerFunctions( 0 );

            // create an property for the ith property type
            mMasterProperties( i ) = new Property( mPropertyTypeList( i ), tActiveDofTypes, tValFunction, tDerFunctions );
        }
    }

//------------------------------------------------------------------------------

//    moris::Cell< Field_Interpolator* > Set::get_IWG_field_interpolators( IWG*                               & aIWG,
//                                                                         moris::Cell< Field_Interpolator* > & aFieldInterpolators )
//    {
//        // ask the IWG for its active dof types
//        moris::Cell< moris::Cell< MSI::Dof_Type > > tIWGActiveDof = aIWG->get_active_dof_types();
//
//        // number of active dof type for the IWG
//        uint tNumOfIWGActiveDof = tIWGActiveDof.size();
//
//        // select associated active interpolators
//        moris::Cell< Field_Interpolator* > tIWGFieldInterpolators( tNumOfIWGActiveDof, nullptr );
//        for( uint i = 0; i < tNumOfIWGActiveDof; i++ )
//        {
//            // find the index of active dof type in the list of element dof type
//            uint tIWGDofIndex = mInterpDofTypeMap( static_cast< int >( tIWGActiveDof( i )( 0 ) ) );
//
//            // select the corresponding interpolator
//            tIWGFieldInterpolators( i ) = aFieldInterpolators( tIWGDofIndex );
//        }
//        return tIWGFieldInterpolators;
//    }

    void Set::create_IWG_set_info()
    {
        // get number of IWGs
        mNumOfIWG = mIWGs.size();

        // set info size
        mIWGNumActiveDof.resize( mNumOfIWG );
        mIWGMasterFieldInterpolators.resize( mNumOfIWG );
        mIWGDofAssemblyMap.resize( mNumOfIWG );

        // loop over the IWGs
        for ( uint iIWG = 0; iIWG < mNumOfIWG; iIWG++ )
        {
            // get IWG
            IWG* tIWG = mIWGs( iIWG );

            // ask the IWG for its residual dof type
            moris::Cell< MSI::Dof_Type > tIWGResidualDof = tIWG->get_residual_dof_type();

            // ask the IWG for its active dof types
            moris::Cell< moris::Cell< MSI::Dof_Type > > tIWGActiveDof = tIWG->get_active_dof_types();

            // get the number of active dof type for the IWG
            mIWGNumActiveDof( iIWG ) = tIWGActiveDof.size();

            // find the index of residual dof type in the list of element dof type
            uint tIWGResDofIndex = mInterpDofTypeMap( static_cast< int >( tIWGResidualDof( 0 ) ), 0 );
            uint startIDof = mInterpDofAssemblyMap( tIWGResDofIndex, 0 );
            uint stopIDof  = mInterpDofAssemblyMap( tIWGResDofIndex, 1 );

            // set size
            mIWGMasterFieldInterpolators( iIWG ).resize( mIWGNumActiveDof( iIWG ), nullptr );
            mIWGDofAssemblyMap( iIWG ).set_size( mIWGNumActiveDof( iIWG ), 4 );

            // select associated active interpolators
            for( uint iDOF = 0; iDOF < mIWGNumActiveDof( iIWG ); iDOF++ )
            {
                // find the index of active dof type in the list of element dof type
                uint tIWGActiveDofIndex = mInterpDofTypeMap( static_cast< int >( tIWGActiveDof( iDOF )( 0 ) ), 0 );
                uint startJDof = mInterpDofAssemblyMap( tIWGActiveDofIndex, 0 );
                uint stopJDof  = mInterpDofAssemblyMap( tIWGActiveDofIndex, 1 );

                // select the corresponding interpolator
                mIWGMasterFieldInterpolators( iIWG )( iDOF ) = mMasterFieldInterpolators( tIWGActiveDofIndex );

                // build the dof assembly map for each IWG
                mIWGDofAssemblyMap( iIWG )( iDOF, 0 ) = startIDof;
                mIWGDofAssemblyMap( iIWG )( iDOF, 1 ) = stopIDof;
                mIWGDofAssemblyMap( iIWG )( iDOF, 2 ) = startJDof;
                mIWGDofAssemblyMap( iIWG )( iDOF, 3 ) = stopJDof;
            }

            // set IWG field interpolators
            tIWG->set_field_interpolators( mIWGMasterFieldInterpolators( iIWG ) );
        }
    }

    void Set::create_IWG_set_info_double()
        {
            // get number of IWGs
            mNumOfIWG = mIWGs.size();

            // set info size
            mIWGNumActiveDof.resize( mNumOfIWG );
            mIWGMasterFieldInterpolators.resize( mNumOfIWG );
            mIWGSlaveFieldInterpolators.resize( mNumOfIWG );
            mIWGDofAssemblyMap.resize( mNumOfIWG );

            // loop over the IWGs
            for ( uint iIWG = 0; iIWG < mNumOfIWG; iIWG++ )
            {
                // get IWG
                IWG* tIWG = mIWGs( iIWG );

                // ask the IWG for its residual dof type
                moris::Cell< MSI::Dof_Type > tIWGResidualDof = tIWG->get_residual_dof_type();

                // ask the IWG for its active dof types
                moris::Cell< moris::Cell< MSI::Dof_Type > > tIWGActiveDof = tIWG->get_active_dof_types();

                // get the number of active dof type for the IWG
                mIWGNumActiveDof( iIWG ) = tIWGActiveDof.size();

                // find the index of residual dof type in the list of element dof type
                uint tIWGResDofIndex = mInterpDofTypeMap( static_cast< int >( tIWGResidualDof( 0 ) ), 0 );
                uint startIDof = mInterpDofAssemblyMap( tIWGResDofIndex, 0 );
                uint stopIDof  = mInterpDofAssemblyMap( tIWGResDofIndex, 1 );

                // set size
                mIWGMasterFieldInterpolators( iIWG ).resize( mIWGNumActiveDof( iIWG ), nullptr );
                mIWGSlaveFieldInterpolators( iIWG ).resize( mIWGNumActiveDof( iIWG ), nullptr );
                mIWGDofAssemblyMap( iIWG ).set_size( mIWGNumActiveDof( iIWG ), 4 );

                // select associated active interpolators
                for( uint iDOF = 0; iDOF < mIWGNumActiveDof( iIWG ); iDOF++ )
                {
                    // find the index of active dof type in the list of element dof type
                    uint tIWGActiveDofIndex = mInterpDofTypeMap( static_cast< int >( tIWGActiveDof( iDOF )( 0 ) ), 0 );
                    uint startJDof = mInterpDofAssemblyMap( tIWGActiveDofIndex, 0 );
                    uint stopJDof  = mInterpDofAssemblyMap( tIWGActiveDofIndex, 1 );

                    // select the corresponding interpolator
                    mIWGMasterFieldInterpolators( iIWG )( iDOF ) = mMasterFieldInterpolators( tIWGActiveDofIndex );
                    mIWGSlaveFieldInterpolators( iIWG )( iDOF )  = mSlaveFieldInterpolators( tIWGActiveDofIndex );

                    // build the dof assembly map for each IWG
                    mIWGDofAssemblyMap( iIWG )( iDOF, 0 ) = startIDof;
                    mIWGDofAssemblyMap( iIWG )( iDOF, 1 ) = stopIDof;
                    mIWGDofAssemblyMap( iIWG )( iDOF, 2 ) = startJDof;
                    mIWGDofAssemblyMap( iIWG )( iDOF, 3 ) = stopJDof;
                }

                // set IWG field interpolators
                tIWG->set_field_interpolators( mIWGMasterFieldInterpolators( iIWG ),
                                               mtk::Master_Slave::MASTER );
                tIWG->set_field_interpolators( mIWGSlaveFieldInterpolators( iIWG ),
                                               mtk::Master_Slave::SLAVE );

            }
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
