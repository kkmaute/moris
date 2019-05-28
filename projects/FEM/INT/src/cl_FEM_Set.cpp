/*
 * cl_FEM_Set.cpp
 *
 *  Created on: Apr 11, 2019
 *      Author: schmidt
 */
#include <iostream>

#include "cl_FEM_Set.hpp"
#include "cl_FEM_Element_Factory.hpp"
#include "cl_MSI_Model_Solver_Interface.hpp"
#include "cl_FEM_Integrator.hpp"   //FEM/INT/src

namespace moris
{
    namespace fem
    {
//------------------------------------------------------------------------------
    Set::Set( moris::Cell< mtk::Cell const * > & aCell,
                                 enum fem::Element_Type      aElementType,
                                 Cell< IWG* >              & aIWGs,
                                 Cell< Node_Base* >        & aNodes) : mMeshElementPointer(aCell),
                                                                       mNodes(aNodes),
                                                                       mIWGs( aIWGs ),
                                                                       mElementType(aElementType)
    {

//    	this->create_unique_dof_type_lists();
//    	this->create_unique_list_of_first_dof_type_of_group();
//    	this->create_dof_type_lists();

        this->create_unique_dof_type_lists();
        this->create_dof_type_lists();


        mEquationObjList.resize( mMeshElementPointer.size(), nullptr);

    	// a factory to create the elements
        fem::Element_Factory tElementFactory;

        for( luint k=0; k < mMeshElementPointer.size(); ++k )
        {

        	// create the element
            //mElements( k ) = tElementFactory.create_element( mElementType,

            // create the element // FIXME replace with mtk::cluster information
            mEquationObjList( k ) = tElementFactory.create_cluster( mElementType,
                                                             mMeshElementPointer( k ),
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
    }

//------------------------------------------------------------------------------

    void Set::delete_pointers()
    {
        // delete the geometry interpolator pointer
        if ( mGeometryInterpolator != nullptr )
        {
            delete mGeometryInterpolator;
        }

        for( auto tFieldInterpolator : mFieldInterpolators )
        {
            delete tFieldInterpolator;
        }
    }

//------------------------------------------------------------------------------

    void Set::finalize( MSI::Model_Solver_Interface * aModelSolverInterface )
    {
        this->delete_pointers();

        if( mMeshElementPointer.size() > 0)
        {
             Interpolation_Rule tGeometryInterpolationRule( mMeshElementPointer( 0 )->get_geometry_type(),       // FIXME change to block information
                                                            Interpolation_Type::LAGRANGE,
                                                            this->get_auto_interpolation_order( mMeshElementPointer( 0 )->get_number_of_vertices(),
                                                                                                mMeshElementPointer( 0 )->get_geometry_type() ),           // FIXME change to block information
                                                            Interpolation_Type::LAGRANGE,
                                                            mtk::Interpolation_Order::LINEAR );

             bool tSpaceSideset = false;
             if (mElementType==fem::Element_Type::SIDESET)
             {
                 tSpaceSideset=true;
             }

             // create the element geometry intepolator
             mGeometryInterpolator = new Geometry_Interpolator( tGeometryInterpolationRule, tSpaceSideset );

            // create the element field interpolators
            this->create_field_interpolators( aModelSolverInterface );

            // create the element dof assembly map
            this->create_dof_assembly_map();

            Integration_Rule* tIntegrationRule;

            if (mElementType==fem::Element_Type::SIDESET)
            {
                mtk::Geometry_Type tSideGeometryType = this->get_block_geometry_interpolator()->get_side_geometry_type();
                enum fem::Integration_Order tSideIntegrationOrder = this->get_auto_integration_order( tSideGeometryType );

                tIntegrationRule = new Integration_Rule( tSideGeometryType,
                                                         Integration_Type::GAUSS,
                                                         tSideIntegrationOrder,
                                                         Integration_Type::GAUSS,
                                                         Integration_Order::BAR_1 );
            }
            else
            {
                enum fem::Integration_Order tIntegrationOrder = this->get_auto_integration_order( mMeshElementPointer( 0 )->get_geometry_type() );

                tIntegrationRule = new Integration_Rule( mMeshElementPointer( 0 )->get_geometry_type(),
                                                         Integration_Type::GAUSS,
                                                         tIntegrationOrder,
                                                         Integration_Type::GAUSS,
                                                         Integration_Order::BAR_1 );
            }

            // create an integrator for the ith IWG
            Integrator tIntegrator( *tIntegrationRule );

            delete tIntegrationRule;

            //get number of integration points
            mNumOfIntegPoints = tIntegrator.get_number_of_points();

            // get integration points
            mSurfRefIntegPoints = tIntegrator.get_points();

            // get integration weights
            mIntegWeights = tIntegrator.get_weights();
        }
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
        Cell<moris::moris_index> tUniqueDofTypeGroupsIndices = unique_index( tInterpDofTypeListBuild );

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
            tDofCounter = tDofCounter + mFieldInterpolators( i )->get_number_of_space_time_coefficients()-1;

            // fill the assembly map with starting dof counter
            mInterpDofAssemblyMap( i, 1 ) = tDofCounter;

            // update dof counter
            tDofCounter = tDofCounter + 1;
        }

        // set mTotalDof
        mTotalDof = tDofCounter;
    }

//------------------------------------------------------------------------------

    mtk::Interpolation_Order Set::get_auto_interpolation_order( const moris::uint aNumVertices,
                                                                          const mtk::Geometry_Type aGeometryType )
    {
        switch( aGeometryType )                                 // FIXME change to block information
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

    void Set::create_field_interpolators(MSI::Model_Solver_Interface * aModelSolverInterface )
    {
        // cell of field interpolators
        mFieldInterpolators.resize( mNumOfInterp, nullptr );

        // loop on the dof type groups and create a field interpolator for each
        for( uint i = 0; i < mNumOfInterp; i++ )
        {
            // get the ith dof type group
            Cell< MSI::Dof_Type > tDofTypeGroup = mInterpDofTypeList( i );

            moris::uint tNumTimeNodes = aModelSolverInterface->get_time_levels_for_type( tDofTypeGroup( 0 ) );

            // create the field interpolation rule for the ith dof type group
            //FIXME: space interpolation based on the mtk::Cell
            //FIXME: time  interpolation set to constant
            Interpolation_Rule tFieldInterpolationRule( mMeshElementPointer( 0 )->get_geometry_type(),           //FIXME
                                                        Interpolation_Type::LAGRANGE,
                                                        this->get_auto_interpolation_order( mMeshElementPointer( 0 )->get_number_of_vertices(),
                                                                                            mMeshElementPointer( 0 )->get_geometry_type()),
                                                        this->get_auto_time_interpolation_type( tNumTimeNodes ),
                                                        // If interpolation type CONSTANT, iInterpolation order is not used
                                                        this->get_auto_interpolation_order( tNumTimeNodes,
                                                                                            mtk::Geometry_Type::LINE ) );

            // get number of field interpolated by the ith field interpolator
            uint tNumOfFields = tDofTypeGroup.size();

            // create an interpolator for the ith dof type group
            mFieldInterpolators( i ) = new Field_Interpolator( tNumOfFields,
                                                               tFieldInterpolationRule,
                                                               mGeometryInterpolator );
        }
    }

//------------------------------------------------------------------------------

    moris::Cell< Field_Interpolator* > Set::get_IWG_field_interpolators ( IWG*                               & aIWG,
                                                                     moris::Cell< Field_Interpolator* > & aFieldInterpolators )
    {
        // ask the IWG for its active dof types
        Cell< Cell< MSI::Dof_Type > > tIWGActiveDof = aIWG->get_active_dof_types();

        // number of active dof type for the IWG
        uint tNumOfIWGActiveDof = tIWGActiveDof.size();

        // select associated active interpolators
        Cell< Field_Interpolator* > tIWGFieldInterpolators( tNumOfIWGActiveDof, nullptr );
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

    void Set::initialize_mJacobianElement()
    {
        if ( !mJacobianExist )
        {
            uint tTotalDof = this->get_total_number_of_dofs();
            mJacobian.set_size( tTotalDof, tTotalDof, 0.0 );

            mJacobianExist = true;
        }
        else
        {
            MORIS_ASSERT( mJacobian.numel() > 0, "Set::initialize_mJacobianElement(): Jacobian not properly initialized.");

            mJacobian.fill( 0.0 );
        }
    }

//------------------------------------------------------------------------------

    void Set::initialize_mResidualElement()
    {
        if ( !mResidualExist )
        {
            mResidual.set_size( this->get_total_number_of_dofs(), 1, 0.0 );

            mResidualExist = true;
        }
        else
        {
            MORIS_ASSERT( mResidual.numel() > 0, "Set::initialize_mJacobianElement(): Residual not properly initialized.");

            mResidual.fill( 0.0 );
        }
    }

//------------------------------------------------------------------------------


    } /* namespace fem */
} /* namespace moris */
