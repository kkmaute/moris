/*
 * cl_FEM_Element_Block.hpp
 *
 *  Created on: Mar 10, 2019
 *      Author: schmidt
 */

#ifndef SRC_FEM_CL_FEM_ELEMENT_BLOCK_HPP_
#define SRC_FEM_CL_FEM_ELEMENT_BLOCK_HPP_

#include "assert.h"
#include "cl_MSI_Equation_Object.hpp"               //FEM/INT/src
#include "cl_FEM_Element_Factory.hpp"               //FEM/INT/src

#include "cl_Communication_Tools.hpp"               //FEM/INT/src

namespace moris
{
    namespace fem
    {

//------------------------------------------------------------------------------
    /**
     * \brief element block class that communicates with the mesh interface
     */
    class Element_Block
    {
    private:
        moris::Cell< mtk::Cell* >     mElementPointer;

        Cell< MSI::Equation_Object* > mElements;

        Geometry_Interpolator       * mGeometryInterpolator;

        moris::Cell< Field_Interpolator* >   mFieldInterpolators;

        // cell of pointers to IWG objects
        moris::Cell< IWG* > mIWGs;

        // map of the element active dof types
        moris::Cell< enum MSI::Dof_Type >    mEqnObjDofTypeList; // List of dof types of this equation obj
        moris::Cell< moris::Cell< DDRMat > > mElemDofTypeList;
        moris::Matrix< DDSMat >              mElemDofTypeMap;
        uint                                 mNumOfElemDofTypes;
        uint                                 mNumOfIWGs;
        moris::Matrix< DDSMat >              mInterpDofTypeMap;
        moris::Cell< Cell< MSI::Dof_Type > > mInterpDofTypeList;
        uint                                 mNumOfInterp;

//------------------------------------------------------------------------------
    public:
//------------------------------------------------------------------------------
        /**
         * constructor
         *
         * @param[ in ]     List of mtk::Cell pointer
         */
        Element_Block( moris::Cell< mtk::Cell* > & aCell,
                       Element_Type                aElementType,
                       Cell< IWG* >              & aIWGs,
                       Cell< Node_Base* >        & aNodes) : mElementPointer(aCell),
                                                             mIWGs( aIWGs )
        {
            if(mElementPointer.size() > 0)
            {
                this->create_dof_type_lists();

                Interpolation_Rule tGeometryInterpolationRule( mElementPointer( 0 )->get_geometry_type(),       // FIXME change to block information
                                                                Interpolation_Type::LAGRANGE,
                                                                this->get_auto_interpolation_order(),           // FIXME change to block information
                                                                Interpolation_Type::LAGRANGE,
                                                                mtk::Interpolation_Order::LINEAR );
\
                bool tSpaceSideset = false;

                if (aElementType==fem::Element_Type::SIDESET)
                {
                    tSpaceSideset=true;
                }

                // create the element geometry intepolator
                mGeometryInterpolator = new Geometry_Interpolator( tGeometryInterpolationRule, tSpaceSideset );

                // create the element field interpolators
                mFieldInterpolators = this->create_field_interpolators( mGeometryInterpolator );

                mElements.resize( mElementPointer.size(), nullptr);

                // a factory to create the elements
                fem::Element_Factory tElementFactory;

                for( luint k=0; k < mElementPointer.size(); ++k )
                {
                    // create the element
                    mElements( k )
                        = tElementFactory.create_element(   aElementType,
                                                            mElementPointer( k ),
                                                            aIWGs,
                                                            aNodes,
                                                            this );
                }
            }
        };

//------------------------------------------------------------------------------

        /**
         * trivial destructor
         */
        ~Element_Block(){};

//------------------------------------------------------------------------------

        void create_dof_type_lists()
        {
            // get the number of IWGs
            mNumOfIWGs = mIWGs.size();

            // set the size of the element active dof type list
            uint tCounter = 0;
            for ( uint i = 0; i < mNumOfIWGs; i++ )
            {
                tCounter = tCounter + mIWGs( i )->get_residual_dof_type().size();
            }
            mEqnObjDofTypeList.resize( tCounter );

            // loop over the IWGs
            tCounter = 0;
            for ( uint i = 0; i < mNumOfIWGs; i++ )
            {
                // get the residual dof type of the ith IWG
                Cell< MSI::Dof_Type > tDofType = mIWGs( i )->get_residual_dof_type();

                for ( uint j = 0; j < tDofType.size(); j++ )
               {
                   // get the residual dof type of the ith IWG
                   mEqnObjDofTypeList( tCounter ) = tDofType( j );
                   tCounter++;
                }
            }

            // use std::unique and std::distance to create a unique list containing all used dof types
            auto last = std::unique( ( mEqnObjDofTypeList.data() ).data(),
                                     ( mEqnObjDofTypeList.data() ).data() + mEqnObjDofTypeList.size() );
            auto pos  = std::distance( ( mEqnObjDofTypeList.data() ).data(), last );
            mEqnObjDofTypeList.resize( pos );

            // create a list of the groups of dof types provided by the IWGs----------------
            // FIXME works as long as the dof type are always grouped in the same way
            moris::Cell< MSI::Dof_Type > tInterpDofTypeListBuild( mNumOfIWGs );

            // loop over the IWGs
            for ( uint i = 0; i < mNumOfIWGs; i++ )
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

        Cell< MSI::Equation_Object* > & get_equation_object_list()
        {
            return mElements;
        };

//------------------------------------------------------------------------------

        moris::Cell< Field_Interpolator* > & get_block_field_interpolator()
        {
            return mFieldInterpolators;
        }

//------------------------------------------------------------------------------

        Geometry_Interpolator * get_block_geometry_interpolator()
        {
            return mGeometryInterpolator;
        }

//------------------------------------------------------------------------------

        /**
         * auto detect full integration scheme
         */
        //FIXME: works for Lagrange only
        mtk::Interpolation_Order get_auto_interpolation_order()
        {
            switch( mElementPointer( 0 )->get_geometry_type() )                                 // FIXME change to block information
            {
                case( mtk::Geometry_Type::LINE ) :
                    switch( mElementPointer( 0 )->get_number_of_vertices() )
                    {
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
                    switch( mElementPointer( 0 )->get_number_of_vertices() )
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
                    switch( mElementPointer( 0 )->get_number_of_vertices() )
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

        Cell< Field_Interpolator* > create_field_interpolators( Geometry_Interpolator* aGeometryInterpolator )
         {
             // cell of field interpolators
             Cell< Field_Interpolator* > tFieldInterpolators( mNumOfInterp, nullptr );

             // loop on the dof type groups and create a field interpolator for each
             for( uint i = 0; i < mNumOfInterp; i++ )
             {
                 // get the ith dof type group
                 Cell< MSI::Dof_Type > tDofTypeGroup = mInterpDofTypeList( i );

                 // create the field interpolation rule for the ith dof type group
                 //FIXME: space interpolation based on the mtk::Cell
                 //FIXME: time  interpolation set to constant
                 Interpolation_Rule tFieldInterpolationRule( mElementPointer( 0 )->get_geometry_type(),           //FIXME
                                                             Interpolation_Type::LAGRANGE,
                                                             this->get_auto_interpolation_order(),
                                                             Interpolation_Type::CONSTANT,
                                                             mtk::Interpolation_Order::CONSTANT );

                 // get number of field interpolated by the ith field interpolator
                 uint tNumOfFields = tDofTypeGroup.size();

                 // create an interpolator for the ith dof type group
                 tFieldInterpolators( i ) = new Field_Interpolator( tNumOfFields,
                                                                    tFieldInterpolationRule,
                                                                    aGeometryInterpolator );
             }
             return tFieldInterpolators;
         }

    };

//------------------------------------------------------------------------------
    } /* namespace fem */
} /* namespace moris */

#endif /* SRC_FEM_CL_FEM_ELEMENT_BLOCK_HPP_ */
