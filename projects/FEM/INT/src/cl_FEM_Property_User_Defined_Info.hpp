/*
 * cl_FEM_Property_User_Defined_Info.hpp
 *
 *  Created on: Aug 23, 2019
 *      Author: noel
 */
#ifndef SRC_FEM_CL_FEM_PROPERTY_USER_DEFINED_INFO_HPP_
#define SRC_FEM_CL_FEM_PROPERTY_USER_DEFINED_INFO_HPP_

#include "typedefs.hpp"                     //MRS/COR/src
#include "cl_Matrix.hpp"                    //LNA/src
#include "cl_Cell.hpp"                      //MRS/CON/src
#include "cl_MSI_Dof_Type_Enums.hpp"        //FEM/MSI/src
#include "cl_FEM_Enums.hpp"                 //FEM/MSI/src
#include <functional>

namespace moris
{
    namespace fem
    {
    class Field_Interpolator;
    typedef std::function< Matrix< DDRMat > ( moris::Cell< Matrix< DDRMat > >         & aCoeff,
                                              moris::Cell< fem::Field_Interpolator* > & aFieldInterpolator ) > PropertyFunc;
//------------------------------------------------------------------------------
        /**
         * Property_User_Defined_Info
         */
        class Property_User_Defined_Info
        {
        protected :

            // properties type list
            moris::Cell< fem::Property_Type > mPropertyTypeList;

            // properties dof type dependency list
            moris::Cell< moris::Cell< moris::Cell< MSI::Dof_Type > > > mPropertyDofList;

            // properties coeff list
            moris::Cell< moris::Cell< Matrix< DDRMat > > > mCoeffList;

            // properties value function list
            moris::Cell< PropertyFunc > mValFuncList;

            // properties derivative functions list
            moris::Cell< moris::Cell< PropertyFunc > > mDerFuncList;

            // property map
            Matrix< DDSMat > mPropertyTypeMap;

//------------------------------------------------------------------------------
        public :

//------------------------------------------------------------------------------
            /**
             * constructor
             */
            Property_User_Defined_Info(){};

            Property_User_Defined_Info( const moris::Cell< fem::Property_Type >                          & aPropertyTypeList,
                                        const moris::Cell< moris::Cell< moris::Cell< MSI::Dof_Type > > > & aPropertyDofList,
                                        const moris::Cell< moris::Cell< Matrix< DDRMat > > >             & aCoeffList,
                                        const moris::Cell< PropertyFunc >                                & aValFuncList,
                                        const moris::Cell< moris::Cell< PropertyFunc > >                 & aDerFuncList )
            : mPropertyTypeList( aPropertyTypeList ),
              mPropertyDofList( aPropertyDofList ),
              mCoeffList( aCoeffList ),
              mValFuncList( aValFuncList ),
              mDerFuncList( aDerFuncList )
            {
                // build the property map
                this->build_property_map();
            };

//------------------------------------------------------------------------------
            /**
             * virtual destructor
             */
             ~Property_User_Defined_Info(){};

//------------------------------------------------------------------------------
            /**
             * returns property type
             */
            const moris::Cell< fem::Property_Type > & get_property_type_list() const
            {
                return mPropertyTypeList;
            };

//------------------------------------------------------------------------------
            /**
             * returns property dof type dependency list
             */
            const moris::Cell< moris::Cell< moris::Cell< MSI::Dof_Type > > > & get_property_dof_type_list() const
            {
                return mPropertyDofList;
            };

//------------------------------------------------------------------------------
            /**
             * returns property coeff list
             */
            const moris::Cell< moris::Cell< Matrix< DDRMat > > > & get_property_coeff_list() const
            {
                return mCoeffList;
            };
//------------------------------------------------------------------------------
            /**
             * returns property value function list
             */
            const moris::Cell< PropertyFunc > & get_property_valFunc_list() const
            {
                return mValFuncList;
            };

//------------------------------------------------------------------------------
            /**
             * returns property derivative function list
             */
            const moris::Cell< moris::Cell< PropertyFunc > > & get_property_derFunc_list() const
            {
                return mDerFuncList;
            };

//------------------------------------------------------------------------------
            /**
             * build a property map
             */
            void build_property_map()
            {
                // get the max enum for the properties to create a property map
                sint tMaxEnum = 0;

                // loop over the property types
                for( uint iProp = 0; iProp < mPropertyTypeList.size(); iProp++ )
                {
                    // check if enum is greater
                    tMaxEnum = std::max( tMaxEnum, static_cast< int >( mPropertyTypeList( iProp ) ) );
                }
                // +1 because 0 based
                tMaxEnum++;

                // set the size of the property map
                mPropertyTypeMap.set_size( tMaxEnum, 1, -1 );

                // loop over the property types
                for( uint iProp = 0; iProp < mPropertyTypeList.size(); iProp++ )
                {
                    // fill the property map
                    mPropertyTypeMap( static_cast< int >( mPropertyTypeList( iProp ) ), 0 ) = iProp;
                }
            };

//------------------------------------------------------------------------------
            /**
             * get a property map
             */
            const Matrix< DDSMat > & get_property_map() const
            {
                return mPropertyTypeMap;
            }
        };

//------------------------------------------------------------------------------
    } /* namespace fem */
} /* namespace moris */
#endif /* SRC_FEM_CL_FEM_PROPERTY_USER_DEFINED_INFO_HPP_ */
