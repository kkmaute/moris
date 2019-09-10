/*
 * cl_FEM_IWG_User_Defined_Info.hpp
 *
 *  Created on: Sep 06, 2019
 *      Author: noel
 */
#ifndef SRC_FEM_CL_FEM_IWG_USER_DEFINED_INFO_HPP_
#define SRC_FEM_CL_FEM_IWG_USER_DEFINED_INFO_HPP_

#include "typedefs.hpp"                     //MRS/COR/src
#include "cl_Matrix.hpp"                    //LNA/src
#include "cl_Cell.hpp"                      //MRS/CON/src
#include "cl_MSI_Dof_Type_Enums.hpp"        //FEM/MSI/src
#include "cl_FEM_Enums.hpp"                 //FEM/MSI/src

namespace moris
{
    namespace fem
    {
//------------------------------------------------------------------------------
        /**
         * IWG_User_Defined_Info
         */
        class IWG_User_Defined_Info
        {
        protected :

            // IWGs type list
            moris::Cell< moris::Cell< fem::IWG_Type > > mIWGTypeList;

            // IWGs master dof type dependency list
            moris::Cell< moris::Cell< moris::Cell< MSI::Dof_Type > > > mMasterDofTypes;
            moris::Cell< moris::Cell< moris::Cell< MSI::Dof_Type > > > mSlaveDofTypes;

            // IWGs master property type dependency list
            moris::Cell< moris::Cell< MSI::Property_Type > > mMasterPropTypes;
            moris::Cell< moris::Cell< MSI::Property_Type > > mSlavePropTypes;

//------------------------------------------------------------------------------
        public :

//------------------------------------------------------------------------------
            /**
             * constructor
             */
            IWG_User_Defined_Info(){};

            IWG_User_Defined_Info( const moris::Cell< moris::Cell< fem::IWG_Type > >                     & aIWGTypeList,
                                   const moris::Cell< moris::Cell< moris::Cell< MSI::Dof_Type > > >      & aMasterDofTypes,
                                   const moris::Cell< moris::Cell< moris::Cell< MSI::Dof_Type > > >      & aSlaveDofTypes,
                                   const moris::Cell< moris::Cell< moris::Cell< MSI::Property_Type > > > & aMasterPropTypes,
                                   const moris::Cell< moris::Cell< moris::Cell< MSI::Property_Type > > > & aSlavePropTypes)
                                 : mIWGTypeList( aIWGTypeList ),
                                   mMasterDofTypes( aMasterDofTypes ),
                                   mSlaveDofTypes( aSlaveDofTypes ),
                                   mMasterPropTypes( aMasterPropTypes ),
                                   mSlavePropTypes( aSlavePropTypes )
            {};

//------------------------------------------------------------------------------
            /**
             * virtual destructor
             */
            ~IWG_User_Defined_Info(){};

//------------------------------------------------------------------------------
    } /* namespace fem */
} /* namespace moris */
#endif /* SRC_FEM_CL_FEM_IWG_USER_DEFINED_INFO_HPP_ */
