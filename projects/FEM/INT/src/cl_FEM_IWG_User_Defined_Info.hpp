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
#include "cl_MTK_Enums.hpp"             //FEM/INT/src

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

            // IWGs residual dof type
            moris::Cell< moris::Cell< moris::Cell< MSI::Dof_Type > > > mResidualDofType;

            // IWGs master dof type dependency list
            moris::Cell< moris::Cell< moris::Cell< moris::Cell< MSI::Dof_Type > > > > mMasterDofTypes;

            // IWGs master property type dependency list
            moris::Cell< moris::Cell< moris::Cell< fem::Property_Type > > > mMasterPropTypes;

            // IWGs slave dof type dependency list
            moris::Cell< moris::Cell< moris::Cell< moris::Cell< MSI::Dof_Type > > > > mSlaveDofTypes;

             // IWGs slave property type dependency list
            moris::Cell< moris::Cell< moris::Cell< fem::Property_Type > > > mSlavePropTypes;

//------------------------------------------------------------------------------
        public :

//------------------------------------------------------------------------------
            /**
             * constructor
             */
            IWG_User_Defined_Info(){};

            IWG_User_Defined_Info( moris::Cell< moris::Cell< fem::IWG_Type > >                                    & aIWGTypeList,
                                   moris::Cell< moris::Cell< moris::Cell< MSI::Dof_Type > > >                     & aResidualDofType,
                                   moris::Cell< moris::Cell< moris::Cell< moris::Cell< MSI::Dof_Type > > > >      & aMasterDofTypes,
                                   moris::Cell< moris::Cell< moris::Cell< fem::Property_Type > > >                & aMasterPropTypes,
                                   moris::Cell< moris::Cell< moris::Cell< moris::Cell< MSI::Dof_Type > > > >      & aSlaveDofTypes,
                                   moris::Cell< moris::Cell< moris::Cell< fem::Property_Type > > >                & aSlavePropTypes )
                                 : mIWGTypeList( aIWGTypeList ),
                                   mResidualDofType( aResidualDofType ),
                                   mMasterDofTypes( aMasterDofTypes ),
                                   mMasterPropTypes( aMasterPropTypes ),
                                   mSlaveDofTypes( aSlaveDofTypes ),
                                   mSlavePropTypes( aSlavePropTypes )
            {};

            IWG_User_Defined_Info( moris::Cell< moris::Cell< fem::IWG_Type > >                                    & aIWGTypeList,
                                   moris::Cell< moris::Cell< moris::Cell< MSI::Dof_Type > > >                     & aResidualDofType,
                                   moris::Cell< moris::Cell< moris::Cell< moris::Cell< MSI::Dof_Type > > > >      & aMasterDofTypes,
                                   moris::Cell< moris::Cell< moris::Cell< fem::Property_Type > > >                & aMasterPropTypes )
                                 : mIWGTypeList( aIWGTypeList ),
                                   mResidualDofType( aResidualDofType ),
                                   mMasterDofTypes( aMasterDofTypes ),
                                   mMasterPropTypes( aMasterPropTypes )
            {};

//------------------------------------------------------------------------------
            /**
             * virtual destructor
             */
            ~IWG_User_Defined_Info(){};

//------------------------------------------------------------------------------
            /**
             * returns IWG type list
             */
           const moris::Cell< moris::Cell< fem::IWG_Type > > & get_IWG_type_list() const
            {
                return mIWGTypeList;
            };

//------------------------------------------------------------------------------
            /**
             * returns IWG residual dof type
             */
            const moris::Cell< moris::Cell< moris::Cell< MSI::Dof_Type > > > & get_residual_dof_type() const
            {
                return mResidualDofType;
            };

//------------------------------------------------------------------------------
            /**
             * returns IWG dof type list
             */
            const moris::Cell< moris::Cell< moris::Cell< moris::Cell< MSI::Dof_Type > > > > & get_dof_type_list( mtk::Master_Slave aIsMaster = mtk::Master_Slave::MASTER ) const
            {
                switch ( aIsMaster )
                {
                    case ( mtk::Master_Slave::MASTER ):
                    {
                        return mMasterDofTypes;
                    }
                    case ( mtk::Master_Slave::SLAVE ):
                    {
                        return mSlaveDofTypes;
                    }
                    default:
                    {
                        MORIS_ERROR( false, "IWG_User_Defined_Info::get_dof_type_list - can only be MASTER or SLAVE.");
                        return mMasterDofTypes;
                    }
                }
            };

//------------------------------------------------------------------------------
            /**
             * returns IWG property type list
             */
            const moris::Cell< moris::Cell< moris::Cell< fem::Property_Type > > > & get_property_type_list( mtk::Master_Slave aIsMaster = mtk::Master_Slave::MASTER ) const
            {
                switch ( aIsMaster )
                {
                    case ( mtk::Master_Slave::MASTER ):
                    {
                        return mMasterPropTypes;
                    }
                    case ( mtk::Master_Slave::SLAVE ):
                    {
                        return mSlavePropTypes;
                    }
                    default:
                    {
                        MORIS_ERROR( false, "IWG_User_Defined_Info::get_property_type_list - can only be MASTER or SLAVE.");
                        return mMasterPropTypes;
                    }
                }
            };

//------------------------------------------------------------------------------
        };
    }/* namespace fem */
} /* namespace moris */
#endif /* SRC_FEM_CL_FEM_IWG_USER_DEFINED_INFO_HPP_ */
