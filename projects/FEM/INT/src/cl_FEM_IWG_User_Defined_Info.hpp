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
#include "cl_MTK_Enums.hpp"                 //FEM/INT/src

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

            // IWG type
            fem::IWG_Type mIWGType;

            // IWG space dimension
            uint mSpaceDim;

            // IWG residual dof type
            moris::Cell< MSI::Dof_Type > mResidualDofType;

            // IWG master dof type dependency list
            moris::Cell< moris::Cell< MSI::Dof_Type > > mMasterDofTypes;

            // IWG master property type dependency list
            moris::Cell< fem::Property_Type > mMasterPropTypes;

            // IWG master constitutive type dependency list
            moris::Cell< fem::Constitutive_Type > mMasterConstitutiveTypes;

            // IWG slave dof type dependency list
            moris::Cell< moris::Cell< MSI::Dof_Type > > mSlaveDofTypes;

            // IWG slave property type dependency list
            moris::Cell< fem::Property_Type > mSlavePropTypes;

            // IWG slave constitutive type dependency list
            moris::Cell< fem::Constitutive_Type > mSlaveConstitutiveTypes;

//------------------------------------------------------------------------------
        public :

//------------------------------------------------------------------------------
            /**
             * constructor
             */
            IWG_User_Defined_Info(){};

            IWG_User_Defined_Info( fem::IWG_Type                               aIWGType,
                                   uint                                        aSpaceDim,
                                   moris::Cell< MSI::Dof_Type >                aResidualDofType,
                                   moris::Cell< moris::Cell< MSI::Dof_Type > > aMasterDofTypes,
                                   moris::Cell< fem::Property_Type >           aMasterPropTypes,
                                   moris::Cell< fem::Constitutive_Type >       aMasterConstitutiveTypes,
                                   moris::Cell< moris::Cell< MSI::Dof_Type > > aSlaveDofTypes,
                                   moris::Cell< fem::Property_Type >           aSlavePropTypes,
                                   moris::Cell< fem::Constitutive_Type >       aSlaveConstitutiveTypes )
                                 : mIWGType( aIWGType ),
                                   mSpaceDim( aSpaceDim ),
                                   mResidualDofType( aResidualDofType ),
                                   mMasterDofTypes( aMasterDofTypes ),
                                   mMasterPropTypes( aMasterPropTypes ),
                                   mMasterConstitutiveTypes( aMasterConstitutiveTypes ),
                                   mSlaveDofTypes( aSlaveDofTypes ),
                                   mSlavePropTypes( aSlavePropTypes ),
                                   mSlaveConstitutiveTypes( aSlaveConstitutiveTypes )
            {};

            IWG_User_Defined_Info( fem::IWG_Type                               aIWGType,
                                   uint                                        aSpaceDim,
                                   moris::Cell< MSI::Dof_Type >                aResidualDofType,
                                   moris::Cell< moris::Cell< MSI::Dof_Type > > aMasterDofTypes,
                                   moris::Cell< fem::Property_Type >           aMasterPropTypes,
                                   moris::Cell< fem::Constitutive_Type >       aMasterConstitutiveTypes )
                                 : mIWGType( aIWGType ),
                                   mSpaceDim( aSpaceDim ),
                                   mResidualDofType( aResidualDofType ),
                                   mMasterDofTypes( aMasterDofTypes ),
                                   mMasterPropTypes( aMasterPropTypes ),
                                   mMasterConstitutiveTypes( aMasterConstitutiveTypes )
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
           fem::IWG_Type get_IWG_type() const
            {
                return mIWGType;
            };

//------------------------------------------------------------------------------
           /**
            * returns IWG space dimension
            */
           uint get_IWG_space_dim() const
           {
               return mSpaceDim;
           };

//------------------------------------------------------------------------------
            /**
             * returns IWG residual dof type
             */
            const moris::Cell< MSI::Dof_Type > & get_residual_dof_type() const
            {
                return mResidualDofType;
            };

//------------------------------------------------------------------------------
            /**
             * returns IWG dof type list
             */
            const moris::Cell< moris::Cell< MSI::Dof_Type > > & get_dof_type_list( mtk::Master_Slave aIsMaster = mtk::Master_Slave::MASTER ) const
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
            const moris::Cell< fem::Property_Type > & get_property_type_list( mtk::Master_Slave aIsMaster = mtk::Master_Slave::MASTER ) const
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
            /**
             * returns IWG constitutive type list
             */
            const moris::Cell< fem::Constitutive_Type > & get_constitutive_type_list( mtk::Master_Slave aIsMaster = mtk::Master_Slave::MASTER ) const
            {
                switch ( aIsMaster )
                {
                    case ( mtk::Master_Slave::MASTER ):
                    {
                        return mMasterConstitutiveTypes;
                    }
                    case ( mtk::Master_Slave::SLAVE ):
                    {
                        return mSlaveConstitutiveTypes;
                    }
                    default:
                    {
                        MORIS_ERROR( false, "IWG_User_Defined_Info::get_constitutive_type_list - can only be MASTER or SLAVE.");
                        return mMasterConstitutiveTypes;
                    }
                }
            };

//------------------------------------------------------------------------------
        };
    }/* namespace fem */
} /* namespace moris */
#endif /* SRC_FEM_CL_FEM_IWG_USER_DEFINED_INFO_HPP_ */
