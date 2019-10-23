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

            // IWG residual dof type
            moris::Cell< MSI::Dof_Type > mResidualDofType;

            // IWG master dof type dependency list
            moris::Cell< moris::Cell< MSI::Dof_Type > > mMasterDofTypes;

            // IWG master dv type dependency list
            moris::Cell< moris::Cell< MSI::Dv_Type > > mMasterDvTypes;

            // IWG master property type dependency list
            moris::Cell< fem::Property_Type > mMasterPropTypes;

            // IWG master constitutive type dependency list
            moris::Cell< fem::Constitutive_Type > mMasterConstitutiveTypes;

            // IWG slave dof type dependency list
            moris::Cell< moris::Cell< MSI::Dof_Type > > mSlaveDofTypes;

            // IWG slave dv type dependency list
            moris::Cell< moris::Cell< MSI::Dv_Type > > mSlaveDvTypes;

            // IWG slave property type dependency list
            moris::Cell< fem::Property_Type > mSlavePropTypes;

            // IWG slave constitutive type dependency list
            moris::Cell< fem::Constitutive_Type > mSlaveConstitutiveTypes;

//------------------------------------------------------------------------------
        public :

//------------------------------------------------------------------------------
            /**
             * trivial constructor
             */
            IWG_User_Defined_Info(){};

            /**
             * constructor with master  and slave without dv
             * @param[ in ] aIWGType                 IWG type
             * @param[ in ] aResidualDofType         group of dof type for residual
             * @param[ in ] aMasterDofTypes          list of group of dof types for master
             * @param[ in ] aMasterPropTypes         list of property types for master
             * @param[ in ] aMasterConstitutiveTypes list of constitutive types for master
             * @param[ in ] aSlaveDofTypes           list of group of dof types for slave
             * @param[ in ] aSlavePropTypes          list of property types for slave
             * @param[ in ] aSlaveConstitutiveTypes  list of constitutive types for slave
             */
            IWG_User_Defined_Info( fem::IWG_Type                               aIWGType,
                                   moris::Cell< MSI::Dof_Type >                aResidualDofType,
                                   moris::Cell< moris::Cell< MSI::Dof_Type > > aMasterDofTypes,
                                   moris::Cell< fem::Property_Type >           aMasterPropTypes,
                                   moris::Cell< fem::Constitutive_Type >       aMasterConstitutiveTypes,
                                   moris::Cell< moris::Cell< MSI::Dof_Type > > aSlaveDofTypes,
                                   moris::Cell< fem::Property_Type >           aSlavePropTypes,
                                   moris::Cell< fem::Constitutive_Type >       aSlaveConstitutiveTypes )
                                 : mIWGType( aIWGType ),
                                   mResidualDofType( aResidualDofType ),
                                   mMasterDofTypes( aMasterDofTypes ),
                                   mMasterPropTypes( aMasterPropTypes ),
                                   mMasterConstitutiveTypes( aMasterConstitutiveTypes ),
                                   mSlaveDofTypes( aSlaveDofTypes ),
                                   mSlavePropTypes( aSlavePropTypes ),
                                   mSlaveConstitutiveTypes( aSlaveConstitutiveTypes )
            {};

            /**
             * constructor with master only without dv
             * @param[ in ] aIWGType                 IWG type
             * @param[ in ] aResidualDofType         group of dof type for residual
             * @param[ in ] aMasterDofTypes          list of group of dof types for master
             * @param[ in ] aMasterPropTypes         list of property types for master
             * @param[ in ] aMasterConstitutiveTypes list of constitutive types for master
             */
            IWG_User_Defined_Info( fem::IWG_Type                               aIWGType,
                                   moris::Cell< MSI::Dof_Type >                aResidualDofType,
                                   moris::Cell< moris::Cell< MSI::Dof_Type > > aMasterDofTypes,
                                   moris::Cell< fem::Property_Type >           aMasterPropTypes,
                                   moris::Cell< fem::Constitutive_Type >       aMasterConstitutiveTypes )
                                 : mIWGType( aIWGType ),
                                   mResidualDofType( aResidualDofType ),
                                   mMasterDofTypes( aMasterDofTypes ),
                                   mMasterPropTypes( aMasterPropTypes ),
                                   mMasterConstitutiveTypes( aMasterConstitutiveTypes )
            {};

            /**
             * constructor with master and slave
             * @param[ in ] aIWGType                 IWG type
             * @param[ in ] aResidualDofType         group of dof type for residual
             * @param[ in ] aMasterDofTypes          list of group of dof types for master
             * @param[ in ] aMasterDvTypes           list of group of dv types for master
             * @param[ in ] aMasterPropTypes         list of property types for master
             * @param[ in ] aMasterConstitutiveTypes list of constitutive types for master
             * @param[ in ] aSlaveDofTypes           list of group of dof types for slave
             * @param[ in ] aSlaveDvTypes            list of group of dv types for slave
             * @param[ in ] aSlavePropTypes          list of property types for slave
             * @param[ in ] aSlaveConstitutiveTypes  list of constitutive types for slave
             */
            IWG_User_Defined_Info( fem::IWG_Type                               aIWGType,
                                   moris::Cell< MSI::Dof_Type >                aResidualDofType,
                                   moris::Cell< moris::Cell< MSI::Dof_Type > > aMasterDofTypes,
                                   moris::Cell< moris::Cell< MSI::Dv_Type > >  aMasterDvTypes,
                                   moris::Cell< fem::Property_Type >           aMasterPropTypes,
                                   moris::Cell< fem::Constitutive_Type >       aMasterConstitutiveTypes,
                                   moris::Cell< moris::Cell< MSI::Dof_Type > > aSlaveDofTypes,
                                   moris::Cell< moris::Cell< MSI::Dv_Type > >  aSlaveDvTypes,
                                   moris::Cell< fem::Property_Type >           aSlavePropTypes,
                                   moris::Cell< fem::Constitutive_Type >       aSlaveConstitutiveTypes )
                                 : mIWGType( aIWGType ),
                                   mResidualDofType( aResidualDofType ),
                                   mMasterDofTypes( aMasterDofTypes ),
                                   mMasterDvTypes( aMasterDvTypes ),
                                   mMasterPropTypes( aMasterPropTypes ),
                                   mMasterConstitutiveTypes( aMasterConstitutiveTypes ),
                                   mSlaveDofTypes( aSlaveDofTypes ),
                                   mSlaveDvTypes( aSlaveDvTypes ),
                                   mSlavePropTypes( aSlavePropTypes ),
                                   mSlaveConstitutiveTypes( aSlaveConstitutiveTypes )
            {};

            /**
             * constructor with master only
             * @param[ in ] aIWGType                 IWG type
             * @param[ in ] aResidualDofType         group of dof type for residual
             * @param[ in ] aMasterDofTypes          list of group of dof types for master
             * @param[ in ] aMasterDvTypes           list of group of dv types for master
             * @param[ in ] aMasterPropTypes         list of property types for master
             * @param[ in ] aMasterConstitutiveTypes list of constitutive types for master
             */
            IWG_User_Defined_Info( fem::IWG_Type                               aIWGType,
                                   moris::Cell< MSI::Dof_Type >                aResidualDofType,
                                   moris::Cell< moris::Cell< MSI::Dof_Type > > aMasterDofTypes,
                                   moris::Cell< moris::Cell< MSI::Dv_Type > >  aMasterDvTypes,
                                   moris::Cell< fem::Property_Type >           aMasterPropTypes,
                                   moris::Cell< fem::Constitutive_Type >       aMasterConstitutiveTypes )
                                 : mIWGType( aIWGType ),
                                   mResidualDofType( aResidualDofType ),
                                   mMasterDofTypes( aMasterDofTypes ),
                                   mMasterDvTypes( aMasterDvTypes ),
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
             * set IWG type
             * @param[ in ] aIWGType IWG type
             */
           void set_IWG_type( fem::IWG_Type aIWGType )
            {
                mIWGType = aIWGType;
            };

//------------------------------------------------------------------------------
            /**
             * return IWG type
             * @param[ out ] mIWGType IWG type
             */
            fem::IWG_Type get_IWG_type() const
            {
                return mIWGType;
            };

//------------------------------------------------------------------------------
            /**
             * set IWG residual dof type
             * @param[ in ] aResidualDofType group of dof types for residual
             */
            void set_residual_dof_type( moris::Cell< MSI::Dof_Type > & aResidualDofType )
            {
                mResidualDofType = aResidualDofType;
            };

//------------------------------------------------------------------------------
            /**
             * return IWG residual dof type
             * @param[ out ] mResidualDofType group of dof types for residual
             */
            const moris::Cell< MSI::Dof_Type > & get_residual_dof_type() const
            {
                return mResidualDofType;
            };

//------------------------------------------------------------------------------
            /**
             * set IWG dof type list
             * @param[ in ] aDofTypes list of group of dof types
             * @param[ in ] aIsMaster enum master or slave
             */
            void set_dof_type_list( moris::Cell< moris::Cell< MSI::Dof_Type > > & aDofTypes,
                                    mtk::Master_Slave                             aIsMaster = mtk::Master_Slave::MASTER )
            {
                switch ( aIsMaster )
                {
                    case ( mtk::Master_Slave::MASTER ):
                    {
                        mMasterDofTypes = aDofTypes;
                        break;
                    }
                    case ( mtk::Master_Slave::SLAVE ):
                    {
                        mSlaveDofTypes = aDofTypes;
                        break;
                    }
                    default:
                    {
                        MORIS_ERROR( false, "IWG_User_Defined_Info::set_dof_type_list - can only be MASTER or SLAVE.");
                        break;
                    }
                }
            };

//------------------------------------------------------------------------------
            /**
             * return IWG dof type list
             * @param[ in ]  aIsMaster enum master or slave
             * @param[ out ] mDofTypes list of group of dof types
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
             * set IWG dv type list
             * @param[ in ] aDvTypes list of group of dv types
             * @param[ in ] aIsMaster enum master or slave
             */
            void set_dv_type_list( moris::Cell< moris::Cell< MSI::Dv_Type > > & aDvTypes,
                                    mtk::Master_Slave                            aIsMaster = mtk::Master_Slave::MASTER )
            {
                switch ( aIsMaster )
                {
                    case ( mtk::Master_Slave::MASTER ):
                    {
                        mMasterDvTypes = aDvTypes;
                        break;
                    }
                    case ( mtk::Master_Slave::SLAVE ):
                    {
                        mSlaveDvTypes = aDvTypes;
                        break;
                    }
                    default:
                    {
                        MORIS_ERROR( false, "IWG_User_Defined_Info::set_dv_type_list - can only be MASTER or SLAVE.");
                        break;
                    }
                }
            };

//------------------------------------------------------------------------------
            /**
             * return IWG dv type list
             * @param[ in ]  aIsMaster enum master or slave
             * @param[ out ] mDvTypes list of group of dv types
             */
            const moris::Cell< moris::Cell< MSI::Dv_Type > > & get_dv_type_list( mtk::Master_Slave aIsMaster = mtk::Master_Slave::MASTER ) const
            {
                switch ( aIsMaster )
                {
                    case ( mtk::Master_Slave::MASTER ):
                    {
                        return mMasterDvTypes;
                    }
                    case ( mtk::Master_Slave::SLAVE ):
                    {
                        return mSlaveDvTypes;
                    }
                    default:
                    {
                        MORIS_ERROR( false, "IWG_User_Defined_Info::get_dv_type_list - can only be MASTER or SLAVE.");
                        return mMasterDvTypes;
                    }
                }
            };

//------------------------------------------------------------------------------
            /**
             * set IWG property type list
             * @param[ in ] aPropTypes list of property types
             * @param[ in ] aIsMaster enum master or slave
             */
            void set_property_type_list( moris::Cell< fem::Property_Type > & aPropTypes,
                                         mtk::Master_Slave                   aIsMaster = mtk::Master_Slave::MASTER )
            {
                switch ( aIsMaster )
                {
                    case ( mtk::Master_Slave::MASTER ):
                    {
                        mMasterPropTypes = aPropTypes;
                        break;
                    }
                    case ( mtk::Master_Slave::SLAVE ):
                    {
                        mSlavePropTypes = aPropTypes;
                        break;
                    }
                    default:
                    {
                        MORIS_ERROR( false, "IWG_User_Defined_Info::set_property_type_list - can only be MASTER or SLAVE.");
                        break;
                    }
                }
            };

//------------------------------------------------------------------------------
            /**
             * return IWG property type list
             * @param[ in ]  aIsMaster  enum master or slave
             * @param[ out ] mPropTypes list of property types
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
             * set IWG constitutive type list
             * @param[ in ] aConstitutiveTypes list of constitutive types
             * @param[ in ] aIsMaster          enum master or slave
             */
            void set_constitutive_type_list( moris::Cell< fem::Constitutive_Type > & aConstitutiveTypes,
                                             mtk::Master_Slave                       aIsMaster = mtk::Master_Slave::MASTER )
            {
                switch ( aIsMaster )
                {
                    case ( mtk::Master_Slave::MASTER ):
                    {
                        mMasterConstitutiveTypes = aConstitutiveTypes;
                        break;
                    }
                    case ( mtk::Master_Slave::SLAVE ):
                    {
                        mSlaveConstitutiveTypes = aConstitutiveTypes;
                        break;
                    }
                    default:
                    {
                        MORIS_ERROR( false, "IWG_User_Defined_Info::set_constitutive_type_list - can only be MASTER or SLAVE.");
                        break;
                    }
                }
            };

//------------------------------------------------------------------------------
            /**
             * return IWG constitutive type list
             * @param[ in ]  aIsMaster          enum master or slave
             * @param[ out ] mConstitutiveTypes list of constitutive types
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

//------------------------------------------------------------------------------

    }/* namespace fem */
} /* namespace moris */
#endif /* SRC_FEM_CL_FEM_IWG_USER_DEFINED_INFO_HPP_ */
