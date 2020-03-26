/*
 * cl_FEM_IWG_Isotropic_Struc_Linear_Bulk.hpp
 *
 *  Created on: Okt 06, 2019
 *      Author: schmidt
 */

#ifndef SRC_FEM_CL_FEM_IWG_ISOTROPIC_STRUC_LINEAR_BULK_HPP_
#define SRC_FEM_CL_FEM_IWG_ISOTROPIC_STRUC_LINEAR_BULK_HPP_

#include <map>

#include "typedefs.hpp"                     //MRS/COR/src
#include "cl_Cell.hpp"                      //MRS/CON/src

#include "cl_Matrix.hpp"                    //LINALG/src
#include "linalg_typedefs.hpp"              //LINALG/src

#include "cl_FEM_IWG.hpp"                   //FEM/INT/src

namespace moris
{
    namespace fem
    {
//------------------------------------------------------------------------------

        class IWG_Isotropic_Struc_Linear_Bulk : public IWG
        {

//------------------------------------------------------------------------------
        public:

            enum class IWG_Property_Type
            {
                LOAD,
                MAX_ENUM
            };

            // Local string to property enum map
            std::map< std::string, IWG_Property_Type > mPropertyMap;

            enum class IWG_Constitutive_Type
            {
                ELAST_LIN_ISO,
                ELAST_LIN_ISO_PRESSURE,
                MAX_ENUM
            };

            // Local string to constitutive enum map
            std::map< std::string, IWG_Constitutive_Type > mConstitutiveMap;

            enum class IWG_Stabilization_Type
            {
                MAX_ENUM
            };

            // Local string to constitutive enum map
            std::map< std::string, IWG_Stabilization_Type > mStabilizationMap;

//------------------------------------------------------------------------------
            /*
             * constructor
             */
            IWG_Isotropic_Struc_Linear_Bulk();

//------------------------------------------------------------------------------
            /**
             * trivial destructor
             */
            ~IWG_Isotropic_Struc_Linear_Bulk(){};

//------------------------------------------------------------------------------
            /**
             * set property
             * @param[ in ] aProperty       a property pointer
             * @param[ in ] aPropertyString a string defining the property
             * @param[ in ] aIsMaster       an enum for master or slave
             */
            void set_property( std::shared_ptr< Property > aProperty,
                               std::string                 aPropertyString,
                               mtk::Master_Slave           aIsMaster = mtk::Master_Slave::MASTER )
            {
                // check that aPropertyString makes sense
                MORIS_ERROR( mPropertyMap.find( aPropertyString ) != mPropertyMap.end(),
                             "IWG_Isotropic_Struc_Linear_Bulk::set_property - Unknown aPropertyString." );

                // check no slave allowed
                MORIS_ERROR( aIsMaster == mtk::Master_Slave::MASTER,
                             "IWG_Isotropic_Struc_Linear_Bulk::set_property - No slave allowed" );

                // set the property in the property cell
                this->get_properties( aIsMaster )( static_cast< uint >( mPropertyMap[ aPropertyString ] ) ) = aProperty;
            }

//------------------------------------------------------------------------------
            /**
             * set constitutive model
             * @param[ in ] aConstitutiveModel  a constitutive model pointer
             * @param[ in ] aConstitutiveString a string defining the constitutive model
             * @param[ in ] aIsMaster           an enum for master or slave
             */
            void set_constitutive_model( std::shared_ptr< Constitutive_Model > aConstitutiveModel,
                                         std::string                           aConstitutiveString,
                                         mtk::Master_Slave                     aIsMaster = mtk::Master_Slave::MASTER )
            {
                // check that aConstitutiveString makes sense
                MORIS_ERROR( mConstitutiveMap.find( aConstitutiveString ) != mConstitutiveMap.end(),
                             "IWG_Isotropic_Struc_Linear_Bulk::set_constitutive_model - Unknown aConstitutiveString." );

                // check no slave allowed
                MORIS_ERROR( aIsMaster == mtk::Master_Slave::MASTER,
                             "IWG_Isotropic_Struc_Linear_Bulk::set_constitutive_model - No slave allowed" );

                // set the constitutive model in the constitutive model cell
                this->get_constitutive_models( aIsMaster )( static_cast< uint >( mConstitutiveMap[ aConstitutiveString ] ) ) = aConstitutiveModel;
            }

//------------------------------------------------------------------------------
            /**
             * compute the residual
             * r =
             * @param[ in ] aResidual residual vector to fill
             */
            void compute_residual( real tWStar );

//------------------------------------------------------------------------------
            /**
             * compute the jacobian
             * j =
             * @param[ in ] aJacobians list of jacobian matrices to fill
             */
            void compute_jacobian( real tWStar );

//------------------------------------------------------------------------------
            /**
             * compute the residual and the jacobian
             * @param[ in ] aJacobians list of jacobian matrices to fill
             * @param[ in ] aResidual  residual vector to fill
             */
            void compute_jacobian_and_residual( real aWStar );

//------------------------------------------------------------------------------
            /**
             * compute the derivative of the residual wrt design variables
             * @param[ in ] aWStar weight associated to the evaluation point
             */
            void compute_dRdp( real aWStar );

//------------------------------------------------------------------------------
        };
//------------------------------------------------------------------------------
    } /* namespace fem */
} /* namespace moris */

#endif /* SRC_FEM_CL_FEM_IWG_ISOTROPIC_STRUC_LINEAR_BULK_HPP_ */
