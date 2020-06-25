/*
 * cl_FEM_IWG_Isotropic_Struc_Linear_Contact_Nitsche.hpp
 *
 *  Created on: Feb 18, 2020
 *      Author: ritzert
 */

#ifndef PROJECTS_FEM_INT_SRC_CL_FEM_IWG_ISOTROPIC_STRUC_LINEAR_CONTACT_NITSCHE_HPP_
#define PROJECTS_FEM_INT_SRC_CL_FEM_IWG_ISOTROPIC_STRUC_LINEAR_CONTACT_NITSCHE_HPP_

#include <map>

#include "typedefs.hpp"                     //MRS/COR/src
#include "cl_Cell.hpp"                      //MRS/CON/src

#include "cl_Matrix.hpp"                    //LINALG/src
#include "linalg_typedefs.hpp"              //LINALG/src

#include "cl_FEM_Field_Interpolator.hpp"    //FEM/INT/src
#include "cl_FEM_IWG.hpp"                   //FEM/INT/src

namespace moris
{
    namespace fem
    {
//------------------------------------------------------------------------------

        class IWG_Isotropic_Struc_Linear_Contact_Nitsche : public IWG
        {

        public:
            
            // sign for symmetric/unsymmetric Nitsche
            sint mBeta = 1;

            enum class IWG_Property_Type
            {
                MATERIAL,
                MAX_ENUM
            };
            // Local string to property enum map
            std::map< std::string, IWG_Property_Type > mPropertyMap;

            enum class IWG_Constitutive_Type
            {
                ELAST_LIN_ISO,
                MAX_ENUM
            };

            // Local string to constitutive enum map
            std::map< std::string, IWG_Constitutive_Type > mConstitutiveMap;

            enum class IWG_Stabilization_Type
                        {
                            PENALTY_CONTACT,
                            MASTER_WEIGHT_INTERFACE,
                            SLAVE_WEIGHT_INTERFACE,
                            MAX_ENUM
                        };

            // Local string to constitutive enum map
            std::map< std::string, IWG_Stabilization_Type > mStabilizationMap;

            //real mStabParam = 1500;

//------------------------------------------------------------------------------
            /*
             * constructor
             */
            IWG_Isotropic_Struc_Linear_Contact_Nitsche(sint aBeta);

//------------------------------------------------------------------------------
            /**
             * trivial destructor
             */
            ~IWG_Isotropic_Struc_Linear_Contact_Nitsche(){};

//------------------------------------------------------------------------------
            /**
             * set property
             * @param[ in ] aProperty       a property pointer
             * @param[ in ] aPropertyString a string defining the property
             * @param[ in ] aIsMaster       an enum for master or slave
             */
            void set_property( std::shared_ptr< Property > aProperty,
                               std::string                 aPropertyString,
                               mtk::Master_Slave           aIsMaster = mtk::Master_Slave::MASTER );

//------------------------------------------------------------------------------
            /**
             * set constitutive model
             * @param[ in ] aConstitutiveModel  a constitutive model pointer
             * @param[ in ] aConstitutiveString a string defining the constitutive model
             * @param[ in ] aIsMaster           an enum for master or slave
             */
            void set_constitutive_model( std::shared_ptr< Constitutive_Model > aConstitutiveModel,
                                         std::string                           aConstitutiveString,
                                         mtk::Master_Slave                     aIsMaster = mtk::Master_Slave::MASTER );

//------------------------------------------------------------------------------
        /**
        * set stabilization parameter
            * @param[ in ] aStabilizationParameter a stabilization parameter pointer
            * @param[ in ] aStabilizationString    a string defining the stabilization parameter
            */
        void set_stabilization_parameter(
                std::shared_ptr< Stabilization_Parameter > aStabilizationParameter,
                std::string                                aStabilizationString );



//------------------------------------------------------------------------------
            /**
             * compute the residual
             * @param[ in ] aResidual cell of residual vectors to fill
             */
            void compute_residual( real aWStar );

//------------------------------------------------------------------------------
            /**
             * compute the jacobian
             * @param[ in ] aJacobians cell of cell of jacobian matrices to fill
             */
            void compute_jacobian( real aWStar );

//------------------------------------------------------------------------------
            /**
             * compute the residual and the jacobian
             * @param[ in ] aJacobians cell of cell of jacobian matrices to fill
             * @param[ in ] aResidual  cell of residual vectors to fill
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





#endif /* PROJECTS_FEM_INT_SRC_CL_FEM_IWG_ISOTROPIC_STRUC_LINEAR_CONTACT_NITSCHE_HPP_ */
