/*
 * cl_FEM_IQI_J_Integral.hpp
 *
 *  Created on: Feb 11, 2020
 *      Author: sonne
 */

#ifndef PROJECTS_FEM_INT_SRC_CL_FEM_IQI_J_INTEGRAL_HPP_
#define PROJECTS_FEM_INT_SRC_CL_FEM_IQI_J_INTEGRAL_HPP_

#include "typedefs.hpp"
#include "cl_Cell.hpp"
#include "cl_Matrix.hpp"
#include "linalg_typedefs.hpp"

#include "cl_FEM_IQI.hpp"


namespace moris
{
    namespace fem
    {
//------------------------------------------------------------------------------

    class IQI_J_Integral : public IQI
    {
        enum class IQI_Property_Type
        {
            MAX_ENUM
        };

        // Local string to property enum map
        std::map< std::string, IQI_Property_Type > mPropertyMap;

        enum class IQI_Constitutive_Type
        {
            ELAST_LIN_ISO,
            MAX_ENUM
        };

        // Local string to constitutive enum map
        std::map< std::string, IQI_Constitutive_Type > mConstitutiveMap;

        enum class IQI_Stabilization_Type
        {
            MAX_ENUM
        };

        // Local string to constitutive enum map
        std::map< std::string, IQI_Stabilization_Type > mStabilizationMap;

    public :
//------------------------------------------------------------------------------
        /*
         * constructor
         */
        IQI_J_Integral();

//------------------------------------------------------------------------------
        /*
         * destructor
         */
        ~IQI_J_Integral(){};

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
            // can only be master
            MORIS_ERROR( aIsMaster == mtk::Master_Slave::MASTER,
                         "IQI::set_property - can only be master." );

            // FIXME check that property type makes sense?

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
            // can only be master
            MORIS_ERROR( aIsMaster == mtk::Master_Slave::MASTER,
                         "IQI::set_constitutive model - can only be master." );

            // FIXME check that constitutive string makes sense?

            // set the constitutive model in the constitutive model cell
            this->get_constitutive_models( aIsMaster )( static_cast< uint >( mConstitutiveMap[ aConstitutiveString ] ) ) = aConstitutiveModel;
        }

//------------------------------------------------------------------------------
        /**
         * set stabilization parameter
         * @param[ in ] aStabilizationParameter a stabilization parameter pointer
         * @param[ in ] aStabilizationString    a string defining the stabilization parameter
         */
        void set_stabilization_parameter( std::shared_ptr< Stabilization_Parameter > aStabilizationParameter,
                                          std::string                                aStabilizationString )
        {
            // FIXME check that stabilization string makes sense?

            // set the stabilization parameter in the stabilization parameter cell
            this->get_stabilization_parameters()( static_cast< uint >( mStabilizationMap[ aStabilizationString ] ) ) = aStabilizationParameter;
        }

//------------------------------------------------------------------------------
       /**
        * compute the quantity of interest
        * @param[ in ] aQI quantity of interest matrix to fill
        */
       void compute_QI( Matrix< DDRMat > & aQI );

//------------------------------------------------------------------------------
       /**
        * compute the derivative of the quantity of interest wrt dof types
        * @param[ in ] adQIdDof derivative of quantity of interest matrix to fill
        */
       void compute_dQIdDof( Matrix< DDRMat > & adQIdDof );

//------------------------------------------------------------------------------
    };
//------------------------------------------------------------------------------
    }   // end fem namespace
}       // end moris namespace




#endif /* PROJECTS_FEM_INT_SRC_CL_FEM_IQI_J_INTEGRAL_HPP_ */
