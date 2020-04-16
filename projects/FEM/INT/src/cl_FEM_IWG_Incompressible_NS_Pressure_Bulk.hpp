/*
 * cl_FEM_IWG_Incompressible_NS_Pressure_Bulk.hpp
 *
 *  Created on: Mar 4, 2020
 *      Author: noel
 */

#ifndef PROJECTS_FEM_INT_SRC_CL_FEM_IWG_INCOMPRESSIBLE_NS_PRESSURE_BULK_HPP_
#define PROJECTS_FEM_INT_SRC_CL_FEM_IWG_INCOMPRESSIBLE_NS_PRESSURE_BULK_HPP_

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

        class IWG_Incompressible_NS_Pressure_Bulk : public IWG
        {

//------------------------------------------------------------------------------
        public:

            // local property enums
            enum class IWG_Property_Type
            {
                DENSITY,
                GRAVITY,
                THERMAL_EXPANSION,
                REF_TEMP,
                MAX_ENUM
            };

            // local string to property enum map
            std::map< std::string, IWG_Property_Type > mPropertyMap;

            // local constitutive enums
            enum class IWG_Constitutive_Type
            {
                INCOMPRESSIBLE_FLUID,
                MAX_ENUM
            };

            // local string to constitutive enum map
            std::map< std::string, IWG_Constitutive_Type > mConstitutiveMap;

            // local stabilization enums
            enum class IWG_Stabilization_Type
            {
                INCOMPRESSIBLE_FLOW,
                MAX_ENUM
            };

            // local string to constitutive enum map
            std::map< std::string, IWG_Stabilization_Type > mStabilizationMap;

//------------------------------------------------------------------------------
            /*
             *  constructor
             */
            IWG_Incompressible_NS_Pressure_Bulk();

//------------------------------------------------------------------------------
            /**
             * trivial destructor
             */
            ~IWG_Incompressible_NS_Pressure_Bulk(){};

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
                             "IWG_Incompressible_NS_Pressure_Bulk::set_property - Unknown aPropertyString." );

                // check no slave allowed
                MORIS_ERROR( aIsMaster == mtk::Master_Slave::MASTER,
                             "IWG_Incompressible_NS_Pressure_Bulk::set_property - No slave allowed." );

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
                             "IWG_Incompressible_NS_Pressure_Bulk::set_constitutive_model - Unknown aConstitutiveString." );

                // check no slave allowed
                MORIS_ERROR( aIsMaster == mtk::Master_Slave::MASTER,
                             "IWG_Incompressible_NS_Pressure_Bulk::set_constitutive_model - No slave allowed." );

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
                // check that aConstitutiveString makes sense
                MORIS_ERROR( mStabilizationMap.find( aStabilizationString ) != mStabilizationMap.end(),
                             "IWG_Incompressible_NS_Pressure_Bulk::set_stabilization_parameter - Unknown aStabilizationString." );

                // set the stabilization parameter in the stabilization parameter cell
                this->get_stabilization_parameters()( static_cast< uint >( mStabilizationMap[ aStabilizationString ] ) ) = aStabilizationParameter;
            }

//------------------------------------------------------------------------------
            /**
             * compute the residual
             * @param[ in ] aWStar weight associated to the evaluation point
             */
            void compute_residual( real aWStar );

//------------------------------------------------------------------------------
            /**
             * compute the jacobian
             * @param[ in ] aWStar weight associated to the evaluation point
             */
            void compute_jacobian( real aWStar );

//------------------------------------------------------------------------------
            /**
             * compute the residual and the jacobian
             * @param[ in ] aWStar weight associated to the evaluation point
             */
            void compute_jacobian_and_residual( real aWStar );

//------------------------------------------------------------------------------
            /**
             * compute the derivative of the residual wrt design variables
             * @param[ in ] aWStar weight associated to the evaluation point
             */
            void compute_dRdp( real aWStar );

        private:
//------------------------------------------------------------------------------
            /**
             * compute the residual strong form
             * @param[ in ] aRM a matrix to fill with RM
             */
            void compute_residual_strong_form( Matrix< DDRMat > & aRM );

//------------------------------------------------------------------------------
            /**
             * compute the residual strong form
             * @param[ in ] aDofTypes a list of dof type wrt which
             *                        the derivative is requested
             * @param[ in ] aJM       a matrix to fill with dRMdDof
             */
            void compute_jacobian_strong_form( moris::Cell< MSI::Dof_Type >   aDofTypes,
                                               Matrix< DDRMat >             & aJM );

//------------------------------------------------------------------------------
            /**
             * compute the term uj vij
             * @param[ in ] aujvij a matrix to fill with uj vij
             */
            void compute_ujvij( Matrix< DDRMat > & aujvij );

//------------------------------------------------------------------------------
            /**
             * compute the term uj vij rm
             * @param[ in ] aujvijrm a matrix to fill with uj vij rm
             * @param[ in ] arm      provided strong form residual
             */
            void compute_ujvijrm( Matrix< DDRMat > & aujvijrm,
                                  Matrix< DDRMat > & arm );

//------------------------------------------------------------------------------
            // FIXME provided directly by the field interpolator?
            /**
             * compute the term dnNdtn
             * @param[ in ] adnNdtn a matrix to fill with dnNdtn
             */
            void compute_dnNdtn( Matrix< DDRMat > & adnNdtn );

//------------------------------------------------------------------------------
        };
//------------------------------------------------------------------------------
    } /* namespace fem */
} /* namespace moris */




#endif /* PROJECTS_FEM_INT_SRC_CL_FEM_IWG_INCOMPRESSIBLE_NS_PRESSURE_BULK_HPP_ */
