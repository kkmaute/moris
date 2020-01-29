/*
 * cl_FEM_SP_Ghost_Virtual_Work.hpp
 *
 *  Created on: Jan 23, 2019
 *  Author: noel
 */

#ifndef SRC_FEM_CL_FEM_SP_GHOST_VIRTUAL_WORK_HPP_
#define SRC_FEM_CL_FEM_SP_GHOST_VIRTUAL_WORK_HPP_

#include <map>

#include "typedefs.hpp"                     //MRS/COR/src
#include "cl_Cell.hpp"                      //MRS/CON/src

#include "cl_Matrix.hpp"                    //LINALG/src
#include "linalg_typedefs.hpp"              //LINALG/src

#include "cl_FEM_Field_Interpolator.hpp"    //FEM/INT/src
#include "cl_FEM_Constitutive_Model.hpp"    //FEM/INT/src
#include "cl_FEM_Stabilization_Parameter.hpp"     //FEM/INT/src

namespace moris
{
    namespace fem
    {
//------------------------------------------------------------------------------

        class SP_Ghost_Virtual_Work : public Stabilization_Parameter
        {

//------------------------------------------------------------------------------
        public:

            moris::real mElementSize;

            enum class SP_Property_Type
            {
                MAX_ENUM
            };

            // Local string to property enum map
            std::map< std::string, SP_Property_Type > mPropertyMap;

            enum class SP_Constitutive_Type
            {
                MAX_ENUM
            };

            // Local string to constitutive enum map
            std::map< std::string, SP_Constitutive_Type > mConstitutiveMap;

//------------------------------------------------------------------------------
            /*
             * constructor
             * Rem: mParameters( 0 ) - gamma penalty parameter
             *      mParameters( 1 ) - interpolation order
             */
            SP_Ghost_Virtual_Work();

//------------------------------------------------------------------------------
            /**
             * trivial destructor
             */
            ~SP_Ghost_Virtual_Work(){};

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
                // FIXME check that constitutive string makes sense?

                // set the constitutive model in the constitutive model cell
                this->get_constitutive_models( aIsMaster )( static_cast< uint >( mConstitutiveMap[ aConstitutiveString ] ) ) = aConstitutiveModel;
            }

//------------------------------------------------------------------------------
            /**
             * evaluate the penalty parameter value
             */
            void eval_SP();

//------------------------------------------------------------------------------
            /**
             * evaluate the penalty parameter derivative wrt to a master dof type
             * @param[ in ] aDofTypes a dof type wrt which the derivative is evaluated
             * dPPdMasterDOF ( 1 x numDerDof )
             */
            void eval_dSPdMasterDOF( const moris::Cell< MSI::Dof_Type > & aDofTypes );

//------------------------------------------------------------------------------
            /**
             * evaluate the penalty parameter derivative wrt to a master dv type
             * @param[ in ] aDvTypes a dv type wrt which the derivative is evaluated
             * dPPdMasterDV ( 1 x numDerDv )
             */
            void eval_dSPdMasterDV( const moris::Cell< MSI::Dv_Type > & aDvTypes )
            {
                MORIS_ERROR( false, "SP_Ghost_Virtual_Work::eval_dSPdMasterDV: not implemented." );
            }

//------------------------------------------------------------------------------
        };
//------------------------------------------------------------------------------
    } /* namespace fem */
} /* namespace moris */

#endif /* SRC_FEM_CL_FEM_SP_GHOST_VIRTUAL_WORK_HPP_ */
