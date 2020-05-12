/*
 * cl_FEM_SP_Convective_Ghost.hpp
 *
 *  Created on: Mar 22, 2020
 *  Author: noel
 */

#ifndef SRC_FEM_CL_FEM_SP_CONVECTIVE_GHOST_HPP_
#define SRC_FEM_CL_FEM_SP_CONVECTIVE_GHOST_HPP_

#include <map>
//MRS/CON/src
#include "typedefs.hpp"
#include "cl_Cell.hpp"
//LINALG/src
#include "cl_Matrix.hpp"
#include "linalg_typedefs.hpp"
//FEM/INT/src
#include "cl_FEM_Field_Interpolator.hpp"
#include "cl_FEM_Constitutive_Model.hpp"
#include "cl_FEM_Stabilization_Parameter.hpp"
#include "cl_FEM_Cluster.hpp"

namespace moris
{
    namespace fem
    {
//------------------------------------------------------------------------------

        class SP_Convective_Ghost : public Stabilization_Parameter
        {

//------------------------------------------------------------------------------
        private:

            // element size
            real mElementSize = 1.0;

        public:

            // property type for the SP
            enum class Property_Type
            {
                DENSITY,  // fluid density
                MAX_ENUM
            };

            // local string to property enum map
            std::map< std::string, Property_Type > mPropertyMap;

            /*
             * Rem: mParameters( 0 ) - alpha_velocity
             */

//------------------------------------------------------------------------------
            /*
             * constructor
             */
            SP_Convective_Ghost();

//------------------------------------------------------------------------------
            /**
             * trivial destructor
             */
            ~SP_Convective_Ghost(){};

//------------------------------------------------------------------------------
            /**
             * reset the cluster measures required for this SP
             */
            void reset_cluster_measures();

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
                             "SP_Viscous_Ghost::set_property - Unknown aPropertyString." );

                // set the property in the property cell
                this->get_properties( aIsMaster )( static_cast< uint >( mPropertyMap[ aPropertyString ] ) ) = aProperty;
            }

//------------------------------------------------------------------------------
            /**
             * evaluate the stabilization parameter value
             */
            void eval_SP();

//------------------------------------------------------------------------------
            /**
             * evaluate the stabilization parameter derivative wrt to a master dof type
             * @param[ in ] aDofTypes a dof type wrt which the derivative is evaluated
             */
            void eval_dSPdMasterDOF( const moris::Cell< MSI::Dof_Type > & aDofTypes );

//------------------------------------------------------------------------------
            /**
             * evaluate the penalty parameter derivative wrt to a master dv type
             * @param[ in ] aDvTypes a dv type wrt which the derivative is evaluated
             */
            void eval_dSPdMasterDV( const moris::Cell< PDV_Type > & aDvTypes )
            {
                MORIS_ERROR( false, "SP_Convective_Ghost::eval_dSPdMasterDV - not implemented." );
            }

//------------------------------------------------------------------------------
        };
//------------------------------------------------------------------------------
    } /* namespace fem */
} /* namespace moris */

#endif /* SRC_FEM_CL_FEM_SP_CONVECTIVE_GHOST_HPP_ */
