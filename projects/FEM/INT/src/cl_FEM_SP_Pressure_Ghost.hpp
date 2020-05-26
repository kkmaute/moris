/*
 * cl_FEM_SP_Pressure_Ghost.hpp
 *
 *  Created on: Mar 21, 2020
 *  Author: noel
 */

#ifndef SRC_FEM_CL_FEM_SP_PRESSURE_GHOST_HPP_
#define SRC_FEM_CL_FEM_SP_PRESSURE_GHOST_HPP_

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
        /*
         * Stabilization parameter for face-oriented pressure ghost-penalty
         * for different flow regimes,
         * i.e., viscous, convective, and transient regime
         * gamma_p = alpha_p * delta_p * h^(2*i)
         * where i       = interpolation order
         *       delta_p = ( ( viscosity / h )
         *                 + ( density * norm_inf( u ) / 6.0 )
         *                 + ( density * h / ( 12.0 * theta * deltat ) ) )^(-1)
         * from Schott et al. (2015)
         */
        class SP_Pressure_Ghost : public Stabilization_Parameter
        {

                //------------------------------------------------------------------------------
            private:

                // element size
                real mElementSize = 1.0;

                // default dof type
                MSI::Dof_Type mMasterDofVelocity = MSI::Dof_Type::VX;

            public:

                // property type for the SP
                enum class Property_Type
                {
                    VISCOSITY,  // fluid viscosity
                    DENSITY,    // fluid density
                    MAX_ENUM
                };

                // local string to property enum map
                std::map< std::string, Property_Type > mPropertyMap;

                /*
                 * Rem: mParameters( 0 ) - alpha_p
                 *      mParamaters( 1 ) - theta
                 */

                //------------------------------------------------------------------------------
                /*
                 * constructor
                 */
                SP_Pressure_Ghost();

                //------------------------------------------------------------------------------
                /**
                 * trivial destructor
                 */
                ~SP_Pressure_Ghost(){};

                //------------------------------------------------------------------------------
                /**
                 * reset the cluster measures required for this SP
                 */
                void reset_cluster_measures();

                //------------------------------------------------------------------------------
                /**
                 * set dof types
                 * @param[ in ] aDofTypes a cell of cell of dof types
                 * @param[ in ] aDofStrings list of strings describing the dof types
                 * @param[ in ] aIsMaster enum for master or slave
                 */
                void set_dof_type_list(
                        moris::Cell< moris::Cell< MSI::Dof_Type > > & aDofTypes,
                        moris::Cell< std::string >                  & aDofStrings,
                        mtk::Master_Slave                             aIsMaster = mtk::Master_Slave::MASTER );

                //------------------------------------------------------------------------------
                /**
                 * set dv types
                 * @param[ in ] aDvTypes   a cell of group of dv types
                 * @param[ in ] aDvStrings list of strings describing the dv types
                 * @param[ in ] aIsMaster enum for master or slave
                 */
                void set_dv_type_list(
                        moris::Cell< moris::Cell< PDV_Type > > & aDvTypes,
                        moris::Cell< std::string >             & aDvStrings,
                        mtk::Master_Slave                        aIsMaster = mtk::Master_Slave::MASTER )
                {
                    Stabilization_Parameter::set_dv_type_list( aDvTypes, aIsMaster );
                }

                //------------------------------------------------------------------------------
                /**
                 * set property
                 * @param[ in ] aProperty       a property pointer
                 * @param[ in ] aPropertyString a string defining the property
                 * @param[ in ] aIsMaster       an enum for master or slave
                 */
                void set_property(
                        std::shared_ptr< Property > aProperty,
                        std::string                 aPropertyString,
                        mtk::Master_Slave           aIsMaster = mtk::Master_Slave::MASTER );

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
                    MORIS_ERROR( false, "SP_Pressure_Ghost::eval_dSPdMasterDV - not implemented." );
                }

                //------------------------------------------------------------------------------
        };
        //------------------------------------------------------------------------------
    } /* namespace fem */
} /* namespace moris */

#endif /* SRC_FEM_CL_FEM_SP_PRESSURE_GHOST_HPP_ */
