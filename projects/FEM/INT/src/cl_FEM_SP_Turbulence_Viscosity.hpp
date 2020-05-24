/*
 * cl_FEM_SP_Turbulence_Viscosity.hpp
 *
 *  Created on: May 19, 2020
 *  Author: noel
 */

#ifndef SRC_FEM_CL_FEM_SP_TURBULENCE_VISCOSITY_HPP_
#define SRC_FEM_CL_FEM_SP_TURBULENCE_VISCOSITY_HPP_

#include <map>
//MRS/COR/src
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

        class SP_Turbulence_Viscosity : public Stabilization_Parameter
        {

                //------------------------------------------------------------------------------
            private :

                // FIXME temp all the constants
                real mCv1 = 7.1;

            public :

                // property type for the SP
                enum class Property_Type
                {
                    VISCOSITY, // fluid viscosity
                    MAX_ENUM
                };

                // local string to property enum map
                std::map< std::string, Property_Type > mPropertyMap;

                //------------------------------------------------------------------------------
                /*
                 * constructor
                 */
                SP_Turbulence_Viscosity();

                //------------------------------------------------------------------------------
                /**
                 * trivial destructor
                 */
                ~SP_Turbulence_Viscosity(){};

                //------------------------------------------------------------------------------
                /**
                 * reset the cluster measures required for this SP
                 */
                void reset_cluster_measures();

                //------------------------------------------------------------------------------
                /**
                 * set function pointers for evaluation
                 */
                void set_function_pointers();

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
                 * evaluate the penalty parameter value
                 */
                void eval_SP();

                //------------------------------------------------------------------------------
                /**
                 * evaluate the penalty parameter derivative wrt to a master dof type
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
                    MORIS_ERROR( false, "SP_Turbulence_Viscosity::eval_dSPdMasterDV - not implemented." );
                }

                //------------------------------------------------------------------------------
            private:

                //------------------------------------------------------------------------------
                /**
                 * compute chi = viscosityDof / viscosityPtop
                 * @param[ out ] chi
                 */
                real compute_chi();

                //------------------------------------------------------------------------------
                /**
                 * compute the derivative of chi wrt to a dof type
                 * @param[ in ] aDofTypes  a list of dof type wrt which
                 *                         the derivative is requested
                 * @param[ in ] adchidu    a matrix to fill with dchidu
                 */
                void compute_dchidu(
                        moris::Cell< MSI::Dof_Type >   aDofTypes,
                        Matrix< DDRMat >             & adchidu );

                //------------------------------------------------------------------------------
                /**
                 * compute fv1 = chi³ / ( chi³ + cv1³)
                 * @param[ out ] fv1
                 */
                real compute_fv1();

                //------------------------------------------------------------------------------
                /**
                 * compute the derivative of fv1 wrt to a dof type
                 * @param[ in ] aDofTypes  a list of dof type wrt which
                 *                         the derivative is requested
                 * @param[ in ] adfv1du    a matrix to fill with dfv1du
                 */
                void compute_dfv1du(
                        moris::Cell< MSI::Dof_Type >   aDofTypes,
                        Matrix< DDRMat >             & adfv1du );

                //------------------------------------------------------------------------------
        };
        //------------------------------------------------------------------------------
    } /* namespace fem */
} /* namespace moris */

#endif /* SRC_FEM_CL_FEM_SP_TURBULENCE_VISCOSITY_HPP_ */
