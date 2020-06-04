/*
 * cl_FEM_SP_Turbulence_Dirichlet_Nitsche.hpp
 *
 *  Created on: Oct 21, 2019
 *  Author: noel
 */

#ifndef SRC_FEM_CL_FEM_SP_TURBULENCE_DIRICHLET_NITSCHE_HPP_
#define SRC_FEM_CL_FEM_SP_TURBULENCE_DIRICHLET_NITSCHE_HPP_

#include <map>

#include "typedefs.hpp"                     //MRS/COR/src
#include "cl_Cell.hpp"                      //MRS/CON/src

#include "cl_Matrix.hpp"                    //LINALG/src
#include "linalg_typedefs.hpp"              //LINALG/src

#include "cl_FEM_Field_Interpolator.hpp"    //FEM/INT/src
#include "cl_FEM_Constitutive_Model.hpp"    //FEM/INT/src
#include "cl_FEM_Stabilization_Parameter.hpp"     //FEM/INT/src
#include "cl_FEM_Cluster.hpp"
namespace moris
{
    namespace fem
    {
        //------------------------------------------------------------------------------

        class SP_Turbulence_Dirichlet_Nitsche : public Stabilization_Parameter
        {

                //------------------------------------------------------------------------------
            private:

                // cluster measures for the SP
                real mElementSize = 1.0;

                // populate the dof map (default)
                MSI::Dof_Type mMasterDofViscosity = MSI::Dof_Type::VISCOSITY;

                // FIXME temp all the constants
                real mSigma = 2.0/3.0;

            public:

                // Property type for the SP
                enum class SP_Property_Type
                {
                    VISCOSITY,
                    MAX_ENUM
                };

                // Local string to property enum map
                std::map< std::string, SP_Property_Type > mPropertyMap;

                //------------------------------------------------------------------------------
                /*
                 * constructor
                 * Rem: mParameters( 0 ) - gamma penalty parameter
                 */
                SP_Turbulence_Dirichlet_Nitsche();

                //------------------------------------------------------------------------------
                /**
                 * trivial destructor
                 */
                ~SP_Turbulence_Dirichlet_Nitsche(){};

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
                void eval_dSPdMasterDV( const moris::Cell< PDV_Type > & aDvTypes )
                {
                    MORIS_ERROR( false, "SP_Turbulence_Dirichlet_Nitsche - eval_dSPdMasterDV: not implemented." );
                }

                //------------------------------------------------------------------------------
        };
        //------------------------------------------------------------------------------
    } /* namespace fem */
} /* namespace moris */

#endif /* SRC_FEM_CL_FEM_SP_TURBULENCE_DIRICHLET_NITSCHE_HPP_ */
