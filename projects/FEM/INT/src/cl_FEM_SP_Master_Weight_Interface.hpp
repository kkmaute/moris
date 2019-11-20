/*
 * cl_FEM_SP_Master_Weight_Interface.hpp
 *
 *  Created on: Nov 18, 2019
 *      Author: noel
 */

#ifndef PROJECTS_FEM_INT_SRC_CL_FEM_SP_MASTER_WEIGHT_INTERFACE_HPP_
#define PROJECTS_FEM_INT_SRC_CL_FEM_SP_MASTER_WEIGHT_INTERFACE_HPP_

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

        class SP_Master_Weight_Interface : public Stabilization_Parameter
        {

//------------------------------------------------------------------------------
        public:
//------------------------------------------------------------------------------
            /*
             * constructor
             */
            SP_Master_Weight_Interface()
            {
                // set the penalty type
                mStabilizationType = fem::Stabilization_Type::MASTER_WEIGHT_INTERFACE;

                // set the list of cluster measures
                mClusterMeasures = { fem::Cluster_Measure::MASTER_VOLUME,
                                     fem::Cluster_Measure::SLAVE_VOLUME };
            };

//------------------------------------------------------------------------------
            /**
             * trivial destructor
             */
            ~SP_Master_Weight_Interface(){};

//------------------------------------------------------------------------------
            /**
             * evaluate the stabilization parameter value
             */
            void eval_SP();

//------------------------------------------------------------------------------
            /**
             * evaluate the stabilization parameter derivative wrt to a master dof type
             * @param[ in ] aDofTypes a dof type wrt which the derivative is evaluated
             * dSPdMasterDOF ( 1 x numDerDof )
             */
            void eval_dSPdMasterDOF( const moris::Cell< MSI::Dof_Type > & aDofTypes );

//------------------------------------------------------------------------------
            /**
             * evaluate the stabilization parameter derivative wrt to a slave dof type
             * @param[ in ] aDofTypes a dof type wrt which the derivative is evaluated
             * dSPdSlaveDOF ( 1 x numDerDof )
             */
             void eval_dSPdSlaveDOF( const moris::Cell< MSI::Dof_Type > & aDofTypes );

//------------------------------------------------------------------------------
            /**
             * evaluate the penalty parameter derivative wrt to a master dv type
             * @param[ in ] aDvTypes a dv type wrt which the derivative is evaluated
             * dPPdMasterDV ( 1 x numDerDv )
             */
            void eval_dSPdMasterDV( const moris::Cell< MSI::Dv_Type > & aDvTypes )
            {
                MORIS_ERROR( false, "SP_Master_Weight_Interface::eval_dSPdMasterDV: not implemented." );
            }

//------------------------------------------------------------------------------
            /**
             * evaluate the penalty parameter derivative wrt to a slave dv type
             * @param[ in ] aDvTypes a dv type wrt which the derivative is evaluated
             * dSPdSlaveDV ( 1 x numDerDv )
             */
             void eval_dSPdSlaveDV( const moris::Cell< MSI::Dv_Type > & aDvTypes )
             {
                 MORIS_ERROR( false, "SP_Master_Weight_Interface::eval_dSPdSlaveDV: not implemented." );
             }
//------------------------------------------------------------------------------
        };
//------------------------------------------------------------------------------
    } /* namespace fem */
} /* namespace moris */

#endif /* PROJECTS_FEM_INT_SRC_CL_FEM_SP_MASTER_WEIGHT_INTERFACE_HPP_ */
