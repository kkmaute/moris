/*
 * cl_FEM_SP_Dirichlet_Nitsche.hpp
 *
 *  Created on: Oct 21, 2019
 *      Author: noel
 */

#ifndef SRC_FEM_CL_FEM_SP_DIRICHLET_NITSCHE_HPP_
#define SRC_FEM_CL_FEM_SP_DIRICHLET_NITSCHE_HPP_

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

        class SP_Dirichlet_Nitsche : public Stabilization_Parameter
        {

//------------------------------------------------------------------------------
        public:
//------------------------------------------------------------------------------
            /*
             * trivial constructor
             */
            SP_Dirichlet_Nitsche()
            {
                // set the penalty type
                mStabilizationType = fem::Penalty_Type::DIRICHLET_NITSCHE;

                // set the list of cluster measures
                mClusterMeasures = { fem::Cluster_Measure::ELEMENT_SIZE };
            };

//------------------------------------------------------------------------------
            /**
             * trivial destructor
             */
            ~SP_Dirichlet_Nitsche(){};

//------------------------------------------------------------------------------
            /**
             * evaluate the penalty parameter value
             */
            void eval_PP();

//------------------------------------------------------------------------------
            /**
             * evaluate the penalty parameter derivative wrt to a master dof type
             * @param[ in ] aDofTypes a dof type wrt which the derivative is evaluated
             * dPPdMasterDOF ( 1 x numDerDof )
             */
            void eval_dPPdMasterDOF( const moris::Cell< MSI::Dof_Type > & aDofTypes );

//------------------------------------------------------------------------------
            /**
             * evaluate the penalty parameter derivative wrt to a master dv type
             * @param[ in ] aDvTypes a dv type wrt which the derivative is evaluated
             * dPPdMasterDV ( 1 x numDerDv )
             */
            void eval_dPPdMasterDV( const moris::Cell< MSI::Dv_Type > & aDvTypes )
            {
                MORIS_ERROR( false, "PP_Dirichlet_Nitsche - eval_dPPdMasterDV: not implemented." );
            }

//------------------------------------------------------------------------------
        };
//------------------------------------------------------------------------------
    } /* namespace fem */
} /* namespace moris */

#endif /* SRC_FEM_CL_FEM_SP_DIRICHLET_NITSCHE_HPP_ */
