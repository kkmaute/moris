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
#include "cl_FEM_Cluster.hpp"     //FEM/INT/src

namespace moris
{
    namespace fem
    {
//------------------------------------------------------------------------------

        class SP_Ghost_Virtual_Work : public Stabilization_Parameter
        {

//------------------------------------------------------------------------------
        private:

            // cluster measures for the SP
            moris::real mElementSize = 1.0;

        public:

//------------------------------------------------------------------------------
            /*
             * constructor
             * Rem: mParameters( 0 ) - gamma penalty parameter
             */
            SP_Ghost_Virtual_Work(){};

//------------------------------------------------------------------------------
            /**
             * trivial destructor
             */
            ~SP_Ghost_Virtual_Work(){};

//------------------------------------------------------------------------------
            /**
             * reset the cluster measures required for this SP
             */
            void reset_cluster_measures();

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
            void eval_dSPdMasterDV( const moris::Cell< GEN_DV > & aDvTypes )
            {
                MORIS_ERROR( false, "SP_Ghost_Virtual_Work::eval_dSPdMasterDV: not implemented." );
            }

//------------------------------------------------------------------------------
        };
//------------------------------------------------------------------------------
    } /* namespace fem */
} /* namespace moris */

#endif /* SRC_FEM_CL_FEM_SP_GHOST_VIRTUAL_WORK_HPP_ */
