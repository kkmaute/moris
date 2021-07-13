/**
 * cl_FEM_IWG_Compressible_NS_Velocity_Dirichlet_Nitsche.hpp
 *
 *  Created on: Mar 17, 2021
 *      Author: wunsch
 */

#ifndef SRC_FEM_CL_FEM_IWG_COMPRESSIBLE_NS_DIRICHLET_NITSCHE_HPP_
#define SRC_FEM_CL_FEM_IWG_COMPRESSIBLE_NS_DIRICHLET_NITSCHE_HPP_

#include <map>
//MRS/COR/src
#include "typedefs.hpp"
#include "cl_Cell.hpp"
//LINALG/src
#include "cl_Matrix.hpp"
#include "linalg_typedefs.hpp"
//FEM/INT/src
#include "cl_FEM_IWG_Compressible_NS_Base.hpp"

namespace moris
{
    namespace fem
    {
        //------------------------------------------------------------------------------

        class IWG_Compressible_NS_Dirichlet_Nitsche : public IWG_Compressible_NS_Base
        {

                //------------------------------------------------------------------------------
            public:

                // sign for symmetric/unsymmetric Nitsche
                sint mBeta;

                // reset flags for storage variables
                bool mSelectMatrixEval = true;
                bool mJumpEval = true;
                bool mJumpDofEval = true;

                bool mTractionEval = true;
                bool mTractionDofEval = true;
                bool mTestTractionEval = true;
                bool mTestTractionDofEval = true;

                bool mUpwindOperatorEval = true;

                // storage for selection matrix spanning all dof types
                Matrix< DDRMat > mSelectMat;

                // storage variables for jump term
                Matrix< DDRMat > mJump;
                Matrix< DDRMat > mdJumpdDOF;

                // cells of matrices containing the traction terms (K*n)
                Matrix< DDRMat > mTraction;
                Matrix< DDRMat > mTractionDOF;

                Matrix< DDRMat > mTestTraction;
                Matrix< DDRMat > mTestTractionDOF;

                // matrices containing the pressure upwind operator and its state variable derivative
                Matrix< DDRMat > mUpwindOperator;
                Matrix< DDRMat > mdUpwindOperatordY;

                // FIXME: only designed for primitive variables right now
                enum class IWG_Property_Type
                {
                    DYNAMIC_VISCOSITY,
                    THERMAL_CONDUCTIVITY,
                    BODY_FORCE,
                    BODY_HEAT_LOAD,
                    PRESCRIBED_DOF_1,
                    PRESCRIBED_VELOCITY,
                    SELECT_VELOCITY,
                    PRESCRIBED_DOF_3,
                    PRESSUREUPWIND,
                    MAX_ENUM
                };

                // local stabilization enums
                enum class IWG_Stabilization_Type
                {
                    NITSCHE_PENALTY_PARAMETER,
                    MAX_ENUM
                };

                //------------------------------------------------------------------------------
                /*
                 *  constructor
                 */
                IWG_Compressible_NS_Dirichlet_Nitsche( sint aBeta );

                //------------------------------------------------------------------------------
                /**
                 * trivial destructor
                 */
                ~IWG_Compressible_NS_Dirichlet_Nitsche(){};

                //------------------------------------------------------------------------------
                /**
                 * reset eval flags specific to this IWG
                 */
                void reset_child_eval_flags();

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

                //------------------------------------------------------------------------------
            
            private:

                //------------------------------------------------------------------------------
                /**
                 * get selection matrix for all dof Types
                 */
                const Matrix< DDRMat > & select_matrix();

                //------------------------------------------------------------------------------
                //------------------------------------------------------------------------------
                /**
                 * assemble the difference between state variables and the prescribed values
                 * into a vector
                 */
                const Matrix< DDRMat > & jump();

                //------------------------------------------------------------------------------
                /**
                 * assemble the DoF derivative of the difference between state variables and 
                 * the prescribed values into a matrix
                 */
                const Matrix< DDRMat > & dJumpdDOF();     

                //------------------------------------------------------------------------------
                //------------------------------------------------------------------------------
                /**
                 * evaluate and get the Traction term ( K_ij * Y_,j * n_i )
                 */
                const Matrix< DDRMat > & Traction();     

                //------------------------------------------------------------------------------
                /**
                 * evaluate and get the Dof-derivative of the Traction term d( K_ij * Y_,j * n_i ) / dDof
                 */
                const Matrix< DDRMat > & dTractiondDOF();    

                //------------------------------------------------------------------------------
                /**
                 * evaluate and get the Test Traction term d( K_ij * Y_,j * n_i )^T / dDof 
                 */
                const Matrix< DDRMat > & TestTraction(); 

                //------------------------------------------------------------------------------
                /**
                 * evaluate and get the Dof-derivative of the Test Traction term d^2( K_ij * Y_,j * n_i )^T / dDof^2 * VR
                 */
                const Matrix< DDRMat > & dTestTractiondDOF( const Matrix< DDRMat > aVL );   

                //------------------------------------------------------------------------------
                /**
                 * evaluate and get the upwind operator
                 */
                const Matrix< DDRMat > & UpwindOperator(); 

                //------------------------------------------------------------------------------
                /**
                 * evaluate and get the state-variable derivative of the upwind operator,
                 * pre-multiplied from the right with a constant vector
                 */
                const Matrix< DDRMat > & dUpwindOperatordY( const Matrix< DDRMat > aVR );

                //------------------------------------------------------------------------------
        };
        //------------------------------------------------------------------------------
    } /* namespace fem */
} /* namespace moris */

#endif /* SRC_FEM_CL_FEM_IWG_COMPRESSIBLE_NS_DIRICHLET_NITSCHE_HPP_ */
