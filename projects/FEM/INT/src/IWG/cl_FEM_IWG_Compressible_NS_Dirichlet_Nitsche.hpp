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
#include "cl_FEM_IWG.hpp"

namespace moris
{
    namespace fem
    {
        //------------------------------------------------------------------------------

        class IWG_Compressible_NS_Dirichlet_Nitsche : public IWG
        {

                //------------------------------------------------------------------------------
            public:

                // sign for symmetric/unsymmetric Nitsche
                sint mBeta;

                // default dof types
                MSI::Dof_Type mDofDensity     = MSI::Dof_Type::RHO;
                MSI::Dof_Type mDofPressure    = MSI::Dof_Type::P;
                MSI::Dof_Type mDofTemperature = MSI::Dof_Type::TEMP;
                MSI::Dof_Type mDofVelocity    = MSI::Dof_Type::VX;

                // flag indicating what variable set is being used
                fem::Variable_Set mVariableSet = fem::Variable_Set::UNDEFINED;

                // save number of spatial dimensions and basis functions per field
                bool mSpaceDimEval = true;
                uint mNumSpaceDims;

                // reset flags for storage variables
                bool mSelectMatrixEval = true;
                bool mTestFunctionsEval = true;
                bool mJumpEval = true;
                bool mJumpDofEval = true;
                bool mFluxAMatEval = true;
                bool mFluxAMatDofEval = true;
                bool mFluxKMatEval = true;

                bool mTractionEval = true;
                bool mTractionDofEval = true;
                bool mTestTractionEval = true;
                bool mTestTractionDofEval = true;

                // storage for selection matrix spanning all dof types
                Matrix< DDRMat > mSelectMat;

                // storage for matrix of test functions
                Matrix< DDRMat > mTestFunctions;

                // storage variables for jump term
                Matrix< DDRMat > mJump;
                Matrix< DDRMat > mdJumpdDOF;

                // cells of matrices containing the A flux matrices and their DoF derivatives
                moris::Cell< Matrix< DDRMat > > mA;
                moris::Cell< moris::Cell< Matrix< DDRMat > > > mADOF;

                // cells of matrices containing the K flux matrices and their DoF derivatives
                moris::Cell< moris::Cell< Matrix< DDRMat > > > mK;
                moris::Cell< Matrix< DDRMat > > mKiji;

                // cells of matrices containing the traction terms (K*n)
                Matrix< DDRMat > mTraction;
                Matrix< DDRMat > mTractionDOF;

                Matrix< DDRMat > mTestTraction;
                Matrix< DDRMat > mTestTractionDOF;

                // FIXME: only designed for primitive variables right now
                enum class IWG_Property_Type
                {
                        PRESCRIBED_DOF_1,
                        PRESCRIBED_VELOCITY,
                        SELECT_VELOCITY,
                        PRESCRIBED_DOF_3,
                        UPWIND,
                        DYNAMIC_VISCOSITY,
                        THERMAL_CONDUCTIVITY,
                        MAX_ENUM
                };

                // local material enums
                enum class IWG_Material_Type
                {
                        FLUID_MM,
                        MAX_ENUM
                };

                // local constitutive enums
                enum class IWG_Constitutive_Type
                {
                        FLUID_CM,
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
                void reset_spec_eval_flags();

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
                 * get number of spatial dimensions
                 * @param[ out ] mNumSpaceDims number of spatial dimensions
                 */
                uint num_space_dims();

                //------------------------------------------------------------------------------
                /**
                 * get selection matrix for all dof Types
                 */
                const Matrix< DDRMat > & select_matrix();

                //------------------------------------------------------------------------------
                /**
                 * get matrix of test functions W for all dof types
                 */
                const Matrix< DDRMat > & test_functions();

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
                /**
                 * evaluate the A flux matrices
                 */
                void eval_A_matrices();    

                //------------------------------------------------------------------------------
                /**
                 * get the A flux matrices
                 * @param[ in ]  aK index
                 * @param[ out ] mA A-matrix
                 */
                const Matrix< DDRMat > & A( const uint aK );   

                //------------------------------------------------------------------------------
                /**
                 * get the K flux matrices
                 * @param[ in ]  aI  first index
                 * @param[ in ]  aJ  second index
                 * @param[ out ] Kij K-matrix
                 */
                const Matrix< DDRMat > & K( const uint aI, const uint aJ );     

                //------------------------------------------------------------------------------
                /**
                 * evaluate the DoF derivatives of the A flux matrices for the Jacobians
                 */
                void eval_A_DOF_matrices();

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
                const Matrix< DDRMat > & dTestTractiondDOF( const Matrix< DDRMat > aVR );             

                //------------------------------------------------------------------------------
        };
        //------------------------------------------------------------------------------
    } /* namespace fem */
} /* namespace moris */

#endif /* SRC_FEM_CL_FEM_IWG_COMPRESSIBLE_NS_DIRICHLET_NITSCHE_HPP_ */
