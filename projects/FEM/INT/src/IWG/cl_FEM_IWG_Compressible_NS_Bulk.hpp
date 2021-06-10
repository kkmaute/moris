/*
 * cl_FEM_IWG_Compressible_NS_Bulk.hpp
 *
 *  Created on: Feb 10, 2021
 *      Author: wunsch
 *
 * IWG for unified NS equations using flux matrices
 * operating on all residual DoF types simultaneously
 */

#ifndef SRC_FEM_CL_FEM_IWG_COMPRESSIBLE_NS_BULK_HPP_
#define SRC_FEM_CL_FEM_IWG_COMPRESSIBLE_NS_BULK_HPP_

#include <map>
#include "typedefs.hpp"                         //MRS/COR/src
#include "cl_Cell.hpp"                          //MRS/CON/src

#include "cl_Matrix.hpp"                        //LINALG/src
#include "linalg_typedefs.hpp"                  //LINALG/src

#include "cl_FEM_IWG_Compressible_NS_Base.hpp"  //FEM/INT/src

namespace moris
{
    namespace fem
    {
        //------------------------------------------------------------------------------

        class IWG_Compressible_NS_Bulk : public IWG_Compressible_NS_Base
        {
                //------------------------------------------------------------------------------

            private:

                // evaluation flags for L-operators
                bool mLYEval = true;
                bool mLWEval = true;
                bool mLDofYEval = true;
                bool mGLSTestFuncEval = true;

                // evaluation flag for the variable mapping operator
                bool mA0invEval = true;
                bool mdA0invdYEval = true;
                
                // evaluation flag for G-operator
                bool mGEval = true;

                // evaluation flags for M-operator and its inverse
                bool mMEval = true;
                bool mMinvEval = true;
                bool mSqrtMinvEval = true;
                bool mdMdYEval = true;

                // evaluation flag for GLS stabilization operator
                bool mTauEval = true; 

                // storage for L-operators
                Matrix< DDRMat > mLY;
                Matrix< DDRMat > mdLdDofY;
                Matrix< DDRMat > mLW;
                Matrix< DDRMat > mdLdDofW;
                Matrix< DDRMat > mGLSTestFunc;
                Matrix< DDRMat > mdGLSTestFuncdDof;

                // storage for the variable mapping operator
                Matrix< DDRMat > mA0inv;
                moris::Cell< Matrix< DDRMat > > mdA0invdY;

                // storage for G-operator
                Matrix< DDRMat > mG;

                // storage for M-operator and its inverse
                Matrix< DDRMat > mM;
                Matrix< DDRMat > mMinv;
                Matrix< DDRMat > mSqrtMinv;
                moris::Cell< Matrix< DDRMat > > mdMdY;

                // storage for the GLS stabiliation operator
                Matrix< DDRMat > mTau;
                Matrix< DDRMat > mdTaudY; 

                //------------------------------------------------------------------------------

            public:

                // local SP enums
                // FIXME: generic SP only for testing purposes for strong form
                enum class IWG_Stabilization_Type
                {
                    GLS,
                    MAX_ENUM
                };

                //------------------------------------------------------------------------------
                /*
                 *  constructor
                 */
                IWG_Compressible_NS_Bulk();

                //------------------------------------------------------------------------------
                /**
                 * trivial destructor
                 */
                ~IWG_Compressible_NS_Bulk(){};

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

            private:

                //------------------------------------------------------------------------------
                /**
                 * compute the L-operator (Navier-Stokes-Operator) 
                 * applied to the state variable vector Y
                 * @param[ out ] mLY  L-operator applied to Y-vector
                 */
                const Matrix< DDRMat > & LY();

                //------------------------------------------------------------------------------
                /**
                 * compute the dof derivative of the L-operator (Navier-Stokes-Operator) 
                 * applied to the state variable vector Y
                 * @param[ out ] mdLdDofY  dof deriv of L-operator applied to Y-vector
                 */
                const Matrix< DDRMat > & dLdDofY();

                //------------------------------------------------------------------------------
                /**
                 * compute the L-operator (Navier-Stokes-Operator) 
                 * applied to the test functions (matrix W)
                 * @param[ out ] mLW  L-operator applied to W-matrix
                 */
                const Matrix< DDRMat > & LW();

                //------------------------------------------------------------------------------
                /**
                 * compute the dof derivative of the L-operator (Navier-Stokes-Operator) 
                 * applied to the test functions (matrix W)
                 * @param[ in ]  aVL       pre-multiplication vector applied from the left
                 * @param[ out ] mdLdDofW  dof deriv of L-operator applied to W-matrix
                 */
                const Matrix< DDRMat > & dLdDofW( const Matrix< DDRMat > & aVL );

                //------------------------------------------------------------------------------
                /**
                 * compute the transpose of the L-operator (Navier-Stokes-Operator) 
                 * applied to the test functions (matrix W)
                 * @param[ out ] mGLSTestFunc  transpose of ( the transpose of the L-operator applied to the W-matrix )
                 */
                const Matrix< DDRMat > & GLSTestFunc();

                //------------------------------------------------------------------------------
                /**
                 * compute the dof derivative of the L-operator (Navier-Stokes-Operator) 
                 * applied to the test functions (matrix W)
                 * @param[ in ]  aVR                pre-multiplication vector applied from the right
                 * @param[ out ] mdGLSTestFuncdDof  dof deriv of transpose of ( the transpose of the L-operator applied to the W-matrix )
                 */
                const Matrix< DDRMat > & dGLSTestFuncdDof( const Matrix< DDRMat > & aVR );

                //------------------------------------------------------------------------------
                /**
                 * get the operator mapping from the state variables to conservative variables Y/U
                 * @param[ out ] mA0inv variale mapping operator mapping
                 */
                const Matrix< DDRMat > & A0inv();

                //------------------------------------------------------------------------------
                /**
                 * get the G operator G_ij = sum_d( dxi_d / dx_i * dxi_d / dx_j )
                 * @param[ out ] mG G operator
                 */
                const Matrix< DDRMat > & G();

                //------------------------------------------------------------------------------
                /**
                 * get the M term for the stabilization operator
                 * @param[ out ] mM M operator
                 */
                const Matrix< DDRMat > & M();

                //------------------------------------------------------------------------------
                /**
                 * get the inverse of the M term for the stabiliztion operator
                 * @param[ out ] mMinv  inverse of the M operator
                 */
                const Matrix< DDRMat > & Minv();

                //------------------------------------------------------------------------------
                /**
                 * get the square root of the inverse of the M term for the stabiliztion operator
                 * @param[ out ] mSqrtMinv  square root of the inverse of the M operator
                 */
                const Matrix< DDRMat > & SqrtMinv();
                
                //------------------------------------------------------------------------------
                /**
                 * get the stabilization operator
                 * @param[ out ] mTau stabilization operator
                 */
                const Matrix< DDRMat > & Tau();

                //------------------------------------------------------------------------------
                /**
                 * get the state variable derivatives of M term
                 * @param[ out ] mdMdY state variable derivatives of M term
                 */
                const Matrix< DDRMat > & dMdY( const uint aYind );

                //------------------------------------------------------------------------------
                /**
                 * get the state variable derivatives of the inverse of the A0 matrix
                 * @param[ out ] mdA0invdY state var derivatives of the inverse of the A0 matrix
                 */
                const Matrix< DDRMat > & dA0invdY( const uint aYind );

                //------------------------------------------------------------------------------
                /**
                 * get the deriv of the stabilization operator wrt to the state variables
                 * @param[ in ]  aVR     a constant pre-multiplication vector acting on the right index of Tau
                 * @param[ out ] dTaudY  deriv of the stabilization operator wrt.to the state variables
                 */
                const Matrix< DDRMat > & dTaudY( const Matrix< DDRMat > aVR );

                //------------------------------------------------------------------------------
                /**
                 * get the deriv of sqrt( inv ( M ) ) operator wrt to the DoFs for debugging purposes
                 * @param[ in ]  aColInd  index of column of sqrt( inv ( M ) ) to be finite difference checked
                 * @param[ out ] aPerturbation  absolute perturbation size of DoFs for FD
                 */
                Matrix< DDRMat > dSqrtMinvdu_FD( const uint aColInd, const real aPerturbation );

                //------------------------------------------------------------------------------

        };
        //------------------------------------------------------------------------------
    } /* namespace fem */
} /* namespace moris */

#endif /* SRC_FEM_CL_FEM_IWG_COMPRESSIBLE_NS_VELOCITY_BULK_HPP_ */
