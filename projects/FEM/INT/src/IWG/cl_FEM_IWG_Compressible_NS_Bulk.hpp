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
#include "typedefs.hpp"                     //MRS/COR/src
#include "cl_Cell.hpp"                      //MRS/CON/src

#include "cl_Matrix.hpp"                    //LINALG/src
#include "linalg_typedefs.hpp"              //LINALG/src

#include "cl_FEM_IWG.hpp"                   //FEM/INT/src

namespace moris
{
    namespace fem
    {
        //------------------------------------------------------------------------------

        class IWG_Compressible_NS_Bulk : public IWG
        {
                //------------------------------------------------------------------------------

            private:

                // default dof types (FIXME: only primitive vars for now)
                MSI::Dof_Type mDofDensity     = MSI::Dof_Type::RHO;
                MSI::Dof_Type mDofPressure    = MSI::Dof_Type::P;
                MSI::Dof_Type mDofVelocity    = MSI::Dof_Type::VX;
                MSI::Dof_Type mDofTemperature = MSI::Dof_Type::TEMP;

                // save number of spatial dimensions and basis functions per field
                bool mSpaceDimEval = true;
                uint mNumSpaceDims;
                bool mNumBasesEval = true;
                uint mNumBasesPerField;

                // evaluation flags for variable and test function sets
                bool mVarSetEval = true;
                bool mTestFuncSetEval = true;

                // evaluation flags for flux matrices
                bool mFluxAMatEval = true;
                bool mFluxADofMatEval = true;
                bool mFluxKMatEval = true;
                bool mFluxKDofMatEval = true;

                // vectors and matrices containing field variables and their spatial derivatives
                Matrix< DDRMat > mY;
                Matrix< DDRMat > mdYdt;
                Matrix< DDRMat > mdYdx;
                Matrix< DDRMat > md2Ydx2;

                // matrices containing test functions and their spatial derivatives
                Matrix< DDRMat > mW;
                Matrix< DDRMat > mdWdt;
                moris::Cell< Matrix< DDRMat > > mdWdx;
                moris::Cell< Matrix< DDRMat > > md2Wdx2;

                // cells of matrices containing the A flux matrices and their DoF derivatives
                moris::Cell< Matrix< DDRMat > > mA;
                moris::Cell< moris::Cell< Matrix< DDRMat > > > mADOF;
                
                // cells of matrices containing the K flux matrices and their DoF derivatives
                moris::Cell< moris::Cell< Matrix< DDRMat > > > mK;
                moris::Cell< moris::Cell< Matrix< DDRMat > > > mKDOF;   

                // multiplication matrices for condensed tensors
                const Matrix< DDRMat > mMultipMat2D = { 
                        { 1.0, 0.0, 0.0 },
                        { 0.0, 1.0, 0.0 },
                        { 0.0, 0.0, 2.0 } };  

                const Matrix< DDRMat > mMultipMat3D = { 
                        { 1.0, 0.0, 0.0, 0.0, 0.0, 0.0 },
                        { 0.0, 1.0, 0.0, 0.0, 0.0, 0.0 },
                        { 0.0, 0.0, 1.0, 0.0, 0.0, 0.0 },
                        { 0.0, 0.0, 0.0, 2.0, 0.0, 0.0 },
                        { 0.0, 0.0, 0.0, 0.0, 2.0, 0.0 },
                        { 0.0, 0.0, 0.0, 0.0, 0.0, 2.0 } };     

                //------------------------------------------------------------------------------

            public:

                // local property enums
                enum class IWG_Property_Type
                {
                    DYNAMIC_VISCOSITY,
                    THERMAL_CONDUCTIVITY,
                    BODY_FORCE,
                    BODY_HEAT_LOAD,
                    MAX_ENUM
                };

                // local constitutive enums
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

                // local SP enums
                // FIXME: generic SP only for testing purposes for strong form
                enum class IWG_Stabilization_Type
                {
                    GENERIC,
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

            private:

                //------------------------------------------------------------------------------
                /**
                 * get number of spatial dimensions
                 * @param[ out ] mNumSpaceDims number of spatial dimensions
                 */
                uint num_space_dims();

                //------------------------------------------------------------------------------
                /**
                 * get number of basis functions per field
                 * @param[ out ] mNumBasesPerField number of basis functions per field
                 */
                uint num_bases();

                //------------------------------------------------------------------------------
                /**
                 * assemble the state variables and their derivatives
                 * into vectors and matrices
                 */
                void assemble_variable_set();

                //------------------------------------------------------------------------------
                /**
                 * get vector of state state variables
                 * @param[ out ] Y  Vector containing state variables
                 */
                const Matrix< DDRMat > & Y();

                //------------------------------------------------------------------------------
                /**
                 * get vector of rate of change of state variables
                 * @param[ out ] dYdt  Vector containing time rate of change of state variables
                 */
                const Matrix< DDRMat > & dYdt();

                //------------------------------------------------------------------------------
                /**
                 * get a matrix of spatial derivative of test functions W for all dof types
                 * @param[ out ] dYdx              spatial derivatives of the state variables
                 */
                const Matrix< DDRMat > & dYdx();

                //------------------------------------------------------------------------------
                /**
                 * get a matrix of spatial derivative of test functions W for all dof types
                 * @param[ out ] d2Ydx2              spatial derivatives of the state variables
                 */
                const Matrix< DDRMat > & d2Ydx2();

                //------------------------------------------------------------------------------
                /**
                 * assemble the DoF derivatives of the state variables and their spatial 
                 * derivatives. These are equal to the test functions
                 */
                void assemble_test_function_set();        

                //------------------------------------------------------------------------------
                /**
                 * get matrix of test functions W for all dof types
                 * @param[ out ] W  Matrix of test functions
                 */
                const Matrix< DDRMat > & W();

                //------------------------------------------------------------------------------
                /**
                 * get of rate of change of test functions
                 * @param[ out ] dWdt  time rate of change of test functions
                 */
                const Matrix< DDRMat > & dWdt();

                //------------------------------------------------------------------------------
                /**
                 * get a matrix of spatial derivative of test functions W for all dof types
                 * @param[ in ]  aSpatialDirection index of the spatial dimensions
                 * @param[ out ] dWdx              spatial derivative of the test functions
                 */
                const Matrix< DDRMat > & dWdx( const uint aSpatialDirection );

                //------------------------------------------------------------------------------
                /**
                 * get a matrix of a second spatial derivative of the test functions W for all dof types
                 * @param[ in ]  aI index of the first spatial dimension for derivative
                 * @param[ in ]  aJ index of the second spatial dimension for derivative
                 * @param[ out ] d2Wdx2   matrix of second spatial derivative of the test functions
                 */
                const Matrix< DDRMat > & d2Wdx2( const uint aI, const uint aJ );

                //------------------------------------------------------------------------------
                /**
                 * compute the residual strong form
                 * @param[ in ] aRM a matrix to fill with RM
                 */
                void compute_residual_strong_form( Matrix< DDRMat > & aRM );

                //------------------------------------------------------------------------------
                /**
                 * compute the residual strong form
                 * @param[ in ] aJM       a matrix to fill with dRMdDof
                 */
                void compute_jacobian_strong_form( Matrix< DDRMat > & aJM );   

                //------------------------------------------------------------------------------
                /**
                 * evaluate the A flux matrices
                 */
                void eval_A_matrices();     

                //------------------------------------------------------------------------------
                /**
                 * evaluate the DoF derivatives of the A flux matrices for the Jacobians
                 */
                void eval_A_DOF_matrices();

                //------------------------------------------------------------------------------
                /**
                 * get the K flux matrices
                 * @param[ in ]  aI  first index
                 * @param[ in ]  aJ  second index
                 * @param[ out ] Kij K-matrix
                 */
                const Matrix< DDRMat > & K( uint aI, uint aJ );     

                //------------------------------------------------------------------------------
                // FIXME provided directly by the field interpolator?
                /**
                 * compute the term dnNdtn
                 * @param[ in ] adnNdtn a matrix to fill with dnNdtn
                 */
                void compute_dnNdtn( Matrix< DDRMat > & adnNdtn );

                //------------------------------------------------------------------------------
                /**
                 * get the multiplication matrix for condensed tensors
                 * @param[ out ] mMultipMat multiplication matrix for condensed tensors
                 */
                const Matrix< DDRMat > & MultipMat();       

        };
        //------------------------------------------------------------------------------
    } /* namespace fem */
} /* namespace moris */

#endif /* SRC_FEM_CL_FEM_IWG_COMPRESSIBLE_NS_VELOCITY_BULK_HPP_ */
