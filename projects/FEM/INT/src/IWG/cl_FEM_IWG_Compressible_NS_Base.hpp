/*
 * cl_FEM_IWG_Compressible_NS_Base.hpp
 *
 *  Created on: Apr 23, 2021
 *      Author: wunsch
 *
 * Base IWG for unified NS equations using flux matrices
 * operating on all residual DoF types simultaneously
 * 
 */

#ifndef SRC_FEM_CL_FEM_IWG_COMPRESSIBLE_NS_BASE_HPP_
#define SRC_FEM_CL_FEM_IWG_COMPRESSIBLE_NS_BASE_HPP_

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

        class IWG_Compressible_NS_Base : public IWG
        {
                //------------------------------------------------------------------------------

            private:

                // default dof types (FIXME: temporary)
                MSI::Dof_Type mDofDensity     = MSI::Dof_Type::RHO;
                MSI::Dof_Type mDofPressure    = MSI::Dof_Type::P;
                MSI::Dof_Type mDofVelocity    = MSI::Dof_Type::VX;
                MSI::Dof_Type mDofTemperature = MSI::Dof_Type::TEMP;

                // Define Variable Set (FIXME: only primitive vars for now, needs to have a set function later)
                fem::Variable_Set mVariableSet = fem::Variable_Set::PRESSURE_PRIMITIVE;
                MSI::Dof_Type mFirstDof   = MSI::Dof_Type::P;
                MSI::Dof_Type mVectorDof  = MSI::Dof_Type::VX;
                MSI::Dof_Type mLastDof    = MSI::Dof_Type::TEMP;
                uint mFirstDofIndex   = 0;
                uint mVectorDofIndex  = 1;
                uint mLastDofINdex    = 2;

                // List of accepted variable sets
                moris::Cell< moris::Cell< MSI::Dof_Type > > mConservativeVars = { 
                        { MSI::Dof_Type::P }, 
                        { MSI::Dof_Type::MX }, 
                        { MSI::Dof_Type::E } };
                moris::Cell< moris::Cell< MSI::Dof_Type > > mDensityPrimitiveVars = { 
                        { MSI::Dof_Type::RHO }, 
                        { MSI::Dof_Type::VX }, 
                        { MSI::Dof_Type::TEMP } };
                moris::Cell< moris::Cell< MSI::Dof_Type > > mPressurePrimitiveVars = { 
                        { MSI::Dof_Type::P }, 
                        { MSI::Dof_Type::VX }, 
                        { MSI::Dof_Type::TEMP } };
                moris::Cell< moris::Cell< MSI::Dof_Type > > mEntropyVars = { 
                        { MSI::Dof_Type::EVP }, 
                        { MSI::Dof_Type::EVX }, 
                        { MSI::Dof_Type::EVT } };

                // evaluation flags for variable set
                bool mYEval = true;
                bool mdYdtEval = true;
                bool mdYdxEval = true;
                bool md2Ydx2Eval = true;

                // evaluation flags for test function sets
                bool mWEval = true;
                bool mdWdtEval = true;
                bool mdWdxEval = true;
                bool md2Wdx2Eval = true;

                // evaluation flags for flux matrices
                bool mAEval = true;
                bool mKEval = true;
                bool mKijiEval = true;

                // vectors and matrices containing field variables and their spatial derivatives
                Matrix< DDRMat > mY;
                Matrix< DDRMat > mdYdt;
                moris::Cell< Matrix< DDRMat > > mdYdx;
                moris::Cell< Matrix< DDRMat > > md2Ydx2;

                // matrices containing test functions and their spatial derivatives
                Matrix< DDRMat > mW;
                Matrix< DDRMat > mdWdt;
                moris::Cell< Matrix< DDRMat > > mdWdx;
                moris::Cell< Matrix< DDRMat > > md2Wdx2;

                // cells of matrices containing the flux matrices 
                moris::Cell< Matrix< DDRMat > > mA;
                moris::Cell< moris::Cell< Matrix< DDRMat > > > mK;
                moris::Cell< Matrix< DDRMat > > mKiji;    

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

                //------------------------------------------------------------------------------
                /*
                 *  constructor
                 */
                IWG_Compressible_NS_Base(){};

                //------------------------------------------------------------------------------
                /**
                 * trivial destructor
                 */
                ~IWG_Compressible_NS_Base(){};

                //------------------------------------------------------------------------------
                /**
                 * reset eval flags specific to this base IWG
                 */
                void reset_spec_eval_flags();

                //------------------------------------------------------------------------------
                /**
                 * reset eval flags specific to the children comp flow IWGs, empty in base class
                 */
                virtual void reset_child_eval_flags(){};

                //------------------------------------------------------------------------------
                /**
                 * compute the residual in children, doesn't do anything in base
                 * @param[ in ] aWStar weight associated to the evaluation point
                 */
                virtual void compute_residual( real aWStar )
                {
                    MORIS_ERROR( false, "IWG_Compressible_NS_Base::compute_residual() - not implemented in comp flow base class" );
                };

                //------------------------------------------------------------------------------
                /**
                 * compute the jacobian in children, doesn't do anything in base
                 * @param[ in ] aWStar weight associated to the evaluation point
                 */
                virtual void compute_jacobian( real aWStar )
                {
                    MORIS_ERROR( false, "IWG_Compressible_NS_Base::compute_jacobian() - not implemented in comp flow base class" );
                };

                //------------------------------------------------------------------------------
                /**
                 * compute the residual and the jacobian in children, doesn't do anything in base
                 * @param[ in ] aWStar weight associated to the evaluation point
                 */
                virtual void compute_jacobian_and_residual( real aWStar )
                {
                    MORIS_ERROR( false, "IWG_Compressible_NS_Base::compute_jacobian_and_residual() - not implemented in comp flow base class" );
                };

                //------------------------------------------------------------------------------
                /**
                 * compute the derivative of the residual wrt design variables in children, doesn't do anything in base
                 * @param[ in ] aWStar weight associated to the evaluation point
                 */
                virtual void compute_dRdp( real aWStar )
                {
                    MORIS_ERROR( false, "IWG_Compressible_NS_Base::compute_jacobian_and_residual() - not implemented in comp flow base class" );
                };

                //------------------------------------------------------------------------------
                /**
                 * get the dof type for the state variables by the index in the standardized state var lists
                 * @param[ out ] PrimaryDofType MSI dof type of the standardized state variable list
                 */
                MSI::Dof_Type get_primary_state_var( uint aIndex )
                {
                    // make sure the index makes sense
                    MORIS_ASSERT( aIndex >= 0 and aIndex <= 2, "get_primary_state_var() - index must be 0, 1, or 2" );

                    // return Dof Type
                    return mPressurePrimitiveVars( aIndex )( 0 );
                };

            //private:

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
                 * @param[ in ]  aJ index of the spatial dimension for derivative
                 * @param[ out ] dYdx  spatial derivatives of the state variables
                 */
                const Matrix< DDRMat > & dYdx( const uint aJ );

                //------------------------------------------------------------------------------
                /**
                 * get a matrix of spatial derivative of test functions W for all dof types
                 * @param[ in ]  aI index of the first spatial dimension for derivative
                 * @param[ in ]  aJ index of the second spatial dimension for derivative
                 * @param[ out ] d2Ydx2  spatial derivatives of the state variables
                 */
                const Matrix< DDRMat > & d2Ydx2( const uint aI, const uint aJ ); 

                //------------------------------------------------------------------------------
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
                 * get the spatial derivatives of the K flux matrices
                 * @param[ in ]  aJ  index
                 * @param[ out ] Kij_i spatial derivatives of the K flux matrices
                 */
                const Matrix< DDRMat > & Kiji ( const uint aJ );   

                //------------------------------------------------------------------------------
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
