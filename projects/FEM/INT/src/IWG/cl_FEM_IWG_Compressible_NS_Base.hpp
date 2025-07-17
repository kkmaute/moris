/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_FEM_IWG_Compressible_NS_Base.hpp
 *
 */

#ifndef SRC_FEM_CL_FEM_IWG_COMPRESSIBLE_NS_BASE_HPP_
#define SRC_FEM_CL_FEM_IWG_COMPRESSIBLE_NS_BASE_HPP_

#include "moris_typedefs.hpp"    //MRS/COR/src
#include "cl_Vector.hpp"         //MRS/CNT/src

#include "cl_Matrix.hpp"          //LINALG/src
#include "linalg_typedefs.hpp"    //LINALG/src

#include "cl_FEM_IWG.hpp"    //FEM/INT/src

namespace moris::fem
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
        fem::Variable_Set mVariableSet    = fem::Variable_Set::PRESSURE_PRIMITIVE;
        MSI::Dof_Type     mFirstDof       = MSI::Dof_Type::P;
        MSI::Dof_Type     mVectorDof      = MSI::Dof_Type::VX;
        MSI::Dof_Type     mLastDof        = MSI::Dof_Type::TEMP;
        uint              mFirstDofIndex  = 0;
        uint              mVectorDofIndex = 1;
        uint              mLastDofIndex   = 2;

        // List of accepted variable sets
        Vector< Vector< MSI::Dof_Type > > mConservativeVars = {
            { MSI::Dof_Type::P },
            { MSI::Dof_Type::MX },
            { MSI::Dof_Type::E }
        };
        Vector< Vector< MSI::Dof_Type > > mDensityPrimitiveVars = {
            { MSI::Dof_Type::RHO },
            { MSI::Dof_Type::VX },
            { MSI::Dof_Type::TEMP }
        };
        Vector< Vector< MSI::Dof_Type > > mPressurePrimitiveVars = {
            { MSI::Dof_Type::P },
            { MSI::Dof_Type::VX },
            { MSI::Dof_Type::TEMP }
        };
        Vector< Vector< MSI::Dof_Type > > mEntropyVars = {
            { MSI::Dof_Type::EVP },
            { MSI::Dof_Type::EVX },
            { MSI::Dof_Type::EVT }
        };

        // evaluation flags for variable set
        bool mYEval      = true;
        bool mdYdtEval   = true;
        bool mdYdxEval   = true;
        bool md2Ydx2Eval = true;

        // evaluation flags for test function sets
        bool mWEval      = true;
        bool mdWdtEval   = true;
        bool mdWdxEval   = true;
        bool md2Wdx2Eval = true;

        bool mWtransEval    = true;
        bool mdWtransdtEval = true;
        bool mdWtransdxEval = true;

        // evaluation flags for flux matrices
        bool mAEval    = true;
        bool mKEval    = true;
        bool mKijiEval = true;

        // evaluation flags for the body load coefficient matrix
        bool mCEval    = true;
        bool mdCdYEval = true;

        // vectors and matrices containing field variables and their spatial derivatives
        Matrix< DDRMat >           mY;
        Matrix< DDRMat >           mdYdt;
        Vector< Matrix< DDRMat > > mdYdx;
        Vector< Matrix< DDRMat > > md2Ydx2;

        // matrices containing test functions and their spatial derivatives
        Matrix< DDRMat >           mW;
        Matrix< DDRMat >           mdWdt;
        Vector< Matrix< DDRMat > > mdWdx;
        Vector< Matrix< DDRMat > > md2Wdx2;

        Matrix< DDRMat >           mWtrans;
        Matrix< DDRMat >           mdWtransdt;
        Vector< Matrix< DDRMat > > mdWtransdx;

        // cells of matrices containing the flux matrices
        Vector< Matrix< DDRMat > >           mA;
        Vector< Vector< Matrix< DDRMat > > > mK;
        Vector< Matrix< DDRMat > >           mKiji;

        // storage vars for the body load coefficient matrix
        Matrix< DDRMat >           mC;
        Matrix< DDRMat >           mdCdYVR;
        Vector< Matrix< DDRMat > > mdCdY;

        // matrix temporarily storing assembly indices
        Matrix< DDSMat > mAssemblyIndices;

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
        ~IWG_Compressible_NS_Base() override{};

        //------------------------------------------------------------------------------
        /**
         * reset eval flags specific to this base IWG
         */
        void reset_spec_eval_flags() override;

        //------------------------------------------------------------------------------
        /**
         * reset eval flags specific to the children comp flow IWGs, empty in base class
         */
        virtual void reset_child_eval_flags() {};

        //------------------------------------------------------------------------------
        /**
         * compute the residual in children, doesn't do anything in base
         * @param[ in ] aWStar weight associated to the evaluation point
         */
        void
        compute_residual( real aWStar ) override
        {
            MORIS_ERROR( false, "IWG_Compressible_NS_Base::compute_residual() - not implemented in comp flow base class" );
        };

        //------------------------------------------------------------------------------
        /**
         * compute the jacobian in children, doesn't do anything in base
         * @param[ in ] aWStar weight associated to the evaluation point
         */
        void
        compute_jacobian( real aWStar ) override
        {
            MORIS_ERROR( false, "IWG_Compressible_NS_Base::compute_jacobian() - not implemented in comp flow base class" );
        };

        //------------------------------------------------------------------------------
        /**
         * compute the residual and the jacobian in children, doesn't do anything in base
         * @param[ in ] aWStar weight associated to the evaluation point
         */
        void
        compute_jacobian_and_residual( real aWStar ) override
        {
            MORIS_ERROR( false, "IWG_Compressible_NS_Base::compute_jacobian_and_residual() - not implemented in comp flow base class" );
        };

        //------------------------------------------------------------------------------
        /**
         * compute the derivative of the residual wrt design variables in children, doesn't do anything in base
         * @param[ in ] aWStar weight associated to the evaluation point
         */
        void
        compute_dRdp( real aWStar ) override
        {
            MORIS_ERROR( false, "IWG_Compressible_NS_Base::compute_jacobian_and_residual() - not implemented in comp flow base class" );
        };

        //------------------------------------------------------------------------------
        /**
         * get the dof type for the state variables by the index in the standardized state var lists
         * @param[ out ] PrimaryDofType MSI dof type of the standardized state variable list
         */
        MSI::Dof_Type
        get_primary_state_var( uint aIndex )
        {
            // make sure the index makes sense
            MORIS_ASSERT( aIndex <= 2, "get_primary_state_var() - index must be 0, 1, or 2" );

            // return Dof Type
            return mPressurePrimitiveVars( aIndex )( 0 );
        };

        // private:

        //------------------------------------------------------------------------------
        /**
         * assemble the standardized element residual in into the set residual
         * @param[ in ]  aStdRes  standardized element residual
         */
        void assemble_residual( const Matrix< DDRMat >& aStdRes );

        //------------------------------------------------------------------------------
        /**
         * assemble the standardized element Jacobian in into the set Jacobian
         * @param[ in ]  aStdJac  standardized element Jacobian
         */
        void assemble_jacobian( const Matrix< DDRMat >& aStdJac );

        //------------------------------------------------------------------------------
        /**
         * get the assembly indices associated with the various dof types for the
         * standardized elemental residual and jacobian
         * @param[ in ]   aDofType  DoF-Type for which the assembly indices are requested
         * @param[ out ]  tIndices  Vector of two entries containing the start and end indices
         */
        Matrix< DDSMat >& get_assembly_indices( const MSI::Dof_Type aDofType );

        //------------------------------------------------------------------------------
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
        const Matrix< DDRMat >& Y();

        //------------------------------------------------------------------------------
        /**
         * get vector of rate of change of state variables
         * @param[ out ] dYdt  Vector containing time rate of change of state variables
         */
        const Matrix< DDRMat >& dYdt();

        //------------------------------------------------------------------------------
        /**
         * get a matrix of spatial derivative of test functions W for all dof types
         * @param[ in ]  aJ index of the spatial dimension for derivative
         * @param[ out ] dYdx  spatial derivatives of the state variables
         */
        const Matrix< DDRMat >& dYdx( const uint aJ );

        //------------------------------------------------------------------------------
        /**
         * get a matrix of spatial derivative of test functions W for all dof types
         * @param[ in ]  aI index of the first spatial dimension for derivative
         * @param[ in ]  aJ index of the second spatial dimension for derivative
         * @param[ out ] d2Ydx2  spatial derivatives of the state variables
         */
        const Matrix< DDRMat >& d2Ydx2( const uint aI, const uint aJ );

        //------------------------------------------------------------------------------
        //------------------------------------------------------------------------------
        /**
         * get matrix of test functions W for all dof types
         * @param[ out ] W  Matrix of test functions
         */
        const Matrix< DDRMat >& W();
        const Matrix< DDRMat >& W_trans();

        //------------------------------------------------------------------------------
        /**
         * get of rate of change of test functions
         * @param[ out ] dWdt  time rate of change of test functions
         */
        const Matrix< DDRMat >& dWdt();
        const Matrix< DDRMat >& dWdt_trans();

        //------------------------------------------------------------------------------
        /**
         * get a matrix of spatial derivative of test functions W for all dof types
         * @param[ in ]  aSpatialDirection index of the spatial dimensions
         * @param[ out ] dWdx              spatial derivative of the test functions
         */
        const Matrix< DDRMat >& dWdx( const uint aSpatialDirection );
        const Matrix< DDRMat >& dWdx_trans( const uint aSpatialDirection );

        //------------------------------------------------------------------------------
        /**
         * get a matrix of a second spatial derivative of the test functions W for all dof types
         * @param[ in ]  aI index of the first spatial dimension for derivative
         * @param[ in ]  aJ index of the second spatial dimension for derivative
         * @param[ out ] d2Wdx2   matrix of second spatial derivative of the test functions
         */
        const Matrix< DDRMat >& d2Wdx2( const uint aI, const uint aJ );

        //------------------------------------------------------------------------------
        //------------------------------------------------------------------------------
        /**
         * get the A flux matrices
         * @param[ in ]  aK index
         * @param[ out ] mA A-matrix
         */
        const Matrix< DDRMat >& A( const uint aK );

        //------------------------------------------------------------------------------
        /**
         * get the K flux matrices
         * @param[ in ]  aI  first index
         * @param[ in ]  aJ  second index
         * @param[ out ] Kij K-matrix
         */
        const Matrix< DDRMat >& K( const uint aI, const uint aJ );

        //------------------------------------------------------------------------------
        /**
         * get the spatial derivatives of the K flux matrices
         * @param[ in ]  aJ  index
         * @param[ out ] Kij_i spatial derivatives of the K flux matrices
         */
        const Matrix< DDRMat >& Kiji( const uint aJ );

        //------------------------------------------------------------------------------
        /**
         * get the coefficient matrix for the body loads
         * @param[ out ] mC coefficient matrix
         */
        const Matrix< DDRMat >& C();

        //------------------------------------------------------------------------------
        /**
         * get the state variable derivative of the coefficient matrix for the body loads
         * pre multiplied with a vector from the right
         * @param[ in ]  aVR      vector for pre-multiplication
         * @param[ out ] mdCdYVR  coefficient matrix
         */
        const Matrix< DDRMat >& dCdY_VR( const Matrix< DDRMat >& aVR );

        //------------------------------------------------------------------------------
        /**
         * get the state variable derivative of the coefficient matrix for the body loads
         * @param[ in ]  aYind  state variable index
         * @param[ out ] mC     coefficient matrix
         */
        const Matrix< DDRMat >& dCdY( const uint aYind );

        //------------------------------------------------------------------------------
        //------------------------------------------------------------------------------
        /**
         * get the multiplication matrix for condensed tensors
         * @param[ out ] mMultipMat multiplication matrix for condensed tensors
         */
        const Matrix< DDRMat >& MultipMat();

        //------------------------------------------------------------------------------
        //------------------------------------------------------------------------------
        /**
         * get elemental residual and jacobian filled with location indices
         * --- for debugging
         */
        Matrix< DDRMat > get_elemental_index_vector();
        Matrix< DDRMat > get_elemental_index_matrix();
    };
    //------------------------------------------------------------------------------
}    // namespace moris::fem

#endif /* SRC_FEM_CL_FEM_IWG_COMPRESSIBLE_NS_VELOCITY_BULK_HPP_ */
