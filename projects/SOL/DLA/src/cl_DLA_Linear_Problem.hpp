/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_DLA_Linear_System.hpp
 *
 */
#ifndef MORIS_DISTLINALG_CL_DLA_LINEAR_SYSTEM_HPP_
#define MORIS_DISTLINALG_CL_DLA_LINEAR_SYSTEM_HPP_

// MORIS header files.
#include "cl_Matrix.hpp"
#include "linalg_typedefs.hpp"
#include "cl_SOL_Enums.hpp"

namespace moris
{
    //--------------------------------------------------------------------------

    namespace sol
    {
        class Dist_Matrix;
        class Dist_Vector;
        class Dist_Map;
        class SOL_Warehouse;
    }    // namespace sol

    //--------------------------------------------------------------------------

    class Solver_Interface;

    //--------------------------------------------------------------------------

    namespace dla
    {
        class Linear_Problem
        {
          private:

          protected:
            sol::Dist_Matrix* mMat            = nullptr;
            sol::Dist_Matrix* mMassMat        = nullptr;
            sol::Dist_Vector* mFreeVectorLHS  = nullptr;
            sol::Dist_Vector* mPointVectorRHS = nullptr;
            sol::Dist_Vector* mPointVectorLHS = nullptr;
            sol::Dist_Vector* mFullVectorLHS  = nullptr;
            sol::Dist_Map*    mMap            = nullptr;
            sol::Dist_Map*    mMapFree        = nullptr;

            //! Pointer to solver intrface
            Solver_Interface* mSolverInterface = nullptr;

            //! Pointer to solver database
            sol::SOL_Warehouse* mSolverWarehouse = nullptr;

            //! Condition number
            moris::real mCondEstimate;

            //! Distributed vector map type
            enum sol::MapType mTplType = sol::MapType::Epetra;

            //! Pointer to RHS Matrix Type
            std::string mRHSMatType;

            //------------------------------------------------------------------

          public:
            //------------------------------------------------------------------    
            
            /**
             * @brief Construct a new Linear_Problem object
             * 
             */
            Linear_Problem(){};

            //------------------------------------------------------------------
            Linear_Problem( Solver_Interface* aInput )
                    : mMat( nullptr )
                    , mMassMat( nullptr )
                    , mFreeVectorLHS( nullptr )
                    , mMap( nullptr )
                    , mSolverInterface( aInput ){};

            //------------------------------------------------------------------

            virtual ~Linear_Problem(){};

            //------------------------------------------------------------------

            // assemble of residual and jacobian
            void assemble_residual_and_jacobian();

            //------------------------------------------------------------------

            // assemble residual only
            void assemble_residual();

            //------------------------------------------------------------------

            // assemble jacobian only
            void assemble_jacobian();

            //------------------------------------------------------------------
            // compute norm of static residual; used by pseudo-time continuation
            real compute_static_residual_norm();

            //------------------------------------------------------------------

            void assemble_staggered_residual_contribution();

            //------------------------------------------------------------------

            void compute_residual_for_adjoint_solve();

            //------------------------------------------------------------------

            Matrix< DDRMat > compute_residual_of_linear_system();

            //------------------------------------------------------------------

            sol::Dist_Vector*
            get_free_solver_LHS()
            {
                return mPointVectorLHS;
            };

            //------------------------------------------------------------------

            void set_free_solver_LHS( sol::Dist_Vector* aFullSolVector );

            //------------------------------------------------------------------

            sol::Dist_Vector* get_full_solver_LHS();

            //------------------------------------------------------------------

            sol::Dist_Vector*
            get_solver_RHS()
            {
                return mPointVectorRHS;
            };

            //------------------------------------------------------------------

            sol::Dist_Matrix*
            get_matrix()
            {
                return mMat;
            };

            //------------------------------------------------------------------

            sol::Dist_Matrix*
            get_mass_matrix()
            {
                return mMassMat;
            };

            //------------------------------------------------------------------

            Solver_Interface*
            get_solver_input() const
            {
                return mSolverInterface;
            };

            //------------------------------------------------------------------

            virtual void get_solution( moris::Matrix< DDRMat >& LHSValues ) = 0;

            //------------------------------------------------------------------

            virtual moris::sint solve_linear_system() = 0;

            //------------------------------------------------------------------

            /**
             * @brief It uses maps and solver interface to construct the linear system using the graph that is cons
             *
             */
            virtual void
            construct_rhs_matrix();

            //------------------------------------------------------------------

            /**
             * @brief this function uses
             *
             */
            void
            assemble_rhs_matrix();

            //------------------------------------------------------------------

            void
            set_rhs_matrix_type( std::string const & aRHSMatType )
            {
                mRHSMatType = aRHSMatType;
            }

            //------------------------------------------------------------------

            std::string const &
            get_rhs_matrix_type( )
            {
                return mRHSMatType;
            }
        };
    }    // namespace dla
}    // namespace moris
#endif /* MORIS_DISTLINALG_CL_DLA_LINEAR_SYSTEM_HPP_ */
