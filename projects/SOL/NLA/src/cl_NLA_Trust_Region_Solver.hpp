/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_NLA_Arc_Length.hpp
 *
 */

#ifndef PROJECTS_SOL_NLA_SRC_CL_NLA_TRUST_REGION_SOLVER_HPP_
#define PROJECTS_SOL_NLA_SRC_CL_NLA_TRUST_REGION_SOLVER_HPP_
#include "cl_SOL_Dist_Matrix.hpp"
#include "moris_typedefs.hpp"
#include "cl_NLA_Nonlinear_Algorithm.hpp"
#include "cl_TSA_Time_Solver_Algorithm.hpp"

namespace moris
{
class Dist_Vector;


namespace tsa
{
class Time_Solver_Algorithm;
}
namespace dla
{
    class Linear_Solver_Algorithm;
}
namespace NLA
{
    struct Trust_Region_Minimizer_Output
    {
            Matrix< DDRMat > mZ;    // CG solve vector (Newton point)
            Matrix< DDRMat > mQ;    // CG solve vector (Cauchy point)
            std::string      mStepType;    // Step type to know if boundary point or interior point or exterior
    };
    class Trust_Region_Solver : public Nonlinear_Algorithm
    {

    public:
        /**
         * @brief Constructor for Arc Length
         *
         */
        Trust_Region_Solver();

        Trust_Region_Solver( const Parameter_List& aParameterlist );

        Trust_Region_Solver( dla::Linear_Solver * aLinSolver );

        ~Trust_Region_Solver();

        /**
         * @brief Call to solve the nonlinear system
         *
         * @param[in] aNonlinearProblem Nonlinear problem
         */

        void solver_nonlinear_system( Nonlinear_Problem * aNonlinearProblem );

        void get_full_solution( moris::Matrix< DDRMat > & LHSValues );

        void get_solution( moris::Matrix< DDRMat > & LHSValues );

        void extract_my_values(
                    const moris::uint&                      aNumIndices,
                    const moris::Matrix< DDSMat >&          aGlobalBlockRows,
                    const moris::uint&                      aBlockRowOffsets,
                    Vector< moris::Matrix< DDRMat > >& LHSValues ) override;

        //void set_my_time_solver_algorithm( std::shared_ptr< tsa::Time_Solver_Algorithm > aMyTimeSolverAlgorithm );

        
        

    private:
        /**
         * @brief Call for solve of linear system
         *
         * @param[in] aIter       Number of iterations
         * @param[in] aHardBreak  Flag for HardBreak
         */
        void solve_linear_system( moris::sint & aIter,
                                  bool        & aHardBreak);

        //------------------------------------------------------------------------------
        //--------------------------- vectors and matrices -----------------------------
        //------------------------------------------------------------------------------
        sol::Dist_Matrix* mJac;

        sol::Dist_Vector* mGlobalRHS;
        //------------------------------------------------------------------------------

        std::shared_ptr< tsa::Time_Solver_Algorithm > mMyTimeSolverAlgorithm = nullptr;

        Trust_Region_Minimizer_Output solver_trust_region_minimizer_system( Matrix< DDRMat > &aX,
                                                                Matrix< DDRMat > &aR,                                                                
                                                                real &aTrSize,
                                                                Nonlinear_Problem *aNonlinearProblem )   ;

        sol::Dist_Vector* preconditioned_project_to_boundary( sol::Dist_Vector* aZ,
                                                              sol::Dist_Vector* aD,
                                                              real          &aTrSize,
                                                              real          &aZz     );
       
        sol::Dist_Vector*   dogleg_step(sol::Dist_Vector* aZ,
                                        sol::Dist_Vector* aQ,
                                        real        &aTrSize);
                                                  
        //void initialize_variables( Nonlinear_Problem *aNonlinearProblem );     
                

    protected:

//        friend class Nonlinear_Problem;

    };
}
}

#endif /* PROJECTS_SOL_NLA_SRC_CL_NLA_ARC_LENGTH_HPP_ */

