/*
 * cl_SOL_Warehouse.hpp
 *
 *  Created on: Okt 6, 2018
 *      Author: schmidt
 */
#ifndef MORIS_DISTLINALG_CL_SOL_WAREHOUSE_HPP_
#define MORIS_DISTLINALG_CL_SOL_WAREHOUSE_HPP_

// MORIS header files.
// MORIS header files.
#include "typedefs.hpp" // CON/src
#include "cl_Cell.hpp"
#include <memory>

#include "cl_MSI_Dof_Type_Enums.hpp"
#include "cl_SOL_Enums.hpp"

#include "cl_NLA_Nonlinear_Solver_Enums.hpp"       //CON/src
#include "cl_TSA_Time_Solver_Enums.hpp"       //CON/src

#include "cl_Param_List.hpp"       //CON/src

//#include "cl_MSI_Dof_Type_Enums.hpp"
//#include "cl_NLA_Nonlinear_Solver.hpp"
//#include "cl_SOL_Dist_Map.hpp"

namespace moris
{
class Solver_Interface;

namespace dla
{
    class Linear_Solver_Algorithm;
    class Linear_Solver;

}
namespace NLA
{
    class Nonlinear_Algorithm;
    class Nonlinear_Solver;
}
namespace tsa
{
    class Time_Solver_Algorithm;
    class Time_Solver;

}
namespace sol
{
    class SOL_Warehouse
    {
    private:
        //! Pointer to the solver interface
        Solver_Interface * mSolverInterface;

        Cell< std::shared_ptr< dla::Linear_Solver_Algorithm > > mLinearSolverAlgorithms;
        Cell< dla::Linear_Solver * >                            mLinearSolvers;

        Cell< std::shared_ptr< NLA::Nonlinear_Algorithm > >     mNonlinearSolverAlgoriths;
        Cell< NLA::Nonlinear_Solver * >                         mNonlinearSolvers;

        Cell< std::shared_ptr< tsa::Time_Solver_Algorithm > >   mTimeSolverAlgorithms;
        Cell< tsa::Time_Solver * >                              mTimeSolvers;

        // Parameterlist for (0) Linear Algorithm (1) Linear Solver (2) nonlinear Algorithm (3) Nonlinear Solver (4) TimeSolver Algorithm (5) Time Solver
        moris::Cell< moris::Cell< moris::ParameterList > >             mParameterlist;


//--------------------------------------------------------------------------------------------------------

        /**
         * @brief Member function called in finalize(). Calculates the downward dependencies based on information received from the nonlinear solver managers
         */
        //void create_solver_manager_dependencies();

//--------------------------------------------------------------------------------------------------------
        /**
         * @brief Member function called in finalize(). Creates all maps far all nonlinear solver managers
         */
        //void create_maps();

//--------------------------------------------------------------------------------------------------------

    public:

//--------------------------------------------------------------------------------------------------------
        /**
         * @brief Constructor.
         *
         * @param[in] aSolverInterface Pointer to the solver interface
         */
        SOL_Warehouse( Solver_Interface * aSolverInterface ) : mSolverInterface( aSolverInterface ){};

        SOL_Warehouse(){};

//--------------------------------------------------------------------------------------------------------

        /**
         * @brief Destructor.
         */
        ~SOL_Warehouse(){};

        void set_parameterlist( moris::Cell< moris::Cell< moris::ParameterList > > aParameterlist )
        {
            mParameterlist = aParameterlist;
        };

//--------------------------------------------------------------------------------------------------------

        /**
         * @brief Memeber function to set the nonliner solver managers. The highest level nonliner solver manager has to be on entry 0
         */
        //void set_nonliner_solver_managers( Nonlinear_Solver * aNonlinerSolverManager );

//--------------------------------------------------------------------------------------------------------

        /**
         * @brief Returns a pointer to the solver interface.
         */
        Solver_Interface * get_solver_interface(){ return mSolverInterface; };

        void set_solver_interface( Solver_Interface * aSolverInterface  )
        {
            mSolverInterface = aSolverInterface;
        };

//--------------------------------------------------------------------------------------------------------

        void initialize()
        {
            this->create_linear_solver_algorithms();

            this->create_linear_solvers();

            this->create_nonlinear_solver_algorithms();

            this->create_nonlinear_solvers();

            this->create_time_solver_algorithms();

            this->create_time_solvers();
        };

        tsa::Time_Solver *  get_main_time_solver()
        {
             return mTimeSolvers( 0 );
        };
//--------------------------------------------------------------------------------------------------------

        void create_linear_solver_algorithms();

        void create_linear_solvers();

        void create_nonlinear_solver_algorithms();

        void create_nonlinear_solvers();

        void create_time_solver_algorithms();

        void create_time_solvers();

//--------------------------------------------------------------------------------------------------------

        /**
         * @brief Finalize call. Calculates dependencies, maps and ships pointers after all information is received.
         */
        //void finalize();

//--------------------------------------------------------------------------------------------------------

        /**
         * @brief Calls finalize() and initializes the solve for the highest system
         */
        //void solve();

//--------------------------------------------------------------------------------------------------------
        /**
          * @brief Returns the index of a certain nonliner solver manager
          *
          * @param[in] aSolverManagerIndex The index of the asking nonlinear solver manager
          * @param[in] aDofTypeListIndex The index of dof type list on the asking nonliner solver manager
          */
        //moris::sint get_nonlinear_solver_manager_index( const moris::sint aSolverManagerIndex,
        //                                                const moris::sint aDofTypeListIndex );

//--------------------------------------------------------------------------------------------------------
        /**
          * @brief Returns the free map for the asking nonlinear solver manager
          *
          * @param[in] aSolverManagerIndex The index of the asking nonlinear solver manager
          */
        //Dist_Map * get_list_of_maps( const moris::sint aSolverManagerIndex );

//--------------------------------------------------------------------------------------------------------
        /**
          * @brief Returns the full map for the asking nonlinear solver manager
          *
          */
        //Dist_Map * get_full_maps(  );

//--------------------------------------------------------------------------------------------------------
        /**
          * @brief Returns the nonlinear solver manager list.
          */
        //moris::Cell< Nonlinear_Solver * > & get_nonliner_solver_manager_list(){ return mListNonlinerSolverManagers; };

//--------------------------------------------------------------------------------------------------------
        /**
          * @brief Returns a pointer to the full vector
          */
        //Dist_Vector * get_full_vector(){ return mFullVector; };

    };




    // creates a parameter list with default inputs
    ParameterList create_linear_algorithm_parameter_list_aztec( )
    {
        ParameterList tLinAlgorithmParameterList;

        enum moris::sol::SolverType tType = moris::sol::SolverType::AZTEC_IMPL;

        tLinAlgorithmParameterList.insert( "Solver_Implementation" , static_cast< uint >( tType ) );

        // ASSIGN DEFAULT PARAMETER VALUES
        // AztecOO User Guide, SAND REPORT, SAND2004-3796, https://trilinos.org/oldsite/packages/aztecoo/AztecOOUserGuide.pdf

        // Determine which solver is used
        //options are: AZ_gmres, AZ_gmres_condnum, AZ_cg, AZ_cg_condnum, AZ_cgs, AZ_tfqmr, AZ_bicgstab
        tLinAlgorithmParameterList.insert( "AZ_solver" ,  INT_MAX );

        // Allowable Aztec solver iterations
        tLinAlgorithmParameterList.insert( "AZ_max_iter", INT_MAX   );

        // Allowable Aztec irelative residual
        tLinAlgorithmParameterList.insert( "rel_residual" , 1e-08 );

        // set Az_conv -convergence criteria
        // options are AZ_r0, AZ_rhs, AZ_Anorm, AZ_noscaled, AZ_sol
        tLinAlgorithmParameterList.insert( "AZ_conv" ,  INT_MAX );

        // set Az_diagnostic parameters
        // Set whether or not diagnostics for every linear iteration are printed or not. options are AZ_all, AZ_none
        tLinAlgorithmParameterList.insert( "AZ_diagnostics" ,  INT_MAX );

        // set AZ_output options
        // options are AZ_all, AZ_none, AZ_warnings, AZ_last, AZ_summary
        tLinAlgorithmParameterList.insert( "AZ_output" ,  INT_MAX );

        // Determines the submatrices factored with the domain decomposition algorithms
        // Option to specify with how many rows from other processors each processor’s local submatrix is augmented.
        tLinAlgorithmParameterList.insert( "AZ_overlap" , INT_MAX );

        // Determines how overlapping subdomain results are combined when different processors have computed different values for the same unknown.
        // Options are AZ_standard, AZ_symmetric
        tLinAlgorithmParameterList.insert( "AZ_type_overlap" , INT_MAX );

        // Determines whether RCM reordering will be done in conjunction with domain decomposition incomplete factorizations.
        // Option to enable (=1) or disable (=0) the Reverse Cuthill–McKee (RCM) algorithm to reorder system equations for smaller bandwidth
        tLinAlgorithmParameterList.insert( "AZ_reorder" , INT_MAX );

        // Use preconditioner from a previous Iterate() call
        // Option are AZ_calc, AZ_recalc, AZ_reuse
        tLinAlgorithmParameterList.insert( "AZ_pre_calc" , INT_MAX );

        // Determines  whether  matrix  factorization  information will be kept after this solve
        // for example for preconditioner_recalculation
        tLinAlgorithmParameterList.insert( "AZ_keep_info" , INT_MAX );

        //--------------------------GMRES specific solver parameters--------------------------------------------------------------------------
        // Set AZ_kspace
        // Krylov subspace size for restarted GMRES
        // Setting mKrylovSpace larger improves the robustness, decreases iteration count, but increases memory consumption. For very difficult problems, set it equal to the maximum number of iterations.
        tLinAlgorithmParameterList.insert( "AZ_kspace" ,INT_MAX );

        // Set AZ_orthog
        //AZ_classic or AZ_modified
        tLinAlgorithmParameterList.insert( "AZ_orthog" , INT_MAX );

        // Set AZ_rthresh
        // Parameter used to modify the relative magnitude of the diagonal entries of the matrix that is used to compute any of the incomplete factorization preconditioners
        tLinAlgorithmParameterList.insert( "AZ_rthresh" , -1.0 );

        // Set AZ_athresh
        //Parameter used to modify the absolute magnitude of the diagonal entries of the matrix that is used to compute any of the incomplete factorization preconditioners
        tLinAlgorithmParameterList.insert( "AZ_athresh" , -1.0 );

        //--------------------------Preconsitioner specific parameters--------------------------------------------------------------------------
        // Determine which preconditioner is used
        // Options are AZ_none, AZ_Jacobi, AZ_sym_GS, AZ_Neumann, AZ_ls, AZ_dom_decomp,
        tLinAlgorithmParameterList.insert( "AZ_precond" ,  INT_MAX );

        // Set preconditioner subdomain solve - direct solve or incomplete
        // Options are AZ_lu, AZ_ilut, , AZ_rilu, AZ_bilu, AZ_icc
        tLinAlgorithmParameterList.insert( "AZ_subdomain_solve" ,  INT_MAX );

        // Set preconditioner polynomial order - polynomial preconditioning, Gauss-Seidel, Jacobi
        tLinAlgorithmParameterList.insert( "AZ_poly_ord" ,  INT_MAX );

        // Set drop tolerance - for LU, ILUT
        tLinAlgorithmParameterList.insert( "AZ_drop" ,  -1.0 );

        // Set level of graph fill in - for ilu(k), icc(k), bilu(k)
        tLinAlgorithmParameterList.insert( "AZ_graph_fill" ,  INT_MAX );

        // Set ilut fill
        tLinAlgorithmParameterList.insert( "AZ_ilut_fill" ,  -1.0 );

        // Set Damping or relaxation parameter used for RILU
        tLinAlgorithmParameterList.insert( "AZ_omega" ,  -1.0 );

        // Set Damping or relaxation parameter used for RILU
        tLinAlgorithmParameterList.insert( "Use_ML_Prec" ,  false );

        // Set Damping or relaxation parameter used for RILU
        tLinAlgorithmParameterList.insert( "ML_reuse" ,  false );

        return tLinAlgorithmParameterList;
    }


    ParameterList create_linear_solver_parameter_list()
    {
        ParameterList tLinSolverParameterList;

        tLinSolverParameterList.insert( "DLA_Linear_solver_algorithms" , std::string("0") );

        // Maximal number of linear solver restarts on fail
        tLinSolverParameterList.insert( "DLA_max_lin_solver_restarts" , 0 );

        // Maximal number of linear solver restarts on fail
        tLinSolverParameterList.insert( "DLA_hard_break" , true );

        // Determines if lin solve should restart on fail
        tLinSolverParameterList.insert( "DLA_rebuild_lin_solver_on_fail" , false );

        return tLinSolverParameterList;
    }

ParameterList create_nonlinear_algorithm_parameter_list()
{
    ParameterList tNonLinAlgorithmParameterList;

    enum moris::NLA::NonlinearSolverType NonlinearSolverType = moris::NLA::NonlinearSolverType::NEWTON_SOLVER;

    tNonLinAlgorithmParameterList.insert( "NLA_Solver_Implementation" , static_cast< uint >( NonlinearSolverType ) );

    tNonLinAlgorithmParameterList.insert( "NLA_Linear_solver" , 0 );

    // Allowable Newton solver iterations
    tNonLinAlgorithmParameterList.insert( "NLA_max_iter", 10 );

    // Allowable Newton solver iterations
    tNonLinAlgorithmParameterList.insert( "NLA_restart", 0 );

    // Allowable Newton irelative residual
//    mParameterListNonlinearSolver.insert( "NLA_rel_residual" , 1e-02 );
    tNonLinAlgorithmParameterList.insert( "NLA_rel_residual" , 1e-08 );

    // Desired total residual norm drop
//    mParameterListNonlinearSolver.insert( "NLA_tot_res_norm_drop" , 1e-02 );
    tNonLinAlgorithmParameterList.insert( "NLA_tot_res_norm_drop" , 1e-08 );

    // Desired total residual norm
//    mParameterListNonlinearSolver.insert( "NLA_tot_res_norm" , 1e-2 );
    tNonLinAlgorithmParameterList.insert( "NLA_tot_res_norm" , 1e-9 );

    // Maximal residual norm drop
//    mParameterListNonlinearSolver.insert( "NLA_max_res_norm_drop" , 1e-2 );
    tNonLinAlgorithmParameterList.insert( "NLA_max_res_norm_drop" , 1e-6 );

    // Maximal number of linear solver restarts on fail
    tNonLinAlgorithmParameterList.insert( "NLA_max_lin_solver_restarts" , 0 );

    // Maximal number of linear solver restarts on fail
    tNonLinAlgorithmParameterList.insert( "NLA_relaxation_parameter" , 1.0 );

    // Maximal number of linear solver restarts on fail
    tNonLinAlgorithmParameterList.insert( "NLA_hard_break" , false );

    // Determines if lin solve should restart on fail
    tNonLinAlgorithmParameterList.insert( "NLA_rebuild_lin_solv_on_fail" , false );

    // Determines if lin solve should restart on fail
    tNonLinAlgorithmParameterList.insert( "NLA_rebuild_jacobian" , true );

    // Determines if newton should restart on fail
    tNonLinAlgorithmParameterList.insert( "NLA_rebuild_nonlin_solv_on_fail" , false );

    // Specifying the number of newton retries
    tNonLinAlgorithmParameterList.insert( "NLA_num_nonlin_rebuild_iterations" , 1 );

    // Determines relaxation multiplier
    tNonLinAlgorithmParameterList.insert( "NLA_relaxation_multiplier_on_fail" , 0.5 );

    // Determines newton maxits multiplier
    tNonLinAlgorithmParameterList.insert( "NLA_maxits_multiplier_on_fail" , 2 );

    return tNonLinAlgorithmParameterList;
}

ParameterList create_nonlinear_solver_parameter_list()
{
    ParameterList tNonLinSolverParameterList;

    enum moris::NLA::NonlinearSolverType NonlinearSolverType = moris::NLA::NonlinearSolverType::NEWTON_SOLVER;

    tNonLinSolverParameterList.insert( "NLA_Solver_Implementation" , static_cast< uint >( NonlinearSolverType ) );

    tNonLinSolverParameterList.insert( "NLA_DofTypes" , std::string("UNDEFINED") );

    tNonLinSolverParameterList.insert( "NLA_Sub_Nonlinear_Solver" , std::string("") );

    tNonLinSolverParameterList.insert( "NLA_Nonlinear_solver_algorithms" , std::string("0") );

    // Maximal number of linear solver restarts on fail
    tNonLinSolverParameterList.insert( "NLA_max_non_lin_solver_restarts" , 0 );

    return tNonLinSolverParameterList;
}

ParameterList create_time_solver_algorithm_parameter_list()
{
    ParameterList tTimeAlgorithmParameterList;

    enum moris::tsa::TimeSolverType tType = tsa::TimeSolverType::MONOLITHIC;

    tTimeAlgorithmParameterList.insert( "TSA_Solver_Implementation" , static_cast< uint >( tType ) );


    tTimeAlgorithmParameterList.insert( "TSA_Nonlinear_solver" , 0 );

    // Number of time steps
    tTimeAlgorithmParameterList.insert( "TSA_Num_Time_Steps", 1 );

    // Time Frame
    tTimeAlgorithmParameterList.insert( "TSA_Time_Frame", 1.0 );

    return tTimeAlgorithmParameterList;
}

ParameterList create_time_solver_parameter_list()
{
    ParameterList tTimeParameterList;

    tTimeParameterList.insert( "TSA_Solver_algorithms" , std::string("0") );

    tTimeParameterList.insert( "TSA_DofTypes" , std::string("UNDEFINED") );

    // Maximal number of linear solver restarts on fail
    tTimeParameterList.insert( "TSA_max_time_solver_restarts" , 0 );

    return tTimeParameterList;
}

// creates a parameter list with default inputs
ParameterList create_linear_algorithm_parameter_list( const enum moris::sol::SolverType aSolverType,
                                                      const uint                        aIndex = 0 )
{
    ParameterList tParameterList;

    switch( aSolverType )
    {
    case ( sol::SolverType::AZTEC_IMPL ):
        return create_linear_algorithm_parameter_list_aztec( );
        break;
    case ( sol::SolverType::AMESOS_IMPL ):
		MORIS_ERROR( false, "No implemented yet" );
        break;

    default:
        MORIS_ERROR( false, "Parameterlist for this solver not implemented yet" );
        break;
    }
    return tParameterList;
}

//    // creates a parameter list with default inputs
//    ParameterList create_nonlinear_algorithm_parameter_list( const enum moris::NLA::NonlinearSolverType aSolverType,
//                                                             const uint                        aIndex = 0 )
//    {
//        ParameterList tParameterList;
//
//        switch( aSolverType )
//        {
//        case ( moris::NLA::NonlinearSolverType::NEWTON_SOLVER ):
//            return create_linear_algorithm_parameter_list_aztec( );
//            break;
//        case ( moris::NLA::NonlinearSolverType::NLBGS_SOLVER ):
//    		MORIS_ERROR( false, "No implemented yet" );
//            break;
//
//        default:
//            MORIS_ERROR( false, "Parameterlist for this solver not implemented yet" );
//            break;
//        }
//
//    return tParameterList;
//}
}}
#endif /* MORIS_DISTLINALG_CL_SOL_WAREHOUSE_HPP_ */

