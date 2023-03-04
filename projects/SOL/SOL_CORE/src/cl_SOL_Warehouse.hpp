/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_SOL_Warehouse.hpp
 *
 */
#ifndef MORIS_DISTLINALG_CL_SOL_WAREHOUSE_HPP_
#define MORIS_DISTLINALG_CL_SOL_WAREHOUSE_HPP_

// MORIS header files.
// MORIS header files.
#include "typedefs.hpp"    // CON/src
#include "cl_Cell.hpp"
#include <memory>

#include "cl_MSI_Dof_Type_Enums.hpp"
#include "cl_SOL_Enums.hpp"

#include "cl_NLA_Nonlinear_Solver_Enums.hpp"    //CON/src
#include "cl_TSA_Time_Solver_Enums.hpp"         //CON/src

#include "cl_Param_List.hpp"    //CON/src

namespace moris
{
    class Solver_Interface;
    class Library_IO;
    namespace dla
    {
        class Linear_Solver_Algorithm;
        class Linear_Solver;

    }    // namespace dla
    namespace NLA
    {
        class Nonlinear_Algorithm;
        class Nonlinear_Solver;
    }    // namespace NLA
    namespace tsa
    {
        class Time_Solver_Algorithm;
        class Time_Solver;
    }    // namespace tsa
    namespace sol
    {
        //--------------------------------------------------------------------------------------------------------

        /**
         * @brief Default function for the output criterion
         *
         * @param aTimeSolver
         * @return true
         * @return false
         */
        inline bool
        Default_Output_Criterion( moris::tsa::Time_Solver* aTimeSolver )
        {
            return true;
        }

        //--------------------------------------------------------------------------------------------------------

        class SOL_Warehouse
        {
          private:
            // Pointer to the solver interface
            moris::Solver_Interface* mSolverInterface;

            // List of solver algorithms and solvers
            Cell< std::shared_ptr< dla::Linear_Solver_Algorithm > > mLinearSolverAlgorithms;
            Cell< dla::Linear_Solver* >                             mLinearSolvers;

            // List of nonlinear solver algorithms and solvers
            Cell< std::shared_ptr< NLA::Nonlinear_Algorithm > > mNonlinearSolverAlgorithms;
            Cell< NLA::Nonlinear_Solver* >                      mNonlinearSolvers;

            // List of time solver algorithms and solvers
            Cell< std::shared_ptr< tsa::Time_Solver_Algorithm > > mTimeSolverAlgorithms;
            Cell< tsa::Time_Solver* >                             mTimeSolvers;

            // Parameter list for (0) Linear Algorithm (1) Linear Solver (2) nonlinear Algorithm (3)
            // Nonlinear Solver (4) TimeSolver Algorithm (5) Time Solver (6) Warehouse
            moris::Cell< moris::Cell< moris::ParameterList > > mParameterlist;

            // pointer to dynamically linked library
            std::shared_ptr< Library_IO > mLibrary = nullptr;

            // TPL type
            enum sol::MapType mTPLType = sol::MapType::Epetra;

            // save operator to matlab string
            std::string mOperatorToMatlab = std::string( "" );

            // save final solution vector to file string
            std::string mSaveFinalSolVecToFile = std::string( "" );

            // load final solution vector from file. This option will skip assembly and solve
            std::string mLoadSolVecFromFile = std::string( "" );

            // save final adjoint vector to file string
            std::string mSaveFinalAdjointVecToFile = std::string( "" );

            // load initial guess solution vector from file
            std::string mFilenameInitialGuess = std::string( "" );

            //--------------------------------------------------------------------------------------------------------

            void create_linear_solver_algorithms();

            //--------------------------------------------------------------------------------------------------------

            void create_linear_solvers();

            //--------------------------------------------------------------------------------------------------------

            void create_nonlinear_solver_algorithms();

            //--------------------------------------------------------------------------------------------------------

            void create_nonlinear_solvers();

            //--------------------------------------------------------------------------------------------------------

            void create_time_solver_algorithms();

            //--------------------------------------------------------------------------------------------------------

            void create_time_solvers();

            //--------------------------------------------------------------------------------------------------------

            void get_default_secondary_dof_types(
                    Cell< Cell< MSI::Dof_Type > >&        aCellOfCellsSecDofTypes,
                    Cell< Cell< MSI::Dof_Type > > const & aCellOfCellDofTypes );

            //--------------------------------------------------------------------------------------------------------

            /**
             * @brief Member function called in finalize(). Calculates the downward dependencies based on information received from the nonlinear solver managers
             */
            // void create_solver_manager_dependencies();

            //--------------------------------------------------------------------------------------------------------
            /**
             * @brief Member function called in finalize(). Creates all maps far all nonlinear solver managers
             */
            // void create_maps();

            //--------------------------------------------------------------------------------------------------------

          public:
            //--------------------------------------------------------------------------------------------------------
            /**
             * @brief Constructor.
             *
             * @param[in] aSolverInterface Pointer to the solver interface
             */
            SOL_Warehouse( moris::Solver_Interface* aSolverInterface )
                    : mSolverInterface( aSolverInterface ){};

            SOL_Warehouse( moris::Solver_Interface* aSolverInterface,
                    std::shared_ptr< Library_IO >   aLibrary )
                    : mSolverInterface( aSolverInterface )
                    , mLibrary( aLibrary ){};

            SOL_Warehouse(){};

            //--------------------------------------------------------------------------------------------------------

            /**
             * @brief Destructor.
             */
            ~SOL_Warehouse();

            //--------------------------------------------------------------------------------------------------------

            void
            set_parameterlist( moris::Cell< moris::Cell< moris::ParameterList > > aParameterlist )
            {
                mParameterlist = aParameterlist;
            };

            //--------------------------------------------------------------------------------------------------------

            /**
             * @brief Returns a pointer to the solver interface.
             */
            // void set_nonlinear_solver_managers( Nonlinear_Solver * aNonlinearSolverManager );

            //--------------------------------------------------------------------------------------------------------

            void
            set_solver_interface( moris::Solver_Interface* aSolverInterface )
            {
                mSolverInterface = aSolverInterface;
            };

            //--------------------------------------------------------------------------------------------------------

            moris::Solver_Interface*
            get_solver_interface()
            {
                return mSolverInterface;
            };

            //--------------------------------------------------------------------------------------------------------

            void initialize();

            //--------------------------------------------------------------------------------------------------------

            tsa::Time_Solver*
            get_main_time_solver()
            {
                return mTimeSolvers( 0 );
            };

            //--------------------------------------------------------------------------------------------------------

            enum sol::MapType
            get_tpl_type()
            {
                return mTPLType;
            }

            //--------------------------------------------------------------------------------------------------------

            const std::string&
            get_output_to_matlab_string()
            {
                return mOperatorToMatlab;
            }

            //--------------------------------------------------------------------------------------------------------

            const std::string&
            get_save_final_sol_vec_to_file()
            {
                return mSaveFinalSolVecToFile;
            }

            //--------------------------------------------------------------------------------------------------------

            const std::string&
            get_save_final_adjoint_vec_to_file()
            {
                return mSaveFinalAdjointVecToFile;
            }

            //--------------------------------------------------------------------------------------------------------

            const std::string&
            get_load_sol_vec_from_file()
            {
                return mLoadSolVecFromFile;
            }

            //-------------------------------------------------------------------------------

            const std::string&
            get_TSA_initial_guess_input_filename()
            {
                return mFilenameInitialGuess;
            }

            //--------------------------------------------------------------------------------------------------------

            void
            set_tpl_type( enum sol::MapType aTPLType )
            {
                mTPLType = aTPLType;
            }

            //--------------------------------------------------------------------------------------------------------

            /**
             * @brief Finalize call. Calculates dependencies, maps and ships pointers after all information is received.
             */
            // void finalize();

            //--------------------------------------------------------------------------------------------------------

            /**
             * @brief Calls finalize() and initializes the solve for the highest system
             */
            // void solve();

            //--------------------------------------------------------------------------------------------------------
            /**
             * @brief Returns the index of a certain nonlinear solver manager
             *
             * @param[in] aSolverManagerIndex The index of the asking nonlinear solver manager
             * @param[in] aDofTypeListIndex The index of dof type list on the asking nonlinear solver manager
             */
            // moris::sint get_nonlinear_solver_manager_index( const moris::sint aSolverManagerIndex,
            //                                                 const moris::sint aDofTypeListIndex );

            //--------------------------------------------------------------------------------------------------------
            /**
             * @brief Returns the free map for the asking nonlinear solver manager
             *
             * @param[in] aSolverManagerIndex The index of the asking nonlinear solver manager
             */
            // Dist_Map * get_list_of_maps( const moris::sint aSolverManagerIndex );

            //--------------------------------------------------------------------------------------------------------
            /**
             * @brief Returns the full map for the asking nonlinear solver manager
             *
             */
            // Dist_Map * get_full_maps(  );

            //--------------------------------------------------------------------------------------------------------
            /**
             * @brief Returns the nonlinear solver manager list.
             */
            // moris::Cell< Nonlinear_Solver * > & get_nonlinear_solver_manager_list(){ return mListNonlinearSolverManagers; };

            //--------------------------------------------------------------------------------------------------------
            /**
             * @brief Returns a pointer to the full vector
             */
            // Dist_Vector * get_full_vector(){ return mFullVector; };
        };


    }    // namespace sol
}    // namespace moris
#endif /* MORIS_DISTLINALG_CL_SOL_WAREHOUSE_HPP_ */
