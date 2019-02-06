/*
 * cl_NLA_Nonlinear_Database.hpp
 *
 *  Created on: Okt 6, 2018
 *      Author: schmidt
 */
#ifndef MORIS_DISTLINALG_CL_NLA_NONLINEAR_DATABASE_HPP_
#define MORIS_DISTLINALG_CL_NLA_NONLINEAR_DATABASE_HPP_

// MORIS header files.
#include "typedefs.hpp" // CON/src
#include "cl_Cell.hpp"
#include <memory>
#include "cl_Param_List.hpp"
#include "cl_MSI_Dof_Type_Enums.hpp"

#include "cl_NLA_Nonlinear_Solver.hpp"

#include "cl_Map_Class.hpp"
#include "cl_Matrix_Vector_Factory.hpp"

namespace moris
{
namespace NLA
{
    class Nonlinear_Problem;
    class Nonlinear_Algorithm;
    class Nonlinear_Database
    {
    private:
        //! Pointer to the solver interface
        Solver_Interface * mSolverInterface;

        //! List containing all nonlinear solver managers
        moris::Cell< Nonlinear_Solver * > mListNonlinerSolverManagers;

        //! List containing the downward dependencies of every nonlinear solver manager
        moris::Cell< moris::Matrix< DDSMat > > mListSolverManagerDepenencies;

        //! List of maps for every nonliner solver manager. The entry corresponds to the nonliner solver manager index. The last map is the pull map
        moris::Cell< Map_Class * > mListOfFreeMaps;

        //! Full Vector
        Dist_Vector * mFullVector = nullptr;

        //!  Counter to count the number of nonlinear solver managers
        moris::uint mCallCounter = 0;

//--------------------------------------------------------------------------------------------------------

        /**
         * @brief Member function called in finalize(). Calculates the downward dependencies based on information received from the nonlinear solver managers
         */
        void create_solver_manager_dependencies();

//--------------------------------------------------------------------------------------------------------
        /**
         * @brief Member function called in finalize(). Creates all maps far all nonlinear solver managers
         */
        void create_maps();

//--------------------------------------------------------------------------------------------------------


//--------------------------------------------------------------------------------------------------------

    protected:

    public:

//--------------------------------------------------------------------------------------------------------
        /**
         * @brief Constructor.
         *
         * @param[in] aSolverInterface Pointer to the solver interface
         */
        Nonlinear_Database( Solver_Interface * aSolverInterface ) : mSolverInterface( aSolverInterface ){};

//--------------------------------------------------------------------------------------------------------

        /**
         * @brief Destructor.
         */
        ~Nonlinear_Database(){};

//--------------------------------------------------------------------------------------------------------

        /**
         * @brief Finalize call. Calculates dependencies, maps and ships pointers after all information is received.
         */
        void finalize();

//--------------------------------------------------------------------------------------------------------

        /**
         * @brief Memeber function to set the nonliner solver managers. The highest level nonliner solver manager has to be on entry 0
         */
        void set_nonliner_solver_managers( Nonlinear_Solver * aNonlinerSolverManager );

//--------------------------------------------------------------------------------------------------------
        /**
         * @brief Calls finalize() and initializes the solve for the highest system
         */
        void solve();

//--------------------------------------------------------------------------------------------------------
        /**
          * @brief Returns the index of a certain nonliner solver manager
          *
          * @param[in] aSolverManagerIndex The index of the asking nonlinear solver manager
          * @param[in] aDofTypeListIndex The index of dof type list on the asking nonliner solver manager
          */
        moris::sint get_nonlinear_solver_manager_index( const moris::sint aSolverManagerIndex,
                                                        const moris::sint aDofTypeListIndex );

//--------------------------------------------------------------------------------------------------------
        /**
          * @brief Returns the free map for the asking nonliner solver manager
          *
          * @param[in] aSolverManagerIndex The index of the asking nonlinear solver manager
          */
        Map_Class * get_list_of_maps( const moris::sint aSolverManagerIndex );

//--------------------------------------------------------------------------------------------------------
        /**
          * @brief Returns a pointer to the solver interface.
          */
        Solver_Interface * get_solver_interface(){ return mSolverInterface; };

//--------------------------------------------------------------------------------------------------------
        /**
          * @brief Returns the nonlinear solver manager list.
          */
        moris::Cell< Nonlinear_Solver * > & get_nonliner_solver_manager_list(){ return mListNonlinerSolverManagers; };

//--------------------------------------------------------------------------------------------------------
        /**
          * @brief Returns a pointer to the full vector
          */
        Dist_Vector * get_full_vector(){ return mFullVector; };

    };
}
}
#endif /* MORIS_DISTLINALG_CL_NLA_NONLINEAR_DATABASE_HPP_ */

