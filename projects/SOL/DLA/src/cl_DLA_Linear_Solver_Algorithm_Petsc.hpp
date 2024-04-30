/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_DLA_Linear_Solver_Algorithm_Petsc.hpp
 *
 */

#pragma once

// MORIS header files.
#include "cl_Matrix.hpp"
#include "linalg_typedefs.hpp"

#include "cl_SOL_Matrix_Vector_Factory.hpp"
#include "cl_SOL_Enums.hpp"

#include "cl_Parameter_List.hpp"    // CON/src
#include "cl_DLA_Linear_Solver_Algorithm.hpp"
#include "cl_DLA_Preconditioner_PETSc.hpp"

namespace moris::dla
{

    //class Preconditioner_PETSc;
    class Linear_Solver_Algorithm_Petsc : public Linear_Solver_Algorithm
    {
      protected:
        Preconditioner_PETSc* mPreconditioner = nullptr;

      public:

        //-----------------------------------------------------------------------------------

        /**
         * @brief Construct a new Linear_Solver_Algorithm_Petsc object
         * 
         * @param aParameterlist 
         */
        Linear_Solver_Algorithm_Petsc( const moris::Parameter_List& aParameterlist )
                : Linear_Solver_Algorithm( aParameterlist ){};

        //-----------------------------------------------------------------------------------

        /**
         * @brief Destroy the Linear_Solver_Algorithm_Petsc object
         * 
         */

        virtual ~Linear_Solver_Algorithm_Petsc(){};

        //-----------------------------------------------------------------------------------

        // cast the preconditioner to the correct type(petsc and assign it to the object)
        
        virtual void
        set_preconditioner( Preconditioner* aPreconditioner ) override
        {
            mPreconditioner = dynamic_cast< Preconditioner_PETSc* >( aPreconditioner );
        }
    };
}    // namespace moris::dla
