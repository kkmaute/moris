/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_SOL_Enums.hpp
 *
 */
#ifndef SRC_DISTLINALG_CL_SOL_ENUMS_HPP_
#define SRC_DISTLINALG_CL_SOL_ENUMS_HPP_

namespace moris
{
    namespace sol
    {
        enum class SolverType
        {
            AZTEC_IMPL,      //< Wrapper around Aztec Solver
            AMESOS_IMPL,     //< Wrapper around Amesos Solver
            BELOS_IMPL,      //< Wrapper around Belos Solver
            PETSC,           //< Wrapper around Petsc Solver
            EIGEN_SOLVER,
            ML,   //< Wrapper around ML Preconditioner as a solver
            END_ENUM
        };

        enum class EigSolMethod
        {
            LINSOL_AMESOS_KLU,
            LINSOL_AMESOS_UMFPACK,
            LINSOL_AMESOS_DSCPACK,
            LINSOL_AMESOS_MUMPS,
            LINSOL_AMESOS_LAPACK,
            LINSOL_AMESOS_SCALAPACK,
            LINSOL_AMESOS_PARDISO,
            END_ENUM
        };

        enum class MapType
        {
            Epetra,    // Indicates the Vector/Matrix/Map type
            Petsc      // Indicates the Vector/Matrix/Map type
        };

        enum class SolverRelaxationType
        {
            Constant,             // Constant relaxation parameter
            InvResNorm,           // Relaxation parameter proportional to inverse of residual norm
            InvResNormAdaptive    // Relaxation parameter proportional to inverse of residual norm with adaptation
        };

        enum class SolverLoadControlType
        {
            Constant,       // Constant load control parameter
            Exponential,    // Exponential growth
            UserDefined     // User defined strategy
        };

        enum class SolverPseudoTimeControlType
        {
            None,                  // No pseudo time step control
            Polynomial,            // Time step index based strategy: polynomial growth
            InvResNorm,            // Residual based strategy
            Hybrid,                // Combined Polynomial and InvResNorm stratgies
            Exponential,           // Time step index based strategy: exponential growth
            SwitchedRelaxation,    // Switched relaxation (based on Ceze and Fidkowski, 2013)
            ResidualDifference,    // Monotonic residual difference method (based on Ceze and Fidkowski, 2013)
            Expur,                 // Exponential with under-relaxation (based on Ceze and Fidkowski, 2013)
            Comsol                 // COMSOL ( see COMSOL_CFDModuleUsersGuide 6.0, page 92, 241)
        };

        // enum for the type of preconditioner
        enum class PreconditionerType
        {
            NONE,
            IFPACK,    // Ifpack
            ML,        // ML
            END_ENUM
        };

        enum class EiegnSolverType
        {
            NONE,
            ANASAZI,    // Anasazi
            SLEPC,      // SLEPc
            END_ENUM
        };
    }    // namespace sol
}    // namespace moris

#endif /* SRC_DISTLINALG_CL_DLA_ENUMS_HPP_ */
