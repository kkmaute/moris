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

#include "cl_Bitset.hpp"
#include "cl_Map.hpp"
#include "fn_enum_macros.hpp"
#include "assert.hpp"
#include "cl_Vector.hpp"

namespace moris::sol
{
    ENUM_MACRO( SolverType,
            AZTEC_IMPL,     //< Wrapper around Aztec Solver
            AMESOS_IMPL,    //< Wrapper around Amesos Solver
            BELOS_IMPL,     //< Wrapper around Belos Solver
            PETSC,          //< Wrapper around Petsc Solver
            EIGEN_SOLVER,
            SLEPC_SOLVER,
            ML,    //< Wrapper around ML Preconditioner as a solver
            END_ENUM )

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

    ENUM_MACRO( MapType,
            Epetra,
            Petsc )

    /**
     * SolverRelaxationType notes
     * Constant: Constant relaxation parameter
     * InvResNorm: Relaxation parameter proportional to inverse of residual norm
     * InvResNormAdaptive: Relaxation parameter proportional to inverse of residual norm with adaptation
     */
    ENUM_MACRO( SolverRelaxationType,
            Constant,
            InvResNorm,
            InvResNormAdaptive )

    /**
     * SolverLoadControlType notes
     * Constant: Constant load control parameter
     * Exponential: Exponential growth
     * UserDefined: User defined strategy
     */
    ENUM_MACRO( SolverLoadControlType,
            Constant,
            Linear,
            Exponential,
            UserDefined )

    ENUM_MACRO( SolverRaytracingStrategy,
        None,                                 // No (re-)ray tracing
        EveryNthIteration,                    // Ray tracing after every Nth newton iteration
        EveryNthLoadStep,                     // Ray tracing after every Nth load step
        EveryNthLoadStepOrNthIteration,       // Ray tracing after every Nth load step or Nth iteration
        ResidualChange,                       // Ray tracing if the change of the residual is below a certain threshold
        MixedNthLoadStepAndResidualChange,    // Uses the nth load steps until full load is reached and then uses the relative residual change
        MixedNthLoadStepAndNthIteration,      // Uses the nth load steps until full load is reached and then uses the nth iteration
    )

    ENUM_MACRO( SolverPseudoTimeControlType,
            None,                  // No pseudo time step control
            Polynomial,            // Time step index based strategy: polynomial growth
            InvResNorm,            // Residual based strategy
            Hybrid,                // Combined Polynomial and InvResNorm stratgies
            Exponential,           // Time step index based strategy: exponential growth
            SwitchedRelaxation,    // Switched relaxation (based on Ceze and Fidkowski, 2013)
            ResidualDifference,    // Monotonic residual difference method (based on Ceze and Fidkowski, 2013)
            Expur,                 // Exponential with under-relaxation (based on Ceze and Fidkowski, 2013)
            Comsol                 // COMSOL ( see COMSOL_CFDModuleUsersGuide 6.0, page 92, 241)
    )

    // enum for the type of preconditioner
    ENUM_MACRO( PreconditionerType,
            NONE,
            IFPACK,    // Ifpack
            ML,        // ML
            PETSC,     // Petsc
            END_ENUM )

    enum class EiegnSolverType
    {
        NONE,
        ANASAZI,    // Anasazi
        SLEPC,      // SLEPc
        END_ENUM
    };

    enum class STTypeSlepc
    {
        NONE,
        STSHELL,
        STSHIFT,
        STSINVERT,
        STCAYLEY,
        STPRECOND,
        STFILTER,
        END_ENUM
    };

    enum class SensitivityAnalysisType
    {
        ADJOINT,
        DIRECT,
        END_ENUM
    };
}    // namespace moris::sol

#endif /* SRC_DISTLINALG_CL_DLA_ENUMS_HPP_ */
