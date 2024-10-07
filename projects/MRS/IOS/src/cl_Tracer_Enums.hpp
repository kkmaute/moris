/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_Tracer_Enums.hpp
 *
 */

#ifndef PROJECTS_MRS_IOS_SRC_CL_TRACER_ENUMS_HPP_
#define PROJECTS_MRS_IOS_SRC_CL_TRACER_ENUMS_HPP_

#include <cstdio>
#include <string>
#include <cstring>

namespace moris
{

    // define Entities ---------------------------------------------------------------------------
    enum class EntityBase {
            LinearSolver,
            LinearProblem,
            NonLinearSolver,
            NonLinearProblem,
            TimeSolver,
            GlobalClock,
            MSI,
            Mesh,
            XTK,
            TEST_CASE,
            OptimizationAlgorithm,
            Unknown
    };

    // copy string
    //inline
    std::string get_enum_str( enum EntityBase aEntityBase );

    enum EntityBase get_entity_base_enum_from_str( const std::string& aEnumString );

    // define Types ---------------------------------------------------------------------------
    enum class EntityType {
            Arbitrary,
            NoType,
            Base,
            Gauss,
            GaussSeidel,
            Amesos,
            Aztec,
            PETSc,
            Monolythic,
            Staggered,
            Newton,
            Arclength,
            GCMMA,
            LBFGS,
            SQP,
            Sweep,
            EquationModel,
            Overall,
            Decomposition,
            Enrichment,
            GhostStabilization,
            Multigrid,
            Unknown
    };

    // copy string
    //inline
    std::string get_enum_str( enum EntityType aEntityType );

    enum EntityType get_entity_type_enum_from_str( const std::string& aEnumString );

    // define Actions ---------------------------------------------------------------------------
    enum class EntityAction {
            Arbitrary,
            Solve,
            Build,
            Assemble,
            Compute,
            Compute_dQIdp_Expl,
            Compute_dQIdp_Impl,
            ComputedQIdpExplImpl,
            Create,
            Evaluate,
            AssembleJacobian,
            AssembleResidual,
            AssembleJacAndRes,
            AssembleRHS,
            Decompose,
            DecomposeRegularHex8,
            DecomposeRegularQuad4,
            DecomposeHierarchyTet4,
            DecomposeHierarchyTri3,
            Enrich,
            Stabilize,
            Visualize,
            Run,
            Unknown
    };

    // copy string
    //inline
    std::string get_enum_str( enum EntityAction aEntityAction );

    enum EntityAction get_entity_action_enum_from_str( const std::string& aEnumString );

    // define Outputs ---------------------------------------------------------------------------
    enum class OutputSpecifier {
            ResidualNorm,
            ResidualDrop,
            SolutionNorm,
            Iteration,
            Count,
            SignIn,
            ElapsedTime,
            Error,
            Step,
            Restart,
            FreeText,
            InfoText,
            DebugText,
            Warning,
            ADVvector,
            Unknown
    };

    // copy string
    //inline
    std::string get_enum_str( enum OutputSpecifier aOutputSpecifier );

    enum OutputSpecifier get_output_spec_enum_from_str( const std::string& aEnumString );

} // namespace moris
#endif /* PROJECTS_MRS_IOS_SRC_CL_TRACER_ENUMS_HPP_ */

