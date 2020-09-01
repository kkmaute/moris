#ifndef PROJECTS_MRS_IOS_SRC_CL_TRACER_ENUMS_HPP_
#define PROJECTS_MRS_IOS_SRC_CL_TRACER_ENUMS_HPP_

#include <cstdio>
#include <string>
#include <cstring>


namespace moris
{

// define Entities ---------------------------------------------------------------------------
enum class EntityBase {
    Unknown,
    LinearSolver,
    LinearProblem,
    NonLinearSolver,
    NonLinearProblem,
    TimeSolver,
    BlockSet,
    GlobalClock,
    MSI,
    Mesh,
    TEST_CASE,
    OptimizationAlgorithm
};

// copy string
//inline
const std::string get_enum_str(enum EntityBase aEntityBase);

enum EntityBase get_entity_base_enum_from_str(std::string aEnumString);

// define Types ---------------------------------------------------------------------------
enum class EntityType {
    Unknown,
    Arbitrary,
    NoType,
    Base,
    Gauss,
    GaussSeidel,
    Amesos,
    Aztec,
    PETSc,
    Monolythic,
    Newton,
    Arclength,
    GCMMA,
    MMA,
    PCG
};

// copy string
//inline
const std::string get_enum_str(enum EntityType aEntityType);

enum EntityType get_entity_type_enum_from_str(std::string aEnumString);

// define Actions ---------------------------------------------------------------------------
enum class EntityAction {
    Unknown,
    Arbitrary,
    Solve,
    Build,
    Assemble,
    Compute,
    Create,
    Evaluate,
    AssembleJacobian,
    AssembleResidual,
    AssembleJacAndRes,
    AssembleRHS,
    Run
};

// copy string
//inline
const std::string get_enum_str(enum EntityAction aEntityAction);

enum EntityAction get_entity_action_enum_from_str(std::string aEnumString);

// define Outputs ---------------------------------------------------------------------------
enum class OutputSpecifier {
    Unknown,
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
    ADVvector
};

// copy string
//inline
const std::string get_enum_str(enum OutputSpecifier aOutputSpecifier);


enum OutputSpecifier get_output_spec_enum_from_str(std::string aEnumString);


} // namespace moris
#endif /* PROJECTS_MRS_IOS_SRC_CL_TRACER_ENUMS_HPP_ */
