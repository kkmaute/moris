#ifndef PROJECTS_MRS_IOS_SRC_CL_TRACER_ENUMS_HPP_
#define PROJECTS_MRS_IOS_SRC_CL_TRACER_ENUMS_HPP_

#include <iostream>
#include <fstream>
#include <cstdio>
#include <string>
#include <cstring>


#include "assert.hpp"


namespace moris
{

// define Entities
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
    Mesh
};

// copy string
inline
const std::string get_enum_str(enum EntityBase aEntityBase)
{
    switch (aEntityBase)
    {
    case EntityBase::Unknown: return "Unknown";
    case EntityBase::LinearSolver: return "LinearSolver";
    case EntityBase::LinearProblem: return "LinearProblem";
    case EntityBase::NonLinearSolver: return "NonLinearSolver";
    case EntityBase::NonLinearProblem: return "NonLinearProblem";
    case EntityBase::TimeSolver: return "TimeSolver";
    case EntityBase::BlockSet: return "BlockSet";
    case EntityBase::GlobalClock: return "GlobalClock";
    case EntityBase::MSI: return "MSI";
    case EntityBase::Mesh: return "Mesh";

    default: return "invalid EntityBase enum provided";
    }
}

// define Types
enum class EntityType {
    Unknown,
    NoType,
    Base,
    Gauss,
    GaussSeidel,
    Amesos,
    Aztec,
    PETSc,
    Monolythic,
    Newton,
    Arclength
};

// copy string
inline
const std::string get_enum_str(enum EntityType aEntityType)
{
    switch (aEntityType)
    {
    case EntityType::Unknown: return "Unknown";
    case EntityType::NoType: return "NoType";
    case EntityType::Base: return "Base";
    case EntityType::Gauss: return "Gauss";
    case EntityType::GaussSeidel: return "GaussSeidel";
    case EntityType::Amesos: return "Amesos";
    case EntityType::Aztec: return "Aztec";
    case EntityType::PETSc: return "PETSc";
    case EntityType::Monolythic: return "Monolythic";
    case EntityType::Newton: return "Newton";
    case EntityType::Arclength: return "Arclength";

    default: return "invalid EntityType enum provided";
    }
}

// define Actions
enum class EntityAction {
    Unknown,
    Solve,
    Build,
    Assemble,
    Compute,
    Create,
    Evaluate,
    AssembleJacobian,
    AssembleResidual,
    AssembleJacAndRes,
    AssembleRHS
};

// copy string
inline
const std::string get_enum_str(enum EntityAction aEntityAction)
{
    switch (aEntityAction)
    {
    case EntityAction::Unknown: return "Unknown";
    case EntityAction::Solve: return "Solve";
    case EntityAction::Build: return "Build";
    case EntityAction::Assemble: return "Assemble";
    case EntityAction::Compute: return "Compute";
    case EntityAction::Create: return "Create";
    case EntityAction::Evaluate: return "Evaluate";
    case EntityAction::AssembleJacobian: return "AssembleJacobian";
    case EntityAction::AssembleResidual: return "AssembleResidual";
    case EntityAction::AssembleJacAndRes: return "AssembleJacAndRes";
    case EntityAction::AssembleRHS: return "AssembleRHS";


    default: return "invalid EntityAction enum provided";
    }
}

// define Outputs
enum class OutputSpecifier {
    Unknown,
    ResidualNorm,
    SolutionNorm,
    Iteration,
    Count,
    Signing,
    Time,
    Error,
    Steps,
    Restart
};

// copy string
inline
const std::string get_enum_str(enum OutputSpecifier aOutputSpecifier)
{
    switch (aOutputSpecifier)
    {
    case OutputSpecifier::Unknown: return "Unknown";
    case OutputSpecifier::ResidualNorm: return "ResidualNorm";
    case OutputSpecifier::SolutionNorm: return "SolutionNorm";
    case OutputSpecifier::Iteration: return "Iteration";
    case OutputSpecifier::Count: return "Count";
    case OutputSpecifier::Signing: return "Signing";
    case OutputSpecifier::Time: return "Time";
    case OutputSpecifier::Error: return "Error";
    case OutputSpecifier::Steps: return "Steps";
    case OutputSpecifier::Restart: return "Restart";


    default: return "invalid OutputSpecifier enum provided";
    }
}


} // namespace moris
#endif /* PROJECTS_MRS_IOS_SRC_CL_TRACER_ENUMS_HPP_ */
