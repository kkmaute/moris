#include "cl_Tracer_Enums.hpp"

#include <cstdio>
#include <string>
#include <cstring>

//#include "assert.hpp"


namespace moris
{


// ------------------------------------------------------------------------------------------ //
// define Entities -------------------------------------------------------------------------- //
// ------------------------------------------------------------------------------------------ //

// copy string
//inline
const std::string get_enum_str(enum EntityBase aEntityBase)
{
    switch (aEntityBase)
    {
    case EntityBase::Unknown:           return "Unknown";
    case EntityBase::LinearSolver:      return "LinearSolver";
    case EntityBase::LinearProblem:     return "LinearProblem";
    case EntityBase::NonLinearSolver:   return "NonLinearSolver";
    case EntityBase::NonLinearProblem:  return "NonLinearProblem";
    case EntityBase::TimeSolver:        return "TimeSolver";
    case EntityBase::BlockSet:          return "BlockSet";
    case EntityBase::GlobalClock:       return "GlobalClock";
    case EntityBase::MSI:               return "MSI";
    case EntityBase::Mesh:              return "Mesh";

    default:
//        MORIS_ASSERT(false, "Invalid EntityBase Enum provided.");
        return "0";
    }
}

enum EntityBase get_entity_base_enum_from_str(std::string aEnumString)
{
    if      (aEnumString == "Unknown")          return EntityBase::Unknown;
    else if (aEnumString == "LinearSolver")     return EntityBase::LinearSolver;
    else if (aEnumString == "LinearProblem")    return EntityBase::LinearProblem;
    else if (aEnumString == "NonLinearSolver")  return EntityBase::NonLinearSolver;
    else if (aEnumString == "NonLinearProblem") return EntityBase::NonLinearProblem;
    else if (aEnumString == "TimeSolver")       return EntityBase::TimeSolver;
    else if (aEnumString == "BlockSet")         return EntityBase::BlockSet;
    else if (aEnumString == "GlobalClock")      return EntityBase::GlobalClock;
    else if (aEnumString == "MSI")              return EntityBase::MSI;
    else if (aEnumString == "Mesh")             return EntityBase::Mesh;

    else
    {
//        MORIS_ASSERT(false, "Invalid EntityBase String provided.");
        return EntityBase::Unknown;
    }
}


// ------------------------------------------------------------------------------------------ //
// define Types ----------------------------------------------------------------------------- //
// ------------------------------------------------------------------------------------------ //

//inline
const std::string get_enum_str(enum EntityType aEntityType)
{
    switch (aEntityType)
    {
    case EntityType::Unknown:       return "Unknown";
    case EntityType::Arbitrary:     return "Arbitrary";
    case EntityType::NoType:        return "NoType";
    case EntityType::Base:          return "Base";
    case EntityType::Gauss:         return "Gauss";
    case EntityType::GaussSeidel:   return "GaussSeidel";
    case EntityType::Amesos:        return "Amesos";
    case EntityType::Aztec:         return "Aztec";
    case EntityType::PETSc:         return "PETSc";
    case EntityType::Monolythic:    return "Monolythic";
    case EntityType::Newton:        return "Newton";
    case EntityType::Arclength:     return "Arclength";

    default:
//        MORIS_ASSERT(false, "Invalid EntityType Enum provided.");
        return "0";
    }
}

enum EntityType get_entity_type_enum_from_str(std::string aEnumString)
{
    if      (aEnumString == "Unknown")      return EntityType::Unknown;
    else if (aEnumString == "Arbitrary")    return EntityType::Arbitrary;
    else if (aEnumString == "NoType")       return EntityType::NoType;
    else if (aEnumString == "Base")         return EntityType::Base;
    else if (aEnumString == "Gauss")        return EntityType::Gauss;
    else if (aEnumString == "GaussSeidel")  return EntityType::GaussSeidel;
    else if (aEnumString == "Amesos")       return EntityType::Amesos;
    else if (aEnumString == "Aztec")        return EntityType::Aztec;
    else if (aEnumString == "PETSc")        return EntityType::PETSc;
    else if (aEnumString == "Monolythic")   return EntityType::Monolythic;
    else if (aEnumString == "Newton")       return EntityType::Newton;
    else if (aEnumString == "Arclength")    return EntityType::Arclength;

    else
    {
//        MORIS_ASSERT(false, "Invalid EntityType String provided.");
        return EntityType::Unknown;
    }
}

// ------------------------------------------------------------------------------------------ //
// define Actions --------------------------------------------------------------------------- //
// ------------------------------------------------------------------------------------------ //

//inline
const std::string get_enum_str(enum EntityAction aEntityAction)
{
    switch (aEntityAction)
    {
    case EntityAction::Unknown:             return "Unknown";
    case EntityAction::Arbitrary:           return "Arbitrary";
    case EntityAction::Solve:               return "Solve";
    case EntityAction::Build:               return "Build";
    case EntityAction::Assemble:            return "Assemble";
    case EntityAction::Compute:             return "Compute";
    case EntityAction::Create:              return "Create";
    case EntityAction::Evaluate:            return "Evaluate";
    case EntityAction::AssembleJacobian:    return "AssembleJacobian";
    case EntityAction::AssembleResidual:    return "AssembleResidual";
    case EntityAction::AssembleJacAndRes:   return "AssembleJacAndRes";
    case EntityAction::AssembleRHS:         return "AssembleRHS";

    default:
//        MORIS_ASSERT(false, "Invalid EntityAction Enum provided.");
        return "0";
    }
}

enum EntityAction get_entity_action_enum_from_str(std::string aEnumString)
{
    if      (aEnumString == "Unknown")            return EntityAction::Unknown;
    else if (aEnumString == "Arbitrary")          return EntityAction::Arbitrary;
    else if (aEnumString == "Solve")              return EntityAction::Solve;
    else if (aEnumString == "Build")              return EntityAction::Build;
    else if (aEnumString == "Assemble")           return EntityAction::Assemble;
    else if (aEnumString == "Compute")            return EntityAction::Compute;
    else if (aEnumString == "Create")             return EntityAction::Create;
    else if (aEnumString == "Evaluate")           return EntityAction::Evaluate;
    else if (aEnumString == "AssembleJacobian")   return EntityAction::AssembleJacobian;
    else if (aEnumString == "AssembleResidual")   return EntityAction::AssembleResidual;
    else if (aEnumString == "AssembleJacAndRes")  return EntityAction::AssembleJacAndRes;
    else if (aEnumString == "AssembleRHS")        return EntityAction::AssembleRHS;

    else
    {
//        MORIS_ASSERT(false, "Invalid EntityAction String provided.");
        return EntityAction::Unknown;
    }
}


// ------------------------------------------------------------------------------------------ //
// define Outputs --------------------------------------------------------------------------- //
// ------------------------------------------------------------------------------------------ //

//inline
const std::string get_enum_str(enum OutputSpecifier aOutputSpecifier)
{
    switch (aOutputSpecifier)
    {
    case OutputSpecifier::Unknown:      return "Unknown";
    case OutputSpecifier::ResidualNorm: return "ResidualNorm";
    case OutputSpecifier::ResidualDrop: return "ResidualDrop";
    case OutputSpecifier::SolutionNorm: return "SolutionNorm";
    case OutputSpecifier::Iteration:    return "Iteration";
    case OutputSpecifier::Count:        return "Count";
    case OutputSpecifier::SignIn:       return "SignIn";
    case OutputSpecifier::ElapsedTime:  return "ElapsedTime";
    case OutputSpecifier::Error:        return "Error";
    case OutputSpecifier::Step:         return "Step";
    case OutputSpecifier::Restart:      return "Restart";
    case OutputSpecifier::FreeText:     return "FreeText";

    default:
//        MORIS_ASSERT(false, "Invalid OutputSpecifier Enum provided.");
        return "0";
    }
}


enum OutputSpecifier get_output_spec_enum_from_str(std::string aEnumString)
{
    if      (aEnumString == "Unknown")      return OutputSpecifier::Unknown;
    else if (aEnumString == "ResidualNorm") return OutputSpecifier::ResidualNorm;
    else if (aEnumString == "ResidualDrop") return OutputSpecifier::ResidualDrop;
    else if (aEnumString == "SolutionNorm") return OutputSpecifier::SolutionNorm;
    else if (aEnumString == "Iteration")    return OutputSpecifier::Iteration;
    else if (aEnumString == "Count")        return OutputSpecifier::Count;
    else if (aEnumString == "SignIn")       return OutputSpecifier::SignIn;
    else if (aEnumString == "ElapsedTime")  return OutputSpecifier::ElapsedTime;
    else if (aEnumString == "Error")        return OutputSpecifier::Error;
    else if (aEnumString == "Step" )        return OutputSpecifier::Step;
    else if (aEnumString == "Restart")      return OutputSpecifier::Restart;
    else if (aEnumString == "FreeText")     return OutputSpecifier::FreeText;
    else if (aEnumString == "InfoText")     return OutputSpecifier::InfoText;
    else if (aEnumString == "DebugText")    return OutputSpecifier::DebugText;
    else if (aEnumString == "Warning")      return OutputSpecifier::Warning;

    else
    {
//        MORIS_ASSERT(false, "Invalid OutputSpecifier String provided.");
        return OutputSpecifier::Unknown;
    }
}

} // namespace moris

