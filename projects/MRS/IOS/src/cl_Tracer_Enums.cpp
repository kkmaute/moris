/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_Tracer_Enums.cpp
 *
 */

#include "cl_Tracer_Enums.hpp"

#include <cstdio>
#include <string>
#include <cstring>

namespace moris
{
    // ------------------------------------------------------------------------------------------ //
    // define Entities -------------------------------------------------------------------------- //
    // ------------------------------------------------------------------------------------------ //

    // copy string
    //inline
    std::string get_enum_str( enum EntityBase aEntityBase )
    {
        switch (aEntityBase)
        {
            case EntityBase::Unknown:                 return "Unknown";
            case EntityBase::LinearSolver:            return "LinearSolver";
            case EntityBase::LinearProblem:           return "LinearProblem";
            case EntityBase::NonLinearSolver:         return "NonLinearSolver";
            case EntityBase::NonLinearProblem:        return "NonLinearProblem";
            case EntityBase::TimeSolver:              return "TimeSolver";
            case EntityBase::GlobalClock:             return "GlobalClock";
            case EntityBase::MSI:                     return "MSI";
            case EntityBase::Mesh:                    return "Mesh";
            case EntityBase::XTK:                     return "XTK";
            case EntityBase::OptimizationAlgorithm:   return "OptimizationAlgorithm";
            case EntityBase::TEST_CASE:               return "TEST_CASE";

            default:
                //        MORIS_ASSERT(false, "Invalid EntityBase Enum provided.");
                return "0";
        }
    }

    enum EntityBase get_entity_base_enum_from_str( const std::string& aEnumString )
    {
        if      (aEnumString == "Unknown")                   return EntityBase::Unknown;
        else if (aEnumString == "LinearSolver")              return EntityBase::LinearSolver;
        else if (aEnumString == "LinearProblem")             return EntityBase::LinearProblem;
        else if (aEnumString == "NonLinearSolver")           return EntityBase::NonLinearSolver;
        else if (aEnumString == "NonLinearProblem")          return EntityBase::NonLinearProblem;
        else if (aEnumString == "TimeSolver")                return EntityBase::TimeSolver;
        else if (aEnumString == "GlobalClock")               return EntityBase::GlobalClock;
        else if (aEnumString == "MSI")                       return EntityBase::MSI;
        else if (aEnumString == "Mesh")                      return EntityBase::Mesh;
        else if (aEnumString == "XTK")                       return EntityBase::XTK;
        else if (aEnumString == "OptimizationAlgorithm")     return EntityBase::OptimizationAlgorithm;
        else if (aEnumString == "TEST_CASE")                 return EntityBase::TEST_CASE;

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
    std::string get_enum_str( enum EntityType aEntityType )
    {
        switch (aEntityType)
        {
            case EntityType::Unknown:            return "Unknown";
            case EntityType::Arbitrary:          return "Arbitrary";
            case EntityType::NoType:             return "NoType";
            case EntityType::Base:               return "Base";
            case EntityType::Gauss:              return "Gauss";
            case EntityType::GaussSeidel:        return "GaussSeidel";
            case EntityType::Amesos:             return "Amesos";
            case EntityType::Aztec:              return "Aztec";
            case EntityType::PETSc:              return "PETSc";
            case EntityType::Monolythic:         return "Monolythic";
            case EntityType::Staggered:          return "Staggered";
            case EntityType::Newton:             return "Newton";
            case EntityType::Arclength:          return "Arclength";
            case EntityType::GCMMA:              return "GCMMA";
            case EntityType::LBFGS:              return "LBFGS";
            case EntityType::SQP:                return "SQP";
            case EntityType::Sweep:              return "Sweep";
            case EntityType::EquationModel:      return "EquationModel";
            case EntityType::Overall:            return "Overall";
            case EntityType::Decomposition:      return "Decomposition";
            case EntityType::Enrichment:         return "Enrichment";
            case EntityType::Multigrid:          return "Multigrid";
            case EntityType::GhostStabilization: return "GhostStabilization";

            default:
                //        MORIS_ASSERT(false, "Invalid EntityType Enum provided.");
                return "0";
        }
    }

    enum EntityType get_entity_type_enum_from_str( const std::string& aEnumString )
    {
        if      (aEnumString == "Unknown")            return EntityType::Unknown;
        else if (aEnumString == "Arbitrary")          return EntityType::Arbitrary;
        else if (aEnumString == "NoType")             return EntityType::NoType;
        else if (aEnumString == "Base")               return EntityType::Base;
        else if (aEnumString == "Gauss")              return EntityType::Gauss;
        else if (aEnumString == "GaussSeidel")        return EntityType::GaussSeidel;
        else if (aEnumString == "Amesos")             return EntityType::Amesos;
        else if (aEnumString == "Aztec")              return EntityType::Aztec;
        else if (aEnumString == "PETSc")              return EntityType::PETSc;
        else if (aEnumString == "Monolythic")         return EntityType::Monolythic;
        else if (aEnumString == "Staggered")          return EntityType::Staggered;
        else if (aEnumString == "Newton")             return EntityType::Newton;
        else if (aEnumString == "Arclength")          return EntityType::Arclength;
        else if (aEnumString == "GCMMA")              return EntityType::GCMMA;
        else if (aEnumString == "SQP")                return EntityType::SQP;
        else if (aEnumString == "LBFGS")              return EntityType::LBFGS;
        else if (aEnumString == "Sweep")              return EntityType::Sweep;
        else if (aEnumString == "EquationModel")      return EntityType::EquationModel;
        else if (aEnumString == "Overall")            return EntityType::Overall;
        else if (aEnumString == "Decomposition")      return EntityType::Decomposition;
        else if (aEnumString == "Enrichment")         return EntityType::Enrichment;
        else if (aEnumString == "Multigrid")          return EntityType::Multigrid;
        else if (aEnumString == "GhostStabilization") return EntityType::GhostStabilization;

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
    std::string get_enum_str( enum EntityAction aEntityAction )
    {
        switch (aEntityAction)
        {
            case EntityAction::Unknown:                return "Unknown";
            case EntityAction::Arbitrary:              return "Arbitrary";
            case EntityAction::Solve:                  return "Solve";
            case EntityAction::Build:                  return "Build";
            case EntityAction::Assemble:               return "Assemble";
            case EntityAction::Compute:                return "Compute";
            case EntityAction::Compute_dQIdp_Expl:     return "Compute_dQIdp_Expl";
            case EntityAction::Compute_dQIdp_Impl:     return "Compute_dQIdp_Impl";
            case EntityAction::ComputedQIdpExplImpl:   return "ComputedQIdpExplImpl";
            case EntityAction::Create:                 return "Create";
            case EntityAction::Evaluate:               return "Evaluate";
            case EntityAction::AssembleJacobian:       return "AssembleJacobian";
            case EntityAction::AssembleResidual:       return "AssembleResidual";
            case EntityAction::AssembleJacAndRes:      return "AssembleJacAndRes";
            case EntityAction::AssembleRHS:            return "AssembleRHS";
            case EntityAction::Enrich:                 return "Enrich";
            case EntityAction::Decompose:              return "Decompose";
            case EntityAction::DecomposeRegularHex8:   return "DecomposeRegularHex8";
            case EntityAction::DecomposeRegularQuad4:  return "DecomposeRegularQuad4";
            case EntityAction::DecomposeHierarchyTet4: return "DecomposeHierarchyTet4";
            case EntityAction::DecomposeHierarchyTri3: return "DecomposeHierarchyTri3";
            case EntityAction::Stabilize:              return "Stabilize";
            case EntityAction::Visualize:              return "Visualize";
            case EntityAction::Run:                    return "Run";

            default:
                //        MORIS_ASSERT(false, "Invalid EntityAction Enum provided.");
                return "0";
        }
    }

    enum EntityAction get_entity_action_enum_from_str( const std::string& aEnumString )
    {
        if      (aEnumString == "Unknown")                 return EntityAction::Unknown;
        else if (aEnumString == "Arbitrary")               return EntityAction::Arbitrary;
        else if (aEnumString == "Solve")                   return EntityAction::Solve;
        else if (aEnumString == "Build")                   return EntityAction::Build;
        else if (aEnumString == "Assemble")                return EntityAction::Assemble;
        else if (aEnumString == "Compute")                 return EntityAction::Compute;
        else if (aEnumString == "Compute_dQIdp_Expl")      return EntityAction::Compute_dQIdp_Expl;
        else if (aEnumString == "Compute_dQIdp_Impl")      return EntityAction::Compute_dQIdp_Impl;
        else if (aEnumString == "ComputedQIdpExplImpl")    return EntityAction::ComputedQIdpExplImpl;
        else if (aEnumString == "Create")                  return EntityAction::Create;
        else if (aEnumString == "Evaluate")                return EntityAction::Evaluate;
        else if (aEnumString == "AssembleJacobian")        return EntityAction::AssembleJacobian;
        else if (aEnumString == "AssembleResidual")        return EntityAction::AssembleResidual;
        else if (aEnumString == "AssembleJacAndRes")       return EntityAction::AssembleJacAndRes;
        else if (aEnumString == "AssembleRHS")             return EntityAction::AssembleRHS;
        else if (aEnumString == "Enrich")                  return EntityAction::Enrich;
        else if (aEnumString == "Decompose")               return EntityAction::Decompose;
        else if (aEnumString == "DecomposeRegularHex8")    return EntityAction::DecomposeRegularHex8;
        else if (aEnumString == "DecomposeRegularQuad4")   return EntityAction::DecomposeRegularQuad4;
        else if (aEnumString == "DecomposeHierarchyTet4")  return EntityAction::DecomposeHierarchyTet4;
        else if (aEnumString == "Stabilize")               return EntityAction::Stabilize;
        else if (aEnumString == "Visualize")               return EntityAction::Visualize;
        else if (aEnumString == "Run")                     return EntityAction::Run;

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
    std::string get_enum_str( enum OutputSpecifier aOutputSpecifier )
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
            case OutputSpecifier::InfoText:     return "InfoText";
            case OutputSpecifier::DebugText:    return "DebugText";
            case OutputSpecifier::Warning:      return "Warning";
            case OutputSpecifier::ADVvector:    return "ADVvector";

            default:
                //        MORIS_ASSERT(false, "Invalid OutputSpecifier Enum provided.");
                return "0";
        }
    }

    enum OutputSpecifier get_output_spec_enum_from_str( const std::string& aEnumString )
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
        else if (aEnumString == "ADVvector")    return OutputSpecifier::ADVvector;

        else
        {
            //        MORIS_ASSERT(false, "Invalid OutputSpecifier String provided.");
            return OutputSpecifier::Unknown;
        }
    }

} // namespace moris

