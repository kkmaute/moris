/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_Equation_Model.cpp
 *
 */

#include "cl_MSI_Design_Variable_Interface.hpp"
#include "cl_MSI_Solver_Interface.hpp"

#include "cl_MSI_Equation_Object.hpp"
#include "cl_MSI_Equation_Set.hpp"
#include "cl_MSI_Equation_Model.hpp"

#include "cl_SOL_Dist_Vector.hpp"
#include "cl_SOL_Dist_Map.hpp"
#include "cl_SOL_Matrix_Vector_Factory.hpp"
#include "cl_SOL_Enums.hpp"

// Logging package
#include "cl_Logger.hpp"
#include "cl_Tracer.hpp"

#include "cl_Communication_Tools.hpp"

namespace moris::MSI
{
    //------------------------------------------------------------------------------

    Equation_Model::~Equation_Model()
    {
        delete mImplicitdQidp;
        delete mExplicitdQidp;
        delete mdQIdpMap;
    }

    //------------------------------------------------------------------------------

    moris::sint
    Equation_Model::get_num_rhs()
    {
        if ( !mIsForwardAnalysis )
        {
            if ( mIsAdjointSensitivityAnalysis )
            {
                mNumSensitivityAnalysisRHS = this->get_requested_IQI_names().size();
            }
            else
            {
                mNumSensitivityAnalysisRHS = mDesignVariableInterface->get_my_local_global_map().n_rows();
            }

            return mNumSensitivityAnalysisRHS;
        }
        else
        {
            return 1;
        }
    }

    //------------------------------------------------------------------------------

    void
    Equation_Model::compute_IQIs()
    {
        // get number of IQI on the model
        uint tNumIQIsOnModel =
                this->get_requested_IQI_names().size();

        // get local number of equation sets
        uint tNumSets = mFemSets.size();

        // loop over local equation sets
        for ( uint tSetIndex = 0; tSetIndex < tNumSets; tSetIndex++ )
        {
            // initialize treated equation set // FIXME????
            mFemSets( tSetIndex )->initialize_set();

            // get number of IQIs on treated equation set
            uint tNumIQIsOnSet = mFemSets( tSetIndex )->get_number_of_requested_IQIs();

            // if some IQI are requested on treated equation set
            if ( tNumIQIsOnSet > 0 )
            {
                // get number of equation objects on treated equation set
                uint tNumEquationObjectOnSet = mFemSets( tSetIndex )->get_num_equation_objects();

                // loop over equation objects on treated equation set
                for ( uint tEquationObjectIndex = 0; tEquationObjectIndex < tNumEquationObjectOnSet; tEquationObjectIndex++ )
                {
                    // compute QI
                    // FIXME this is elemental right now??
                    mFemSets( tSetIndex )->get_equation_object_list()( tEquationObjectIndex )->compute_QI();

                    // loop over IQIs on model
                    for ( uint tIQIIndex = 0; tIQIIndex < tNumIQIsOnModel; tIQIIndex++ )
                    {
                        // if the IQI vector is empty (due to basis extension), continue
                        if ( mFemSets( tSetIndex )->get_QI().size() == 0 )
                        {
                            continue;
                        }
                        // assemble QI values into global vector
                        mGlobalIQIVal( tIQIIndex ) += mFemSets( tSetIndex )->get_QI()( tIQIIndex );
                    }
                }
                // free memory on treated equation set
                mFemSets( tSetIndex )->free_matrix_memory();
            }
        }

        // Normalization
        if ( gLogger.mIteration == 0 )
        {
            this->normalize_IQIs();
        }
    }

    //------------------------------------------------------------------------------

    void
    Equation_Model::initialize_explicit_and_implicit_dQIdp()
    {
        // create map for dQIdpMap
        moris::sol::Matrix_Vector_Factory tMatFactory( sol::MapType::Epetra );

        delete mdQIdpMap;

        mdQIdpMap = tMatFactory.create_map(
                mDesignVariableInterface->get_my_local_global_map() );

        // get number of RHS for either adjoint or direct method
        uint tNumRHMS = this->get_num_rhs();

        // create vector for dQIdp implicit and explicit contributions
        delete mImplicitdQidp;
        delete mExplicitdQidp;

        mImplicitdQidp = tMatFactory.create_vector( mdQIdpMap, tNumRHMS );
        mExplicitdQidp = tMatFactory.create_vector( mdQIdpMap, tNumRHMS );

        // fill dQIdp implicit/explicit vectors with zeros
        mExplicitdQidp->vec_put_scalar( 0.0 );
        mImplicitdQidp->vec_put_scalar( 0.0 );
    }

    //------------------------------------------------------------------------------

    void
    Equation_Model::compute_explicit_and_implicit_dQIdp()
    {
        // Trace this function
        Tracer tTracer( "MSI", "Equation Model", "Compute dQI/dp explicit and implicit" );

        // get local number of equation sets
        uint tNumSets = mFemSets.size();

        // loop over local equation sets
        for ( uint tSetIndex = 0; tSetIndex < tNumSets; tSetIndex++ )
        {
            // get access to the equation objects
            Vector< MSI::Equation_Object* >& tEqnObjList = mFemSets( tSetIndex )->get_equation_object_list();

            // get number of equation object on treated equation set
            uint tNumEqObjOnSet = tEqnObjList.size();

            // initialize treated equation set
            mFemSets( tSetIndex )->initialize_set();

            // log which set is currently being treated
            if ( tNumEqObjOnSet > 0 )
            {
                MORIS_LOG_SPEC( "Process FEM set", mFemSets( tSetIndex )->get_set_name() );
            }

            // loop over equation objects on treated equation set
            for ( uint tEqObjIndex = 0; tEqObjIndex < tNumEqObjOnSet; tEqObjIndex++ )
            {
                // compute dQIdp implicit
                tEqnObjList( tEqObjIndex )->compute_dQIdp_explicit_implicit();
            }

            // free memory on treated equation set
            mFemSets( tSetIndex )->free_matrix_memory();
        }

        // global assembly to switch entries to the right processor
        mExplicitdQidp->vector_global_assembly();

        // global assembly to switch entries to the right processor
        mImplicitdQidp->vector_global_assembly();
    }

    //------------------------------------------------------------------------------

    //    void
    //    Equation_Model::compute_implicit_dQIdp()
    //    {
    //        // Trace this function
    //        Tracer tTracer( "MSI", "EquationModel", "Compute_dQIdp_Impl" );
    //
    //        // get local number of equation sets
    //        uint tNumSets = mFemSets.size();
    //
    //        // loop over local equation sets
    //        for ( uint tSetIndex = 0; tSetIndex < tNumSets; tSetIndex++ )
    //        {
    //            // get number of equation object on treated equation set
    //            uint tNumEquationObjectOnSet =
    //                    mFemSets( tSetIndex )->get_num_equation_objects();
    //
    //            // initialize treated equation set //FIXME????
    //            mFemSets( tSetIndex )->initialize_set();
    //
    //            // loop over equation objects on treated equation set
    //            for ( uint tEquationObjectIndex = 0; tEquationObjectIndex < tNumEquationObjectOnSet; tEquationObjectIndex++ )
    //            {
    //                // compute dQIdp implicit
    //                mFemSets( tSetIndex )->get_equation_object_list()( tEquationObjectIndex )->compute_dQIdp_implicit();
    //            }
    //            // free memory on treated equation set
    //            mFemSets( tSetIndex )->free_matrix_memory();
    //        }
    //
    //        // global assembly to switch entries to the right processor
    //        mImplicitdQidp->vector_global_assembly();
    //    }
    //
    //    //------------------------------------------------------------------------------
    //
    //    void
    //    Equation_Model::compute_explicit_dQIdp()
    //    {
    //        // Trace this function
    //        Tracer tTracer( "MSI", "EquationModel", "Compute_dQIdp_Expl" );
    //
    //        // get local number of equation sets
    //        uint tNumSets = mFemSets.size();
    //
    //        // loop over local equation sets
    //        for ( uint iSet = 0; iSet < tNumSets; iSet++ )
    //        {
    //            // get number of equation objects on treated equation set
    //            uint tNumEquationObjectOnSet =
    //                    mFemSets( iSet )->get_num_equation_objects();
    //
    //            // initialize treated equation set //FIXME????
    //            mFemSets( iSet )->initialize_set();
    //
    //            // if some IQI are requested on treated equation set
    //            if ( mFemSets( iSet )->get_number_of_requested_IQIs() > 0 )
    //            {
    //                // loop over equation objects on treated equation set
    //                for ( uint iEqObj = 0; iEqObj < tNumEquationObjectOnSet; iEqObj++ )
    //                {
    //                    // compute dQIdp explicit
    //                    mFemSets( iSet )->get_equation_object_list()( iEqObj )->compute_dQIdp_explicit();
    //                }
    //                // free memory on treated equation set
    //                mFemSets( iSet )->free_matrix_memory();
    //            }
    //        }
    //
    //        // global assembly to switch entries to the right processor
    //        mExplicitdQidp->vector_global_assembly();
    //    }

    //-------------------------------------------------------------------------------------------------

    sol::Dist_Vector*
    Equation_Model::get_dQIdp()
    {
        // create map object
        moris::sol::Matrix_Vector_Factory tMatFactory( sol::MapType::Epetra );

        // get number of RHD
        uint tNumRHMS = this->get_num_rhs();

        // create vector for dQIdp
        delete mdQIdp;
        mdQIdp = tMatFactory.create_vector( mdQIdpMap, tNumRHMS );

        // fill vector with zero
        mdQIdp->vec_put_scalar( 0.0 );

        // add explicit contribution to dQIdp
        mdQIdp->vec_plus_vec( 1.0, *mExplicitdQidp, 1.0 );

        // add implicit contribution to dQIdp
        mdQIdp->vec_plus_vec( 1.0, *mImplicitdQidp, 1.0 );

        // return dQIdp
        return mdQIdp;
    }

}    // namespace moris::MSI
