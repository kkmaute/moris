/*
 * cl_Equation_Model.cpp
 *
 *  Created on: May 08, 2020
 *      Author: schmidt
 */

#include "cl_MSI_Design_Variable_Interface.hpp"
#include "cl_MSI_Solver_Interface.hpp"

#include "cl_MSI_Equation_Object.hpp"
#include "cl_MSI_Equation_Set.hpp"
#include "cl_MSI_Equation_Model.hpp"

#include "cl_SOL_Dist_Vector.hpp"
#include "cl_Matrix_Vector_Factory.hpp"
#include "cl_SOL_Enums.hpp"

namespace moris
{
    namespace MSI
    {
        //------------------------------------------------------------------------------

        moris::sint Equation_Model::get_num_rhs( )
        {
            if( !mIsForwardAnalysis )
            {
                mNumSensitivityAnalysisRHS = this->get_requested_IQI_names().size();

                MORIS_ASSERT( mNumSensitivityAnalysisRHS > 0,
                        "MSI::Equation_Model::get_num_rhs(), num rhs not set for sensitivity analysis");

                return mNumSensitivityAnalysisRHS;
            }
            else
            {
                return 1;
            }
        }

        //------------------------------------------------------------------------------

        moris::Cell< moris::Matrix< DDRMat > > Equation_Model::compute_IQIs()
        {
            // Get local number of elements
            moris::uint tNumSets = mFemSets.size();

            moris::uint tNumIQIs = this->get_requested_IQI_names().size();

            moris::Cell< moris::Matrix< DDRMat > > tGloablIQIVal( tNumIQIs );

            for( auto & tQI : tGloablIQIVal )
            {
                // set size for the QI value
                // FIXME assumed scalar
                tQI.set_size( 1, 1, 0.0 );
            }

            // Loop over all local elements to build matrix graph
            for ( moris::uint Ii=0; Ii < tNumSets; Ii++ )
            {
                if( mFemSets( Ii )->get_element_type() == fem::Element_Type::BULK )
                {
                    moris::uint tNumEquationObjectOnSet = mFemSets( Ii )->get_num_equation_objects();

                    mFemSets( Ii )->initialize_set( true );   //FIXME

                    for ( moris::uint Ik=0; Ik < tNumEquationObjectOnSet; Ik++ )
                    {
                        // FIXME this is elemental right now
                        mFemSets( Ii )->get_equation_object_list()( Ik )->compute_QI();

                        for( moris::uint Ij = 0; Ij < tNumIQIs; Ij++ )
                        {
                            tGloablIQIVal( Ij )( 0 ) += mFemSets( Ii )->get_QI()( Ij )( 0 );
                        }
                    }
                    //this->free_block_memory( Ii );
                }
            }
            return tGloablIQIVal;
        }

        //------------------------------------------------------------------------------

        void Equation_Model::compute_implicit_dQIdp()
        {
            // create map object
            moris::Matrix_Vector_Factory tMatFactory( sol::MapType::Epetra );
            mdQiduMap = tMatFactory.create_map( mDesignVariableInterface->get_my_local_global_map() );

            uint tNumRHMS = this->get_num_rhs();

            // full vector and prev full vector
            mImplicitdQidp = tMatFactory.create_vector( mdQiduMap, tNumRHMS );

            mImplicitdQidp->vec_put_scalar( 0.0 );

            // Get local number of elements
            moris::uint tNumSets = mFemSets.size();

            // Loop over all local elements to build matrix graph
            for ( moris::uint Ii=0; Ii < tNumSets; Ii++ )
            {
                moris::uint tNumEquationObjectOnSet = mFemSets( Ii )->get_num_equation_objects();

                mFemSets( Ii )->initialize_set( true );   //FIXME

                for ( moris::uint Ik=0; Ik < tNumEquationObjectOnSet; Ik++ )
                {
                    mFemSets( Ii )->get_equation_object_list()( Ik )->compute_dQIdp_implicit();
                }

                //this->free_block_memory( Ii );
            }

            // global assembly to switch entries to the right processor
            mImplicitdQidp->vector_global_asembly();

            //mImplicitdQidp->print();
        }

        //------------------------------------------------------------------------------

        void Equation_Model::compute_explicit_dQIdp()
        {
            // create map object
            moris::Matrix_Vector_Factory tMatFactory( sol::MapType::Epetra );

            // FIXME create map only once. either implicit or explicit
            mdQiduMap = tMatFactory.create_map( mDesignVariableInterface->get_my_local_global_map() );

            // get the number of rhs
            uint tNumRHMS = this->get_num_rhs();

            // create full vector for dQIdp explicit
            mExplicitdQidp = tMatFactory.create_vector( mdQiduMap, tNumRHMS );
            mExplicitdQidp->vec_put_scalar( 0.0 );

            // get local number of equation sets
            moris::uint tNumSets = mFemSets.size();

            // loop over the local equation sets
            for ( moris::uint iSet = 0; iSet < tNumSets; iSet++ )
            {
                // get number of equation object on the treated equation set
                moris::uint tNumEquationObjectOnSet =
                        mFemSets( iSet )->get_num_equation_objects();

                // initialize the treated equation set
                mFemSets( iSet )->initialize_set( true );

                // loop over the equation objects
                for ( moris::uint iEqObj = 0; iEqObj < tNumEquationObjectOnSet; iEqObj++ )
                {
                    // if some IQI are requested on the treated equation set
                    if( mFemSets( iSet )->get_number_of_requested_IQIs() > 0 )
                    {
                        // compute dQIdp explicit
                        mFemSets( iSet )->get_equation_object_list()( iEqObj )->compute_dQIdp_explicit();
                    }
                }
                //this->free_block_memory( iSet );
            }

            // global assembly to switch entries to the right processor
            mExplicitdQidp->vector_global_asembly();
        }

        //-------------------------------------------------------------------------------------------------

        sol::Dist_Vector * Equation_Model::get_dQidu()
        {
            // create map object
            moris::Matrix_Vector_Factory tMatFactory( sol::MapType::Epetra );

            // get
            uint tNumRHMS = this->get_num_rhs();

            // full vector and prev full vector
            mQidu = tMatFactory.create_vector( mdQiduMap, tNumRHMS );

            mQidu->vec_put_scalar( 0.0 );

            //            mExplicitdQidp->print();
            //            mImplicitdQidp->print();

            mQidu->vec_plus_vec( 1.0, *mExplicitdQidp, 1.0 );
            mQidu->vec_plus_vec( 1.0, *mImplicitdQidp, 1.0 );

            //mQidu->print();

            return mQidu;
        }

        //-------------------------------------------------------------------------------------------------

    }/* end_namespace_msi */
}/* end_namespace_moris */
