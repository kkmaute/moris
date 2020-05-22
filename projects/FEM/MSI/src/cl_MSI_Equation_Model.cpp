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
//            mSolutionVector->print();
//             mSensitivitySolutionVector->print();

            // create map object
            moris::Matrix_Vector_Factory tMatFactory( sol::MapType::Epetra );
            mdQiduMap = tMatFactory.create_map( mDesignVariableInterface->get_my_local_global_map() );

            uint tNumRHMS = this->get_num_rhs();

            // full vector and prev full vector
            mImplicitdQidu = tMatFactory.create_vector( mdQiduMap, tNumRHMS );

            mImplicitdQidu->vec_put_scalar( 0.0 );

            // Get local number of elements
            moris::uint tNumSets = mFemSets.size();

            // Loop over all local elements to build matrix graph
            for ( moris::uint Ii=0; Ii < tNumSets; Ii++ )
            {
                std::cout<<"Set "<<Ii<<std::endl;

                moris::uint tNumEquationObjectOnSet = mFemSets( Ii )->get_num_equation_objects();

                mFemSets( Ii )->initialize_set( true );   //FIXME

                for ( moris::uint Ik=0; Ik < tNumEquationObjectOnSet; Ik++ )
                {
                    mFemSets( Ii )->get_equation_object_list()( Ik )->compute_dQIdp_implicit();
                }

                //this->free_block_memory( Ii );
            }

            // global assembly to switch entries to the right processor
            mImplicitdQidu->vector_global_asembly();

            //mImplicitdQidu->print();

            sleep( 5 );
        }

//------------------------------------------------------------------------------
        void Equation_Model::compute_explicit_dQIdp()
        {


            // create map object
            moris::Matrix_Vector_Factory tMatFactory( sol::MapType::Epetra );
            // FIXME create map only once. eiteher implicit or explicit
            mdQiduMap = tMatFactory.create_map( mDesignVariableInterface->get_my_local_global_map() );

            uint tNumRHMS = this->get_num_rhs();

            // full vector and prev full vector
            mExplicitdQidu = tMatFactory.create_vector( mdQiduMap, tNumRHMS );

            mExplicitdQidu->vec_put_scalar( 0.0 );

            // Get local number of elements
            moris::uint tNumSets = mFemSets.size();

            // Loop over all local elements to build matrix graph
            for ( moris::uint Ii=0; Ii < tNumSets; Ii++ )
            {
                moris::uint tNumEquationObjectOnSet = mFemSets( Ii )->get_num_equation_objects();

                std::cout<<"Set "<<Ii<<std::endl;

                mFemSets( Ii )->initialize_set( true );

                for ( moris::uint Ik=0; Ik < tNumEquationObjectOnSet; Ik++ )
                {
                    mFemSets( Ii )->get_equation_object_list()( Ik )->compute_dQIdp_explicit();
                }

                //this->free_block_memory( Ii );
            }

            // global assembly to switch entries to the right processor
            mExplicitdQidu->vector_global_asembly();
        }

//-------------------------------------------------------------------------------------------------

        sol::Dist_Vector * Equation_Model::get_dQidu()
        {
            moris::Matrix_Vector_Factory tMatFactory( sol::MapType::Epetra );

            uint tNumRHMS = this->get_num_rhs();

            // full vector and prev full vector
            mQidu = tMatFactory.create_vector( mdQiduMap, tNumRHMS );

            mQidu->vec_put_scalar( 0.0 );

//            mExplicitdQidu->print();
//            mImplicitdQidu->print();

            mQidu->vec_plus_vec( 1.0, *mExplicitdQidu, 1.0 );
            mQidu->vec_plus_vec( 1.0, *mImplicitdQidu, 1.0 );

//            mQidu->print();

            return mQidu;
        }

//-------------------------------------------------------------------------------------------------

    }/* end_namespace_msi */
}/* end_namespace_moris */
