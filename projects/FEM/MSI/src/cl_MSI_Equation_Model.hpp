/*
 * cl_MSI_Model.hpp
 *
 *  Created on: Aug 22, 2018
 *      Author: schmidt
 */

#ifndef PROJECTS_FEM_MDL_SRC_CL_MSI_MODEL_HPP_
#define PROJECTS_FEM_MDL_SRC_CL_MSI_MODEL_HPP_

#include "typedefs.hpp"                       //MRS/COR/src
#include "cl_Cell.hpp"                        //MRS/CON/src

#include "cl_Matrix.hpp"
#include "linalg_typedefs.hpp"

namespace moris
{
    class Dist_Vector;
//------------------------------------------------------------------------------

    namespace MSI
    {
        class Model_Solver_Interface;
        class MSI_Solver_Interface;
        class Equation_Set;
        class Equation_Object;
        class Design_Variable_Interface;
        enum class Dof_Type;
//------------------------------------------------------------------------------

        class Equation_Model
        {
        protected:

            // list of FEM sets
            moris::Cell< MSI::Equation_Set * > mFemSets;

            // list of FEM clusters
            moris::Cell< MSI::Equation_Object* > mFemClusters;

            // map from mesh set indices to fem set indices
            //map< moris_index, moris_index >   mMeshSetToFemSetMap;
            map< std::pair< moris_index, bool >, moris_index > mMeshSetToFemSetMap;

            // distributed solution vectors for current and previous time slabs
            Dist_Vector * mSolutionVector     = nullptr;
            Dist_Vector * mPrevSolutionVector = nullptr;
            Dist_Vector * mSensitivitySolutionVector = nullptr;

            Dist_Vector * mdQidu = nullptr;

            // matrices for current and previous time slabs
            Matrix< DDRMat > mTime;
            Matrix< DDRMat > mPrevTime;

            MSI::Design_Variable_Interface * mDesignVariableInterface = nullptr;

            bool mIsForwardAnalysis = true;

            moris::sint mNumSensitivityAnalysisRHS = -1;

//------------------------------------------------------------------------------
            // Dummy Variables
            moris::Cell< std::string > mDummy;

//------------------------------------------------------------------------------
        public:
//------------------------------------------------------------------------------
            /**
             * constructor
             */
            Equation_Model(){};

//------------------------------------------------------------------------------
            /**
             * destructor
             */
            virtual ~Equation_Model(){};

//------------------------------------------------------------------------------
            /**
             * get equation sets for test
             */
            moris::sint get_num_rhs( )
            {
                if( !mIsForwardAnalysis )
                {
                    mNumSensitivityAnalysisRHS = this->get_requested_IQI_names().size();

                    MORIS_ASSERT( mNumSensitivityAnalysisRHS <= 0, "MSI::Equation_Model::get_num_rhs(), num rhs not set for sensitivity analysis");

                    return mNumSensitivityAnalysisRHS;
                }
                else
                {
                    return 1;
                }
            };

//------------------------------------------------------------------------------
            /**
             * get equation sets for test
             */
            moris::Cell< MSI::Equation_Set * > & get_equation_sets( )
            {
                return mFemSets;
            };

//------------------------------------------------------------------------------
            /**
             * get equation objects
             */
            moris::Cell< MSI::Equation_Object * > & get_equation_objects( )
            {
                return mFemClusters;
            };

//------------------------------------------------------------------------------
            /**
             * MTK set to fem set index map
             */
            //map< moris_index, moris_index > & get_mesh_set_to_fem_set_index_map( )
            map< std::pair< moris_index, bool >, moris_index > & get_mesh_set_to_fem_set_index_map( )
            {
                return mMeshSetToFemSetMap;
            };

//------------------------------------------------------------------------------
            /**
             * finalize the equation sets
             * @param[ in ] aModelSolverInterface pointer to a model solver interface
             */
            virtual void finalize_equation_sets( MSI::Model_Solver_Interface * aModelSolverInterface )
            {
                MORIS_ERROR( false, "Equation_Model::finalize_equation_sets - not implemented for base class." );
            }

//------------------------------------------------------------------------------
            /**
             * set solution vector
             * @param[ in ] aSolutionVector distributed solution vector
             */
            void set_solution_vector( Dist_Vector * aSolutionVector )
            {
                mSolutionVector = aSolutionVector;
            }

//------------------------------------------------------------------------------
            /**
             * get solution vector
             * @param[ out ] aSolutionVector distributed solution vector
             */
            Dist_Vector * get_solution_vector()
            {
                return mSolutionVector;
            }

//------------------------------------------------------------------------------
            /**
             * set previous solution vector
             * @param[ in ] aSolutionVector previous distributed solution vector
             */
            void set_previous_solution_vector( Dist_Vector * aSolutionVector )
            {
                mPrevSolutionVector = aSolutionVector;
            }

//------------------------------------------------------------------------------
            /**
             * get previous solution vector
             * @param[ out ] aSolutionVector previous distributed solution vector
             */
            Dist_Vector * get_previous_solution_vector()
            {
                return mPrevSolutionVector;
            }

//------------------------------------------------------------------------------
            /**
             * set sensitivity solution vector
             * @param[ in ] aSensitivitySolutionVector distributed solution vector for sensitivity
             */
            void set_adjoint_solution_vector( Dist_Vector * aSolutionVector )
            {
                mSensitivitySolutionVector = aSolutionVector;
            }

//------------------------------------------------------------------------------
            /**
             * get adjoint solution vector
             * @param[ out ] aSolutionVector adjoint distributed solution vector
             */
            Dist_Vector * get_adjoint_solution_vector()
            {
                return mSensitivitySolutionVector;
            }

//------------------------------------------------------------------------------
            /**
             * set time for current time slab
             * @param[ in ] aTime matrix for time in current time slab
             */
            void set_time( Matrix< DDRMat > & aTime )
            {
                mTime = aTime;
            }

//------------------------------------------------------------------------------
            /**
             * get time for current time slab
             * @param[ out ] mTime matrix for time in current time slab
             */
            Matrix< DDRMat > & get_time()
            {
                return mTime;
            }

//------------------------------------------------------------------------------
            /**
             * set time for previous time slab
             * @param[ in ] aPrevTime matrix for time in previous time slab
             */
            void set_previous_time( Matrix< DDRMat > & aPrevTime )
            {
                mPrevTime = aPrevTime;
            }

//------------------------------------------------------------------------------

            /**
             * get time for previous time slab
             * @param[ out ] mPrevTime matrix for time in previous time slab
             */
            Matrix< DDRMat > & get_previous_time()
            {
                return mPrevTime;
            }

//------------------------------------------------------------------------------
            /**
             * set pointer to design variable interface
             * @param[ out ] aDesignVariableInterface pointer to design variable interface
             */
            void set_design_variable_interface( MSI::Design_Variable_Interface * aDesignVariableInterface )
            {
                mDesignVariableInterface = aDesignVariableInterface;
            }


//------------------------------------------------------------------------------
            /**
             * set requested IQI names
             * @param[ in ] aRequestedIQINames List of requested IQI names
             */
            virtual void set_requested_IQI_names( const moris::Cell< std::string > & aRequestedIQINames )
            {
                MORIS_ERROR( false, "Equation_Model::set_requested_IQI_names - not implemented for base class." );
            }

//------------------------------------------------------------------------------
            /**
             * get requested IQI names
             */
            virtual const moris::Cell< std::string > & get_requested_IQI_names()
            {
                MORIS_ERROR( false, "Equation_Model::get_requested_IQI_names - not implemented for base class." );
                return mDummy;
            }

//------------------------------------------------------------------------------
            /**
             * indicated that this equation model is used for the sensitivity analysis
             */
            void set_is_sensitivity_analysis()
            {
                mIsForwardAnalysis = false;
            };

//------------------------------------------------------------------------------
            /**
             * indicated that this equation model is used for the forward analysis
             */
            void set_is_forward_analysis()
            {
                mIsForwardAnalysis = true;
            };

//------------------------------------------------------------------------------
            /**
             * returns if this is the a forward analysis
             * @param[ out ] mIsForwardAnalysis
             */
            bool get_is_forward_analysis()
            {
                return mIsForwardAnalysis;
            };

//------------------------------------------------------------------------------
            /**
             * retruns if this is the a forward analysis
             */
            void compute_dQIdp()
            {
//                // Get local number of elements
//                moris::uint tNumSets = mFemSets.size();
//
//                moris::uint tNumRHS = this->get_num_rhs();
//
//                // Loop over all local elements to build matrix graph
//                for ( moris::uint Ii=0; Ii < tNumSets; Ii++ )
//                {
//                    moris::uint tNumEquationObjectOnSet = mFemSets( Ii )->get_num_equation_objects();
//
//                    this->initialize_set( Ii, true );
//
//                    for ( moris::uint Ik=0; Ik < tNumEquationObjectOnSet; Ik++ )
//                    {
//                        mFemSets( Ii )->get_equation_object_list()( Ik )->compute_dQIdp();
//                    }
//
//                    this->free_block_memory( Ii );
//                }
            };

//------------------------------------------------------------------------------
        };
//------------------------------------------------------------------------------
    } /* namespace MSI */
} /* namespace moris */


#endif /* PROJECTS_FEM_MDL_SRC_CL_MSI_MODEL_HPP_ */
