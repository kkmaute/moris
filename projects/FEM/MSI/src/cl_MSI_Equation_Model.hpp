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

            // matrices for current and previous time slabs
            Matrix< DDRMat > mTime;
            Matrix< DDRMat > mPrevTime;

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
        };
//------------------------------------------------------------------------------
    } /* namespace MSI */
} /* namespace moris */


#endif /* PROJECTS_FEM_MDL_SRC_CL_MSI_MODEL_HPP_ */
