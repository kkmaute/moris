/*
 * cl_MDL_Model.hpp
 *
 *  Created on: Aug 22, 2018
 *      Author: messe
 */

#ifndef PROJECTS_FEM_MDL_SRC_CL_MDL_MODEL_HPP_
#define PROJECTS_FEM_MDL_SRC_CL_MDL_MODEL_HPP_

#include "typedefs.hpp"                       //MRS/COR/src
#include "cl_Cell.hpp"                        //MRS/CON/src

#include "cl_Matrix.hpp"
#include "linalg_typedefs.hpp"
#include "cl_MTK_Mesh.hpp"                    //MTK/src
#include "cl_MTK_Enums.hpp"

//#include "cl_MSI_Model_Solver_Interface.hpp"  //FEM/MSI/src
#include "cl_MSI_Equation_Object.hpp"
namespace moris
{
//------------------------------------------------------------------------------

    namespace fem
    {
        class IWG;
        class Node_Base;
        enum class IWG_Type;
        enum class BC_Type;
    }

    namespace dla
    {
        class Linear_Solver;
        class Linear_Solver_Algorithm;
    }

    namespace NLA
	{
    	class Nonlinear_Algorithm;
    	class Nonlinear_Problem;
    	class Nonlinear_Solver;
	}

    namespace MSI
    {
        class Model_Solver_Interface;
        class MSI_Solver_Interface;
    }
    namespace mdl
    {
//------------------------------------------------------------------------------

        class Model
        {
            // pointer to reference mesh
            mtk::Mesh                       * mMesh;
            Cell< fem::Node_Base* >           mNodes;
            Cell< MSI::Equation_Object* >     mElements;
            Cell< Cell< fem::IWG* > >         mIWGs;

            Cell< fem::IWG* >         mIWGs1;

            // by default, this value is set to the order of the
            // Lagrange modes
            moris::uint                       mDofOrder = 0;

            MSI::Model_Solver_Interface                   * mModelSolverInterface;
            MSI::MSI_Solver_Interface                     * mSolverInterface;
            NLA::Nonlinear_Problem                        * mNonlinearProblem;
            NLA::Nonlinear_Solver                         * mNonlinearSolver;
            std::shared_ptr< NLA::Nonlinear_Algorithm >     mNonlinearSolverAlgorithm;
            std::shared_ptr< dla::Linear_Solver_Algorithm > mLinearSolverAlgorithm;
            dla::Linear_Solver                            * mLinSolver;

            // fixme: maybe introduce a cell of maps for different orders?
            map< moris_id, moris_index >      mCoefficientsMap;
            Matrix< DDUMat >                  mAdofMap;

            Matrix< DDRMat> mSolHMR;

//------------------------------------------------------------------------------
        public:
//------------------------------------------------------------------------------

           /**
            * simple constructor
            * @param[ in ] aMesh  Mesh for this problem
            * @param[ in ] aIWG   Integrant Weak form of Governing Equation
            */
//            Model(       mtk::Mesh   * aMesh,
//                   const uint          aBSplineOrder,
//                   Cell< Cell< fem::IWG_Type > >aIWGTypeList );

            Model(       mtk::Mesh *                   aMesh,
                   const uint                          aBSplineOrder,
                         Cell< Cell< fem::IWG_Type > > aIWGTypeList,
                         Cell< moris_index >           aSidesetList,
                         Cell< fem::BC_Type >          aSidesetBCTypeList );
//------------------------------------------------------------------------------

            Matrix< DDRMat> &
            get_mSolHMR( )
            {
                return mSolHMR;
            };

            void
            set_dof_order( const uint aOrder );

//------------------------------------------------------------------------------

            uint
            get_lagrange_order_from_mesh();

//------------------------------------------------------------------------------

            void
            set_weak_bcs( const Matrix< DDRMat > & aWeakBCs );

//------------------------------------------------------------------------------

            void
            set_weak_bcs_from_nodal_field( moris_index aFieldIndex );

//------------------------------------------------------------------------------

            void
            solve( Matrix< DDRMat > & aSolution );

//------------------------------------------------------------------------------

           /**
            * destructor
            */
           ~Model();

//------------------------------------------------------------------------------

           real
           compute_integration_error(
                               real (*aFunction)( const Matrix< DDRMat > & aPoint ) );

//------------------------------------------------------------------------------

           real
           compute_element_average( const uint aElementIndex );

//------------------------------------------------------------------------------

           void output_solution( const std::string & aFilePath );

//------------------------------------------------------------------------------

           void
           output_solution_nils_HACK( const std::string & aFilePath, mtk::Mesh * aMesh );

//------------------------------------------------------------------------------
        };
//------------------------------------------------------------------------------
    } /* namespace mdl */
} /* namespace moris */


#endif /* PROJECTS_FEM_MDL_SRC_CL_MDL_MODEL_HPP_ */
