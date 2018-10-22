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
//#include "cl_MSI_Model_Solver_Interface.hpp"  //FEM/MSI/src
#include "cl_MSI_Equation_Object.hpp"
namespace moris
{
//------------------------------------------------------------------------------

    namespace fem
    {
        class IWG;
        class Node_Base;
    }

    namespace NLA
	{
    	class Nonlinear_Solver;
    	class Nonlinear_Problem;
	}

    namespace mdl
    {
//------------------------------------------------------------------------------

        class Model
        {

            Cell< fem::Node_Base* >       mNodes;
            Cell< MSI::Equation_Object* > mElements;

//------------------------------------------------------------------------------
        public:
//------------------------------------------------------------------------------

           /**
            * simple constructor
            * @param[ in ] aMesh  Mesh for this problem
            * @param[ in ] aIWG   Integrant Weak form of Governing Equation
            */
//           Model(
//        		   NLA::Nonlinear_Problem * aNonlinearProblem,
//        		   std::shared_ptr< NLA::Nonlinear_Solver > aSolver,
//                   mtk::Mesh         * aMesh,
//                   fem::IWG          & aIWG,
//                   const Matrix< DDRMat > & aWeakBCs,
//                         Matrix< DDRMat > & aDOFs );

            Model(
                    mtk::Mesh         * aMesh,
                    fem::IWG          & aIWG,
                    const Matrix< DDRMat > & aWeakBCs,
                          Matrix< DDRMat > & aDOFs );

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
        };
//------------------------------------------------------------------------------
    } /* namespace mdl */
} /* namespace moris */


#endif /* PROJECTS_FEM_MDL_SRC_CL_MDL_MODEL_HPP_ */
