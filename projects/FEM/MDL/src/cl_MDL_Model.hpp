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

#include "cl_Mat.hpp"                         // LNA/src
#include "cl_MTK_Mesh.hpp"                    //MTK/src
//#include "cl_MSI_Model_Solver_Interface.hpp"  //FEM/MSI/src
#include "cl_FEM_IWG.hpp"                     //FEM/INT/src
#include "cl_MSI_Node.hpp"
#include "cl_MSI_Equation_Object.hpp"
namespace moris
{
    namespace mdl
    {
//------------------------------------------------------------------------------

        class Model
        {
//------------------------------------------------------------------------------
        public:
//------------------------------------------------------------------------------




           /**
            * simple constructor
            * @param[ in ] aMesh  Mesh for this problem
            * @param[ in ] aIWG   Integrant Weak form of Governing Equation
            */
           Model(
                   mtk::Mesh         & aMesh,
                   fem::IWG          & aIWG,
                   const Mat< real > & aWeakBCs,
                         Mat< real > & aDOFs );

//------------------------------------------------------------------------------

           /**
            * destructor
            */
           ~Model();

        };
//------------------------------------------------------------------------------
    } /* namespace mdl */
} /* namespace moris */


#endif /* PROJECTS_FEM_MDL_SRC_CL_MDL_MODEL_HPP_ */
