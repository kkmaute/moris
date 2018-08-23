/*
 * cl_MDL_Model.hpp
 *
 *  Created on: Aug 22, 2018
 *      Author: messe
 */

#ifndef PROJECTS_FEM_MDL_SRC_CL_MDL_MODEL_HPP_
#define PROJECTS_FEM_MDL_SRC_CL_MDL_MODEL_HPP_

#include "typedefs.hpp"                     //MRS/COR/src
#include "cl_Cell.hpp"                      //MRS/CON/src

#include "cl_Mat.hpp"                       // LNA/src
#include "cl_MTK_Mesh.hpp"                  //MTK/src
#include "cl_Model_Solver_Interface.hpp"    //FEM/MSI/src
#include "cl_FEM_IWG.hpp"                   //FEM/INT/src
#include "cl_MSI_Node.hpp"                  //FEM/INT/src

namespace moris
{
    namespace mdl
    {
//------------------------------------------------------------------------------

        class Model
        {
           //! result vector
           Mat< real > & mResult;

           //! list of equation objects
           moris::Cell< MSI::Equation_Object * > mEquationObjects;

           //! list of nodes ( may be moved somewhere else )
           moris::Cell< MSI::Node * > mNodes;
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
                   const Mat< real > & aInput,
                         Mat< real > & aResult );

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
