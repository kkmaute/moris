/*
 * cl_FEM_Property_Temp_Dirichlet.hpp
 *
 *  Created on: Jul 12, 2019
 *      Author: noel
 */
#ifndef SRC_FEM_CL_FEM_PROPERTY_TEMP_DIRICHLET_HPP_
#define SRC_FEM_CL_FEM_PROPERTY_TEMP_DIRICHLET_HPP_

#include "typedefs.hpp"                     //MRS/COR/src
#include "cl_Cell.hpp"                      //MRS/CON/src
#include "cl_Matrix.hpp"                    //LINALG/src
#include "linalg_typedefs.hpp"              //LINALG/src
#include "cl_FEM_Property.hpp"                   //FEM/INT/src

namespace moris
{
    namespace fem
    {
//------------------------------------------------------------------------------

        class Property_Temp_Dirichlet : public Property
        {


//------------------------------------------------------------------------------
        public:
//------------------------------------------------------------------------------
            //constructor
            Property_Temp_Dirichlet();

//------------------------------------------------------------------------------
            // trivial destructor
            ~Property_Temp_Dirichlet(){};

//------------------------------------------------------------------------------

            void val_coeff( Matrix< DDRMat > & aCoeff );

//------------------------------------------------------------------------------

            void val( Matrix< DDRMat > & aCoeff,
                      Matrix< DDRMat > & aSpacePhysCoords,
                      Matrix< DDRMat > & aTimePhysCoords )
            {
                MORIS_ERROR(false, "Property_Temp_Dirichlet::val - This function is not implemented for this property type");
            }

        private:

        };
//------------------------------------------------------------------------------
    } /* namespace fem */
} /* namespace moris */

#endif /* SRC_FEM_CL_FEM_PROPERTY_TEMP_DIRICHLET_HPP_ */
