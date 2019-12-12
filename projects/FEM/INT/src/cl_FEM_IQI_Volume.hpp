/*
 * cl_FEM_IQI_Volume.hpp
 *
 *  Created on: Nov 20, 2019
 *      Author: noel
 */

#ifndef PROJECTS_FEM_INT_SRC_CL_FEM_IQI_VOLUME_HPP_
#define PROJECTS_FEM_INT_SRC_CL_FEM_IQI_VOLUME_HPP_

#include "typedefs.hpp"                     //MRS/COR/src
#include "cl_Cell.hpp"                      //MRS/CON/src

#include "cl_Matrix.hpp"                    //LINALG/src
#include "linalg_typedefs.hpp"              //LINALG/src

#include "cl_FEM_IQI.hpp"                   //FEM/INT/src

namespace moris
{
    namespace fem
    {
//------------------------------------------------------------------------------

        class IQI_Volume : public IQI
        {
//------------------------------------------------------------------------------
            public:
//------------------------------------------------------------------------------
            /*
             *  constructor
             */
            IQI_Volume(){};

//------------------------------------------------------------------------------
            /**
             * trivial destructor
             */
            ~IQI_Volume(){};

//------------------------------------------------------------------------------
            /**
             * compute the quantity of interest
             * @param[ in ] aQI quantity of interest matrix to fill
             */
            void compute_QI( Matrix< DDRMat > & aQI );

//------------------------------------------------------------------------------
        };
    }/* end namespace fem */
} /* end namespace moris */

#endif /* PROJECTS_FEM_INT_SRC_CL_FEM_IQI_VOLUME_HPP_ */
