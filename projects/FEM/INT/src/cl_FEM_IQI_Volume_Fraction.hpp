/*
 * cl_FEM_IQI_Volume_Fraction.hpp
 *
 *  Created on: Nov 20, 2019
 *      Author: schmidt
 */

#ifndef PROJECTS_FEM_INT_SRC_CL_FEM_IQI_VOLUME_FRACTION_HPP_
#define PROJECTS_FEM_INT_SRC_CL_FEM_IQI_VOLUME_FRACTION_HPP_

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

        class IQI_Volume_Fraction : public IQI
        {
//------------------------------------------------------------------------------
            public:
//------------------------------------------------------------------------------
            /*
             *  constructor
             */
            IQI_Volume_Fraction()
            {
                // set IQI type
                mIQIType = vis::Output_Type::VOLUME_FRACTION;
            };

//------------------------------------------------------------------------------
            /**
             * trivial destructor
             */
            ~IQI_Volume_Fraction(){};

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
