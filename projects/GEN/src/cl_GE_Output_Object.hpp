/*
 * cl_GE_Output_Object.hpp
 *
 *  Created on: May 1, 2019
 *      Author: sonne
 */

#ifndef PROJECTS_GEN_SRC_CL_GE_OUTPUT_OBJECT_HPP_
#define PROJECTS_GEN_SRC_CL_GE_OUTPUT_OBJECT_HPP_

#include "cl_GE_Core.hpp"

#include "cl_MTK_Mapper.hpp"

namespace moris
{
namespace ge
{
    class Output_Object
    {

    public:
        Output_Object()
        {};

        Output_Object( mapper::Mapper & aMyMapper ) :
                     mMyMapper(&aMyMapper)
        {};

        Output_Object( Matrix< DDRMat > & aNodalADVs, Matrix< DDRMat > & aFieldVals ) :
                     mNodalADVs(aNodalADVs),
                     mFieldVals(aFieldVals)
        {};

        ~Output_Object(){};

        //------------------------------------------------------------------------------
        /*
         * @brief returns the nodal field values
         */
        Matrix< DDRMat >
        get_field_vals()
        {
            return mFieldVals;
        }
        //------------------------------------------------------------------------------
        /*
         * @brief returns the nodal ADVs
         */
        Matrix< DDRMat >
        get_nodal_advs()
        {
            return mNodalADVs;
        }

//------------------------------------------------------------------------------
    private:
        Matrix< DDRMat > mNodalADVs;
        Matrix< DDRMat > mFieldVals;

        mapper::Mapper*     mMyMapper = nullptr;

//------------------------------------------------------------------------------
    protected:

    };

} //end ge namespace
} //end moris namespace


#endif /* PROJECTS_GEN_SRC_CL_GE_OUTPUT_OBJECT_HPP_ */
