/*
 * cl_GE_Output_Object.hpp
 *
 *  Created on: May 1, 2019
 *      Author: sonne
 */

#ifndef PROJECTS_GEN_SRC_CL_GE_OUTPUT_OBJECT_HPP_
#define PROJECTS_GEN_SRC_CL_GE_OUTPUT_OBJECT_HPP_

namespace moris
{
namespace ge
{
    class Output_Object
    {

    public:
        Output_Object(){};

        Output_Object( Matrix< DDRMat > &aNodalADVs ) : mNodalADVs(aNodalADVs)
        {

        };
        ~Output_Object(){};

    private:
        Matrix< DDRMat > mNodalADVs;

    protected:

    };

} //end ge namespace
} //end moris namespace


#endif /* PROJECTS_GEN_SRC_CL_GE_OUTPUT_OBJECT_HPP_ */
