/*
 * cl_GEN_Field.hpp
 *
 *  Created on: Nov 5, 2019
 *      Author: sonne
 */

#ifndef PROJECTS_GEN_SRC_NEW_FIELD_CL_GEN_FIELD_HPP_
#define PROJECTS_GEN_SRC_NEW_FIELD_CL_GEN_FIELD_HPP_

namespace moris{
namespace ge{

/*
 * @brief the geometry engine performs operations on fields which are not geometry (e.g. density field)
 *
 * this is the base class for said fields
 */

class GEN_Field
{
public:
    GEN_Field(  )
    {

    }

    ~GEN_Field(  )
    {

    }
    //------------------------------------------------------------------------------
    virtual real eval_function( Matrix< DDRMat > aParam )
    {
        MORIS_ASSERT( false, "GEN_Field::eval_function() - field not created, returning 0" );
        return 0;
    }

};

}   // end ge namespace
}   // end moris namespace



#endif /* PROJECTS_GEN_SRC_NEW_FIELD_CL_GEN_FIELD_HPP_ */
