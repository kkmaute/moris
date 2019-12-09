/*
 * cl_GEN_Property.hpp
 *
 *  Created on: Dec 4, 2019
 *      Author: sonne
 */

#ifndef PROJECTS_GEN_SRC_PROPERTY_CL_GEN_PROPERTY_HPP_
#define PROJECTS_GEN_SRC_PROPERTY_CL_GEN_PROPERTY_HPP_

// FEM includes
#include "cl_FEM_Field_Interpolator.hpp"

namespace moris
{
namespace ge
{
typedef std::function< Matrix< DDRMat > ( moris::Cell< Matrix< DDRMat > >         & aCoeff,
                                          moris::Cell< fem::Field_Interpolator* > & aDofFI,
                                          moris::Cell< fem::Field_Interpolator* > & aDvFI,
                                          fem::Geometry_Interpolator              * aGeometryInterpolator ) > PropertyFunc;
//------------------------------------------------------------------------------
class GEN_Property
{
public:
    GEN_Property(){};

    ~GEN_Property(){};
    //------------------------------------------------------------------------------
    /**
     * set val function
     * @param[ in ] aValFunction function for property evaluation
     */
    void set_val_function( PropertyFunc aValFunction )
    {
        mValFunction = aValFunction;
    };
    //------------------------------------------------------------------------------


    //------------------------------------------------------------------------------
private:

    //------------------------------------------------------------------------------
protected:

    //------------------------------------------------------------------------------
private:    // member data
    PropertyFunc mValFunction = nullptr;
};

}   // end ge namespace
}   // end moris namespace



#endif /* PROJECTS_GEN_SRC_PROPERTY_CL_GEN_PROPERTY_HPP_ */
