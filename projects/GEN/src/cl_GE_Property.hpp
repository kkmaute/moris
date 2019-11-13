/*
 * cl_GE_Property.hpp
 *
 *  Created on: Sep 25, 2019
 *      Author: sonne
 */

#ifndef PROJECTS_GEN_SRC_CL_GE_PROPERTY_HPP_
#define PROJECTS_GEN_SRC_CL_GE_PROPERTY_HPP_

#include "typedefs.hpp"                     //MRS/COR/src
#include "cl_Matrix.hpp"                    //LNA/src
#include "cl_Cell.hpp"

// FEM includes
#include "cl_FEM_Field_Interpolator.hpp"    //FEM/INT/src

namespace moris{
namespace ge{

typedef std::function< moris::real ( const Matrix< DDRMat >         & aPoint,
                                     const moris::Cell< moris::real > aProps ) > PropertyFunction;

typedef std::function< Matrix< DDRMat > ( const Matrix< DDRMat >    & aPoint,
                                          const moris::Cell< real >   aProps ) > PropertyDerFunc;
//------------------------------------------------------------------------------
class Property
{
public:
    //------------------------------------------------------------------------------
    //constructor to set both the analytical function and the derivative(s)
    Property( PropertyFunction aValFunc,
              PropertyDerFunc  aDerFuncs,
              moris::Cell< moris::real > aParams )
    : mValFunction(aValFunc),
      mDerivativeFunction(aDerFuncs),
      mParams(aParams)
    {
    }

    // constructor to set only the analytical function
    Property( PropertyFunction aValFunc,
              moris::Cell< moris::real > aParams )
    : mValFunction(aValFunc),
      mParams(aParams)
    {
    }

    ~Property()
    {}

    //------------------------------------------------------------------------------
    bool is_derivative_set(  )
    {
        bool tSet = false;
        if ( mDerivativeFunction == nullptr )
        {
            tSet = false;
        }
        else
        {
            tSet = true;
        }
        return tSet;
    }
    //------------------------------------------------------------------------------
    void set_derivative_functions( PropertyDerFunc aDerFuncs )
    {
        mDerivativeFunction = aDerFuncs;
    }
    //------------------------------------------------------------------------------
    void set_params( moris::Cell< real > aParams )
    {
        mParams = aParams;
    }
    //------------------------------------------------------------------------------
    moris::real get_field_val_at_coordinate( const Matrix< DDRMat >  & aPoint )
    {
        return mValFunction( aPoint, mParams );
    };
    //------------------------------------------------------------------------------
    Matrix< DDRMat > get_sensitivity_dphi_dp_at_coordinate( const Matrix< DDRMat >  & aPoint )
    {
        MORIS_ASSERT( mDerivativeFunction != nullptr, "ge::Property::get_sensitivity_dphi_dp_at_coordinate() - the derivative function is not set " );
        return mDerivativeFunction( aPoint, mParams );
    };

    //------------------------------------------------------------------------------
private:
    //------------------------------------------------------------------------------
    // note: all DOF types are ADVs

    // value function
    PropertyFunction mValFunction = nullptr;

    // derivative functions
    PropertyDerFunc mDerivativeFunction = nullptr;

    // field interpolators
    fem::Field_Interpolator* mFieldInterpolators = nullptr;

    // geometry interpolator
    fem::Geometry_Interpolator* mGeometryInterpolator = nullptr;

    // parameters
    moris::Cell< moris::real > mParams;     // constant values to be used in the functions

    //------------------------------------------------------------------------------
protected:
    //------------------------------------------------------------------------------

    //------------------------------------------------------------------------------
};

}   // end ge namespace
}   // end moris namespace


#endif /* PROJECTS_GEN_SRC_CL_GE_PROPERTY_HPP_ */
