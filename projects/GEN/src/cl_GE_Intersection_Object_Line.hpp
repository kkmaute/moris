/*
 * cl_GE_Intersection_Object_Line.hpp
 *
 *  Created on: May 23, 2019
 *      Author: sonne
 */

#ifndef PROJECTS_GEN_SRC_CL_GE_INTERSECTION_OBJECT_LINE_HPP_
#define PROJECTS_GEN_SRC_CL_GE_INTERSECTION_OBJECT_LINE_HPP_

// GE includes
#include "cl_GE_Intersection_Object.hpp"

// MTK includes
#include "cl_MTK_Enums.hpp"

namespace moris
{
namespace ge
{

class Intersection_Object_Line : public Intersection_Object
{
    /*
     * @brief The intersection object is passed to the geometry engine and asked if there is an intersection. If there is, the information gets stored in the
     *        intersection object as member data.
     *
     *        The line type in space and time is defaulted to:
     *                                  Interpolation_Type: Legrange
     *                                  Interpolation_Rule: Linear
     */
public:
    Intersection_Object_Line(  )
    {
        mtk::Geometry_Type       tGeomType        = mtk::Geometry_Type::LINE;
        fem::Interpolation_Type  tGeomInterpType  = fem::Interpolation_Type::LAGRANGE;
        mtk::Interpolation_Order tGeomInterpOrder = mtk::Interpolation_Order::LINEAR;
        fem::Interpolation_Type  tTimeInterpType  = fem::Interpolation_Type::LAGRANGE;
        mtk::Interpolation_Order tTimeInterpOrder = mtk::Interpolation_Order::LINEAR;

        fem::Interpolation_Rule tGeomInterpRule( tGeomType, tGeomInterpType, tGeomInterpOrder, tTimeInterpType, tTimeInterpOrder );
        mMyGeomInterp = new fem::Geometry_Interpolator( tGeomInterpRule );

        mMyIntersectionFlag = false; // initially assume there are no intersection points
    };
    ~Intersection_Object_Line()
    {
        delete mMyGeomInterp;
        delete mMyFieldInterp;
    };

    //------------------------------------------------------------------------------
    void is_intersected()
    {
        mMyIntersectionFlag = true;
    }
    //******************************* set functions ********************************
    //------------------------------------------------------------------------------
    void set_coords_and_param_point( std::shared_ptr< Geometry > & aGeomPointer,
                                     Matrix<DDRMat> const & aGlobCoords,
                                     Matrix<DDRMat> const & aTimeCoords,
                                     Matrix<DDRMat> const & aFieldVals,     // LS field vals need to come from the geometry engine
                                     Matrix<DDRMat> const & aParamPoint = Matrix< DDRMat >({{-1},{0}}) )  // start at the beginning of the line(local coordinate xi=-1) and at t=0
    {
        fem::Interpolation_Rule tLevelSetInterpRule( this->get_my_geom_type(),
                                                     aGeomPointer->get_my_space_interpolation_type(),
                                                     aGeomPointer->get_my_space_interpolation_order(),
                                                     aGeomPointer->get_my_time_interpolation_type(),
                                                     aGeomPointer->get_my_time_interpolation_order() );
        //fixme this is a default to just look at one LS field at a time, need to update this to look at all the LS fields,
        //      additionally, the geometry pointer information needs to be accessed in a better way
        uint tNumFields = 1;
        mMyFieldInterp = new fem::Field_Interpolator( tNumFields, tLevelSetInterpRule, mMyGeomInterp );

        mMyParamPoint   = aParamPoint; // set the parameter point where the interpolation begins
        mMyGlobalCoords = aGlobCoords;
        mMyTimeCoords   = aTimeCoords;
        mMyFieldVals    = aFieldVals;

        mMyGeomInterp->set_coeff( mMyGlobalCoords, mMyTimeCoords );
        mMyFieldInterp->set_coeff( mMyFieldVals );
    }
    //------------------------------------------------------------------------------
    void set_intersection_point( Matrix< DDRMat > aPoint )
    {
        mMyIntersectionPoints = aPoint;  // currently only set up for one intersection point
    }

    //******************************* get functions ********************************
    //------------------------------------------------------------------------------
    mtk::Geometry_Type
    get_my_geom_type()
    {
        return mtk::Geometry_Type::LINE;
    }
    //------------------------------------------------------------------------------
    fem::Geometry_Interpolator*
    get_my_geom_interp()
    {
        return mMyGeomInterp;
    }
    //------------------------------------------------------------------------------
    fem::Field_Interpolator*
    get_my_field_interp()
    {
        return mMyFieldInterp;
    }
    //------------------------------------------------------------------------------
    Matrix< DDRMat >
    get_my_param_point()
    {
        return mMyParamPoint;
    }
    //------------------------------------------------------------------------------
    Matrix< DDRMat >
    get_my_global_coord()
    {
        return mMyGlobalCoords;
    }
    //------------------------------------------------------------------------------
    Matrix< DDRMat >
    get_my_time_coord()
    {
        return mMyTimeCoords;
    }
    //------------------------------------------------------------------------------
    Matrix< DDRMat >
    get_my_field_vals()
    {
        return mMyFieldVals;
    }
    //------------------------------------------------------------------------------
    bool
    get_intersection_flag()
    {
        return mMyIntersectionFlag;
    }
    //------------------------------------------------------------------------------
    Matrix< DDRMat >
    get_intersection_point()
    {
        return mMyIntersectionPoints;    // currently only set for a single intersection point
    }

//------------------------------------------------------------------------------
private:
    //------------------------------------------------------------------------------

    fem::Geometry_Interpolator* mMyGeomInterp  = nullptr;
    fem::Field_Interpolator*    mMyFieldInterp = nullptr;

    Matrix<DDRMat> mMyGlobalCoords;
    Matrix<DDRMat> mMyTimeCoords;
    Matrix<DDRMat> mMyFieldVals;
    Matrix<DDRMat> mMyParamPoint;

    bool mMyIntersectionFlag;

    Matrix< DDRMat > mMyIntersectionPoints;
//------------------------------------------------------------------------------
protected:

};

}// end ge namespace
}// end moris namespace
#endif /* PROJECTS_GEN_SRC_CL_GE_INTERSECTION_OBJECT_LINE_HPP_ */
