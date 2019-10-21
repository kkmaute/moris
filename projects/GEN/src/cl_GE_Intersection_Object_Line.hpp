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
     *                                  Interpolation_Type: Lagrange
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

        mMyIntersectionFlag = false; // initially assume it is not intersected
    };
    ~Intersection_Object_Line()
    {
        delete mMyGeomInterp;
        delete mMyFieldInterp;
    };
    //------------------------------------------------------------------------------
    void flag_as_intersected( )
    {
        mMyIntersectionFlag = true;
    }
    //------------------------------------------------------------------------------
    moris_index compute_intersection( )
    {//fixme what if there is no intersection?
        mMyIntersectionPoints.push_back({{(mMyFieldVals(0) + mMyFieldVals(1))/(mMyFieldVals(0) - mMyFieldVals(1))}});

        //        mMyIntersFieldSensVal = Cell<Matrix<DDRMat>>(mMyIntersectionPoints.size()-1, Matrix<DDRMat>(2,3));
        mMyIntersFieldSensVal.push_back( {{0},{0}} );
        return mMyIntersectionPoints.size()-1;
    }
    //------------------------------------------------------------------------------
    Matrix< DDRMat > get_intersection_point_global_coord( moris_index aMyIndex )
    {
        Matrix<moris::DDRMat> tN = {{0.5 * (1-mMyIntersectionPoints(aMyIndex)(0))},     // 0.5*(1-zeta)
                                    {0.5 * (1+mMyIntersectionPoints(aMyIndex)(0))}};    // 0.5*(1+zeta)
        return trans(tN)*this->get_my_global_coord();
    }
    //------------------------------------------------------------------------------
    /*
     * @brief - compute the sensitivity at an intersection point
     *
     * @param[in] aMyIndex - index to the intersection point from the cell of intersection points
     */
    void compute_intersection_sensitivity( moris_index aMyIndex )
    {
        MORIS_ASSERT( mMyIntersFieldSensVal.size() == mMyIntersectionPoints.size(),"ge::Intersection_Object_Line::compute_intersection_sensitivity() - sensitivity information not allocated");
//        MORIS_ASSERT(mMyIntersSensInd.size() == mMyIntersectionPoints.size(),"Sensitivity information not allocated");

        // Local Coordinate of edge
        moris::real tZeta = mMyIntersectionPoints(aMyIndex)(0);

        moris::real const & tPhiA = mMyFieldVals(0);
        moris::real const & tPhiB = mMyFieldVals(1);

        // linear shape function along edge
//        Matrix<moris::DDRMat> tN = {{0.5 * (1-tZeta)},
//                                    {0.5 * (1+tZeta)}};

        // Derivatives of shape functions wrt global coordinate
        Matrix<moris::DDRMat> tdN_dzeta = {{-0.5},
                                           { 0.5}};

        // denominator (used throughout)
        moris::real tDenom = 1/std::pow((tPhiA-tPhiB),2);

        // partial derivatives of the interface change wrt change
        moris::real tdzetagamma_dphiA = -2*tPhiB*tDenom;
        moris::real tdzetagamma_dphiB =  2*tPhiA*tDenom;

        // derivative of x_gamma wrt p (of edge node 0)
        Matrix<moris::DDRMat> tdxgamma_dphiA = trans(tdN_dzeta)*this->get_my_global_coord()*tdzetagamma_dphiA;
        Matrix<moris::DDRMat> tdxgamma_dphiB = trans(tdN_dzeta)*this->get_my_global_coord()*tdzetagamma_dphiB;

        // Compute dx/dphi
        mMyIntersFieldSensVal(aMyIndex).set_size(2,this->get_my_global_coord().n_cols());

        MORIS_ASSERT( this->get_my_global_coord().n_cols() == 2 || this->get_my_global_coord().n_cols() == 3, "ge::Intersection_Object_Line::compute_intersection_sensitivity() - invalid global coordinate vector (check that this is a row vector?) " );

        mMyIntersFieldSensVal(aMyIndex).get_row(0) = tdxgamma_dphiA.get_row(0);
        mMyIntersFieldSensVal(aMyIndex).get_row(1) = tdxgamma_dphiB.get_row(0);
        //fixme: add field value indexes for sparse storage
    }
    //------------------------------------------------------------------------------
    //******************************* set functions ********************************
    //------------------------------------------------------------------------------
    void set_coords_and_param_point( std::shared_ptr< Geometry > & aGeomPointer,
                                     Matrix<DDRMat> const & aGlobCoords,
                                     Matrix<DDRMat> const & aTimeCoords,
                                     Matrix<DDRMat> const & aFieldVals,     // LS field vals need to come from the geometry engine?
                                     Matrix<DDRMat> const & aParamPoint = Matrix< DDRMat >({{-1},
                                                                                            {0}}) )  // start at the beginning of the line(local coordinate xi=-1) and at t=0
    {
        fem::Interpolation_Rule tLevelSetInterpRule( this->get_my_geom_type(),
                                                     aGeomPointer->get_my_space_interpolation_type(),
                                                     aGeomPointer->get_my_space_interpolation_order(),
                                                     aGeomPointer->get_my_time_interpolation_type(),
                                                     aGeomPointer->get_my_time_interpolation_order() );
        //fixme the geometry pointer information needs to be accessed in a better way
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
    void set_intersection( Matrix< DDRMat > aPoint,
                           moris_index      aMyIndex )
    {
        mMyIntersectionPoints( aMyIndex ) = aPoint;
    }
    //------------------------------------------------------------------------------
    //******************************* get functions ********************************
    //------------------------------------------------------------------------------
    mtk::Geometry_Type get_my_geom_type()
    {
        return mtk::Geometry_Type::LINE;
    }
    //------------------------------------------------------------------------------
    fem::Geometry_Interpolator* get_my_geom_interp()
    {
        return mMyGeomInterp;
    }
    //------------------------------------------------------------------------------
    fem::Field_Interpolator* get_my_field_interp()
    {
        return mMyFieldInterp;
    }
    //------------------------------------------------------------------------------
    Matrix< DDRMat > get_my_param_point()
    {
        return mMyParamPoint;
    }
    //------------------------------------------------------------------------------
    Matrix< DDRMat > get_my_global_coord()
    {
        return mMyGlobalCoords;
    }
    //------------------------------------------------------------------------------
    Matrix< DDRMat > get_my_time_coord()
    {
        return mMyTimeCoords;
    }
    //------------------------------------------------------------------------------
    Matrix< DDRMat > get_my_field_vals()
    {
        return mMyFieldVals;
    }
    //------------------------------------------------------------------------------
    bool get_intersection_flag()
    {
        return mMyIntersectionFlag;
    }
    //------------------------------------------------------------------------------
    Matrix< DDRMat > get_intersection_point_local_coord( moris_index aMyIndex )
    {
        return mMyIntersectionPoints(aMyIndex);
    }
    //------------------------------------------------------------------------------
    uint get_num_intersection_point( )
    {
        return mMyIntersectionPoints.size();
    }
    //------------------------------------------------------------------------------
    Matrix< DDRMat > get_field_sensitivity_vals( moris_index aMyIndex )
    {
        MORIS_ASSERT( std::abs(mMyIntersFieldSensVal( aMyIndex )( 0,0 )) >= 0.0, "ge::Intersection_Object_Line::get_field_sensitivity_vals() - value at given index is not set " );
        return mMyIntersFieldSensVal( aMyIndex );
    }

    //***************** get functions for asserts/debugs ***************************
    //------------------------------------------------------------------------------
    moris::Cell< Matrix< DDRMat > > get_field_sens_vals_cell()
    {
        return mMyIntersFieldSensVal;
    }

//------------------------------------------------------------------------------
private:
    //------------------------------------------------------------------------------

    fem::Geometry_Interpolator* mMyGeomInterp  = nullptr;
    fem::Field_Interpolator*    mMyFieldInterp = nullptr;

    Matrix<DDRMat> mMyGlobalCoords;
    Matrix<DDRMat> mMyTimeCoords;
    Matrix<DDRMat> mMyFieldVals;    // phi
    Matrix<DDRMat> mMyParamPoint;

    bool mMyIntersectionFlag;

    moris::Cell< Matrix< DDRMat > >   mMyIntersectionPoints;
    moris::Cell< Matrix< DDRMat > >   mMyIntersFieldSensVal;    // dxgamma/dphi
    moris::Cell< Matrix< DDRMat > >   mMyIntersSensVal;
    moris::Cell< Matrix< IndexMat > > mMyIntersSensInd;
//------------------------------------------------------------------------------
protected:

};

}// end ge namespace
}// end moris namespace
#endif /* PROJECTS_GEN_SRC_CL_GE_INTERSECTION_OBJECT_LINE_HPP_ */
