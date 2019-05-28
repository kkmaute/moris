/*
 * cl_GE_Intersection_Object.hpp
 *
 *  Created on: May 22, 2019
 *      Author: sonne
 */

#ifndef PROJECTS_GEN_SRC_CL_GE_INTERSECTION_OBJECT_HPP_
#define PROJECTS_GEN_SRC_CL_GE_INTERSECTION_OBJECT_HPP_

// GE includes
#include "cl_GE_Enums.hpp"

// FEM includes
#include "cl_FEM_Enums.hpp"
#include "cl_FEM_Interpolation_Rule.hpp"
#include "cl_FEM_Geometry_Interpolator.hpp"
#include "cl_FEM_Field_Interpolator.hpp"

// MTK includes
#include "cl_MTK_Enums.hpp"

namespace moris
{
namespace ge
{
/*
 * @brief Define a simple geometric primitive, such as a line.
 *        A primitive can be passed into the geometry engine to ask if there is an intersection with it and
 *        a known geometric representation.
 */
class Intersection_Object
    {
    public:

        Intersection_Object()
        {
        };

        ~Intersection_Object(){};

        //------------------------------------------------------------------------------
        /*
         * @brief set the global coordinates and the time coordinates for the interpolators
         */
        virtual void set_coords_and_param_point( std::shared_ptr< Geometry > & aGeomPointer,
                                                 Matrix<DDRMat> const & aGlobCoords,
                                                 Matrix<DDRMat> const & aTimeCoords,
                                                 Matrix<DDRMat> const & aFieldVals,     // LS field vals need to come from the geometry engine
                                                 Matrix<DDRMat> const & aParamPoint = Matrix< DDRMat >({{-1},{0}}) )  // start at the beginning of the line(local coordinate xi=-1) and at t=0        {
        {
            MORIS_ERROR( false, "ge::Intersection_Object::set_coords_and_param_point(): not implemented" );
        }
        //******************************* get functions ********************************
        //------------------------------------------------------------------------------
        virtual mtk::Geometry_Type
        get_my_geom_type()
        {
            MORIS_ERROR( false, "ge::Intersection_Object::get_my_geom_type(): not implemented " );
            return mtk::Geometry_Type::LINE;
        }
        //------------------------------------------------------------------------------
        virtual fem::Geometry_Interpolator*
        get_my_geom_interp()
        {
            MORIS_ERROR( false, "ge::Intersection_Object::get_my_geom_interp(): not implemented " );
            return mDummyGeomInterp;
        }
        //------------------------------------------------------------------------------
        virtual fem::Field_Interpolator*
        get_my_field_interp()
        {
            MORIS_ERROR( false, "ge::Intersection_Object::get_my_field_interp(): not implemented " );
            return mDummmyFieldInterp;
        }
        //------------------------------------------------------------------------------
        virtual Matrix<DDRMat>
        get_my_param_point()
        {
            MORIS_ERROR( false, "ge::Intersection_Object::get_my_param_point(): not implemented " );
            return mDummyParamPoint;
        }
        //------------------------------------------------------------------------------
        virtual Matrix<DDRMat>
        get_my_global_coord()
        {
            MORIS_ERROR( false, "ge::Intersection_Object::get_my_global_coord(): not implemented " );
            return mDummyGlobalCoords;
        }
        //------------------------------------------------------------------------------
        virtual Matrix<DDRMat>
        get_my_time_coord()
        {
            MORIS_ERROR( false, "ge::Intersection_Object::get_my_time_coord(): not implemented " );
            return mDummyTimeCoords;
        }
        //------------------------------------------------------------------------------
        virtual Matrix<DDRMat>
        get_my_field_vals()
        {
            MORIS_ERROR( false, "ge::Intersection_Object::get_my_field_vals(): not implemented " );
            return mDummyFieldVals;
        }
//------------------------------------------------------------------------------
    private:
        // dummy member variables
        fem::Geometry_Interpolator* mDummyGeomInterp  = nullptr;
        fem::Field_Interpolator*    mDummmyFieldInterp = nullptr;

        Matrix<DDRMat> mDummyGlobalCoords;
        Matrix<DDRMat> mDummyTimeCoords;
        Matrix<DDRMat> mDummyFieldVals;
        Matrix<DDRMat> mDummyParamPoint;
//------------------------------------------------------------------------------
    protected:

    };
} /* namespace gen */
} /* namespace moris */



#endif /* PROJECTS_GEN_SRC_CL_GE_INTERSECTION_OBJECT_HPP_ */
