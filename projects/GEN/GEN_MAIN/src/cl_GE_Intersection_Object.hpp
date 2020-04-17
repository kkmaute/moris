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
class Geometry_Analytic;
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

        virtual moris_index compute_intersection(  )
        {
            MORIS_ASSERT( false, "ge::Intersection_Object::compute_intersection() - not implemented " );
            return 0;
        }
        virtual void compute_intersection_sensitivity( moris_index aMyIndex )
        {
            MORIS_ASSERT( false, "ge::Intersection_Object::compute_intersection_sensitivity() - not implemented " );
        }
        //------------------------------------------------------------------------------
        /*
         * @brief set the global coordinates and the time coordinates for the interpolators
         */
        virtual void set_coords_and_param_point( std::shared_ptr< Geometry_Analytic > & aGeomPointer,
                                                 Matrix<DDRMat> const & aGlobCoords,
                                                 Matrix<DDRMat> const & aTimeCoords,
                                                 Matrix<DDRMat> const & aFieldVals,     // LS field vals need to come from the geometry engine
                                                 Matrix<DDRMat> const & aParamPoint = Matrix< DDRMat >({{-1},{0}}) )  // start at the beginning of the line(local coordinate xi=-1) and at t=0        {
        {
            MORIS_ASSERT( false, "ge::Intersection_Object::set_coords_and_param_point() - not implemented" );
        }
        //------------------------------------------------------------------------------
        virtual void set_intersection( Matrix< DDRMat > aPoint )
        {
            MORIS_ASSERT( false, "ge::Intersection_Object::set_intersection() - not implemented " );
        }
        //******************************* get functions ********************************
        //------------------------------------------------------------------------------
        virtual mtk::Geometry_Type get_my_geom_type()
        {
            MORIS_ASSERT( false, "ge::Intersection_Object::get_my_geom_type() - not implemented " );
            return mtk::Geometry_Type::LINE;
        }
        //------------------------------------------------------------------------------
        virtual fem::Geometry_Interpolator* get_my_geom_interp()
        {
            MORIS_ASSERT( false, "ge::Intersection_Object::get_my_geom_interp() - not implemented " );
            return mDummyGeomInterp;
        }
        //------------------------------------------------------------------------------
        virtual fem::Field_Interpolator* get_my_field_interp()
        {
            MORIS_ASSERT( false, "ge::Intersection_Object::get_my_field_interp() - not implemented " );
            return mDummmyFieldInterp;
        }
        //------------------------------------------------------------------------------
        virtual Matrix<DDRMat> get_my_param_point()
        {
            MORIS_ASSERT( false, "ge::Intersection_Object::get_my_param_point() - not implemented " );
            return mDummyMat;
        }
        //------------------------------------------------------------------------------
        virtual Matrix<DDRMat> get_my_global_coord()
        {
            MORIS_ASSERT( false, "ge::Intersection_Object::get_my_global_coord() - not implemented " );
            return mDummyMat;
        }
        //------------------------------------------------------------------------------
        virtual Matrix<DDRMat> get_my_time_coord()
        {
            MORIS_ASSERT( false, "ge::Intersection_Object::get_my_time_coord() - not implemented " );
            return mDummyMat;
        }
        //------------------------------------------------------------------------------
        virtual Matrix<DDRMat> get_my_field_vals()
        {
            MORIS_ASSERT( false, "ge::Intersection_Object::get_my_field_vals() - not implemented " );
            return mDummyMat;
        }
        //------------------------------------------------------------------------------
        virtual bool get_intersection_flag()
        {
            MORIS_ASSERT( false, "ge::Intersection_Object::get_intersection_flag() - not implemented " );
            return false;
        }
        //------------------------------------------------------------------------------
        virtual Matrix< DDRMat > get_intersection_point_local_coord( moris_index aMyIndex )
        {
            return mDummyMat;
        }
        //------------------------------------------------------------------------------
        virtual Matrix< DDRMat > get_intersection_point_global_coord( moris_index aMyIndex )
        {
            MORIS_ASSERT( false, "ge::Intersection_Object::get_intersection_global_coord() - not implemented " );
            return mDummyMat;
        }
        //------------------------------------------------------------------------------
        virtual uint get_num_intersection_point( )
        {
            MORIS_ASSERT( false, "ge::Intersection_Object::get_num_intersection_point() - not implemented " );
            return 0;
        }
        //------------------------------------------------------------------------------
        virtual Matrix< DDRMat > get_field_sensitivity_vals( moris_index aMyIndex )
        {
            MORIS_ASSERT( false, "ge::Intersection_Object::get_field_sensitivity_vals() - not implemented " );
            return mDummyMat;
        }

        //***************** get functions for asserts/debugs ***************************
        //------------------------------------------------------------------------------
        virtual moris::Cell< Matrix< DDRMat > > get_field_sens_vals_cell()
        {
            MORIS_ASSERT( false, "ge::Intersection_Object::get_field_sens_vals_cell() - not implemented " );
            return mDummyCell;
        }
//------------------------------------------------------------------------------
    private:
        // dummy member variables
        fem::Geometry_Interpolator* mDummyGeomInterp  = nullptr;
        fem::Field_Interpolator*    mDummmyFieldInterp = nullptr;

        Matrix<DDRMat> mDummyMat;
        Cell< Matrix<DDRMat> > mDummyCell;
//------------------------------------------------------------------------------
    protected:

    };
} /* namespace gen */
} /* namespace moris */



#endif /* PROJECTS_GEN_SRC_CL_GE_INTERSECTION_OBJECT_HPP_ */
