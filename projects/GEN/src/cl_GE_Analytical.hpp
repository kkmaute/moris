/*
 * cl_GE_Analytical.hpp
 *
 *  Created on: Jan 4, 2019
 *      Author: sonne
 */

#ifndef PROJECTS_GEN_SRC_CL_GE_ANALYTICAL_HPP_
#define PROJECTS_GEN_SRC_CL_GE_ANALYTICAL_HPP_

// GE includes
#include "cl_GE_Geometry.hpp"

// other includes
#include "assert.hpp"
//------------------------------------------------------------------------------

namespace moris{
    namespace ge{

		class Analytical : public Geometry
		{
		public:
		    /*
		     * function needs to be set with set_analytical_function because class is constructed through the factory
		     */

			Analytical(){};

			~Analytical(){};

			/*
			 * *****************************************************************************
			 * pass in a user-defined function
			 * *****************************************************************************
			 */

            void
            set_analytical_function( real ( *funcPointer )( const Matrix< DDRMat > & aCoordinate, Cell< real > aConst ) )
            {
                mFuncAnalytic = funcPointer;
            };

            //------------------------------------------------------------------------------
            /* used when the function needs multiple cells as inputs
             * ( currently only multi_cylinder_function() requires this )
             */
            void
            set_analytical_function( real ( *funcPointer )( const Matrix< DDRMat >        & aCoordinates,
                                                                  Cell<Cell<moris::real>> & aCenter,
                                                                  Cell<moris::real>       & aRadius,
                                                                  Cell<moris::real>       & aLength,
                                                                  Cell<Cell<moris::real>> & aAxis ) )
            {
                mFuncAnalyticExtra = funcPointer;
            };

            //------------------------------------------------------------------------------

            void
            set_analytical_function_dphi_dx( Matrix< DDRMat > ( *funcPointer )( const Matrix< DDRMat > & aCoordinate, Cell< real > aConst ) )
            {
                mFuncAnalyticDphiDx = funcPointer;
            };

			//------------------------------------------------------------------------------
            /*
             * *****************************************************************************
             * use function from list of standardized functions
             * *****************************************************************************
             */
            void
            set_analytical_function( type aGeomType )
            {
                switch(aGeomType)
                {
                case( type::SPHERE ):
                        this->set_analytical_function( sphere_function );
                        break;
                case( type::PLANE ):
                        this->set_analytical_function( plane_function );
                        break;
                case( type::SPIRAL ):
                        this->set_analytical_function( spiral_function );
                        break;
                case( type::GYROID ):
                        this->set_analytical_function( gyroid_function );
                        break;
                case( type::COMPOSITE_FIBER ):
                        this->set_analytical_function( composite_fiber_function );
                        break;
                case( type::COMPOSITE_FIBER_WAVE_1 ):
                        this->set_analytical_function( composite_fiber_wave_1_function );
                        break;
                case( type::COMPOSITE_FIBER_WAVE_2 ):
                        this->set_analytical_function( composite_fiber_wave_2_function );
                        break;
                case( type::COMPOSITE_FIBER_WAVE_3 ):
                        this->set_analytical_function( composite_fiber_wave_3_function );
                        break;
                case( type::COMPOSITE_FIBER_STRAIGHT_1 ):
                        this->set_analytical_function( composite_fiber_straight_1_function );
                        break;
                case( type::COMPOSITE_FIBER_STRAIGHT_2 ):
                        this->set_analytical_function( composite_fiber_straight_2_function );
                        break;
                case( type::COMPOSITE_FIBER_STRAIGHT_3 ):
                        this->set_analytical_function( composite_fiber_straight_3_function );
                        break;
                case( type::MULTI_CYLINDER ):
                        this->set_analytical_function( multi_cylinder_function );
                        break;
                default:
                        MORIS_ERROR(false, "cl_GE_Analytical::set_analytical_function() please choose from standardized functions or pass in your own function");
                        break;
                }

            };

            //------------------------------------------------------------------------------
            void
            set_analytical_function_dphi_dx( type aGeomType )
            {
                switch(aGeomType)
                {
                case(type::SPHERE):
                        this->set_analytical_function_dphi_dx( sphere_function_dphi_dx );
                        break;
                case(type::PLANE):
                        MORIS_ERROR(false, "dphi/dx currently not set for plane_function");
                        break;
                case( type::SPIRAL ):
                        MORIS_ERROR(false, "dphi/dx currently not set for spiral_function");
                        break;
                case( type::GYROID ):
                        MORIS_ERROR(false, "dphi/dx currently not set for gyroid_function");
                        break;
                case( type::COMPOSITE_FIBER ):
                        MORIS_ERROR(false, "dphi/dx currently not set for composite_fiber_function");
                        break;
                case( type::COMPOSITE_FIBER_WAVE_1 ):
                        MORIS_ERROR(false, "dphi/dx currently not set for composite_fiber_wave_1_function");
                        break;
                case( type::COMPOSITE_FIBER_WAVE_2 ):
                        MORIS_ERROR(false, "dphi/dx currently not set for composite_fiber_wave_2_function");
                        break;
                case( type::COMPOSITE_FIBER_WAVE_3 ):
                        MORIS_ERROR(false, "dphi/dx currently not set for composite_fiber_wave_3_function");
                        break;
                case( type::COMPOSITE_FIBER_STRAIGHT_1 ):
                        MORIS_ERROR(false, "dphi/dx currently not set for composite_fiber_straight_1_function");
                        break;
                case( type::COMPOSITE_FIBER_STRAIGHT_2 ):
                        MORIS_ERROR(false, "dphi/dx currently not set for composite_fiber_straight_2_function");
                        break;
                case( type::COMPOSITE_FIBER_STRAIGHT_3 ):
                        MORIS_ERROR(false, "dphi/dx currently not set for composite_fiber_straight_3_function");
                        break;
                case( type::MULTI_CYLINDER ):
                        MORIS_ERROR(false, "dphi/dx currently not set for mutli_cylinder_function");
                        break;
                default:
                        MORIS_ERROR(false, "cl_GE_Analytical::set_analytical_function() please choose from standardized functions or pass in your own function");
                        break;
                }

            };

            //------------------------------------------------------------------------------

			real
			get_field_val_at_coordinate( const Matrix< DDRMat >  & aPoint,
			                                   moris::Cell< real > aConst )
			{
				return mFuncAnalytic( aPoint, aConst );
			};

            //------------------------------------------------------------------------------
			/* used when the function does not require any constants as inputs ( e.g. gyroid_function() ) */
            real
            get_field_val_at_coordinate( const Matrix< DDRMat >  & aPoint )
            {
                moris::Cell< real > aConst(1);
                aConst(1) = 0.0;
                return mFuncAnalytic( aPoint, aConst );
            };

            //------------------------------------------------------------------------------
            /* used when the function requires multiple cell inputs ( e.g. multi_cylinder_function )
             * ---update/add defaults if needed---
             */
            real
            get_field_val_at_coordinate( const Matrix< DDRMat >        & aCoordinates,
                                               Cell<Cell<moris::real>> & aCenter,
                                               Cell<moris::real>       & aRadius,
                                               Cell<moris::real>       & aLength,
                                               Cell<Cell<moris::real>> & aAxis )
            {
                return mFuncAnalyticExtra( aCoordinates, aCenter, aRadius, aLength, aAxis );
            };
            //------------------------------------------------------------------------------

			Matrix< DDRMat >
			get_sensitivity_dphi_dp_at_coordinate( const Matrix< DDRMat >  & aPoint,
			                                             moris::Cell< real > aConst)
            {
			    return mFuncAnalyticDphiDx( aPoint, aConst );
            };

			//------------------------------------------------------------------------------


        private:

            real ( *mFuncAnalytic )( const Matrix< DDRMat > & aPoint, moris::Cell< real> aConstant );

            Matrix< DDRMat > ( *mFuncAnalyticDphiDx)( const Matrix< DDRMat > & aPoint, moris::Cell< real > aConstant );

            /*
             * this is needed for tricky functions which require multiple cell inputs (e.g. multi_cylinder_function)
             * ---update/add default values if needed---
             */
            real ( *mFuncAnalyticExtra )( const Matrix< DDRMat >        & aCoordinates,
                                                Cell<Cell<moris::real>> & aCenter,
                                                Cell<moris::real>       & aRadius,
                                                Cell<moris::real>       & aLength,
                                                Cell<Cell<moris::real>> & aAxis );
            //------------------------------------------------------------------------------


        protected:

		};
	} /* namespace gen */
} /* namespace moris */



#endif /* PROJECTS_GEN_SRC_CL_GE_ANALYTICAL_HPP_ */
