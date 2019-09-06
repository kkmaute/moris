/*
 * cl_GE_Analytic.hpp
 *
 *  Created on: Jan 4, 2019
 *      Author: sonne
 */

#ifndef PROJECTS_GEN_SRC_CL_GE_ANALYTIC_HPP_
#define PROJECTS_GEN_SRC_CL_GE_ANALYTIC_HPP_

// GE includes
#include "cl_GE_Geometry.hpp"

// other includes
#include "assert.hpp"
//------------------------------------------------------------------------------

namespace moris{
    namespace ge{

		class Analytic : public Geometry
		{
		public:
		    /*
		     * function needs to be set with set_analytical_function because class is constructed through the factory
		     */
			Analytic()
		    {
	            mMySpaceInterpType  = fem::Interpolation_Type::LAGRANGE;
	            mMySpaceInterpOrder = mtk::Interpolation_Order::LINEAR;
	            mMyTimeInterpType   = fem::Interpolation_Type::CONSTANT;
	            mMyTimeInterpOrder  = mtk::Interpolation_Order::CONSTANT;
		    };

			~Analytic()
			{
			};
			//------------------------------------------------------------------------------
			/*
			 * @brief function to report if geometry representation type
			 */
		    enum GeomType
		    get_geom_type() const
		    {
		        return GeomType::ANALYTIC;
		    }
		    //------------------------------------------------------------------------------
		    bool
		    check_if_function_is_set()
		    {
//		        MORIS_ASSERT(mFuncAnalytic != nullptr,"ge::GE_Analytic::check_if_function_is_set(): analytic function not set ");
		        bool tBool;
		        if ( mFuncAnalytic != nullptr )
		        {
		            tBool = true;
		        }
		        else
		        {
		            tBool = false;
		        }
		        return tBool;
		    }
		    //------------------------------------------------------------------------------
		    bool
		    check_if_sensitivity_function_is_set()
		    {
//		        MORIS_ASSERT(mFuncAnalyticDphiDx != nullptr,"ge::GE_Analytic::check_if_sensitivity_function_is_set(): analytic sensitivity function not set ");
                bool tBool;
                if ( mFuncAnalyticDphiDx != nullptr )
                {
                    tBool = true;
                }
                else
                {
                    tBool = false;
                }
                return tBool;
		    }
		    //------------------------------------------------------------------------------
		    void
		    set_my_constants( moris::Cell< real > aMyConstants )
		    {
		        mMyConstants = aMyConstants;
		    }
		    //------------------------------------------------------------------------------
	        void
	        set_my_mesh(mtk::Mesh_Manager* aMyMesh)
	        {
	            mMyMesh = aMyMesh;
	        }
	        //------------------------------------------------------------------------------
	        void
	        set_my_interpolation_rules( fem::Interpolation_Type  aSpaceInterpType,
	                                    mtk::Interpolation_Order aSpaceInterpOrder,
	                                    fem::Interpolation_Type  aTimeInterpType  = fem::Interpolation_Type::CONSTANT,
	                                    mtk::Interpolation_Order aTimeInterpOrder = mtk::Interpolation_Order::CONSTANT )
	        {
	            mMySpaceInterpType  = aSpaceInterpType;
	            mMySpaceInterpOrder = aSpaceInterpOrder;
	            mMyTimeInterpType   = aTimeInterpType;
	            mMyTimeInterpOrder  = aTimeInterpOrder;
	        }
	        //------------------------------------------------------------------------------
	        fem::Interpolation_Type
	        get_my_space_interpolation_type()
	        {
	            return mMySpaceInterpType;
	        }
            //------------------------------------------------------------------------------
            mtk::Interpolation_Order
            get_my_space_interpolation_order()
            {
                return mMySpaceInterpOrder;
            }
            //------------------------------------------------------------------------------
            fem::Interpolation_Type
            get_my_time_interpolation_type()
            {
                return mMyTimeInterpType;
            }
            //------------------------------------------------------------------------------
            mtk::Interpolation_Order
            get_my_time_interpolation_order()
            {
                return mMyTimeInterpOrder;
            }
		    //------------------------------------------------------------------------------
	        mtk::Mesh_Manager*
	        get_my_mesh()
	        {
	            MORIS_ASSERT( mMyMesh != nullptr, "ge::Geometry::get_my_mesh(): the associated mesh has not been set" );
	            return mMyMesh;
	        }
			/*
			 * *****************************************************************************
			 * pass in a user-defined function
			 * *****************************************************************************
			 */
//------------------------------------------------------------------------------
			/*
			 * @brief sets the analytic function \phi for the analytic geometry class
			 *
			 * @param[in] *funcPointer - pointer to the user defined function
			 *
			 */
            void
            set_analytical_function( real ( *funcPointer )( const Matrix< DDRMat > & aCoordinate, Cell< real > aConst ) )
            {
                mFuncAnalytic = funcPointer;
            };

            //------------------------------------------------------------------------------
            /*
             * @brief sets the analytic function \phi for the analytic geometry class
             * note: this type is used when the function needs multiple cells as inputs
             * ( currently only multi_cylinder_function() requires this )
             *
             * @param[in] *funcPointer - pointer to the user defined function
             *
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
            /*
             * @brief sets the sensitivity function d\phi/dx for the analytic geometry class
             *
             * @param[in] *funcPointer - pointer to the user defined sensitivity function
             *
             */
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
//------------------------------------------------------------------------------
            /*
             * @brief sets the function \phi for the analytic geometry class
             *
             * @param[in] aGeomType - enum which names a function from the list of known standardized functions
             *
             */
            void
            set_analytical_function( AnalyticType aGeomType )
            {
                switch(aGeomType)
                {
                case(AnalyticType::CIRCLE):
                {
                    MORIS_ASSERT( mMyConstants.size() == 3, "Analytic::set_analytical_function(): incorrect size for constants list (needs to be 3 constants for the circle function: 1) x_c, 2) y_c, 3) r " );
                    MORIS_ASSERT( mMyConstants(2) > 0, "Geometry_Library::circle_function(): radius must be > 0" );
                    this->set_analytical_function( circle_function );
                    break;
                }
                case( AnalyticType::COMPOSITE_FIBER ):
                {
                    MORIS_ASSERT( mMyConstants.size() == 10, "Analytic::set_analytical_function(): incorrect size for constants list " );
                    MORIS_ASSERT( mMyConstants(0) > 0, "Analytic::set_analytical_function(): fiber radius must be > 0 " );
                    this->set_analytical_function( composite_fiber_function );
                    break;
                }
                case( AnalyticType::COMPOSITE_FIBER_STRAIGHT_1 ):
                {
                    MORIS_ASSERT( mMyConstants.size() == 10, "Analytic::set_analytical_function(): incorrect size for constants list " );
                    MORIS_ASSERT( mMyConstants(0) > 0, "Analytic::set_analytical_function(): fiber radius must be > 0 " );
                    this->set_analytical_function( composite_fiber_straight_1_function );
                    break;
                }
                case( AnalyticType::COMPOSITE_FIBER_STRAIGHT_2 ):
                {
                    MORIS_ASSERT( mMyConstants.size() == 10, "Analytic::set_analytical_function(): incorrect size for constants list " );
                    MORIS_ASSERT( mMyConstants(0) > 0, "Analytic::set_analytical_function(): fiber radius must be > 0 " );
                    this->set_analytical_function( composite_fiber_straight_2_function );
                    break;
                }
                case( AnalyticType::COMPOSITE_FIBER_STRAIGHT_3 ):
                {
                    MORIS_ASSERT( mMyConstants.size() == 10, "Analytic::set_analytical_function(): incorrect size for constants list " );
                    MORIS_ASSERT( mMyConstants(0) > 0, "Analytic::set_analytical_function(): fiber radius must be > 0 " );
                    this->set_analytical_function( composite_fiber_straight_3_function );
                    break;
                }
                case( AnalyticType::COMPOSITE_FIBER_WAVE_1 ):
                {
                    MORIS_ASSERT( mMyConstants.size() == 10, "Analytic::set_analytical_function(): incorrect size for constants list " );
                    MORIS_ASSERT( mMyConstants(0) > 0, "Analytic::set_analytical_function(): fiber radius must be > 0 " );
                    this->set_analytical_function( composite_fiber_wave_1_function );
                    break;
                }
                case( AnalyticType::COMPOSITE_FIBER_WAVE_2 ):
                {
                    MORIS_ASSERT( mMyConstants.size() == 10, "Analytic::set_analytical_function(): incorrect size for constants list " );
                    MORIS_ASSERT( mMyConstants(0) > 0, "Analytic::set_analytical_function(): fiber radius must be > 0 " );
                    this->set_analytical_function( composite_fiber_wave_2_function );
                    break;
                }
                case( AnalyticType::COMPOSITE_FIBER_WAVE_3 ):
                {
                    MORIS_ASSERT( mMyConstants.size() == 10, "Analytic::set_analytical_function(): incorrect size for constants list " );
                    MORIS_ASSERT( mMyConstants(0) > 0, "Analytic::set_analytical_function(): fiber radius must be > 0 " );
                    this->set_analytical_function( composite_fiber_wave_3_function );
                    break;
                }
                case( AnalyticType::GYROID ):
                {
                    MORIS_ASSERT( mMyConstants.size() == 0, "Analytic::set_analytical_function(): gyroid function does not require a list of constants " );
                    this->set_analytical_function( gyroid_function );
                    break;
                }
                case( AnalyticType::MULTI_CYLINDER ):
                {
                        this->set_analytical_function( multi_cylinder_function );
                        break;
                }
                case( AnalyticType::PLANE ):
                {
                    MORIS_ASSERT( mMyConstants.size() == 6, "Analytic::set_analytical_function(): incorrect size for constants list " );
                        this->set_analytical_function( plane_function );
                        break;
                }
                case( AnalyticType::SPHERE ):
                {
                    MORIS_ASSERT( mMyConstants.size() == 4, "Analytic::set_analytical_function(): incorrect size for constants list (needs to be 4 constants for the sphere function: 1) x_c, 2) y_c, 3) z_c, 4) r " );
                    MORIS_ASSERT( mMyConstants(3) > 0, "Geometry_Library::sphere_function(): radius must be > 0" );
                    this->set_analytical_function( sphere_function );
                    break;
                }
                case( AnalyticType::SPIRAL ):
                {
                    MORIS_ASSERT( mMyConstants.size() == 10, "Analytic::set_analytical_function(): incorrect size for constants list " );
                    MORIS_ASSERT( mMyConstants(1) > 0, "Analytic::set_analytical_function(): fiber radius must be > 0 " );
                    this->set_analytical_function( spiral_function );
                    break;
                }
                default:
                {
                        MORIS_ERROR(false, "cl_GE_Analytical::set_analytical_function() please choose from standardized functions or pass in your own function");
                        break;
                }

                }
            };

            //------------------------------------------------------------------------------
            /*
             * @brief sets the sensitivity function d\phi/dx for the analytic geometry class
             *
             * @param[in] aGeomType - enum which names a sensitivity function from the list of known standardized functions
             *
             */
            void
            set_analytical_function_dphi_dx( AnalyticType aGeomType )
            {
                switch(aGeomType)
                {
                case( AnalyticType::CIRCLE ):
                {
                    MORIS_ASSERT( mMyConstants.size() == 3, "Analytic::set_analytical_function(): incorrect size for constants list, be sure to set the constants before setting the function itself (needs to be three constants for the circle function: 1) x_c, 2) y_c, 3) r " );
                    MORIS_ASSERT( mMyConstants(2) > 0, "Geometry_Library::circle_function(): radius must be > 0" );
                    this->set_analytical_function_dphi_dx( circle_function_dphi_dx );
                    break;
                }
                case( AnalyticType::COMPOSITE_FIBER ):
                {
                    MORIS_ERROR(false, "dphi/dx currently not set for composite_fiber_function");
                    break;
                }
                case( AnalyticType::COMPOSITE_FIBER_STRAIGHT_1 ):
                {
                    MORIS_ERROR(false, "dphi/dx currently not set for composite_fiber_straight_1_function");
                    break;
                }
                case( AnalyticType::COMPOSITE_FIBER_STRAIGHT_2 ):
                {
                    MORIS_ERROR(false, "dphi/dx currently not set for composite_fiber_straight_2_function");
                    break;
                }
                case( AnalyticType::COMPOSITE_FIBER_STRAIGHT_3 ):
                {
                    MORIS_ERROR(false, "dphi/dx currently not set for composite_fiber_straight_3_function");
                    break;
                }
                case( AnalyticType::COMPOSITE_FIBER_WAVE_1 ):
                {
                        MORIS_ERROR(false, "dphi/dx currently not set for composite_fiber_wave_1_function");
                        break;
                }
                case( AnalyticType::COMPOSITE_FIBER_WAVE_2 ):
                {
                        MORIS_ERROR(false, "dphi/dx currently not set for composite_fiber_wave_2_function");
                        break;
                }
                case( AnalyticType::COMPOSITE_FIBER_WAVE_3 ):
                {
                        MORIS_ERROR(false, "dphi/dx currently not set for composite_fiber_wave_3_function");
                        break;
                }
                case( AnalyticType::GYROID ):
                {
                        MORIS_ERROR(false, "dphi/dx currently not set for gyroid_function");
                        break;
                }
                case( AnalyticType::MULTI_CYLINDER ):
                {
                        MORIS_ERROR(false, "dphi/dx currently not set for mutli_cylinder_function");
                        break;
                }
                case( AnalyticType::PLANE ):
                {
                        MORIS_ERROR(false, "dphi/dx currently not set for plane_function");
                        break;
                }
                case( AnalyticType::SPIRAL ):
                {
                        MORIS_ERROR(false, "dphi/dx currently not set for spiral_function");
                        break;
                }
                case( AnalyticType::SPHERE ):

                {
                    MORIS_ASSERT( mMyConstants.size() == 4, "Analytic::set_analytical_function(): incorrect size for constants list (needs to be 4 constants for the circle function: 1) x_c, 2) y_c, 3) z_c, 4) r " );
                    MORIS_ASSERT( mMyConstants(3) > 0, "Geometry_Library::circle_function(): radius must be > 0" );
                    this->set_analytical_function_dphi_dx( sphere_function_dphi_dx );
                    break;
                }
                default:
                {
                        MORIS_ERROR(false, "cl_GE_Analytical::set_analytical_function() please choose from standardized functions or pass in your own function");
                        break;
                }

                }
            };

            //------------------------------------------------------------------------------
            /*
             * @brief determines the value of the function at a specified coordinate
             *
             * @param[in] aPoint - point for value to be determined at
             * @oaram[in] aConst - cell of real values which describe functional constants (such as the radius of a sphere)
             *
             * @param[out] real - function value at specified coordinate
             */
			real
			get_field_val_at_coordinate( const Matrix< DDRMat >  & aPoint )
			{
				return mFuncAnalytic( aPoint, mMyConstants );
			};

            //------------------------------------------------------------------------------
            /*
             * @brief determines the value of the function at a specified coordinate;
             * used when the function requires multiple cell inputs ( e.g. multi_cylinder_function )
             *
             * @param[in] aPoint - point for value to be determined at
             * @oaram[in] aConst - cell of real values which describe functional constants (such as the radius of a sphere)
             *
             * @param[out] real - function value at specified coordinate
             *
             */
			//fixme update/add default argument values if needed
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
            /*
             * @brief determines the sensitivity value of the function at a specified coordinate
             *
             * @param[in] aPoint - point for value to be determined at
             * @oaram[in] aConst - cell of real values which describe functional constants (such as the radius of a sphere)
             *
             * @param[out] real - function sensitivty value at specified coordinate
             */
			Matrix< DDRMat >
			get_sensitivity_dphi_dp_at_coordinate( const Matrix< DDRMat >  & aPoint )
            {
			    return mFuncAnalyticDphiDx( aPoint, mMyConstants );
            };


			//------------------------------------------------------------------------------

        private:

			moris::Cell< real > mMyConstants;
			mtk::Mesh_Manager*  mMyMesh = nullptr;

            real ( *mFuncAnalytic )( const Matrix< DDRMat > & aPoint, moris::Cell< real> aConstant ) = nullptr;
            Matrix< DDRMat > ( *mFuncAnalyticDphiDx)( const Matrix< DDRMat > & aPoint, moris::Cell< real > aConstant ) = nullptr;
            /*
             * this is needed for tricky functions which require multiple cell inputs (e.g. multi_cylinder_function)
             */
            real ( *mFuncAnalyticExtra )( const Matrix< DDRMat >        & aCoordinates,
                                                Cell<Cell<moris::real>> & aCenter,
                                                Cell<moris::real>       & aRadius,
                                                Cell<moris::real>       & aLength,
                                                Cell<Cell<moris::real>> & aAxis ) = nullptr;

            fem::Interpolation_Type  mMySpaceInterpType;
            mtk::Interpolation_Order mMySpaceInterpOrder;
            fem::Interpolation_Type  mMyTimeInterpType;
            mtk::Interpolation_Order mMyTimeInterpOrder;
            //------------------------------------------------------------------------------

        protected:

		};
	} /* namespace gen */
} /* namespace moris */



#endif /* PROJECTS_GEN_SRC_CL_GE_ANALYTIC_HPP_ */
