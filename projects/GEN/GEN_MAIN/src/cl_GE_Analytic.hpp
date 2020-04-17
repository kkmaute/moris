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
#include "cl_GE_Property.hpp"

#include <functional>

// other includes
#include "assert.hpp"
//------------------------------------------------------------------------------

namespace moris{
    namespace ge{

//    //------------------------------------------------------------------------------
//    typedef std::function< moris::real ( const Matrix< DDRMat >    & aPoint,
//                                         const moris::Cell< real >   aConstant ) > AnalyticFunction;
//    //------------------------------------------------------------------------------
////    typedef std::function< moris::real ( const Matrix< DDRMat >         & aCoordinates,
////                                               Cell<Cell<moris::real>>  & aCenter,
////                                               Cell<moris::real>        & aRadius,
////                                               Cell<moris::real>        & aLength,
////                                               Cell<Cell<moris::real>>  & aAxis ) > AnalyticFunction;       // special case of the analytic function
//    //------------------------------------------------------------------------------
//    typedef std::function< Matrix< DDRMat > ( const Matrix< DDRMat >    & aPoint,
//                                              const moris::Cell< real >   aConstant ) > AnalyticFunctionDphiDp;
    //------------------------------------------------------------------------------
		class Analytic : public Geometry_Analytic
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
		    check_if_function_is_set( const moris_index aSubIndex )
		    {
		        bool tBool;
		        if ( mPropertyList.size()-1 == (uint)aSubIndex )
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
		    check_if_sensitivity_function_is_set( const moris_index aSubIndex )
		    {
                bool tBool;
                if ( mPropertyList(aSubIndex).is_derivative_set() )
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
		    set_constants( moris::Cell< real > aConstants,
		                          moris_index  aSubIndex )
		    {
		        mPropertyList( aSubIndex ).set_params( aConstants );
		    }
		    //------------------------------------------------------------------------------
	        void
	        set_my_mesh(std::shared_ptr< moris::hmr::Mesh > aMesh)
	        {
	            mMesh_HMR = aMesh;
	        }
            //------------------------------------------------------------------------------
            void
            set_my_mesh(mtk::Mesh_Manager* aMesh)
            {
                mMesh = aMesh;
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
	            MORIS_ASSERT( mMesh != nullptr, "ge::Geometry::get_my_mesh(): the associated mesh has not been set" );
	            return mMesh;
	        }

	        std::shared_ptr< moris::hmr::Mesh >
            get_my_mesh_HMR()
            {
                MORIS_ASSERT( mMesh_HMR != nullptr, "ge::Geometry::get_my_mesh(): the associated mesh has not been set" );
                return mMesh_HMR;
            }
			/*
			 * *****************************************************************************
			 * pass in a user-defined function
			 * *****************************************************************************
			 */
//------------------------------------------------------------------------------
            /*
             * @brief sets the analytic function \phi and the derivative(s) \partial \phi / \partial p
             *
             * @param[in] aFunc01 - pointer to the user defined function
             * @param[in] aFunc02 - pointer to the user defined derivative(s)
             * @param[in] aConstants - list of constants for the function to use
             *
             * @param[out] index to the sub-type
             *
             */
            moris_index
            set_analytical_function_and_dphi_dp( std::function< moris::real ( const Matrix< DDRMat >    & aPoint,
                                                                              const moris::Cell< real >   aConstant ) > aFunc01,
                                                 std::function< Matrix< DDRMat > ( const Matrix< DDRMat >    & aPoint,
                                                                                   const moris::Cell< real >   aConstant ) > aFunc02,
                                                 moris::Cell< real > aConstants )
            {
                ge::Property tNewProperty( aFunc01,
                                           aFunc02,
                                           aConstants );
                mPropertyList.push_back( tNewProperty );

                return mPropertyList.size()-1;
            };

			/*
			 * @brief sets the analytic function \phi for the analytic geometry class
			 *
			 * @param[in] *funcPointer - pointer to the user defined function
			 * @param[in] aConstants - list of constants for the function to use
			 *
			 * @param[out] index to the sub-type
			 *
			 */
            moris_index
            set_analytical_function( real ( *funcPointer )( const Matrix< DDRMat > &  aCoordinate,
                                                                  Cell< real >        aConst ),
                                                                  moris::Cell< real > aConstants = { { 0 } } )
            {
                ge::Property tNewProperty( funcPointer,
                                           aConstants );
                mPropertyList.push_back( tNewProperty );

                return mPropertyList.size()-1;
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
             * @brief sets the sensitivity function \partial \phi / \partial p for the analytic geometry class
             *
             * @param[in] *funcPointer - pointer to the user defined sensitivity function
             * @param[in] aFuncIndex   - index of Property to set \partial \phi / \partial p for
             *
             */
            void
            set_analytical_function_dphi_dp( Matrix< DDRMat > ( *funcPointer )( const Matrix< DDRMat > & aCoordinate, Cell< real > aConst ),
                                             moris_index aFuncIndex = 0 )
            {
                mPropertyList( aFuncIndex ).set_derivative_functions( funcPointer );
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
             * @param[in] aConstants - constants needed by the function
             *
             * @param[out] index of the function in the list of functions
             *
             */
            moris_index
            set_analytical_function( AnalyticType aGeomType, moris::Cell< real > aConstants )
            {
                moris_index tFuncIndex = 0;
                switch(aGeomType)
                {
                case(AnalyticType::CIRCLE):
                {
                    MORIS_ASSERT( aConstants.size() == 3, "Analytic::set_analytical_function(): incorrect size for constants list (needs to be 3 constants for the circle function: 1) r, 2) x_c, 3) y_c " );
                    MORIS_ASSERT( aConstants(0) > 0, "Geometry_Library::circle_function(): radius must be > 0" );

                    tFuncIndex = this->set_analytical_function( circle_function, aConstants );
                                 mPropertyList( tFuncIndex ).set_derivative_functions( circle_function_dphi_dp );
                    return tFuncIndex;
                    break;
                }
                case( AnalyticType::COMPOSITE_FIBER ):
                {
                    MORIS_ASSERT( aConstants.size() == 10, "Analytic::set_analytical_function() - incorrect size for constants list " );
                    MORIS_ASSERT( aConstants(0) > 0, "Analytic::set_analytical_function(): fiber radius must be > 0 " );
                    tFuncIndex = this->set_analytical_function( composite_fiber_function, aConstants );
                    return tFuncIndex;
                    break;
                }
                case( AnalyticType::COMPOSITE_FIBER_STRAIGHT_1 ):
                {
                    MORIS_ASSERT( aConstants.size() == 10, "Analytic::set_analytical_function() - incorrect size for constants list " );
                    MORIS_ASSERT( aConstants(0) > 0, "Analytic::set_analytical_function(): fiber radius must be > 0 " );
                    tFuncIndex = this->set_analytical_function( composite_fiber_straight_1_function, aConstants );
                    return tFuncIndex;
                    break;
                }
                case( AnalyticType::COMPOSITE_FIBER_STRAIGHT_2 ):
                {
                    MORIS_ASSERT( aConstants.size() == 10, "Analytic::set_analytical_function() - incorrect size for constants list " );
                    MORIS_ASSERT( aConstants(0) > 0, "Analytic::set_analytical_function(): fiber radius must be > 0 " );
                    tFuncIndex = this->set_analytical_function( composite_fiber_straight_2_function, aConstants );
                    return tFuncIndex;
                    break;
                }
                case( AnalyticType::COMPOSITE_FIBER_STRAIGHT_3 ):
                {
                    MORIS_ASSERT( aConstants.size() == 10, "Analytic::set_analytical_function() - incorrect size for constants list " );
                    MORIS_ASSERT( aConstants(0) > 0, "Analytic::set_analytical_function(): fiber radius must be > 0 " );
                    tFuncIndex = this->set_analytical_function( composite_fiber_straight_3_function, aConstants );
                    return tFuncIndex;
                    break;
                }
                case( AnalyticType::COMPOSITE_FIBER_WAVE_1 ):
                {
                    MORIS_ASSERT( aConstants.size() == 10, "Analytic::set_analytical_function() - incorrect size for constants list " );
                    MORIS_ASSERT( aConstants(0) > 0, "Analytic::set_analytical_function(): fiber radius must be > 0 " );
                    tFuncIndex = this->set_analytical_function( composite_fiber_wave_1_function, aConstants );
                    return tFuncIndex;
                    break;
                }
                case( AnalyticType::COMPOSITE_FIBER_WAVE_2 ):
                {
                    MORIS_ASSERT( aConstants.size() == 10, "Analytic::set_analytical_function() - incorrect size for constants list " );
                    MORIS_ASSERT( aConstants(0) > 0, "Analytic::set_analytical_function(): fiber radius must be > 0 " );
                    tFuncIndex = this->set_analytical_function( composite_fiber_wave_2_function, aConstants );
                    return tFuncIndex;
                    break;
                }
                case( AnalyticType::COMPOSITE_FIBER_WAVE_3 ):
                {
                    MORIS_ASSERT( aConstants.size() == 10, "Analytic::set_analytical_function() - incorrect size for constants list " );
                    MORIS_ASSERT( aConstants(0) > 0, "Analytic::set_analytical_function(): fiber radius must be > 0 " );
                    tFuncIndex = this->set_analytical_function( composite_fiber_wave_3_function, aConstants );
                    return tFuncIndex;
                    break;
                }
                case( AnalyticType::GYROID ):
                {
                    MORIS_ASSERT( aConstants.size() == 0, "Analytic::set_analytical_function() - gyroid function does not require a list of constants " );
                    tFuncIndex = this->set_analytical_function( gyroid_function, aConstants );
                    return tFuncIndex;
                    break;
                }
                case( AnalyticType::MULTI_CYLINDER ):
                {
                    // TODO make this guy play nice with the rest of the implementation
                    MORIS_ASSERT( false, "Analytic::set_analytical_function() -  MULTI_CYLINDER is not currently implemented " );
                    this->set_analytical_function( multi_cylinder_function );
                    return 0;
                    break;
                }
                case( AnalyticType::PLANE ):
                {
                    MORIS_ASSERT( aConstants.size() == 6, "Analytic::set_analytical_function() - incorrect size for constants list " );
                    tFuncIndex = this->set_analytical_function( plane_function, aConstants );
                    return tFuncIndex;
                    break;
                }
                case( AnalyticType::SPHERE ):
                {
                    MORIS_ASSERT( aConstants.size() == 4, "Analytic::set_analytical_function() - incorrect size for constants list (needs to be 4 constants for the sphere function: 1) r, 2) x_c, 3) y_c, 4) z_c " );
                    MORIS_ASSERT( aConstants(0) > 0, "Geometry_Library::sphere_function() - radius must be > 0" );

                    tFuncIndex = this->set_analytical_function( sphere_function, aConstants );
                                 mPropertyList( tFuncIndex ).set_derivative_functions( sphere_function_dphi_dp );
                    return tFuncIndex;
                    break;
                }
                case( AnalyticType::SPIRAL ):
                {
                    MORIS_ASSERT( aConstants.size() == 10, "Analytic::set_analytical_function() - incorrect size for constants list " );
                    MORIS_ASSERT( aConstants(1) > 0, "Analytic::set_analytical_function() - fiber radius must be > 0 " );
                    tFuncIndex = this->set_analytical_function( spiral_function, aConstants );
                    return tFuncIndex;
                    break;
                }
                default:
                {
                        MORIS_ERROR(false, "cl_GE_Analytical::set_analytical_function() - please choose from standardized functions or pass in your own function");
                        return tFuncIndex;
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
            set_analytical_function_dphi_dp( AnalyticType aGeomType, moris_index aFuncIndex )
            {
                switch(aGeomType)
                {
                case( AnalyticType::CIRCLE ):
                {
                    mPropertyList( aFuncIndex ).set_derivative_functions( circle_function_dphi_dp );
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
                    mPropertyList( aFuncIndex ).set_derivative_functions( sphere_function_dphi_dp );
//                    this->set_analytical_function_dphi_dp( sphere_function_dphi_dp );
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
             *index
             * @param[in] aPoint - point for value to be determined at
             * @oaram[in] aSubIndex - index to the sub type which is to be evaluated
             *
             * @param[out] real - function value at specified coordinate
             */
			real
			get_field_val_at_coordinate( const Matrix< DDRMat >  & aPoint,
			                             const moris_index aSubIndex = 0 )
			{
			    return mPropertyList( aSubIndex ).get_field_val_at_coordinate( aPoint );
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
			//TODO update/add default argument values if needed
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
             * @oaram[in] aSubIndex - index to the sub-type which is to be evaluated
             *
             * @param[out] real - function sensitivity value at specified coordinate
             */
			Matrix< DDRMat >
			get_sensitivity_dphi_dp_at_coordinate( const Matrix< DDRMat >  & aPoint,
			                                       const moris_index aSubIndex = 0 )
            {
			    return mPropertyList( aSubIndex ).get_sensitivity_dphi_dp_at_coordinate( aPoint );
            };
			//------------------------------------------------------------------------------
			uint get_number_of_sub_types()
			{
			    return mPropertyList.size();
			};
			//------------------------------------------------------------------------------

        private:
			mtk::Mesh_Manager*  mMesh = nullptr;

			std::shared_ptr< moris::hmr::Mesh >  mMesh_HMR = nullptr;               //FIXME delete this one

            fem::Interpolation_Type  mMySpaceInterpType;
            mtk::Interpolation_Order mMySpaceInterpOrder;
            fem::Interpolation_Type  mMyTimeInterpType;
            mtk::Interpolation_Order mMyTimeInterpOrder;

            moris::Cell< Property > mPropertyList;
            /*
             * this is needed for tricky functions which require multiple cell inputs (e.g. multi_cylinder_function)
             */
            real ( *mFuncAnalyticExtra )( const Matrix< DDRMat >        & aCoordinates,
                                                Cell<Cell<moris::real>> & aCenter,
                                                Cell<moris::real>       & aRadius,
                                                Cell<moris::real>       & aLength,
                                                Cell<Cell<moris::real>> & aAxis ) = nullptr;
            //------------------------------------------------------------------------------

        protected:

		};
	} /* namespace gen */
} /* namespace moris */



#endif /* PROJECTS_GEN_SRC_CL_GE_ANALYTIC_HPP_ */
