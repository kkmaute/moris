/*
 * cl_GE_Geometry.hpp
 *
 *  Created on: Dec 28, 2018
 *      Author: sonne
 */

#ifndef PROJECTS_MTK_GE_SRC_CL_GE_BASE_HPP_
#define PROJECTS_MTK_GE_SRC_CL_GE_BASE_HPP_

// GE includes
#include "cl_GE_Element.hpp"
#include "cl_GE_Enums.hpp"
#include "cl_GE_Interface.hpp"
#include "cl_GE_Geometry_Library.hpp"
//------------------------------------------------------------------------------
// linalg includes
#include "cl_Matrix.hpp"
#include "linalg_typedefs.hpp"
#include "fn_norm.hpp"
//------------------------------------------------------------------------------
// other includes
#include <cmath>
//------------------------------------------------------------------------------

namespace moris
{
namespace ge
{
	class Geometry
	{

	public:
		Geometry(){};

		~Geometry(){};
		/*
		 * *****************************************************************************
		 * ************************ ANALYTIC GEOMETRY FUNCTIONS ************************
		 * *****************************************************************************
		 */
        //------------------------------------------------------------------------------
		virtual
		void set_analytical_function( real ( *mFuncAnalytic )( const Matrix< DDRMat > & aPoint, moris::Cell< real> aConstant ) )
		{
			std::cout<<"set_analytical_function(): please specify your own analytic function"<<std::endl;
		};
		//------------------------------------------------------------------------------
		virtual
        void set_analytical_function( real ( *funcPointer )( const Matrix< DDRMat >        & aCoordinates,
                                                                   Cell<Cell<moris::real>> & aCenter,
                                                                   Cell<moris::real>       & aRadius,
                                                                   Cell<moris::real>       & aLength,
                                                                   Cell<Cell<moris::real>> & aAxis ) )
        {
		    std::cout<<"set_analytical_function(): please choose a valid function"<<std::endl;
        };
        //------------------------------------------------------------------------------
        virtual
        void set_analytical_function( type aGeomType )
        {
            std::cout<<"set_analytical_function(): please choose a valid function"<<std::endl;
        };

        //------------------------------------------------------------------------------
        virtual
        void set_analytical_function_dphi_dx( Matrix< DDRMat > ( *mFuncAnalyticDphiDx )( const Matrix< DDRMat > & aPoint, Cell< real > aConst ) )
        {
            std::cout<<"set_analytical_function_dphi_dx(): please specify your own analytic function dphi/dx"<<std::endl;
        };

        //------------------------------------------------------------------------------
        virtual
        void set_analytical_function_dphi_dx( type aGeomType )
        {
            std::cout<<"set_analytical_function_dphi_dx(): please choose a valid dphi/dx function"<<std::endl;
        };

        //------------------------------------------------------------------------------
        virtual
        real get_field_val_at_coordinate( const Matrix< DDRMat > & aPoint, moris::Cell< real > aConst = 0.0 )
        {
            std::cout<<"get_field_val_at_coordinate(): function not implemented"<<std::endl;
            return 0.0;
        };

        //------------------------------------------------------------------------------

        virtual
        Matrix< DDRMat >
        get_sensitivity_dphi_dp_at_coordinate( const Matrix< DDRMat >  & aPoint,
                                                     moris::Cell< real > aConst)
        {
            Matrix< DDRMat > tSensitivityDxDp(4, 3, 0.0);
            return tSensitivityDxDp;
        };

        //------------------------------------------------------------------------------
        /*
         * *****************************************************************************
         * ************************ DISCRETE GEOMETRY FUNCTIONS ************************
         * *****************************************************************************
         */
        //------------------------------------------------------------------------------




	//------------------------------------------------------------------------------
    private:

	//------------------------------------------------------------------------------
    protected:

	};


//------------------------------------------------------------------------------

} /* namespace gen */
} /* namespace moris */


#endif /* PROJECTS_GEN_SRC_CL_GE_BASE_HPP_ */
