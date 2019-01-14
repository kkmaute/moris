/*
 * cl_GE_Analytical.hpp
 *
 *  Created on: Jan 4, 2019
 *      Author: sonne
 */

#ifndef PROJECTS_GEN_SRC_CL_GE_ANALYTICAL_HPP_
#define PROJECTS_GEN_SRC_CL_GE_ANALYTICAL_HPP_

namespace moris{
	namespace ge{

		class Analytical : public Geometry
		{
		private:

	        real ( *mFunctionAnalytical )( const Matrix< DDRMat > & aPoint, moris::Cell< real> aConstant );

		protected:

		public:
			Analytical(){};
			~Analytical(){};
            //------------------------------------------------------------------------------

			real get_val_at_vertex( const Matrix< DDRMat > & aPoint, moris::Cell< real > aConst )
			{
				return mFunctionAnalytical( aPoint, aConst );
			};

            //------------------------------------------------------------------------------

		    void set_analytical_function( real ( *funcPointer )( const Matrix< DDRMat > & aPoint, Cell< real > aConst ) )
			{
		    	mFunctionAnalytical = funcPointer;
			};

		    //------------------------------------------------------------------------------


		};
	} /* namespace gen */
} /* namespace moris */



#endif /* PROJECTS_GEN_SRC_CL_GE_ANALYTICAL_HPP_ */
