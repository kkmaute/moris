/*
 * cl_FEM_Interpolator.hpp
 *
 *  Created on: Jul 13, 2018
 *      Author: messe
 */

#ifndef SRC_FEM_CL_FEM_INTERPOLATOR_HPP_
#define SRC_FEM_CL_FEM_INTERPOLATOR_HPP_

#include "cl_FEM_Interpolation_Function_Base.hpp"
#include "typedefs.hpp" //MRS/COR/src
#include "cl_Matrix.hpp" //LNA/src
#include "linalg_typedefs.hpp" //LNA/src

#include "cl_FEM_Interpolation_Rule.hpp" //FEM/INT/src
#include "cl_FEM_Interpolation_Matrix.hpp" //FEM/INT/src
#include "cl_FEM_Geometry_Interpolator.hpp" //FEM/INT/src
#include "cl_FEM_Integrator.hpp" //FEM/INT/src

namespace moris
{
    namespace fem
    {
//------------------------------------------------------------------------------

        // forward declaration of element
        class Element;

//------------------------------------------------------------------------------

        class Interpolator
        {
            //! how many fields are to be interpolated
            const uint mNumberOfFields;

            //! pointer to space interpolation object
            Interpolation_Function_Base * mSpaceInterpolation = nullptr;

            //! pointer to time interpolation object
            Interpolation_Function_Base * mTimeInterpolation  = nullptr;

            //! pointer to space time interpolation object
            Interpolation_Function_Base * mSpaceTimeInterpolation  = nullptr;

            //! pointer to function that creates matrices
            Interpolation_Function_Base * mMatrixCreator = nullptr;

            //! geometry interpolation object
            Geometry_Interpolator       * mGeometryInterpolator = nullptr;

            //! pointer to integrator object (must be pointer because may be not assigned
            //                                  in alternative constructor )
            Integrator                  * mIntegrator = nullptr;

            //! container for node coordinates
            Matrix< DDRMat > mNodeCoords;

            //! container for integration points
            Matrix< DDRMat > mIntegrationPoints;

            //! container for integration weights
            Matrix< DDRMat > mIntegrationWeights;

            //! flag telling if element is isoparametric
            bool mIsoparametricFlag;

            //! transposed of Geometry Jacobi
            Matrix< DDRMat > mJt;

            //! last point for which mJt was evaluated
            Matrix< DDRMat > mLastPointJt;

            //! dmatrix for the geometry ( needs destructor )
            Interpolation_Matrix * mGN = nullptr;

            //! derivative matrix for the geometry ( needs destructor )
            Interpolation_Matrix * mGdNdXi = nullptr;

//------------------------------------------------------------------------------
        public:
//------------------------------------------------------------------------------

            /**
             * constructor with two interpolation rules and integrator
             *
             * @param[ in ] aSpaceInterpolationRule   pointer to space interpolation rule
             * @param[ in ] aTimeInterpolationRule    pointer to time interpolation rule
             * @param[ in ] aGeometryInterpolator     pointer to geometry interpolation object
             * @param[ in ] aIntegrator               pointer to integration object
             *
             */
            Interpolator(
                          Element               * aElement,
                    const uint                  & aNumberOfFields,
                    const Interpolation_Rule    & aFieldInterpolationRule,
                    const Interpolation_Rule    & aGeometryInterpolationRule,
                    const Integration_Rule      & aIntegrationRule );

//------------------------------------------------------------------------------

            /**
             * default constructor
             */
            ~Interpolator();

//------------------------------------------------------------------------------


            /**
             * a constructor for the interpolation matrix object
             *
             * @param[ in ] derivative in space
             * @param[ in ] derivative in time  ( currently unused )
             * @param[ in ] type of coeffs  :
             *               0 : Interpolated values
             *               1 : shape function matrix
             *               2 : geometry jacobian
             */
            Interpolation_Matrix *
            create_matrix( const uint & aDerivativeInSpace,
                    	   const uint & aDerivativeInTime  );

//------------------------------------------------------------------------------

            /**
             * returns the number of DOFs of this interpolation
             */
            uint
            get_number_of_dofs();

//------------------------------------------------------------------------------

            /**
             * returns the number of points of the integrator
             */
            uint
            get_number_of_integration_points();

//------------------------------------------------------------------------------

            /**
             * returns the coordinates of an integration point
             */
            Matrix< DDRMat >
            get_point( const uint & aPoint )
            {
                  return mIntegrationPoints.get_column( aPoint );
            }

//------------------------------------------------------------------------------

            const Matrix< DDRMat > &
            get_node_coords() const
            {
                return mNodeCoords;
            }

//------------------------------------------------------------------------------

            /**
             * returns the integration weight of this point
             */
            real
            get_integration_weight( const uint & aPoint );

//------------------------------------------------------------------------------

            /**
             * returns the determinant of the geometry Jacobian by index
             */
            real
            get_det_J( const uint & aPoint );

//------------------------------------------------------------------------------

            /**
             * returns the determinant of the geometry Jacobian by point
             */
            real
            get_det_J( const Matrix< DDRMat > & aPoint );

//------------------------------------------------------------------------------

            /**
             * evaluate the geometry coordinates of a point
             */
            Matrix< DDRMat >
            eval_geometry_coords( const Matrix< DDRMat > & aPoint );

            Matrix< DDRMat >
            eval_geometry_coords( const uint & aPoint );

//------------------------------------------------------------------------------

            void
            eval_N( Interpolation_Matrix 		& aMatrix,
                    const Matrix< DDRMat >    	& aPoint );

//------------------------------------------------------------------------------
            void
            eval_dNdx( Interpolation_Matrix 	& aMatrix,
                       const Matrix< DDRMat >   & aPoint );

//------------------------------------------------------------------------------
        };

//------------------------------------------------------------------------------

        // free function called by interpolation matrix
        void
        interpolator_eval_N(
                Interpolator          	* aInterpolator,
                Interpolation_Matrix  	* aMatrix,
                const Matrix< DDRMat >  & aPoint )
        {
            aInterpolator->eval_N( *aMatrix, aPoint );
        }

//------------------------------------------------------------------------------

        // free function called by interpolation matrix
        void
        interpolator_eval_dNdx(
                Interpolator          	* aInterpolator,
                Interpolation_Matrix  	* aMatrix,
                const Matrix< DDRMat >  & aPoint )
        {
            aInterpolator->eval_dNdx( *aMatrix, aPoint );
        }


//------------------------------------------------------------------------------
    } /* namespace fem */
} /* namespace moris */

#endif /* SRC_FEM_CL_FEM_INTERPOLATOR_HPP_ */
