/*
 * cl_FEM_IWG.hpp
 *
 *  Created on: Aug 13, 2018
 *      Author: messe
 */
#ifndef SRC_FEM_CL_FEM_IWG_HPP_
#define SRC_FEM_CL_FEM_IWG_HPP_

#include "typedefs.hpp"                     //MRS/COR/src
#include "linalg_typedefs.hpp"              //MRS/COR/src
#include "cl_Cell.hpp"                      //MRS/CON/src
#include "cl_Matrix.hpp"                    //LNA/src
#include "cl_FEM_Field_Interpolator.hpp"    //FEM/INT/src
#include "cl_MSI_Dof_Type_Enums.hpp"        //FEM/MSI/src


namespace moris
{
    namespace fem
    {
//------------------------------------------------------------------------------
        //class Element;
        /**
         * Integrand of Weak Form of Governing Equations
         */
        class IWG
        {
        protected :

//            // element pointer
//            Element* mElement = nullptr;

            // nodal weak BCs
            Matrix< DDRMat > mNodalWeakBCs;

            // normal
            Matrix< DDRMat > mNormal;

            // residual dof type
            Cell< MSI::Dof_Type > mResidualDofType;

            // active dof types
            Cell< Cell< MSI::Dof_Type > > mActiveDofTypes;

            // FIXME temporary until other way
            uint mSpaceDim = 3;

//------------------------------------------------------------------------------
        public :

//------------------------------------------------------------------------------
            /**
             * constructor
             */
            IWG(){};

//------------------------------------------------------------------------------
            /**
             * virtual destructor
             */
            virtual ~IWG(){};

//------------------------------------------------------------------------------

//            void set_element_pointer( Element * aElement )
//            {
//                mElement = aElement;
//            }

            void set_nodal_weak_bcs( Matrix< DDRMat > & aNodalWeakBCs )
            {
                mNodalWeakBCs = aNodalWeakBCs;
            }

            void set_normal( Matrix< DDRMat > & aNormal )
            {
                mNormal = aNormal;
            }

//------------------------------------------------------------------------------

            /**
             * creates the required interpolation matrices as member variables
             */
            virtual void create_matrices( Field_Interpolator * aInterpolator )
            {
                MORIS_ASSERT( 0, "IWG::create_matrices - not implemented. " );
            }

//------------------------------------------------------------------------------

//            /**
//             * deletes the matrices
//             */
//            virtual void delete_matrices() = 0;

//------------------------------------------------------------------------------

//            /**
//             * returns a cell with the dof types, assuming that all nodes
//             * have the same type
//             */
//            virtual Cell< MSI::Dof_Type > get_dof_types() = 0;

//------------------------------------------------------------------------------
            /**
             * returns a dof type for the residual
             */
            Cell< MSI::Dof_Type > get_residual_dof_type() const
            {
                return mResidualDofType;
            };

//------------------------------------------------------------------------------
            /**
             * returns a cell of dof types used to evaluate the residual
             * and the jacobian
             */
            Cell< Cell< MSI::Dof_Type > > get_active_dof_types() const
            {
                return mActiveDofTypes;
            };

//------------------------------------------------------------------------------
            /**
             * evaluates the residual
             */
//            virtual void
//            compute_residual(       Matrix< DDRMat > & aResidual,
//                              const Matrix< DDRMat > & aNodalDOF,
//                              const Matrix< DDRMat > & aNodalWeakBC )
//            {
//                MORIS_ERROR( false, "This function does nothing" );
//             }

            virtual void compute_residual( Matrix< DDRMat >            & aResidual,
                                           Cell< Field_Interpolator* > & aFieldInterpolators)
            {
                MORIS_ERROR( false, "IWG::compute_residual - This function does nothing. " );
            }

//------------------------------------------------------------------------------
            /**
             * evaluates the Jacobian
             */
//            virtual void
//            compute_jacobian(       Matrix< DDRMat > & aJacobian,
//                              const Matrix< DDRMat > & aNodalDOF,
//                              const Matrix< DDRMat > & aNodalWeakBC )
//            {
//                MORIS_ERROR( false, "This function does nothing" );
//            }

            virtual void compute_jacobian( Cell< Matrix< DDRMat > >    & aJacobians,
                                           Cell< Field_Interpolator* > & aFieldInterpolators)
            {
                MORIS_ERROR( false, "IWG::compute_jacobian - This function does nothing. " );
            }
//------------------------------------------------------------------------------
//
//            /**
//             * evaluates the Jacobian
//             */
//            virtual void
//            compute_jacobian_and_residual(
//                    Matrix< DDRMat >       & aJacobian,
//                    Matrix< DDRMat >       & aResidual,
//                    const Matrix< DDRMat > & aNodalDOF,
//                    const Matrix< DDRMat > & aNodalWeakBC,
//                    const uint             & aPointIndex )
//            {
//                MORIS_ERROR( false, "This function does nothing" );
//            }
//
            /**
             * evaluates the residual and the Jacobian
             */
//            virtual void
//            compute_jacobian_and_residual(       Matrix< DDRMat > & aJacobian,
//                                                 Matrix< DDRMat > & aResidual,
//                                           const Matrix< DDRMat > & aNodalDOF,
//                                           const Matrix< DDRMat > & aNodalWeakBC )
//            {
//                MORIS_ERROR( false, "This function does nothing" );
//            }
            virtual void compute_jacobian_and_residual( Cell< Matrix< DDRMat > >    & aJacobians,
                                                        Matrix< DDRMat >            & aResidual,
                                                        Cell< Field_Interpolator* > & aFieldInterpolators )
            {
                MORIS_ERROR( false, " IWG::compute_jacobian_and_residual - This function does nothing. " );
            }
//------------------------------------------------------------------------------
//
//            virtual real
//            compute_integration_error(
//                    const Matrix< DDRMat > & aNodalDOF,
//                    real (*aFunction)( const Matrix< DDRMat > & aPoint ) ,
//                    const uint        & aPointIndex )
//            {
//                MORIS_ERROR( false, "This function does nothing" );
//                return 0.0;
//            }

//------------------------------------------------------------------------------
//
//            real
//            interpolate_scalar_at_point(
//                    const Matrix< DDRMat > & aNodalWeakBC,
//                    const uint             & aPointIndex )
//            {
//                MORIS_ERROR( false, "This function does nothing" );
//                return 0.0;
//            }
//
//        };

        };

//------------------------------------------------------------------------------
    } /* namespace fem */
} /* namespace moris */
#endif /* SRC_FEM_CL_FEM_IWG_HPP_ */
