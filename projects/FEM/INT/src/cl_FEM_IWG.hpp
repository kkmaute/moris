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
#include "cl_FEM_Enums.hpp"        //FEM/MSI/src


namespace moris
{
    namespace fem
    {
//------------------------------------------------------------------------------
        /**
         * Integrand of Weak Form of Governing Equations
         */
        class IWG
        {
        protected :

            // nodal weak BCs
            Matrix< DDRMat > mNodalWeakBCs;

            // normal
            Matrix< DDRMat > mNormal;

            // residual dof type
            moris::Cell< MSI::Dof_Type > mResidualDofType;

            // active dof types
            moris::Cell< moris::Cell< MSI::Dof_Type > > mActiveDofTypes;

            // active model parameter type
            moris::Cell< fem::Mp_Type > mActiveMpTypes;

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
            /**
             * set nodal weak BCs
             * @param[ in ] aNodalWeakBCs matrix with nodal values
             */
            void set_nodal_weak_bcs( Matrix< DDRMat > & aNodalWeakBCs )
            {
                mNodalWeakBCs = aNodalWeakBCs;
            }

            /**
             * set normal
             * @param[ in ] aNormal normal vector
             */
            void set_normal( Matrix< DDRMat > & aNormal )
            {
                mNormal = aNormal;
            }

//------------------------------------------------------------------------------
            /**
             * returns a dof type for the residual
             */
            moris::Cell< MSI::Dof_Type > get_residual_dof_type() const
            {
                return mResidualDofType;
            };

//------------------------------------------------------------------------------
            /**
             * returns a cell of dof types used to evaluate the residual
             * and the jacobian
             */
            moris::Cell< moris::Cell< MSI::Dof_Type > > get_active_dof_types() const
            {
                return mActiveDofTypes;
            };

//------------------------------------------------------------------------------
            /**
             * returns a cell of mp types used to evaluate the residual
             * and the jacobian
             */
            moris::Cell< fem::Mp_Type > get_active_mp_types() const
            {
                return mActiveMpTypes;
            };

//------------------------------------------------------------------------------
            /**
             * evaluates the residual
             */
            virtual void compute_residual( Matrix< DDRMat >                   & aResidual,
                                           moris::Cell< Field_Interpolator* > & aFieldInterpolators)
            {
                MORIS_ERROR( false, "IWG::compute_residual - This function does nothing. " );
            }

            virtual void compute_residual( Matrix< DDRMat >                   & aResidual,
                                           moris::Cell< Field_Interpolator* > & aLeftFieldInterpolators,
                                           moris::Cell< Field_Interpolator* > & aRightFieldInterpolators)
            {
                MORIS_ERROR( false, "IWG::compute_residual - This function does nothing. " );
            }

//------------------------------------------------------------------------------
            /**
             * evaluates the Jacobian
             */
            virtual void compute_jacobian( moris::Cell< Matrix< DDRMat > >    & aJacobians,
                                           moris::Cell< Field_Interpolator* > & aFieldInterpolators)
            {
                MORIS_ERROR( false, "IWG::compute_jacobian - This function does nothing. " );
            }

            virtual void compute_jacobian( moris::Cell< Matrix< DDRMat > >    & aJacobians,
                                           moris::Cell< Field_Interpolator* > & aLeftFieldInterpolators,
                                           moris::Cell< Field_Interpolator* > & aRightFieldInterpolators )
            {
                MORIS_ERROR( false, "IWG::compute_jacobian - This function does nothing. " );
            }
//------------------------------------------------------------------------------
            /**
             * evaluates the residual and the Jacobian
             */
            virtual void compute_jacobian_and_residual( moris::Cell< Matrix< DDRMat > >    & aJacobians,
                                                        Matrix< DDRMat >                   & aResidual,
                                                        moris::Cell< Field_Interpolator* > & aFieldInterpolators )
            {
                MORIS_ERROR( false, " IWG::compute_jacobian_and_residual - This function does nothing. " );
            }

            virtual void compute_jacobian_and_residual( moris::Cell< Matrix< DDRMat > >    & aJacobians,
                                                        Matrix< DDRMat >                   & aResidual,
                                                        moris::Cell< Field_Interpolator* > & aLeftFieldInterpolators,
                                                        moris::Cell< Field_Interpolator* > & aRightFieldInterpolators)
            {
                MORIS_ERROR( false, " IWG::compute_jacobian_and_residual - This function does nothing. " );
            }
//------------------------------------------------------------------------------
//
//            virtual real compute_integration_error( const Matrix< DDRMat > & aNodalDOF,
//                                                    real (*aFunction)( const Matrix< DDRMat > & aPoint ) ,
//                                                    const uint        & aPointIndex )
//            {
//                MORIS_ERROR( false, "This function does nothing" );
//                return 0.0;
//            }

//------------------------------------------------------------------------------
//
//            real interpolate_scalar_at_point( const Matrix< DDRMat > & aNodalWeakBC,
//                                              const uint             & aPointIndex )
//            {
//                MORIS_ERROR( false, "This function does nothing" );
//                return 0.0;
//            }


        };
//------------------------------------------------------------------------------
    } /* namespace fem */
} /* namespace moris */
#endif /* SRC_FEM_CL_FEM_IWG_HPP_ */
