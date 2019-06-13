/*
 * cl_FEM_Cluster.hpp
 *
 *  Created on: Apr 11, 2019
 *      Author: schmidt
 */
#include <iostream>

#include "cl_FEM_Element.hpp" //FEM/INT/src
#include "cl_FEM_Cluster.hpp"   //FEM/INT/src

namespace moris
{
    namespace fem
    {
//------------------------------------------------------------------------------

//------------------------------------------------------------------------------

    Cluster::~Cluster()
    {
      mElements.clear();
    }

//------------------------------------------------------------------------------

        void Cluster::compute_jacobian()
        {
            // if double side cluster
            if( mElementType == fem::Element_Type::DOUBLE_SIDESET )
            {
                // set the IP geometry interpolator physical space and time coefficients for the left interpolation cell
                mSet->get_left_IP_geometry_interpolator()->set_space_coeff( mLeftInterpolationCell->get_vertex_coords() );
                mSet->get_left_IP_geometry_interpolator()->set_time_coeff( this->mTime );

//                // set the IP geometry interpolator param space and time coefficients for the left interpolation cell
//                // fixme param coeff from cluster
//                mSet->get_left_IP_geometry_interpolator()->set_param_coeff();

                // set the IP geometry interpolator physical space and time coefficients for the right interpolation cell
                mSet->get_right_IP_geometry_interpolator()->set_space_coeff( mRightInterpolationCell->get_vertex_coords() );
                mSet->get_right_IP_geometry_interpolator()->set_time_coeff( this->mTime );

//                // set the IP geometry interpolator param space and time coefficients for the right interpolation cell
//                // default implementation
//                mSet->get_right_IP_geometry_interpolator()->set_param_coeff();

                //Fixme do this only once
                this->get_my_pdof_values();

                // init the jacobian //fixme still ok?
                mSet->initialize_mJacobian();

                // set the field interpolators coefficients
                this->set_field_interpolators_coefficients_double();

                // loop over the elements
                for ( uint Ik = 0; Ik < mElements.size(); Ik++ )
                {
                    // compute the jacobian for the element
                    mElements( Ik )->compute_jacobian();
                }
            }
            // if cell or side cluster
            else
            {
                // set the IP geometry interpolator physical space and time coefficients
                mSet->get_IP_geometry_interpolator()->set_space_coeff( mInterpolationCell->get_vertex_coords() );
                mSet->get_IP_geometry_interpolator()->set_time_coeff( this->mTime );

//                // set the IP geometry interpolator param space and time coefficients
//                // default implementation
//                mSet->get_IP_geometry_interpolator()->set_param_coeff();

                //Fixme do this only once
                this->get_my_pdof_values();

                // init the jacobian
                mSet->initialize_mJacobian();

                // set the field interpolators coefficients
                this->set_field_interpolators_coefficients();

                // loop over the elements
                for ( uint Ik = 0; Ik < mElements.size(); Ik++ )
                {
                    // compute the jacobian for the element
                    mElements( Ik )->compute_jacobian();
                }
            }
        }

//------------------------------------------------------------------------------

        void Cluster::compute_residual()
        {
            // if double side cluster
            if( mElementType == fem::Element_Type::DOUBLE_SIDESET )
            {
                // set the IP geometry interpolator physical space and time coefficients for the left interpolation cell
                mSet->get_left_IP_geometry_interpolator()->set_space_coeff( mLeftInterpolationCell->get_vertex_coords() );
                mSet->get_left_IP_geometry_interpolator()->set_time_coeff( this->mTime );

//                // set the IP geometry interpolator param space and time coefficients for the left interpolation cell
//                // fixme param coeff from cluster
//                mSet->get_left_IP_geometry_interpolator()->set_param_coeff();

                // set the IP geometry interpolator physical space and time coefficients for the right interpolation cell
                mSet->get_right_IP_geometry_interpolator()->set_space_coeff( mRightInterpolationCell->get_vertex_coords() );
                mSet->get_right_IP_geometry_interpolator()->set_time_coeff( this->mTime );

//                // set the IP geometry interpolator param space and time coefficients for the right interpolation cell
//                // fixme param coeff from cluster
//                mSet->get_right_IP_geometry_interpolator()->set_param_coeff();

                //Fixme do this only once
                this->get_my_pdof_values();

                // init the residual
                mSet->initialize_mResidual();

                // set the field interpolators coefficients
                this->set_field_interpolators_coefficients_double();

                // loop over the elements
                for ( uint Ik = 0; Ik < mElements.size(); Ik++ )
                {
                    // compute the jacobian for the element
                    mElements( Ik )->compute_residual();
                }
            }
            // if cell ore side cluster
            else
            {
                // set the IP geometry interpolator physical space and time coefficients
                mSet->get_IP_geometry_interpolator()->set_space_coeff( mInterpolationCell->get_vertex_coords() );
                mSet->get_IP_geometry_interpolator()->set_time_coeff( this->mTime );

//                // set the IP geometry interpolator param space and time coefficients
//                // default implementation
//                mSet->get_IP_geometry_interpolator()->set_param_coeff();

                // fixme do this only once
                this->get_my_pdof_values();

                // init the residual
                mSet->initialize_mResidual();

                // set the field interpolators coefficients
                this->set_field_interpolators_coefficients();

                // loop over the elements
                for ( uint Ik = 0; Ik < mElements.size(); Ik++ )
                {
                    // compute the residual for the element
                    mElements( Ik )->compute_residual();
                }
            }
        }

//------------------------------------------------------------------------------

//------------------------------------------------------------------------------
    } /* namespace fem */
} /* namespace moris */
