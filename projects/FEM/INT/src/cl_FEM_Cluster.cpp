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

        void Cluster::compute_jacobian()
        {
            // set the IP geometry interpolator physical space and time coefficients
            mSet->get_block_IP_geometry_interpolator()->set_space_coeff( mInterpolationCell->get_vertex_coords() );
            mSet->get_block_IP_geometry_interpolator()->set_time_coeff( this->mTime );

            // set the IP geometry interpolator param space and time coefficients
            // fixme param coeff from cluster
            mSet->get_block_IP_geometry_interpolator()->set_param_coeff();

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

//------------------------------------------------------------------------------

        void Cluster::compute_residual()
        {
            // set the IP geometry interpolator physical space and time coefficients
            mSet->get_block_IP_geometry_interpolator()->set_space_coeff( mInterpolationCell->get_vertex_coords() );
            mSet->get_block_IP_geometry_interpolator()->set_time_coeff( this->mTime );

            // set the IP geometry interpolator param space and time coefficients
            // fixme param coeff from cluster
            mSet->get_block_IP_geometry_interpolator()->set_param_coeff();

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

//------------------------------------------------------------------------------

//------------------------------------------------------------------------------
    } /* namespace fem */
} /* namespace moris */
