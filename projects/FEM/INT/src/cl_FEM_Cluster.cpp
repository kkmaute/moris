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
            for( auto tElements : mElements )
            {
                delete tElements;
            }
            mElements.clear();
        }

//------------------------------------------------------------------------------
        void Cluster::compute_jacobian()
         {

            // set the IP geometry interpolator physical space and time coefficients for the master interpolation cell
            mSet->get_IP_geometry_interpolator( mtk::Master_Slave::MASTER )->set_space_coeff( mMasterInterpolationCell->get_vertex_coords() );
            mSet->get_IP_geometry_interpolator( mtk::Master_Slave::MASTER )->set_time_coeff( this->mTime );

            // if double side cluster
             if( mElementType == fem::Element_Type::DOUBLE_SIDESET )
             {
                 // set the IP geometry interpolator physical space and time coefficients for the slave interpolation cell
                 mSet->get_IP_geometry_interpolator( mtk::Master_Slave::SLAVE )->set_space_coeff( mSlaveInterpolationCell->get_vertex_coords() );
                 mSet->get_IP_geometry_interpolator( mtk::Master_Slave::SLAVE )->set_time_coeff( this->mTime );
             }

             //Fixme do this only once
             this->get_my_pdof_values();

             // init the jacobian //fixme still ok?
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
            // set the IP geometry interpolator physical space and time coefficients for the master interpolation cell
            mSet->get_IP_geometry_interpolator( mtk::Master_Slave::MASTER )->set_space_coeff( mMasterInterpolationCell->get_vertex_coords() );
            mSet->get_IP_geometry_interpolator( mtk::Master_Slave::MASTER )->set_time_coeff( this->mTime );

            // if double side cluster
            if( mElementType == fem::Element_Type::DOUBLE_SIDESET )
            {
                // set the IP geometry interpolator physical space and time coefficients for the slave interpolation cell
                mSet->get_IP_geometry_interpolator( mtk::Master_Slave::SLAVE )->set_space_coeff( mSlaveInterpolationCell->get_vertex_coords() );
                mSet->get_IP_geometry_interpolator( mtk::Master_Slave::SLAVE )->set_time_coeff( this->mTime );
            }

            //Fixme do this only once
            this->get_my_pdof_values();

            // init the residual
            mSet->initialize_mResidual();

            // set the field interpolators coefficients
            this->set_field_interpolators_coefficients();

            // loop over the elements
            for ( uint Ik = 0; Ik < mElements.size(); Ik++ )
            {
                // compute the jacobian for the element
                mElements( Ik )->compute_residual();
            }
        }

//------------------------------------------------------------------------------

//------------------------------------------------------------------------------
    } /* namespace fem */
} /* namespace moris */
