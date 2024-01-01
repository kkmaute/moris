/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_FEM_SP_Crosswind.cpp
 *
 */

#include "cl_FEM_SP_Crosswind.hpp"
#include "cl_FEM_Cluster.hpp"
#include "cl_FEM_Field_Interpolator_Manager.hpp"
// LINALG/src
#include "fn_norm.hpp"
#include "fn_dot.hpp"
#include "fn_eye.hpp"

namespace moris
{
    namespace fem
    {

        //------------------------------------------------------------------------------

        SP_Crosswind::SP_Crosswind(){}

        //------------------------------------------------------------------------------

        void
        SP_Crosswind::set_parameters( moris::Vector< Matrix< DDRMat > > aParameters )
        {
            // set mParameters
            mParameters = aParameters;

            // get number of parameters
            uint tParamSize = aParameters.size();

            // check for proper size of constant function parameters
            MORIS_ERROR( tParamSize > 0 && tParamSize < 3,
                    "SP_Crosswind::set_parameters - 1 or 2 constant parameters need to be set." );

            // check for proper parameter type; here just a scalar
            MORIS_ERROR( aParameters( 0 ).numel() == 1,
                    "SP_Crosswind::set_parameters - 1st parameter is not a scalar but a vector." );

            // check for proper parameter value
            MORIS_ERROR( mParameters( 0 )( 0 ) > 0.0,
                    "SP_Crosswind::set_parameters - C parameter needs to be larger than zero." );

            // set C
            mC = mParameters( 0 )( 0 );

            // if more than one parameter
            if( tParamSize > 1 )
            {
                // check for proper parameter type; here just a scalar
                MORIS_ERROR( aParameters( 1 ).numel() == 1,
                        "SP_Crosswind::set_parameters - 2nd parameter is not a scalar but a vector." );

                // set epsilon
                mEpsilon = mParameters( 1 )( 0 );
            }
        }

        //------------------------------------------------------------------------------

        void
        SP_Crosswind::eval_SP()
        {
            // evaluate the metric tensor Gij = sum_k dxi_k/dx_i dxi_k/dx_j
            const Matrix< DDRMat > & tG = //
                    mLeaderFIManager->get_IP_geometry_interpolator()->metric_tensor();

            // get flattened G to row vector
            Matrix< DDRMat > tFlatG = vectorize( tG );

            // get trace of G
            real tGijGij  = dot( tFlatG, tFlatG );

            // check against zero
            real tDeno = std::max( std::pow( tGijGij, 0.25 ), mEpsilon );

            // check against zero
            real tH = 1.0 / tDeno ;

            // compute stabilization parameter and set tolerance
            mPPVal = { { 0.5 * mC * tH, mEpsilon } };
        }

        //------------------------------------------------------------------------------

        void
        SP_Crosswind::eval_dSPdLeaderDOF(
                const moris::Vector< MSI::Dof_Type >& aDofTypes )
        {
            // get the dof type index
            const uint tDofIndex = mLeaderGlobalDofTypeMap( static_cast< uint >( aDofTypes( 0 ) ) );

            // get the dof type FI
            Field_Interpolator* tFIDer =
                    mLeaderFIManager->get_field_interpolators_for_type( aDofTypes( 0 ) );

            // get the number of coefficients
            uint tNumCoeff = tFIDer->get_number_of_space_time_coefficients();

            // set size for dSPdLeaderDof
            mdPPdLeaderDof( tDofIndex ).set_size( 2, tNumCoeff, 0.0 );
        }

        //------------------------------------------------------------------------------
    } /* namespace fem */
} /* namespace moris */

