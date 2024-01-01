/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_FEM_IQI_Volume.cpp
 *
 */

#include "cl_FEM_IQI_Volume.hpp"
#include "cl_FEM_Set.hpp"
#include "cl_FEM_Field_Interpolator_Manager.hpp"

namespace moris
{
    namespace fem
    {
        //------------------------------------------------------------------------------

        IQI_Volume::IQI_Volume()
        {
            // set size for the property pointer cell
            mLeaderProp.resize( static_cast< uint >( IQI_Property_Type::MAX_ENUM ), nullptr );

            // populate the property map
            mPropertyMap[ "Density" ] = static_cast< uint >( IQI_Property_Type::DENSITY );
        }

        //------------------------------------------------------------------------------

        void IQI_Volume::compute_QI( Matrix< DDRMat > & aQI )
        {
            // check the point is inside the bounded box
            if ( !this->is_within_box_bounds() )
            {
                return;
            }

            // get density property
            std::shared_ptr< Property > & tPropDensity =
                    mLeaderProp( static_cast< uint >( IQI_Property_Type::DENSITY ) );

            // if density property
            if ( tPropDensity != nullptr )
            {
                // evaluate the density
                aQI = tPropDensity->val();
            }
            else
            {
                // set density to 1
                aQI = {{ 1.0 }};
            }
        }

        //------------------------------------------------------------------------------

        void IQI_Volume::compute_QI( real aWStar )
        {
            // check the point is inside the bounded box
            if ( !this->is_within_box_bounds() )
            {
                return;
            }

            // get index for QI
            sint tQIIndex = mSet->get_QI_assembly_index( mName );

            // get density property
            std::shared_ptr< Property > & tPropDensity =
                    mLeaderProp( static_cast< uint >( IQI_Property_Type::DENSITY ) );

            // if density property
            if ( tPropDensity != nullptr )
            {
                // evaluate the density
                mSet->get_QI()( tQIIndex ) += aWStar * ( tPropDensity->val() );
            }
            else
            {
                // set density to 1
                mSet->get_QI()( tQIIndex ) += aWStar;
            }
        }

        //------------------------------------------------------------------------------

        void IQI_Volume::compute_dQIdu( real aWStar )
        {
            // check the point is inside the bounded box
            if ( !this->is_within_box_bounds() )
            {
                return;
            }

            // get the column index to assemble in residual
            sint tQIIndex = mSet->get_QI_assembly_index( mName );

            // get density property
            const std::shared_ptr< Property > & tPropDensity =
                    mLeaderProp( static_cast< uint >( IQI_Property_Type::DENSITY ) );

            // get the number of leader dof type dependencies
            uint tNumDofDependencies = mRequestedLeaderGlobalDofTypes.size();

            // compute dQIdu for indirect dof dependencies
            for( uint iDof = 0; iDof < tNumDofDependencies; iDof++ )
            {
                // get the treated dof type
                Vector< MSI::Dof_Type > & tDofType = mRequestedLeaderGlobalDofTypes( iDof );

                // get leader index for residual dof type, indices for assembly
                uint tLeaderDofIndex      = mSet->get_dof_index_for_type( tDofType( 0 ), mtk::Leader_Follower::LEADER );
                uint tLeaderDepStartIndex = mSet->get_res_dof_assembly_map()( tLeaderDofIndex )( 0, 0 );
                uint tLeaderDepStopIndex  = mSet->get_res_dof_assembly_map()( tLeaderDofIndex )( 0, 1 );

                // Dof dependency
                if ( tPropDensity != nullptr && tPropDensity->check_dof_dependency( tDofType ) )
                {
                    // compute dQIdu
                    mSet->get_residual()( tQIIndex )(
                            { tLeaderDepStartIndex, tLeaderDepStopIndex },
                            { 0, 0 } ) += aWStar * ( tPropDensity->dPropdDOF( tDofType ) );
                }
            }
        }

        //------------------------------------------------------------------------------

        void IQI_Volume::compute_dQIdu(
                moris::Vector< MSI::Dof_Type > & aDofType,
                Matrix< DDRMat >             & adQIdu )
        {
            // check the point is inside the bounded box
            if ( !this->is_within_box_bounds() )
            {
                return;
            }

            // get density property
            std::shared_ptr< Property > & tPropDensity =
                    mLeaderProp( static_cast< uint >( IQI_Property_Type::DENSITY ) );

            // Dof dependency
            if ( tPropDensity != nullptr && tPropDensity->check_dof_dependency( aDofType ) )
            {
                // compute dQIdu
                adQIdu = tPropDensity->dPropdDOF( aDofType );
            }
        }

        //------------------------------------------------------------------------------

        bool IQI_Volume::is_within_box_bounds()
        {
                // check if the box bounds are empty then skip
                if ( mParameters.empty() ) 
                {
                        return true;
                }

                //if the box bounds are not empty then check if it is inside the box
                else
                {
                        // get the coordinate 
                        const Matrix<DDRMat> & tGaussPoint = mLeaderFIManager->get_IG_geometry_interpolator()->valx();

                        //check if the calculation point coordinates are more then lower corner of the box 
                        bool tLowerBound = std::equal(mParameters(0).begin(), mParameters(0).end(), tGaussPoint.begin(),tGaussPoint.end() , 
                        [](real aA, real aB) -> bool { return aA < aB ;} );
                        
                        //check if the calculation point coordinates are less then upper corner of the box 
                        bool tUpperBound = std::equal(tGaussPoint.begin(),tGaussPoint.end(), mParameters(1).begin(), mParameters(1).end(),
                         [](real aA, real aB) -> bool { return aA < aB ;} );
                        
                        //combine the two bounds that satisfy both
                        return tUpperBound and tLowerBound; 
                }
        }

    }/* end namespace fem */
}/* end namespace moris */

