/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_FEM_IQI_Strain_Energy.cpp
 *
 */

#include "cl_FEM_Set.hpp"
#include "cl_FEM_Field_Interpolator_Manager.hpp"
#include "cl_FEM_IQI_Strain_Energy.hpp"

namespace moris
{
    namespace fem
    {
        //------------------------------------------------------------------------------

        IQI_Strain_Energy::IQI_Strain_Energy()
        {
            // set fem IQI type
            mFEMIQIType = fem::IQI_Type::STRAIN_ENERGY;
            init_property( "Bedding", IWG_Property_Type::BEDDING );
            init_property( "Thickness", IWG_Property_Type::THICKNESS );
            init_constitutive_model( "Elast", IQI_Constitutive_Type::ELAST );
        }

        //------------------------------------------------------------------------------

        void IQI_Strain_Energy::compute_QI( Matrix< DDRMat > &aQI )
        {
            // check the point is inside the bounded box
            if ( !this->is_within_box_bounds() )
            {
                return;
            }

            // get the elasticity CM
           std::shared_ptr< Constitutive_Model > const &tCMElasticity = get_leader_constitutive_model(IQI_Constitutive_Type::ELAST);

            // get bedding property
            const std::shared_ptr< Property > &tPropBedding = get_leader_property(IWG_Property_Type::BEDDING);

            // evaluate the QI
            aQI = 0.5 * trans( tCMElasticity->flux() ) * tCMElasticity->strain();

            // if bedding
            if ( tPropBedding != nullptr )
            {
                // get field interpolator for displacements
                Field_Interpolator *tDisplacementFI = get_leader_fi_manager()->get_field_interpolators_for_type( MSI::Dof_Type::UX );

                // check that field interpolator exists
                MORIS_ASSERT( tDisplacementFI != nullptr,
                        "IQI_Strain_Energy::compute_QI - field interpolator for dof type UX does not exist." );

                // compute body load contribution
                aQI += 0.5 * trans( tDisplacementFI->val() ) * tDisplacementFI->val() * tPropBedding->val();
            }
        }

        //------------------------------------------------------------------------------

        void IQI_Strain_Energy::compute_QI( real aWStar )
        {
            // check the point is inside the bounded box
            if ( !this->is_within_box_bounds() )
            {
                return;
            }

            // get index for QI
            sint tQIIndex = mSet->get_QI_assembly_index( get_name() );

            // get the elasticity CM
           std::shared_ptr< Constitutive_Model > const &tCMElasticity = get_leader_constitutive_model(IQI_Constitutive_Type::ELAST);

            // get bedding property
            const std::shared_ptr< Property > &tPropBedding = get_leader_property(IWG_Property_Type::BEDDING);

            // get thickness property
            const std::shared_ptr< Property > &tPropThickness = get_leader_property(IWG_Property_Type::THICKNESS);

            // multiplying aWStar by user defined thickness (2*pi*r for axisymmetric)
            aWStar *= ( tPropThickness != nullptr ) ? tPropThickness->val()( 0 ) : 1;

            // evaluate the QI
            mSet->get_QI()( tQIIndex ) += aWStar * ( 0.5 * trans( tCMElasticity->flux() ) * tCMElasticity->strain() );

            // if bedding
            if ( tPropBedding != nullptr )
            {
                // get field interpolator for displacements
                Field_Interpolator *tDisplacementFI = get_leader_fi_manager()->get_field_interpolators_for_type( MSI::Dof_Type::UX );

                // check that field interpolator exists
                MORIS_ASSERT( tDisplacementFI != nullptr,
                        "IQI_Strain_Energy::compute_QI - field interpolator for dof type UX does not exist." );

                // compute body load contribution
                mSet->get_QI()( tQIIndex ) += aWStar * ( 0.5 * trans( tDisplacementFI->val() ) * tDisplacementFI->val() * tPropBedding->val() );
            }
        }

        //------------------------------------------------------------------------------

        void IQI_Strain_Energy::compute_dQIdu( real aWStar )
        {
            // check the point is inside the bounded box
            if ( !this->is_within_box_bounds() )
            {
                return;
            }

            // get the column index to assemble in residual
            sint tQIIndex = mSet->get_QI_assembly_index( get_name() );

            // get the elasticity CM
           std::shared_ptr< Constitutive_Model > const &tCMElasticity = get_leader_constitutive_model(IQI_Constitutive_Type::ELAST);

            // get bedding property
            const std::shared_ptr< Property > &tPropBedding = get_leader_property(IWG_Property_Type::BEDDING);

            // get thickness property
            const std::shared_ptr< Property > &tPropThickness = get_leader_property(IWG_Property_Type::THICKNESS);

            // multiplying aWStar by user defined thickness (2*pi*r for axisymmetric)
            aWStar *= ( tPropThickness != nullptr ) ? tPropThickness->val()( 0 ) : 1;

            // get the number of leader dof type dependencies
            uint tNumDofDependencies = get_requested_leader_dof_types().size();

            // compute dQIdu for indirect dof dependencies
            for ( uint iDof = 0; iDof < tNumDofDependencies; iDof++ )
            {
                // get the treated dof type
                Vector< MSI::Dof_Type > const &tDofType = get_requested_leader_dof_types()( iDof );

                // get leader index for residual dof type, indices for assembly
                uint tLeaderDofIndex      = mSet->get_dof_index_for_type( tDofType( 0 ), mtk::Leader_Follower::LEADER );
                uint tLeaderDepStartIndex = mSet->get_res_dof_assembly_map()( tLeaderDofIndex )( 0, 0 );
                uint tLeaderDepStopIndex  = mSet->get_res_dof_assembly_map()( tLeaderDofIndex )( 0, 1 );

                // if elasticity CM depends on dof type
                if ( tCMElasticity->check_dof_dependency( tDofType ) )
                {
                    // compute dQIdu
                    mSet->get_residual()( tQIIndex )(
                            { tLeaderDepStartIndex, tLeaderDepStopIndex } ) += aWStar * 0.5 * ( trans( tCMElasticity->dFluxdDOF( tDofType ) ) * tCMElasticity->strain() + trans( tCMElasticity->dStraindDOF( tDofType ) ) * tCMElasticity->flux() );
                }

                // if bedding
                if ( tPropBedding != nullptr )
                {
                    // get field interpolator for displacements
                    Field_Interpolator *tDisplacementFI = get_leader_fi_manager()->get_field_interpolators_for_type( MSI::Dof_Type::UX );

                    // compute bedding contribution - displacements
                    if ( tDofType( 0 ) == MSI::Dof_Type::UX )
                    {
                        mSet->get_residual()( tQIIndex )(
                                { tLeaderDepStartIndex, tLeaderDepStopIndex } ) += aWStar * ( trans( tDisplacementFI->N() ) * tDisplacementFI->val() * tPropBedding->val() );
                    }

                    if ( tPropBedding->check_dof_dependency( tDofType ) )
                    {
                        mSet->get_residual()( tQIIndex )(
                                { tLeaderDepStartIndex, tLeaderDepStopIndex } ) += aWStar * ( 0.5 * trans( tDisplacementFI->val() ) * tDisplacementFI->val() * tPropBedding->dPropdDOF( tDofType ) );
                    }
                }
            }
        }

        //------------------------------------------------------------------------------

        void IQI_Strain_Energy::compute_dQIdu(
                Vector< MSI::Dof_Type > &aDofType,
                Matrix< DDRMat >        &adQIdu )
        {
            // check the point is inside the bounded box
            if ( !this->is_within_box_bounds() )
            {
                return;
            }

            // get the elasticity CM
           std::shared_ptr< Constitutive_Model > const &tCMElasticity = get_leader_constitutive_model(IQI_Constitutive_Type::ELAST);

            // get bedding property
            const std::shared_ptr< Property > &tPropBedding = get_leader_property(IWG_Property_Type::BEDDING);

            // initialize derivative
            adQIdu.fill( 0.0 );

            // if elasticity CM depends on dof type
            if ( tCMElasticity->check_dof_dependency( aDofType ) )
            {
                // compute dQIdu
                adQIdu = 0.5 * ( trans( tCMElasticity->dFluxdDOF( aDofType ) ) * tCMElasticity->strain() + trans( tCMElasticity->dStraindDOF( aDofType ) ) * tCMElasticity->flux() );
            }

            // if bedding
            if ( tPropBedding != nullptr )
            {
                // get field interpolator for displacements
                Field_Interpolator *tDisplacementFI = get_leader_fi_manager()->get_field_interpolators_for_type( MSI::Dof_Type::UX );

                // compute bedding contribution - displacements
                if ( aDofType( 0 ) == MSI::Dof_Type::UX )
                {
                    adQIdu += trans( tDisplacementFI->N() ) * tDisplacementFI->val() * tPropBedding->val();
                }

                if ( tPropBedding->check_dof_dependency( aDofType ) )
                {
                    adQIdu += 0.5 * trans( tDisplacementFI->val() ) * tDisplacementFI->val() * tPropBedding->dPropdDOF( aDofType );
                }
            }
        }

        //------------------------------------------------------------------------------

        bool IQI_Strain_Energy::is_within_box_bounds()
        {
            Vector< Matrix< DDRMat > > tParameters = get_parameters();

            // check if the box bounds are empty then skip
            if ( tParameters.empty() )
            {
                return true;
            }

            // if the box bounds are not empty then check if it is inside the box
            else
            {
                // get the coordinate
                const Matrix< DDRMat > &tGaussPoint = get_leader_fi_manager()->get_IG_geometry_interpolator()->valx();

                // check if the calculation point coordinates are more then lower corner of the box
                bool tLowerBound = std::equal( tParameters( 0 ).begin(), tParameters( 0 ).end(), tGaussPoint.begin(), tGaussPoint.end(), []( real aA, real aB ) -> bool { return aA < aB; } );

                // check if the calculation point coordinates are less then upper corner of the box
                bool tUpperBound = std::equal( tGaussPoint.begin(), tGaussPoint.end(), tParameters( 1 ).begin(), tParameters( 1 ).end(), []( real aA, real aB ) -> bool { return aA < aB; } );

                // combine the two bounds that satisfy both
                return tUpperBound and tLowerBound;
            }
        }

        //------------------------------------------------------------------------------
    }    // namespace fem
}    // namespace moris
