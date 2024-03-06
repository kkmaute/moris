/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_FEM_IQI.cpp
 *
 */

#include "cl_FEM_IQI.hpp"
#include "cl_FEM_Set.hpp"
#include "cl_FEM_Cluster.hpp"
#include "cl_FEM_Field_Interpolator_Manager.hpp"

#include "fn_norm.hpp"
#include "fn_min.hpp"
#include "fn_max.hpp"

namespace moris
{
    namespace fem
    {

        void
        IQI::set_function_pointers()
        {
            // switch on element type
            switch ( mSet->get_element_type() )
            {
                case fem::Element_Type::BULK:
                {
                    m_compute_dQIdu_FD          = &IQI::select_dQIdu_FD;
                    m_compute_dQIdp_FD_material = &IQI::select_dQIdp_FD_material;
                    m_compute_dQIdp_FD_geometry = &IQI::select_dQIdp_FD_geometry_bulk;
                    break;
                }
                case fem::Element_Type::SIDESET:
                {
                    m_compute_dQIdu_FD          = &IQI::select_dQIdu_FD;
                    m_compute_dQIdp_FD_material = &IQI::select_dQIdp_FD_material;
                    m_compute_dQIdp_FD_geometry = &IQI::select_dQIdp_FD_geometry_sideset;
                    break;
                }
                case fem::Element_Type::TIME_SIDESET:
                case fem::Element_Type::TIME_BOUNDARY:
                {
                    m_compute_dQIdu_FD          = &IQI::select_dQIdu_FD;
                    m_compute_dQIdp_FD_material = &IQI::select_dQIdp_FD_material;
                    m_compute_dQIdp_FD_geometry = &IQI::select_dQIdp_FD_geometry_bulk;
                    break;
                }
                case fem::Element_Type::DOUBLE_SIDESET:
                case fem::Element_Type::NONCONFORMAL_SIDESET:
                {
                    m_compute_dQIdu_FD          = &IQI::select_dQIdu_FD;
                    m_compute_dQIdp_FD_material = &IQI::select_dQIdp_FD_material_double;
                    m_compute_dQIdp_FD_geometry = &IQI::select_dQIdp_FD_geometry_double;
                    break;
                }
                default:
                    MORIS_ERROR( false, "IQI::set_function_pointers - unknown element type." );
            }

            // switch on perturbation strategy
            switch ( mSet->get_perturbation_strategy() )
            {
                case fem::Perturbation_Type::RELATIVE:
                {
                    m_build_perturbation_size = &IQI::build_perturbation_size_relative;
                    break;
                }
                case fem::Perturbation_Type::ABSOLUTE:
                {
                    m_build_perturbation_size = &IQI::build_perturbation_size_absolute;
                    break;
                }
                default:
                    MORIS_ERROR( false, "IQI::set_function_pointers - unknown perturbation type." );
            }
        }


        //------------------------------------------------------------------------------

        void
        IQI::set_reference_value( real aReferenceValue )
        {
            if ( not mNormalized )
            {
                mReferenceValue = aReferenceValue;
                mNormalized     = true;
            }
        }

        //------------------------------------------------------------------------------


        void
        IQI::select_dQIdu_FD(
                real               aWStar,
                real               aPerturbation,
                fem::FDScheme_Type aFDSchemeType )
        {
            // get the column index to assemble in residual
            sint tQIIndex = mSet->get_QI_assembly_index( get_name() );

            // store QI value
            Matrix< DDRMat > tQIStore = mSet->get_QI()( tQIIndex );

            // get the FD scheme info
            Vector< Vector< real > > tFDScheme;
            fd_scheme( aFDSchemeType, tFDScheme );
            uint tNumFDPoints = tFDScheme( 0 ).size();

            // get leader number of dof types
            uint tNumLeaderDofType = get_requested_leader_dof_types().size();

            // reset the QI
            mSet->get_QI()( tQIIndex ).fill( 0.0 );

            // compute the QI
            this->compute_QI( aWStar );

            // store QI value
            Matrix< DDRMat > tQI = mSet->get_QI()( tQIIndex );

            // loop over the IWG dof types
            for ( uint iFI = 0; iFI < tNumLeaderDofType; iFI++ )
            {
                // init dof counter
                uint tDofCounter = 0;

                // get the dof type
Vector< MSI::Dof_Type > const &tDofType = get_requested_leader_dof_types()( iFI );

                // get leader index for residual dof type, indices for assembly
                uint tLeaderDofIndex      = mSet->get_dof_index_for_type( tDofType( 0 ), mtk::Leader_Follower::LEADER );
                uint tLeaderDepStartIndex = mSet->get_res_dof_assembly_map()( tLeaderDofIndex )( 0, 0 );

                // get field interpolator for dependency dof type
                Field_Interpolator* tFI = get_leader_fi_manager()->get_field_interpolators_for_type( tDofType( 0 ) );

                // get number of leader FI bases and fields
                uint tDerNumBases  = tFI->get_number_of_space_time_bases();
                uint tDerNumFields = tFI->get_number_of_fields();

                // coefficients for dof type wrt which derivative is computed
                Matrix< DDRMat > tCoeff = tFI->get_coeff();

                // loop over the coefficient column
                for ( uint iCoeffCol = 0; iCoeffCol < tDerNumFields; iCoeffCol++ )
                {
                    // loop over the coefficient row
                    for ( uint iCoeffRow = 0; iCoeffRow < tDerNumBases; iCoeffRow++ )
                    {
                        // compute the perturbation absolute value
                        real tDeltaH = aPerturbation * tCoeff( iCoeffRow, iCoeffCol );

                        // check that perturbation is not zero
                        if ( std::abs( tDeltaH ) < 1e-12 )
                        {
                            tDeltaH = aPerturbation;
                        }

                        // set starting point for FD
                        uint tStartPoint = 0;

                        // if backward or forward fd
                        if ( ( aFDSchemeType == fem::FDScheme_Type::POINT_1_BACKWARD ) ||    //
                                ( aFDSchemeType == fem::FDScheme_Type::POINT_1_FORWARD ) )
                        {
                            // add unperturbed QI contribution to dQIdu
                            mSet->get_residual()( tQIIndex )(
                                    { tLeaderDepStartIndex + tDofCounter, tLeaderDepStartIndex + tDofCounter },
                                    { 0, 0 } ) +=
                                    tFDScheme( 1 )( 0 ) * tQI /    //
                                    ( tFDScheme( 2 )( 0 ) * tDeltaH );

                            // skip first point in FD
                            tStartPoint = 1;
                        }

                        // loop over the points for FD
                        for ( uint iPoint = tStartPoint; iPoint < tNumFDPoints; iPoint++ )
                        {
                            // reset the perturbed coefficients
                            Matrix< DDRMat > tCoeffPert = tCoeff;

                            // perturb the coefficient
                            tCoeffPert( iCoeffRow, iCoeffCol ) += tFDScheme( 0 )( iPoint ) * tDeltaH;

                            // set the perturbed coefficients to FI
                            tFI->set_coeff( tCoeffPert );
                            tFI->reset_eval_flags();    // not useful

                            // reset properties, CM and SP for IWG
                            this->reset_eval_flags();

                            // reset the QI
                            mSet->get_QI()( tQIIndex ).fill( 0.0 );

                            // compute the QI
                            this->compute_QI( aWStar );

                            // assemble the dQIdu
                            mSet->get_residual()( tQIIndex )(
                                    { tLeaderDepStartIndex + tDofCounter, tLeaderDepStartIndex + tDofCounter },
                                    { 0, 0 } ) +=
                                    tFDScheme( 1 )( iPoint ) *      //
                                    mSet->get_QI()( tQIIndex ) /    //
                                    ( tFDScheme( 2 )( 0 ) * tDeltaH );
                        }
                        // update dof counter
                        tDofCounter++;
                    }
                }
                // reset the coefficients values
                tFI->set_coeff( tCoeff );
            }

            // get follower number of dof types
            uint tNumFollowerDofType = get_requested_follower_dof_types().size();

            // loop over the follower dof types
            for ( uint iFI = 0; iFI < tNumFollowerDofType; iFI++ )
            {
                // init dof counter
                uint tDofCounter = 0;

                // get the dof type
                Vector< MSI::Dof_Type > tDofType = get_requested_follower_dof_types()( iFI );

                // get the index for the dof type
                sint tFollowerDepDofIndex   = mSet->get_dof_index_for_type( tDofType( 0 ), mtk::Leader_Follower::FOLLOWER );
                uint tFollowerDepStartIndex = mSet->get_jac_dof_assembly_map()( tFollowerDepDofIndex )( tFollowerDepDofIndex, 0 );

                // get follower dependency field interpolator
                Field_Interpolator* tFI = get_follower_fi_manager()->get_field_interpolators_for_type( tDofType( 0 ) );

                // get number of leader FI bases and fields
                uint tDerNumBases  = tFI->get_number_of_space_time_bases();
                uint tDerNumFields = tFI->get_number_of_fields();

                // coefficients for dof type wrt which derivative is computed
                Matrix< DDRMat > tCoeff = tFI->get_coeff();

                // loop over the coefficients columns
                for ( uint iCoeffCol = 0; iCoeffCol < tDerNumFields; iCoeffCol++ )
                {
                    // loop over the coefficients rows
                    for ( uint iCoeffRow = 0; iCoeffRow < tDerNumBases; iCoeffRow++ )
                    {
                        // compute the perturbation absolute value
                        real tDeltaH = aPerturbation * tCoeff( iCoeffRow, iCoeffCol );

                        // check that perturbation is not zero
                        if ( std::abs( tDeltaH ) < 1e-12 )
                        {
                            tDeltaH = aPerturbation;
                        }

                        // set starting point for FD
                        uint tStartPoint = 0;

                        // if backward or forward fd
                        if ( ( aFDSchemeType == fem::FDScheme_Type::POINT_1_BACKWARD ) ||    //
                                ( aFDSchemeType == fem::FDScheme_Type::POINT_1_FORWARD ) )
                        {
                            // add unperturbed QI contribution to dQIdu
                            mSet->get_residual()( tQIIndex )(
                                    { tFollowerDepStartIndex + tDofCounter, tFollowerDepStartIndex + tDofCounter },
                                    { 0, 0 } ) +=
                                    tFDScheme( 1 )( 0 ) * tQI /    //
                                    ( tFDScheme( 2 )( 0 ) * tDeltaH );

                            // skip first point in FD
                            tStartPoint = 1;
                        }

                        // loop over the points for FD
                        for ( uint iPoint = tStartPoint; iPoint < tNumFDPoints; iPoint++ )
                        {
                            // reset the perturbed coefficients
                            Matrix< DDRMat > tCoeffPert = tCoeff;

                            // perturb the coefficient
                            tCoeffPert( iCoeffRow, iCoeffCol ) += tFDScheme( 0 )( iPoint ) * tDeltaH;

                            // set the perturbed coefficients to FI
                            tFI->set_coeff( tCoeffPert );

                            // reset properties, CM and SP for IWG
                            this->reset_eval_flags();
                            tFI->reset_eval_flags();    // not useful

                            // reset the QI
                            mSet->get_QI()( tQIIndex ).fill( 0.0 );

                            // compute the QI
                            this->compute_QI( aWStar );

                            // assemble the dQIdu
                            mSet->get_residual()( tQIIndex )(
                                    { tFollowerDepStartIndex + tDofCounter, tFollowerDepStartIndex + tDofCounter },
                                    { 0, 0 } ) +=
                                    tFDScheme( 1 )( iPoint ) *      //
                                    mSet->get_QI()( tQIIndex ) /    //
                                    ( tFDScheme( 2 )( 0 ) * tDeltaH );
                        }
                        // update dof counter
                        tDofCounter++;
                    }
                }
                // reset the coefficients values
                tFI->set_coeff( tCoeff );
            }

            // reset QI value
            mSet->get_QI()( tQIIndex ) = tQIStore;
        }

        //------------------------------------------------------------------------------

        bool
        IQI::check_dQIdu_FD(
                real               aWStar,
                real               aPerturbation,
                real               aEpsilon,
                Matrix< DDRMat >&  adQIdu,
                Matrix< DDRMat >&  adQIduFD,
                bool               aErrorPrint,
                fem::FDScheme_Type aFDSchemeType )
        {
            // get the column index to assemble in residual
            sint tQIIndex = mSet->get_QI_assembly_index( get_name() );

            // compute dQIdu with IQI
            this->compute_dQIdu( aWStar );
            adQIdu = mSet->get_residual()( tQIIndex );

            // reset dQIdu
            mSet->get_residual()( tQIIndex ).fill( 0.0 );

            // compute dQIdu by FD
            this->compute_dQIdu_FD( aWStar, aPerturbation, aFDSchemeType );
            adQIduFD = mSet->get_residual()( tQIIndex );

            // define a boolean for check
            bool tCheckdQIdu = true;

            // check if adQIdu and adQIduFD have the same size
            tCheckdQIdu = tCheckdQIdu && ( adQIdu.n_rows() == adQIduFD.n_rows() );
            tCheckdQIdu = tCheckdQIdu && ( adQIdu.n_cols() == adQIduFD.n_cols() );

            // check that matrices to compare have same size
            MORIS_ERROR(
                    ( adQIdu.n_rows() == adQIduFD.n_rows() ) &&    //
                            ( adQIdu.n_cols() == adQIduFD.n_cols() ),
                    "IQI::check_dQIdu - matrices to check do not share same dimensions." );

            // define a real for absolute difference
            real tAbsolute = 0.0;

            // define a real for relative difference
            real tRelative = 0.0;

            // loop over the rows
            for ( uint iRow = 0; iRow < adQIdu.n_rows(); iRow++ )
            {
                // loop over the columns
                for ( uint jCol = 0; jCol < adQIdu.n_cols(); jCol++ )
                {
                    // get absolute difference
                    tAbsolute = std::abs( adQIduFD( iRow, jCol ) - adQIdu( iRow, jCol ) );

                    // get relative difference
                    tRelative = std::abs( ( adQIduFD( iRow, jCol ) - adQIdu( iRow, jCol ) ) / adQIduFD( iRow, jCol ) );

                    // update check value
                    tCheckdQIdu = tCheckdQIdu && ( ( tAbsolute < aEpsilon ) || ( tRelative < aEpsilon ) );

                    // debug print
                    if ( ( ( tAbsolute < aEpsilon ) || ( tRelative < aEpsilon ) ) == false )
                    {
                        if ( aErrorPrint )
                        {
                            std::cout << "iRow " << iRow << " - jCol " << jCol << "\n"
                                      << std::flush;
                            std::cout << "adQIdu( iRow, jCol )    " << std::setprecision( 12 ) << adQIdu( iRow, jCol ) << "\n"
                                      << std::flush;
                            std::cout << "adQIduFD( iRow, jCol )  " << std::setprecision( 12 ) << adQIduFD( iRow, jCol ) << "\n"
                                      << std::flush;
                            std::cout << "Absolute difference " << tAbsolute << "\n"
                                      << std::flush;
                            std::cout << "Relative difference " << tRelative << "\n"
                                      << std::flush;
                        }
                    }
                }
            }
            // return bool
            return tCheckdQIdu;
        }

        //------------------------------------------------------------------------------

        void
        IQI::select_dQIdp_FD_material(
                moris::real        aWStar,
                moris::real        aPerturbation,
                fem::FDScheme_Type aFDSchemeType )
        {
            // get the column index to assemble in residual
            sint tIQIAssemblyIndex = mSet->get_QI_assembly_index( get_name() );

            // store QI value
            Matrix< DDRMat > tQIStore = mSet->get_QI()( tIQIAssemblyIndex );

            // get the FD scheme info
            Vector< Vector< real > > tFDScheme;
            fd_scheme( aFDSchemeType, tFDScheme );
            uint tNumFDPoints = tFDScheme( 0 ).size();

            // get the requested ip pdv types
            Vector< Vector< gen::PDV_Type > > tRequestedPdvTypes;
            mSet->get_ip_dv_types_for_set( tRequestedPdvTypes );

            // get number of requested dv types
            uint tNumPdvType = tRequestedPdvTypes.size();

            // reset the QI
            mSet->get_QI()( tIQIAssemblyIndex ).fill( 0.0 );

            // compute the QI
            this->compute_QI( aWStar );

            // store QI value
            Matrix< DDRMat > tQI = mSet->get_QI()( tIQIAssemblyIndex );

            // loop over the pdv types
            for ( uint iFI = 0; iFI < tNumPdvType; iFI++ )
            {
                // get the FI for the dv type
                Field_Interpolator* tFI = get_leader_fi_manager()->get_field_interpolators_for_type( tRequestedPdvTypes( iFI )( 0 ) );

                // get number of leader FI bases and fields
                uint tDerNumBases  = tFI->get_number_of_space_time_bases();
                uint tDerNumFields = tFI->get_number_of_fields();

                // coefficients for dof type wrt which derivative is computed
                Matrix< DDRMat > tCoeff = tFI->get_coeff();

                // init pdv coeff counter
                uint tPdvCoeffCounter = 0;

                // loop over the pdv coefficient column
                for ( uint iCoeffCol = 0; iCoeffCol < tDerNumFields; iCoeffCol++ )
                {
                    // loop over the pdv coefficient row
                    for ( uint iCoeffRow = 0; iCoeffRow < tDerNumBases; iCoeffRow++ )
                    {
                        // compute the perturbation absolute value
                        real tDeltaH = aPerturbation * tCoeff( iCoeffRow, iCoeffCol );

                        // check that perturbation is not zero
                        if ( std::abs( tDeltaH ) < 1e-12 )
                        {
                            tDeltaH = aPerturbation;
                        }

                        // get the pdv index for assembly
                        uint tPdvAssemblyIndex = mSet->get_mat_pdv_assembly_map()( iFI )( 0, 0 ) + tPdvCoeffCounter;

                        // set starting point for FD
                        uint tStartPoint = 0;

                        // if backward or forward fd
                        if ( ( aFDSchemeType == fem::FDScheme_Type::POINT_1_BACKWARD ) ||    //
                                ( aFDSchemeType == fem::FDScheme_Type::POINT_1_FORWARD ) )
                        {
                            // add unperturbed QI contribution to dQIdp
                            mSet->get_dqidpmat()( tIQIAssemblyIndex )( tPdvAssemblyIndex ) +=
                                    tFDScheme( 1 )( 0 ) * tQI( 0 ) / ( tFDScheme( 2 )( 0 ) * tDeltaH );

                            // skip first point in FD
                            tStartPoint = 1;
                        }

                        // loop over the points for FD
                        for ( uint iPoint = tStartPoint; iPoint < tNumFDPoints; iPoint++ )
                        {
                            // reset the perturbed coefficients
                            Matrix< DDRMat > tCoeffPert = tCoeff;

                            // perturb the coefficient
                            tCoeffPert( iCoeffRow, iCoeffCol ) += tFDScheme( 0 )( iPoint ) * tDeltaH;

                            // set the perturbed coefficients to FI
                            tFI->set_coeff( tCoeffPert );

                            // reset properties, CM and SP for IWG
                            this->reset_eval_flags();

                            // reset the QI
                            mSet->get_QI()( tIQIAssemblyIndex ).fill( 0.0 );

                            // compute the QI
                            this->compute_QI( aWStar );

                            // assemble the jacobian
                            mSet->get_dqidpmat()( tIQIAssemblyIndex )( tPdvAssemblyIndex ) +=
                                    tFDScheme( 1 )( iPoint ) *                    //
                                    mSet->get_QI()( tIQIAssemblyIndex )( 0 ) /    //
                                    ( tFDScheme( 2 )( 0 ) * tDeltaH );
                        }
                        // update the pdv coefficient counter
                        tPdvCoeffCounter++;
                    }
                }
                // reset the coefficients values
                tFI->set_coeff( tCoeff );
            }

            // reset QI value
            mSet->get_QI()( tIQIAssemblyIndex ) = tQIStore;

            // check for nan, infinity
            MORIS_ASSERT( isfinite( mSet->get_dqidpmat()( tIQIAssemblyIndex ) ),
                    "IQI::compute_dQIdp_FD_material - dQIdp contains NAN or INF, exiting!" );
        }

        //------------------------------------------------------------------------------

        real
        IQI::build_perturbation_size(
                const real& aPerturbation,
                const real& aCoefficientToPerturb,
                const real& aMaxPerturbation,
                const real& aTolerance )
        {
            return ( this->*m_build_perturbation_size )(
                    aPerturbation,
                    aCoefficientToPerturb,
                    aMaxPerturbation,
                    aTolerance );
        }

        //------------------------------------------------------------------------------

        real
        IQI::build_perturbation_size_relative(
                const real& aPerturbation,
                const real& aCoefficientToPerturb,
                const real& aMaxPerturbation,
                const real& aTolerance )
        {
            // check that maximum perturbation size is larger than tolerance
            MORIS_ASSERT( aMaxPerturbation >= aTolerance,
                    "IQI::build_perturbation_size_relative - maximum perturbation size is smaller than tolerance: max = %e  tol = %e\n",
                    aMaxPerturbation,
                    aTolerance );

            // determine actual tolerance (only useful when above assert inactive)
            const real tActualTol = std::min( aMaxPerturbation, aTolerance );

            // compute the perturbation value using fraction of maximum allowable perturbation
            real tDeltaH = std::abs( aPerturbation * aMaxPerturbation );

            // compute perturbation such that it is not smaller than tolerance
            // and not larger than maximum value
            tDeltaH = std::max( std::min( tDeltaH, aMaxPerturbation ), tActualTol );

            return tDeltaH;
        }

        //------------------------------------------------------------------------------

        real
        IQI::build_perturbation_size_absolute(
                const real& aPerturbation,
                const real& aCoefficientToPerturb,
                const real& aMaxPerturbation,
                const real& aTolerance )
        {
            // check that maximum perturbation size is larger than tolerance
            MORIS_ASSERT( aMaxPerturbation >= aTolerance,
                    "IQI::build_perturbation_size_absolute - maximum perturbation size is smaller than tolerance: max = %e  tol = %e\n",
                    aMaxPerturbation,
                    aTolerance );

            // determine actual tolerance (only useful when above assert inactive)
            const real tActualTol = std::min( aMaxPerturbation, aTolerance );

            // check that absolute value of perturbation is not smaller than tolerance
            // and not larger than maximum value
            return std::max( std::min( std::abs( aPerturbation ), aMaxPerturbation ), tActualTol );
        }

        //------------------------------------------------------------------------------

        real
        IQI::check_ig_coordinates_inside_ip_element(
                const real&         aPerturbation,
                const real&         aCoefficientToPerturb,
                const uint&         aSpatialDirection,
                fem::FDScheme_Type& aUsedFDScheme )
        {
            // FIXME: only works for rectangular IP elements
            // FIXME: only works for forward, backward, central, not for higher as 5-point FD

            // get the IP element geometry interpolator
            Geometry_Interpolator* tIPGI =
                    mSet->get_field_interpolator_manager()->get_IP_geometry_interpolator();

            // IP element max/min
            real tMaxIP = max( tIPGI->get_space_coeff().get_column( aSpatialDirection ) );    // get maximum values of coordinates of IP nodes
            real tMinIP = min( tIPGI->get_space_coeff().get_column( aSpatialDirection ) );    // get minimum values of coordinates of IP nodes

            // get maximum possible perturbation
            real tMaxPerturb = ( tMaxIP - tMinIP ) / 3.0;

            // compute the perturbation value
            real tDeltaH = build_perturbation_size(
                    aPerturbation,
                    aCoefficientToPerturb,
                    tMaxPerturb,
                    mToleranceFD );

            // check that IG node coordinate is consistent with minimum and maximum IP coordinates
            MORIS_ASSERT(
                    tMaxIP >= aCoefficientToPerturb - tDeltaH &&    //
                            tMinIP <= aCoefficientToPerturb + tDeltaH,
                    "ERROR: IG coordinates are outside IP element: dim: %d  minIP: %e  maxIP: %e  cordIG: %e  \n",
                    aSpatialDirection,
                    tMinIP,
                    tMaxIP,
                    aCoefficientToPerturb );

            // check point location
            if ( aCoefficientToPerturb + tDeltaH >= tMaxIP )
            {
                aUsedFDScheme = fem::FDScheme_Type::POINT_1_BACKWARD;

                // check for correctness of perturbation size for backward FD
                MORIS_ASSERT( tDeltaH < aCoefficientToPerturb - tMinIP,
                        "ERROR: backward perturbation size exceed limits of interpolation element:\n"
                        "dim: %d  minIP: %e  maxIP: %e  cordIG: %e  maxPert: %e  delta: %e  precPert: %e\n.",
                        aSpatialDirection,
                        tMinIP,
                        tMaxIP,
                        aCoefficientToPerturb,
                        tMaxPerturb,
                        tDeltaH,
                        aPerturbation );
            }
            else
            {
                if ( aCoefficientToPerturb - tDeltaH <= tMinIP )
                {
                    aUsedFDScheme = fem::FDScheme_Type::POINT_1_FORWARD;

                    // check for correctness of perturbation size for backward FD
                    MORIS_ASSERT( tDeltaH < tMaxIP - aCoefficientToPerturb,
                            "ERROR: forward perturbation size exceeds limits of interpolation element: dim: %d  minIP: %e  maxIP: %e  cordIG: %e  maxPert: %e  delta: %e  precPert: %e\n.",
                            aSpatialDirection,
                            tMinIP,
                            tMaxIP,
                            aCoefficientToPerturb,
                            tMaxPerturb,
                            tDeltaH,
                            aPerturbation );
                }
                else
                {
                    // check for correctness of perturbation size for central FD
                    MORIS_ASSERT(
                            tDeltaH < tMaxIP - aCoefficientToPerturb &&    //
                                    tDeltaH < aCoefficientToPerturb - tMinIP,
                            "ERROR: central perturbation size exceed limits of interpolation element: dim: %d  minIP: %e  maxIP: %e  cordIG: %e  maxPert: %e  delta: %e  precPert: %e\n.",
                            aSpatialDirection,
                            tMinIP,
                            tMaxIP,
                            aCoefficientToPerturb,
                            tMaxPerturb,
                            tDeltaH,
                            aPerturbation );
                }
            }

            return tDeltaH;
        }

        //------------------------------------------------------------------------------

        void
        IQI::select_dQIdp_FD_geometry_bulk(
                moris::real                   aWStar,
                moris::real                   aPerturbation,
                fem::FDScheme_Type            aFDSchemeType,
                Matrix< DDSMat >&             aGeoLocalAssembly,
                Vector< Matrix< IndexMat > >& aVertexIndices )
        {
            // get the IQI index
            uint tIQIAssemblyIndex = mSet->get_QI_assembly_index( get_name() );

            // store QI value
            Matrix< DDRMat > tQIStore = mSet->get_QI()( tIQIAssemblyIndex );

            // get the GI for the IP and IG element considered
            Geometry_Interpolator* tIPGI =
                    mSet->get_field_interpolator_manager()->get_IP_geometry_interpolator();
            Geometry_Interpolator* tIGGI =
                    mSet->get_field_interpolator_manager()->get_IG_geometry_interpolator();

            // get number of leader GI bases and space dimensions
            uint tDerNumBases      = tIGGI->get_number_of_space_bases();
            uint tDerNumDimensions = tIPGI->get_number_of_space_dimensions();

            // coefficients for dv type wrt which derivative is computed
            Matrix< DDRMat > tCoeff      = tIGGI->get_space_coeff();
            Matrix< DDRMat > tParamCoeff = tIGGI->get_space_param_coeff();
            Matrix< DDRMat > tEvaluationPoint;
            tIGGI->get_space_time( tEvaluationPoint );
            real tGPWeight = aWStar / tIGGI->det_J();

            // reset the QI
            mSet->get_QI()( tIQIAssemblyIndex ).fill( 0.0 );

            // compute the QI
            this->compute_QI( aWStar );

            // store QI value
            Matrix< DDRMat > tQI = mSet->get_QI()( tIQIAssemblyIndex );

            // init FD scheme
            Vector< Vector< real > > tFDScheme;

            // init perturbation
            real tDeltaH = 0.0;

            // loop over the spatial directions/loop on pdv type
            for ( uint iCoeffCol = 0; iCoeffCol < tDerNumDimensions; iCoeffCol++ )
            {
                // loop over the IG nodes/loop one nodes
                for ( uint iCoeffRow = 0; iCoeffRow < tDerNumBases; iCoeffRow++ )
                {
                    // get the geometry pdv assembly index
                    sint tPdvAssemblyIndex = aGeoLocalAssembly( iCoeffRow, iCoeffCol );

                    // if pdv is active
                    if ( tPdvAssemblyIndex != -1 )
                    {
                        // check point location and define perturbation size and FD scheme accordingly
                        fem::FDScheme_Type tUsedFDScheme = aFDSchemeType;
                        tDeltaH                          = this->check_ig_coordinates_inside_ip_element(
                                aPerturbation,
                                tCoeff( iCoeffRow, iCoeffCol ),
                                iCoeffCol,
                                tUsedFDScheme );

                        // finalize FD scheme
                        fd_scheme( tUsedFDScheme, tFDScheme );
                        uint tNumFDPoints = tFDScheme( 0 ).size();

                        // set starting point for FD
                        uint tStartPoint = 0;

                        // if backward or forward fd
                        if ( ( tUsedFDScheme == fem::FDScheme_Type::POINT_1_BACKWARD ) ||    //
                                ( tUsedFDScheme == fem::FDScheme_Type::POINT_1_FORWARD ) )
                        {
                            // add unperturbed QI contribution to dQIdp
                            mSet->get_dqidpgeo()( tIQIAssemblyIndex )( tPdvAssemblyIndex ) +=
                                    tFDScheme( 1 )( 0 ) * tQI( 0 ) / ( tFDScheme( 2 )( 0 ) * tDeltaH );

                            // skip first point in FD
                            tStartPoint = 1;
                        }

                        // loop over point of FD scheme
                        for ( uint iPoint = tStartPoint; iPoint < tNumFDPoints; iPoint++ )
                        {
                            // reset the perturbed coefficients
                            Matrix< DDRMat > tCoeffPert = tCoeff;

                            // perturb the coefficient
                            tCoeffPert( iCoeffRow, iCoeffCol ) += tFDScheme( 0 )( iPoint ) * tDeltaH;

                            // setting the perturbed coefficients
                            tIGGI->set_space_coeff( tCoeffPert );

                            // update local coordinates
                            Matrix< DDRMat > tXCoords  = tCoeffPert.get_row( iCoeffRow );
                            Matrix< DDRMat > tXiCoords = tParamCoeff.get_row( iCoeffRow );

                            tIPGI->update_parametric_coordinates( tXCoords, tXiCoords );

                            Matrix< DDRMat > tParamCoeffPert     = tParamCoeff;
                            tParamCoeffPert.get_row( iCoeffRow ) = tXiCoords.matrix_data();

                            tIGGI->set_space_param_coeff( tParamCoeffPert );

                            // set evaluation point for interpolators (FIs and GIs)
                            mSet->get_field_interpolator_manager()->    //
                                    set_space_time_from_local_IG_point( tEvaluationPoint );

                            // reset properties, CM and SP for IWG
                            this->reset_eval_flags();

                            // reset the QI
                            mSet->get_QI()( tIQIAssemblyIndex ).fill( 0.0 );

                            // compute the QI
                            real tWStarPert = tGPWeight * tIGGI->det_J();
                            this->compute_QI( tWStarPert );

                            // evaluate dQIdpGeo
                            mSet->get_dqidpgeo()( tIQIAssemblyIndex )( tPdvAssemblyIndex ) +=
                                    tFDScheme( 1 )( iPoint ) *                    //
                                    mSet->get_QI()( tIQIAssemblyIndex )( 0 ) /    //
                                    ( tFDScheme( 2 )( 0 ) * tDeltaH );
                        }
                    }
                }
                // reset the coefficients values
                tIGGI->set_space_coeff( tCoeff );
                tIGGI->set_space_param_coeff( tParamCoeff );
                mSet->get_field_interpolator_manager()->set_space_time_from_local_IG_point( tEvaluationPoint );
            }

            // reset QI value
            mSet->get_QI()( tIQIAssemblyIndex ) = tQIStore;

            // if active cluster measure on IQI
            if ( is_active_cluster_measure() )
            {
                // add their contribution to dQIdp
                this->add_cluster_measure_dQIdp_FD_geometry(
                        aWStar,
                        aPerturbation,
                        aFDSchemeType );
            }

            // check for nan, infinity
            MORIS_ASSERT( isfinite( mSet->get_dqidpgeo()( tIQIAssemblyIndex ) ),
                    "IQI::compute_dQIdp_FD_geometry - dQIdp contains NAN or INF, exiting!" );
        }

        //------------------------------------------------------------------------------

        void
        IQI::select_dQIdp_FD_geometry_sideset(
                moris::real                   aWStar,
                moris::real                   aPerturbation,
                fem::FDScheme_Type            aFDSchemeType,
                Matrix< DDSMat >&             aGeoLocalAssembly,
                Vector< Matrix< IndexMat > >& aVertexIndices )
        {
            // get the IQI index
            uint tIQIAssemblyIndex = mSet->get_QI_assembly_index( get_name() );

            // store QI value
            Matrix< DDRMat > tQIStore = mSet->get_QI()( tIQIAssemblyIndex );

            // get the GI for the IP and IG element considered
            Geometry_Interpolator* tIPGI =
                    mSet->get_field_interpolator_manager()->get_IP_geometry_interpolator();
            Geometry_Interpolator* tIGGI =
                    mSet->get_field_interpolator_manager()->get_IG_geometry_interpolator();

            // store unperturbed xyz coordinates
            Matrix< DDRMat > tCoeff = tIGGI->get_space_coeff();

            // store unperturbed local coordinates
            Matrix< DDRMat > tParamCoeff = tIGGI->get_space_param_coeff();

            // store unperturbed evaluation point
            Matrix< DDRMat > tEvaluationPoint;
            tIGGI->get_space_time( tEvaluationPoint );

            // store unperturbed evaluation point weight
            real tGPWeight = aWStar / tIGGI->det_J();

            // store unperturbed normal
            Matrix< DDRMat > tNormal;
            tIGGI->get_normal( tNormal );

            // reset the QI
            mSet->get_QI()( tIQIAssemblyIndex ).fill( 0.0 );

            // compute the QI
            this->compute_QI( aWStar );

            // store QI value
            Matrix< DDRMat > tQI = mSet->get_QI()( tIQIAssemblyIndex );

            // init FD scheme
            Vector< Vector< real > > tFDScheme;

            // init perturbation
            real tDeltaH = 0.0;

            // get number of leader GI bases and space dimensions
            uint tDerNumBases      = tIGGI->get_number_of_space_bases();
            uint tDerNumDimensions = tIPGI->get_number_of_space_dimensions();

            // loop over the spatial directions/loop on pdv type
            for ( uint iCoeffCol = 0; iCoeffCol < tDerNumDimensions; iCoeffCol++ )
            {
                // loop over the IG nodes/loop one nodes
                for ( uint iCoeffRow = 0; iCoeffRow < tDerNumBases; iCoeffRow++ )
                {
                    // get the geometry pdv assembly index
                    sint tPdvAssemblyIndex = aGeoLocalAssembly( iCoeffRow, iCoeffCol );

                    // if pdv is active
                    if ( tPdvAssemblyIndex != -1 )
                    {
                        // check point location and define perturbation size and FD scheme accordingly
                        fem::FDScheme_Type tUsedFDScheme = aFDSchemeType;
                        tDeltaH                          = this->check_ig_coordinates_inside_ip_element(
                                aPerturbation,
                                tCoeff( iCoeffRow, iCoeffCol ),
                                iCoeffCol,
                                tUsedFDScheme );

                        // finalize FD scheme
                        fd_scheme( tUsedFDScheme, tFDScheme );
                        uint tNumFDPoints = tFDScheme( 0 ).size();

                        // set starting point for FD
                        uint tStartPoint = 0;

                        // if backward or forward fd
                        if ( ( tUsedFDScheme == fem::FDScheme_Type::POINT_1_BACKWARD ) ||    //
                                ( tUsedFDScheme == fem::FDScheme_Type::POINT_1_FORWARD ) )
                        {
                            // add unperturbed QI contribution to dQIdp
                            mSet->get_dqidpgeo()( tIQIAssemblyIndex )( tPdvAssemblyIndex ) +=
                                    tFDScheme( 1 )( 0 ) * tQI( 0 ) / ( tFDScheme( 2 )( 0 ) * tDeltaH );

                            // skip first point in FD
                            tStartPoint = 1;
                        }

                        // loop over point of FD scheme
                        for ( uint iPoint = tStartPoint; iPoint < tNumFDPoints; iPoint++ )
                        {
                            // reset the perturbed coefficients
                            Matrix< DDRMat > tCoeffPert = tCoeff;

                            // perturb the coefficient
                            tCoeffPert( iCoeffRow, iCoeffCol ) += tFDScheme( 0 )( iPoint ) * tDeltaH;

                            // setting the perturbed coefficients
                            tIGGI->set_space_coeff( tCoeffPert );

                            // update local coordinates
                            Matrix< DDRMat > tXCoords  = tCoeffPert.get_row( iCoeffRow );
                            Matrix< DDRMat > tXiCoords = tParamCoeff.get_row( iCoeffRow );

                            tIPGI->update_parametric_coordinates( tXCoords, tXiCoords );

                            Matrix< DDRMat > tParamCoeffPert = tParamCoeff;

                            tParamCoeffPert.get_row( iCoeffRow ) = tXiCoords.matrix_data();

                            tIGGI->set_space_param_coeff( tParamCoeffPert );

                            // set evaluation point for interpolators (FIs and GIs)
                            mSet->get_field_interpolator_manager()->    //
                                    set_space_time_from_local_IG_point( tEvaluationPoint );

                            // reset the normal
                            Matrix< DDRMat > tNormalPert;
                            tIGGI->get_normal( tNormalPert );
                            this->set_normal( tNormalPert );

                            // reset properties, CM and SP for IWG
                            this->reset_eval_flags();

                            // reset the QI
                            mSet->get_QI()( tIQIAssemblyIndex ).fill( 0.0 );

                            // compute the QI
                            real tWStarPert = tGPWeight * tIGGI->det_J();
                            this->compute_QI( tWStarPert );

                            // evaluate dQIdpGeo
                            mSet->get_dqidpgeo()( tIQIAssemblyIndex )( tPdvAssemblyIndex ) +=
                                    tFDScheme( 1 )( iPoint ) * mSet->get_QI()( tIQIAssemblyIndex )( 0 ) /    //
                                    ( tFDScheme( 2 )( 0 ) * tDeltaH );
                        }
                    }
                }
                // reset xyz values
                tIGGI->set_space_coeff( tCoeff );

                // reset local coordinates
                tIGGI->set_space_param_coeff( tParamCoeff );

                // reset evaluation point
                mSet->get_field_interpolator_manager()->set_space_time_from_local_IG_point( tEvaluationPoint );

                // reset normal
                this->set_normal( tNormal );
            }

            // reset QI value
            mSet->get_QI()( tIQIAssemblyIndex ) = tQIStore;

            // if active cluster measure on IQI
            if ( is_active_cluster_measure() )
            {
                // add their contribution to dQIdp
                this->add_cluster_measure_dQIdp_FD_geometry(
                        aWStar,
                        aPerturbation,
                        aFDSchemeType );
            }

            // check for nan, infinity
            MORIS_ASSERT( isfinite( mSet->get_dqidpgeo()( tIQIAssemblyIndex ) ),
                    "IQI::compute_dQIdp_FD_geometry - dQIdp contains NAN or INF, exiting!" );
        }

        //------------------------------------------------------------------------------

        void
        IQI::add_cluster_measure_dQIdp_FD_geometry(
                moris::real        aWStar,
                moris::real        aPerturbation,
                fem::FDScheme_Type aFDSchemeType )
        {
            // get the IQI index
            uint tIQIAssemblyIndex = mSet->get_QI_assembly_index( get_name() );

            // store QI value
            Matrix< DDRMat > tQIStore = mSet->get_QI()( tIQIAssemblyIndex );

            // reset the QI
            mSet->get_QI()( tIQIAssemblyIndex ).fill( 0.0 );

            // compute the QI
            this->compute_QI( aWStar );

            // store QI value
            Matrix< DDRMat > tQI = mSet->get_QI()( tIQIAssemblyIndex );

            // init FD scheme
            Vector< Vector< real > > tFDScheme;

            // initialize perturbation of cluster measure
            real tDeltaCM = 0.0;

            // loop over the cluster measures
            for ( auto const& [ _, tClusterMeasure ] : mCluster->get_cluster_measures() )
            {
                // evaluate the perturbation of cluster measure
                tDeltaCM = this->build_perturbation_size(
                        aPerturbation,
                        tClusterMeasure->val()( 0 ),
                        std::max( tClusterMeasure->val()( 0 ), mToleranceFD ),
                        mToleranceFD );

                // create FD scheme
                fd_scheme( aFDSchemeType, tFDScheme );
                uint tNumFDPoints = tFDScheme( 0 ).size();

                // set starting point for FD
                uint tStartPoint = 0;

                // if backward or forward fd
                if ( ( aFDSchemeType == fem::FDScheme_Type::POINT_1_BACKWARD ) ||    //
                        ( aFDSchemeType == fem::FDScheme_Type::POINT_1_FORWARD ) )
                {
                    // add unperturbed QI contribution to dQIdp
                    mSet->get_dqidpgeo()( tIQIAssemblyIndex ) +=
                            tFDScheme( 1 )( 0 ) * tQI( 0 ) *    //
                            tClusterMeasure->dMEAdPDV() / ( tFDScheme( 2 )( 0 ) * tDeltaCM );

                    // skip first point in FD
                    tStartPoint = 1;
                }

                // loop over point of FD scheme
                for ( uint iPoint = tStartPoint; iPoint < tNumFDPoints; iPoint++ )
                {
                    // perturb the cluster measure
                    tClusterMeasure->perturb_cluster_measure( tFDScheme( 0 )( iPoint ) * tDeltaCM );

                    // reset properties, CM and SP for IWG
                    this->reset_eval_flags();

                    // reset the QI
                    mSet->get_QI()( tIQIAssemblyIndex ).fill( 0.0 );

                    // compute the QI
                    this->compute_QI( aWStar );

                    // evaluate dQIdpGeo
                    mSet->get_dqidpgeo()( tIQIAssemblyIndex ) +=
                            tFDScheme( 1 )( iPoint ) *                                                  //
                            mSet->get_QI()( tIQIAssemblyIndex )( 0 ) * tClusterMeasure->dMEAdPDV() /    //
                            ( tFDScheme( 2 )( 0 ) * tDeltaCM );

                    // reset cluster measures
                    mCluster->reset_cluster_measure();
                }
            }

            // reset QI value
            mSet->get_QI()( tIQIAssemblyIndex ) = tQIStore;

            // check for nan, infinity
            MORIS_ASSERT( isfinite( mSet->get_dqidpgeo()( tIQIAssemblyIndex ) ),
                    "IQI::compute_dQIdp_FD_geometry - dQIdp contains NAN or INF, exiting!" );
        }

        //------------------------------------------------------------------------------
    }    // namespace fem
}    // namespace moris
