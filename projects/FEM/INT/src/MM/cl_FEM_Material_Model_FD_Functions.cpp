/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_FEM_Material_Model_FD_Functions.cpp
 *
 */

#include "cl_FEM_Material_Model.hpp"
#include "cl_FEM_Set.hpp"
#include "cl_FEM_Field_Interpolator_Manager.hpp"

namespace moris
{
    namespace fem
    {
        //------------------------------------------------------------------------------
        // FINITE DIFFERENCE FUNCTIONS
        //------------------------------------------------------------------------------

        void Material_Model::eval_EintDOF_FD(
                const Vector< MSI::Dof_Type > & aDofTypes,
                Matrix< DDRMat >                   & aEintDOF_FD,
                real                                 aPerturbation,
                fem::FDScheme_Type                   aFDSchemeType)
        {
            // get the FD scheme info
            Vector< Vector< real > > tFDScheme;
            fd_scheme( aFDSchemeType, tFDScheme );
            uint tNumPoints = tFDScheme( 0 ).size();

            // get the field interpolator for type
            Field_Interpolator* tFI = mFIManager->get_field_interpolators_for_type( aDofTypes( 0 ) );

            // get number of coefficients, fields and bases for the considered FI
            uint tDerNumDof    = tFI->get_number_of_space_time_coefficients();
            uint tDerNumBases  = tFI->get_number_of_space_time_bases();
            uint tDerNumFields = tFI->get_number_of_fields();

            // evaluate unperturbed internal energy
            Matrix< DDRMat > tEint = this->Eint();

            // set size for derivative
            aEintDOF_FD.set_size( tEint.n_rows(), tDerNumDof, 0.0 );

            // coefficients for dof type wrt which derivative is computed
            Matrix< DDRMat > tCoeff = tFI->get_coeff();

            // initialize dof counter
            uint tDofCounter = 0;

            // loop over coefficients columns
            for( uint iCoeffCol = 0; iCoeffCol < tDerNumFields; iCoeffCol++ )
            {
                // loop over coefficients rows
                for( uint iCoeffRow = 0; iCoeffRow < tDerNumBases; iCoeffRow++ )
                {
                    // compute the perturbation absolute value
                    real tDeltaH = aPerturbation * tCoeff( iCoeffRow, iCoeffCol );

                    // check that perturbation is not zero
                    if( std::abs( tDeltaH ) < 1e-12 )
                    {
                        tDeltaH = aPerturbation;
                    }

                    // set starting point for FD
                    uint tStartPoint = 0;

                    // if backward or forward add unperturbed contribution
                    if( ( aFDSchemeType == fem::FDScheme_Type::POINT_1_BACKWARD ) ||
                            ( aFDSchemeType == fem::FDScheme_Type::POINT_1_FORWARD ) )
                    {
                        // add unperturbed flux contribution to dfluxdu
                        aEintDOF_FD.get_column( tDofCounter ) +=
                                tFDScheme( 1 )( 0 ) * tEint /
                                ( tFDScheme( 2 )( 0 ) * tDeltaH );

                        // skip first point in FD
                        tStartPoint = 1;
                    }

                    // loop over the points for FD
                    for( uint iPoint = tStartPoint; iPoint < tNumPoints; iPoint++ )
                    {
                        // reset the perturbed coefficients
                        Matrix< DDRMat > tCoeffPert = tCoeff;

                        // perturb the coefficient
                        tCoeffPert( iCoeffRow, iCoeffCol ) += tFDScheme( 0 )( iPoint ) * tDeltaH;

                        // set the perturbed coefficients to FI
                        tFI->set_coeff( tCoeffPert );

                        // reset properties
                        this->reset_eval_flags();

                        // assemble the jacobian
                        aEintDOF_FD.get_column( tDofCounter ) +=
                                        tFDScheme( 1 )( iPoint ) * this->Eint() /
                                        ( tFDScheme( 2 )( 0 ) * tDeltaH );
                    }
                    // update dof counter
                    tDofCounter++;
                }
            }
            // reset the coefficients values
            tFI->set_coeff( tCoeff );
        }

        //------------------------------------------------------------------------------

        void Material_Model::eval_EintDotDOF_FD(
                const Vector< MSI::Dof_Type > & aDofTypes,
                Matrix< DDRMat >                   & aEintDotDOF_FD,
                real                                 aPerturbation,
                fem::FDScheme_Type                   aFDSchemeType)
        {
            // get the FD scheme info
            Vector< Vector< real > > tFDScheme;
            fd_scheme( aFDSchemeType, tFDScheme );
            uint tNumPoints = tFDScheme( 0 ).size();

            // get the field interpolator for type
            Field_Interpolator* tFI = mFIManager->get_field_interpolators_for_type( aDofTypes( 0 ) );

            // get number of coefficients, fields and bases for the considered FI
            uint tDerNumDof    = tFI->get_number_of_space_time_coefficients();
            uint tDerNumBases  = tFI->get_number_of_space_time_bases();
            uint tDerNumFields = tFI->get_number_of_fields();

            // evaluate unperturbed internal energy
            Matrix< DDRMat > tEintDot = this->EintDot();

            // set size for derivative
            aEintDotDOF_FD.set_size( tEintDot.n_rows(), tDerNumDof, 0.0 );

            // coefficients for dof type wrt which derivative is computed
            Matrix< DDRMat > tCoeff = tFI->get_coeff();

            // initialize dof counter
            uint tDofCounter = 0;

            // loop over coefficients columns
            for( uint iCoeffCol = 0; iCoeffCol < tDerNumFields; iCoeffCol++ )
            {
                // loop over coefficients rows
                for( uint iCoeffRow = 0; iCoeffRow < tDerNumBases; iCoeffRow++ )
                {
                    // compute the perturbation absolute value
                    real tDeltaH = aPerturbation * tCoeff( iCoeffRow, iCoeffCol );

                    // check that perturbation is not zero
                    if( std::abs( tDeltaH ) < 1e-12 )
                    {
                        tDeltaH = aPerturbation;
                    }

                    // set starting point for FD
                    uint tStartPoint = 0;

                    // if backward or forward add unperturbed contribution
                    if( ( aFDSchemeType == fem::FDScheme_Type::POINT_1_BACKWARD ) ||
                            ( aFDSchemeType == fem::FDScheme_Type::POINT_1_FORWARD ) )
                    {
                        // add unperturbed flux contribution to dfluxdu
                        aEintDotDOF_FD.get_column( tDofCounter ) +=
                                tFDScheme( 1 )( 0 ) * tEintDot /
                                ( tFDScheme( 2 )( 0 ) * tDeltaH );

                        // skip first point in FD
                        tStartPoint = 1;
                    }

                    // loop over the points for FD
                    for( uint iPoint = tStartPoint; iPoint < tNumPoints; iPoint++ )
                    {
                        // reset the perturbed coefficients
                        Matrix< DDRMat > tCoeffPert = tCoeff;

                        // perturb the coefficient
                        tCoeffPert( iCoeffRow, iCoeffCol ) += tFDScheme( 0 )( iPoint ) * tDeltaH;

                        // set the perturbed coefficients to FI
                        tFI->set_coeff( tCoeffPert );

                        // reset properties
                        this->reset_eval_flags();

                        // assemble the jacobian
                        aEintDotDOF_FD.get_column( tDofCounter ) +=
                                        tFDScheme( 1 )( iPoint ) * this->EintDot() /
                                        ( tFDScheme( 2 )( 0 ) * tDeltaH );
                    }
                    // update dof counter
                    tDofCounter++;
                }
            }
            // reset the coefficients values
            tFI->set_coeff( tCoeff );
        }

        //------------------------------------------------------------------------------

        void Material_Model::eval_dnEintdxnDOF_FD(
                const Vector< MSI::Dof_Type > & aDofTypes,
                Matrix< DDRMat >                   & adnEintdxnDOF_FD,
                real                                 aPerturbation,
                uint                                 aOrder,
                fem::FDScheme_Type                   aFDSchemeType)
        {
            // get the FD scheme info
            Vector< Vector< real > > tFDScheme;
            fd_scheme( aFDSchemeType, tFDScheme );
            uint tNumPoints = tFDScheme( 0 ).size();

            // get the field interpolator for type
            Field_Interpolator* tFI = mFIManager->get_field_interpolators_for_type( aDofTypes( 0 ) );

            // get number of coefficients, fields and bases for the considered FI
            uint tDerNumDof    = tFI->get_number_of_space_time_coefficients();
            uint tDerNumBases  = tFI->get_number_of_space_time_bases();
            uint tDerNumFields = tFI->get_number_of_fields();

            // evaluate unperturbed internal energy
            Matrix< DDRMat > tdnEintdxn = this->dnEintdxn( aOrder );

            // set size for derivative
            adnEintdxnDOF_FD.set_size( tdnEintdxn.n_rows(), tDerNumDof, 0.0 );

            // coefficients for dof type wrt which derivative is computed
            Matrix< DDRMat > tCoeff = tFI->get_coeff();

            // initialize dof counter
            uint tDofCounter = 0;

            // loop over coefficients columns
            for( uint iCoeffCol = 0; iCoeffCol < tDerNumFields; iCoeffCol++ )
            {
                // loop over coefficients rows
                for( uint iCoeffRow = 0; iCoeffRow < tDerNumBases; iCoeffRow++ )
                {
                    // compute the perturbation absolute value
                    real tDeltaH = aPerturbation * tCoeff( iCoeffRow, iCoeffCol );

                    // check that perturbation is not zero
                    if( std::abs( tDeltaH ) < 1e-12 )
                    {
                        tDeltaH = aPerturbation;
                    }

                    // set starting point for FD
                    uint tStartPoint = 0;

                    // if backward or forward add unperturbed contribution
                    if( ( aFDSchemeType == fem::FDScheme_Type::POINT_1_BACKWARD ) ||
                            ( aFDSchemeType == fem::FDScheme_Type::POINT_1_FORWARD ) )
                    {
                        // add unperturbed flux contribution to dfluxdu
                        adnEintdxnDOF_FD.get_column( tDofCounter ) +=
                                tFDScheme( 1 )( 0 ) * tdnEintdxn /
                                ( tFDScheme( 2 )( 0 ) * tDeltaH );

                        // skip first point in FD
                        tStartPoint = 1;
                    }

                    // loop over the points for FD
                    for( uint iPoint = tStartPoint; iPoint < tNumPoints; iPoint++ )
                    {
                        // reset the perturbed coefficients
                        Matrix< DDRMat > tCoeffPert = tCoeff;

                        // perturb the coefficient
                        tCoeffPert( iCoeffRow, iCoeffCol ) += tFDScheme( 0 )( iPoint ) * tDeltaH;

                        // set the perturbed coefficients to FI
                        tFI->set_coeff( tCoeffPert );

                        // reset properties
                        this->reset_eval_flags();

                        // assemble the jacobian
                        adnEintdxnDOF_FD.get_column( tDofCounter ) +=
                                        tFDScheme( 1 )( iPoint ) * this->dnEintdxn( aOrder ) /
                                        ( tFDScheme( 2 )( 0 ) * tDeltaH );
                    }
                    // update dof counter
                    tDofCounter++;
                }
            }
            // reset the coefficients values
            tFI->set_coeff( tCoeff );
        }

        //------------------------------------------------------------------------------
        //------------------------------------------------------------------------------

        void Material_Model::eval_TDvarDOF_FD(
                const Vector< MSI::Dof_Type > & aDofTypes,
                Matrix< DDRMat >                   & aTDvarDOF_FD,
                real                                 aPerturbation,
                MSI::Dof_Type                        aTDvar,
                fem::FDScheme_Type                   aFDSchemeType )
        {
            // get the FD scheme info
            Vector< Vector< real > > tFDScheme;
            fd_scheme( aFDSchemeType, tFDScheme );
            uint tNumPoints = tFDScheme( 0 ).size();

            // get the field interpolator for type
            Field_Interpolator* tFI = mFIManager->get_field_interpolators_for_type( aDofTypes( 0 ) );

            // get number of coefficients, fields and bases for the considered FI
            uint tDerNumDof    = tFI->get_number_of_space_time_coefficients();
            uint tDerNumBases  = tFI->get_number_of_space_time_bases();
            uint tDerNumFields = tFI->get_number_of_fields();

            // get dof index for thermodynamic variable
            uint tTDvarIndex = static_cast< uint >( aTDvar );
            Matrix< DDRMat > tTDvarVal;

            // evaluate unperturbed variable
            switch ( tTDvarIndex )
            {
                // Density DoF type
                case static_cast< int >( MSI::Dof_Type::RHO ):
                {
                    tTDvarVal = this->density();
                    break;
                }

                // Pressure DoF Type
                case static_cast< int >( MSI::Dof_Type::P ):
                {
                    tTDvarVal = this->pressure();
                    break;
                }

                // Temperature DoF Type
                case static_cast< int >( MSI::Dof_Type::TEMP ):
                {
                    tTDvarVal = this->temperature();
                    break;
                }

                default:
                {
                    // throw error
                    MORIS_ERROR( false, "Material_Model::eval_TDvarDOF_FD - Only thermondynamic state variables (rho,p,T) supported." );
                    break;
                }
            }

            // set size for derivative
            aTDvarDOF_FD.set_size( tTDvarVal.n_rows(), tDerNumDof, 0.0 );

            // coefficients for dof type wrt which derivative is computed
            Matrix< DDRMat > tCoeff = tFI->get_coeff();

            // initialize dof counter
            uint tDofCounter = 0;

            // loop over coefficients columns
            for( uint iCoeffCol = 0; iCoeffCol < tDerNumFields; iCoeffCol++ )
            {
                // loop over coefficients rows
                for( uint iCoeffRow = 0; iCoeffRow < tDerNumBases; iCoeffRow++ )
                {
                    // compute the perturbation absolute value
                    real tDeltaH = aPerturbation * tCoeff( iCoeffRow, iCoeffCol );

                    // check that perturbation is not zero
                    if( std::abs( tDeltaH ) < 1e-12 )
                    {
                        tDeltaH = aPerturbation;
                    }

                    // set starting point for FD
                    uint tStartPoint = 0;

                    // if backward or forward add unperturbed contribution
                    if( ( aFDSchemeType == fem::FDScheme_Type::POINT_1_BACKWARD ) ||
                            ( aFDSchemeType == fem::FDScheme_Type::POINT_1_FORWARD ) )
                    {
                        // add unperturbed flux contribution to dfluxdu
                        aTDvarDOF_FD.get_column( tDofCounter ) +=
                                tFDScheme( 1 )( 0 ) * tTDvarVal /
                                ( tFDScheme( 2 )( 0 ) * tDeltaH );

                        // skip first point in FD
                        tStartPoint = 1;
                    }

                    // loop over the points for FD
                    for( uint iPoint = tStartPoint; iPoint < tNumPoints; iPoint++ )
                    {
                        // reset the perturbed coefficients
                        Matrix< DDRMat > tCoeffPert = tCoeff;

                        // perturb the coefficient
                        tCoeffPert( iCoeffRow, iCoeffCol ) += tFDScheme( 0 )( iPoint ) * tDeltaH;

                        // set the perturbed coefficients to FI
                        tFI->set_coeff( tCoeffPert );

                        // reset properties
                        this->reset_eval_flags();

                        // evaluate perturbed variable
                        Matrix< DDRMat > tTDvarValPert;
                        switch ( tTDvarIndex )
                        {
                            // Density DoF type
                            case static_cast< int >( MSI::Dof_Type::RHO ):
                            {
                                tTDvarValPert = this->density();
                                break;
                            }

                            // Pressure DoF Type
                            case static_cast< int >( MSI::Dof_Type::P ):
                            {
                                tTDvarValPert = this->pressure();
                                break;
                            }

                            // Temperature DoF Type
                            case static_cast< int >( MSI::Dof_Type::TEMP ):
                            {
                                tTDvarValPert = this->temperature();
                                break;
                            }

                            default:
                            {
                                // throw error
                                MORIS_ERROR( false, "Material_Model::eval_TDvarDOF_FD - Only thermondynamic state variables (rho,p,T) supported." );
                                break;
                            }
                        }

                        // assemble the jacobian
                        aTDvarDOF_FD.get_column( tDofCounter ) +=
                                        tFDScheme( 1 )( iPoint ) * tTDvarValPert /
                                        ( tFDScheme( 2 )( 0 ) * tDeltaH );
                    }
                    // update dof counter
                    tDofCounter++;
                }
            }
            // reset the coefficients values
            tFI->set_coeff( tCoeff );
        }

        //------------------------------------------------------------------------------

        void Material_Model::eval_TDvarDotDOF_FD(
                const Vector< MSI::Dof_Type > & aDofTypes,
                Matrix< DDRMat >                   & aTDvarDotDOF_FD,
                real                                 aPerturbation,
                MSI::Dof_Type                        aTDvar,
                fem::FDScheme_Type                   aFDSchemeType )
        {
            // get the FD scheme info
            Vector< Vector< real > > tFDScheme;
            fd_scheme( aFDSchemeType, tFDScheme );
            uint tNumPoints = tFDScheme( 0 ).size();

            // get the field interpolator for type
            Field_Interpolator* tFI = mFIManager->get_field_interpolators_for_type( aDofTypes( 0 ) );

            // get number of coefficients, fields and bases for the considered FI
            uint tDerNumDof    = tFI->get_number_of_space_time_coefficients();
            uint tDerNumBases  = tFI->get_number_of_space_time_bases();
            uint tDerNumFields = tFI->get_number_of_fields();

            // get dof index for thermodynamic variable
            uint tTDvarIndex = static_cast< uint >( aTDvar );
            Matrix< DDRMat > tTDvarDotVal;

            // evaluate unperturbed variable
            switch ( tTDvarIndex )
            {
                // Density DoF type
                case static_cast< int >( MSI::Dof_Type::RHO ):
                {
                    tTDvarDotVal = this->DensityDot();
                    break;
                }

                // Pressure DoF Type
                case static_cast< int >( MSI::Dof_Type::P ):
                {
                    tTDvarDotVal = this->PressureDot();
                    break;
                }

                // Temperature DoF Type
                case static_cast< int >( MSI::Dof_Type::TEMP ):
                {
                    tTDvarDotVal = this->TemperatureDot();
                    break;
                }

                default:
                {
                    // throw error
                    MORIS_ERROR( false, "Material_Model::eval_TDvarDOF_FD - Only thermondynamic state variables (rho,p,T) supported." );
                    break;
                }
            }

            // set size for derivative
            aTDvarDotDOF_FD.set_size( tTDvarDotVal.n_rows(), tDerNumDof, 0.0 );

            // coefficients for dof type wrt which derivative is computed
            Matrix< DDRMat > tCoeff = tFI->get_coeff();

            // initialize dof counter
            uint tDofCounter = 0;

            // loop over coefficients columns
            for( uint iCoeffCol = 0; iCoeffCol < tDerNumFields; iCoeffCol++ )
            {
                // loop over coefficients rows
                for( uint iCoeffRow = 0; iCoeffRow < tDerNumBases; iCoeffRow++ )
                {
                    // compute the perturbation absolute value
                    real tDeltaH = aPerturbation * tCoeff( iCoeffRow, iCoeffCol );

                    // check that perturbation is not zero
                    if( std::abs( tDeltaH ) < 1e-12 )
                    {
                        tDeltaH = aPerturbation;
                    }

                    // set starting point for FD
                    uint tStartPoint = 0;

                    // if backward or forward add unperturbed contribution
                    if( ( aFDSchemeType == fem::FDScheme_Type::POINT_1_BACKWARD ) ||
                            ( aFDSchemeType == fem::FDScheme_Type::POINT_1_FORWARD ) )
                    {
                        // add unperturbed flux contribution to dfluxdu
                        aTDvarDotDOF_FD.get_column( tDofCounter ) +=
                                tFDScheme( 1 )( 0 ) * tTDvarDotVal /
                                ( tFDScheme( 2 )( 0 ) * tDeltaH );

                        // skip first point in FD
                        tStartPoint = 1;
                    }

                    // loop over the points for FD
                    for( uint iPoint = tStartPoint; iPoint < tNumPoints; iPoint++ )
                    {
                        // reset the perturbed coefficients
                        Matrix< DDRMat > tCoeffPert = tCoeff;

                        // perturb the coefficient
                        tCoeffPert( iCoeffRow, iCoeffCol ) += tFDScheme( 0 )( iPoint ) * tDeltaH;

                        // set the perturbed coefficients to FI
                        tFI->set_coeff( tCoeffPert );

                        // reset properties
                        this->reset_eval_flags();

                        // evaluate perturbed variable
                        Matrix< DDRMat > tTDvarDotValPert;
                        switch ( tTDvarIndex )
                        {
                            // Density DoF type
                            case static_cast< int >( MSI::Dof_Type::RHO ):
                            {
                                tTDvarDotValPert = this->DensityDot();
                                break;
                            }

                            // Pressure DoF Type
                            case static_cast< int >( MSI::Dof_Type::P ):
                            {
                                tTDvarDotValPert = this->PressureDot();
                                break;
                            }

                            // Temperature DoF Type
                            case static_cast< int >( MSI::Dof_Type::TEMP ):
                            {
                                tTDvarDotValPert = this->TemperatureDot();
                                break;
                            }

                            default:
                            {
                                // throw error
                                MORIS_ERROR( false, "Material_Model::eval_TDvarDOF_FD - Only thermondynamic state variables (rho,p,T) supported." );
                                break;
                            }
                        }

                        // assemble the jacobian
                        aTDvarDotDOF_FD.get_column( tDofCounter ) +=
                                        tFDScheme( 1 )( iPoint ) * tTDvarDotValPert /
                                        ( tFDScheme( 2 )( 0 ) * tDeltaH );
                    }
                    // update dof counter
                    tDofCounter++;
                }
            }
            // reset the coefficients values
            tFI->set_coeff( tCoeff );
        }

        //------------------------------------------------------------------------------

        void Material_Model::eval_dnTDvardxnDOF_FD(
                const Vector< MSI::Dof_Type > & aDofTypes,
                Matrix< DDRMat >                   & adnTDvardxnDOF_FD,
                real                                 aPerturbation,
                MSI::Dof_Type                        aTDvar,
                uint                                 aOrder,
                fem::FDScheme_Type                   aFDSchemeType )
        {
            // get the FD scheme info
            Vector< Vector< real > > tFDScheme;
            fd_scheme( aFDSchemeType, tFDScheme );
            uint tNumPoints = tFDScheme( 0 ).size();

            // get the field interpolator for type
            Field_Interpolator* tFI = mFIManager->get_field_interpolators_for_type( aDofTypes( 0 ) );

            // get number of coefficients, fields and bases for the considered FI
            uint tDerNumDof    = tFI->get_number_of_space_time_coefficients();
            uint tDerNumBases  = tFI->get_number_of_space_time_bases();
            uint tDerNumFields = tFI->get_number_of_fields();

            // get dof index for thermodynamic variable
            uint tTDvarIndex = static_cast< uint >( aTDvar );
            Matrix< DDRMat > tdnTDvardxnVal;

            // evaluate unperturbed variable
            switch ( tTDvarIndex )
            {
                // Density DoF type
                case static_cast< int >( MSI::Dof_Type::RHO ):
                {
                    tdnTDvardxnVal = this->dnDensitydxn( aOrder );
                    break;
                }

                // Pressure DoF Type
                case static_cast< int >( MSI::Dof_Type::P ):
                {
                    tdnTDvardxnVal = this->dnPressuredxn( aOrder );
                    break;
                }

                // Temperature DoF Type
                case static_cast< int >( MSI::Dof_Type::TEMP ):
                {
                    tdnTDvardxnVal = this->dnTemperaturedxn( aOrder );
                    break;
                }

                default:
                {
                    // throw error
                    MORIS_ERROR( false, "Material_Model::eval_TDvarDOF_FD - Only thermondynamic state variables (rho,p,T) supported." );
                    break;
                }
            }

            // set size for derivative
            adnTDvardxnDOF_FD.set_size( tdnTDvardxnVal.n_rows(), tDerNumDof, 0.0 );

            // coefficients for dof type wrt which derivative is computed
            Matrix< DDRMat > tCoeff = tFI->get_coeff();

            // initialize dof counter
            uint tDofCounter = 0;

            // loop over coefficients columns
            for( uint iCoeffCol = 0; iCoeffCol < tDerNumFields; iCoeffCol++ )
            {
                // loop over coefficients rows
                for( uint iCoeffRow = 0; iCoeffRow < tDerNumBases; iCoeffRow++ )
                {
                    // compute the perturbation absolute value
                    real tDeltaH = aPerturbation * tCoeff( iCoeffRow, iCoeffCol );

                    // check that perturbation is not zero
                    if( std::abs( tDeltaH ) < 1e-12 )
                    {
                        tDeltaH = aPerturbation;
                    }

                    // set starting point for FD
                    uint tStartPoint = 0;

                    // if backward or forward add unperturbed contribution
                    if( ( aFDSchemeType == fem::FDScheme_Type::POINT_1_BACKWARD ) ||
                            ( aFDSchemeType == fem::FDScheme_Type::POINT_1_FORWARD ) )
                    {
                        // add unperturbed flux contribution to dfluxdu
                        adnTDvardxnDOF_FD.get_column( tDofCounter ) +=
                                tFDScheme( 1 )( 0 ) * tdnTDvardxnVal /
                                ( tFDScheme( 2 )( 0 ) * tDeltaH );

                        // skip first point in FD
                        tStartPoint = 1;
                    }

                    // loop over the points for FD
                    for( uint iPoint = tStartPoint; iPoint < tNumPoints; iPoint++ )
                    {
                        // reset the perturbed coefficients
                        Matrix< DDRMat > tCoeffPert = tCoeff;

                        // perturb the coefficient
                        tCoeffPert( iCoeffRow, iCoeffCol ) += tFDScheme( 0 )( iPoint ) * tDeltaH;

                        // set the perturbed coefficients to FI
                        tFI->set_coeff( tCoeffPert );

                        // reset properties
                        this->reset_eval_flags();

                        // evaluate perturbed variable
                        Matrix< DDRMat > tdnTDvardxnValPert;
                        switch ( tTDvarIndex )
                        {
                            // Density DoF type
                            case static_cast< int >( MSI::Dof_Type::RHO ):
                            {
                                tdnTDvardxnValPert = this->dnDensitydxn( aOrder );
                                break;
                            }

                            // Pressure DoF Type
                            case static_cast< int >( MSI::Dof_Type::P ):
                            {
                                tdnTDvardxnValPert = this->dnPressuredxn( aOrder );
                                break;
                            }

                            // Temperature DoF Type
                            case static_cast< int >( MSI::Dof_Type::TEMP ):
                            {
                                tdnTDvardxnValPert = this->dnTemperaturedxn( aOrder );
                                break;
                            }

                            default:
                            {
                                // throw error
                                MORIS_ERROR( false, "Material_Model::eval_TDvarDOF_FD - Only thermondynamic state variables (rho,p,T) supported." );
                                break;
                            }
                        }

                        // assemble the jacobian
                        adnTDvardxnDOF_FD.get_column( tDofCounter ) +=
                                        tFDScheme( 1 )( iPoint ) * tdnTDvardxnValPert /
                                        ( tFDScheme( 2 )( 0 ) * tDeltaH );
                    }
                    // update dof counter
                    tDofCounter++;
                }
            }
            // reset the coefficients values
            tFI->set_coeff( tCoeff );
        }

        //------------------------------------------------------------------------------

        void Material_Model::eval_QuantityDOF_FD(
                const Vector< MSI::Dof_Type > & aDofTypes,
                Matrix< DDRMat >                   & aQuantityDOF_FD,
                std::string                          aQuantityString,
                real                                 aPerturbation,
                fem::FDScheme_Type                   aFDSchemeType )
        {
            // set function pointer based on string
            const Matrix< DDRMat > &  ( Material_Model:: * t_get_Quantity )();
            if ( aQuantityString == "Cp" )
            {
                t_get_Quantity = &Material_Model::Cp;
            }
            else if ( aQuantityString == "Cv" )
            {
                t_get_Quantity = &Material_Model::Cv;
            }
            else if ( aQuantityString == "Gamma" )
            {
                t_get_Quantity = &Material_Model::Gamma;
            }
            else if ( aQuantityString == "AlphaP" )
            {
                t_get_Quantity = &Material_Model::AlphaP;
            }
            else if ( aQuantityString == "BetaT" )
            {
                t_get_Quantity = &Material_Model::BetaT;
            }
            else
            {
                MORIS_ERROR( false, "Material_Model::eval_QuantityDOF_FD - Unknown aQuantityString." );
            }

            // get the FD scheme info
            Vector< Vector< real > > tFDScheme;
            fd_scheme( aFDSchemeType, tFDScheme );
            uint tNumPoints = tFDScheme( 0 ).size();

            // get the field interpolator for type
            Field_Interpolator* tFI = mFIManager->get_field_interpolators_for_type( aDofTypes( 0 ) );

            // get number of coefficients, fields and bases for the considered FI
            uint tDerNumDof    = tFI->get_number_of_space_time_coefficients();
            uint tDerNumBases  = tFI->get_number_of_space_time_bases();
            uint tDerNumFields = tFI->get_number_of_fields();

            // evaluate unperturbed internal energy
            Matrix< DDRMat > tQuantity = ( this->*t_get_Quantity )();

            // set size for derivative
            aQuantityDOF_FD.set_size( tQuantity.n_rows(), tDerNumDof, 0.0 );

            // coefficients for dof type wrt which derivative is computed
            Matrix< DDRMat > tCoeff = tFI->get_coeff();

            // initialize dof counter
            uint tDofCounter = 0;

            // loop over coefficients columns
            for( uint iCoeffCol = 0; iCoeffCol < tDerNumFields; iCoeffCol++ )
            {
                // loop over coefficients rows
                for( uint iCoeffRow = 0; iCoeffRow < tDerNumBases; iCoeffRow++ )
                {
                    // compute the perturbation absolute value
                    real tDeltaH = aPerturbation * tCoeff( iCoeffRow, iCoeffCol );

                    // check that perturbation is not zero
                    if( std::abs( tDeltaH ) < 1e-12 )
                    {
                        tDeltaH = aPerturbation;
                    }

                    // set starting point for FD
                    uint tStartPoint = 0;

                    // if backward or forward add unperturbed contribution
                    if( ( aFDSchemeType == fem::FDScheme_Type::POINT_1_BACKWARD ) ||
                            ( aFDSchemeType == fem::FDScheme_Type::POINT_1_FORWARD ) )
                    {
                        // add unperturbed flux contribution to dfluxdu
                        aQuantityDOF_FD.get_column( tDofCounter ) +=
                                tFDScheme( 1 )( 0 ) * tQuantity /
                                ( tFDScheme( 2 )( 0 ) * tDeltaH );

                        // skip first point in FD
                        tStartPoint = 1;
                    }

                    // loop over the points for FD
                    for( uint iPoint = tStartPoint; iPoint < tNumPoints; iPoint++ )
                    {
                        // reset the perturbed coefficients
                        Matrix< DDRMat > tCoeffPert = tCoeff;

                        // perturb the coefficient
                        tCoeffPert( iCoeffRow, iCoeffCol ) += tFDScheme( 0 )( iPoint ) * tDeltaH;

                        // set the perturbed coefficients to FI
                        tFI->set_coeff( tCoeffPert );

                        // reset properties
                        this->reset_eval_flags();

                        // assemble the jacobian
                        aQuantityDOF_FD.get_column( tDofCounter ) +=
                                        tFDScheme( 1 )( iPoint ) * ( this->*t_get_Quantity )() /
                                        ( tFDScheme( 2 )( 0 ) * tDeltaH );
                    }
                    // update dof counter
                    tDofCounter++;
                }
            }
            // reset the coefficients values
            tFI->set_coeff( tCoeff );
        }

        //------------------------------------------------------------------------------

    }/* end_fem_namespace */
}/* end_moris_namespace */

