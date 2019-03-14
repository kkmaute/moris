#include <iostream>
#include "cl_FEM_Element_Sideset.hpp" //FEM/INT/src
#include "cl_FEM_Integrator.hpp"      //FEM/INT/src

namespace moris
{
    namespace fem
    {
//------------------------------------------------------------------------------

        Element_Sideset::Element_Sideset( mtk::Cell                 * aCell,
                                          moris::Cell< IWG* >       & aIWGs,
                                          moris::Cell< Node_Base* > & aNodes )
                                        : Element( aCell, aIWGs, aNodes )
        {
        }

//------------------------------------------------------------------------------

        Element_Sideset::~Element_Sideset()
        {
        }

//------------------------------------------------------------------------------

        void Element_Sideset::compute_jacobian()
        {
            MORIS_ERROR( false, " Element_Sideset::compute_jacobian - not implemented. ");
        }

//------------------------------------------------------------------------------

        void Element_Sideset::compute_residual()
        {
            MORIS_ERROR( false, " Element_Sideset::compute_residual - not implemented. ");
        }

//------------------------------------------------------------------------------

        void Element_Sideset::compute_jacobian_and_residual()
        {
            MORIS_ERROR( false, " Element_Sideset::compute_jacobian_and_residual - not implemented. ");
        }

//------------------------------------------------------------------------------

        moris::Cell< fem::Field_Interpolator* >
        Element_Sideset::create_element_field_interpolators
        ( fem::Geometry_Interpolator* aGeometryInterpolator )
        {
            // cell of field interpolators
            Cell< Field_Interpolator* > tFieldInterpolators( mNumOfInterp, nullptr );

            // loop on the dof type groups and create a field interpolator for each
            for( uint i = 0; i < mNumOfInterp; i++ )
            {
                // get the ith dof type group
                Cell< MSI::Dof_Type > tDofTypeGroup = mInterpDofTypeList( i );

                // create the field interpolation rule for the ith dof type group
                //FIXME: space interpolation based on the mtk::Cell
                //FIXME: time  interpolation set to constant
                Interpolation_Rule tFieldInterpolationRule( mCell->get_geometry_type(),
                                                            Interpolation_Type::LAGRANGE,
                                                            this->get_auto_interpolation_order(),
                                                            Interpolation_Type::CONSTANT,
                                                            mtk::Interpolation_Order::CONSTANT );

                // get number of field interpolated by the ith field interpolator
                uint tNumOfFields = tDofTypeGroup.size();

                // create an interpolator for the ith dof type group
                tFieldInterpolators( i ) = new Field_Interpolator( tNumOfFields,
                                                                   tFieldInterpolationRule,
                                                                   aGeometryInterpolator );
            }
            return tFieldInterpolators;
        }

//------------------------------------------------------------------------------

        void
        Element_Sideset::set_element_field_interpolators_coefficients
        ( moris::Cell< Field_Interpolator* > & aFieldInterpolators )
        {
            // loop on the dof types
            for( uint i = 0; i < mNumOfInterp; i++ )
            {
                // get the ith dof type group
                Cell< MSI::Dof_Type > tDofTypeGroup = mInterpDofTypeList( i );

                //FIXME:forced coefficients
                // get the pdof values for the ith dof type group
                Matrix< DDRMat > tCoeff;
                this->get_my_pdof_values( tDofTypeGroup, tCoeff );

                // set the field coefficients
                aFieldInterpolators( i )->set_coeff( tCoeff );
            }
        }

//------------------------------------------------------------------------------
        void
        Element_Sideset::initialize_mJacobianElement_and_mResidualElement
        ( moris::Cell< Field_Interpolator* > & aFieldInterpolators )
        {
            mJacobianElement.resize( mNumOfInterp * mNumOfInterp );
            mResidualElement.resize( mNumOfInterp );

            uint tTotalDof = 0;
            for( uint i = 0; i < mNumOfInterp; i++ )
            {
                // get number of pdofs for the ith dof type
                uint tNumOfDofi = aFieldInterpolators( i )->get_number_of_space_time_coefficients();

                // get total number of dof
                tTotalDof = tTotalDof + tNumOfDofi;

                // set mResidualElement size
                mResidualElement( i ).set_size( tNumOfDofi, 1, 0.0 );

                for( uint j = 0; j < mNumOfInterp; j++ )
                {
                    // get number of pdofs for the ith dof type
                    uint tNumOfDofj = aFieldInterpolators( j )->get_number_of_space_time_coefficients();

                    // set mResidualElement size
                    mJacobianElement( i * mNumOfInterp + j ).set_size( tNumOfDofi, tNumOfDofj, 0.0 );
                }
            }

            mJacobian.set_size( tTotalDof, tTotalDof, 0.0 );
            mResidual.set_size( tTotalDof, 1, 0.0 );
        }


//------------------------------------------------------------------------------
        moris::Cell< Field_Interpolator* >
        Element_Sideset::get_IWG_field_interpolators( IWG*                               & aIWG,
                                                      moris::Cell< Field_Interpolator* > & aFieldInterpolators )
        {
            // ask the IWG for its active dof types
            Cell< Cell< MSI::Dof_Type > > tIWGActiveDof = aIWG->get_active_dof_types();

            // number of active dof type for the IWG
            uint tNumOfIWGActiveDof = tIWGActiveDof.size();

            // select associated active interpolators
            Cell< Field_Interpolator* > tIWGFieldInterpolators( tNumOfIWGActiveDof, nullptr );
            for( uint i = 0; i < tNumOfIWGActiveDof; i++ )
            {
                // find the index of active dof type in the list of element dof type
                uint tIWGDofIndex = mInterpDofTypeMap( static_cast< int >( tIWGActiveDof( i )( 0 ) ) );

                // select the corresponding interpolator
                tIWGFieldInterpolators( i ) = aFieldInterpolators( tIWGDofIndex );
            }
            return tIWGFieldInterpolators;
        }

//------------------------------------------------------------------------------

    } /* namespace fem */
} /* namespace moris */
