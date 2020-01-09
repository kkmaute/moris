/*
 * cl_FEM_Field_Interpolator_Manager.hpp
 *
 *  Created on: Apr 20, 2019
 *      Author: schmidt
 */

#ifndef SRC_FEM_CL_FEM_FIELD_INTERPLATOR_MANAGER_HPP_
#define SRC_FEM_CL_FEM_FIELD_INTERPLATOR_MANAGER_HPP_

#include "assert.h"
#include <cmath>

#include "typedefs.hpp"                     //MRS/COR/src
#include "cl_Matrix.hpp"
#include "linalg_typedefs.hpp"
#include "cl_Cell.hpp"
#include "cl_MTK_Cell.hpp"                  //MTK/src
#include "cl_MSI_Equation_Object.hpp"       //FEM/MSI/src
#include "cl_MSI_Equation_Set.hpp"          //FEM/MSI/src
#include "cl_MSI_Model_Solver_Interface.hpp"          //FEM/MSI/src
#include "cl_FEM_Enums.hpp"                 //FEM/INT/src
#include "cl_FEM_Node.hpp"                  //FEM/INT/src
#include "cl_FEM_IWG.hpp"                   //FEM/INT/src
#include "cl_FEM_Geometry_Interpolator.hpp" //FEM/INT/src
#include "cl_FEM_Field_Interpolator.hpp"    //FEM/INT/src
#include "cl_FEM_Integrator.hpp"            //FEM/INT/src
#include "cl_FEM_Element_Factory.hpp"       //FEM/INT/src
#include "cl_FEM_Set.hpp"                   //FEM/INT/src

namespace moris
{
    namespace fem
    {
    class Set;
//------------------------------------------------------------------------------
    /**
     * Field interpolator manager
     */
    class Field_Interpolator_Manager
    {

    private:

        // dof type list for the FI manager
        const moris::Cell< moris::Cell< enum MSI::Dof_Type > > & mDofTypes;

        // pointer to the equation set
        MSI::Equation_Set * mEquationSet = nullptr;

        // enum for master or slave
        mtk::Master_Slave mIsMaster;

        // dof type map
        moris::Matrix< DDSMat > mDofTypeMap;

        // list of field intepolators
        moris::Cell< Field_Interpolator* > mFI;

        // maximum number of field interpolators
        moris::uint mMaxNumDofFI;

        // dof type list for the FI manager
        const moris::Cell< moris::Cell< enum MSI::Dv_Type > > mDvTypes;

        // dof type map
        moris::Matrix< DDSMat > mDvTypeMap;

        // list of field intepolators
        moris::Cell< Field_Interpolator* > mDvFI;

        // maximum number of field interpolators
        moris::uint mMaxNumDvFI;

//------------------------------------------------------------------------------
    public:
//------------------------------------------------------------------------------
        /**
         * constructor
         * @param[ in ] aDofTypes    a list of group of dof types
         * @param[ in ] aEquationSet a pointer to the corresponding equation set
         */
        Field_Interpolator_Manager( const moris::Cell< moris::Cell< enum MSI::Dof_Type > > & aDofTypes,
                                          MSI::Equation_Set                                * aEquationSet,
                                          mtk::Master_Slave                                  aIsMaster = mtk::Master_Slave::MASTER )
        : mDofTypes( aDofTypes ),
          mEquationSet( aEquationSet ),
          mIsMaster( aIsMaster )
        {
            // set the dof type map
            mDofTypeMap = mEquationSet->get_dof_type_map( aIsMaster );

            // maximum number of dof field interpolators
            mMaxNumDofFI =  mEquationSet->get_num_unique_dof_types();

            // FIXME default
            mMaxNumDvFI = 0;
        };

        Field_Interpolator_Manager( const moris::Cell< moris::Cell< enum MSI::Dof_Type > > & aDofTypes,
                                    const moris::Cell< moris::Cell< enum MSI::Dv_Type > >  & aDvTypes,
                                          MSI::Equation_Set                                * aEquationSet,
                                          mtk::Master_Slave                                  aIsMaster = mtk::Master_Slave::MASTER )
        : mDofTypes( aDofTypes ),
          mEquationSet( aEquationSet ),
          mIsMaster( aIsMaster ),
          mDvTypes( aDvTypes )
        {
            // set the dof type map
            mDofTypeMap = mEquationSet->get_dof_type_map( aIsMaster );

            // maximum number of dof field interpolators
            mMaxNumDofFI =  mEquationSet->get_num_unique_dof_types();

            // FIXME maximum number of dv field interpolators
            mMaxNumDvFI =  mDvTypes.size();

        };

//------------------------------------------------------------------------------
        /**
         * constructor
         * @param[ in ] aDofTypes a list of group of dof types
         * @param[ in ] aEquationSet a pointer to the corresponding equation set
         * @param[ in ] aModelSolverInterface a pointer to the corresponding model solver interface
         */
        Field_Interpolator_Manager( const moris::Cell< moris::Cell< enum MSI::Dof_Type > > & aDofTypes,
                                          MSI::Equation_Set                                * aEquationSet,
                                          MSI::Model_Solver_Interface                      * aModelSolverInterface,
                                          mtk::Master_Slave                                  aIsMaster = mtk::Master_Slave::MASTER )
        : mDofTypes( aDofTypes ),
          mEquationSet( aEquationSet ),
          mIsMaster( aIsMaster )
        {
            // set the dof type map
            mDofTypeMap = mEquationSet->get_dof_type_map( aIsMaster );

            // maximum number of dof field interpolators
            mMaxNumDofFI = mEquationSet->get_num_unique_dof_types();

            // FIXME default
            mMaxNumDvFI = 0;
        };

//------------------------------------------------------------------------------
        /**
         * destructor
         */
        ~Field_Interpolator_Manager()
        {
            // delete pointers on the FI manager
            this->delete_pointers();
        };

//------------------------------------------------------------------------------
        /**
         * delete pointers
         */
        void delete_pointers()
        {
            // delete the field interpolator pointers
            for( Field_Interpolator* tFI : mFI )
            {
                delete tFI;
            }
            mFI.clear();
        }

//------------------------------------------------------------------------------
        /**
         * create the field interpolator for the FI manager
         * @param[ in ] aModelSolverInterface a pointer to the corresponding model solver interface
         */
        void create_field_interpolators( MSI::Model_Solver_Interface * aModelSolverInterface )
        {
            // dof field interpolators------------------------------------------

            // set the size of the cell of field interpolators
            mFI.resize( mMaxNumDofFI, nullptr );

            // loop over the dof type groups
            for( uint iDof = 0; iDof < mDofTypes.size(); iDof++ )
            {
                // get the number of time level for the dof type group
                uint tNumTimeNodes = aModelSolverInterface->get_time_levels_for_type( mDofTypes( iDof )( 0 ) );

                // get the set index for the dof type group
                uint tDofIndex = mEquationSet->get_dof_index_for_type_1( mDofTypes( iDof )( 0 ), mIsMaster );

                // create the field interpolation rule for the dof type group
                Interpolation_Rule tFieldInterpolationRule( static_cast< Set* >( mEquationSet )->mIPGeometryType,
                                                            Interpolation_Type::LAGRANGE,
                                                            static_cast< Set* >( mEquationSet )->mIPSpaceInterpolationOrder,
                                                            static_cast< Set* >( mEquationSet )->get_auto_time_interpolation_type( tNumTimeNodes ), // fixme
                                                            // If interpolation type CONSTANT, iInterpolation order is not used
                                                            static_cast< Set* >( mEquationSet )->get_auto_interpolation_order( tNumTimeNodes, mtk::Geometry_Type::LINE ) ); //fixme

                // check if the fiedl interpolator was created previously
                MORIS_ASSERT( mFI( tDofIndex ) == nullptr, "Field_Interpolator_Manager::create_field_interpolators - Field interpolator was created previously" );

                // create a field interpolator for the dof type group
                mFI( tDofIndex ) = new Field_Interpolator( mDofTypes( iDof ).size(),
                                                           tFieldInterpolationRule,
                                                           static_cast< Set* >( mEquationSet )->get_IP_geometry_interpolator( mIsMaster ),
                                                           mDofTypes( iDof ) );
            }

            // dv field interpolators------------------------------------------

            // set the size of the cell of field interpolators
            mDvFI.resize( mMaxNumDvFI, nullptr );

            // loop over the dv type groups
            for( uint iDv = 0; iDv < mDvTypes.size(); iDv++ )
            {
                // get the number of time level for the dv type group
                // FIXME where do we get this info
                uint tNumTimeNodes = 1;

                // get the set index for the dv type group
                uint tDvIndex = mEquationSet->get_dv_index_for_type_1( mDvTypes( iDv )( 0 ), mIsMaster );

                // create the field interpolation rule for the dv type group
                Interpolation_Rule tFieldInterpolationRule( static_cast< Set* >( mEquationSet )->mIPGeometryType,
                                                            Interpolation_Type::LAGRANGE,
                                                            static_cast< Set* >( mEquationSet )->mIPSpaceInterpolationOrder,
                                                            static_cast< Set* >( mEquationSet )->get_auto_time_interpolation_type( tNumTimeNodes ), // fixme
                                                            // If interpolation type CONSTANT, iInterpolation order is not used
                                                            static_cast< Set* >( mEquationSet )->get_auto_interpolation_order( tNumTimeNodes, mtk::Geometry_Type::LINE ) ); //fixme

                // check if the field interpolator was created previously
                MORIS_ASSERT( mDvFI( tDvIndex ) == nullptr,
                              "Field_Interpolator_Manager::create_field_interpolators - Field interpolator was created previously." );

                // create a field interpolator for the dof type group
                mFI( tDvIndex ) = new Field_Interpolator( mDvTypes( iDv ).size(),
                                                          tFieldInterpolationRule,
                                                          static_cast< Set* >( mEquationSet )->get_IP_geometry_interpolator( mIsMaster ),
                                                          mDvTypes( iDv ) );
            }
        };

//------------------------------------------------------------------------------
        /**
         * get the maximum number of dof field interpolators on the manager
         * @param[ out ] mMaxNumDofFI the maximum number of dof FI on the manager
         */
        moris::uint get_max_num_field_interpolators()
        {
            return mMaxNumDofFI;
        }

//------------------------------------------------------------------------------
        /**
         * get the maximum number of dv field interpolators on the manager
         * @param[ out ] mMaxNumDvFI the maximum number of dv FI on the manager
         */
        moris::uint get_max_num_dv_field_interpolators()
        {
            return mMaxNumDvFI;
        }

//------------------------------------------------------------------------------
        /**
         * get the field interpolator for a given dof type
         * @param[ in ] aDofType a dof type enum
         */
        Field_Interpolator * get_field_interpolators_for_type( enum MSI::Dof_Type aDofType )
        {
            // check of the equation set pointer was set for the FI manager
            MORIS_ASSERT( mEquationSet != nullptr, "Field_Interpolator_Manager::get_field_interpolators_for_type - Equation Set pointer not set");

            // get the set index for the requested dof type
            sint tDofIndex = mEquationSet->get_dof_index_for_type_1( aDofType, mIsMaster );

            // if the index was set for the equation set
            if( tDofIndex != -1 )
            {
                // check if the FI exists for the FI manager
                MORIS_ASSERT( (sint)mFI.size() > tDofIndex,
                              "Field_Interpolator_Manager::get_field_interpolators_for_type - field interpolator does not exist" );

                // return the FI
                return mFI( tDofIndex );
            }
            else
            {
                return nullptr;
            }
        }

//------------------------------------------------------------------------------
        /**
         * get the field interpolator for a given dv type
         * @param[ in ] aDvType a dv type enum
         */
        Field_Interpolator * get_field_interpolators_for_type( enum MSI::Dv_Type aDvType )
        {
            // get the set index for the requested dv type
            sint tDvIndex =  mEquationSet->get_dv_index_for_type_1( aDvType, mIsMaster );

            // if the index was set for the equation set
            if( tDvIndex != -1 )
            {
                // check if the FI exists for the FI manager
                MORIS_ASSERT( (sint)mDvFI.size() > tDvIndex,
                              "Field_Interpolator_Manager::get_field_interpolators_for_type - field interpolator does not exist" );

                // return the FI
                return mDvFI( tDvIndex );
            }
            else
            {
                return nullptr;
            }
        }

//------------------------------------------------------------------------------
        /**
         * set an evaluation point in space and time
         * @param[ in ] aParamPoint coordinates of an evaluation point
         */
        void set_space_time( Matrix< DDRMat > & aParamPoint )
        {
            // loop over the dof field interpolators
            for ( uint iDofFI = 0; iDofFI < mDofTypes.size(); iDofFI++ )
            {
                // get the set index for the dof type
                sint tDofIndex = mDofTypeMap( static_cast< uint >( mDofTypes( iDofFI )( 0 ) ) );

                // set the evaluation point
                mFI( tDofIndex )->set_space_time( aParamPoint );
            }

            // loop over the dv field interpolators
            for ( uint iDvFI = 0; iDvFI < mDvTypes.size(); iDvFI++ )
            {
                // get the set index for the dv type
                sint tDvIndex = mDvTypeMap( static_cast< uint >( mDvTypes( iDvFI )( 0 ) ) );

                // set the evaluation point
                mDvFI( tDvIndex )->set_space_time( aParamPoint );
            }
        }

//------------------------------------------------------------------------------
        /**
         * set coefficients for field interpolator with specific dof type
         * @param[ in ] aDofType a dof type for which the coeff will be set
         * @param[ in ] aCoeff   coefficients to be set
         */
        void set_coeff_for_type( enum MSI::Dof_Type   aDofType,
                                 Matrix< DDRMat >   & aCoeff )
        {
            // get field interpolator for dof type and set coefficients
            this->get_field_interpolators_for_type( aDofType )->set_coeff( aCoeff );
        }

//------------------------------------------------------------------------------
        /**
         * set coefficients for field interpolator with specific dof type
         * @param[ in ] aDofType a dof type for which the coeff will be set
         * @param[ in ] aCoeff   coefficients to be set
         */
        void set_coeff_for_type( enum MSI::Dv_Type    aDvType,
                                 Matrix< DDRMat >   & aCoeff )
        {
            // get field interpolator for dof type and set coefficients
            this->get_field_interpolators_for_type( aDvType )->set_coeff( aCoeff );
        }

//------------------------------------------------------------------------------
    };

//------------------------------------------------------------------------------
    } /* namespace fem */
} /* namespace moris */

#endif /* SRC_FEM_CL_FEM_FIELD_INTERPLATOR_MANAGER_HPP_ */
