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
     * \brief element class that communicates with the mesh interface
     */
    class Field_Interpolator_Manager
    {

    private:

        const moris::Cell< moris::Cell< enum MSI::Dof_Type > > & mMasterDofTypes;
        const moris::Cell< moris::Cell< enum MSI::Dof_Type > > & mSlaveDofTypes;

        //! list of master field interpolator pointers
        moris::Cell< Field_Interpolator* >  mMasterFI;

        //! list of slave field interpolator pointers
        moris::Cell< Field_Interpolator* >  mSlaveFI;

        MSI::Equation_Set * mEquationSet = nullptr;

        moris::uint mMaxNumFieldInterpolators;

//------------------------------------------------------------------------------
    public:
//------------------------------------------------------------------------------

        Field_Interpolator_Manager( const moris::Cell< moris::Cell< enum MSI::Dof_Type > > & aMasterDofTypes,
                                    const moris::Cell< moris::Cell< enum MSI::Dof_Type > > & aSlaveDofTypes,
                                          MSI::Equation_Set                                * aEquationSet ) : mMasterDofTypes( aMasterDofTypes ),
                                                                                                              mSlaveDofTypes( aSlaveDofTypes ),
                                                                                                              mEquationSet( aEquationSet ),
                                                                                                              mMaxNumFieldInterpolators( mEquationSet->get_num_unique_dof_types() )
        {

        };
//------------------------------------------------------------------------------

        Field_Interpolator_Manager( const moris::Cell< moris::Cell< enum MSI::Dof_Type > > & aMasterDofTypes,
                                    const moris::Cell< moris::Cell< enum MSI::Dof_Type > > & aSlaveDofTypes,
                                          MSI::Equation_Set                                * aEquationSet,
                                          MSI::Model_Solver_Interface                      * aModelSolverInterface ) : mMasterDofTypes( aMasterDofTypes ),
                                                                                                                       mSlaveDofTypes( aSlaveDofTypes ),
                                                                                                                       mEquationSet( aEquationSet ),
                                                                                                                       mMaxNumFieldInterpolators( mEquationSet->get_num_unique_dof_types() )
        {

        };

//------------------------------------------------------------------------------

        ~Field_Interpolator_Manager()
        {
            this->delete_pointers();
        };

//------------------------------------------------------------------------------

        void delete_pointers()
        {
            // delete the master field interpolator pointers
            for( Field_Interpolator* tMasterFieldInterpolator : mMasterFI )
            {
                delete tMasterFieldInterpolator;
            }
            mMasterFI.clear();

            // delete the slave field interpolator pointers
            for( Field_Interpolator* tSlaveFieldInterpolator : mSlaveFI )
            {
                delete tSlaveFieldInterpolator;
            }
            mSlaveFI.clear();
        }

//------------------------------------------------------------------------------

        void create_field_interpolators( MSI::Model_Solver_Interface * aModelSolverInterface )
        {
            // set the size of the cell of field interpolators
            mMasterFI.resize( mMaxNumFieldInterpolators, nullptr );
            mSlaveFI .resize( mMaxNumFieldInterpolators, nullptr );

            // loop over the master dof type groups
            for( uint iDOF = 0; iDOF < mMasterDofTypes.size(); iDOF++ )
            {
                // get the number of time level for the dof type group
                uint tNumTimeNodes = aModelSolverInterface->get_time_levels_for_type( mMasterDofTypes( iDOF )( 0 ) );

                uint tDofIndex = mEquationSet->get_dof_index_for_type_1( mMasterDofTypes( iDOF )( 0 ), mtk::Master_Slave::MASTER );

                // create the field interpolation rule for the ith dof type group
                Interpolation_Rule tFieldInterpolationRule( static_cast< Set* >( mEquationSet )->mIPGeometryType,
                                                            Interpolation_Type::LAGRANGE,
                                                            static_cast< Set* >( mEquationSet )->mIPSpaceInterpolationOrder,
                                                            static_cast< Set* >( mEquationSet )->get_auto_time_interpolation_type( tNumTimeNodes ), // fixme
                                                            // If interpolation type CONSTANT, iInterpolation order is not used
                                                            static_cast< Set* >( mEquationSet )->get_auto_interpolation_order( tNumTimeNodes, mtk::Geometry_Type::LINE ) ); //fixme

                MORIS_ASSERT( mMasterFI( tDofIndex ) == nullptr, "Field_Interpolator_Manager(), Field interpolator was created previously" );

                // create a field interpolator for the dof type group
                mMasterFI( tDofIndex ) = new Field_Interpolator( mMasterDofTypes( iDOF ).size(),
                                                                 tFieldInterpolationRule,
                                                                 static_cast< Set* >( mEquationSet )->mMasterIPGeometryInterpolator,
                                                                 mMasterDofTypes( iDOF ) );
            }

            // loop over the master dof type groups
            for( uint iDOF = 0; iDOF < mSlaveDofTypes.size(); iDOF++ )
            {
                // get the number of time level for the dof type group
                uint tNumTimeNodes = aModelSolverInterface->get_time_levels_for_type( mSlaveDofTypes( iDOF )( 0 ) );

                uint tDofIndex = mEquationSet->get_dof_index_for_type_1( mSlaveDofTypes( iDOF )( 0 ), mtk::Master_Slave::SLAVE );

                // create the field interpolation rule for the ith dof type group
                Interpolation_Rule tFieldInterpolationRule( static_cast< Set* > (mEquationSet )->mIPGeometryType,
                                                            Interpolation_Type::LAGRANGE,
                                                            static_cast< Set* >( mEquationSet )->mIPSpaceInterpolationOrder,
                                                            static_cast< Set* >( mEquationSet )->get_auto_time_interpolation_type( tNumTimeNodes ), // fixme
                                                            // If interpolation type CONSTANT, iInterpolation order is not used
                                                            static_cast< Set* >( mEquationSet )->get_auto_interpolation_order( tNumTimeNodes, mtk::Geometry_Type::LINE ) ); //fixme

//                MORIS_ASSERT( mSlaveFI( tDofIndex ) == nullptr, "Field_Interpolator_Manager(), Field interpolator was created previously" );

                // create a field interpolator for the dof type group
                mSlaveFI( tDofIndex ) = new Field_Interpolator( mSlaveDofTypes( iDOF ).size(),
                                                                tFieldInterpolationRule,
                                                                static_cast< Set* >( mEquationSet )->mSlaveIPGeometryInterpolator,
                                                                mSlaveDofTypes( iDOF ) );
            }
        };

//------------------------------------------------------------------------------

        moris::uint get_max_num_field_interpolators()
        {
            return mMaxNumFieldInterpolators;
        }

//------------------------------------------------------------------------------

        Field_Interpolator * get_field_interpolators_for_type( enum MSI::Dof_Type aDofType,
                                                                    mtk::Master_Slave  aIsMaster = mtk::Master_Slave::MASTER)
        {
            MORIS_ASSERT( mEquationSet != nullptr, "get_field_interpolators_for_type(), Equation Set pointer not set");

            sint tDofIndex = mEquationSet->get_dof_index_for_type_1( aDofType, aIsMaster );

            if( tDofIndex != -1 )
            {

                switch ( aIsMaster )
                {
                case ( mtk::Master_Slave::MASTER ):
                {
                    MORIS_ASSERT( (sint)mMasterFI.size() > tDofIndex, "Field_Interpolator_Manager::get_field_interpolators_for_type(), field interpolator does not exist" );
//                    MORIS_ASSERT( mMasterFI( tDofIndex ) != nullptr, "Field_Interpolator_Manager::get_field_interpolators_for_type(), field interpolator does not exist" );

                    return mMasterFI( tDofIndex );
                }
                case( mtk::Master_Slave::SLAVE ):
                {
                    MORIS_ASSERT( (sint)mSlaveFI.size() > tDofIndex, "Field_Interpolator_Manager::get_field_interpolators_for_type(), field interpolator does not exist" );
//                    MORIS_ASSERT( mSlaveFI( tDofIndex ) != nullptr, "Field_Interpolator_Manager::get_field_interpolators_for_type(), field interpolator does not exist" );

                    return mSlaveFI( tDofIndex );
                }
                default:
                {
                    MORIS_ERROR(false, "Set::get_dof_type_map - can only be MASTER or SLAVE");
                    return 0;
                }
                }
            }
            else
            {
                return nullptr;
            }
        }

        //------------------------------------------------------------------------------

        void set_coeff()
        {

        }


        //------------------------------------------------------------------------------
    };

    //------------------------------------------------------------------------------
    } /* namespace fem */
} /* namespace moris */

#endif /* SRC_FEM_CL_FEM_FIELD_INTERPLATOR_MANAGER_HPP_ */
