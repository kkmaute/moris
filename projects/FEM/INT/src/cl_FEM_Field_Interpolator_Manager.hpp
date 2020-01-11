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

            std::cout<<"222"<<std::endl;
            // maximum number of dof field interpolators
            mMaxNumDofFI =  mEquationSet->get_num_unique_dof_types();
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
        };

//------------------------------------------------------------------------------
        /**
         * get the maximum number of field interpolators on the manager
         * @param[ out ] mMaxNumFieldInterpolators the maximum number of FI on the manager
         */
        moris::uint get_max_num_field_interpolators()
        {
            return mMaxNumDofFI;
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
        void set_coeff()
        {

        }

//------------------------------------------------------------------------------
    };

//------------------------------------------------------------------------------
    } /* namespace fem */
} /* namespace moris */

#endif /* SRC_FEM_CL_FEM_FIELD_INTERPLATOR_MANAGER_HPP_ */
