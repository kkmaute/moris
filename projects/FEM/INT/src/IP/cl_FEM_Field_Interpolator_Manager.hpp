/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_FEM_Field_Interpolator_Manager.hpp
 *
 */

#ifndef SRC_FEM_CL_FEM_FIELD_INTERPLATOR_MANAGER_HPP_
#define SRC_FEM_CL_FEM_FIELD_INTERPLATOR_MANAGER_HPP_

#include "assert.h"
#include <cmath>

#include "moris_typedefs.hpp"    //MRS/COR/src
#include "cl_Matrix.hpp"
#include "linalg_typedefs.hpp"
#include "cl_Cell.hpp"
#include "cl_MTK_Cell.hpp"                      //MTK/src
#include "cl_MTK_Enums.hpp"                     //MTK/src
#include "cl_MTK_Enums.hpp"                    //MTK/src
#include "cl_MSI_Equation_Object.hpp"           //FEM/MSI/src
#include "cl_MSI_Equation_Set.hpp"              //FEM/MSI/src
#include "cl_MSI_Model_Solver_Interface.hpp"    //FEM/MSI/src
#include "cl_FEM_Enums.hpp"                     //FEM/INT/src
#include "cl_FEM_Node.hpp"                      //FEM/INT/src
#include "cl_FEM_IWG.hpp"                       //FEM/INT/src
#include "cl_FEM_Geometry_Interpolator.hpp"     //FEM/INT/src
#include "cl_FEM_Field_Interpolator.hpp"        //FEM/INT/src
#include "cl_MTK_Integrator.hpp"                //FEM/INT/src
#include "cl_FEM_Element_Factory.hpp"           //FEM/INT/src

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
            const moris::Cell< moris::Cell< enum MSI::Dof_Type > >& mDofTypes;

            // pointer to the equation set
            MSI::Equation_Set* mEquationSet = nullptr;

            // enum for leader or follower
            mtk::Leader_Follower mIsLeader;

            // dof type map
            moris::Matrix< DDSMat > mDofTypeMap;

            // list of field interpolators
            moris::Cell< Field_Interpolator* > mFI;

            // maximum number of field interpolators for all dof types
            moris::uint mMaxNumDofFI;

            // number of solution sets
            moris::uint mNumSolutionSets = 1;

            // dof type list for the FI manager
            const moris::Cell< moris::Cell< enum ge::PDV_Type > > mDvTypes;

            // dof type map
            moris::Matrix< DDSMat > mDvTypeMap;

            // list of field interpolators
            moris::Cell< Field_Interpolator* > mDvFI;

            // maximum number of field interpolators
            moris::uint mMaxNumDvFI;

            // field type list for the FI manager
            const moris::Cell< moris::Cell< mtk::Field_Type > > mFieldTypes;

            // field type map
            moris::Matrix< DDSMat > mFieldTypeMap;

            // list of field interpolators
            moris::Cell< Field_Interpolator* > mFieldFI;

            // maximum number of mtk::field field interpolators
            moris::uint mMaxNumFieldFI;

            // pointer to geometry interpolator for IP element
            Geometry_Interpolator* mIPGeometryInterpolator = nullptr;

            // pointer to geometry interpolator for IG element
            Geometry_Interpolator* mIGGeometryInterpolator = nullptr;

            // flag to indicate that geometry interpolators are owned
            bool mGeometryInterpolatorOwned = false;

            // cell shape used for interpolation
            mtk::CellShape mIGCellShape = mtk::CellShape::GENERAL;
            mtk::CellShape mIPCellShape = mtk::CellShape::GENERAL;

            // point to previous time slab field manager
            Field_Interpolator_Manager* mFieldManagerPrevious = nullptr;

            //------------------------------------------------------------------------------

          public:
            //------------------------------------------------------------------------------
            /**
             * constructor
             * @param[ in ] aDofTypes    a list of group of dof types
             * @param[ in ] aEquationSet a pointer to the corresponding equation set
             * @param[ in ] aIsLeader    enum for leader or follower
             */
            Field_Interpolator_Manager(
                    const moris::Cell< moris::Cell< enum MSI::Dof_Type > >& aDofTypes,
                    MSI::Equation_Set*                                      aEquationSet,
                    mtk::Leader_Follower                                       aIsLeader = mtk::Leader_Follower::LEADER );

            /**
             * constructor
             * @param[ in ] aDofTypes    a list of group of dof types
             * @param[ in ] aDvTypes     a list of group of dv types
             * @param[ in ] aEquationSet a pointer to the corresponding equation set
             * @param[ in ] aIsLeader    enum for leader or follower
             */
            Field_Interpolator_Manager(
                    const moris::Cell< moris::Cell< enum MSI::Dof_Type > >&   aDofTypes,
                    const moris::Cell< moris::Cell< enum ge::PDV_Type > >&        aDvTypes,
                    const moris::Cell< moris::Cell< enum mtk::Field_Type > >& aFieldTypes,
                    MSI::Equation_Set*                                        aEquationSet,
                    mtk::Leader_Follower                                         aIsLeader = mtk::Leader_Follower::LEADER );

            /**
             * constructor
             * @param[ in ] aDofTypes a list of group of dof types
             * @param[ in ] aEquationSet a pointer to the corresponding equation set
             * @param[ in ] aModelSolverInterface a pointer to the corresponding model solver interface
             * @param[ in ] aIsLeader    enum for leader or follower
             */
            Field_Interpolator_Manager(
                    const moris::Cell< moris::Cell< enum MSI::Dof_Type > >& aDofTypes,
                    MSI::Equation_Set*                                      aEquationSet,
                    MSI::Model_Solver_Interface*                            aModelSolverInterface,
                    mtk::Leader_Follower                                       aIsLeader = mtk::Leader_Follower::LEADER );

            //------------------------------------------------------------------------------
            /**
             * destructor
             */
            ~Field_Interpolator_Manager();

            //------------------------------------------------------------------------------
            /**
             * delete pointers
             */
            void delete_pointers();

            //------------------------------------------------------------------------------
            /**
             * create the field interpolator for the FI manager
             * @param[ in ] aModelSolverInterface  pointer to the corresponding model solver interface
             * @param[ in ] tNumSolutionSets       number of solutions sets (e.g. eigen vectors)
             */
            void create_field_interpolators(
                    MSI::Model_Solver_Interface* aModelSolverInterface,
                    uint                         tNumSolutionSets = 1 );

            //------------------------------------------------------------------------------
            /**
             * create IP and IG geometry interpolator for the FI manager
             */
            void create_geometry_interpolators();

            //------------------------------------------------------------------------------
            /**
             * get IP geometry interpolator pointer
             */
            Geometry_Interpolator*
            get_IP_geometry_interpolator()
            {
                return mIPGeometryInterpolator;
            }

            //------------------------------------------------------------------------------
            /**
             * get IG geometry interpolator pointer
             */
            Geometry_Interpolator*
            get_IG_geometry_interpolator()
            {
                return mIGGeometryInterpolator;
            }

            //------------------------------------------------------------------------------
            /**
             * get the dof field interpolators
             */
            moris::Cell< Field_Interpolator* >&
            get_dof_field_interpolators()
            {
                return mFI;
            }

            //------------------------------------------------------------------------------
            /**
             * get the dv field interpolators
             */
            moris::Cell< Field_Interpolator* >&
            get_dv_field_interpolators()
            {
                return mDvFI;
            }

            //------------------------------------------------------------------------------
            /**
             * get the maximum number of dof field interpolators on the manager
             * @param[ out ] mMaxNumDofFI the maximum number of dof FI on the manager
             */
            moris::uint
            get_max_num_field_interpolators()
            {
                return mMaxNumDofFI;
            }

            //------------------------------------------------------------------------------
            /**
             * get the number of solution sets
             * @param[ out ] mMaxNumDofFI the maximum number of dof FI on the manager
             */
            moris::uint
            get_num_set_field_interpolators()
            {
                return mNumSolutionSets;
            }

            //------------------------------------------------------------------------------
            /**
             * get the maximum number of dv field interpolators on the manager
             * @param[ out ] mMaxNumDvFI the maximum number of dv FI on the manager
             */
            moris::uint
            get_max_num_dv_field_interpolators()
            {
                return mMaxNumDvFI;
            }

            //------------------------------------------------------------------------------
            /**
             * get the field interpolator for a given dof type
             * @param[ in ] aDofType a dof type enum
             */
            Field_Interpolator* get_field_interpolators_for_type(
                    const enum MSI::Dof_Type aDofType,
                    const uint               aSolutionSetIndex = 0 );

            //------------------------------------------------------------------------------
            /**
             * get the field interpolator for a given dv type
             * @param[ in ] aDvType a dv type enum
             */
            Field_Interpolator* get_field_interpolators_for_type( enum ge::PDV_Type aDvType );

            //------------------------------------------------------------------------------
            /**
             * get the field interpolator for a given field type
             * @param[ in ] aFieldType a field type enum
             */
            Field_Interpolator* get_field_interpolators_for_type( mtk::Field_Type aFieldType );

            //------------------------------------------------------------------------------
            /**
             * set an evaluation point in space and time
             * @param[ in ] aParamPoint coordinates of an evaluation point
             */
            void set_space_time( const Matrix< DDRMat >& aParamPoint );

            //------------------------------------------------------------------------------
            /**
             * set an evaluation point in space and time
             * @param[ in ] aParamPoint coordinates of an evaluation point
             */
            void set_space_time_from_local_IG_point( const Matrix< DDRMat >& aLocalParamPoint );

            //------------------------------------------------------------------------------
            /**
             * set coefficients for field interpolator with specific dof type
             * @param[ in ] aDofType a dof type for which the coeff will be set
             * @param[ in ] aCoeff   coefficients to be set
             */
            void set_coeff_for_type(
                    const enum MSI::Dof_Type aDofType,
                    const Matrix< DDRMat >&  aCoeff,
                    const uint               aSolutionSetIndex = 0 );

            //------------------------------------------------------------------------------
            /**
             * set coefficients for field interpolator with specific dv type
             * @param[ in ] aDvType a dv type for which the coeff will be set
             * @param[ in ] aCoeff   coefficients to be set
             */
            void set_coeff_for_type(
                    enum ge::PDV_Type           aDvType,
                    const Matrix< DDRMat >& aCoeff );

            //------------------------------------------------------------------------------
            /**
             * set coefficients for field interpolator with specific field type
             * @param[ in ] aFieldType a field type for which the coeff will be set
             * @param[ in ] aCoeff     coefficients to be set
             */
            void set_coeff_for_type(
                    mtk::Field_Type    aFieldType,
                    const Matrix< DDRMat >& aCoeff );

            //------------------------------------------------------------------------------
            /**
             * sets the cell shape used for interpolation
             * @param[ in ] aCellShape cell shape, ie, rectangular, straight, general
             */
            void set_IG_cell_shape( mtk::CellShape aCellShape );

            //------------------------------------------------------------------------------
            /**
             * sets the cell shape used for interpolation
             * @param[ in ] aCellShape cell shape, ie, rectangular, straight, general
             */
            void set_IP_cell_shape( mtk::CellShape aCellShape );

            //------------------------------------------------------------------------------
            /**
             * sets the field manager for the previous time step
             */
            void
            set_field_interpolator_manager_previous( Field_Interpolator_Manager* aFieldManagerPrevious )
            {
                mFieldManagerPrevious = aFieldManagerPrevious;
            }

            //------------------------------------------------------------------------------
            /**
             * gets the field manager for the previous time step
             */
            Field_Interpolator_Manager*
            get_field_interpolator_manager_previous()
            {
                return mFieldManagerPrevious;
            }
            //------------------------------------------------------------------------------
        };

        //------------------------------------------------------------------------------
    } /* namespace fem */
} /* namespace moris */

#endif /* SRC_FEM_CL_FEM_FIELD_INTERPLATOR_MANAGER_HPP_ */
