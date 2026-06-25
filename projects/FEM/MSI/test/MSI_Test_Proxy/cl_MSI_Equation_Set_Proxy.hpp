/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_MSI_Equation_Set_Proxy.hpp
 *
 */

#ifndef SRC_FEM_CL_MSI_EQUATION_BLOCK_HPP_
#define SRC_FEM_CL_MSI_EQUATION_BLOCK_HPP_

#include "assert.h"
#include "cl_Communication_Tools.hpp"    //FEM/INT/src
#include "cl_Map.hpp"
#include "cl_MTK_Enums.hpp"
#include "cl_FEM_Enums.hpp"
#include "GEN_Data_Types.hpp"
#include "cl_MSI_Equation_Set.hpp"

#include <unordered_set>
#include <unordered_map>

namespace moris
{
    namespace mtk
    {
        class Cell;
        class Set;
    }    // namespace mtk
    namespace fem
    {
        class Field;
    }

    namespace vis
    {
        enum class Output_Type;
        enum class Field_Type;
    }    // namespace vis

    namespace MSI
    {
        class Model_Solver_Interface;
        class Design_Variable_Interface;
        class Equation_Model;
        enum class Dof_Type;

        //------------------------------------------------------------------------------
        /**
         * \brief Equation Set class
         */
        class Equation_Set_Proxy : public Equation_Set
        {
          public:
            //------------------------------------------------------------------------------
            /**
             * trivial constructor
             */
            Equation_Set_Proxy() {};

            //------------------------------------------------------------------------------
            /**
             * trivial destructor
             */
            virtual ~Equation_Set_Proxy() {};

            //-------------------------------------------------------------------------------------------------
            /**
             * initialize set
             */
            virtual void
            initialize_set(
                    const bool                      aIsStaggered            = false,
                    const fem::Time_Continuity_Flag aTimeContinuityOnlyFlag = fem::Time_Continuity_Flag::DEFAULT ) final
            {
                MORIS_ERROR( false, "Equation_Set_Proxy::initialize_set - not implemented for virtual member function" );
            }

            virtual fem::Element_Type get_element_type() const final
            {
                MORIS_ERROR( false, "Equation_Set_Proxy::get_element_type - not implemented for virtual member function" );
                return fem::Element_Type::BULK;
            }

            /**
             * \brief Get the displacement for every node in the set.
             * \return A map from the node index to the displacement vector.
             */
            virtual std::unordered_map< moris_index, Vector< real > > get_nodal_displacements( const std::unordered_set< moris_index >& aRequestedNodes ) final
            {
                MORIS_ERROR( false, "Equation_Set_Proxy::get_nodal_displacements - not implemented for virtual member function" );
                return {};
            };

            virtual void update() final
            {
                MORIS_ERROR( false, "Equation_Set_Proxy::update - not implemented for virtual member function" );
            };

            //-------------------------------------------------------------------------------------------------
            /**
             * free memory
             */
            virtual void
            free_memory() final
            {
                MORIS_ERROR( false, "Equation_Set_Proxy::free_memory - not implemented for virtual member function" );
            }

            //-------------------------------------------------------------------------------------------------
            /**
             * finalize
             */
            virtual void
            finalize( MSI::Model_Solver_Interface* aModelSolverInterface ) final
            {
                MORIS_ERROR( false, "Equation_Set_Proxy::finalize - not implemented for msi base class." );
            }

            //-------------------------------------------------------------------------------------------------
            /**
             * set GEN/MSI interface
             * @param[ in ] aDesignVariableInterface a GEN/MSI interface pointer
             */
            // virtual void
            // set_dv_interface( MSI::Design_Variable_Interface* aDesignVariableInterface ) final
            // {
            //     MORIS_ERROR( false, "Equation_Set_Proxy::set_dv_interface - not implemented for msi base class." );
            // }

            /**
             * set visualization set
             * @param[ in ] aMeshIndex
             * @param[ in ] aVisMeshSet
             * @param[ in ] aOnlyPrimaryCells
             */
            virtual void
            set_visualization_set(
                    const uint       aMeshIndex,
                    moris::mtk::Set* aVisMeshSet,
                    const bool       aOnlyPrimaryCells ) final
            {
                MORIS_ASSERT( false, "Equation_Set_Proxy::set_visualization_set(), not implemented for base class" );
            }

            //------------------------------------------------------------------------------

            //-----------------------------------------------------------------------------------------
            /**
             * get number of requested IQI for SA on set
             * @param[ out ] uint number of requested IQI for SA on set
             */
            virtual uint
            get_number_of_requested_IQIs() final
            {
                MORIS_ASSERT( false, "Equation_Set_Proxy::get_number_of_requested_IQIs(), not implemented for base class" );
                return 0;
            }

            //------------------------------------------------------------------------------
            /**
             * compute quantity of interest global
             * @param[ in ] aMeshIndex   mesh index to defined IG mesh to use
             * @param[ in ] aFieldValues matrix to be filled with QI global values
             * @param[ in ] aQINames     list of QI names to compute
             */
            virtual void
            compute_quantity_of_interest_global(
                    const uint                   aMeshIndex,
                    Matrix< DDRMat >*            aFieldValues,
                    const Vector< std::string >& aQINames ) final
            {
                MORIS_ASSERT( false, "Equation_Set_Proxy::compute_quantity_of_interest_global - not implemented for base class." );
            }

            //------------------------------------------------------------------------------
            /**
             * compute quantity of interest nodal
             * @param[ in ] aMeshIndex   mesh index to defined IG mesh to use
             * @param[ in ] aFieldValues matrix to be filled with QI nodal values
             * @param[ in ] aQINames     list of QI names to compute
             */
            virtual void
            compute_quantity_of_interest_nodal(
                    const uint                   aMeshIndex,
                    Matrix< DDRMat >*            aFieldValues,
                    const Vector< std::string >& aQINames ) final
            {
                MORIS_ASSERT( false, "Equation_Set_Proxy::compute_quantity_of_interest_nodal - not implemented for base class." );
            }

            //------------------------------------------------------------------------------
            /**
             * compute quantity of interest elemental
             * @param[ in ] aMeshIndex          mesh index to defined IG mesh to use
             * @param[ in ] aFieldValues        matrix to be filled with QI elemental values
             * @param[ in ] aQINames            list of QI names to compute
             * @param[ in ] aOutputAverageValue whether the value is an average on the element, or the integrated quantity on the element
             */
            virtual void
            compute_quantity_of_interest_elemental(
                    const uint                   aMeshIndex,
                    Matrix< DDRMat >*            aFieldValues,
                    const Vector< std::string >& aQINames,
                    const bool                   aOutputAverageValue = true ) final
            {
                MORIS_ASSERT( false, "Equation_Set_Proxy::compute_quantity_of_interest_elemental - not implemented for base class." );
            }

            //------------------------------------------------------------------------------

            virtual void
            populate_fields(
                    Vector< std::shared_ptr< fem::Field > >& aFieldToPopulate,
                    Vector< std::string > const &            aFieldIQINames ) final
            {
                MORIS_ERROR( false, "populate_fields(), no child implementation." );
            }

            //------------------------------------------------------------------------------

            virtual std::string
            get_set_name() final
            {

                MORIS_ERROR( false, "get_set_name(), not implemented for base class." );
                return "";
            }

            //------------------------------------------------------------------------------
        };
        //------------------------------------------------------------------------------
    }    // namespace MSI
} /* namespace moris */

#endif /* SRC_FEM_CL_MSI_Equation_Set_Proxy_HPP_ */
