/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_Submodule_Parameter_Lists.hpp
 *
 */

#pragma once

#include "cl_Parameter_List.hpp"
#include "GEN_Data_Types.hpp"
#include "cl_SOL_Enums.hpp"
#include "fn_PRM_OPT_Parameters.hpp"

namespace moris
{
    class Submodule_Parameter_Lists
    {
      private:
        std::string              mType;
        Vector< Parameter_List > mParameterLists;

      public:
        /**
         * Constructor, optionally pass in an already-constructed vector of parameter lists.
         *
         * @param aParameterLists Vector of constructed parameter lists
         */
        explicit Submodule_Parameter_Lists( std::string aType );

        /**
         * Gets the type of this submodule parameter list.
         *
         * @return Submodule type
         */
        const std::string& get_type();

        /**
         * Gets the number of parameter lists being stored.
         *
         * @return Parameter list vector size.
         */
        [[nodiscard]] uint size() const;

        /**
         * Adds a new parameter list to be stored in this class.
         *
         * @param aParameterList Parameter list to be pushed back into the parameter list vector
         */
        void add_parameter_list( const Parameter_List& aParameterList );

        /**
         * Erases the parameter list at the given index.
         *
         * @param aParameterListIndex Index of the parameter list
         */
        void erase( uint aParameterListIndex );

        /**
         * Adds a new parameter list to be stored in this class. The type of parameter list depends on the stored submodule type.
         * If this submodule does not support adding a new parameter list, an error will be thrown.
         */
        void add_parameter_list();

        /**
         * Adds a new optimization algorithm parameter list to be stored in this class.
         * Will throw an error if this is not an optimization algorithm submodule.
         *
         * @param aOptimizationAlgorithmType Type of optimization algorithm
         */
        void add_parameter_list( opt::Optimization_Algorithm_Type aOptimziationAlgorithmType );

        /**
         * Adds a new geometry parameter list to be stored in this class.
         * Will throw an error if this is not a geometry submodule.
         *
         * @param aGeometryType Type of geometry
         */
        //void add_parameter_list( uint aGeometryType, gen::Field_Type aFieldType );

        /**
         * Adds a new GEN property parameter list to be stored in this class.
         * Will throw an error if this is not a GEN property submodule.
         *
         * @param aFieldType Type of field
         */
        void add_parameter_list( gen::Field_Type aFieldType );

        /**
         * Adds a new linear algorithm parameter list to be stored in this class.
         * Will throw an error if this is not a linear algorithm submodule.
         *
         * @param aSolverType Type of linear solver
         */
        void add_parameter_list( sol::SolverType aSolverType );

        /**
         * Adds a new preconditioner parameter list to be stored in this class.
         * Will throw an error if this is not a preconditioner submodule.
         *
         * @param aPreconditionerType Type of preconditioner
         */
        void add_parameter_list( sol::PreconditionerType aPreconditionerType );

        /**
         * Appends another submodule parameter list to this submodule parameter list.
         *
         * @param aSubmoduleToAppend Submodule parameter lists to append
         */
        void append( const Submodule_Parameter_Lists& aSubmoduleToAppend );

        /**
         * Gets if no parameter lists are currently being stored.
         *
         * @return If the vector of parameter lists is empty
         */
        bool empty();

        /**
         * Sets a parameter to a value in the most recently added parameter list, if it exists.
         *
         * @param aName Parameter name
         * @param aValue Parameter value
         */
        template< typename T >
        void set(
                const std::string& aName,
                const T&           aValue )
        {
            // Set in last parameter list
            mParameterLists.back().set( aName, aValue );
        }

        /**
         * Sets a parameter to a moris vector in the most recently added parameter list, using a parameter pack
         *
         * @param aName Parameter name
         * @param aFirstValue First value to put into the vector
         * @param aSecondValue Second value to put into the vector
         * @param aMoreValues Parameter pack of more values
         */
        template< typename T, typename... Arg_Types >
        void set(
                const std::string& aName,
                T                  aFirstValue,
                T                  aSecondValue,
                Arg_Types...       aMoreValues )
        {
            // Set in last parameter list
            mParameterLists.back().set( aName, aFirstValue, aSecondValue, aMoreValues... );
        }

        /**
         * Inserts a design variable into the most recently added parameter list.
         *
         * @param aName Design variable name
         * @param aDesignVariable Design variable parameters
         */
        void insert(
                const std::string&     aName,
                const Design_Variable& aDesignVariable );

        /**
         * Accessor operator to get the parameter list stored at a given index.
         *
         * @param aParameterListIndex Parameter list index
         * @return Reference to stored parameter list
         */
        Parameter_List& operator()( uint aParameterListIndex );

        /**
         * Const accessor operator to get the parameter list stored at a given index.
         *
         * @param aParameterListIndex Parameter list index
         * @return Const reference to stored parameter list
         */
        const Parameter_List& operator()( uint aParameterListIndex ) const;

        /**
         * Gets a begin() iterator through the underlying vector of parameter lists.
         *
         * @return Beginning of the parameter list vector
         */
        [[nodiscard]] auto begin()->decltype( mParameterLists.begin() );

        /**
         * Gets an end() iterator through the underlying vector of parameter lists.
         *
         * @return End of the parameter list vector
         */
        [[nodiscard]] auto end()->decltype( mParameterLists.end() );

      private:
        /**
         * This function errors out if the given type does not match the stored type for adding parameter lists
         */
        void check_submodule_type( const std::string& aExpectedType );
    };
}
