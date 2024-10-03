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
        Submodule_Parameter_Lists( const Vector< Parameter_List >& aParameterLists = {} );

        /**
         * Sets the type of the submodule parameter list. Can only be done once.
         *
         * @param aType Type of parameter list, as a string
         */
        void set_type( const std::string& aType );

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
    };
}
