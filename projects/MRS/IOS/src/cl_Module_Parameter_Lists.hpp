/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_Module_Parameter_Lists.hpp
 *
 */

#pragma once

#include "cl_Submodule_Parameter_Lists.hpp"

namespace moris
{
    class Module_Parameter_Lists
    {
      private:
        Parameter_List_Type                 mParameterListType;
        Vector< Submodule_Parameter_Lists > mSubmoduleParameterLists;
        uint                                mCurrentSubmoduleIndex = 0;

      public:
        /**
         * Module parameter list constructor. Correctly resizes submodule parameter lists.
         *
         * @param aParameterListType Type of parameter list, to determine how to resize the parameter lists
         */
        explicit Module_Parameter_Lists( Parameter_List_Type aParameterListType );

        /**
         * Gets the number of parameter lists being stored.
         *
         * @return Parameter list vector size.
         */
        [[nodiscard]] uint size() const;

        /**
         * Clears this module parameter list, indicating it is not to be used by MORIS.
         */
        void clear();

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
            mSubmoduleParameterLists( mCurrentSubmoduleIndex ).set( aName, aValue );
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
            mSubmoduleParameterLists( mCurrentSubmoduleIndex ).set( aName, aFirstValue, aSecondValue, aMoreValues... );
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
         * Access operator to get the submodule parameter list stored at a given index.
         *
         * @tparam T Index type, may be a uint or an appropriate enum
         * @param aSubmoduleType Submodule type enum (or uint)
         * @return Submodule parameter lists
         */
        template< typename T >
        Submodule_Parameter_Lists& operator()( T aSubmoduleType )
        {
            // Set new current index
            mCurrentSubmoduleIndex = index( aSubmoduleType, std::is_enum< T >() );

            // Return submodule parameter list
            return mSubmoduleParameterLists( mCurrentSubmoduleIndex );
        }

        /**
         * Const accessor operator to get the submodule parameter list stored at a given index.
         *
         * @param aParameterListIndex Parameter list index
         * @return Const reference to stored parameter list
         */
        template< typename T >
        const Submodule_Parameter_Lists& operator()( T aSubmoduleType ) const
        {
            return mSubmoduleParameterLists( index( aSubmoduleType, std::is_enum< T >() ) );
        }

        /**
         * Gets a begin() iterator through the underlying vector of parameter lists.
         *
         * @return Beginning of the parameter list vector
         */
        [[nodiscard]] auto begin()->decltype( mSubmoduleParameterLists.begin() );

        /**
         * Gets an end() iterator through the underlying vector of parameter lists.
         *
         * @return End of the parameter list vector
         */
        [[nodiscard]] auto end()->decltype( mSubmoduleParameterLists.end() );

        /**
         * Should only be used for old input files, with old FEM parameter list setup
         */
        void hack_for_legacy_fem()
        {
            mSubmoduleParameterLists.resize( 8 );
        }

      private:

        /**
         * Gets a submodule index from an enum.
         *
         * @tparam T Enum type
         * @param SubmoduleType Submodule enum
         * @return Submodule index
         */
        template< typename T >
        static uint index( T aSubmoduleType, std::true_type )
        {
            return static_cast< uint >( aSubmoduleType );
        }

        /**
         * Overload for getting a submodule index, if the index is already a uint.
         *
         * @tparam T Uint type
         * @param SubmoduleType Submodule index
         * @return Submodule index
         */
        template< typename T >
        static uint index( T aSubmoduleType, std::false_type )
        {
            return aSubmoduleType;
        }
    };
}
