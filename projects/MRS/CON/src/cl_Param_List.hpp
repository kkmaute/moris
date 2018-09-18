#ifndef MORIS_CONTAINERS_CL_PARAM_LIST_HPP_
#define MORIS_CONTAINERS_CL_PARAM_LIST_HPP_

// c++ header files.
#include <cstring>
#include <map>

// Third-party header files.
#include <boost/variant.hpp>

// MORIS library header files.
#include "assert.hpp"
#include "core.hpp"
#include "ios.hpp"

namespace moris
{
    namespace containers
    {
        struct strcmp
        {
            bool operator()( char const *a, char const *b )
            {
                return std::strcmp(a, b) < 0;
            }
        };

    }   // namespace containers
}       // namespace moris

namespace moris
{
    /**
     * The parameter list class.
     *
     * A parameter list can be created as follows:
     * @include CON/src/cl_param_list/cl_param_list_insert.inc
     */
    template< typename Variant >
    class Param_List
    {
    private:

        std::map< const char*, Variant, moris::containers::strcmp > mParamMap;

    public:

        /**
         * Constructor
         */
        Param_List() = default;

        /**
         * Destructor
         */
        ~Param_List() = default;

        /**
         * @brief Extends an existing container by the element inserted.
         *
         * @param[in] aKey Key corresponding to the mapped value that
         *            needs to be accessed
         * @param[in] aVal Value corresponding to aKey
         */
        void insert( const char* aKey, Variant aVal )
        {
            mParamMap.insert( { aKey, aVal } );
        }

        /**
         * @brief Sets an element to a value if it exists, otherwise insets it
         *
         * @param[in] aKey Key corresponding to the mapped value that
         *            needs to be accessed
         * @param[in] aVal Value corresponding to aKey
         */
        void
        set( const char* aKey, Variant aVal )
        {
            auto it = mParamMap.find( aKey );

            if( it ==  mParamMap.end() )
            {
                mParamMap.insert( { aKey, aVal } );
            }
            else
            {
                it->second  = aVal;
            }
        }

        /**
         * @brief Access operator.
         *
         * @param[in] aKey Key corresponding to the mapped value that
         *            needs to be accessed
         *
         * The mapped values can be accessed directly using this
         * operator. If the input parameter matches the key of an
         * element in the container, the function returns a reference to
         * the mapped value else an error is encountered.
         *
         * @include CON/src/cl_param_list/cl_param_list_access.inc
         */
        Variant &
        operator()( const char* aKey )
        {
            auto it = mParamMap.find( aKey );

            // check if parameter exists
            if( it == mParamMap.end() )
            {
                MORIS_LOG_ERROR << "The requested parameter does not exist.\n";
                moris::assert::error( "In cl_Param_List.hpp" );
            }

            return it->second;
        }

        /**
         * @brief Determine the underlying type of the entry.
         *
         * @param[in] aKey Key corresponding to the mapped value that
         *            needs to be accessed
         *
         * @return The index of the type as ordered in the variant entry
         *         of mParamMap
         */
        moris::sint which( const char* aKey )
        {
            auto it = mParamMap.find( aKey );

            // check if parameter exists
            if( it == mParamMap.end() )
            {
                MORIS_LOG_ERROR << "The requested parameter does not exist.\n";
                moris::assert::error( "In cl_Param_List.hpp" );
            }

            return ( it->second ).which();
        }

        /**
         * @brief Get the value corresponding to a given key.
         *
         * @param[in] aKey Key corresponding to the mapped value that
         *            needs to be accessed
         *
         * @return The value corresponding to aKey.
         */
        template< typename Key >
        const Key & get( const char* aKey )
        {
            auto it = mParamMap.find( aKey );

            // check if parameter exists
            if( it == mParamMap.end() )
            {
                MORIS_LOG_ERROR << "The requested parameter does not exist.\n";
                moris::assert::error( "In cl_Param_List.hpp" );
            }

            // check if the requested type is correct
            if( boost::get< Key >( &( it->second ) ) == nullptr )
            {
                MORIS_LOG_ERROR << "cl_Param_List: Incorrect type for requested parameter.\n";
                moris::assert::error( "In cl_Param_List.hpp" );
            }

            return boost::get< Key >( it->second );
        }

        /**
         * @brief Get the beginning iterator of the underlying map.
         *
         * @return An iterator pointing to the beginning of mParamMap.
         */
        auto
        begin()
        -> decltype( mParamMap.begin() )
        {
            return mParamMap.begin();
        }

        /**
         * @brief Get the end iterator of the underlying map.
         *
         * @return An iterator pointing to the end of mParamMap.
         */
        auto
        end()
        -> decltype( mParamMap.end() )
        {
            return mParamMap.end();
        }
    };
}

#endif /* MORIS_CONTAINERS_CL_PARAM_LIST_HPP_ */
