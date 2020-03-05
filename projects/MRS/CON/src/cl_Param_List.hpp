#ifndef MORIS_CONTAINERS_CL_PARAM_LIST_HPP_
#define MORIS_CONTAINERS_CL_PARAM_LIST_HPP_

// c++ header files.
#include <cstring>
#include <map>
#include <string>

// Third-party header files.
#include <boost/variant.hpp>

// MORIS library header files.
#include "cl_Matrix.hpp"
#include "assert.hpp"
#include "core.hpp"
#include "ios.hpp"


namespace moris
{
    namespace containers
    {
        struct strcmp
        {
            bool operator()( const std::string & a, const std::string & b )
            {
                return std::strcmp( a.c_str(), b.c_str() ) < 0;
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

        std::map< std::string, Variant, moris::containers::strcmp > mParamMap;

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
        void insert( const std::string & aKey, Variant aVal )
        {
            std::string tName(aVal.type().name());
            if (!tName.compare("PKc"))
            {
                std::string tVal(boost::get<const char*>(aVal));
                aVal = tVal;
            }
            mParamMap.insert( { aKey, aVal } );
        }

        /**
         * Removes the specified key and its associated value from the map
         *
         * @param aKey the key to be erased
         */
        void erase(const std::string & aKey)
        {
            mParamMap.erase(aKey);
        }

        /**
         * @brief Sets an element to a value if it exists, otherwise an error is thrown
         *
         * @param[in] aKey Key corresponding to the mapped value that
         *            needs to be accessed
         * @param[in] aVal Value corresponding to aKey
         */
        void
        set( const std::string & aKey, Variant aVal )
        {
            auto it = mParamMap.find( aKey );

            if( it ==  mParamMap.end() )
            {

                // create error message
                std::string tError =  "The requested parameter '" + aKey + "' can not be set because it does not exist.\n";

                // throw error
                MORIS_ERROR( false, tError.c_str() );

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
        operator()( const std::string & aKey )
        {
            auto it = mParamMap.find( aKey );

            // check if parameter exists
            if( it == mParamMap.end() )
            {
                // create error message
                std::string tError =  "The requested parameter '" + aKey + "' does not exist.\n";

                // throw error
                MORIS_ERROR( false, tError.c_str() );
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
        moris::sint which( const std::string & aKey )
        {
            auto it = mParamMap.find( aKey );

            // check if parameter exists
            if( it == mParamMap.end() )
            {
                // create error message
                std::string tError =  "The requested parameter '" + aKey + "' does not exist.\n";

                // throw error
                MORIS_ERROR( false, tError.c_str() );
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
        const Key & get( const std::string & aKey )
        {
            auto it = mParamMap.find( aKey );

            // check if parameter exists
            if( it == mParamMap.end() )
            {
                // create error message
                std::string tError =  "The requested parameter '" + aKey + "' does not exist.\n";

                // throw error
                MORIS_ERROR( false, tError.c_str() );
            }

            // check if the requested type is correct
            if( boost::get< Key >( &( it->second ) ) == nullptr )
            {
                // create error message
                std::string tError =  "The parameter '" + aKey + "' was requested with an incorrect type.\n";
                MORIS_ERROR( false, tError.c_str() );
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


    //datatype for parameter lists
    typedef boost::variant< bool, sint, real, const char*, std::string, uint, std::pair< std::string, std::string > > ParameterListTypes;
    typedef Param_List< ParameterListTypes > ParameterList;

}

#endif /* MORIS_CONTAINERS_CL_PARAM_LIST_HPP_ */
