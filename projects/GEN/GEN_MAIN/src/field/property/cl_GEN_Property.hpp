#ifndef MORIS_CL_GEN_PROPERTY_HPP_
#define MORIS_CL_GEN_PROPERTY_HPP_

#include "cl_GEN_Field.hpp"
#include "st_GEN_Property_Parameters.hpp"

namespace moris
{
    namespace ge
    {
        class Property : virtual public Field
        {
        private:
            Property_Parameters mParameters;

        public:

            /**
             * Constructor for property, needs to know about other fields that it depends on.
             *
             * @param aFieldDependencies This property's dependencies
             */
            Property(Property_Parameters aParameters);

            /**
             * Copy constructor
             *
             * @param aProperty Property to copy
             */
            Property(std::shared_ptr<Property> aProperty);

        };
    }
}

#endif /* MORIS_CL_Property_HPP_ */
