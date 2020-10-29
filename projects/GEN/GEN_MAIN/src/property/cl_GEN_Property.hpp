#ifndef MORIS_CL_GEN_PROPERTY_HPP_
#define MORIS_CL_GEN_PROPERTY_HPP_

#include "cl_GEN_Field.hpp"

namespace moris
{
    namespace ge
    {
        class Property : virtual public Field
        {

        protected:
            Cell<std::shared_ptr<Field>> mFieldDependencies;

        public:

            /**
             * Constructor for property, needs to know about other fields that it depends on.
             *
             * @param aFieldDependencies This property's dependencies
             */
            Property(Cell<std::shared_ptr<Field>> aFieldDependencies = {});

        };
    }
}

#endif /* MORIS_CL_Property_HPP_ */
