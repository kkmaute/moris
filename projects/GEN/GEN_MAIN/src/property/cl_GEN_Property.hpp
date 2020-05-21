#ifndef MORIS_CL_GEN_PROPERTY_HPP_
#define MORIS_CL_GEN_PROPERTY_HPP_

#include "cl_GEN_Field_Base.hpp"

namespace moris
{
    namespace ge
    {
        class Property : virtual public Field
        {

        protected:
            Cell<std::shared_ptr<Property>> mPropertyDependencies;

        public:

            /**
             * Constructor for property, needs to know about other properties that it depends on
             *
             * @param aPropertyDependencies This property's dependencies
             */
            Property(Cell<std::shared_ptr<Property>> aPropertyDependencies = Cell<std::shared_ptr<Property>>(0));

        };

    }   // end ge namespace
}   // end moris namespace

#endif /* MORIS_CL_Property_HPP_ */
