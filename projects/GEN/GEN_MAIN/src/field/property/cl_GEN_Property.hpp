/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_GEN_Property.hpp
 *
 */

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
            Cell<std::string> mDependencies;

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

            /**
             * Updates the dependencies of the property based on the given fields which the property may depend on
             * (fields may have been mapped/updated).
             *
             * @param aAllUpdatedFields All fields (this property will take the ones it needs)
             */
            void update_dependencies(Cell<std::shared_ptr<Field>> aAllUpdatedFields);

            /**
             * Gets the PDV type that this property defines.
             *
             * @return PDV type
             */
            PDV_Type get_pdv_type();

            /**
             * Gets if this property's PDV type is defined on the interpolation mesh.
             *
             * @return mesh type switch
             */
            bool is_interpolation_pdv();

            /**
             * Gets the mesh set indices where this property's PDV is defined.
             *
             * @return Mesh set indices
             */
            Matrix<DDUMat> get_pdv_mesh_set_indices();

            /**
             * Gets the mesh set names where this property's PDV is defined.
             *
             * @return Mesh set names
             */
            Cell<std::string> get_pdv_mesh_set_names();

        private:

            /**
             * Sets the dependencies of this property after they have been found by update_dependencies(). By default
             * does nothing.
             *
             * @param aDependencyFields Fields that this property depends on.
             */
            virtual void set_dependencies(Cell<std::shared_ptr<Field>> aDependencyFields);

        };
    }
}

#endif /* MORIS_CL_Property_HPP_ */

