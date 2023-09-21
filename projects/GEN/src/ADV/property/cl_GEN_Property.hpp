/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_GEN_Property.hpp
 *
 */

#pragma once

#include "cl_GEN_Design_Field.hpp"
#include "GEN_Data_Types.hpp"

namespace moris::ge
{
    /**
     * This struct contains additional parameters that are used by properties.
     */
    struct Property_Parameters : public Field_Parameters
    {
        Cell<std::string> mDependencyNames = {};  //! Names of the dependencies of this property
        PDV_Type mPDVType = PDV_Type::UNDEFINED;  //! The type of PDV that this property will be assigned to
        bool mInterpolationPDV = true;            //! If the PDV is defined on the interpolation mesh (always true for now)
        Matrix<DDUMat> mPDVMeshSetIndices = {{}}; //! Mesh set indices for assigning PDVs
        Cell<std::string> mPDVMeshSetNames = {};  //! Mesh set names for assigning PDVs
    };

    class Property : public Design_Field
    {
    private:
        Property_Parameters mParameters;

    public:

        /**
         * Constructor taking in a field pointer and a set of parameters.
         *
         * @param aField Field for computing nodal values
         * @param aParameters Field parameters
         */
        Property(
              std::shared_ptr< Field > aField,
              Property_Parameters      aParameters = {} );

        /**
         * Updates the dependencies of the property based on the given fields which the property may depend on
         * (fields may have been mapped/updated).
         *
         * @param aAllUpdatedFields All fields (this property will take the ones it needs)
         */
        void update_dependencies( Cell< std::shared_ptr< Design_Field > > aAllUpdatedFields );

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
