/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_GEN_Property.cpp
 *
 */

#include "cl_GEN_Property.hpp"

namespace moris
{
    namespace ge
    {

        //--------------------------------------------------------------------------------------------------------------

        Property::Property(Property_Parameters aParameters)
                : mParameters(aParameters)
        {
        }

        //--------------------------------------------------------------------------------------------------------------

        Property::Property(std::shared_ptr<Property> aProperty)
                : mParameters(aProperty->mParameters)
        {
        }

        //--------------------------------------------------------------------------------------------------------------

        void Property::update_dependencies(Vector<std::shared_ptr<Field>> aAllUpdatedFields)
        {
            // Set up dependency fields
            uint tNumDependencies = mParameters.mDependencyNames.size();
            Vector<std::shared_ptr<Field>> tDependencyFields(tNumDependencies);

            // Grab dependencies
            for (uint tDependencyIndex = 0; tDependencyIndex < tNumDependencies; tDependencyIndex++)
            {
                for (uint tFieldIndex = 0; tFieldIndex < aAllUpdatedFields.size(); tFieldIndex++)
                {
                    if (aAllUpdatedFields(tFieldIndex)->get_name() == mParameters.mDependencyNames(tDependencyIndex))
                    {
                        tDependencyFields(tDependencyIndex) = aAllUpdatedFields(tFieldIndex);
                    }
                }
            }

            // Set dependencies
            this->set_dependencies(tDependencyFields);
        }

        //--------------------------------------------------------------------------------------------------------------

        PDV_Type Property::get_pdv_type()
        {
            return mParameters.mPDVType;
        }

        //--------------------------------------------------------------------------------------------------------------

        bool Property::is_interpolation_pdv()
        {
            return mParameters.mInterpolationPDV;
        }

        //--------------------------------------------------------------------------------------------------------------

        Matrix<DDUMat> Property::get_pdv_mesh_set_indices()
        {
            return mParameters.mPDVMeshSetIndices;
        }

        //--------------------------------------------------------------------------------------------------------------

        Vector<std::string> Property::get_pdv_mesh_set_names()
        {
            return mParameters.mPDVMeshSetNames;
        }

        //--------------------------------------------------------------------------------------------------------------

        void Property::set_dependencies(Vector<std::shared_ptr<Field>> aDependencyFields)
        {
        }

        //--------------------------------------------------------------------------------------------------------------

    }
}

