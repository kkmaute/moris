//
// Created by christopherson on 5/19/20.
//

#include "fn_GEN_create_properties.hpp"
#include "fn_GEN_create_property.hpp"
#include "fn_Parsing_Tools.hpp"

namespace moris
{
    namespace ge
    {
        Cell<std::shared_ptr<Property>> create_properties( Cell<ParameterList> aPropertyParameterLists,
                                                           Matrix<DDRMat>& aADVs,
                                                           std::shared_ptr<moris::Library_IO> aLibrary )
        {
            // Initialize
            uint tNumProperties = aPropertyParameterLists.size();
            Cell<std::shared_ptr<Property>> tProperties(tNumProperties);
            Cell<std::string> tPropertyNames(tNumProperties);
            Cell<Cell<std::string>> tPropertyDependencies(tNumProperties);
            Cell<Cell<std::shared_ptr<Property>>> tNeededProperties(tNumProperties);
            
            // Fill names, dependencies
            for (uint tPropertyIndex = 0; tPropertyIndex < tNumProperties; tPropertyIndex++)
            {
                tPropertyNames(tPropertyIndex) = aPropertyParameterLists(tPropertyIndex).get<std::string>("name");
                tPropertyDependencies(tPropertyIndex) = 
                        string_to_cell<std::string>(aPropertyParameterLists(tPropertyIndex).get<std::string>("dependencies"));
                tNeededProperties(tPropertyIndex).resize(tPropertyDependencies(tPropertyIndex).size());
            }
            
            // Build based on dependencies (this is not optimally efficient, but doesn't need to be)
            bool tBuild;
            uint tNumPropertiesLeft = tNumProperties;
            uint tLoopCount = 0;
            while (tNumPropertiesLeft > 0)
            {
                for (uint tBuildPropertyIndex = 0; tBuildPropertyIndex < tNumProperties; tBuildPropertyIndex++)
                {
                    // Check if property needs to be built
                    if (tProperties(tBuildPropertyIndex) == nullptr)
                    {
                        tBuild = true;
                        
                        // Check if dependencies are built
                        if (tPropertyDependencies(tBuildPropertyIndex).size() > 0)
                        {
                            for (uint tDependencyIndex = 0; tDependencyIndex < tPropertyDependencies(tBuildPropertyIndex).size(); tDependencyIndex++)
                            {
                                for (uint tCheckPropertyIndex = 0; tCheckPropertyIndex < tNumProperties; tCheckPropertyIndex++)
                                {
                                    if (tPropertyDependencies(tBuildPropertyIndex)(tDependencyIndex) == tPropertyNames(tCheckPropertyIndex)
                                    && tProperties(tCheckPropertyIndex) == nullptr)
                                    {
                                        tBuild = false;
                                    }
                                }
                            }
                        }

                        // Build
                        if (tBuild)
                        {
                            create_property(aPropertyParameterLists(tBuildPropertyIndex), aADVs, tNeededProperties(tBuildPropertyIndex), aLibrary);
                            tNumPropertiesLeft--;
                        }
                    }
                }
                tLoopCount++;
                MORIS_ERROR(tLoopCount <= tNumProperties, "In fn_GEN_create_properties, a circular property dependency was detected. Exiting.");
            }
            
            return tProperties;
        }
    }
}
