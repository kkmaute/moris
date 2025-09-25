#pragma once
#include "moris_typedefs.hpp"

namespace moris
{
    /**
     * Struct to locate a parameter given four inices already given in Moris_Tree_Widget_Item.
     */
    struct Parameter_Locator {
        uint module;
        uint submodule;
        uint list;
        uint parameter;
    };
}