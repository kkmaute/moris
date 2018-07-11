#ifndef MORIS_IOS_FN_LOGGER_HPP_
#define MORIS_IOS_FN_LOGGER_HPP_

// C++ header files.
#include <string>

namespace moris
{
namespace ios
{
    /**
     * Dash separator.
     *
     * Defined here to have uniform length across the source code.
     * @return
     */
    inline
    std::string
    dash_separator()
    {
        return "-------------------------------------------------------------------------------";
    }

}    // namespace ios
}    // namespace moris

#endif    /* MORIS_IOS_FN_LOGGER_HPP_ */
