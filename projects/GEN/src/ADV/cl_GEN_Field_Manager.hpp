/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_GEN_Field_Manager.hpp
 *
 */

#include "cl_GEN_ADV_Manager.hpp"
#include "cl_GEN_Field.hpp"

namespace moris::ge
{
    class Field_Manager
    {
      private:
        Cell< std::shared_ptr< Field > > mFields;

        
    };
}    // namespace moris::ge
