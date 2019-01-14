/*
 * fn_HMR_Exec_initialize_fields.hpp
 *
 *  Created on: Nov 13, 2018
 *      Author: messe
 */

#ifndef PROJECTS_HMR_SRC_FN_HMR_EXEC_INITIALIZE_FIELDS_HPP_
#define PROJECTS_HMR_SRC_FN_HMR_EXEC_INITIALIZE_FIELDS_HPP_

#include <string>
#include <memory>

#include "../../../HMR/src/cl_HMR.hpp"
#include "../../../HMR/src/cl_HMR_Arguments.hpp"
#include "../../../HMR/src/cl_HMR_Field.hpp"
#include "../../../HMR/src/cl_HMR_Paramfile.hpp"
#include "../../../HMR/src/HMR_Globals.hpp"
#include "../../../HMR/src/HMR_Tools.hpp"
#include "cl_Cell.hpp"
#include "typedefs.hpp"

namespace moris
{
    namespace hmr
    {
// -----------------------------------------------------------------------------

        void
        initialize_fields(
                const Arguments                  & aArguments,
                const Paramfile                  & aParamfile,
                HMR                              * aHMR,
                Cell< std::shared_ptr< Field > > & aFields )
        {
            // reset field container
            aFields.clear();

            // get number of fields from parameter file
            uint tNumberOfFields = aParamfile.get_number_of_fields();

            // loop over all requested fields
            for( uint f=0; f<tNumberOfFields; ++f )
            {
                // create field pointer
                aFields.push_back(
                        aHMR->create_field(  aParamfile.get_field_params( f ) ) );
            }

        }

// -----------------------------------------------------------------------------
    } /* namespace hmr */
} /* namespace moris */


#endif /* PROJECTS_HMR_SRC_FN_HMR_EXEC_INITIALIZE_FIELDS_HPP_ */
