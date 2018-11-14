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

#include "cl_Cell.hpp"
#include "typedefs.hpp"
#include "HMR_Globals.hpp"
#include "HMR_Tools.hpp"
#include "cl_HMR.hpp"
#include "cl_HMR_Arguments.hpp"
#include "cl_HMR_Field.hpp"

namespace moris
{
    namespace hmr
    {
// -----------------------------------------------------------------------------

        void
        initialize_fields(
                const Arguments                  & aArguments,
                Cell< ParameterList >            & aFieldParameters,
                HMR                              * aHMR,
                Cell< std::shared_ptr< Field > > & aFields )
        {
            // list with files
            ParameterList tFileList;

            // load the passed files from the parameter list
            load_file_list_from_xml(
                    aArguments.get_parameter_path(), tFileList );

            // path for hdf input
            // std::string tHDF5Path = tFileList.get< std::string >("input_field_database");

            // path to input exodus file
            std::string tInputExodus = tFileList.get< std::string >("input_exodus");

            // cell that contains field specific parameters
            Cell< ParameterList > tFieldParams;

            // initialize field parameters
            (
                    aArguments.get_parameter_path(),
                    tFieldParams );


        }

// -----------------------------------------------------------------------------
    } /* namespace hmr */
} /* namespace moris */


#endif /* PROJECTS_HMR_SRC_FN_HMR_EXEC_INITIALIZE_FIELDS_HPP_ */
