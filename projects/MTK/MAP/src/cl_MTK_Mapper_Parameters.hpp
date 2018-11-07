/*
 * cl_MTK_Mapper_Parameters.hpp
 *
 *  Created on: Nov 5, 2018
 *      Author: messe
 */

#ifndef PROJECTS_MTK_MAP_SRC_CL_MTK_MAPPER_PARAMETERS_HPP_
#define PROJECTS_MTK_MAP_SRC_CL_MTK_MAPPER_PARAMETERS_HPP_

#include <string>
#include <algorithm>

#include "typedefs.hpp"
#include "cl_Cell.hpp"
#include "cl_Mesh_Enums.hpp"
#include "cl_XML_Parser.hpp"
#include "cl_Param_List.hpp"

//datatype for field paramater list
typedef Param_List< boost::variant< std::string > > ParameterList;

namespace moris
{
    namespace mapper
    {
// -----------------------------------------------------------------------------

       ParameterList
       create_parameter_list()
       {
           aParameterList.insert( "input_hmr_db", "" );
           aParameterList.insert( "output_hmr_db", "" );
           aParameterList.insert( "output_exodus","" );
       }

// -----------------------------------------------------------------------------

       void
       load_parameters_from_xml(
               const std::string          & aPath,
               ParameterList & aParams )
       {
           // clear output cell
           aParams.clear();

           // create a parser object for this settings path
           XML_Parser tParser( aPath );

           // get number of fields from parser
           uint tNumberOfFields = tParser.count_keys_in_subtree( "moris.mapper", "parameters" );

           // create an empty parameter list
           ParameterList tParams =  create_parameter_list();

           Cell< std::string > tFirst;
           Cell< std::string > tSecond;

           tParser.get_keys_from_subtree( "moris.mapper", "field", f, tFirst, tSecond );

           // copy key to settings struct
           for( uint k=0; k<tFirst.size(); ++k )
           {
               // get key
               std::string tKey = tFirst( k );

               if ( tKey == "input_hmr_db" )
               {
                   aSettings.set( "input_hmr_db", std::string(  tSecond( k ) ) );
               }
               else if ( tKey == "output_hmr_db" )
               {
                   aSettings( f ).set( "output_hmr_db", std::string(  tSecond( k ) ) );
               }
               else if ( tKey == "output_exodus" )
               {
                   aSettings( f ).set( "output_exodus", std::string(  tSecond( k ) ) );
               }
           }

       }
// -----------------------------------------------------------------------------
    } /* namespace mapper */
} /* namespace moris */


#endif /* PROJECTS_MTK_MAP_SRC_CL_MTK_MAPPER_PARAMETERS_HPP_ */
