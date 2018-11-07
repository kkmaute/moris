/*
 * cl_MTK_Mapper_Field_Parameters.hpp
 *
 *  Created on: Nov 5, 2018
 *      Author: messe
 */

#ifndef PROJECTS_MTK_MAP_SRC_CL_MTK_MAPPER_FIELD_PARAMETERS_HPP_
#define PROJECTS_MTK_MAP_SRC_CL_MTK_MAPPER_FIELD_PARAMETERS_HPP_


#include <string>
#include <algorithm>

#include "typedefs.hpp"
#include "cl_Cell.hpp"
#include "cl_Mesh_Enums.hpp"
#include "cl_XML_Parser.hpp"
#include "cl_Param_List.hpp"

//datatype for field paramater list
typedef Param_List< boost::variant< std::string, EntityRank  > > FieldParameterList;

namespace moris
{
    namespace mapper
    {
// -----------------------------------------------------------------------------

       FieldParameterList
       create_field_parameter_list()
       {
           ParameterList aParameterList;

           aParameterList.insert( "label", "" );
           aParameterList.insert( "input_hdf5", "" );
           aParameterList.insert( "output_hdf5", "" );
           aParameterList.insert( "input_rank", EntityRank::INVALID );
           aParameterList.insert( "output_rank", EntityRank::INVALID );

           return FieldParameterList;
       }

// -----------------------------------------------------------------------------

       void
       load_field_parameters_from_xml(
               const std::string          & aPath,
               Cell< FieldParameterList > & aParams )
       {
           // clear output cell
           aParams.clear();

           // create a parser object for this settings path
           XML_Parser tParser( aPath );

           // get number of fields from parser
           uint tNumberOfFields = tParser.count_keys_in_subtree( "moris.mapper", "field" );

           // create an empty parameter list
           FieldParameterList tParams =  create_fields_parameter_list();

           // initialize output vector
           aParams.resize( tNumberOfFields, tParams );

           // loop over all fields
           for(  uint f=0; f<tNumberOfFields; ++f )
           {
               Cell< std::string > tFirst;
               Cell< std::string > tSecond;

               tParser.get_keys_from_subtree( "moris.mapper", "field", f, tFirst, tSecond );

               // copy key to settings struct
               for( uint k=0; k<tFirst.size(); ++k )
               {
                   // get key
                   std::string tKey = tFirst( k );

                   if ( tKey == "label" )
                   {
                       aSettings( f ).set( "label", std::string(  tSecond( k ) ) );
                   }
                   else if ( tKey == "input_hdf5" )
                   {
                       aSettings( f ).set( "input_hdf5", tSecond( k ) );

                   }
                   else if ( tKey == "output_hdf5" )
                   {
                       aSettings( f ).set( "output_hdf5", tSecond( k ) );

                   }
                   else if ( tKey == "input_rank" )
                   {
                       aSettings( f ).set( "input_rank", string_to_entity_rank( tSecond( k ) ) );

                   }
                   else if ( tKey == "output_rank" )
                   {
                       aSettings( f ).set( "output_rank", string_to_entity_rank( tSecond( k ) ) );

                   }
               }
           }
       }

// -----------------------------------------------------------------------------

       EntityRank
       string_to_entity_rank( const std::string & aString )
       {
           // convert string to lower case
           std::string tString( aString );
           std::transform( data.begin(), data.end(), data.begin(), ::tolower );

           // test string and convert it to corresponding rank
           switch ( tString )
           {
               case( "node" ) :
               {
                   return EntityRank::NODE;
                   break;
               }
               case( "edge" ) :
               {
                   return EntityRank::EDGE;
                   break;
               }
               case( "face" ) :
               {
                   return EntityRank::FACE;
                   break;
               }
               case( "element" ) :
               {
                   return EntityRank::ELEMENT;
                   break;
               }
               case( "bspline_1" ) :
               {
                   return EntityRank::BSPLINE_1;
                   break;
               }
               case( "bspline_2" ) :
               {
                   return EntityRank::BSPLINE_2;
                   break;
               }
               case( "bspline_3" ) :
               {
                   return EntityRank::BSPLINE_3;
                   break;
               }
               default :
               {
                   return EntityRank::INVALID;
                   break;
               }
           }
       }

// -----------------------------------------------------------------------------

    } /* namespace mapper */
} /* namespace moris */

#endif /* PROJECTS_MTK_MAP_SRC_CL_MTK_MAPPER_FIELD_PARAMETERS_HPP_ */
