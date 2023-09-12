/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_SDF_Object.hpp
 *
 */

#ifndef PROJECTS_GEN_SDF_SRC_CL_SDF_OBJECT_HPP_
#define PROJECTS_GEN_SDF_SRC_CL_SDF_OBJECT_HPP_

#include <cl_SDF_Triangle_Vertex.hpp>
#include <string>

#include "typedefs.hpp"
#include "cl_Cell.hpp"
#include "cl_SDF_Triangle.hpp"

namespace moris
{
    namespace sdf
    {
//-------------------------------------------------------------------------------
        class Object
        {
            const real                        mMeshHighPass=1e-9;
            moris::Cell< Triangle_Vertex * >  mVertices;
            moris::Cell< Triangle * >         mTriangles;

            Matrix< DDRMat > mOffsets;

//-------------------------------------------------------------------------------
        public:
//-------------------------------------------------------------------------------

            Object ( const std::string & aFilePath,
                     Matrix< DDRMat >    aOffsets = { {0 ,0 , 0} } );

//-------------------------------------------------------------------------------

            ~Object();

//-------------------------------------------------------------------------------

            moris::Cell< Triangle * > &
            get_triangles()
            {
                return mTriangles;
            }

//-------------------------------------------------------------------------------

            moris::Cell< Triangle_Vertex * > &
            get_vertices()
            {
                return mVertices;
            }

//-------------------------------------------------------------------------------
// MTK
//-------------------------------------------------------------------------------

            Matrix< IndexMat >
            get_nodes_connected_to_element_loc_inds
                ( moris_index aElementIndex ) const;

//-------------------------------------------------------------------------------
        private:
//-------------------------------------------------------------------------------

            /**
             * loads an ascii file and creates vertex and triangle objects
             */
            void
            load_from_object_file( const std::string& aFilePath );

//-------------------------------------------------------------------------------

            /**
             * loads an ascii file and creates vertex and triangle objects
             */
            void
            load_from_stl_file( const std::string& aFilePath );

//-------------------------------------------------------------------------------

            /**
             * loads an ASCII file into a buffer of strings.
             * Called through construction.
             */
            void
            load_ascii_to_buffer( const std::string& aFilePath,
                    moris::Cell<std::string>& aBuffer);

//-------------------------------------------------------------------------------
        };
//-------------------------------------------------------------------------------
    } /* namespace sdf */
} /* namespace moris */

#endif /* PROJECTS_GEN_SDF_SRC_CL_SDF_OBJECT_HPP_ */

