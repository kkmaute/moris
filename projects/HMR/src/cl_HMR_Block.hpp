/*
 * cl_HMR_Block.hpp
 *
 *  Created on: Jul 24, 2018
 *      Author: messe
 */

#ifndef SRC_HMR_CL_HMR_BLOCK_HPP_
#define SRC_HMR_CL_HMR_BLOCK_HPP_


#include <string>
#include "typedefs.hpp" //COR/src
#include "cl_Map.hpp"
#include "cl_MTK_Vertex.hpp" //MTK/src
#include "cl_MTK_Cell.hpp" //MTK/src
#include "cl_MTK_Block.hpp" //MTK/src
namespace moris
{
    namespace hmr
    {
//-----------------------------------------------------------------------------
        class Lagrange_Mesh_Base;
//-----------------------------------------------------------------------------

        class Block  : public moris::mtk::Block
        {
            // pointer to lagrange mesh
            Lagrange_Mesh_Base* mMesh;

            //! ID of this block
            const luint mID;

            //! label telling something about the block
            std::string mLabel;

//------------------------------------------------------------------------------
        public:
//------------------------------------------------------------------------------

            /**
             * Block constructor
             */
            Block( Lagrange_Mesh_Base*  aMesh, luint aID ) :
                mMesh( aMesh ),
                mID( aID )
            {

            }

//------------------------------------------------------------------------------

            /**
             * destructor
             */
            ~Block(){};

//-----------------------------------------------------------------------------

            /**
             * return a label that describes the block
             */
            std::string
            get_label() const ;

//------------------------------------------------------------------------------

            /**
             * sets the name of a block
             */
            void
            set_label( const std::string & aLabel ) ;
//------------------------------------------------------------------------------

            /**
             * returns the number of nodes on this mesh
             */
            uint
            get_number_of_vertices() const ;

//------------------------------------------------------------------------------

            /**
             * returns the pointer to a node
             */
            mtk::Vertex *
            get_vertex_by_index( const moris_index & aIndex ) ;

            const mtk::Vertex *
            get_vertex_by_index( const moris_index & aIndex ) const ;

//------------------------------------------------------------------------------

            /**
             * returns the number of elements on the Lagrange mesh
             */
            uint
            get_number_of_cells() const ;

//------------------------------------------------------------------------------

            /**
             * returns a pointer to the element object on the mesh
             */
            mtk::Cell *
            get_cell_by_index( const moris_index & aIndex ) ;

            const mtk::Cell *
            get_cell_by_index( const moris_index & aIndex )  const;

//------------------------------------------------------------------------------

            /**
             * returns the ID of the block, which in this case, is identical
             * to the order of the mesh
             */
            moris_id
            get_id() const ;

//------------------------------------------------------------------------------

            sint
            get_number_of_adofs_used_by_proc() const;

//------------------------------------------------------------------------------

            void
            get_adof_map( map< moris_id, moris_index > & aAdofMap ) const;

 //------------------------------------------------------------------------------

            Lagrange_Mesh_Base*
            get_lagrange_mesh()
            {
                return mMesh;
            }

//------------------------------------------------------------------------------

            uint
            get_interpolation_order() const;


//------------------------------------------------------------------------------
        };

//------------------------------------------------------------------------------
    } /* namespace hmr */
} /* namespace moris */

#endif /* SRC_HMR_CL_HMR_BLOCK_HPP_ */
