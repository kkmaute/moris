#ifndef SRC_XTK_CONNECTIVITY_HPP_
#define SRC_XTK_CONNECTIVITY_HPP_

// MORIS library header files.
//#include "cl_Mesh_Connectivity.hpp"

// MORIS project header files.
#include "core.hpp"
#include "cl_Mat.hpp" // LNA/src
#include "assert.hpp"
#include "cl_Cell.hpp" // CON/src

namespace moris
{
namespace xtk
{
    class Connectivity
    {
    public:

        /**
         * Constructor.
         *
         * This is a very simplistic wrapper around vectors of vectors, to serve as a temporary interface for enrichment strategy work.
         *
         * Examples of use:
         * @include /xtk/cl_xtk_Connectivity1.inc
         */
        Connectivity();

        /**
         * Destructor.
         */
        ~Connectivity();

        /**
         * Set Connectivity.
         *
         * @param aConn       A vector of vectors connectivity.
         */
//        void
//        set_connectivity(moris::Cell<moris::Cell<moris::lint>> & aConn);


        void
        set_connectivity(moris::Mat<moris::uint> & aConn);

        /**
         * Get number of entries in outer vector
         *
         * @return Number of entries in outer vector.
         */
        moris::lint
        num_entities();

        /**
         * Get number of entries in outer vector
         *
         * @return Number of entries in outer vector.
         */
        moris::lint
        num_entities() const;

        /**
         * Get number of entries in inner vector
         *
         * @param arow       Outer vector index.
         * @return Number of entries in inner vector.
         */
        moris::lint
        num_conn( moris::lint & arow);

        /**
         * Get number of entries in inner vector
         *
         * @param arow       Outer vector index.
         * @return Number of entries in inner vector.
         */
        moris::lint
        num_conn( moris::lint & arow) const;

        moris::uint
        num_conn( moris::uint & arow);



        /**
         * Accessor operator for entry acol inner vector of arow outer vector
         *
         * @param arow       Outer vector index.
         * @param acol       Inner vector index.
         * @return Entry value.
         */
        moris::lint
        operator()(
                moris::lint const & arow,
                moris::lint const & acol);

        /**
         * Constant accessor operator for entry acol inner vector of arow outer vector
         *
         * @param arow       Outer vector index.
         * @param acol       Inner vector index.
         * @return Entry value.
         */
        moris::lint const
        operator()(
                moris::lint const & arow,
                moris::lint const & acol) const;

        moris::lint
        operator()(
                moris::uint const & arow,
                moris::uint const & acol);

        moris::lint
        get_node(moris::uint    aRow,
                 moris::uint    aCol);

        moris::Cell<moris::uint>
        get_row(moris::uint & aRow);






    private:

        moris::Cell<moris::Cell<moris::uint>> mConn;
        moris::Mat<moris::uint> mConn2;
    };

}   // namespace mesh
}   // namespace moris

//void
//test_func();
#endif  /* SRC_XTK_CONNECTIVITY_HPP_ */
