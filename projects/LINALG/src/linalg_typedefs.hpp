/*
 * LINALG_typedefs.hpp
 *
 *  Created on: Aug 24, 2018
 *      Author: doble
 */

#ifndef PROJECTS_LINALG_SRC_LINALG_TYPEDEFS_HPP_
#define PROJECTS_LINALG_SRC_LINALG_TYPEDEFS_HPP_


#include "typedefs.hpp"
namespace moris
{
/*
 * Each typedef for matrix type has a string of letters which indicate the type of matrix it is.
 * Letters are used here to cut down on the length of the typedef which will be used throughout MORIS
 *
 *
 * DD - Dense Dynamic
 * ST - size_t
 * L  - lint
 * S  - sint
 * U  - uint
 * C - Comparison Type
 *  */

struct
MDynamic {};


#ifdef MORIS_USE_EIGEN
#include "Eigen/Dense"
typedef bool ncomp;     // native type of compare operators
typedef uint nint; // native integer type
typedef Eigen::Matrix<real,        Eigen::Dynamic, Eigen::Dynamic>  DDRMat; // Dense dynamic Real Mat
typedef Eigen::Matrix<size_t,      Eigen::Dynamic, Eigen::Dynamic>  DDSTMat; // Dense dynamic size_t Mat
typedef Eigen::Matrix<lint,        Eigen::Dynamic, Eigen::Dynamic>  DDLMat;     // Dense dynamic lint Mat
typedef Eigen::Matrix<luint,       Eigen::Dynamic, Eigen::Dynamic>  DDLUMat;     // Dense dynamic lint Mat
typedef Eigen::Matrix<sint,        Eigen::Dynamic, Eigen::Dynamic>  DDSMat; // Dense dynamic sint  Mat
typedef Eigen::Matrix<uint,        Eigen::Dynamic, Eigen::Dynamic>  DDUMat; // Dense dynamic uint  Mat
typedef Eigen::Matrix<cplx,        Eigen::Dynamic, Eigen::Dynamic>  DDCMat; // Dense dynamic cmplx Mat
typedef Eigen::Matrix<ncomp,       Eigen::Dynamic, Eigen::Dynamic>  DDBMat; // Dense dynamic bool Mat
typedef Eigen::Matrix<nint,        Eigen::Dynamic, Eigen::Dynamic>  DDNIMat; // Densie Dynamic Native Integer Matrix
typedef Eigen::Matrix<moris_id,    Eigen::Dynamic, Eigen::Dynamic>  IdMat; // Id Matrix
typedef Eigen::Matrix<moris_index, Eigen::Dynamic, Eigen::Dynamic>  IndexMat; // Index Matrix
typedef Eigen::Matrix<real,                3,              3>  F33RMat; // Fixed 3x3 Real Mat
typedef Eigen::Matrix<real,                3,              1>  F31RMat; // Fixed 3x1 Real Mat
typedef Eigen::Matrix<uint,                3,              1>  F31UMat; // Fixed 3x1 uint Mat


#else
#include <armadillo>
typedef arma::uword             ncomp;     // native type of compare operators
typedef arma::uword             nint;     // native type of compare operators
typedef arma::Mat< real >       DDRMat;  // Dense dynamic Real Mat
typedef arma::Mat<size_t>       DDSTMat; // Dense dynamic size_t Mat
typedef arma::Mat< lint >       DDLMat;  // Dense dynamic lint Mat
typedef arma::Mat< luint >      DDLUMat;  // Dense dynamic lint Mat
typedef arma::Mat< sint >       DDSMat;  // Dense sint size_t Mat
typedef arma::Mat< uint >       DDUMat;  // Dense uint size_t Mat
typedef arma::Mat< cplx >       DDCMat;  // Dense dynamic cmplx Mat
typedef arma::Mat< ncomp >      DDBMat;  // Dense dynamic native Mat (type that comes out of <,>,== operators
typedef arma::Mat< nint >       DDNIMat;  // Dense dynamic native int Mat
typedef arma::Mat< moris_id >   IdMat; // Id Matrix
typedef arma::Mat< moris_index> Integer_Matrix; // Index Matrix
typedef arma::Mat< real > F33RMat; // Fixed 3x3 Real Mat (for arma this is the same as DDRMat)
typedef arma::Mat< real > F31RMat; // Fixed 3x1 Real Mat (for arma this is the same as DDRMat)
typedef arma::Mat< uint > F31UMat; // Fixed 3x1 Uint Mat (for arma this is the same as DDRMat)
#endif

}
#ifdef MORIS_USE_EIGEN
namespace xtk
{
    typedef Eigen::Matrix<double,   Eigen::Dynamic, Eigen::Dynamic> Default_Matrix_Real;
    typedef Eigen::Matrix<size_t, Eigen::Dynamic, Eigen::Dynamic> Default_Matrix_Integer;
}
#else
namespace xtk
{
    typedef arma::Mat< moris::real > Default_Matrix_Real;
    typedef arma::Mat<moris::size_t> Default_Matrix_Integer;
}
#endif

/*
 * Methods to deduce a matrix type
 */
namespace moris
{
template< typename Data_Type, typename N_Cols, typename N_Rows >
struct
deduce_matrix
{

};

template <>
struct
deduce_matrix<real,MDynamic,MDynamic>
{
    typedef DDRMat Matrix_Type;
};

template <>
struct
deduce_matrix<cplx,MDynamic,MDynamic>
{
    typedef DDCMat Matrix_Type;
};



}



#endif /* PROJECTS_LINALG_SRC_LINALG_TYPEDEFS_HPP_ */
