#ifndef SRC_LINALG_CL_TENSOR_HPP_
#define SRC_LINALG_CL_TENSOR_HPP_


//C++ libraries
#include <stdlib.h>

// MORIS library header files.
#include "typedefs.hpp" // COR/src
#include "assert.hpp"
#include "cl_Mat.hpp" // LNA/src
#include "cl_TensorMapCreator.hpp" // LNA/src
#include "op_ostream.hpp" // LNA/src

// class forward declaration and specializations
namespace moris
{

/**
 * @brief Templated Tensor class.
 *
 * The Tensor class is specialized for each order, since each order
 * of tensors behaves slightly differently. Currently, the only common
 * functionality to all tensors is the copy constructors.
 *
 * All template arguments, the type, order, dimension and a boolean flag
 * for symmetry, are mandatory. For more on the Tensor, see the page:
 * @ref TensorClass.
 *
 * A 2nd order, 3D symmetric tensor of reals would be declared as:
 * @include LNA/src/cl_Tensor/cl_Tensor_23sym.inc
 */
template< typename T, int Order, int Dim, bool Sym >
class Tensor;
}


template< typename T, int Order, int Dim, bool Sym >
class moris::Tensor
{

protected:
    moris::Mat< T > mMat;                              /// underlying storage matrix
    static const moris::Mat< moris::uint > mTensorMap; /// Voigt notation map

public:

    /**
     * Tensor default destructor
     */
    ~Tensor() = default;

    // -------------------------------------------------------------------------

    /**
     * Tensor copy constructor
     *
     * @param[in] aTensor Given tensor to be copied
     */
    template< typename A >
    Tensor(
            moris::Tensor< A, Order, Dim, Sym > const & aTensor)
    {
        mMat   = aTensor.data();
    }

    // -------------------------------------------------------------------------

    /**
     * Wrapper constructor
     *
     * Allows a Tensor to be constructed from a moris::Mat. There are checks for symmetric
     * and specific order tensors. aMat should be given in Voigt notation.
     *
     * @param[in] aMat Canditate storage matrix
     */
    Tensor(
            moris::Mat< T> const & aMat )
    {
#ifndef NDEBUG
        if (
             !( (Order == 1 &&         aMat.n_rows() == Dim)                                             ||
                (Order == 2 &&  Sym && aMat.n_rows() == Dim*(Dim+1)/2 && aMat.n_cols() == 1)             ||
                (Order == 2 && !Sym && aMat.n_rows() == Dim*Dim       && aMat.n_cols() == 1)             ||
                (Order == 3 &&  Sym && aMat.n_rows() == Dim*(Dim+1)/2 && aMat.n_cols() == Dim)           ||
                (Order == 3 && !Sym && aMat.n_rows() == Dim*Dim       && aMat.n_cols() == Dim)           ||
                (Order == 4 &&  Sym && aMat.n_rows() == Dim*(Dim+1)/2 && aMat.n_cols() == Dim*(Dim+1)/2) ||
                (Order == 4 && !Sym && aMat.n_rows() == Dim*Dim       && aMat.n_cols() == Dim*Dim)
              )
        )
        {
            MORIS_ASSERT( false, "Tensor order, dimension and symmetry does not match size of given matrix size");
        }
#endif
        mMat = aMat; // after some checks
    }

    // -------------------------------------------------------------------------

    /**
     * Default construction
     *
     * Initializes the underlying storage matrix based on the order,
     * dimension and symmetry of the tensor; also defines the map based
     * on the given dimension and symmetry of the tensor.
     */
    template< bool M = Sym, typename std::enable_if<M>::type* = nullptr>
    Tensor()
    {
        switch( Order )
        {
        case 1:
            mMat.set_size(Dim, 1); // 1st order tensor cannot be symmetric
            break;
        case 2:
            mMat.set_size(Dim*(Dim+1)/2, 1);
            break;
        case 3:
            mMat.set_size(Dim*(Dim+1)/2, Dim);
            break;
        case 4:
            mMat.set_size(Dim*(Dim+1)/2, Dim*(Dim+1)/2);
            break;
        default:
            MORIS_ASSERT( false, "Tensor of given order is not currently supported!");
        }
    }

    // -------------------------------------------------------------------------

    /*
     * Specialization for unsymmetric tensor (see above)
     */
    template< bool M = Sym, typename std::enable_if<!M>::type* = nullptr>
    Tensor()
    {
        switch( Order )
        {
        case 1:
            mMat.set_size(Dim, 1);
            break;
        case 2:
            mMat.set_size(Dim*Dim, 1);
            break;
        case 3:
            mMat.set_size(Dim*Dim, Dim);
            break;
        case 4:
            mMat.set_size(Dim*Dim, Dim*Dim);
            break;
        default:
            MORIS_ASSERT( false, "Tensor of given order is not currently supported!");
        }
    }

    // -------------------------------------------------------------------------

    /**
     * Assignment operator.
     *
     * @param[in] aTensor Tensor to be assigned.
     *
     * @return Assignment.
     */
    const moris::Tensor< T, Order, Dim, Sym > &
    operator=(
            moris::Tensor< T, Order,Dim, Sym > const & aTensor )
    {
        mMat = aTensor.data();

        return *this;
    }

    // -------------------------------------------------------------------------

    template< typename A >
    const moris::Tensor< T, Order, Dim, Sym > &
    operator=(
            A const & X )
    {
        mMat = X;

        return *this;
    }

    // -------------------------------------------------------------------------

    /**
     * Returns a constant reference to the underlying matrix.
     *
     * @return Reference to the underlying matrix.
     */
    moris::Mat< T > const &
    data() const
    {
        return mMat;
    }

    // -------------------------------------------------------------------------

    /**
     * @return Order of the Tensor (typically: 2, 3 or 4)
     */
    moris::size_t
    order() const
    {
        return Order;
    }

    // -------------------------------------------------------------------------

    /**
     * @return Dimension of the Tensor
     */
    moris::size_t
    dim() const
    {
        return Dim;
    }

    // -------------------------------------------------------------------------

    /**
     * Allows external map access
     *
     * Const keeps the tensor map protected from external changes.
     *
     * @param[in] ind0 Desired row on map
     * @param[in] ind1 Desired column on map
     *
     * @return Corresponding Voigt index
     */
    moris::uint
    map(
            moris::size_t const & ind0 ,
            moris::size_t const & ind1 ) const
    {
        return mTensorMap(ind0,ind1);
    }

    // -------------------------------------------------------------------------

    /**
     * Gives mTensorMap(ind0,ind1)one-index access to the underlying matrix
     * directly.
     *
     * @param[in] ind0 Single index
     */
    auto
    operator[](
            moris::size_t const & ind0 )
    ->decltype( mMat(ind0 ) )
    {
        return mMat(ind0 ); // makes use of moris::Mat's single index access
    }

    /**
     * const version of operator[] tensors.
     */
    auto
    operator[](
            moris::size_t const & ind0 ) const
    ->decltype( mMat(ind0 ) )
    {
        return mMat(ind0 ); // makes use of moris::Mat's single index access
    }

    // -------------------------------------------------------------------------

    /**
     * Gives one-index access to a 1st-order Tensor.
     *
     * @param[in] ind0 Single index
     *
     * This allows the user to do something like val = A(0) or A(1) = 5.
     * Mapping is trivial for 1st order tensors. Tensors of 1st order are only
     * wrappers.
     */
    template< bool M = (Order == 1), typename std::enable_if<M>::type* = nullptr>
    auto
    operator()(
            moris::size_t const & ind0 )
    ->decltype( mMat(ind0 ) )
    {
        return mMat(ind0 ); // makes use of moris::Mat's single index access
    }

    // -------------------------------------------------------------------------

    /**
     * const version of operator() for 1st order tensors.
     */
    template< bool M = (Order == 1), typename std::enable_if<M>::type* = nullptr>
    auto
    operator()(
            moris::size_t const & ind0 ) const
    ->decltype( mMat(ind0 ) )
    {
        return mMat(ind0 ); // makes use of moris::Mat's single index access
    }

    // -------------------------------------------------------------------------

    /**
     * Gives two-index access to a 2nd-order Tensor.
     *
     * @param[in] ind0 Row index
     * @param[in] ind1 Column index
     *
     * This allows the user to do something like val = A(0,0) or A(1,2) = 5.
     * The tensor automatically figures out which component of the underlying storage
     * matrix this corresponds to (i.e. the user need not know that the (1,2) component
     * of the tensor is stored in the (5) position in the storage matrix.
     */
    template< bool M = (Order == 2), typename std::enable_if<M>::type* = nullptr>
    auto
    operator()(
            moris::size_t const & ind0 ,
            moris::size_t const & ind1 )
    ->decltype( mMat(this->map(ind0,ind1),0) )
    {
        return mMat(this->map(ind0,ind1),0);
    }

    // -------------------------------------------------------------------------

    /**
     * const version of operator() for 2nd order tensors.
     */
    template< bool M = (Order == 2), typename std::enable_if<M>::type* = nullptr>
    auto
    operator()(
            moris::size_t const & ind0 ,
            moris::size_t const & ind1 ) const
    ->decltype( mMat(mTensorMap(ind0,ind1),0) )
    {
        return mMat(mTensorMap(ind0,ind1),0);
    }

    // -------------------------------------------------------------------------

    /**
     * Gives three-index access to a 3rd-order Tensor.
     *
     * @param[in] ind0 Row index
     * @param[in] ind1 Column index
     * @param[in] ind2 Slice index
     *
     * This allows the user to do something like val = A(0,0,0) or A(1,2,2) = 5.
     * The tensor automatically figures out which component of the underlying storage
     * matrix this corresponds to (i.e. the user need not know that the (1,2) component
     * of the tensor is stored in the (5) position in the storage matrix.
     */
    template< bool M = (Order == 3), typename std::enable_if<M>::type* = nullptr>
    auto
    operator()(
            moris::size_t const & ind0,
            moris::size_t const & ind1,
            moris::size_t const & ind2 )
    ->decltype( mMat(mTensorMap(ind0, ind1), ind2 ) )
    {
        return mMat(mTensorMap(ind0, ind1), ind2 );
    }

    // -------------------------------------------------------------------------

    /**
     * const version of operator() for 3rd order tensors.
     */
    template< bool M = (Order == 3), typename std::enable_if<M>::type* = nullptr>
    auto
    operator()(
            moris::size_t const & ind0,
            moris::size_t const & ind1,
            moris::size_t const & ind2 ) const
    ->decltype( mMat(mTensorMap(ind0, ind1), ind2 ) )
    {
        return mMat(mTensorMap(ind0, ind1), ind2 );
    }

    // -------------------------------------------------------------------------

    /**
     * Gives four-index access to a 4th-order Tensor.
     *
     * @param[in] ind0 Row index
     * @param[in] ind1 Column index
     * @param[in] ind2 Slice index
     * @param[in] ind3 Second slice index
     *
     * This allows the user to do something like val = A(0,0,1,2) or A(1,2,0,1) = 5.
     * The tensor automatically figures out which component of the underlying storage
     * matrix this corresponds to (i.e. the user need not know that the (0,0,1,1) component
     * of the tensor is stored in the (0,1) position in the storage matrix.
     */
    template< bool M = (Order == 4), typename std::enable_if<M>::type* = nullptr>
    auto
    operator()(
            moris::size_t const & ind0,
            moris::size_t const & ind1,
            moris::size_t const & ind2,
            moris::size_t const & ind3 )
    ->decltype( mMat(                           // the auto + decltype allows
                     mTensorMap(ind0, ind1),    // for both writing into the tensor
                     mTensorMap(ind2, ind3) ) ) // and accessing the value
    {
        return mMat( mTensorMap(ind0, ind1), mTensorMap(ind2, ind3) );
    }

    // -------------------------------------------------------------------------

    /**
     * const version of operator() for 4th order tensors.
     */
    template< bool M = (Order == 4), typename std::enable_if<M>::type* = nullptr>
    auto
    operator()(
            moris::size_t const & ind0,
            moris::size_t const & ind1,
            moris::size_t const & ind2,
            moris::size_t const & ind3 ) const
    ->decltype( mMat(                           // the auto + decltype allows
                     mTensorMap(ind0, ind1),    // for both writing into the tensor
                     mTensorMap(ind2, ind3) ) ) // and accessing the value
    {
        return mMat( mTensorMap(ind0, ind1), mTensorMap(ind2, ind3) );
    }

}; //ends template class declaration

// Tensor Map initialization
template< typename T, int Order, int Dim, bool Sym >
const moris::Mat< moris::uint > moris::Tensor<T, Order, Dim, Sym>::mTensorMap = moris::TensorMapCreator<Dim, Sym>::makeMap();

#endif /* SRC_LINALG_CL_TENSOR_HPP_ */
