
template< typename T >
moris::Eigen_Mat< T >::Eigen_Mat()
    : moris::Base_Eigen_Mat< Eigen::Matrix< T, Eigen::Dynamic, Eigen::Dynamic > >()
{
}

// ----------------------------------------------------------------------------

template< typename T >
moris::Eigen_Mat< T >::Eigen_Mat(
        moris::Eigen_Mat< T > const & mat )
    : moris::Base_Eigen_Mat< Eigen::Matrix< T, Eigen::Dynamic, Eigen::Dynamic  > >()
{
    this->mMat = mat.data();
}

// ----------------------------------------------------------------------------

template< typename T >
template< typename A >
moris::Eigen_Mat< T >::Eigen_Mat(
        A const & X )
    : moris::Base_Eigen_Mat< Eigen::Matrix< T, Eigen::Dynamic, Eigen::Dynamic  > >()
{
    this->mMat = X;
}

// ----------------------------------------------------------------------------

template< typename T >
moris::Eigen_Mat< T >::Eigen_Mat(
        moris::size_t const & i_elems,
        moris::size_t const & j_elems )
    : moris::Base_Eigen_Mat< Eigen::Matrix< T, Eigen::Dynamic, Eigen::Dynamic > >(i_elems, j_elems)
{
}

// ----------------------------------------------------------------------------

template< typename T >
moris::Eigen_Mat< T >::Eigen_Mat(
            T*                  & array,
            moris::size_t const & i_elems,
            moris::size_t const & j_elems )
    : moris::Base_Eigen_Mat< Eigen::Matrix< T, Eigen::Dynamic, Eigen::Dynamic > >()
{
    this->mMat = Eigen::Map< Eigen::Matrix< T, Eigen::Dynamic, Eigen::Dynamic > >(array, i_elems, j_elems);
}

// ----------------------------------------------------------------------------

template< typename T >
moris::Eigen_Mat< T >::Eigen_Mat(
        std::initializer_list< std::initializer_list< T > > list )
    : moris::Base_Eigen_Mat< Eigen::Matrix< T, Eigen::Dynamic, Eigen::Dynamic > >( list )
{
}

// ----------------------------------------------------------------------------

template< typename T >
moris::Eigen_Mat< T >::~Eigen_Mat()
{
}

// ----------------------------------------------------------------------------

template< typename T >
auto
moris::Eigen_Mat< T >::operator()(
        moris::size_t const & i_index,
        moris::size_t const & j_index )
-> decltype( this->mMat( i_index, j_index ) )
{
    return this->mMat( i_index, j_index );
}

// ----------------------------------------------------------------------------

template< typename T >
auto
moris::Eigen_Mat< T >::operator()(
        moris::size_t const & i_index,
        moris::size_t const & j_index ) const
-> decltype( this->mMat( i_index, j_index ) )
{
    return this->mMat( i_index, j_index );
}

// ----------------------------------------------------------------------------

template< typename T >
auto
moris::Eigen_Mat< T >::operator()(
        moris::size_t const & i_index )
-> decltype( this->mMat( i_index ) )
{
    return this->mMat( i_index );
}

// ----------------------------------------------------------------------------

template< typename T >
auto
moris::Eigen_Mat< T >::operator()(
        moris::size_t const & i_index ) const
-> decltype( this->mMat( i_index ) )
{
    return this->mMat( i_index );
}

// ----------------------------------------------------------------------------

template< typename T >
auto
moris::Eigen_Mat< T >::operator()(
        std::pair< moris::size_t, moris::size_t > const & i,
        std::pair< moris::size_t, moris::size_t > const & j )
-> decltype( this->mMat.block( i.first, j.first, i.second-i.first+1,j.second-j.first+1 ) )
{
    return this->mMat.block( i.first, j.first, i.second-i.first+1,j.second-j.first+1 );
}

// ----------------------------------------------------------------------------

template< typename T >
void
moris::Eigen_Mat< T >::fill(
        T const & aVal )
{
    this->mMat.fill( aVal );
}

// ----------------------------------------------------------------------------

template< typename T >
auto
moris::Eigen_Mat<T>::eye(
            moris::size_t const & aNumElems )
-> decltype ( this->mMat.setIdentity(aNumElems, aNumElems) )
{
    return this->mMat.setIdentity(aNumElems, aNumElems);
}

// ----------------------------------------------------------------------------

template< typename T >
auto
moris::Eigen_Mat< T >::resize(
        moris::size_t const & aNumRows,
        moris::size_t const & aNumCols )
-> decltype( this->mMat.conservativeResize( aNumRows, aNumCols ) )
{
    return this->mMat.conservativeResize( aNumRows, aNumCols );
}

// ----------------------------------------------------------------------------

template< typename T >
auto
moris::Eigen_Mat< T >::copy_size(
        moris::Eigen_Mat< T > const & aMat)
-> decltype( this->mMat.resizeLike( aMat.data() ) )
{
    return this->mMat.resizeLike( aMat.data() );
}

// ----------------------------------------------------------------------------

template< typename T >
auto
moris::Eigen_Mat< T >::min( ) const
-> decltype( this->mMat.minCoeff( ) )
{
    return this->mMat.minCoeff( );
}

// ----------------------------------------------------------------------------

template< typename T >
auto
moris::Eigen_Mat< T >::min(
        moris::uint & aRowIndex,
        moris::uint & aColIndex ) const
-> decltype( this->mMat.minCoeff( ) ) 
{
    moris::uint row;
    moris::uint col;
    auto val = this->mMat.minCoeff( &row, &col );

    aRowIndex = row;
    aColIndex = col;
    return val;
}

// ----------------------------------------------------------------------------

template< typename T >
auto
moris::Eigen_Mat< T >::max( ) const
-> decltype( this->mMat.maxCoeff( ) )
{
    return this->mMat.maxCoeff( );
}

// ----------------------------------------------------------------------------

template< typename T >
auto
moris::Eigen_Mat< T >::max(
        moris::uint & aRowIndex,
        moris::uint & aColIndex ) const
-> decltype( this->mMat.maxCoeff( ) )
{
    moris::uint row;
    moris::uint col;
    auto val = this->mMat.maxCoeff( &row, &col );

    aRowIndex = row;
    aColIndex = col;
    return val;
}

// ----------------------------------------------------------------------------

template< typename T >
moris::real
moris::Eigen_Mat< T >::norm()
{
    return this->mMat.norm();
}
