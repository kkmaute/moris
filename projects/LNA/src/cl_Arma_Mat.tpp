
// ----------------------------------------------------------------------------

template< typename T >
moris::Arma_Mat< T >::Arma_Mat()
    : moris::Base_Arma_Mat< arma::Mat< T > >()
{
}

// ----------------------------------------------------------------------------

template< typename T >
moris::Arma_Mat< T >::Arma_Mat(
        moris::Arma_Mat< T > const & mat )
    : moris::Base_Arma_Mat< arma::Mat< T > >()
{
    this->mMat = mat.data();
}

// ----------------------------------------------------------------------------

template< typename T >
template< typename A >
moris::Arma_Mat< T >::Arma_Mat(
        A const & X )
    : moris::Base_Arma_Mat< arma::Mat< T > >()
{
    this->mMat = X;
}

// ----------------------------------------------------------------------------

template< typename T >
moris::Arma_Mat< T >::Arma_Mat(
        moris::size_t const & i_elems,
        moris::size_t const & j_elems )
    : moris::Base_Arma_Mat< arma::Mat< T > >(i_elems, j_elems)
{
}

// ----------------------------------------------------------------------------

template< typename T >
moris::Arma_Mat< T >::Arma_Mat(
        T*                  & array,
        moris::size_t const & i_elems,
        moris::size_t const & j_elems )
    : moris::Base_Arma_Mat< arma::Mat< T > >()
{
    this->mMat = arma::Mat< T >(array, i_elems, j_elems);
}

// ----------------------------------------------------------------------------

template< typename T >
moris::Arma_Mat< T >::Arma_Mat(
        std::initializer_list< std::initializer_list< T > > list )
    : moris::Base_Arma_Mat< arma::Mat< T > >( list )
{
}

// ----------------------------------------------------------------------------

template< typename T >
moris::Arma_Mat< T >::~Arma_Mat()
{
}

// ----------------------------------------------------------------------------

template< typename T >
auto
moris::Arma_Mat< T >::operator()(
        moris::size_t const & i_index,
        moris::size_t const & j_index )
-> decltype( this->mMat( i_index, j_index ) )
{
    return this->mMat( i_index, j_index );
}

// ----------------------------------------------------------------------------

template< typename T >
auto
moris::Arma_Mat< T >::operator()(
        moris::size_t const & i_index,
        moris::size_t const & j_index ) const
-> decltype( this->mMat( i_index, j_index ) )
{
    return this->mMat( i_index, j_index );
}

// ----------------------------------------------------------------------------

template< typename T >
auto
moris::Arma_Mat< T >::operator()(
        moris::size_t const & i_index )
-> decltype( this->mMat( i_index ) )
{
    return this->mMat( i_index );
}

// ----------------------------------------------------------------------------

template< typename T >
auto
moris::Arma_Mat< T >::operator()(
        moris::size_t const & i_index ) const
-> decltype( this->mMat( i_index ) )
{
    return this->mMat( i_index );
}

// ----------------------------------------------------------------------------

template< typename T >
auto
moris::Arma_Mat< T >::operator()(
        std::pair< moris::size_t, moris::size_t > const & i,
        std::pair< moris::size_t, moris::size_t > const & j )
-> decltype( this->mMat( arma::span( i.first, i.second ), arma::span( j.first, j.second ) ) )
{
    return this->mMat( arma::span( i.first, i.second ), arma::span( j.first, j.second ) );
}

// ----------------------------------------------------------------------------

template< typename T >
void
moris::Arma_Mat< T >::fill(
        T const & aVal )
{
    this->mMat.fill( aVal );
}

// ----------------------------------------------------------------------------

template< typename T >
auto
moris::Arma_Mat< T >::eye(
            moris::size_t const & aNumElems )
-> decltype ( this->mMat.eye(aNumElems, aNumElems) )
{
    return this->mMat.eye(aNumElems, aNumElems);
}

// ----------------------------------------------------------------------------

template< typename T >
auto
moris::Arma_Mat< T >::resize(
        moris::size_t const & aNumRows,
        moris::size_t const & aNumCols )
-> decltype( this->mMat.resize( aNumRows, aNumCols ) )
{
    return this->mMat.resize( aNumRows, aNumCols );
}

// ----------------------------------------------------------------------------

template< typename T >
auto
moris::Arma_Mat< T >::copy_size(
        moris::Arma_Mat< T > const & aMat )
-> decltype( this->mMat.copy_size( aMat.data() ) )
{
    return this->mMat.copy_size( aMat.data() );
}

// ----------------------------------------------------------------------------

template< typename T >
auto
moris::Arma_Mat< T >::min( ) const
-> decltype( this->mMat.min( ) )
{
    return this->mMat.min( );
}

// ----------------------------------------------------------------------------

template< typename T >
auto
moris::Arma_Mat< T >::min(
        moris::uint & aRowIndex,
        moris::uint & aColIndex ) const
-> decltype( this->mMat.min( (arma::uword&) aRowIndex, (arma::uword&)aColIndex ) )
{
    arma::uword rowIndex;
    arma::uword colIndex;

    auto retValue = this->mMat.min( rowIndex, colIndex );

    aRowIndex = rowIndex;
    aColIndex = colIndex;

    return retValue;
}

// ----------------------------------------------------------------------------

template< typename T >
auto
moris::Arma_Mat< T >::max( ) const
-> decltype( this->mMat.max( ) )
{
    return this->mMat.max( );
}

// ----------------------------------------------------------------------------

template< typename T >
auto
moris::Arma_Mat< T >::max(
        moris::uint & aRowIndex,
        moris::uint & aColIndex ) const
-> decltype( this->mMat.max( (arma::uword&) aRowIndex, (arma::uword&) aColIndex ) )
{
    arma::uword rowIndex;
    arma::uword colIndex;

    auto retValue =  this->mMat.max( rowIndex, colIndex );

    aRowIndex = rowIndex;
    aColIndex = colIndex;

    return retValue;
}

// ----------------------------------------------------------------------------

