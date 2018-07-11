
// -----------------------------------------------------------------------------

template< typename T >
moris::Arma_Sp_Mat< T >::Arma_Sp_Mat()
    : moris::Base_Arma_Mat< arma::SpMat< T > >()
{
}

// -----------------------------------------------------------------------------

template< typename T >
moris::Arma_Sp_Mat< T >::Arma_Sp_Mat(
        moris::Arma_Sp_Mat< T > const & mat )
    : moris::Base_Arma_Mat< arma::SpMat< T > >()
{
    this->mMat = mat.data();
}

// -----------------------------------------------------------------------------

// FIXME: need to resolve issue with arma::word
template< typename T >
template< typename A >
moris::Arma_Sp_Mat< T >::Arma_Sp_Mat(
        A const & X )
    : moris::Base_Arma_Mat< arma::SpMat< T > >()
{
    this->mMat=X;
}

// -----------------------------------------------------------------------------

template< typename T >
moris::Arma_Sp_Mat< T >::Arma_Sp_Mat(
        moris::size_t const & i_elems,
        moris::size_t const & j_elems )
    : moris::Base_Arma_Mat< arma::SpMat< T > >(i_elems, j_elems)
{
}

// ----------------------------------------------------------------------------

template< typename T >
moris::Arma_Sp_Mat< T >::Arma_Sp_Mat(
        moris::Mat< moris::uint > const & aRowInd,
        moris::Mat< moris::uint > const & aColInd,
        moris::Mat< T >           const & aValues )
    : moris::Base_Arma_Mat< arma::SpMat< T > >()
{
    arma::Mat< moris::uint > tLocations( 2, aRowInd.numel() );
    arma::Mat< T > tValues( aValues.numel(), 1 );

    if ( !moris::isrow( aRowInd ) )
    {
        tLocations.row( 0 ) = moris::trans( aRowInd );
    }
    else
    {
        tLocations.row( 0 ) = aRowInd.data();
    }

    if ( !moris::isrow( aColInd ) )
    {
        tLocations.row( 1 ) = moris::trans( aColInd );
    }
    else
    {
        tLocations.row( 1 ) = aColInd.data();
    }

    if ( moris::isrow( aValues ) )
    {
        tValues = moris::trans( aValues );
    }
    else
    {
        tValues = aValues.data();
    }

    arma::Mat< arma::uword > tLoc( 2, aRowInd.numel() );

    for (moris::size_t i=0;i<aRowInd.numel();++i)
       for (moris::size_t j=0;j<2;++j)
          tLoc(j,i)=tLocations(j,i);

  
    this->mMat = arma::SpMat<T>( tLoc, tValues );
}

// ----------------------------------------------------------------------------

template< typename T >
moris::Arma_Sp_Mat< T >::Arma_Sp_Mat(
        moris::Mat< moris::uint > const & aRowInd,
        moris::Mat< moris::uint > const & aColInd,
        moris::Mat< T >           const & aValues,
        moris::size_t             const & aIElems,
        moris::size_t             const & aJElems )
    : moris::Base_Arma_Mat< arma::SpMat< T > >()
{
     
     arma::Mat< moris::uint> tLocations( 2, aRowInd.numel() );
     arma::Mat< T > tValues( aValues.numel(), 1 );

     if ( !moris::isrow( aRowInd ) )
     {
         tLocations.row( 0 ) = moris::trans( aRowInd );
     }
     else
     {
         tLocations.row( 0 ) = aRowInd.data();
     }

     if ( !moris::isrow( aColInd ) )
     {
         tLocations.row( 1 ) = moris::trans( aColInd );
     }
     else
     {
         tLocations.row( 1 ) = aColInd.data();
     }

     if ( moris::isrow( aValues ) )
     {
         tValues = moris::trans( aValues );
     }
     else
     {
         tValues = aValues.data();
     }

     arma::Mat< arma::uword > tLoc( 2, aRowInd.numel() );

     for (moris::size_t i=0;i<aRowInd.numel();++i)
        for (moris::size_t j=0;j<2;++j)
           tLoc(j,i)=tLocations(j,i);

     const arma::uword tIElem = aIElems;
     const arma::uword tJElem = aJElems;
     const bool add_values = true;

     this->mMat = arma::SpMat<T>( add_values, tLoc, tValues, tIElem, tJElem );
}

// -----------------------------------------------------------------------------

template< typename T >
moris::Arma_Sp_Mat< T >::~Arma_Sp_Mat()
{
}

// ----------------------------------------------------------------------------

template< typename T >
auto
moris::Arma_Sp_Mat< T >::operator()(
        moris::size_t const & i_index,
        moris::size_t const & j_index )
-> decltype( this->mMat( i_index, j_index ) )
{
    return this->mMat( i_index, j_index );
}

// ----------------------------------------------------------------------------

template< typename T >
auto
moris::Arma_Sp_Mat< T >::operator()(
        moris::size_t const & i_index,
        moris::size_t const & j_index ) const
-> decltype( this->mMat( i_index, j_index ) )
{
    return this->mMat( i_index, j_index );
}

// ----------------------------------------------------------------------------

template< typename T >
moris::uint
moris::Arma_Sp_Mat< T >::get_nnz()
{
    return this->mMat.n_nonzero;
}

// ----------------------------------------------------------------------------

template< typename T >
void
moris::Arma_Sp_Mat< T >::clear_sparsity()
{
    this->mMat.zeros();
    return;
}

// ----------------------------------------------------------------------------

template< typename T >
void
moris::Arma_Sp_Mat< T >::compress()
{
    return;
}
