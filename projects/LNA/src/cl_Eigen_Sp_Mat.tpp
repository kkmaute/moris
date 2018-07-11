
template< typename T >
moris::Eigen_Sp_Mat< T >::Eigen_Sp_Mat()
    : moris::Base_Eigen_Mat< Eigen::SparseMatrix< T > >()
{
}

// ----------------------------------------------------------------------------

template< typename T >
moris::Eigen_Sp_Mat< T >::Eigen_Sp_Mat(
        moris::Eigen_Sp_Mat< T > const & mat )
    : moris::Base_Eigen_Mat< Eigen::SparseMatrix< T > >()
{
    this->mMat = mat.data();
}

// ----------------------------------------------------------------------------

template< typename T >
template< typename A >
moris::Eigen_Sp_Mat< T >::Eigen_Sp_Mat(
        A const & X )
    : moris::Base_Eigen_Mat< Eigen::SparseMatrix< T > >()
{
    this->mMat = X;
}

// ----------------------------------------------------------------------------

template< typename T >
moris::Eigen_Sp_Mat< T >::Eigen_Sp_Mat(
        moris::size_t const & i_elems,
        moris::size_t const & j_elems )
    : moris::Base_Eigen_Mat< Eigen::SparseMatrix< T > >(i_elems, j_elems)
{
}

// ----------------------------------------------------------------------------

template< typename T >
moris::Eigen_Sp_Mat< T >::Eigen_Sp_Mat(
        moris::Mat< moris::uint > const & aRowInd,
        moris::Mat< moris::uint > const & aColInd,
        moris::Mat< T >           const & aValues )
    : moris::Base_Eigen_Mat< Eigen::SparseMatrix< T > >()
{
    Eigen::Matrix< moris::uint, Eigen::Dynamic, Eigen::Dynamic > tLocations( aRowInd.numel(), 2 );
    Eigen::Matrix< T, Eigen::Dynamic, Eigen::Dynamic > tValues( aValues.numel(), 1 );

    if ( moris::isrow( aRowInd ) )
    {
        tLocations.col( 0 ) = moris::trans( aRowInd );
    }
    else
    {
        tLocations.col( 0 ) = aRowInd.data();
    }

    if ( moris::isrow( aColInd ) )
    {
        tLocations.col( 1 ) = moris::trans( aColInd );
    }
    else
    {
        tLocations.col( 1 ) = aColInd.data();
    }

    if ( moris::isrow( aValues ) )
    {
        tValues = moris::trans( aValues );
    }
    else
    {
        tValues = aValues.data();
    }

    typedef Eigen::Triplet< T > tT;
    std::vector< tT > tripletList;
    tripletList.reserve( aValues.numel());

    for ( moris::size_t ie = 0; ie < (moris::size_t)aValues.numel(); ++ie )
    {
        tripletList.push_back( tT(tLocations( ie, 0 ), tLocations( ie, 1 ), tValues( ie, 0 ) ) );
    }

    Eigen::SparseMatrix< T > tMat(tLocations.col( 0 ).maxCoeff() + 1, tLocations.col( 1 ).maxCoeff() + 1);
    tMat.setFromTriplets( tripletList.begin(), tripletList.end() );

    this->mMat = tMat;
}


// ----------------------------------------------------------------------------

template< typename T >
moris::Eigen_Sp_Mat< T >::Eigen_Sp_Mat(
        moris::Mat< moris::uint > const & aRowInd,
        moris::Mat< moris::uint > const & aColInd,
        moris::Mat< T >           const & aValues,
        moris::size_t                  const & aIElems,
        moris::size_t                  const & aJElems )
    : moris::Base_Eigen_Mat< Eigen::SparseMatrix< T > >( aIElems, aJElems )
{
    Eigen::Matrix< moris::uint, Eigen::Dynamic, Eigen::Dynamic > tLocations( aRowInd.numel(), 2 );
     Eigen::Matrix< T, Eigen::Dynamic, Eigen::Dynamic > tValues( aValues.numel(), 1 );

     if ( moris::isrow( aRowInd ) )
     {
         tLocations.col( 0 ) = moris::trans( aRowInd );
     }
     else
     {
         tLocations.col( 0 ) = aRowInd.data();
     }

     if ( moris::isrow( aColInd ) )
     {
         tLocations.col( 1 ) = moris::trans( aColInd );
     }
     else
     {
         tLocations.col( 1 ) = aColInd.data();
     }

     if ( moris::isrow( aValues ) )
     {
         tValues = moris::trans( aValues );
     }
     else
     {
         tValues = aValues.data();
     }

     typedef Eigen::Triplet< T > tT;
     std::vector< tT > tripletList;
     tripletList.reserve( aValues.numel());

     for ( moris::size_t ie = 0; ie < (moris::size_t)aValues.numel(); ++ie )
     {
         tripletList.push_back( tT(tLocations( ie, 0 ), tLocations( ie, 1 ), tValues( ie, 0 ) ) );
     }

     this->mMat.setFromTriplets( tripletList.begin(), tripletList.end() );
}

// ----------------------------------------------------------------------------

template< typename T >
moris::Eigen_Sp_Mat< T >::~Eigen_Sp_Mat()
{
}

// ----------------------------------------------------------------------------

template< typename T >
auto
moris::Eigen_Sp_Mat< T >::operator()(
        moris::size_t const & i_index,
        moris::size_t const & j_index )
-> decltype( this->mMat.coeffRef( i_index, j_index ) )
{
    return this->mMat.coeffRef( i_index, j_index );
}

// ----------------------------------------------------------------------------

template< typename T >
auto
moris::Eigen_Sp_Mat< T >::operator()(
        moris::size_t const & i_index,
        moris::size_t const & j_index ) const
-> decltype( this->mMat.coeff( i_index, j_index ) )
{
    return this->mMat.coeff( i_index, j_index );
}

// ----------------------------------------------------------------------------

template< typename T >
moris::uint
moris::Eigen_Sp_Mat< T >::get_nnz()
{
    return this->mMat.nonZeros();
}

// ----------------------------------------------------------------------------

template< typename T >
void
moris::Eigen_Sp_Mat< T >::clear_sparsity()
{
    this->mMat.setZero();
    return;
}

// ----------------------------------------------------------------------------

template< typename T >
void
moris::Eigen_Sp_Mat< T >::compress()
{
    this->mMat.data().squeeze();
    this->mMat.makeCompressed();
    return;
}
