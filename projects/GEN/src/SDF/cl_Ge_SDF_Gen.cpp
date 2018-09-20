/*
 * ge_SDF_Gen.cpp
 *
 *  Created on: Mar 6, 2018
 *      Author: messe
 */

#include "cl_Ge_SDF_Gen.hpp"
#include "fn_unique.hpp"

// -----------------------------------------------------------------------------

ge::SDF_Gen::SDF_Gen(
        const moris::database                    &aBackgroundMesh, // !< a wrapper of the mesh we are working with
        const moris::Cell<std::string>           &aFilePaths,      // !< file paths
        moris::Cell< moris::Matrix< moris::DDRMat > >  &aObjectSDFs,      // !< individual SDFs for each object
        moris::Cell< moris::BoostBitset >        &aObjectSDFFlags      // !< individual SDFs for each object
        ):
        mCore(), mMeshData(aBackgroundMesh), mNumberOfObjects(aFilePaths.size())
{
    // reset objects
    mObjects.clear();
    mObjects.reserve(mNumberOfObjects);
    aObjectSDFs.clear();
    aObjectSDFs.reserve(mNumberOfObjects);
    moris::Matrix< moris::DDRMat > tEmptyMatrix;

    aObjectSDFFlags.clear();
    aObjectSDFFlags.reserve(mNumberOfObjects);
    moris::BoostBitset tEmptyBitset;

    // create objects
    for( moris::uint k=0; k < mNumberOfObjects; ++k)
    {

        // we load the triangles from the file
        ge::SDF_Triangle_File tTriangleFile( aFilePaths( k ) );

        // create an empty SDF matrix in output vector
        aObjectSDFs.push_back(tEmptyMatrix);

        // add empty bitset to Cell
        aObjectSDFFlags.push_back(tEmptyBitset);

        // create the data object
        ge::SDF_Data tData(aObjectSDFs( k ),
                           aObjectSDFFlags( k ),
                           tTriangleFile.get_triangles(),
                           tTriangleFile.get_node_coords());

        // save object file in list
        mObjects.push_back(tData);

    }

}

// -----------------------------------------------------------------------------

void
ge::SDF_Gen::set_candidate_search_depth(const moris::uint aCandidateSearchDepth)
{
    for( moris::uint k=0; k < mNumberOfObjects; ++k)
    {
        mObjects( k ).mSettings.mCandidateSearchDepth = aCandidateSearchDepth;
    }
}

// -----------------------------------------------------------------------------

void
ge::SDF_Gen::calculate_raycast()
{
    // update mesh data
    mMeshData.update();

    // calculate raycast for each object
    for( moris::uint k=0; k < mNumberOfObjects; ++k)
    {
        mCore.calculate_raycast(mMeshData, mObjects( k ));
        std::cout << " SDF : Inside Nodes in Object " << k << " " << mObjects( k ).get_number_of_inside_nodes() << std::endl;
    }
}

// -----------------------------------------------------------------------------

void
ge::SDF_Gen::calculate_raycast_and_first_sdf()
{
    // update mesh data
    mMeshData.update();

    // calculate SDF of first entity
    mCore.calculate_raycast_and_sdf(mMeshData, mObjects( 0 ));

    // only raycast rest
    if (mNumberOfObjects > 1)
        {
        // calculate raycast for each object
        for( moris::uint k=1; k < mNumberOfObjects; ++k)
        {
            mCore.calculate_raycast(mMeshData, mObjects( k ));
        }
    }
}

// -----------------------------------------------------------------------------

void
ge::SDF_Gen::calculate_raycast_and_sdf()
{
    // update mesh data
    mMeshData.update();

    // loop over all objects
    for( moris::uint k=0; k < mNumberOfObjects; ++k)
    {
        // calculate ray cast and sdf for this object
        mCore.calculate_raycast_and_sdf( mMeshData, mObjects( k ) );

        // create vtk file for debugging
        //std::cout << "writing SDF VTK debugging file" << std::endl;
        //std::string tFilePath = "sdf_" +  std::to_string(k) + ".vtk";
        //mCore.save_to_vtk( mMeshData, mObjects( k ), tFilePath );
    }
}

// -----------------------------------------------------------------------------

void
ge::SDF_Gen::calculate_raycast_and_sdf( const uint aObject )
{
    // update mesh data
    mMeshData.update();

    // calculate ray cast and sdf for this object
    mCore.calculate_raycast_and_sdf( mMeshData, mObjects( aObject ) );
}
// -----------------------------------------------------------------------------
/**
* @brief Returns moris::Mat<uint> containing all elements at the surface
*
*/
moris::Matrix< moris::DDUMat >
ge::SDF_Gen::get_elements_at_surface()
{

    // count max size of output matrix
    moris::uint tMemory = 0;
    for( moris::uint k=0; k < mNumberOfObjects; ++k )
    {
        tMemory += mObjects( k ).mLocalElementsAtSurface.length();
    }

    // output matrix
    moris::Matrix< moris::DDUMat > tElementsAtSurface(tMemory, 1);

    // counter
    moris::uint tCount = 0;
    // loop over all objects
    for( moris::uint k=0; k < mNumberOfObjects; ++k)
    {
        for( moris::uint i=0; i<mObjects( k ).mLocalElementsAtSurface.length(); ++i )
        {
            // save global element ID in vector
            tElementsAtSurface(tCount) = mMeshData.get_global_element_id(
                    mObjects( k ).mLocalElementsAtSurface( i ));
            ++tCount;

        }
    }

    // make result unique
    moris::Matrix< moris::DDUMat > tUnique;
    moris::unique(tElementsAtSurface,tUnique);
    return tUnique;
}

// -----------------------------------------------------------------------------

moris::Matrix< moris::DDUMat >
ge::SDF_Gen::get_elements_at_surface( const moris::uint aObject ) const
{
    // get number of element in surface
    moris::uint tNumberOfElements = mObjects( aObject ).mLocalElementsAtSurface.length();

    // initialize memory for output
    moris::Matrix< moris::DDUMat > aOutput( tNumberOfElements, 1 );

    // loop over all elements
    for ( moris::uint e=0; e<tNumberOfElements; ++e)
    {
        aOutput( e ) = mMeshData.get_global_element_id(
                mObjects( aObject ).mLocalElementsAtSurface( e ) );
    }

    std::fprintf( stdout, "SDF Object %u : Number of elements at surface %u.\n", aObject, tNumberOfElements );
    return  aOutput;
}

// -----------------------------------------------------------------------------

moris::Matrix< moris::DDUMat >
ge::SDF_Gen::get_elements_in_volume( const moris::uint aObject ) const
{
    // get number of element in surface
    moris::uint tNumberOfElements = mObjects( aObject ).mLocalElementsInVolume.length();

    // initialize memory for output
    moris::Matrix< moris::DDUMat > aOutput( tNumberOfElements, 1 );

    // loop over all elements
    for ( moris::uint e=0; e<tNumberOfElements; ++e)
    {
        aOutput( e ) = mMeshData.get_global_element_id(
                mObjects( aObject ).mLocalElementsInVolume( e ) );
    }

    std::fprintf( stdout, "SDF Object %u : Number of elements in volume %u.\n", aObject, tNumberOfElements );

    return  aOutput;
}

// -----------------------------------------------------------------------------

moris::Matrix< moris::DDRMat >
ge::SDF_Gen::get_inside_outside_sign_from_raycast( const moris::uint aObject ) const
{

    // get number of nodes
    moris::uint tNumberOfNodes = mMeshData.get_number_of_nodes();

    // assign memory for output
    moris::Matrix< moris::DDRMat > aSigns( tNumberOfNodes, 1, 1 );

    // loop over all nodes
    for ( uint i=0; i<tNumberOfNodes; ++i )
    {
        // test if node is inside or outside
        if ( mObjects( aObject ).mLocalNodeInsideFlags.test( i ) )
        {
            // flip sign of this node
            aSigns( i ) = -1;
        }
    }

    // return values as output
    return aSigns;
}

// -----------------------------------------------------------------------------

moris::Matrix< moris::DDUMat >
ge::SDF_Gen::get_elements_at_surface_for_hmr_ref()
{
    if ( mNumberOfObjects > 1)
    {
        // get maximum element id
        uint tMaxElemID = mMeshData.get_max_element_id() ;
        // create Boost Bitset
        moris::BoostBitset tElements( tMaxElemID + 1 );

        // loop over all objects (including first one)
        for( moris::uint k=0; k < mNumberOfObjects; ++k)
        {
            for( moris::uint i=0; i<mObjects( k ).mLocalElementsAtSurface.length(); ++i)
            {

                // activate all surface elements
                tElements.set(
                        mMeshData.get_global_element_id( mObjects( k ).mLocalElementsAtSurface( i )) );

            }
        }

        for( moris::uint k=1; k < mNumberOfObjects; ++k)
        {
            for( moris::uint i=0; i<mObjects( k ).mLocalElementsInVolume.length(); ++i)
            {
                // deactivate all volume elements
                tElements.reset(
                        mMeshData.get_global_element_id( mObjects( k ).mLocalElementsInVolume( i )) );

            }
        }

        return ge::BitsetToMat( tElements );
    }
    else
    {
        // return empty matrix
        moris::Matrix< moris::DDUMat > tEmptyMatrix;
        return tEmptyMatrix;
    }
}

// -----------------------------------------------------------------------------
moris::Matrix< moris::DDUMat >
ge::SDF_Gen::get_elements_in_volume(){
    // count max size of output matrix
    moris::uint tMemory = 0;
    for( moris::uint k=0; k < mNumberOfObjects; ++k)
    {
        tMemory += mObjects( k ).mLocalElementsInVolume.length();
    }

    // output matrix
    moris::Matrix< moris::DDUMat > tElementsInVolume( tMemory, 1 );

    // counter
    moris::uint tCount = 0;

    // loop over all objects
    for( moris::uint k=0; k < mNumberOfObjects; ++k)
    {
        for( moris::uint i=0; i<mObjects( k ).mLocalElementsInVolume.length(); ++i )
        {
            // save global element ID in vector
            tElementsInVolume(tCount) = mMeshData.get_global_element_id(
                    mObjects( k ).mLocalElementsInVolume( i ));
            ++tCount;

        }
    }

    // make result unique
    moris::Matrix< moris::DDUMat > tUnique;
    moris::unique(tElementsInVolume,tUnique);
    return tUnique;
}

// -----------------------------------------------------------------------------
/* moris::Matrix< moris::DDUMat >
ge::SDF_Gen::get_elements_in_volume_for_HMR_ref(){
    if ( mNumberOfObjects > 1)
    {
        // count max size of output matrix
        moris::uint tMemory = 0;
        for( moris::uint k=0; k < mNumberOfObjects; ++k)
        {
            tMemory += mObjects( k ).mLocalElementsInVolume.length();
        }

        // output matrix
        moris::Matrix< moris::DDUMat > tElementsInVolume(tMemory, 1);

        // counter
        moris::uint tCount = 0;

        // loop over all objects
        for( moris::uint k=0; k < mNumberOfObjects; ++k)
        {
            for( moris::uint i=0; i<mObjects( k ).mLocalElementsInVolume.length(); ++i)
            {
                // save global element ID in vector
                tElementsInVolume(tCount) = mMeshData.get_global_element_id(
                        mObjects( k ).mLocalElementsInVolume( i ));
                ++tCount;

            }
        }

        // make result unique
        return moris::unique(tElementsInVolume);
    }
    else
    {
        // return empty matrix
        moris::Matrix< moris::DDUMat > tEmptyMatrix;
        return tEmptyMatrix;
    }
} */

// -----------------------------------------------------------------------------

moris::Matrix< moris::DDUMat >
ge::SDF_Gen::get_node_ids(){
    return mMeshData.get_node_ids();
}

// -----------------------------------------------------------------------------

moris::Matrix< moris::DDRMat >
ge::SDF_Gen::get_node_signs()
{
    moris::uint tNumberOfNodes = mMeshData.get_number_of_nodes();

    // initialize output and assume all positive
    moris::Matrix< moris::DDRMat > aMat(tNumberOfNodes , 1, 1);

    // loop over all objects
    for ( uint k=0; k<mNumberOfObjects; ++k )
    {
        // loop over all nodes
        for( uint i=0; i<tNumberOfNodes; ++i )
        {
            // check if this node is inside in this object
            if ( mObjects( k ).mLocalNodeInsideFlags.test( i ) )
            {
                // flip sign of this node
                aMat( i ) = -1;
            }
        }
    }

    return aMat;
}

// -----------------------------------------------------------------------------

