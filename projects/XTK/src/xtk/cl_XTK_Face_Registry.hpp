/*
 * cl_XTK_Face_Registry.hpp
 *
 *  Created on: Feb 5, 2018
 *      Author: ktdoble
 */

#ifndef XTK_SRC_XTK_CL_XTK_FACE_REGISTRY_HPP_
#define XTK_SRC_XTK_CL_XTK_FACE_REGISTRY_HPP_


// XTKL: Linear Algebra Includes
#include "../linalg/cl_XTK_Matrix_Base_Utilities.hpp"
#include "linalg/cl_XTK_Matrix_Base.hpp"
#include "tools/fn_bubble_sort.hpp"

namespace xtk
{

// Forward declaration of child mesh class
template<typename Real, typename Integer, typename Real_Matrix, typename Integer_Matrix> class Child_Mesh;

/*
 * This class tracks the creation of new faces and returns a element to face connectivity
 */
template<typename Real, typename Integer, typename Real_Matrix, typename Integer_Matrix>
class Face_Registry
{
public:
    /*
     * Constructor for face registry which needs to include ancestry
     */
    Face_Registry(Integer aMaxNumNewElements,
                  Integer aFacesPerElement,
                  moris::Matrix< moris::IndexMat > const &  aFacetoNodeConnectivity,
                  moris::Matrix< moris::IndexMat > const &  aFaceToElementConnectivity,
                  moris::Matrix< moris::IndexMat > &  aFaceParentIndices,
                  moris::Matrix< Integer_Matrix > &  aFaceParentRanks):
                      mDummyValue(std::numeric_limits<moris::moris_index>::max()),
                      mModificationOpen(true)
    {
        Integer tNumRows  = aFacetoNodeConnectivity.n_rows();
        Integer tNumCols  = aFacetoNodeConnectivity.n_cols();
        Integer tMaxSize  = tNumRows+(aMaxNumNewElements*aFacesPerElement);

        // Mark the first place where a face is appended and initialize a counter for the total number of registered faces
        mFirstAppendedFace   = tNumRows;
        mNumRegisteredFaces  = tNumRows;

        // Allocate Space then make a copy of the face to node connectivity
        mFaceToNode = moris::Matrix< moris::IndexMat >(tMaxSize,tNumCols,0);
        conservative_copy(aFacetoNodeConnectivity,mFaceToNode);

        // Allocate Space then make a copy of the face ancestry connectivity
        mFaceAncestryIndices = moris::Matrix< moris::IndexMat >(1,tMaxSize,0);
        conservative_copy(aFaceParentIndices, mFaceAncestryIndices);
        mFaceAncestryRanks = moris::Matrix< Integer_Matrix >(1,tMaxSize,0);
        conservative_copy(aFaceParentRanks, mFaceAncestryRanks);
        mReplacedInheritanceMarker = moris::Matrix< Integer_Matrix >(1,tMaxSize,0);


        // Initialize a face replacement marker
        mReplacedFaceMarker  = moris::Matrix< Integer_Matrix >(1,tNumRows,0);

        //Sort the face to node connectivity
        //TODO: REMOVE with new method of determining face
        row_bubble_sort(mFaceToNode);

        // Face to element connectivity, has a maximum of number of elements
        // each face can be connected to 2 elements
        Integer tMaxFaceToElement = aFaceToElementConnectivity.n_rows() + (aMaxNumNewElements*aFacesPerElement);
        mFaceToElement = moris::Matrix< moris::IndexMat >(tMaxFaceToElement,2);
        mFaceToElement.fill(mDummyValue);
        conservative_copy(aFaceToElementConnectivity, mFaceToElement);

    }


    // Constructor for a face registry used to generate face indices,no ancestry information, also not used for modifying the face connectivity
    Face_Registry(Integer aFacesPerElement,
                  moris::Matrix< moris::IndexMat > const &  aFacetoNodeConnectivity,
                  moris::Matrix< moris::IndexMat > const &  aFaceToElementConnectivity):
                      mDummyValue(std::numeric_limits<moris::moris_index>::max()),
                      mModificationOpen(false)
    {
        Integer tNumRows  = aFacetoNodeConnectivity.n_rows();
        Integer tNumCols  = aFacetoNodeConnectivity.n_cols();
        Integer tMaxSize  = tNumRows;

        // Mark the first place where a face is appended and initialize a counter for the total number of registered faces
        mFirstAppendedFace   = tNumRows;
        mNumRegisteredFaces  = tNumRows;

        // Allocate Space then make a copy of the face to node connectivity
        mFaceToNode = moris::Matrix< moris::IndexMat >(tMaxSize,tNumCols,0);
        conservative_copy(aFacetoNodeConnectivity,mFaceToNode);

        // Initialize a face replacement marker
        mReplacedFaceMarker  = moris::Matrix< Integer_Matrix >(1,tNumRows,0);

        //Sort the face to node connectivity
        row_bubble_sort(mFaceToNode);

        // Face to element connectivity, has a maximum of number of elements
        // each face can be connected to 2 elements
        Integer tMaxFaceToElement = aFaceToElementConnectivity.n_rows();
        mFaceToElement = moris::Matrix< moris::IndexMat >(tMaxFaceToElement,2);
        mFaceToElement.fill(mDummyValue);
        conservative_copy(aFaceToElementConnectivity, mFaceToElement);

    }

    // Dummy value
    moris::moris_index const mDummyValue;

    // Returns the face indices of a provide list of face indices
    moris::Matrix< moris::IndexMat >
    get_face_indices(moris::Matrix< moris::IndexMat > & aFaceToNodeConnectivity,
                     bool aSort = false)
    {
        // Initialize information
        Integer tNumFaces = aFaceToNodeConnectivity.n_rows();
        moris::Matrix< moris::IndexMat > tFaceIndices(1,tNumFaces);

        // Sort the rows in ascending order (Note this changes the input matrix but does not change the row order)
        if(aSort){row_bubble_sort(aFaceToNodeConnectivity);}

        // Loop over all faces and determine if its equivalent to any in mFaceToNodes
        for(Integer i = 0; i<tNumFaces; i++)
        {
            (tFaceIndices)(0,i) = get_face_index(i,aFaceToNodeConnectivity);
        }

        return tFaceIndices;
    }

    // single face version of the above
    Integer
    get_face_indices(moris::Matrix< moris::IndexMat > & aFaceToNodeConnectivity,
                     Integer aRowIndex,
                     bool aSort = false)
    {
        // Sort the rows in ascending order (Note this changes the input matrix but does not change the row order)
        if(aSort){row_bubble_sort(aFaceToNodeConnectivity);}

        return get_face_index(aRowIndex,aFaceToNodeConnectivity);
    }

    /**
     * Replaces a face in the face to node list at location aFaceIndexToReplace,
     * with the face located in row aFaceIndexInFaceNodeConn of aFaceToNodeConnectivity.
     */
    void replace_face( Integer const & aFaceIndexToReplace,
                       Integer const & aFaceIndexInFaceNodeConn,
                       moris::Matrix< moris::IndexMat > const & aFaceToNodeConnectivity)
    {
        XTK_ASSERT(aFaceIndexToReplace<mReplacedFaceMarker.n_cols(),"Replacing a newly added face is not allowed");
        XTK_ASSERT(mReplacedFaceMarker(0,aFaceIndexToReplace) == 0,"Face has already been replaced, choose another face to replace");
        replace_row(aFaceIndexInFaceNodeConn,aFaceToNodeConnectivity,aFaceIndexToReplace,mFaceToNode);
        fill_row(mDummyValue,aFaceIndexToReplace,mFaceToElement);
        mReplacedFaceMarker(0,aFaceIndexToReplace) = 1;

    }

    /*
     * Appends the face index to the first available face index
     */
    Integer append_face(Integer const & aFaceIndexInFaceNodeConn,
                        moris::Matrix< moris::IndexMat >  const & aFaceToNodeConnectivity)
    {
        Integer tNewFaceIndex = mNumRegisteredFaces;
        replace_row(aFaceIndexInFaceNodeConn, aFaceToNodeConnectivity, tNewFaceIndex, mFaceToNode);
        mNumRegisteredFaces++;
        return tNewFaceIndex;
    }

    /*
     * Modifies the face to element connectivity by telling the face registry this element is attached to these faces
     */
    void
    set_face_to_element(moris::Matrix< moris::IndexMat >  const & aElementIndices,
                        moris::Matrix< moris::IndexMat >  const & aElementToFaceIndices)
    {

        Integer tNumElems  = aElementIndices.n_cols();
        Integer tNumFaces  = aElementToFaceIndices.n_cols();
        Integer tNumCols   = mFaceToElement.n_cols();
        Integer tFaceIndex = 0;

        // Loop over elements
        for(Integer iE = 0; iE<tNumElems; iE++)
        {
            // Loop over faces and add face to element connectivity
            for(Integer iF = 0; iF<tNumFaces; iF++)
            {
                // Face Index iF of element iE (for clarity later)
                tFaceIndex = aElementToFaceIndices(iE,iF);


                // Loop over the two entries in the
                for(Integer j = 0; j<tNumCols; j++)
                {
                    // If the entry at face index does not have a value add one
                    if((mFaceToElement)(tFaceIndex,j) == mDummyValue)
                    {
                        (mFaceToElement)(tFaceIndex,j) = aElementIndices(0,iE);
                        break;
                    }
                }
            }
        }

    }

    /*
     * Reset face to element connectivity for a given element index.
     * If the elements that a face are connected to changes then this resets the entry corresponding to the element index
     * This does not change the inheritance but does indicate that the face has been replaced
     */
    void
    reset_face_to_element(moris::moris_index aFaceIndexToReset,
                          moris::moris_index aElementToReplace)
    {


        for(Integer i = 0; i<mFaceToElement.n_cols(); i++)
        {
            if( (mFaceToElement)(aFaceIndexToReset,i) == aElementToReplace)
            {
                (mFaceToElement)(aFaceIndexToReset,i) = mDummyValue;
            }
        }

        mReplacedFaceMarker(0,aFaceIndexToReset) = 1;

    }


    void set_face_ancestry(moris::moris_index aFaceIndex,
                           moris::moris_index aParentIndex,
                           Integer            aParentRank,
                           bool               aCheckDuplicates = false)
    {
        // check for a duplicate face change
        if(aCheckDuplicates)
        {
            // If the face inheritance has been replaced
            // make sure the rank and index are consistent
             if(has_inheritance_been_replaced(aFaceIndex))
             {

                 if(mFaceAncestryRanks(0,aFaceIndex) != aParentRank || mFaceAncestryIndices(0,aFaceIndex) != aParentIndex)
                 {
                     std::cout<<"Existing Parent Rank: " <<mFaceAncestryRanks(0,aFaceIndex)<< " New Parent Rank: "<< aParentRank<<std::endl;
                     std::cout<<"Existing Parent Index: " <<mFaceAncestryIndices(0,aFaceIndex)<< " New Parent Index: "<< aParentIndex<<std::endl;
                 }
                 XTK_ASSERT(mFaceAncestryRanks(0,aFaceIndex) == aParentRank,"Inheritance Rank mismatch on a given face");
                 XTK_ASSERT(mFaceAncestryIndices(0,aFaceIndex) == aParentIndex,"Inheritance Index mismatch on a given face");
             }
        }

        mFaceAncestryRanks(0,aFaceIndex)   = aParentRank;
        mFaceAncestryIndices(0,aFaceIndex) = aParentIndex;

        mReplacedInheritanceMarker(0,aFaceIndex) = 1;
    }


    /*
     * Ask if a face has been replaced, returns true if the face has already been replaced.
     */
    bool has_face_been_replaced(Integer aFaceIndexToReplace)
    {
        // if the parent has been replaced
        if (mReplacedFaceMarker(0,aFaceIndexToReplace) == 1 )
        {
            return true;
        }

        return false;
    }

    /*
     * Checks whether a face's inheritance has been replaced already
     */
    bool has_inheritance_been_replaced(Integer aFaceIndexToReplace)
    {
        if (mReplacedInheritanceMarker(0,aFaceIndexToReplace) == 1 )
         {
            return true;
         }

        return false;
    }

    /*
     * Returns whether the provided face index has been created during this face registry cycle or existed at construction time.
     * Returns true if the face was appended
     */
    bool is_appended_face(Integer aFaceIndex)
    {
        bool IsAppended = false;

        if(aFaceIndex>=mFirstAppendedFace)
        {
            IsAppended = true;
        }

        return IsAppended;
    }

    /*
     *  Returns the final face to node connectivity
     *  Note: this should be called after end_modification cycle
     */
    moris::Matrix< moris::IndexMat > & get_face_to_node()
        {
        return mFaceToNode;
        }

    /*
     * Returns the face ancestry parent indices
     */
    moris::Matrix< moris::IndexMat >  & get_face_inheritance_indices()
        {
        return mFaceAncestryIndices;
        }

    /*
     * Return the face ancestry parent ranks
     */
    moris::Matrix< Integer_Matrix >  & get_face_inheritance_ranks()
    {
        return mFaceAncestryRanks;
    }

    moris::Matrix< moris::IndexMat > &  get_face_to_element()
    {
        return mFaceToElement;
    }


    /*
     * Ends the modification cycle,
     * After this call, no modifications should be made via calling append_face,  or replace face
     * because all excess size has been removed from member variables
     */
    void end_modification()
    {
        Integer tNumCols = mFaceToNode.n_cols();
        mFaceToNode.resize(mNumRegisteredFaces,tNumCols);
        mFaceAncestryIndices.resize(1,mNumRegisteredFaces);
        mFaceAncestryRanks.resize(1,mNumRegisteredFaces);
        mFaceToElement.resize(mNumRegisteredFaces,2);


        // Close the modification cycle (this should only be called here)
        mModificationOpen = false;
    }



private:
    bool mModificationOpen;
    Integer mFirstAppendedFace;
    Integer mNumRegisteredFaces;
    moris::Matrix< Integer_Matrix > mReplacedFaceMarker;
    moris::Matrix< Integer_Matrix > mFaceToElementCounter;
    moris::Matrix< Integer_Matrix > mReplacedInheritanceMarker;
    moris::Matrix< moris::IndexMat > mFaceToNode;
    moris::Matrix< moris::IndexMat > mFaceAncestryIndices;
    moris::Matrix< Integer_Matrix > mFaceAncestryRanks;
    moris::Matrix< moris::IndexMat > mFaceToElement;

    /*
     * Returns face index if it exists otherwise it returns the numerical limit
     */
    moris::moris_index
    get_face_index( Integer const & aRowIndex,
                    moris::Matrix< moris::IndexMat > & aFaceToNodeConnectivity)
    {
        return face_exists(aRowIndex,aFaceToNodeConnectivity);
    }

    /*
     * Checks whether a face exists and if it does returns its index
     */
    moris::moris_index
    face_exists(Integer const & aRowIndex,
                moris::Matrix< moris::IndexMat > & aFaceToNodeConnectivity)
    {
        moris::moris_index  tFaceIndex = mDummyValue;

        bool tExists = false;
        for(size_t i = 0; i<mFaceToNode.n_rows(); i++)
        {

            tExists =  row_equal(aRowIndex,aFaceToNodeConnectivity,i,mFaceToNode);

            if(tExists)
            {
                tFaceIndex = i;
                break;
            }
        }

        return tFaceIndex;

    }
};

}
#endif /* XTK_SRC_XTK_CL_XTK_FACE_REGISTRY_HPP_ */
