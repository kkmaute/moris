/*
 * cl_Entity_Tracker.cpp
 *
 *  Created on: Apr 3, 2017
 *      Author: doble
 */


#include"cl_Entity_Tracker.hpp"

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// Entity Tracker Object
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

moris::xtk::EntityTracker::EntityTracker(enum EntityRank aEntityRanktoTrack,
                                         enum EntityRank aChildEntityRank,
                                         moris::uint     aNumEntitiestoTrack,
                                         moris::uint     aNumChildrenAllowed):
                 mUseMarker(aNumEntitiestoTrack,1,0),
                 mEntityTrackerInfo(aNumEntitiestoTrack,aNumChildrenAllowed*3,-1),
                 mRequestIndexTracker(aNumEntitiestoTrack,1,UINT_MAX),
                 mReqCounter(0)
{
    mChildrenAllowed = aNumChildrenAllowed;
}

//-----------------------------------------------------------------------------

moris::xtk::EntityTracker::~EntityTracker()
{

}

//-----------------------------------------------------------------------------

void
moris::xtk::EntityTracker::mark_entity_as_used(moris::uint aEntityIndex)
{
    mUseMarker(aEntityIndex) = true;
}

//-----------------------------------------------------------------------------

void
moris::xtk::EntityTracker::set_child_entity_lcl_index(moris::uint aParentEntityIndex,
                                                      moris::uint aSecondaryIndex,
                                                      moris::uint aChildEntityIndex)
{

    moris::uint tLoc = get_child_location(aParentEntityIndex,aSecondaryIndex);

//    MORIS_LOG_INFO<< "Setting lcl ind";
//    MORIS_LOG_INFO<< "index : "<<tLoc + 2*mChildrenAllowed;
    mEntityTrackerInfo(aParentEntityIndex,tLoc + 2*mChildrenAllowed) = aChildEntityIndex;
}

//-----------------------------------------------------------------------------

void
moris::xtk::EntityTracker::set_child_entity_glb_id(moris::uint aParentEntityIndex,
                                                   moris::uint aSecondaryIndex,
                                                   moris::uint aChildEntityId)
{

    moris::uint tLoc = get_child_location(aParentEntityIndex,aSecondaryIndex);

//    MORIS_LOG_INFO<< "Setting glb Id";
    mEntityTrackerInfo(aParentEntityIndex,tLoc + mChildrenAllowed) = aChildEntityId;
}

//-----------------------------------------------------------------------------

bool
moris::xtk::EntityTracker::is_parent_entity_used(moris::uint aEntityIndex)
{
    return mUseMarker(aEntityIndex);
}

moris::Cell<moris::uint*>
moris::xtk::EntityTracker::is_parent_entity_used(moris::uint aEntityIndex,
                                                 moris::uint aSecondaryIndex)
{
    MORIS_ASSERT(mChildrenAllowed!=1,"If only one child is allowed then secondary index is not needed. Use other is_parent_entity_used(aEntityIndex)");

    // Intialize as null pointers
    moris::Cell<moris::uint*> tAnswer(3);
    tAnswer(0) = NULL;

    bool used = false;
    for(moris::uint i = 0; i<mUseMarker(aEntityIndex);i++)
    {
        if(mEntityTrackerInfo(aEntityIndex,i)==aSecondaryIndex)
        {
            tAnswer(0) = & mEntityTrackerInfo(aEntityIndex,i +   mChildrenAllowed);
            tAnswer(1) = & mEntityTrackerInfo(aEntityIndex,i +   mChildrenAllowed); // Id pointer
            tAnswer(2) = & mEntityTrackerInfo(aEntityIndex,i + 2*mChildrenAllowed); // Index pointer

            used =true;
        }
    }

    // If used was never marked true the entity has not been used
    if(used == false)
    {
        mEntityTrackerInfo(aEntityIndex,mUseMarker(aEntityIndex)) = aSecondaryIndex;
        tAnswer(1) = & mEntityTrackerInfo(aEntityIndex,mUseMarker(aEntityIndex) +   mChildrenAllowed);
        tAnswer(2) = & mEntityTrackerInfo(aEntityIndex,mUseMarker(aEntityIndex) + 2*mChildrenAllowed);
        mUseMarker(aEntityIndex)++;
    }

    return tAnswer;
}

//-----------------------------------------------------------------------------

moris::uint*
moris::xtk::EntityTracker::get_index_pointer(moris::uint aParentEntityIndex)
{
    return & mEntityTrackerInfo(aParentEntityIndex,2);
}

//-----------------------------------------------------------------------------

moris::uint*
moris::xtk::EntityTracker::get_id_pointer(moris::uint aParentEntityIndex)
{
    return & mEntityTrackerInfo(aParentEntityIndex,1);
}

//-----------------------------------------------------------------------------

moris::uint
moris::xtk::EntityTracker::get_request_index_from_entity_index(moris::uint aEntityIndex)
{
    if(mRequestIndexTracker(aEntityIndex,0) == UINT_MAX)
    {
        mRequestIndexTracker(aEntityIndex,0) = mReqCounter;
        mUseMarker(aEntityIndex) = true;
        mReqCounter++;
    }

    return mRequestIndexTracker(aEntityIndex);
}

//-----------------------------------------------------------------------------

moris::uint
moris::xtk::EntityTracker::get_num_children_allowed()
{
    return mChildrenAllowed;
}


//-----------------------------------------------------------------------------
void
moris::xtk::EntityTracker::print()
{
    for(moris::uint i = 0; i<mEntityTrackerInfo.n_rows();i++)
    {
        //MORIS_LOG_INFO<<i<<"| "<< mEntityTrackerInfo.row(i);
    }
}

moris::uint
moris::xtk::EntityTracker::get_child_location(moris::uint  aParentIndex,
                                              moris::uint  aSecondaryIndex)
{
    moris::uint tLoc = UINT_MAX;
    MORIS_ASSERT(aParentIndex<mEntityTrackerInfo.n_rows(),"Attempted to access entity outside of entity tracker bounds.");
    for(moris::uint i = 0; i<mChildrenAllowed; i++)
    {
        if(mEntityTrackerInfo(aParentIndex,i) == aSecondaryIndex)
        {
            tLoc = i;
        }
    }
    if(tLoc==UINT_MAX)
    {
        MORIS_LOG_INFO<<"Parent Index: "<< aParentIndex;
        MORIS_LOG_INFO<<"Secondary Index: "<< aSecondaryIndex;
        //MORIS_LOG_INFO<<"Entity Tracker: "<< mEntityTrackerInfo.row(aParentIndex);
    }
    MORIS_ASSERT(tLoc!=UINT_MAX,"get_child_location did not find a child on the provided parent entity with the given secondary index");
    return tLoc;
}

