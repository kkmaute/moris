/*
 * cl_Mesh_Element_Bucket.cpp
 *
 *  Created on: Sep 1, 2017
 *      Author: doble
 */

#include "cl_Mesh_Element_Bucket.hpp" // STK/src/Hierarchical
#include "assert.hpp"
namespace moris
{
    Element_Bucket::Element_Bucket()
    {

    }

    Element_Bucket::Element_Bucket(uint const & aMaxEntitiesInBucket,
            uint const & aNumberOfSecondaryEntitiesPerPrimaryEntity,
            uint const & aMaxNumberOfParts,
            uint const & aMaxStringLength)
    {
        mEntitiesConnectivity.reserve(aMaxEntitiesInBucket*aNumberOfSecondaryEntitiesPerPrimaryEntity);
        mPartNames.reserve(aMaxNumberOfParts*aMaxStringLength);
        mEntityIds.reserve(aMaxEntitiesInBucket);

    }

    Element_Bucket::~Element_Bucket()
    {

    }

    void Element_Bucket::add_entities(Cell< Mat< uint > > const & aEntitiesToAdd)
    {
        mEntitiesConnectivity.append(aEntitiesToAdd);
    }

    void Element_Bucket::add_entity(Mat< uint > const & aEntityToAdd )
    {
        mEntitiesConnectivity.push_back(aEntityToAdd);
    }
    void Element_Bucket::add_part_names(Cell< std::string > const & aPartNames)
    {
        mPartNames.append(aPartNames);
    }

    /*
     * Singular Form of the above
     */
    void Element_Bucket::add_part_name(std::string const & aPartName)
    {
        mPartNames.push_back(aPartName);
    }
    void Element_Bucket::add_entity_ids(Cell< uint > const & aEntityIds)
    {
        mEntityIds.append(aEntityIds);
    }

    void Element_Bucket::add_entity_id(uint const & aEntityId)
    {
        mEntityIds.push_back(aEntityId);
    }

    bool Element_Bucket::has_entities() const
    {
        if(this->get_num_entities_in_bucket() == 0)
        {
            return false;
        }
        else
        {
            return true;
        }

    }
    uint Element_Bucket::get_num_entities_in_bucket() const
    {
        return mEntitiesConnectivity.size();
    }
    uint Element_Bucket::get_num_parts() const
    {
        return mPartNames.size();
    }
    Cell< std::string > const & Element_Bucket::get_part_names() const
    {
        return mPartNames;
    }
    Mat< uint > const & Element_Bucket::get_entity( uint const & aBucketEntityIndex)
    {
        MORIS_ASSERT( get_num_entities_in_bucket() > aBucketEntityIndex, "Request index is out of bounds in element bucket");
        return  mEntitiesConnectivity(aBucketEntityIndex);
    }

    uint const & Element_Bucket::get_entity_id(uint const & aBucketEntityIndex)
    {
        MORIS_ASSERT( get_num_entities_in_bucket() > aBucketEntityIndex, "Request index is out of bounds in element bucket");
        return mEntityIds(aBucketEntityIndex);
    }
}

