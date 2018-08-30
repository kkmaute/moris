/*
 * cl_Mesh_Bucket.hpp
 *
 *  Created on: Aug 17, 2017
 *      Author: ktdoble
 */

#ifndef SRC_MESH_CL_MESH_BUCKET_HPP_
#define SRC_MESH_CL_MESH_BUCKET_HPP_

using namespace xtk;
namespace mesh
{
template<typename Integer, typename Integer_Matrix>
class Bucket
{
public:
    Bucket()
{

}

    Bucket(Integer const & aMaxEntitiesInBucket,
           Integer const & aNumSubEntities,
           Integer const & aMaxNumberOfParts,
           Integer const & aStringLength):
               mBucketPhaseIndex(std::numeric_limits<Integer>::max())

    {
        mEntitiesConnectivity.reserve(aMaxEntitiesInBucket*aNumSubEntities);
        mPartNames.reserve(aMaxNumberOfParts*aStringLength);
        mEntityIds.reserve(aMaxEntitiesInBucket);
    }


    Bucket(Cell<std::string> const & aPartNames):
        mPartNames(aPartNames)
    {

    }

    void set_bucket_phase_index(Integer const & aBucketPhase)
    {
        mBucketPhaseIndex = aBucketPhase;
    }

    void add_entity(moris::Matrix<Integer, Integer_Matrix> & aEntityToAdd)
    {
        mEntitiesConnectivity.push_back(aEntityToAdd.copy());
    }

    void add_entities(Cell<moris::Matrix<Integer, Integer_Matrix>> const & aEntitiesToAdd)
    {
        mEntitiesConnectivity.append(aEntitiesToAdd);
    }

    void add_part_names(Cell<std::string> const & aPartNames)

    {
        mPartNames.append(aPartNames);
    }

    void add_part_name(std::string const & aPartName)
    {
        mPartNames.push_back(aPartName);
    }

    void add_entity_ids(Cell<Integer> const & aEntityIds)
    {
        mEntityIds.append(aEntityIds);
    }

    void add_entity_ids(moris::Matrix<Integer, Integer_Matrix> const & aEntityIds)
    {

        Integer tNumIds = aEntityIds.n_cols();

        for(Integer i = 0; i<tNumIds; i++)
        {
            mEntityIds.push_back(aEntityIds(0,i));
        }
    }


    void add_entity_ids(Integer const & aEntityId)
    {

//        if(mEntityIds.capacity()-mEntityIds.size() == 0 )
//        {
//            XTK_ERROR<<"Warning: Dynamic Allocation in add_entity_id";
//        }

        mEntityIds.push_back(aEntityId);
    }

    Integer get_num_entities_in_bucket() const
    {
        return mEntitiesConnectivity.size();
    }

    Integer get_num_parts() const
    {
        return mPartNames.size();
    }

    Cell<std::string> const & get_part_names() const
    {
        return mPartNames;
    }


    moris::Matrix<Integer, Integer_Matrix> const & get_entity(Integer aBucketLocalEntityIndex) const
    {
        return mEntitiesConnectivity(aBucketLocalEntityIndex);
    }


    Integer get_entity_id(Integer aBucketLocalEntityIndex) const
    {
        return mEntityIds(aBucketLocalEntityIndex);
    }

    bool bucket_has_entities() const
    {
        if(mEntitiesConnectivity.size()!=0)
        {
            return true;
        }

        else
            return false;
    }

    void print() const
    {
        std::cout<<" Number of Parts: " << mPartNames.size()<< std::endl;

        std::cout<<" Part Names: ";
        for(Integer iPart = 0; iPart<mPartNames.size(); iPart++)
        {
            std::cout<< mPartNames(iPart)<<"  ";
        }
        std::cout<<std::endl;
    }

private:
    Integer mBucketPhaseIndex;
    Cell<std::string> mPartNames;
    Cell<Integer> mEntityIds;
    Cell<moris::Matrix<Integer, Integer_Matrix>> mEntitiesConnectivity;



};
}


#endif /* SRC_MESH_CL_MESH_BUCKET_HPP_ */
