/*
 * cl_Mesh_Element_Bucket.hpp
 *
 *  Created on: Sep 1, 2017
 *      Author: doble
 */

#ifndef SRC_MESH_CL_MESH_ELEMENT_BUCKET_HPP_
#define SRC_MESH_CL_MESH_ELEMENT_BUCKET_HPP_
#include "cl_Mat.hpp" // LNA/src
#include "cl_Cell.hpp" // CON/src

namespace moris
{
     class Element_Bucket
     {
     public:
         Element_Bucket();

         Element_Bucket(uint const & aMaxEntitiesInBucket,
                        uint const & aNumberOfSecondaryEntitiesPerPrimaryEntity,
                        uint const & aMaxNumberOfParts,
                        uint const & aMaxStringLength);

         ~Element_Bucket();

         /**
          *  Add node to element connectivity
          *  @param[in] aEntitiesToAdd -  element to node connectivity
          **/
         void add_entities(Cell< Mat< uint > > const & aEntitiesToAdd);

         /*
          * Single element version of above
          */
         void add_entity(Mat< uint > const & aEntityToAdd );

         /*
          * Add part names associated with element bucket
          */
         void add_part_names(Cell< std::string > const & aPartNames);

         /*
          * Singular Form of the above
          */
         void add_part_name(std::string const & aPartName);

         /*
          * Add element ids
          */
         void add_entity_ids(Cell< uint > const & aEntityIds);

         /*
          * Singular form of above
          */
         void add_entity_id(uint const & aEntityId);

         /*
          * Returns whether or not the bucket has anything in it
          */
         bool has_entities() const;

         /*
          * Returns the number of entities in the bucket
          */
         uint get_num_entities_in_bucket() const;

         /*
          * Returns the number of parts associated with the bucket
          */
         uint get_num_parts() const;

         /*
          * Returns the part names as a const reference
          */
         Cell< std::string > const & get_part_names() const;

         /*
          * Access an entity in the bucket using the index that it is located at in the element bucket
          */
         Mat< uint > const & get_entity( uint const & aBucketEntityIndex);

         /*
          * Returns the entity id
          */
         uint const & get_entity_id(uint const & aBucketEntityIndex);
     private:
         Cell< std::string > mPartNames;
         Cell< uint > mEntityIds;
         Cell<Mat< uint > > mEntitiesConnectivity;
     };
}
#endif /* SRC_MESH_CL_MESH_ELEMENT_BUCKET_HPP_ */
