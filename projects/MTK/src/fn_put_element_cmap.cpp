

namespace moris
{
  namespace mtk
  {
    /**
    */
    void add_element_cmap_to_exodus(std::string aExodusFile)
     {
       // Open the Exodus file
       // TODO: SEE SYNTAX IN STK_CLASSIC

       // Get the sides on processor boundaries
       stk::mesh::EntityVector tSidesOnProcBoundaries;
       stk::mesh::get_selected_entities(mMtkMeshMetaData->globally_shared_part(),
                                        mMtkMeshBulkData->buckets(mMtkMeshMetaData->side_rank()),
                                        tSidesOnProcBoundaries);

       // Number of sides on mesh along processor boundaries.
       uint tNumSidesOnProcBoundaries = tSidesOnProcBoundaries.size();

       // Allocate vector for element ids attached to sides and side ordinals
       // attached to the (note: we limit the side to be attached to one element along the boundary)
       Matrix< IdMat > tElementIdsOnBoundaries(1,tNumSidesOnProcBoundaries);
       Matrix< IdMat > tSideOrdinalsOnBoundaries(1,tNumSidesOnProcBoundaries);
       Matrix< IdMat > tSideSharedProc(1,tNumSidesOnProcBoundaries);

       // Allocate vector for elements connected to a given faces
       Matrix< IdMat > tElementToFace(1,2);
       Matrix< IdMat > tSharedProcessorsOfFace(1,2);
       uint            tProcRank = par_rank();

       // Populate the above vectors
       for(uint iS = 0; iS<tNumSidesOnProcBoundaries; iS++)
       {
          // Get the side id
          moris_id tSideId = mMtkMeshBulkData->identifier(tSidesOnProcBoundaries[i]);

          // Get the elements attached to this faces
          tElementToFace = this->get_entities_connected_to_entities_glob_id(EntityRank::FACE,EntityRank::ELEMENT,tSideId);

          // Figure out which other processors share this sides
          tSharedProcessorsOfFace = this->get_procs_whom_share_face();
          MORIS_ASSERT(tSharedProcessorsOfFace.numel() == 2,"For the current implementation, a face can only be shared by 2 processors");
          for( uint  iP = 0; iP<tSharedProcessorsOfFace.numel(); iP++)
          {
              if(tSharedProcessorsOfFace(iP) != tProcRank)
              {
                tSideSharedProc(iS) = tSharedProcessorsOfFace(iP);
                break;
              }
          } // shared proc loop

          // iterate through elements connected to faces
          for(uint iE = 0; iE<tElementToFace.numel(); iE++)
          {
              // get the element id from the vector
              moris_id tElementId = tElementToFace(iE);

              // check whether the element is in the aura
              if(!this->is_aura_element(tElementId))
              {
                // figure out the side ordinal of this element where tSideId is
                moris_id tSideOrdinal = this->get_facet_ordinal_from_element_and_facet_id_glob_ids(tElementId,tSideId);

                // Add to vectors
                tElementIdsOnBoundaries(iS) = tElementId;
                tSideOrdinalsOnBoundaries(iS) = tSiderOrdinal;
              }
          } // Element loop
       } // side loop

       // now that data has been assembled, call ex_put_element_cmap()....
     }
  }
}
