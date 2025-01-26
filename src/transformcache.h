#ifndef TRANSFORMCACHEH
#define TRANSFORMCACHEH

#include "transform.h"

class TransformCache {
public:
  TransformCache() : hashTable(512), hashTableOccupancy(0) {}
  
  // TransformCache Public Methods
  Transform* Lookup(const Transform &t);
  void Clear();
private:
  void Insert(std::shared_ptr<Transform> tNew);
  void Grow();
  static uint64_t Hash(const Transform &t);
  // TransformCache Private Data
  std::vector<std::shared_ptr<Transform> > hashTable;
  unsigned int hashTableOccupancy;
  // MemoryArena arena;
};


#endif