#ifndef TRANSFORMCACHEH
#define TRANSFORMCACHEH

#include "../math/transform.h"

class TransformCache {
public:
  TransformCache() : hashTable(512), hashTableOccupancy(0) {}
  
  // TransformCache Public Methods
  Transform* Lookup(const Transform &t);
  void Clear();
private:
  void Grow();
  static size_t FindSlot(const Transform &t,
                        const std::vector<std::shared_ptr<Transform>> &table);
  static uint64_t Hash(const Transform &t);
  // TransformCache Private Data
  std::vector<std::shared_ptr<Transform> > hashTable;
  unsigned int hashTableOccupancy;
  // MemoryArena arena;
};


#endif
