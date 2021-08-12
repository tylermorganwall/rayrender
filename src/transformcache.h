#ifndef TRANSFORMCACHEH
#define TRANSFORMCACHEH

#include "transform.h"

class TransformCache {
public:
  TransformCache() : hashTable(512), hashTableOccupancy(0) {}
  
  // TransformCache Public Methods
  std::shared_ptr<Transform> Lookup(const Transform &t) {
    int offset = Hash(t) & (hashTable.size() - 1);
    int step = 1;
    while (true) {
      // Keep looking until we find the Transform or determine that
      // it's not present.
      if (!hashTable[offset] || *hashTable[offset] == t)
        break;
      // Advance using quadratic probing.
      offset = (offset + step * step) & (hashTable.size() - 1);
      ++step;
    }
    std::shared_ptr<Transform> tCached = hashTable[offset];
    if (tCached) {
    } else {
      tCached = std::make_shared<Transform>();
      *tCached = t;
      Insert(tCached);
    }
    return tCached;
  }
  
  void Clear() {
    // transformCacheBytes += arena.TotalAllocated() + hashTable.size() * sizeof(Transform *);
    hashTable.clear();
    hashTable.resize(512);
    hashTableOccupancy = 0;
    // arena.Reset();
  }
  
private:
  void Insert(std::shared_ptr<Transform> tNew);
  void Grow();
  
  static uint64_t Hash(const Transform &t) {
    const char *ptr = (const char *)(&t.GetMatrix());
    size_t size = sizeof(Matrix4x4);
    uint64_t hash = 14695981039346656037ull;
    while (size > 0) {
      hash ^= *ptr;
      hash *= 1099511628211ull;
      ++ptr;
      --size;
    }
    return hash;
  }
  
  // TransformCache Private Data
  std::vector<std::shared_ptr<Transform> > hashTable;
  int hashTableOccupancy;
  // MemoryArena arena;
};


#endif