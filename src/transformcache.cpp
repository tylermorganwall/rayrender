#include "transformcache.h"


void TransformCache::Insert(std::shared_ptr<Transform> tNew) {
  if (++hashTableOccupancy == hashTable.size() / 2) {
    Grow();
  }
  
  int baseOffset = Hash(*tNew) & (hashTable.size() - 1);
  for (int nProbes = 0;; ++nProbes) {
    // Quadratic probing.
    int offset = (baseOffset + nProbes/2 + nProbes*nProbes/2) & (hashTable.size() - 1);
    if (hashTable[offset] == nullptr) {
      hashTable[offset] = tNew;
      return;
    }
  }
}

void TransformCache::Grow() {
  std::vector<std::shared_ptr<Transform> > newTable(2 * hashTable.size());
  // LOG(INFO) << "Growing transform cache hash table to " << newTable.size();
  
  // Insert current elements into newTable.
  for (std::shared_ptr<Transform>&  tEntry : hashTable) {
    if (!tEntry) {
      continue;
    }
    
    int baseOffset = Hash(*tEntry.get()) & (hashTable.size() - 1);
    for (int nProbes = 0;; ++nProbes) {
      // Quadratic probing.
      int offset = (baseOffset + nProbes/2 + nProbes*nProbes/2) & (hashTable.size() - 1);
      if (newTable[offset] == nullptr) {
        newTable[offset] = tEntry;
        break;
      }
    }
  }
  std::swap(hashTable, newTable);
}

Transform* TransformCache::Lookup(const Transform &t) {
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
  Transform* tCached = hashTable[offset].get();
  std::shared_ptr<Transform> tCached_new;
  if (tCached) {
    //Do nothing
  } else {
    std::shared_ptr<Transform> tCached_new = std::make_shared<Transform>(t);
    Insert(tCached_new);
    tCached = tCached_new.get();
  }
  return tCached;
}

void TransformCache::Clear() {
  // transformCacheBytes += arena.TotalAllocated() + hashTable.size() * sizeof(Transform *);
  hashTable.clear();
  hashTable.resize(512);
  hashTableOccupancy = 0;
  // arena.Reset();
}

uint64_t TransformCache::Hash(const Transform &t) {
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