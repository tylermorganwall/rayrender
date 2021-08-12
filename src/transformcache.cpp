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