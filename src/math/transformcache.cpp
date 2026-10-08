#include "../math/transformcache.h"

size_t TransformCache::FindSlot(
    const Transform &t, const std::vector<std::shared_ptr<Transform>> &table) {
  const size_t mask = table.size() - 1;
  size_t offset = Hash(t) & mask;
  // Triangular probing visits every slot in a power-of-two table. Lookup,
  // insertion and growth must use exactly the same sequence and table mask.
  for (size_t step = 1; table[offset] && *table[offset] != t; ++step)
    offset = (offset + step) & mask;
  return offset;
}

void TransformCache::Grow() {
  std::vector<std::shared_ptr<Transform>> enlarged(2 * hashTable.size());
  for (const auto &entry : hashTable)
    if (entry) enlarged[FindSlot(*entry, enlarged)] = entry;
  hashTable.swap(enlarged);
}

Transform *TransformCache::Lookup(const Transform &t) {
  size_t offset = FindSlot(t, hashTable);
  if (!hashTable[offset]) {
    if (hashTableOccupancy + 1 >= hashTable.size() / 2) {
      Grow();
      offset = FindSlot(t, hashTable);
    }
    hashTable[offset] = std::make_shared<Transform>(t);
    ++hashTableOccupancy;
  }
  return hashTable[offset].get();
}

void TransformCache::Clear() {
  hashTable.clear();
  hashTable.resize(512);
  hashTableOccupancy = 0;
}

uint64_t TransformCache::Hash(const Transform &t) {
  uint64_t hash = 14695981039346656037ull;
  for (int row = 0; row < 4; ++row) {
    for (int column = 0; column < 4; ++column) {
      // Matrix equality treats signed zeros as equal, so their hashes must
      // agree as well. Other finite values retain their bit representation.
      const Float original = t.GetMatrix().m[row][column];
      const Float value = original == 0 ? Float(0) : original;
      const auto *bytes = reinterpret_cast<const unsigned char *>(&value);
      for (size_t i = 0; i < sizeof(Float); ++i) {
        hash ^= bytes[i];
        hash *= 1099511628211ull;
      }
    }
  }
  return hash;
}
