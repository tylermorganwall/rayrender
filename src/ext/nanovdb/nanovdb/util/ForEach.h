// Copyright Contributors to the OpenVDB Project
// SPDX-License-Identifier: Apache-2.0
// Modified for rayrender: run CPU work on the calling thread so exceptions
// propagate to R and NanoVDB does not create an independent worker pool.

/*!
    \file nanovdb/util/ForEach.h

    \author Ken Museth

    \date August 24, 2020

    \brief Serial range traversal for rayrender's NanoVDB loading and validation
*/

#ifndef NANOVDB_UTIL_FOREACH_H_HAS_BEEN_INCLUDED
#define NANOVDB_UTIL_FOREACH_H_HAS_BEEN_INCLUDED

#include <nanovdb/util/Range.h>// for Range1D

namespace nanovdb {

namespace util {

/// @brief Visit a nonempty range on the calling thread.
///
/// @param range Range, CoordBBox, tbb::blocked_range, blocked_range2D, or blocked_range3D.
/// @param func functor with the signature [](const RangeT&){...},
///
/// @code
///     std::vector<int> array(100);
///     auto func = [&array](auto &r){for (auto i=r.begin(); i!=r.end(); ++i) array[i]=i;};
///     forEach(array, func);
/// @endcode
template <typename RangeT, typename FuncT>
inline void forEach(RangeT range, const FuncT &func)
{
    if (range.empty()) return;
    func(range);
}

/// @brief Simple wrapper for the function defined above
template <typename FuncT>
inline void forEach(size_t begin, size_t end, size_t grainSize, const FuncT& func)
{
    forEach(Range1D(begin, end, grainSize), func);
}

/// @brief Simple wrapper for the function defined above, which works with std::containers
template <template<typename...> class ContainerT, typename... T, typename FuncT>
inline void forEach(const ContainerT<T...> &c, const FuncT& func)
{
    forEach(Range1D(0, c.size(), 1), func);
}

/// @brief Simple wrapper for the function defined above, which works with std::containers
template <template<typename...> class ContainerT, typename... T, typename FuncT>
inline void forEach(const ContainerT<T...> &c, size_t grainSize, const FuncT& func)
{
    forEach(Range1D(0, c.size(), grainSize), func);
}

}// namespace util

/// @brief Simple wrapper for the function defined above
template <typename FuncT>
[[deprecated("Use nanovdb::util::forEach instead")]]
inline void forEach(size_t begin, size_t end, size_t grainSize, const FuncT& func)
{
    util::forEach(util::Range1D(begin, end, grainSize), func);
}

/// @brief Simple wrapper for the function defined above, which works with std::containers
template <template<typename...> class ContainerT, typename... T, typename FuncT>
[[deprecated("Use nanovdb::util::forEach instead")]]
inline void forEach(const ContainerT<T...> &c, const FuncT& func)
{
    util::forEach(util::Range1D(0, c.size(), 1), func);
}

/// @brief Simple wrapper for the function defined above, which works with std::containers
template <template<typename...> class ContainerT, typename... T, typename FuncT>
[[deprecated("Use nanovdb::util::forEach instead")]]
inline void forEach(const ContainerT<T...> &c, size_t grainSize, const FuncT& func)
{
    util::forEach(util::Range1D(0, c.size(), grainSize), func);
}

}// namespace nanovdb

#endif // NANOVDB_UTIL_FOREACH_H_HAS_BEEN_INCLUDED
