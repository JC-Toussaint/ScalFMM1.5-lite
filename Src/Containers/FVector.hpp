// See LICENCE file at project root
#ifndef FVECTOR_HPP
#define FVECTOR_HPP

#include <vector>

// Compatibility with the code written for ScalFMM 1.5 (e.g. getIndexes() of the
// indexed containers, which now returns a std::vector): use std::vector directly.
template <class ObjectType>
using FVector = std::vector<ObjectType>;

#endif // FVECTOR_HPP
