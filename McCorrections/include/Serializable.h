#pragma once

// The archives must be listed before any boost/serialization header.
// Otherwise, in some cases the export macros trigger compilation errors.
// #include "CondFormats/Serialization/interface/Archive.h"

#include <boost/serialization/access.hpp>

#include <boost/serialization/string.hpp>
#include <boost/serialization/vector.hpp>
#include <boost/serialization/list.hpp>
#include <boost/serialization/set.hpp>
#include <boost/serialization/map.hpp>
#include <boost/serialization/bitset.hpp>
#include <boost/serialization/shared_ptr.hpp>
#include <boost/serialization/unordered_map.hpp>

#include "AnalysisTools/Core/include/exception.h"

// We cannot include Equal.h here since it is C++11
namespace cond_standalone {
namespace serialization {
    template <typename CondSerializationT, typename Enabled = void>
    struct access;
}
}

// Marks a class/struct as serializable Conditions.
// It must be used in the end of the class/struct, to avoid
// changing the default access specifier.
// Note: the serialization code generator script depends on
//       the implementation of the macro.
#define COND_SERIALIZABLE_STANDALONE \
    private: \
        template <class Archive> void serialize(Archive & ar, const unsigned int version); \
        template <typename CondSerializationT, typename Enabled> friend struct cond_standalone::serialization::access; \
        friend class boost::serialization::access;

