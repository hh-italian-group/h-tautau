/*! Definition of EventId class which represent unique CMS event identifier.
This file is part of https://github.com/hh-italian-group/h-tautau. */

#pragma once

#include <limits>
#include <iostream>

namespace analysis {

struct EventId{
    unsigned runId;
    unsigned lumiBlock;
    unsigned eventId;
    static const EventId& Undef_event() {
        static const EventId undef_event;
        return undef_event;
    }

    EventId() : runId(std::numeric_limits<UInt_t>::max()), lumiBlock(std::numeric_limits<UInt_t>::max()),
                eventId(std::numeric_limits<UInt_t>::max()){}

    EventId(unsigned _runId, unsigned _lumiBlock, unsigned _eventId) : runId(_runId), lumiBlock(_lumiBlock),
                eventId(_eventId){}

    bool operator == (const EventId& other) const
    {
        return !(*this != other);
    }

    bool operator != (const EventId& other) const
    {
        return runId != other.runId || lumiBlock != other.lumiBlock || eventId != other.eventId;
    }

    bool operator < (const EventId& other) const
    {
        if(runId < other.runId) return true;
        if(runId > other.runId) return false;
        if(lumiBlock < other.lumiBlock) return true;
        if(lumiBlock > other.lumiBlock) return false;
        return eventId < other.eventId;
    }
};

} // analysis

inline std::ostream& operator <<(std::ostream& s, const analysis::EventId& event)
{
    s << "run = " << event.runId << ", lumi = " << event.lumiBlock << ", evt = " << event.eventId;
    return s;
}
