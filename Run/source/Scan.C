/*! ROOT macro to scan analysis tree file for the specified event.
This file is part of https://github.com/hh-italian-group/h-tautau. */


void Scan(UInt_t event_id, const char* file_name)
{
    TFile f(file_name, "READ");
    std::ostringstream ss;
    ss << "EventId==" << event_id;
    const TSQLResult* result = events->Query("EventId", ss.str().c_str());
    if(result->GetRowCount() > 0)
        std::cerr << "found in " << file_name << std::endl;
}

