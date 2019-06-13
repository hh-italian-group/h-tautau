/*! Merge multiple root files into a single file splitting by channels.
This file is part of https://github.com/hh-italian-group/h-tautau. */

#include "AnalysisTools/Core/include/RootFilesMerger.h"
#include "AnalysisTools/Core/include/EventIdentifier.h"
#include "AnalysisTools/Run/include/program_main.h"
#include "h-tautau/Core/include/EventTuple.h"
#include "h-tautau/Core/include/SummaryTuple.h"
#include "h-tautau/Core/include/AnalysisTypes.h"

struct Arguments {
    run::Argument<std::string> output{"output", "output root file"};
    run::Argument<std::vector<std::string>> input_dirs{"input-dir", "input directory"};
    run::Argument<std::string> file_name_pattern{"file-name-pattern", "regex expression to match file names",
                                                 "^.*\\.root$"};
    run::Argument<std::string> exclude_list{"exclude-list", "comma separated list of files to exclude", ""};
    run::Argument<std::string> exclude_dir_list{"exclude-dir-list",
                                                "comma separated list of directories to exclude", ""};
    run::Argument<unsigned> n_threads{"n-threads", "number of threads", 1};
};

class TupleMerger : public analysis::RootFilesMerger {
public:
    TupleMerger(const Arguments& args) :
        RootFilesMerger(args.output(), args.input_dirs(), args.file_name_pattern(), args.exclude_list(),
                        args.exclude_dir_list(), args.n_threads(), ROOT::kLZ4, 5),
        output_summaryTuple("summary", output_file.get(), false),
        output_expressTuple("all_events", output_file.get(), false),
        n_total_duplicates(0)
    {
    }

    void Run()
    {
        Process(false, false);

        std::cout << "All file has been merged. Number of files = " << input_files.size() << std::endl;

        for(const auto& iter : output_eventTuple){
          iter.second->Write();
          std::cout << ". Number of output entries = " << iter.second->GetEntries() << std::endl;
        }
        output_summaryTuple.Write();
        output_expressTuple.Write();

        std::cout << ". Total number of duplicated entires = " << n_total_duplicates << "." << std::endl;
    }

  private:
      virtual void ProcessFile(const std::string& /*file_name*/, const std::shared_ptr<TFile>& file) override
      {
          ntuple::EventTuple eventTuple("events", file.get(), true);
          size_t n_duplicates = 0;
          for(const auto& event : eventTuple) {
              std::cout << "In eventTuple" << std::endl;
              const analysis::Channel channel = static_cast<analysis::Channel>(event.channelId);
              const analysis::EventIdentifier eventId(event.run,event.lumi,event.evt);
              if(processed_entries[channel].count(eventId)) {
                  ++n_duplicates;
                  continue;
              }
              processed_entries[channel].insert(eventId);
              std::cout << "Insert event" << std::endl;
              if(!output_eventTuple[channel])
                output_eventTuple[channel] = ntuple::CreateEventTuple(analysis::ToString(channel), file.get(), false,
                                                    ntuple::TreeState::Full);
              output_eventTuple[channel]->Fill();
              std::cout << "Filled tuple" << std::endl;
          }
          n_total_duplicates += n_duplicates;

          ntuple::SummaryTuple input_summaryTuple("summary", file.get(), true);
          for(const ntuple::ProdSummary& summary : input_summaryTuple) {
              output_summaryTuple() = summary;
              output_summaryTuple.Fill();
          }

          ntuple::ExpressTuple input_expressTuple("all_events", file.get(), true);
          for(const ntuple::ExpressEvent& express : input_expressTuple) {
              output_expressTuple() = express;
              output_expressTuple.Fill();
          }

          std::cout << "\tn_entries = " << eventTuple.GetEntries() << ", n_duplicates = " << n_duplicates << ".\n";
      }

  private:
      std::map<analysis::Channel, std::shared_ptr<ntuple::EventTuple>> output_eventTuple;
      ntuple::SummaryTuple output_summaryTuple;
      ntuple::ExpressTuple output_expressTuple;
      std::map<analysis::Channel, std::set<analysis::EventIdentifier>> processed_entries;
      size_t n_total_duplicates;
};

PROGRAM_MAIN(TupleMerger, Arguments)
