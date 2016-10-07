// simple_analyzer.cxx
#include <boost/format.hpp>
#include <thread>

#include "AnalysisTools/Core/include/RootExt.h"
#include "AnalysisTools/Run/include/program_main.h" // definition of wrappers for the program main and program arguments.
#include "h-tautau/Analysis/include/SyncTupleHTT.h"
#include "h-tautau/Analysis/include/EventTuple.h"

struct Arguments { // list of all program arguments
  REQ_ARG(std::string, input_file); // required argument "input_file"
  REQ_ARG(std::string, tree_name); // required argument "Name_Tree"
  REQ_ARG(std::string, output_file); // required argument "output_file"
};
namespace analysis{
class SyncTreeProducer { // simple analyzer definition
  public:
   using Event = ntuple::Event;
   using EventTuple = ntuple::EventTuple;
   using SyncEvent = htt_sync::SyncEvent;
   using SyncTuple = htt_sync::SyncTuple;
    SyncTreeProducer(const Arguments& _args) : args(_args)
  {
    // Analyzer initialization (e.g. open input/output files, parse configs...)
  }
    void Run()
    {
      // analyzer code
      std::cout << boost::format("Processing input file '%1%' into output file '%2%' with Tree Name = %3%.\n")
	% args.input_file() % args.output_file() % args.tree_name();

      auto originalFile = root_ext::OpenRootFile(args.input_file());
      auto outputFile = root_ext::CreateRootFile(args.output_file());
      std::shared_ptr<EventTuple> originalTuple(new EventTuple(args.tree_name(), originalFile.get(), true,
{ "lhe_n_partons", "lhe_HT" }));
     SyncTuple sync(args.tree_name(), outputFile.get(), false );
      const Long64_t n_entries = originalTuple->GetEntries();
      for(Long64_t current_entry = 0; current_entry < n_entries; ++current_entry) {
          originalTuple->GetEntry(current_entry);
	  const Event& event = originalTuple->data();
          if(event.p4_1.Pt()<20 /*|| event.p4_2.Pt() < 30*/) continue;
          sync().run = event.run;
          sync().lumi = event.lumi;
          sync().evt = event.evt;
          sync().npv = event.npv;
          sync().npu = event.npu;
          sync().pt_1 = event.p4_1.Pt(); 
          sync().phi_1 = event.p4_1.Phi(); 
          sync().eta_1 = event.p4_1.Eta(); 
          sync().m_1 = event.p4_1.mass(); 
          sync().q_1 = event.q_1;
          sync().d0_1 = event.d0_1; 
          sync().dZ_1 = event.dZ_1; 
          sync().mt_1 = event.mt_1;
          sync().pfmt_1 = event.pfmt_1;
          sync().puppimt_1 = event.puppimt_1;
          sync().iso_1   =  event.iso_1;
          sync().id_e_mva_nt_loose_1 = event.id_e_mva_nt_loose_1;
          sync().gen_match_1 = event.gen_match_1;
                    

 
          sync.Fill(); 
      }
      sync.Write();
    }
  private:
    Arguments args;
    void ReadThread(const std::string& treeName, const std::string& originalFileName){



    }

};
}
PROGRAM_MAIN(analysis::SyncTreeProducer, Arguments) // definition of the main program function
