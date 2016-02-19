/*! Print stack with specified name superimposing several files.
This file is part of https://github.com/hh-italian-group/h-tautau. */

#include <TTree.h>
#include <iostream>
#include <fstream>
#include <sstream>

#include "Analysis/include/FlatAnalyzerDataCollection.h"
#include "AnalysisTools/Core/include/ConfigReader.h"
#include "AnalysisTools/Print/include/RootPrintToPdf.h"

struct LoopOptions {
    bool eventCategory, eventSubCategory, eventRegion, eventEnergyScale;
    LoopOptions() : eventCategory(true), eventSubCategory(true), eventRegion(true), eventEnergyScale(true) {}
};

static const std::map<analysis::EventCategory, std::string> eventCategoryMap =
          { { analysis::EventCategory::Inclusive, "Inclusive" },
//            { analysis::EventCategory::TwoJets_Inclusive, "2jet-Inclusive" },
            { analysis::EventCategory::TwoJets_Inclusive, "inclusive" },
            { analysis::EventCategory::TwoJets_ZeroBtag, "2jet-0tag" },
            { analysis::EventCategory::TwoJets_OneBtag, "2jet-1tag"},
            { analysis::EventCategory::TwoJets_TwoBtag, "2jet-2tag" },
            { analysis::EventCategory::TwoJets_ZeroLooseBtag, "2jet-0Loosebtag" },
            { analysis::EventCategory::TwoJets_OneLooseBtag, "2jet-1Loosebtag" },
            { analysis::EventCategory::TwoJets_TwoLooseBtag, "2jet-2Loosebtag" },
            { analysis::EventCategory::TwoJets_AtLeastOneBtag, "2jets_at_least_1btag" },
          { analysis::EventCategory::TwoJets_AtLeastOneLooseBtag, "2jets_at_least_1Loosebtag" }};

class Print_Stack {
public:
    Print_Stack(const std::string& source_cfg, const std::string& anaDataFileName, const std::string& output_path,
                const std::string& channel_name, const std::string& id_selection, const std::string& _hist_name,
                const std::string& signal_list, bool _is_blind = false, bool _draw_ratio = true,
                bool _draw_bkg_errors = false)
        : anaDataReader(anaDataFileName), hist_name(_hist_name), is_blind(_is_blind), draw_ratio(_draw_ratio),
          draw_bkg_errors(_draw_bkg_errors)
    {
        gROOT->SetMustClean(kFALSE);

        std::istringstream s_channel(channel_name);
        s_channel >> channel;

        dataCategories = std::shared_ptr<analysis::DataCategoryCollection>(
                    new analysis::DataCategoryCollection(source_cfg, signal_list, channel));

        ParseIdSelection(id_selection);

        const std::string output_file_name = GenerateOutputFileName(output_path);

        printer = std::shared_ptr<root_ext::PdfPrinter>(new root_ext::PdfPrinter(output_file_name));
    }

    void Run()
    {
        for(analysis::EventCategory eventCategory : analysis::AllEventCategories) {
            if(!loop_options.eventCategory && eventCategory != meta_id.eventCategory) continue;
            for(analysis::EventSubCategory eventSubCategory : analysis::AllEventSubCategories) {
                if(!loop_options.eventSubCategory && eventSubCategory != meta_id.eventSubCategory) continue;
                for(analysis::EventRegion eventRegion : analysis::AllEventRegions) {
                    if(!loop_options.eventRegion && eventRegion != meta_id.eventRegion) continue;
                    for(analysis::EventEnergyScale eventEnergyScale : analysis::AllEventEnergyScales) {
                        if(!loop_options.eventEnergyScale && eventEnergyScale != meta_id.eventEnergyScale) continue;

                        const analysis::FlatAnalyzerDataMetaId_noName id(eventCategory, eventSubCategory, eventRegion,
                                                                         eventEnergyScale);
                        PrintStackedPlots(id, eventCategory);
                    }
                }
            }
        }
    }

private:
    void PrintStackedPlots(const analysis::FlatAnalyzerDataMetaId_noName& id, analysis::EventCategory eventCategory)
    {
        std::ostringstream ss_title;
        ss_title << id << ": " << hist_name;

        analysis::StackedPlotDescriptor stackDescriptor(ss_title.str(), false,
                                                        analysis::detail::ChannelNameMapLatex.at(channel),
                                                        eventCategoryMap.at(eventCategory),
                                                        draw_ratio, draw_bkg_errors);

        for(const analysis::DataCategory* category : dataCategories->GetAllCategories()) {
            if(!category->draw) continue;

            const analysis::FlatAnalyzerDataId full_id = id.MakeId(category->name);
            const auto histogram = anaDataReader.GetHistogram<TH1D>(full_id, channel, hist_name);
            if(!histogram) continue;

            if(category->IsSignal() && id.eventCategory != analysis::EventCategory::TwoJets_Inclusive)
                stackDescriptor.AddSignalHistogram(*histogram, category->title, category->color,
                                                   category->draw_sf);
            else if(category->IsBackground())
                stackDescriptor.AddBackgroundHistogram(*histogram, category->title, category->color);
            else if(category->IsData())
                stackDescriptor.AddDataHistogram(*histogram, category->title, is_blind,
                                                 GetBlindRegion(id.eventSubCategory, hist_name));
        }

        printer->PrintStack(stackDescriptor);
    }

    void ParseIdSelection(const std::string& id_selection)
    {
        static const char list_separator = '/';
        static const std::string wildcard = "*";
        static const std::string error_message =
                "Bad selection format. Expected eventCategory/eventSubCategory/eventRegion/eventEnergyScale.";
        const auto id_vector = analysis::ConfigReader::ParseOrderedParameterList(id_selection, true, list_separator);
        if(id_vector.size() != 4)
            throw analysis::exception(error_message);
        try {
            if(id_vector.at(0) != wildcard) {
                Parse(id_vector.at(0), meta_id.eventCategory);
                loop_options.eventCategory = false;
            }
            if(id_vector.at(1) != wildcard) {
                Parse(id_vector.at(1), meta_id.eventSubCategory);
                loop_options.eventSubCategory = false;
            }
            if(id_vector.at(2) != wildcard) {
                Parse(id_vector.at(2), meta_id.eventRegion);
                loop_options.eventRegion = false;
            }
            if(id_vector.at(3) != wildcard) {
                Parse(id_vector.at(3), meta_id.eventEnergyScale);
                loop_options.eventEnergyScale = false;
            }
        } catch(analysis::exception& e) {
            throw analysis::exception(e.message()) << "\n" << error_message;
        }
    }

    template<typename Type>
    static void Parse(const std::string& str, Type& t)
    {
        std::istringstream ss(str);
        ss >> t;
    }

    std::string GenerateOutputFileName(const std::string& output_path) const
    {
        static const std::string split = "_";
        std::ostringstream ss;
        ss << output_path;
        if(output_path.size() && output_path.at(output_path.size() - 1) != '/')
            ss << "/";
        ss << channel << split << hist_name;
        if(!loop_options.eventCategory)
            ss << split << meta_id.eventCategory;
        if(!loop_options.eventSubCategory)
            ss << split << meta_id.eventSubCategory;
        if(!loop_options.eventRegion)
            ss << split << meta_id.eventRegion;
        if(!loop_options.eventEnergyScale)
            ss << split << meta_id.eventEnergyScale;
        ss << ".pdf";
        return ss.str();
    }

    static const std::vector< std::pair<double, double> >& GetBlindRegion(analysis::EventSubCategory subCategory,
                                                                          const std::string& hist_name)
    {
        using analysis::EventSubCategory;
        using analysis::FlatAnalyzerData;
        using analysis::FlatAnalyzerData_semileptonic;

        static const std::vector< std::vector< std::pair<double, double> > > blindingRegions = {
            { { std::numeric_limits<double>::max(), std::numeric_limits<double>::lowest() } },
            { { 100, 150 } },
            { { 200, 400 } },
            { { 100, 150 }, { 450, 500 }, { 800, 850 }, { 1150, 1200 }, { 1500, 1550 } }
        };

        static const std::map<std::string, size_t> histogramsToBlind = {
            { FlatAnalyzerData_semileptonic::m_sv_Name(), 1 }, { FlatAnalyzerData::m_vis_Name(), 1 },
            { FlatAnalyzerData_semileptonic::m_bb_Name(), 1 }, { FlatAnalyzerData::m_ttbb_Name(), 2 },
            { FlatAnalyzerData::m_ttbb_kinfit_Name(), 2 }, { FlatAnalyzerData::m_bb_slice_Name(), 3 }
        };

        static const std::set<EventSubCategory> sidebandSubCategories = {
            EventSubCategory::OutsideMassWindow, EventSubCategory::KinematicFitConvergedOutsideMassWindow
        };

        const auto findRegionId = [&]() -> size_t {
            if(sidebandSubCategories.count(subCategory) || !histogramsToBlind.count(hist_name))
                return 0;
            return histogramsToBlind.at(hist_name);
        };

        const size_t regionId = findRegionId();
        if(regionId >= blindingRegions.size())
            throw analysis::exception("Bad blinding region index = ") << regionId;
        return blindingRegions.at(regionId);
    }

private:
    std::shared_ptr<analysis::DataCategoryCollection> dataCategories;
    analysis::FlatAnalyzerDataCollectionReader anaDataReader;
    std::shared_ptr<root_ext::PdfPrinter> printer;
    analysis::Channel channel;
    analysis::FlatAnalyzerDataMetaId_noName meta_id;
    LoopOptions loop_options;
    std::string hist_name;
    bool is_blind, draw_ratio, draw_bkg_errors;
};


namespace make_tools {
template<typename T>
struct Factory;

template<>
struct Factory<Print_Stack> {
    static Print_Stack* Make(int argc, char *argv[])
    {
        if(argc != 11)
            throw std::runtime_error("Invalid number of command line arguments.");

        int n = 0;
        const std::string source_cfg = argv[++n];
        const std::string anaDataFileName = argv[++n];
        const std::string output_path = argv[++n];
        const std::string channel_name = argv[++n];
        const std::string id_selection = argv[++n];
        const std::string hist_name = argv[++n];
        const std::string signal_list = argv[++n];
        const bool is_blind = ReadParameter<bool>(argv[++n], "is_blind");
        const bool draw_ratio = ReadParameter<bool>(argv[++n], "draw_ratio");
        const bool draw_bkg_errors = ReadParameter<bool>(argv[++n], "draw_bkg_errors");

        return new Print_Stack(source_cfg, anaDataFileName, output_path, channel_name, id_selection, hist_name,
                               signal_list, is_blind, draw_ratio, draw_bkg_errors);
    }

private:
    template<typename Param>
    static Param ReadParameter(const std::string& value, const std::string& param_name)
    {
        try {
            char c;
            Param param;
            std::istringstream ss_param(value);
            ss_param.exceptions(std::istream::failbit | std::istream::badbit);
            ss_param >> std::boolalpha;
            ss_param >> c;
            if(c != '@')
                throw std::istream::failure("Expected '@' before the parameter value.");
            ss_param >> param;
            return param;
        } catch(std::istream::failure& e) {
            std::ostringstream ss;
            ss << "Unable to parse parameter '" << param_name << "' from string '" << value << "'.\n" << e.what();
            throw std::runtime_error(ss.str());
        }
    }
};
} // make_tools

