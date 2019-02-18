/*! Classes that represent analysis objects at the tuple level.
This file is part of https://github.com/hh-italian-group/h-tautau. */

#pragma once

#include "AnalysisTools/Core/include/Tools.h"
#include "AnalysisTools/Core/include/AnalysisMath.h"
#include "AnalysisTypes.h"
#include "EventTuple.h"
#include "DiscriminatorIdResults.h"

namespace ntuple {

class TupleObject {
public:
    using MetType = analysis::MetType;
    using DiscriminatorWP = analysis::DiscriminatorWP;
    using DiscriminatorIdResults = analysis::DiscriminatorIdResults;
    using DiscriminatorResult = float;
    using Integer = int;
    using RealNumber = float;

    TupleObject(const ntuple::Event& _event);

protected:
    const Event* event;
};

class TupleLepton : public TupleObject {
public:
    TupleLepton(const ntuple::Event& _event, size_t _leg_id);
    const LorentzVectorM& p4() const;
    Integer charge() const;
    RealNumber dxy() const;
    RealNumber dz() const;
    RealNumber iso() const;
    Integer gen_match() const;

protected:
    size_t leg_id;
};

class TupleElectron : public TupleLepton {
public:
    explicit TupleElectron(const ntuple::Event& _event, size_t _leg_id = 1);
    DiscriminatorIdResults iso_wp() const { return DiscriminatorIdResults(); }
};

class TupleMuon : public TupleLepton {
public:
    explicit TupleMuon(const ntuple::Event& _event, size_t _leg_id = 1);
    DiscriminatorIdResults iso_wp() const { return DiscriminatorIdResults(); }
};

class TupleTau : public TupleLepton {
public:
    using IdKey = uint32_t;
    using TupleLepton::TupleLepton;
    using ValueKeyPair = std::pair<std::string, IdKey>;

public:
//     static ValueKeyPair GetNameKeyPair(const std::string& discriminator)
//     {
//         static std::map<std::string, IdKey> hashes;
//         auto iter = hashes.find(discriminator);
//         if(iter == hashes.end()) {
//             const IdKey key = analysis::tools::hash(discriminator);
//             iter = hashes.emplace(discriminator, key).first;
//         }
//         return *iter;
//     }
//
//     DiscriminatorResult tauID(const std::string& discriminator) const
//     {
//         const IdKey key = GetNameKeyPair(discriminator).second;
//         return _tauID(key, discriminator);
//     }
//
//     bool tauID(IdKey key, DiscriminatorResult& result) const
//     {
//         if(!tauIds.size()) {
//             const auto& keys = leg_id == 1 ? event->tauId_flags_1 :event->tauId_flags_2;
//             const auto& values = leg_id == 1 ? event->tauId_flags_1 :event->tauId_flags_2;
//             if(keys.size() != values.size())
//                 throw analysis::exception("Invalid tauID data");
//             for(size_t n = 0; n < keys.size(); ++n)
//                 tauIds[keys.at(n)] = values.at(n);
//         }
//         auto result_iter = tauIds.find(key);
//         bool has_result = result_iter != tauIds.end();
//         if(has_result)
//             result = result_iter->second;
//         return has_result;
//     }
//
    // DiscriminatorResult againstElectronMVA6(DiscriminatorWP wp) const
    // {
    //     std::ostringstream ss_name;
    //     ss_name << "againstElectron" << wp << "MVA6";
    //     return tauID(ss_name.str());
    // }
    //
    // DiscriminatorResult againstMuon3(DiscriminatorWP wp) const
    // {
    //     std::ostringstream ss_name;
    //     ss_name << "againstMuon" << wp << "3";
    //     return tauID(ss_name.str());
    // }
//
//     DiscriminatorResult byIsolationMVAraw(bool use_new_dm = false, bool use_lifetime = true) const
//     {
//         using IsoKey = std::tuple<bool, bool>;
//         static std::map<IsoKey, ValueKeyPair> keys;
//         auto iso_key = std::make_tuple(use_new_dm, use_lifetime);
//         auto iter = keys.find(iso_key);
//         if(iter == keys.end()) {
//             const std::string discriminator = GetByIsolationName(use_new_dm, use_lifetime, true);
//             iter = keys.emplace(iso_key, GetNameKeyPair(discriminator)).first;
//         }
//         return _tauID(iter->second.second, iter->second.first);
//     }
//
//     bool byIsolationMVA(DiscriminatorWP wp, bool use_new_dm = false, bool use_lifetime = true) const
//     {
//         using IsoKey = std::tuple<DiscriminatorWP, bool, bool>;
//         static std::map<IsoKey, ValueKeyPair> keys;
//         auto iso_key = std::make_tuple(wp, use_new_dm, use_lifetime);
//         auto iter = keys.find(iso_key);
//         if(iter == keys.end()) {
//             const std::string discriminator = GetByIsolationName(use_new_dm, use_lifetime, false, wp);
//             iter = keys.emplace(iso_key, GetNameKeyPair(discriminator)).first;
//         }
//         return _tauID(iter->second.second, iter->second.first) > 0.5;
//     }
//
//     DiscriminatorIdResults iso_wp() const
//     {
//         static const std::vector<DiscriminatorWP> available_wp = {
//             DiscriminatorWP::VLoose, DiscriminatorWP::Loose, DiscriminatorWP::Medium, DiscriminatorWP::Tight,
//             DiscriminatorWP::VTight, DiscriminatorWP::VVTight
//         };
//         DiscriminatorIdResults id_results;
//         for(auto wp : available_wp)
//             id_results.SetResult(wp, byIsolationMVA(wp));
//         return id_results;
//     }
//
// private:
//     DiscriminatorResult _tauID(IdKey key, const std::string& discriminator) const
//     {
//         DiscriminatorResult result;
//         if(!tauID(key, result))
//             throw analysis::exception("TauID discriminator '%1%' not found.") % discriminator;
//         return result;
//     }
//
//     static std::string GetByIsolationName(bool use_new_dm, bool use_lifetime, bool raw,
//                                           DiscriminatorWP wp = DiscriminatorWP::Medium)
//     {
//         const std::string dm_str = use_new_dm ? "new" : "old";
//         const std::string lt_str = use_lifetime ? "w" : "wo";
//         std::ostringstream ss_name;
//         ss_name << "by";
//         if(!raw)
//             ss_name << wp;
//         ss_name << "IsolationMVArun2v1DB" << dm_str << "DM" << lt_str << "LT";
//         if(raw)
//             ss_name << "raw";
//         return ss_name.str();
//     }

private:
    mutable std::map<IdKey, DiscriminatorResult> tauIds;
};

class TupleJet : public TupleObject {
public:
    TupleJet(const ntuple::Event& _event, size_t _jet_id);
    const LorentzVectorE& p4() const;
    bool PassPuId(DiscriminatorWP wp) const;
    DiscriminatorResult csv() const;
    DiscriminatorResult deepcsv() const;
    DiscriminatorResult deepFlavour() const;
    Integer hadronFlavour() const;
    RealNumber rawf() const;
    RealNumber resolution() const;
    ULong64_t triggerFilterMatch() const;
    size_t jet_index() const;
private:
    size_t jet_id;
};

class TupleSubJet : public TupleObject {
public:
    TupleSubJet(const ntuple::Event& _event, size_t _jet_id);
    const LorentzVectorE& p4() const;

private:
    size_t jet_id;
};

class TupleFatJet : public TupleObject {
public:
    enum class MassType { Pruned, Filtered, Trimmed, SoftDrop };

    TupleFatJet(const ntuple::Event& _event, size_t _jet_id);
    const LorentzVectorE& p4() const;
    float m(MassType massType) const;
    DiscriminatorResult jettiness(size_t tau_index) const;
    const std::vector<TupleSubJet>& subJets() const;

private:
    size_t jet_id;
    std::vector<TupleSubJet> sub_jets;
};

class TupleMet : public TupleObject {
public:
    using CovMatrix = analysis::SquareMatrix<2>;
    TupleMet(const ntuple::Event& _event, MetType _met_type);
    MetType type() const;
    const LorentzVectorM& p4() const;
    const CovMatrix& cov() const;
    RealNumber pt() const;
    RealNumber phi() const;

private:
    MetType met_type;
};

} // namespace ntuple
