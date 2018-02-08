/*! Classes that represent analysis objects at the tuple level.
This file is part of https://github.com/hh-italian-group/h-tautau. */

#pragma once

#include "AnalysisTools/Core/include/Tools.h"
#include "AnalysisMath.h"
#include "AnalysisTypes.h"
#include "EventTuple.h"

namespace ntuple {

class TupleObject {
public:
    using MetType = analysis::MetType;
    using DiscriminatorWP = analysis::DiscriminatorWP;
    using DiscriminatorResult = float;
    using Integer = int;
    using RealNumber = float;

    TupleObject(const ntuple::Event& _event) : event(&_event) {}

protected:
    const Event* event;
};

class TupleLepton : public TupleObject {
public:    
    TupleLepton(const ntuple::Event& _event, size_t _leg_id)
        : TupleObject(_event), leg_id(_leg_id)
    {
        if(leg_id < 1 || leg_id > 2)
            throw analysis::exception("Invalid leg id = %1%.") % leg_id;
    }

    const LorentzVectorM& p4() const { return leg_id == 1 ? event->p4_1 : event->p4_2; }
    Integer charge() const { return leg_id == 1 ? event->q_1 : event->q_2; }
    RealNumber dxy() const { return leg_id == 1 ? event->dxy_1 : event->dxy_2; }
    RealNumber dz() const { return leg_id == 1 ? event->dz_1 : event->dz_2; }
    RealNumber iso() const { return leg_id == 1 ? event->iso_1 : event->iso_2; }
    Integer gen_match() const { return leg_id == 1 ? event->gen_match_1 : event->gen_match_2; }

protected:
    size_t leg_id;
};

class TupleElectron : public TupleLepton {
public:
    explicit TupleElectron(const ntuple::Event& _event, size_t _leg_id = 1) : TupleLepton(_event, _leg_id) {}
};

class TupleMuon : public TupleLepton {
public:
    explicit TupleMuon(const ntuple::Event& _event, size_t _leg_id = 1) : TupleLepton(_event, _leg_id) {}
};

class TupleTau : public TupleLepton {
public:
    using IdKey = uint32_t;
    using TupleLepton::TupleLepton;
    using ValueKeyPair = std::pair<std::string, IdKey>;

public:
    static ValueKeyPair GetNameKeyPair(const std::string& discriminator)
    {
        static std::map<std::string, IdKey> hashes;
        auto iter = hashes.find(discriminator);
        if(iter == hashes.end()) {
            const IdKey key = analysis::tools::hash(discriminator);
            iter = hashes.emplace(discriminator, key).first;
        }
        return *iter;
    }

    DiscriminatorResult tauID(const std::string& discriminator) const
    {
        const IdKey key = GetNameKeyPair(discriminator).second;
        return _tauID(key, discriminator);
    }

    bool tauID(IdKey key, DiscriminatorResult& result) const
    {
        if(!tauIds.size()) {
            const auto& keys = leg_id == 1 ? event->tauId_keys_1 :event->tauId_keys_2;
            const auto& values = leg_id == 1 ? event->tauId_values_1 :event->tauId_values_2;
            if(keys.size() != values.size())
                throw analysis::exception("Invalid tauID data");
            for(size_t n = 0; n < keys.size(); ++n)
                tauIds[keys.at(n)] = values.at(n);
        }
        auto result_iter = tauIds.find(key);
        bool has_result = result_iter != tauIds.end();
        if(has_result)
            result = result_iter->second;
        return has_result;
    }

    DiscriminatorResult againstElectronMVA6(DiscriminatorWP wp) const
    {
        std::ostringstream ss_name;
        ss_name << "againstElectron" << wp << "MVA6";
        return tauID(ss_name.str());
    }

    DiscriminatorResult againstMuon3(DiscriminatorWP wp) const
    {
        std::ostringstream ss_name;
        ss_name << "againstMuon" << wp << "3";
        return tauID(ss_name.str());
    }

    DiscriminatorResult byIsolationMVAraw(bool use_new_dm = false, bool use_lifetime = true) const
    {
        using IsoKey = std::tuple<bool, bool>;
        static std::map<IsoKey, ValueKeyPair> keys;
        auto iso_key = std::make_tuple(use_new_dm, use_lifetime);
        auto iter = keys.find(iso_key);
        if(iter == keys.end()) {
            const std::string discriminator = GetByIsolationName(use_new_dm, use_lifetime, true);
            iter = keys.emplace(iso_key, GetNameKeyPair(discriminator)).first;
        }
        return _tauID(iter->second.second, iter->second.first);
    }

    bool byIsolationMVA(DiscriminatorWP wp, bool use_new_dm = false, bool use_lifetime = true) const
    {
        using IsoKey = std::tuple<DiscriminatorWP, bool, bool>;
        static std::map<IsoKey, ValueKeyPair> keys;
        auto iso_key = std::make_tuple(wp, use_new_dm, use_lifetime);
        auto iter = keys.find(iso_key);
        if(iter == keys.end()) {
            const std::string discriminator = GetByIsolationName(use_new_dm, use_lifetime, false, wp);
            iter = keys.emplace(iso_key, GetNameKeyPair(discriminator)).first;
        }
        return _tauID(iter->second.second, iter->second.first) > 0.5;
    }

private:
    DiscriminatorResult _tauID(IdKey key, const std::string& discriminator) const
    {
        DiscriminatorResult result;
        if(!tauID(key, result))
            throw analysis::exception("TauID discriminator '%1%' not found.") % discriminator;
        return result;
    }

    static std::string GetByIsolationName(bool use_new_dm, bool use_lifetime, bool raw,
                                          DiscriminatorWP wp = DiscriminatorWP::Medium)
    {
        const std::string dm_str = use_new_dm ? "new" : "old";
        const std::string lt_str = use_lifetime ? "w" : "wo";
        std::ostringstream ss_name;
        ss_name << "by";
        if(!raw)
            ss_name << wp;
        ss_name << "IsolationMVArun2v1DB" << dm_str << "DM" << lt_str << "LT";
        if(raw)
            ss_name << "raw";
        return ss_name.str();
    }

private:
    mutable std::map<IdKey, DiscriminatorResult> tauIds;
};

class TupleJet : public TupleObject {
public:
    TupleJet(const ntuple::Event& _event, size_t _jet_id)
        : TupleObject(_event), jet_id(_jet_id)
    {
        if(jet_id >= event->jets_p4.size())
            throw analysis::exception("Jet id = %1% is out of range.") % jet_id;
    }

    const LorentzVectorE& p4() const { return event->jets_p4.at(jet_id); }
    DiscriminatorResult mva() const { return event->jets_mva.at(jet_id); }
    DiscriminatorResult csv() const { return event->jets_csv.at(jet_id); }
    DiscriminatorResult deepcsv() const { return event->jets_deepCsv_b.at(jet_id); }
    Integer partonFlavour() const { return event->jets_partonFlavour.at(jet_id); }
    Integer hadronFlavour() const { return event->jets_hadronFlavour.at(jet_id); }
    RealNumber rawf() const { return event->jets_rawf.at(jet_id); }

private:
    size_t jet_id;
};

class TupleSubJet : public TupleObject {
public:
    TupleSubJet(const ntuple::Event& _event, size_t _jet_id)
        : TupleObject(_event), jet_id(_jet_id)
    {
        if(jet_id >= event->subJets_p4.size())
            throw analysis::exception("Fat sub-jet id = %1% is out of range.") % jet_id;
    }

    const LorentzVectorE& p4() const { return event->subJets_p4.at(jet_id); }
    DiscriminatorResult csv() const { return event->subJets_csv.at(jet_id); }

private:
    size_t jet_id;
};

class TupleFatJet : public TupleObject {
public:
    enum class MassType { Pruned, Filtered, Trimmed, SoftDrop };

    TupleFatJet(const ntuple::Event& _event, size_t _jet_id)
        : TupleObject(_event), jet_id(_jet_id)
    {
        if(jet_id >= event->fatJets_p4.size())
            throw analysis::exception("Fat jet id = %1% is out of range.") % jet_id;

        for(size_t n = 0; n < event->subJets_p4.size(); ++n) {
            if(event->subJets_parentIndex.at(n) == jet_id)
                sub_jets.push_back(TupleSubJet(_event, n));
        }
    }

    const LorentzVectorE& p4() const { return event->fatJets_p4.at(jet_id); }
    DiscriminatorResult csv() const { return event->fatJets_csv.at(jet_id); }

    float m(MassType massType) const
    {
        if(massType == MassType::Pruned) return event->fatJets_m_pruned.at(jet_id);
        if(massType == MassType::SoftDrop) return event->fatJets_m_softDrop.at(jet_id);
        throw analysis::exception("Unsupported fat jet mass type");
    }

    DiscriminatorResult n_subjettiness(size_t tau_index) const
    {
        if(tau_index == 1) return event->fatJets_n_subjettiness_tau1.at(jet_id);
        if(tau_index == 2) return event->fatJets_n_subjettiness_tau2.at(jet_id);
        if(tau_index == 3) return event->fatJets_n_subjettiness_tau3.at(jet_id);
        throw analysis::exception("Unsupported tau index = %1% for fat jet subjettiness.") % tau_index;
    }

    const std::vector<TupleSubJet>& subJets() const { return sub_jets; }

private:
    size_t jet_id;
    std::vector<TupleSubJet> sub_jets;
};

class TupleMet : public TupleObject {
public:
    using CovMatrix = analysis::SquareMatrix<2>;
    TupleMet(const ntuple::Event& _event, MetType _met_type)
        : TupleObject(_event), met_type(_met_type)
    {
        static const std::set<MetType> supported_types = { MetType::PF, MetType::MVA, MetType::PUPPI };
        if(!supported_types.count(met_type))
            throw analysis::exception("Unsupported met type.");
    }

    MetType type() const { return met_type; }

    const LorentzVectorM& p4() const
    {
        return event->pfMET_p4;
    }

    const CovMatrix& cov() const
    {
        return event->pfMET_cov;
    }

    RealNumber pt() const { return p4().pt(); }
    RealNumber phi() const { return p4().phi(); }

private:
    MetType met_type;
};

} // namespace ntuple
