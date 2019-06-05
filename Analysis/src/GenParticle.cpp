/*! Definiton of GenParticle class.
This file is part of https://github.com/hh-italian-group/h-tautau. */

#include "h-tautau/Analysis/include/GenParticle.h"

namespace analysis {

GenParticle::GenParticle(const ntuple::Event& events, size_t n)
    : index(n)
{
    if(n >= events.genParticles_p4.size())
        throw std::runtime_error("events index is out of range.");

    pdg = events.genParticles_pdg.at(n);
    status = events.genParticles_status.at(n);
    genStatusFlags = GenStatusFlags(events.genParticles_statusFlags.at(n));
    momentum = events.genParticles_p4.at(n);
    vertex = events.genParticles_vertex.at(n);
}

GenEvent::GenEvent(const ntuple::Event& event)
{
    genParticles.reserve(event.genParticles_p4.size());
    for(size_t n = 0; n < event.genParticles_p4.size(); ++n)
        genParticles.emplace_back(event,n);

    for (size_t relIndex = 0; relIndex < event.genParticles_rel_mIndex.size(); relIndex++) {
        int motherIndex = event.genParticles_rel_mIndex.at(relIndex);
        int daughterIndex = event.genParticles_rel_pIndex.at(relIndex);

        GenParticle* mother = &genParticles.at(static_cast<size_t>(motherIndex));
        GenParticle* daughter = &genParticles.at(static_cast<size_t>(daughterIndex));

        mother->daughters.push_back(daughter);
        daughter->mothers.push_back(mother);
    }

    for(const GenParticle& genParticle : genParticles ){
        if (!genParticle.mothers.size())
            primaryParticles.insert(&genParticle);


         if (genParticle.genStatusFlags.isPrompt())
            particleCodeMap[std::abs(genParticle.pdg)].insert(&genParticle);
    }
}

GenParticleSet GenEvent::GetParticles(int particle_pgd, bool requireIsLastCopy) const
{
    GenParticleSet results;
    const ParticleCodeMap::const_iterator code_iter = particleCodeMap.find(particle_pgd);

    if (code_iter != particleCodeMap.end()){
        for (const GenParticle* particle : code_iter->second){
            if(requireIsLastCopy && !particle->genStatusFlags.isLastCopy()) continue;
            results.insert(particle);
        }
    }
    return results;
}

void GenEvent::GetTypesParticles(const std::set<particles::ParticleCode::ParticleType>& type_names, const GenParticle* mother, GenParticleSet& result) const
{
    if(mother->genStatusFlags.isPrompt() && mother->genStatusFlags.isLastCopy() && particle_types.count(mother->pdg)
        && type_names.count(particle_types.at(mother->pdg))) {
//        if(shared && areParented(mother,mother->daughters.at(0)) && areParented(mother,mother->daughters.at(1)) )
        result.insert(mother);
    } else {
        for(const GenParticle* daughter : mother->daughters)
            GetTypesParticles(type_names, daughter, result);
    }
}

bool GenEvent::areParented(const GenParticle* daughter, const GenParticle* mother) const
{
    for (size_t i = 0; i < daughter->mothers.size(); i++) {
        const GenParticle* mother_ = daughter->mothers.at(i);
        if(mother_ == mother)
            return true;
        if(areParented(mother_,mother))
            return true;
    }
    return false;
}

std::string& GenEvent::GetParticleName(int pdgId) const
{
    auto iter = particle_names.find(pdgId);
    if(iter == particle_names.end()) throw exception("Name not found for particle with pdgId = %1%") % pdgId;
    return iter->second;
}

void GenEvent::intializeNames(const std::string& fileName)
{
    particle_names.clear();
    particle_types.clear();
    particle_charge.clear();
    std::ifstream f(fileName.c_str());
    if(!f.is_open() )
        throw analysis::exception("Unable to read the configuration file ");

    while(f.good()) {
        std::string line;
        std::getline(f, line);
        if(line.empty() || line[0] == '#')
            continue;
        auto pdg_name_type = SplitValueList(line, false, ",");
        if(pdg_name_type.size() < 2)
            throw exception("Invalid particle definition (1): '%1%'") % line;

        const int pdgId = Parse<int>(pdg_name_type.at(0));
        const std::string& name = pdg_name_type.at(1);
        if(particle_names.count(pdgId))
            throw exception("Duplicated definition of particle with pdgId = %1%") % pdgId;
        particle_names[pdgId] = name;

        if(pdg_name_type.size() > 2){

            const std::string& type_str = pdg_name_type.at(2);

            std::map<std::string, particles::ParticleCode::ParticleType> ParticleTypeMap;
            if(type_str == "baryon")
                ParticleTypeMap[type_str] = particles::ParticleCode::ParticleType::baryon;
            else if(type_str == "meson")
                ParticleTypeMap[type_str] = particles::ParticleCode::ParticleType::meson;

            if(ParticleTypeMap.count(type_str)){
                if(particle_types.count(pdgId))
                    throw exception("Duplicated definition of particle with pdgId = %1%") % pdgId;
                particle_types[pdgId] = ParticleTypeMap.at(type_str);
            }
        }

        if(pdg_name_type.size() == 4){
            const int& charge = Parse<int>(pdg_name_type.at(3));
            if(particle_charge.count(pdgId))
                throw exception("Duplicated definition of particle with pdgId = %1%") % pdgId;
            particle_charge[pdgId] = charge;

        }
    }
}

void GenEvent::PrintChain(const GenParticle* particle, const std::string& pre) const //, unsigned iteration = 0)
{
    const int pdgParticle = particle->pdg;
    const auto particleName = GetParticleName(pdgParticle);
    const int particleStatus = particle->status;
    const LorentzVectorM genParticle_momentum = particle->momentum;
    const std::bitset<15> genStatusFlags_ = particle->genStatusFlags.flags ;
    // for (unsigned n = 0; n < iteration; ++n)
    //     std::cout << "  ";
    auto mothers_index = particle->mothers.size() > 0 ?  particle->mothers.at(0)->index : 0;
    std::cout << particleName                              << " <" << pdgParticle
              << "> pt=" << genParticle_momentum.Pt()      << " eta=" << genParticle_momentum.Eta()
              << " phi=" << genParticle_momentum.Phi()     << " E=" << genParticle_momentum.E()
              << " m=" << genParticle_momentum.M()         << " index=" << particle->index
              << " mother_index=" << mothers_index;
    if(particle->mothers.size() > 0){
        for(size_t index_mother = 1; index_mother < particle->mothers.size(); ++index_mother)
            std::cout  << "," << particle->mothers.at(index_mother)->index;
    }

    std::cout << " vertex=" << particle->vertex << " status=" << particleStatus
              << " statusFlags=" << genStatusFlags_ << std::endl;

    for(unsigned n = 0; n < particle->daughters.size(); ++n) {
        const GenParticle* daughter = particle->daughters.at(n);
        std::cout << pre << "+-> ";
        const char pre_first = n == particle->daughters.size() -1 ? ' ' : '|';
        const std::string pre_d = pre + pre_first + "   ";
        PrintChain(daughter, pre_d);
    }
}

void GenEvent::Print() const
{
    for (const GenParticle* particle : primaryParticles) {
        GenEvent::PrintChain(particle, "");
    }
}

void GenEvent::FindFinalStateDaughters(const GenParticle& particle, std::set<const GenParticle*>& daughters,
                             const std::set<int>& pdg_to_exclude)
{
    if(!particle.daughters.size()) {
        const int abs_pdg = std::abs(particle.pdg);
        if(!pdg_to_exclude.count(abs_pdg))
            daughters.insert(&particle);
    } else {
        for(const auto& daughter : particle.daughters)
            FindFinalStateDaughters(*daughter, daughters, pdg_to_exclude);
    }
}


LorentzVectorM GenEvent::GetFinalStateMomentum(const GenParticle& particle, std::vector<const GenParticle*>& visible_daughters,
                                   bool excludeInvisible, bool excludeLightLeptons)
{
    using pair = std::pair<bool, bool>;
    static const std::set<int> empty = {};

    static const std::map<pair, const std::set<int>*> to_exclude {
        { pair(false, false), &empty }, { pair(true, false), &particles::neutrinos },
        { pair(false, true), &particles::light_leptons }, { pair(true, true), &particles::light_and_invisible },
    };

    std::set<const GenParticle*> daughters_set;

    FindFinalStateDaughters(particle, daughters_set, *to_exclude.at(pair(excludeInvisible, false)));
    visible_daughters.clear();
    visible_daughters.insert(visible_daughters.begin(), daughters_set.begin(), daughters_set.end());

    LorentzVectorM p4;
    for(auto daughter : visible_daughters) {
        if(excludeLightLeptons && particles::light_leptons.count(std::abs(daughter->pdg))
            && daughter->genStatusFlags.isDirectTauDecayProduct()) continue;
        p4 += daughter->momentum;
    }
    return p4;
}

std::map<int, std::string> GenEvent::particleNames() {return particle_names; }

std::map<int, int> GenEvent::particleCharge() {return particle_charge; }


std::map<int, std::string> GenEvent::particle_names;
std::map<int, particles::ParticleCode::ParticleType> GenEvent::particle_types;
std::map<int,int> GenEvent::particle_charge;

} //analysis
