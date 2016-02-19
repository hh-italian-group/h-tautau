/*! Enumaration of the PDG MC particle codes.
This file is part of https://github.com/hh-italian-group/h-tautau. */

#pragma once

#include <cmath>
#include <string>
#include <map>
#include <set>
#include <stdexcept>
#include <sstream>
#include <vector>
#include <iostream>
#include <limits>

namespace particles {
class ParticleCode {
private:
    typedef std::map<int, std::string> CodeToNameMap;
    typedef std::map<std::string, int> NameToCodeMap;
    static CodeToNameMap& Names()
    {
        static CodeToNameMap names;
        return names;
    }

    static NameToCodeMap& Codes()
    {
        static NameToCodeMap codes;
        return codes;
    }

    static bool& Freeze()
    {
        static bool freeze = false;
        return freeze;
    }

public:
    ParticleCode() : raw_code(0) {}
    ParticleCode(int _raw_code)
        : raw_code(_raw_code) {
        if(Names().find(raw_code) == Names().end()) {
            std::ostringstream ss;
            ss << "Unknown particle id " << raw_code << ".";
            std::cerr << ss.str() << std::endl;
            //throw std::runtime_error(ss.str());
            raw_code = 0;
        }
    }
    ParticleCode(const std::string& name)
    {
        const NameToCodeMap::const_iterator iter = Codes().find(name);
        if(iter == Codes().end()) {
            std::ostringstream ss;
            ss << "Unknown particle name '" << name << "'.";
            throw std::runtime_error(ss.str());
        }
        raw_code = iter->second;
    }
    ParticleCode(int _raw_code, const std::string& name)
        : raw_code(_raw_code) {
        if(Freeze())
            throw std::runtime_error("New particle code can not be defined.");
        if(Names().find(raw_code) != Names().end()) {
            std::ostringstream ss;
            ss << "Particle id " << raw_code << " is already defined.";
            throw std::runtime_error(ss.str());
        }
        if(Codes().find(name) != Codes().end()) {
            std::ostringstream ss;
            ss << "Particle name '" << name << "' is already taken.";
            throw std::runtime_error(ss.str());
        }

        Names()[raw_code] = name;
        Codes()[name] = raw_code;
        if(raw_code == std::numeric_limits<int>::max())
            Freeze() = true;
    }

    int RawCode() const { return raw_code; }
    const std::string& Name() const { return Names().find(raw_code)->second; }
    bool operator<(const ParticleCode& other) const { return raw_code < other.raw_code; }
    bool operator>(const ParticleCode& other) const { return raw_code > other.raw_code; }
    bool operator==(const ParticleCode& other) const { return raw_code == other.raw_code; }
    bool operator!=(const ParticleCode& other) const { return raw_code != other.raw_code; }
private:
    int raw_code;
};

#define PARTICLE(name, code) \
    static const ParticleCode name(code, #name)

PARTICLE(d, 1);
PARTICLE(u, 2);
PARTICLE(s, 3);
PARTICLE(c, 4);
PARTICLE(b, 5);
PARTICLE(t, 6);
PARTICLE(b_prime, 7);
PARTICLE(t_prime, 8);
PARTICLE(e, 11);
PARTICLE(nu_e, 12);
PARTICLE(mu, 13);
PARTICLE(nu_mu, 14);
PARTICLE(tau, 15);
PARTICLE(nu_tau, 16);
PARTICLE(tau_prime, 17);
PARTICLE(nu_tau_prime, 18);
PARTICLE(g, 21);
PARTICLE(gamma, 22);
PARTICLE(pi_zero, 111);
PARTICLE(rho_zero, 113);
PARTICLE(K_zero_L, 130);
PARTICLE(pi, 211);
PARTICLE(rho, 213);
PARTICLE(eta, 221);
PARTICLE(omega_782, 223);
PARTICLE(K_zero_S, 310);
PARTICLE(K_zero, 311);
PARTICLE(K, 321);
PARTICLE(K_star, 323);
PARTICLE(p, 2212);
PARTICLE(n, 2112);
PARTICLE(uu1, 2203);
PARTICLE(sigma_star_0, 3214);
PARTICLE(sigma_minus, 3112);
PARTICLE(MC_internal_92, 92);
PARTICLE(Delta_plusplus, 2224);
PARTICLE(Delta_plus, 2214);
PARTICLE(Delta_zero, 2114);
PARTICLE(Delta_minus, 1114);
PARTICLE(Higgs, 25);
PARTICLE(K_star_0, 313);
PARTICLE(B_0, 511);
PARTICLE(Lambda, 3122);
PARTICLE(eta_prime, 331);
PARTICLE(ud0, 2101);
PARTICLE(xi_minus, 3312);
PARTICLE(phi_1020, 333);
PARTICLE(sigma_plus, 3222);
PARTICLE(sigma_0, 3212);
PARTICLE(ud1, 2103);
PARTICLE(D_star_0, 423);
PARTICLE(D_0, 421);
PARTICLE(sigma_star_minus, 3114);
PARTICLE(sigma_star_plus, 3224);
PARTICLE(D_plus, 411);
PARTICLE(a1_plus, 20213);
PARTICLE(D_star_plus, 413);
PARTICLE(B_star_0, 513);
PARTICLE(xi_0, 3322);
PARTICLE(B_star_plus, 523);
PARTICLE(Lambda_b_zero, 5122);
PARTICLE(Lambda_c_plus, 4122);
PARTICLE(D_s_plus, 431);
PARTICLE(D_s_star_plus, 433);
PARTICLE(xi_star_0, 3324);
PARTICLE(omega_minus, 3334);
PARTICLE(xi_star_minus, 3314);
PARTICLE(sigma_c_star_0, 4114);
PARTICLE(B_plus, 521);
PARTICLE(MC_internal_91, 91);
PARTICLE(sigma_c_star_plus, 4214);
PARTICLE(sigma_c_plus, 4212);
PARTICLE(xi_c_prime_0, 4312);
PARTICLE(D_s2_star_plus, 435);
PARTICLE(sigma_c_star_plusplus, 4224);
PARTICLE(K1_plus, 10323);
PARTICLE(K1_0, 20313);
PARTICLE(B_s_star_0, 533);
PARTICLE(B_s_0, 531);
PARTICLE(sigma_b_star_plus, 5224);
PARTICLE(xi_c_zero, 4132);
PARTICLE(xi_b_minus, 5132);
PARTICLE(sigma_b_0, 5212);
PARTICLE(sigma_b_star_minus, 5114);
PARTICLE(sigma_c_plusplus, 4222);
PARTICLE(xi_b_0, 5232);
PARTICLE(xi_b_star_0, 5324);
PARTICLE(D1_0_H, 20423);
PARTICLE(D2_star_plus, 415);
PARTICLE(xi_c_plus, 4232);
PARTICLE(D1_plus_H, 20413);
PARTICLE(sigma_b_star_0, 5214);
PARTICLE(xi_c_star_zero, 4314);
PARTICLE(D0_star_plus, 10411);
PARTICLE(D1_plus, 10413);
PARTICLE(D1_0, 10423);
PARTICLE(J_psi, 443);
PARTICLE(D_s0_star_plus, 10431);
PARTICLE(sigma_c_0, 4112);
PARTICLE(D2_star_0, 425);
PARTICLE(D_s1_plus_2536, 10433);
PARTICLE(f0, 10221);
PARTICLE(sigma_b_plus, 5222);
PARTICLE(xi_c_star_plus, 4324);
PARTICLE(p_diffr_plus, 9902210);
PARTICLE(deuteron, 1000010020);
PARTICLE(eta_c_1S, 441);
PARTICLE(dd1, 1103);
PARTICLE(su0, 3201);
PARTICLE(su1, 3203);
PARTICLE(D0_star_0, 10421);
PARTICLE(xi_b_prime_minus, 5312);
PARTICLE(chi_c1_1P, 20443);
PARTICLE(sd0, 3101);
PARTICLE(D_s1_plus_2460,20433);
PARTICLE(sigma_b_minus, 5112);
PARTICLE(junction, 88);
PARTICLE(W_plus, 24);
PARTICLE(Z, 23);
PARTICLE(xi_c_plus_prime, 4322);
PARTICLE(upsilon_1s, 553);
PARTICLE(B_c_plus, 541);
PARTICLE(B_c_star_plus, 543);
PARTICLE(xi_b_prime_0, 5322);
PARTICLE(xi_b_minus_star, 5314);
PARTICLE(omega_b_minus_star, 5334);
PARTICLE(omega_b_minus, 5332);
PARTICLE(eta_b_1S, 551);
PARTICLE(omega_c_star_0, 4334);
PARTICLE(omega_c_0, 4332);
PARTICLE(sd_1, 3103);
PARTICLE(MSSM_H, 35);
PARTICLE(MSSM_A, 36);

PARTICLE(pi_zero_1300, 100111);
PARTICLE(rho_zero_1450, 100113);
PARTICLE(pi_plus_1300, 100211);
PARTICLE(rho_plus_1450, 100213);
PARTICLE(eta_1295, 100221);
PARTICLE(omega_1420, 100223);
PARTICLE(K_zero_1460, 100311);
PARTICLE(K_zero_star_1410, 100313);
PARTICLE(K_plus_1460, 100321);
PARTICLE(K_plus_star_1410, 100323);
PARTICLE(eta_1475, 100331);
PARTICLE(phi_1680, 100333);
PARTICLE(eta_c_2S, 100441);
PARTICLE(psi_2S, 100443);
PARTICLE(a_zero_zero_1450, 10111);
PARTICLE(b_one_zero_1235, 10113);
PARTICLE(pi_two_zero_1670, 10115);
PARTICLE(a_zero_plus_1450, 10211);
PARTICLE(b_one_plus_1235, 10213);
PARTICLE(pi_two_plus_1670, 10215);
PARTICLE(h_one_1170, 10223);
PARTICLE(eta_two_1645, 10225);
PARTICLE(K_zero_star_zero_1430, 10311);
PARTICLE(K_one_zero_1270, 10313);
PARTICLE(K_two_zero_1770, 10315);
PARTICLE(K_zero_star_plus_1430, 10321);
PARTICLE(K_two_plus_1770, 10325);
PARTICLE(f_zero_1710, 10331);
PARTICLE(h_one_1380, 10333);
PARTICLE(eta_two_1870, 10335);
PARTICLE(chi_c0_1P, 10441);
PARTICLE(B_zero_star_zero, 10511);
PARTICLE(B_one_zero_L, 10513);
PARTICLE(B_zero_star_plus, 10521);
PARTICLE(B_one_plus_L, 10523);
PARTICLE(B_s0_star_zero, 10531);
PARTICLE(B_s1_zero_L, 10533);
PARTICLE(a_two_zero_1320, 115);
PARTICLE(rho_three_zero_1690, 117);
PARTICLE(P11_N_zero_1440, 12112);
PARTICLE(P11_N_plus_1440, 12212);
PARTICLE(sigma_minus_1660, 13112);
PARTICLE(lambda_1404, 13122);
PARTICLE(sigma_zero_1660, 13212);
PARTICLE(sigma_plus_1660, 13222);
PARTICLE(xi_minus_1690, 13312);
PARTICLE(xi_zero_1690, 13322);
PARTICLE(lambda_c_2593, 14122);
PARTICLE(lambda_b1_zero, 15122);
PARTICLE(a_one_zero_1260, 20113);
PARTICLE(f_one_1285, 20223);
PARTICLE(K_two_zero_1820, 20315);
PARTICLE(K_one_plus_1400, 20323);
PARTICLE(K_two_plus_1820, 20325);
PARTICLE(f_one_1420, 20333);
PARTICLE(B_one_zero_H, 20513);
PARTICLE(B_one_plus_H, 20523);
PARTICLE(B_s1_zero_H, 20533);
PARTICLE(a_two_plus_1320, 215);
PARTICLE(rho_three_plus_1690, 217);
PARTICLE(S11_N_zero_1535, 22112);
PARTICLE(S11_N_plus_1535, 22212);
PARTICLE(f_two_1270, 225);
PARTICLE(omega_three_1670, 227);
PARTICLE(sigma_minus_1750, 23112);
PARTICLE(lambda_zero_1600, 23122);
PARTICLE(sigma_zero_1750, 23212);
PARTICLE(sigma_plus_1750, 23222);
PARTICLE(xi_minus_1950, 23312);
PARTICLE(xi_zero_1950, 23322);
PARTICLE(rho_zero_1700, 30113);
PARTICLE(rho_plus_1700, 30213);
PARTICLE(omega_1650, 30223);
PARTICLE(K_star_zero_1680, 30313);
PARTICLE(K_star_plus_1680, 30323);
PARTICLE(lambda_zero_1520, 3124);
PARTICLE(K_two_star_zero, 315);
PARTICLE(K_three_star_zero_1780, 317);
PARTICLE(K_two_star_plus, 325);
PARTICLE(K_three_star_plus_1780, 327);
PARTICLE(lambda_zero_1670, 33122);
PARTICLE(f_two_prime_1525, 335);
PARTICLE(phi_three_1850, 337);
PARTICLE(chi_c2_1P, 445);
PARTICLE(B_two_star_zero, 515);
PARTICLE(B_two_star_plus, 525);
PARTICLE(B_s2_star_zero, 535);
PARTICLE(MC_internal_use_81, 81);
PARTICLE(MC_internal_use_82, 82);
PARTICLE(a_zero_zero_980, 9000111);
PARTICLE(a_zero_plus_980, 9000211);
PARTICLE(f_zero_600, 9000221);
PARTICLE(f_zero_980, 9010221);
PARTICLE(eta_1405, 9020221);
PARTICLE(f_zero_1500, 9030221);
PARTICLE(psi_3770, 30443);
PARTICLE(xi_c_zero_2790, 14312);
PARTICLE(xi_c_plus_2790, 14322);
PARTICLE(xi_b1_zero, 15322);
PARTICLE(h_c_1P, 10443);
PARTICLE(lambda_c_plus_2625, 4124);
PARTICLE(xi_b1_minus, 15312);

PARTICLE(NONEXISTENT, std::numeric_limits<int>::max());

#undef PARTICLE

static const std::set< particles::ParticleCode > neutrinos = { nu_e, nu_mu, nu_tau };


enum Status {
    FinalStateParticle = 1, Decayed_or_fragmented = 2, HardInteractionProduct = 3,
    Generator_dependet_41 = 41, Generator_dependet_51 = 51
};

enum ParticleType {
    Particle = 1, AntiParticle = -1
};

namespace detail {
template<typename T>
struct NameCollection;
} // detail

template<typename T>
class NameProvider {
public:
    typedef std::map<int, std::string> NameMap;

    static const std::string& Name(const T& value)
    {
        const NameMap::const_iterator iter = detail::NameCollection<T>::Names().find(value);
        Check(iter, (int) value);
        return iter->second;
    }

    static T Convert(int id)
    {
        const NameMap::const_iterator iter = detail::NameCollection<T>::Names().find(id);
        Check(iter, id);
        return (T) id;
    }

    static T Convert(const std::string& name)
    {
        NameMap::const_iterator iter = detail::NameCollection<T>::Names().begin();
        for(; iter != detail::NameCollection<T>::Names().end(); ++iter) {
            if(iter->second == name)
                break;
        }
        Check(iter, name);
        return (T) iter->first;
    }

private:
    NameProvider() {}

    template<typename Value>
    static void Check(const NameMap::const_iterator& iter, const Value& value)
    {
        if(iter == detail::NameCollection<T>::Names().end())
        {
            std::ostringstream ss;
            ss << "Unknown " << detail::NameCollection<T>::TypeName() << " '" << value << "'.";
            throw std::runtime_error(ss.str());
        }
    }
};

namespace detail {
template<>
struct NameCollection<Status> {
    static const char* TypeName() { return "particle status"; }
    static const NameProvider<Status>::NameMap& Names()
    {
        static NameProvider<Status>::NameMap names;
        if(!names.size())
        {
            names[FinalStateParticle] = "FinalStateParticle";
            names[Decayed_or_fragmented] = "Decayed_or_fragmented";
            names[HardInteractionProduct] = "HardInteractionProduct";
            names[Generator_dependet_41] = "Generator_dependet_41";
            names[Generator_dependet_51] = "Generator_dependet_51";
        }
        return names;
    }
};

template<>
struct NameCollection<ParticleType> {
    static const char* TypeName() { return "particle type"; }
    static const NameProvider<ParticleType>::NameMap& Names()
    {
        static NameProvider<ParticleType>::NameMap names;
        if(!names.size())
        {
            names[Particle] = "";
            names[AntiParticle] = "anti-";
        }
        return names;
    }
};
} // detail

struct PdgParticle {
    ParticleCode Code;
    ParticleType Type;
    PdgParticle() : Code(NONEXISTENT), Type(Particle) {}
    PdgParticle(int id) : Code(std::abs(id)), Type(id >= 0 ? Particle : AntiParticle) { }
    std::string Name() const { return NameProvider<ParticleType>::Name(Type) + Code.Name(); }
    bool operator<(const PdgParticle& other) const
    {
        if(Code < other.Code) return true;
        if(Code > other.Code) return false;
        return Type < other.Type;
    }

    bool operator!=(const PdgParticle& other) const
    {
        return Code != other.Code || Type != other.Type;
    }

    bool operator==(const PdgParticle& other) const
    {
        return Code == other.Code && Type == other.Type;
    }

    int ToInteger() const
    {
        return Type == Particle ? Code.RawCode() : - Code.RawCode();
    }
};

typedef std::vector<particles::ParticleCode> ParticleCodes;
typedef std::vector<ParticleCodes> ParticleCodes2D;

std::ostream& operator<<(std::ostream& s, const Status& particleStatus){
    s << particles::NameProvider<particles::Status>::Name(particleStatus);
    return s;
}

std::ostream& operator<<(std::ostream& s, const ParticleCode& code){
    s << code.Name();
    return s;
}

std::ostream& operator<<(std::ostream& s, const PdgParticle& particle){
    s << particle.Name();
    return s;
}

} // particles
