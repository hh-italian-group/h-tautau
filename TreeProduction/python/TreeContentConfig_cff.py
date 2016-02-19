##  Configuration file that imports all other configurations related with ROOT-tuple production.
## This file is part of https://github.com/hh-italian-group/h-tautau.

import FWCore.ParameterSet.Config as cms

from HHbbTauTau.TreeProduction.EventBlock_cfi import eventBlock
from HHbbTauTau.TreeProduction.VertexBlock_cfi import vertexBlock
from HHbbTauTau.TreeProduction.JetBlock_cfi import jetBlock
from HHbbTauTau.TreeProduction.ElectronBlock_cfi import electronBlock
from HHbbTauTau.TreeProduction.METBlock_cfi import metBlock
from HHbbTauTau.TreeProduction.MuonBlock_cfi import muonBlock
from HHbbTauTau.TreeProduction.TauBlock_cfi import tauBlock
from HHbbTauTau.TreeProduction.PFCandBlock_cfi import pfCandBlock
from HHbbTauTau.TreeProduction.GenParticleBlock_cfi import genParticleBlock
from HHbbTauTau.TreeProduction.GenEventBlock_cfi import genEventBlock
from HHbbTauTau.TreeProduction.GenJetBlock_cfi import genJetBlock
from HHbbTauTau.TreeProduction.GenMETBlock_cfi import genMETBlock
from HHbbTauTau.TreeProduction.TriggerBlock_cfi import triggerBlock
from HHbbTauTau.TreeProduction.TriggerObjectBlock_cfi import triggerObjectBlock

