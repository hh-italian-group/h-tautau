## Configuration file that defines the producer of ROOT-tuple for gen-particles.
## This file is part of https://github.com/hh-italian-group/h-tautau.

import FWCore.ParameterSet.Config as cms

genParticleBlock = cms.EDAnalyzer("GenParticleBlock",
  verbosity      = cms.int32(0),
  genParticleSrc = cms.InputTag('genParticles')
)
