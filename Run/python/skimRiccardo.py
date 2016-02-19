## Definition of the function to skim events from CERN sync tuples.
## This file is part of https://github.com/hh-italian-group/h-tautau.

def skimmingIso(tree):
  if ( tree.byCombinedIsolationDeltaBetaCorrRaw3Hits_1<1 and tree.byCombinedIsolationDeltaBetaCorrRaw3Hits_2 < 1
       and tree.againstElectronLooseMVA3_2 ) : return True

  else : return False
