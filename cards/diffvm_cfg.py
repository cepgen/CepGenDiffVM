import Config.Core as cepgen
#from Config.pythia8_cff import pythia8 as hadroniser
from Config.PDG_cfi import PDG
from Config.logger_cfi import logger
from Config.generator_cfi import generator as _gen
from Integrators.vegas_cfi import vegas as integrator


DiffVMPDG = cepgen.Parameters(
    Jpsi      = Parameters(pdgid = 443),
    Upsilon1S = Parameters(pdgid = 553),
    proton    = Parameters(pdgid = 2212),
    neutron   = Parameters(pdgid = 2112),
    Upsilon2S = Parameters(pdgid = 100553),
    Upsilon3S = Parameters(pdgid = 200553),
)

logger.enabledModules += (
#    'DiffVM.*',
#    'EPA.*',
)

import Config.DiffVM_cff as diffvm

process = cepgen.Module('diffvm',
    processParameters = diffvm.defaultProcessParameters.clone(
        vmFlavour = DiffVMPDG.Upsilon1S,
        #protonMode = diffvm.BeamMode.StandardFragmentation,
    ),
    inKinematics = cepgen.Parameters(
        pdgIds = (PDG.electron, PDG.proton),
        pz = (27.55, 820.),
    ),
    outKinematics = cepgen.Parameters(
        w = (20.,),
        #invmass = ( 0., 6. ),
        mx = (1.07, 1000.),
    ),
)

generator = _gen.clone(
    numEvents = 25000,
    printEvery = 5000,
)
