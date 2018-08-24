import Config.Core as cepgen
#from Config.pythia8_cff import pythia8 as hadroniser
from Config.logger_cfi import logger
from Config.generator_cfi import generator as _gen
from Integrators.vegas_cfi import vegas as integrator


logger.enabledModules += (
#    'DiffVM.*',
#    'EPA.*',
)

import Config.DiffVM_cff as diffvm

process = cepgen.Module('diffvm',
    processParameters = diffvm.defaultProcessParameters.clone(
        #vmFlavour = 23, # Z
        #protonMode = diffvm.BeamMode.StandardFragmentation,
    ),
    inKinematics = cepgen.Parameters(
        pdgIds = (11, 2212),
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
