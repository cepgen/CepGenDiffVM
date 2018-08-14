import Config.Core as cepgen
#from Config.pythia8_cff import pythia8 as hadroniser
from Config.logger_cfi import logger
from Config.generator_cfi import generator as _gen
from Integrators.vegas_cfi import vegas as integrator


logger.enabledModules += (
    'DiffVM',
)

process = cepgen.Module('diffvm',
    mode = cepgen.ProcessMode.ElasticElastic,
    inKinematics = cepgen.Parameters(
        pz = (6500., 6500.),
        structureFunctions = cepgen.StructureFunctions.LUXlike,
    ),
    outKinematics = cepgen.Parameters(
        pair = 13,
        pt = (25.,),
        energy = (0.,),
        eta = (-2.5, 2.5),
        mx = (1.07, 1000.),
    ),
)

generator = _gen.clone(
    numEvents = 100000,
    printEvery = 10000,
)

