from Config.containers_cff import Parameters
from Config.ProcessMode_cff import ProcessMode


class BeamMode:
    '''Beam particles treatment mode'''
    GluonFragmentation    = -1
    Elastic               =  0
    StandardFragmentation =  1
    NucleonPionsDecay     =  2

class PhotonMode:
    '''Photon generation mode'''
    Fixed    = -1
    InvK     =  0
    WWA      =  1
    ABTSmith =  2
    AandS    =  3

PDG = Parameters(
    omega782  = Parameters(pdgid = 223),
    Jpsi      = Parameters(pdgid = 443),
    h1380_1   = Parameters(pdgid = 10333),
    phi1680   = Parameters(pdgid = 100333),
    Upsilon1S = Parameters(pdgid = 553),
    Upsilon2S = Parameters(pdgid = 100553),
    Upsilon3S = Parameters(pdgid = 200553),
)

defaultProcessParameters = Parameters(
    mode = ProcessMode.ElasticElastic,
    vmFlavour = PDG.Jpsi,
    photonMode = PhotonMode.WWA,
    protonMode = BeamMode.Elastic,
)

