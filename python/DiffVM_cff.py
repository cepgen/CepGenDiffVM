from containers_cfi import Module, Parameters
from ProcessMode_cfi import ProcessMode
from PDG_cfi import PDG

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

defaultProcessParameters = Parameters(
    mode = ProcessMode.ElasticElastic,
    vmFlavour = PDG.Jpsi,
    photonMode = PhotonMode.WWA,
    protonMode = BeamMode.Elastic,
)

