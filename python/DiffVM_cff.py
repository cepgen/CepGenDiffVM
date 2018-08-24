import Config.Core as cepgen

class BeamMode:
    GluonFragmentation    = -1
    Elastic               =  0
    StandardFragmentation =  1
    NucleonPionsDecay     =  2

class PhotonMode:
    Fixed    = -1
    InvK     =  0
    WWA      =  1
    ABTSmith =  2
    AandS    =  3

defaultProcessParameters = cepgen.Parameters(
    mode = cepgen.ProcessMode.ElasticElastic,
    vmFlavour = 443, # J/psi
    #vmFlavour = 23, # Z
    photonMode = PhotonMode.WWA,
    protonMode = BeamMode.Elastic,
)

