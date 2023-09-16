from .containers_cff import Parameters


class EPAMode:
    WWA        = 1
    T          = 2
    TL         = 3
    TLinPframe = 4

EPA = Parameters(
    mode = EPAMode.WWA,
    yRange = (0., 1.),
    dyRange = (0., 1.),
)

