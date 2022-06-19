class NotReachable(Exception):
    """Raised when a pass is not reachable on a time window"""


class PropagationError(Exception):
    """Raised when a calculation issue is found"""


class CelestrakError(Exception):
    """ Raised when TLEs can't be downloaded from Celestrak """


class PassAlgorithmError(Exception):
    """ Raised when pass finding algorithm fails """