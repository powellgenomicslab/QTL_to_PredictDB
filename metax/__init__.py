__version__ = "0.5.7"
from metax.WeightDBUtilities import WeightDBUtilities
from metax.PrediXcanFormatUtilities import PrediXcanFormatUtilities
from metax.ThousandGenomesUtilities import ThousandGenomesUtilities
from metax.Logging import Logging
from metax.Utilities import Utilities
from metax.Formats import Formats

def exitIf(doExit, Exception, msg):
    if doExit:
        raise Exception(msg)
