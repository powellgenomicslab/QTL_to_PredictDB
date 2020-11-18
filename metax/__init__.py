__version__ = "0.5.7"
import metax.WeightDBUtilities as WeightDBUtilities
import metax.PrediXcanFormatUtilities as PrediXcanFormatUtilities
import metax.ThousandGenomesUtilities as ThousandGenomesUtilities
import metax.Logging as Logging
import metax.Utilities as Utilities
import metax.Formats as Formats

def exitIf(doExit, Exception, msg):
    if doExit:
        raise Exception(msg)
