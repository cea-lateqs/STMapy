#Set PyQt4 to API level 2
#Import this module as first thing

import sip
API_NAMES = ["QString", "QVariant"]
API_VERSION = 2
for name in API_NAMES:
    sip.setapi(name, API_VERSION)