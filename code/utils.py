import os
import re
import xml.etree.ElementTree as ET

# helper function to load physicell simulation parameters from config.xml
def loadPhysiCellSettings(dirpath):
    tree = ET.parse(dirpath+"/config.xml")
    root = tree.getroot()
    return(root)

# fetches a single parameter value from the physicell settings ET root
def loadParamValue(root, paramname):
    val = None
    if ('.' in paramname):
        k = paramname.split('.')
        uep = root
        for idx in range(len(k)):
            uep = uep.find('.//' + k[idx])  # unique entry point (uep) into xml
        val = uep.text
    return(val)