#!/usr/bin/env python
"""
This module mimics some of the functionality of
WsExample.sh which can be found at:
http://cdaweb.gsfc.nasa.gov/WebServices/REST/WsExample.sh

This module requires the following external utilities:
1) curl
2) xmllint ??? (This one may not actually be necessary)

Some of the docstrings on the functions here are not
exactly correct as I just borrowed them from the comments
in WsExample.sh,  but some of the API was changed slightly.

In any event,  this interface should be a little less
verbose than the original.

I think the only functions that I will actually use are:
1) getTextData
2) downloadResults

@TODO:  Create the xml files in a more sane manner than
        as simple text strings (Yuck!)
@TODO:  See if I can replace the subprocess calls to curl
        with urllib/urllib2
"""
from __future__ import annotations

import os
import shlex
from subprocess import PIPE, Popen, call

try:
    from urllib import urlretrieve
except ImportError:
    from urllib.request import urlretrieve
import datetime
import xml.sax


###################################################################
####                    Parsing helpers                        ####
###################################################################
class FileNameHandler(xml.sax.handler.ContentHandler, list):
    """
    Parse the "filenames" (actually URLS) from the returned
    results file
    """
    def __init__(self):
        list.__init__(self)
        self.isName = False
    def startElement(self, name, attr):
        if name == "Name":
            self.isName = True
    def endElement(self, name):
        if name == "Name":
            self.isName = False
    def characters(self, data):
        if self.isName:
            self.append(data)


if call(['which', 'xmllint'], stderr=PIPE, stdout=PIPE) == 0:
    _has_xmllint = True
else:
    _has_xmllint = False

_mypid = str(os.getpid())
userAgent = "curlWsExample"
traceOption = ''
Debug = False

class DownloadError(Exception):
    pass

def downloadResults(results, directory=os.curdir, prefix='', keepgoing=False):
    """
    Given a list of URLs,  download each.
    keyword args:
       keepgoing: If True,  do not raise a DownloadError if a file does not download
       directory: Where to store the data.
       prefix:    Like tempfile's prefix argument.
    """
    if not hasattr(results, '__iter__'):
        results = (results, )
    if not os.path.exists(directory):
        os.makedirs(directory)

    out = []
    for r in results:
        out1 = os.path.join(directory, prefix+os.path.basename(r))
        try:
            lfname = urlretrieve(r, out1)
            out.append(out1)
        except Exception:
            if keepgoing:
                pass
            else:
                raise DownloadError(f"Could not download file: {r}")
    return out


def datetimeToCDAWebTimeString(obj):
    """
    Convert a datetime object to a CDAWebTimeString
    """
    try:
        out = obj.strftime("%Y-%m-%dT%H:%M:%S.")
        if hasattr(obj, 'microseconds'):
            out += str(obj.microseconds)[:3]+"Z"
        else:
            out += '000Z'
        return out
    except Exception:
        return obj #Lets hope it is a string...



def getResultFilenames(resultFile):
    """
    Extracts the result filenames (//FileDescription/Name/text()) from a
    DataResult file and returns a list of the results

    Arguments:
      resultFile name of file containing a DataResult XML document
    """
    parser = xml.sax.make_parser()
    handler = FileNameHandler()
    parser.setContentHandler(handler)
    parser.parse(resultFile)
    return handler

def printResultFilenames(resultFile):
    for l in getResultFilenames(resultFile):
        print(l)

def appendDatasetRequest(fileToAppend, info):
    """
    Creates a DatasetRequest element from the given DatasetId and
    VariableName values and appends it to the given (DataRequest XML)
    file.

    Arguments:
       fileToAppend: name of file to append to
       info: DatasetId and VariableName values separated by a '/'
       """
    linfo = info.split('/')
    datasetID = linfo[0]
    variableName = linfo[1]
    f = open(fileToAppend, 'a')
    #I wonder if multiple variables can be requested here
    # by adding multiple <VariableName>...</VariableName> tags
    f.write("""
    <DatasetRequest>
      <DatasetId>%s</DatasetId>
      <VariableName>%s</VariableName>
    </DatasetRequest>
"""[1:]%(datasetID, variableName))
    f.close()

###################################################################
## beautify only works if we have xmllint ... however,  I        ##
## do not think that it's use in this script is strictly         ##
## necessary ...                                                 ##
###################################################################
if _has_xmllint:
    def beautify(stream, stdout=PIPE):
        proc = Popen(['xmllint', '--format', '-'], stdin=stream,
                     stdout=stdout)
        if stdout is PIPE:
            return proc.stdout
        return None
else:
    def beautify(stream, stdout=None):
        return stream
###################################################################

def getWadl(endpointURL):
    """
    Makes an HTTP OPTIONS request and echos the results to standard-out.

    Arguments:
       endpointURL: endpoint URL
    """
    endpointURL = endpointURL+'/'
    f = open('cdas.wadl', 'w')
    proc1 = Popen(shlex.split('curl --user-agent %s --silent --request OPTIONS "%s"'
                              ''%(userAgent, endpointURL)), stdout=PIPE)
    f.write(beautify(proc1.stdout).read())
    f.close()
    return open('cdas.wadl')


def getDataviews(endpointURL):
    """
    Gets (HTTP GET) the dataviews and echos the results to standard-out.

    Arguments:
       endpointURL: endpoint URL
    """
    endpointURL = "%s/dataviews"%endpointURL
    return beautify(Popen(shlex.split('curl --user-agent %s --silent "%s"'
                    ''%(userAgent, endpointURL)), stdout=PIPE).stdout)

def getObservatoryGroups(endpointURL, dataview, instrumentTypes=None):
    """
    # Gets (HTTP GET) the specified dataview's observatory groups.
    # The result is displayed to standard-out.
    #
    # Arguments:
    #   $1 endpoint URL
    #   $2 dataview
    #   $3 instrumentTypes optional list of instrument types
    """
    endpointURL = "%s/dataviews/%s/observatoryGroups?"%(endpointURL, dataview)
    if instrumentTypes is not None:
        for iType in instrumentTypes:
            endpointURL += 'instruementType=%s&'%iType
    #The next 2 lines replace the sed command
    if endpointURL.endswith('?'):
        endpontURL = endpointURL[:-1]
    if endpointURL.endswith('&'):
        endpontURL = endpointURL[:-1]

    return beautify(Popen(shlex.split('curl --user-agent %s --globoff --silent '
                                      '%s "%s"'%(userAgent, traceOption, endpointURL)),
                    stdout=PIPE).stdout)


def HTTPGET(endpointURL, dataview, key):
    endpointURL = "%s/dataviews/%s/%s"%(endpointURL, dataview, key)
    return beautify(
                    Popen(shlex.split('curl --user-agent %s --globoff --silent '
                                      '%s "%s"'%(userAgent, traceOption, endpointURL)),
                          stdout=PIPE).stdout)

def getInstrumentTypes(endpointURL, dataview):
    """
    # Gets (HTTP GET) the specified dataview's instrument types.
    # The result is displayed to standard-out.
    #
    # Arguments:
    #   $1 endpoint URL
    #   $2 dataview
    """
    return HTTPGET(endpointURL, dataview, 'instrumentTypes')

def getInstruments(endpointURL, dataview):
    """
    # Gets (HTTP GET) the specified dataview's instruments     .
    # The result is displayed to standard-out.
    #
    # Arguments:
    #   $1 endpoint URL
    #   $2 dataview
    """
    return HTTPGET(endpointURL, dataview, 'instruments')

def getObservatories(endpointURL, dataview):
    """
    # Gets (HTTP GET) the specified dataview's instruments.
    # The result is returned as a file object
    #
    # Arguments:
    #   $1 endpoint URL
    #   $2 dataview
    """
    return HTTPGET(endpointURL, dataview, 'observatories')

def getDatasets(endpointURL, dataview):
    """
    # Gets (HTTP GET) the specified dataview's instruments.
    # The result is returned as a file object
    #
    # Arguments:
    #   $1 endpoint URL
    #   $2 dataview
    """
    return HTTPGET(endpointURL, dataview, 'datasets')

lastInventoryURL = None
lastInventoryLM = None

def getInventory(endpointURL, dataview, dataset):
    """
    # Gets (HTTP GET) the specified dataset's inventory.
    # The result is displayed to standard-out.
    #
    # Arguments:
    #   $1 endpoint URL
    #   $2 dataview
    #   $3 dataset
    """
    global lastInventoryURL
    global lastInventoryLM

    endpointURL = "%s/dataviews/%s/datasets/%s/inventory"%(endpointURL, dataview, dataset)
    resultFile = '/tmp/getInventory%s.xml'%(_mypid)
    if endpointURL == lastInventoryURL:
        proc = Popen(shlex.split('curl --user-agent %s %s --header '
                                 '"If-Modified-Since: %s" --silent --include '
                                 '"%s"'%(userAgent, traceOption, lastInventoryLM,
                                         endpointURL)),
                     stdout=PIPE)
        proc.wait()
        return proc.stdout
    args = shlex.split('curl --user-agent %s %s --globoff --silent '
                       '--include --output %s "%s"'
                       ''%(userAgent, traceOption, resultFile, endpointURL))
    proc = Popen(args)
    proc.wait()
    lastInventoryURL = endpointURL
    s = "awk '/Last-Modified/{printf \"%s\\n\", substr($0,  16)}' " + resultFile
    proc = Popen(shlex.split(s), stdout=PIPE)
    lastInventoryLM = proc.stdout.read()


    out = beautify(Popen(shlex.split('tail -n +8 %s'%resultFile), stdout=PIPE).stdout)

    os.remove(resultFile)
    return out


def getVariables(endpointURL, dataview, dataset):
    """
    # Gets (HTTP GET) the specified dataset's variables.
    # The result is displayed to standard-out.
    #
    # Arguments:
    #   $1 endpoint URL
    #   $2 dataview
    #   $3 dataset
    """
    endpointURL = ("%s/dataviews/%s/datasets/%s/variables"
                   ""%(endpointURL, dataview, dataset))
    return beautify(Popen(shlex.split('curl --user-agent %s %s --globoff '
                                      '--silent "%s"'
                                      ''%(userAgent, traceOption, endpointURL)),
                    stdout=PIPE).stdout)


#I have a feeling that I am going to change the API on this one...
def doPost(endpointURL, postDataFile, outFile):
    """
    # Performs an HTTP POST with the specified data.
    #
    # Arguments:
    #   $1 endpoint URL
    #   $2 name of file containing POST data
    #   $3 name of file to store results in
    """
    if Debug:
        print(beautify(open(postDataFile)).read())
    s = ("curl --user-agent %s %s "%(userAgent, traceOption) +
         "--globoff --silent " +
         r'--header "Content-Type: application/xml" ' +
         r'--header "Accept: application/xml" ' +
         r'--data-ascii "@%s" '%postDataFile +
         '--output %s '%outFile+
         '"%s"'%endpointURL)
    proc = Popen(shlex.split(s))
    proc.wait()

def getTextData(endpointURL, dataview, starttime, endtime,
                datasetID, compression=('Uncompressed', )):
    """
    Returns a list of filenames
    """
    starttime = datetimeToCDAWebTimeString(starttime)
    endtime = datetimeToCDAWebTimeString(endtime)
    postData = "/tmp/getDataPost%s.xml"%_mypid
    f = open(postData, 'w')
    f.write("""
<?xml version="1.0" encoding="UTF-8" standalone="yes"?>
<DataRequest xmlns="http://cdaweb.gsfc.nasa.gov/schema">
  <TextRequest>
    <TimeInterval>
      <Start>XXX</Start>
      <End>YYY</End>
    </TimeInterval>
""".replace('XXX', starttime).replace('YYY', endtime)[1:])
    f.close()
    # if not hasattr(datasetID, '__iter__'):  # this broke in py3k
    if not isinstance(datasetID, (list, tuple)):
        datasetID = (datasetID, )
    for d in datasetID:
        appendDatasetRequest(postData, d)

    # if not hasattr(compression, '__iter__'):  # this broke in py3k
    if not isinstance(compression, (list, tuple)):
        compression = (compression, )

    s = "".join("\n    <Compression>%s</Compression>"%comp for comp in compression)

    f = open(postData, 'a')
    f.write((s+"""
  </TextRequest>
</DataRequest>
""")[1:])
    f.close()
    resultFile = "/tmp/getTextData%s.xml"%_mypid
    endpointURL = "%s/dataviews/%s/datasets"%(endpointURL, dataview)
    doPost(endpointURL, postData, resultFile)
    out = getResultFilenames(resultFile)
    os.remove(postData)
    os.remove(resultFile)
    return out


if __name__ == "__main__":
    endpoint = 'http://cdaweb.gsfc.nasa.gov/WS/cdasr/1'
    dataview = 'sp_phys'
    start = "2005-01-01T00:00:00.000Z"
    end   = "2005-01-02T00:00:00.000Z"  # pylint: disable=bad-whitespace

    start = "2000-09-16T00:00:00.000Z"
    end   = "2000-09-16T02:00:00.000Z"  # pylint: disable=bad-whitespace
    start = datetime.datetime(year=2000, month=9, day=16)
    end = datetime.datetime(year=2000, month=9, day=16, hour=2)
    datasetVar = "AC_H1_MFI/Magnitude"
    dataset = datasetVar.split('/')[0]
    print(getWadl(endpoint).read())
    print('='*80)
    print(getDataviews(endpoint).read())
    print('='*80)
    instrumentTypes = "Magnetic%20Fields%20(space)"
    print(getObservatoryGroups(endpoint, dataview, instrumentTypes).read())
    print('='*80)
    print(getInstrumentTypes(endpoint, dataview).read())
    print('='*80)
    print(getInstruments(endpoint, dataview).read())
    print('='*80)
    print(getObservatories(endpoint, dataview).read())
    print('='*80)
    print(getDatasets(endpoint, dataview).read())
    print('='*80)
    print(getInventory(endpoint, dataview, dataset).read())
    print('='*80)
    print(getInventory(endpoint, dataview, dataset).read())
    print('='*80)
    print(getVariables(endpoint, dataview, dataset).read())
    print('='*80)
    print(getTextData(endpoint, dataview, start, end, datasetVar))
    datasetVars = "AC_H0_MFI/BGSEc"
    print(downloadResults(getTextData(endpoint, dataview, start, end,
                                      datasetVars)))
