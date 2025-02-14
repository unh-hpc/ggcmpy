from __future__ import absolute_import, print_function
import xml.sax.handler as handler
import xml.sax

from . import *


class CDA_XML(object):
    """
    This is the base class for all the data
    parsed out of a CDA_XML file,
    provided so something can do "isinstance"
    to check the object.
    """
    pass

def CDA_XML_factory(expectedAttr, NameAttr='_Name'):
    if 'Name' in expectedAttr:
        expectedAttr.remove('Name')
    if NameAttr != "_Name":
        expectedAttr.append(NameAttr)

    class cls(CDA_XML):
        def __init__(self):
            self._TYPE = (self.__class__.__name__).replace("CDA_XML_", "")
            self.attr = ['Name'] + expectedAttr
            for a in expectedAttr:
                setattr(self, a, None)
            self.Name = None

        @property
        def TYPE(self):
            return self._TYPE

        @property
        def Name(self):
            return getattr(self, NameAttr)

        @Name.setter
        def Name(self, value):
            setattr(self, NameAttr, value)

        def __str__(self):
            s = (str(self.Name) + "\n\t" +
                 "\n\t".join([(v + ': ' + str(getattr(self, v))) for v in self.attr[1:]]))
            return s

        def __setattr__(self, attr, val):
            #Attribute doesn't exist
            try:
                n = getattr(self, attr)
            except AttributeError:
                object.__setattr__(self, attr, val)
                return

            #Attribute is None,  reset
            if n is None:
                object.__setattr__(self, attr, val)
            else:
                #Attribute is already a list,  append
                if isinstance(n, list):
                    n.append(val)
                else:
                    #Attribute was set once,  now append to list
                    object.__setattr__(self, attr, [n, val])
            return

    return cls

class CDA_XML_Variable(CDA_XML_factory(["ShortDescription", "LongDescription"])):
    pass
class CDA_XML_Dataview(CDA_XML_factory(["EndpointAddress", "Title", "Subtitle",
                                        "UnderConstruction", "NoticeUrl"], NameAttr="Id")):
    pass
class CDA_XML_ObservatoryGroup(CDA_XML_factory(["ObservatoryId"])):
    pass
class CDA_XML_InstrumentType(CDA_XML_factory([])):
    pass
class CDA_XML_Instrument(CDA_XML_factory(["ShortDescription", "LongDescription"])):
    pass
class CDA_XML_Observatory(CDA_XML_factory(["ShortDescription", "LongDescription"])):
    pass

class CDA_XML_Dataset(CDA_XML_factory(["Observatory", "Instrument", "ObservatoryGroup",
                                       "InstrumentType", "Label", "PiName",
                                       "PiAffiliation", "Notes"], NameAttr='Id')):
    pass

class MyHandler(handler.ContentHandler, dict):
    subclasses = {"Variable":CDA_XML_Variable,
                  "Dataview":CDA_XML_Dataview,
                  "ObservatoryGroup":CDA_XML_ObservatoryGroup,
                  "InstrumentType":CDA_XML_InstrumentType,
                  "Instrument":CDA_XML_Instrument,
                  "Observatory":CDA_XML_Observatory,
                  "Dataset":CDA_XML_Dataset,
                  }
    def __init__(self):
        dict.__init__(self)
        self.majorElement = None
        self.Element = None

    def startElement(self, name, attributes):
        if name.endswith("Description") and name[:-11] in self.subclasses:
            self.majorElement = self.subclasses[name[:-11]]()
            self.Element = None
            assert name[:-11] == self.majorElement.TYPE
        else:
            self.Element = (name, attributes)

    def characters(self, data):
        if str(data) == data:
            data = str(data)
        if self.majorElement and self.Element:
            setattr(self.majorElement, self.Element[0], data)

    def endElement(self, name):
        if self.majorElement and name[:-11] == self.majorElement.TYPE:
            try:
                self[self.majorElement.Name] = self.majorElement
            except Exception:
                import sys
                sys.stderr.write("Invalid major Element:%s\nIGNORING\n"%str(self.majorElement))
            self.majorElement = None
            self.Element = None
        else:
            self.Element = None

    def __str__(self):
        return "\n".join([str(s) for s in sorted(self.values(), key=lambda x: x.Name)])


def pyXMLparse(func, *args, **kwargs):
    """
    factor a little.  This function handles
    the parsing of the XML ...
    """
    parser = xml.sax.make_parser()
    handler = MyHandler()
    parser.setContentHandler(handler)
    parser.parse(func(*args, **kwargs))
    return handler


def pyGetVariables(endpointURL, dataview, dataset):
    return pyXMLparse(getVariables, endpointURL, dataview, dataset)

def pyGetDataviews(endpointURL):
    return pyXMLparse(getDataviews, endpointURL)

def pyGetObservatoryGroups(endpointURL, dataview, **kwargs):
    return pyXMLparse(getObservatoryGroups, endpointURL, dataview, **kwargs)

def pyGetInstrumentTypes(endpointURL, dataview):
    return pyXMLparse(getInstrumentTypes, endpointURL, dataview)

def pyGetInstruments(endpointURL, dataview):
    return pyXMLparse(getInstruments, endpointURL, dataview)

def pyGetObservatories(endpointURL, dataview):
    return pyXMLparse(getObservatories, endpointURL, dataview)

def pyGetDatasets(endpointURL, dataview):
#    return getDatasets(endpointURL, dataview).read()
    return pyXMLparse(getDatasets, endpointURL, dataview)

##########################################################################
##########################################################################
##########################################################################
##########################################################################
if __name__ == "__main__":
    end = 'http://cdaweb.gsfc.nasa.gov/WS/cdasr/1'
    dv = 'sp_phys'
    dset = 'GE_K0_MGF'
    #dset='WI_H0_MFI'
    print(pyGetVariables(end, dv, dset))
    #exit(1) #This is nice for inspecting to try to figure out
    #        #which variables are available for a particular dataset.
    print()
    print(pyGetVariables(end, dv, dset))
    print()
    print(pyGetDataviews(end))
    print()
    print(pyGetObservatoryGroups(end, dv))
    print()
    print(pyGetInstruments(end, dv))
    print()
    print(pyGetObservatories(end, dv))
    print()
    print(pyGetDatasets(end, dv))
    print()
    print(getDatasets(end, dv).read())

#    pyGetDataviews(end)
