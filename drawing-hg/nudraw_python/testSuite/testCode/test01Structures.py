# Nudraw tests
# Conrad Steenberg <conrad.steenberg@caltech.edu>
# Nov 20, 2008

import unittest
import random
import sys
import os
import ConfigParser
import Utils


plist=os.path.dirname(__file__).split(os.sep)[:-2]
nudraw_path=Utils.join_path(plist)
nudraw_exe=Utils.join_path(plist+["nudraw.py"])

data_path=Utils.join_path(plist+["testSuite","testData"])
result_path=Utils.join_path(plist+["testSuite","testResults"])

class TestStructures(unittest.TestCase):
  def __init__(self, methodName="runTest", configobj=None):
    self.methodName=methodName
    if methodName!="runTest":
      setattr(self,methodName, self.runTest)
    unittest.TestCase.__init__(self, methodName)
    self.configobj=configobj
    #~ print "nudraw_path =",nudraw_path

    self.data_path=Utils.join_path([data_path, 'structure'])
    #~ print "data_path =", self.data_path
    #~ self.data_file=

  #-------------------------------------------------------------------------------
  def setUp(self):
    self.test_path=Utils.join_path([result_path, 'structure', self.methodName])
    Utils.clearData(self.test_path)
    self.createPaths('structure', self.methodName)

  def runTest(self):
    #~ sys.stderr.write(self.prefix+"...")
    structure="...(((...)))..." # TODO replace string
    prefix="structure"
    options=self.configobj.options(self.methodName)
    options.sort()
    sys.stderr.write(self.configobj.get(self.methodName, 'comment')+'...')
    for o in options:
      if o=='comment': continue
      sys.stderr.write(o)
      structure=self.configobj.get(self.methodName,o)
      sequence, domains, domainnames = Utils.generateSequence(structure)
      self.jsonout=Utils.join_path([self.test_path, prefix+'_%s_%%s.json'%o])
      self.pngfile2d=Utils.join_path([self.test_path, prefix+'_%s_%%s.png'%o])
      self.svgfile2d=Utils.join_path([self.test_path, prefix+'_%s_%%s.svg'%o])
      tmpout=Utils.join_path([self.test_path, prefix+'_%s_%%s.out'%o])
      tmperr=Utils.join_path([self.test_path, prefix+'_%s_%%s.err'%o])
      args=[ "",
             "--colorbar=1,0,1,0,0,0","--filecounter2d=0,1,2,3,4,5,6",
             "--colorbarspace=1,1,1,1,1,1,1 ", "--colorbaseprob=1,0,1,0,0,0,0",
             "--drawbases=0,1,1,0,0,1,0", "--colorstrands=0,0,0,0,1,0,0",
             "--drawbasenumbers=1,1,1,1,1,0,1", "--colorbaseid=0,0,0,0,1,1,1",
             "--baseidbar=0,0,0,0,1,1,1 ",
             "--drawbaseticks=0,0,0,0,0,0,1 ",
             "--colordomains=0,0,0,0,0,1,0 ",
             "--labeldomains=0,0,0,0,0,1,0 ",
             "--jsonout="+self.jsonout,
             "--pngfile2d="+self.pngfile2d,
             "--svgfile2d="+self.svgfile2d,
             "--structure=%s"%structure,
             "--sequence=%s"%sequence,
             "--domains=%s"%domains,
             "--domaincolors=ff0000,ff9900,00ff00,00ff99,0000ff,9900ff",
             "--domainnames=a,a*,b,b*,c,c*",
             "--stdout=%s"%tmpout,
             "--stderr=%s"%tmperr]

      #~ print reduce(lambda x,y:x+" "+y, args)
      result = Utils.execTest(nudraw_exe, args)
      if result:
        sys.stderr.write("(ok)...")
      else:
        sys.stderr.write("(er)...")
      assert(result == True)

  #-------------------------------------------------------------------------------
  def createPaths(self, path_suffix, prefix):

    try:
      os.makedirs(self.test_path)
    except OSError, v:
      if v.errno!=17: raise


def suiteInit(configobj, sections=[]):
    global config
    if configobj == None:
      configobj = config
    tests=[]

    tests+=map(lambda x: TestStructures(x, configobj), sections)
    return unittest.TestSuite(tests)

# Read config file
suite=None
data_file=Utils.join_path([data_path, 'structure','teststrands.ini'])
config = ConfigParser.SafeConfigParser()
config.read([data_file])
sections=config.sections()
sections.sort()
#~ print sections



if __name__ == '__main__':
  pass
  #~ def_suite=suite()
  #~ def_suite.run()
else:
  suite=suiteInit(config, sections)
  #~ print suite
