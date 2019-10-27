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
print "nudraw_exe=",nudraw_exe

data_path=Utils.join_path(plist+["testSuite","testData"])
result_path=Utils.join_path(plist+["testSuite","testResults"])

render_exe=Utils.join_path(plist[:-1]+["RenderGL","JsonRenderer"])

class TestStructures25D(unittest.TestCase):
  def __init__(self, methodName="runTest", configobj=None, clearData=True, material="dna"):
    self.methodName=methodName
    self.material=material
    if methodName!="runTest":
      setattr(self,methodName, self.runTest)
    self.clearData=clearData
    unittest.TestCase.__init__(self, methodName)
    self.configobj=configobj
    #~ print "nudraw_path =",nudraw_path

    self.data_path=Utils.join_path([data_path, 'structure2.5D'])
    #~ print "data_path =", self.data_path
    #~ self.data_file=

  #-------------------------------------------------------------------------------
  def setUp(self):
    self.test_path=Utils.join_path([result_path, 'structure3D', self.methodName])
    if self.clearData:
      Utils.clearData(self.test_path)
    self.createPaths('structure3D', self.methodName)

  def runTest(self):
    #~ sys.stderr.write(self.prefix+"...")
    structure="...(((...)))..." # TODO replace string
    prefix="structure"
    options=self.configobj.options(self.methodName)
    options.sort()
    sys.stderr.write(self.configobj.get(self.methodName, 'comment')+" - %s:"%self.material)
    for o in options:
      if o=='comment': continue
      sys.stderr.write('...'+o+".(l.")
      structure=self.configobj.get(self.methodName,o)
      sequence, domains, domainnames = Utils.generateSequence(structure)
      self.json3d=Utils.join_path([self.test_path, prefix+'_%s_%s_3d.json'%(self.material,o)])
      self.json25d=Utils.join_path([self.test_path, prefix+'_%s_%s_25d.json'%(self.material,o)])
      self.png3d=Utils.join_path([self.test_path, prefix+'_%s_%s_3d.png'%(self.material,o)])

      tmpout=Utils.join_path([self.test_path, prefix+'_%s_%s.out'%(self.material,o)])
      tmperr=Utils.join_path([self.test_path, prefix+'_%s_%s.err'%(self.material,o)])
      args=[ "",
             "--json3d="+self.json3d,
             "--json25d="+self.json25d,
             "--structure=%s"%structure,
             "--sequence=%s"%sequence,
             "--domainnames=a,a*,b,b*,c,c*",
             "--domains=%s"%domains,
             "--colordomains",
             "--material=%s"%self.material,
             "--stdout=%s"%tmpout,
             "--stderr=%s"%tmperr]

      #~ print reduce(lambda x,y:x+" "+y, args)
      result = Utils.execTest(nudraw_exe, args)

      if result:
        sys.stderr.write("ok)")
      else:
        sys.stderr.write("err)")
      assert(result == True)
      continue
      # Now Render PNG images
      sys.stderr.write(".(r.")
      args=[ "",
         "--jsonfile="+self.json3d,
         "--pngfile=%s"%self.png3d,
         "--pngaa=2",
         "--width=%s"%1024,
         "--height=%s"%1024]
      result = Utils.execTest(render_exe, args)
      if result:
        sys.stderr.write("ok)")
      else:
        sys.stderr.write("err)")
      assert(result == True)
      break

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
  clearData=True
  tests=[]
  for s in sections:
    tests.append(TestStructures25D(s, configobj, clearData=True, material="dna"))
    clearData=False
    tests.append(TestStructures25D(s, configobj, clearData=False, material="rna"))
  return unittest.TestSuite(tests)

# Read config file
suite=None
data_file=Utils.join_path([data_path, 'structure','teststrands.ini'])
config = ConfigParser.SafeConfigParser()
config.read([data_file])
sections=config.sections()
sections.sort()

if __name__ == '__main__':
  pass
  #~ def_suite=suite()
  #~ def_suite.run()
else:
  suite=suiteInit(config, sections)
  #~ print suite
