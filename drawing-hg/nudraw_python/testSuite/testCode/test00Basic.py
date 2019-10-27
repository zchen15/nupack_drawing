# Nudraw tests
# Conrad Steenberg <conrad.steenberg@caltech.edu>
# Nov 20, 2008

import unittest
import random
import sys
import os
import Utils

plist=os.path.dirname(__file__).split(os.sep)[:-2]
nudraw_path=Utils.join_path(plist)
nudraw_exe=Utils.join_path(plist+["nudraw.py"])

data_path=Utils.join_path(plist+["testSuite","testData"])
result_path=Utils.join_path(plist+["testSuite","testResults"])

class TestBasic(unittest.TestCase):
  def __init__(self, methodName="MfeBasic", prefix="844_0", clearPath=True):
    self.prefix=prefix
    self.methodName=methodName
    self.clearPath=clearPath
    unittest.TestCase.__init__(self, methodName)
    self.data_path=Utils.join_path([data_path, 'basic'])

  #-------------------------------------------------------------------------------
  def setUp(self):
    if self.methodName.find('MfeBasic')>=0:
      self.test_path=Utils.join_path([result_path, 'basic', self.methodName])
      if self.clearPath:
        Utils.clearData(self.test_path)
      self.createPathsBasic('basic', self.prefix)

  #-------------------------------------------------------------------------------
  def createPathsMulti(self, path_suffix, prefix):
    self.test_path=Utils.join_path([result_path, path_suffix, prefix])
    self.mfefile=Utils.join_path([data_path, path_suffix, prefix+'.mfe'])
    self.probfile=Utils.join_path([data_path, path_suffix, prefix+'.ppairs'])
    self.jsonout=Utils.join_path([data_path, path_suffix, prefix+'.json'])
    self.pngfile2d=Utils.join_path([self.test_path, prefix+'_%s.png'])
    self.svgfile2d=Utils.join_path([self.test_path, prefix+'_%s.svg'])

    try:
      os.makedirs(self.test_path)
    except OSError, v:
      if v.errno!=17: raise

  #-------------------------------------------------------------------------------
  def createPathsBasic(self, path_suffix, prefix):
    self.mfefile=Utils.join_path([data_path, path_suffix, prefix+'.mfe'])
    self.probfile=Utils.join_path([data_path, path_suffix, prefix+'.ppairs'])
    self.jsonout=Utils.join_path([data_path, path_suffix, prefix+'.json'])
    self.pngfile2d=Utils.join_path([self.test_path, prefix+'.png'])
    self.svgfile2d=Utils.join_path([self.test_path, prefix+'.svg'])

    try:
      os.makedirs(self.test_path)
    except OSError, v:
      if v.errno!=17: raise

  #-------------------------------------------------------------------------------
  def testHelp(self):
    args=[ "", "--help"]
    return Utils.execTest(nudraw_exe, args)

  #-------------------------------------------------------------------------------
  def MfeBasic(self):
    sys.stderr.write(self.prefix)
    self.createPathsBasic('basic', self.prefix)
    tmpout=Utils.join_path([self.test_path, self.prefix+'.out'])
    tmperr=Utils.join_path([self.test_path, self.prefix+'.err'])
    args=[ "",
           "--colorbar",
           "--drawbasenumbers",
           "--jsonout="+self.jsonout,
           "--pngfile2d="+self.pngfile2d,
           "--mfefile="+self.mfefile,
           "--stdout=%s"%tmpout,
           "--stderr=%s"%tmperr]
    try:
      status=os.stat(self.probfile)
      args+="--probfile="+self.probfile,
    except OSError:
      pass

    result = Utils.execTest(nudraw_exe, args)

    if result:
      sys.stderr.write("(ok)...")
    else:
      sys.stderr.write("(er)...")
    assert(result == True)


  #-------------------------------------------------------------------------------
  def MfeFull(self):
    sys.stderr.write(self.prefix)
    self.createPathsMulti('basic', self.prefix)
    Utils.clearData(self.test_path)
    tmpout=Utils.join_path([self.test_path, self.prefix+'.out'])
    tmperr=Utils.join_path([self.test_path, self.prefix+'.err'])
    args=[ "",
           "--colorbar=1,0,1,0,0,0","--filecounter2d=0,1,2,3,4,5,6",
           "--colorbarspace=1,1,1,1,1,1,1 ", "--colorbaseprob=1,0,1,0,0,0,0",
           "--drawbases=0,1,1,0,0,1,0", "--colorstrands=0,0,0,0,0,0,0",
           "--drawbasenumbers=1,1,1,1,1,1,1", "--colorbaseid=0,0,0,0,1,1,1",
           "--baseidbar=0,0,0,0,1,1,1 ",
           "--drawbaseticks=0,0,0,0,0,0,1 ",
           "--jsonout="+self.jsonout,
           "--pngfile2d="+self.pngfile2d,
           "--svgfile2d="+self.svgfile2d,
           "--probfile="+self.probfile,
           "--mfefile="+self.mfefile,
           "--stdout=%s"%tmpout,
           "--stderr=%s"%tmperr]

    print reduce(lambda x,y:x+" "+y, args)
    result = Utils.execTest(nudraw_exe, args)
    if result:
      sys.stderr.write("(ok)...")
    else:
      sys.stderr.write("(er)...")
    assert(result == True)

def suiteInitBasic(prefixes=['832_0', '844_0']):
    test_names = ['MfeBasic']
    clearPath=True
    tests=[]
    for p in prefixes:
      for n in test_names:
        tests.append(TestBasic(n, p, clearPath=clearPath))
        clearPath=False
    return tests

def suiteInitMulti(prefixes=['997_0']):
    test_names = ['MfeFull']
    tests=[]
    for p in prefixes:
      tests+=map(lambda n: TestBasic(n, p), test_names)
    return tests


suite=unittest.TestSuite()
suite.addTests(suiteInitBasic())
suite.addTests(suiteInitMulti())
