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
  def __init__(self, methodName="MfeBasic", prefix="844_0", clearPath=True, pair=(1,1)):
    self.prefix=prefix
    self.pair=pair
    self.methodName=methodName
    self.clearPath=clearPath
    unittest.TestCase.__init__(self, methodName)
    self.data_path=Utils.join_path([data_path, 'complex'])

  #-------------------------------------------------------------------------------
  def setUp(self):
    if self.methodName.find('MfeBasic')>=0:
      self.test_path=Utils.join_path([result_path, 'complex', self.methodName])
      if self.clearPath:
        Utils.clearData(self.test_path)
      self.createPathsBasic('basic', self.prefix)

  #-------------------------------------------------------------------------------
  def createPathsMulti(self, path_suffix, prefix):
    pair_prefix="%s_%s_%s"%(prefix, self.pair[0], self.pair[1])
    self.test_path=Utils.join_path([result_path, path_suffix, prefix])
    self.mfefile=Utils.join_path([data_path, path_suffix, prefix+'.ocx-mfe'])
    self.probfile=Utils.join_path([data_path, path_suffix, prefix+'.ocx-ppairs'])
    self.keyfile=Utils.join_path([data_path, path_suffix, prefix+'.ocx-key'])
    self.jsonout=Utils.join_path([data_path, path_suffix, pair_prefix+'.json'])
    self.pngfile2d=Utils.join_path([self.test_path, pair_prefix+'_%s.png'])
    self.svgfile2d=Utils.join_path([self.test_path, pair_prefix+'_%s.svg'])

    try:
      os.makedirs(self.test_path)
    except OSError, v:
      if v.errno!=17: raise

  #-------------------------------------------------------------------------------
  def createPathsBasic(self, path_suffix, prefix):
    pair_prefix="%s_%s_%s"%(prefix, self.pair[0], self.pair[1])
    self.mfefile=Utils.join_path([data_path, path_suffix, prefix+'.ocx-mfe'])
    self.probfile=Utils.join_path([data_path, path_suffix, prefix+'.ocx-ppairs'])
    self.keyfile=Utils.join_path([data_path, path_suffix, prefix+'.ocx-key'])
    self.jsonout=Utils.join_path([data_path, path_suffix, pair_prefix+'.json'])
    self.pngfile2d=Utils.join_path([self.test_path, pair_prefix+'.png'])
    self.svgfile2d=Utils.join_path([self.test_path, pair_prefix+'.svg'])

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
    self.createPathsBasic('complex', self.prefix)
    tmpout=Utils.join_path([self.test_path, self.prefix+'.out'])
    tmperr=Utils.join_path([self.test_path, self.prefix+'.err'])
    args=[ "",
           "--colorbar",
           "--drawbasenumbers",
           "--jsonout="+self.jsonout,
           "--pngfile2d="+self.pngfile2d,
           "--mfefile="+self.mfefile,
           "--keyfile="+self.keyfile,
           "--stdout=%s"%tmpout,
           "--stderr=%s"%tmperr,
           "--complex_permutation=%s,%s"%self.pair]
    try:
      status=os.stat(self.probfile)
      args+="--probfile="+self.probfile,
    except OSError:
      pass

    #~ print reduce(lambda x,y:x+' '+y, args)
    result = Utils.execTest(nudraw_exe, args)

    if result:
      sys.stderr.write("(ok)...")
    else:
      sys.stderr.write("(er)...")
    assert(result == True)


  #-------------------------------------------------------------------------------
  def MfeFull(self):
    sys.stderr.write(self.prefix)
    self.createPathsMulti('complex', self.prefix)
    if self.clearPath:
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
           "--keyfile="+self.keyfile,
           "--complex_permutation=%s,%s"%self.pair,
           "--stdout=%s"%tmpout,
           "--stderr=%s"%tmperr]
    #~ print reduce(lambda x,y:x+' '+y, args)
    result = Utils.execTest(nudraw_exe, args)
    if result:
      sys.stderr.write("(ok)...")
    else:
      sys.stderr.write("(er)...")
    assert(result == True)

def suiteInitBasic(sets=[('1387',(1,1), (2,1))]):
  tests=[]
  clearPath=True
  for s in sets:
    for p in s[1:]:
      #~ print "appendieng (%s, %s)"%(s[0],p)
      tests.append(TestBasic('MfeBasic', prefix=s[0],pair=p, clearPath=clearPath))
      clearPath=False
  return tests

def suiteInitMulti(sets=[('1387',(1,1), (2,1))]):
  test_names = ['MfeFull']
  tests=[]
  for s in sets:
    clearPath=True
    for p in s[1:]:
      #~ print "appending (%s, %s)"%(s[0],p)
      tests.append(TestBasic('MfeFull', prefix=s[0],pair=p, clearPath=clearPath))
      clearPath=False
  return tests


suite=None

if __name__ == '__main__':
  def_suite=suite()
  def_suite.run()
else:
  suite=unittest.TestSuite()
  suite.addTests(suiteInitBasic())
  suite.addTests(suiteInitMulti())
