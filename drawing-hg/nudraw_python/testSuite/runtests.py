#!/usr/bin/env python
# Run test suite
# 
# Conrad Steenberg <conrad.steenberg@caltech.edu>
# Nov 25, 2008

import sys
import os
import unittest
import testCode
import re
import pdb

#-------------------------------------------------------------------------------
# Options processing
def join_path(plist):
  return reduce(lambda x,y:x+os.sep+y, plist)

plist=os.path.dirname(testCode.__file__).split(os.sep)[:-2]
nudraw_path=testCode.Utils.join_path(plist)
sys.path.append(nudraw_path)
import Options
import testOptions

option_vars=Options.parse(__file__, sys.argv,
  testOptions.def_option_vars, testOptions.option_vars_doc)
cmd_options=Options.ComplexOptions(option_vars)

testexp=re.compile('test\d{2}.*')

#-------------------------------------------------------------------------------
def getSuites():
  codeModules=filter(testexp.match, dir(testCode))
  codeModules.sort()
  return codeModules  

#-------------------------------------------------------------------------------
def listSuites(codeModules):
  map(lambda x: sys.stdout.write(x+"\n"),codeModules)

#-------------------------------------------------------------------------------
# Run test suites
if __name__=='__main__':

  # Determine suites available
  testSuites=getSuites()

  # List suites
  if cmd_options.list:
    listSuites(testSuites)
    sys.exit(0)

  # Use user supplied or complete list of suites
  if not cmd_options.suites:
    runSuites=testSuites
  else:
    runSuites=cmd_options.suites

  # Warn user if a specified suite is not available
  for mod_info in runSuites:
    mod_name = mod_info.split(":")[0]
    mod=getattr(testCode,mod_name)
    if not mod:
      print "Error: could not find suite named %s"%mod_name
      sys.exit(0)

  # Run suites
  for mod_spec in runSuites:
    mod_info = mod_spec.split(":")
    mod_name = mod_info[0]
    mod=getattr(testCode, mod_name)
    if len(mod_info)<2:
      suite=getattr(mod, "suite")
    else:
      sections = mod_info[1].split("/")
      suite = mod.suiteInit(None, sections)
    if suite and not cmd_options.dryrun:
      unittest.TextTestRunner(verbosity=2).run(suite)


