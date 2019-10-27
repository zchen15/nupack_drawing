# Parse Configuration options
# Conrad D. Steenberg <conrad.steenberg@caltech.edu>
# Sep 24, 2008


#~ -freeenergybuffer -drawbases -colorbarbuffer -svgfile 792_0_2_mfe.svg
#~ -freeenergy -energypreamble "Free energy of secondary structure:" -prob
#~ 792_0.ppairs 792_0.mfe &

import sys
import re
import string
import json2 as json

# ------------------------------------------------------------------------------
optionexp=re.compile("^-{1,2}([^ \t\n\r\f\v=]+)\=?(.*)")
class IntList(list):
  pass

  def __getattribute__(self, attr):
    print "IntList: getting attribute %s"%attr
    return super(list, self).__getattribute__(attr)

  def __getitem__(self, i):
    if 0<=i<len(self):
      return list.__getitem__(self,i)
    else:
      return False

class StringList(list):
  def __getattribute__(self, attr):
    print "StringList: getting attribute %s"%attr
    return super(list, self).__getattribute__(attr)

  def __getitem__(self, i):
    if 0<=i<len(self):
      return list.__getitem__(self,i)
    else:
      return False


globaloptions={}

def get_globaloption(containerclass, item):
  global globaloptions
  gid="%x"%(id(containerclass)&0xffff)
  return globaloptions[gid][item]

def get_globaloptions(containerclass):
  global globaloptions
  gid="%x"%(id(containerclass)&0xffff)
  return globaloptions[gid]


def set_globaloption(containerclass, vars):
  global globaloptions
  gid="%x"%(id(containerclass)&0xffff)
  if not globaloptions.has_key(gid):
    globaloptions[gid]={}
  globaloptions[gid].update(vars)

def set_globaloption_item(containerclass, key, val):
  global globaloptions
  gid="%x"%(id(containerclass)&0xffff)
  globaloptions[gid]=vars


class ComplexOptions(dict):
  def __init__(self, default_vars, changed_vars={}):
    set_globaloption(self, default_vars)
    set_globaloption(self, changed_vars)

  def __getattribute__(self, attr):
    return get_globaloption(self, attr)

  def __setattribute__(self, attr, val):
    return set_globaloption_item(self, attr, val)

  def __getitem__(self, attr):
    val= get_globaloption(self, attr)
    return val

  def __str__(self):
    retval=""
    for key, val in get_globaloptions(self).items():
      retval+="%s: %s\n"%(key,val)
    return retval


def subst(name, elem):
  try:
    return name%elem
  except TypeError:
    return name


# ------------------------------------------------------------------------------
def print_usage(execfile, vars, vars_doc):
  print("Usage:\n python %s --<option>=<value>\n"%execfile+\
  "   or \n python %s --<switch>\n"%execfile)
  print("Switches:")
  keys=vars.keys()
  keys.sort()

  maxlen=0
  for key in keys:
    if len(key)>maxlen: maxlen=len(key)
  for key in keys:
    if vars_doc.has_key(key):
      if type(vars[key])==str:
        continue
      elif type(vars[key])==int:
        if vars[key]:
          value="True"
        else:
          value="False"
      elif type(vars[key])==IntList:
        value=[]
        for opt in vars[key]:
          if opt>0:
            value.append("True")
          else:
            value.append("False")
        value=str(value)
      else:
        continue

      print("%s : %s\n%s %s"%(string.rjust(key,maxlen),value,
            " "*maxlen,vars_doc[key]))

  print("\nOptions:")
  for key in keys:
    if vars_doc.has_key(key):
      if type(vars[key])==str:
        value=vars[key]
      elif type(vars[key])==StringList:
        value=str(vars[key])
      else:
        continue
      print("%s : %s\n%s %s"%(string.rjust(key,maxlen),value,
            " "*maxlen,vars_doc[key]))

# ------------------------------------------------------------------------------
def parse(execfile, argv, vars, vars_doc):
  global optionexp
  printhelp=False
  for i in range(len(argv)):
    #~ print argv[i]
    res_o=optionexp.match(argv[i])
    if not res_o: continue
    res=res_o.groups()
    try:
      if res[1]:
        if vars.has_key(res[0]) and type(vars[res[0]])!=int:
          if type(vars[res[0]])==str:
            vars[res[0]]=res[1]
          if type(vars[res[0]])==StringList:
            vars[res[0]]=StringList(res[1].split(","))
          if type(vars[res[0]])==IntList:
            vars[res[0]]=IntList(map(lambda n: n>0,map(int,res[1].split(","))))
            #~ print "%s=Intlist%s"%(res[0], vars[res[0]])
        else:
          raise ValueError("*** Error: No such option '%s'"%res[0])
      else:
        if string.lower(res[0])=='help' or string.lower(res[0])=='h':
          #~ print "printhelp!"
          printhelp=True
          continue
        elif vars.has_key(res[0]):
          if type(vars[res[0]])==type(1):
            vars[res[0]]=True
          elif type(vars[res[0]])==IntList:
            vars[res[0]]=IntList([True])
          else: raise ValueError("*** Error: Option '%s' takes an argument"%res[0])

        else:
          if string.find(res[0],'enable-')==0 and \
            vars.has_key(res[0][7:]):
              if type(vars[res[0]])==type(1):
                vars[res[0][7:]]=True
              elif type(vars[res[0]])==IntList:
                vars[res[0][7:]]=IntList([True])
          elif string.find(res[0],'disable-')==0 and \
            vars.has_key(res[0][8:]):
              if type(vars[res[0]])==type(1):
                vars[res[0][7:]]=False
              elif type(vars[res[0]])==IntList:
                vars[res[0][7:]]=IntList([False])
          else: raise ValueError("*** Error: No such switch '%s'"%res[0])
    except ValueError,v:
      print v.args[0]
      sys.exit(0)
    except:
      sys.exit(0)
  #~ print "printhelp=%s"%printhelp
  if printhelp:
    print_usage(execfile, vars, vars_doc)
    sys.exit(0)

  return vars

# ------------------------------------------------------------------------------
def parse_json(execfile, argv, vars, vars_doc, stream):
  global optionexp
  printhelp=False

  # Read JSON lines until an empty line is found
  json_string=""
  while 1:
    line=stream.readline()
    sys.stderr.write(line)
    if len(line)==0 or line[0]=='\n': break
    json_string+=line

  if len(line)==0 and len(json_string)==0:
    return None
  #~ try:
  json_dict=json.read(json_string)
  #~ except json.ReadException, v:
    #~ return {}

  for res in json_dict["options"].items():
    try:
      if res[1]:
        if vars.has_key(res[0]) and type(vars[res[0]])!=int:
          if type(vars[res[0]])==str:
            vars[res[0]]=res[1]
          if type(vars[res[0]])==StringList:
            vars[res[0]]=StringList(res[1].split(","))
          if type(vars[res[0]])==IntList:
            vars[res[0]]=IntList(map(lambda n: n>0,map(int,res[1].split(","))))
            #~ print "%s=Intlist%s"%(res[0], vars[res[0]])
        else:
          raise ValueError("*** Error: No such option '%s'"%res[0])
      else:
        if string.lower(res[0])=='help' or string.lower(res[0])=='h':
          #~ print "printhelp!"
          printhelp=True
          continue
        elif vars.has_key(res[0]):
          if type(vars[res[0]])==type(1):
            vars[res[0]]=True
          elif type(vars[res[0]])==IntList:
            vars[res[0]]=IntList([True])
          else: raise ValueError("*** Error: Option '%s' takes an argument"%res[0])

        else:
          if string.find(res[0],'enable-')==0 and \
            vars.has_key(res[0][7:]):
              if type(vars[res[0]])==type(1):
                vars[res[0][7:]]=True
              elif type(vars[res[0]])==IntList:
                vars[res[0][7:]]=IntList([True])
          elif string.find(res[0],'disable-')==0 and \
            vars.has_key(res[0][8:]):
              if type(vars[res[0]])==type(1):
                vars[res[0][7:]]=False
              elif type(vars[res[0]])==IntList:
                vars[res[0][7:]]=IntList([False])
          else: raise ValueError("*** Error: No such switch '%s'"%res[0])
    except ValueError,v:
      print v.args[0]
      sys.exit(0)
    except:
      raise
      sys.exit(0)
  #~ print "printhelp=%s"%printhelp
  if printhelp:
    print_usage(execfile, vars, vars_doc)
    sys.exit(0)

  json_dict["options"]=vars
  return json_dict


if __name__=="__main__":
  def_option_vars={\
        "material"        : "rna",
        "svgfile2d"       : "",
        "pngfile2d"       : "",
        "filecounter2d"   : StringList(),
        "dotfile"         : "",
        "mfefile"         : "",
        "probfile"        : "",
        "jsonout"         : "",
        "jsonin"          : "",
        "scenefile3d"     : "",
        "show3d"          : 0,
        "show2d"          : 0,
        "drawbases"       : IntList(),
        "colorbases"      : IntList(),
        "colorbar"        : IntList(),
        "colorbarspace"   : IntList(),
        "colorstrands"    : IntList(),
  }

  option_vars_doc={ \
        "material"        : "Either rna or dna",
        "svgfile2d"       : "Name of the 2D svg file to be produced",
        "pngfile2d"       : "Name of the 2D PNG file to be produced",
        "dotfile"         : "Name of the Graphviz file. In this file loops are nodes and helices are edges",
        "mfefile"         : "Input file to read structure information from",
        "probfile"        : "Input file to read probability information from",
        "filecounter2d"   : "Comma separated list of strings to substitute in export filenames",
        "jsonout"         : "Name of json file to be created, skip if not provided",
        "jsonin"          : "Name of json file to be read, skip if not provided",
        "scenefile3d"     : "3D Scene file export for Tachyon raytracer",
        "show3d"          : "Show 3D OpenGL model",
        "show2d"          : "Show 2D layout",
        "drawbases"       : "Draw base letters if available from mfe file",
        "colorbases"      : "Shade bases according to probfile",
        "colorbar"        : "Draw probability colorbar",
        "colorbarspace"   : "Reserve space for colorbar, but draw only if colorbar option set",
        "colorstrands"    : "Draw strands in different colors"
  }


  newvars=parse_options(__file__, sys.argv, def_option_vars, option_vars_doc)
  print_usage(__file__, def_option_vars, option_vars_doc)
  testoptions=ComplexOptions(newvars)
  print "testoptions.material=%s"%testoptions.material
