# Parse mfe file
# Conrad D. Steenberg <conrad.steenberg@caltech.edu>
# Sep 24, 2008

import sys
import re
import json2 as json
import time

# Expression for extracting parameters from comment lines
cexp=re.compile('%\s+(.*):\s*(\S+).*')
param_exp=re.compile('%\s+(\w+)\s+param[ea]ters.*', flags=re.I)
sexp=re.compile('%\s+%+')
nexp=re.compile('(\d+).*')
mfexp=re.compile('([\+\-\d\.eE]+|nan|inf).*')
pairexp=re.compile('(\d+)\s+(\d+).*')
ppairexp=re.compile('(\d+)\s+(\d+)\s+([\deE\.\+-]+|nan|inf).*')
# ppairexp_pseudo=re.compile('(\d+)\s+(\d+)\s+([\deE\.\+-]*)\s+([\deE\.\+-]*)\s+([\deE\.\+-]+|nan|inf)')
ppairexp_pseudo=re.compile('fuck this code')

compid_exp=re.compile("% id sequence")
cparam_exp=re.compile('%\s+param[ae]ters.*:\s*(\w+).*', flags=re.I)
seqid_exp=re.compile("%\s+(\d+)\s*(\w+)")
compname_exp=re.compile("%\s+(complex\d+-order\d+).*")
compname_exp2=re.compile("complex(\d+)-order(\d+)")
temp_exp=re.compile("%\s+T\s+=\s+([\S]+)")

#-------------------------------------------------------------------------------
# try opening a file, sleeping and repeating a few times if it's not found
def delay_open(fname, mode='r'):
  i=0
  delay=0.1
  f=None
  while f==None and i<5:
    try:
      f=open(fname, mode)
    except IOError,e:
      if e.errno == 2:
        i+=1
        time.sleep(delay)
        sys.stderr.write("Sleeping %1.2f\n"%delay)
        delay*=2
      else:
        break
  if f == None:
    sys.stderr.write("Could not open file %s\n"%fname)
  return f

#-------------------------------------------------------------------------------
class Structure:
  def __init__(self, seq, temp, mfe, struc, pairs, complex_id=-1, order=0):
    self.seq=seq
    self.temp=temp
    self.mfe=mfe
    self.struc=struc
    self.pairs=pairs
    self.probs=[]
    self.complex_id=complex_id
    self.order=order

  def set_probs(self, probs):
    self.probs=probs

  def __str__(self):
    s=\
    "seq:   %s\n"%self.seq+\
    "temp:  %f\n"%self.temp+\
    "mfe:   %f\n"%self.mfe+\
    "struc: %s\n"%self.struc+\
    "pairs:"
    for pair in self.pairs:
      s=s+"\n%3d -> %3d"%(pair[0], pair[1])
    s=s+"\nprobs:\n"
    for basep in self.probs:
      s=s+"%s\n"%(basep)
      for i in range(1,len(basep)-1,2):
        s=s+"  %s -> %s\n"%(basep[i],basep[i+1])

    s=s+"\n"
    return s

#-------------------------------------------------------------------------------
class Complex:
  def __init__(self, seqs):
    self.seqs=seqs
    self.nseqs=len(seqs)
    self.names=[]
    self.structs=[]
    self.params=""
    self.complex_perm=[]

  def add_struct(self, name, comb, params, temp, mfe, struc, pairs):
    self.params=params
    seq=""
    self.complex_perm.append(comb[:2])
    #~ print "self.complex_perm=%s"%(self.complex_perm)
    #~ print "seqs=\n%s"%repr(self.seqs)
    #~ print "combs=\n%s"%repr(comb)
    for i in range(2,len(comb)):
      c=comb[i]
      #~ print "i=%d"%i
      seq=seq+self.seqs[c]
      if i<len(comb)-1:
        seq=seq+'+'
    #~ print seq
    self.structs.append(Structure(seq, temp, mfe, struc, pairs))
    self.names.append(name)

  def __str__(self):
    s=""
    for i in range(len(self.structs)):
      s=s+\
      "%s:\n"%self.names[i] + \
      "-"*(len(self.names[i])+1) + "\n"+\
      str(self.structs[i])
    s+="params: %s\n"%(self.params)
    return s
#-------------------------------------------------------------------------------
def pairsFromStruc(struc):
  sp=-1 # stack pointer
  stack=range(len(struc)/2) # preallocate space
  pairs=[]
  base_index = 0
  for i in xrange(len(struc)):
    if struc[i]=='(':
      sp+=1
      stack[sp]=base_index
      base_index += 1
    elif struc[i]==')':
      pairs.append([stack[sp],base_index])
      sp-=1
      base_index += 1
    elif struc[i]=='.':
      base_index += 1

  return pairs

#-------------------------------------------------------------------------------
def combFromSeq(seq):
  seqs=seq.split('+')
  comb={}
  seqOrd=[]
  lst=[]
  i=0
  for s in seqs:
    c=comb.get(s)
    if c==None:
      comb[s]=i
      lst.append(i)
      seqOrd.append(s)
      i+=1
    else:
      lst.append(c)


  return seqOrd,lst

#-------------------------------------------------------------------------------
def parseMfe(fname):
  structs=[]
  f=delay_open(fname)
  lines=f.readlines() # Read the file into memory, should be small enough
  f.close()
  info={}
  for li in xrange(len(lines)):
    #~ print lines[li]
    mo=cexp.match(lines[li])
    if mo:
      g=mo.groups()
      if g and len(g)==2:
        info[g[0]]=g[1]
        #~ print g
        if info.has_key("Sequence") and \
          (info.has_key("Temperature") or info.has_key("Temperature (C)")):
          if not info.has_key("Temperature"):
            info["Temperature"]=info["Temperature (C)"]
          break
    mo=param_exp.match(lines[li])
    if mo:
      #~ sys.stderr.write(lines[li]+'\n')
      g=mo.groups()
      if g and len(g)==1:
        info["Parameters"]=g[0]

  #~ print info
  foundstruct=False
  li=li+1
  while li < len(lines):
    line=lines[li]
    #~ sys.stdout.write(line)

    if not foundstruct: # We're still looking for a comment line
      mo=sexp.match(line)
      if mo:
        foundstruct=True
        li=li+1
        continue
    else:
      mo=nexp.match(line) # Number or bases
      if mo:
        g=mo.groups()[0]
        numbases=int(g)

        #~ print "numbases=%d"%numbases

      li=li+1
      if li>=len(lines): break
      line=lines[li]
      mo=mfexp.match(line)
      if mo:
        g=mo.groups()
        mfe=float(g[0])
        #~ print "mfe=%f"%mfe

      li=li+1
      line=lines[li]
      struc=line.strip()
      #~ print "struc='%s'"%struc

      li=li+1
      pairs=[]

      while li<len(lines):
        line=lines[li]
        #~ sys.stdout.write(line)
        mo=pairexp.match(line)
        if mo:
          g=mo.groups()
          #~ print g
          pairs.append((int(g[0])-1, int(g[1])-1))

        li=li+1

      complex=Complex([info["Sequence"]])
      params=info.get("Parameters", False)
      #~ sys.stderr.write("info=%s"%(info))
      if not params:
        params=info.get("Paramaters", False)
      #~ sys.stderr.write("Params=%s\n"%params)
      complex.add_struct("simplex", [0,0,0], params, float(info["Temperature"]), mfe, struc, pairs)
      #~ structs.append(Structure(info["Sequence"], float(info["Temperature"]), mfe, struc, pairs))
      #~ print complex
      foundstruc=False # Continue parsing file

    li=li+1

  return complex

#-------------------------------------------------------------------------------
def parseKeyComplexes(fname, complexnum, permutation):
  structs=[]
  f=delay_open(fname)
  info={}
  line = f.readline()
  while line != "":
    #~ print lines[li]
    mo=cexp.match(line)
    if mo:
      g=mo.groups()
      if g and len(g)==2:
        info[g[0]]=g[1]
        #~ print g
        if info.has_key("Temperature") or info.has_key("Temperature (C)"):
          if not info.has_key("Temperature"):
            info["Temperature"]=info["Temperature (C)"]
          break
    line = f.readline()
  #~ print info
  foundseq=False
  foundstruct=False
  seqs=[]
  combs=[]

  # Look for the sequences
  while line != "":

    if not foundseq: # We're still looking for a comment line
      mo=compid_exp.match(line)
      if mo:
        foundseq=True
        continue
    else:
      mo=seqid_exp.match(line) # Sequence
      if mo:
        #~ print "Found sequence"
        #~ print mo.groups()
        g=mo.groups()
        seqs.append(g[1])
      else:
        mo=temp_exp.match(line)
        if mo:
          break
    line = f.readline()

  line = f.readline()

  while line != "":
    comb=line.split()
    cur_comb = map(lambda x:x-1,map(int,comb))
    if cur_comb[0] == complexnum and cur_comb[1] == permutation:
        combs.append(cur_comb)
    line = f.readline()

  f.close()

  return seqs, combs

#-------------------------------------------------------------------------------
def parseMfeComplexes(fname, seqs, combs):
  f=delay_open(fname)
  #~ sys.stderr.write("fname=%s"%(fname))
  # It is quite possible this will be too large....
  info={}
  line = f.readline()

  while line != "":
    #~ print lines[li]
    mo=cexp.match(line)
    if mo:
      g=mo.groups()
      if g and len(g)==2:
        info[g[0]]=g[1]
        #~ print g
        if info.has_key("Temperature") or info.has_key("Temperature (C)"):
          if not info.has_key("Temperature"):
            info["Temperature"]=info["Temperature (C)"]
          break
    mo=cparam_exp.match(line)
    if mo:
      g=mo.groups()
      if g and len(g)==1:
        info["Parameters"]=g[0]
    line = f.readline()

  foundstruct=False
  line

  complex=Complex(seqs)
  line = f.readline()
    # Look for the structure
  for comb in combs:
    while line != "":
      #~ sys.stdout.write(line)

      if not foundstruct: # We're still looking for a comment line
        mo=compname_exp.match(line)
        if mo:
          compname=mo.groups()[0]
          mo2=compname_exp2.match(compname)
          g2 = mo2.groups([1,2])
          comp_id = int(g2[0]) - 1
          perm_id = int(g2[1]) - 1
          if comp_id == comb[0] and perm_id == comb[1]:
            foundstruct=True
          line = f.readline()
          continue
      else:
        mo=nexp.match(line) # Number or bases
        if mo:
          g=mo.groups()[0]
          #~ print g
          numbases=int(g)

          #~ print "numbases=%d"%numbases

        line = f.readline()
        if line == "": break
        mo=mfexp.match(line)
        if mo:
          g=mo.groups()
          #~ print g
          mfe=float(g[0])
          #~ print "mfe=%f"%mfe

        line = f.readline()
        struc=line.strip()
        #~ print "struc='%s'"%struc

        line = f.readline()
        pairs=[]

        while line != "":
          #~ sys.stdout.write(line)
          mo=pairexp.match(line)
          if mo:
            g=mo.groups()
            #~ print g
            pairs.append((int(g[0])-1, int(g[1])-1))
          else:
            mo=sexp.match(line)
            if mo: break

          line = f.readline()
        params=info.get("Parameters", False)
        #~ sys.stderr.write("info=%s\n"%(info))
        if not params:
          params=info.get("Paramaters", False)
        #~ sys.stderr.write("Params=%s\n"%params)
        complex.add_struct(compname, comb, params, float(info["Temperature"]), mfe, struc, pairs)

        foundstruct=False # Continue parsing file
        break

      line = f.readline()
  f.close()
  #~ print complex
  return complex

#-------------------------------------------------------------------------------
def parsePairs(fname, complexes, add_info=False):
  f=delay_open(fname)
  info={}
  si=0
  nested=False
  line = f.readline()

  #~ print "Parsing ppairs:"
  while line != "":
    mo=cexp.match(line)
    if mo:
      g=mo.groups()
      if g and len(g)==2:
        info[g[0]]=g[1]
        if info.has_key("Sequence") and \
          (info.has_key("Temperature") or info.has_key("Temperature (C)")) and \
          (info.has_key("Free") or info.has_key("Ensemble free energy")):
          if not info.has_key("Temperature"):
            info["Temperature"]=info["Temperature (C)"]
          if not info.has_key("Free"):
            info["Free"]=info["Ensemble free energy"]
          break
    line = f.readline()

  foundstruct=False
  line = f.readline()
  while line != "":

    if not foundstruct: # We're still looking for a comment line
      mo=nexp.match(line) # Number or bases
      if mo:
        g=mo.groups()[0]
        numbases=int(g)
        foundstruct=True
        probs=[]
        line = f.readline()
        #~ print "numbases=%d"%numbases
        continue

    else:
      if line == "": break

      while line != "":
        #~ sys.stdout.write(line)

        mo=ppairexp.match(line)
        mo_pseudo=ppairexp_pseudo.match(line)
        if mo_pseudo:
          g=mo_pseudo.groups()
          base=int(g[0])-1
          if len(probs)==0 or probs[-1][0]!=base:
            probs.append([base]) # format: [[base, (pairbase, prob), (pairbase, prob)],[base...]
          #~ print g
          fl=map(float,g[2:])
          cbase=int(g[1])-1
          probs[-1].append((cbase,fl[0], fl[1], fl[2]))
          #~ if f[1]>f[2] and cbase<numbases-1 and f[1]>0.5:
            #~ sys.stdout.write(line)
            #~ print "Possibly nested: %d -> %d: %f"%(base, cbase, f[1])

        elif mo:
          g=mo.groups()
          base=int(g[0])-1
          if len(probs)==0 or probs[-1][0]!=base:
            probs.append([base]) # format: [[base, (pairbase, prob), (pairbase, prob)],[base...]
          #~ print g
          probs[-1].append((int(g[1])-1,float(g[2])))
          #~ print probs[-1]
        else:
          break
        line = f.readline()
      #~ for p in probs:
        #~ print p
      complexes.structs[si].set_probs(probs)
      if add_info:
        complexes.structs[si].temp=float(info["Temperature"])
        complexes.structs[si].mfe=float(info["Free"])

      foundstruc=False # Continue parsing file
      nested=complexes.structs[si].struc.find('{')>=0
      si=si+1

    line = f.readline()
  f.close()
  return nested

#-------------------------------------------------------------------------------
def parsePairsComplexes(fname, complexes):
  f=delay_open(fname)
  info={}
  line = f.readline()

  #~ print "Parsing ppairs:"
  while line != "":
    mo=cexp.match(line)
    if mo:
      g=mo.groups()
      if g and len(g)==2:
        info[g[0]]=g[1]
        #~ print g
        if info.has_key("Temperature") or info.has_key("Temperature (C)"):
          if not info.has_key("Temperature"):
            info["Temperature"]=info["Temperature (C)"]
          break
    line = f.readline()

  #~ print info
  foundstruct=False
  line = f.readline()
  for c in range(len(complexes.complex_perm)):
    probs=[]
    while line != "":
      #~ sys.stdout.write(line)

      if not foundstruct: # We're still looking for a comment line
        mo=compname_exp.match(line)
        if mo:
          compname=mo.groups()[0]
          mo2=compname_exp2.match(compname)
          g2 = mo2.groups([1,2])
          comp_id = int(g2[0]) - 1
          perm_id = int(g2[1]) - 1
          if comp_id == complexes.complex_perm[c][0] and perm_id == complexes.complex_perm[c][1]:
            foundstruct=True
          line = f.readline()
          continue
      else:
        mo=nexp.match(line) # Number or bases
        if mo:
          g=mo.groups()[0]
          #~ print g
          numbases=int(g)

          #~ print "numbases=%d"%numbases

        line = f.readline()
        if line == "": break

        while line != "":
          #~ sys.stdout.write(line)

          mo=ppairexp.match(line)
          if mo:
            g=mo.groups()
            base=int(g[0])-1
            if len(probs)==0 or probs[-1][0]!=base:
              probs.append([base]) # format: [[base, (pairbase, prob), (pairbase, prob)],[base...]
            #~ print g
            probs[-1].append((int(g[1])-1,float(g[2])))
            #~ print probs[-1]
          else:
            break
          line = f.readline()
        complexes.structs[c].set_probs(probs)
        #~ print c,complexes.structs[c]
        foundstruct=False # Continue parsing file
        break # Look for next block

      line = f.readline()
  f.close()

#-------------------------------------------------------------------------------
def parseDesign(fname):
  f=delay_open(fname)
  lines=f.readlines() # Read the file into memory, should be small enough
  f.close()
  info={}

  info['structure']=lines[0]
  info['sequence']=lines[1]
  prob=[]
  li=2
  while li<len(lines):
    prob.append(float(line[li]))
  info['prob']=prob

  return info

def parseDesignJson(fname):
  f=delay_open(fname)
  data=f.read()
  info=json.read(data)
  return info

if __name__=='__main__':
  structs=parseMfe("example.mfe")
  parsePairs("example.ppairs", structs)
