#!/usr/bin/env python
# -*- coding: iso-8859-15 -*-
#global nh nb rhc rbc rhelix rdh rcap ncap thcap dzb bturn brise dthb dsb dz dth ...
#       cmap rdye ddye ndye majorgroove minorgroove dthgroove strutangle strutlength ...
#       strutrise proptwist inclination extend5prime extend3prime extendcoil rcaphead;


# ------------------------------------------------------------------------------
# Process command-line options
import sys
import Options
import NudrawOptions

option_vars=Options.parse(__file__, sys.argv,
  NudrawOptions.def_option_vars, NudrawOptions.option_vars_doc)
cmd_options=Options.ComplexOptions(option_vars)

# Redirect output if needed
if cmd_options.stdout:
  sys.stdout=open(cmd_options.stdout, "w")

if cmd_options.stderr:
  sys.stderr=open(cmd_options.stderr, "w")

if cmd_options.quiet:
  sys.stdout=open("/dev/null", "w")
  sys.stderr=open("/dev/null", "w")

# ------------------------------------------------------------------------------
import time
starttime_begin=time.time()
print 'Importing libraries'
import os
import math
import warnings
warnings.simplefilter("ignore",DeprecationWarning)

import numpy
import numpy.linalg
import scipy
import string
import traceback
import re

if cmd_options.show2d or cmd_options.svgfile2d or \
    cmd_options.pngfile2d or cmd_options.dotfile or \
    cmd_options.jsonout or \
    cmd_options.show25d or cmd_options.svgfile25d or \
    cmd_options.pngfile25d or \
    cmd_options.json25d or cmd_options.jsonin25d or \
    cmd_options.optionsjson:
  import matplotlib

  if len(cmd_options.backend)>0:
    matplotlib.use(cmd_options.backend)
    matplotlib.rcParams['svg.fonttype'] = 'svgfont'
  elif cmd_options.show2d or cmd_options.show25d or cmd_options.interactive: # Non-interactive only
    matplotlib.use("Qt4Agg")
  else:
    matplotlib.use("agg")
  import pylab
  import Export
else:
  pylab=None
import copy
#from Export import cmap, color_names, colors, icolors

try:
  import Coil
except:
  Coil=None
import Utils
import Layout
import pdb
import json
import JsonLayout
import Geometric
from Classes import *
import Newton
from spline import spl, polyfit, findpoint, circlefit
import Utils3D

#-------------------------------------------------------------------------------
def toc(starttime):
  endtime=time.time()
  # print "Elapsed time is %2.6f\n"%(endtime-starttime)

toc(starttime_begin)

#-------------------------------------------------------------------------------
# Global variables
#

import Globals as gg

#-------------------------------------------------------------------------------
# Start main program

global_count=1
while global_count>0:
  global_count-=1

  try:
    if cmd_options.optionsjson:
      print "Parsing JSON parameters from stdin"
      sys.stderr.write("status:READY\n")
      option_vars_json = Options.parse_json(__file__, sys.argv,
        option_vars, NudrawOptions.option_vars_doc, sys.stdin)
      if type(option_vars_json)==dict and len(option_vars_json)==0:
        global_count=1
        sys.stderr.write("status:ERROR parsing input\n")
        continue
      if not option_vars_json: # EOF reached
        sys.stderr.write("status:ERROR parsing input\n")
        continue
      try:
        cmds=option_vars_json["commands"] #{"commands": ["exit"]}
        if cmds[0]=="exit":
          break
      except KeyError:
        pass

      cmd_options=Options.ComplexOptions(option_vars,
                  option_vars_json["options"])
      if cmd_options.debug: print cmd_options
      global_count=1 # repeat outer loop once more

    sys.stderr.write("status:BUSY\n")
    starttime_begin=time.time()
    seq="" # Sequence
    s=""   # Structure

    if cmd_options.mfefile:
      print "Parsing MFE and/or Prob files"
      starttime=time.time()
      if cmd_options.jsonin or cmd_options.structure or cmd_options.jsonin25d:
        raise RuntimeError("Multiple input options specified (mfefile, jsonin or structure)")
      if not cmd_options.complex_permutation:
        complex=Utils.parseMfe(cmd_options.mfefile)
        #~ print "complex=%d"%len(complex.names)
        complexnum=0
        permutation=0
        titlebottom="Free energy of secondary structure: %2.2f kcal/mol"%complex.structs[complexnum].mfe
        titletop="MFE structure at %2.1f C"%complex.structs[complexnum].temp
        if complex.params.find('DNA')>=0:
          material=2
        else:
          material=1

      else:
        if not cmd_options.keyfile: raise RuntimeError("Need keyfile if complexes specified")

        complexnum=int(cmd_options.complex_permutation[0])-1
        permutation=int(cmd_options.complex_permutation[1])-1
        seqs,combs=Utils.parseKeyComplexes(cmd_options.keyfile,complexnum,permutation)

        complex=Utils.parseMfeComplexes(cmd_options.mfefile, seqs, combs)
        titlebottom="Free energy of secondary structure: %2.2f kcal/mol"%complex.structs[0].mfe
        titletop="MFE structure at %2.1f C"%complex.structs[0].temp

      print "Done parsing MFE and/or Prob files"
      toc(starttime)
    elif cmd_options.designjson:
      design_info=Utils.parseDesignJson(cmd_options.designjson)
      s=design_info['structure']
      seq=design_info['sequence']
      if 'energy' in design_info:
        titlebottom=design_info['energypreamble']+" %2.2f kcal/mol"%float(design_info['energy'])
      elif 'target free energy' in design_info:
        titlebottom=design_info['energypreamble']+" %2.2f kcal/mol"%float(design_info['target free energy'])
      else:
        titlebottom="Normalized ensemble defect: "+" %2.2f"% (100*float(design_info['normalized ensemble defect']))
        titlebottom += '%'
      if design_info.get("temperature (C)"):
        titletop="Target structure at %s"%design_info["temperature (C)"]
      complex=False
    else:
      if cmd_options.structure:
        s=cmd_options.structure
      if cmd_options.sequence:
        seq=cmd_options.sequence

      titlebottom=""
      titletop=""
      complex=False
      if len(s)>0 and len(seq)>0 and cmd_options.probfile:
        seqOrd,lst=Utils.combFromSeq(seq)
        complex=Utils.Complex(seqOrd)
        pairs=Utils.pairsFromStruc(s)
        comb=[0,0]+lst
        # Dummy params
        params=""
        temp=37.0
        mfe=-200.0
        complexnum=0
        permutation=0
        complex.add_struct("One", comb, params, temp, mfe, s, pairs)
        Utils.parsePairs(cmd_options.probfile, complex, add_info=True)

    if cmd_options.toptitle and cmd_options.toptitle!="-":
      titletop=cmd_options.toptitle
    if cmd_options.bottomtitle  and cmd_options.bottomtitle!="-":
      titlebottom=cmd_options.bottomtitle

    # Alternative structures read from JSON file
    json_base=[]
    json_loop=[]
    json_helix=[]

    # Read 2D JSON file
    if cmd_options.jsonin:
      if cmd_options.jsonin25d:
         raise RuntimeError("Cannot load 2D and 2.5D JSON file!")
      json_data, s, seq = JsonLayout.loadJson(fname=cmd_options.jsonin)
      s=s.replace(' ','+') # Somewhere Javascript converts +'s to spaces
      seq=seq.replace(' ','+')

    # READ 2.5D JSON file
    if cmd_options.jsonin25d:
      if cmd_options.jsonin:
        raise RuntimeError("Cannot load 2D and 2.5D JSON file!")
      json_data, s, seq = JsonLayout.loadJson25d(fname=cmd_options.jsonin25d)
      s=s.replace(' ','+') # Somewhere Javascript converts +'s to spaces
      seq=seq.replace(' ','+')

    if len(s)==0 and not complex:
        raise RuntimeError("No structure specified")

    if cmd_options.domaincolors:
      colspec_re=re.compile('[a-fA-F0-9]{6}')
      dc=cmd_options.domaincolors
      gg.domainColors = [0] * len(dc)
      gg.domainCMap = [0] * len(dc)
      for i in range(len(dc)):
        if dc[i]=='n': continue
        if not colspec_re.match(dc[i]):continue
        gg.domainColors[i]=dc[i]
        gg.domainCMap[i]=map(lambda j: float(int(dc[i][j*2:(j+1)*2],16))/255.0, range(3))


    #-------------------------------------------------------------------------------
    planar=True
    colorscheme = 1 # 1: by strands, 2: by coil or helix
    smartrot = True  # 1: rotate leaves of tree for aesthetics
    randomendcoils = False  # make coils random on either side of nicks
    stacknickedhelices = True # stack nicked helices
    if not cmd_options.jsonin and not cmd_options.jsonin25d:
      if cmd_options.material=="rna":
        material = 1 # 1: A-DNA, 2: B-DNA
      elif cmd_options.material=="dna":
        material = 2
      else:
        raise RuntimeError("Material must be one of rna or dna")
    else:
      material=json_data["properties"]["material"]

    gg.init_globals(material, cmd_options.show3d or cmd_options.tachyon3d or
       cmd_options.json3d or cmd_options.show25d or cmd_options.json25d or cmd_options.pov3d or
       cmd_options.pngfile25d or cmd_options.svgfile25d or cmd_options.jsonin25d)
    eye=[numpy.identity(0), numpy.identity(1),numpy.identity(2), numpy.identity(3)]
    I3=eye[3]

    #-------------------------------------------------------------------------------
    starttime=time.time()
    print 'Logical Definition of Loops'

    # Assign probabilities from probfile.
    if not complex and len(s)>0 and len(seq)>0 and cmd_options.probfile:
      seqOrd,lst=Utils.combFromSeq(seq)
      complex=Utils.Complex(seqOrd)
      pairs=Utils.pairsFromStruc(s)
      comb=[0,0]+lst
      # Dummy params
      params=""
      complexnum=0
      temp=-200
      mfe=-200
      permutation=0

      complex.add_struct("One", comb, params, temp, mfe, s, pairs)
      Utils.parsePairs(cmd_options.probfile, complex, add_info=True)

    if complex:
      ncombs=len(complex.structs)
    else:
      ncombs=1

    pseudo=False

    for com in range(ncombs):
      if complex:
        #~ print complex.complex_perm[com], complexnum, permutation
        if not (complex.complex_perm[com]==[complexnum,permutation]):
          continue
        #~ else:
          #~ print "rendering complex.complex_perm=%s"%complex.complex_perm
        s=complex.structs[com].struc
        seq=string.upper(complex.structs[com].seq)
        if cmd_options.probfile and cmd_options.mfefile:
          if cmd_options.complex_permutation:
            Utils.parsePairsComplexes(cmd_options.probfile, complex)
          else:
            pseudo=Utils.parsePairs(cmd_options.probfile, complex)
            if pseudo:
              s=len(s)*Layout.nopair
              if cmd_options.debug:
                print "*** PSEUDOKNOT DETECTED!!!"

          if not cmd_options.bottomtitle:
            titlebottom="Free energy of secondary structure: %2.2f kcal/mol"%complex.structs[com].mfe
          if not cmd_options.toptitle:
            titletop="MFE structure at %2.1f C"%complex.structs[com].temp
        if complex.params.find('DNA')>=0:
          material=2
        elif complex.params.find('RNA')>=0:
          material=1
        gg.init_globals(material, cmd_options.show3d or cmd_options.tachyon3d or
           cmd_options.json3d or cmd_options.json25d or cmd_options.pov3d)

      # set up logical representation of strand with nick, base and strand
      # numbers
      #

      # Our indices begin at 0, not 1!
      ibase = 0
      iunit = 0
      istrand = 0
      unit=[]
      base=[]
      loop=[]
      helix=[]
      coil=[]

      # Set up unit structure
      nunit, nbase, nstrand = Layout.setup_units(s, seq, unit, base, cmd_options.debug)

      # Set up loop structure
      nloop = Layout.setup_loops(unit, loop, cmd_options.debug)

      toc(starttime)

      print("Geometric Definition of Loops")

      Geometric.calc(gg.rdh, gg.dzb, loop, base, stacknickedhelices, gg.dsb,
                    show3d=cmd_options.show3d, debug=cmd_options.debug,
                    circle=(cmd_options.unpairedcircle or pseudo))

      if pseudo:
        Layout.setup_pseudopairs(unit, base, complex.structs[0])

      #-----------------------------------------------------------------------------
      toc(starttime)
      print("Assign Letters and Probabilites")
      totalbases=len(base)

      # Assign letters
      if len(seq)>0:
        ibase=0
        iseq=0
        seq=seq.upper()

        # Make sure correct letters are used
        if material==1: # RNA
          ttable=string.maketrans('tT','uU')
        elif material==2: # DNA
          ttable=string.maketrans('uU', 'tT')
        seq=seq.translate(ttable)

        oneseq=reduce(lambda x,y:x+y, seq.split('+'))
        while ibase<len(base) and iseq<len(oneseq):
          if oneseq[iseq] in ['A', 'C', 'G', 'T', 'U']:
            base[ibase].type=oneseq[iseq]
          iseq=iseq+1
          ibase=ibase+1

      # Handle domain numbers and ids
      if cmd_options.domains:
        onedom=reduce(lambda x,y:x+y,cmd_options.domains.split('+'))
        onedomid=reduce(lambda x,y:x+y,cmd_options.domainids.split('+'))
        splitInt=filter(None, onedom.split(","))
        splitIntId=filter(None, onedomid.split(","))

        ibase=0
        iseq=0

        if len(splitInt)==len(base): # We got a comma-seaprated list of ints
          while ibase<len(base) and iseq<len(splitInt):
            base[ibase].domain=int(splitInt[iseq])
            iseq=iseq+1
            ibase=ibase+1
        else:
          mo=re.search('([A-Z\+]+)',cmd_options.domains)
          if not mo or not mo.groups() or len(mo.groups()[0])!=len(cmd_options.domains):
            raise RuntimeError("Domain spec %s is invalid"%cmd_options.domains)
          while ibase<len(base) and iseq<len(onedom):
            base[ibase].domain=ord(onedom[iseq])-ord('A')
            iseq=iseq+1
            ibase=ibase+1

        if len(splitIntId)==0:
            splitIntId = splitInt

        if len(splitIntId)==len(base): # We got a comma-seaprated list of ints
          ibase=0
          iseq=0
          while ibase<len(base) and iseq<len(splitIntId):
            base[ibase].domainid=int(splitIntId[iseq])
            iseq=iseq+1
            ibase=ibase+1
        else:
          raise RuntimeError("Domain IDs must be specified as a comma-separated list")

      #-----------------------------------------------------------------------------
      # Assign probabilities
      if complex :
        ui=0  # unit index (index of the current unit?)
        base_index = 0
        bpi=0 # essentially the line of the ppairs file we are considering
        spi=0 # the next index of the target structure specification we are considering
        st=complex.structs[com]
        restarted=False
        pairs_index = 0
        st.pairs.sort()
        new_pairs = [(x,totalbases) for x in range(0,totalbases)]
        for new_index in range(0,totalbases):
          if pairs_index < len(st.pairs) and new_index == st.pairs[pairs_index][0]:
            paired_index = st.pairs[pairs_index][1]
            new_pairs[new_index] = (new_index,paired_index)
            new_pairs[paired_index] = (paired_index,new_index)
            pairs_index += 1

        sl=len(st.pairs)
        st.probs.sort()

        for base_index in range(0,totalbases):
          while bpi < len(st.probs) and st.probs[bpi][0] < base_index:
            bpi += 1
          #~ print

          while ui<len(unit) and (unit[ui].base < base_index or type(unit[ui].type)!=str) :
            ui=ui+1
          while bpi < len(st.probs) and st.probs[bpi][0] == base_index and unit[ui].prob < 0:
            for pair in st.probs[bpi][1:]:       # pair = baseind, probability
              if new_pairs[base_index][1]==pair[0]: # Use this as the probability to shade base with
                unit[ui].prob=pair[1]
                if unit[ui].pair < len(unit) and unit[ui].pair >= 0:
                  unit[unit[ui].pair].prob=pair[1]
            bpi += 1

          ui += 1

        # If we have a probability source, set all unset probs to 0.0
        if len(st.probs)>0:
          ui=0
          while ui<len(unit):
            if unit[ui].prob<0:
              unit[ui].prob=0.0
            ui+=1


      #-----------------------------------------------------------------------------
      # Assign probabilities from a designjson file
      if cmd_options.designjson :
        ui=0  # unit index
        bi=0  # base
        if 'nucleotide probabilities' in design_info:
          prob=design_info['nucleotide probabilities']
        elif 'nucleotide defects' in design_info:
          prob=design_info['nucleotide defects']
        else:
          prob=design_info['prob']
        while bi < len(base) and ui < len(unit):
          while ui<len(unit) and not (unit[ui].base==bi and type(unit[ui].type)==str): # Advance unit index
            ui=ui+1
            if ui >= len(unit):
              ui=0

          unit[ui].prob=prob[bi]
          bi+=1
          ui+=1

      #-----------------------------------------------------------------------------
      # Assign probabilities from a simpleprobfile
      if cmd_options.simpleprobfile:
        ui=0  # unit index
        bi=0  # base
        fo=open(cmd_options.simpleprobfile)
        lines=fo.readlines()
        fo.close()
        while bi < len(base) and ui < len(unit):
          while ui<len(unit) and not (unit[ui].base==bi and type(unit[ui].type)==str): # Advance unit index
            ui=ui+1
            if ui >= len(unit):
              ui=0

          if cmd_options.debug:
            print "Setting prob(%d)=%f"%(bi, float(lines[bi]))
          unit[ui].prob=float(lines[bi])
          bi+=1
          ui+=1


      #-------------------------------------------------------------------------
      if cmd_options.targetfile:
        tgf=open(cmd_options.targetfile,"w")
        tgf.write("% Target file created by nudraw.py")
        tgf.write("%d\n"%len(base))
        tgf.write("%s\n"%structure)
        bi=0

        while ui<len(unit):
          if unit[ui].base>=0 and type(unit[ui].type)==str:
            if unit[ui].pair>=0:
              tgf.write("%d \t%d\n"%(unit[ui].base+1, unit[unit[ui].pair].base+1))

      #-------------------------------------------------------------------------
      if cmd_options.debugprobfile:
        dbf=open(cmd_options.debugprobfile,"w")
        dbf.write('{"file" : "%s",\n'%cmd_options.probfile+
                  ' "desc" :"Probabilities extracted from ppairs file",\n'+
                  ' "pairs": [\n')
        paired_probs=[]
        unpaired_probs=[]
        desc=""
        ui=0

        while ui<len(unit):
          if unit[ui].base>=0 and type(unit[ui].type)==str:
            if unit[ui].pair>=0:
              new_entry=[unit[ui].base+1, unit[unit[ui].pair].base+1, unit[ui].prob]
              paired_probs.append(new_entry)
            else:
              new_entry=[unit[ui].base+1, -1, unit[ui].prob]
              unpaired_probs.append(new_entry)
            dbf.write(str(new_entry)+",\n")
          ui+=1
        dbf.write(']}\n')
        dbf.close()

      #-------------------------------------------------------------------------
      if cmd_options.jsonin and json_data:
        # Update base coordinates from JSON file
        update_prob = not cmd_options.probfile and not cmd_options.simpleprobfile
        JsonLayout.updateBases(json_data, base, unit, update_prob)

        # Update loop coordinates from JSON file
        JsonLayout.updateLoops(json_data, loop, base)


      # Plot loop structures
      if pylab and cmd_options.show2d and cmd_options.debug:
        pylab.clf()
        toc(starttime)
        print("Plotting loop structures")

        for iloop in range(nloop):
           nside = loop[iloop].nside #  number of stems in loop
           for jside in range(nside):
              xc = loop[iloop].geo[jside].xc
              nc = loop[iloop].geo[jside].nc
              pylab.plot([xc[0]],[xc[1]],'rs')

           th = numpy.linspace(0,2*math.pi,num=50)

           xcirc = loop[iloop].center[0] + loop[iloop].radius*scipy.cos(th)
           ycirc = loop[iloop].center[1] + loop[iloop].radius*scipy.sin(th)
           pylab.plot(xcirc,ycirc,'g-')
           pylab.text(loop[iloop].center[0],loop[iloop].center[1],'L%d'%(iloop))

        for ibase in range(nbase):
            pylab.text(base[ibase].x[0], base[ibase].x[1], '%d'%(ibase))
            pylab.plot([base[ibase].x[0]], [base[ibase].x[1]], 'mo')

      #-----------------------------------------------------------------------------
      toc(starttime)
      #
      # define helices with links to previous loop number and side number
      # and next loop number and side number
      # define coil segments with link to relevant loop number and side number
      #
      #
      print('Define Loops')
      starttime=time.time()

      iloop = 1
      ihelix = 0
      while iloop < nloop:
      #
      #   define helices
      #   each helix has a pointer to the previous loop before the helix
      #   and the next loop after the helix as well as pointers to the
      #   appropriate side in these neighboring loops: these pointers
      #   are used to link helix ends to coil ends below
      #
          helix.append(helixClass())
          helixlength = 1
          firstloop = iloop

          l=loop[iloop]
          if cmd_options.debug:
            print "\nL%d:\nsidenbase=%s\nsideunit=%s"%(iloop,l.sidenbase, l.sideunit)
            print "sidebase=%s\ncoilnum=%s"%(l.sidebase, l.coilnum)
            print "toloop=%s\ntohelix=%s"%(l.toloop, l.tohelix)
            print "center=%s\nradius=%s"%(l.center, l.radius)
            print "strand=%s"%(l.strand)


          helixobj=helix[-1]
          while (loop[iloop].nside == 2 and \
                 loop[iloop].sidenbase[0] == 0 and \
                 loop[iloop].sidenbase[1] == 0 and \
                 loop[iloop].nick[1,0] != -3):

            l=loop[iloop]
            if cmd_options.debug:
              print "\n**L%d:\nsidenbase=%s\nsideunit=%s"%(iloop,l.sidenbase, l.sideunit)
              print "sidebase=%s\ncoilnum=%s"%(l.sidebase, l.coilnum)
              print "toloop=%s\ntohelix=%s"%(l.toloop, l.tohelix)
              print "center=%s\nradius=%s"%(l.center, l.radius)
              print "strand=%s"%(l.strand)

            helixlength = helixlength + 1
            helixobj.bases.append([l.sidebase[0],l.sidebase[1]+1])

            iloop = iloop + 1 # because loops in a helix are always numbered consecutively

          helix[ihelix].npair = helixlength
          l=loop[iloop]

          sn=0
          while sn<len(l.sidebase)-1 and l.sidebase[sn+1]>=0:
            sn+=1

          if cmd_options.debug:
            print "\n**L%d:\nsidenbase=%s\nsideunit=%s"%(iloop,l.sidenbase, l.sideunit)
            print "sidebase=%s\ncoilnum=%s"%(l.sidebase, l.coilnum)
            print "toloop=%s\ntohelix=%s"%(l.toloop, l.tohelix)
            print "center=%s\nradius=%s"%(l.center, l.radius)
            print "strand=%s"%(l.strand)

          helixobj.bases.append([l.sidebase[0], l.sidebase[sn]+l.sidenbase[sn]+1])

          firstside = loop[firstloop].nside-1
          prevloop = loop[firstloop].toloop[firstside]
          prevside = loop[firstloop].toside[firstside]
          helix[ihelix].prevloop  = prevloop
          helix[ihelix].prevside  = prevside
          helix[ihelix].xc     =  loop[firstloop].geo[firstside].xc.copy() # coord of helix axis start
          helix[ihelix].nc     = -loop[firstloop].geo[firstside].nc.copy() # helix axis at start (into helix)
          helix[ihelix].thc    = 180./math.pi*(math.pi - gg.strutangle)
                  #rotate so base-pair struct is parallel to plane of loop it borders
                  # reversed sign of thc relative to nupack6 to get first
                  # basepair strut into the plane of the loop it exits from
                  # (works for A-DNA or B-DNA)

          if not smartrot:
              helix[ihelix].thc = helix[ihelix].thc  + 90. + \
              180./math.pi*numpy.arctan2(helix[ihelix].nc[1],helix[ihelix].nc[0])

          # rotate so that projection from canonical
          # u_z orientation creates helix that joins the loop
          # with the same rotation regardless of helix axis

          #
          # make pointers from loop sides to helix numbers: from root to
          # leaves (don't need pointers from loops to helix entering the loop
          # -- but could add below using nextloop, nextside
          #
          loop[prevloop].tohelix[prevside] = ihelix

          if helixlength == 1: #special treatment for isolated base pairs
              helix[ihelix].nextloop = firstloop
              helix[ihelix].nextside = firstside
          else:                # regular helix with two or more stacked pairs
              lastloop  = iloop - 1
              lastside = 0
              helix[ihelix].nextloop  = loop[lastloop].toloop[lastside]
              helix[ihelix].nextside  = loop[lastloop].toside[lastside]


          if stacknickedhelices and loop[prevloop].nickedhelix:
             helix[ihelix].thc = helix[ihelix-1].thc + \
             180./math.pi*gg.dthb*helix[ihelix-1].npair

          # this takes advantage of fact that helices connected by a nick are consecutive in the helix array
          # pass on rotation info from the previous helix (can't do this by
          # defining reference angles for the nicked loops because that would
          # pass on non-2D orientations to the subsequent loops)

          ihelix = ihelix + 1
          iloop = iloop + 1

      nhelix = ihelix


      #-----------------------------------------------------------------------------
      #
      #   define coils
      #
      toc(starttime)
      print('Define Coils')
      starttime=time.time()

      iloop = 0
      icoil = 0

      while iloop < nloop:
          l=loop[iloop]
          if cmd_options.debug:
            print "\nL%d:\nsidenbase=%s\nsideunit=%s"%(iloop,l.sidenbase, l.sideunit)
            print "sidebase=%s\ncoilnum=%s"%(l.sidebase, l.coilnum)
            print "toloop=%s\ntohelix=%s"%(l.toloop, l.tohelix)
            print "center=%s\nradius=%s"%(l.center, l.radius)
            print "strand=%s\nnside=%d"%(l.strand, l.nside)
            print "nick=%s"%(l.nick)

          if not (loop[iloop].nside == 2 and loop[iloop].sidenbase[0] == 0 and \
                loop[iloop].sidenbase[1] == 0) or loop[iloop].nick[1,0] == -3 \
                or loop[iloop].nick[0,0] == -5 or loop[iloop].tohelix[0]>=0: # if not part of helix
              for jside in range(loop[iloop].nside):
                  coil.append(coilClass(loop[iloop].nside))
                  coil[icoil].nbase = loop[iloop].sidenbase[jside] + 2 # including end pts
                  coil[icoil].loopnum = iloop # check
                  coil[icoil].loopside = jside # check
                  coil[icoil].center = loop[iloop].center
                  loop[iloop].coilnum[jside] = icoil


                  coil[icoil].strand = loop[iloop].strand[jside]
                  coil[icoil].nick = loop[iloop].nick[:,jside]
                  if loop[iloop].nside==1:
                    coil[icoil].nick[1] = loop[iloop].nick[1,jside+1]
                  if coil[icoil].nick[0] == -5:
                      coil[icoil].tohelix[0,:] = numpy.array([-1,-5],numpy.int32) # could change zero to istrand if need info later
                                                      # or just not use tohelix for strand ends and rely on coil.nick
                  if coil[icoil].nick[1] == -3:
                      coil[icoil].tohelix[1,:] = numpy.array([-1, -3], numpy.int32)

                  if loop[iloop].nick[0,jside] == -5:
                      startval = 0
                  else:
                      startval = 1

                  if loop[iloop].nick[1,jside] == -3:
                      stopval = loop[iloop].sidenbase[jside] + 1
                  else:
                      stopval = loop[iloop].sidenbase[jside]

                  jcount = 0
                  coil[icoil].setSplineSize(stopval-startval+1)
                  for jval in range(startval,stopval+1): # store spline pts for coil, including end point if a nick
                      jbase = jval + loop[iloop].sidebase[jside]
                      coil[icoil].xsplinepts[:,jcount] = base[jbase].x
                      coil[icoil].bases.append(jbase)
                      jcount = jcount + 1

                  coil[icoil].nsplinepts = jcount
                  icoil = icoil + 1

          iloop = iloop + 1

      ncoil = icoil

      #
      # drawing reality checks
      if cmd_options.debug and pylab and cmd_options.show2d:
        for icoil in range(ncoil): # draw coil bases for reality check
            if coil[icoil].nsplinepts>0:
                for ipt in range(coil[icoil].nsplinepts):
                    pylab.plot([coil[icoil].xsplinepts[0,ipt]], [coil[icoil].xsplinepts[1,ipt]],'mv')
                pylab.text(coil[icoil].xsplinepts[0,0]+3.,coil[icoil].xsplinepts[1,0],'C%d'%icoil)

        for ihelix in range(nhelix):
           pylab.plot([helix[ihelix].xc[0]], [helix[ihelix].xc[1]],'cs')
           pylab.text(helix[ihelix].xc[0]+2, helix[ihelix].xc[1],'H%d'%ihelix, color="red")

      #-----------------------------------------------------------------------------
      #
      # make linked list between helix strand ends and coil strand ends
      # loop of helices
      # follow pointers helix -> end loops -> loop sides -> coil
      #
      # helix definitions in drawing program:
      # a chain runs from 5' -> 3' (indices 1 and 2, respectively, in drawing program)
      # b chain rums from 3' -> 5' (indices 1 and 2, respectively, in drawing program)
      #
      # helix definitions in automatic structure generation program:
      # a1 = 1, a2 = 2, b1 = 3, b2 = 4
      #
      for ihelix in range(nhelix):
          prevloop  = helix[ihelix].prevloop
          side1     = helix[ihelix].prevside
          ic1       = loop[prevloop].coilnum[side1]

          side3     = side1 + 1
          nside     = loop[prevloop].nside
          if side3 >= nside:
              side3 = 0

          ic3       = loop[prevloop].coilnum[side3]

          nextloop  = helix[ihelix].nextloop
          side4     = helix[ihelix].nextside
          ic4       = loop[nextloop].coilnum[side4]

          side2     = side4 + 1
          nside     = loop[nextloop].nside
          if side2  >= nside:
              side2 = 0

          ic2       = loop[nextloop].coilnum[side2]

          # coil number and position
          helix[ihelix].tocoil = numpy.array([[ic1, 1], [ic2, 0], [ic3, 0], [ic4, 1]])
          # strand number for chains a and b
          helix[ihelix].strand = numpy.array([coil[ic1].strand, coil[ic3].strand])

          coil[ic1].tohelix[1,:] = numpy.array([ihelix, 0]) # helix number and position (1->4)
          coil[ic2].tohelix[0,:] = numpy.array([ihelix, 1])
          coil[ic3].tohelix[0,:] = numpy.array([ihelix, 2])
          coil[ic4].tohelix[1,:] = numpy.array([ihelix, 3])

          #
          # if coil is at strand end and has no length then just draw cap on
          # relevant helix strand end and don't draw coil
          if coil[ic1].nick[0] == -5 and coil[ic1].nbase == 1:
            helix[ihelix].cap[0] = -5

          if coil[ic2].nick[1] == -3 and coil[ic2].nbase == 1:
            helix[ihelix].cap[1] = -3

          if coil[ic3].nick[1] == -3 and coil[ic3].nbase == 1:
            helix[ihelix].cap[2] = -3

          if coil[ic4].nick[0] == -5 and coil[ic4].nbase == 1:
            helix[ihelix].cap[3] = -5
      toc(starttime)
      print('Export figures')
      starttime=time.time()

    #-------------------------------------------------------------------------------
    # Export 2D layout and pictures
      if cmd_options.jsonout or cmd_options.svgfile2d or \
        cmd_options.dotfile or cmd_options.pngfile2d or cmd_options.show2d:

        if cmd_options.numberinterval in ['auto','0', '-1']: # Auto line numbering
          intok=[1.0,2.0,5.0,10.0,15.0,25.0,40.0,50.0]
          realintv=intv= max(1.0,math.ceil(math.fabs(math.log(float(nbase)/4)))*2)
          mindiff=1000
          for ok in intok:
            newdiff=math.fabs(intv-ok)
            if newdiff>mindiff: break
            if newdiff<mindiff:
              mindiff=newdiff
              realintv=ok

          numberinterval=realintv
        else:
          numberinterval=int(cmd_options.numberinterval)

        # Export with and without multiple file options
        if cmd_options.filecounter2d:
          for i in range(len(cmd_options.filecounter2d)):
            exvars={}
            for opt in ["svgfile2d", "pngfile2d"]:
              exvars[opt]=Options.subst(cmd_options[opt],cmd_options.filecounter2d[i])
            if i==0:
              for opt in ["jsonout","dotfile"]:
                exvars[opt]=Options.subst(cmd_options[opt],cmd_options.filecounter2d)
            else:
              exvars["jsonout"]=""
              exvars["dotfile"]=""

            for opt in ["drawbases", "colorbar", "colorbarspace", "colorbaseprob",
                        "colorbaseid", "baseidbar", "colorstrands", "drawbasenumbers",
                        "drawbaseticks", "colordomains", "labeldomains", "domainkey"]:
              exvars[opt]=cmd_options[opt][i]

            tmp_top = titletop
            if cmd_options.hidetoptitle and cmd_options.hidetoptitle[i]:
              tmp_top = ""

            tmp_bottom = titlebottom
            if cmd_options.hidebottomtitle and cmd_options.hidebottomtitle[i]:
              tmp_bottom = ""

            exopt=Options.ComplexOptions(exvars)


            Export.TwoD(base, unit, helix, loop, seq, s, gg.icolors, gg.dsb, material,
              titlebottom=tmp_bottom, titletop=tmp_top,
              jname=exopt.jsonout, svgname=exopt.svgfile2d,
              dotname=exopt.dotfile, pngname=exopt.pngfile2d,
              drawbases=exopt.drawbases, colorbar=exopt.colorbar,
              colorbarspace=exopt.colorbarspace, colorbaseprob=exopt.colorbaseprob,
              colorstrands=exopt.colorstrands, colorbarlabel=cmd_options.legend_label,
              show2d=cmd_options.show2d, labeldomains=exopt.labeldomains,
              drawbasenumbers=exopt.drawbasenumbers, drawbaseticks=exopt.drawbaseticks,
              baseidbar=exopt.baseidbar, colordomains=exopt.colordomains,
              numberinterval=numberinterval, colorbaseid=exopt.colorbaseid,
              pseudo=pseudo, debug=cmd_options.debug, domainnames=cmd_options.domainnames,
              domainkey=exopt.domainkey, dpi=int(cmd_options.dpi),
              skipfirstnumber=cmd_options.skipfirstnumber, unpairedcircle=cmd_options.unpairedcircle,
              figwidth=int(cmd_options.width), figheight=int(cmd_options.height),
              noscale=cmd_options.noscale, bbwidth=float(cmd_options.bbwidth))

        else:
          # Decide whether to show titles
          if cmd_options.hidetoptitle:
            tmp_top=""
          else:
            tmp_top=titletop
          if cmd_options.hidebottomtitle:
            tmp_bottom=""
          else:
            tmp_bottom=titlebottom

          Export.TwoD(base, unit, helix, loop, seq, s, gg.icolors, gg.dsb, material,
            titlebottom=tmp_bottom, titletop=tmp_top,
            jname=cmd_options.jsonout, svgname=cmd_options.svgfile2d,
            dotname=cmd_options.dotfile, pngname=cmd_options.pngfile2d,
            drawbases=cmd_options.drawbases, colorbar=cmd_options.colorbar,
            colorbarspace=cmd_options.colorbarspace, colorbaseprob=cmd_options.colorbaseprob,
            colorstrands=cmd_options.colorstrands, colorbarlabel=cmd_options.legend_label,
            show2d=cmd_options.show2d, colordomains=cmd_options.colordomains,
            labeldomains=cmd_options.labeldomains, domainnames=cmd_options.domainnames,
            drawbasenumbers=cmd_options.drawbasenumbers,
            baseidbar=cmd_options.baseidbar, drawbaseticks=cmd_options.drawbaseticks,
            numberinterval=numberinterval, colorbaseid=cmd_options.colorbaseid,
            pseudo=pseudo, debug=cmd_options.debug, unpairedcircle=cmd_options.unpairedcircle,
            domainkey=cmd_options.domainkey, dpi=int(cmd_options.dpi),
            skipfirstnumber=cmd_options.skipfirstnumber,
            figwidth=int(cmd_options.width), figheight=int(cmd_options.height),
            noscale=cmd_options.noscale, bbwidth=float(cmd_options.bbwidth))

      toc(starttime)
      if pylab and cmd_options.show2d and not cmd_options.show25d:
        pylab.show()
      if not (cmd_options.show3d or cmd_options.tachyon3d or cmd_options.json3d \
        or cmd_options.show25d or cmd_options.svgfile25d or cmd_options.pngfile25d \
        or cmd_options.json25d or cmd_options.jsonin25d or cmd_options.pov3d):
        toc(starttime_begin)
        #~ sys.stderr.write("status:DONE\n")
        continue

      #-----------------------------------------------------------------------------
      starttime=time.time()

      print('\nPerform Leaf to Root Aesthetic Rotations')

      #
      # need pointers from each loopside to first base to store coordinates: done
      # need pointers from coil to loopside to retreive base coordinates: done
      # calculate base coordinates base.x: done
      # for each coil take internal spline points from stored base coordinates: done
      # make linked list from coil ends to helix strand ends and back: done
      #
      # need to fix drawhelix so that can draw a helix with a single base pair
      #

      #
      # introduce rotations so that the loop following a helix is in the plane
      # rotated to match the end of the helix
      #
      # actually, decided that for visual clarity it's desirable to keep all the
      # loops in the plane. So, for the loop exiting each helix, rotate loop on
      # helix axis either zero or pi radians depending on exit spin of helix
      #
      # implementation: start with iloop = nloop (the external loop -- with
      # rotation reference angle = 0) so all helixes exiting that loop have
      # refangle = 0.
      #
      # for each loop, loop over exiting helices making a stack of
      # loops at the end of these helices, and assigning accumlated
      # refangles and reforigins to these loops
      #
      # continue until nothing in the stack
      #
      # then loop over coils and helices and use pointers to loops to do transformations
      # based on refangles and reforigins
      # for coils, refloop is coil(icoil).loopnum
      # for helices, refloop is helix(ihelix).prevloop
      #
      # it is critical to do these transformations in the reverse order that the
      # transformations were encountered by the above process: this way we work
      # from the leaf of the tree to the trunk, so that the reference origin and
      # axis is valid for each transformation at the time it is applied
      #

      globalangle = 0.0                    # can rotate around z axis to change orientation (in radians)
      loop[0].ref.axis[:,0]   = numpy.array([0, 0, 0]) # rotations work for z axis but not for arbitrary axes
      loop[0].ref.num    = 0         # number of reference rotations

      nstack = 1
      loopstack=[0]
      while nstack > 0:
        prevloop = loopstack[0]         # always work from the front of the stack
        loopstack = loopstack[1:nstack] # shift stack entries
        nstack = nstack - 1

        nside = loop[prevloop].nside

        if prevloop == 0:
            maxside = nside
        else:
            maxside = nside-1

        for iside in range (maxside):            # loop over existing helices
          if loop[prevloop].nick[1,iside] != -3: # don't process nicks
            ihelix = loop[prevloop].tohelix[iside]
            if ihelix==-1: continue
            nextloop = helix[ihelix].nextloop

            #loop[nextloop].ref = loop[prevloop].ref #inherit rotation information from next higher loop

            loop[nextloop].ref.num=loop[prevloop].ref.num
            nrefnext = loop[nextloop].ref.num+1

            ex=nrefnext-loop[prevloop].ref.angle.size+1
            #~ if ex>0:
              #~ loop[nextloop].ref.angle=numpy.hstack([loop[prevloop].ref.angle,numpy.zeros(ex)])
              #~ loop[nextloop].ref.origin=numpy.hstack([loop[prevloop].ref.origin, numpy.zeros([3,ex])])
              #~ loop[nextloop].ref.axis=numpy.hstack([loop[prevloop].ref.axis, numpy.zeros([3,ex])])

            #~ else:
            loop[nextloop].ref.angle=loop[prevloop].ref.angle.copy()
            loop[nextloop].ref.origin=loop[prevloop].ref.origin.copy()
            loop[nextloop].ref.axis=loop[prevloop].ref.axis.copy()


            # add new rotation from intervening helix: angle= thc(2)-thc[1]
            # origin = xc(:,2) axis = nc(:,2)
            a=helix[ihelix].npair
            b=helix[ihelix].nc
            c=helix[ihelix].xc

            xc2 =  gg.dzb*(helix[ihelix].npair - 1)*helix[ihelix].nc + \
              helix[ihelix].xc
            binangle = gg.dthb*(helix[ihelix].npair - 1) \
                     + helix[ihelix].thc*math.pi/180. - (math.pi - gg.strutangle)
                     # second and third terms are for nicked helices
           # because need to know rotation since last planar entry point to
           # a helix so pass on rotation infor acquired from previous
           # helices (the third term just cancels out the initial entry
           # angle of first helix since the binning below was developed
           # based on the first term alone)
           #
           # there was a bug in binangle in nudraw6 because helix.thc was
           # in degrees instead of radians so binning was bogus

            if planar:
              modangle = binangle%(2*math.pi)

              if stacknickedhelices and loop[nextloop].nickedhelix:
                   binangle = 0.0 # no flipping if in a nick (just use rotation defined earlier and stored in helix.thc
              else:
                if material == 2:  # B-DNA
                  if modangle > math.pi/4. and modangle <= 5*math.pi/4: # bin rotations based on exit angle from helix -- keep the loops in the plane for clarity
                      binangle = math.pi
                  else:
                      binangle = 0.

                elif material == 1: # A-DNA
                  if modangle > math.pi/4. and modangle <= 5*math.pi/4: # bin rotations based on exit angle from helix -- keep the loops in the plane for clarity
                      binangle = math.pi
                  else:
                      binangle = 0.

            loop[nextloop].ref.angle[nrefnext]    = binangle
            loop[nextloop].ref.origin[:,nrefnext] = xc2.copy()
            loop[nextloop].ref.axis[:,nrefnext]   = helix[ihelix].nc.copy() # nc2 = nc1
            loop[nextloop].ref.num                = nrefnext # This is an idex, not a number!!
            nstack = nstack + 1
            loopstack.append(nextloop)

      if smartrot:
        #
        # loop over helices and do transformations on xc1, nc1, thc1
        # it's really important to loop from leaves to trunk of tree so origins and
        # axes remain valid during each step
        #
        for ihelix in range(nhelix):
          irefloop = helix[ihelix].prevloop
          nref = loop[irefloop].ref.num
          for iref in range(nref,-1,-1):  # leaf to trunk
            refangle  = loop[irefloop].ref.angle[iref]
            reforigin = loop[irefloop].ref.origin[:,iref]
            refaxis   = numpy.matrix(loop[irefloop].ref.axis[:,iref]).T

            a0=-refangle
            a1=scipy.cos(-refangle)*I3  \
                   + (1-scipy.cos(-refangle))*refaxis*refaxis.T
            a2=scipy.sin(-refangle) *\
                   numpy.matrix([[0,           refaxis[2], -refaxis[1]], \
                                [-refaxis[2], 0,           refaxis[0]], \
                                [refaxis[1], -refaxis[0],  0]])
            rmat    = a1+a2

            nc1 = helix[ihelix].nc.copy()
            xc1 = helix[ihelix].xc.copy()
            thc1 = helix[ihelix].thc
            nc1 = numpy.dot(rmat,nc1)
            xc1 = numpy.dot(rmat,(xc1 - reforigin)) + reforigin

            thc1 = thc1 + refangle*180/math.pi #spin helix on axis so exit angle stays same relative to that of entering helix
            helix[ihelix].nc = numpy.array(nc1).flatten()
            helix[ihelix].xc = numpy.array(xc1).flatten()
            helix[ihelix].thc  = thc1

          helix[ihelix].thc  = helix[ihelix].thc \
            + 90. + 180./math.pi*(scipy.arctan2(helix[ihelix].nc[1],helix[ihelix].nc[0]) - globalangle)
        # rotate so that projection from canonical
        # u_z orientation creates helix that joins the loop
        # with the same rotation regardless of
        # helix axis -- this has nothing to do with the coils since they start
        # out defined in the x-y plane (so no project issue from u_z)

        #
        # loop over coils and do transformations on xsplinespts: loop from leaves
        # to trunk so origins and axes remain valid as you work
        #
        #~print "*** rotate"
        for icoil in range(ncoil):
          irefloop = coil[icoil].loopnum
          nref = loop[irefloop].ref.num
          nsplinepts = coil[icoil].nsplinepts
          for iref in range(nref,-1,-1): #leaf to trunk
            refangle  = loop[irefloop].ref.angle[iref]
            reforigin = loop[irefloop].ref.origin[:,iref]
            refaxis   = numpy.array(loop[irefloop].ref.axis[:,iref])
            #~ print"coil[%d].refaxis=%s"%(icoil, refaxis)
            refaxisMatrix=numpy.matrix(refaxis)
            a= scipy.cos(-refangle)*I3
            b=(1.-scipy.cos(-refangle))*(refaxisMatrix.T*refaxis)
            c=scipy.sin(-refangle)* numpy.array([[0.,        refaxis[2], -refaxis[1]], \
                                                      [-refaxis[2], 0.,          refaxis[0]], \
                                                       [refaxis[1], -refaxis[0], 0.]])

            rmat    = a+b+c
            #~ print "nsplinepts=%d"%(nsplinepts)
            for ipt in range(nsplinepts):
              xsplinepts = coil[icoil].xsplinepts[:,ipt]
              xsplinepts = numpy.dot(rmat,(xsplinepts - reforigin)) + reforigin
              #~ print "xsplinepts=%s"%(xsplinepts)
              #~ if loop[irefloop].tohelix[0]<0 and nsplinepts>1:
                #~ print "moving coil"
                #~ xsplinepts+=refaxis*(gg.dzb/(nsplinepts/3.0))
              coil[icoil].xsplinepts[:,ipt] = xsplinepts

            center=coil[icoil].center
            coil[icoil].center=numpy.array((numpy.dot(rmat,(center - reforigin)) + reforigin)).flatten()
            #~ if loop[irefloop].tohelix[0]<0 and nsplinepts>1:
              #~ coil[icoil].center+=refaxis*(gg.dzb/(nsplinepts/3.0))


      toc (starttime)
      starttime=time.time()

      #-----------------------------------------------------------------------------
      if cmd_options.jsonin25d and json_data:
        # Update base coordinates from JSON file
        update_prob = not cmd_options.probfile and not cmd_options.simpleprobfile
        JsonLayout.updateBases25d(json_data, base, unit, update_prob)

        # Update loop coordinates from JSON file
        JsonLayout.updateHelix25d(json_data, loop, base, helix, gg)

        JsonLayout.updateCoils25d(json_data, loop, base, coil)
        JsonLayout.updateCoils25dSweep(json_data, loop, base, coil)
        #~ continue

      #-----------------------------------------------------------------------------
      print('\nRefine Helices, Render helix bases')

      # loop over helices
      # draw each helix starting from xc,nc
      # save output for strand locations and orientations at ends of helix
      # helix vector is into helix at xc1 and out of helix at xc2

      mycoil=Coil.DnaCoil()

      for ihelix in range(nhelix):
          npair   = helix[ihelix].npair
          xc1     = helix[ihelix].xc  # start of helix axis
          nc1     = helix[ihelix].nc  # direction of helix axis
          thc1    = helix[ihelix].thc # rotation around helix axis in degrees

          dzc1    = 0    # translation along helix axis from start of helix axis
          if colorscheme == 1:
             helixshades =  numpy.array([helix[ihelix].strand[0],helix[ihelix].strand[1]])
          elif colorscheme == 2:
             helixshades = numpy.array([2, 3])

          shades = numpy.array([helixshades[0], helixshades[1], 17])  # colors of chain a, chain b, base pairs
          render = numpy.array([1, 1, 1]).T # render chain a, chain b, base pair struts

          acap = numpy.array([helix[ihelix].cap[0], helix[ihelix].cap[1]])  # cap 5', 3' end of a chain
          bcap = numpy.array([helix[ihelix].cap[3], helix[ihelix].cap[2]])  # cap 5', 3' end of b chain

          adye = numpy.array([0, 0])  # dye 5', 3' end of a chain
          bdye = numpy.array([0, 0])  # dye 5', 3' end of b chain
          #pdb.set_trace()
         #xc2, nc2, thc2, x1a, n1a, x2a, n2a, x1b, n1b, x2b, n2b, nhtot, xc_a, xc_b
          xc2, nc2, thc2, x1a, n1a, x2a, n2a, x1b, n1b, x2b, n2b, nhtot, xc_a, xc_b= \
            Utils3D.drawhelix(base, helix[ihelix], npair, xc1, nc1, thc1, dzc1, shades, render,
            acap, bcap, adye, bdye, mycoil, gg, helpers=False)
          helix[ihelix].xc = xc2
          helix[ihelix].nc = nc2
          helix[ihelix].points=[xc_a, xc_b]
          #~ print "helix[%d].nc=%s, tocoil=\n%s"%(ihelix, nc2, helix[ihelix].tocoil)


          helix[ihelix].thc  = thc2
          helix[ihelix].x = numpy.array([x1a, x2a, x1b, x2b]) # Needs to be transposed, beware!
          helix[ihelix].n = numpy.array([n1a, n2a, n1b, n2b])

      # Add more points to xsplinepts
      #~ for icoil in range(len(coil)):
        #~ if coil[icoil].xsplinepts.shape[1]>0:
          #~ print coil[icoil].xsplinepts.shape, coil[icoil].xsplinepts
          #~ loopcenter=coil[icoil].center
          #~ loopradius=loop[coil[icoil].loopnum].radius
          #~ xr, yr, zr, dxr, dyr, dzr=circlefit(coil[icoil].xsplinepts, loopcenter,
              #~ loopradius, len(coil[icoil].xsplinepts)*3, mycoil, gg.colors)
          #~ xc_r = numpy.array([xr, yr, zr])
          #~ print xc_r.shape, xc_r
          #~ coil[icoil].xsplinepts=xc_r
          #~ coil[icoil].nsplinepts=xc_r.shape[1]



      # define end points and vectors for coils
      # vectors are into the coil
      for icoil in range(len(coil)):
          if coil[icoil].nbase > 1:
            for icoilend in range(2):
              ihelix    = coil[icoil].tohelix[icoilend,0]
              ihelixend = coil[icoil].tohelix[icoilend,1]

              if ihelix >= 0:
                x = helix[ihelix].x[ihelixend,:] # transposed from matlab code
                n = helix[ihelix].n[ihelixend,:]
                bctype = 1 # bc based on first derivatives

              else: # coil ends at a nick
                nsplinepts = coil[icoil].nsplinepts

                if ihelixend == -5:
                  if randomendcoils>0:
                    nchoose = int(scipy.floor(0.5*nsplinepts))
                    x = coil[icoil].xsplinepts[:,nchoose:]

                  else:
                    x = numpy.array([coil[icoil].xsplinepts[:,0]]).T

                  n = numpy.array([0, 0, 0]) # set second derivative to zero
                  bctype = 2 # bc based on second derivative
                elif ihelixend == -3:
                  if randomendcoils>0:
                    nchoose = int(scipy.floor(0.5*nsplinepts))
                    x = numpy.array([coil[icoil].xsplinepts[:,nchoose]]).T

                  else:
                    x = numpy.array([coil[icoil].xsplinepts[:,coil[icoil].nsplinepts-1]]).T

                  n = numpy.array([0, 0, 0]) # set second derivative to zero
                  bctype = 2    # bc based on second derivative
              # helix strand normals at strand ends point out of the strands

              if icoilend==0:
                coil[icoil].x       = x
              else:
                coil[icoil].x=numpy.hstack([coil[icoil].x,x])

              coil[icoil].n[:,icoilend]       = n # normals pointing into coil
              coil[icoil].bctype[icoilend]    = bctype

            coil[icoil].th[:,0] = 0 # starting rotation for bases around coil (I think)

            l=coil[icoil].x.shape[1]
            x1          = coil[icoil].x[:,:l/2]
            x2          = coil[icoil].x[:,l/2:]

            if coil[icoil].nsplinepts > 0:
              xsplinepts = coil[icoil].xsplinepts
              np=True
              if coil[icoil].tohelix[0,1] == -5:
                if randomendcoils>0:
                  xp      = numpy.hstack([x1[0], x2[0]])                  # points to spline between start and finish
                  yp      = numpy.hstack([x1[1], x2[1]])
                  zp      = numpy.hstack([x1[2], x2[2]])
                  np=False
                else:
                  xp      = numpy.hstack([xsplinepts[0,:], x2[0]])        # points to spline between start and finish
                  yp      = numpy.hstack([xsplinepts[1,:], x2[1]])
                  zp      = numpy.hstack([xsplinepts[2,:], x2[2]])

              elif coil[icoil].tohelix[1,1] == -3:
                if randomendcoils>0:
                  xp      = numpy.hstack([x1[0], x2[0]])                  # points to spline between start and finish
                  yp      = numpy.hstack([x1[1], x2[1]])
                  zp      = numpy.hstack([x1[2], x2[2]])
                  np=False

                else:
                  xp      = numpy.hstack([x1[0], xsplinepts[0,:]])        # points to spline between start and finish
                  yp      = numpy.hstack([x1[1], xsplinepts[1,:]])
                  zp      = numpy.hstack([x1[2], xsplinepts[2,:]])

              else:
                  xp      = numpy.hstack([x1[0], xsplinepts[0,:], x2[0]]) # points to spline between start and finish
                  yp      = numpy.hstack([x1[1], xsplinepts[1,:], x2[1]])
                  zp      = numpy.hstack([x1[2], xsplinepts[2,:], x2[2]])

            else:
              xp      = numpy.hstack([x1[0], x2[0]])                      # points to spline between start and finish
              yp      = numpy.hstack([x1[1], x2[1]])
              zp      = numpy.hstack([x1[2], x2[2]])
              np=False
            if np:
              coil[icoil].xsplinepts = numpy.array([xp, yp, zp])

            coil[icoil].nsplinepts = xp.size

      toc(starttime)
      starttime=time.time()

      #-----------------------------------------------------------------------------
      print('\nRefine Coils')

      #
      #   refine coils
      #
      for icoil in range(len(coil)):
        thiscoil=coil[icoil]
        if thiscoil.nbase > 1:
          tohelix=thiscoil.tohelix
          nbases  = coil[icoil].nbase    # need to add one only to length of coils not at end of strand
          x1      = coil[icoil].x[:,0]   # start of coil
          n1      = coil[icoil].n[:,0]   # starting coil vector (into coil)
          th1     = coil[icoil].th[:,0]  # starting rotation around n1
          x2      = coil[icoil].x[:,-1]  # end of coil
          n2      = coil[icoil].n[:,1]   # ending coiled vector (into coil)
          xp      = coil[icoil].xsplinepts[0,:]
          yp      = coil[icoil].xsplinepts[1,:]
          zp      = coil[icoil].xsplinepts[2,:]
          bctype  = coil[icoil].bctype   # 1: first derivative bc, 2: 2nd derivative bc
          if colorscheme == 1:
              coilshade = coil[icoil].strand
          elif colorscheme == 2:
              coilshade = 5

          shades  = coilshade*numpy.ones(nbases+1)
          shades[-1]=17

          render  = numpy.array([1, coil[icoil].nick[0], 1, coil[icoil].nick[1]])  # render coil, first base, middle bases, last base
          ccap    = numpy.array([coil[icoil].nick[0], coil[icoil].nick[1]])        # render cap at 5', 3' ends of coil

          cdye    = numpy.array([0, 0])           # render cap at 5', 3' ends of coil
          direction = 1          # 1: 5' end at x1, 2: 5' end at x2

          if randomendcoils and (coil[icoil].nick[0] or coil[icoil].nick[1]):
              helical = 1        # make coil helical around spline
          else:
              helical = 0


          random = 1             # randomize angle increments in helix and also base pair orientations
                                 # if helical=0, still randomizes base pair rotations around chain
          seed = 10              # random seed (change to get different random chain)
          # arcratio             # adjust spline points to get arcratio \approx 1 for correct chain length
          loopcenter=coil[icoil].center
          loopradius=loop[coil[icoil].loopnum].radius

          use_start=(tohelix[0,0]>=0)
          use_end=(tohelix[1,1]>=0)


          coil[icoil].arcratio, xc = Utils3D.drawcoil(base, coil[icoil], nbases, x1,
            n1, th1, x2, n2, xp, yp, zp, bctype, helical, random, shades, render,ccap, cdye,
            direction, seed,mycoil, loopcenter, loopradius, coil[icoil].tohelix.min(),
            gg, use_start, use_end, helpers=False)
          coil[icoil].points=xc

      # Display strands using computed coil and helix points
      toc(starttime)
      starttime=time.time()
      print('\nRender Backbone')

      splines=[]
      i=0
      helpers=False
      if helpers:
        for icoil in coil:
          if icoil.points.shape[0]>1:
            mycoil.addPolyCylinder(icoil.points, colors=cmap[i], radius=0.5)
            mycoil.addSphere(numpy.array(icoil.points[0]).flatten(), colors=colors["white"], radius=rhc)
            mycoil.addSphere(numpy.array(icoil.points[-1]).flatten(), colors=colors["red"], radius=rhc)
          elif icoil.points.shape[0]==1:
            #~ print "small coil points=%s"%repr(icoil.points[0])
            mycoil.addSphere(numpy.array(icoil.points[0]).flatten(), colors=cmap[i], radius=gg.rhc)
          i=i+1

      #-----------------------------------------------------------------------------

      ihelix=0
      strandsdone={}
      while len(strandsdone)< nstrand:
        points=None
        #~ print "points1:None"
        bindices=None
        #~ strand_coils=[]
        #~ strand_colors=[]
        i=0 # coil index
        j=0 # helix index
        # Now find start of strand
        foundcoil=False
        while i < len(coil):
          for e in [0,1]:
            #~ print "loop[%d].tohelix=\n%s\nstrand=%d\n"%(i,loop[i].tohelix,loop[i].strand[e])
            #~ print "unit[%d]=%s\nprevloop=%d"
            if coil[i].tohelix[e,1]==-5 and not strandsdone.has_key(coil[i].strand):
            #~ if loop[i].tohelix[e,1]==-5 and not strandsdone.has_key(loop[i].strand[e]):
              foundcoil=True
              break
          if foundcoil:
            break
          i=i+1

        #~ print "i=%d"%i

        if not foundcoil:
          strand=0
        else:
          strand=coil[i].strand
        #~ strand_color=len(gg.cmap)-2
        strand_color=strand%len(gg.cmap)

        foundhelix=False
        k=0 # start of helix, see later
        e=0 # coil end: 0=start 1=end
        # Now follow coil, helix pairs
        m=0 # Segment counter, not used
        strandbases=0
        helpers=False
        nhelix=len(helix)

        ordered_coil=[]
        while 1:
          cindex=(ihelix*2)%len(gg.cmap)
          m=m+1
          # Should we start a new strand?
          added=False
          if i<ncoil: # process this coil
            if coil[i].points!=None and points==None:
              points=coil[i].points
              added=True
              bindices=Utils3D.stackBases(i, coil, None, base, unit, cmd_options, gg)

              if helpers:
                mycoil.addPolyCylinder(points, colors=gg.cmap[ihelix], radius=gg.rhc)
              added=True


            elif coil[i].points!=None and coil[i].points.shape[0]>0:
              if helpers:
                mycoil.addPolyCylinder(coil[i].points, colors=gg.cmap[ihelix], radius=gg.rhc)
              if e==0:
                if Utils3D.dist(coil[i].points[0], points[-1])<1e-5:
                  ind=1
                else:
                  ind=0
                newpoints=coil[i].points[ind:]
                coil[i].x[:,0]=points[-1]
              else:
                if Utils3D.dist(coil[i].points[-1], points[-1])<1e-5:
                  ind=-2
                else:
                  ind=-1
                newpoints=coil[i].points[::ind] # Reverse order
                coil[i].x[:,-1]=points[-1]

              points=numpy.vstack([points,newpoints])
              bindices=Utils3D.stackBases(i,coil, bindices, base, unit, cmd_options, gg)

              added=True

            if added:
              ordered_coil.append(coil[i])
              lp=loop[coil[i].loopnum]
              if lp.nick.min()<-1: addnum=1
              else: addnum=0
              coil[i].start_base=lp.sidebase[coil[i].loopside]
              strandbases=strandbases+lp.sidenbase[coil[i].loopside]+addnum

            # Draw begin/end points
            if coil[i].points!=None and coil[i].points.shape[0]>0 and helpers:
              mycoil.addSphere(numpy.array(coil[i].points[0]).flatten(), color=gg.colors["white"], radius=gg.rhc)
              mycoil.addSphere(numpy.array(coil[i].points[-1]).flatten(), color=gg.colors["red"], radius=gg.rhc)

            j=coil[i].tohelix[(e+1)%2,0] # The index of the helix at the end of the coil
            k=coil[i].tohelix[(e+1)%2,1] # The side of the helix 0=start, 1=end
                                   # 2=start of second, 3=end of second
            l=k/2                  # Index of coil 0 or 1

            if k<0:
              strandsdone[coil[i].strand]=coil[i]
              break # We found a nick
            cindex=cindex+1


          if not added and helix[j].added[0] and helix[j].added[1]:
            #~ print "Skipping helix %d"%j
            j=j+1
            continue
          l=k/2                  # Index of coil 0 or 1
          if j<nhelix:
            added=False
            if helpers:
              mycoil.addSphere(helix[j].points[l][0], color=gg.colors["yellow"], radius=gg.rhc)

            if bindices==None:
              bindices=numpy.array([], dtype=int)
              if helpers:
                mycoil.addSphere(helix[j].points[l][-1], color=gg.colors["orange"], radius=gg.rhc)

            if k==0 or k==2:
              newpoints=helix[j].points[l]
              new_bindices=numpy.array(helix[j].bases, dtype=int)[:,l]
            else:
              newpoints=helix[j].points[l][::-1] # reverse order
              new_bindices=numpy.array(helix[j].bases, dtype=int)[:,l][::-1]
            if points!=None:
              if Utils3D.dist(points[-1],newpoints[0])<1e-5:
                ind=1
              else:
                ind=0
              points=numpy.vstack([points,newpoints[ind:]])
              if e==0:
                coil[i].x[:,-1]=newpoints[0]
              else:
                coil[i].x[:,0]=newpoints[0]
              added=True
            else:
              points=newpoints
              added=True

            bindices=numpy.hstack([bindices, new_bindices])

            if added:
              strandbases=strandbases+helix[j].npair
              helix[j].added[k%2]=True

            if k==0 or k==2:
              i=helix[j].tocoil[k+1,0]
              e=helix[j].tocoil[k+1,1]
            else:
              i=helix[j].tocoil[k-1,0]
              e=helix[j].tocoil[k-1,1]

        helpers=False
        # Draw piecewise polyCylinder
        if helpers:
          mycoil.addPolyCylinder(points, colors=gg.cmap[strand_color], radius=0.25)

        # Fit spline through points
        spoints=points.shape[0]
        sp=numpy.zeros(spoints)
        sp[0] = 0.

        xp=numpy.array(points[:,0]).flatten()
        yp=numpy.array(points[:,1]).flatten()
        zp=numpy.array(points[:,2]).flatten()
        for i in range(1,spoints):
            sp[i] = sp[i-1] + numpy.sqrt((xp[i]-xp[i-1])**2 + (yp[i]-yp[i-1])**2 + \
                    (zp[i]-zp[i-1])**2)

        xbc1 = xbc2 = numpy.zeros(3)

        ss = numpy.linspace(0.,1.,num=spoints*4)
        xs,ys,zs, dxs, dys, dzs=spl("",sp,[xp, yp, zp], ss ,xbc1, xbc2, 0.,0., order=4, stiff=20.0)

        #~ except TypeError, v:
          #~ xs, ys, zs, dxs, dys, dzs=polyfit([xp, yp, zp], order=2, n=nhtot[0])
        #~ xs,ys,zs, dxs, dys, dzs=spl("",sp,[xp, yp, zp], ss,xbc1,xbc2,0.,0., order=2, stiff=10.0)

        points2=numpy.array([xs,ys,zs]).T
        ind_per_base=points2.shape[0]/strandbases
        if cmd_options.colordomains:
          points2_colors,texinds=Utils3D.findPoints(points2, ind_per_base, bindices, base, gg,
            debug=cmd_options.debug)
          #Utils3D.findPointsInt(points2, bindices, base, gg)
          start_color=points2_colors[texinds[0]]
          end_color=points2_colors[texinds[-1]]
        else:
          points2_colors=start_color=end_color=gg.cmap[strand_color]
          texinds=None

        points2_deriv=numpy.array([dxs,dys,dzs]).T

        # Draw the strand
        if cmd_options.debug: nsides=10
        else: nsides=25

        mycoil.addPolyCylinder(points2, colors=points2_colors,
          radius=gg.rhc, nsides=nsides, texinds=texinds)

        # Draw 5' end
        helix_space = numpy.sqrt(gg.dzb**2 + (gg.rdh * numpy.sin(gg.dthb))**2) - (1.5*gg.rhc + 2.0)
        Utils3D.drawcap(points2[::-1], helix_space/2, gg.rhc, gg.dzb, mycoil,
          start_color, arrow=False)

        # Draw 3' end
        #~ if cmd_options.tachyon3d:
          #~ mycoil.addSphere(numpy.array(points2[-1]).flatten(),
            #~ color=end_color, radius=gg.rhc*1.5)
        #~ else:
        Utils3D.drawcap(points2, helix_space/2, gg.rhc, gg.dzb, mycoil,
           end_color)
          #~ pass

        toc(starttime)
        starttime=time.time()
        print('\nRender Coil Bases')
        imax=0
        if True:
          for icoil in range(len(ordered_coil)):
            thiscoil=ordered_coil[icoil]
            if thiscoil.nbase > 1 and not thiscoil.drawn_bases:
              tohelix=thiscoil.tohelix
              nbases  = thiscoil.nbase    # need to add one only to length of coils not at end of strand
              x1      = thiscoil.x[:,0]   # start of coil
              n1      = thiscoil.n[:,0]   # starting coil vector (into coil)
              th1     = thiscoil.th[:,0]  # starting rotation around n1
              x2      = thiscoil.x[:,-1]   # end of coil
              n2      = thiscoil.n[:,1]   # ending coiled vector (into coil)
              xp      = thiscoil.xsplinepts[0,:]
              yp      = thiscoil.xsplinepts[1,:]
              zp      = thiscoil.xsplinepts[2,:]
              bctype  = thiscoil.bctype   # 1: first derivative bc, 2: 2nd derivative bc
              if colorscheme == 1:
                  coilshade = thiscoil.strand
              elif colorscheme == 2:
                  coilshade = 5

              shades  = coilshade*numpy.ones(nbases+1)
              shades[-1]=17
              #shades  = numpy.array([coilshade*numpy.ones([nbases,1]), 17])      # shades for coil and bases

              render  = numpy.array([1, thiscoil.nick[0], 1, thiscoil.nick[1]])  # render coil, first base, middle bases, last base
              #~ print "thiscoil.nick=%s"%(thiscoil.nick)
              ccap    = numpy.array([thiscoil.nick[0], thiscoil.nick[1]])           # render cap at 5', 3' ends of coil

              cdye    = numpy.array([0, 0])           # render cap at 5', 3' ends of coil
              direction = 1          # 1: 5' end at x1, 2: 5' end at x2

              if randomendcoils and (thiscoil.nick[0] or thiscoil.nick[1]):
                  helical = 1        # make coil helical around spline
              else:
                  helical = 0

              #helical = 0
              random = 1             # randomize angle increments in helix and also base pair orientations
                                      # if helical=0, still randomizes base pair rotations around chain
              seed = 10              # random seed (change to get different random chain)
              # arcratio              # adjust spline points to get arcratio \approx 1 for correct chain length

              if icoil==0:
                si=thiscoil.start_base*ind_per_base
              else:
                si=ordered_coil[icoil-1].end_index + (thiscoil.start_base -
                (ordered_coil[icoil-1].start_base + ordered_coil[icoil-1].nbase))*ind_per_base

              ei=int(si + (nbases+3)*ind_per_base)
              if si>=points2.shape[0]:
                ei=points2.shape[0]
                si=int(ei - (nbases+3)*ind_per_base)


              use_start=(tohelix[0,0]>=0)
              use_end=(tohelix[1,1]>=0)

              helpers=False
              thiscoil.end_index=Utils3D.drawcoil_bases(thiscoil, nbases, th1,
                                  shades, render, mycoil, points2,
                                  points2_deriv, x1, x2, si, ei, gg, use_start=use_start,
                                  use_end=use_end, ipb=ind_per_base, random=True,
                                  helpers=helpers)
              thiscoil.drawn_bases=True
              #~ print "tohelix=\n%s"%(thiscoil.tohelix)

      toc(starttime)

      if cmd_options.tachyon3d:
        mycoil.export_tachyon(cmd_options.tachyon3d)

      if cmd_options.pov3d:
        mycoil.export_povray(cmd_options.pov3d)


      if cmd_options.json3d:
        mycoil.export_json(cmd_options.json3d)

      #~ print "cmd_options.show25d =",cmd_options.show25d
      if cmd_options.json25d or cmd_options.pngfile25d or cmd_options.svgfile25d \
        or cmd_options.show25d:
        if cmd_options.numberinterval in ['auto','0', '-1']: # Auto line numbering
          intv= max(1.0,math.ceil(math.fabs(math.log(float(nbase)/4)))*2)
          numberinterval=intv
        else:
          numberinterval=int(cmd_options.numberinterval)
        #~ print "numberinterval=%d"%numberinterval
        Export.ThreeD(base, unit, helix, loop, coil, seq, s, gg.icolors, gg.dsb,
          material, titlebottom=titlebottom, titletop=titletop,
          jname=cmd_options.json25d, svgname=cmd_options.svgfile25d,
          dotname=cmd_options.dotfile, pngname=cmd_options.pngfile25d,
          drawbases=cmd_options.drawbases, colorbar=cmd_options.colorbar,
          colorbarspace=cmd_options.colorbarspace, colorbaseprob=cmd_options.colorbaseprob,
          colorstrands=cmd_options.colorstrands, show25d=cmd_options.show25d,
          drawbasenumbers=cmd_options.drawbasenumbers, baseidbar=cmd_options.baseidbar,
          numberinterval=numberinterval, colorbaseid=cmd_options.colorbaseid,
          debug=cmd_options.debug)

      toc(starttime)
      if pylab and cmd_options.show25d:
        pylab.show()

      if not cmd_options.show3d:
        print "Total Processing"
        toc(starttime_begin)
        continue

      if cmd_options.show3d:
        mycoil.Display()
  except Exception,v:
    if cmd_options.optionsjson:
      global_count=1
      traceback.print_exc(10, file=sys.stderr)

    sys.stderr.write("Exception caught:\n%s\nParameters:\n%s"%(v, cmd_options))
    sys.stderr.write("status:ERROR\n")
    if cmd_options.optionsjson:
      continue
    else:
      raise
  sys.stderr.write("status:DONE\n\n")
