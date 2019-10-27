# Logical Setup and Layout
# Conrad Steenberg <conrad.steenberg@caltech.edu>
# Nov 20, 2008

import numpy
import scipy
from Classes import *
import Newton

leftpair = '('
nopair = '.'
rightpair = ')'
nick = '+'

#-------------------------------------------------------------------------------
# Set up unit array describing all the structure
def setup_units(s, seq, unit, base, debug=False):
  ibase = 0
  iunit = 0
  istrand = 0
  nsymbol=len(s)
  #newtype=1, newbase=1, newstrand=1, baseInfo=None):
  unit.append(unitClass(-5, ibase, istrand))
  iunit = iunit + 1
  for isymbol in range(0,nsymbol):
      if s[isymbol] == nick:
          unit.append(unitClass(-3, ibase-1, istrand))
          iunit               = iunit + 1
          istrand             = istrand + 1
          unit.append(unitClass(-5, ibase, istrand))
          iunit               = iunit + 1
      else:
          if len(seq)==len(s):
            bc=baseClass(iunit, seq[isymbol])
          else:
            bc=baseClass(iunit)
          unit.append(unitClass(s[isymbol], ibase, istrand, baseInfo=bc))
          ibase               = ibase + 1
          iunit               = iunit + 1
          base.append(bc)
  unit.append(unitClass(-3, ibase-1, istrand))
  nunit = iunit
  nbase = ibase
  nstrand = istrand + 1

  #
  # identify partners for all logical pairs (strand ends or base pairs)
  # used to traverse loops
  #
  istack = 0
  stack=range(nunit)
  iprev = nunit
  for iunit in range (0,nunit):
      if unit[iunit].type     == -3:
          unit[iprev].next    = iunit
          iprev               = iunit
      elif unit[iunit].type == -5:
          unit[iunit].pair   = iprev
          unit[iprev].pair   = iunit
          unit[iprev].next    = iunit
          iprev               = iunit
      elif unit[iunit].type == leftpair:
          istack              = istack + 1
          stack[istack]       = iunit
          unit[iprev].next    = iunit
          iprev               = iunit # previous paired base along strand
      elif unit[iunit].type == rightpair:
          unit[iunit].pair   = stack[istack]
          unit[stack[istack]].pair  = iunit
          unit[iprev].next    = iunit
          iprev               = iunit
          istack              = istack - 1
      u=unit[iunit]
      if debug:
        print "\nu[%d]:\ntype=%s\nbase=%s"%(iunit,u.type, u.base+1)
        print "pair=%d\nnext=%s"%(u.pair, u.next)
        print "prevloop=%d\nprevside=%s"%(u.prevloop, u.prevside)
  return nunit, nbase, nstrand

def testbase(u, fb):
  if u.base==fb and type(u.type)==str:
    return True
  return False

#-------------------------------------------------------------------------------
def findbase(unit, fb, gs):
  units=[]
  i=gs-1
  j=gs
  k=gs+1
  c=0
  while i>=0 or k<len(unit):
    if c==0:
      if testbase(unit[j], fb):
        return j
    if i>=0 and testbase(unit[i], fb):
        return i
    if k<len(unit) and testbase(unit[k], fb):
        return k
    k+=1
    i-=1
    c+=1
  return -1

#-------------------------------------------------------------------------------
# set up ppairs for pseudoknots
def setup_pseudopairs(unit, base, st):
  for pair in st.pairs:
    #~ print "%3d -> %3d"%pair
    ui0=base[pair[0]].iunit
    ui1=base[pair[1]].iunit

    unit[ui0].pair=ui1
    unit[ui1].pair=ui0
  return

#-------------------------------------------------------------------------------
# traverse all loops starting from a left base pair
# (except last exterior loop, which starts from a right base)
#
# save pointer to previous loop number and loop side
# at left base of pair at point of entry to new loop
#
# starting from a left base in a pair (ibase) and traversing a loop in the
# clockwise direction, jside = 1 indexes the first encountered single-stranded
# region and the first encountered stem. For any given loop, jside = nside
# indexes the last single-stranded region and the last stem (of which ibase
# is part of the closing pair)
#
# nicks are treated logically as pairs between the base at the 3' end of one strand
# and the 5' end of the next strand (by the variable base(ibase).next and
# base(ibase).ipair) but there is extra logic in the variable
# base(ibase).nick to keep track of the fact that they actually represent a
# nick
#
# the first loop starts at ibase = 1 with iside = 1 representing the first
# encountered single stranded region and the subsequent helix

def setup_loops(unit, loop, debug=False):
  iunit = 0 # use i index for moving along chain to find left unit in pair to start each loop
  iloop = 0

  # Get jstartMax
  jsideMax=0
  while (iunit == 0) ^ (iloop > 0):
    iu=unit[iunit]
    if debug:
      print "%d: type=%s, base=%d, next=%d, pair=%d"%(iunit, iu.type, iu.base, iu.next, iu.pair)
    if unit[iunit].type == leftpair or iunit == 0: # special case for first exterior loop
      jstart  = iunit # use j index for traversing loop
      jside   = 0
      junit   = jstart
      loop.append(loopClass())
      if debug:
        print "   junit=%d, jstart=%d, jside=%d next=%d"%(junit, jstart, jside, unit[junit].next)

      while (junit == jstart) ^ (jside > 0):
          jnext = unit[junit].next
          jside = jside + 1
          if jside>len(unit):
            raise RuntimeError("Possibly malformed structure, exiting")
          if jside>jsideMax:
            jsideMax=jside
          junit = unit[jnext].pair
          if debug:
            print "   junit=%d, jstart=%d, jside=%d next=%d"%(junit, jstart, jside, unit[junit].next)


      loop[iloop].nside = jside
      if debug:
        print "loop[%d].nside=%d"%(iloop,jside)
      iloop = iloop + 1

    iunit = unit[iunit].next

  for myloop in loop:
    myloop.setSize(jsideMax)

  iunit = 0 # use i index for moving along chain to find left unit in pair to start each loop
  iloop = 0

  while (iunit == 0) ^ (iloop > 0):
      if debug:
        print "\n1.unit[%d].type=%s, base=%d"%(iunit,unit[iunit].type, unit[iunit].base)
      if unit[iunit].type == leftpair or iunit == 0: # special case for first exterior loop
          jstart  = iunit # use j index for traversing loop
          junit   = jstart
          jside   = 0

          while (junit == jstart) ^ (jside > 0):
              jnext                           = unit[junit].next
              if debug:
                print "\n2.unit[%d].type=%s jside=%d jnext=%d"%(jnext,unit[jnext].type,jside,jnext)

              loop[iloop].sideunit[jside]     = junit            # unit right before single stranded region for each side
              loop[iloop].sidebase[jside]     = unit[junit].base # base right before single stranded region for each side
              loop[iloop].sidenbase[jside]    = unit[jnext].base - unit[junit].base - 1

              if unit[junit].type             == -5:
                  loop[iloop].nick[:,jside]   = numpy.array([-5.,0.])
                  if unit[junit].next         == -1:
                    loop[iloop].nick[:,jside]   = numpy.array([0.,-3.])
              elif unit[jnext].type           == -3:
                  loop[iloop].nick[:,jside]   = numpy.array([0.,-3.])

              else:
                  loop[iloop].nick[:,jside] = numpy.zeros(2)
              if debug:
                print "  loop[%d].nick=\n%s"%(iloop, loop[iloop].nick)
                print "  loop[%d].sidenbase=\n%s"%(iloop, loop[iloop].sidenbase)
              loop[iloop].strand[jside] = unit[junit].strand

              if unit[jnext].type  == leftpair:
                  unit[jnext].prevloop = iloop # use left base in pair to store prev loop
                  unit[jnext].prevside = jside # use left base in pair to store prev side

              jside = jside + 1
              if jside>jsideMax:
                jsideMax=jside

              junit = unit[jnext].pair
          if debug:
            print "loop[%d].nick=\n%s"%(iloop,loop[iloop].nick)
            print "loop[%d].sidenbase=\n%s"%(iloop, loop[iloop].sidenbase)
          iloop = iloop + 1

      iunit = unit[iunit].next

  nloop = iloop


  #
  # make linked list of loops and sides using saved pointers
  #
  for iloop in range(nloop):
      #~ print "setting loop[%d].setref(%d)"%(iloop, jsideMax*2)
      loop[iloop].setRef(jsideMax*2)

  for iloop in range(1,nloop):
      if debug:
        sidebase=map(lambda x: unit[x].base, loop[iloop].sideunit)
        print "%d: sideunit=%s sidebase=%s sidenbase=%s"%(iloop, loop[iloop].sideunit, sidebase, loop[iloop].sidenbase)
      iunit = loop[iloop].sideunit[0]
      nside = loop[iloop].nside
      prevloop = unit[iunit].prevloop
      prevside = unit[iunit].prevside
      loop[iloop].toloop[nside-1] = prevloop # link to previous loop
      loop[iloop].toside[nside-1] = prevside # link to previous side
      if debug:
        print "%d: prevloop=%d, prevside=%d"%(iloop, prevloop, prevside)
      loop[prevloop].toloop[prevside] = iloop  # link from previous loop
      loop[prevloop].toside[prevside] = nside-1  # link from previous side

      #
      # make linked list of sides within each loop (just for convenience when
      # traversing from side = 1 to side = nside or back)
      #

      for iside in range(nside):
          prevside = iside - 1
          if prevside < 0:
              prevside = nside-1

          nextside = iside + 1
          if nextside >= nside:
              nextside = 0

          loop[iloop].prevside[iside] = prevside
          loop[iloop].nextside[iside] = nextside
           #[iloop iside nside0 prevside nextside]
  for iloop in xrange(nloop):
      if debug:
        sidebase=map(lambda x: unit[x].base, loop[iloop].sideunit)
        print "%d: sideunit=%s sidebase=%s"%(iloop, loop[iloop].sideunit, sidebase)

  return nloop


