
#
# 2D file exporter for json, svg, png and dot
#
# Conrad Steenberg <conrad.steenberg@caltech.edu>
# 12 Oct 2008

import numpy
import scipy
import math
import Utils
import json2 as json
import pdb
import matplotlib

from matplotlib import pyplot
import matplotlib as mpl
from matplotlib.ticker import NullLocator
import pylab

import Globals
from Globals import cmap, color_names, colors, icolors, colorid, \
      domainColors, getDomainColors

#-------------------------------------------------------------------------------
# letter offsets in base radius units
letter_xoffset={
  'A':  -0.01,
  'C': -0.05,
  'G': -0.05,
  'U': -0.07,
  'T': 0.01
}

letter_yoffset={
  'A':  0.0,
  'C':  -0.025,
  'G':  -0.0,
  'U':  -0.04,
  'T':  -0.07
}

newver='0.99'
realver=matplotlib.__version__

try:
  arc = matplotlib.patches.Arc
  ssc=arc.set_solid_capstyle
  isNewMatplotlib=True
except:
  isNewMatplotlib=False



#-------------------------------------------------------------------------------
def get_angle(center, point):
  """ Returns the angle between the line point -> center and the line y==0
  """

  # cos TH=(point.x-center.x)
  #           point-center
  l=point-center
  th = scipy.arccos((point[0] - center[0])/ numpy.sqrt(numpy.vdot(l,l)))

  if th<0:
    th=scipy.pi+th
  if center[1]>point[1]:
    th=-th

  th=th*180/scipy.pi
  return th

#-------------------------------------------------------------------------------
def get_angle2(center, point):
  """ Returns the angle between the line point -> center and the line y==0
  """

  # cos TH=(point.x-center.x)
  #           point-center
  l=point-center
  th = scipy.arccos((point[0] - center[0])/ numpy.sqrt(numpy.vdot(l,l)))

  if th<0:
    th=scipy.pi+th
  if center[1]>point[1]:
    th=-th

  th=th*180/scipy.pi

  if th<0.0:
    th=360.0+th

  return th


#-------------------------------------------------------------------------------
def arrow_head(P1, P2, s):
  """Calculates points needed to draw arrow head, given a vector
  s  : segment length
  rdh: base radius
  """
  rdh=s/3.
  v1=P2-P1
  nv1=v1/numpy.linalg.norm(v1) #Normalized vector

  v2=P2+nv1*rdh*1.85
  #~ v2=P2+nv1*rdh*1.15
  pv=numpy.array([-nv1[1],nv1[0]])
  a1=P2+rdh*pv*1.0
  a2=P2-rdh*pv*1.0

  b1=a1+rdh*nv1*0.85
  b2=a2+rdh*nv1*0.85
  v3=P2+nv1*rdh*1.0

  #~ return a1,a2,v2
  return b1, v2, b2, a2, v3, a1

#-------------------------------------------------------------------------------
def arrow_coil(unit, iunit, center, r, th, dth, s, ax, debug=None):
  """Get vector for arrow at the end of a coil
  """
  if debug:
    print "----- Arrow"
    print "th=%f, dth=%f"%(th,dth)

  th2=scipy.arctan(s/r)
  r2=s/scipy.sin(th2)
  if debug: print "r=%f r2=%f"%(r,r2)
  th3=th-dth/math.fabs(dth)*th2
  if debug:  print "th3=%f"%th3
  P1=unit[iunit].baseInfo.x[:2]
  P2=center[:2]+numpy.array([r2*scipy.cos(th3), r2*scipy.sin(th3)])

  return P1,P2

#-------------------------------------------------------------------------------
def arrow_helix(pts, s):
  """Get vector for arrow at the end of a helix
  """
  P2 = pts[-1]
  v1 = P2 - pts[0]
  nv1 = v1/numpy.linalg.norm(v1) # Normalized vector

  v2 = P2 + nv1*s

  return numpy.array([P2, v2])[:,:2]

#-------------------------------------------------------------------------------
def arrow_helix2(pts, s):
  """Get vector for arrow at the end of a helix
  """
  P2 = pts[-1]
  v1 = P2 - pts[0]
  nv1 = v1/numpy.linalg.norm(v1) # Normalized vector

  v2 = pts[-2] + nv1*s

  return pts[-2], v2

#-------------------------------------------------------------------------------
def arrow_helix_perp(pts, s):
  """Get vector for arrow at the end of a helix
  """
  P2 = pts[-1]
  v1 = P2 - pts[0]

  nv1 = v1/numpy.linalg.norm(v1) # Normalized vector
  nv2=nv1[::-1]
  nv2[0]=-nv2[0]

  nv2=nv1[::-1]*-1 # Quick 'n Dirty rotate by 90deg

  v2 = P2 + nv2*s

  return P2, v2

#-------------------------------------------------------------------------------
def dist(v1,v2):
  return numpy.sqrt((v2[0]-v1[0])**2 + (v2[1]-v1[1])**2 + \
              (v2[2]-v1[2])**2)

#-------------------------------------------------------------------------------
def draw_arrow(arp, strand, col, dsb_in, bbonewidth, strutwidth, maxy_pad, miny_pad, ax):

  a1,a2,a3,a4,a5,a6=arrow_head(arp[0],arp[1], dsb_in)
  A1=numpy.vstack([arp[0], a5])
  A2=numpy.vstack([a1,a2,a3,a4,a5,a6])
  arrow_desc={"strand"   : strand,
              "color"     : col,
              "shaftx"    : map(float,-A1[:,0]+maxy_pad+miny_pad),
              "shafty"    : map(float,-A1[:,1]+maxy_pad+miny_pad),
              "arrowheadx": map(float,-A2[:,0]+maxy_pad+miny_pad),
              "arrowheady": map(float,-A2[:,1]+maxy_pad+miny_pad)
            }
  if ax:
    poly1=matplotlib.lines.Line2D([arp[0][0], a5[0]], [arp[0][1], a5[1]],
      lw=bbonewidth, color=col, solid_capstyle='round', zorder=1)

    ax.add_artist(poly1)

    poly2=matplotlib.patches.Polygon(A2, lw=strutwidth/2, fill=True,
      fc=col, color=col)
    ax.add_patch(poly2)

  return arrow_desc

#-------------------------------------------------------------------------------
def get_orient_angle(x1, x2, radius, center):

  ang1=get_angle2(center, x1)
  ang2=get_angle2(center, x2)

  # Determine sense/direction of loop
  if ang1>ang2:
    if ang1-ang2>180.0:
      orient=1
      angle=360.-(ang1-ang2)
    else:
      orient=-1
      angle=ang2-ang1
  else:
    if ang2-ang1<180.0:
      orient=1
      angle=ang2-ang1
    else:
      orient=-1
      angle=ang1-(360.- ang2)
  return orient, angle, ang1, ang2

#-------------------------------------------------------------------------------
def is_nicked_helix(l, j, unit):
  if l.sidenbase[j]==0:
    u=l.sideunit[j]
    nbases=filter(lambda x: x>-1, l.sidenbase)

    if unit[u].pair>=0 and unit[u+1].pair>=0 and len(nbases)<2:
      return True
  return False

#-------------------------------------------------------------------------------
# Draw a possibly 2-colored line segment denoting the domain
def drawHelixLines(base, baseInfo, hnb, ax, lw):
  domains = []
  segments = []
  for i in xrange(hnb-1):
    pos=[base[baseInfo[i]].x]
    if base[baseInfo[i]].domain!=base[baseInfo[i+1]].domain:
      pos.append((base[baseInfo[i]].x+base[baseInfo[i+1]].x)/2)
    pos.append(base[baseInfo[i+1]].x)

    for j in range(len(pos)-1):
      domains.append(base[baseInfo[i+j]].domain)
      segments.append((pos[j], pos[j+1]))

  final_segs = []
  final_doms = []

  dstart = domains[0]
  pstart = segments[0]
  istart = 0

  for i in xrange(1, len(domains)):
    dseg = segments[i][1] - segments[i][0]
    dprev = segments[i-1][1] - segments[i-1][0]

    diff = dseg - dprev

    if abs(diff[0]) > 1e-5 or abs(diff[1]) > 1e-5 or domains[i] != domains[i-1]:
      final_doms.append(domains[i-1])
      final_segs.append((pstart[0], segments[i-1][1]))
      pstart = segments[i]
      dstart = domains[i]
      istart = i

  final_doms.append(domains[-1])
  final_segs.append((pstart[0], segments[-1][1]))

  for dom, seg in zip(final_doms, final_segs):
    col="#%s"%getDomainColors(dom)
    line=matplotlib.lines.Line2D([seg[0][0], seg[1][0]],
            [seg[0][1], seg[1][1]],
            lw=lw, zorder=1, solid_capstyle='butt', color=col)
    ax.add_artist(line)

  col="#%s"%getDomainColors(final_doms[0])
  startpoint = (final_segs[0][1] + final_segs[0][0]) / 2
  endpoint = final_segs[0][0]
  startline=matplotlib.lines.Line2D([startpoint[0], endpoint[0]],
          [startpoint[1], endpoint[1]], lw=lw, zorder=1, solid_capstyle='round',
          color=col)

  ax.add_artist(startline)


  col="#%s"%getDomainColors(final_doms[-1])
  startpoint = (final_segs[-1][1] + final_segs[-1][0]) / 2
  endpoint = final_segs[-1][1]
  endline=matplotlib.lines.Line2D([startpoint[0], endpoint[0]],
          [startpoint[1], endpoint[1]], lw=lw, zorder=1, solid_capstyle='round',
          color=col)

  ax.add_artist(endline)

def drawUnPairedDomains(unit, base, ax, lw):
  # Construct the bases array of arrays
  strands=[]
  for u in unit:
    if u.type == - 5:
      strands.append([])
    elif u.type != - 3:
      strands[-1].append(u.base)

  for baseInfo in strands:
    drawHelixLines(base, baseInfo, len(baseInfo), ax, lw)


#-------------------------------------------------------------------------------
# Draw a possibly 2-colored arc segment denoting the domain
def  drawSegment(base, unit, l, nbasej, segment, lw, ax, debug=False, drawArc=True):
  begin_ind=segment[0]
  end_ind=segment[1]

  if debug:
    print "drawSegment base %d -> %d"%(begin_ind, end_ind)
  pos=[base[begin_ind].x]
  if base[begin_ind].domain!=base[end_ind].domain: # add midpoint
    if debug: print "adding midpoint"
    pos.append((base[begin_ind].x+base[end_ind].x)/2)
  pos.append(base[end_ind].x)

  for i in range(len(pos)-1):

    ang1=get_angle(l.center, pos[i])
    ang2=get_angle(l.center, pos[i+1])

    # Correct winding direction for arc call
    if ang1<ang2:
      begin_ang=ang1
      end_ang=ang2
    else:
      begin_ang=ang2
      end_ang=ang1

    col="#%s"%getDomainColors(base[segment[i]].domain)
    if not is_nicked_helix(l, nbasej, unit) and drawArc:
      if debug: print "Drawing arc from %f to %f, nbasej=%d"%(ang2, ang1, nbasej)
      if isNewMatplotlib:
        arc = matplotlib.patches.Arc((float(l.center[0]), float(l.center[1])),
          float(l.radius*2), float(l.radius*2), 0., float(ang2.real),
          float(ang1.real), lw=lw, solid_capstyle="butt",
          ec=col, fill=False, zorder=1)
      else:
        arc = matplotlib.patches.Arc((float(l.center[0]), float(l.center[1])),
          float(l.radius*2), float(l.radius*2), 0., float(ang2.real),
          float(ang1.real), lw=lw, ec=col, fill=False, zorder=1)
      ax.add_patch(arc)
    else:
      line=matplotlib.lines.Line2D([pos[i][0],pos[i+1][0]],
          [pos[i][1],pos[i+1][1]], lw=lw,
          color=col, zorder=1)
      ax.add_artist(line)

#-------------------------------------------------------------------------------
# bbox=[[x[0], x[1]], [y[0], y[1]]]
def evenlySpaceBases(unit, base, bbox, dXin = 0.0, debug = False):
  # Construct the bases array of arrays
  strands=[]
  for u in unit:
    if u.type == - 5:
      strands.append([])
    elif u.type != - 3:
      strands[-1].append(u.base)

  nDivs = len(strands)+1
  dY = (bbox[1][1] - bbox[1][0])/nDivs


  # Calculate dX from the longest strand
  maxLen = len(strands[0])
  for si in range(1, nDivs - 1):
    ul=len(strands[si])
    if ul > maxLen: maxLen = ul

  boxWidth = bbox[0][1] - bbox[0][0]
  if dXin ==0:
    dX = boxWidth/(maxLen + 2)
  else:
    dX = dXin

  for si in range(0, nDivs - 1):
    y = (si+1) * dY + bbox[1][0]
    xi = 0
    xOffset = bbox[0][0] + (boxWidth - dX*len(strands[si]))/2

    for bi in strands[si]:
      base[bi].x = numpy.array([dX * xi + xOffset, y, 0])
      xi+=1

  if debug:
    for bi in range(len(base)):
      print "base[%s].x = %s"%(bi, base[bi].x)

# end

#-------------------------------------------------------------------------------
def interpDomainIndices(bases, domainlabels):
  # Move domain labels to center of domains
  domainlabels.append([len(bases),"A"]) # We'll remove this next
  for i in xrange(len(domainlabels)-1):
    domainlabels[i][0]=(domainlabels[i][0]+domainlabels[i+1][0]-1)/2

  # Domain lookup map
  return dict(domainlabels)

#-------------------------------------------------------------------------------
# Overloaded 2D export function for
# - SVG
# - PNG
# - json
# - dot
def TwoD(base, unit, helix, loop, seq, dotparens, icolors, dsb, material,
          jname="", svgname="",
          pngname="", dotname="", dpi=100., colorstrands=False, drawbases=True,
          show2d=False, colorbaseprob=True, colorbaseid=False, colorbar=True,
          baseidbar=False, colorbarspace=False, titlebottom="", titletop="",
          drawbaseticks=False,  drawbasenumbers=True, numberinterval=2,
          pseudo=False, debug=False, skipfirstnumber=False, domainkey=False,
          colordomains=False, labeldomains=False, domainnames=None,
          figwidth=5.12, figheight=5.12, colorbarlabel="Equilibrium Probability",
          unpairedcircle=False, noscale=False, bbwidth=24):

  global global_cmap, pyplot
  # Set up canvas if SVG or PNG file to be produced
  if debug:
    deltab=0
  else:
    deltab=1

  if len(svgname)>0 or len(pngname)>0 or show2d or jname or dotname:
    matplotlib.rcParams['ytick.labelsize'] = 8
    fig=pyplot.figure(dpi=dpi,figsize=(float(figwidth)/dpi, float(figheight)/dpi))
    xaxis_start=0.075
    if colorbarspace or colorbar or baseidbar:
      xaxis_start=0.0

    if not (titletop or titlebottom or colorbarspace or colorbar or \
            baseidbar or domainkey):
      ax=fig.add_axes([0.0, 0.0, 1.0, 1.0])
    elif not domainkey:
      ax=fig.add_axes([xaxis_start, 0.075, 0.875, 0.875]) # Square area, fix for non colorbar
    else:
      ax=fig.add_axes([xaxis_start, 0.1, 0.85, 0.85])     # Square area, but smaller
    ax.set_frame_on(False)

    # Turn off minor ticking altogether
    ax.set_xticks([])
    ax.set_xticks([],minor=True)
    ax.set_yticks([])
    ax.set_yticks([],minor=True)
  else: ax=False

  # First do bases
  bases=[]
  arrows=[]
  strandmax=-1
  strands=[]
  domainlabels=[]
  processdomains=(colordomains or labeldomains)

  # For creating dot files
  dotlines=[]
  loopdeps={}
  helixdeps={}
  dotloops={}

  # Set up options for how bases will be drawn
  if drawbases and not drawbaseticks:
    baseradius=dsb/5
  else:
    baseradius=dsb/6

  endcaplen=baseradius+dsb/10
  arrowlen=dsb/3
  baseTicklen=Globals.rdh


  # Find mins/maxes for picture
  minx=[100000, 100000]
  maxx=[-100000, -100000]

  for i in xrange(len(unit)):
    u=unit[i]
    if u.baseInfo:
      for j in [0,1]:
        if u.baseInfo.x[j]<minx[j]: minx[j]=u.baseInfo.x[j]
        if u.baseInfo.x[j]>maxx[j]: maxx[j]=u.baseInfo.x[j]

  # Now pad the smaller edge
  dx=[0,0]

  for j in [0,1]:
    dx[j]=maxx[j]-minx[j]

  if len(helix)==0 and not unpairedcircle and not pseudo:
    drawArc=False
  else:
    drawArc=True

  pad=max(dx)/50+baseradius*2+dsb/3
  if not titletop and not titlebottom:
    x=[minx[0]-pad,maxx[0]+pad]
    y=[minx[1]-pad,maxx[1]+pad]
  else:
    x=[minx[0]-pad*1.3,maxx[0]+pad*1.3]
    y=[minx[1]-pad,maxx[1]+pad*1.6]

  if not drawArc:
    x[1]+=dsb/4
  bbox=[[x[0], x[1]], [y[0], y[1]]]
  j=0

  if dx[0]>dx[1]:
    j=1

  ddx=dx[(j+1)%2]-dx[j]
  bbox[j][0]=bbox[j][0]-ddx/2
  bbox[j][1]=bbox[j][1]+ddx/2

  if ax:
    if len(seq)==0:
      drawbases=False

    if titletop:
      pyplot.text((maxx[0]+minx[0])/2, bbox[1][1], titletop,
        family="sans-serif", ha="center", size=11)

    if titlebottom:
      pyplot.text((maxx[0]+minx[0])/2, bbox[1][0]-dx[1]/15, titlebottom,
        family="sans-serif", ha="center", size=11)

    ax.set_xlim(bbox[0])
    ax.set_ylim(bbox[1])

    miny_pad=bbox[1][0] # min y
    maxy_pad=bbox[1][1] # max y

    rec=matplotlib.patches.Rectangle([bbox[0][0], bbox[1][0]],
      bbox[0][1]-bbox[0][0], bbox[1][1]-bbox[1][0], ec="#0000ff", fill=False, lw=0.0)
    ax.add_patch(rec)

    # Set up line widths, dsb is the coil arc length
    black="#000000"

    box_width=rec.get_width()
    box_height=rec.get_height()
    if noscale:
        scalefac = 30
        fig_width = box_width / scalefac
        fig_height = box_height / scalefac
        fig.set_size_inches(fig_width, fig_height)
    else:
        fig_width,fig_height=fig.get_size_inches() # These should be the same as given earlier

        if box_width>box_height:
            scalefac=box_width/fig_width
        else:
            scalefac=box_height/fig_height

    bbonewidth = bbwidth / scalefac
    strutwidth = bbwidth / 2 / scalefac
    numsize=float(max(6.0, bbonewidth*2)) # In points
    numsize=float(min(numsize, 12.0)) # In points
    labelsize=numsize*1.5
    numtmp=numsize/72.72*scalefac
    numoffset_raw=max(baseradius*1.4, dsb/2*scipy.log10(float(len(base)))*scipy.log10(scalefac))
    digits=scipy.ceil(scipy.log10(len(base)))
    numoffset=min(numoffset_raw+numtmp*digits, baseradius+numtmp*digits)
    numoffset=max(numoffset,baseradius*2.25)

    ticklen=numoffset/2.5
    if drawbases:
      ticklen=numoffset/2

    numwidth=numwidthHelix=numwidthCoil=strutwidth/2
    if drawbaseticks:
      numwidthCoil=0.0

    tickwidth=strutwidth*1.5

    if drawbaseticks:
      numoffsetCoil=max(numoffset, baseTicklen+numtmp)
    else:
      numoffsetCoil=numoffset
    numoffsetHelix=numoffset

    labeloffsetHelix=numoffsetHelix*1.3
    labeloffsetCoil=numoffsetCoil*1.3

    # Set up options for how bases will be drawn
    if drawbaseticks:
      clw=0.0
    else:
      clw=strutwidth

    if colorbar or colorbaseprob:
      cmap = mpl.cm.jet
      norm = mpl.colors.Normalize(vmin=0.0, vmax=1.0)
      mappable=mpl.cm.ScalarMappable(norm,cmap)

    if colorbar: # Draw a colorbar
      ax1 = fig.add_axes([0.875, 0.1, 0.0175, 0.8])
      ax1.set_frame_on(False)
      ax1.yaxis.tick_right()
      ax1.set_axisbelow(True)

      ax1.xaxis.set_visible(False)
      locator=matplotlib.ticker.MultipleLocator(0.2)
      ax1.yaxis.set_major_locator(locator)

      cb1 = mpl.colorbar.ColorbarBase(ax1, cmap=cmap, norm=norm,
              orientation='vertical', drawedges=False, ticks=None)
      cb1.outline.set_linewidth(0.0)
      rec2=matplotlib.patches.Rectangle([0, 0], 1,1,
        ec="#ffffff", fill=False, lw=0.0)
      ax1.add_patch(rec2)

      ax1.autoscale_view(tight=True)
      ax1.yaxis.set_ticks([0.0, 0.2, 0.4, 0.6, 0.8, 1.0], minor=False)
      ax1.yaxis.set_ticklabels(["0.0", "0.2", "0.4", "0.6", "0.8", "1.0"],
                                minor=False, size=9)
      ax1.yaxis.set_label_text(colorbarlabel, size=11,
                                rotation=-90.0)


    if domainkey: # Draw a base ID legend
      strucdomain_names={}
      dn_maxlen=0
      dn_maxlen_struct=0

      for bi in xrange(len(base)):
        dn=domainnames[base[bi].domain]
        if (len(dn))> dn_maxlen: dn_maxlen=len(dn)
        strucdomain_names[dn]=base[bi].domain # bi of domain start

      sd={}
      for key in strucdomain_names.keys():
        if (len(key))> dn_maxlen_struct: dn_maxlen_struct=len(key)
        if key[-1]!='*' and not sd.get(key):
          sd[key]=1
        if key[-1]=='*' and not sd.get(key[:-1]):
          sd[key[:-1]]=1

      sd=sd.keys()
      sd.sort()

      strucdn_map = {} # Use this to look up short domain names
      if dn_maxlen>4:
        for dni in xrange(len(sd)):
          domain1 = strucdomain_names.get(sd[dni],None)
          if domain1!=None:
            strucdn_map[domain1]=str(dni+1)
          domain2=strucdomain_names.get(sd[dni]+"*",None)
          if domain2!=None:
            strucdn_map[domain2]=str(dni+1)+"*"
      else:
        for dni in xrange(len(sd)):
          domain1 = strucdomain_names.get(sd[dni],None)
          if domain1!=None:
            strucdn_map[domain1]=sd[dni]
          domain2=strucdomain_names.get(sd[dni]+"*",None)
          if domain2!=None:
            strucdn_map[domain2]=sd[dni]+"*"

      ax3 = fig.add_axes([0.0, 0.00, 1.00, 1.00])
      ax3.set_frame_on(False)
      ax3.yaxis.tick_right()
      ax3.set_axisbelow(True)

      ax3.xaxis.set_visible(False)
      ax3.yaxis.set_visible(False)

      ax3.set_xlim([0, 1])
      ax3.set_ylim([0, 1])

      dmy=0.03
      dmx=0.025*dn_maxlen_struct * dpi/100.0
      rows=2
      cols=3

      if dn_maxlen>4:
        ordered_dnmap={}
        for dni in xrange(len(sd)):
          dn=sd[dni]
          domain=strucdomain_names.get(dn, None)
          if domain==None:
            domain=strucdomain_names.get(dn+'*', None)
          ordered_dnmap[int(domain)]=dn
        ordered_indices=ordered_dnmap.keys()
        ordered_indices.sort()

        for dni in xrange(len(sd)):
          idom = ordered_indices[dni]
          dn=ordered_dnmap[idom]
          x=0.025+(dni/rows)*dmx
          y=0.10-(dni%rows)*dmy
          ax3.text(x, y,
            "%d: %s"%(dni+1, sd[dni]),
            color="#000000", va="center", ha="left", stretch="semi-condensed",
            zorder=3, family="sans-serif", size=9.0)

    else: # Draw a base ID legend
      strucdn_map = {}
      for bi in xrange(len(base)):
        dni=base[bi].domain
        strucdn_map[dni]=domainnames[dni]


    if baseidbar: # Draw a base ID legend
      ax2 = fig.add_axes([0.9, 0.1, 0.35, 0.8])
      ax2.set_frame_on(False)
      ax2.yaxis.tick_right()
      ax2.set_axisbelow(True)

      ax2.xaxis.set_visible(False)
      ax2.yaxis.set_visible(False)
      ax2.autoscale_view(True)
      ax2.set_aspect(1.,adjustable="box")
      i=0
      dy=0.125
      letters=['A', 'C', 'G','T']
      if material==1:
        letters[-1]='U'

      if len(seq)==0 or seq.find('N')>=0:
        letters.append('?')

      for bi in letters:
        fc=colorid[bi]
        lw=0.0
        #~ if bi=='?':
          #~ lw=strutwidth/2
        dot = matplotlib.patches.Circle([0.05,0.9-i*dy], radius=0.035,
          fill=True, zorder=2, ec="#000000", fc=colorid[bi], lw=lw)
        ax2.add_patch(dot)
        pyplot.text(0.15, 0.9-i*dy, bi, color="#000000", va="center", ha="center",
          zorder=3, family="sans-serif", size=11.0)
        i=i+1

  ddx=[dx[0]-dx[1], dx[1]-dx[0]]
  if ddx[0]<0:
    j=0
    i=1
  else:
    j=1
    i=0
  minx[j]=minx[j]-ddx[i]/2
  maxx[j]=maxx[j]+ddx[i]/2


  # Draw bases -----------------------------------------------------------------
  if not drawArc:
    evenlySpaceBases(unit, base, bbox, dsb, debug)

  dotlines.append("digraph G {\n")
  domainids=[]

  for i in xrange(len(unit)):
    u=unit[i]
    if debug:
      print "\nu[%d]:\ntype=%s\nbase=%s"%(i,u.type, u.base)
      print "pair=%d\nnext=%s"%(u.pair, u.next)
      print "prevloop=%d\nprevside=%s"%(u.prevloop, u.prevside)

    if u.baseInfo:
      pair=-1
      if u.pair>=0:
        pair=int(unit[u.pair].base)
      if u.next>=0:
        next=int(unit[u.next].base)
      btype=0
      if i>0 and unit[i-1].type==-5:
        btype=-5

      if i+1<len(unit) and unit[i+1].type==-3:
        btype=-3

      bx=float(u.baseInfo.x[0])
      by=float(maxy_pad-u.baseInfo.x[1]+miny_pad)

      desc= {
          "id":"%s"%u.base,
          "pair":pair,
          "loc":(bx,by),
          "strand":int(u.strand),
          "structype":u.type,
          "type":btype,
      }

      if u.baseInfo.domain!=None:
        desc["domain"]=u.baseInfo.domain
      if processdomains and (btype==-5 or len(domainlabels)==0 or \
                            domainlabels[-1][1]!=u.baseInfo.domain or \
                            domainids[-1][1]!=u.baseInfo.domainid or \
                            u.base==0):
        domainlabels.append([u.base, u.baseInfo.domain])
        domainids.append([u.base, u.baseInfo.domainid])

      bases.append(desc)
      if ax:
        if drawbases:
          fillcolor="#ffffff"
          if len(u.baseInfo.type)>0:
            btext=matplotlib.text.Text(u.baseInfo.x[0]+
              letter_xoffset.get(u.baseInfo.type, 0.0),
              u.baseInfo.x[1]+letter_yoffset.get(u.baseInfo.type, 0.0),
              text=u.baseInfo.type, color="#000000", va="center", ha="center",
              zorder=3, family="sans-serif", size=float(bbonewidth*3))
            ax.add_artist(btext)
        else:
          if colorbaseid:
            fillcolor=colorid['?']
          else:
            fillcolor="#000000"

        ec="#000000"
        if colorbaseprob and u.prob>=0:
          ec=matplotlib.colors.rgb2hex(mappable.to_rgba(u.prob)[:3])
          if not drawbases:
            fillcolor=ec
        elif colorbaseid:
          tempcol=colorid.get(u.baseInfo.type)
          if tempcol:
            if drawbases:
              ec=tempcol
            else:
              ec=fillcolor=tempcol
          else:
            if not drawbases:
              ec=fillcolor
            else:
              ec=colorid['?']

        desc["color"]=ec
        desc["probability"]=u.prob

        if u.type in ['.', '(', ')']:
          if drawbases:
            u.baseInfo.color=ec
          else:
            u.baseInfo.color=fillcolor
          if not drawbaseticks or drawbases:
            dot = matplotlib.patches.Circle(u.baseInfo.x[:2], radius=baseradius,
              fill=True, label=u.type, ec=ec, lw=clw, zorder=2, fc=fillcolor)
            ax.add_patch(dot)

        if pseudo:
          for cu in [u.pair]:
            if cu<0: continue
            pts=numpy.array([[base[u.base].x[0],base[unit[cu].base].x[0]],
              [base[u.base].x[1],base[unit[cu].base].x[1]]]).T
            line=matplotlib.patches.Polygon(pts, zorder=1,
              ec=col, lw=strutwidth, fill=False)
            ax.add_patch(line)

    if u.strand>strandmax:
      strandmax=u.strand
      if colorstrands:
        col=icolors[u.strand]
      else:
        col="000000"
      strands.append({"id":int(u.strand), "color":col})

  # Start loops and helices
  helices=[]
  hloops=[]
  loops=[]
  lastloop=-1
  parentloops=[]
  helices_update=[]

  lasthelix=-1

  i=0
  loopsdone={}
  helicesdone={}
  loopstodo={}
  loopmap={} # Maps our loop ID to JSON loop ID
  lc=0

  # Move domain labels to center of domains, create domain lookup map
  domainmap=interpDomainIndices(bases, domainlabels)
  if debug:
    print "drawArc=",drawArc
    print "domainmap=",domainmap
    print "domainlabels=",domainlabels

  # Iterate through loops
  for li in xrange(len(loop)):
    loopstodo[li]=loop[li]
    if loop[li].coilnum[0]>-1:
      loopmap[li]=lc
      lc=lc+1
  todo=[i] # current work queue
  tlilist=[] # Special case work queue
  newstrand=True
  deltaH=0.0
  strnd=0
  while len(loopsdone) < len(loop):
    if len(todo)==0:
      todo.append(loopstodo.keys()[0]) # Ugh...
      newstrand=True
    i=todo.pop()
    if not loopsdone.has_key(i):
      l=loop[i]
      if debug:
        print "\nL%d:\nsidenbase=%s\nsideunit=%s"%(i,l.sidenbase, l.sideunit)
        print "sidebase=%s\ncoilnum=%s"%(l.sidebase, l.coilnum)
        print "toloop=%s\ntohelix=%s"%(l.toloop, l.tohelix)
        print "center=%s\nradius=%s"%(l.center, l.radius)
        print "strand=%s"%(l.strand)


      if l.coilnum[0]>=0 or l.tohelix[0]>=0: # This is a hairpin/multiloop
        if debug:
          pylab.text(l.center[0],l.center[1],'L%d,%d'%(len(loops),i), size=numsize)

        # Determine bases and segments
        lastloop=lastloop+1
        strand=l.strand

        if l.tohelix[0]>=0: # This loop has children
          parentloops.append(lastloop)

        # Dot file creation
        if not dotloops.has_key(lastloop):
          dotloops[lastloop]={"size": float(l.radius), "deps": [], "esizes":[]}

        # End dotfile creation
        baseInfo=[]
        segments=[]
        units=[]
        nbasej=[]

        for j in xrange(len(l.sidenbase)):
          bnum=l.sidenbase[j]
          try:
            btype=int(unit[l.sideunit[j]].type)
          except ValueError:
            btype=0
          if bnum<0: continue
          unit_id=l.sideunit[j]

          strand_id=l.strand[j]
          base_id=l.sidebase[j]

          oldk=-1
          donespecial=False
          segunits=[]

          for k in xrange(bnum+2):
            addbase=False

            if k==0 or k==bnum+1:
              ind=k%bnum
              try:
                btype=int(unit[l.sideunit[j]+k+ind].type)
              except ValueError:
                btype=0

              if btype==-5 and (l.tohelix[j]>-1 or len(loop)==1 or l.toloop[j]>-1):
                addbase=True

              if  (k>0 and l.toloop[j]==-1 and l.tohelix[j]==-1):
                addbase=True

            if (k>0 and k<bnum+1):
              addbase=True

            if addbase:
                baseInfo.append(int(base_id))
                if not loopdeps.has_key(int(base_id)):
                  loopdeps[base_id]={}
                segunits.append(unit_id)

            if k>0:
              nu=unit[l.sideunit[j]]
              strnd=strand_id
              if processdomains: domain=base[base_id].domain
              else: domain=None
              segments.append([int(base_id-1), int(base_id), int(strnd), domain])
              if debug:
                print "Added segment: %s"%segments[-1]
              nbasej.append(j)
              if not loopdeps.has_key(int(base_id-1)):
                loopdeps[base_id-1]={}

              loopdeps[int(base_id-1)][int(base_id)]=1
            oldk=k
            base_id=base_id+1
            unit_id=unit_id+1

          units.append(segunits)

        # Determine child helices
        childHelices=[]
        for ch in l.tohelix:
          if ch>=0: # Found a child helix
            childHelices.append(int(ch))
        if newstrand:
          parentHelix=-1
        elif tlilist and tlilist[-1][0]==i:
          parentHelix=lasthelix+1
        else:
          parentHelix=lasthelix

        desc= {
            "id": "%d"%lastloop,
            "bases":baseInfo,
            "parentHelix":parentHelix,
            "childHelices": childHelices,
            "radius":float(l.radius),
            "center":(float(l.center[0]), float(maxy_pad-l.center[1] + miny_pad)),
            "segments":segments,
            "loopDirection":"0",
        }
        loops.append(desc)

        # Draw arc
        if debug: print segments
        numseg=len(segments)
        if numseg>0:
          begin_ind=segments[0][0]
          strnd=segments[0][2]
        if debug: print "units=\n%s"%units
        icoil=0


        for s in range(numseg):
          end_ind=segments[s][1]

        if colordomains and not drawArc:
              drawUnPairedDomains(unit, base, ax, bbonewidth)

        for s in range(numseg):
          # Count nuber of bases
          end_ind=segments[s][1]
          domain=segments[s][3]

          if colordomains:
            if drawArc:
              drawSegment(base, unit, l, nbasej[s], segments[s], bbonewidth, ax, debug, drawArc)

          if numseg==0 or s==numseg-1 or segments[s+1][2]!=strnd or \
              segments[s][1]!=segments[s+1][0]:
            ang1=get_angle(l.center, base[begin_ind].x)
            ang2=get_angle(l.center, base[end_ind].x)
            if ang1<ang2:
              begin_ang=ang1
              end_ang=ang2
            else:
              begin_ang=ang2
              end_ang=ang1

            if ax and not colordomains:
              if colorstrands: col="#%s"%icolors[strnd]
              else: col=black
              if not is_nicked_helix(l, nbasej[s], unit) and drawArc:
                if debug: print "Drawing arc1 from %f to %f, nbasej=%d"%(ang2, ang1, nbasej[s])
                if isNewMatplotlib:
                  arc = matplotlib.patches.Arc((float(l.center[0]), float(l.center[1])),
                    float(l.radius*2), float(l.radius*2), 0., float(ang2.real),
                    float(ang1.real), lw=bbonewidth, solid_capstyle="round",
                    ec=col, fill=False, zorder=1)
                else:
                  arc = matplotlib.patches.Arc((float(l.center[0]), float(l.center[1])),
                    float(l.radius*2), float(l.radius*2), 0., float(ang2.real),
                    float(ang1.real), lw=bbonewidth, ec=col, fill=False, zorder=1)
                ax.add_patch(arc)
              else:
                if debug:
                  print "Nicked helix: drawing line from base=%d to base=%d"%(begin_ind, end_ind)
                  print [base[begin_ind].x[0], base[end_ind].x[0]], [base[begin_ind].x[1], base[end_ind].x[1]]
                line=matplotlib.lines.Line2D([base[begin_ind].x[0],base[end_ind].x[0]],
                    [base[begin_ind].x[1],base[end_ind].x[1]], lw=bbonewidth,
                    color=col, zorder=1)
                ax.add_artist(line)

            if s<numseg-1:
              begin_ind=segments[s+1][0]
              strnd=segments[s+1][2]

            # Draw arrow
            drawArrow=False
            if debug:
              print "units[%d]=%s"%(icoil,units[icoil])
            if units[icoil]:
              ii=units[icoil][0]
              for ih in range(len(l.sideunit)):
                if l.sideunit[ih]<0: continue
                if l.tohelix[ih]<0: drawArrow=True # No helix but there is a coil

            if drawArrow:
              if debug: print "drawArrow = True ----->"
              ei=ii+len(units[icoil])-1
              eb=unit[ei].base
              bb=unit[ii].base

              while ei-ii<len(units[icoil]) and unit[ei].base<=eb \
                and ei-eb<len(units[icoil]) and type(unit[ei+1].type)==str:
                ei=ei+1
                if debug: print "ei =",ei

              arp=False
              endcap=False
              # Test for arrows
              offset=numoffsetCoil
              if units[icoil] and unit[ii].type==-3:
                if drawArc:
                  delta_ang=scipy.pi/180*(end_ang-begin_ang)/len(units[icoil])/3
                  arp=arrow_coil(unit, ii, l.center, l.radius,
                    scipy.pi/180*ang1, -delta_ang, dsb/3, ax)
                else:
                  arp=arrow_helix([base[eb].x, base[eb+1].x], dsb/3)
                if colordomains: col="#%s"%getDomainColors(base[segments[s][0]].domain)
              elif units[icoil] and unit[ei+1].type==-3 and (unit[ei].pair<0 or pseudo):
                if drawArc:
                  delta_ang=scipy.pi/180*(end_ang-begin_ang)/len(units[icoil])/3
                  arp=arrow_coil(unit, ei, l.center, l.radius,
                    scipy.pi/180*ang2, delta_ang, dsb/3, ax)
                else:
                  arp=(base[eb+1].x[:2], base[eb+1].x[:2]+numpy.array([dsb/3,0]))

                if colordomains: col="#%s"%getDomainColors(base[segments[s][1]].domain)

              # We have an arrow vector
              if type(arp)!=type(False):
                if debug: print "<----- drawing arrow"
                arrows.append(draw_arrow(arp, unit[units[icoil][0]].strand, col,
                  dsb*1.1, bbonewidth, strutwidth, maxy_pad, miny_pad, ax))

              # Test for endcaps
              if drawbaseticks:
                endcap=False
                if unit[ii].type==-5:
                  if drawArc:
                    delta_ang=scipy.pi/180*(end_ang-begin_ang)/len(units[icoil])/5
                    endcap=arrow_coil(unit, ii+1, l.center, l.radius,
                      scipy.pi/180*ang1, -delta_ang, endcaplen, ax)
                  else:
                    endcap=[base[bb].x, base[bb].x-numpy.array([endcaplen,0,0])]

                  if colordomains: col="#%s"%getDomainColors(base[bb].domain)
                elif units[icoil] and unit[ei+1].type==-5 and (unit[ei].pair<0 or pseudo):
                  if drawArc:
                    delta_ang=scipy.pi/180*(end_ang-begin_ang)/len(units[icoil])/5
                    endcap=arrow_coil(unit, ei, l.center, l.radius,
                      scipy.pi/180*ang2, -delta_ang, endcaplen, ax)
                  else:
                    endcap=arrow_helix([base[eb-1].x, base[eb].x], endcaplen)

                  if colordomains: col="#%s"%getDomainColors(base[eb].domain)

                # We have an endcap vector
                if type(endcap)!=type(False):
                  if debug:
                    print "<----- drawing endcap"
                    print "endcap="
                    print endcap
                  line=matplotlib.lines.Line2D([endcap[0][0],endcap[1][0]],
                    [endcap[0][1],endcap[1][1]], lw=bbonewidth, zorder=1,
                    solid_capstyle='round', color=col)
                  ax.add_artist(line)

            icoil=icoil+1

        # Draw base numbers
        if drawbasenumbers or processdomains:
          for bn in baseInfo:
            if ((bn==0 and skipfirstnumber) or not (bn+deltab)%numberinterval==0) \
              and not domainmap.has_key(bn):
              continue

            if domainmap.has_key(bn) and labeldomains:
              labeltext=strucdn_map[domainmap[bn]]
              fontsize=labelsize
              offset=labeloffsetHelix
              linewidth=0
            else:
              if not drawbasenumbers: continue
              labeltext="%d"%(bn+deltab)
              fontsize=numsize
              offset=numoffsetCoil
              linewidth=numwidthCoil
            if drawArc:
              tarp=arrow_helix([l.center, base[bn].x], offset)
            else:
              tarp=[base[bn].x, base[bn].x + numpy.array([0,offset,0])]

            ax.annotate(labeltext, xy=base[bn].x[:2], zorder=1,
              xytext=(tarp[1][0],tarp[1][1]), fontsize=fontsize, ha="center", va="center",
              arrowprops=dict(ec="#777777", shrinkA=2,shrinkB=2,
              arrowstyle="-", lw=linewidth, fill=False))

        # Draw base ticks
        if drawbaseticks:
          for bn in baseInfo:
            if drawArc:
              tarp=arrow_helix([l.center, base[bn].x], baseTicklen)
            else:
              tarp=[base[bn].x, base[bn].x + numpy.array([0,baseTicklen,0])]
            line=matplotlib.lines.Line2D([tarp[0][0],tarp[1][0]],
              [tarp[0][1],tarp[1][1]], lw=tickwidth, zorder=0,
              solid_capstyle='round', color=base[bn].color)
            if debug:
              print "drawing line",bn, "tickwidth=",tickwidth, "color=",base[bn].color
            ax.add_artist(line)

        loopsdone[i]=l
        del loopstodo[i]
        newstrand=False
        for t in l.toloop[::-1]:
          if t>0: todo.append(t)

# Special case of a helix without a loop
        tlilist_new=[]
        sp_hi=len(l.toloop)-1

        for tli in l.toloop[::-1]:
          if tli<0:
            sp_hi-=1
            continue
          sl=loop[tli] # Special loop
          #~ pdb.set_trace()
          if debug:
            print "i=%i"%i
            print "\n**L%d:\nsidenbase=%s\nsideunit=%s"%(tli,sl.sidenbase, sl.sideunit)
            print "sidebase=%s\ncoilnum=%s"%(sl.sidebase, sl.coilnum)
            print "toloop=%s\ntohelix=%s"%(sl.toloop, sl.tohelix)
            print "center=%s\nradius=%s"%(sl.center, sl.radius)
            print "strand=%s"%(sl.strand)
            print "%d: unit.base=%d sidebase=%d"%(tli,unit[sl.sideunit[0]].base,sl.sidebase[0])

          if sl.coilnum[0]>=0 and unit[sl.sideunit[0]].base==sl.sidebase[0] \
            and l.tohelix[sp_hi]>=0 \
            and not helicesdone.has_key(l.tohelix[sp_hi]):
            dohelix=tlilist_new.append((tli,l.tohelix[sp_hi],i))
          sp_hi-=1

        for ind_i in range(len(tlilist)-1,-1,-1):
          if tlilist and tlilist[ind_i][0]==i:
            tli, hnum, parent_i=tlilist[ind_i]
            del tlilist[ind_i]
            helicesdone[hnum]=True
            sl=loop[tli] # Special loop

            hloops.append((tli,sl))
            ci=tli
            child=loop[ci]

            if ci<0 or child.coilnum[0]>=0 and debug: # We came to the end of a helix
              pass

            lasthelix=lasthelix+1
            if debug:
              print "lasthelix=%d"%lasthelix
              print "parentloops=%s"%(parentloops)
              if parentloops:
                print "parentloops=%s childHelices=%s"%(parentloops,
                loops[parentloops[-1]]["childHelices"])
              else:
                print "parentloops=%s"

            while parentloops and not lasthelix in loops[parentloops[-1]]["childHelices"]:
              parentloops.pop()

            baseInfo=[]
            unitInfo=[]
            oldh=None
            angle=0.0

            if debug:
              print "hloops=%s"%repr(hloops)
            hlength=0.0
            for hindex in xrange(len(hloops)):
              ind,h=hloops[hindex]

              baseInfo.append(int(h.sidebase[0]))
              unitInfo.append(int(h.sideunit[0]))
              bi=h.sidebase[0]

            # One special case
            iside=0
            while 1:
              if iside+1>=len(h.sidebase) or h.sidebase[iside+1]<0: break
              iside=iside+1

            bi2=h.sidebase[iside]+h.sidenbase[iside]+1
            baseInfo.append(int(bi2))
            uis=int(h.sideunit[iside]+h.sidenbase[iside]+1)
            while unit[uis].type not in ['(',')']:
              uis+=1
            unitInfo.append(uis)

            # Draw helix lines to canvas
            if ax:
              strnd=unit[unitInfo[0]].strand
              if colorstrands: col="#%s"%icolors[strnd]
              else: col=black

              pts=numpy.array([[base[baseInfo[0]].x[0],base[baseInfo[1]].x[0]],
                [base[baseInfo[0]].x[1],base[baseInfo[1]].x[1]]]).T
              #~ if not drawbaseticks:
                #~ line=matplotlib.patches.Polygon(pts, zorder=1,
                  #~ ec=col, lw=strutwidth, fill=False)
                #~ ax.add_patch(line)

              # Draw arrow if needed
              pts_vec=[pts[::-1],pts]

              for ui in range(len(unitInfo)):
                if debug:
                  print "unit_index=%d, type=%s, base=%d"%(unitInfo[ui],unit[unitInfo[ui]+1].type, unit[unitInfo[ui]].base)
                if unit[unitInfo[ui]+1].type==-3:
                  arp=arrow_helix_perp(pts_vec[ui], dsb/3)
                  astrnd=unit[unitInfo[ui]].strand
                  if colorstrands: acol="#%s"%icolors[astrnd]
                  else: acol=black

                  if colordomains:
                    acol="#%s"%getDomainColors(base[unit[unitInfo[ui]].base].domain)

                  if debug:
                    print "drawing arrow:\n%s, %d -> %d"%(pts_vec[ui], baseInfo[0], baseInfo[1])

                  if arp:
                    arrows.append(draw_arrow(arp, unit[unitInfo[ui]+1].strand, acol,
                    dsb*1.1, bbonewidth, strutwidth, maxy_pad, miny_pad, ax))

              if drawbasenumbers:
                bnList=baseInfo[:2]
                conds=[not (baseInfo[0]==0 and skipfirstnumber) \
                      and (baseInfo[0]+deltab)%numberinterval==0, \
                      not (baseInfo[1]+deltab==0  and skipfirstnumber) \
                      and (baseInfo[1]+deltab)%numberinterval==0]
                tarp_pts=[pts[::-1], pts]

                for i in range(len(bnList)):
                  bn=bnList[i]
                  if not conds[i] and not domainmap.has_key(bn):
                    continue

                  if False and domainmap.has_key(bn):
                    labeltext=strucdn_map[domainmap[bn]]
                    fontsize=labelsize
                    offset=labeloffsetHelix
                    linewidth=0
                  else:
                    if not drawbasenumbers: continue
                    labeltext="%d"%(bn+deltab)
                    fontsize=numsize
                    offset=numoffsetHelix
                    linewidth=numwidthHelix

                  tarp=arrow_helix(tarp_pts[i], offset)

                  ax.annotate(labeltext, xy=base[bn].x[:2], zorder=0,
                    xytext=tarp[1], fontsize=fontsize, ha="center", va="center",
                    arrowprops=dict(ec="#777777", shrinkA=2, shrinkB=2,
                    arrowstyle="-", lw=linewidth, fill=False))

              # Draw base ticks
              tarp=arrow_helix2(pts[::-1], baseTicklen)
              for i in [0, 1]:
                x=base[baseInfo[i]].x
                if drawbaseticks:
                  ladderColor=base[baseInfo[i]].color
                  ladderWidth=tickwidth
                elif colorstrands:
                  strnd=unit[unitInfo[i]].strand
                  ladderColor="#%s"%icolors[strnd]
                  ladderWidth=tickwidth
                else:
                  ladderColor=col
                  ladderWidth=strutwidth

                line=matplotlib.lines.Line2D([x[0],tarp[1][0]],
                  [x[1],tarp[1][1]], lw=ladderWidth, zorder=0,
                  solid_capstyle='butt', color=ladderColor)
                ax.add_artist(line)

              # Calculate angle
              deltaX=(base[bi2].x[0]-base[bi].x[0])
              deltaY=(base[bi2].x[1]-base[bi].x[1])

              deltaH=scipy.sqrt(deltaX**2 + deltaY**2)

              thP=scipy.arccos(deltaX/deltaH)
              if (deltaY>0):
                  angle=thP+scipy.pi/2
              else:
                angle=2.5*scipy.pi-thP

              if angle>2*scipy.pi:
                angle-=2*scipy.pi
              if debug:
                print "angle3=%f"%(angle*180/scipy.pi)
                print "angle4=%f"%((2*scipy.pi-float(angle))*180/scipy.pi)

            if newstrand:
              parentLoop=-1
            else:
              try:
                parentLoop=loopmap[parent_i]
              except:
                parentLoop=-1

            childLoop=loopmap[ci]
            desc= {
                "id": "%d"%lasthelix,
                "bases":baseInfo,
                "parentLoop":parentLoop,
                "childLoop":childLoop,
                "strand0":int(unit[h.sideunit[0]].strand),
                "strand1":int(unit[h.sideunit[1]].strand),
                "angle":2*scipy.pi-float(angle),
                "length":float(hlength),
                "flips":"-1"
            }

            bl=len(baseInfo)/2

            for bi in xrange(bl):
              if not helixdeps.has_key(baseInfo[bi]):
                helixdeps[baseInfo[bi]]={}
              if not helixdeps.has_key(baseInfo[bi+bl]):
                helixdeps[baseInfo[bi+bl]]={}
              if bi<bl-1:
                helixdeps[baseInfo[bi]][baseInfo[bi+1]]=1
                helixdeps[baseInfo[bi+bl]][baseInfo[bi+bl+1]]=1
              helixdeps[baseInfo[bi]][baseInfo[bi+bl]]=1

            helices.append(desc)
            newstrand=False

            hloops=[]
            if not tli in todo: todo.append(tli)
            loopstodo[tli]=sl
        tlilist=tlilist_new+tlilist # Prepend the new special loops

# End Special case

      elif l.coilnum[0]<0: # We found a helix ----------------------------------
        hloops.append((i,l))
        ci=l.toloop[0]
        child=loop[ci]

        if ci<0 or child.coilnum[0]>=0: # We came to the end of a helix
          lasthelix=lasthelix+1

          while parentloops and not lasthelix in loops[parentloops[-1]]["childHelices"]:
            parentloops.pop()

          baseInfo=[]
          unitInfo=[]
          oldh=None
          angle=0.0

          hlength=0.0
          for hindex in xrange(len(hloops)):
            ind,h=hloops[hindex]
            loopsdone[ind]=h
            del loopstodo[ind]

            baseInfo.append(int(h.sidebase[0]))
            unitInfo.append(int(h.sideunit[0]))
            if hindex>0:
              bi=h.sidebase[0] # index of unit
              oldbi=oldh.sidebase[0]
              deltaH=scipy.sqrt((base[bi].x[0]-base[oldbi].x[0])**2 + \
                      (base[bi].x[1]-base[oldbi].x[1])**2)

              hlength=hlength+deltaH
            oldh=h

          # Two special cases
          baseInfo.append(int(oldh.sidebase[0]+1))
          unitInfo.append(int(oldh.sideunit[0]+1))
          hlength=hlength+deltaH
          ind,h=hloops[0]
          baseInfo.append(int(h.sidebase[1]+1))
          unitInfo.append(int(h.sideunit[1]+1))

          # And now the second series of bases
          for ind,h in hloops:
            baseInfo.append(int(h.sidebase[1]))
            unitInfo.append(int(h.sideunit[1]))
            oldh=h

          # Draw helix lines to canvas
          hnb=len(baseInfo)/2 # helix_nbases
          strnd=unit[unitInfo[0]].strand
          if colorstrands: col="#%s"%icolors[strnd]
          else: col=black

          # Draw arrow
          for s,e in [[0,hnb-1],[-1,hnb]]:
            P1=base[baseInfo[s]].x[:2]
            P2=base[baseInfo[e]].x[:2]

            if colorstrands: col="#%s"%icolors[strnd]
            else: col=black

            if debug:
              print "Drawing arrow/endcap[%d]:"%(baseInfo[s]),unit[unitInfo[s]-1]
              print "                    [%d]:"%(unitInfo[e]+1),unit[unitInfo[e]+1]
              print "                     %d -> %d :"%(baseInfo[s],baseInfo[e])

            if unit[unitInfo[e]+1].type==-3:
              arp=arrow_helix([P1, P2], dsb/3)
              astrnd=unit[unitInfo[e]].strand
              if colorstrands: acol="#%s"%icolors[astrnd]
              else: acol=black

              if colordomains:
                acol="#%s"%getDomainColors(base[baseInfo[e]].domain)
              arrows.append(draw_arrow(arp, astrnd, acol, dsb*1.1, bbonewidth,
                strutwidth, maxy_pad, miny_pad, ax))
            if unit[unitInfo[s]-1].type==-5 and drawbaseticks:
              if unitInfo[s]-1>0 and unit[unitInfo[s]-2].type==-3:
                elen=baseradius/2
              else:
                elen=endcaplen
              arp=arrow_helix([P2, P1], elen)

              astrnd=unit[unitInfo[s]].strand
              if colorstrands: acol="#%s"%icolors[astrnd]
              else: acol=black

              if colordomains:
                acol="#%s"%getDomainColors(base[baseInfo[s]].domain)

              line=matplotlib.lines.Line2D([arp[0][0],arp[1][0]],
                [arp[0][1],arp[1][1]], lw=bbonewidth, zorder=1,
                solid_capstyle='round', color=acol)
              ax.add_artist(line)

          # Set parent and child loops
          childLoop=loopmap[hloops[-1][1].toloop[0]]
          parentLoop=loopmap[hloops[0][1].toloop[1]]


          if ax:
            if not colordomains:
              line=matplotlib.lines.Line2D([base[baseInfo[0]].x[0],base[baseInfo[hnb-1]].x[0]],
                      [base[baseInfo[0]].x[1],base[baseInfo[hnb-1]].x[1]],
                      lw=bbonewidth, zorder=1,
                      solid_capstyle='butt', color=col)
              ax.add_artist(line)

            strnd=unit[unitInfo[hnb]].strand
            if colorstrands: col="#%s"%icolors[strnd]
            else: col=black

            if not colordomains:
              line=matplotlib.lines.Line2D([base[baseInfo[hnb]].x[0],base[baseInfo[hnb*2-1]].x[0]],
                      [base[baseInfo[hnb]].x[1],base[baseInfo[hnb*2-1]].x[1]],
                      lw=bbonewidth, zorder=1,
                      solid_capstyle='butt', color=col)
              ax.add_artist(line)

            if colordomains:
              drawHelixLines(base, baseInfo[:hnb], hnb, ax, bbonewidth)
              drawHelixLines(base, baseInfo[hnb:], hnb, ax, bbonewidth)

            # Draw struts
            for bi in range(hnb):
              pts=numpy.array([[base[baseInfo[bi]].x[0],base[baseInfo[hnb+bi]].x[0]],
                [base[baseInfo[bi]].x[1],base[baseInfo[hnb+bi]].x[1]]])
              pts_old=pts

              #~ if not drawbaseticks:
                #~ line=matplotlib.patches.Polygon(pts.T,
                  #~ ec=col, lw=strutwidth, fill=False, zorder=0)
                #~ ax.add_patch(line)

              # Move number away from loop
              if hnb>=2 and (bi==0 or bi==hnb-1) :
                if bi==0 and parentLoop>=0 and loops[parentLoop]["bases"] and \
                    unit[unitInfo[e]].type!=-3:
                  x1=(base[baseInfo[bi]].x*4 + base[baseInfo[bi+1]].x*1)/5
                  x2=(base[baseInfo[hnb+bi]].x*4 + base[baseInfo[hnb+bi+1]].x*1)/5
                  pts=numpy.array([[x1[0],x2[0]], [x1[1],x2[1]]])
                elif bi==hnb-1 and childLoop>=0:# and len(loops)>childLoop+1 and loops[childLoop]["bases"]:
                  x1=(base[baseInfo[bi]].x*4 + base[baseInfo[bi-1]].x*1)/5
                  x2=(base[baseInfo[hnb+bi]].x*4 + base[baseInfo[hnb+bi-1]].x*1)/5
                  pts=numpy.array([[x1[0],x2[0]], [x1[1],x2[1]]])

              if drawbasenumbers  or processdomains:
                bnList=[baseInfo[bi], baseInfo[hnb+bi]]
                conds=[not (baseInfo[bi]==0 and skipfirstnumber) and \
                  (baseInfo[bi]+deltab)%numberinterval==0, \
                      not (baseInfo[hnb+bi]==0 and skipfirstnumber) \
                and (baseInfo[hnb+bi]+deltab)%numberinterval==0]
                tarp_pts=[pts.T[::-1],pts.T]

                for i in range(len(bnList)):
                  bn=bnList[i]
                  if not conds[i] and not domainmap.has_key(bn):
                    continue

                  if domainmap.has_key(bn) and labeldomains:
                    labeltext=strucdn_map[domainmap[bn]]
                    fontsize=labelsize
                    offset=labeloffsetHelix
                    linewidth=0
                  else:
                    if not drawbasenumbers: continue
                    labeltext="%d"%(bn+deltab)
                    fontsize=numsize
                    offset=numoffsetHelix
                    linewidth=numwidthHelix

                  tarp=arrow_helix(tarp_pts[i], offset)
                  ax.annotate(labeltext, xy=base[bn].x[:2], zorder=0,
                    xytext=tarp[1], fontsize=fontsize, ha="center", va="center",
                    arrowprops=dict(ec="#777777", shrinkA=2, shrinkB=2,
                    arrowstyle="-", lw=linewidth, fill=False))


              # Draw base ticks
              x1=base[baseInfo[bi]].x
              x2=base[baseInfo[hnb+bi]].x
              pts=numpy.array([[x2[0],x1[0]],  [x2[1],x1[1]] ])

              tarp=arrow_helix2(pts.T[::-1], baseTicklen)
              for i in [0, 1]:
                x=base[baseInfo[hnb*i+bi]].x
                if drawbaseticks:
                  ladderColor=base[baseInfo[hnb*i+bi]].color
                  ladderWidth=tickwidth
                elif colorstrands:
                  strnd=unit[unitInfo[hnb*i]].strand
                  ladderColor="#%s"%icolors[strnd]
                  ladderWidth=tickwidth
                else:
                  ladderColor=col
                  ladderWidth=strutwidth

                line=matplotlib.lines.Line2D([x[0],tarp[1][0]],
                  [x[1],tarp[1][1]], lw=ladderWidth, zorder=0,
                  solid_capstyle='butt', color=ladderColor)
                ax.add_artist(line)

          childLoop=loopmap[hloops[-1][1].toloop[0]]
          if childLoop>=0:
            x1=hloops[-1][1].center
            x2=loop[hloops[-1][1].toloop[0]].center

            deltaH=scipy.sqrt((x2[0]-x1[0])**2 + \
              (x2[1]-x1[1])**2)

            angle=scipy.arccos((x2[0]-x1[0])/deltaH)

            if (x2[1]<x1[1]):angle=2*scipy.pi-angle

          desc= {
              "id": "%d"%lasthelix,
              "bases":baseInfo,
              "parentLoop":parentLoop,
              "childLoop":childLoop,
              "strand0":int(unit[h.sideunit[0]].strand),
              "strand1":int(unit[h.sideunit[1]].strand),
              "angle":2*scipy.pi-float(angle),
              "length":float(hlength),
              "flips":"-1"
          }

          bl=len(baseInfo)/2
          if parentLoop>=0:
            dotloops[parentLoop]["deps"].append(childLoop)
            dotloops[parentLoop]["esizes"].append(hlength)

          for bi in xrange(bl):
            if not helixdeps.has_key(baseInfo[bi]):
              helixdeps[baseInfo[bi]]={}
            if not helixdeps.has_key(baseInfo[bi+bl]):
              helixdeps[baseInfo[bi+bl]]={}
            if bi<bl-1:
              helixdeps[baseInfo[bi]][baseInfo[bi+1]]=1
              helixdeps[baseInfo[bi+bl]][baseInfo[bi+bl+1]]=1
            helixdeps[baseInfo[bi]][baseInfo[bi+bl]]=1

          helices.append(desc)
          newstrand=False

          hloops=[]
        if ci>=0:
          todo.append(ci)

  # Create the final dictionary
  properties={"is25d":"0",
              "baseradius":"%f"%float(baseradius),
              "backbonewidth":"%f"%float(bbonewidth/72.72*dpi),
              "strutwidth":"%f"%float(strutwidth/72.72*dpi),
              "material": material,
              }

  if ax:
    td=ax.transData
    bbox=ax.get_window_extent()
    if hasattr(td,"inverted"):
      tdi=td.inverted()
      ll=tdi.transform(numpy.array([bbox.xmin,bbox.ymin]))
      ur=tdi.transform(numpy.array([bbox.xmax,bbox.ymax]))
      properties["ll_fig"]=[float(ll[0]), float(ll[1])]
      properties["ur_fig"]=[float(ur[0]), float(ur[1])]
      properties["ypad"]=[float(miny_pad), float(maxy_pad)]

      tda = td.get_affine().get_matrix()

      properties["affine"]=(float(tda[0,0]), float(tda[1,0]), float(tda[0,1]),
                            float(tda[1,1]), float(tda[0,2]), float(tda[1,2]))
      for bd in bases: # Array of dicts
        bx,by=bd["loc"]
        pl=td.transform(numpy.array([bx,by]))
        bd["pixloc"]=[float(pl[0]),float(pl[1])]
    else:
      ll=td.inverse_xy_tup((bbox.xmin(),bbox.ymin()))
      ur=td.inverse_xy_tup((bbox.xmax(),bbox.ymax()))
      properties["ll_fig"]=[float(ll[0]), float(ll[1])]
      properties["ur_fig"]=[float(ur[0]), float(ur[1])]
      properties["ypad"]=[float(miny_pad), float(maxy_pad)]

      afn=td.as_vec6_val()

      properties["affine"]=(float(afn[0]), float(afn[2]), float(afn[1]),
                            float(afn[3]), float(afn[4]), float(afn[5]))
      for bd in bases: # Array of dicts
        bx,by=bd["loc"]
        pl=td.xy_tup((bx,by))
        bd["pixloc"]=[float(pl[0]),float(pl[1])]

  struc={ "properties"  : properties,
          "bases"       : bases,
          "strands"     : strands,
          "helices"     : helices,
          "loops"       : loops,
          "anchor"      : [{ "base":"0", "type":"2", "elementID":"0" }],
          "arrows"      : arrows,
          "sequence"    : seq,
          "dotparens"   : dotparens
        }

  pyobject={"struc" : struc}

  if len(jname)>0:
    jfile=open(jname,"w")
    jfile.write(json.write_debug(pyobject))

  if dotname:
    mfac=25.4
    dotlines.append("nodesep=%f;\noverlap=false;\nmode=KK\n"%(deltaH/mfac))
    for origin,info in dotloops.items():
      r=info["size"]/mfac
      dotlines.append("%d "%origin +
      "[shape=circle, fixedsize=true, width=%f, height=%f ]\n"%(r,r))
      for i in xrange(len(info["deps"])):
        rt=info["deps"][i]
        hl=info["esizes"][i]/mfac
        dotlines.append("%d -> %d [arrowhead=none, dir=none, len=%f weight=10]\n"%(origin, rt, hl))

    dotlines.append("}\n")

    df=open(dotname,"w")
    map(df.write, dotlines)
    df.close()


  if len(svgname)>0:
    fig.savefig(svgname, type="svg", dpi=100)
  if len(pngname)>0:
    fig.savefig(pngname, type="png")

#-------------------------------------------------------------------------------
# Overloaded 3D export function
def ThreeD(base, unit, helix, loop, coil, seq, dotparens, icolors, dsb, material,
          jname="", svgname="",
          pngname="", dotname="", dpi=100., colorstrands=False, drawbases=True,
          show25d=False, colorbaseprob=True, colorbaseid=False, colorbar=True,
          baseidbar=False, colorbarspace=False, titlebottom="", titletop="",
          drawbasenumbers=True, numberinterval=2, debug=False):

  global global_cmap, pyplot
  if debug:
    deltab=0
  else:
    deltab=1
  #~ pdb.set_trace()
  # Set up canvas if SVG or PNG file to be produced
  if len(svgname)>0 or len(pngname)>0 or show25d or jname or dotname:
    fig=pyplot.figure(dpi=dpi,figsize=(5.12,5.12))
    xaxis_start=0.075
    if colorbarspace or colorbar or baseidbar:
      xaxis_start=0.0
    ax=fig.add_axes([xaxis_start, 0.075, 0.875, 0.875]) # Square area, fix for non colorbar
    ax.set_frame_on(False)

    # Turn off minor ticking altogether
    ax.set_xticks([])
    ax.set_xticks([],minor=True)
    ax.set_yticks([])
    ax.set_yticks([],minor=True)

    #~ ax.set_aspect(1.,adjustable="box")
    #~ ax.autoscale_view(tight=True, scalex=False, scaley=False)

  else: ax=False

  # First do bases
  bases=[]
  arrows=[]
  strandmax=-1
  strands=[]

  # For creating dot files
  dotlines=[]
  loopdeps={}
  helixdeps={}
  dotloops={}

  # Set up options for how bases will be drawn
  if drawbases:
    baseradius=dsb/5
  else:
    baseradius=dsb/7


  # Find mins/maxes for picture
  minx=[100000, 100000]
  maxx=[-100000, -100000]

  for i in xrange(len(unit)):
    u=unit[i]
    if u.baseInfo:
      for j in [0,1]:
        if u.baseInfo.x3[j]<minx[j]: minx[j]=u.baseInfo.x3[j]
        if u.baseInfo.x3[j]>maxx[j]: maxx[j]=u.baseInfo.x3[j]


  # Now pad the smaller edge
  dx=[0,0]

  for j in [0,1]:
    dx[j]=maxx[j]-minx[j]

  pad=max(dx)/50+baseradius*2+dsb/3
  x=[minx[0]-pad,maxx[0]+pad]
  y=[minx[1]-pad,maxx[1]+pad]

  bbox=[[x[0], x[1]], [y[0], y[1]]]
  j=0
  if dx[0]>dx[1]:
    j=1
  ddx=dx[(j+1)%2]-dx[j]

  bbox[j][0]=bbox[j][0]-ddx/2
  bbox[j][1]=bbox[j][1]+ddx/2

  if ax:
    if len(seq)==0:
      drawbases=False

    if titletop:
      txt=matplotlib.text.Text(x=(maxx[0]+minx[0])/2,y=bbox[1][1],
        text=titletop, family="sans-serif", ha="center", size="14")
      ax.add_artist(txt)
    if titlebottom:
      txt=matplotlib.text.Text(x=(maxx[0]+minx[0])/1.9,y=bbox[1][0]-dx[1]/11,
        text=titlebottom, family="sans-serif", ha="center", size="14")
      ax.add_artist(txt)
      #~ txt.set_x(0.5)
      #~ txt.set_y(0.05)


    ax.set_xlim(bbox[0])
    ax.set_ylim(bbox[1])

    miny_pad=bbox[1][0] # min y
    maxy_pad=bbox[1][1] # max y
    rec=matplotlib.patches.Rectangle([bbox[0][0], bbox[1][0]],
      bbox[0][1]-bbox[0][0], bbox[1][1]-bbox[1][0], ec="#0000ff", fill=False, lw=0.0)
    ax.add_patch(rec)

    # Set up line widths, dsb is the coil arc length
    black="#000000"
    fig_width,fig_height=fig.get_size_inches() # These should be the same as given earlier
    box_width=rec.get_width()
    box_height=rec.get_height()
    if box_width>box_height:
      scalefac=box_width/fig_width
    else:
      scalefac=box_height/fig_height

    bbonewidth=48/scalefac
    strutwidth=24/scalefac
    numsize=float(max(6.0, bbonewidth*2)) # In points
    numsize=float(min(numsize, 12.0)) # In points
    numtmp=numsize/72.72*scalefac
    numoffset=max(baseradius*1.4, dsb/2*scipy.log10(float(len(base)))*scipy.log10(scalefac))
    digits=scipy.ceil(scipy.log10(len(base)))
    numoffset=min(numoffset, baseradius+numtmp*digits)
    ticklen=numoffset/2.5
    if drawbases:
      ticklen=numoffset/2

    # Set up options for how bases will be drawn
    if drawbases:
      if colorbaseprob or colorbaseid:
        clw=strutwidth
      else:
        clw=0.0
    else:
      clw=0.0

    if colorbar or colorbaseprob:
      cmap = mpl.cm.jet
      norm = mpl.colors.Normalize(vmin=0.0, vmax=1.0)
      mappable=mpl.cm.ScalarMappable(norm,cmap)

    if colorbar: # Draw a colorbar
      ax1 = fig.add_axes([0.9, 0.1, 0.0175, 0.8])
      ax1.set_frame_on(False)
      ax1.yaxis.tick_right()
      ax1.set_axisbelow(True)

      ax1.xaxis.set_visible(False)
      #~ ax1.yaxis.set_visible(False)
      locator=matplotlib.ticker.MultipleLocator(0.2)
      ax1.yaxis.set_major_locator(locator)

      cb1 = mpl.colorbar.ColorbarBase(ax1, cmap=cmap, norm=norm,
              orientation='vertical', drawedges=False, ticks=None)
      cb1.outline.set_linewidth(0.0)
      rec2=matplotlib.patches.Rectangle([0, 0], 1,1,
        ec="#ffffff", fill=False, lw=0.0)
      ax1.add_patch(rec2)

      ax1.autoscale_view(tight=True)
      ax1.yaxis.set_ticks([0.0, 0.2, 0.4, 0.6, 0.8, 1.0], minor=False)
      ax1.yaxis.set_ticklabels(["0.0", "0.2", "0.4", "0.6", "0.8", "1.0"], minor=False)

    if baseidbar: # Draw a colorbar
      ax2 = fig.add_axes([0.9, 0.1, 0.35, 0.8])
      ax2.set_frame_on(False)
      ax2.yaxis.tick_right()
      ax2.set_axisbelow(True)

      ax2.xaxis.set_visible(False)
      ax2.yaxis.set_visible(False)
      ax2.autoscale_view(True)
      ax2.set_aspect(1.,adjustable="box")
      i=0
      dy=0.125
      for bi in ['A', 'C', 'G', 'T']:
        dot = matplotlib.patches.Circle([0.05,0.9-i*dy], radius=0.035,
          fill=True, zorder=2, ec="#000000", fc=colorid[bi], lw=0.0)
        ax2.add_patch(dot)
        btext=matplotlib.text.Text(0.15, 0.9-i*dy,
          text=bi, color="#000000", va="center", ha="center",
          zorder=3, family="sans-serif", size=11.0)
        ax2.add_artist(btext)

        i=i+1

  ddx=[dx[0]-dx[1], dx[1]-dx[0]]
  if ddx[0]<0:
    j=0
    i=1
  else:
    j=1
    i=0
  minx[j]=minx[j]-ddx[i]/2
  maxx[j]=maxx[j]+ddx[i]/2


  dotlines.append("digraph G {\n")

  for i in xrange(len(unit)):
    u=unit[i]
    if False and debug:
      print "\nu[%d]:\ntype=%s\nbase=%s"%(i,u.type, u.base)
      print "pair=%d\nnext=%s"%(u.pair, u.next)
      print "prevloop=%d\nprevside=%s"%(u.prevloop, u.prevside)

    if u.baseInfo:
      pair=-1
      if u.pair>=0:
        pair=int(unit[u.pair].base)
      if u.next>=0:
        next=int(unit[u.next].base)
      btype=0
      if i>0 and unit[i-1].type==-5:
        btype=-5

      if i+1<len(unit) and unit[i+1].type==-3:
        btype=-3

      bx=float(u.baseInfo.x3[0])
      by=float(maxy_pad-u.baseInfo.x3[1]+miny_pad)
      bz=float(u.baseInfo.x3[2])

      desc= {
          "id":"%s"%u.base,
          "pair":pair,
          "loc":(bx, by, bz),
          "strand":int(u.strand),
          "structype":u.type,
          "type":btype,
      }
      bases.append(desc)
      if ax:
        if drawbases:
          fillcolor="#ffffff"
          if len(u.baseInfo.type)>0:
            btext=matplotlib.text.Text(u.baseInfo.x3[0],
            u.baseInfo.x3[1],
            text=u.baseInfo.type, color="#000000", va="center", ha="center",
            zorder=3, family="sans-serif", size=float(bbonewidth*3))
            ax.add_artist(btext)
        else:
          fillcolor="#CCCCCC"
        ec="#000000"
        if colorbaseprob and u.prob>=0:
          ec=matplotlib.colors.rgb2hex(mappable.to_rgba(u.prob)[:3])
          if not drawbases:
            fillcolor=ec
        elif colorbaseid and u.baseInfo.type:
          ec=colorid[u.baseInfo.type]
          if not drawbases:
            fillcolor=ec
        desc["color"]=ec
        desc["probability"]=u.prob

        dot = matplotlib.patches.Circle(u.baseInfo.x3[:2], radius=baseradius,
          fill=True, label=u.type, ec=ec, lw=clw, zorder=2, fc=fillcolor)
        ax.add_patch(dot)

    if u.strand>strandmax:
      strandmax=u.strand
      if colorstrands:
        col=icolors[u.strand]
      else:
        col="000000"
      strands.append({"id":int(u.strand), "color":col})

  # Start loops and helices
  helices=[]
  hloops=[]
  loops=[]
  lastloop=-1
  lastloopi=-1
  parentloops=[]
  helices_update=[]
  lasthelix=-1

  i=0
  loopsdone={}
  loopstodo={}
  helicesdone={}
  loopmap={} # Maps our loop ID to JSON loop ID
  lc=0
  for li in xrange(len(loop)):
    loopstodo[li]=loop[li]
    if loop[li].coilnum[0]>-1:
      loopmap[li]=lc
      lc=lc+1
  todo=[i] # current work queue
  tlilist=[] # special helix work queue
  newstrand=True
  deltaH=0.0
  strnd=0
  while len(loopsdone) < len(loop):
    if len(todo)==0:
      todo.append(loopstodo.keys()[0]) # Ugh...
      newstrand=True

    i=todo.pop()
    if not loopsdone.has_key(i):
      l=loop[i]
      if debug:
        print "\nL%d:\nsidenbase=%s\nsideunit=%s"%(i,l.sidenbase, l.sideunit)
        print "sidebase=%s\ncoilnum=%s"%(l.sidebase, l.coilnum)
        print "toloop=%s\ntohelix=%s"%(l.toloop, l.tohelix)
        print "center=%s\nradius=%s"%(l.center, l.radius)
        print "strand=%s"%(l.strand)


      if l.coilnum[0]>=0 or l.tohelix[0]>=0: # This is a hairpin/multiloop
        # Determine bases and segments
        lastloop=lastloop+1
        lastloopi=i
        strand=l.strand
        if l.tohelix[0]>=0: # This loop has children
          parentloops.append(lastloop)

        # Dot file creation
        if not dotloops.has_key(lastloop):
          dotloops[lastloop]={"size": float(l.radius), "deps": [], "esizes":[]}

        # End dotfile creation
        baseInfo=[]
        segments=[]
        units=[]

        for j in xrange(len(l.sidenbase)):
          bnum=l.sidenbase[j]
          try:
            btype=int(unit[l.sideunit[j]].type)
          except ValueError:
            btype=0
          if bnum<0: continue
          unit_id=l.sideunit[j]
          strand_id=l.strand[j]
          base_id=l.sidebase[j]

          oldk=-1
          donespecial=False
          segunits=[]
          for k in xrange(bnum+2):
            addbase=False

            if k==0 or k==bnum+1:
              ind=k%bnum
              try:
                btype=int(unit[l.sideunit[j]+k+ind].type)
              except ValueError:
                btype=0
              if btype<0 or (k>0 and l.toloop[j]==-1 and l.tohelix[j]==-1):
                addbase=True
            if (k>0 and k<bnum+1):
              addbase=True
            if addbase:
                baseInfo.append(int(base_id))
                if not loopdeps.has_key(int(base_id)):
                  loopdeps[base_id]={}
                segunits.append(unit_id)


            if k>0:
              strnd=unit[l.sideunit[j]].strand
              segments.append([int(base_id-1), int(base_id), int(strnd)])
              if not loopdeps.has_key(int(base_id-1)):
                loopdeps[base_id-1]={}

              loopdeps[int(base_id-1)][int(base_id)]=1
            oldk=k
            base_id=base_id+1
            unit_id=unit_id+1
          units.append(segunits)

        # Determine child helices
        childHelices=[]
        for ch in l.tohelix:
          if ch>=0: # Found a child helix
            childHelices.append(int(ch))
        if newstrand:
          parentHelix=-1
        elif tlilist and tlilist[-1][0]==i:
          parentHelix=lasthelix+1
        else:
          parentHelix=lasthelix

        # Use coil center info
        if l.coilnum[0]>=0:
          coilcenter=coil[l.coilnum[0]].center
        else:
          coilcenter=l.center
        if debug:
          pylab.text(coilcenter[0],coilcenter[1],'L%d,%d'%(len(loops),i), size=numsize)
          pylab.plot([coilcenter[0]],[coilcenter[1]],'rs')

        desc= {
            "id": "%d"%lastloop,
            "bases":baseInfo,
            "parentHelix":parentHelix,
            "childHelices": childHelices,
            "radius":float(l.radius),
            "center":(float(coilcenter[0]), float(maxy_pad-coilcenter[1] + miny_pad), float(coilcenter[2])),
            "segments":segments,
            "loopDirection":"0",
            "coils": filter(lambda x: x>-1, map(int, l.coilnum)),
            "loop": int(i)
        }
        loops.append(desc)

        # Draw arcs
        if debug: print segments
        numseg=len(segments)
        if numseg>0:
          begin_ind=segments[0][0]
          strnd=segments[0][2]
        if debug: print "numseg=%d, units=\n%s"%(numseg, units)
        icoil=0
        angle=0. # Delta theta for each segment
        orient=-2

        tseg=[]
        seg_orient=[]
        seg_flip=[]
        si=0
        for s in range(numseg):
          if si==0:
            tseg.append([])
          si+=1
          tseg[-1].append(segments[s])
          if numseg==0 or s==numseg-1 or segments[s][2]!= segments[s+1][2]or \
              segments[s][1]!=segments[s+1][0]:
            si=0

        for s in range(len(tseg)): # Number of "coil"s
          ltsg=len(tseg[s])
          if ltsg==1:
            begin_ind=tseg[s][0][0]
            end_ind=tseg[s][0][1]
            orient, angle, ang1, ang2=\
              get_orient_angle(base[begin_ind].x3, base[end_ind].x3,  l.radius, coilcenter)

            if orient==1:
              begin_ang=ang1
              end_ang=ang2
            elif orient==-1:
              begin_ang=ang2
              end_ang=ang1
            flip=0
            orient2=orient
          else:
            begin_ind1=tseg[s][0][0]
            end_ind1=tseg[s][0][1]
            orient1, angle1, ang1, ang2=\
              get_orient_angle(base[begin_ind1].x3, base[end_ind1].x3,  l.radius, coilcenter)

            begin_ind2=end_ind1
            end_ind2=tseg[s][1][1]
            orient2, angle2, ang3, ang4=\
              get_orient_angle(base[begin_ind2].x3, base[end_ind2].x3,  l.radius, coilcenter)

            end_ind3=tseg[s][-1][1]
            ang5=get_angle2(coilcenter, base[end_ind3].x3)

            dst=numpy.zeros(3)
            dst[0]=dist(base[begin_ind1].x3, base[end_ind1].x3) # 7->8
            dst[1]=dist(base[begin_ind2].x3, base[end_ind2].x3) # 8->9
            dst[2]=dist(base[end_ind3].x3, base[end_ind1].x3) # 20->8

            if math.fabs((dst[0]-dst[1])/dst[1])>0.3: #We have a flip
              flip=1
              begin_ind=end_ind3
              end_ind=begin_ind1
              begin_ang=ang5
              end_ang=ang1
            else:
              flip=0
              begin_ind=begin_ind1
              end_ind=end_ind3
              begin_ang=ang1
              end_ang=ang5

            if angle2>0 and orient2>0:
              tmp=end_ang
              end_ang=begin_ang
              begin_ang=tmp

          seg_flip.append(flip)
          seg_orient.append(orient2)
          if ax:
            if colorstrands: col="#%s"%icolors[strnd]
            else: col=black
            if debug: print "Drawing arc from %f to %f\nangle=%s"%(begin_ang, end_ang, angle)
            arc = matplotlib.patches.Arc((float(coilcenter[0]), float(coilcenter[1])),
              float(l.radius*2), float(l.radius*2), 0., float(end_ang.real),
              float(begin_ang.real), lw=bbonewidth,
              ec=col, fill=False, zorder=1)
            ax.add_patch(arc)

          # Draw arrow
          if units[icoil]:
            ii=units[icoil][0]

            ei=ii+len(units[icoil])-1
            eb=unit[ei].base
            while ei-ii<len(units[icoil])+1 and unit[ei].base<eb \
              and ei-eb<len(units[icoil])+1:
              ei=ei+1
              if debug: print "ei =",ei

            delta_ang=scipy.pi/180*(end_ang-begin_ang)/len(units[icoil])/3

            arp=False
            if units[icoil] and unit[ii].type==-3:
              arp=arrow_coil(unit, ii, coilcenter, l.radius,
                scipy.pi/180*begin_ang, -delta_ang, dsb/3, ax)
            elif units[icoil] and unit[ei+1].type==-3:
              arp=arrow_coil(unit, ei, coilcenter, l.radius,
                scipy.pi/180*end_ang, delta_ang, dsb/3, ax)

            if arp: # We have an arrow vector
              arrows.append(draw_arrow(arp, unit[units[icoil][0]].strand, col,
                dsb*1.1, bbonewidth, strutwidth, maxy_pad, miny_pad, ax))

          icoil=icoil+1
          angles=[]
          orient=-2

        desc["sense"]=orient
        desc["seg_flip"]=seg_flip
        desc["seg_orient"]=seg_orient

        # Draw base numbers
        if drawbasenumbers:
          for bn in baseInfo:
            if bn>0 and not (bn+deltab)%numberinterval==0: continue
            tarp=arrow_helix([coilcenter,base[bn].x3], numoffset)
            dsc="%d"%(bn+deltab)
            ax.annotate(dsc, xy=base[bn].x[:2], zorder=1,
              xytext=(tarp[1][0],tarp[1][1]), fontsize=numsize, ha="center", va="center",
              arrowprops=dict(ec="#777777", shrinkA=2,shrinkB=2,
              arrowstyle="-", lw=strutwidth/2, fill=False))

        loopsdone[i]=l

        del loopstodo[i]
        newstrand=False
        for t in l.toloop[::-1]:
          if t>0: todo.append(t)

# Special case of a helix without a loop
        tlilist_new=[]
        sp_hi=len(l.toloop)-1
        for tli in l.toloop[::-1]:
          if tli<0:
            sp_hi-=1
            continue

          sl=loop[tli] # Special loop
          #~ pdb.set_trace()
          if debug:
            print "\n**L%d:\nsidenbase=%s\nsideunit=%s"%(tli,sl.sidenbase, sl.sideunit)
            print "sidebase=%s\ncoilnum=%s"%(sl.sidebase, sl.coilnum)
            print "toloop=%s\ntohelix=%s"%(sl.toloop, sl.tohelix)
            print "center=%s\nradius=%s"%(sl.center, sl.radius)
            print "strand=%s"%(sl.strand)
            print "%d: unit.base=%d sidebase=%d"%(tli,unit[sl.sideunit[0]].base,sl.sidebase[0])

          if sl.coilnum[0]>=0 and unit[sl.sideunit[0]].base==sl.sidebase[0] \
            and l.tohelix[sp_hi]>=0 \
            and not helicesdone.has_key(l.tohelix[sp_hi]):
            dohelix=tlilist_new.append((tli,l.tohelix[sp_hi],i))
          sp_hi-=1



        if tlilist and tlilist[-1][0]==i:
          tli, hnum, parent_i=tlilist.pop()
          helicesdone[hnum]=True
          sl=loop[tli] # Special loop

          hloops.append((tli,sl))
          ci=tli
          child=loop[ci]

          if ci<0 or child.coilnum[0]>=0 and debug: # We came to the end of a helix
            print "found helix in loop[%d], child_index=%d"%(i,ci)

          lasthelix=lasthelix+1

          while parentloops and not lasthelix in loops[parentloops[-1]]["childHelices"]:
            parentloops.pop()

          baseInfo=[]
          unitInfo=[]
          oldh=None
          angle=0.0

          hlength=0.0
          for hindex in xrange(len(hloops)):
            ind,h=hloops[hindex]

            baseInfo.append(int(h.sidebase[0]))
            unitInfo.append(int(h.sideunit[0]))
            bi=h.sidebase[0]

          # One special case
          iside=0
          while 1:
            if iside+1>=len(h.sidebase) or h.sidebase[iside+1]<0: break
            iside=iside+1

          bi2=h.sidebase[iside]+h.sidenbase[iside]+1
          baseInfo.append(int(bi2))
          unitInfo.append(int(h.sideunit[iside]+h.sidenbase[iside]+1))

          # Draw helix lines to canvas

          if ax:
            strnd=unit[unitInfo[0]].strand
            if colorstrands: col="#%s"%icolors[strnd]
            else: col=black
            pts=numpy.array([[base[baseInfo[0]].x3[0],base[baseInfo[1]].x3[0]],
              [base[baseInfo[0]].x3[1],base[baseInfo[1]].x3[1]]]).T
            line=matplotlib.patches.Polygon(pts, zorder=1,
              ec=col, lw=strutwidth, fill=False)
            ax.add_patch(line)

            # Draw arrow if needed
            pts_vec=[pts[::-1],pts]
            for ui in range(len(unitInfo)):
              if unit[unitInfo[ui]+1].type==-3:
                arp=arrow_helix_perp(pts_vec[ui], numoffset)
                if arp:
                  arrows.append(draw_arrow(arp, unit[unitInfo[ui]+1].strand, col,
                  dsb*1.1, bbonewidth, strutwidth, maxy_pad, miny_pad, ax))

            if drawbasenumbers and (baseInfo[0]==0 or (baseInfo[0]+deltab)%numberinterval==0):
              tarp=arrow_helix(pts[::-1], numoffset)
              if debug:
                dsc="%d, H%d"%(baseInfo[0]+deltab, len(helices))
              else:
                dsc="%d"%(baseInfo[0]+deltab)
              ax.annotate(dsc, xy=base[baseInfo[0]].x[:2], zorder=1,
                xytext=(tarp[1][0],tarp[1][1]), fontsize=numsize, ha="center", va="center",
                arrowprops=dict(ec="#777777", shrinkA=2,shrinkB=2,
                arrowstyle="-", lw=strutwidth/2, fill=False))

            if drawbasenumbers and (baseInfo[1]+deltab==0 or (baseInfo[1]+deltab)%numberinterval==0):
              tarp=arrow_helix(pts, numoffset)
              if debug:
                dsc="%d, H%d"%(baseInfo[1]+deltab, len(helices))
              else:
                dsc="%d"%(baseInfo[1]+deltab)
              ax.annotate(dsc, xy=base[baseInfo[1]].x[:2], zorder=1,
                xytext=(tarp[1][0],tarp[1][1]), fontsize=numsize, ha="center", va="center",
                arrowprops=dict(ec="#777777", shrinkA=0.075,shrinkB= 0.075,
                arrowstyle="-", lw=strutwidth/2, fill=False))


          bx1=base[baseInfo[0]].x3
          bx2=base[baseInfo[1]].x3
          get_angle=False

          if ci>-1:
            if child.coilnum[0]>=0:
              coilcenter=coil[child.coilnum[0]].center
            else:
              coilcenter=child.center

            x1=(bx1+bx2)/2
            x2=coilcenter
            get_angle=True
          elif lastloopi>-1:

            parent=loop[lastloopi]
            if parent.coilnum[0]>=0:
              coilcenter=coil[parent.coilnum[0]].center
            else:
              coilcenter=parent.center
            x2=(bx1+bx2)/2
            x1=coilcenter
            get_angle=True

          if get_angle:
            if debug:
              print "Calculate angle:"
              print "x1=\n%s"%x1
              print "x2=\n%s"%x2
              ptsc=numpy.array([[x1[0],x2[0]], [x1[1],x2[1]]]).T
              print ptsc
              linec=matplotlib.patches.Polygon(ptsc, zorder=1,
                ec="#ff0000", fc="#ff0000", lw=bbonewidth, fill=True)
              ax.add_patch(linec)

            deltaH=scipy.sqrt((x2[0]-x1[0])**2 + \
              (x2[1]-x1[1])**2)

            angle=scipy.arccos((x2[0]-x1[0])/deltaH)
            if (x2[1]<x1[1]):angle=2*scipy.pi-angle

          if newstrand:
            parentLoop=-1
          else:
            try:
              parentLoop=loopmap[parent_i]
            except:
              parentLoop=-1

          childLoop=loopmap[ci]

          desc= {
              "id": "%d"%lasthelix,
              "bases":baseInfo,
              "parentLoop":parentLoop,
              "childLoop":childLoop,
              "strand0":int(unit[h.sideunit[0]].strand),
              "strand1":int(unit[h.sideunit[1]].strand),
              "angle":2*scipy.pi-float(angle),
              "length":float(hlength),
              "flips":"-1"
          }

          bl=len(baseInfo)/2
          if len(parentloops)>0:
            dotloops[parentloops[-1]]["deps"].append(childLoop)
            dotloops[parentloops[-1]]["esizes"].append(hlength)

          for bi in xrange(bl):
            if not helixdeps.has_key(baseInfo[bi]):
              helixdeps[baseInfo[bi]]={}
            if not helixdeps.has_key(baseInfo[bi+bl]):
              helixdeps[baseInfo[bi+bl]]={}
            if bi<bl-1:
              helixdeps[baseInfo[bi]][baseInfo[bi+1]]=1
              helixdeps[baseInfo[bi+bl]][baseInfo[bi+bl+1]]=1
            helixdeps[baseInfo[bi]][baseInfo[bi+bl]]=1

          helices.append(desc)
          newstrand=False

          hloops=[]
          if not tli in todo: todo.append(tli) # Do this loop again
          loopstodo[tli]=sl
        tlilist=tlilist_new+tlilist # Prepend the new special loops

# End Special case
      elif l.coilnum[0]<0: # We found a helix ----------------------------------
        hloops.append((i,l))
        ci=l.toloop[0]
        child=loop[ci]
        if ci<0 or child.coilnum[0]>=0: # We came to the end of a helix
          lasthelix=lasthelix+1

          while parentloops and not lasthelix in loops[parentloops[-1]]["childHelices"]:
            parentloops.pop()

          baseInfo=[]
          unitInfo=[]
          oldh=None
          angle=0.0

          hlength=0.0
          for hindex in xrange(len(hloops)):
            ind,h=hloops[hindex]
            loopsdone[ind]=h
            del loopstodo[ind]

            baseInfo.append(int(h.sidebase[0]))
            unitInfo.append(int(h.sideunit[0]))
            if hindex>0:
              bi=h.sidebase[0] # index of unit
              oldbi=oldh.sidebase[0]
              deltaH=scipy.sqrt((base[bi].x3[0]-base[oldbi].x3[0])**2 + \
                      (base[bi].x3[1]-base[oldbi].x3[1])**2)

              hlength=hlength+deltaH
            oldh=h

          # Two special cases
          baseInfo.append(int(oldh.sidebase[0]+1))
          unitInfo.append(int(oldh.sideunit[0]+1))
          hlength=hlength+deltaH
          ind,h=hloops[0]
          baseInfo.append(int(h.sidebase[1]+1))
          unitInfo.append(int(h.sideunit[1]+1))

          # And now the second series of bases
          for ind,h in hloops:
            baseInfo.append(int(h.sidebase[1]))
            unitInfo.append(int(h.sideunit[1]))
            oldh=h

          # Draw helix lines to canvas

          hnb=len(baseInfo)/2 # helix_nbases
          strnd=unit[unitInfo[0]].strand
          if colorstrands: col="#%s"%icolors[strnd]
          else: col=black

          # Draw arrow
          for s,e in [[0,hnb-1],[-1,hnb]]:
            if unit[unitInfo[e]+1].type==-3:
              P1=base[baseInfo[s]].x3[:2]
              P2=base[baseInfo[e]].x3[:2]

              arp=arrow_helix([P1, P2], dsb/4)

              astrnd=unit[unitInfo[e]].strand
              arrows.append(draw_arrow(arp, astrnd, col, dsb*1.1, bbonewidth,
                strutwidth, maxy_pad, miny_pad, ax))


          # Set parent and child loops
          childLoop=loopmap[hloops[-1][1].toloop[0]]
          parentLoop=loopmap[hloops[0][1].toloop[1]]

          if ax:
            # Draw struts
            for bi in range(hnb):
              pts=numpy.array([[base[baseInfo[bi]].x3[0],base[baseInfo[hnb+bi]].x3[0]],
                [base[baseInfo[bi]].x3[1],base[baseInfo[hnb+bi]].x3[1]]])
              pts_old=pts

              line=matplotlib.patches.Polygon(pts.T,
                ec="#000000", lw=strutwidth, fill=False, zorder=0)
              ax.add_patch(line)

              # Move number away from loop
              if hnb>=2 and (bi==0 or bi==hnb-1):
                if bi==0 and parentLoop>=0 and loops[parentLoop]["bases"]:
                  x1=(base[baseInfo[bi]].x3*3 + base[baseInfo[bi+1]].x3*2)/5
                  x2=(base[baseInfo[hnb+bi]].x3*3 + base[baseInfo[hnb+bi+1]].x3*2)/5
                  pts=numpy.array([[x1[0],x2[0]], [x1[1],x2[1]]])
                elif bi==hnb-1 and childLoop>=0:# and len(loops)>childLoop+1 and loops[childLoop]["bases"]:
                  x1=(base[baseInfo[bi]].x3*3 + base[baseInfo[bi-1]].x3*2)/5
                  x2=(base[baseInfo[hnb+bi]].x3*3 + base[baseInfo[hnb+bi-1]].x3*2)/5
                  pts=numpy.array([[x1[0],x2[0]], [x1[1],x2[1]]])


              # Draw base numbers
              if drawbasenumbers and (baseInfo[bi]==0 or (baseInfo[bi]+deltab)%numberinterval==0):
                tarp=arrow_helix(pts.T[::-1], numoffset)
                if debug:
                  dsc="%d, H%d"%(baseInfo[bi]+deltab, len(helices))
                else:
                  dsc="%d"%(baseInfo[bi]+deltab)

                ax.annotate(dsc, xy=base[baseInfo[bi]].x[:2], zorder=1,
                  xytext=(tarp[1][0],tarp[1][1]), fontsize=numsize, ha="center", va="center",
                  arrowprops=dict(ec="#777777", shrinkA=0.075,shrinkB= 0.075,
                  arrowstyle="-", lw=strutwidth/2, fill=False))

                #~ txt=matplotlib.text.Text(x=arp[1][0],y=arp[1][1], size=numsize,
                  #~ text=dsc, family="sans-serif", ha="center", va="center")
                #~ ax.add_artist(txt)
                #~ # Draw line:
                #~ arp=arrow_helix2([pts_old[:,0],arp[1]], ticklen)
                #~ line=matplotlib.patches.Polygon(arp, zorder=1,
                  #~ ec="#777777", lw=strutwidth/2, fill=False)
                #~ ax.add_patch(line)

              if drawbasenumbers and (baseInfo[hnb+bi]+deltab==0 or (baseInfo[hnb+bi]+deltab)%numberinterval==0):
                tarp=arrow_helix(pts.T, numoffset)
                if debug:
                  dsc="%d, H%d"%(baseInfo[hnb+bi]+deltab, len(helices))
                else:
                  dsc="%d"%(baseInfo[hnb+bi]+deltab)
                ax.annotate(dsc, xy=base[baseInfo[hnb+bi]].x[:2], zorder=1,
                  xytext=(tarp[1][0],tarp[1][1]), fontsize=numsize, ha="center", va="center",
                  arrowprops=dict(ec="#777777", shrinkA=0.075,shrinkB= 0.075,
                  arrowstyle="-", lw=strutwidth/2, fill=False))

                #~ txt=matplotlib.text.Text(x=arp[1][0],y=arp[1][1], size=numsize,
                  #~ text="%d"%(baseInfo[hnb+bi]+deltab), family="sans-serif", ha="center", va="center")
                #~ ax.add_artist(txt)
                #~ # Draw line:
#~
                #~ arp=arrow_helix2([pts_old[:,1],arp[1]], ticklen)
                #~ line=matplotlib.patches.Polygon(arp, zorder=1,
                  #~ ec="#777777", lw=strutwidth/2, fill=False)
                #~ ax.add_patch(line)

          bx1=base[baseInfo[0]].x3
          bx2=base[baseInfo[hnb]].x3
          get_angle=False

          if ci>-1:
            if child.coilnum[0]>=0:
              coilcenter=coil[child.coilnum[0]].center
            else:
              coilcenter=child.center

            x1=(bx1+bx2)/2
            x2=coilcenter
            get_angle=True
          elif lastloopi>-1:

            parent=loop[lastloopi]
            if parent.coilnum[0]>=0:
              coilcenter=coil[parent.coilnum[0]].center
            else:
              coilcenter=parent.center
            x2=(bx1+bx2)/2
            x1=coilcenter
            get_angle=True

          if get_angle:
            if debug:
              print "Calculate angle:"
              print "x1=\n%s"%x1
              print "x2=\n%s"%x2
              ptsc=numpy.array([[x1[0],x2[0]], [x1[1],x2[1]]]).T
              print ptsc
              linec=matplotlib.patches.Polygon(ptsc, zorder=1,
                ec="#ff0000", fc="#ff0000", lw=bbonewidth, fill=True)
              ax.add_patch(linec)

            deltaH=scipy.sqrt((x2[0]-x1[0])**2 + \
              (x2[1]-x1[1])**2)

            angle=scipy.arccos((x2[0]-x1[0])/deltaH)
            if (x2[1]<x1[1]):angle=2*scipy.pi-angle

          desc= {
              "id": "%d"%lasthelix,
              "bases":baseInfo,
              "parentLoop":parentLoop,
              "childLoop":childLoop,
              "strand0":int(unit[h.sideunit[0]].strand),
              "strand1":int(unit[h.sideunit[1]].strand),
              "angle":2*scipy.pi-float(angle),
              "length":float(hlength),
              "flips":"-1"
          }

          bl=len(baseInfo)/2
          if parentLoop>=0:
            dotloops[parentLoop]["deps"].append(childLoop)
            dotloops[parentLoop]["esizes"].append(hlength)

          for bi in xrange(bl):
            if not helixdeps.has_key(baseInfo[bi]):
              helixdeps[baseInfo[bi]]={}
            if not helixdeps.has_key(baseInfo[bi+bl]):
              helixdeps[baseInfo[bi+bl]]={}
            if bi<bl-1:
              helixdeps[baseInfo[bi]][baseInfo[bi+1]]=1
              helixdeps[baseInfo[bi+bl]][baseInfo[bi+bl+1]]=1
            helixdeps[baseInfo[bi]][baseInfo[bi+bl]]=1

          helices.append(desc)
          newstrand=False

          hloops=[]
        if ci>=0:
          todo.append(ci)


  # Create the final dictionary
  properties={"is25d"        :"1",
              "baseradius"   :"%f"%float(baseradius),
              "backbonewidth":"%f"%float(bbonewidth/72.72*dpi),
              "strutwidth"   :"%f"%float(strutwidth/72.72*dpi),
              "strutlength"  : float(2*Globals.rdh),
              "material"     : material,
              }
  if ax:
    td=ax.transData
    bbox=ax.get_window_extent()
    if hasattr(td,"inverted"):
      tdi=td.inverted()
      ll=tdi.transform(numpy.array([bbox.xmin,bbox.ymin]))
      ur=tdi.transform(numpy.array([bbox.xmax,bbox.ymax]))
      properties["ll_fig"]=[float(ll[0]), float(ll[1])]
      properties["ur_fig"]=[float(ur[0]), float(ur[1])]
      properties["ypad"]=[float(miny_pad), float(maxy_pad)]

      tda = td.get_affine().get_matrix()

      properties["affine"]=(float(tda[0,0]), float(tda[1,0]), float(tda[0,1]),
                            float(tda[1,1]), float(tda[0,2]), float(tda[1,2]))

      for bd in bases: # Array of dicts
        bx,by,bz=bd["loc"]
        pl=td.transform(numpy.array([bx,by]))
        bd["pixloc"]=[float(pl[0]), float(pl[1])]
    else:
      ll=td.inverse_xy_tup((bbox.xmin(),bbox.ymin()))
      ur=td.inverse_xy_tup((bbox.xmax(),bbox.ymax()))
      properties["ll_fig"]=[float(ll[0]), float(ll[1])]
      properties["ur_fig"]=[float(ur[0]), float(ur[1])]
      properties["ypad"]=[float(miny_pad), float(maxy_pad)]

      afn=td.as_vec6_val()
      properties["affine"]=(float(afn[0]), float(afn[2]), float(afn[1]),
                            float(afn[3]), float(afn[4]), float(afn[5]))
      for bd in bases: # Array of dicts
        bx=bd["loc"][0]
        by=bd["loc"][1]
        pl=td.xy_tup((bx,by))
        bd["pixloc"]=[float(pl[0]),float(pl[1])]

  struc={ "properties"  : properties,
          "bases"       : bases,
          "strands"     : strands,
          "helices"     : helices,
          "loops"       : loops,
          "anchor"      : [{ "base":"0", "type":"2", "elementID":"0" }],
          "arrows"      : arrows,
          "sequence"    : seq,
          "dotparens"   : dotparens
        }

  pyobject={"struc" : struc}

  if len(jname)>0:
    jfile=open(jname,"w")
    jfile.write(json.write_debug(pyobject))

  if dotname:
    mfac=25.4
    dotlines.append("nodesep=%f;\noverlap=false;\nmode=KK\n"%(deltaH/mfac))
    for origin,info in dotloops.items():
      r=info["size"]/mfac
      dotlines.append("%d "%origin +
      "[shape=circle, fixedsize=true, width=%f, height=%f ]\n"%(r,r))
      for i in xrange(len(info["deps"])):
        rt=info["deps"][i]
        hl=info["esizes"][i]/mfac
        dotlines.append("%d -> %d [arrowhead=none, dir=none, len=%f weight=10]\n"%(origin, rt, hl))

    dotlines.append("}\n")

    df=open(dotname,"w")
    map(df.write, dotlines)
    df.close()


  if len(svgname)>0:
    fig.savefig(svgname, type="svg")

  if len(pngname)>0:
    fig.savefig(pngname, type="png")


if __name__=='__main__':
  load_matplotlib()
