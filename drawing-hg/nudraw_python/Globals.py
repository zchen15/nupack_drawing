# Global parameters
# Conrad D. Steenberg <conrad.steenberg@caltech.edu>
# Nov 20, 2008

import math
import Newton
import copy

nh=[0,0]
nb=[0,0]
rbc=[0,0]
nh[0]   = 20 + 1        # points along chain per base %20 default
nh[1]   = 30 + 1        # points around chain         % 30 default
nb[0]   = 1  + 1        # points along base pair
nb[1]   = 30 + 1        # points around base pair   % 30 default



rhc = None
rhelix = None
rdh = None
rcap = None
ncap = None
thcap = None
dzb = None
bturn = None
brise = None
dthb = None
dsb = None
dz = None
dth = None
majorgroove = None
minorgroove = None
dthgroove = None
strutangle = None
strutlength = None
strutrise = None
proptwist = None
inclination = None
extend5prime = None
extend3prime = None
extendcoil = None
rcaphead = None
ndye = None
ddye = None
rdye = None
cmap = None
colors = None
color_names = None
icolors = None

ipp=0.0139

# Colormap
cmap = [[ 109./255,  216./255,  45./255 ],     # light green 0
          [1,       .2857,     0       ],      # orange 2
          [164./255, 0,        0       ],      # red 3
          [181./255, 145./255,  209./255],     # lavender 1
          [0./255,   100./255,  0./255  ],     # dark green 4
          [1.,       .9375,    0       ],      # yellow 5
          [0.,        0.9286,  1       ],      # turquoise 6
          [0.,        0.,      0.8571  ],      # dark blue 7
          [142./255, 73./255,   5./255  ],     # brown 8
          [28./255,  206./255,  40./255 ],     # med green 9
          [255./255, 160./255,  204./255],     # pink 10
          [0.4,      0.4,      0.4     ],      # dark gray 11
          [0./255,   0./255,   0./255   ],     # black 12
          [79./255,  0./255,   147./255 ],     # purple 13
          [229./255, 0./255,   153./255 ],     # magenta 14
          [242./255, 206./255, 204./255 ],     # tan 15
          [.7778,   .7778,    .7778    ],      # pale gray 16
          [68./255,  102./255, 255./255 ],     # cornflowerblue 17
          [135./255, 206./255, 235./255 ],     # lightblue 18
          [0,        0,       128./255 ],      # navyblue 19
          [34./255,  139./255, 34./255  ],     # forest green 20
          [124./255, 252./255, 0.       ],     # lawn green 21
          [255./255, 204./255, 102./255 ],     # light orange 22
          [96./255,  114./255, 175./255 ],     # shadyblue 23
          [169./255, 212./255, 113./255 ],     # shadygreen 24
          [24/255.,  73/255., 1.0],            # skyblue
          [34/255.,  139/255., 34/255.],       # leaf green
          [0./255.,  120./255., 80./255.],     # sea green
          [.58,   .58,    .58    ]             # med gray 16
          #[1., 1., 1. ]                        # white
        ]
global_cmap=cmap
color_names=["light_green", "orange", "red", "lavender","dark_green", "yellow",
        "turquoise", "dark_blue", "brown", "med_green", "pink", "dark_gray",
        "black", "purple", "magenta", "tan", "light_gray", "cornflowerblue",
        "light_blue", "navy_blue", "forest_green", "lawn_green", "light_orange",
        "shady_blue", "shady_green", "sky_blue", "leaf_green", "sea_green",
        "med_gray", "white"]
colors={}
icolors={}
for i in range(len(cmap)):
  colors[color_names[i]]=cmap[i]
  icolors[i]="%02x%02x%02x"%(int(cmap[i][0]*255), \
                          int(cmap[i][1]*255), int(cmap[i][2]*255))

#Domain colors - extend this to all colors?
domainCMap=[ \
            [128./255,             208./255,              79./255],              #(med green -- NUPACK)
            [255./255,             153./255,               0./255],               #(med orange -- NUPACK)
            [ 68./255,             102./255,             255./255],             #(med blue -- NUPACK)
            [237./255,              31./255,              36./255],              #(bright red)
            [145./255,              74./255,             156./255],             #(magenta)
            [159./255,              95./255,              47./255],              #(med brown)
            #~ Dark colors
            [  1./255,             105./255,              55./255],              #(dark green)
            [241./255,              90./255,              36./255],              #(dark orange)
            [ 33./255,              28./255,              92./255],              #(dark blue)
            [133./255,               0./255,              14./255],              #(dark red)
            [ 77./255,              19./255,             110./255],             #(dark purple)
            [ 76./255,              61./255,              44./255],              #(dark brown)
            #~ Light colors
            [202./255,             224./255,             145./255],             #(pale green)
            [255./255,             218./255,               0./255],               #(yellow)
            [128./255,             242./255,             255./255],             #(pale blue)
            [247./255,             169./255,             168./255],             #(pink)
            [195./255,             156./255,             255./255],             #(lavender)
            [210./255,             180./255,             140./255]              #(tan)
          ]



domainColors={} #=copy.deepcopy(icolors)
for i in range(len(domainCMap)):
  domainColors[i]="%02x%02x%02x"%(int(domainCMap[i][0]*255), \
                          int(domainCMap[i][1]*255), int(domainCMap[i][2]*255))


colorid={"A": colors["leaf_green"],
         "C": colors["sky_blue"],
         "G": colors["black"],
         "T": colors["red"],
         "U": colors["red"],
         ".": colors["black"],
         "(": colors["yellow"],
         ")": colors["orange"],
         "?": "#999999"
        }

def getDomainColors(ind):
  return domainColors[ind%len(domainColors)]

def init_globals(material, show3d=False):
  global nh, nb, rhc, rbc, rhelix, rdh, rcap, ncap, thcap, dzb, bturn, brise, \
  dthb, dsb, dz, dth, majorgroove, minorgroove, dthgroove, strutangle, \
  strutlength, strutrise, proptwist, inclination, extend5prime, extend3prime, \
  extendcoil, rcaphead, ndye, ddye, rdye, cmap, colors, color_names, icolors

  golden=1.61803399
  e=2.71828183
  pi=math.pi
  cbrt2=math.pow(2.,0.2)
  if show3d:
    #~ print "show3d=%s"%show3d
    rhc     = 1.5           # radius of helix chain
    rbc[0]  = .3            # major radius of base pair strut (elliptical: .5 default)
    rbc[1]  = .25            # minor radius of base pair strut (elliptical: .2 default)
    if material == 1:      # A DNA (RNA/RNA duplex or RNA/DNA duplex)
      rhelix     = 13.0       # radius of double helix to outside of chain axis
      dzb     = 2.6          # stacking height per base along helix
      bturn   = 10.0         # bases per turn in the helix (10 is default value)
      majorgroove = .45*dzb*bturn          # height of major groove along helix axis in angstroms (less than .5 for narrow, deep major groove)
      inclination = -19*math.pi/180  # angle of strut relative x axis (spinning around y axis)
                                              # should tilt down from 5' (start
                                              # of A chain) to 3' (start of B
                                              # chain) strand
      proptwist = 11.8*math.pi/180   # prop twist around base pair axis in ccw direction
    elif material == 2: # B DNA
      rhelix     = 10.          # radius of double helix to outside of chain axis
      dzb     = 3.4             # stacking height per base along helix
      bturn   = 10.5            # bases per turn in the helix (10.5 is modern value)
      majorgroove = 22./34.*(dzb*bturn)          # height of major groove along helix axis in angstroms (22 is default, 17 for symmetric debugging)
                                  # take ratio of 22/34 and multiply by
                                  # modern value of brise
      inclination = 1.2*math.pi/180    # angle of strut relative x axis (spinning around y axis)
                                              # should tilt up from 5' (start
                                              # of A chain) to 3' (start of B
                                              # chain) (1.2 for B DNA)
      proptwist = 11.4*math.pi/180     # prop twist around base pair axis in ccw direction %11.4
  else:
    rhc     = 0.5              # radius of helix chain
    rbc[0]  = .5/3             # major radius of base pair strut (elliptical: .5 default)
    rbc[1]  = .25/3            # minor radius of base pair strut (elliptical: .2 default)

    if material == 1:      # A DNA (RNA/RNA duplex or RNA/DNA duplex)
      rhelix     = 3.0     # radius of double helix to outside of chain axis
      dzb     = rhelix     # stacking height per base along helix
      bturn   = 10.0       # bases per turn in the helix (10 is default value)
      majorgroove = .45*dzb*bturn          # height of major groove along helix axis in angstroms (less than .5 for narrow, deep major groove)
      inclination = -19*math.pi/180  # angle of strut relative x axis (spinning around y axis)
                                              # should tilt down from 5' (start
                                              # of A chain) to 3' (start of B
                                              # chain) strand
      proptwist = 11.8*math.pi/180   # prop twist around base pair axis in ccw direction
    elif material == 2: # B DNA
      rhelix     = 2.6          # radius of double helix to outside of chain axis
      dzb     = 3.0             # stacking height per base along helix
      bturn   = 9.5            # bases per turn in the helix (10.5 is modern value)
      majorgroove = 22./34.*(dzb*bturn)          # height of major groove along helix axis in angstroms (22 is default, 17 for symmetric debugging)
                                  # take ratio of 22/34 and multiply by
                                  # modern value of brise
      inclination = 1.2*math.pi/180    # angle of strut relative x axis (spinning around y axis)
                                              # should tilt up from 5' (start
                                              # of A chain) to 3' (start of B
                                              # chain) (1.2 for B DNA)
      proptwist = 11.4*math.pi/180     # prop twist around base pair axis in ccw direction %11.4

  extend5prime = 1*rbc[0]   # extension of chain at end so that strut doesn't hit cap
  extend3prime = 3*rbc[0]   # extension of chain at end so that strust doesn't hit cap
                             # if majorgroove=minorgroove then works with
                             # both = rbc[0]
  extendcoil = 1*rbc[0]  # for coil can be same at each end becauses no major/minor groove issue

  rdh     = rhelix - rhc     # radius of helix to center of chain
  rcap    = .5               # radius of capping corner
  ncap    = 30               # number of points along cap corners as function of r
  thcap   = 55*math.pi/180        # degrees of cone at 3' end of chain
  rcaphead = 1.25*rhc + rcap # radius of arrowhead at 3' end of chain

  brise   = dzb*bturn        # rise along helix axis per full turn
  dthb    = 2*math.pi/bturn       # stacking twist per base along helix
  dsb     = cbrt2*math.sqrt(math.pow(rdh*dthb,2) + math.pow(dzb,2)) # arc length along chain for one base
  dz      = dzb/(nh[0]-1)    # patch increment along chain
  dth     = 2*math.pi/(nh[1]-1)   # patch increment around chain
  rdye    = 2*rhc            # radius of dye molecule
  ddye    = math.sqrt(math.pow(rdye,2)  - math.pow(rhc,2))    # displacement of dye center along chain axis from last base
  ndye    = 50               # resolution of dye sphere
  minorgroove = brise - majorgroove
  #~print "minorgroove=%f majorgroove=%f"%(minorgroove, majorgroove)
  dthgroove,strutrise,strutlength,strutangle = \
    Newton.groove(inclination, rdh, brise, minorgroove)

  #~ gg.cmap = [[ 109./255,  216./255,  45./255 ],     # light green 0
            #~ [1,       .2857,     0       ],     # orange 2
            #~ [164./255, 0,        0       ],     # red 3
            #~ [181./255, 145./255,  209./255],     # lavender 1
            #~ [0./255,   100./255,  0./255  ],     # dark green 4
            #~ [1.,       .9375,    0       ],     # yellow 5
            #~ [0.,        0.9286,  1       ],     # turquoise 6
            #~ [0.,        0.,      0.8571  ],     # dark blue 7
            #~ [142./255, 73./255,   5./255  ],     # brown 8
            #~ [28./255,  206./255,  40./255 ],     # med green 9
            #~ [255./255, 160./255,  204./255],     # pink 10
            #~ [0.4,      0.4,      0.4     ],     # dark gray 11
            #~ [0./255,   0./255,   0./255   ],     # black 12
            #~ [79./255,  0./255,   147./255 ],     # purple 13
            #~ [229./255, 0./255,   153./255 ],     # magenta 14
            #~ [242./255, 206./255, 204./255 ],     # tan 15
            #~ [.7778,   .7778,    .7778    ],     # pale gray 16
            #~ [68./255,  102./255, 255./255 ],   # cornflowerblue 17
            #~ [135./255, 206./255, 235./255 ],    #lightblue 18
            #~ [0,        0,       128./255 ],                # navyblue 19
            #~ [34./255,  139./255, 34./255  ],     # forest green 20
            #~ [124./255, 252./255, 0.       ],      # lawn green 21
            #~ [255./255, 204./255, 102./255 ],  #light orange 22
            #~ [96./255,  114./255, 175./255 ],   # shadyblue 23
            #~ [169./255, 212./255, 113./255 ],  # shadygreen 24
            #~ [1., 1., 1. ]  # white
          #~ ]
#~
  #~ gg.color_names=["light_green", "orange", "red", "lavender","dark_green", "yellow",
          #~ "turquoise", "dark_blue", "brown", "med_green", "pink", "dark_gray",
          #~ "black", "purple", "magenta", "tan", "light_gray", "cornflowerblue",
          #~ "light_blue", "navy_blue", "forest_green", "lawn green", "light_orange",
          #~ "shady_blue", "shady_green", "white"]
  #~ colors={}
  #~ icolors={}
  #~ for i in range(len(cmap)):
    #~ colors[color_names[i]]=cmap[i]
    #~ icolors[i]="%02x%02x%02x"%(int(cmap[i][0]*255), \
                            #~ int(cmap[i][1]*255), int(cmap[i][2]*255))

