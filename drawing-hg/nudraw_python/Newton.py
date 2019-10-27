import math
import numpy
import scipy

#-------------------------------------------------------------------------------
def loop(nside, sidesecant_in, ds, nds):
  sidesecant=sidesecant_in[:nside]
  r = sidesecant.max()/2;     # max of stem and nick gaps
  r = max([r,ds/2]) + .001;   # max of ssDNA gap
  f = nds*scipy.arcsin(ds/(2*r)) + numpy.sum(scipy.arcsin(sidesecant/(2*r))) \
     - math.pi;
  while f > 1e-12: # use Newton's method to find radius on which bases fall for this loop (default: 3 steps)
    #~ dfdr        = -1/(2*math.pow(r,2))*nds*ds/math.sqrt(1 - math.pow(ds/(2*r),2)) \
                  #~ -1/(2*math.pow(r,2))*numpy.sum(sidesecant/scipy.sqrt(1 - \
                  #~ scipy.power(sidesecant/(2*r),2) ) )

    a=-1/(2*math.pow(r,2))*nds*ds/math.sqrt(1 - math.pow(ds/(2*r),2))
    b=sidesecant/scipy.sqrt(1 - scipy.power(sidesecant/(2*r),2) )
    c=-1/(2*math.pow(r,2))*numpy.sum(b)
    dfdr=a+c

    r       = r - f/dfdr;
    f = nds*scipy.arcsin(ds/(2*r)) + numpy.sum(scipy.arcsin(sidesecant/(2*r))) \
    - math.pi;

  sideangle = 2*scipy.arcsin(sidesecant/(2*r));
  ang_ds = 2*scipy.arcsin(ds/(2*r))
  return r, sideangle, ang_ds, f
#new = [r sideangle(1:nside) ang_ds]


#-------------------------------------------------------------------------------
def groove(inclination, rdh, brise, minorgroove):
    #~print "newtongroove\n"
    #global strutrise, dthgroove
    dthgroove = math.pi
    #~print "inclination=%f, rdh=%f, brise=%f, minorgroove=%f"%(inclination,rdh,brise,minorgroove)
    sa = math.sin(dthgroove - math.pi)
    ca = math.cos(dthgroove - math.pi)
    strutrise = rdh*math.tan(inclination)*math.sqrt(math.pow(sa,2) +\
                math.pow(1 + ca,2)) # rise of strut along its length
    f = dthgroove/(2*math.pi) * brise - strutrise - minorgroove
    #~print "dthgroove=%s\nsa=%s\nca=%s\nstrutrise=%s\nf=%s"%(dthgroove, sa,ca,strutrise,f)
    for iter in range(3): # use Newton's method to find radius on which bases fall for this loop (default: 3 steps)
        dstrutdth = rdh*math.tan(inclination)* 0.5*math.pow(math.pow(sa, 2) + \
                    math.pow(1+ca, 2), -0.5) \
            *(2*sa*ca - 2*(1 + ca)*sa)
        dfdth = brise/(2*math.pi)  - dstrutdth

        dthgroove       = dthgroove - f/dfdth
        #~print "dstrutdth=%s dfdth=%s dthgroove=%s"%(dstrutdth, dfdth, dthgroove)
        sa = math.sin(dthgroove - math.pi)
        ca = math.cos(dthgroove - math.pi)
        strutrise = rdh*math.tan(inclination)*math.sqrt(pow(sa,2) + pow(1 + ca,2))
        f = dthgroove/(2*math.pi) * brise - strutrise - minorgroove
        #~print 'sa=%f ca=%f strutrise=%f f=%f'%(sa, ca, strutrise, f)

    strutlength = rdh/math.cos(inclination) * math.sqrt(pow(sa,2) + pow(1 + ca,2)) #length of strut
    strutangle = math.pi + scipy.arctan2(rdh*sa,rdh*(1+ca)) # angle of strut relative to x axis
                    # with origin at one end of strut and rotating
                    # around z axis
    return dthgroove, strutrise, strutlength, strutangle
