# Python code for determining the rgb values bases on a color map
#
# I.e., input is an intensity x (on a scale from 0 to 1)
# and output is a tuple of the (red,green,blue) values
# for the colormap.  red, green, and blue all range from 0 to 255.
#
# Low colors are darker blue and higher colors are red.  Based on
# Matlab jet colormap
#
# The colormap is defined as follows:
#  
#

def getRGB(x):
       
       if x < 0:
              print 'Warning: x must satisfy 0 <= x <= 1, using x = 0'
              x = 0

       if x > 1:
              print 'Warning: x must satisfy 0 <= x <= 1, using x = 1'
              x = 1
              
       
       # Red channel
       if x <= 0.375:
              red = 0.0
       elif x <= 0.625:
              red = 4.0*x - 1.5
       elif x <= 0.875:
              red = 1.0
       else: # x > 0.875
              red = -4.0*x + 4.5

       # Green channel
       if x <= 0.125 or x >= 0.875:
              green = 0.0
       elif x <= 0.375:
              green = 4*x - 0.5
       elif x <= 0.625:
              green = 1.0
       else: # 0.625 < x < 0.875
              green = -4.0*x + 3.5

       # Blue channel
       if x <= 0.125:
              blue = 4.0*x + 0.5
       elif x <= 0.375:
              blue = 1.0
       elif x <= 0.625:
              blue = -4.0*x + 2.5
       else: # 0.626 < x <= 1
              blue = 0.0

       return (int(red*255), int(green*255), int(blue*255))

