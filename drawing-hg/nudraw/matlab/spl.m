%
% nudraw
% 
% A program for visualizing DNA Structures and Devices
% Copyright 2004
% Niles A. Pierce
% California Institute of Technology
% Version 1.0 coded 2004 (April 9,10,13,14)
%
%
function [ys,dys,ddys] = spl(SPLINE,xp,yp,xs,ybc1,ybc2,ksbc1,ksbc2)
%
%
%

if SPLINE == 3
      sp = spapi(4,xp,yp');
elseif SPLINE == 4
      sp = spapi(5,xp,yp');
elseif SPLINE == 5
      sp = spapi(6,xp,yp');
elseif SPLINE == 6
      sp = spapi(7,xp,yp');
elseif SPLINE == 7
      sp = spapi(8,xp,yp');
elseif SPLINE == 9
      sp = spapi(10,xp,yp');
end

if SPLINE >=3

      dsp = fnder(sp);
      ddsp = fnder(dsp);

      ys = fnval(sp,xs);
      dys = fnval(dsp,xs);
      ddys = fnval(ddsp,xs);

elseif SPLINE == 0
%
%
%     ******************************************************************
%     *                                                                *
%     *  cubic spline: returns y,y',y''                                *
%     *                                                                *
%     ******************************************************************
%
%     *****************************************************************
%     * spline unkowns are the second derivatives at the np control
%     * points, found by equating the first derivatives at control
%     * points.
%     * np-1 = the number of segments. the points go from 1->np 
%     * sd: diagonal 1->np, su: upper 1->np-1, sl: lower 2->np,
%     * sb: rhs 1->np, sm: second deriv solution 1->np
%     *****************************************************************
%
      np = length(xp);
      ns = length(xs);
%
%     define tridiagonal spline system
%
      for i=2:np
      dx(i) = xp(i) - xp(i-1);
      dy(i) = yp(i) - yp(i-1);
      end
%
      for i=2:np-1
      sd(i) = (dx(i) + dx(i+1))/3.0;
      sl(i) = dx(i)/6.0;
      su(i) = dx(i+1)/6.0;
      sb(i) = dy(i+1)/dx(i+1) - dy(i)/dx(i);
      end
%
%     set spline bc's
%
      if ksbc1 == 2              % set y''(x_0)=ybc1
        sd(1)    = 1.0;
        su(1)    = 0.0;
        sb(1)    = ybc1;
      elseif ksbc1 == 1          % set y'(x_0) = ybc1
        sd(1)    = dx(2)/3.0;
        su(1)    = dx(2)/6.0;
        sb(1)    = dy(2)/dx(2) - ybc1;
      elseif ksbc1 == 3       % not a knot condition 
        sd(1)    = -1.0/dx(2);          % (y''' continuous 1st and last interior nodes)
        su(1)    = 1.0/dx(2) + 1.0/dx(3);
        suu      = -1.0/dx(3);          % second upper diagonal in top row
        sb(1)    = 0.0;
%                               
        sd(1)    = sd(1) - sl(2)*suu/su(2);    % preprocess to remove second 
        su(1)    = su(1) - sd(2)*suu/su(2);    % upper diagonal to make tridiagonal
        sb(1)    = sb(1) - sb(2)*suu/su(2);
      else
      disp('no valid spline bc: quit spline')
      return
      end
      
      if ksbc2 == 2              % set y''(x_np) = ybc2
        sd(np)   = 1.0;
        sl(np)   = 0.0;
        sb(np)   = ybc2;
      elseif ksbc2 == 1          % set y'(x_np) = ybc2
        sd(np)   = dx(np)/3.0;
        sl(np)   = dx(np)/6.0;
        sb(np)   = ybc2 - dy(np)/dx(np);
      elseif ksbc2 == 3       % not a knot condition 	
        sd(np)   = -1.0/dx(np);
        sl(np)   = 1.0/dx(np-1) + 1.0/dx(np);
        sll      = -1.0/dx(np-1);                     % second lower diagonal in bottom row
        sb(np)   = 0.0;
%                               
        sd(np)    = sd(np) - su(np-1)*sll/sl(np-1);    % preprocess to remove second 
        sl(np)    = sl(np) - sd(np-1)*sll/sl(np-1);    % upper diagonal to make tridiagonal
        sb(np)    = sb(np) - sb(np-1)*sll/sl(np-1);
      else
      disp('no valid spline bc: quit spline')
      return
      end
      
      INVERT=1;
      if INVERT == 1
%
%       form matrix and solve with matlab
%      
        smat = diag(sd) + diag(su(1:np-1),1) + diag(sl(2:np),-1);
        sm = smat\sb';
      else
%
%       solve tridiagonal system for second derivatives
%
        for i=2:np
          sd(i) = sd(i) - su(i-1)*sl(i)/sd(i-1);
          sb(i) = sb(i) - sb(i-1)*sl(i)/sd(i-1);
        end
%
        sm(np) = sb(np)/sd(np);
%
        for i=np-1:-1:1
          sm(i) = (sb(i) - su(i)*sm(i+1))/sd(i);
        end
      end
%
%     solve for points along the spline
%
      for j=1:ns
%
        i = 2;
        while xs(j) > xp(i)
          i = i+1;
        end
%        disp(xp(i-1));disp(xs(j));disp(xp(i))
%
        dxp    = xp(i) - xs(j);
        dxm    = xs(j) - xp(i-1);
%
%       y
%
        ys(j)  = sm(i-1)*dxp^3/(6.0*dx(i)) ...                   
               + sm(i)*dxm^3/(6.0*dx(i))   ...   
               + dxm*(yp(i)   - dx(i)^2*sm(i)/6.0)/dx(i) ...
               + dxp*(yp(i-1) - dx(i)^2*sm(i-1)/6.0)/dx(i);   
%
%       dydx
%   
        dys(j) = -.5*sm(i-1)*dxp^2/dx(i)   ...     
               +  .5*sm(i)*dxm^2/dx(i)      ...
               + (yp(i) - yp(i-1))/dx(i)    ...
               - dx(i)*(sm(i) - sm(i-1))/6.0;
%
%       ddyddx
%   
        ddys(j) = (sm(i-1)*dxp + sm(i)*dxm)/dx(i);
      end
%
elseif SPLINE == 1
      np = length(xp);
      ns = length(xs);
%
%     define tridiagonal spline system
%
      for i=2:np
      dx(i) = xp(i) - xp(i-1);
      dy(i) = yp(i) - yp(i-1);
      end
      for j=1:ns
%
        i = 2;
        while xs(j) > xp(i)
          i = i+1;
        end
%
%       y
%
        ys(j)  = yp(i-1) + dy(i)/dx(i) * (xs(j) - xp(i-1));  
%
%       dydx
%   
        dys(j) = dy(i)/dx(i) ;
%
%       ddyddx
%   
        ddys(j) = 0;
      end
end
