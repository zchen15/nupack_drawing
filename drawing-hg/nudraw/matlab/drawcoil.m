%
% nudraw
% 
% A program for visualizing DNA Structures and Devices
% Copyright 2004
% Niles A. Pierce
% California Institute of Technology
% Version 1.0 coded 2004 (April 9,10,13,14)

function [arcratio] = drawcoil(nbase,xc1,nc1,thc1,xc2,nc2,xp,yp,zp,bctype,helical,random,shades,render,ccap,cdye,direction,seed)
% nbase 	        % number of bases in helix
% xc1(3)     	    % start coil position
% nc1(3)            % start coil vector (into coil)
% thc1 		        % start coil rotation around nc1
% xc2(3)     	    % end coil position
% nc2(3)            % end coil vector (into coil)
% xp,yp,zp          % points to be splined between end points
% bctype(2)         % 5' and 3' ends (1: first deriv bc, 2: second deriv bc)
% helical           % make a helix around the spline
% random            % make the helix angular progression and base orientations random 
                    % in absence of helix, this still randomizes base
                    % orientations around spline
% shades(2)         % colormap indices for coil, base pairs
% render(4)         % render coil, first base, middle bases, last base 
% ccap(2)           % cap 5', 3' end of coil
% cdye(2)           % dye 5', 3' end of coil
% direction         % 1: 5' end at x1, 2: 5' end at x2
% seed              % set seed if random


global nh nb rhc rbc rhelix rdh rcap ncap thcap dzb bturn brise dthb dsb dz dth ...
       cmap rdye ddye ndye majorgroove minorgroove dthgroove strutangle strutlength ...
       strutrise proptwist inclination extend5prime extend3prime extendcoil rcaphead;

   
if nc1
    nc1 = nc1/norm(nc1);  % normalize vectors (can be all zero if represent second derivatives)
end
if nc2
    nc2 = nc2/norm(nc2);
end

randn('state',seed);

nhtot(1) = (nh(1)-1)*(nbase-1) + 1;  % total points along chain
nhtot(2) = nh(2);   % total points around chain




SPLINE = 0;

helpers = 0;  % turn off or on diagnostic spline plots



ksbc1 = bctype(1);
ksbc2 = bctype(2);
sp(1) = 0;
for i=2:length(xp)
    sp(i) = sp(i-1) + sqrt((xp(i)-xp(i-1)).^2 + (yp(i)-yp(i-1)).^2 + (zp(i)-zp(i-1)).^2);
end

ss = linspace(0,sp(end),nhtot(1));

xbc1 = nc1(1);
xbc2 = -nc2(1);
[xs,dxs,ddxs] = spl(SPLINE,sp,xp,ss,xbc1,xbc2,ksbc1,ksbc2);

ybc1 = nc1(2);
ybc2 = -nc2(2);
[ys,dys,ddys] = spl(SPLINE,sp,yp,ss,ybc1,ybc2,ksbc1,ksbc2);

zbc1 = nc1(3);
zbc2 = -nc2(3);
[zs,dzs,ddzs] = spl(SPLINE,sp,zp,ss,zbc1,zbc2,ksbc1,ksbc2);

arc(1) = 0;
for i=2:length(ss)
    arc(i) = arc(i-1) + sqrt((xs(i)-xs(i-1)).^2 + (ys(i)-ys(i-1)).^2 + (zs(i)-zs(i-1)).^2);
end

if helpers
  plot3(xp,yp,zp,'go',xs,ys,zs,'b-'); hold on;
  plot3([0 nc1(1)] + xc1(1),[0 nc1(2)] + xc1(2), [0 nc1(3)] + xc1(3),'g-','linewidth',2)
  plot3([0 nc2(1)] + xc2(1),[0 nc2(2)] + xc2(2), [0 nc2(3)] + xc2(3),'g-','linewidth',2)
end

xc_a = [xs;ys;zs];

thc_a(1) = thc1*pi/180;
for i=2:nbase
   if random
     dthb_rnd   = dthb + randn(1)*sqrt(pi/8);
       thc_a(i)      = thc_a(i-1) + dthb_rnd;
   else 
       thc_a(i)      = thc_a(i-1) + dthb;           
     end
 end

if helical  % make the coil helical
%     
   xph(1) = xp(1);
   yph(1) = yp(1);
   zph(1) = zp(1);
%   
   xph(nbase) = xp(end);
   yph(nbase) = yp(end);
   zph(nbase) = zp(end);

   
   j = 1;
 
   zc_a  = linspace(0,arc(end),nbase);
   for i=2:nbase-1 %leave first and last base as is
       while arc(j) < zc_a(i) 
         j = j+1;    
       end
       frac = (zc_a(i) - arc(j-1))/(arc(j) - arc(j-1));
       xi = (1-frac)*xs(j-1) + frac*xs(j);
       yi = (1-frac)*ys(j-1) + frac*ys(j);
       zi = (1-frac)*zs(j-1) + frac*zs(j);
       dxi = (1-frac)*dxs(j-1) + frac*dxs(j);
       dyi = (1-frac)*dys(j-1) + frac*dys(j);
       dzi = (1-frac)*dzs(j-1) + frac*dzs(j);
       
       nci = [dxi dyi dzi];
       nci = nci/norm(nci); 
       % 2) rotate helix axis to vector u_n \equiv nci(3)
       % axis starts out as u_z
       % rotation is around vector u_rot = u_z x u_n
       sinth_rot = sqrt(nci(1)^2 + nci(2)^2);
       costh_rot = nci(3); 
       if sinth_rot
         u_rot = [-nci(2); nci(1); 0];
         u_rot = u_rot/norm(u_rot); % make unit vectors
         th_rot = 180/pi * atan2(sinth_rot,costh_rot);
         % rotate(h,u_rot,th_rot,[0 0 0]); % "rotate" won't work properly if xc1 \neq 0
         % (since intrinsic function won't do translation)

         % th_rot needs to be reversed compared to value for using
         % matlab intrinsic function "rotate"
         rmat    = cos(-th_rot*pi/180)*eye(3)  ...
                 + (1-cos(-th_rot*pi/180))*u_rot*u_rot' ...
                 +    sin(-th_rot*pi/180) *[0                 u_rot(3)           -u_rot(2); ...
                                           -u_rot(3)          0                   u_rot(1); ...
                                            u_rot(2)         -u_rot(1)            0];                        
      else
         rmat  = eye(3);    
      end

      
      xtmp = zeros(3,1);
      if i==2 | i == nbase-1  % move out from center of helix to radius gradually
         xtmp(1)  = rdh*cos(thc_a(i)) * 1/3;
         xtmp(2)  = rdh*sin(thc_a(i)) * 1/3;
      elseif i==3 | i == nbase-2
         xtmp(1)  = rdh*cos(thc_a(i)) * 2/3;
         xtmp(2)  = rdh*sin(thc_a(i)) * 2/3;
        
      else
         xtmp(1)  = rdh*cos(thc_a(i));
         xtmp(2)  = rdh*sin(thc_a(i));      
      end
      xtmp(3)  = 0;
      xtmp       = rmat*xtmp;
      xph(i)  = xtmp(1) + xi;
      yph(i)  = xtmp(2) + yi;
      zph(i)  = xtmp(3) + zi;
   end
   sph(1) = 0;  % base arc length on new (quasi)-helical points for each base
   for i=2:length(xph)
       sph(i) = sph(i-1) + sqrt((xph(i)-xph(i-1)).^2 + (yph(i)-yph(i-1)).^2 + (zph(i)-zph(i-1)).^2);
   end
   ssh = linspace(0,sph(end),nhtot(1));  % spline (quasi)-helical version
   [xs,dxs,ddxs] = spl(SPLINE,sph,xph,ssh,xbc1,xbc2,ksbc1,ksbc2);
   [ys,dys,ddys] = spl(SPLINE,sph,yph,ssh,ybc1,ybc2,ksbc1,ksbc2);
   [zs,dzs,ddzs] = spl(SPLINE,sph,zph,ssh,zbc1,zbc2,ksbc1,ksbc2);
   xc_a = [xs;ys;zs];
   if helpers
     plot3(xph,yph,zph,'ro',xs,ys,zs,'m-'); hold on;
   end
        
end  % helical
 
   
%
% now equidistribute points along the spline
%
%     
xpeq(1) = xp(1);
ypeq(1) = yp(1);
zpeq(1) = zp(1);
%   
xpeq(nbase) = xp(end);
ypeq(nbase) = yp(end);
zpeq(nbase) = zp(end);

for m = 1:1  % iterate on arc length to get equidistribution of points along arc
  arceq(1) = 0;  % arc length on (quasi)-helical spline
    for i=2:length(xs)
    arceq(i) = arceq(i-1) + sqrt((xs(i)-xs(i-1)).^2 + (ys(i)-ys(i-1)).^2 + (zs(i)-zs(i-1)).^2);
    end
    
    speq = linspace(0,arceq(end),nbase);
    j = 1;
    for i=2:nbase-1 %leave first and last base as is
       while arceq(j) < speq(i)
         j = j+1;    
       end
       frac = (speq(i) - arceq(j-1))/(arceq(j) - arceq(j-1));
       xpeq(i) = (1-frac)*xs(j-1) + frac*xs(j);
       ypeq(i) = (1-frac)*ys(j-1) + frac*ys(j);
       zpeq(i) = (1-frac)*zs(j-1) + frac*zs(j);
    end
    sseq = linspace(0,speq(end),nhtot(1));  % spline (quasi)-helical version
   [xs,dxs,ddxs] = spl(SPLINE,speq,xpeq,sseq,xbc1,xbc2,ksbc1,ksbc2);
   [ys,dys,ddys] = spl(SPLINE,speq,ypeq,sseq,ybc1,ybc2,ksbc1,ksbc2);
   [zs,dzs,ddzs] = spl(SPLINE,speq,zpeq,sseq,zbc1,zbc2,ksbc1,ksbc2);
%    figure(1)
%    if m==1
%        clf
%    end
%    for i=2:nbase
%       dsp(i) = speq(i) - speq(i-1);    
%    end
%    subplot(5,1,m)
%    hist(dsp,40)
%    figure(2)
end
xc_a = [xs;ys;zs];
if helpers
   plot3(xpeq,ypeq,zpeq,'ro',xs,ys,zs,'m-'); hold on;
end

arcratio = (sseq(end)/(nbase-1))/dsb;   
    
[xdot ydot zdot] = sphere(ndye);    

if helpers
   for i=1:nbase
      surface(2*xdot+xpeq(i),2*ydot+ypeq(i),2*zdot+zpeq(i),0*zdot+8)
   end
end

% define base pair
%
xb = zeros(nb);
yb = zeros(nb);
zb = zeros(nb);
dm1 = linspace(0,2/3 * rdh,nb(1));
dm2 = linspace(0,2*pi,nb(2)); 
for i=1:nb(1)
 xb(i,:) = dm1(i);  
 yb(i,:) = rbc(1)*cos(dm2); 
 zb(i,:) = rbc(2)*sin(dm2); 
end


% base pair struts

for j=1:nbase
   i = 1 + (j-1)*(nh(1)-1);
   
   % 0) normalize chain axis vector
   ncj = [dxs(i); dys(i); dzs(i)];
   ncj = ncj/norm(ncj); 
   % 2) rotate helix axis to vector u_n \equiv ncj(3)
   % axis starts out as u_z
   % rotation is around vector u_rot = u_z x u_n
   sinth_rot = sqrt(ncj(1)^2 + ncj(2)^2);
   costh_rot = ncj(3); 
   if sinth_rot
      u_rot = [-ncj(2); ncj(1); 0];
      u_rot = u_rot/norm(u_rot); % make unit vectors
      th_rot = 180/pi * atan2(sinth_rot,costh_rot);
      % rotate(h,u_rot,th_rot,[0 0 0]); % "rotate" won't work properly if xc1 \neq 0
      % (since intrinsic function won't do translation)

      % th_rot needs to be reversed compared to value for using
      % matlab intrinsic function "rotate"
      rmat    = cos(-th_rot*pi/180)*eye(3)  ...
              + (1-cos(-th_rot*pi/180))*u_rot*u_rot' ...
              +    sin(-th_rot*pi/180) *[0                 u_rot(3)           -u_rot(2); ...
                                        -u_rot(3)          0                   u_rot(1); ...
                                         u_rot(2)         -u_rot(1)            0];                        
   else
      rmat  = eye(3);    
   end
   
   rmat1     = cos(thc_a(j))*eye(3) + (1-cos(thc_a(j)))*ncj*ncj' ...  
             + sin(thc_a(j))*[0 ncj(3) -ncj(2); -ncj(3) 0 ncj(1); ncj(2) -ncj(1) 0]; % rotate around chain tangent
   for k=1:size(xb,1)
      for m = 1:size(xb,2)
         xtmp = rmat1*rmat*[xb(k,m); yb(k,m); zb(k,m)];
         xbtmp(k,m) = xtmp(1,:) + xpeq(j); 
         ybtmp(k,m) = xtmp(2,:) + ypeq(j);
         zbtmp(k,m) = xtmp(3,:) + zpeq(j);  %translate to base location
      end
   end

   % define base pair cap
   endtype = 0;
   rmat0 = cos(pi/2)*eye(3) + (1-cos(pi/2))*[0 0 0; 0 -1 0; 0 0 0] + sin(pi/2)*[0 0 1; 0 0 0; -1 0 0]; 
   [xcap,ycap,zcap] =  drawcap(0,dm1(end),rmat0,endtype,0);  % last argument not used for this case
   
   % cap
   for k=1:size(xcap,1)
      for m=1:size(xcap,2)
      xtmp = rmat1*rmat*[xcap(k,m); ycap(k,m); zcap(k,m)]; % rotate cap to end of rotated base and follow rotation around chain tangent
      xcap(k,m) = xtmp(1) + xpeq(j);
      ycap(k,m) = xtmp(2) + ypeq(j); 
      zcap(k,m) = xtmp(3) + zpeq(j);
    end
   end

   if (j==1 & render(2)) | (j==nbase & render(4)) | (j > 1 & j < nbase & render(3))
     surface(xbtmp,ybtmp,zbtmp,0*zbtmp+shades(end),'CDataMapping','direct'); hold on;
     surface(xcap,ycap,zcap,0*zcap+shades(end),'CDataMapping','direct'); hold on;
   end

end


clear xtmp;
x_a = zeros(nhtot);
y_a = zeros(nhtot);
z_a = zeros(nhtot);

%
% rotate cross-section of chain to be normal to axis of splined chain
%
th = 0:dth:2*pi;
for i=1:nhtot(1)
   % 0) normalize chain axis vector
   nci = [dxs(i); dys(i); dzs(i)];
   nci = nci/norm(nci); 
   % 2) rotate helix axis to vector u_n \equiv nci(3)
   % axis starts out as u_z
   % rotation is around vector u_rot = u_z x u_n
   sinth_rot = sqrt(nci(1)^2 + nci(2)^2);
   costh_rot = nci(3); 
   if sinth_rot
      u_rot = [-nci(2); nci(1); 0];
      u_rot = u_rot/norm(u_rot); % make unit vectors
      th_rot = 180/pi * atan2(sinth_rot,costh_rot);
      % rotate(h,u_rot,th_rot,[0 0 0]); % "rotate" won't work properly if xc1 \neq 0
      % (since intrinsic function won't do translation)

      % th_rot needs to be reversed compared to value for using
      % matlab intrinsic function "rotate"
      rmat    = cos(-th_rot*pi/180)*eye(3)  ...
              + (1-cos(-th_rot*pi/180))*u_rot*u_rot' ...
              +    sin(-th_rot*pi/180) *[0                 u_rot(3)           -u_rot(2); ...
                                        -u_rot(3)          0                   u_rot(1); ...
                                         u_rot(2)         -u_rot(1)            0];                        
   else
      rmat  = eye(3);    
   end
   
   xtmp(1,:) = rhc*cos(th);
   xtmp(2,:) = rhc*sin(th);
   xtmp(3,:) = 0;  
   xtmp = rmat*xtmp;
   x_a(i,:) = xtmp(1,:) + xc_a(1,i); % translate rotated coords 
   y_a(i,:) = xtmp(2,:) + xc_a(2,i);
   z_a(i,:) = xtmp(3,:) + xc_a(3,i); 
   
   if i == 1  % rotation for end caps in reverse direction at beginning of chain
      if (sinth_rot)
          rmat    = cos(-th_rot*pi/180 + pi)*eye(3)  ...
              + (1-cos(-th_rot*pi/180 + pi))*u_rot*u_rot' ...
              +    sin(-th_rot*pi/180 + pi) *[0                 u_rot(3)      -u_rot(2); ...
                                        -u_rot(3)          0                   u_rot(1); ...
                                         u_rot(2)         -u_rot(1)            0];  
      else
        rmat  = [-1 0 0;0 1 0; 0 0 -1];                          
      end
      if direction == 1
        endtype = 5;
        [x5_a,y5_a,z5_a] =  drawcap(i,xc_a,rmat,endtype,extendcoil);
      elseif direction == 2
        endtype = 3;
        [x3_a,y3_a,z3_a] =  drawcap(i,xc_a,rmat,endtype,extendcoil);
      end
   elseif i==nhtot(1) % rotation for end caps in same direction at end of chain
      if direction == 1
        endtype = 3;
        [x3_a,y3_a,z3_a] =  drawcap(i,xc_a,rmat,endtype,extendcoil);
      elseif direction == 2
        endtype = 5;
        [x5_a,y5_a,z5_a] =  drawcap(i,xc_a,rmat,endtype,extendcoil);        
      end
   end
end


if render(1)
   %surface(x_a,y_a,z_a, 0*z_a   +shades(1),'CDataMapping','direct'); hold on;
   
    range(1) = 1;
    range(2) = (nh(1)-1)/2 + 1;
    for ibase=3:nbase
        range(ibase) = range(ibase-1) + (nh(1)-1);    
    end
    range(nbase+1) = range(nbase) + (nh(1)-1)/2;

    for ibase=1:nbase
        ival = range(ibase):range(ibase+1);
        surface(x_a(ival,:),y_a(ival,:),z_a(ival,:), 0*z_a(ival,:) +shades(ibase),'CDataMapping','direct'); hold on;
    end

   if ccap(1)
      if direction == 1
          ishade = 1;
      elseif direction == 2
          ishade = nbase;
      end
      surface(x5_a,y5_a,z5_a,0*z5_a+shades(ishade),'CDataMapping','direct');
   end
   if ccap(2)
       if direction == 1
          ishade = nbase;
      elseif direction == 2
          ishade = 1;
      end
      surface(x3_a,y3_a,z3_a,0*z3_a+shades(ishade),'CDataMapping','direct');
   end
   if (cdye(1) & direction == 1) | (cdye(2) & direction == 2)
     xtmp = xc1 - nc1*(extendcoil + ddye);
     surface(xtmp(1)+rdye*xdot,xtmp(2)+rdye*ydot,xtmp(3)+rdye*zdot,0*zdot+shades(1),'CDataMapping','direct');       
   end
   if (cdye(2) & direction == 1) | (cdye(1) & direction == 2)
     xtmp = xc2 - nc2*(extendcoil + ddye);
     surface(xtmp(1)+rdye*xdot,xtmp(2)+rdye*ydot,xtmp(3)+rdye*zdot,0*zdot+shades(nbase),'CDataMapping','direct');              
   end
   
   

end
xlabel('x')
ylabel('y')
zlabel('z')


