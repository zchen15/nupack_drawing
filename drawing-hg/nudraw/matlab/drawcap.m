%
% nudraw
% 
% A program for visualizing DNA Structures and Devices
% Copyright 2004
% Niles A. Pierce
% California Institute of Technology
% Version 1.0 coded 2004 (April 9,10,13,14)

function [xcap_a,ycap_a,zcap_a] =  cap(i,xc_a,rmat_a,endtype,extend)
% i                 index along chain patches
% xc_a(3)           axis of helix (z is same for both chains before rotation)
% rcap              radius of end cap corner
% rbc(2)            radii of base strut
% th                angle array around the chain
% ncap              number of points going over radius of cap
% thcap              angle of conical cap at 3' end of chain
% rhc               radius of helix chain
% rmat_a            rotation matrix for cap to match angle of chain tube
% endtype           end of chain
% extend            amount to extend chain at end of strand before cap
                    % differs for 3' and 5' ends in helix and in coil

global nh nb rhc rbc rhelix rdh rcap ncap thcap dzb bturn brise dthb dsb dz dth ...
       cmap rdye ddye ndye majorgroove minorgroove dthgroove strutangle strutlength ...
       strutrise proptwist inclination extend5prime extend3prime extendcoil rcaphead;

   
th = 0:dth:2*pi;


if endtype ==  5
%if endtype ==  3 % quick fix for polarity issue
   xtmp(1,:) = rhc*cos(th);
   xtmp(2,:) = rhc*sin(th);
   xtmp(3,:) = 0;
   xtmp = rmat_a*xtmp;
   xcap_a(1,:) = xtmp(1,:) + xc_a(1,i);
   ycap_a(1,:) = xtmp(2,:) + xc_a(2,i);
   zcap_a(1,:) = xtmp(3,:) + xc_a(3,i);

   for j=2:ncap-1
      dcap = pi/2 * (j-2)/(ncap-3);
      xtmp(1,:) = ((rhc - rcap) + rcap*cos(dcap))*cos(th);
      xtmp(2,:) = ((rhc - rcap) + rcap*cos(dcap))*sin(th);
      xtmp(3,:) = rcap*sin(dcap) + extend; % offset from last base by rbc(1) 
      xtmp = rmat_a*xtmp;
      xcap_a(j,:) = xtmp(1,:) + xc_a(1,i);
      ycap_a(j,:) = xtmp(2,:) + xc_a(2,i);
      zcap_a(j,:) = xtmp(3,:) + xc_a(3,i);
   end
   xtmp(1,:) = 0*cos(th);
   xtmp(2,:) = 0*sin(th);
   xtmp(3,:) = rcap + extend; %*sin(phi);
   xtmp = rmat_a*xtmp;
   xcap_a(ncap,:) = xtmp(1,:) + xc_a(1,i);
   ycap_a(ncap,:) = xtmp(2,:) + xc_a(2,i);
   zcap_a(ncap,:) = xtmp(3,:) + xc_a(3,i);   
% elseif endtype == 3
%       xtmp(1,:) = rhc*cos(th);
%       xtmp(2,:) = rhc*sin(th);
%       xtmp(3,:) = 0;  
%       xtmp = rmat_a*xtmp;
%       xcap_a(1,:) = xtmp(1,:) + xc_a(1,i);
%       ycap_a(1,:) = xtmp(2,:) + xc_a(2,i);
%       zcap_a(1,:) = xtmp(3,:) + xc_a(3,i);
%    for j=2:ncap-1
%       dcap = (pi/2 - thcap) * (j-2)/(ncap-3);
%       xtmp(1,:) = ((rhc - rcap) + rcap*cos(dcap))*cos(th);
%       xtmp(2,:) = ((rhc - rcap) + rcap*cos(dcap))*sin(th);
%       xtmp(3,:) = rcap*sin(dcap) + extend; %offset from last base by rbc(1) 
%       xtmp = rmat_a*xtmp;
%       xcap_a(j,:) = xtmp(1,:) + xc_a(1,i);
%       ycap_a(j,:) = xtmp(2,:) + xc_a(2,i);
%       zcap_a(j,:) = xtmp(3,:) + xc_a(3,i);
%    end
%    for j = ncap:2*ncap
%       dcap = thcap * ((2*ncap)-j)/ncap;  % goes from thcap to zero
%       xtmp(1,:) = rcap*sin(dcap)*cos(th);
%       xtmp(2,:) = rcap*sin(dcap)*sin(th);
%       xtmp(3,:) = (rhc - rcap)*tan(thcap) + rcap*cos(dcap) + extend; 
%       xtmp = rmat_a*xtmp;
%       xcap_a(j,:) = xtmp(1,:) + xc_a(1,i);
%       ycap_a(j,:) = xtmp(2,:) + xc_a(2,i);
%       zcap_a(j,:) = xtmp(3,:) + xc_a(3,i);
%    end  
elseif endtype == 3
%elseif endtype == 5 % quick fix
      xtmp(1,:) = rhc*cos(th);
      xtmp(2,:) = rhc*sin(th);
      xtmp(3,:) = 0;  
      xtmp = rmat_a*xtmp;
      xcap_a(1,:) = xtmp(1,:) + xc_a(1,i);
      ycap_a(1,:) = xtmp(2,:) + xc_a(2,i);
      zcap_a(1,:) = xtmp(3,:) + xc_a(3,i);
      
      xtmp(1,:) = rhc*cos(th);
      xtmp(2,:) = rhc*sin(th);
      xtmp(3,:) = extend;  
      xtmp = rmat_a*xtmp;
      xcap_a(2,:) = xtmp(1,:) + xc_a(1,i);
      ycap_a(2,:) = xtmp(2,:) + xc_a(2,i);
      zcap_a(2,:) = xtmp(3,:) + xc_a(3,i);
      
   for j=3:ncap-1
      dcap = (pi - thcap) * (j-3)/(ncap-4) - pi/2;
      xtmp(1,:) = ((rcaphead - rcap) + rcap*cos(dcap))*cos(th);
      xtmp(2,:) = ((rcaphead - rcap) + rcap*cos(dcap))*sin(th);
      xtmp(3,:) = rcap + rcap*sin(dcap) + extend; %offset from last base by rbc(1) 
      xtmp = rmat_a*xtmp;
      xcap_a(j,:) = xtmp(1,:) + xc_a(1,i);
      ycap_a(j,:) = xtmp(2,:) + xc_a(2,i);
      zcap_a(j,:) = xtmp(3,:) + xc_a(3,i);
   end
   for j = ncap:2*ncap
      dcap = thcap * ((2*ncap)-j)/ncap;  % goes from thcap to zero
      xtmp(1,:) = rcap*sin(dcap)*cos(th);
      xtmp(2,:) = rcap*sin(dcap)*sin(th);
      xtmp(3,:) = rcap + (rcaphead - rcap)*tan(thcap) + rcap*cos(dcap) + extend; 
      xtmp = rmat_a*xtmp;
      xcap_a(j,:) = xtmp(1,:) + xc_a(1,i);
      ycap_a(j,:) = xtmp(2,:) + xc_a(2,i);
      zcap_a(j,:) = xtmp(3,:) + xc_a(3,i);
   end  
elseif endtype ==  0  % cap for base pair
   rx = rbc(2);
   ry = rbc(1);
   dcap = linspace(pi/2,0,ncap);
   for j=1:ncap
      xtmp(1,:) = rx*cos(dcap(j))*cos(th);
      xtmp(2,:) = ry*cos(dcap(j))*sin(th);
      xtmp(3,:) = rcap*sin(dcap(j)); % *sin(phi);  % radius starts at edge of last pair 
      xtmp = rmat_a*xtmp;
      xcap_a(j,:) = xtmp(1,:) + xc_a;
      ycap_a(j,:) = xtmp(2,:);
      zcap_a(j,:) = xtmp(3,:);
   end   
end

