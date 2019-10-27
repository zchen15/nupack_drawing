%
% nudraw
% 
% A program for visualizing DNA Structures and Devices
% Copyright 2004
% Niles A. Pierce
% California Institute of Technology
% Version 1.0 coded 2004 (April 9,10,13,14)

function [xc2,nc2,thc2,x1a,n1a,x2a,n2a,x1b,n1b,x2b,n2b] = drawhelix(nbase,xc1,nc1,thc1,dzc1,shades,render,acap,bcap,adye,bdye)
% nbase 	        % number of bases in helix
% xc1(3)     	    % start helix axis position
% nc1(3)            % start helix axis vector
% thc1 		        % start helix rotation around axis
% dzc1              % translation along helix axis from start of helix axis
% shades(3)         % colormap indices for chain a, chain b, base pairs
% render(3)         % render chain a, chain b, base pair struts
% acap(2)           % cap 5', 3' end of a chain
% bcap(2)           % cap 5', 3' end of b chain
% adye(2)           % dye 5', 3' end of a chain
% bdye(2)           % dye 5', 3' end of b chain

global nh nb rhc rbc rhelix rdh rcap ncap thcap dzb bturn brise dthb dsb dz dth ...
       cmap rdye ddye ndye majorgroove minorgroove dthgroove strutangle strutlength ...
       strutrise proptwist inclination extend5prime extend3prime extendcoil rcaphead;
   
nhtot(1) = (nh(1)-1)*(nbase-1) + 1;  % total points along chain
nhtot(2) = nh(2);   % total points around chain

helpers = 0;

% 0) normalize target helix axis vector
nc1 = nc1/norm(nc1); 

dz = dzb/(nh(1)-1);
dth = 2*pi/(nh(2)-1);
x_a = zeros(nhtot);
y_a = zeros(nhtot);
z_a = zeros(nhtot);
x_b = zeros(nhtot);
y_b = zeros(nhtot);
z_b = zeros(nhtot);

th = 0:dth:2*pi;

thc_a      = 0:dz*dthb/dzb:dthb*(nbase-1);
xc_a(1,:)  = rdh*cos(thc_a);
xc_a(2,:)  = rdh*sin(thc_a);
xc_a(3,:)  = (0:dz:dzb*(nbase-1)) - .5*strutrise; % move start of a chain so rise due to inclination of base pair is centered
                        % on helix origin

thc_b      = dthgroove + (0:dz*dthb/dzb:dthb*(nbase-1));
xc_b(1,:)  = rdh*cos(thc_b);
xc_b(2,:)  = rdh*sin(thc_b);
xc_b(3,:)  = (0:dz:dzb*(nbase-1)) + .5*strutrise;  % move start of chain b so rise due to inclination of base pair is 
                        % centered on the origin in the z direction
%
% define backbone
%

phi = pi/2 - atan(dzb/(rdh*dthb));
for i=1:nhtot(1)
   xtmp(1,:) = rhc*cos(th);
   xtmp(2,:) = rhc*sin(th);
   xtmp(3,:) = 0;  
   ca = cos(thc_a(i));  
   sa = sin(thc_a(i));  % rotate circles so that cross section is circular normal to chain axis instead of helix axis
   rmat_a = cos(phi)*eye(3) ...  % see p. 58 of "Dynamics": crandall, karnopp, kurtz and pridmore-brown
          + (1-cos(phi))*[ca^2 ca*sa 0;ca*sa sa^2 0; 0 0 0] ...
          + sin(phi)*[0 0 -sa; 0 0 ca; sa -ca 0];
   xtmp = rmat_a*xtmp;
   x_a(i,:) = xtmp(1,:) + xc_a(1,i); % translate
   y_a(i,:) = xtmp(2,:) + xc_a(2,i);
   z_a(i,:) = xtmp(3,:) + xc_a(3,i); 
   
   
   xtmp(1,:) = rhc*cos(th);
   xtmp(2,:) = rhc*sin(th);
   xtmp(3,:) = 0;  
   cb = cos(thc_b(i));
   sb = sin(thc_b(i));
   rmat_b = cos(phi)*eye(3) ...  % see p. 58 of "Dynamics": crandall, karnopp, kurtz and pridmore-brown
          + (1-cos(phi))*[cb^2 cb*sb 0;cb*sb sb^2 0; 0 0 0] ...
          + sin(phi)*[0 0 -sb; 0 0 cb; sb -cb 0];
   xtmp = rmat_b*xtmp;
   x_b(i,:) = xtmp(1,:) + xc_b(1,i); % translate
   y_b(i,:) = xtmp(2,:) + xc_b(2,i);
   z_b(i,:) = xtmp(3,:) + xc_b(3,i);   
  
   if i == nhtot(1) % must go before i==1 case so don't overwrite rmat_a and rmat_b
     endtype = 3;
     [x3_a,y3_a,z3_a] =  drawcap(i,xc_a,rmat_a,endtype,extend3prime);
     endtype = 5;
     [x5_b,y5_b,z5_b] =  drawcap(i,xc_b,rmat_b,endtype,extend5prime);
   end
   if i == 1
      rmat_a = cos(phi+pi)*eye(3) ...  % see p. 58 of "Dynamics": crandall, karnopp, kurtz and pridmore-brown
             + (1-cos(phi+pi))*[ca^2 ca*sa 0;ca*sa sa^2 0; 0 0 0] ...
             + sin(phi+pi)*[0 0 -sa; 0 0 ca; sa -ca 0];
      rmat_b = cos(phi+pi)*eye(3) ...  % see p. 58 of "Dynamics": crandall, karnopp, kurtz and pridmore-brown
             + (1-cos(phi+pi))*[cb^2 cb*sb 0;cb*sb sb^2 0; 0 0 0] ...
             + sin(phi+pi)*[0 0 -sb; 0 0 cb; sb -cb 0];
     endtype = 5;
     [x5_a,y5_a,z5_a] =  drawcap(i,xc_a,rmat_a,endtype,extend5prime);
     endtype = 3;
     [x3_b,y3_b,z3_b] =  drawcap(i,xc_b,rmat_b,endtype,extend3prime);
   end
   

end

% if rotating with "rotate" intrinsic (but can't do translation)
% h(1)     = surface(x_a,y_a,z_a, 0*z_a   +shades(1),'CDataMapping','direct'); hold on;
% h(end+1) = surface(x5_a,y5_a,z5_a,0*z5_a+shades(1),'CDataMapping','direct');
% h(end+1) = surface(x3_a,y3_a,z3_a,0*z3_a+shades(1),'CDataMapping','direct');
% h(end+1) = surface(x_b,y_b,z_b, 0*z_b   +shades(2),'CDataMapping','direct'); 
% h(end+1) = surface(x5_b,y5_b,z5_b,0*z5_a+shades(2),'CDataMapping','direct');
% h(end+1) = surface(x3_b,y3_b,z3_b,0*z3_a+shades(2),'CDataMapping','direct');


% define canonical base pair strut
%
xb = zeros(nb);
yb = zeros(nb);
zb = zeros(nb);




dm1 = linspace(0,strutlength,nb(1));
%dm1 = linspace(-rdh,rdh,nb(1));
dm2 = linspace(0,2*pi,nb(2)); 
for i=1:nb(1) 
 xb(i,:) = dm1(i);  
 yb(i,:) = rbc(1)*cos(dm2); 
 zb(i,:) = rbc(2)*sin(dm2); 
end

%rotate propeller twist along the axis of the strut
for k=1:nb(1)
   twistang = (dm1(k) - strutlength/2)/strutlength * proptwist;
   ca = cos(twistang);
   sa = sin(twistang);  %rotate around x axis in ccw direction
   rmat0 = ca*eye(3) + (1-ca)*[1 0 0; 0 0 0; 0 0 0] + sa*[0 0 0; 0 0 1; 0 -1 0];                
   for m = 1:nb(2)
      xtmp = rmat0*[xb(k,m); yb(k,m); zb(k,m)];
      xb(k,m) = xtmp(1); 
      yb(k,m) = xtmp(2);
      zb(k,m) = xtmp(3);  
   end
end

%rotate inclination around the y axis
ca = cos(inclination);   
sa = sin(inclination);   
rmat0 = ca*eye(3) + (1-ca)*[0 0 0; 0 1 0; 0 0 0] + sa*[0 0 -1; 0 0 0; 1 0 0];
for k=1:nb(1)           
   for m = 1:nb(2)
      xtmp = rmat0*[xb(k,m); yb(k,m); zb(k,m)];
      xb(k,m) = xtmp(1); 
      yb(k,m) = xtmp(2);
      zb(k,m) = xtmp(3);  
   end
end


%rotate and then translate canonical strut
ca = cos(2*pi-strutangle);   % correct for fact that rotation is CW and strutangle is CCW
sa = sin(2*pi-strutangle);   % rotate around z axis
rmat0 = ca*eye(3) + (1-ca)*[0 0 0; 0 0 0; 0 0 1] + sa*[0 1 0; -1 0 0; 0 0 0];
for k=1:nb(1)
   for m = 1:nb(2)
      xtmp = rmat0*[xb(k,m); yb(k,m); zb(k,m)];
      xb(k,m) = xtmp(1) + rdh; % translate to canonical location  
      yb(k,m) = xtmp(2);
      zb(k,m) = xtmp(3);  %translate to proper base height
   end
end


% if rotating with "rotate" intrinsic, but can't do translation
% for j=1:nbase
%     i = 1 + (j-1)*(nh(1)-1);
%    dbase = sqrt((xc_b(1,i)-xc_a(1,i))^2 + (xc_b(2,i)-xc_a(2,i))^2);
%    nx = (xc_b(1,i)-xc_a(1,i))/dbase;
%    ny = (xc_b(2,i)-xc_a(2,i))/dbase;
%    nz = 0;
%    xbtmp = xc_a(1,i) + dbase*nx*xb - ny*yb; 
%    ybtmp = xc_a(2,i) + dbase*ny*xb + nx*yb;
%    zbtmp = xc_a(3,i) + zb;
%    h(end+1) = surface(xbtmp,ybtmp,zbtmp,0*zbtmp+58,'CDataMapping','direct'); hold on;
% end

% start chain information
x1a  = xc_a(:,1);   % start helix chain position and vector 5'->3' chain
n1a  = [0; -rdh*dthb/sqrt(dzb^2 + (rdh*dthb)^2); -dzb/sqrt(dzb^2 + (rdh*dthb)^2)];
x1b  = xc_b(:,1);   % start helix chain position and vector 3'->5' chain

ca   = cos(2*pi-dthgroove); %rotation matrix is cw but groove angle is ccw
sa   = sin(2*pi-dthgroove);  % rotate around z axis 
rmat = ca*eye(3) + (1-ca)*[0 0 0; 0 0 0; 0 0 1] + sa*[0 1 0; -1 0 0; 0 0 0];
n1b  = rmat*n1a;

%n1b  = [0;  rdh*dthb/sqrt(dzb^2 + (rdh*dthb)^2); -dzb/sqrt(dzb^2 + (rdh*dthb)^2)];

% end chain information
xc2     = [0; 0; dzb*(nbase-1)];        % end helix axis
nc2     = nc1;  	                    % end helix axis vector
thc2    = dthb*(nbase-1)*180/pi;	    % end helix rotation around axis
rmat    = [cos(thc2*pi/180) -sin(thc2*pi/180) 0; sin(thc2*pi/180) cos(thc2*pi/180) 0; 0 0 1];  % rotate around z axis
x2a     = xc_a(:,end);                  % end helix chain position 5'->3' chain
n2a     = -rmat*n1a;                    % end helix chain vector 5'->3' chain
x2b     = xc_b(:,end);                  % end helix chain position 3'->5' chain	
n2b     = -rmat*n1b;                    % end helix chain vector 3'->5' chain

n1a = n1a/norm(n1a);  % normalize normal vectors
n1b = n1b/norm(n1b);
n2a = n2a/norm(n2a);
n2b = n2b/norm(n2b);

% convenient to rotate surface using matlab function rotate
% however, still need to keep track of end positions and vectors using
% rotation matrices, hence, might be more consistent just to explicitly
% compute rotation matrix and do everything manually
%
% actually, would be useful reference check to keep moving surfaces using 
% "rotate" and move end info manually, unfortunately, since there is no
% "translate" equivalent for translation, have to translate manually before
% "rotate" and this makes it messy to rotate the end points since have
% to change origin of rotation for them
%
% decided just to do everything manually in the end
%
    

% 1) first rotate around z axis amount thc1
% rotate(h,[0 0 1],thc1,[0 0 0]);  % "rotate" won't work properly if xc1 \neq  (since "rotate" won't do translation)
% (sign of sin terms seems reversed to me...????)
rmat1   = [cos(thc1*pi/180) -sin(thc1*pi/180) 0; sin(thc1*pi/180) cos(thc1*pi/180) 0; 0 0 1];  % rotate around z axis

% 2) rotate helix axis to vector u_n \equiv nc1(3)
% axis starts out as u_z
% rotation is around vector u_rot = u_z x u_n
sinth_rot = sqrt(nc1(1)^2 + nc1(2)^2);
costh_rot = nc1(3); 
if sinth_rot
   u_rot = [-nc1(2); nc1(1); 0];
   u_rot = u_rot/norm(u_rot); % make unit vectors
   th_rot = 180/pi * atan2(sinth_rot,costh_rot);
   % rotate(h,u_rot,th_rot,[0 0 0]); % "rotate" won't work properly if xc1 \neq 0
   % (since intrinsic function won't do translation)

   % th_rot needs to be reversed compared to value for using
   % matlab intrinsic function "rotate"
   rmat2   = cos(-th_rot*pi/180)*eye(3)  ...
           + (1-cos(-th_rot*pi/180))*u_rot*u_rot' ...
           +    sin(-th_rot*pi/180) *[0                 u_rot(3)           -u_rot(2); ...
                                     -u_rot(3)          0                   u_rot(1); ...
                                      u_rot(2)         -u_rot(1)            0];                        
elseif costh_rot == -1  % need special case for u_n = [0; 0; -1]
   u_rot = [0; 1; 0];
   th_rot = 180;
   rmat2   = cos(-th_rot*pi/180)*eye(3)  ...
           + (1-cos(-th_rot*pi/180))*u_rot*u_rot' ...
           +    sin(-th_rot*pi/180) *[0                 u_rot(3)           -u_rot(2); ...
                                     -u_rot(3)          0                   u_rot(1); ...
                                      u_rot(2)         -u_rot(1)            0];                        
else % special case for u_n = [0; 0; 1]
   rmat2 = eye(3);    
end
% 3) then translate the helix

% chains
for i=1:nhtot(1)
   for j=1:nhtot(2)
       xtmp = rmat2*rmat1*[x_a(i,j); y_a(i,j); z_a(i,j)]  + xc1 + nc1*dzc1;
       x_a(i,j) = xtmp(1);
       y_a(i,j) = xtmp(2);
       z_a(i,j) = xtmp(3);
       xtmp = rmat2*rmat1*[x_b(i,j); y_b(i,j); z_b(i,j)]  + xc1 + nc1*dzc1;
       x_b(i,j) = xtmp(1);
       y_b(i,j) = xtmp(2);
       z_b(i,j) = xtmp(3);
  end
end
% 5 caps
for i=1:size(x5_a,1)
   for j=1:size(x5_a,2)
      xtmp = rmat2*rmat1*[x5_a(i,j); y5_a(i,j); z5_a(i,j)] + xc1 + nc1*dzc1;
      x5_a(i,j) = xtmp(1);
      y5_a(i,j) = xtmp(2); 
      z5_a(i,j) = xtmp(3);
      xtmp = rmat2*rmat1*[x5_b(i,j); y5_b(i,j); z5_b(i,j)] + xc1 + nc1*dzc1;
      x5_b(i,j) = xtmp(1);
      y5_b(i,j) = xtmp(2); 
      z5_b(i,j) = xtmp(3);
  end
end
% 3 caps
for i=1:size(x3_a,1)
   for j=1:size(x3_a,2)
      xtmp = rmat2*rmat1*[x3_a(i,j); y3_a(i,j); z3_a(i,j)] + xc1 + nc1*dzc1;
      x3_a(i,j) = xtmp(1);
      y3_a(i,j) = xtmp(2); 
      z3_a(i,j) = xtmp(3); 
      xtmp = rmat2*rmat1*[x3_b(i,j); y3_b(i,j); z3_b(i,j)] + xc1 + nc1*dzc1;
      x3_b(i,j) = xtmp(1); 
      y3_b(i,j) = xtmp(2); 
      z3_b(i,j) = xtmp(3); 
   end
end

% base pair struts
for j=1:nbase
   i = 1 + (j-1)*(nh(1)-1);
   ca = cos(-thc_a(i));
   sa = sin(-thc_a(i));   % rotate around axis
   rmat0 = ca*eye(3) + (1-ca)*[0 0 0; 0 0 0; 0 0 1] + sa*[0 1 0; -1 0 0; 0 0 0];
   for k=1:size(xb,1)
      for m = 1:size(xb,2)
         xtmp = rmat0*[xb(k,m); yb(k,m); zb(k,m)];
         xbtmp(k,m) = xtmp(1,:); 
         ybtmp(k,m) = xtmp(2,:);
         zbtmp(k,m) = xtmp(3,:) + xc_a(3,i);  %translate to proper base height
     end
   end
   for k=1:size(xbtmp,1)
       for m=1:size(xbtmp,2)
           xtmp = rmat2*rmat1*[xbtmp(k,m); ybtmp(k,m); zbtmp(k,m)] + xc1 + nc1*dzc1;
           xbtmp(k,m) = xtmp(1);
           ybtmp(k,m) = xtmp(2);
           zbtmp(k,m) = xtmp(3);
      end
  end
  if render(3)
     surface(xbtmp,ybtmp,zbtmp,0*zbtmp+shades(3),'CDataMapping','direct'); hold on;
  end
end


x1a     = rmat2*rmat1*x1a + xc1 + nc1*dzc1; % displace along axis a distance dzc1
n1a     = rmat2*rmat1*n1a;
x1b     = rmat2*rmat1*x1b + xc1 + nc1*dzc1;
n1b     = rmat2*rmat1*n1b;
x2a     = rmat2*rmat1*x2a + xc1 + nc1*dzc1;
n2a     = rmat2*rmat1*n2a;
x2b     = rmat2*rmat1*x2b + xc1 + nc1*dzc1;
n2b     = rmat2*rmat1*n2b;
xc2     =       rmat2*xc2 + xc1 + nc1*dzc1;
nc2     =       nc2;
thc2    = thc2 + thc1;


[xs ys zs] = sphere(ndye);


if render(1) 
   if nhtot(1)>1
   surface(x_a,y_a,z_a, 0*z_a   +shades(1),'CDataMapping','direct'); %hold on;
   end
   hold on;
   if acap(1)
     surface(x5_a,y5_a,z5_a,0*z5_a+shades(1),'CDataMapping','direct');
   end 
   if acap(2)
      surface(x3_a,y3_a,z3_a,0*z3_a+shades(1),'CDataMapping','direct');
   end
   if adye(1)
     xtmp = x1a + n1a*(extend5prime + ddye); % need to add cap in addition to dye to get extension drawn out to dye location
     surface(xtmp(1)+rdye*xs,xtmp(2)+rdye*ys,xtmp(3)+rdye*zs,0*zs+shades(1),'CDataMapping','direct');       
   end
   if adye(2)
     xtmp = x2a + n2a*(extend3prime + ddye);
     surface(xtmp(1)+rdye*xs,xtmp(2)+rdye*ys,xtmp(3)+rdye*zs,0*zs+shades(1),'CDataMapping','direct');              
   end
   if helpers
      plot3([0 5*n1a(1)] + x1a(1),[0 5*n1a(2)] + x1a(2),[0 5*n1a(3)] + x1a(3),'b-','linewidth',2);
      plot3([0 5*n2a(1)] + x2a(1),[0 5*n2a(2)] + x2a(2),[0 5*n2a(3)] + x2a(3),'b-','linewidth',2);
      surf(2*xs+x1a(1),2*ys+x1a(2),2*zs+x1a(3))
      surf(2*xs+x2a(1),2*ys+x2a(2),2*zs+x2a(3))
   end
end

if render(2)
    if nhtot(1) > 1
        surface(x_b,y_b,z_b, 0*z_b   +shades(2),'CDataMapping','direct'); %hold on;
    end
    hold on;
   if bcap(1)
      surface(x5_b,y5_b,z5_b,0*z5_b+shades(2),'CDataMapping','direct');
   end
   if bcap(2)
      surface(x3_b,y3_b,z3_b,0*z3_b+shades(2),'CDataMapping','direct');
   end
   if bdye(1)
     xtmp = x2b + n2b*(extend5prime + ddye);
     surface(xtmp(1)+rdye*xs,xtmp(2)+rdye*ys,xtmp(3)+rdye*zs,0*zs+shades(2),'CDataMapping','direct');       
   end
   if bdye(2)
     xtmp = x1b + n1b*(extend3prime + ddye);
     surface(xtmp(1)+rdye*xs,xtmp(2)+rdye*ys,xtmp(3)+rdye*zs,0*zs+shades(2),'CDataMapping','direct');              
   end
   if helpers
      plot3([0 5*n1b(1)] + x1b(1),[0 5*n1b(2)] + x1b(2),[0 5*n1b(3)] + x1b(3),'r-','linewidth',2);
      plot3([0 5*n2b(1)] + x2b(1),[0 5*n2b(2)] + x2b(2),[0 5*n2b(3)] + x2b(3),'r-','linewidth',2);
      surf(2*xs+x1b(1),2*ys+x1b(2),2*zs+x1b(3))
      surf(2*xs+x2b(1),2*ys+x2b(2),2*zs+x2b(3))
  end
end

%surf(2*xs+xc2(1),2*ys+xc2(2),2*zs+xc2(3))



%colorbar


xlabel('x')
ylabel('y')
zlabel('z')
