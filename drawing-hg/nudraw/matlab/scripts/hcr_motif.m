%
% nudraw
% 
% A program for visualizing DNA Structures and Devices
% Copyright 2004
% Niles A. Pierce
% California Institute of Technology
% Version 1.0 coded 2004 (April 9,10,13,14)



clear 
clf


global nh nb rhc rbc rhelix rdh rcap ncap thcap dzb bturn brise dthb dsb dz dth ...
       cmap rdye ddye ndye majorgroove minorgroove dthgroove strutangle strutlength ...
       strutrise proptwist inclination extend5prime extend3prime extendcoil rcaphead;

     global red  pink orange yellow lightgreen turquoise darkblue lavender brown ...
    mediumgreen darkgreen darkgray black  purple  magenta  tan  palegray cornflowerblue ...
    lightblue navyblue forestgreen lawngreen shadyblue shadygreen;

material = 2;  % 1: A-DNA, 2: B-DNA
globals(material);


shade_jnk = 2;
shade_base = darkgray;
shade_t1 = 6;
shade_t2 = 1;  
shade_t3 = 1;
shade_t4 = 9;  
shade_t5 = 10; 
shade_t6 = 8; 
shade_w1 = 3;
shade_w2 = 15; 
shade_a1 = 4; 
shade_a2 = 2;
shade_d1 = 13;


% % helix 1
nbase   = 18;
xc1     = [0; 0; 0]; % start of helix axis
nc1     = [0; 0; 1]; % direction of helix axis
thc1    = -150;      % rotation around helix axis in degrees
dzc1    = 0;    % translation along helix axis from start of helix axis
shades = [shadyblue; shadygreen; shade_base];  % colors of chain a, chain b, base pairs
render = [1; 1; 1]; % render chain a, chain b, base pair struts
acap = [1 0];  % cap 5', 3' end of a chain
bcap = [0 1];  % cap 5', 3' end of b chain
adye = [0 0];  % dye 5', 3' end of a chain
bdye = [0 0];  % dye 5', 3' end of b chain
[xc2,nc2,thc2,x1a,n1a,x2a,n2a,x1b,n1b,x2b,n2b] = drawhelix(nbase,xc1,nc1,thc1,dzc1,shades,render,acap,bcap,adye,bdye);


nbase = 8;
x1 = x2a;        % start of coil
n1 = n2a;        % starting coil vector (into coil)
th1 = 90;               % starting rotation around nc1
x2 = x2b;       % end of coil
n2 = n2b;       % ending coild vector (into coil)
xp = [x1(1) .5*(x1(1)+x2(1)) x2(1)];   % points to spline between start and finish
yp = [x1(2) .5*(x1(2)+x2(2)) x2(2)];
zp = [x1(3) .5*(x1(3)+x2(3))+20 x2(3)];
bctype = [1 1];
shades = [shadyblue*ones(nbase-1,1); shadygreen; shade_base];      % shades for coil and bases shades  
render = [1; 0; 1; 0];  % render coil, first base (at x1), middle bases, last base (at x2)
ccap = [0 0];           % render cap at 5', 3' ends of coil
cdye = [0 0];           % render cap at 5', 3' ends of coil
direction = 1;          % 1: 5' end at x1, 2: 5' end at x2
helical = 0;            % make coil helical around spline
random = 0;             % randomize angle increments in helix and also base pair orientations
                        % if helical=0, still randomizes base pair rotations around chain
seed = 11;              % random seed (change to get different random chain)
% arcratio              % adjust spline points to get arcratio \approx 1 for correct chain length
[arcratio_t2] = drawcoil(nbase,x1,n1,th1,x2,n2,xp,yp,zp,bctype,helical,random,shades,render,ccap,cdye,direction,seed)

nbase = 7;
x1 = x1b + [5; 0; -17];        % start of coil
n1 = [-1; -1; 1];        % starting coil vector (into coil)
th1 = 160;               % starting rotation around nc1
x2 = x1b;       % end of coil
n2 = n1b;       % ending coild vector (into coil)
xp = [x1(1) x2(1)];   % points to spline between start and finish
yp = [x1(2) x2(2)];
zp = [x1(3) x2(3)];
bctype = [1 1];
shades = [shadygreen*ones(nbase,1); shade_base];      % shades for coil and bases shades  
render = [1; 1; 1; 0];  % render coil, first base (at x1), middle bases, last base (at x2)
ccap = [0 1];           % render cap at 5', 3' ends of coil
cdye = [0 0];           % render cap at 5', 3' ends of coil
direction = 2;          % 1: 5' end at x1, 2: 5' end at x2
helical = 1;            % make coil helical around spline
random = 1;             % randomize angle increments in helix and also base pair orientations
                        % if helical=0, still randomizes base pair rotations around chain
seed = 12;              % random seed (change to get different random chain)
% arcratio              % adjust spline points to get arcratio \approx 1 for correct chain length
[arcratio_t2] = drawcoil(nbase,x1,n1,th1,x2,n2,xp,yp,zp,bctype,helical,random,shades,render,ccap,cdye,direction,seed)


xback = [-20 20; -20 20];
yback = -[80 80; 80 80];
zback = [-25 -25; 100 100];
cback = [brown brown; brown brown];
surface(xback,yback,zback,cback,'CDataMapping','direct');


%
% rendering options
%
%set(gcf,'Position',[400 100 1050 750])
colormap(cmap);
axis off
%view(2)
view([0 1 0])
%axis([-10 310 -20 22 -10 230])
%axis([0 40 -50 50 -20 20])
axis equal
%axis image
%view(-53,8)

%view(0,0)
shading interp
lightangle(0,25)
lightangle(180,0)
set(gcf,'Renderer','zbuffer')
set(findobj(gca,'type','surface'),...
    'FaceLighting','phong',...
    'AmbientStrength',.5,'DiffuseStrength',.8,...
    'SpecularStrength',1.1,'SpecularExponent',25,...
    'BackFaceLighting','lit')

% print options
% print -dtiff -r400 walker3.tif
