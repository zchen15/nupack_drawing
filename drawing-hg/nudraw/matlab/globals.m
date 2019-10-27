%
% nudraw
% 
% A program for visualizing DNA Structures and Devices
% Copyright 2004
% Niles A. Pierce
% California Institute of Technology
% Version 1.0 coded 2004 (April 9,10,13,14)

function globals(material);
global nh nb rhc rbc rhelix rdh rcap ncap thcap dzb bturn brise dthb dsb dz dth ...
       cmap rdye ddye ndye majorgroove minorgroove dthgroove strutangle strutlength ...
       strutrise proptwist inclination extend5prime extend3prime extendcoil rcaphead;
   
   global red  pink orange yellow lightgreen turquoise darkblue lavender brown ...
    mediumgreen darkgreen darkgray black  purple  magenta  tan  palegray cornflowerblue ...
    lightblue navyblue forestgreen lawngreen lightorange shadyblue shadygreen;
   
if material == 1 % A DNA (RNA/RNA duplex or RNA/DNA duplex)
    rhelix     = 13;               % radius of double helix to outside of chain axis 
    dzb     = 2.6;		        % stacking height per base along helix
    bturn   = 11;               % bases per turn in the helix (10 is default value)   
    majorgroove = .45*dzb*bturn;          % height of major groove along helix axis in angstroms (less than .5 for narrow, deep major groove) 
    inclination = -19*pi/180;  % angle of strut relative x axis (spinning around y axis)
                                            % should tilt down from 5' (start
                                            % of A chain) to 3' (start of B
                                            % chain) strand 
    proptwist = 11.8*pi/180;   % prop twist around base pair axis in ccw direction       
elseif material == 2 % B DNA
    rhelix     = 10;               % radius of double helix to outside of chain axis 
    dzb     = 3.4;		        % stacking height per base along helix
    bturn   = 10.5;               % bases per turn in the helix (10.5 is modern value)   
    majorgroove = 22/34*(dzb*bturn);          % height of major groove along helix axis in angstroms (22 is default, 17 for symmetric debugging)
                                % take ratio of 22/34 and multiply by
                                % modern value of brise
    inclination = 1.2*pi/180;  % angle of strut relative x axis (spinning around y axis)
                                            % should tilt up from 5' (start
                                            % of A chain) to 3' (start of B
                                            % chain) (1.2 for B DNA)
    proptwist = 11.4*pi/180;   % prop twist around base pair axis in ccw direction %11.4
end

nh(1)   = 20 + 1;	        % points along chain per base %20 default
nh(2)   = 30 + 1;	        % points around chain         % 30 default
nb(1)   = 1  + 1; 	        % points along base pair
nb(2)   = 30 + 1; 	        % points around base pair   % 30 default
rhc     = 1.5;		        % radius of helix chain
rbc(1)  = .5;		        % major radius of base pair strut (elliptical: .5 default)
rbc(2)  = .2;               % minor radius of base pair strut (elliptical: .2 default)

rdh     = rhelix - rhc;     % radius of helix to center of chain
rcap    = .5;               % radius of capping corner
ncap    = 30;               % number of points along cap corners as function of r
thcap   = 55*pi/180;        % degrees of cone at 3' end of chain
rcaphead = 1.25*rhc + rcap; % radius of arrowhead at 3' end of chain

brise   = dzb*bturn;        % rise along helix axis per full turn
dthb    = 2*pi/bturn;       % stacking twist per base along helix  
dsb     = sqrt((rdh*dthb)^2 + dzb^2); % arc length along chain for one base
dz      = dzb/(nh(1)-1);    % patch increment along chain
dth     = 2*pi/(nh(2)-1);   % patch increment around chain
rdye    = 2*rhc;            % radius of dye molecule
ddye    = sqrt(rdye^2  - rhc^2);    % displacement of dye center along chain axis from last base
ndye    = 50;               % resolution of dye sphere
minorgroove = brise - majorgroove; 

extend5prime = 1*rbc(1);   % extension of chain at end so that strut doesn't hit cap
extend3prime = 3*rbc(1);   % extension of chain at end so that strust doesn't hit cap
                           % if majorgroove=minorgroove then works with
                           % both = rbc(1)
extendcoil = 1*rbc(1);  % for coil can be same at each end becauses no major/minor groove issue

%dthgroove = 2*pi*majorgroove/brise; % major groove in radians (should be > pi, with 3' trailing 5' at each slice)
%sa = sin(dthgroove-pi);
%ca = cos(dthgroove-pi);
%strutangle = pi + atan2(rdh*sa,rdh*(1+ca)); % angle of strut with origin at one end of strut
%strutlength = rdh*sqrt(sa^2 + (1+ca)^2);    % length of strut

                                            
[dthgroove,strutrise,strutlength,strutangle] = newtongroove(inclination,rdh,brise,minorgroove);
% calculate dthgroove such that dthgroove/(2*pi) * brise  - strutrise =
% majorgroove (solve via newton iteration)
% for B DNA, dthgroove > pi
% for A DNA, dthgroove < pi (major groove is smaller than minor groove)

%dthgroove = 1.5*pi;             


cmap = [109/255 216/255  45/255;    % light green
%    127/255 255/255  89/255;    % light green
        181/255 145/255 209/255;    % lavender
        1   .2857       0;    % orange
        164/255       0       0;    % red
        0/255 100/255   0/255;    % dark green
              1   .9375       0;    % yellow
              0  0.9286       1;    % turquoise
              0       0   .8571;    % dark blue
        142/255  73/255   5/255;    % brown
         28/255 206/255  40/255;    % med green
        255/255 160/255 204/255;    % pink
         0.4  0.4  0.4;    % dark gray
          0/255   0/255   0/255;    % black
         79/255   0/255 147/255;    % purple
        229/255   0/255 153/255;    % magenta
        242/255 206/255 204/255;    % tan
          .7778   .7778   .7778;    % pale gray
         % 100/255 149/255 237/255;  % cornflowerblue
           68/255 102/255 255/255;  % cornflowerblue

          135/255 206/255 235/255;   %lightblue
          0 0 128/255;               % navyblue
          34/255 139/255 34/255;    % forest green
           124/255 252/255 0 ;      % lawn green
           255/255 204/255 102/255; %light orange
           96/255 114/255 175/255;  % shadyblue
           169/255 212/255 113/255; % shadygreen
           
          ];    
   
lightgreen = 1;
lavender = 2;
orange = 3;
red = 4;
darkgreen = 5;
yellow = 6;
turquoise = 7;
darkblue = 8;
brown = 9;
mediumgreen = 10;
pink = 11;
darkgray = 12;
black = 13;
purple = 14;
magenta = 15;
tan = 16;
palegray = 17;
cornflowerblue = 18;
lightblue = 19;
navyblue = 20;
forestgreen = 21;
lawngreen = 22;
lightorange = 23;
shadyblue = 24;
shadygreen = 25;
