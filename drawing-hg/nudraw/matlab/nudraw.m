clear

global nh nb rhc rbc rhelix rdh rcap ncap thcap dzb bturn brise dthb dsb dz dth ...
       cmap rdye ddye ndye majorgroove minorgroove dthgroove strutangle strutlength ...
       strutrise proptwist inclination extend5prime extend3prime extendcoil rcaphead;

% known issues
%
% 1) nicks are now treated geometrically as if there was a helix at that
% location but instead we leave a gap of the same width. There is a
% separate "side" on either flank of the "nick" so there is no a 1-to-1
% relationship between sides of loops and coils (for those loops that are
% not part of a helix)
%
% 2) need to fix drawhelix so it can handle a single base-pair helix: done
%
% 3) introduce rotations to geometry at the end of each helix: done
%
% 4) 
% the default helix orientation is along u_z in drawhelix.m
% as a result, when it is rotated into the x-y plane at different points
% around a multiloop, the strands don't all enter the multiloop with the
% same orientations (as they would if the default helix axis was along u_x
% or u_y) -- correct for this by spinning around the axis of the helix a
% different amount depending on theta in x-y plane
%
% 5) should I change the helix width in loops to be strutwidth? probably not 
%
% 6) need to fix Newton solver so can handle hairpins of size three: done
%
% 7) need to look into sizing of loops depending on whether they flip or
% not to reduce the amount of stretching (and depending on how close the 
% ends of the entering helix is to the plane) *****
% 
% 8) loop special cases
%    a) need to make it so exterior loops can be of size zero or other
%    small number: done
%    b) need to make it so exterior loops between two helices stack but 3
%    or more look like multiloop with nick !!!
%    c) need to make it so 1-base bulges stack (based on energy model
%    treatment including the stacking energy) !!!
% 9) redid the loop ordering so iloop=1 corresponds to the loop that starts
% with side = 1 at ibase = 1: done
%
% 10) define blunt exterior loop geometry so that two disjoint helices appear
% in a sensible way (instead of currently on top of eachother) Note that the exterior
% loop geometry is not used for anything at the end of leaf helices that have blunt
% ends
%

 s = '..(((...(((.....)))..(((....))).)))....'; % classic test with multiloop
%s = '((((((((((((((|))))))))))))))';
% s = '..(((...(((.....)))(((....))).)))....(.(....).)...((.....))...'; %
% extended test with large multistem external loop

s = '..((((((....((((((....))))))((((((....))))))....))))))..';
s = '..(((((((((....((((((....))))))((((((....))))))....)))))))))..';

%s = '...(((...((((..(((.(((|((((....)))))))))).(((.....)))...))))..(((......)))..)))...((((..(((....)))..))))...';
%s = '...((((((...((((((...(((((((((((....))))))))))).((((((.....))))))...))))))...((((((......))))))...))))))...((((((((..(((((((((((....)))))))))))..))))))))...';
%s = '...((((((((...)))).))))..';

% s = '(|)((((|))))'; % test one-pair helix

% test case
s = '(((((((((...((((((..((((((((|))))))))|..))))))(((((((((((((..((((((.|.))))))...(((((((((((.|)))))))))))..)))))))))))))...)))))))))....((((((((((((((((..((((((((((..|...))))))))))((((((((((((((((...|..)))))....)))))))))))((((((((((((...|..))))))))))))..))))))))))))))))...((((((((...)))))))).(((((((..|...)))))))'; 
%s = '(((((((((...((((((..((((((((...))))))))..))))))(((((((((((((..((((((....))))))...(((((((((((....)))))))))))..)))))))))))))...)))))))))....((((((((((((((((..((((((((((.....))))))))))((((((((((((((((.....)))))....)))))))))))((((((((((((.....))))))))))))..))))))))))))))))...((((((((...)))))))).(((((((.....)))))))'; 
%
% the planar geometry for this output can be tested using 
%"diff geocheck geocheck_safe" with randomendcoils = 0;

%s = '...((((((((((((((.(((((((...............))))))).)))))))(((((((....))))))).....(((((((((...(((((((((........(((((((((.....((((((((.......))))))))))))))))))))))))))......)))))))))...(((((((((((......)))))))))))..(((((((((((((((.........))))))))))))))))))))))(((((((((((((.....)))))))))))))(((((((((.......((((((((.......)))))))).(((((((((((((((.......((((((((((((((((((((....)))))))))))..........((((((((...((((((((((.............))))))))))(((((((......)))))))((((((((.....)))))))).......(((((((((((....)))))))))))..))))))))...)))))))))))))))))))))))))))))))))..((((((((((((.(((((((((((....)))))))))))..(((((((((((.(((((((((((.........(((((((((....................))))))))))))))))))))...)))))))))))..(((((((((((((((........))))))))))))))).(((((((((((....((((((((((.....))))))))))...)))))))))))))))))))))))....((((((((..((((((((.......))))))))...((((((((((((....))))))))))))(((((((((....))))))))).(((((((((....))))))))).((((((((..........((((((((((........))))))))))))))))))))))))))...(((((((((....)))))))))..';
%s = '...(..((((..((((((((....((((........((((((.(((.(....(((((.((((((......)))))))))))..).)))))).)))..((((.((((((..(((.(((((.(((.(((((........((((.((((((......)))))).))))........((((((..((((......))))..))))))........(((((((((.....((((...(((.((((((...........(((.....)))............)))))).))).))))..)))))).))).))))).))).)))))............))).))))))))))..)))).))).)))).)..))))..)((((((...((((((((((((......((((.(((.((((.((...)).))))..(((((((.(((((......)))))......((((((((.....))).))))).))))))).((...((((.(((((((......(((...((((((((((((.((((((((((((.(..(....(((..(((...)))..)))....)..).))))))))..........(((((.........)))))..))))))))))))))))...)))((((.((((..........((((((..(((((((.......((((.....(((((..(((((...((((((..(((((((..((((((((...........((((.....))))...........)))))))).....((((((...)))))))))))))...))))))))))))))))...))))..(((((((.....))))))).((((((...((((.......((((((((((((..(((..(((((.....((((....).))).....)))))..)))))))..))).)))))...))))..))))))....)))))))..)))))).)))).)))).)))))))))))..)).)))))))....))))))......))))))....................(((((((..(((((((((((..((.(((((....))))).))..))).))..(((....)))...)))))).......((((((......((((((((.((((((.............((((...(((.....((((.......))))..((....)).)))...)))).(((((....)))))......)))))).)))))))).....((((((....))))))......))))))....))).)))).....(((((..(((((((..((((................))))..(((.(((((.((((......(.((((....)))))..))))..((((((((((((((.........))))))))))...))))..)))))..)))..............))))))).))))).))))))...'; 
%  s = '...(((((((((((...(((...(((....))))))...)))))))))))......';
%s = '...((((((((((((((((...((((((..((((((.....))))))..))))))...((((((((((((((((....))))))))))))))))...))))))))))))))))......';

% s = '.....((((((...(((((((((((....)))))))))))...((((((...........))))))...))))))..';
%
%
%s = '((((((((((((((((....(((((((((((((((|)))))))))))))))...(((((((((((((|)))))))))))))...))))))))))))))))';
%s = '((((|((((((((((((((((((..................|.................)))))|)))))).......|..............((((((|((((..|..))))))))))..)))))))))))';
% s = '..(((((((((...((((((..((((((((.....))))))))....))))))...(((((((((((((..((((((.....))))))...(((((((((((....)))))))))))..)))))))))))))...)))))))))....'; %((((((((((((((((..((((((((((.....))))))))))((((((((((((((((.....)))))....)))))))))))((((((((((((.....))))))))))))..))))))))))))))))...((((((((.....)))))))).(((((((.....)))))))..'; 
%s = '..(((....)))(((...)))..';
%s = '(((|)))(((...)))';
%s = '(((((((((((((((((|(((((((((((((((|(((((((((((((((((........................|))))))))))))))|))))))))|)))))))))))))))))))))))))))';
%s = '((((.(((((...((((((...)))))))))))(((|)))))))';
%s = '(((|((|((|((((...(((((..................))))))))))))))))';
%s = '((((|))))';
%s = '((((|((|((|((((((((((((|((|((((|((((((((((|)))))))))))))))))|)))))))))))))))))))';
%s = '((((((((((((((((|((((((((((((((((((((((((((((((((((((|((((((((((((((((((..................|))))))))))))))))))))))))))))))))))))|((((((((((((((((|..........))))))))))))))))))))))))))))))))))))))))))))))))))..........';
%s = '......((((((((((((((((((...|...))))))))))))))))))((((((((((((((((((.|.....))))))))))))))))))((((((((((((((((((...|...))))))))))))))))))';
%s = '((((((((((((((((|((((((((((((((((((((((((((((((((((((|((((((((((((((((((..................|))))))))))))))))))))))))))))))))))))|((((((((((((((((|..........))))))))))))))))))))))))))))))))))))))))))))))))))..........';
%s = '((((((((((((((((((((((((|((((((((((((((((((((((((........................|))))))))))))))))))))))))))))))))))))))))))))))))';
%s = '(((((((((((((((....)))))))))))))))';

tic
disp(sprintf('Logical Definition of Loops'))

colorscheme = 1; % 1: by strands, 2: by coil or helix
smartrot = 1;  % 1: rotate leaves of tree for aesthetics (should always be on)
planar = 1; % 1: keep structure in the plane, 0: let loops rotate out of plane
randomendcoils = 1;  % make coils random on either side of nicks
stacknickedhelices = 1; % stack nicked helices
material = 2; % 1: A-DNA, 2: B-DNA
globals(material);

leftpair = '(';
nopair = '.';
rightpair = ')';
nick = '|';

nsymbol = length(s);

%
% set up logical representation of strand with nick, base and strand
% numbers
%
ibase = 1;
iunit = 1;
istrand = 1;
unit(iunit).type    = -5;
unit(iunit).base    = ibase;
unit(iunit).strand  = istrand;
iunit = iunit + 1;
for isymbol = 1:nsymbol
    if s(isymbol) == nick
        unit(iunit).type    = -3; % 3' end
        unit(iunit).base    = ibase - 1;
        unit(iunit).strand  = istrand;
        iunit               = iunit + 1;
        istrand             = istrand + 1;
        unit(iunit).type    = -5; % 5' end
        unit(iunit).base    = ibase;
        unit(iunit).strand  = istrand;
        iunit               = iunit + 1;        
    else
        unit(iunit).type    = s(isymbol);
        unit(iunit).base    = ibase;
        unit(iunit).strand  = istrand;
        ibase               = ibase + 1;
        iunit               = iunit + 1;
    end
end
unit(iunit).type    = -3;
unit(iunit).base    = ibase - 1;
unit(iunit).strand  = istrand;
nunit = iunit;
nbase = ibase - 1;
nstrand = istrand;


%
% identify partners for all logical pairs (strand ends or base pairs)
% used to traverse loops
%
istack = 0;
iprev = nunit;
for iunit = 1:nunit
    if unit(iunit).type     == -3 
        unit(iprev).next    = iunit;
        iprev               = iunit;
    elseif unit(iunit).type == -5 
        unit(iunit).pair   = iprev;
        unit(iprev).pair   = iunit;
        unit(iprev).next    = iunit;
        iprev               = iunit;
    elseif unit(iunit).type == leftpair    
        istack              = istack + 1;
        stack(istack)       = iunit;
        unit(iprev).next    = iunit;
        iprev               = iunit; % previous paired base along strand 
    elseif unit(iunit).type == rightpair   
        unit(iunit).pair   = stack(istack);
        unit(stack(istack)).pair  = iunit;        
        unit(iprev).next    = iunit;
        iprev               = iunit;
        istack              = istack - 1;
    end  
end

% for iunit = 1:nunit
%    unit(iunit)
% end

%
% traverse all loops starting from a left base pair 
% (except last exterior loop, which starts from a right base)
%
% save pointer to previous loop number and loop side
% at left base of pair at point of entry to new loop
%
% starting from a left base in a pair (ibase) and traversing a loop in the
% clockwise direction, jside = 1 indexes the first encountered single-stranded
% region and the first encountered stem. For any given loop, jside = nside 
% indexes the last single-stranded region and the last stem (of which ibase
% is part of the closing pair)
%
% nicks are treated logically as pairs between the base at the 3' end of one strand
% and the 5' end of the next strand (by the variable base(ibase).next and
% base(ibase).ipair) but there is extra logic in the variable
% base(ibase).nick to keep track of the fact that they actually represent a
% nick
%
% the first loop starts at ibase = 1 with iside = 1 representing the first
% encountered single stranded region and the subsequent helix

iunit = 1; % use i index for moving along chain to find left unit in pair to start each loop
iloop = 1;
while xor(iunit == 1,iloop > 1) 
    if unit(iunit).type == leftpair | iunit == 1 % special case for first exterior loop
        jstart  = iunit; % use j index for traversing loop
        junit   = jstart;
        jside   = 1;
        while xor(junit == jstart,jside > 1)
            jnext                           = unit(junit).next;
            loop(iloop).sideunit(jside)     = junit;            % unit right before single stranded region for each side 
            loop(iloop).sidebase(jside)     = unit(junit).base; % base right before single stranded region for each side
            loop(iloop).sidenbase(jside)   = unit(jnext).base - unit(junit).base - 1;
            
            if unit(junit).type             == -5
                loop(iloop).nick(:,jside)   = [-5; 0];
            elseif unit(jnext).type         == -3
                loop(iloop).nick(:,jside)   = [0; -3];
            else
                loop(iloop).nick(:,jside)   = [0; 0];
            end
            loop(iloop).strand(jside) = unit(junit).strand;

            if unit(jnext).type  == leftpair 
                unit(jnext).prevloop = iloop; % use left base in pair to store prev loop 
                unit(jnext).prevside = jside; % use left base in pair to store prev side
            end
            jside = jside + 1;
            junit = unit(jnext).pair;
        end
        loop(iloop).nside = jside-1;         
        iloop = iloop + 1;
    end
    iunit = unit(iunit).next;
end
nloop = iloop-1;



%
% make linked list of loops and sides using saved pointers
%
for iloop = 2:nloop
    iunit = loop(iloop).sideunit(1);
    nside = loop(iloop).nside;
    prevloop = unit(iunit).prevloop;
    prevside = unit(iunit).prevside;
    loop(iloop).toloop(nside) = prevloop; % link to previous loop
    loop(iloop).toside(nside) = prevside; % link to previous side
    loop(prevloop).toloop(prevside) = iloop;  % link from previous loop
    loop(prevloop).toside(prevside) = nside;  % link from previous side
   
    %
    % make linked list of sides within each loop (just for convenience when 
    % traversing from side = 1 to side = nside or back)
    %
    for iside = 1:nside
        prevside = iside - 1;
        if prevside < 1
            prevside = nside;
        end
        nextside = iside + 1; 
        if nextside > nside
            nextside = 1;
        end
        loop(iloop).prevside(iside) = prevside;
        loop(iloop).nextside(iside) = nextside;
         %[iloop iside nside;0 prevside nextside]
    end
end



toc
tic
disp(sprintf('\nGeometric Definition of Loops'))

%
% define loop geometries in terms of center point and direction of all
% stems on 
% each loop
%
% convention: normal points out from current loop at each base pair 
%

nside = 2;
sidelength(1:2) = 2*rdh;
ds  = dzb;
nds = 2;    
[stacked.radius stacked.sideangle stacked.ang_ds] = newtonloop(nside,sidelength,ds,nds); % get helix parameters

nicklength = 2*rdh;
safe_count = 1;
for iloop = 1:nloop 
    loop(iloop).nickedhelix = 0; % put in zeros for all loops and overwrite for nicked helices below
    nside       = loop(iloop).nside; % number of stems in loop
   


    is_stacked      = nside == 2 && loop(iloop).sidenbase(1) == 0 && loop(iloop).sidenbase(2) == 0;
    is_smallhairpin = nside == 1 && loop(iloop).sidenbase(1) == 3;   % hairpin with 3 bases in loop
    is_blunthelix   = nside == 2 && loop(iloop).sidenbase(1) == -1 && loop(iloop).sidenbase(2) == -1;  % exterior loop of size zero
    is_nickedhelix  = nside == 3 && sum(loop(iloop).sidenbase) == -2;
    
   
    
    for iside = 1:nside
       sidenbase(iside) = loop(iloop).sidenbase(iside); 
    end
    
    if is_stacked
        r           = stacked.radius;
        ang_ds      = stacked.ang_ds;     
        sideangle   = stacked.sideangle;
        sidelength  = [2*rdh 2*rdh];
    elseif is_smallhairpin
        r = rdh;
        ang_ds = pi/4; % for hairpin of length three array bases on half circle centered at center of paired base
        sideangle(1) = pi; % angle between ends of helix
        sidelength(1) = 2*rdh;
    elseif is_blunthelix % external loop of size zero
        r = rdh;
        ang_ds = 0;    % no angle between bases along strand
        sideangle(1) = pi;   % angle between strand ends
        sideangle(2) = pi;
        sidelength   = [2*rdh 2*rdh];
    elseif stacknickedhelices && is_nickedhelix
         r           = stacked.radius;
         ang_ds      = stacked.ang_ds;     
         if loop(iloop).sidenbase(1) == -1
            sideangle   = [(pi-stacked.sideangle(2)) stacked.sideangle(1) stacked.sideangle(2) ];
            sidelength  = [dzb 2*rdh 2*rdh];   
         elseif loop(iloop).sidenbase(2) == -1
            sideangle   = [stacked.sideangle(1) (pi-stacked.sideangle(2)) stacked.sideangle(2) ];
            sidelength  = [2*rdh dzb 2*rdh];     
         end   
         loop(iloop).nickedhelix = 1;
    else
        nds         = 0; % number of gaps 
        for iside = 1:nside
            nds                 = nds + loop(iloop).sidenbase(iside) + 1;
            if loop(iloop).nick(2,iside) == -3    % space for nick
                loop(iloop).geo(iside).length = nicklength;
            else                            % space for stem
                loop(iloop).geo(iside).length = 2*rdh;
            end 
            sidelength(iside)   = loop(iloop).geo(iside).length;
        end
        ds      = dsb; %use arc length between bases in single-stranded regions
         
        [r sideangle ang_ds f] = newtonloop(nside,sidelength,ds,nds);
    end
    

    if iloop == 1 % taking starting point from global origin
        isideref = 1;
        loop(iloop).geo(isideref).xc = [0; 0; 0]; % starting position for secondary structure drawing
        loop(iloop).geo(isideref).nc = [0; 1; 0]; % starting direction for secondary structure drawing
        jsidestart = 2;
        jsidestop = nside;
    else % take starting point from previous loop       
        isideref = nside;
        prevloop = loop(iloop).toloop(isideref);
        prevside = loop(iloop).toside(isideref);
        loop(iloop).geo(isideref).xc =  loop(prevloop).geo(prevside).xc;
        loop(iloop).geo(isideref).nc = -loop(prevloop).geo(prevside).nc;
        jsidestart = 1;
        jsidestop = nside-1;
    end
    
    % assign stem position and direction for jside stem (by
    % convention this follows the jside single-stranded region)
    %
    xc = loop(iloop).geo(isideref).xc;
    nc = loop(iloop).geo(isideref).nc;
    
    loopcenter = xc - nc*sqrt(r^2 - (sidelength(isideref)/2)^2);
    loop(iloop).center = loopcenter;
    loop(iloop).radius = r;
    jprev = isideref;
    for jside = jsidestart:jsidestop   
        jgap = loop(iloop).sidenbase(jside) + 1;
        jang    = ang_ds*jgap + .5*(sideangle(jside) + sideangle(jprev));
        
        ca      = cos(jang);   % rotation matrix is CW around z axis
        sa      = sin(jang); 
        rmat    = ca*eye(3) + (1-ca)*[0 0 0; 0 0 0; 0 0 1] + sa*[0 1 0; -1 0 0; 0 0 0];
        nc      = rmat*nc;
        xc      = loopcenter + nc*sqrt(r^2 - (sidelength(jside)/2)^2);
        loop(iloop).geo(jside).nc = nc;
        loop(iloop).geo(jside).xc = xc;
        jprev   = jside;
        
        stem_safe(safe_count,:) = [iloop jside xc' nc']; 
        safe_count = safe_count + 1;
    end    

%      xc_ref = loop(iloop).geo(isideref).xc;
%      nc_ref = loop(iloop).geo(isideref).nc;
%      center = xc_ref - nc_ref*sqrt(r^2 - rdh^2);
%      if jsidestop >= jsidestart
%      [xcstore ncstore] = placestems(jsidestart,jsidestop,xc_ref,nc_ref,center,r,ang_ds,sideangle,sidenbase);
%      end  
%        
%      loop(iloop).center = center;
%      loop(iloop).radius = r;
%      for jside = jsidestart:jsidestop  
%         loop(iloop).geo(jside).nc = ncstore(:,jside);
%         loop(iloop).geo(jside).xc = xcstore(:,jside);
%         stem_safe(safe_count,:) = [iloop jside xcstore(:,jside)' ncstore(:,jside)']; 
%         safe_count = safe_count + 1;
%      end
     
     
    %
    % assign base positions for jside single-stranded region
    %
    jprev = nside;
    for jside = 1:nside
        jbasestart = loop(iloop).sidebase(jside);
        jbasestop = loop(iloop).sidebase(jside) + loop(iloop).sidenbase(jside);
        if loop(iloop).nick(2,jside) == -3  
            jbasestop = jbasestop + 1; % handle base at nick that would not otherwise get assigned an x value
        end

        nc = loop(iloop).geo(jprev).nc;
        ang = .5*sideangle(jprev);
        for jbase = jbasestart:jbasestop %loop over paired base and following single-stranded region               
            ca      = cos(ang);   % rotate CW by ang
            sa      = sin(ang);   % rotation matrix is CW around z axis
            rmat    = ca*eye(3) + (1-ca)*[0 0 0; 0 0 0; 0 0 1] + sa*[0 1 0; -1 0 0; 0 0 0];
            ncstep = rmat*nc;
            base(jbase).x = loop(iloop).center + ncstep*r;
            ang = ang + ang_ds;
        end
        jprev = jside;
    end
end    

figure(1)
clf
hold on;

for iloop = 1:nloop
   nside = loop(iloop).nside; % number of stems in loop
   for jside = 1:nside 
       xc = loop(iloop).geo(jside).xc;
       nc = loop(iloop).geo(jside).nc;
       plot(xc(1),xc(2),'.')
   end
   th = linspace(0,2*pi,20);
   xcirc = loop(iloop).center(1) + loop(iloop).radius*cos(th);
   ycirc = loop(iloop).center(2) + loop(iloop).radius*sin(th);
   plot(xcirc,ycirc,'g-')
   text(loop(iloop).center(1),loop(iloop).center(2),sprintf('L%d',iloop));
end

for ibase = 1:nbase
    text(base(ibase).x(1),base(ibase).x(2),sprintf('%d',ibase));
    base_safe(ibase,:) = [ibase base(ibase).x'];
end
%plot(base.x(1),base.x(2),'r*')

save geocheck stem_safe base_safe -ascii

axis('equal')

toc
tic
disp(sprintf('\nDefine Helices and Coils'))

%
% define helices with links to previous loop number and side number
% and next loop number and side number
% define coil segments with link to relevant loop number and side number
%
%
iloop = 2;
ihelix = 1;
while iloop <= nloop
%
%   define helices
%   each helix has a pointer to the previous loop before the helix
%   and the next loop after the helix as well as pointers to the
%   appropriate side in these neighboring loops: these pointers
%   are used to link helix ends to coil ends below
%
    helixlength = 1;  
    firstloop = iloop;
    while (loop(iloop).nside == 2 && loop(iloop).sidenbase(1) == 0 && loop(iloop).sidenbase(2) == 0 ...
           && loop(iloop).nick(2,1) ~= -3)
        helixlength = helixlength + 1;
        iloop = iloop + 1; % because loops in a helix are always numbered consecutively
    end 
    
    
        helix(ihelix).npair = helixlength;
        
        firstside = loop(firstloop).nside; 
        prevloop = loop(firstloop).toloop(firstside);
        prevside = loop(firstloop).toside(firstside);
        helix(ihelix).prevloop  = prevloop;
        helix(ihelix).prevside  = prevside;
        helix(ihelix).xc     =  loop(firstloop).geo(firstside).xc; % coord of helix axis start
        helix(ihelix).nc     = -loop(firstloop).geo(firstside).nc; % helix axis at start (into helix)
        helix(ihelix).thc    = 180/pi*(pi - strutangle) ; %rotate so base-pair struct is parallel to plane of loop it borders
                % reversed sign of thc relative to nupack6 to get first
                % basepair strut into the plane of the loop it exits from
                % (works for A-DNA or B-DNA)
   
        
        if ~smartrot
            helix(ihelix).thc = helix(ihelix).thc  + 90 + 180/pi*atan2(helix(ihelix).nc(2),helix(ihelix).nc(1)); 
        end
                                % rotate so that projection from canonical 
                                % u_z orientation creates helix that joins the loop
                                % with the same rotation regardless of helix axis
                                
        %
        % make pointers from loop sides to helix numbers: from root to
        % leaves (don't need pointers from loops to helix entering the loop
        % -- but could add below using nextloop, nextside
        %
        loop(prevloop).tohelix(prevside) = ihelix;                                         
        
        if helixlength == 1 %special treatment for isolated base pairs
            helix(ihelix).nextloop = firstloop;
            helix(ihelix).nextside = firstside;
        else                % regular helics with two or more stacked pairs
            lastloop  = iloop - 1;
            lastside = 1;
            helix(ihelix).nextloop  = loop(lastloop).toloop(lastside);
            helix(ihelix).nextside  = loop(lastloop).toside(lastside);
        end
        
        if stacknickedhelices && loop(prevloop).nickedhelix
           helix(ihelix).thc = helix(ihelix-1).thc + 180/pi*dthb*helix(ihelix-1).npair;
        end % this takes advantage of fact that helices connected by a nick are consecutive in the helix array
        % pass on rotation info from the previous helix (can't do this by
        % defining reference angles for the nicked loops because that would
        % pass on non-2D orientations to the subsequent loops)
    
        ihelix = ihelix + 1;
    iloop = iloop + 1; 
end
nhelix = ihelix - 1; % substract last increment 

%
%   define coils
%
iloop = 1;
icoil = 1;
while iloop <= nloop
    if ~(loop(iloop).nside == 2 && loop(iloop).sidenbase(1) == 0 && loop(iloop).sidenbase(2) == 0) ...
        | loop(iloop).nick(2,1) == -3 % if not part of helix
        for jside = 1:loop(iloop).nside
            coil(icoil).nbase = loop(iloop).sidenbase(jside) + 2; % including end pts
            coil(icoil).loopnum = iloop; % check
            coil(icoil).loopside = jside; % check
            loop(iloop).coilnum(jside) = icoil; 
                
            coil(icoil).strand = loop(iloop).strand(jside);
            coil(icoil).nick = loop(iloop).nick(:,jside);
            if coil(icoil).nick(1) == -5
                coil(icoil).tohelix(1,:) = [0 -5]; % could change zero to istrand if need info later
            end                                    % or just not use tohelix for strand ends and rely on coil.nick
            if coil(icoil).nick(2) == -3
                coil(icoil).tohelix(2,:) = [0 -3];
            end
                
            if loop(iloop).nick(1,jside) == -5
                startval = 0;
            else
                startval = 1;
            end
            if loop(iloop).nick(2,jside) == -3
                stopval = loop(iloop).sidenbase(jside) + 1;
            else
                stopval = loop(iloop).sidenbase(jside);
            end

            jcount = 1;
            for jval = startval:stopval % store spline pts for coil, including end point if a nick
                jbase = jval + loop(iloop).sidebase(jside);
                coil(icoil).xsplinepts(:,jcount) = base(jbase).x;
                jcount = jcount + 1;
            end
            coil(icoil).nsplinepts = jcount-1;             
            icoil = icoil + 1;        
        end
    end
    iloop = iloop + 1; 
end
ncoil = icoil - 1;
%
% drawing reality checks
%
for icoil = 1:ncoil % draw coil bases for reality check
    if coil(icoil).nsplinepts
        for ipt = 1:coil(icoil).nsplinepts
            plot(coil(icoil).xsplinepts(1,ipt),coil(icoil).xsplinepts(2,ipt),'m*')
        end
        text(coil(icoil).xsplinepts(1,1)+3,coil(icoil).xsplinepts(2,1),sprintf('c%d',icoil));
    end
end

for ihelix = 1:nhelix
   plot(helix(ihelix).xc(1),helix(ihelix).xc(2),'c+')
   text(helix(ihelix).xc(1)+2,helix(ihelix).xc(2),sprintf('H%d',ihelix)) 
end


%
% make linked list between helix strand ends and coil strand ends
% loop of helices
% follow pointers helix -> end loops -> loop sides -> coil
% 
% helix definitions in drawing program:
% a chain runs from 5' -> 3' (indices 1 and 2, respectively, in drawing program)
% b chain rums from 3' -> 5' (indices 1 and 2, respectively, in drawing program)
%
% helix definitions in automatic structure generation program:
% a1 = 1, a2 = 2, b1 = 3, b2 = 4
%
for ihelix = 1:nhelix
    prevloop  = helix(ihelix).prevloop;
    side1     = helix(ihelix).prevside;
    ic1       = loop(prevloop).coilnum(side1);
    
    side3     = side1 + 1;
    nside     = loop(prevloop).nside;
    if side3 > nside
        side3 = 1;
    end  
    ic3       = loop(prevloop).coilnum(side3);
   
    nextloop  = helix(ihelix).nextloop;
    side4     = helix(ihelix).nextside;
    ic4       = loop(nextloop).coilnum(side4);
    
    side2     = side4 + 1;
    nside     = loop(nextloop).nside;
    if side2  > nside
        side2 = 1;
    end
    ic2       = loop(nextloop).coilnum(side2);
%    
     helix(ihelix).tocoil = [ic1 2; ic2 1; ic3 1; ic4 2]; % coil number and position
     helix(ihelix).strand = [coil(ic1).strand coil(ic3).strand]; % strand number for chains a and b
%     helix(ihelix).strand = [astrand bstrand]; % strand number for chains a and b
     
     coil(ic1).tohelix(2,:) = [ihelix 1]; % helix number and position (1->4)
     coil(ic2).tohelix(1,:) = [ihelix 2];
     coil(ic3).tohelix(1,:) = [ihelix 3];
     coil(ic4).tohelix(2,:) = [ihelix 4];
     %
     % if coil is at strand end and has no length then just draw cap on
     % relevant helix strand end and don't draw coil
     %
     helix(ihelix).cap = [0; 0; 0; 0];
     if coil(ic1).nick(1) == -5 && coil(ic1).nbase == 1 
        helix(ihelix).cap(1) = -5;
     end
     if coil(ic2).nick(2) == -3 && coil(ic2).nbase == 1
        helix(ihelix).cap(2) = -3;
     end
     if coil(ic3).nick(2) == -3 && coil(ic3).nbase == 1
        helix(ihelix).cap(3) = -3;
     end
     if coil(ic4).nick(1) == -5 && coil(ic4).nbase == 1
        helix(ihelix).cap(4) = -5;
     end
        
end

toc
tic
disp(sprintf('\nPerform Leaf to Root Aesthetic Rotations'))

%
% need pointers from each loopside to first base to store coordinates: done
% need pointers from coil to loopside to retreive base coordinates: done
% calculate base coordinates base.x: done
% for each coil take internal spline points from stored base coordinates: done
% make linked list from coil ends to helix strand ends and back: done
%
% need to fix drawhelix so that can draw a helix with a single base pair
%

%
% introduce rotations so that the loop following a helix is in the plane
% rotated to match the end of the helix
%
% actually, decided that for visual clarity it's desirable to keep all the
% loops in the plane. So, for the loop exiting each helix, rotate loop on
% helix axis either zero or pi radians depending on exit spin of helix
% 
% implementation: start with iloop = nloop (the external loop -- with
% rotation reference angle = 0) so all helixes exiting that loop have
% refangle = 0.  
%
% for each loop, loop over exiting helices making a stack of 
% loops at the end of these helices, and assigning accumlated
% refangles and reforigins to these loops
%
% continue until nothing in the stack
%
% then loop over coils and helices and use pointers to loops to do transformations 
% based on refangles and reforigins 
% for coils, refloop is coil(icoil).loopnum
% for helices, refloop is helix(ihelix).prevloop
%
% it is critical to do these transformations in the reverse order that the
% transformations were encountered by the above process: this way we work
% from the leaf of the tree to the trunk, so that the reference origin and
% axis is valid for each transformation at the time it is applied
%

globalangle = 0;                    % can rotate around z axis to change orientation (in radians)
loop(1).ref.angle  = globalangle;    
loop(1).ref.origin = [0; 0; 0]; % keep the origin the same as helix(1).xc(:,1)
loop(1).ref.axis   = [0; 0; 1]; % rotations work for z axis but not for arbitrary axes
loop(1).ref.num    = 1;         % number of reference rotations

nstack = 1;
loopstack(nstack) = 1;
while nstack > 0
    prevloop = loopstack(1); % always work from the front of the stack
    loopstack = loopstack(2:nstack); % shift stack entries
    nstack = nstack - 1;
    
    nside = loop(prevloop).nside;   
    
    if prevloop == 1; 
        maxside = nside; 
    else
        maxside = nside-1; 
    end
    
    for iside = 1:maxside % loop over exiting helices
        if loop(prevloop).nick(2,iside) ~= -3 % don't process nicks 
            ihelix = loop(prevloop).tohelix(iside);        
            nextloop = helix(ihelix).nextloop; 
            loop(nextloop).ref = loop(prevloop).ref; %inherit rotation information from next higher loop
            % add new rotation from intervening helix: angle= thc(2)-thc(1)
            % origin = xc(:,2); axis = nc(:,2)
            nrefnext = loop(nextloop).ref.num + 1;        
            xc2      =  dzb*(helix(ihelix).npair - 1)*helix(ihelix).nc(:,1) + helix(ihelix).xc(:,1);  
            binangle = dthb*(helix(ihelix).npair - 1) ...
                     + helix(ihelix).thc*pi/180 - (pi - strutangle); % second and third terms are for nicked helices
           % because need to know rotation since last planar entry point to
           % a helix so pass on rotation infor acquired from previous
           % helices (the third term just cancels out the initial entry
           % angle of first helix since the binning below was developed
           % based on the first term alone)
           %
           % there was a bug in binangle in nudraw6 because helix.thc was
           % in degrees instead of radians so binning was bogus
           
           if planar
               
            modangle = mod(binangle,2*pi);
        
             if stacknickedhelices && loop(nextloop).nickedhelix 
                   binangle = 0; % no flipping if in a nick (just use rotation defined earlier and stored in helix.thc
             else
               %%%if ~(loop(nextloop).nside == 1)  && loop(nextloop).sidenbase == 3)
                if material == 2  % B-DNA
                    if modangle >= 0 && modangle < pi % bin rotations based on exit angle from helix -- keep the loops in the plane for clarity
                        binangle = pi;
                    else
                        binangle = 0;
                    end
                elseif material == 1 % A-DNA
                    if modangle > pi/8 && modangle <= 9*pi/8 % bin rotations based on exit angle from helix -- keep the loops in the plane for clarity
                        binangle = pi;
                    else
                        binangle = 0;
                    end                
                end
             %%%end % test
             end
           end % end planar
            
            
            loop(nextloop).ref.angle(nrefnext)    = binangle;
            loop(nextloop).ref.origin(:,nrefnext) = xc2;
            loop(nextloop).ref.axis(:,nrefnext)   = helix(ihelix).nc(:,1); % nc2 = nc1
            loop(nextloop).ref.num                = nrefnext;
            nstack = nstack + 1;
            loopstack(nstack) = nextloop;
        end
    end
end


if smartrot
%
% loop over helices and do transformations on xc1, nc1, thc1
% it's really important to loop from leaves to trunk of tree so origins and
% axes remain valid during each step
%
 for ihelix = 1:nhelix
   irefloop = helix(ihelix).prevloop;
   nref = loop(irefloop).ref.num;
   for iref = nref:-1:1  % leaf to trunk
       refangle  = loop(irefloop).ref.angle(iref);
       reforigin = loop(irefloop).ref.origin(:,iref);
       refaxis   = loop(irefloop).ref.axis(:,iref);
       
   
       % this rotation matrix is clockwise so put in minus sign on angle to get ccw
       rmat    = cos(-refangle)*eye(3)  ...
               + (1-cos(-refangle))*refaxis*refaxis' ...
               +    sin(-refangle) *[0                 refaxis(3)           -refaxis(2); ...
                                    -refaxis(3)          0                   refaxis(1); ...
                                     refaxis(2)         -refaxis(1)            0];   

       nc1 = helix(ihelix).nc(:,1);
       xc1 = helix(ihelix).xc(:,1);
       thc1 = helix(ihelix).thc(1);
       nc1 = rmat*nc1;
       xc1 = rmat*(xc1 - reforigin) + reforigin;
       thc1 = thc1 + refangle*180/pi; %spin helix on axis so exit angle stays same relative to that of entering helix 
       helix(ihelix).nc(:,1) = nc1;
       helix(ihelix).xc(:,1) = xc1;
       helix(ihelix).thc(1)  = thc1;
   end
   
   %un_dot_uz = [0 0 1]*helix(ihelix).nc(:,1); % dot product between z axis and helix axis 
   helix(ihelix).thc(1)  = helix(ihelix).thc(1)  ... 
                          + 90 + 180/pi*(atan2(helix(ihelix).nc(2),helix(ihelix).nc(1)) - globalangle); 
  % rotate so that projection from canonical 
  % u_z orientation creates helix that joins the loop
  % with the same rotation regardless of
  % helix axis -- this has nothing to do with the coils since they start
  % out defined in the x-y plane (so no project issue from u_z)
 end


%
% loop over coils and do transformations on xsplinespts: loop from leaves
% to trunk so origins and axes remain valid as you work
%
for icoil = 1:ncoil
   irefloop = coil(icoil).loopnum;
   nref = loop(irefloop).ref.num;
   nsplinepts = coil(icoil).nsplinepts; 
   for iref = nref:-1:1; %leaf to trunk
       refangle  = loop(irefloop).ref.angle(iref);
       reforigin = loop(irefloop).ref.origin(:,iref);
       refaxis   = loop(irefloop).ref.axis(:,iref);
   
       % this rotation matrix is clockwise so put in minus sign to get ccw
       rmat    = cos(-refangle)*eye(3)  ...
               + (1-cos(-refangle))*refaxis*refaxis' ...
               +    sin(-refangle) *[0                  refaxis(3)           -refaxis(2); ...
                                     -refaxis(3)          0                   refaxis(1); ...
                                      refaxis(2)         -refaxis(1)            0]; 
                                      
                                     
       if nsplinepts > 0        
          for ipt = 1:nsplinepts
              xsplinepts = coil(icoil).xsplinepts(:,ipt);
              xsplinepts = rmat*(xsplinepts - reforigin) + reforigin;
              coil(icoil).xsplinepts(:,ipt) = xsplinepts;
          end
       end
   end  
end
end % if smartrot

toc
tic
disp(sprintf('\nRender Helices'))

figure(2)
clf

% loop over helices 
% draw each helix starting from xc,nc
% save output for strand locations and orientations at ends of helix
% helix vector is into helix at xc1 and out of helix at xc2

for ihelix = 1:nhelix
    npair   = helix(ihelix).npair;
    xc1     = helix(ihelix).xc(:,1);  % start of helix axis
    nc1     = helix(ihelix).nc(:,1);  % direction of helix axis
    thc1    = helix(ihelix).thc(:,1); % rotation around helix axis in degrees
    dzc1    = 0;    % translation along helix axis from start of helix axis
    if colorscheme == 1
       helixshades =  [helix(ihelix).strand(1); helix(ihelix).strand(2)];
    elseif colorscheme == 2
       helixshades = [2 3];
    end
    
    
    shades = [helixshades(1) helixshades(2) 17];  % colors of chain a, chain b, base pairs
    render = [1; 1; 1]; % render chain a, chain b, base pair struts
    
    acap = [helix(ihelix).cap(1) helix(ihelix).cap(2)];  % cap 5', 3' end of a chain
    bcap = [helix(ihelix).cap(4) helix(ihelix).cap(3)];  % cap 5', 3' end of b chain
    
    adye = [0 0];  % dye 5', 3' end of a chain
    bdye = [0 0];  % dye 5', 3' end of b chain
    [xc2,nc2,thc2,x1a,n1a,x2a,n2a,x1b,n1b,x2b,n2b] = drawhelix(npair,xc1,nc1,thc1,dzc1,shades,render,acap,bcap,adye,bdye);
    helix(ihelix).xc(:,2) = xc2;
    helix(ihelix).nc(:,2) = nc2;
    helix(ihelix).thc(2)  = thc2;
    helix(ihelix).x = [x1a, x2a, x1b, x2b];
    helix(ihelix).n = [n1a, n2a, n1b, n2b];
end


% define end points and vectors for coils
% vectors are into the coil

for icoil = 1:ncoil
    if coil(icoil).nbase > 1
        
    for icoilend = 1:2
        ihelix    = coil(icoil).tohelix(icoilend,1);
        ihelixend = coil(icoil).tohelix(icoilend,2);
        if ihelix > 0
          x = helix(ihelix).x(:,ihelixend);
          n = helix(ihelix).n(:,ihelixend); 
          bctype = 1; % bc based on first derivatives
        else % coil ends at a nick 
          nsplinepts = coil(icoil).nsplinepts;
          if ihelixend == -5
             if randomendcoils
                nchoose = floor(1/2*nsplinepts) + 1;
                x = coil(icoil).xsplinepts(:,nchoose);
             else
                x = coil(icoil).xsplinepts(:,1);
             end
             n = [0; 0; 0]; % set second derivative to zero
             bctype = 2; % bc based on second derivative
          elseif ihelixend == -3
             if randomendcoils
                nchoose = floor(1/2*nsplinepts) + 1;
                x = coil(icoil).xsplinepts(:,nchoose); 
             else
                x = coil(icoil).xsplinepts(:,coil(icoil).nsplinepts);  
             end
             n = [0; 0; 0]; % set second derivative to zero
             bctype = 2;    % bc based on second derivative
          end
        end
        % helix strand normals at strand ends point out of the strands
        coil(icoil).x(:,icoilend)       = x;
        coil(icoil).n(:,icoilend)       = n; % normals pointing into coil
        coil(icoil).bctype(icoilend)    = bctype;
    end
    coil(icoil).th(:,1) = 0; % starting rotation for bases around coil (I think)
    
    
    x1          = coil(icoil).x(:,1);
    x2          = coil(icoil).x(:,2);
    if coil(icoil).nsplinepts > 0 
        xsplinepts = coil(icoil).xsplinepts;
        if coil(icoil).tohelix(1,2) == -5
            if randomendcoils
                xp      = [x1(1) x2(1)];        % points to spline between start and finish
                yp      = [x1(2) x2(2)];
                zp      = [x1(3) x2(3)];  
            else
                xp      = [xsplinepts(1,:) x2(1)];        % points to spline between start and finish
                yp      = [xsplinepts(2,:) x2(2)];
                zp      = [xsplinepts(3,:) x2(3)]; 
            end  
        elseif coil(icoil).tohelix(2,2) == -3
            if randomendcoils
                xp      = [x1(1) x2(1)];        % points to spline between start and finish
                yp      = [x1(2) x2(2)];
                zp      = [x1(3) x2(3)];
            else
                xp      = [x1(1) xsplinepts(1,:)];        % points to spline between start and finish
                yp      = [x1(2) xsplinepts(2,:)];
                zp      = [x1(3) xsplinepts(3,:)];
            end
        else
            xp      = [x1(1) xsplinepts(1,:) x2(1)];        % points to spline between start and finish
            yp      = [x1(2) xsplinepts(2,:) x2(2)];
            zp      = [x1(3) xsplinepts(3,:) x2(3)];               
        end
    else
        xp      = [x1(1) x2(1)];        % points to spline between start and finish
        yp      = [x1(2) x2(2)];
        zp      = [x1(3) x2(3)];    
    end  
    coil(icoil).xsplinepts = [xp; yp; zp];
    end
end 

toc
tic
disp(sprintf('\nRender Coils'))
   
%
%   draw coils
%
for icoil = 1:ncoil
    if coil(icoil).nbase > 1        
    nbases  = coil(icoil).nbase;    % need to add one only to length of coils not at end of strand
    x1      = coil(icoil).x(:,1);   % start of coil 
    n1      = coil(icoil).n(:,1);   % starting coil vector (into coil)
    th1     = coil(icoil).th(:,1);  % starting rotation around n1
    x2      = coil(icoil).x(:,2);   % end of coil
    n2      = coil(icoil).n(:,2);   % ending coiled vector (into coil)
    xp      = coil(icoil).xsplinepts(1,:);
    yp      = coil(icoil).xsplinepts(2,:);
    zp      = coil(icoil).xsplinepts(3,:);
    bctype  = coil(icoil).bctype;   % 1: first derivative bc, 2: 2nd derivative bc
    if colorscheme == 1
        coilshade = coil(icoil).strand;
    elseif colorscheme == 2
        coilshade = 5; 
    end
    shades  = [coilshade*ones(nbases,1); 17];      % shades for coil and bases  
   
    render  = [1; coil(icoil).nick(1); 1; coil(icoil).nick(2)];  % render coil, first base, middle bases, last base
    ccap    = [coil(icoil).nick(1) coil(icoil).nick(2)];           % render cap at 5', 3' ends of coil
    
    cdye    = [0 0];           % render cap at 5', 3' ends of coil
    direction = 1;          % 1: 5' end at x1, 2: 5' end at x2
    
    if randomendcoils && (coil(icoil).nick(1) | coil(icoil).nick(2))
        helical = 1;        % make coil helical around spline
    else
        helical = 0;
    end
    
    %helical = 0;            
    random = 1;             % randomize angle increments in helix and also base pair orientations
                            % if helical=0, still randomizes base pair rotations around chain
    seed = 10;              % random seed (change to get different random chain)
    % arcratio              % adjust spline points to get arcratio \approx 1 for correct chain length
    coil(icoil).arcratio = drawcoil(nbases,x1,n1,th1,x2,n2,xp,yp,zp,bctype,helical,random,shades,render,ccap,cdye,direction,seed);
    end
end


%
% rendering options
%

%set(gcf,'Position',[400 100 1050 750])
colormap(cmap);
axis off
%axis([-10 310 -20 22 -10 230])
axis equal
%axis image
%view(-53,8)
view(2)
shading interp
lightangle(0,25)
lightangle(180,0)
set(gcf,'Renderer','zbuffer')
set(findobj(gca,'type','surface'),...
    'FaceLighting','phong',...
    'AmbientStrength',.5,'DiffuseStrength',.8,...
    'SpecularStrength',1.1,'SpecularExponent',25,...
    'BackFaceLighting','lit')

toc

% loop over coils
% take end information from helix strands bcs at either end
%

