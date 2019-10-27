function [r sideangle ang_ds f] = newtonloop(nside,sidesecant,ds,nds)

global nh nb rhc rbc rhelix rdh rcap ncap thcap dzb bturn brise dthb dsb dz dth ...
       cmap rdye ddye ndye majorgroove minorgroove dthgroove strutangle strutlength ...
       strutrise proptwist inclination extend5prime extend3prime extendcoil rcaphead;

%         r = nds*ds;
%         for iside = 1:nside
%             r = r + sidesecant(iside);
%         end
%         r = r/(2*pi);  % first approx is for large loops
%         r = max(r,1.05*rdh); % second approx is for stacked bases    
%                    
%         temp_ds   = ds/2 * (r^2 - (ds/2)^2)^-.5;
%         ang_ds    = 2*atan(temp_ds);
%         f         = nds*ang_ds - 2*pi;
%         for iside = 1:nside
%             sidetemp(iside)     = sidesecant(iside)/2 * (r^2 - (sidesecant(iside)/2)^2)^-.5;
%             sideangle(iside)    = 2*atan(sidetemp(iside));
%             f                   = f + sideangle(iside);
%         end
% 
%         for iter=1:5 % use Newton's method to find radius on which bases fall for this loop (default: 3 steps)
%             dtemp_ds    = -1/2 * ds/2 * (r^2 - (ds/2)^2)^-1.5 * 2*r;
%             dfdr        = nds*2/(1 + temp_ds^2) * dtemp_ds;
%             for iside = 1:nside
%                 dtemp    = -1/2 * sidesecant(iside)/2 * (r^2 - (sidesecant(iside)/2)^2)^-1.5 * 2*r;
%                 dfdr     = dfdr + 2/(1 + sidetemp(iside)^2) * dtemp;
%             end
% 
%             r       = r - f/dfdr;
%             temp_ds = ds/2 * (r^2 - (ds/2)^2)^-.5;
%             ang_ds  = 2*atan(temp_ds);
%             f       = nds*ang_ds - 2*pi;
%             for iside = 1:nside
%                 sidetemp(iside)     = sidesecant(iside)/2 * (r^2 - (sidesecant(iside)/2)^2)^-.5;
%                 sideangle(iside)    = 2*atan(sidetemp(iside));
%                 f                   = f + sideangle(iside);
%             end
%         end
%    
% 
% original = [r sideangle(1:nside) ang_ds]

    
       r = max(sidesecant(1:nside)/2); % max of stem and nick gaps
       r = max(r,ds/2) + .001;   % max of ssDNA gap
       f = nds*asin(ds/(2*r)) + sum(asin(sidesecant(1:nside)/(2*r))) - pi;
        while f > 1e-12 % use Newton's method to find radius on which bases fall for this loop (default: 3 steps)
            dfdr        = -1/(2*r^2)*nds*ds/sqrt(1 - (ds/(2*r))^2) ...
                          -1/(2*r^2)*sum(sidesecant(1:nside)./sqrt(1 - (sidesecant(1:nside)/(2*r)).^2));

            r       = r - f/dfdr;
            f = nds*asin(ds/(2*r)) + sum(asin(sidesecant(1:nside)/(2*r))) - pi;
        end
        sideangle(1:nside) = 2*asin(sidesecant(1:nside)/(2*r));
        ang_ds = 2*asin(ds/(2*r));

%new = [r sideangle(1:nside) ang_ds]
    