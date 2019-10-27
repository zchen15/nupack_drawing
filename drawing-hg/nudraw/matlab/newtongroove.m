function [dthgroove,strutrise,strutlength,strutangle] = newtongroove(inclination,rdh,brise,minorgroove)

dthgroove = pi;
sa = sin(dthgroove - pi);
ca = cos(dthgroove - pi);
strutrise = rdh*tan(inclination)*sqrt(sa^2 + (1 + ca)^2); % rise of strut along its length
f = dthgroove/(2*pi) * brise - strutrise - minorgroove; 

        for iter=1:3 % use Newton's method to find radius on which bases fall for this loop (default: 3 steps)
            dstrutdth = rdh*tan(inclination)* .5*(sa^2 + (1+ca)^2)^-.5 ...
                *(2*sa*ca - 2*(1 + ca)*sa);
            dfdth = brise/(2*pi)  - dstrutdth;

            dthgroove       = dthgroove - f/dfdth;
            
            sa = sin(dthgroove - pi);
            ca = cos(dthgroove - pi);
            strutrise = rdh*tan(inclination)*sqrt(sa^2 + (1 + ca)^2);
            f = dthgroove/(2*pi) * brise - strutrise - minorgroove;
        end
        
        strutlength = rdh/cos(inclination) * sqrt(sa^2 + (1 + ca)^2); %length of strut
        strutangle = pi + atan2(rdh*sa,rdh*(1+ca)); % angle of strut relative to x axis 
                        % with origin at one end of strut and rotating
                        % around z axis