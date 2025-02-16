function value = bilinearInterpolate(flt, frt, frb, flb, weights)

% bilinear interpolation formula where 
% flt = top-left value 
% frt = top-right value 
% frb = bottom-right value 
% flb = bottom-left value 

valuet = weights.wl*flt + weights.wr*frt;
valueb = weights.wl*flb + weights.wr*frb;
value  = weights.wt*valuet + weights.wb*valueb;

end