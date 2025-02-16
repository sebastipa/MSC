function [ids, weights, iClose, jClose] = findBilinaerInterpParam(xquery, yquery, xs, ys, dx, dy, nx, ny)

% IMPORTANT: assumes uniformly spaced mesh to speed up the algorithm

iClose = ceil((xquery-xs)/dx);
jClose = ceil((yquery-ys)/dy);

ids.il     = iClose;
ids.ir     = iClose + 1;
ids.jb     = jClose;
ids.jt     = jClose + 1;

% would go out of right boundary
if(ids.il>=nx)
    ids.ir     = nx;
    ids.il     = nx;
    weights.wl = 0.5;
    weights.wr = 0.5;
end

% would go out of left boundary
if(ids.ir<=1)
    ids.il     = 1;
    ids.ir     = 1;
    weights.wl = 0.5;
    weights.wr = 0.5;
end

% would go out of top boundary
if(ids.jb>=ny)
    ids.jt     = ny;
    ids.jb     = ny;
    weights.wt = 0.5;
    weights.wb = 0.5;
end

% would go out of bottom boundary
if(ids.jt<=1)
    ids.jb     = 1;
    ids.jt     = 1;
    weights.wt = 0.5;
    weights.wb = 0.5;
end

if(ids.il~=ids.ir)
    xr         = xs + (ids.ir-1)*dx;
    xl         = xs + (ids.il-1)*dx;
    weights.wl = (xr - xquery)/dx;
    weights.wr = (xquery - xl)/dx;
else
    weights.wl = 0.5;
    weights.wr = 0.5;
end

if(ids.jt~=ids.jb)
    yt         = ys + (ids.jt-1)*dy;
    yb         = ys + (ids.jb-1)*dy;
    weights.wb = (yt - yquery)/dy;
    weights.wt = (yquery - yb)/dy;
else
    weights.wt = 0.5;
    weights.wb = 0.5;
end

end