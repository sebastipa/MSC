function meso  = setWindFarmBufferOnDomain(meso, farm)

fprintf('\nCalculating wind farm buffers...\n');

nx = meso.nx;
ny = meso.ny;

xmin = max(min(farm.x), min(meso.x));
ymin = max(min(farm.y), min(meso.y));
xmax = min(max(farm.x), max(meso.x));
ymax = min(max(farm.y), max(meso.y));

imin = 0;
imax = 0;
jmin = 0;
jmax = 0;

dmin_min = 1e20;
dmin_max = 1e20;

for i=1:nx

    dmin = abs(xmin - meso.x(i));
    dmax = abs(xmax - meso.x(i));

    if(dmin < dmin_min)
        dmin_min = dmin;
        imin     = i;
    end

    if(dmax < dmin_max)
        dmin_max = dmax;
        imax     = i;
    end
end

dmin_min = 1e20;
dmin_max = 1e20;

for j=1:ny

    dmin = abs(ymin - meso.y(j));
    dmax = abs(ymax - meso.y(j));

    if(dmin < dmin_min)
        dmin_min = dmin;
        jmin     = j;
    end

    if(dmax < dmin_max)
        dmin_max = dmax;
        jmax     = j;
    end
end

meso.farmBuf.imin = imin;
meso.farmBuf.imax = imax;
meso.farmBuf.jmin = jmin;
meso.farmBuf.jmax = jmax;

end