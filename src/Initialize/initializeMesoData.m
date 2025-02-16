function meso = initializeMesoData(sol, abl)

fprintf('\nInitializing meso-scale data...\n');

fid = fopen('./input/settings.msm', 'r'); 

if (fid<=0)
    application = mfilename;
    error([application ':io'],'Unable to open settings.msm file');
end 

while ~feof(fid)
    tline        = fgetl(fid);
    keyword      = sscanf(tline,'%s',[1 1]);
    tline        = tline(length(keyword)+1:end);
    
    intflgValue  = sscanf(tline,'%d',[1 Inf]);
    scalarValue  = sscanf(tline,'%f',[1 Inf]);
    vectorValue  = sscanf(tline,'%f %f',[1 Inf]);
    word         = sscanf(tline,'%s',[1 1]);

    switch keyword
        case 'xsMeso',  xs                 = scalarValue * 1000;
        case 'ysMeso',  ys                 = scalarValue * 1000;
        case 'xeMeso',  xe                 = scalarValue * 1000;
        case 'yeMeso',  ye                 = scalarValue * 1000;
        case 'dxMeso',  dx                 = scalarValue;
        case 'dyMeso',  dy                 = scalarValue;
    end
end

fclose(fid);

% check input 
if(xs>=xe)
    application = mfilename;
    error([application ':io'],'mesoscale domain xStartMeso cannot be bigger than xEndMeso');
end

if(ys>=ye)
    application = mfilename;
    error([application ':io'],'mesoscale domain yStartMeso cannot be bigger than yEndMeso');
end

% create mesoscale domain
% meso.Lx = xe-xs;
% meso.Ly = ye-ys;
% nxdes = ceil(meso.Lx / dx) + 1;
% nydes = ceil(meso.Ly / dx) + 1;
% nxexp = ceil(log2(nxdes));
% nyexp = ceil(log2(nydes));
% meso.nx = 2^nxexp; 
% meso.ny = 2^nyexp; 
% meso.dx = meso.Lx / (meso.nx-1);
% meso.dy = meso.Ly / (meso.ny-1);
% meso.xs = xs;
% meso.xe = xe;
% meso.ys = ys;
% meso.ye = ye;

meso.nx = ceil((xe-xs) / dx) + 1;
meso.ny = ceil((ye-ys) / dy) + 1;

% make nx and ny even
if(mod(meso.nx,2)==1), meso.nx = meso.nx+1; end
if(mod(meso.ny,2)==1), meso.ny = meso.ny+1; end

meso.dx = dx;
meso.dy = dy;
meso.xs = xs;
meso.xe = xs + (meso.nx-1)*meso.dx;
meso.ys = ys;
meso.ye = ys + (meso.ny-1)*meso.dy;
meso.Lx = xe-xs;
meso.Ly = ye-ys;

% physical space 
meso.x       = zeros(meso.nx,1);
meso.y       = zeros(meso.ny,1);
meso.xx      = zeros(meso.nx,meso.ny);
meso.yy      = zeros(meso.nx,meso.ny);

% fourier space 
meso.dk      = 2*pi / (meso.nx*meso.dx);
meso.dl      = 2*pi / (meso.ny*meso.dy);
meso.k       = meso.dk*(-meso.nx/2:meso.nx/2-1) + meso.dk/10; % offset to avoid singular matrix 
meso.l       = meso.dl*(-meso.ny/2:meso.ny/2-1) + meso.dl/20; % offset to avoid singular matrix
%meso.k       = meso.dk*(-meso.nx/2:meso.nx/2-1); 
%meso.l       = meso.dl*(-meso.ny/2:meso.ny/2-1); 

% create arrays 
for i=1:meso.nx
    meso.x(i) = xs + (i-1)*meso.dx;
end

for j=1:meso.ny
    meso.y(j) = ys + (j-1)*meso.dy;    
end

for i=1:meso.nx
    for j=1:meso.ny
        meso.xx(i,j) = xs + (i-1)*meso.dx;
        meso.yy(i,j) = ys + (j-1)*meso.dy; 
    end
end

% solution of the 1st 3LM run 
meso.u1        = zeros(meso.nx,meso.ny);
meso.u2        = zeros(meso.nx,meso.ny);
meso.v1        = zeros(meso.nx,meso.ny);
meso.v2        = zeros(meso.nx,meso.ny);
meso.p         = zeros(meso.nx,meso.ny);
meso.eta       = zeros(meso.nx,meso.ny);

% solution of the 2nd 3LM run (background velocity due to pressure only)
meso.ur1  = zeros(meso.nx,meso.ny);
meso.vr1  = zeros(meso.nx,meso.ny);

% hub background velocity on the mesoscale domain (pressure effects)
[ubk, vbk]  = mesoscaleAvg2Hub(meso.ur1, meso.vr1, abl);
meso.ubk    = ubk; 
meso.vbk    = vbk;

% hub background velocity on the mesoscale domain (pressure and deep array effects)
meso.ubk_da = ubk; 
meso.vbk_da = vbk; 

% body force 
meso.fx     = zeros(meso.nx,meso.ny);
meso.fy     = zeros(meso.nx,meso.ny);
meso.f      = zeros(meso.nx,meso.ny);
meso.bfMask = zeros(meso.nx,meso.ny);

% turbine spacing 
meso.sx     = -ones(meso.nx,meso.ny);
meso.sy     = -ones(meso.nx,meso.ny);
meso.ctf    = zeros(meso.nx,meso.ny);
meso.tclose = zeros(meso.nx,meso.ny);

% dummy fields for debugging 
meso.placeHldr1 = zeros(meso.nx,meso.ny);
meso.placeHldr2 = zeros(meso.nx,meso.ny);

% allocate wind farm buffer 
meso.farmBuf.imin = 1;
meso.farmBuf.imax = 1;
meso.farmBuf.jmin = 1;
meso.farmBuf.jmax = 1;

end