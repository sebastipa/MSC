function [udeep, vdeep, uStar_hi, uStar_lo, deltaFarms, deltaWakes, z0_hi] = topDownModelPoint(xs, ys, x, y, ubak, vbak, limit, cft, uStarFree, z0, H, L, zh, D)

% xs,ys are the coordinates of the query point 
% x and y are the coordinates of the mesh assumed cartesian and uniform  
% ubak and vbak are the backgroun wind components to move along streamlines
% limit is the upstream x limit after which there are no more turbines 
% cft is the wind farm thrust coefficient field 
% uStar is the freestream friction velocity 
% z0 is the freestream roughness length 
% H is the freestream ABL height 
% zh is the average hub height 
% D is the average wind turbine diameter 
x0 = x(1);
y0 = y(1);
dx = x(2) - x(1);
dy = y(2) - y(1);
nx = length(x);
ny = length(y);

% compute local uStar 
uStar = bilinearInterp(x0,y0,dx,dy,nx,ny,uStarFree,xs,ys);

% compute local freestream 
uFree          = topDownModel(uStar, z0, H, L, 0, zh, D, -1, -1);

deficit        = 0;
vdeep          = 0;
uStar_hi       = uStar;
uStar_lo       = uStar;
deltaFarms     = zh + D/2;
deltaWakes     = zh + D/2;
z0_hi          = z0;

% find the number of wind farms upstream
dx             = x(2) - x(1);
dy             = y(2) - y(1);
ds             = min(dx,dy);

nfarms         = 0;
farm_cfts      = 0;
farm_uStarFree = 0;
inFarm_dist    = 0;
inWake_dist    = 0;

xq             = xs - 1.0; % move a bit downstream for the detection to work 
yq             = ys;

% it is of interest only the portion of the line upstream the point, so we
% loop until xp < limit. 

% initialize wake distance to 0
wake_dist = 0;

while(xq>limit)

    cftq = bilinearInterp(x0,y0,dx,dy,nx,ny,cft,xq,yq);

    % detected the end of a wind farm
    if(cftq>0)
        
        nfarms              = nfarms + 1;
        inFarm_dist(nfarms) = 0;         % will compute inside the next while loop 
        inWake_dist(nfarms) = wake_dist; % already computed: we are entering a wind farm from behind 
        navg                = 0;
        cftavg              = 0;
        
        % loop until the end of the wind farm or until stopid. The end of
        % the wind farm is when the body force goes back to zero
        while(cftq>0 && xq>limit)
            cftavg = cftavg + cftq;
            navg   = navg + 1;
            
            % go upstream 
            ux    = -1.0*bilinearInterp(x0,y0,dx,dy,nx,ny,ubak,xq,yq);
            uy    = -1.0*bilinearInterp(x0,y0,dx,dy,nx,ny,vbak,xq,yq);
            um    = sqrt(ux^2 + uy^2);

            xq    = xq + ds*ux/um;
            yq    = yq + ds*uy/um;
            
            inFarm_dist(nfarms) = inFarm_dist(nfarms) + ds;

            cftq  = bilinearInterp(x0,y0,dx,dy,nx,ny,cft,xq,yq);
        end

        % save freestream velocity for weighted superposition
        farm_uStarFree(nfarms) = bilinearInterp(x0,y0,dx,dy,nx,ny,uStarFree,xq,yq);

        if(navg~=0)
            farm_cfts(nfarms)    = cftavg / navg;
        else
            farm_cfts(nfarms)    = 0;
        end
    end

    % go upstream 
    ux    = -1.0*bilinearInterp(x0,y0,dx,dy,nx,ny,ubak,xq,yq);
    uy    = -1.0*bilinearInterp(x0,y0,dx,dy,nx,ny,vbak,xq,yq);
    um    = sqrt(ux^2 + uy^2);

    xq    = xq + ds*ux/um;
    yq    = yq + ds*uy/um;

    wake_dist = wake_dist + ds;
end

% if there are farms upstream this point evaluate
if(nfarms>0)

    uStar_hi_isolated_p  = zeros(nfarms,1);
    uStar_lo_isolated_p  = zeros(nfarms,1);
    deltaFarm_isolated_p = zeros(nfarms,1);
    deltaWake_isolated_p = zeros(nfarms,1);
    z0_hi_isolated_p     = zeros(nfarms,1);
        
    for f=1:nfarms
       
        % last farm upwind to the query point. Query point is the true point 
        wakeDist_f = inWake_dist(f);
        farmDist_f = inFarm_dist(f);
        cft_f      = farm_cfts(f);

        % run top-down model for this chunk of line and update characteristic stresses
        [u,                          ...
            uStar_hi_isolated_p(f),  ...
            uStar_lo_isolated_p(f),  ...
            deltaFarm_isolated_p(f), ...
            deltaWake_isolated_p(f), ...
            z0_hi_isolated_p(f),     ...
        ]...
        = topDownModel(uStar, z0, H, L, cft_f, zh, D, farmDist_f, wakeDist_f);

        uFree_up = farm_uStarFree(f)/0.4*(log(zh/z0) + 4.7*zh/L);
        deficit  = deficit + (uFree_up/uFree)*(uFree - u);
    end

    % do superposition 
    uStar_hi   = max(uStar_hi_isolated_p);
    uStar_lo   = min(uStar_lo_isolated_p);
    
    z0_hi      = max(z0_hi_isolated_p(f));
    deltaFarms = max(deltaFarm_isolated_p);
    deltaWakes = max(deltaWake_isolated_p);
end

ux         = bilinearInterp(x0,y0,dx,dy,nx,ny,ubak,xs,ys);
uy         = bilinearInterp(x0,y0,dx,dy,nx,ny,vbak,xs,ys);
um         = sqrt(ux^2 + uy^2);

udeep      = (uFree - deficit)*ux/um;
vdeep      = (uFree - deficit)*uy/um;

end