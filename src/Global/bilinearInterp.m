function value = bilinearInterp(xs,ys,dx,dy,nx,ny,f,xquery,yquery)

    % IMPORTANT: assumes uniformly spaced mesh to speed up the algorithm

    iClose = ceil((xquery-xs)/dx);
    jClose = ceil((yquery-ys)/dy);
    
    il     = iClose;
    ir     = iClose + 1;
    jb     = jClose;
    jt     = jClose + 1;
    
    % would go out of right boundary
    if(il>=nx)
        ir     = nx;
        il     = nx;
    end
    
    % would go out of left boundary
    if(ir<=1)
        il     = 1;
        ir     = 1;
    end
    
    % would go out of top boundary
    if(jb>=ny)
        jt     = ny;
        jb     = ny;
    end
    
    % would go out of bottom boundary
    if(jt<=1)
        jb     = 1;
        jt     = 1;
    end
    
    if(il~=ir)
        xr         = xs + (ir-1)*dx;
        xl         = xs + (il-1)*dx;
        wl = (xr - xquery)/dx;
        wr = (xquery - xl)/dx;
    else
        wl = 0.5;
        wr = 0.5;
    end
    
    if(jt~=jb)
        yt         = ys + (jt-1)*dy;
        yb         = ys + (jb-1)*dy;
        wb = (yt - yquery)/dy;
        wt = (yquery - yb)/dy;
    else
        wt = 0.5;
        wb = 0.5;
    end

    valuet = wl*f(il,jt) + wr*f(ir,jt);
    valueb = wl*f(il,jb) + wr*f(ir,jb);
    value  = wt*valuet + wb*valueb;

end