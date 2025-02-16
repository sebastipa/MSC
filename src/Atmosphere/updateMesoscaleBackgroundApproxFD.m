function meso = updateMesoscaleBackgroundApproxFD(meso, sol, abl)
%%
fprintf('3LMR calculating background velocity...\n');

tstart = tic;
fprintf(' > Assembled matrix in ');

% get some info on the system
spparms('spumoni',2)

% Performs incomplete LU factorization of the linear system matrix
nx = meso.nx;
ny = meso.ny;
dx = meso.dx;
dy = meso.dy;

% layer velocities
U1 = abl.U1;
V1 = abl.V1;
U2 = abl.U2;
V2 = abl.V2;
DU = U2 - U1;
DV = V2 - V1;
fc  = abl.fc;
rho = abl.rho;
nut1 = abl.nut1;

H1       = abl.H1;
c_coeff  = abl.tauMag1 / (U1^2 + V1^2) / sqrt(U1^2 + V1^2);
d_coeff = abl.tauMag2 / (DU^2 + DV^2) / sqrt(DU^2 + DV^2);
C11      = ( 2*U1^2 + V1^2 ) * c_coeff;
C12      = ( U1 * V1       ) * c_coeff;
C21      = ( V1 * U1       ) * c_coeff;
C22      = ( 2*V1^2 + U1^2 ) * c_coeff;

D11      = ( 2*DU^2 + DV^2 ) * d_coeff;

mtid_i   = zeros(12*nx*ny,1);
mtid_j   = zeros(12*nx*ny,1);
values   = zeros(12*nx*ny,1);
b        = zeros(2*nx*ny,1);
id       = 0;

for q=1:nx
    for j=1:ny

        [L,R,B,T]   = getStencil(q,j,nx,ny);

        uID         = getID(q,j,ny);
        uID_Left    = getID(L,j,ny);
        uID_Right   = getID(R,j,ny);
        uID_Bottom  = getID(q,B,ny);
        uID_Top     = getID(q,T,ny);

        vID         = getID(q,j,ny) + nx*ny;
        vID_Left    = getID(L,j,ny) + nx*ny;
        vID_Right   = getID(R,j,ny) + nx*ny;
        vID_Bottom  = getID(q,B,ny) + nx*ny;
        vID_Top     = getID(q,T,ny) + nx*ny;

        % matrix columns to set (10)
        mtid_j(id+1 ) = uID;
        mtid_j(id+2 ) = uID_Left;
        mtid_j(id+3 ) = uID_Right;
        mtid_j(id+4 ) = uID_Bottom;
        mtid_j(id+5 ) = uID_Top;

        mtid_j(id+6 ) = vID;
        mtid_j(id+7 ) = vID_Left;
        mtid_j(id+8 ) = vID_Right;
        mtid_j(id+9 ) = vID_Bottom;
        mtid_j(id+10) = vID_Top;

        % matri rows to set (2)
        mtid_i(id+1:id+5)  = uID;
        mtid_i(id+6:id+10) = vID;
        
        % i j
        values(id+1) = values(id+1) + (2*nut1/dx^2 + 2*nut1/dy^2 + C11/H1);
        % i-1 j
        values(id+2) = values(id+2) + (-0.5*U1/dx - nut1/dx^2);
        % i+1 j
        values(id+3) = values(id+3) + ( 0.5*U1/dx - nut1/dx^2);
        % i j-1
        values(id+4) = values(id+4) + (-0.5*V1/dy - nut1/dy^2);
        % i j+1
        values(id+5) = values(id+5) + ( 0.5*V1/dy - nut1/dy^2);
        
        % i j
        values(id+6) = values(id+6) + (2*nut1/dx^2 + 2*nut1/dy^2 + C22/H1);
        % i-1 j
        values(id+7) = values(id+7) + (-0.5*U1/dx - nut1/dx^2);
        % i+1 j
        values(id+8) = values(id+8) + ( 0.5*U1/dx - nut1/dx^2);
        % i j-1
        values(id+9) = values(id+9) + (-0.5*V1/dy - nut1/dy^2);
        % i j+1
        values(id+10)= values(id+10)+ ( 0.5*V1/dy - nut1/dy^2);

        % coupling terms 
        values(id+11) = values(id+11) + ( - fc + C12/H1);
        mtid_i(id+11) = uID;
        mtid_j(id+11) = vID;

        values(id+12) = values(id+12) + ( fc + C21/H1);
        mtid_i(id+12) = vID;
        mtid_j(id+12) = uID;

        % compute pressure gradient 
        dpdx = (meso.p(R,j) - meso.p(L,j))/(2*rho*dx);
        dpdy = (meso.p(q,T) - meso.p(q,B))/(2*rho*dy);

        b(uID) = -dpdx;
        b(vID) = -dpdy;

        id = id + 12;

        % set stress field to placeholder 
        meso.placeHldr1(q,j) = D11*(meso.u2(q,j) - meso.u1(q,j));
        meso.placeHldr2(q,j) = -C11* meso.u1(q,j);
    end
end

Amat        = sparse(mtid_i,mtid_j,values);

fprintf('%.2f s\n',toc(tstart));

tstart = tic;
fprintf(' > Solved in ');

s = Amat\b;

for q=1:nx
    for j=1:ny
        
        uID = getID(q,j,ny);
        vID = getID(q,j,ny) + nx*ny;

        meso.ur1(q,j) = s(uID);
        meso.vr1(q,j) = s(vID);
    end
end

fprintf('%.2f s\n',toc(tstart));

% update background hub velocity with pressure effect
[ubk, vbk] = mesoscaleAvg2Hub(meso.ur1, meso.vr1, abl);
meso.ubk   = ubk; 
meso.vbk   = vbk; 

if(sol.deepArray && sol.localShearStress)
    
    fprintf('3LMR calculating background velocity with shear stress...\n');
    tstart = tic;
    fprintf(' > Assembled matrix in ');

    mtid_i   = zeros(12*nx*ny,1);
    mtid_j   = zeros(12*nx*ny,1);
    values   = zeros(12*nx*ny,1);
    b        = zeros(2*nx*ny,1);
    id       = 0;
    
    for q=1:nx
        for j=1:ny
    
            [L,R,B,T]   = getStencil(q,j,nx,ny);
    
            uID         = getID(q,j,ny);
            uID_Left    = getID(L,j,ny);
            uID_Right   = getID(R,j,ny);
            uID_Bottom  = getID(q,B,ny);
            uID_Top     = getID(q,T,ny);
    
            vID         = getID(q,j,ny) + nx*ny;
            vID_Left    = getID(L,j,ny) + nx*ny;
            vID_Right   = getID(R,j,ny) + nx*ny;
            vID_Bottom  = getID(q,B,ny) + nx*ny;
            vID_Top     = getID(q,T,ny) + nx*ny;
    
            % matrix columns to set (10)
            mtid_j(id+1 ) = uID;
            mtid_j(id+2 ) = uID_Left;
            mtid_j(id+3 ) = uID_Right;
            mtid_j(id+4 ) = uID_Bottom;
            mtid_j(id+5 ) = uID_Top;
    
            mtid_j(id+6 ) = vID;
            mtid_j(id+7 ) = vID_Left;
            mtid_j(id+8 ) = vID_Right;
            mtid_j(id+9 ) = vID_Bottom;
            mtid_j(id+10) = vID_Top;
    
            % matri rows to set (2)
            mtid_i(id+1:id+5)  = uID;
            mtid_i(id+6:id+10) = vID;

            if(abs(meso.taux1(q,j)) > 1e-10 )
                c11 = 0.0;
                c22 = 0.0;
                c12 = 0.0;
                c21 = 0.0;
            else
                c11 = C11/H1;
                c22 = C22/H1;
                c12 = C12/H1;
                c21 = C21/H1;
            end
            
            % i j
            values(id+1) = values(id+1) + (2*nut1/dx^2 + 2*nut1/dy^2 + c11);
            % i-1 j
            values(id+2) = values(id+2) + (-0.5*U1/dx - nut1/dx^2);
            % i+1 j
            values(id+3) = values(id+3) + ( 0.5*U1/dx - nut1/dx^2);
            % i j-1
            values(id+4) = values(id+4) + (-0.5*V1/dy - nut1/dy^2);
            % i j+1
            values(id+5) = values(id+5) + ( 0.5*V1/dy - nut1/dy^2);
            
            % i j
            values(id+6) = values(id+6) + (2*nut1/dx^2 + 2*nut1/dy^2 + c22);
            % i-1 j
            values(id+7) = values(id+7) + (-0.5*U1/dx - nut1/dx^2);
            % i+1 j
            values(id+8) = values(id+8) + ( 0.5*U1/dx - nut1/dx^2);
            % i j-1
            values(id+9) = values(id+9) + (-0.5*V1/dy - nut1/dy^2);
            % i j+1
            values(id+10)= values(id+10)+ ( 0.5*V1/dy - nut1/dy^2);

            % coupling terms 
            values(id+11) = values(id+11) + ( - fc + c12);
            mtid_i(id+11) = uID;
            mtid_j(id+11) = vID;
    
            values(id+12) = values(id+12) + ( fc + c21);
            mtid_i(id+12) = vID;
            mtid_j(id+12) = uID;
    
            % compute pressure gradient 
            dpdx = (meso.p(R,j) - meso.p(L,j))/(2*rho*dx);
            dpdy = (meso.p(q,T) - meso.p(q,B))/(2*rho*dy);
    
            b(uID) = -dpdx + ( - meso.taux2(q,j) + meso.taux1(q,j)) / H1;
            b(vID) = -dpdy + ( - meso.tauy2(q,j) + meso.tauy1(q,j)) / H1;
    
            id = id + 12;
        end
    end
    
    Amat        = sparse(mtid_i,mtid_j,values);
    
    fprintf('%.2f s\n',toc(tstart));
    
    tstart = tic;
    fprintf(' > Solved in ');
    
    s = Amat\b;
    
    for q=1:nx
        for j=1:ny
            
            uID = getID(q,j,ny);
            vID = getID(q,j,ny) + nx*ny;
    
            meso.ur1_da(q,j) = s(uID);
            meso.vr1_da(q,j) = s(vID);
        end
    end
    
    fprintf('%.2f s\n',toc(tstart));

    % update background hub velocity with pressure and deep array effect
    [ubk_da, vbk_da] = mesoscaleAvg2Hub(meso.ur1_da, meso.vr1_da, abl);
    meso.ubk_da      = ubk_da; 
    meso.vbk_da      = vbk_da; 

end 

end

function [L,R,B,T] = getStencil(i,j,nx,ny)
    
    L = i - 1;
    R = i + 1;
    B = j - 1;
    T = j + 1;

    if(i==1); L = 1; end
    if(j==1); B = 1; end
    if(i==nx); R = nx; end
    if(j==ny); T = ny; end
end

function [id] = getID(i,j,ny)
    
    id = j + (i-1)*ny;
end
