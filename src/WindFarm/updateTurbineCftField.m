function meso = updateTurbineCftField(meso, farm)

nx = meso.nx;
ny = meso.ny;

for i=1:nx
    for j=1:ny
        if(meso.sx(i,j) ~= -1 && meso.sy(i,j) ~= -1)
            meso.ctf(i,j) = (pi*farm.turbines(meso.tclose(i,j)).Ct*farm.turbines(meso.tclose(i,j)).D^2) / (4*meso.sx(i,j)*meso.sy(i,j));
        else
            meso.ctf(i,j) = 0.0;
        end
    end
end

end