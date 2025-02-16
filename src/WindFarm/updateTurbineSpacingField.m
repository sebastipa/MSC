function meso = updateTurbineSpacingField(meso, farm)

fprintf('\nCalculating turbine spacing field...\n');

% search radius 
refDist = 1200;

% visibility cone half angle 
alpha   = 45;

% edge detection cone half angle 
gamma   = 80;

for i=meso.farmBuf.imin:meso.farmBuf.imax
    for j=meso.farmBuf.jmin:meso.farmBuf.jmax

        meso.tclose(i,j) = -1;
        minDist          = 1e20;
        avgSx            = 0;
        avgSy            = 0;
        avgNx            = 0;
        avgNy            = 0;

        cellOutside = 0;
        
        % find spacing relative to this wind turbine 
        for t = 1:farm.Nt

            % edge detection 
            hasTopTrbs = 0;
            hasBotTrbs = 0;
            hasLftTrbs = 0;
            hasRgtTrbs = 0;

            xRef = farm.turbines(t).rotor(1);
            yRef = farm.turbines(t).rotor(2);

            % turbine-cell distance 
            dist_t2c_x = meso.x(i) - xRef;
            dist_t2c_y = meso.y(j) - yRef;
            angle_t2c_x = abs(atan2d(dist_t2c_y,dist_t2c_x));
            angle_t2c_y = abs(atan2d(dist_t2c_x,dist_t2c_y));

            pointDist = sqrt(dist_t2c_x^2 + dist_t2c_y^2);

            if(pointDist > refDist)
                continue;
            else

                % closest wind turbine to this cell 
                if(pointDist < minDist)
                    minDist = pointDist;
                    meso.tclose(i,j) = t;
                end
    
                minDist_x = 1e20;
                minDist_y = 1e20;
    
                % find the spacing relative to the ref turbine 
                for tt = 1:farm.Nt
                    if(tt ~= t)
    
                        dist_t2t_x = farm.turbines(tt).rotor(1) - xRef;
                        dist_t2t_y = farm.turbines(tt).rotor(2) - yRef;
                        dist   = sqrt(dist_t2t_x^2 + dist_t2t_y^2);

                        if(dist > 2*refDist)
                            continue;
                        end
    
                        % turbines are considered if angles are between 
                        % smaller than abs(+-30) 
                        % greater than abs(+-150)
                        angle_t2t_x = abs(atan2d(dist_t2t_y,dist_t2t_x));
                        angle_t2t_y = abs(atan2d(dist_t2t_x,dist_t2t_y));
    
                        % downstream visibility cone 
                        if(angle_t2t_x < alpha) 
                            if(abs(dist_t2t_x) < minDist_x)
                                minDist_x = abs(dist_t2t_x);
                            end
                            hasRgtTrbs = 1;
                        end

                        % upstream visibility cone 
                        if(angle_t2t_x > 180-alpha)
                            if(abs(dist_t2t_x) < minDist_x)
                                minDist_x = abs(dist_t2t_x);
                            end
                            hasLftTrbs = 1;
                        end
    
                        % top visibility cone 
                        if(angle_t2t_y < alpha) 
                            if(abs(dist_t2t_y) < minDist_y)
                                minDist_y = abs(dist_t2t_y);
                            end
                            hasTopTrbs = 1;
                        end

                        % bottom visibility cone 
                        if(angle_t2t_y > 180-alpha)
                            if(abs(dist_t2t_y) < minDist_y)
                                minDist_y = abs(dist_t2t_y);
                            end
                            hasBotTrbs = 1;
                        end
                    end
                end

                % left edge  
                if(hasLftTrbs==0 && angle_t2c_x > 180 - gamma)
                    cellOutside = 1;
                end

                % right edge  
                if(hasRgtTrbs==0 && angle_t2c_x < gamma)
                    cellOutside = 1;
                end

                % top edge  
                if(hasTopTrbs==0 && angle_t2c_y < gamma)
                    cellOutside = 1;
                end

                % bottom edge  
                if(hasBotTrbs==0 && angle_t2c_y > 180 - gamma)
                    cellOutside = 1;
                end
           
                % average spacings from turbines in the radius 
                if(minDist_y < 1e20)
                    avgSy = avgSy + minDist_y;
                    avgNy = avgNy + 1;
                end
                if(minDist_x < 1e20)
                    avgSx = avgSx + minDist_x;
                    avgNx = avgNx + 1;
                end
            end    
        end

        % average spacing of at this cell
        if(avgNx > 0 && avgNy > 0 && cellOutside == 0)
            meso.sx(i,j) = avgSx / avgNx;
            meso.sy(i,j) = avgSy / avgNy;
        end
    end
end

end