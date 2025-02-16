function avgD = averageTurbineDiameter(farm)

avgD = 0;
for t=1:farm.Nt
    avgD = avgD + farm.turbines(t).D;
end
avgD = avgD / farm.Nt;

end