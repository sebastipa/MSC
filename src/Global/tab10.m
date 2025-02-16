function [color] = tab10(number)

% Matplotlib tab10 colormap 

colors = ...
    {...
        "#1F77B4",...
        "#FF7f0E",...
        "#2CA02C",...
        "#D62728",...
        "#9467BD",...
        "#8c564B",...
        "#E377C2",...
        "#7F7F7F",...
        "#BCBd22",...
        "#17BECF",...
    };

if(number>10)
    fprintf("Error: number exceeds number of defined colors in tab10 palette\n");
else
    color = colors{number};
end

end