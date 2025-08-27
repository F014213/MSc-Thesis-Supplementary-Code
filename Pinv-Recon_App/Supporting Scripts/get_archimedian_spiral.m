%sofia pearson
%21/05/2025
%get_archimedian_spiral

%growth_factor = 0.001; % for the Archimedian spiral (how close together
%arms are)

%k_max = 1; %max radius of the spiral trajectory
%if rad_trajectory =1, then the kspace will range from -1 to 1

%d_theta = 0.1; % angle step size


function [Kx, Ky] = get_archimedian_spiral(growth_factor, k_max, d_theta)
    
    %initialise variables
    Kx = 0;
    Ky = 0;
    theta = 0;

    while Kx(end)<k_max && Kx(end)>((-1)*k_max)
        theta = theta + d_theta;
        r = growth_factor * theta;
        x = r .* cos(theta);
        y = r .* sin(theta);
        Kx = [Kx, x];
        Ky = [Ky, y];
    end
  
    % max(Kx)
    % max(Ky)
    % min(Kx)
    % min(Ky)
    
    % Plot
    % figure;
    % plot(Kx, Ky), axis equal, title('Archimedean Spiral');
end