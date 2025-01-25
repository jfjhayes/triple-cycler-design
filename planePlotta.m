% Debugging Planetary Orbits, Flyby Points, and Unified Trajectory in AU (XY and XZ planes)
clear; clc; close all;

% Inputs - change based on fmincon output 
% INBOUND BELOW
% start_date = 2.460863657910894e+06; % Julian date of trajectory start
% relative_dates = [526.0622, 381.3189, 2.0040e+03, 373.5517];
% end_date = start_date + 365.25 * 12.8; % Julian date of trajectory end

% OUTBOUND BELOW
start_date = 2.461841570758452e+06; % Julian date of trajectory start
relative_dates = [274.7570, 339.7034, 477.2836, 157.5392];
end_date = start_date + 365.25 * 12.8; % Julian date of trajectory end

mu_sun = 1.32712440018e11; % Gravitational parameter of the Sun

% Conversion factor: km to AU
km_to_au = 1 / 149597870.7; % 1 AU = 149,597,870.7 km

% Compute flyby dates
flyby_dates = [start_date];
for i = 1:length(relative_dates)
    flyby_dates = [flyby_dates, flyby_dates(end) + relative_dates(i)];
end
flyby_dates = [flyby_dates, end_date];

% Flyby planets - change based on fmincon
% INBOUND BELOW
%flyby_planets = [4, 3, 2, 3, 4, 4]; % Planet IDs

% OUTBOUND BELOW
flyby_planets = [3, 3, 4, 3, 2, 3]; % Planet IDs

% Define line styles for planets using containers.Map
planet_styles = containers.Map({'2', '3', '4'}, {'g-', 'b-', 'r-'}); % Venus, Earth, Mars

% Define custom colors and markers for each flyby
flyby_colors = {'black', "#0072BD", "#D95319", "#EDB120", "#A2142F", 'black'}; % Example colors
flyby_markers = {'diamond', 'o', 'o', 'o', 'o', 'square'}; % Example markers

%% XY Plane Plot
figure;
hold on; grid on; axis equal;
xlabel('$X$ [AU]', 'Interpreter', 'latex'); 
ylabel('$Y$ [AU]', 'Interpreter', 'latex');
%title('Planetary Orbits, Flyby Points, and Trajectory (XY Plane)', 'Interpreter', 'latex');

% Crop the axes to ±1.75 AU
%xlim([-1.75, 1.75]);
%ylim([-1.75, 1.75]);

% Plot the Sun
plot(0, 0, 'yo', 'MarkerSize', 12, 'MarkerFaceColor', 'yellow', 'DisplayName', 'Sun');

% Add a single legend entry for the trajectory
plot(NaN, NaN, 'k-', 'DisplayName', 'Trajectory'); % Invisible line for legend entry

% Plot each planet's orbit in the XY plane
unique_planets = unique(flyby_planets);
for planet_id = unique_planets
    time_steps = linspace(flyby_dates(1), flyby_dates(end), 1000);
    planet_positions = zeros(3, length(time_steps));
    for t = 1:length(time_steps)
        [~, r, ~, ~] = planet_elements_and_sv_jd(planet_id, time_steps(t));
        planet_positions(:, t) = r * km_to_au; % Convert to AU
    end
    style = planet_styles(num2str(planet_id)); % Retrieve style for the planet
    plot(planet_positions(1, :), planet_positions(2, :), style, ...
         'DisplayName', sprintf('%s', planet_name(planet_id)));
end

% Plot flyby points and trajectory segments in the XY plane
start_positions = [];
for i = 1:length(flyby_dates)
    [~, r_flyby, ~, ~] = planet_elements_and_sv_jd(flyby_planets(i), flyby_dates(i));
    r_flyby_au = r_flyby * km_to_au;
    start_positions = [start_positions; r_flyby];

    color = flyby_colors{i};
    marker = flyby_markers{i};
    plot(r_flyby_au(1), r_flyby_au(2), marker, 'MarkerSize', 6, ...
         'MarkerFaceColor', color, 'Color', color, ...
         'DisplayName', sprintf('Flyby %d (%s)', i, planet_name(flyby_planets(i))));
end

for i = 1:length(flyby_dates)-1
    r1 = start_positions(i, :);
    [~, r2, ~, ~] = planet_elements_and_sv_jd(flyby_planets(i+1), flyby_dates(i+1));
    tof = (flyby_dates(i+1) - flyby_dates(i)) * 86400;

    num_revolutions = 0;
    if i== 3
        [~, ~, ~, error, v1, ~, ~, ~] = lambertMR(r1, r2, tof-580000, mu_sun, 0, num_revolutions, 0, 0);
    else
        [~, ~, ~, error, v1, ~, ~, ~] = lambertMR(r1, r2, tof, mu_sun, 0, 0, 0, 0);
    end

    if error ~= 0
        warning('Lambert solver failed for segment %d-%d.', i, i+1);
        continue;
    end

    [t, rv] = ode45(@(t, rv) two_body_ode(t, rv, mu_sun), ...
                    linspace(0, tof, 10000), [r1'; v1']);

    trajectory_au = rv(:, 1:3) * km_to_au;
    plot(trajectory_au(:, 1), trajectory_au(:, 2), 'k-', 'HandleVisibility', 'off');
end

legend('Location', 'northeastoutside', 'Interpreter', 'latex');
legend show;

% Save XY Plane plot as EPS
saveas(gcf, 'outbound_trajectory_xy_plane.eps', 'epsc');

%% XZ Plane Plot
figure;
hold on; grid on; axis equal;
xlabel('$X$ [AU]', 'Interpreter', 'latex'); 
ylabel('$Z$ [AU]', 'Interpreter', 'latex');
%title('Planetary Orbits, Flyby Points, and Trajectory (XZ Plane)', 'Interpreter', 'latex');

% Crop the axes to ±1.75 AU
%xlim([-1.75, 1.75]);
%ylim([-1.75, 1.75]);

% Plot the Sun
plot(0, 0, 'yo', 'MarkerSize', 12, 'MarkerFaceColor', 'yellow', 'DisplayName', 'Sun');

% Add a single legend entry for the trajectory
plot(NaN, NaN, 'k-', 'DisplayName', 'Trajectory'); % Invisible line for legend entry

% Plot each planet's orbit in the XZ plane
for planet_id = unique_planets
    time_steps = linspace(flyby_dates(1), flyby_dates(end), 1000);
    planet_positions = zeros(3, length(time_steps));
    for t = 1:length(time_steps)
        [~, r, ~, ~] = planet_elements_and_sv_jd(planet_id, time_steps(t));
        planet_positions(:, t) = r * km_to_au; % Convert to AU
    end
    style = planet_styles(num2str(planet_id)); % Retrieve style for the planet
    plot(planet_positions(1, :), planet_positions(3, :), style, ...
         'DisplayName', sprintf('%s', planet_name(planet_id)));
end

% Plot flyby points and trajectory segments in the XZ plane
for i = 1:length(flyby_dates)
    [~, r_flyby, ~, ~] = planet_elements_and_sv_jd(flyby_planets(i), flyby_dates(i));
    r_flyby_au = r_flyby * km_to_au;

    color = flyby_colors{i};
    marker = flyby_markers{i};
    plot(r_flyby_au(1), r_flyby_au(3), marker, 'MarkerSize', 6, ...
         'MarkerFaceColor', color, 'Color', color, ...
         'DisplayName', sprintf('Flyby %d (%s)', i, planet_name(flyby_planets(i))));
end

for i = 1:length(flyby_dates)-1
    r1 = start_positions(i, :);
    [~, r2, ~, ~] = planet_elements_and_sv_jd(flyby_planets(i+1), flyby_dates(i+1));
    tof = (flyby_dates(i+1) - flyby_dates(i)) * 86400;

    if i == length(flyby_dates) - 1
        num_revolutions = 0; % Final leg with multiple revolutions
        [~, ~, ~, error, v1, ~, ~, ~] = lambertMR(r1, r2, tof+300000, mu_sun, 0, num_revolutions, 0, 0);
    else
        [~, ~, ~, error, v1, ~, ~, ~] = lambertMR(r1, r2, tof, mu_sun, 0, 0, 0, 0);
    end

    if error ~= 0
        warning('Lambert solver failed for segment %d-%d.', i, i+1);
        continue;
    end

    if i == length(flyby_dates) - 1
        [t, rv] = ode45(@(t, rv) two_body_ode(t, rv, mu_sun), ...
                        linspace(0, tof, 10000), [r1'; v1']);
    else
        [t, rv] = ode45(@(t, rv) two_body_ode(t, rv, mu_sun), ...
                        linspace(0, tof, 10000), [r1'; v1']);
    end

    trajectory_au = rv(:, 1:3) * km_to_au;
    plot(trajectory_au(:, 1), trajectory_au(:, 3), 'k-', 'HandleVisibility', 'off');
end

legend('Location', 'northeastoutside', 'Interpreter', 'latex');
legend show;

% Save XZ Plane plot as EPS
saveas(gcf, 'outbound_trajectory_xz_plane.eps', 'epsc');

% Helper function for planet names
function name = planet_name(id)
    switch id
        case 2
            name = 'Venus';
        case 3
            name = 'Earth';
        case 4
            name = 'Mars';
        otherwise
            name = 'Unknown';
    end
end

% Two-body dynamics function
function drdt = two_body_ode(~, rv, mu)
    r = rv(1:3); % Position vector (km)
    v = rv(4:6); % Velocity vector (km/s)
    r_norm = norm(r);
    drdt = [v; -mu * r / r_norm^3]; % Velocity and acceleration
end
