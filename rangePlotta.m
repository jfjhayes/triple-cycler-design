% Plot range from Sun for Cycler and Planets over time
figure;
hold on; grid on;

% LaTeX axes labels
xlabel('Calendar Time [days since start]', 'Interpreter', 'latex');
ylabel('Range from Sun [AU]', 'Interpreter', 'latex');

% Compute range for planets
time_steps = linspace(flyby_dates(1), flyby_dates(end), 1000); % Time steps for range computation
planet_ranges = struct();

% Planets: Venus, Earth, Mars
for planet_id = [2, 3, 4]
    planet_ranges.(sprintf('Planet%d', planet_id)) = zeros(1, length(time_steps));
    for t = 1:length(time_steps)
        [~, r, ~, ~] = planet_elements_and_sv_jd(planet_id, time_steps(t));
        planet_ranges.(sprintf('Planet%d', planet_id))(t) = norm(r * km_to_au); % Compute range in AU
    end
end

% Plot ranges for Venus, Earth, and Mars
% plot((time_steps - flyby_dates(1)), planet_ranges.Planet2, 'g-', 'DisplayName', 'Venus');
% plot((time_steps - flyby_dates(1)), planet_ranges.Planet3, 'b-', 'DisplayName', 'Earth');
% plot((time_steps - flyby_dates(1)), planet_ranges.Planet4, 'r-', 'DisplayName', 'Mars');

cycler_range = []; % To store range from Sun
cycler_time = [];  % To store time in days
for i = 1:length(flyby_dates)-1
    r1 = start_positions(i, :);
    [~, r2, ~, ~] = planet_elements_and_sv_jd(flyby_planets(i+1), flyby_dates(i+1));
    tof = (flyby_dates(i+1) - flyby_dates(i)) * 86400;

    if i == length(flyby_dates) - 1
        num_revolutions = 0; % Final leg with multiple revolutions
        [~, ~, ~, error, v1, ~, ~, ~] = lambertMR(r1, r2, tof-300000, mu_sun, 0, num_revolutions, 0, 0);
    else
        [~, ~, ~, error, v1, ~, ~, ~] = lambertMR(r1, r2, tof, mu_sun, 0, 0, 0, 0);
    end

    if error ~= 0
        warning('Lambert solver failed for segment %d-%d.', i, i+1);
        continue;
    end

    [t, rv] = ode45(@(t, rv) two_body_ode(t, rv, mu_sun), ...
                    linspace(0, tof, 1000), [r1'; v1']);
    trajectory_au = rv(:, 1:3) * km_to_au; % Convert to AU
    ranges = vecnorm(trajectory_au, 2, 2); % Compute range from Sun
    times = linspace(flyby_dates(i), flyby_dates(i+1), length(ranges)); % Time for this segment

    % Concatenate results
    cycler_range = [cycler_range; ranges]; % Append ranges
    cycler_time = [cycler_time, times];    % Append times (note the use of comma here)
end

yticks([1, 1.5, 2, 2.5, 3, 3.5, 4, 5, 6, 7, 8, 9, 10]); % Example: Set linear ticks for narrow ranges
yticklabels({'1', '1.5', '2', '2.5', '3', '3.5', '4', '5', '6', '7', '8', '9', '10'}); % Labels in linear format
set(gca, 'YScale', 'log');

% Convert Julian dates to datetime format
cycler_datetime = datetime(cycler_time, 'convertfrom', 'juliandate');
planet_datetime = datetime(time_steps, 'convertfrom', 'juliandate');

% Plot cycler range
plot(cycler_datetime, cycler_range, 'k-', 'DisplayName', 'Cycler Trajectory');

% Plot planetary ranges
plot(planet_datetime, planet_ranges.Planet2, 'g-', 'DisplayName', 'Venus');
plot(planet_datetime, planet_ranges.Planet3, 'b-', 'DisplayName', 'Earth');
plot(planet_datetime, planet_ranges.Planet4, 'r-', 'DisplayName', 'Mars');

% Adjust axis limits
%xlim([min(cycler_datetime), max(cycler_datetime)]);
%ylim([1e-2, max([cycler_range; planet_ranges.Planet2'; planet_ranges.Planet3'; planet_ranges.Planet4'])]);

% Add legend with LaTeX interpreter
legend('Location', 'northeast', 'Interpreter', 'latex');
legend show;

% Add labels and grid
xlabel('Calendar Date', 'Interpreter', 'latex');
ylabel('Range from Sun [AU] (log scale)', 'Interpreter', 'latex');
grid on;

% Save plot as EPS
saveas(gcf, 'outbound_range.eps', 'epsc');

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