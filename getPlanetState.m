function [r, v, jd] = getPlanetState(date, planet_id)
    [~, r, v, jd] = planet_elements_and_sv(planet_id, ...
        date(1), date(2), date(3), date(4), date(5), date(6));
end