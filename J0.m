function j0 = J0(year, month, day)
% -------------------------------------------------------------------------
% Function: J0
% Purpose:
% This function computes the Julian Day Number (J0) at 0 UT for a given 
% Gregorian calendar date, from Curtis. It is valid for any year 
% between 1900 and 2100.
%
% Inputs:
%   year  - Gregorian year (valid range: 1901–2099)
%   month - Month of the year (1–12)
%   day   - Day of the month (1–31)
%
% Outputs:
%   j0 - Julian Day Number at 0 UT (Universal Time)
%
% Dependencies:
%   None.
% -------------------------------------------------------------------------

% Compute the Julian Day Number using the standard algorithm
j0 = 367 * year - fix(7 * (year + fix((month + 9) / 12)) / 4) ...
     + fix(275 * month / 9) + day + 1721013.5;
end

