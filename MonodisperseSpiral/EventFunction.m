% Define the combined event function
function [value, isterminal, direction] = EventFunction(t, y,phim)
    value(1) = y(1);%y(1) is phi, the particle volume fraction and y(2) is sigma
    isterminal(1) = 1; % Halt integration when the event is triggered
    direction(1) = -1; % Detect events only when y is decreasing
    value(2) = y(1)-phim;
    isterminal(2) = 1; % Halt integration when the event is triggered
    direction(2) = 1; 
    % Detect events only when y is increasing and larger than phim, which
    % is the max particle volume fraction
end
