
clearvars, clc
% Input the initial pressure, volume, temperature, number of gas moles
iPoint = 1;
P(iPoint) = input('Enter initial pressure P0(Pa): ');
V(iPoint) = input('Enter initial volume V0(m^3): ');
nMoles    = input('Enter number of gas moles(moles): ');
R = 8.314; % Ideal gas constant
T(iPoint) = P(iPoint)*V(iPoint)/(nMoles*R);

% Select type of path(isobar, isochore or isotherm) or quit
iPoint = iPoint + 1;
fprintf('For leg #%g \n', iPoint-1);
PathType = menu(sprintf('Leg %g: Select next path',iPoint-1), ...
    'Isobar(Constant P)', 'Isochore (Constant V)', 'Isotherm (Select new V)', 'QUIT');

% If the next path leg is an isobar (Constant P)
if (PathType == 1)
    close(gcf);  % Close the figure window
    % Determine the new volume, pressure and temperature
    V(iPoint) = input('Enter new volume: ');
    P(iPoint) = P(iPoint-1);
    T(iPoint) = P(iPoint)*V(iPoint)/(nMoles*R);
    % Compute the work on done an isobar
    Work = P(iPoint)*( V(iPoint) - V(iPoint-1))
    % Plot data
    plot(V,P);

% If the next path legs is an isochore(Constant V)
elseif(PathType == 2)
    close(gcf); % Close the figure window
    % Determine the new volume, pressure and temperature
    V(iPoint) = V(iPoint-1);
    P(iPoint) = input('Enter new pressure: ');
    T(iPoint) = P(iPoint)*V(iPoint)/(nMoles*R);
    % Compute the work on done an isochore
    Work = 0
    % Plot data
    plot(V,P);
 
    
 % If the next path legs is an isothermal(Constant T)   
elseif(PathType == 3)
    close(gcf); % Close the figure window
    % Determine the new volume, pressure and temperature
    T(iPoint) = T(iPoint-1);
    V(iPoint) = input('Enter new volume: ');
    P(iPoint) = (nMoles*R*T(iPoint))/V(iPoint);
    % Compute the work on done an isothermal
    Work = nMoles*R*T(iPoint)*log(V(iPoint)/V(iPoint-1))
    % Plot data
    x = linspace(V(iPoint-1), V(iPoint), 10000);
    y = (nMoles*R*T(iPoint))./x;
    plot(x,y);

end

