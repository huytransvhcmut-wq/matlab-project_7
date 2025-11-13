
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
    'Isobar(Constant P)', 'Isochore (Constant V)', 'Isotherm (Select new V)', 'Adiabatic', 'QUIT');

% If the next path leg is an isobar (Constant P)
if (PathType == 1)
    close(gcf);  % Close the figure window
    % Determine the new volume, pressure and temperature
    V(iPoint) = input('Enter new volume: ');
    P(iPoint) = P(iPoint-1);
    T(iPoint) = P(iPoint)*V(iPoint)/(nMoles*R);
    % Compute the work on done an isobar
    Work = P(iPoint)*( V(iPoint) - V(iPoint-1));
    fprintf('Work done on this process is %.2f J\n', Work);
    % Plot data (P,V)
    plot(V,P);
    xlabel('Volume (m^3)');
    ylabel('Pressure (Pa)');
    title('P-V diagram');

% If the next path legs is an isochore(Constant V)
elseif(PathType == 2)
    close(gcf); % Close the figure window
    % Determine the new volume, pressure and temperature
    V(iPoint) = V(iPoint-1);
    P(iPoint) = input('Enter new pressure: ');
    T(iPoint) = P(iPoint)*V(iPoint)/(nMoles*R);
    % Compute the work on done an isochore
    Work = 0;
    fprintf('Work done on this process is 0 J\n');
    % Plot data (P,V)
    plot(V,P);
    xlabel('Volume (m^3)');
    ylabel('Pressure (Pa)');
    title('P-V diagram');
 
    
 % If the next path legs is an isothermal(Constant T)   
elseif(PathType == 3)
    close(gcf); % Close the figure window
    % Determine the new volume, pressure and temperature
    T(iPoint) = T(iPoint-1);
    V(iPoint) = input('Enter new volume: ');
    P(iPoint) = (nMoles*R*T(iPoint))/V(iPoint);
    % Compute the work on done an isothermal
    Work = nMoles*R*T(iPoint)*log(V(iPoint)/V(iPoint-1));
    fprintf('Work done on this process is %.2f J\n',Work);
    % Plot data (P,V)
    x = linspace(V(iPoint-1), V(iPoint), 10000);
    y = (nMoles*R*T(iPoint))./x;
    plot(x,y);
    xlabel('Volume (m^3)');
    ylabel('Pressure (Pa)');
    title('P-V diagram');

    
 % If the next path legs is an adiabatic (Q = 0)
elseif(PathType == 4)
     close(gcf); % Close the figure window
     % Input the adiabatic index
     a = input('Enter the Poisson constant: ');
     % Determine the new volume, pressure and temperature
     V(iPoint) = input('Enter new volume: ');
     P(iPoint) = P(iPoint-1)*((V(iPoint-1)/V(iPoint))^a);
     T(iPoint) = (P(iPoint)*V(iPoint))/(nMoles*R);
     % Compute the work on done an adiabatic
     Work = (P(iPoint)*V(iPoint) - P(iPoint-1)*V(iPoint-1))/(1-a);
     fprintf('Work done on this process is %.2f J\n', Work);
     % Plot data (P,V)
     x = linspace(V(iPoint-1), V(iPoint), 10000);
     y = P(iPoint-1)*((V(iPoint-1)./x).^a);
     plot(x,y);
     xlabel('Volume (m^3)');
     ylabel('Pressure (Pa)');
     title('P-V diagram');
   
% Handle quitting the program
elseif (PathType == 5)
    disp('Program terminated by user.');
    return; % Ends the script

end

