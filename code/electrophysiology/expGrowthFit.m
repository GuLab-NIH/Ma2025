% Fitting 1-phase exponential decay/growth function
%
% by June Hoan Kim, December 30, 2025

function fitResult = expGrowthFit(data)
    x = data(:,1);
    y = data(:,2);
    
    %% Define exponential growth model: y = A + B*exp(-k*x)
    expModel = fittype('A + B*exp(-k*x)', ...
        'independent', 'x', ...
        'coefficients', {'A','B','k'});
    
    %% Set initial conditions to improve fitting stability
    A0 = mean(y(end-3:end));   % approximate asymptote using the last few values
    B0 = y(1) - A0;            % offset from starting value to asymptote
    k0 = 1e-5;                 % small initial k because x values are large
    
    startPoints = [A0, B0, k0];
    
    %% Perform curve fitting
    try
        [fitParams, gof] = fit(x, y, expModel, 'Start', startPoints);
    catch
        warning('Fit failed at index %d, skipping...', i);
        fitResult.Params.k = NaN;
        fitResult.gof = NaN;
        return;
    end
    fitResult = struct();
    fitResult.Params = fitParams;
    fitResult.gof = gof;
    %% Display results
    disp(fitParams);
    disp(gof);
    
    %% Plot raw data and fitted curve
    % figure; hold on;
    % h = plot(fitParams, x, y);  % let cfit/plot handle data & fit together
    % set(h(1), 'Marker', 'o', 'MarkerFaceColor', 'k'); % data
    % set(h(2), 'Color', 'r', 'LineWidth', 2);          % fitted curve
    % 
    % xlabel('x');
    % ylabel('y');
    % title('Exponential Growth Fit:  y = A + B e^{-kx}');
    % legend('Data', 'Exponential Fit');
    % grid on;
end