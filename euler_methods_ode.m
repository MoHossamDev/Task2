function euler_methods_ode()
    % Manual input of values
    t0 = input('Enter the initial time t0: ');
    tf = input('Enter the final time tf: ');
    y0 = input('Enter the initial value y0: ');
    h = input('Enter the step size h: ');

    % Create the time vector
    t = t0:h:tf;

    % Display options to the user
    fprintf('Choose the numerical method you want to use:\n');
    fprintf('1 - Forward Euler\n');
    fprintf('2 - Modified Euler\n');
    fprintf('3 - Backward Euler\n');
    fprintf('4 - Runge-Kutta 2nd order\n');
    fprintf('5 - Runge-Kutta 4th order\n');
    fprintf('6 - Adams-Bashforth 2-step\n');

    choice = input('Enter the desired method number: ');

    % Call the selected method and display the result
    switch choice
        case 1
            y = forward_euler(t, y0, h);
            method_name = 'Forward Euler';
        case 2
            y = modified_euler(t, y0, h);
            method_name = 'Modified Euler';
        case 3
            y = backward_euler(t, y0, h);
            method_name = 'Backward Euler';
        case 4
            y = runge_kutta2(t, y0, h);
            method_name = 'Runge-Kutta 2nd Order';
        case 5
            y = runge_kutta4(t, y0, h);
            method_name = 'Runge-Kutta 4th Order';
        case 6
            y_FE = forward_euler(t, y0, h); % Need initial value for AB2
            y = adams_bashforth2(t, y0, h, y_FE);
            method_name = 'Adams-Bashforth 2-step';
        otherwise
            error('Invalid choice.');
    end

    % Plot the results
    figure;
    plot(t, y, 'b-', 'LineWidth', 1.5);
    xlabel('Time t'); ylabel('y(t)');
    title(['Numerical Solution using ', method_name]);
    grid on;
end

function y_FE = forward_euler(t, y0, h)
    n = length(t);
    y_FE = zeros(1, n);
    y_FE(1) = y0;
    for i = 1:n-1
        y_FE(i+1) = y_FE(i) + h * (-50*y_FE(i) + sin(t(i)));
    end
end

function y_ME = modified_euler(t, y0, h)
    n = length(t);
    y_ME = zeros(1, n);
    y_ME(1) = y0;
    for i = 1:n-1
        y_predictor = y_ME(i) + h * (-50*y_ME(i) + sin(t(i)));
        y_ME(i+1) = y_ME(i) + (h/2) * ((-50*y_ME(i) + sin(t(i))) + (-50*y_predictor + sin(t(i+1))));
    end
end

function y_BE = backward_euler(t, y0, h)
    n = length(t);
    y_BE = zeros(1, n);
    y_BE(1) = y0;
    for i = 1:n-1
        y_BE(i+1) = (y_BE(i) + h * sin(t(i+1))) / (1 + 50*h);
    end
end

function y_RK2 = runge_kutta2(t, y0, h)
    n = length(t);
    y_RK2 = zeros(1, n);
    y_RK2(1) = y0;
    for i = 1:n-1
        k1 = -50*y_RK2(i) + sin(t(i));
        k2 = -50*(y_RK2(i) + h*k1/2) + sin(t(i) + h/2);
        y_RK2(i+1) = y_RK2(i) + h*k2;
    end
end

function y_RK4 = runge_kutta4(t, y0, h)
    n = length(t);
    y_RK4 = zeros(1, n);
    y_RK4(1) = y0;
    for i = 1:n-1
        k1 = -50*y_RK4(i) + sin(t(i));
        k2 = -50*(y_RK4(i) + h*k1/2) + sin(t(i) + h/2);
        k3 = -50*(y_RK4(i) + h*k2/2) + sin(t(i) + h/2);
        k4 = -50*(y_RK4(i) + h*k3) + sin(t(i) + h);
        y_RK4(i+1) = y_RK4(i) + (h/6)*(k1 + 2*k2 + 2*k3 + k4);
    end
end

function y_AB2 = adams_bashforth2(t, y0, h, y_FE)
    n = length(t);
    y_AB2 = zeros(1, n);
    y_AB2(1) = y0;
    y_AB2(2) = y_FE(2); % Use Forward Euler for the first step
    for i = 2:n-1
        y_AB2(i+1) = y_AB2(i) + (h/2) * (3*(-50*y_AB2(i) + sin(t(i))) - (-50*y_AB2(i-1) + sin(t(i-1))));
    end
end
