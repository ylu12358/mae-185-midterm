function theta = machAngle(field, fieldname, xx, yy, dx, dy, M)
    [nx, ~] = size(field);

    % Establish threshold difference for which to estimate shock line
    % Tresholds distinct for first and last points as flow physics are more
    % dynamic near leading edge
    threshold = 0.01;
    end_threshold = 0.05;
    pt1 = 10;

    % Compute the percentage difference for first to last row
    diff_first = diff(field(pt1, 2:end), 1, 2) ./ field(pt1, 2:end-1);
    diff_end = diff(field(end, 2:end), 1, 2) ./ field(end, 2:end-1);
    
    % Find the first element that meets the threshold
    col_first = find(abs(diff_first) > threshold, 1, 'last');
    col_end = find(abs(diff_end) > end_threshold, 1, 'last');

    % Estimate slopes between points, get angle
    theta = atan((col_end-col_first)*dy/(dx*(nx-10)))*180/pi;

    % Verify slope of line found by plotting field
    pcolor(xx,yy,field);
    hold on;
    % Plot line using computed data
    plot([xx(pt1,col_first) xx(nx,col_end)], ...
        [yy(pt1,col_first) yy(nx,col_end)], ...
        'LineWidth',5,'Color',[1 0.1 0.1]);
    shading interp; 
    axis equal tight;
    % Set labels and title
    xlabel('x'); 
    ylabel('y');
    title(['Mach angle = ', num2str(theta), 'º verification using ', ...
        fieldname, ' for Ma = ', num2str(M), ' (arcsin(1/Ma) = ', ...
        num2str(asin(1/M)*180/pi), 'º)'])
end