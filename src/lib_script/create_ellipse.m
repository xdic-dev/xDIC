
function [xyFit,center] = create_ellipse(CxFit,CyFit,Rx,Ry,theta_radians)% Negate the angle because the function gives you the negative of the actual angle for some reason.
    thetaDegrees = -rad2deg(theta_radians); % Negate and convert to degrees.
    % First make the ellipse that is NOT rotated.
    numPoints = 100;
    angles = linspace(0, 360, numPoints);
    xFit = Rx * cosd(angles);
    yFit = Ry * sind(angles);
    % Compute the rotation matrix.  Ref: https://en.wikipedia.org/wiki/Rotation_matrix
    rotationMatrix = [cosd(thetaDegrees), -sind(thetaDegrees); sind(thetaDegrees), cosd(thetaDegrees)];
    xy = [xFit(:), yFit(:)]; % Load x and y into N-by-2 matrix.
    % Do the rotation by doing a matrix multiply of the data by the rotation matrix.
    xyRotated = xy * rotationMatrix;
    % Extract back out the rotated x and y.
    xFit = xyRotated(:, 1);
    yFit = xyRotated(:, 2);
    % Now translate/offset to put the center at where it should be.
    xyFit = [xFit + CxFit,yFit + CyFit]; 
    center = [CxFit,CyFit]; 
end