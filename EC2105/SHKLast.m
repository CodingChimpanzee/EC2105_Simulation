a = 10; b = 10; c = 10;
syms x y z;
syms m n;
TopV = 1000;

constant = 4*TopV/(a*b);
k3 = pi*(sqrt((power(m, 2)/power(a, 2))+(power(n, 2)/power(b, 2))));
integ = int(sin(m*pi*x/a), x, 0, a);
value = constant*int(integ*sin(n*pi*y/b), y, 0, b);
coeff = value/sinh(k3*c);
% Get the k3, coefficient

sumX = 0;
for m = 1 : 100
    sumX = sumX + sin(m*pi*x/a);
end

sumY = 0;
for n = 1 : 100
    sumY = sumY + sin(n*pi*y/b);
end

sumZ = 0;
for m = 1 : 100
    for n = 1 : 100
        sumZ = sumZ + sin(z*pi*(sqrt((power(m, 2)/power(a, 2))+(power(n, 2)/power(b, 2)))));
    end
    disp("done");
end

L = sumX*sumY*sumZ;     % Laplace equation solved

V = zeros(a, b, c);      % 3D Voltage potential Matrix 

TopV = 1000;    % I will let Top voltage potential V_0 = 1000
Bot1V = 0;      % Grounded
Sid1V = 0;      % Grounded
Sid2V = 0;      % Grounded
Sid3V = 0;      % Grounded
Sid4V = 0;      % Grounded

% Boundary Conditions
V(:,:,c) = TopV;
V(:,:,1) = Bot1V;
V(:,1,:) = Sid1V;
V(1,:,:) = Sid2V;
V(:,b,:) = Sid3V;
V(a,:,:) = Sid4V;
% Set all initial voltage level of the box


%xAxis = linspace(1,10,10);
%yAxis = linspace(1,10,10);
%zAxis = linspace(1,10,10);
%[X, Y, Z] = meshgrid(xAxis,yAxis,zAxis);
%[U, K, W];
%quiver3(X, Y, Z, U, K, W);   % Plotting quiver(flux lines) in cartesian coordinates
%hold on;
%contour(X, Y, Z);   % Drawing equipotential lines
%axis equal;
%hold off;

