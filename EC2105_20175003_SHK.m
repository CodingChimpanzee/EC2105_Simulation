% 20175003 강성현
% MATLAB R2021a version code
% Rectangular box a, b, c
% side/bottom grounded
% top isolated, potential V_0
% Potential distribution? => equpotential lines, flux lines are necessary

a = 10;
b = 10;
c = 10;     % Total Range of x, y, z axis (i.e. a, b, c. I will let this as 100)
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

X = 2 : a-1;
Y = 2 : b-1;
Z = 2 : c-1;

% Using Laplace equation grad(grad(V)) = 0
% Solution of 3D Laplace equation in cartesian coordinates is
% (c1cos(k1x) + c2sin(k1x))*(c3cos(k2y) + c4sin(k2y))*(c5cosh(k3z) +
% c6sinh(k3z)), and by BC,
% c1, c3, c5 = 0, k1 = m*pi/a, k2 = n*pi/b, 
% k3 = sqrt(k1^2 + k2^2), Let coeff = c2*c4*c6
% hence sigma_sigma_(coeff)*sin(k1)*sin(k2)*sinh(k3), putting value V_0,
% Then, the value = coeff*sinh(k3) will be below

syms m n;
syms x y z;
constant = 4*TopV/(a*b);
k3 = pi*(sqrt((power(m, 2)/power(a, 2))+(power(n, 2)/power(b, 2))));
integ = int(sin(m*pi*x/a), x, 0, a);
value = constant*int(integ*sin(n*pi*y/b), y, 0, b);
coeff = value/sinh(k3*c);
% Get the coefficient

sum = V;

M = 10;        % Number of iterations
N = 10;        % Number of iterations
for p = 1 : M
    for q = 1 : N
        sumX = sin(p*pi*X/a);
        sumY = sin(q*pi*Y/b);
        sumZ = sinh(k3*Z);
    end
end
sum = sum*coeff;
% In this part, I followed Laplace solution from
% https://www.youtube.com/watch?v=QopKseQYnyw&ab_channel=Vidya-mitra

% Conclusion within z condition
divider = (coeff*(sumX.*sumY)).^(-1);
Z = asinh(TopV*divider)*(1/(k3));

quiver3(x, y, z, Z);   % Plotting quiver(flux lines) in cartesian coordinates
hold on;
contour();   % Drawing equipotential lines
axis equal;
hold off;
