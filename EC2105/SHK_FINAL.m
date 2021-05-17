
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

%V = 0;
%for m = 1 : 10
%    for n = 1 : 10
%        V = V + coeff*sin(m*pi*x/a)*sin(n*pi*y/b)*sinh(z*k3);
%    end
%end
% Voltage output of laplace equation

% dividerx = (coeff*(sin(m*pi*x/a)*sinh(z*k3)))^(-1);
% x = asind(V*dividerx)*(a/(m*pi));
% get the separable x equation

% dividery = (coeff*(sin(n*pi*y/b)*sinh(z*k3)))^(-1);
% y = asind(V*dividery)*(b/(n*pi));
% get the separable y equation

% dividerz = (coeff*(sin(m*pi*x/a)*sin(n*pi*y/b)))^(-1);
% z = asinh(V*dividerz)*(1/(k3));
% get the separable z equation

X = 1 : 10;
Y = 1 : 10;
Z = 1 : 10;
V = 0;
for m = 1 : 10
    for n = 1 : 10
        V = V + coeff*sin(m*pi*X(m)/a)*sin(n*pi*Y(n)/b).*sinh(Z.*k3);
    end
end

disp(size(V));
%[X,Y,Z] = meshgrid(0:1:a, 0:1:b, 0:1:c);
%[U,K,W] = surfnorm(V);

%quiver3(X, Y, Z, U, K, W);   % Plotting quiver(flux lines) in cartesian coordinates
%hold on;
%contour(X, Y, Z);   % Drawing equipotential lines
%axis equal;
%hold off;



