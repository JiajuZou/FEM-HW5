% clear the memory and the screen
clear all; clc;

% exact solution
exact = @(x) sin(x);

% problem definition
f = @(x) sin(x);
g = sin(1);
h = -cos(0);

for n_el = 2:2:16
% generate my mesh
% n_el = 5;
hh = 1 / n_el;
x_coor = 0 : hh : 1;

% IEN
IEN = zeros(2, n_el);

for ee = 1 : n_el
    IEN(1,ee) = ee;
    IEN(2,ee) = ee+1;
end

% ID
n_pt = n_el + 1; % number of points
ID = 1 : n_pt;
ID(end) = 0;

n_eq = n_pt - 1; % number of equations

% generate the quadrature rule
n_int = 30;
[xi, weight] = Gauss(n_int, -1, 1);

for ee = 1 : n_el

    x_ele = zeros(2,1);
    for aa = 1 : 2
        x_ele(aa) = x_coor(IEN(aa,ee)); % A = IEN(a,e)
    end

    for l = 1 : n_int
        dx_dxi = 0.0;
        x_l = 0.0;
        u_h = 0.0;
        u = 0.0;
        for aa = 1 : 2
            dx_dxi = dx_dxi + x_ele(aa) * PolyShape(aa, xi(l), 1);
            x_l = x_l + x_ele(aa) * PolyShape(aa, xi(l), 0);
         
            u_h = u_h + exact(x_ele(aa)) * PolyShape(aa, xi(aa), 0);
            u = exact(x_l);
            
        end
        dxi_dx = 1.0 / dx_dxi;
        
        Error = 0.0;
        Error = Error + weight(l) * (u_h - u)^2 * dx_dxi;
        
    end    
end

plot(log(hh),log(Error),'o-');
hold on
end





