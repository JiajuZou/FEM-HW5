% clear the memory and the screen
clear; clc;

% exact solution
exact = @(x) sin(x);

% problem definition
f = @(x) sin(x);
g = sin(1);
h = -cos(0);

% generate my mesh
n_el = 100; %the number of elements
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
n_int = 2;
[xi, weight] = Gauss(n_int, -1, 1);

K = zeros(n_eq, n_eq); % allocate the global stiffness matrix
F = zeros(n_eq, 1);    % allocate the global load vector

% Assembly of K and F
for ee = 1 : n_el

    k_e = zeros(2,2); %2*2 k matrix
    f_e = zeros(2,1);

    x_ele = zeros(2,1);
    for aa = 1 : 2
        x_ele(aa) = x_coor(IEN(aa,ee)); % A = IEN(a,e)
    end

    for l = 1 : n_int
        dx_dxi = 0.0;
        x_l = 0.0;
        for aa = 1 : 2
            dx_dxi = dx_dxi + x_ele(aa) * PolyShape(aa, xi(l), 1);
            x_l = x_l + x_ele(aa) * PolyShape(aa, xi(l), 0);
        end
        dxi_dx = 1.0 / dx_dxi;

        for aa = 1 : 2
            for bb = 1 : 2
                k_e(aa,bb) = k_e(aa,bb) + weight(l) * PolyShape(aa, xi(l), 1) * PolyShape(bb, xi(l), 1) * dxi_dx;
            end
        end

        for aa = 1 : 2
            f_e(aa) = f_e(aa) + weight(l) * PolyShape(aa, xi(l), 0) * f(x_l) * dx_dxi;
        end

    end

    % Now we need to put element k and f into global K and F
    for aa = 1 : 2
        AA = IEN(aa,ee);
        PP = ID(AA);
        if PP > 0
            F(PP) = F(PP) + f_e(aa);
            for bb = 1 : 2
                BB = IEN(bb,ee);
                QQ = ID(BB);
                if QQ > 0
                    K(PP,QQ) = K(PP,QQ) + k_e(aa,bb);
                else
                    F(PP) = F(PP) - k_e(aa,bb) * g;
                end
            end
        end
    end

    if ee == 1
        F(ID(IEN(1,ee))) = F(ID(IEN(1,ee))) + h;
    end
end

% Now we have K and F
% Solve Kd = F
uh = K \ F;

d = [uh; g];


%% Bonus
%Define the parameters of gmres
restart = 10000;
maxit = 10000;

% Set the value of tol
tol_values = [1e-2, 1e-4, 1e-6];

% Create array to store the solutions
solutions_gmres = cell(length(tol_values), 1);
Comparasion = cell(length(tol_values), 1);

% Solve the linear system using gmres with different tols
for i = 1:length(tol_values)
    tol = tol_values(i);
    x = gmres(K, F, restart, tol, maxit);
    solutions_gmres{i} = x;
    Comparasion{i}= x-uh;
end

Differrence = Comparasion{2}-Comparasion{1};
Differrence2 = Comparasion{3}-Comparasion{2};
%Both of them are zero








