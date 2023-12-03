% clear the memory and the screen
clear all; clc;

%Problem defination
f = @(x) sin(x);
g = sin(1);
h = -cos(0);

% exact solution
exact = @(x) sin(x);
exact_dx = @(x) cos(x);

%solutions with different quadrature points
n_quad = 6;
solutions = cell(n_quad, 1);


    n_el = 2; % generate mesh
    hh = 1 / n_el / 3;
    x_coor = 0 : hh : 1;

    % IEN
    IEN = zeros(4, n_el);
    for ee = 1 : n_el
        IEN(1,ee) = 3*ee - 2;
        IEN(2,ee) = 3*ee - 1;
        IEN(3,ee) = 3*ee;
        IEN(4,ee) = 3*ee + 1;
    end

    % ID
    n_pt = 3*n_el + 1; % number of points
    ID = 1 : n_pt;
    ID(end) = 0;

    n_eq = n_pt - 1; % number of equations

    % generate the quadrature rule to caculate the integral of one element
    
    for n_int = 1:1:n_quad

        [xi, weight] = Gauss(n_int, -1, 1);
    
        % Assembly of K and F and solve Kd=F
        K = zeros(n_eq,n_eq);
        F = zeros(n_eq,1);
    
        for ee = 1 : n_el

            k_e = zeros(4,4); 
            f_e = zeros(4,1);

            x_ele = zeros(4,1);
            for aa = 1 : 4
                x_ele(aa) = x_coor(IEN(aa,ee)); % A = IEN(a,e)
            end

        for l = 1 : n_int
            dx_dxi = 0.0;
            x_l = 0.0;
            for aa = 1 : 4
                dx_dxi = dx_dxi + x_ele(aa) * PolyShape(aa, xi(l), 1);
                x_l = x_l + x_ele(aa) * PolyShape(aa, xi(l), 0);
            end
                dxi_dx = 1.0 / dx_dxi;

            for aa = 1 : 4
                for bb = 1 : 4
                    k_e(aa,bb) = k_e(aa,bb) + weight(l) * PolyShape(aa, xi(l), 1) * PolyShape(bb, xi(l), 1) * dxi_dx;
                end
            end

            for aa = 1 : 4
                f_e(aa) = f_e(aa) + weight(l) * PolyShape(aa, xi(l), 0) * f(x_l) * dx_dxi;
            end

    end

    % Now we need to put element k and f into global K and F
    for aa = 1 : 4
        AA = IEN(aa,ee);
        PP = ID(AA);
        if PP > 0
            F(PP) = F(PP) + f_e(aa);
            for bb = 1 : 4
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
    
   %将不同数目quadrature points获得的结果进行存储
    solutions{n_int} = d;
    
    end






