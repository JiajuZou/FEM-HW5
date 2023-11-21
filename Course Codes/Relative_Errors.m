% clear the memory and the screen
clear all; clc;

% exact solution
exact = @(x) sin(x);
exact_square = @(x) sin(x).^2;
exact_dx = @(x) cos(x);
exact_dx_square = @(x) cos(x).^2;
% integral the exact solution and its derivitive
result_L2_down = sqrt(integral(exact_square, 0, 1)); %分母部分
result_H1_down = sqrt(integral(exact_dx_square, 0, 1)); %分母部分

% creat two tables to store the results
resultTable_L2 = table();
resultTable_H1 = table();

for n_el = 2:2:16 % generate different mesh
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

    % generate the quadrature rule to caculate the integral of one element
    n_int = 30;
    [xi, weight] = Gauss(n_int, -1, 1);
    
    Error_L2 = 0.0;
    for ee = 1 : n_el
        
        x_ele = zeros(2,1); % one element with two points
        for aa = 1 : 2
            x_ele(aa) = x_coor(IEN(aa,ee)); % A = IEN(a,e)
        end
       
        AA = 0.0; 
        for l = 1 : n_int
            u_h = 0.0;
            x_l = 0.0;
            dx_dxi = 0.0;
            for aa = 1 : 2
                u_h = u_h + exact(x_ele(aa)) * PolyShape(aa, xi(aa), 0);
                x_l = x_l + x_ele(aa) * PolyShape(aa, xi(l), 0);
                dx_dxi = dx_dxi + x_ele(aa) * PolyShape(aa, xi(l), 1);
            end
            u_exact = exact(x_l);
            AA = AA + weight(l) * (u_h - u_exact)^2 * dx_dxi;
        end
        
        Error_L2 = Error_L2 + AA;
    end
    
    Error_Final_L2 = Error_L2^0.5 / result_L2_down;
    %Store the results into the table and plot them
    resultTable_L2 = [resultTable_L2; table(hh, Error_Final_L2)];
    plot(log10(resultTable_L2.hh),log10(resultTable_L2.Error_Final_L2),'o-');
    hold on
    xlabel('lg(hh)');
    ylabel('lg(Error L2)');
    title('Plot of Error L2 vs. Mesh Size');

end








