% clear the memory and the screen
clear all; clc;

% exact solution
exact = @(x) sin(x);
% 创建一个空的表格
resultTable = table();

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

    for ee = 1 : n_el
        x_ele = zeros(2,1); % one element with two points
        
        for aa = 1 : 2
            x_ele(aa) = x_coor(IEN(aa,ee)); % A = IEN(a,e)
        end
        
        for l = 1 : n_int
            dx_dxi = 0.0;
            x_l = 0.0;
            u_h = 0.0;
            u = 0.0;
            Error = 0.0;
            for aa = 1 : 2
                dx_dxi = dx_dxi + x_ele(aa) * PolyShape(aa, xi(l), 1);
                x_l = x_l + x_ele(aa) * PolyShape(aa, xi(l), 0);
                u_h = u_h + exact(x_ele(aa)) * PolyShape(aa, xi(aa), 0);
                u = exact(x_l);
            end
        
            Error = Error + weight(l) * (u_h - u)^2 * dx_dxi;
            %还没除以分母
            
        end    
    end
    
    %Store the results into one table and plot them
    resultTable = [resultTable; table(hh, Error)];
    plot(log10(resultTable.hh),log10(resultTable.Error),'o-');
    hold on
    xlabel('lg(hh)');
    ylabel('lg(Error)');
    title('Plot of Error vs. Mesh Size');
end








