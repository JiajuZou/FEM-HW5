% clear the memory and the screen
clear all; clc;

%Problem defination
f = @(x) sin(x);
g = sin(1);
h = -cos(0);

% exact solution
exact = @(x) sin(x);
exact_square = @(x) sin(x).^2;
exact_dx = @(x) cos(x);
exact_dx_square = @(x) cos(x).^2;
% integral the exact solution and its derivitive
Error_L2_down = sqrt(integral(exact_square, 0, 1)); % L2分母部分
Error_H1_down = sqrt(integral(exact_dx_square, 0, 1)); % H1分母部分

% creat two tables to store the results
resultTable_L2 = table();
resultTable_H1 = table();


for n_el = 2:2:16 % generate different mesh
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
    n_int = 30;
    [xi, weight] = Gauss(n_int, -1, 1);
    
    
   %% Assembly of K and F and solve Kd=F
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


    %% Calculate the Error L2 and H1
    Error_L2 = 0.0;
    Error_H1 = 0.0;
    for ee = 1 : n_el
        
        x_ele = zeros(4,1); % one element with four points
        for aa = 1 : 4
            x_ele(aa) = x_coor(IEN(aa,ee)); % A = IEN(a,e)
        end
       
        AA = 0.0; 
        BB = 0.0;
        for l = 1 : n_int
            u_h = 0.0;
            u_h_dx = 0.0;
            x_l = 0.0;
            dx_dxi = 0.0;
            
            for aa = 1 : 4
                dx_dxi = dx_dxi + x_ele(aa) * PolyShape(aa, xi(l), 1);
            end
            
            for aa = 1 : 4
%                 u_h = u_h + exact(x_ele(aa)) * PolyShape(aa, xi(l), 0); 
                u_h = u_h + d(IEN(aa,ee)) * PolyShape(aa, xi(l), 0); 
                
                
                u_h_dx = u_h_dx + exact(x_ele(aa)) * PolyShape(aa, xi(l), 1) / dx_dxi;
      
                x_l = x_l + x_ele(aa) * PolyShape(aa, xi(l), 0);
                 
            end
            u_exact = exact(x_l);
            AA = AA + weight(l) * (u_h - u_exact)^2 * dx_dxi;
            % Error_L2 = Error_L2 + weight(l) * (u_h - u_exact)^2 * dx_dxi;
            % We can obtain the same results with the above equation only.
            % AA is just here for helping understanding the loop.
            
            u_exact_dx = exact_dx(x_l);
            BB = BB + weight(l) * (u_h_dx - u_exact_dx)^2 * dx_dxi;
            
        end
        
        Error_L2 = Error_L2 + AA;
        Error_H1 = Error_H1 + BB;
    end
    
    Error_Final_L2 = Error_L2^0.5 / Error_L2_down;
    Error_Final_H1 = Error_H1^0.5 / Error_H1_down;
    %Store the results into the table
    resultTable_L2 = [resultTable_L2; table(hh, Error_Final_L2)];
    resultTable_H1 = [resultTable_H1; table(hh, Error_Final_H1)];
    
end

    %Plot Error L2 and H1
    figure
    plot(log(resultTable_L2.hh),log(resultTable_L2.Error_Final_L2),'o-')
    hold on
    xlabel('log(hh)');
    ylabel('log(Error L2)');
    title('Plot of Error L2 vs. Mesh Size - Cubic');
    
    figure
    plot(log(resultTable_H1.hh),log(resultTable_H1.Error_Final_H1),'o-')
    hold on
    xlabel('log(hh)');
    ylabel('log(Error H1)');
    title('Plot of Error H1 vs. Mesh Size - Cubic');








