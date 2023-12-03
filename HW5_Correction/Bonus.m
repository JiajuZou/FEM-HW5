% clear the memory and the screen
clear; clc;

% exact solution
exact = @(x) sin(x);

% problem definition
f = @(x) sin(x);
g = sin(1);
h = -cos(0);

n_ele = [2 4 8 16 32 64 128 256 512 1024];%依次增大单元数目进行计算
ErrorCell = cell(length(n_ele),1); %用于存储不同单元数目下计算的误差
for n = 1:length(n_ele)
% generate my mesh
n_el = n_ele(n);
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

%% Bonus
%Define the parameters of gmres
restart = 10000;
maxit = 10000;

% Set the value of tol
tol_values = [1e-2, 1e-4, 1e-6];

uh1 = gmres(K,F,restart,tol_values(1),maxit);
uh2 = gmres(K,F,restart,tol_values(2),maxit);
uh3 = gmres(K,F,restart,tol_values(3),maxit);

%计算不同单元数目下的误差并储存
solution = [uh uh1 uh2 uh3;g g g g];
for i = 1:4
   ErrorCell{n}(i,:) = solution(:,i)' - (exact(x_coor));
end

end

%通过成倍地增加单元个数进行计算，发现：
%在单元数目较少时，采用不同tol计算得到的结果与直接计算得到的结果相差不大

%逐渐增大单元数目，采用10e-2 tol计算得到的结果与直接计算结果相比出现明显偏差
%此时采用10e-4和10e-6 tol计算得到的结果与直接计算结果相比误差更小

%继续增大单元数目，采用10e-4 tol计算得到的结果也出现和直接计算结果相比出现明显偏差
%此时直接计算结果是这4种结果里面误差最小的

%出现上述情况的原因：
%指定restart = maxit = 10000，总迭代次数不超过其两者之和
%指定容差分别为10e-2,10e-4,10e-6
%在单元数目较小时，直接计算和不同容差下的gmres方法计算结果相同
%逐渐增大单元数目，容差较小的gmres方法在不超过总迭代次数下结果最先出现偏差，因为给定的容差为10e-2
%继续增大单元数目，在不超过总迭代次数的情况下，最小容差情况10e-6计算所得的结果也会次于直接计算结果
%说明在此单元数目下，给定的容差和迭代次数无法满足计算要求，此时可以通过继续缩小容差tol的值来实现更高的计算精度









