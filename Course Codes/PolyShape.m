function val = PolyShape(a, xi, der)

if a == 1
    if der == 0
        val = 0.5 * (xi-1) * xi;
    elseif der == 1
        val = 0.5 * (2 * xi-1);
    end
elseif a == 2
    if der == 0
        val = 1 - xi * xi;
    elseif der == 1
        val = -2*xi;
    else
    end
elseif a == 3
    if der == 0
        val = 0.5 * (1+xi) * xi;
    elseif der == 1
        val = 0.5 * (2 * xi+1);
    else
    end
end