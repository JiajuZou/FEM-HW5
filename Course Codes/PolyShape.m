function val = PolyShape(a, xi, der)

if a == 1
    if der == 0
        val = (xi^3-xi^2-1/9*xi+1/9)/(-16/9);
    elseif der == 1
        val = (3*xi^2-2*xi-1/9)/(-16/9);
    end
elseif a == 2
    if der == 0
        val = (xi^3-1/3*xi^2-xi+1/3)/(16/27);
    elseif der == 1
        val = (3*xi^2-2/3*xi-1)/(16/27);
    else
    end
elseif a == 3
    if der == 0
        val = (xi^3+1/3*xi^2-xi-1/3)/(-16/27);
    elseif der == 1
        val = (3*xi^2+2/3*xi-1)/(-16/27);
    else
    end
elseif a == 4
    if der == 0
        val = (xi^3+xi^2-1/9*xi-1/9)/(16/9);
    elseif der == 1
        val = (3*xi^2+2*xi-1/9)/(16/9);
    else
    end
end