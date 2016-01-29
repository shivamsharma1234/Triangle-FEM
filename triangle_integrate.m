function [ value ] = triangle_integrate(f,degree)
    syms x y z;
    f(x, y, z) = matlabFunction(f(x, y, z));
    value = 0;

%    
% 
%     try
%         x_degree = length(sym2poly(f(x, 1))) - 1;
%         y_degree = length(sym2poly(f(1, y))) - 1;
% 
%         degree = max(x_degree, y_degree);
%     catch
%         degree = 2;
%     end
% 
%     if degree == 0
%         degree = degree + 1;
%     end
% 

    if (degree == 1)
        W(1) = 1;
        x_values(1, :) = [1/3 1/3 1/3];
    end
    
    if (degree == 2)
        W(1) = 1/3;
        W(2) = 1/3;
        W(3) = 1/3;
        
        x_values(1, :) = [1/2 0 1/2];
        x_values(2, :) = [1/2 1/2 0];
        x_values(3, :) = [0 1/2 1/2];
    end

    for i=1:length(W)
        value = value + 0.5*W(i)*f(x_values(i, 1), x_values(i, 2), x_values(i, 3));
    end

end
