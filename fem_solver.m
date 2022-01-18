% Aga Patro
% Finite Elements Method Project 
% Nb. 1: Heat Transport Equation


% main function which solves our problem 
function fem_solver
    % here you can change input data (number of elements)
    n = 10;
    
    BMatrix = zeros(n+1, n+1);
    for i=0:n-1
        for j=0:n
           if ( abs(i-j) > 1) % if so, integrals will be equal 0 anyway, so i do not count them
               continue 
           end
           from = countLowerLimit(i, j, n);
           to = countUpperLimit(i, j, n);
           % insert counted B(u, v) in right place
           BMatrix(i+1, j+1) = countLeftIntegral(determinatePrimElement(i, n), determinatePrimElement(j, n), from, to)...
               - elementsMultiplication(determinateElement(i, n), determinateElement(j, n));
        end
    end

    % matrix of L(v) functions
    LMatrix = zeros(n+1, 1);
    
    for i=0:n-1
        LMatrix(i+1) = countRight(determinateElement(i, n));
    end

    BMatrix(n+1, n+1) = 1;
    Results = linsolve(BMatrix, LMatrix);
    
    % printing matrixs and results
    disp([BMatrix LMatrix]);  
    disp(Results);

    % making plot 
    X = [0:0.1:2];
    Y = zeros(1, length(X));

    for i=0:length(X)-1
        for j=0:length(Results)-1
            element = determinateElement(j, n);
            Y(i+1) = Y(i+1) + Results(j+1) .* element(X(i+1));
        end
    end

    plot(X, Y, 'm'), xlabel('x'), ylabel('f(x)'), title("equation graph")

end


% functions that help main function solve problem:


% returns integral counted by Legenfre-Gauss Quadrature
function integration = integrate(uPrim_vPrim, from, to)
    integration = ( ...
        (to-from) ./2 .*( ...
        uPrim_vPrim((to-from) ./2 .* (1/sqrt(3)) + (from+to)./2) ...
        + uPrim_vPrim((from-to) ./2 .* (-1/sqrt(3)) + (from+to)./2)...
        )...
        );
end


% returns u(0)*v(0) 
function product = elementsMultiplication(u, v)
    product = u(0) .* v(0);
end


%returns right side of equation
function l = countRight(v)
    l = -20 .* v(0);
end


% returns counted integral from left side of equation 
function integralB = countLeftIntegral(uPrim, vPrim, from, to)
    integralB = integrate (@(x) (uPrim(x) .* vPrim(x)), from, to);
end


% returns i-element function 
function element = determinateElement(i, n)
    element = @(x)(max( 1- abs( (x - i.*2/n) ./2.*n), 0));
end


% returns derivative from i-element function
function elementPrim = determinatePrimElement(i, n)
    elementPrim = @(x) (n./2 .* (2.*(i-1)/n <= x).*(x<2.*i/n) .*(0 <= x) + (-n)./2 .*(2.*i/n <= x) .* (x < 2.*(i+1) / n).*(x <= 2));
end


% returns upper integration limit <to>
function upperLimit = countUpperLimit(i, j, n)
    if (abs(i-j) == 1)
        upperLimit = 2 .* min(1, max(i, j)/n);
    else
        upperLimit = 2 .* min(1, (i+1)/n);
    end
end


% returns lower integration limit <from>
function lowerLimit = countLowerLimit(i, j, n)
    if (abs(i-j) == 1)
         lowerLimit = 2 .* max(0, min(i, j)/n);
    else
         lowerLimit = 2 .* max(0, (i-1)/n);
    end 
end



