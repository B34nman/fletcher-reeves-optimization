warning off;

x = [-1.2;1];
A = 100;

function fletcher_reeves(x, A)

    storage = zeros(2, 6000); %storage for path plot
    c = 0.01;
    rho = .5;
    f=func(x, A);
    g=grad(x, A);
    k = 0; % k = # iterations
    funcEval=1;	% funcEval = # function eval.	


    while 

        %this is our step size calculation
        pk = -g; % steepest descent direction
        a = 1;
        newf = func(x + a*pk, A);
        funcEval = funcEval+1;
        while (newf > f + c*a*g'*pk)
            a = a*rho;
            newf = func(x + a*pk, A);
            funcEval = funcEval+1;
        end

    end
    

end

function y = func(x, A)
    y = A*(x(1)^2 - x(2))^2 + (x(1)-1)^2;
end

function y = grad(x, A)
    y(1) = A*(2*(x(1)^2-x(2))*2*x(1)) + 2*(x(1)-1);
    y(2) = A*(-2*(x(1)^2-x(2)));
    y = y';
end

function y = hessian(x, A)
    y(1,1) = (12*A) * x(1)^2 - (4*A) * x(2) + 2;
    y(1,2) = -(4*A) * x(1);
    y(2,1) = -(4*A) * x(1);
    y(2,2) = (2*A);
end