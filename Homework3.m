%warning off;
format long;

x = [-1.2;1];
A = 100;

fletcher_reeves(x, A);

function fletcher_reeves(x, A)

    %storage = zeros(2, 6000); %storage for path plot
    c = 0.01;
    rho = .5;
    f=func(x, A);
    g=grad(x, A);
    pk = -g; % steepest descent direction
    k = 0; % k = # iterations
    funcEval=1;	% funcEval = # function eval.	


    while g ~= 0

        %this is our step size calculation
        a = 1;
        newf = func(x + a*pk, A);
        funcEval = funcEval+1;
        while (newf > f + c*a*g'*pk)
            a = a*rho;
            newf = func(x + a*pk, A);
            funcEval = funcEval+1;
        end

        x = x + a*pk; % gradient descent
        gNew=grad(x, A); %gk+1 = grad(f(x+1))
        
        %beta evaluation
        top = dot(transpose(gNew), gNew);
        bot = dot(transpose(g), g);
        beta = top/bot;

        %update g value
        g = gNew;

        %new pk value
        pk = -1*g + beta*pk;
        if (dot(transpose(pk), g) > 0)
            pk = -g;
        end

        %iterations
        k = k + 1;

        

    end
    
    disp(x);

end

function y = func(x, A)
    y = A*(x(1)^2 - x(2))^2 + (x(1)-1)^2;
end

function y = grad(x, A)
    y(1) = A*(2*(x(1)^2-x(2))*2*x(1)) + 2*(x(1)-1);
    y(2) = A*(-2*(x(1)^2-x(2)));
    y = y';
end