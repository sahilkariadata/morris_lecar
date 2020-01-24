function [x,conv,J] = MySolve(f,x0,df,tol,maxit)
if maxit ~= 0 %if maxit is zero the convergence hasn't worked
    J = df(x0);
    x1 = x0 - J\f(x0);
    if norm(abs(x1-x0)) && norm(f(x0)) < tol
        conv = 1;
        x = x1;
        return
    else
        maxit = maxit - 1; %uses recursion by decreasing maxit
        [x,conv,J] = MySolve(f,x1,df,tol,maxit); %uses newest approximation
    end
else
    conv = 0;
    x = 0;
    J = 0; 
    return
end
end