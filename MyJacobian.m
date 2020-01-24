function df = MyJacobian(f,x,h)
len = length(x);
xs = x;
xa = x;
%% carries out centred difference to evaluate f'(x)
for i = 1:len
    xs(i) = xs(i) + h; 
    xa(i) = xa(i) - h;
    df(:,i) = (f(xs) - f(xa))/(2*h); 
    xs(i) = x(i);
    xa(i) = x(i);
end
end