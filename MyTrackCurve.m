function ylist = MyTrackCurve(userf,userdf,y0,ytan,varargin)
%% uses mysetoptions for optional variable in function
default = {'s',[1e-6,0.01,2],'nmax',200,'limit',0};
options = MySetOptions(default,varargin);
s = options.s;
nmax = options.nmax;
limit = options.limit;
%% initialises variables used
smin = s(1);
smax = s(3);
s = s(2);
d = size(y0);
l = d(1,2);
ypred = y0(:,l); % + s*ytan (initial step value 0)
ylist(:,1) = y0;
i=1; % initialises for loop to work

while  i < nmax
    %% checks variable is in limit to prevent unnecessary computations
    if limit ~= 0
        if ylist(3,i) > limit
            ylist = ylist(:,1:i-1);
            return;
        end
    end
    %% solves for ylist
    F = @(y)[userf(y);(ytan.')*(y-ypred)];
    DF = @(y)[userdf(y);ytan.']; %Defines the functions
    [ylist(:,i+1),conv,J] = MySolve(F,ylist(:,i),DF,(1e-6),100);
    %% checks convergence for variable step-size
    if conv == 0
        sold = s;
        s = s/2; 
        if s < smin
            break;
        end
        ypred = ylist(:,i) + s*ytan;
    else
        %% step has converged so assigns variable values for next step
        s = min([1.2*s,smax]);
        [~,n] = size(J);
        z1 = zeros(n,1);
        z1(n,1) = 1;
        Z = (J^-1)*z1;
        ytan = (Z/norm(Z,Inf))*(sign((Z.')*ytan));
        ypred = ylist(:,i+1) + s*ytan;
        i = i+1; %goes next iteration when successful
    end
end
end
