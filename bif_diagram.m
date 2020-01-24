load('mlvars')
h = 1e-6;
nmax = 200;
y0 = [0;0;0.04];
ytan = [-0.1;0;1]; %combine v, w and I
limit = 50;

minv = @(v) 0.5*(1+tanh((v-ml.V1)/ml.V2));
winv = @(v) 0.5*(1+tanh((v-ml.V3)/ml.V4));
lambda = @(v) ml.phi*cosh((v-ml.V3)/(2*ml.V4));
dvdt = @(v,w,I) (1/ml.C)*(ml.gL*(ml.VL-v) + ml.gK*w.*(ml.VK-v) + ml.gCa*minv(v).*(ml.VCa-v) + I);
dwdt = @(v,w) lambda(v).*(winv(v) - w);
dydt = @(y) [dvdt(y(1),y(2),y(3));dwdt(y(1),y(2))];
df = @(y) MyJacobian(dydt,y,h);
ylist = MyTrackCurve(dydt,df,y0,ytan,'s',[1e-6,0.001,0.005],'nmax',nmax,'limit',limit);

xlist = ylist(1:2,:);
I0 = ylist(3,:);
l = length(I0);
sink = zeros(2,l);
saddle = zeros(2,l);
source = zeros(2,l);

dydt2 = @(y,I) [dvdt(y(1),y(2),I);dwdt(y(1),y(2))];
for i = 1:l
    dydt3 = @(y) dydt2(y,I0(i)); %fixes I0 values to retrieve eigenvectors/eigenvalues for each one
    J = @(y)MyJacobian(dydt3,y,h);
    [eV(2*i-1:2*i,:),eD(2*i-1:2*i,:)] = eig(J(xlist(:,i))); % gets eigenvectors
    test(:,i) = (eig(J(xlist(:,i)))); %eigenvalues to test for sinks,saddles and sources
    e = sum(real(test(:,i))>0);
    if  e == 0
        sink(:,i) = xlist(:,i);
    elseif e == 2
        source(:,i) = xlist(:,i);
    elseif e == 1
        saddle(:,i) = xlist(:,i);
    end
end

%% turns zeroes into NaN's so they aren't plotted
sinks = sink(1,:);
sinks(sinks==0) = nan;
saddles = saddle(1,:);
saddles(saddles==0) = nan;
sources = source(1,:);
sources(sources==0) = nan;
%% plots the equilibria in terms of type
plot(I0,sinks,'.',I0,saddles,'x',I0,sources,'s')
xlim([ylist(3,1) ylist(3,l)])
legend({'sinks','saddles','sources'},'Location','southeast')
xlim([0 0.1])
ylim([-0.5 0.2])
xlabel('I0')
ylabel('V')