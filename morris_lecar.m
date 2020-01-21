%% initialised variables
ml.I = 0.075;
ml.C = 1;
ml.VK = -0.7;
ml.VL = -0.5;
ml.VCa = 1;
ml.gK = 2;
ml.gL = 0.5;
ml.gCa = 1.33;
ml.V1 = -0.01;
ml.V2 = 0.15;
ml.V3 = 0.1;
ml.V4 = 0.145;
ml.phi = 1/3;
%% initial conditions
v0 = -60;
w0 = 0;
y0 = [v0;w0];
%% solve the DE's using ode45
t0 = 0;
tend = 300;
tspan = [t0 tend];

options = odeset('Abstol',1e-6,'RelTol',1e-6);
[t,y] = ode45(@mlsolve,tspan,y0,options,ml);
%% plot both v and w against time
figure(1);
plot(t,y(:,1),'k');
title('Solution of the Morris-Lecar model');
xlabel('Time (ms)');
ylabel('Voltage (mV)');

figure(2);
plot(t,y(:,2),'r');
title('Solution of the Morris-Lecar model');
xlabel('Time (ms)');
ylabel('W');
%% set up equations for phase plane portrait
minv = @(v) 0.5*(1+tanh((v-ml.V1)/ml.V2));
winv = @(v) 0.5*(1+tanh((v-ml.V3)/ml.V4));
lambda = @(v) ml.phi*cosh((v-ml.V3)/(2*ml.V4));
dvdt = @(v,w) (1/ml.C)*(ml.gL*(ml.VL-v) + ml.gK*w.*(ml.VK-v) + ml.gCa*minv(v).*(ml.VCa-v) + ml.I);
dwdt = @(v,w) lambda(v).*(winv(v) - w);
[X,Y] = meshgrid(-0.8:0.1:0.8);
DV = dvdt(X,Y);
DW = dwdt(X,Y);

figure(3);
plot(y(100:end,1),y(100:end,2),'b');
title('attempted phase plane portrait at I = 0.075');
xlim([-0.5 0.3]);
ylim([-0.1 0.8]);
hold on
quiver(X,Y,DV,DW)
hold on 
fcontour(dvdt,[-1 1 -1 1],'-r')
hold on
fcontour(dwdt,[-1 1 -1 1],'-g')
hold off