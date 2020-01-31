function dydt = mlsolve(t,y,ml)

minv = @(v) 0.5*(1+tanh((v-ml.V1)/ml.V2)); %computes ion channel functions
winv = @(v) 0.5*(1+tanh((v-ml.V3)/ml.V4));
lambda = @(v) ml.phi*cosh((v-ml.V3)/(2*ml.V4));

dvdt = (1/ml.C)*(ml.gL*(ml.VL-y(1)) + ml.gK*y(2)*(ml.VK-y(1)) + ml.gCa*minv(y(1))*(ml.VCa - y(1)) + ml.I);
dwdt = lambda(y(1))*(winv(y(1)) - y(2));
dydt = [dvdt;dwdt]; %outputs gradient of dV/dt and dW/dt 
end
