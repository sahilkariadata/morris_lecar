data = importdata('pdata.dat');
data(data(:,5)~=4,:) = [];
I = data(:,1);
spo = data;
upo = data;
spo(spo(:,4)~=3,:) = nan;
upo(upo(:,4)~=4,:) = nan;
spo = spo(:,2);
upo = upo(:,2);

plot(I,spo,'b*',I,upo,'r*')
xlabel('I')
ylabel('Period')
legend('Stable Periodic Orbits','Unstable Periodic Orbits')
