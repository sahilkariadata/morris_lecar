data = importdata('pdata.dat'); %imports periodic orbit data from xppaut .dat file
data(data(:,5)~=4,:) = []; %gets rid of all data not periodic orbits
I = data(:,1); %applied current for each orbit
spo = data; 
upo = data;
spo(spo(:,4)~=3,:) = nan; %singles out the stable orbits
upo(upo(:,4)~=4,:) = nan; %singles out the unstable orbits
spo = spo(:,2);
upo = upo(:,2); %period values

plot(I,spo,'b*',I,upo,'r*') %plots the stable and unstable orbits
hold on
line([0.072932 0.072932], [0 3500000],'LineStyle','--'); %plots a dashed line at I = Ic
%closer approximation I_c = 0.07293189486804816 as xppaut rounds upon conversion to .dat file
xlabel('I (uA/cm^2)')
ylabel('Period(ms)')
%ylim([0 3.5])
% title('Plot of period of orbits against applied current I')
legend('Stable Periodic Orbits','Unstable Periodic Orbits')
set(gcf,'color','w');
set(gca,'fontsize',20);
