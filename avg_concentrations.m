clear all
close all

dd = 'C:\Users\log24\Dropbox\autologous3\';
pos = load([dd 'dat\cellcoordinates.csv']); % um
D = load([dd 'dat\collective-transport3d.dat']);
x = D(:,1); % um
y = D(:,2); % um
z = D(:,3); % um
c = D(:,6); % mol/m^3

a = 10; % um
[N,d] = size(pos);
f = 0.0001;
for i = 2:N
    x0 = pos(i,1);
    y0 = pos(i,2);
    z0 = pos(i,3);
    r = sqrt((x-x0).^2+(y-y0).^2+(z-z0).^2);
    j = find(abs(r-a)/a < f);
    cavg(i-1) = mean(c(j));
    numpts(i-1) = length(j);
    theta(i-1) = acos(z0/sqrt(x0^2+y0^2+z0^2));
    phi(i-1) = atan2(y0,x0)+pi;
end

cavg'
% numpts'
% theta'
% phi'

Acol = sum(cavg.*cos(theta))/sum(cavg)

figure(1); clf
scatter3(pos(2:N,1),pos(2:N,2),pos(2:N,3),100,cavg)
colorbar
xlabel('x')
ylabel('y')
zlabel('z')