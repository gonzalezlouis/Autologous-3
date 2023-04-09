clear all

% parameters
a = 10; % um
beta_ = .05;
alpha_ = .75;
gamma = beta_/(1+alpha_);
epsilon = .07;
kappa = 1e-2;

% variables
N = 100;
theta_vec = linspace(0,pi,N);
phi_vec = linspace(0,2*pi,N);
theta = theta_vec'*ones(1,N);
%phi = ones(N,1)*phi_vec;

% parameters/functions for flow
w = 1 + 1/kappa - exp(1/kappa)*expint(1/kappa)/kappa^2;

% with flow, at surface, c(theta,phi) at r = a
c0 = gamma/a^3; % 1/um^3
c0_surf = c0*(1e6)^3/6e23*ones(N,N); % surface of cell, moles/m^3

% with flow, at surface, c(theta,phi) at r = a
chi1 = gamma/2*(alpha_/(1+alpha_) - 1 ...
    + cos(theta)*((1-alpha_)*w/(alpha_+2) + w)/4);
c1 = chi1/a^3;
c = c0 + epsilon*c1; % 1/um^3
c_surf = c*(1e6)^3/6e23; % surface of cell, moles/m^3

figure(1); clf
subplot(2,2,1)
imagesc(theta_vec,phi_vec,c0_surf')
colorbar
xlabel('\theta')
ylabel('\phi')
title('c_0 (mole/m^3)')
set(gca,'ydir','normal')

subplot(2,2,2)
imagesc(theta_vec,phi_vec,c_surf')
colorbar
xlabel('\theta')
ylabel('\phi')
title('c (mole/m^3)')
set(gca,'ydir','normal')

subplot(2,2,3)
imagesc(theta_vec,phi_vec,(c0_surf - c_surf)'/c0_surf(1))
colorbar
xlabel('\theta')
ylabel('\phi')
title('(c_0 - c)/c_0 (mole/m^3)')
set(gca,'ydir','normal')

% numerically integrate to find anisotropy (one factor of cos is from
% Jacobian; other (in numerator only) is from definition of A)
A_num = sum(sum(c_surf.*sin(theta).*cos(theta)))...
    /sum(sum(c_surf.*sin(theta)))

% anisotropy, theory 
A = epsilon*w/8/(2+alpha_)

% COMSOL
dd = 'C:\Users\log24\Dropbox\autologous3\dat\';
D = load([dd 'brinkman-transport3d.dat']); % 'cleaned.dat' has isosurface 10.00001 mum while 'cleaned10.dat' has isosurface 10 mum 
x = D(:,1);
y = D(:,2);
z = D(:,3);
c = D(:,6);
N = length(c);

r = sqrt(x.^2+y.^2+z.^2);
theta = acos(z./r);
phi = atan2(y,x)+pi;

% tile 9 times
theta_ = [];
phi_ = [];
for i = -1:1:1
    for j = -1:1:1
        theta_ = [theta_; theta + i*pi];
        phi_ = [phi_; phi + j*2*pi];
    end
end

% get voronoi vertices, cells, and cell areas of tiled data
[v_,C_] = voronoin([theta_ phi_]);
figure(3); clf
voronoi(theta_,phi_)
a_ = zeros(length(C_),1);
for i = 1:length(C_)
    v1_ = v_(C_{i},1); 
    v2_ = v_(C_{i},2);
    a_(i) = polyarea(v1_,v2_);
end

% untile cell area (take just central square)
a = a_(4*N+1:5*N);

figure(2); clf
subplot(2,2,1)
scatter3(x,y,z,15,c)
colorbar

subplot(2,2,2)
hist(r)

subplot(2,2,3)
hist(theta)

subplot(2,2,4)
hist(phi)

figure(1)
subplot(2,2,4)
scatter(theta,phi,15,c)
colorbar
ylim([0 2*pi])
xlim([0 pi])
xlabel('\theta')
ylabel('\phi')
title('c (mole/m^3)')
set(gca,'ydir','normal')
box on

% numerically integrate to find anisotropy
badA_com = sum(c.*sin(theta).*cos(theta))...
    /sum(c.*sin(theta))

goodA_com = sum(a.*c.*sin(theta).*cos(theta))...
    /sum(a.*c.*sin(theta))

%New Method
A1=scatteredInterpolant(phi,theta,c); %Interpolate the data
s=0;
%Do the integral with all points in between the given ones are now
%interpolated
for i=0:0.1:pi
    for j=0:0.1:2*pi
        s=s+ cos(i)*sin(i)*A1(j,i);
    end
end
 
k=0;
for i=0:0.1:pi
    for j=0:0.1:2*pi
        k=k+sin(i)*A1(j,i);
    end
end
 
newA_com=s/k
