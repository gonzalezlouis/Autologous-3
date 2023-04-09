% Monte Carlo scheme: multiple cells
clear all

N = 5;
M = 1;

T = 1e4;
a = 10; % um
D = 145; % um^2/s
v0 = 0.5; % um/s
nu = 3200/3600; % 1/s
kappa = .01;
alpha_ = 0; %(sqrt(17)-1)/4;

chibar = chi_calc(1,0,a,D,0,nu,kappa,alpha_);
% 
rstar = 10;
ustar = -(1/(2*rstar))*10^(4);


% lambda = 10;
lambda = sqrt(-ustar)*rstar;
psi1 = 10;
psi2 = -2*ustar*rstar;

% % Balanced
% lambda = sqrt(-ustar)*rstar;
% psi1 = 10;
% psi2 = -2*ustar*rstar;

% % collective only
% lambda = sqrt(-ustar)*rstar;
% psi1 = 10;
% psi2 = -2*ustar*rstar;

% % individual only
% lambda = 10;
% psi1 = lambda;
% psi2 = -2*ustar*rstar;

d = [ 1  0  0
    -1  0  0
    0  1  0
    0 -1  0
    0  0  1
    0  0 -1];
s_save = zeros(T,1);
l_save = zeros(T,1);
U_save = zeros(T,1);
dU_save = zeros(T,N);
lnn_save = zeros(T,1);

for m = 1:M
    %m
    
    % initialize
    r = [zeros(N,1) zeros(N,1) (0:N-1)'];
    
    for t = 1:T

        
        % record trajectory (m = 1)
        if m == 1
            r_save(t,:,:) = r;
        end
        
        % record cell-averaged and trial-averaged
        %   cell-cell spacing at each time
        s = 0;
        for n = 1:N
            for n_ = n+1:N
                s = s + sqrt(sum((r(n,:)-r(n_,:)).^2));
            end
        end
        s = s/(N*(N-1)/2);
        s_save(t) = ((m-1)*s_save(t) + s)/m;
        
        % record distance of c-o-m from origin at each time
        rcom = mean(r,1);
        l = sqrt(sum(rcom.^2));
        l_save(t) = ((m-1)*l_save(t) + l)/m;
                
        for n = 1:N
            Un = 0;
            for n_ = 1:N
                 if n_ ~= n % for all other cells
                    
                    % calculate concentration
                    dr = r(n,:)-r(n_,:);
                    rho = sqrt(sum((dr).^2));
                    if rho == 0; rho = eps; end
                    theta = acos(dr(3)/rho);
                    chi_n_ = chi_calc(rho,theta,a,D,v0,nu,kappa,alpha_);
                    chi_n_(chi_n_ < 0) = 0; % correct for limitations of perturbation
                    
                    % update potential
                    Un = Un + (lambda/rho)^2 - psi2*(chi_n_/chibar);
                 end
            end

            Uold = Un - psi1;
            
            Upn(n) = Uold;


            
            % choose direction at random
            R = rand;
            j = 0;
            while R > j/6
                j = j + 1;
            end
            r_ = r(n,:) + d(j,:);
            
            % calculate new pseudo-potential
            Up = 0;
            for n_ = 1:N
                  if n_ ~= n % for all other cells
                    
                    % calculate concentration
                    drp = r(n,:)-r(n_,:) + d(j,:);
                    rho = sqrt(sum((drp).^2));
                    if rho == 0; rho = eps; end
                    theta = acos(drp(3)/rho);
                    chi_n_ = chi_calc(rho,theta,a,D,v0,nu,kappa,alpha_);
                    chi_n_(chi_n_ < 0) = 0; % correct for limitations of perturbation

                    % update potential
                    Up = Up + (lambda/rho)^2 - psi2*(chi_n_/chibar);
                 end
            end


            drj = d(j,:);
% 
            rho = sum(((drj).^2));
            if rho == 0; rho = eps; end
            theta = acos(drj(3)/rho);

            c = chi_calc(rho,theta,a,D,v0,nu,kappa,alpha_);

%             %collective chemotaxis terms
%             for n_ = 1:N
%                 
%                 dr_ = r_ - r(n_,:);
%                 rho = sqrt(sum((dr_).^2));
%                 if rho == 0; rho = eps; end
%                 theta = acos(dr_(3)/rho);
%                 c = chi_calc(rho,theta,a,D,v0,nu,kappa,alpha_);
%                 c(c < 0) = 0;
%                 
%             end

            for i = 1:6
                chi(i) = 0;
                for n_ = 1:N
                    dri = r(n,:) - r(n_,:) + d(i,:);
                    rho = sqrt(sum((dri).^2));
                    if rho == 0; rho = eps; end
                    theta = acos(dri(3)/rho);
                    chicalc = chi_calc(rho,theta,a,D,v0,nu,kappa,alpha_);
                    chicalc(chicalc < 0) = 0;
                    chi(i) = chi(i) + chicalc;

                    cbar = mean(chi);
                    Unew = Up - psi1*(c/cbar);
                   
                end
            end

            

            dU = Unew - Uold;

            if dU < 0
                r(n,:) = r_;
            else
                R = rand;
                if R < exp(-dU)
                    r(n,:) = r_;
                end
            end

            % record cell-averaged and trial-averaged potential at each time
            U_save(t) = ((m-1)*U_save(t) + mean(Upn))/m;

%             for i = 1:N
% 
%                 if i ~= n
%                     lnn(t) = (m-1)*(1/N * min(norm(r(n,:)-r(i,:))))/m;
%                 end
%             end

            lnn = 0;

            for n_ = 1:N
                if n ~= n_
                    lnn = lnn + min(norm(r(n,:)-r(n_,:)));
                end
            end

            lnn = lnn/(N*(N-1)/2);
            lnn_save(t) = ((m-1)*lnn_save(t) + lnn)/m;
        end
    end
end

% plot
figure(1); clf
% h = subplot(2,2,1);


plot3(r_save(:,:,3),r_save(:,:,2),r_save(:,:,1),'LineWidth',3)
plot = gcf;
ax = gca;
view(2)
xlabel('z','Fontsize',30)
ylabel('y','Fontsize',30)
ax.XAxis.FontSize = 20;
ax.YAxis.FontSize = 20;
exportgraphics(plot,'../fig/ratio10000zyplot.png','Resolution',300)
% myAxes=findobj(h,'Type','Axes');
% baseFileName = sprintf('%dtrajectory%i.png', psi2/psi1, N);
% fullFileName = fullfile('../fig/', baseFileName);
% exportgraphics(myAxes, fullFileName);

% 
% tvec = [10 100];
% dvec = sqrt(tvec);
% t = 1:T;
% 
% zavg = sum(r_save(:,:,3),2)/N;
% speed = polyfit(t,zavg,1)
% 
% subplot(2,2,2);
% % loglog(1:T,s_save,tvec,dvec)
% % xlabel('time steps')
% % ylabel('avg cell-cell spacing')
% loglog(1:T, lnn_save, tvec, dvec)
% xlabel('time steps')
% ylabel('nearest neighbor spacing')
% lnnbar = mean(lnn_save);
% rho = 1/(lnnbar^3)
% % myAxes=findobj(h,'Type','Axes');
% % baseFileName = sprintf('%dnearestneighbor%i.jpg', psi2/psi1, N);
% % fullFileName = fullfile('../dat/nearestneighborplots/', baseFileName);
% % exportgraphics(myAxes, fullFileName);
% 
% 
% subplot(2,2,3)
% loglog(1:T,l_save,tvec,dvec)
% xlabel('time steps')
% ylabel('avg c.o.m. distance')
% 
% % subplot(2,2,4)
% % plot(1:T,U_save);
% % hold on
% % yline(0)
% % xlabel('time steps')
% % ylabel('average potential')