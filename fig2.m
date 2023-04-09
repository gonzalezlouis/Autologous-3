clear all

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
a = 10^-5; %m
n = 100;
eps = 0.07;
eps3 = 3*eps;
w = 2;
alpha = 0.75;
alpha0 = 0;
width = 3*10^-3; %m
rho = logspace(6,16,100);
number = floor(linspace(1,600,20)');
y = floor(logspace(log10(10),log10(600),20)');
vol = (3000*2000*100)*10^-8;

N = 30*ones(1,22)';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

H = readtable('../dat/celldensitysimulationcomsol.xlsx');
G = readtable('../dat/averagecollectiveanisotropy.xlsx');
J = readtable('../dat/comsolclusteranisotropy-individual.xlsx');
K = readtable('../dat/comsolclusteranisotropy-average.xlsx');

collectivenumber_ = table2array(G(:,1));
collectiveanisotropy = table2array(G(:,7));
collectiveerror = table2array(G(:,8));
comsol_number_ = table2array(H(:,1));
comsol_anisotropy = table2array(H(:,7));
comsol_error = table2array(H(:,8));
comsol_anisotropyalpha0001 = table2array(H(:,31));
comsol_erroralpha0001 = table2array(H(:,32));
clusteranisotropylambda = table2array(J(:,1));
clusteranisotropyindividual = table2array(J(:,7));
clusteranisotropyindividualerror = table2array(J(:,8));
clusteranisotropyaveragelambda = table2array(K(:,1));
clusteranisotropyaverage = table2array(K(:,7));
clusteranisotropyaverageerror = table2array(K(:,8));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clustercollectiverho = 1e18 * 30 * (4/3 * pi .* clusteranisotropyaveragelambda.^3).^(-1);
comsol_number = ((comsol_number_ ./ (3000*2000*100)) .* 10^9); % mm^-3;
collective_number = 1e9 * (collectivenumber_ ./ (3000*2000*100)); % mm^-3

clusterrho = 1e18 * (4/3 * pi .* clusteranisotropylambda.^3).^(-1);

A = (eps3*w/(8*(2 + alpha))) .* (eps3*a)./ ... 
    (eps3.*a + 4*pi*width*(1 + alpha).*rho.*(a^3));

A_3alpha0 = (eps3*w/(8*(2 + alpha0))) .* (eps3*a)./ ... 
    (eps3.*a + 4*pi*width*(1 + alpha0).*rho.*(a^3));

A_col =  4/3 * pi *vol^(1/3) * eps * a^3 .* rho.^(4/3);

A_colavg = 4/3 * pi * 30^(1/3) * 3 * eps * a^3 * rho;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fakerho1ind = logspace(0.2,1.2,10); 
y1 = linspace(0.035,0.035,10);
fakerho2ind = 8*logspace(1.2,1.8,10);
y2 = 6*logspace(-2.5,-3,10);
fakerho1col = logspace(2,2.65,10);
y3 = 5*logspace(-2,-1.2,10);

figure(1); clf
% plot(rho,A_3alpha0,'k','DisplayName','analytic A_{ind}, \alpha = 0');
% hold on
% plot(rho,A_col,'r','DisplayName','analytic A_{col}, \alpha = 0');
% hold on
errorbar(collective_number(7:end),collectiveanisotropy(7:end),collectiveerror(7:end),'ro', ...
    'DisplayName','comsol A_{col,box}, \alpha = 0.0001',LineWidth=2,MarkerSize=10)
hold on
errorbar(comsol_number,comsol_anisotropyalpha0001, comsol_erroralpha0001, 'bo', ...
    'DisplayName','comsol A_{ind,box}, \alpha = 0.0001',LineWidth=2,MarkerSize=10)
hold on
plot(fakerho1ind,y1,'b-',LineWidth=2)
hold on
text(fakerho1ind(4),y1(1)+0.02,'~\rho^0','FontSize',40,'Color','b')
hold on
plot(fakerho2ind,y2,'b-',LineWidth=2)
hold on
text(fakerho2ind(4),y2(1)+0.001,'~\rho^{-1}','FontSize',40,'Color','b')
hold on
plot(fakerho1col,y3,'r-',LineWidth=2)
hold on
text(fakerho1col(4),y3(3) + 0.00001,'~\rho^{4/3}','FontSize',40,'Color','r')
hold on
text(2,0.6,'Individual, A_{I}','FontSize',40,'Color','b')
hold on
text(1.9,0.3,'Collective, A_{C}','FontSize',40,'Color','r')
hold on
[xl,xt] = xlin1(5.57*10^1,'\rho_c', 10^-3, 2*10^-2, 4*10^-3);
set(gca, 'XScale', 'log', 'YScale', 'log');
ax = gca;
ax.XAxis.FontSize = 35;
ax.YAxis.FontSize = 35;
xlabel('Cell density, $\rho$ [cells/mm$^{3}$]','Interpreter','Latex', ... 
        'fontsize',40);
ylabel('Anisotropy, A','Interpreter','Latex','fontsize',40)
set(gcf,'Position',[100 100 720 720])

function [hl,ht] = xlin1(x,txt,ylo,yhi,ytxt) 
% Documentation: 
%   x       = x-Position
%   txt     = Text String
%   ylo     = Low y-Value (Start)
%   yhi     = High y-Value (End)
%   ytxt    = Text Starting Position
hold on
hl = plot([x x],[ylo yhi],'DisplayName',txt, 'LineWidth',2,'LineStyle','--','Color','b');
ht = text(4.7*10^1,ytxt, txt, 'Horiz','right', 'Vert','top', 'Rotation',0,'FontSize',40,'Color','b');
hold off
end

function [hl,ht] = xlin2(x,txt,ylo,yhi,ytxt) 
% Documentation: 
%   x       = x-Position
%   txt     = Text String
%   ylo     = Low y-Value (Start)
%   yhi     = High y-Value (End)
%   ytxt    = Text Starting Position
hold on
hl = plot([x x],[ylo yhi],'DisplayName',txt, 'LineWidth',2,'LineStyle','--','Color','r');
ht = text(x,ytxt, txt, 'Horiz','left', 'Vert','top', 'Rotation',270,'FontSize',20,'Color','r');
hold off
end
