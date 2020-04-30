clear all %#ok<CLALL>

addpath('./bsweep_code')
savepath = '../Data/stable_bsweep_bandstructure/miniBand/';


% Bandstructure inputs
filename = 'TabEfiner.mat';

%Ef = -135;
% Ef = -150;
Ef = -188.5; % meV
nband = 3;
Kinit = [0.012220543885626;-0.026];
%Kinit = [-0.010114671119369;-0.019825000000000];

phi = 0*pi/6;
xi = 1; % Valley index +/- 1
h = 0.0025;
%h =  0.01;
delta = 10^-8;
nstep = 1000;


%% Generate closed caustic loop

% For positive fields
mbcp = miniBandCaustic;
mbcp.loadmesh(filename);
mbcp.Ef = Ef;
mbcp.nband = nband;
mbcp.phi = phi;
mbcp.xi = xi;
mbcp.B = 1;
mbcp.h =  h;
mbcp.delta = delta;
mbcp.nstep = nstep;
mbcp.plotBands(2,3)


Kcut = mbcp.FermiSurfaceArea();

areaBZ = polyarea(mbcp.mBZ(1,:),mbcp.mBZ(2,:));
x = mbcp.mBZ(1,:);
y = mbcp.mBZ(2,:);

%%
% Closed Fermi Surface
% areaFS = polyarea(Kcut{nband,1}(1,:),Kcut{nband,1}(2,:));
% area = 2*areaBZ-areaFS;
% e = 1.602E-19; % Charge of electron
% C = 5.51E-4; % F/m^2
% n = (14/12.85)^2*area*1E16/pi^2;
% 
% V = e*n*1E4/C;
% display(V)

%%
% Open Fermi Surface
% for j = 1:size(Kcut,2)
%     k = dsearchn(mbcp.mBZ',Kcut{nband,j}(:,end)');
%     %Kcut{nband,j} = [Kcut{nband,j}(:,1:end), mbcp.mBZ(:,k),Kcut{nband,j}(:,1)];
%     Kcut{nband,j} = [Kcut{nband,j}(:,1:end), mbcp.mBZ(:,8-2*j),Kcut{nband,j}(:,1)];
%     figure(3); hold on
%     plot(Kcut{nband,j}(1,:),Kcut{nband,j}(2,:))
%     [x,y] = poly2cw(x,y);
%     [x,y] = polybool('subtraction',x,y,Kcut{nband,j}(1,:),Kcut{nband,j}(2,:));
%     figure(4); clf; hold on;
%     plot(x,y)
% end

j = 1;
Kcut{nband,j} = [Kcut{nband,j}(:,1:end), mbcp.mBZ(:,4),Kcut{nband,j}(:,1)];
figure(3); hold on
plot(Kcut{nband,j}(1,:),Kcut{nband,j}(2,:))
[x,y] = poly2cw(x,y);
[x,y] = polybool('subtraction',x,y,Kcut{nband,j}(1,:),Kcut{nband,j}(2,:));
figure(4); clf; hold on;
plot(x,y)

j = 2;
Kcut{nband,j} = [Kcut{nband,j}(:,1:end), mbcp.mBZ(:,6),Kcut{nband,j}(:,1)];
figure(3); hold on
plot(Kcut{nband,j}(1,:),Kcut{nband,j}(2,:))
[x,y] = poly2cw(x,y);
[x,y] = polybool('subtraction',x,y,Kcut{nband,j}(1,:),Kcut{nband,j}(2,:));
figure(4); clf; hold on;
plot(x,y)

j = 3;
Kcut{nband,j} = [Kcut{nband,j}(:,1:end), mbcp.mBZ(:,2),Kcut{nband,j}(:,1)];
figure(3); hold on
plot(Kcut{nband,j}(1,:),Kcut{nband,j}(2,:))
[x,y] = poly2cw(x,y);
[x,y] = polybool('subtraction',x,y,Kcut{nband,j}(1,:),Kcut{nband,j}(2,:));
figure(4); clf; hold on;
plot(x,y)



[x,y] = poly2cw(x,y);
figure(4); clf; hold on;
plot(x,y)
area = polyarea(x,y); % A^-2

area = 2*areaBZ-area;

e = 1.602E-19; % Charge of electron
C = 5.51E-4; % F/m^2
n = (14/12.85)^2*area*1E16/pi^2;

V = e*n*1E4/C;
display(V)

%%
n = (14/12.85)^2*areaBZ*1E16/pi^2;

VsDP = e*n*1E4/C;
display(VsDP)





