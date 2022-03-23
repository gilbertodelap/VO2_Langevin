%% test kinetic energy


addpath(genpath('~/Dropbox/xdscode/'))
% addpath('~/Documents/CDWs/SwissFEL_Feb2021/remodelGL/')
addpath(genpath('~/Dropbox/include/matlab/xspde_matlab/XSPDE_code/'));

cd ~/Documents/VO2/simu/2TMD/dynamics/Jan2022/

%% define paths

% addpath(genpath('~/Dropbox/xdscode/'))
% addpath('~/Documents/CDWs/SwissFEL_Feb2021/remodelGL/')
addpath(genpath('~/Documents/DATA/2022/Jan2022_model/xspde_matlab/XSPDE_code/'));
addpath(genpath('~/Documents/DATA/2020/VO2_2015_data/xdscode/'))
% cd ~/Documents/VO2/simu/2TMD/dynamics/Jan2022/
cd ~/Documents/DATA/2022/Jan2022_model/untitled_folder/untitled_folder/

%% parameters to simulate the dynamics in 2D %%%%%%%%%%%%%%%%%%%%%%%%%%%
Lx = 90;
Ly = 60;
T = 3;
Nx = Lx;
Ny = Ly;
Nt = 100;
%
clear opts videofile
opts.ranges = [T, Lx, Ly];
opts.npoints = [Nt, Nx, Ny];
opts.nfields = 4;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% the potential paramters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% V(psi) = 0.5 a psi^2 + 0.25 b psi^4;
a       = 1.6;              % psi^2 term [eV/Angstrom^2]
b       = 25.6;             % psi^4 term [eV/Angstrom^4]
m       = 1.4e-25;          % mass [kg]
decay   = 7.69e-25;         % decay const. [kg/ps]
Vmin    = -0.25*a^2/b;      % potential energy at minimum [eV]
psi0    = sqrt(a/b);        % position of the minimum [Angstroms]
freq    = sqrt(2*1.602e-19*1e20*a/m)/pi/2*1e-12; % 
decay/m;  % 5.5 THz
% opts.gam = decay/m/opts.a;  %  [THz/THz^2 = ps ]

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% setup the options structure
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

opts.tauf = 1.2;
opts.xi0 = 0.2*Lx;       % pump penetration depth [nm]
opts.xi0xray = 1*Lx;     % xray penetration depth [nm]
opts.kT =1.1116e-04*100;        % kT (units?)
opts.cohlengths2 = 1e-2*[1.2^2,1.2^2];  % nm^2
opts.cohlengths2 = 20.6e-2*[1.2^2,1.2^2];  % nm^2
opts.a = (2*pi*3)^2/2;  % should be a/m
% opts.a = 303; %THz^2 from 1.6eV/A^2/(Mass of vanadium)
opts.gam = 0.03;        % [ps]

dt = opts.ranges(1)/opts.npoints(1);
dx = opts.ranges(2)/opts.npoints(2);
dy = opts.ranges(3)/opts.npoints(3);
t = 0:dt:opts.ranges(1);  t=t(1:end-1);
x = 0:dx:opts.ranges(2);  x=x(1:end-1);
y = 0:dy:opts.ranges(3);  y=y(1:end-1);
xs{1} = x(:);
xs{2} = y(:);
[xx,yy]=ndgrid(x,y);

% define the Laplacian operator:
opts.Laplacian = make_Laplacian_operator(opts);

% setup the staggered initial condition 
X00 = zeros(Nx, Ny);
for ii=1:opts.npoints(2)
    for jj=1:opts.npoints(3)
        X00(ii,jj) = (-1)^(ii+jj);
    end
end
X0 = sqrt(1)*[  X00(:)', zeros(1,prod(opts.npoints(2:end))),...
        X00(:)', zeros(1,prod(opts.npoints(2:end)))];
    

opts.r = @(xs,t) - 1;                           % constant quadratic coefficient
opts.sig = @(xs,t) sqrt(2*opts.gam*opts.kT);   % (2*opts.gam*opts.kT) % constant temperature


%% thermalize
[F,G] = make_update_function_VO2_tdeptemp(opts);
dd = integrate_with_XSPDE(X0, F, G, opts);
X = permute(dd, [2, 3, 1, 4]);
X0 = X(:,:,end,:);
X0 = X0(:)';

%% compute the kinetic energy
[F,G] = make_update_function_VO2_tdeptemp(opts);
dd = integrate_with_XSPDE(X0, F, G, opts);
X = permute(dd, [2, 3, 1, 4]);
X0 = X(:,:,end,:);
X0 = X0(:)';

%%
% psi0=0.25;
% m = 1.4e-25; %kg
kgTHz2_to_eVperAngstroms = 6.242e22;

% p  = 1/~a[THz^2] d~psi/dt[THz]
% psi0 opts.a * p [L/T]
p2 = mean(reshape(mean(X(:,:,:,[2,4]).^2,4),[],1));  % p^2 [ps^2]
% std(reshape(sum(X(:,:,:,[2,4]).^2,4),[],1))
KEunits = m*psi0^2*opts.a^2;
KEunits = KEunits*kgTHz2_to_eVperAngstroms;
% was KE = 18.3[eV/ps^2]*p2[ps^2] % KE = (psi0*a)^2*p^2/m = (0.4 eV/A)^2/1.4e-25kg * p2[ps^2]

fprintf('\nKinetic energy (per DOF) = %f meV\n', 1e3*KEunits*p2/2)
fprintf('Temperature = %f K\n', 2*11605*KEunits*p2)

% compute potential energy (normalized units)
x2 = mean(reshape(sum(X(:,:,:,[1,3]).^2,4)/2,[],1));
x4 = mean(reshape(sum(X(:,:,:,[1,3]).^4,4)/2,[],1));
Vpot = 0.25*a*psi0^2*(-2*x2 + x4);
fprintf('Potential energy = %f meV\n\n', 1e3*kgTHz2_to_eVperAngstroms*0.25*opts.a*m*psi0^2*Vpot)
% 
% def: p = \dot{~psi}/~a 
% def: P = m dot{psi} = a psi0 1/~a dot{\~psi} = a psi0 p
% KE = (a psi0)^2 p^2/m = 18.3 eV/ps^2 p^2

% compute velocities directly (finite differences)
kg_Angstrom_per_ps2_to_ev = 6.262e22;
v2 = psi0*sqrt(mean(reshape(mean(diff(X(:,:,:,[1,3])/dt,[],3).^2,4),[],1)));
fprintf('Average velocity = %f Angstroms/ps\n', v2)
fprintf('KE from average velocity = %f meV\n\n', 0.5*1e3*v2*m*kg_Angstrom_per_ps2_to_ev)
% Wolfram alpha gives kT = m <v>^2 -> v ~ 2.18 Angstroms/ps



%% 2D potential in relative units
V2D = @(x,y) 0.5*0.025*(-2*(x.^2 + y.^2) + x.^4+y.^4)
xx = linspace(-1.5, 1.5, 50);
yy = linspace(-1.5, 1.5, 50);

figure(15);clf
imagesc(xx,yy,V2D(xx,yy'));
axis image; colorbar
figure(16); clf
plot(xx, V2D(xx,xx))
