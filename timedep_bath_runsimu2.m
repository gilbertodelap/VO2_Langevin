
%% define paths

% addpath(genpath('~/Dropbox/xdscode/'))
% addpath('~/Documents/CDWs/SwissFEL_Feb2021/remodelGL/')
% addpath(genpath('~/Dropbox/include/matlab/xspde_matlab/XSPDE_code/'));

% cd ~/Documents/VO2/simu/2TMD/dynamics/Jan2022/
cd ~/Documents/DATA/2022/VO2_Langevin/
addpath(genpath('./xspde_matlab/XSPDE_code/'));
addpath(genpath('~/Documents/xdscode-master/'))

%% parameters to simulate the dynamics in 2D %%%%%%%%%%%%%%%%%%%%%%%%%%%
Lx = 90*4;
Ly = 60;
T = 6;
Nx = Lx;
Ny = Ly;
Nt = 600;
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
freq    = sqrt(2*1.602e-19*1e20*a/m)/pi/2*1e-12; % [THz]
% decay/m;  % 5.5 THz
% opts.gam = decay/m/opts.a;  %  [THz/THz^2 = ps ]

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% setup the options structure
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

opts.a = (2*pi*3)^2/2;  % should be a/m
% opts.cohlengths2 = 15.6e-2*[1.2^2,1.2^2];  % nm^2
opts.cohlengths2 = 20.6e-2*[1.2^2,1.2^2];  % nm^2
% opts.cohlengths2 = 1e-2*[1.2^2,1.2^2];  % nm^2

opts.gam = 0.03;        % [ps]
% opts.tauf = 1.2;        % [ps]
opts.tauf = 1.2;        % [ps]
% opts.xi0 = 0.2*Lx;       % pump penetration depth [nm] 
opts.xi0 =.1* Lx;       % pump penetration depth [nm] 
% opts.xi0xray = 1*Lx;     % xray penetration depth [nm]
opts.xi0xray = 1e8;     % xray penetration depth [nm]
opts.kT =1e-04*70;        % kT 
opts.kT =1e-04*20;        % kT 
opts.kT =1e-04*0;   

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
X00 = ones(Nx, Ny);
% for ii=1:opts.npoints(2)
%     for jj=1:opts.npoints(3)
%         X00(ii,jj) = (-1)^(ii+jj);
%     end
% end
X0 = [  X00(:)', zeros(1,prod(opts.npoints(2:end))),...
        X00(:)', zeros(1,prod(opts.npoints(2:end)))];
    
opts.r = @(xs,t) - 1;                           % constant quadratic coefficient
opts.sig = @(xs,t) sqrt(2*opts.gam*opts.kT);    % constant temperature

    
%% run over the ensemble of initial conditions
Nensembles = 1;
% setup k-space
nq = 200;
Qs = linspace(6, 10, nq)'*[1,1]; 
% Qs2 = linspace(6, 10, nq)'*[1,-1];
Qs2 = linspace(6, 10, nq)'*[1,2]; 
Qs = [Qs; Qs2]';
Fall = zeros(size(Qs,2), opts.npoints(1), Nensembles);

% for a 2D grid in Q-space use this:
%{
Nensembles = 1;
[qqx,qqy] = ndgrid(linspace(6, 10, 30));
Qs = cat(3, qqx, qqy);
Qs = permute(Qs, [3, 1,2]);
Qs = reshape(Qs, [], prod(size(qqx)));
Fall = zeros(size(Qs,2), opts.npoints(1), Nensembles);
%}

%% this option sets the equilibrium at t < 0
opts0 = opts;
opts0.r = @(xs,t) - 1;                      % constant quadratic coefficient
opts0.sig = @(xs,t) sqrt(2*opts0.gam*opts0.kT);     % constant temperature

prefact_r = 100;
prefact_sig = .1;

% this option set the potential (r param.) and temperature (sig) for t>0
% rT      = @(t) exp(-(t)/opts.tauf).*nheaviside(t);
rT      = @(t) (1/opts.tauf).*exp(-(t)/opts.tauf).*nheaviside(t);
% rT      = @(t) exp(-z/xi0)*exp(-(t)/opts.tauf).*nheaviside(t).*errf(0);

% opts.r = @(xs,t) repmat(prefact_r*rT(t-1).*exp(-xs{1}/opts.xi0) + opts0.r(0,0), [length(xs{2}),1] );
% 
% opts.r = @(xs,t) repmat(prefact_r*rT(t-1).*exp(-(xs{1}/(2*opts.xi0)).^2) + opts0.r(0,0), [length(xs{2}),1] );

% play with definition of gaussian
% opts.r = @(xs,t) repmat(prefact_r*rT(t-1).*exp(-(xs{1}/(2*opts.xi0)).^2) + opts0.r(0,0), [length(xs{2}),1] );
% opts.sig = @(xs,t) repmat(sqrt(2*opts.gam*(exp(-xs{1}/opts.xi0).*prefact_sig*rT(t-1)...
%                     +opts.kT)), [length(xs{2}),1] );

opts.r = @(xs,t) repmat(prefact_r*rT(t-1).*((1/(2*opts.xi0))).*exp(-(xs{1}/(2*opts.xi0)).^2) + opts0.r(0,0), [length(xs{2}),1] );
opts.sig = @(xs,t) repmat(sqrt(2*opts.gam.*(1/(opts.xi0)).*(exp(-xs{1}/opts.xi0).*prefact_sig*rT(t-1)...
                    +opts.kT)), [length(xs{2}),1] );         
       
%%
for niter = 1:Nensembles
niter
%% thermalize
[F,G] = make_update_function_VO2_tdeptemp(opts0);
dd = integrate_with_XSPDE(X0, F, G, opts0);
X = permute(dd, [2, 3, 1, 4]);
X0 = X(:,:,end,:);
X0 = X0(:)';

%% run the dynamics with time dependent parameters
[F,G] = make_update_function_VO2_tdeptemp(opts);
dd = integrate_with_XSPDE(X0, F, G, opts);
X = permute(dd, [2, 3, 1, 4]);

%% build 2D square crystal model and calculate diffraction pattern
a0 = 5.78;
at = 1*[1 0; 0 1];

% generate the locations of the dimerized lattice:
rs = zeros(Nx,Ny,2);
atten = zeros(Nx,Ny);
staggered_phase = zeros(Nx,Ny);
for i=1:Nx
    for j=1:Ny
        rs(i,j,:) = (i-1)*at(:,1)+(j-1)*at(:,2);
        atten(i,j) = exp(-i*(0.1*a0)/opts.xi0xray);
        staggered_phase(i,j) = (-1)^(i+j);
    end
end

% the positions are slightly perturbed by the computed displacements:
X2=reshape(staggered_phase.*X, Nx*Ny, Nt, 4);

Xcr = zeros(Nx*Ny, Nt, 2);
rs = reshape(rs, Nx*Ny, 2);
delta = psi0/a0;

for j=1:Nt
    Xcr(:,j,:) = delta*[X2(:, j, 1), X2(:, j, 3)] + rs;
end

%save the diffraction pattern
% Fall(:,:,niter) = calc_Fhkl_from_trajectory_1d(Qs, Xcr);
Fhkl = zeros(size(Qs,2), Nt);
qim = 2*pi*1i*Qs;

tic
for j=1:Nt
    r1 = squeeze(Xcr(:, j, :));
    phases = r1(:,:)*qim;
    Fall(:,j,niter) = sum(atten(:).*exp(phases), 1);
%     Fall(:,j,niter) = sum(1.*exp(phases), 1);
end
toc

end
%%
% F2all = abs(Fall).^2;
% F2mean = mean(F2all,3);
% F2max = (F2mean-mean(F2mean(:,t<1),2));
% F2max = (F2mean./mean(F2mean(:,t<.5),2));


F2all = abs(Fall).^2;
F2mean = mean(F2all,3);
F2max = (F2mean-mean(F2mean(:,t<.5),2));
F2max = F2max./max(F2mean,[],2);


% F2max = F2mean./mean(F2mean(:,t<.5),2); % this is correct, but almost fully kills the intensity
%% Plot the averaged structure factors
%
figure(10); clf
idx = [1, 25, 50, 75, 100,150,200];
idxQs = 1;
offsets = [0, 1, 2.7, 3 ,4.8,5.4, 6]';
plot(t-1, F2max(idx,:) + offsets)
hold on
plot(t-1, 3*mean(F2max(idx(1):idx(end),:),1)+6.5,'b')
% xlim([-.2, 3])
legend([num2str([num2str(Qs(idxQs,idx)')]); "m"])
xlabel('time')
ylabel('Relative Intensity (arb)')
xlim([-0.1, 4.5])
grid on 
grid minor

%% Plot the averaged structure factors
%
figure(13); clf
idx = [1, 26, 50, 75, 100,150,200]+200;
offsets = [0, 1, 2.7, 3 ,4.8,5.4, 6]';
plot(t-1, F2max(idx,:) + offsets)
hold on
plot(t-1, 3*mean(F2max(idx(1):idx(end),:),1)+6.5,'b')
% xlim([-.2, 3])
legend([num2str([num2str(Qs(1,idx)')]); "m"])
xlabel('time')
ylabel('Relative Intensity (arb)')
xlim([-0.1, 4.5])
grid on 
grid minor

%% Plot the averaged structure factors some
%
figure(10); clf
idxQs = 1;
idx = [1, 50, 100,150,200];
offsets = [0, 1 ,2,3,4]';
idx = [1];
offsets = [0]';
plot(t, F2max(idx,:) + offsets)
% plot(t, F2max(idx,:)/(F2max(idx,10)) + offsets)
hold on
F2difuse = mean(F2mean(F2mean(:,12) < 1e2,:),1);
plot(t, F2difuse./F2difuse(end))
% plot(t, 3*mean(F2max(idx(1):idx(end),:),1)+2,'b')
% xlim([-.2, 3])
legend([num2str([num2str(Qs(1,idx)')]); "m"])
xlabel('time')
ylabel('Relative Intensity (arb)')
xlim([0.9, 5.5])
xlim([-0.1, 5.5])
grid on 
grid minor
%% For a 2D section of Q-space use this 
% 
% figure(1); 
% figpos = get(figure(1),'position');
% close(figure(1));
% figure(1);
% set(figure(1), 'position', figpos)
% 
% F2mean = reshape(F2mean, [size(qqx), Nt]);
% plotdatacube(t, (F2mean),[],[0, 10],[0,0]); 
% % plotdatacube(t, log10(F2mean),[],[0, 10],[0,0]); 
% % F2max = reshape(F2max, [size(qqx), Nt]);
% % plotdatacube(t, F2max,[],[0, 10],[0,0]); 
% caxis auto


%%



%% save simu output
% outputfilename=sprintf('output/simT=%1g_tdep.mat',opts.kT);
% save(outputfilename, 'X','opts','F2all')

save_folder = '../VO2_Langevin_sims';
% myfile = sprintf('simT=%1g_prefr=%1g_tdep.mat',opts.kT,prefact_r)
myfile = sprintf('simT=%1g_gam=%1g_prefr=%1g_prefsig=%1g_xi0=%1g_tdep.mat',opts.kT,opts.gam,prefact_r,prefact_sig,opts.xi0)
save(strcat(save_folder,myfile), 'F2all','X00','X','Nx','Ny','Nt', 'opts0','opts','t','Qs')


videofile=sprintf('movies/movT=%1.g_tdep',opts.kT)

%% Plot the diffuse only
% figure(9); hold all
% plot(t, mean(F2max,1))
% xlabel('time')
% ylabel('Relative Intensity (arb)')
% xlim([0.8, 3])

% 
F2difuse = mean(F2mean(F2mean(:,12) < 1e2,:),1);
figure(8); hold all
plot(t, F2difuse./F2difuse(end))
xlabel('time')
ylabel('Relative Intensity (arb)')
xlim([0.8, 3])

