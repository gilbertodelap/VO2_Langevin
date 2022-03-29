%% save data as .mat file
% save_folder = './saved_sims/';
% myfile = 'test1.mat';
% save(strcat(save_folder,myfile), 'F2all','X00','X','Nx','Ny','Nt', 'opts0','opts')
cd ~/Documents/DATA/2022/VO2_Langevin/
addpath(genpath('./xspde_matlab/XSPDE_code/'));
addpath(genpath('~/Documents/xdscode-master/'))

%%
save_folder = '../VO2_Langevin_sims/';
% myfile = 'simT=0.002_prefr=3.2_xi0=0.1_tdep.mat';
load(strcat(save_folder,myfile))
%% plot lattice distortion field
figure(1); close(figure(1)); figure(1)
figpos = get(figure(1),'position');
close(figure(1));
figure(1);
set(figure(1), 'position', figpos)
X2=reshape(X00.*X, Nx,Ny, Nt, 4);

plotdatacube(t, angle(X2(:,:,:, 1)+1i*X2(:,:,:, 3)),[],[-pi,pi],[0,0]); colormap hsv
% plotdatacube(t, (X(:,:,:, 1)),[],[0,1.2],[0,0]); 
axis image
% caxis('auto')
% caxis([-3 3])
% ylim([0, 100])

tfsavedistortion = exist('videofile','var');


%%
% plot distribution in REAL space
tfsavemovie = exist('videofile','var');
tfsavemovie = 0;
if tfsavemovie
    vid = VideoWriter(videofile,'MPEG-4');
    vid.Quality = 50;
    open(vid);
end
% [xx,yy]=ndgrid(1:Nx,1:Ny);
% ss=xx.^2+yy.^2;
figure(11); clf; 
axes
set(gca, 'xlimmode','manual',...
'ylimmode','manual',...
'zlimmode','manual',...
'climmode','manual',...
'alimmode','manual');

X2=reshape(X00.*X, Nx*Ny, Nt, 4);
xmin=min(min(X2(:,:,1)));
xmax=max(max(X2(:,:,1)));
ymin=min(min(X2(:,:,3)));
ymax=max(max(X2(:,:,3)));
% p1=scatter(X2(:,1,1), X2(:,1,3),15,'fill','markerfacecolor',[0.5, 0.2 0.9]);
p1=scatter(X2(:,1,1), X2(:,1,3),15,'fill','markerfacecolor',[0.5, 0.2 0.9]);
% p1=scatter(X2(:,1,1), X2(:,1,3),15,xx(:),'fill');
hold on
title_ = title(sprintf('t = %2.2f ps', t(1)));

xlabel('\Deltax_1')
ylabel('\Deltax_2')
axis image 
xlim([xmin xmax]); ylim([ymin ymax]); 
% axis(ax)
for j=1:1:Nt
%%
    set(p1,'xdata', X2(:,j,1), 'ydata', X2(:,j,3))
    title_.String = sprintf('t = %2.2f ps', t(j));
    drawnow; pause(0.02);
    if tfsavemovie
        pause(0.05)
        frame = getframe(gcf);
        writeVideo(vid, frame); 
    end
end
if tfsavemovie
    close(vid);
end



%% plot the crystal dynamics
% % filename = 'crystal2D_64x64_16bit.gif';
% figure(19); clf; 
% hold on
% set(gca, 'xlimmode','manual',...
% 'ylimmode','manual',...
% 'zlimmode','manual',...
% 'climmode','manual',...
% 'alimmode','manual');
% 
% % define the range of coordinates to plot
% xmin=min(min(X2(:,:,1)))-1;
% xmax=max(max(X2(:,:,1)))+22*1;
% ymin=min(min(X2(:,:,3)))-1;
% ymax=max(max(X2(:,:,3)))+22*1;
% p2=scatter(Xcr(:,1,1), Xcr(:,1,2),150,'fill','markerfacecolor',[0.6, 0.6 0.6]);
% p1=plot(Xcr(:,2,1), Xcr(:,2,2),'.', 'color',[0.5, 0.2 0.9],'markersize',40);
% 
% xlabel('x1')
% ylabel('x2')
% xlim([xmin xmax/2]); ylim([ymin ymax/2]); 
% axis image 
% axis([-.2, 22.1*1, -0.2, 22.1*1])
% 
% for j=1:2:Nt %+1
% %%
%     %%
% %     j = 1;
% %     r1 = [phs(:).*xx(:, j, 1), phs(:).*xx(:, j, 3)] + reshape(rs, N,2);
% %     set(p1,'xdata', r1(:,1), 'ydata', r1(:,2))
%     set(p1,'xdata', Xcr(:,j,1), 'ydata', Xcr(:,j,2))
%     drawnow; 
% %     
% %     % Capture the plot as an image 
% %     frame = getframe(2); 
% %     im = frame2im(frame); 
% %     [imind,cm] = rgb2ind(im,16); 
% % %    % Write to the GIF File 
% %     if j == 2 
% %       imwrite(imind,cm,filename,'gif', 'Loopcount',inf,'delaytime',0.002); 
% %     else 
% %       imwrite(imind,cm,filename,'gif','WriteMode','append','delaytime',0.002); 
% %     end 
% 
% end
% % fprintf('Done.\n')


%% load the relevant runs
% save_folder = './saved_sims/';
save_folder = "../VO2_Langevin_sims/";

% initialize what to plot
simT_plot = 'all'; gam_plot = 'all'; prefr_plot = 'all'; 
prefsig_plot = 'all'; prefcohl2_plot = 'all'; xi0_plot = 'all';

% simT_plot = [0.002];
simT_plot = [0.011116];
% simT_plot = [0.033348];

% simT_plot = [0.002];

% xi0_plot = [.7,.9]*360;
% xi0_plot = 252;
xi0_plot = [.1]*360;

prefcohl2_plot = [1];

prefr_plot = [200];

prefsig_plot = [0.1];
% runnames = everything in folder
runname_list = string({dir(strcat(save_folder,'*.mat')).name})';

%
runnames = {};
list_iter = 1;
for myindex = 1:size(runname_list)
    flag = 1; 
    myname = runname_list(myindex);
    myname_split = split(myname,'_');
    
    if ~strcmp(simT_plot,"all")
        valcheck = split(myname,'simT=');
        valcheck = split(valcheck(2),'_');
        valcheck = str2num(valcheck(1));
        if ~ismember(valcheck,simT_plot) % check if it is not there
            flag = 0;
        end
    end
    
    if ~strcmp(gam_plot,"all")
        valcheck = split(myname,'gam=');
        valcheck = split(valcheck(2),'_');
        valcheck = str2num(valcheck(1));
        if ~ismember(valcheck,gam_plot) % check if it is not there
            flag = 0;
        end
    end

    if ~strcmp(prefr_plot,"all")
        valcheck = split(myname,'prefr=');
        valcheck = split(valcheck(2),'_');
        valcheck = str2num(valcheck(1));
        if ~ismember(valcheck,prefr_plot) % check if it is not there
            flag = 0;
        end
    end

    if ~strcmp(prefsig_plot,"all")
        valcheck = split(myname,'prefsig=');
        valcheck = split(valcheck(2),'_');
        valcheck = str2num(valcheck(1));
        if ~ismember(valcheck,prefsig_plot) % check if it is not there
            flag = 0;
        end
    end

    if ~strcmp(prefcohl2_plot,"all")
        valcheck = split(myname,'prefcohl2=');
        valcheck = split(valcheck(2),'_');
        valcheck = str2num(valcheck(1));
        if ~ismember(valcheck,prefcohl2_plot) % check if it is not there
            flag = 0;
        end
    end

    if ~strcmp(xi0_plot,"all")
        valcheck = split(myname,'xi0=');
        valcheck = split(valcheck(2),'_');
        valcheck = str2num(valcheck(1));
        if ~ismember(valcheck,xi0_plot) % check if it is not there
            flag = 0;
        end
    end
%     flag
%     myname
    if flag == 1
        runnames{list_iter} = myname;
        list_iter = list_iter + 1; 
    end
end
runnames = string(runnames)';

myfile = runnames(1)
load(strcat(save_folder,myfile))
%% plot same peak diffeerent runs
% on demand
% save_folder = "/media/gilberto/data2/DATA/2022/VO2_sims/saved_sims/";
% runnames = [
%     "simT=0_gam=0.03_prefr=100_prefsig=0_xi0=36_tdep.mat"
%     "simT=0_gam=0.03_prefr=200_prefsig=0_xi0=36_tdep.mat"
% ];.

figure(20); clf
legend_info = cell(size(runnames));
for run = 1:size(runnames)
    
    myfile = runnames(run)
    load(strcat(save_folder,myfile))
    
    F2mean = mean(F2all,3);
    F2max = F2mean./mean(F2mean(:,t<.5),2);
    
    % plot peak at indicated Q
    idxQs = 1;
    whichQs = [.5,.5]';
    idx1 = find(Qs'==whichQs',1);
    offsets = [0]';
    plot(t, F2max(idx1,:) + offsets, 'LineWidth',4)

    hold on
    whichQs = [1,1]';
    idx2 = find(Qs'==whichQs',1);
    offsets = [0]';
    plot(t, F2max(idx2,:) + offsets, 'LineWidth',4)
%     plot(t, F2max(idx,:)/(F2max(idx,10)) + offsets)

    hold on
    F2difuse = mean(F2mean(F2mean(:,12) < .5e2,:),1);
    plot(t, F2difuse./F2difuse(end)+1,'LineWidth',4)
%     plot(t, 3*mean(F2max(idx(1):idx(end),:),1)+2,'b')
    % xlim([-.2, 3])
    legend_info{run} = [strcat('Q= ', num2str(Qs(1,idx1)'),' xi0= ',  num2str(opts.xi0))];
    legend_info{run} = [strcat('Q= ', num2str(Qs(1,idx2)'),' xi0= ',  num2str(opts.xi0))];


    hold on
end

legend(legend_info)

xlabel('time')
ylabel('Relative Intensity (arb)')
xlim([0.9, 5.5])
xlim([-0.1, 5.5])
ylim([0,2.1])
grid on 
grid minor
set(gca,'FontSize',15)



%% plot different peak same run
% on demand
% save_folder = "/media/gilberto/data2/DATA/2022/VO2_sims/saved_sims/";
% runnames = [
%     "simT=0_gam=0.03_prefr=200_prefsig=0_xi0=36_tdep.mat"
% ];

figure(20); clf
legend_info = cell(size(runnames));
for run = 1:size(runnames)
    
    myfile = runnames(run)
    load(strcat(save_folder,myfile))
    
    F2mean = mean(F2all,3);
    F2max = F2mean./mean(F2mean(:,t<.5),2);

    idxQs = 1;
    qcellsiz = 19;
    idx = [12,12+qcellsiz] %+200 % for other direction of Qs;
    offsets = [0,0]';

%     idx = [2,4,6,8,10,12,14,16,18]+19 %+200 % for other direction of Qs;
%     offsets = [0,0,0,0,0,0,0,0,0]';
%     idx = [2,4,6,8]+19 %+200 % for other direction of Qs;
%     offsets = [0,0,0,0]';

%     idx = [2,4,6,8,10]+19+8 %+200 % for other direction of Qs;
%     offsets = [0,0,0,0,0]';

    plot(t, F2max(idx,:) + offsets)
    % plot(t, F2max(idx,:)/(F2max(idx,10)) + offsets)
    
%     F2difuse = mean(F2mean(F2mean(:,12) < 1e2,:),1);
    % plot(t, F2difuse./F2difuse(end))
    % plot(t, 3*mean(F2max(idx(1):idx(end),:),1)+2,'b')
    % xlim([-.2, 3])
    legend_info{run} = [strcat('Q= ', num2str(Qs(1,idx)'),' xi0= ',  num2str(opts.xi0))];

    hold on
end

legend(legend_info)

xlabel('time')
ylabel('Relative Intensity (arb)')
xlim([0.9, 5.5])
xlim([-0.1, 5.5])
grid on 
grid minor

%% plot intensity in reciprocal space
% figure(20); clf
% % ii=4000;
% idx = [ 600, 800,  1000, 1500, 2000,4000];
% % idx = [1, 600, 1500, 3000, 4000];
% plot(log(F2mean(:,idx)))
% legend(num2str(t(idx)',2))
% xlabel('q')
% ylabel('Log(I)')
% title(['Te = ', num2str(opts.Te)])
% 
% 
% 
% %% Plot the diffuse from saved sims
% 
% figure(8); clf
% 
% % thefiles = dir('output/*.mat');
% thefiles = dir('output/*tdep.mat');
% kTs = [];
% hold all
% 
% for ii=1:length(thefiles)
%     load([thefiles(ii).folder,'/', thefiles(ii).name],'F2all','opts');
%     t = 0:opts.ranges(1)/opts.npoints(1):opts.ranges(1);  t=t(1:end-1);
%     kTs = [kTs, opts.kT];
%     F2mean = mean(F2all,3);
%     F2difuse = mean(F2mean(F2mean(:,12) < 1e2,:),1);
%     plot(t, F2difuse./F2difuse(end))
% end
% xlabel('time')
% ylabel('Relative Intensity (arb)')
% xlim([0.8, 3])
% legend(num2str(1e3*kTs'))
% 
% 
% %% Plot the Bragg peaks from saved sims
% figure(7); clf
% 
% hold all
% for ii=1:length(thefiles)
%     load([thefiles(ii).folder,'/', thefiles(ii).name],'F2all','opts');
%     F2mean = mean(F2all,3);
%     F2Bragg = F2mean(1,:);
%     plot(t, F2Bragg./F2Bragg(5))
% end
% 
% xlabel('time')
% ylabel('Relative Intensity (arb)')
% xlim([0.8, 3])
% legend(num2str(1e3*kTs'))
