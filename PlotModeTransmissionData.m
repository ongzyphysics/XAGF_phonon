function PlotModeTransmissionData(DataFilesDir,wvec,Left,Right,LeftPhonon,RightPhonon,PhononData)
% Plots the phonon dispersion of the left and right leads along with the transmission coefficient 
% of the individual modes in the left and right leads

symb_size = 5;
symb_type = '.';
close all;
nwmax = length(wvec);

filename_t_R = 'Plot_Right_Transmission';
filename_t_L = 'Plot_Left_Transmission';
filename_r_L = 'Plot_Left_Reflection';
filename_xi = 'Plot_Transmission';
format_t_R = 'png';
format_t_L = 'png';
format_r_L = 'png';
format_xi = 'png';

% ----------------------------------------------------------------
% Plot right lead phonon dispersion (transmission)
% ----------------------------------------------------------------
set(figure(1),'visible','off');
set(gcf,'position',[400 150 800 800])
colormap jet;
caxis([0 1]);

subplot(1,2,1);
[kphon, wphon, vphon] = GetPhononDispersion(Right);
plot(kphon/(2*pi)*Right.a_long,real(wphon),'k-','linewidth',0.5);
ax1 = gca; 
set(ax1,'fontsize',16);

subplot(1,2,2);
% tphon = [PhononData(:).Xi_mode]./[RightPhonon(:).Xi_mode]; tphon(eq(round([RightPhonon(:).Xi_mode]),0)) = 0;
tphon = [PhononData(:).Xi_mode];
area([0; wvec],[1 real(tphon)],'facecolor',[1 1 1]*0.95);
hold on;
plot([0; wvec],[1 real(tphon)],'k-','linewidth',1);
ax2 = gca; 
set(ax2,'fontsize',16);


for nw = 1:1:nwmax
    w = wvec(nw);
    q_L_temp = [LeftPhonon(nw).VecQ_minus_adv];
    q_R_temp = [RightPhonon(nw).VecQ_plus];
    t_L_temp = [PhononData(nw).tt_L];
    t_R_temp = [PhononData(nw).tt_R];
    % q_L_temp = -[LeftPhonon(nw).QR_adv];
    % q_R_temp = [RightPhonon(nw).QR];
    % t_L_temp = [PhononData(nw).tt_L];
    % t_R_temp = [PhononData(nw).tt_R];

    t_L_temp = t_L_temp(gt(abs(q_L_temp),0));
    t_R_temp = t_R_temp(gt(abs(q_R_temp),0));    
    q_L_temp = q_L_temp(gt(abs(q_L_temp),0));
    q_R_temp = q_R_temp(gt(abs(q_R_temp),0));
             
    if gt(numel(q_R_temp),0)
        w_temp = w*ones(size(q_R_temp));
        % s_temp = [];
        s_temp = symb_size;

        subplot(1,2,1);
        hold on;
        scatter(q_R_temp/(2*pi),w_temp,s_temp,t_R_temp,symb_type,'linewidth',2);        

        subplot(1,2,2);
        hold on;
        scatter(w_temp,t_R_temp,s_temp,t_R_temp,symb_type,'linewidth',2);                
    end
end

wlim = [0 max(wvec)];

subplot(1,2,1);
box on; caxis([0 1]);
xlim([-1 1]*0.5); ylim(wlim);
ax1 = gca;
 
subplot(1,2,2);
xlim(wlim);
box on; caxis([0 1]);
view(90,-90);
set(colorbar('north'),'fontsize',16);
set(1,'position',[400 150 800 800]);

cd(DataFilesDir);
saveas(1,filename_t_R,format_t_R);
cd('..');

fprintf(1,'\t  <%s> \n', [filename_t_R '.' format_t_R]);


%%

% ----------------------------------------------------------------
% Plot left lead phonon dispersion (transmission)
% ----------------------------------------------------------------
set(figure(2),'visible','off');
set(gcf,'position',[400 150 800 800])
colormap jet;
caxis([0 1]);

subplot(1,2,1);
[kphon, wphon, vphon] = GetPhononDispersion(Left);
plot(kphon/(2*pi)*Left.a_long,real(wphon),'k-','linewidth',0.5);
ax1 = gca; 
set(ax1,'fontsize',16);

subplot(1,2,2);
% tphon = [PhononData(:).Xi_mode]./[LeftPhonon(:).Xi_mode]; tphon(eq(round([LeftPhonon(:).Xi_mode]),0)) = 0;
tphon = [PhononData(:).Xi_mode];
area([0; wvec],[1 real(tphon)],'facecolor',[1 1 1]*0.95);
hold on;
plot([0; wvec],[1 real(tphon)],'k-','linewidth',1);
ax2 = gca; 
set(ax2,'fontsize',16);


for nw = 1:1:nwmax
    w = wvec(nw);
    q_L_temp = [LeftPhonon(nw).VecQ_minus_adv];
    q_R_temp = [RightPhonon(nw).VecQ_plus];
    t_L_temp = [PhononData(nw).tt_L];
    t_R_temp = [PhononData(nw).tt_R];
    % q_L_temp = -[LeftPhonon(nw).QR_adv];
    % q_R_temp = [RightPhonon(nw).QR];
    % t_L_temp = [PhononData(nw).tt_L];
    % t_R_temp = [PhononData(nw).tt_R];

    t_L_temp = t_L_temp(gt(abs(q_L_temp),0));
    t_R_temp = t_R_temp(gt(abs(q_R_temp),0));    
    q_L_temp = q_L_temp(gt(abs(q_L_temp),0));
    q_R_temp = q_R_temp(gt(abs(q_R_temp),0));
             
    if gt(numel(q_L_temp),0)
        w_temp = w*ones(size(q_L_temp));
        % s_temp = [];
        s_temp = symb_size;

        subplot(1,2,1);
        hold on;
        scatter(q_L_temp/(2*pi),w_temp,s_temp,t_L_temp,symb_type,'linewidth',2);        

        subplot(1,2,2);
        hold on;
        scatter(w_temp,t_L_temp,s_temp,t_L_temp,symb_type,'linewidth',2);                
    end
end

wlim = [0 max(wvec)];

subplot(1,2,1);
box on; caxis([0 1]);
xlim([-1 1]*0.5); ylim(wlim);
ax1 = gca;
 
subplot(1,2,2);
xlim(wlim);
box on; caxis([0 1]);
view(90,-90);
set(colorbar('north'),'fontsize',16);
set(2,'position',[400 150 800 800]);

cd(DataFilesDir);
saveas(2,filename_t_L,format_t_L);
cd('..');

fprintf(1,'\t  <%s> \n', [filename_t_L '.' format_t_L]);

%%

% ----------------------------------------------------------------
% Plot left lead phonon dispersion (reflection)
% ----------------------------------------------------------------
set(figure(20),'visible','off');
set(gcf,'position',[400 150 800 800])
colormap jet;
caxis([0 1]);

subplot(1,2,1);
[kphon, wphon, vphon] = GetPhononDispersion(Left);
plot(kphon/(2*pi)*Left.a_long,real(wphon),'k-','linewidth',0.5);
ax1 = gca; 
set(ax1,'fontsize',16);

subplot(1,2,2);
% tphon = [PhononData(:).Xi_mode]./[LeftPhonon(:).Xi_mode]; tphon(eq(round([LeftPhonon(:).Xi_mode]),0)) = 0;
tphon = [PhononData(:).AntiXi_mode];
area([0; wvec],[1 real(tphon)],'facecolor',[1 1 1]*0.95);
hold on;
plot([0; wvec],[1 real(tphon)],'k-','linewidth',1);
ax2 = gca; 
set(ax2,'fontsize',16);


for nw = 1:1:nwmax
    w = wvec(nw);
    q_L_temp = [LeftPhonon(nw).VecQ_minus];
    t_L_temp = [PhononData(nw).rr_L];
    t_L_temp = t_L_temp(gt(abs(q_L_temp),0));
    q_L_temp = q_L_temp(gt(abs(q_L_temp),0));

             
    if gt(numel(q_L_temp),0)
        w_temp = w*ones(size(q_L_temp));
        % s_temp = [];
        s_temp = symb_size;

        subplot(1,2,1);
        hold on;
        scatter(q_L_temp/(2*pi),w_temp,s_temp,t_L_temp,symb_type,'linewidth',2);        

        subplot(1,2,2);

        hold on;
        scatter(w_temp,t_L_temp,s_temp,t_L_temp,symb_type,'linewidth',2);                
    end
end

wlim = [0 max(wvec)];

subplot(1,2,1);
box on; caxis([0 1]);
xlim([-1 1]*0.5); ylim(wlim);
ax1 = gca;
 
subplot(1,2,2);
xlim(wlim);
box on; caxis([0 1]);
view(90,-90);
set(colorbar('north'),'fontsize',16);
set(2,'position',[400 150 800 800]);

cd(DataFilesDir);
saveas(20,filename_r_L,format_r_L);
cd('..');

fprintf(1,'\t  <%s> \n', [filename_r_L '.' format_r_L]);

%%

% ----------------------------------------------------------------
% Plot total transmission
% ----------------------------------------------------------------
set(figure(3),'visible','off');
set(gcf,'position',[400 150 800 600])

tphon = [PhononData(:).Xi_negf];
tphon2 = [PhononData(:).Xi_mode];
area([0; wvec],[1 real(tphon)],'facecolor',[1 1 1]*0.95);
hold on;
plot([0; wvec],[1 real(tphon)],'k-','linewidth',1);
plot([0; wvec],[1 real(tphon2)],'ko','linewidth',1);
plot([0; wvec],[1 real([LeftPhonon(:).Xi_mode])],'r--','linewidth',2);
plot([0; wvec],[1 real([RightPhonon(:).Xi_mode])],'b--','linewidth',2);
ax3 = gca; 
set(ax3,'fontsize',16);

wlim = [0 max(wvec)];

xlim(wlim);
box on; caxis([0 1]);

cd(DataFilesDir);
saveas(3,filename_xi,format_xi);
cd('..');

fprintf(1,'\t  <%s> \n', [filename_xi '.' format_xi]);


%%

close all;
