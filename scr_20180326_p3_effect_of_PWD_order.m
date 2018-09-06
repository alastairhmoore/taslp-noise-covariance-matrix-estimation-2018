exp_id = '20180326_x3';

f_target_seq = [100 220 470 1000 2200 4700] ;
pwd_order_seq{1} = [1 2 3 4 14];
pwd_order_seq{2} = [1 2 3 4 14];
pwd_order_seq{3} = [1 2 3 4 14];
pwd_order_seq{4} = [1 2 3 4 5 6 7 8 14];
pwd_order_seq{5} = [1 2 3 4 5 6 7 8 10 12 14];
pwd_order_seq{6} = [1 2 3 4 5 6 7 8 10 12 14];
% f_target_seq = [1000 2200 4700] ;
% pwd_order_seq{1} = [1 2 3 4 5];
% pwd_order_seq{2} = [1 2 3 4 5 6 7 8 10 12 14];
% pwd_order_seq{3} = [1 2 3 4 5 6 7 8 10 12 14];

sphHarmOrdDiffuseGt = 2;
sphHarmOrdDiffuseEst = 2;


nFreq = length(f_target_seq);
collated_mse = cell(nFreq,1);

for ifreqseq = nFreq:-1:1
    nPwdOrd = length(pwd_order_seq{ifreqseq});
    collated_mse{ifreqseq} = zeros(nPwdOrd,1);
    f_target = f_target_seq(ifreqseq);
    for iordseq = 1:nPwdOrd
        sphHarmOrdPWDModel = pwd_order_seq{ifreqseq}(iordseq)

        scr_id = sprintf('%s_gt_o%d_%d_Hz_pwd_ord_%d',...
            exp_id,sphHarmOrdDiffuseGt,...
            f_target,sphHarmOrdPWDModel);
        
        
        load(sprintf('%s/dat_%s_summary.mat',data_dir,scr_id),...
    'nr_mat','mse_mat',...
    'legstr');

       collated_mse{ifreqseq}(iordseq) = mean(mse_mat(:,1,2),1); 
    end
end

%%
nColors = 3;
v_colors = get_voicebox_thermliny_linecolors(nColors);
styles = {'-','--'}

fig_width = column_width
fig_height = 3/4*fig_width


fh = figure;
setFigureSize(fh,[fig_width fig_height],'centimeters');
lh = [];
legstr = {};
legend_order = [nFreq:-1:1];
for ifreqseq = 1:nFreq%:-1:1
    lh(ifreqseq) = plot(pwd_order_seq{ifreqseq},collated_mse{ifreqseq},...
        'Marker','s',...
        'linewidth',1.5);
    legstr{ifreqseq} = strcat(num2str(f_target_seq(ifreqseq)),' Hz');
    hold all;
end

%cycle through colors and styles
for ii=1:nColors
    set(lh(ii:3:end),'Color',v_colors(ii,:),'MarkerFaceColor',v_colors(ii,:));
end
set(lh(1:nColors),'LineStyle',styles{1});
set(lh(nColors+[1:nColors]),'LineStyle',styles{2},'MarkerFaceColor',[1 1 1]);
% 
% xlabel('PWD order');
% ylabel('Norm of covariance matrix estimation error [dB]')
xlabel('$N_{H}$','Interpreter','latex')
ylabel('$\bar{\mathcal{E}}$ \textsf{[dB]}','Interpreter','latex')


legh = legend(lh(legend_order),legstr(legend_order));
title(legh,'Frequency')

alltext = findall(fh,'-property','FontSize');%findall(gcf,'Type','text');
set(alltext,'FontSize',0.8*font_size);
set(gca,'linewidth',1.1)

ahm_print_to_pdf(fh,fullfile(graphics_dir,...
    sprintf('fig_%s_effect_of_PWD_ord_per_freq.pdf',exp_id)))
