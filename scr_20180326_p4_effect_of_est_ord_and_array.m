exp_id = '20180326_x4';
sphHarmOrdPWDModel = 15;
sphHarmOrdDiffuseGt = [];                       % variable under test
sphHarmOrdDiffuseMax = 4;
pose_measurement_noise_var_rad = deg2rad(1);
f_target = [];                                  % variable under test
ema = [];                                       % variable under test

methods_to_run = [ 1 ... % gt
    11 ... % EWLS ord 1 0.99999
    12 ... % EWLS ord 2 0.99999
    13 ... % EWLS ord 3 0.99999
    14 ... % EWLS ord 4 0.99999
    ];

f_target_seq = [2200];
sph_ord_gt_seq = [1 2 3 4];
mic_arrays = {@RigidSphere_20171120_01_horiz_up_20deg_circle4,...
    @RigidSphere_20171102_02_horiz_up_20deg_circle16};
mic_lab = {'circ4','circ16'};
mic_lab_plot = {'$Q=4$ microphones','$Q=16$ microphones'};
ylims = {[0 1.2],[0 12]}

nFreq = length(f_target_seq);
comb_mse_mat = cell(length(mic_lab),length(f_target_seq));
comb_nr_mat = cell(length(mic_lab),length(f_target_seq));

for imic = 1:length(mic_lab)
    for iftarget = 1:length(f_target_seq)
        f_target = f_target_seq(iftarget);
        
        comb_mse_mat{imic,iftarget}=[];
        comb_nr_mat{imic,iftarget}=[];
        
        for isphord = 1:length(sph_ord_gt_seq)
            
            sphHarmOrdDiffuseGt = sph_ord_gt_seq(isphord);
            
            scr_id = sprintf('%s_gt_o%d_%d_Hz_%s_array',...
                exp_id,sphHarmOrdDiffuseGt,...
                f_target,mic_lab{imic});
            
            load(sprintf('%s/dat_%s_summary.mat',data_dir,scr_id))
            
            comb_mse_mat{imic,iftarget}(:,:,isphord) = squeeze(mse_mat);
            comb_nr_mat{imic,iftarget}(:,:,isphord) = squeeze(nr_mat);
        end
    end
end

%% make a single legend
legstr = {'White','Sph iso',...
    'RS: 5e-2','RS: 1e-2','RS: 5e-3','RS: 1e-3',...
    'EWLS: 1','EWLS: 2',...
    'EWLS: 3','EWLS: 4'};

fig_width = column_width;
fig_height = 1/8*fig_width;
highlightwidth = 1.6;
figure;
fh = gcf;
setFigureSize(fh,[fig_width fig_height],'centimeters');

nCol = 5;
nRow = 2;
ax_width = 1/nCol;
ax_height = 1;
set(gcf,'NextPlot','add')
for iax = 1:nCol
    lax(iax) = axes('position',[(iax-1)*ax_width 0 ax_width ax_height]);
    h{iax} = bar(lax(iax),nan(10,length(legstr)));
    if ismember(iax,[2 3])
        set(h{iax},'LineWidth',highlightwidth)
    end
    colormap(v_colormap('v_thermliny'));
    legend(lax(iax),...
        h{iax}((iax-1)*nRow + (1:nRow)),...
        legstr((iax-1)*nRow + (1:nRow)),...
        'position',lax(iax).Position,...
        'interpreter','none',...
        'orientation','vertical',...
        'color','none',...
        'edgecolor','none');
end
set(lax,'visible','off')
alltext = findall(gcf,'-property','FontSize');%findall(gcf,'Type','text');
set(alltext,'FontSize',0.8*font_size);
set(lax,'linewidth',1.1)
ahm_print_to_pdf(fh,fullfile(graphics_dir,...
    sprintf('fig_%s_summary_effect_of_order_legend.pdf',...
    exp_id)));
%bug there is a black line?

%% plot the results
fig_width = column_width;
fig_height = 5/8 * fig_width;
for imic = 1:length(mic_lab)
    for iftarget = 1:length(f_target_seq)
        f_target = f_target_seq(iftarget);
        
        mean_mse_mat = permute(mean(comb_mse_mat{imic,iftarget},1),[3 2 1]);
        mean_nr_mat = -permute(mean(comb_nr_mat{imic,iftarget},1),[3 2 1]);

        excess_noisewrt_gt = repmat(mean_nr_mat(:,1),1,size(mean_nr_mat,2)-1)-mean_nr_mat(:,2:end)
        
        figure;
        fh = gcf;
        setFigureSize(fh,[fig_width fig_height],'centimeters');
        
        ax2 = gca;
        colormap(v_colormap('v_thermliny'))
        h2 = bar(ax2,excess_noisewrt_gt);
        set(h2(3:6),'LineWidth',highlightwidth)
        set(ax2,'ylim',ylims{imic})
        ylabel('$\Delta\gamma$ \textrm{[dB]}','Interpreter','latex')
        xlabel('True order of noise field, $N_{s}$','interpreter','latex')
        
        %% add text overlay
        ax2.YTickMode = 'manual';
        ax2.YTickLabelMode = 'manual';
        baseline = ax2.YLim(2);
        ax2.YLim = ax2.YLim .* 1.3;
        topline = ax2.YLim(2);
        text_middle = baseline+2*(topline-baseline)/3;
        for ientry = 1:size(mean_nr_mat(:,1),1)
            text(ientry,text_middle,...
                sprintf('%2.2f',mean_nr_mat(ientry,1)),...
                'HorizontalAlignment','center','VerticalAlignment','middle');
        end
        
        this_xlims = ax2.XLim;
        right_anchor = this_xlims(1) - 0.007*diff(this_xlims);
        text(right_anchor,text_middle,'$\gamma$ [dB]:',...
            'interpreter','latex',...
            'HorizontalAlignment','right','VerticalAlignment','middle');
        
        title(mic_lab_plot{imic},...
            'HorizontalAlignment','center',...
            'FontWeight','normal',...
            'Interpreter','latex')
        
        %% export
        alltext = findall(gcf,'-property','FontSize');
        set(alltext,'FontSize',0.8*font_size);
        alltext = findall(gcf,'Type','text')
        for ih = 1:length(alltext)
            if strcmp(alltext(ih).Interpreter,'latex')
                set(alltext(ih),'FontSize',font_size);
            end
        end
        set([ax2],'linewidth',1.1)
        ahm_print_to_pdf(fh,fullfile(graphics_dir,...
            sprintf('fig_%s_summary_effect_of_order_mic_%s_freq_%d_Hz.pdf',...
            exp_id,mic_lab{imic},f_target_seq(iftarget))))
        
    end
end


