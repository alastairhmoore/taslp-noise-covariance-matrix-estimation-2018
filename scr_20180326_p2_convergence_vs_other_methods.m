% collate the results for effect of diffuse field SH order - model vs gt
exp_id = '20180326_x2';

nFields = 20;
f_target = 2200;
sphHarmOrdDiffuseGt = 2;
sphHarmOrdDiffuseEst = 2;
nSH = (sphHarmOrdDiffuseEst+1)^2;


nMethods = 13;

cov_mse = cell(nMethods-1,nFields);

            scr_id = sprintf('%s_gt_o%d_%d_Hz',...
                exp_id,sphHarmOrdDiffuseGt,...
                f_target);


            load(sprintf('%s/dat_%s.mat',data_dir,scr_id))
for imethod = 1:nMethods
            for ifield = 1:nFields
                cov_mse{imethod,ifield} = est(ifield,1,imethod+1).est_cov_mean_sq_err;
            end
end

%% average error
cum_error = zeros(size(cov_mse{1,1},1),nMethods);
for imethod = 1:nMethods
            for ifield = 1:nFields
                cum_error(:,imethod) = cum_error(:,imethod) + cov_mse{imethod,ifield};
            end
end
avg_error_db = 10*log10(cum_error./nFields);

%%
for iplot = 1:2
fig_width = column_width
fig_height = fig_width

nPoses = 5;
nIterPerPose = 250;
xlims = [0 nPoses*nIterPerPose];
ylims = [-40 10];

switch iplot
    case 1
        fig_lab = 'rs'
        methods_to_plot = [1 2 3:7]
legend_labels = {'White','Sph iso',...
    'RS: 1e-1','RS: 5e-2','RS: 1e-2','RS: 5e-3','RS: 1e-3'};
line_colors = [0 0 0;0 0 0;...
    get_voicebox_thermliny_linecolors(5)];
line_styles = {'-.',':',...
    '--','-','--','-','--'};
line_widths = [1 1,...
    1.5 * ones(1,5)];
    case 2
        fig_lab = 'ewls'
        methods_to_plot = [1 2 8:13]
        legend_labels = {'White','Sph iso',...
            'EWLS: 1e-1','EWLS: 1e-2','EWLS: 1e-3',...
            'EWLS: 1e-4','EWLS: 1e-5','EWLS: 1e-6'}
        line_colors = [0 0 0;0 0 0;...
            get_voicebox_thermliny_linecolors(3);
            get_voicebox_thermliny_linecolors(3)];
        line_styles = {'-.',':',...
            '--','--','--',...
            '-','-','-'};
        line_widths = [1 1,...
             1.5 * ones(1,3),...
             1.5 * ones(1,3)];
end

figure;
fh = gcf;
setFigureSize(fh,[fig_width fig_height],'centimeters');
clear lineh
for imethod = 1:length(methods_to_plot)
    lineh(imethod) = plot(avg_error_db(:,methods_to_plot(imethod)),...
        'color',line_colors(imethod,:),...
        'linestyle',line_styles{imethod},...
        'linewidth',line_widths(imethod));
    hold all
end



scr_20171206_05_frame_index_xlabel;
set(gca,'xtick',(1:nPoses)*nIterPerPose)
ylabel('$\mathcal{E}$ \textsf{[dB]}','Interpreter','latex')

xlim(xlims);
ylim(ylims);


lh = legend(legend_labels);

set(lh,'Location','SouthWest','Orientation','Vertical');
set(lh,'AutoUpdate','off')
fprintf('avg_error_db_at_end\n')
legend_labels
avg_error_db(end,methods_to_plot)

drawnow
pos = get(gca,'position');
set(gca,'pos',[pos(1)-0.03,pos(2)-0.03,pos(3)+0.06,pos(4)+0.03])


%reorder so we can see the model
uistack(lineh(1),'top')
uistack(lineh(2),'top')
uistack(lineh(end),'top')

ph = patch([950 1250 1250 950],ylims([1 1 2 2]),[0.9128    0.9666    0.8950]);
ph.LineStyle='none';
uistack(ph,'bottom')

alltext = findall(fh,'-property','FontSize');%findall(gcf,'Type','text');
set(alltext,'FontSize',0.8*font_size);
allax = findobj(fh,'type','axes');
set(allax,'linewidth',1.1)
set(fh,'Renderer','painters')

%% overlay arrows - requires no more adjustment of axis size/scale
frame_centers = (0:nPoses-1)*nIterPerPose + nIterPerPose/2;
yaw_angle = [0 30 60 90 0];
my_scale = 0.15;
plot_scale = get(gca,'DataAspectRatio');
clear xvals yvals
for ii=1:length(frame_centers)
    xvals(ii,:) = frame_centers(ii) + [0 my_scale*plot_scale(1)*cosd(yaw_angle(ii)+90)];
    yvals(ii,:) = 5 + [0 my_scale*plot_scale(2)*sind(yaw_angle(ii)+90)];
end
ax_pos = get(gca,'position');
xlims = get(gca,'xlim');
ylims = get(gca,'ylim');

xvals = ax_pos(1) + (xvals-xlims(1))./diff(xlims) * ax_pos(3);
yvals = ax_pos(2) + (yvals-ylims(1))./diff(ylims) * ax_pos(4);
for ii=1:length(frame_centers)
    annotation('arrow',xvals(ii,:),yvals(ii,:),...
        'linewidth',1.1);
end

ahm_print_to_pdf(fh,fullfile(graphics_dir,sprintf('fig_%s_convergence_compared_%s.pdf',exp_id,fig_lab)))

end