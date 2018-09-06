% collate the results for effect of diffuse field SH order - model vs gt
exp_id = '20180326_x1';

nFields = 20;
nPoseSeq = 3;
f_target = 2200;
sphHarmOrdDiffuseGt = 2;
sphHarmOrdDiffuseEst = 2;
nSH = (sphHarmOrdDiffuseEst+1)^2;
cov_mse = cell(nPoseSeq,nFields);
fnm_hist_all = cell(nPoseSeq,nFields);
gt_fnm = zeros(nFields,nSH);
for iposeseq = 1:nPoseSeq
    
    scr_id = sprintf('%s_gt_o%d_%d_Hz_poseseq_%d',...
        exp_id,sphHarmOrdDiffuseGt,...
        f_target,iposeseq);
    
    
    load(sprintf('%s/dat_%s.mat',data_dir,scr_id))
    
    for ifield = 1:nFields
        cov_mse{iposeseq,ifield} = est(ifield,1,end).est_cov_mean_sq_err;
        fnm_hist_all{iposeseq,ifield} = est(ifield,1,end).model_params.fnm_est_hist(:,1:nSH);
        %fnm_hist_all{iposeseq,ifield} = conj(est(ifield,1,end).model_params.fnm_est_hist(:,1:nSH));
        if iposeseq==1
            %gt_fnm(ifield,:) = 1/(4*pi) * saved_gt.fnm{ifield}.'; %[1,9]
            gt_fnm(ifield,:) = saved_gt.fnm{ifield}.'; %[1,9]
        end
    end
end

%%
pose_lab = {'yaw only','constrained','unconstrained'};
line_style = {'-',':','--'};
line_color = get_voicebox_thermliny_linecolors(3);
line_width = 1.1;

fig_width = column_width;
fig_height = 3/4*fig_width;

fh=figure;
setFigureSize(fh,[fig_width fig_height],'centimeters');

final_cov_error=zeros(1,3);
lh=[];
pose_seq_lab = {'yaw only','constrained rotation','free rotation'};
for iposeseq = 1:nPoseSeq
    
    mean_cov_error = zeros(size(est(1,1,end).est_cov_mean_sq_err));
    for ifield = 1:nFields
        mean_cov_error = mean_cov_error + cov_mse{iposeseq,ifield};
    end
    mean_cov_error = mean_cov_error/nFields;
    final_cov_error(iposeseq) = mean_cov_error(end);
    lh(iposeseq)=plot(10*log10(mean_cov_error),...
        'LineStyle',line_style{iposeseq},...
        'Color',line_color(iposeseq,:),...
        'LineWidth',1.5);
    hold all;
    
end
ylabel('$\mathcal{E}$ \textrm{[dB]}','Interpreter','latex')
scr_20171206_05_frame_index_xlabel
legend(pose_lab,'Location','NorthEast')
%logx
alltext = findall(fh,'-property','FontSize');%findall(gcf,'Type','text');
set(alltext,'FontSize',0.8*font_size);
set(gca,'linewidth',1.1)
set(fh,'Renderer','painters')
ahm_print_to_pdf(fh,fullfile(graphics_dir,sprintf('fig_%s_effect_of_pose_seq_noise_cov.pdf',exp_id)))

fprintf('final_cov_error')
10*log10(final_cov_error)

%%
fig_width = text_width;
fig_height = 1/2*fig_width;

sh_lab = {'p=1\quad(n=0,m=0)',...
    'p=2\quad(n=1,m=-1)',...
    'p=3\quad(n=1,m=0)',...
    'p=4\quad(n=1,m=1)',...
    'p=5\quad(n=2,m=-2)',...
    'p=6\quad(n=2,m=-1)',...
    'p=7\quad(n=2,m=0)',...
    'p=8\quad(n=2,m=1)',...
    'p=9\quad(n=2,m=2)'};
spi = [3, 7,8,9, 11,12,13,14,15];
spi_ylab = [1,2,5];
spi_xlab = [5:9];
nSH = 9;
markers = {'o','x','+'};

xsp_scale = 1.2;
ysp_scale = 1.2;

xlims = [0 5000];
ylims = [-40 10];

final_error = zeros(nSH,nPoseSeq);
fh=figure;
setFigureSize(fh,[fig_width fig_height],'centimeters');
for iSH = 1:9
    subplot(3,5,spi(iSH))
    pos = get(gca,'position');
    new_pos = pos;
    
    new_width = xsp_scale * pos(3);
    diff_width = new_width-pos(3);
    new_pos(1) = pos(1)-diff_width/2;
    new_pos(3) = new_width;
    
    new_height = ysp_scale * pos(4);
    diff_heigth = new_height-pos(4);
    new_pos(2) = pos(2)-diff_heigth/2;
    new_pos(4) = new_height;
    
    
    set(gca,'position',new_pos);
    for iposeseq = 1:nPoseSeq
        av_error = zeros(size(fnm_hist_all{1,1}(:,1)));
        %reference = 1/(4*pi) * saved_gt.fnm{ifield}(1);
        reference = saved_gt.fnm{ifield}(1);

        for ifield = 1:nFields
            %target = 1/(4*pi) * saved_gt.fnm{ifield}(iSH);
            target = saved_gt.fnm{ifield}(iSH);

            av_error = av_error + ...
                abs(target-fnm_hist_all{iposeseq,ifield}(:,iSH))./abs(reference);
        end
        av_error = av_error/nFields;
        final_error(iSH,iposeseq) = av_error(end);
        plot(20*log10(av_error),...
            'LineStyle',line_style{iposeseq},...
            'Color',line_color(iposeseq,:),...
            'LineWidth',1.5);
        hold all;
    end
    text(4800,5,sprintf('$%s$',sh_lab{iSH}),...
        'Interpreter','latex',...
        'HorizontalAlignment','right',...
        'VerticalAlignment','middle');
    xlim(xlims);
    ylim(ylims);
    
    if ismember(iSH,spi_ylab)
        ylabel('$\epsilon_{p}$ \textrm[dB]','Interpreter','latex');
    else
        set(gca,'yticklabel',[])
    end
    if ismember(iSH,spi_xlab)
        scr_20171206_05_frame_index_xlabel;
    else
        set(gca,'xticklabel',[])
    end
end
allax = findobj(fh,'type','axes');

axl=subplot(3,5,5);
for iposeseq=1:nPoseSeq
    plot(NaN(10,1),...
        'LineStyle',line_style{iposeseq},...
        'Color',line_color(iposeseq,:),...
        'LineWidth',1.5);
    hold all
end
legh= legend(pose_lab);
legh.Position = axl.Position;
axl.Visible = 'off';

if 0
    xlims = cell2mat(get(allax,'xlim'));
    ylims = cell2mat(get(allax,'ylim'));
    set(allax,'xlim',[min(xlims(:,1)),max(xlims(:,2))]);
    set(allax,'ylim',[min(ylims(:,1)),max(ylims(:,2))]);
    setappdata(fh,'lp',linkprop(allax,{'xlim','ylim'}))
end


alltext = findall(fh,'-property','FontSize');
set(alltext,'FontSize',0.8*font_size);
set(allax,'linewidth',1.1)
set(fh,'Renderer','painters')
ahm_print_to_pdf(fh,fullfile(graphics_dir,...
    sprintf('fig_%s_effect_of_pose_seq_sh_est.pdf',exp_id)))

fprintf('Summary\nYaw only\tConstrained\tUnconstrained\tDelta Constrained\tDelta Unconstrained\n');
final_error_db = 20*log10(final_error);
[final_error_db final_error_db(:,[2 3])-final_error_db(:,1)]
fprintf('Deg 0 only')
idc = [1 3 7];
[final_error_db(idc,:) final_error_db(idc,[2 3])-final_error_db(idc,1)]
