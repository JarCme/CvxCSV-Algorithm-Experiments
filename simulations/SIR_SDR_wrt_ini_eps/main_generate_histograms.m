% clearvars -except ISR_picked SNR_picked SDR_picked params_to_get_vals methods;
clear;clc;close all;


result_path = ['..',filesep,'cvxCSV_task',filesep,'outputFiles',filesep];

if (exist(result_path,'dir') ~= 7)
    disp('unzipping precomputed results')
    unzip(['..',filesep,'cvxCSV_task',filesep,'outputFiles.zip'],['..',filesep,'cvxCSV_task',filesep]);
    disp('done');
end

result_path = [result_path,'cfg_ini_model_lap_(dot)_m'];

if ~exist('ISR_picked','var')
    [ISR_picked, SNR_picked, SDR_picked, params_to_get_vals, methods] =collectData(result_path, 'ini_pert');
end


data_met = SNR_picked;

f =figure(2);

set(f,'position',[420,295,1500,860]);

spacing = 0.02;
margin_left = 0.1;


edges = -55:2:55;
axs = [];
t_means = [];
t_names = [];
title_hs = [];
plot_idx = 1;
for m_idx = 1:length(methods)
    m = methods(m_idx);
    for ini_idx = 1:length(params_to_get_vals)
        ini = params_to_get_vals(ini_idx);
       
        ax = subaxis(length(methods),length(params_to_get_vals),plot_idx, 'Spacing', spacing, 'Padding', 0, 'MarginLeft', margin_left);
        
        axs = [axs,ax];
        histogram(10*log10(data_met(:,m_idx,ini_idx)),edges);
        xline(0,'--')

        if ini_idx == 1
            name = char(m);
            switch char(m)
                case 'ini_lcmp'
                    name = 'init. LCMP';
                case 'full'
                    name = 'CvxCSV alg.';
                case 'bogice'
                    name = 'BOGICE alg.';
                case 'fastdiva'
                    name = 'FastDIVA alg.';
            end
    

            t_name = text(-0.55,0.5,name,'Units','normalized', ...
                                        'HorizontalAlignment','center', ...
                                        'FontWeight','bold', ...
                                        'FontSize',16,...
                                        'Rotation',90);
            t_names = [t_names,t_name];

    
            ylabel('Count');

        else
            yticks([]);
        end
        title_mean = mean(10*log10(data_met(:,m_idx,ini_idx)));
        if m_idx ==1
            title_h = title(['||\epsilon||_2= ',num2str(ini)]);
            title_hs = [title_hs, title_h];
        end
        t_obj = text(-50,1300,{['mean: '], num2str(title_mean,'%.2f')},'FontWeight','bold');
        t_means = [t_means, t_obj];
        if m_idx == length(methods)
            xlabel('SIR [dB]');
        else
            xticks([]);
        end

        plot_idx = plot_idx+1;
    end
end

linkaxes(axs);
set(axs,'FontWeight','bold')
set(axs,'FontSize',12)

xli = xlim;
yli = ylim;

set(t_means,'Position',[xli(1)+(xli(2)-xli(1))/2*0.15, yli(2)-(yli(2)-yli(1))*0.15, 0 ])


set(f, 'Renderer', 'painters');

saveas(f,'SIR.eps','epsc')



data_met = SDR_picked;

f2 =figure(3);

set(f2,'position',[420,295,1500,860]);

spacing = 0.02;
margin_left = 0.1;


edges = 0:2:105;
axs = [];
t_means = [];
t_names = [];
plot_idx = 1;
for m_idx = 1:length(methods)
    m = methods(m_idx);
    for ini_idx = 1:length(params_to_get_vals)
        ini = params_to_get_vals(ini_idx);
        % ax = subplot(length(methods),length(params_to_get_vals),plot_idx);
        ax = subaxis(length(methods),length(params_to_get_vals),plot_idx, 'Spacing', spacing, 'Padding', 0, 'MarginLeft', margin_left);
        % set(ax,'LooseInset', [0,0,0,0]);
        axs = [axs,ax];
        histogram(10*log10(data_met(:,m_idx,ini_idx)),edges);
        xline(0,'--')
        % legend;
        % ylim([0,150])
        if ini_idx == 1
            
            name = char(m);
            switch char(m)
                case 'ini_lcmp'
                    name = 'init. LCMP';
                case 'full'
                    name = 'CvxCSV alg.';
                case 'bogice'
                    name = 'BOGICE alg.';
                case 'fastdiva'
                    name = 'FastDIVA alg.';
            end
            t_name = text(-0.55,0.5,name,'Units','normalized', ...
                                        'HorizontalAlignment','center', ...
                                        'FontWeight','bold', ...
                                        'FontSize',16,...
                                        'Rotation',90);
            t_names = [t_names,t_name];

    
            ylabel('Count');

        else
            yticks([]);
        end
        title_mean = mean(10*log10(data_met(:,m_idx,ini_idx)));
        if m_idx ==1
            title(['||\epsilon||_2= ',num2str(ini)]);
        % else
        %     title(['mean: ',num2str(title_mean), ' dB']);
        end
        t_obj = text(-50,1300,{['mean: '], num2str(title_mean,'%.2f')},'FontWeight','bold');
        t_means = [t_means, t_obj];
        if m_idx == length(methods)
            xlabel('SDR [dB]');
        else
            xticks([]);
        end
        plot_idx = plot_idx+1;
    end
end

linkaxes(axs);
set(axs,'FontWeight','bold')
set(axs,'FontSize',12)

xli = xlim;
yli = ylim;

set(t_means,'Position',[xli(2)-(xli(2)-xli(1))/2*0.65, yli(2)-(yli(2)-yli(1))*0.15, 0 ])

set(f2, 'Renderer', 'painters');
saveas(f2,'SDR.eps','epsc')





