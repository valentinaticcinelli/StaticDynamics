%% file format
folder = 'C:\Users\ticcinel\Dropbox\Dynamic_exp_QUEST\data\';
formatSpec = '%f %f %f %f %f %f %f %f %f %f %s %s %s %f\n';
VaribleName = {'Trial', 'Resp', 'ACC', 'rt', 'percent', ...
    'cRect1', 'cRect2', 'cRect3', 'cRect4' 'randseq', ...
    'Condition', 'Expression', 'Stim_id', 'randseed'};
Resp_exp = {'Fear', 'Anger', 'Disgust', 'Happiness', 'Sadness', 'Surprise'};
Resp_id  = {'Male1', 'Female1', 'Male2', 'Female2', 'Male3', 'Female3', 'Male4', 'Female4'};
Name_id  = {'Antoine', 'Sophie', 'Didier', 'Fanny', 'Gregoire', 'Helene', 'Joseph', 'Katia'};
%% list subject name
filename_all = dir([folder, '*_g*.csv']);
subjnametemp = [];
for ifile = 1:length(filename_all)
    temp = strsplit(filename_all(ifile).name,'_');
    subjnametemp = [subjnametemp; temp(1:3)];
end
[a,b] = unique(subjnametemp(:,1));
Infotbl = cell2table(subjnametemp(b, :), 'VariableNames', {'sbjname', 'gender', 'age'});
Ns = size(Infotbl, 1);
Infotbl.gender = nominal(Infotbl.gender, {'female', 'male'});
agevec = cell2mat(Infotbl.age);
Infotbl.age = str2num(agevec(:,2:end));
%% plotting parameters
nboot = 1000;
nbin = 100;
q_x1 = linspace(0, 1, 1000);
q_x2 = linspace(0, 1, nbin);
scrsz  = get(0,'ScreenSize');% get screen size for output display
cc = [1, .5, 0; 0, .5, 1];

%% check each subject

for is = 15:Ns
    subjectname = Infotbl.sbjname{is};
    
    %% Expression with both uniform and sampled noise
    filename = dir([folder,subjectname,'*_DyExpthr.csv']);
    if isempty(filename)
        warning(['Subject ', subjectname, ' has not done the expression task.'])
    else
        ds = table;
        for file = 1:length(filename)
            dstmp = readtable([folder,filename(file).name],...
                'Delimiter',' ', 'ReadVariableNames', false, 'Format', formatSpec);
            dstmp.Properties.VariableNames = VaribleName;
            ds = [ds;dstmp];
        end
        
        Cond = unique(ds.Condition);
        Nc = length(Cond);
        Exp = unique(ds.Expression);
        Nexp = length(Exp);
        
        f1 = figure('Numbertitle', 'off', 'Name', subjectname, 'Position', [1 1 scrsz(3)-100 scrsz(4)-100]); hold on
        for ie = 1:Nexp
            subplot(2, 3, ie)
            l2 = [];
            for ic = 1:Nc
                idx = strcmp(ds.Expression, Exp{ie}) & strcmp(ds.Condition, Cond{ic});
                x1 = ds.percent(idx);
                y1 = ds.ACC(idx);
                
                % Bayesian update testing display psychometric function
                
                [b,dev,stats] = glmfit(x1, y1, 'binomial', 'link', 'probit');
                yfit = glmval(b, q_x1, 'probit');
                
                % bootstrapping
                ybootfit = zeros(nboot, nbin);
                for ib = 1:nboot
                    randidx = randi(length(x1), 1, length(x1));
                    [b,dev,stats] = glmfit(x1(randidx), y1(randidx), 'binomial', 'link', 'probit');
                    ybootfit(ib, :) = glmval(b, q_x2, 'probit');
                end
                
                plot(x1, y1+(ic-1.5)*.05, 'o', 'color', cc(ic,:)); % observed data
                hold on
                l2(ic) = plot(q_x1, yfit, 'LineWidth', 2, 'color', cc(ic,:));
                plot(q_x2, prctile(ybootfit, [2.5, 97.5]), '--', 'LineWidth', 1.5, 'color', cc(ic,:))
                xlim([-0.05,1.05])
                ylim([-0.1,1.1])
                title([Exp{ie}])
            end
        end
        subplot(2, 3, ie)
        xlabel('Presented Signal')
        ylabel('Response')
        legend(l2, Cond)
        saveas(f1, ['./sbj_check_fig/', subjectname, '_Exp.png'])
    end
    %% Identity
    filename = dir([folder,subjectname,'*_DyIdthr.csv']);
    if isempty(filename)
        warning(['Subject ', subjectname, ' has not done the Identity task.'])
    else
        ds = table;
        for file = 1:length(filename)
            dstmp = readtable([folder,filename(file).name],...
                'Delimiter',' ', 'ReadVariableNames', false, 'Format', formatSpec);
            dstmp.Properties.VariableNames = VaribleName;
            ds = [ds;dstmp];
        end
        
        Cond = unique(ds.Condition);
        Nc = length(Cond);
        Id = unique(ds.Stim_id);
        Nid = length(Id);
        
        f2 = figure('Numbertitle', 'off', 'Name', subjectname, 'Position', [1 1 scrsz(3)-100 scrsz(4)-100]); hold on
        for jd = 1:Nid
            subplot(2, 4, jd)
            l2 = [];
            for ic = 1:Nc
                idx = strcmp(ds.Stim_id, Id{jd}) & strcmp(ds.Condition, Cond{ic});
                x1 = ds.percent(idx);
                y1 = ds.ACC(idx);
                
                % Bayesian update testing display psychometric function
                
                [b,dev,stats] = glmfit(x1, y1, 'binomial', 'link', 'probit');
                yfit = glmval(b, q_x1, 'probit');
                
                % bootstrapping
                ybootfit = zeros(nboot, nbin);
                for ib = 1:nboot
                    randidx = randi(length(x1), 1, length(x1));
                    [b,dev,stats] = glmfit(x1(randidx), y1(randidx), 'binomial', 'link', 'probit');
                    ybootfit(ib, :) = glmval(b, q_x2, 'probit');
                end
                
                plot(x1, y1+(ic-1.5)*.05, 'o', 'color', cc(ic,:)); % observed data
                hold on
                l2(ic) = plot(q_x1, yfit, 'LineWidth', 2, 'color', cc(ic,:));
                plot(q_x2, prctile(ybootfit, [2.5, 97.5]), '--', 'LineWidth', 1.5, 'color', cc(ic,:))
                xlim([-0.05,1.05])
                ylim([-0.1,1.1])
                title([Name_id(strcmp(Resp_id,Id{jd}))])
            end
        end
        subplot(2, 4, jd)
        xlabel('Presented Signal')
        ylabel('Response')
        legend(l2, Cond)
        saveas(f2, ['./sbj_check_fig/', subjectname, '_Id.png'])
    end
    %pause
    close all
end


%% save table
dsall = table;
for is = 1:Ns
    subjectname = Infotbl.sbjname{is};
    % Expression
    filename = dir([folder,subjectname,'*_DyExpthr.csv']);
    if isempty(filename)
        warning(['Subject ', subjectname, ' has not done the expression task.'])
    else
        ds = table;
        for file = 1:length(filename)
            dstmp = readtable([folder,filename(file).name],...
                'Delimiter',' ', 'ReadVariableNames', false, 'Format', formatSpec);
            dstmp.Properties.VariableNames = VaribleName;
            temp = strsplit(filename(file).name,'_');
            dstmp.randtype = repmat(temp(4), [size(dstmp, 1), 1]);
            dstmp.tasktype = repmat({'Expression'}, [size(dstmp, 1), 1]);
            ds = [ds;dstmp];
        end
        ds.subjectname = repmat({subjectname}, [size(ds, 1), 1]);
        dsall = [dsall; ds];
    end
    % Identity
    filename = dir([folder,subjectname,'*_DyIdthr.csv']);
    if isempty(filename)
        warning(['Subject ', subjectname, ' has not done the Identity task.'])
    else
        ds = table;
        for file = 1:length(filename)
            dstmp = readtable([folder,filename(file).name],...
                'Delimiter',' ', 'ReadVariableNames', false, 'Format', formatSpec);
            dstmp.Properties.VariableNames = VaribleName;
            temp = strsplit(filename(file).name,'_');
            dstmp.randtype = repmat(temp(4), [size(dstmp, 1), 1]);
            dstmp.tasktype = repmat({'Identity'}, [size(dstmp, 1), 1]);
            ds = [ds;dstmp];
        end
        ds.subjectname = repmat({subjectname}, [size(ds, 1), 1]);
        dsall = [dsall; ds];
    end
end
t = datetime('now','TimeZone','local','Format','yMMMd');
writetable(dsall, ['Rawdata', char(t), '.csv'])


%% Quick summary
dsall = readtable('Rawdata2018May30.csv');
sbjname = unique(dsall.subjectname);
Ns = length(sbjname);
Cond = unique(dsall.Condition);
Nc = length(Cond);
Exp = unique(dsall.Expression);
Nexp = length(Exp);
Id = unique(dsall.Stim_id);
Nid = length(Id);
Task = unique(dsall.tasktype);
%
Exp_fitted = nan(Ns, Nexp, Nc, 2);
Id_fitted = nan(Ns, Nid, Nc, 2);
for is = 1:Ns
    %% Expression
    index1 = strcmp(dsall.subjectname, sbjname(is)) & ...
        strcmp(dsall.tasktype, 'Expression');
    if sum(index1)==0
        warning(['Subject ', sbjname{is}, ' has not done the expression task.'])
    else
        ds = dsall(index1, :);
        for ie = 1:Nexp
            for ic = 1:Nc
                idx = strcmp(ds.Expression, Exp{ie}) & strcmp(ds.Condition, Cond{ic});
                x1 = ds.percent(idx);
                y1 = ds.ACC(idx);
                % Bayesian update testing display psychometric function
                [b,dev,stats] = glmfit(x1, y1, 'binomial', 'link', 'probit');
                Exp_fitted(is, ie, ic, :) = b;
            end
        end
    end
    %% Identity
    index2 = strcmp(dsall.subjectname, sbjname(is)) & ...
        strcmp(dsall.tasktype, 'Identity');
    if sum(index2)==0
        warning(['Subject ', sbjname{is}, ' has not done the identity task.'])
    else
        ds = dsall(index2, :);
        for jd = 1:Nid
            for ic = 1:Nc
                idx = strcmp(ds.Stim_id, Id{jd}) & strcmp(ds.Condition, Cond{ic});
                x1 = ds.percent(idx);
                y1 = ds.ACC(idx);
                % Bayesian update testing display psychometric function
                [b,dev,stats] = glmfit(x1, y1, 'binomial', 'link', 'probit');
                Id_fitted(is, jd, ic, :) = b;
            end
        end
    end
end


%% Plotting Expression
trimm_perc = 15;
f1 = figure('Numbertitle', 'off', 'Name', 'Expression Task', 'Position', [1 1 scrsz(3) scrsz(4)]); hold on
for ie = 1:Nexp
    subplot(2, 3, ie); hold on
    l2 = [];
    for ic = 1:Nc
        bvect = squeeze(Exp_fitted(:, ie, ic, :));
        b = trimmean(bvect, trimm_perc, 'round', 1)';
        yfit = glmval(b, q_x1, 'probit');
        
        % bootstrapping
        ybootfit = zeros(nboot, nbin);
        for ib = 1:nboot
            randidx = randi(Ns, 1, Ns);
            b_boot = trimmean(bvect(randidx, :), trimm_perc, 'round', 1)';
            ybootfit(ib, :) = glmval(b_boot, q_x2, 'probit');
        end
        
        l2(ic) = plot(q_x1, yfit, 'LineWidth', 2, 'color', cc(ic,:));
        plot(q_x2, prctile(ybootfit, [2.5, 97.5]), '--', 'LineWidth', 1.5, 'color', cc(ic,:))
        xlim([-0.05,1.05])
        ylim([-0.1,1.1])
    end
    plot(q_x2, ones(size(q_x2))/Nexp, 'k')
    title([Exp{ie}])
end
subplot(2, 3, ie)
xlabel('Presented Signal')
ylabel('Response')
legend(l2, Cond)
saveas(f1, ['./Summary_stats/Expressions.png'])

f2 = figure('Numbertitle', 'off', 'Name', 'Expression Task - overall', 'Position', [1 1 scrsz(3)/2 scrsz(4)/2]); hold on
l2 = [];
for ic = 1:Nc
    bvect = squeeze(nanmean(Exp_fitted(:, :, ic, :), 2));
    b = trimmean(bvect, trimm_perc, 'round', 1)';
    yfit = glmval(b, q_x1, 'probit');
    
    % bootstrapping
    ybootfit = zeros(nboot, nbin);
    for ib = 1:nboot
        randidx = randi(Ns, 1, Ns);
        b_boot = trimmean(bvect(randidx, :), trimm_perc, 'round', 1)';
        ybootfit(ib, :) = glmval(b_boot, q_x2, 'probit');
    end
    
    l2(ic) = plot(q_x1, yfit, 'LineWidth', 2, 'color', cc(ic,:));
    plot(q_x2, prctile(ybootfit, [2.5, 97.5]), '--', 'LineWidth', 1.5, 'color', cc(ic,:))
    xlim([-0.05,1.05])
    ylim([-0.1,1.1])
end
plot(q_x2, ones(size(q_x2))/Nexp, 'k')
title('Expression task')
xlabel('Presented Signal')
ylabel('Response')
legend(l2, Cond)
saveas(f2, ['./Summary_stats/Expression_overall.png'])

% Plotting Identity
f1 = figure('Numbertitle', 'off', 'Name', 'Identity Task', 'Position', [1 1 scrsz(3) scrsz(4)]); hold on
for jd = 1:Nid
    subplot(2, 4, jd); hold on
    l2 = [];
    for ic = 1:Nc
        bvect = squeeze(Id_fitted(:, jd, ic, :));
        b = trimmean(bvect, trimm_perc, 'round', 1)';
        yfit = glmval(b, q_x1, 'probit');
        
        % bootstrapping
        ybootfit = zeros(nboot, nbin);
        for ib = 1:nboot
            randidx = randi(Ns, 1, Ns);
            b_boot = trimmean(bvect(randidx, :), trimm_perc, 'round', 1)';
            ybootfit(ib, :) = glmval(b_boot, q_x2, 'probit');
        end
        
        l2(ic) = plot(q_x1, yfit, 'LineWidth', 2, 'color', cc(ic,:));
        plot(q_x2, prctile(ybootfit, [2.5, 97.5]), '--', 'LineWidth', 1.5, 'color', cc(ic,:))
        xlim([-0.05,1.05])
        ylim([-0.1,1.1])
    end
    plot(q_x2, ones(size(q_x2))/Nid, 'k')
    title([Name_id(strcmp(Resp_id,Id{jd}))])
end
subplot(2, 4, jd)
xlabel('Presented Signal')
ylabel('Response')
legend(l2, Cond)
saveas(f1, ['./Summary_stats/Identities.png'])

f2 = figure('Numbertitle', 'off', 'Name', 'Identity Task - overall', 'Position', [1 1 scrsz(3)/2 scrsz(4)/2]); hold on
l2 = [];
for ic = 1:Nc
    bvect = squeeze(nanmean(Id_fitted(:, :, ic, :), 2));
    b = trimmean(bvect, trimm_perc, 'round', 1)';
    yfit = glmval(b, q_x1, 'probit');
    
    % bootstrapping
    ybootfit = zeros(nboot, nbin);
    for ib = 1:nboot
        randidx = randi(Ns, 1, Ns);
        b_boot = trimmean(bvect(randidx, :), trimm_perc, 'round', 1)';
        ybootfit(ib, :) = glmval(b_boot, q_x2, 'probit');
    end
    
    l2(ic) = plot(q_x1, yfit, 'LineWidth', 2, 'color', cc(ic,:));
    plot(q_x2, prctile(ybootfit, [2.5, 97.5]), '--', 'LineWidth', 1.5, 'color', cc(ic,:))
    xlim([-0.05,1.05])
    ylim([-0.1,1.1])
end
plot(q_x2, ones(size(q_x2))/Nid, 'k')
title('Identity task')
xlabel('Presented Signal')
ylabel('Response')
legend(l2, Cond)
saveas(f2, ['./Summary_stats/Identity_overall.png'])