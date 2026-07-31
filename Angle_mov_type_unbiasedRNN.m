function Angle_mov_type_unbiasedRNN


folderM1 = '.\Output_files\TrainedRNNs\M1_unbiased';

folderSMA = '.\Output_files\TrainedRNNs\SMA_unbiased';

files = dir([folderSMA '\*.mat']);

figure 
subplot(2,2,1)
hold on
angles_random=control_angles(length(files));

plot_fancy_errorbars(1,angles_random,[0.5 0.5 0.5])

compute_angle(folderM1,2)

compute_angle(folderSMA,3)

ylim([0,180])
xlim([0.5 3.5])
xticks([1 2 3])
xticklabels({'Random','M1','SMA'})
ylabel('Angle between rhythmic and discrete movement input vectors [\circ]')
box off

%% Show example networks
% select animal and condition to evaluate in both example RNN
animal = 'Cousteau';
i_dir = 1;
i_pos = 1;
area = 'SMA';

subplot(4,3,11)
rec_number = 4;
plot_dynamics_different(animal,rec_number,area,i_dir,i_pos)

subplot(4,3,12)
rec_number = 19;
plot_dynamics_different(animal,rec_number,area,i_dir,i_pos)



area = 'M1';
subplot(4,3,9)
rec_number = 1;
plot_dynamics_different(animal,rec_number,area,i_dir,i_pos)

subplot(4,3,8)
rec_number = 8;
plot_dynamics_different(animal,rec_number,area,i_dir,i_pos)


end

function theta=control_angles(N)

theta=nan(N,1);
for i=1:N
    %create 2 random vectors
    %u=normrnd(0, 0.3,50,1);
    %v=normrnd(0, 0.3,50,1);

    u=randn(50,1);
    v=randn(50,1);

    % measure angle
    theta(i)=angle_between_vectors(u,v);
end

end

function theta=angle_between_vectors(u,v)

theta=subspace(u(:),v(:))*180/pi;
%CosTheta = max(min(dot(u, v) / (norm(u) * norm(v)), 1), -1);
%theta = real(acos(CosTheta))*180/pi;
end

function compute_angle(folder,i)
files = dir([folder '\*.mat']);

Nfiles = length(files);
theta=nan(Nfiles,1);
for k = 1:Nfiles
    fileName = [folder '\' files(k).name];

    load(fileName, 'Btask');
    u=Btask(:,1);

    v=Btask(:,2);

    theta(k)=angle_between_vectors(u,v);
end

hold on 
plot_fancy_errorbars(i,theta,[0.9 0.1 0.1])
plot(i,median(theta),'o')
end

function plot_dynamics_different(animal,i,area,i_dir,i_pos)
load(['.\Output_files\RNNs_Inputs\' area '_' animal '_different.mat'],'exec','idx_current_cycle')
training=[0,1,2,3,16,17,18,19]+1;
test = 5:16 ;
exec=[exec(training,:);exec(test,:)];
idx_current_cycle=[idx_current_cycle(training,:);idx_current_cycle(test,:)];

if strcmp(area,'SMA')==1
    output='M1';
else
    output='EMG';
end
filename = ['Output_files\TrainedRNNs\' area '_unbiased\Trained_' output '_Hyp_unbiased' animal '_' num2str(i) '.mat'];

info = load_RNN_info(filename);
[scores,trials_idx,~,~,~,~,explained] = Eval_RNN_all_conditions(info.Input,info.Output,info.Test_input,info.Test_Outputs,info.net_params,exec,info.idx_conditions_train,info.idx_conditions_test,0);

[idx_dir,idx_pos,~,idx_dist]=find_idx_conds(trials_idx,info.idx_conds_all,idx_current_cycle);

discrete = find(idx_dir == i_dir & idx_pos == i_pos & idx_dist == 0.5);
rhythmic = find(idx_dir == i_dir & idx_pos== i_pos & idx_dist == 7);
Ndist = numel(unique(idx_dist));
colours=plasma(Ndist);


hold on
plot3(scores(rhythmic(1),1),scores(rhythmic(1),2),scores(rhythmic(1),3),'o','MarkerFaceColor',colours(end,:),'MarkerEdgeColor',colours(end,:))
plot3(scores(rhythmic,1),scores(rhythmic,2),scores(rhythmic,3),'Color',colours(end,:))

plot3(scores(discrete,1),scores(discrete,2),scores(discrete,3),'Color',colours(1,:))
plot3(scores(discrete(1),1),scores(discrete(1),2),scores(discrete(1),3),'o','MarkerFaceColor',colours(1,:),'MarkerEdgeColor',colours(end,:))

load(filename, 'Btask');
theta=angle_between_vectors(Btask(:,1),Btask(:,2));
title(['Recording number = ' num2str(i) ' angle input = ' num2str(theta,3) ' VE = ' num2str(sum(explained(1:3)))])


end
