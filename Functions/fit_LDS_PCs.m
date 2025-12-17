function [Max_eigen,Exec,idx_dir2,idx_duration2,r2_all,all_eigen,All_A]=fit_LDS_PCs(score,idx_dir,idx_duration,variance,t_from,t_upto,threshold,threshold_rec,do_plot,binsl)
%% fit_LDS_PCs fits the neural trajectories in the selected recording to an LDS
%
%% INPUTS
%
% session: session name e.g.'MM_S2_raw.mat'
%
% Area: Area to be analysed e.g. 'M1'
%
% OUTPUTS
%
% max_eigen: maximum eigenvalue for each segment.
%
%
% Exec: binary number indicating if that timebin corresponds to movement preparation
% or execution
%
%
% Andrea Colins Rodriguez
% 05/05/2023
%
% 09/05/2023
% % Note: for matlab imag(nan)=0; that's silly! so we need to especify a
% % "complex nan" as complex(nan,nan) for specific cases. Do not take the
% % average of the entire complex number, take the separate parts


ms=1000; %to transform from S to ms
%threshold_rec=10;

ndim=max(find(cumsum(variance)>=threshold,1,'First'),3);


traj_length=(t_upto-t_from)*ms;

Ndir=numel(unique(idx_dir));
Nbins=numel(unique(idx_duration));
t_upto=ones(Nbins,1).*t_upto;
colour_dir=hsv(Ndir);

%prep_eig_traj=nan(abs(t_from)*ms,2,Ndir*Nbins);
Max_eigen=nan(size(idx_dir,1),1);
Exec=nan(size(idx_dir,1),1);
idx_dir2=nan(size(idx_dir,1),1);
idx_duration2=nan(size(idx_dir,1),1);
r2_all=nan(size(idx_dir,1),1);
all_eigen=cell(Nbins,1);
%i_traj=1;

if do_plot
    figure
end

% recenter all trajectories to region of recurrence
%score=recentre_activity(score,idx_dir,idx_duration);
%[score]=geometrical_center_activity(score(:,1:ndim),idx_dir,idx_duration);
%Distance=pdist(score(:,1:ndim));
%recurrence=squareform(Distance);
%imagesc(recurrence)
%limit=prctile(Distance,threshold_rec);
do_video=0;

if do_video
figure
end

counter_cond=1;

%imagesc(recurrence<limit)
for i_dir=1:Ndir
    
    Period_prep=nan(Nbins,1);
    Period_exec=nan(Nbins,1);
    
    for i_dur=1:Nbins
        
        idx_this=idx_dir==i_dir & idx_duration==i_dur;
        This_score=score(idx_this,1:1+ndim-1);%-mean(score(idx_this,1:ndim),1);

        % Define slength according to the distance
        
        nTpoints=size(This_score,1);
        eigsA=nan(ndim,nTpoints);
        All_A=nan(ndim,ndim,nTpoints);
        r2=nan(nTpoints,1);

        % constant lenght
        if isempty(binsl)
        bins=round(nTpoints/4);
        else
            bins=binsl;
        end
        %bins=125;
        sLength=nan(nTpoints,1);
        sLength(bins+1:end-bins-1)=bins;

        %sLength=calculate_length(This_score,limit);
        %sLength(51:end-51)=50;
        t=(1:nTpoints)+t_from*ms;
        if do_video
         counter=1;

        end
        % recenter

        %This_score=This_score-mean(This_score,1);
        %This_score=This_score-mean(This_score([end],:),1);

        parfor i_time=1:nTpoints
            if ~isnan(sLength(i_time))
                %[i_time sLength(i_time)]
                [eigsA(:,i_time),r2(i_time),r1,All_A(:,:,i_time)]=LDS(This_score(i_time-sLength(i_time):i_time+sLength(i_time),:));
                % debugging
                if do_video && i_dir==1 && i_dur==2

                            % video russo
%                            plot(This_score(
% :,1),This_score(:,2),'Color','k','LineWidth',2)
%                             hold on
%                             plot(This_score(1,1),This_score(1,2),'ok','MarkerFaceColor','k')
%                             plot(r1(:,1),r1(:,2),'Color',[0.8 0.2 0.2],'LineWidth',2)
%                             plot(r1(sLength(i_time),1),r1(sLength(i_time),2),'o','Color',[0.8 0.2 0.2],'MarkerFaceColor',[0.8 0.2 0.2])
%                             hold off

                            %xlim([-0.3000 0.6])
                            %ylim([-0.3 0.3032])
%                             title(['Time to mov onset = ' num2str(t(i_time)) ' [ms]'])
%                             xlabel('PC 1')
%                             ylabel('PC 2')
%                             pause(0.01)
                            %print(['./Output_files/fit_LDS_example/im_' num2str(counter)],'-dpng')
 

                            %% create a video for debugging
%                             subplot(5,1,1:4)
%                             plot3(This_score(:,1),This_score(:,2),This_score(:,3),'Color',colour_dir(i_dir,:),'LineWidth',2)
%                             hold on
%                             plot3(This_score(1,1),This_score(1,2),This_score(1,3),'o','Color',colour_dir(i_dir,:),'MarkerFaceColor',colour_dir(i_dir,:))
% 
%                             plot3(0,0,0,'ok','MarkerFaceColor','k')
%                             plot3(r1(:,1),r1(:,2),r1(:,3),'Color',[0.5 0.5 0.5],'LineWidth',2)
%                             plot3(r1(sLength(i_time),1),r1(sLength(i_time),2),r1(sLength(i_time),3),'o','Color',[0.5 0.5 0.5],'MarkerFaceColor',[0.5 0.5 0.5])
%                             hold off
%                             %view(-149.7,50)
% %                             xlim([-0.01 0.03])
% %                             ylim([-0.01 0.01])
% %                             zlim([-0.01 0.01])
%                             title(['Time to mov onset = ' num2str(t(i_time)) ' [ms]'])
%                             xlabel('PC 1')
%                             ylabel('PC 2')
%                             zlabel('PC 3')
%                             pause(0.01)
%                             print(['./Output_files/fit_LDS_example/im_' num2str(counter)],'-dpng')
%                             counter=counter+1;
                 end
           end
        end
        
        max_eig=max(eigsA,[],1);
        all_eigen{counter_cond}=eigsA;
         
         counter_cond=counter_cond+1;

         % ignore all ill fitted segments
        max_eig(r2<0.3)=complex(nan,nan);
        
        % if the imaginary part is zero, then replace by nan
        
        % Note: for matlab imag(nan)=0;
        % that's silly! so we change it for nan+nan*i
         
        idx_imag=imag(max_eig)==0;
        max_eig(idx_imag)=complex(real(max_eig(idx_imag)),nan);
        
        
        max_eig_re=real(max_eig);
        max_eig_imag=imag(max_eig);
        
       
        Mexec=double(t>0);
        Mexec(t>t_upto(i_dur)*ms)=2;

        Period_prep(i_dur)=2*pi/mean(max_eig_imag(Mexec==0),'omitnan');
        Period_exec(i_dur)=2*pi/mean(max_eig_imag(Mexec==1),'omitnan');
        %prep_eig_traj(:,:,i_traj)=[max_eig_re(~Mexec)', max_eig_imag(~Mexec)'];
        
        %i_traj=i_traj+1;
        
        % store outputs
        Max_eigen(idx_this)=max_eig;
        Exec(idx_this)=Mexec;
        idx_dir2(idx_this)=idx_dir(idx_this);
        idx_duration2(idx_this)=idx_duration(idx_this);
        r2_all(idx_this)=r2;
        
        %% PLOTS
        if do_video
            subplot(5,1,5)
            plot(t,movmedian(max_eig_re,50,'omitnan'),'Color',colour_dir(i_dir,:)./sqrt(i_dur))
            hold on
        end

        if do_plot
            kpoints=1;% ms ~ filter lenght

            subplot(4,4,1:3)
            plot(t,movmedian(max_eig_re,kpoints,'omitnan'),'.','Color',colour_dir(i_dir,:)./sqrt(i_dur))
            hold on
            
            subplot(4,4,5:7)
            plot(t,movmedian(max_eig_imag,kpoints,'omitnan'),'.','Color',colour_dir(i_dir,:)./sqrt(i_dur))
            hold on
            
            subplot(4,4,9:11)
            plot(t,r2,'Color',colour_dir(i_dir,:)./sqrt(i_dur))
            hold on
            
            subplot(4,4,13:15)
            plot(t,sLength*2,'Color',colour_dir(i_dir,:)./sqrt(i_dur))
            hold on
            
            subplot(4,4,4)
            %plot(max_eig_re(~Mexec),max_eig_imag(~Mexec),'.-','Color',colour_dir(i_dir,:)/sqrt(i_dur))
            plot(mean(max_eig_re(Mexec==0),'omitnan'),mean(max_eig_imag(Mexec==0),'omitnan'),'.','Color',colour_dir(i_dir,:)./sqrt(i_dur))
            hold on
            
            subplot(4,4,8)
            %plot(max_eig_re(Mexec),max_eig_imag(Mexec),'.-','Color',colour_dir(i_dir,:)/sqrt(i_dur))
            plot(mean(max_eig_re(Mexec==1),'omitnan'),mean(max_eig_imag(Mexec==1),'omitnan'),'.','Color',colour_dir(i_dir,:)./sqrt(i_dur))
            hold on
%             figure
%             
%             plot3(score(idx_dir==i_dir & idx_duration==i_dur,1),score(idx_dir==i_dir & idx_duration==i_dur,2),score(idx_dir==i_dir & idx_duration==i_dur,3),'Color',colour_dir(i_dir,:))

        end
        
        
    end
    if do_plot
        subplot(4,4,12)
        plot((t_upto-t_from)*ms,Period_prep,'.','Color',colour_dir(i_dir,:))
        hold on
        
        subplot(4,4,16)
        plot((t_upto-t_from)*ms,Period_exec,'.','Color',colour_dir(i_dir,:))
        hold on
    end
    
    
end

if do_plot
    subplot(4,4,1:3)
    title(['Real eigenvalue, ndim = ' num2str(ndim)])
    box off
    plot([t(1)  t(end)],[0 0],'k','LineWidth',2)
    plot([0  0],[-0.01 0.01],'k','LineWidth',2)
    plot([t_upto(end) t_upto(end)]*ms,[-0.01 0.01],'k','LineWidth',2)
    xlim([t(1)  max(traj_length)])
    
    subplot(4,4,5:7)
    title('Imag eigenvalue')
    box off
    xlim([t(1)  max(traj_length)])
    
    subplot(4,4,9:11)
    title('R')
    box off
    xlim([t(1)  max(traj_length)])
    ylim([0 1])
    
    subplot(4,4,13:15)
    box off
    xlabel('Time from movement onset [ms]')
    title('Segment Length [ms]')
    xlim([t(1)  max(traj_length)])
    
    subplot(4,4,4)
    title('Mean eigs prep')
    xlabel('Real')
    ylabel('Im')
    plot([0 0],[0 0.02],'k')
    
    box off
    
    subplot(4,4,8)
    title('Mean eigs exec')
    xlabel('Real')
    ylabel('Im')
    
    plot([0 0],[0 0.02],'k')
    box off
    
    
    subplot(4,4,12)
    title('Pred traj period prep')
    xlabel('Trajectory length')
    ylabel('Predicted rotation period [ms]')
    box off
    plot(traj_length,traj_length,'k')
    ylim([0 max(traj_length)*1.1])
    
    subplot(4,4,16)
    title('Pred traj period Exec')
    xlabel('Trajectory length')
    ylabel('Predicted rotation period [ms]')
    box off
    plot(traj_length,traj_length,'k')
    ylim([0 max(traj_length)*1.1])
 
end

Exec=logical(Exec);
end