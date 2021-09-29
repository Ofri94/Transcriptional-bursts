A=[0.15:0.1:1.05];%act concentration in nanoMolar
theta_on=[];sig_on=[];
theta_off=[];sig_off=[];
On=[];

for act=A
    N=6;%total number of binding sites
    states=0:N;%available binding states in the system
    kon=0.008*act;%in Hz*nM.
    koff_active=0.004;%for active state. rate of single dissociation, assuming only 1 dna strand. in Hz.
    koff_inactive=12;%for active state.
    
    rate_Mon=0.0017;%rate for transitioning to active state at 0 occupancy
    rate_Moff=0.005;%rate for transitioning to inactive state at 0 occupancy
    m=1;%memory of system's initial condition
    
    window_dur=240;%in sec
    
    sample_size=1;
    Occ_P=zeros(sample_size,length(states));
    
    for itr2=1:sample_size
        
        n=0;%initial condition - number of sites currently occupied
        t=0;%initial condition - time in sec
        occ=n;%occupancy vector - for convenience
        M=m;
        Mtimes=[M];%timing the states
        while t(end)<360000
            
             current_Mon=rate_Mon;%M appropriate to the occupancy
            current_Moff=rate_Moff/(n+1);%M appropriate to the occupancy
            k_add=(N-n)*kon;%rate constant for adding a TF. depends only on activator concentration
            
            if M==1
                k_dec=n*koff_active;%rate constant for TF dissociation
                Rtot=k_add+k_dec+current_Moff;%Total rate of event
                
                p_add=k_add/Rtot;%probability of adding a TF
                p_dec=k_dec/Rtot;%probability of removing a TF
            elseif M==0
                k_dec=n*koff_inactive;%rate constant for TF dissociation
                Rtot=k_add+k_dec+current_Mon;%Total rate of event
                
                p_add=k_add/Rtot;%probability of adding a TF
                p_dec=k_dec/Rtot;%probability of removing a TF
            end
            
            num=rand;%rand chooses a number between 0-1 from a uniform distribution
            
            if num<=p_add %updating occupancy (creating a vec).
                n=n+1;
            elseif num>p_add && num<=(p_add+p_dec)
                n=n-1;
            elseif num>(p_add+p_dec)
                M=double(~M);%changing mode
            end
            
            t=[t, t(end) + exprnd(1/Rtot)];%creating a vec with updated times
            Mtimes=[Mtimes, M];
            occ=[occ, n];%adding the current occupancy to the occupancy vec
        end
%         OCC=movmean_time(occ,t,window_dur);
        
%         figure;
% %         subplot(1,3,1)
%         plot(t./60,occ)%show time in Minutes
%         title('Transcription rate as a function of time')
%         xlabel('Time [min]');ylabel('Transcription rate [a.u.]')
        
%         figure;
% %         subplot(1,3,2)
%         plot(t./60,OCC)%show time in minutes
%         title('Transcription rate as a function of time')
%         xlabel('Time [min]');ylabel('Transcription rate [a.u.]')
        
        
        %     figure;
        %     plot(t./60,movmean(occ,10))
        %     hold on
        %     plot(t./60, Mtimes)
        
        
        dur=diff(t);%discrete duartaion of modes
        
        M_ind=logical([abs(diff([Mtimes, 2]))]);%indicating locations where the system's state changes.
        %concatanating 2 to make sure that the last durstion is accounted for
        t_states=t(M_ind);%the times at which the system's state changed
        states_dur=diff([0 t_states]);%duration vec for the states. combining on and off states. cocncatanating 0 to fix diff's offset
        states_dur=states_dur./60;%switching units to minutes
        
        state_ind=m:2:length(states_dur);%getting indication to asign state duration for each state. based on the intermittence of the system.
        
        if state_ind(1)==0
            if state_ind(end)==length(states_dur)
                state_ind(1)=[];
                off_dur=states_dur(state_ind-1);%noticing off-state durations
                on_dur=states_dur(state_ind);%noticing on-state durations
            else
                off_dur=states_dur(state_ind+1);%noticing off-state durations
                state_ind(1)=[];%ussing first term in row above
                on_dur=states_dur(state_ind);%noticing on-state durations
            end
        elseif state_ind(1)==1
            if state_ind(end)==length(states_dur)
                on_dur=states_dur(state_ind);%noticing on-state durations
                state_ind(1)=[];%ussing first term in row above
                off_dur=states_dur(state_ind-1);%noticing off-state durations
            else
                on_dur=states_dur(state_ind);%noticing on-state durations
                off_dur=states_dur(state_ind+1);%noticing off-state durations
            end
        else
            error(['state ' num2str(m) ' does not exsit'])
        end
        
        
        %     figure;
        %     subplot(1,3,1)
        %     h=histogram(on_dur,'Normalization','probability');
        %     title({'On durations', ['total time spent - ' num2str(sum(on_dur)./60) ' Hours']})
        %     ylabel('Fraction of epoches');xlabel('Time [min]');
        %     xlim([0, 100]);ylim([0 1])
        %     subplot(1,3,2)
        %     h=histogram(off_dur,'Normalization','probability');
        %     title({'Off durations', ['total time spent - ' num2str(sum(off_dur)./60) ' Hours']})
        %     ylabel('Fraction of epoches');xlabel('Time [min]');
        %     xlim([0, 100]);ylim([0 1])
        %     subplot(1,3,3)
        %     h=histogram(states_dur,'Normalization','probability');
        %     title(['total simulation time - ' num2str((sum(on_dur)+sum(off_dur))./60) ' Hours'])
        %     ylabel('Fraction of epoches');xlabel('Time [min]');
        %     xlim([0, 100]);ylim([0 1])
        
        occ_prob=zeros(1,N+1);%allocating space
        
        theta_on=[theta_on, mean(on_dur)];%average on-burst duration for given activator concentration
        theta_off=[theta_off, mean(off_dur)];%average off-burst duration for given activator concentration
        sig_on=[sig_on, std(on_dur)/sqrt(length(on_dur))];%error of on-burst duration for given activator concentration
        sig_off=[sig_off, std(off_dur)/sqrt(length(off_dur))];%error of off-burst duration for given activator concentration
        
        for itr=states
            occ_prob(itr+1)=sum(dur(occ(1:end-1)==itr))/t(end);%caculating probability of oocurance by temporal duration of a state
        end
        
%                     figure;
%                     bar(states,occ_prob)
%                     title('Occurence of occupancy')
%                     xlabel('Occupancy');ylabel('Occurence')
%                     ylim([0 1])
        
        
        Occ_P(itr2,:)=occ_prob;%concatanating runs
        
        
        th=3;%threshold for active transcription
        occOn_loc=occ>=th;%indices for active-transcription momments
        occ_ind=logical([abs(diff([occOn_loc, 2]))]);%indicating locations where the system's transcriptional activity changes.
        %concatanating 2 to make sure that the last durstion is accounted for
        t_activity=t(occ_ind);%the times at which the system's transcriptional activity changed
        activity_dur=diff([0 t_activity]);%duration vec for the activity. combining on and off states. cocncatanating 0 to fix diff's offset
        activity_dur=activity_dur./60;%switching units to minutes
        
        activity_ind=1:2:length(activity_dur);
        
        if occ(1)<3
            if activity_ind(end)==length(activity_dur)
                
                T_off=activity_dur(activity_ind);%noticing Transcription-off durations
                activity_ind(1)=[];
                T_on=activity_dur(activity_ind-1);%noticing Transcription-on durations
            else
                T_off=activity_dur(activity_ind);%noticing Transcription-off durations
                T_on=activity_dur(activity_ind+1);%noticing Transcription-on durations
            end
        elseif occ(1)>=3
            if activity_ind(end)==length(activity_dur)
                T_on=activity_dur(activity_ind);%noticing Transcription-on durations
                activity_ind(1)=[];%ussing first term in row above
                T_off=activity_dur(activity_ind-1);%noticing Transcription-off durations
            else
                T_on=activity_dur(activity_ind);%noticing Transcription-on durations
                T_off=activity_dur(activity_ind+1);%noticing Transcription-off durations
            end
        else
            error(['state ' num2str(m) ' does not exsit'])
        end
        
    end
    
    mean_op=mean(Occ_P);
    std_op=std(Occ_P)./sqrt(sample_size);
    
    
    
    
    
    %     figure;
    %     subplot(1,3,1)
    %     h=histogram(T_on,'Normalization','probability');
    %     title({'Transcription-on durations', ['total time spent - ' num2str(sum(T_on)./60) ' Hours']})
    %     ylabel('Fraction of epoches');xlabel('Time [min]');
    %     xlim([0, 100]);ylim([0 1])
    %     subplot(1,3,2)
    %     h=histogram(T_off,'Normalization','probability');
    %     title({'Transcription-off durations', ['total time spent - ' num2str(sum(T_off)./60) ' Hours']})
    %     ylabel('Fraction of epoches');xlabel('Time [min]');
    %     xlim([0, 100]);ylim([0 1])
    %     subplot(1,3,3)
    %     h=histogram(activity_dur,'Normalization','probability');
    %     title({'Combined', ['total simulation time - ' num2str((sum(T_on)+sum(T_off))./60) ' Hours']})
    %     ylabel('Fraction of epoches');xlabel('Time [min]');
    %     xlim([0, 100]);ylim([0 1])
    
    % figure;
    % bar(states,mean_op)
    %
    % hold on
    %
    % er = errorbar(states,mean_op,std_op,std_op);
    % er.Color = [0 0 0];
    % er.LineStyle = 'none';
       save(['Stochastich cooperative new, A conc. = ' num2str(act)])
end
% figure;
% subplot(1,2,1)
% plot(A,theta_on)
% title('Time constant as afunction of activator concentration');
% xlabel('Activator conentration [nM]');ylabel('Time [min]');
% hold on
% tube = fill([A,fliplr(A)], [(theta_on-sig_on),fliplr(theta_on+sig_on)], [0 1 0]);%plotting a tube of the error around the values
% set(tube, 'edgecolor', 'none');
% set(tube, 'FaceAlpha', 0.3);
% hold off
% subplot(1,2,2)
% plot(A,theta_off)
% title('Time constant as afunction of activator concentration');
% xlabel('Activator conentration [nM]');ylabel('Time [min]');
% hold on
% tube = fill([A,fliplr(A)], [theta_off-sig_off,fliplr(theta_off+sig_off)], [0 1 0]);
% set(tube, 'edgecolor', 'none');
% set(tube, 'FaceAlpha', 0.3);
% hold off