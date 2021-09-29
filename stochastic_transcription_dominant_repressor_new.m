A=[5];%act concentration in nanoMolar
R=[0.1 1 2 5 10 25 100];%rep concentration in nanoMolar
theta_on=[];sig_on=[];
theta_off=[];sig_off=[];
On=[];
omega=1.5;%the factor by which a repressor multiplies the off rate

for rep=R
for act=A
    
    N=6;%total number of binding sites
    states=0:N;%available binding states in the system
    kon=0.008;
    kon_a=kon*act;%activator's binding rate in Hz*nM.
    kon_r=kon*rep;%repressor binding rate
    koff_active_a=0.004;%for activator in active state. rate of single dissociation, assuming only 1 dna strand. in Hz.
    koff_inactive_a=12;%for activator in inactive state.
    
    koff_active_r=koff_active_a;%for repressor in active state. rate of single dissociation, assuming only 1 dna strand. in Hz.
    koff_inactive_r=koff_inactive_a;%for repressor in inactive state.
    
    
    rate_Mon=0.0017;%rate for transitioning to active state at 0 occupancy
    rate_Moff=0.005;%rate for transitioning to inactive state at 0 occupancy
    m=1;%memory of system's initial condition
    
    window_dur=240;%in sec
    
    sample_size=1;
    Occ_P=zeros(sample_size,length(states));
    
    for itr2=1:sample_size
        
        n_a=0;%initial condition - number of sites currently occupied by activator
        n_r=0;%number of sites currently occupied by repressor
        t=0;%initial condition - time in sec
        occ_a=n_a;%occupancy vector - for convenience
        occ_r=n_r;
        M=m;
        Mtimes=[M];%timing the states
        while t(end)<360000
            n=n_a+n_r;%total occupany
            current_Mon=rate_Mon;%M appropriate to the occupancy. dependant on total occupancy
            current_Moff=(rate_Moff/(n+1))*(omega*n_r+1);%M appropriate to the occupancy. dependant negatively on total occupancy and positivly on rep occ
            k_add_a=(N-n)*kon_a;%rate constant for adding an activator.
            k_add_r=(N-n)*kon_r;%rate constant for adding a repressor.
            
            if M==1
                k_dec_a=n_a*koff_active_a;%rate constant for activator dissociation
                k_dec_r=n_r*koff_active_r;%rate constant for repressor dissociation
                
                Rtot=k_add_a+k_dec_a+k_add_r+k_dec_r+current_Moff;%Total rate of event
                
                p_add_a=k_add_a/Rtot;%probability of adding an activator
                p_add_r=k_add_r/Rtot;%probability of adding a repressor
                p_dec_a=k_dec_a/Rtot;%probability of removing an activator
                p_dec_r=k_dec_r/Rtot;%probability of removing a repressor
                
            elseif M==0
                k_dec_a=n_a*koff_inactive_a;%rate constant for activator dissociation
                k_dec_r=n_r*koff_inactive_r;%rate constant for repressor dissociation
                
                Rtot=k_add_a+k_dec_a+k_add_r+k_dec_r+current_Mon;%Total rate of event
                
                p_add_a=k_add_a/Rtot;%probability of adding an activator
                p_add_r=k_add_r/Rtot;%probability of adding a repressor
                p_dec_a=k_dec_a/Rtot;%probability of removing an activator
                p_dec_r=k_dec_r/Rtot;%probability of removing a repressor
            end
            
            num=rand;%rand chooses a number between 0-1 from a uniform distribution
            
            if num<=p_add_a %increasing activator
                n_a=n_a+1;
            elseif num>p_add_a && num<=(p_add_a+p_add_r)%increasing repressor
                n_r=n_r+1;
            elseif num>(p_add_a+p_add_r) && num<=(p_add_a+p_add_r+p_dec_a)%decrease activator
                n_a=n_a-1;
               elseif num>(p_add_a+p_add_r+p_dec_a) && num<=(p_add_a+p_add_r+p_dec_a+p_dec_r)%decrease repressor 
                   n_r=n_r-1;
                elseif num>(p_add_a+p_add_r+p_dec_a+p_dec_r)
                M=double(~M);%changing mode
            end
            
            t=[t, t(end) + exprnd(1/Rtot)];%creating a vec with updated times
            Mtimes=[Mtimes, M];
            occ_a=[occ_a, n_a];%adding the current occupancy to the occupancy vec
            occ_r=[occ_r, n_r];
        end
%         OCC=movmean_time(occ,t,window_dur); 
        occ=occ_a+occ_r;%total promoter's occupancy
        
%         figure;
%          subplot(1,3,1)
%         plot(t./60,occ_a)%show time in Minutes
%         title('Activator occupancy as a function of time')
%         xlabel('Time [min]');ylabel('Occupancy')
%         
% %         figure;
%          subplot(1,3,2)
%         plot(t./60,occ_r)%show time in minutes
%         title('Repressor occupancy as a function of time')
%         xlabel('Time [min]');ylabel('Occupancy')
%         
%        
%         subplot(1,3,3)
%         plot(t./60,occ)%show time in minutes
%         title('Total occupancy as a function of time')
%         xlabel('Time [min]');ylabel('Occupancy')
        
%             figure;
%             plot(t./60,movmean(occ,10))
%             hold on
%             plot(t./60, Mtimes)
        
        
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
            occ_prob(itr+1)=sum(dur(occ_a(1:end-1)==itr))/t(end);%caculating probability of oocurance by temporal duration of a state
        end
        
%                     figure;
%                     bar(states,occ_prob)
%                     title('Occurence of occupancy')
%                     xlabel('Occupancy');ylabel('Occurence')
%                     ylim([0 1])
%         
        
        Occ_P(itr2,:)=occ_prob;%concatanating runs
        
        
        th=3;%threshold for active transcription
        occOn_loc=occ_a>=th;%indices for active-transcription momments
        occ_ind=logical([abs(diff([occOn_loc, 2]))]);%indicating locations where the system's transcriptional activity changes.
        %concatanating 2 to make sure that the last durstion is accounted for
        t_activity=t(occ_ind);%the times at which the system's transcriptional activity changed
        activity_dur=diff([0 t_activity]);%duration vec for the activity. combining on and off states. cocncatanating 0 to fix diff's offset
        activity_dur=activity_dur./60;%switching units to minutes
        
        activity_ind=1:2:length(activity_dur);
        
        if occ_a(1)<3
            if activity_ind(end)==length(activity_dur)
                
                T_off=activity_dur(activity_ind);%noticing Transcription-off durations
                activity_ind(1)=[];
                T_on=activity_dur(activity_ind-1);%noticing Transcription-on durations
            else
                T_off=activity_dur(activity_ind);%noticing Transcription-off durations
                T_on=activity_dur(activity_ind+1);%noticing Transcription-on durations
            end
        elseif occ_a(1)>=3
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
      save(['Stochastich T, dominant rep new, A conc. = ' num2str(act)  ' , R coc. = ' num2str(rep)])
end
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