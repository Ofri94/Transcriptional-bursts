N=3;%number of SPS sites
states=0:2*N;%available binding states in the system

loc=zeros(2,N);%creating location index for occupancy
% n=randi(2*N,1);%start with a random occupancy
n=0;

loc(1:n)=1;%setting initial condition of location

A=[0.15:0.1:1.05];%act concentration in nanoMolar

theta_on=[];sig_on=[];
theta_off=[];sig_off=[];

C=20;%cooperativity factor


for act=A
    kon=0.008*act;%in Hz*nM.
    koff_active=0.004;%for active state. rate of single dissociation, assuming only 1 dna strand. in Hz.
    koff_inactive=12;%for active state.
    
    rate_Mon=0.0017;%rate for transitioning to active state at 0 occupancy
    rate_Moff=0.005;%rate for transitioning to inactive state at 0 occupancy
    m=1;%memory of system's initial condition
    
    
    sample_size=1;
    Occ_P=zeros(sample_size,length(states));
    
    for itr2=1:sample_size%for statistical analysis
        
        t=0;%initial condition - time in sec
        occ=n;%occupancy vector - for convenience;
        M=m;
        Mtimes=[M];%timing the states
        while t(end)<360000
            
            current_Mon=rate_Mon;%M appropriate to the occupancy
            current_Moff=rate_Moff/(n+1);%M appropriate to the occupancy
            
            k_add=(2*N-n)*kon;%rate constant for adding a TF. depends only on activator concentration
            
            full_sps=find(sum(loc)==2);%location of full sps sites
            
            
            if M==1
                k_dec_star=(2*length(full_sps)*koff_active)/C;%rate constant for TF dissociation when cooperatively bound
                k_dec=(n-2*length(full_sps))*koff_active;%rate constant for TF dissociation when bound non-cooperativly
                
                Rtot=k_add+k_dec_star+k_dec+current_Moff;%Total rate of event
                
                p_add=k_add/Rtot;%probability of adding a TF
                p_dec_star=k_dec_star/Rtot;%probability of removing a TF which is cooperativly bound
                p_dec=k_dec/Rtot;%probability of removing a TF which is non-cooperativly bound
            elseif M==0
                k_dec_star=2*length(full_sps)*koff_inactive/C;%rate constant for TF dissociation when cooperatively bound
                k_dec=(n-2*length(full_sps))*koff_inactive;%rate constant for TF dissociation when bound non-cooperativly
                
                Rtot=k_add+k_dec_star+k_dec+current_Mon;%Total rate of event
                
                p_add=k_add/Rtot;%probability of adding a TF
                p_dec_star=k_dec_star/Rtot;%probability of removing a TF which is cooperativly bound
                p_dec=k_dec/Rtot;%probability of removing a TF which is non-cooperativly bound
            end
            
            num=rand;%rand chooses a number between 0-1 from a uniform distribution
            
            if num<=p_add %updating occupancy (creating a vec).
                n=n+1;
                f=find(loc==0);%locations of empty sites
                loc(f(randi(length(f),1)))=1;%adding a TF randomly
                
            elseif num>p_add && num<=(p_add+p_dec_star)
                n=n-1;
                
                chosen=full_sps(randi(length(full_sps),1));%chosen sps site to take TF from
                loc(1,chosen)=0;%removing TF from SPS site
                
            elseif num>(p_add+p_dec_star) && num<=(p_add+p_dec_star+p_dec)
                n=n-1;
                half_sps=find(sum(loc)==1);%location of half full sps sites
                
                loc(:,half_sps(randi(length(half_sps),1)))=0;%removing a TF from a random half-full site
                
            elseif num>(p_add+p_dec_star+p_dec)
                M=double(~M);%changing mode
            end
            
            t=[t, t(end) + exprnd(1/Rtot)];%creating a vec with updated times
            Mtimes=[Mtimes, M];
            occ=[occ, n];%adding the current occupancy to the occupancy vec
        end
        
    end
    dur=diff(t);%discrete duartaion of modes
    occ_prob=zeros(1,N+1);%allocating space
    for itr=states
        occ_prob(itr+1)=sum(dur(occ(1:end-1)==itr))/t(end);%caculating probability of oocurance from the temporal duration of a state
    end
    
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
    
 save(['Stochastich cooperativity analysis new, A conc. = ' num2str(act)])
end

% figure;
% plot(t./60,occ)%show time in Minutes
% title('Transcription rate as a function of time')
% xlabel('Time [min]');ylabel('Transcription rate [a.u.]')
% 
% figure;
% bar(states,occ_prob)
% title('Occurence of occupancy')
% xlabel('Occupancy');ylabel('Occurence')
% ylim([0 1])