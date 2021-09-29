function [smoothed] = movmean_time(occ,t,window_dur)
%this function applies a moving average filter over time rather than a fixed
%window of data-points.
%Importantly!! the time window is only ahead of the point. 
%more specifically - the nth point will be averaged over the n-->n+loc dots,
%where loc is the location of the final point to average across
%in gillespi simulations, the passage of time
%between two points of "measurements" is not fixed, therefore this function
%will average the signal over time, rather than time-points
%occ is the occupancy vector
%t is the time vector
%window_dur is a scalar stating the window length is units matching t
smoothed=[];%output variable
t_diff=diff(t);

for ind=1:length(t_diff)
  window=t_diff(ind); 
  loc=ind+1;
while window<window_dur
    if loc<length(occ)
        window=window+t_diff(loc);
        loc=loc+1;
    else
        break
    end
end
current=occ(ind:(loc-1));
smoothed=[smoothed, mean(current)];
end

smoothed=[smoothed, occ(end)];
end
