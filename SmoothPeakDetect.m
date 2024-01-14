% Detect whether a given array at the peak index is a locally smooth peak
% (i.e. monotonic decay and increase on the left and right side
% respectively, within the checking_range), can also handle detecting max
% (i.e. leftr mono inc, right mono dec)
function is_smooth_peak=SmoothPeakDetect(array,max_or_min,peak_idx,checking_range)

%calculate differences between element up to the end
left_df=diff(array(max(1,peak_idx-checking_range):peak_idx));
right_df=diff(array(peak_idx:min(length(array),peak_idx+checking_range)));

switch max_or_min
    case 'max'
        if all(left_df>=0) && all(right_df<=0)
            is_smooth_peak=1;
        else
            is_smooth_peak=0;
        end
    case 'min'
        if all(left_df<=0) && all(right_df>=0)
            is_smooth_peak=1;
        else
            is_smooth_peak=0;
        end
end


end