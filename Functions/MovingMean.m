function [mm]=MovingMean(A,wl,t_step)
    %A is the vector of protein values in time
    %wl is the window length which the mean should be taken
    %t_step is the current time step (only calculates one time point at a
    %time (use for loop))
    
    if mod(wl,2)==1
        error('wl must be even')
    end
    
    elems=size(A,2);
    st=t_step-wl/2;
    en=t_step+wl/2;
    
    if st<1
        st=1;
    end
    
    if en>elems
        en=elems;
    end
    
    mm=mean(A(:,st:en),2);   
    
end


%%Old function
% function [mm]=MovingMean(A,wl,t_step)
%     %A is the vector of protein values in time
%     %wl is the window length which the mean should be taken
%     
%     if t_step<=wl
%         mm=mean(A(:,1:t_step),2); %Moving mean value
%     else
%         mm=mean(A(:,t_step-wl:t_step),2);
%     end
%     
% end