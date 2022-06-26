function pidstate=pid_loop(pidstate)

%pidstate - A software PID loop with good integerator windup protection and
%slew rate limits.
%
% This function is intended for low bandwidth computer-in-the-loop control.
% If you are aiming for maximum performance I recomend either an anlog
% control loop or a FPGA. If you want better than pid performance look
% elsewhere at some adaptive controllers.
% 
%
% Syntax:  [output1,output2] = function_name(input1,input2,input3)
%
% Inputs:
%   pidstate.meas 
%   pidstate.setpt          - desired setpoint
%   pidstate.k_int=         -integerator gain in [out]/(sec*[pidstate.erroror])
%   pidstate.k_prop=        -proportional gain with unity output scaling
%   pidstate.outlims        -output limits
%   pidstate.aw_thresh_range -how far away from the edge aw starts (full range 0-1)
%   pidstate.int_lim        -integerator limits 
%   pidstate.slew_lim       -output slew limit in multiples of output range/sec

%
% Outputs:
%    pidstate.ctr_output - float,actuator control  between pidstate.outlims
%    pidstate.integrator     - double used to initalize integerator;
% 
% Example: 
% pidstate.ctr_output=0;
% pidstate.setpt=0.65; %max val set pt
% pidstate.k_int=-2e-6;
% pidstate.k_prop=1e-8;
% pidstate.outlims=[-1 1];
% pidstate.aw_thresh_range=0.05; %how far away from the edge AW starts 
% pidstate.int_lim=3;        %limits on the integerator term
% pidstate.slew_lim=1e-5;    %slew rate limit


% Bugs,Improvements,Ideas
% input external timing
% output is not continious with proportional control gain change
% Derivative term not implemented
% this would probably be better written as a class
% slew lim always activates on first loop

% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
% See also: none

% Author: Bryce Henson
% Work address
% email: bryce.henson@live.com
% Last revision: 2018-08-23


%------------- BEGIN CODE --------------
if ~isfield(pidstate,'initalize')
    warning('call pid_loop with the pidstate.initalize=true to set things up')
    pidstate.initalize=true;
end

if pidstate.initalize
    if ~isfield(pidstate,'outlims')
       error('requires output limits')
    end 
    if ~isfield(pidstate,'verbose')
        pidstate.verbose=0;
    end   
    if ~isfield(pidstate,'ctr_output')
        warning('no inital output specified using middle of output limits')
        pidstate.ctr_output=mean(pidstate.outlims);
    end
    if ~isfield(pidstate,'slew_lim')
        warning('no slew rate limit specified setting as inf')
        pidstate.slew_lim=inf;
    end
    if ~isfield(pidstate,'dout_lim')
        warning('no output delta specified setting as inf')
        pidstate.dout_lim=inf;
    end
    if ~isfield(pidstate,'int_lim')
        warning('no integerator limits specified, setting as twice the output range')
        pidstate.int_lim=mean(pidstate.outlims)+[-1,1]*2*range(pidstate.outlims);
        fprintf('int_lim=[%s]\n',sprintf('%d ', pidstate.int_lim))
    end
    if ~isfield(pidstate,'aw_thresh_range')
        %default to 1e-3 of the output range
        pidstate.aw_thresh_range=range(pidstate.outlims)*1e-3;
    end
    if ~isfield(pidstate,'time')
        pidstate.time=posixtime(datetime('now'));
    end
    
    pidstate.aw=1;
    pidstate.meas=pidstate.setpt;
    pidstate.integrator=pidstate.ctr_output;
    pidstate.ctr_prev=pidstate.ctr_output;
    pidstate.initalize=false;

end

old_time=pidstate.time;
if isfield(pidstate,'set_time')
    time_ctr_start=pidstate.set_time;
else
    time_ctr_start=posixtime(datetime('now'));
end
pidstate.loop_time=(time_ctr_start-old_time);
if pidstate.verbose>2,fprintf('feedback bandwidth %.2f\n',1/pidstate.loop_time), end

pidstate.ctr_prev=pidstate.ctr_output;

logist=@(x) 1./(1+exp(-x));
%scale the anti windup range to a unit output range
scaled_aw_range=pidstate.aw_thresh_range/range(pidstate.outlims);
%needs to be defined for di==0
aw_fun_range=@(x,y) (logist((x-scaled_aw_range)*10/scaled_aw_range))*(y<=0)...
    +(1-logist((x-1+scaled_aw_range)*10/scaled_aw_range))*(y>0);

% xvals=linspace(0,1,1e4);
% plot(xvals,aw_fun_range(xvals,1))
% hold on
% plot(xvals,aw_fun_range(xvals,-1))
% hold off

if pidstate.verbose>1
    fprintf('previous control output %.3f\n',pidstate.ctr_prev)
end

pidstate.error=pidstate.meas-pidstate.setpt;
if pidstate.verbose>1,fprintf('error %.3f\n',pidstate.error), end
di=pidstate.k_int*pidstate.loop_time*pidstate.error;
scaled_last_control=(pidstate.ctr_prev-pidstate.outlims(1))/range(pidstate.outlims);
if pidstate.verbose>1,fprintf('scaled last control %.3f\n',scaled_last_control), end
pidstate.aw=aw_fun_range(scaled_last_control,di);
if pidstate.verbose>1, fprintf('anti windup %.3f\n',pidstate.aw), end
di=di*pidstate.aw; %actuator range anti windup
pidstate.integrator=pidstate.integrator+di;%+aw_slew*k_aw_rate;
pidstate.integrator=min(max(pidstate.int_lim(1),pidstate.integrator),pidstate.int_lim(2));
pidstate.ctr_output=pidstate.integrator+pidstate.k_prop*pidstate.error;%-aw_rate*k_aw_rate);

pidstate.delt_out=(pidstate.ctr_prev-pidstate.ctr_output);
pidstate.slew=pidstate.delt_out/(pidstate.loop_time);
if pidstate.verbose>1,fprintf('desired slew %.3f\n',pidstate.slew), end


%combine the delta output step limits and slew limits as an effective slew limit
combined_delt_lim=min(pidstate.slew_lim*pidstate.loop_time,pidstate.dout_lim);
%fprintf('slew before corr %f\n',pidstate.slew)
if abs(pidstate.delt_out)>combined_delt_lim
    %fprintf('pre mod int %f\n',pidstate.integrator)
    pidstate.integrator=pidstate.ctr_prev-combined_delt_lim*sign(pidstate.slew)-pidstate.k_prop*pidstate.error;
    %fprintf('slew mod int %f\n',pidstate.integrator)
    pidstate.aw_slew=1; %not a anti windup should just call slew lim
else
    pidstate.aw_slew=0;
end

%recalc output
%integrator=integrator+aw_slew;
%pidstate.integrator=min(max(-pidstate.int_lim,pidstate.integrator),pidstate.int_lim);
pidstate.ctr_output=pidstate.integrator+pidstate.k_prop*pidstate.error;%-aw_rate*k_aw_rate);
pidstate.ctr_output=min(max(pidstate.outlims(1),pidstate.ctr_output),pidstate.outlims(2));
pidstate.slew=(pidstate.ctr_prev-pidstate.ctr_output)/(pidstate.loop_time);
if pidstate.verbose>1 && pidstate.aw_slew
    fprintf('clamped slew %.3f\n',pidstate.slew)
end
if pidstate.verbose>1 &&pidstate.aw_slew
    fprintf('output %.3f\n',pidstate.ctr_output)
end
%fprintf('slew after corr %f\n',pidstate.slew)


pidstate.time=posixtime(datetime('now'));



if pidstate.verbose>1
    fprintf('loop time %.3f\n',pidstate.loop_time)
end
end