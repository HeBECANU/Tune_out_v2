function motion_solution=motion_in_harmonic_w_gaussian(opts)
% simulate the motion of a classical particle in a (an)harmonic trap combined with a guasian potential
% used to calculate the nonlinearity of the polarizability measurment method used in the tune out
% and the optimal oscillation amp/ beam size
% required inputs
%   opts.trap_derivs or opts.trap_freq
%   opts.masss
% optional inputs
%   opts.gauss_amp
%   opts.gauss_radius % beam radius
%   opts.damping_coef (aka viscous damping coefficient) or damping_time (assuming harmonic osc.)


verbose=1;


if isfield(opts,'trap_freq')
    if isfield(opts,'trap_derivs')
        error('if omega is specified cannot use trap derivs')
    end
    opts.trap_derivs=[0,0,opts.mass*(2*pi*opts.trap_freq)^2];
end
% calculate for the ode step size
trap_omega= sqrt(opts.trap_derivs(3)/opts.mass);


if ~isfield(opts,'damping_time') && ~isfield(opts,'damping_coef')
    opts.damping_coef =0;
elseif isfield(opts,'damping_time') && isfield(opts,'damping_coef')
    error('input damping_time or damping_coef')
elseif isfield(opts,'damping_time') && ~isfield(opts,'damping_coef')
    if verbose>2 && (...
                    numel( opts.trap_derivs)>3 || ...
                    (isfield(opts,'gauss_amp') && ~isnan(opts.gauss_amp) && abs(opts.gauss_amp)>0) ...
                    )
        warning('damping coef calculation from damping time will not be correct with anharmonic terms')
    end
    damping_ratio=1/(trap_omega*opts.damping_time); % damping ratio is usually denoted as \zeta
    % damp_coef=zeta*2*sqrt( m* k(spring_const) )
    % omega=sqrt(k/m)
    % k= omega^2*m
    % damp_coef=zeta*2*sqrt( m^2* omega^2 )
    % damp_coef=zeta*2*m* omega
    opts.damping_coef=damping_ratio*2*trap_omega*opts.mass;
end


% set up the params for the state update function
sim_st=[];
sim_st.trap_derivs=opts.trap_derivs;
sim_st.mass=opts.mass;
sim_st.damping_coef =opts.damping_coef ;
if isfield(opts,'gauss_amp') && ~isnan(opts.gauss_amp)
     sim_st.gauss_amp=opts.gauss_amp;
     sim_st.gauss_radius=opts.gauss_radius;
     if ~isfield(opts,'gauss_displacement')
         sim_st.gauss_displacement=0;
     else
        sim_st.gauss_displacement=opts.gauss_displacement;
     end
else
    sim_st.gauss_amp=nan;
end
    
% if the gaussian is displaced we will find the new trap minima pos
% before simulation of the ocillation
initial_state=opts.inital_state;
if isfield(opts,'find_new_min') && opts.find_new_min ...
        && ~isnan(sim_st.gauss_amp) && abs(sim_st.gauss_amp)>0 ...
        && abs(sim_st.gauss_displacement)>0 ...

    %%
    % for testing, use a random amplidude and displacemnt
%     sim_st.gauss_displacement=(rand(1)-0.5)*4*sim_st.gauss_radius;
%     sim_st.gauss_amp=1e-29*(rand(1)-0.5)*10
    %
    % this is very rough but is good enough to start fminsearch on
    start_min_pt=-sign(sim_st.gauss_amp)*sim_st.gauss_displacement*3; %
    
    fminopt = optimset('TolX',1e-11,'TolFun',abs(sim_st.gauss_amp)/1000);                                          
    trap_min=fminsearch(@(x) potential_val(x,sim_st),start_min_pt,fminopt);                                 
    %trap_min=fzero(@(x) potential_force(x,sim_st) ,start_min_pt);
    if isempty(trap_min)
        error('fminbnd did not return a value')
    end
    
    if isfield(opts,'plot_find_new_min') && opts.plot_find_new_min
        %%
        global const
        stfig('finding new trap minimum');
        clf
        xscale=1e6;
        yscale_u=1e6/const.kb;
        yscale_f=1e3/const.kb;
        subplot(1,2,1)
        xsamp=linspace(-1,1,1e3)*4*sim_st.gauss_radius;
        xsamp=col_vec(xsamp);
        fsamp=potential_force(xsamp,sim_st);
        plot(xsamp*xscale,fsamp*yscale_f)
        hold on
        plot(trap_min*xscale,potential_force(trap_min,sim_st)*yscale_f,'bx');
        hold off
        xlabel('position ($\mu$m)')
        ylabel('force/$k_{b}$ (pK/$\mu$m)')
        subplot(1,2,2)
        hold on
        plot(xsamp*xscale,potential_val(xsamp,sim_st)*yscale_u)
        plot(start_min_pt*xscale,potential_val(start_min_pt,sim_st)*yscale_u,'rx')
        plot(trap_min*xscale,potential_val(trap_min,sim_st)*yscale_u,'bx')
        hold off
        xlabel('position ($\mu$m)')
        ylabel('potential/$k_{b}$ ($\mu$K)')
        
        txt = {sprintf('minima pos %.3f $\\mu$m',trap_min*1e6),...
               sprintf('beam pot %.3f pK',1e12*sim_st.gauss_amp/const.kb),...
               sprintf('beam pos %.3f $\\mu$m',sim_st.gauss_displacement*1e6),...
               sprintf('beam radius %.3f $\\mu$m',sim_st.gauss_radius*1e6)};
        text(0.5,0.5,txt,'Units','normalized','HorizontalAlignment','center')
        pause(1e-6)
    end
    initial_state(1)=trap_min;

end


% 
ode_opts = odeset('RelTol',1e-9,'AbsTol',1e-12,'Stats','off','InitialStep',1e-3/trap_omega,'MaxStep',0.1/trap_omega,'Vectorized',1); %
motion_solution = ode45(@(t,X) particle_cord_deriv(t,X,sim_st),opts.tlims,opts.inital_state,ode_opts);
%ode113



end



function dx_dt=particle_cord_deriv(t,pos_vel,sim_st)
% the particle has 2 cordinates the fist is position and the second is the velocity
% required arguements
% sim_st.sim_st
% sim_st.mass
% sim_st.gauss_sigma
% sim_st.gauss_amp
% sim_st.damping_coef

dx_dt=nan*pos_vel;

dx_dt(1,:)=pos_vel(2,:);
pforce=potential_force(pos_vel(1,:),sim_st);

force=pforce- sim_st.damping_coef*pos_vel(2,:);
dx_dt(2,:) = force./sim_st.mass;

end

% the force from the potential is calculated in a sep function so we can acess it in the main function
% to find the new trap minimum with a displaced gaussian
function pforce=potential_force(pos,sim_st)
%pos=col_vec(pos);
pot_grad=deriv_taylor_series(pos,sim_st.trap_derivs,0,1);
if ~isnan(sim_st.gauss_amp) && sim_st.gauss_amp~=0 
    %pot_grad=pot_grad+...
    %    gaussian_function_1d(pos_vel(1,:),sim_st.gauss_sigma,0,sim_st.gauss_amp,0,'norm','amp','derivative',1);
    % the above is too slow because gaussian_function_1d has a parse on the input
    % we defie our own low level version
    pot_grad=pot_grad+...
        gaussian_1_deriv(pos,...
                        sim_st.gauss_radius/2,...
                        sim_st.gauss_displacement,...
                        sim_st.gauss_amp);
end
    
pforce=-pot_grad;
end

function pval=potential_val(pos,sim_st)
%pos=col_vec(pos);
pval=deriv_taylor_series(pos,sim_st.trap_derivs,0,0);
if ~isnan(sim_st.gauss_amp) && sim_st.gauss_amp~=0 
    %pot_grad=pot_grad+...
    %    gaussian_function_1d(pos_vel(1,:),sim_st.gauss_sigma,0,sim_st.gauss_amp,0,'norm','amp','derivative',1);
    % the above is too slow because gaussian_function_1d has a parse on the input
    % we defie our own low level version
    pval=pval+...
        gaussian_val(pos,...
                        sim_st.gauss_radius/2,...
                        sim_st.gauss_displacement,...
                        sim_st.gauss_amp);
end
    
end


function y=gaussian_1_deriv(x,sigma,mu,amp)
% take the first derivative of a gaussian
% used because gaussian_function_1d takes too long to parse inputs
 y = -amp*( (x-mu)./ (sigma.^2) ).*exp(-(1/2)*((x-mu)./sigma).^2);
end


function y=gaussian_val(x,sigma,mu,amp)
% take the first derivative of a gaussian
% used because gaussian_function_1d takes too long to parse inputs
 y = amp.*exp(-(1/2)*((x-mu)./sigma).^2);
end
