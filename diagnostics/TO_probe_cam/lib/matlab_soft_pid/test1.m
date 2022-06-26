%testing the controler against a simulated plant-actuator model


pidstate=[];
pidstate.initalize=1;
pidstate.ctr_output=0;
pidstate.setpt=1;
pidstate.verbose=0;
pidstate.k_int=-1e1;
pidstate.k_prop=-2e1;
pidstate.outlims=[-1 1]*1e5;
pidstate.aw_thresh_range=0.01; %how far away from the edge AW starts 
pidstate.int_lim=500;        %limits on the integerator term
pidstate.slew_lim=1e6;    %

sim.run.time=0;
sim.run.itt=0;
sim.set.tmax=10;
sim.set.dt=1e-2;
sim.set.pid_poll_time=1e-2;
sim.set.pid_calltime=1e-3;
pidstate.time=sim.run.time; %only for simulation for realtime(ish) do not set

sim.run.plant_temp=2;
sim.set.plant_heater_strength=1;
sim.set.plant_thermal_mass=1;
sim.run.plant_load=0;
sim.set.plant_drift=1e-2;

iimax=ceil(sim.set.tmax/sim.set.dt);
sim.history.plant_temp=nan(1,iimax);
sim.history.plant_load=nan(1,iimax);
sim.history.time=nan(1,iimax);
sim.history.ctr=nan(1,iimax);


fprintf('itt %06i:%06i',iimax,0)
while sim.run.itt<iimax
    if mod(sim.run.itt,100)==0; fprintf('\b\b\b\b\b\b%06i',sim.run.itt);end
    sim.run.time=sim.run.time+sim.set.dt;
    sim.run.itt=sim.run.itt+1;
    if sim.run.time>pidstate.time+sim.set.pid_poll_time
        pidstate.meas=sim.run.plant_temp;
        pidstate.set_time=sim.run.time;
        pidstate=pid_loop(pidstate);
        pidstate.time=sim.run.time+sim.set.pid_calltime;
    end
    sim.run.plant_load=sim.run.plant_load+sim.set.plant_drift*(rand(1)-0.5)*2;
    sim.run.plant_temp=sim.run.plant_temp+(1/sim.set.plant_thermal_mass)*sim.set.dt*...
        (sim.set.plant_heater_strength*pidstate.ctr_output-sim.run.plant_load);

    sim.history.plant_temp(sim.run.itt)=sim.run.plant_temp;
    sim.history.plant_load(sim.run.itt)=sim.run.plant_load;
    sim.history.ctr(sim.run.itt)=pidstate.ctr_output;
    sim.history.time(sim.run.itt)=sim.run.time;

end
fprintf('\n')


%%
figure(1)
subplot(3,1,1)
plot(sim.history.time,sim.history.plant_temp)
subplot(3,1,2)
plot(sim.history.time,sim.history.ctr)
subplot(3,1,3)
plot(sim.history.time,sim.history.plant_load)


%%

fprintf('itt %06i:%06i',iimax,0)
fprintf('\b\b\b\b\b\b%06i',1)
