function vzxy_out=txy_to_vel(txy_in,out_time,gravity,fall_distance)
% gravity is assumed to be in -ve z

if size(txy_in,2)~=3
    error('wrong txy size')
end 
if ~isscalar(out_time)
    error('time_offset should be scalar')
end 
if ~isscalar(gravity)
    error('gravity should be scalar')
end 

gravity=-abs(gravity);


%convert the position data into velocity

% the z/time axis is slightly complicated because
% z=z(0) + vz(0)t+1/2 g t^2
% want to solve for v(0) given x(t_fall) the fall distance and gravity(in -ve x) , say bec pos is x=0 at t=0
% and is released at t=0
% -fall_dist=vz(0)t+1/2 g t^2
% vz(0)=-fall_dist/t  -1/2 g t 

% lets compare this to what i have calculated in the past
% fall_time=sqrt(2*fall_dist/g)
% vel_fall=sqrt(2*fall_dist*g)
% vz(0)~(t-fall_time)*vel_fall/fall_time
% vz(0)~(t-sqrt(2*fall_dist/g))*sqrt(2*fall_dist*g)/sqrt(2*fall_dist/g)
% vz(0)~(t-sqrt(2*fall_dist/g))*g
% vz(0)~ g t -sqrt(2*g*fall_dist)

% initalize the output array
vzxy_out=txy_in*nan;
% now we dont have the fall time directly, we can use the offset from the output time (generaly atom laser
% pulse time)

fall_time=txy_in(:,1)-out_time;

vel_z=-fall_distance./fall_time-(1/2)*gravity*fall_time;

% then we want to convert the x,y data into velocity using the fall time
% bacause the TOF changes a little for each count we can correct for this

vzxy_out(:,[2,3])=txy_in(:,[2,3])./repmat(fall_time,1,2);
vzxy_out(:,1)=vel_z;



end