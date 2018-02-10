%%
%% Optimization based code for flying Pelican
%% The optimization satisfied the velocity constraints
%% It tries to minimize the error of cost function
%%

%%
%% First go to terminal and run following
%% cd rotors_ws
%% source devel/setup.bash
%% roslaunch rotors_gazebo mav_hovering_example.launch mav_name:=pelican
%% Only after running that, you should run this code.
%% Don't forget to bring a sphere into the gazebo simulation.
%% Use gazebo top bar to get a sphere, on the left column, 
%% go to models, unit_sphere_0 and turn off the gravity factor!
%%

%%
%% The code integrates position controller and uses a velocity controller.
%% The code assumes that the starting position of drone is (0,0,0)
%% The updated version will be more generic in nature!
%% This code uses Bhattacharyya distance as a metric to avoid Obstalcle
%% Here, both, the obstacle as well as the robot(drone) have gaussian uncertainties
%%

clc;
clear all;
close all;


%% Drone is supposed to reach (dest_x, dest_y, dest_z) from (0,0,0) in
%% 'n' timesteps with each timestep of duration del_t
%% Version 2 (for higher values of n)
n = 30;      % Number of timesteps
del_t = 1.0; % Time duration for which some command is executed.
tot_time = n*del_t; %% Total time of flight

%% defining drone uncertainty
drone_sig = [0.8, 0.0, 0.0; 0.0, 0.7, 0.0; 0.0, 0.0, 0.8];

%% Defining obstacle trajectory
%% For now we will assume that obstacle is moving along a straing line. it's starting from (0,0,0) and It's trajectory will be determined based on the destination of the drone

%% Drone Destination variables
dest_x = 15;
dest_y = 15;
dest_z = 15;

%% Getting linearization points
%% These are the points around which the obstacle avoidance constraints will be linearized
for i = 1:n
   x(i) = 50; 
   y(i) = 50;
   z(i) = 50;
%    x(i) = dest_x*(i-1)/n;
%    y(i) = dest_y*(i-1)/n;
%    z(i) = dest_z*(i-1)/n;
end


%%%%%%%%%%%%%% Obstacle stuff starts here!!!! %%%%%%%%%%%%%%%%%
%% Here, the drone is moving from (0,0,0) to (dest_x, dest_y, dest_z) in n timesteps
%% The obstacle will reach from start_obs to end_obs in n timesteps 
Ob_start=[0, dest_y/2, dest_z];
Ob_end=[dest_x, dest_y/2, 0];
%Ob_start=[0,3,0];
%Ob_end=[dest_x, dest_y+5, dest_z];
%% Calculating velocity of the obstacle in all the directions
vel_obs = (Ob_end - Ob_start)/tot_time;

%% Defining Obstacle path location at each timestep
for i = 1:n
  Ox(i) = Ob_start(1) + del_t*i*vel_obs(1);
  Oy(i) = Ob_start(2) + del_t*i*vel_obs(2);
  Oz(i) = Ob_start(3) + del_t*i*vel_obs(3);
  %Ox(i) = dest_x/2;
  %Oy(i) = dest_y/2;
  %Oz(i) = 1;
end


%% Defining obstalce uncertainty
obs_sig = [1.0, 0.0, 0.0; 0.0, 0.4, 0.0; 0.0, 0.0, 0.7];;

%% Criterion threshold
crit_thresh =  0.0001;

%% The minimum confidence by which both, the drone and obstacle should be avoided
confidence = 0.90;

load Linear_approx.mat
%%load Linear_approx_table.mat
%% Finding U1 threshold for the above confidence
P_u1 = Overlap(int16(confidence*length(Conf)))/2;
CP_u1 = 1 - P_u1; %% Cumulative probability of U1
u1_thresh = norminv(CP_u1);

%%%%%%%%%Obstacle stuff ends here %%%%%%%%%%%%%%%

%% Difference between subsequent velocities
del_Vx = 0.02;
del_Vy = 0.02;
del_Vz = 0.02;

%% Initializing linearization point for ts
for i = 1:n
	td(i) = 0.3;
end

l = 0;
no_of_iter = 10;

while l<no_of_iter
    
    cvx_begin
    
    variables Vx(n) Vy(n) Vz(n) t(n) %% Velocities
    %variables Px(n) Py(n) Pz(n) %% Drone locations at each timestep
    
    for i = 1:n
        Px(i)  =  sum(Vx(1:i))*del_t;
        Py(i)  =  sum(Vy(1:i))*del_t;
        Pz(i)  =  sum(Vz(1:i))*del_t;
    end
    
    %% Cost function
    %%minimise((sum(Vx)*del_t - dest_x)^2 + (sum(Vy)*del_t - dest_y)^2 + (sum(Vz)*del_t - dest_z)^2)
    minimise((Px(n) - dest_x)^2 + (Py(n) - dest_y)^2 + (Pz(n) - dest_z)^2)
    
    %%Obstacle avoindace constrains
    %% original constraint is as follows
    %% (Px(i) - Ox(i))^2/(sig_x^2) + (Py(i) - Oy(i))^2/(sig_y^2) + (Pz(i) - Oz(i))^2/(sig_z^2) >= R^2
    %% Since the above constraint is quadratic in nature, we will linearize it around point x(i), y(i), z(i)
    %% f(Px(i),Py(i),Pz(i)) = (Px(i) - Ox(i))^2 + (Py(i) - Oy(i))^2 + (Pz(i) - Oz(i))^2 linearized around x(i), y(i), z(i) will look like following
    %% (x(i) - Ox(i))^2 + (y(i) - Oy(i))^2 + (z(i) - Oz(i))^2 + 2(x(i) - Ox(i))*(Px(i) -  Ox(i)) + 2(y(i) - Oy(i))*(Py(i) -  Oy(i)) + 2(y(i) - Oy(i))*(Py(i) -  Oy(i))
    %% Constraint will be (x(i) - Ox(i))^2 + (y(i) - Oy(i))^2 + (z(i) - Oz(i))^2 + 2(x(i) - Ox(i))*(Px(i) -  Ox(i)) + 2(y(i) - Oy(i))*(Py(i) -  Oy(i)) + 2(z(i) - Oz(i))*(Pz(i) -  Oz(i)) > R^2
    if(l~=0)
    	for i = 1:n
    	   %% Linearized criteria for minimizing overlap through linear approximation
    	   -crit_thresh <= linconstraintconstFunc(obs_sig, drone_sig, Ox(i), Oy(i), Oz(i), x(i), y(i), z(i), td(i)) +  ...
    	   					(Px(i) - x(i))*linconstraintPxFunc(obs_sig, drone_sig, Ox(i), Oy(i), Oz(i), x(i), y(i), z(i), td(i)) +  ...
    	   						(Py(i) - y(i))*linconstraintPyFunc(obs_sig, drone_sig, Ox(i), Oy(i), Oz(i), x(i), y(i), z(i), td(i)) +  ...
    	   							(Pz(i) - z(i))*linconstraintPzFunc(obs_sig, drone_sig, Ox(i), Oy(i), Oz(i), x(i), y(i), z(i), td(i)) +  ...
    	   								(t(i) - td(i))*linconstrainttFunc(obs_sig, drone_sig, Ox(i), Oy(i), Oz(i), x(i), y(i), z(i), td(i))  <=  crit_thresh;
    	   
    	   %% Making U1 greater than certain number
    	   linU1constFunc(obs_sig, drone_sig, Ox(i), Oy(i), Oz(i), x(i), y(i), z(i), td(i)) +  ...
	   	   (Px(i) - x(i))*linU1PxFunc(obs_sig, drone_sig, Ox(i), Oy(i), Oz(i), x(i), y(i), z(i), td(i)) +  ...
    	   			(Py(i) - y(i))*linU1PyFunc(obs_sig, drone_sig, Ox(i), Oy(i), Oz(i), x(i), y(i), z(i), td(i)) +  ...
    	   				(Pz(i) - z(i))*linU1PzFunc(obs_sig, drone_sig, Ox(i), Oy(i), Oz(i), x(i), y(i), z(i), td(i)) +  ...
    	   					(t(i) - td(i))*linU1tFunc(obs_sig, drone_sig, Ox(i), Oy(i), Oz(i), x(i), y(i), z(i), td(i)) >= u1_thresh;
    	   
    	end
    end
    %%%%%%%% Constraints %%%%%%%%
    Vx >= 0;
    Vy >= 0;
    Vz >= 0;
    
    %% Constraint on t
    0 <= t <= 1;
    
    %% The velocity during initialization shouldn't be much higher
    0 <= Vx(1) <= 1.0;
    0 <= Vy(1) <= 1.0;
    0 <= Vz(1) <= 1.0;
    
    %% Ensuring that subsequent X vel commands do not have difference more than del_Vx
    0 <= Vx(2:n) - Vx(1:n-1) <= del_Vx;
    
    %% Ensuring that subsequent Y vel commands do not have difference more than del_Vy
    0 <= Vy(2:n) - Vy(1:n-1) <= del_Vy;
    
    %% Ensuring that subsequent Z vel commands do not have difference more than del_Vz
    0 <= Vz(2:n) - Vz(1:n-1) <= del_Vz;
    
    %%%%%%%% Constraints definition over!! %%%%%%%%%
    cvx_end
    %%	Defining path of the obstacle based on first optimization
    %if l ==1
    %	    Ob_middle=[Px(n/2), Py(n/2), Pz(n/2)];
    %	    Ob_end = 2*Ob_middle - Ob_start; %% This will be determined based on first velocity solution of the drone
    %end   
    %% Updating the linearization points!
    x = Px;
    y = Py;
    z = Pz;
    %cvx_optval
    if(l ~= no_of_iter-1)
        clear Vx Vy Vz Px Py Pz
    end	
    if(l>1)
    	td = t;
    	%confidence = confidence + 0.025;
    	% Finding U1 threshold for the above confidence
	P_u1 = Overlap(int16(confidence*length(Conf)))/2;
	CP_u1 = 1 - P_u1; %% Cumulative probability of U1
	u1_thresh = norminv(CP_u1);
    end
    
    l = l + 1;
    %R = R + 1
end %% End of while loop

for i = 1:n
	overlap_area = Gaussian_overlap_func([Px(i); Py(i); Pz(i)], drone_sig,  [Ox(i); Oy(i); Oz(i)], obs_sig);
	%% Getting confidence interval on a scale of 0-1 for above overlap_area
		flag = true;
		overlap_lower = overlap_area - 0.00001;
		overlap_upper = overlap_area + 0.00001;

		while(flag)
	 	
			conf_indices = find(Overlap >= overlap_lower & Overlap <= overlap_upper);
			s_conf = size(conf_indices);
			if s_conf(2) == 0
				overlap_lower = overlap_lower - 0.000001;
				overlap_upper = overlap_upper + 0.000001;
			else
				flag = false;
			end
		end

		%% This is the confidence interval by which drone will avoid the obstacle
		conf_interval(i) = Conf(conf_indices(1));
end


