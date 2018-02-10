
clc;
clear all;
close all;


%% Drone is supposed to reach (dest_x, dest_y, dest_z) from (0,0,0) in
%% 'n' timesteps with each timestep of duration del_t
%% Version 2 (for higher values of n)
n = 40;      % Number of timesteps
del_t = 1.0; % Time duration for which some command is executed.
tot_time = n*del_t; %% Total time of flight

%% defining drone uncertainty
drone_sig = [0.4, 0.0, 0.0; 0.0, 0.5, 0.0; 0.0, 0.0, 0.4];

%% Defining obstacle trajectory
%% For now we will assume that obstacle is moving along a straing line. it's starting from (0,0,0) and It's trajectory will be determined based on the destination of the drone

%% Drone Destination variables
dest_x = 12;
dest_y = 12;
dest_z = 12;

%%%%%%%%%%%%%% Obstacle stuff starts here!!!! %%%%%%%%%%%%%%%%%
%% Here, the drone is moving from (0,0,0) to (dest_x, dest_y, dest_z) in n timesteps
%% The obstacle will reach from start_obs to end_obs in n timesteps 
%%Ob_start=[0, dest_y/2, dest_z];
%%Ob_end=[dest_x, dest_y/2, 0];
Ob_start=[0,2.2,0];
Ob_end=[dest_x, dest_y+2.2, dest_z];
%% Calculating velocity of the obstacle in all the directions
vel_obs = (Ob_end - Ob_start)/tot_time;

%% Defining Obstacle path location at each timestep
for i = 1:n
  Ox(i) = Ob_start(1) + del_t*i*vel_obs(1);
  Oy(i) = Ob_start(2) + del_t*i*vel_obs(2);
  Oz(i) = Ob_start(3) + del_t*i*vel_obs(3);
end


%% Defining obstalce uncertainty
obs_sig = [0.6, 0.0, 0.0; 0.0, 0.4, 0.0; 0.0, 0.0, 0.7];;

%% Criterion threshold
crit_thresh =  0.00001;

%% The minimum confidence by which both, the drone and obstacle should be avoided
confidence_low = 0.30;
confidence_upper = 0.90;

confidence_e = 0.90;

load Linear_approx.mat
%% Finding U1 threshold for the above confidence
P_u1_low = Overlap(int16(confidence_low*length(Conf)))/2;
CP_u1_low = 1 - P_u1_low; %% Cumulative probability of U1
u1_thresh_low = norminv(CP_u1_low);

%% Finding U1 threshold for the above confidence
P_u1_upper = Overlap(int16(confidence_upper*length(Conf)))/2;
CP_u1_upper = 1 - P_u1_upper; %% Cumulative probability of U1
u1_thresh_upper = norminv(CP_u1_upper);

%%%%%%%%%Obstacle stuff ends here %%%%%%%%%%%%%%%

%% Difference between subsequent velocities
del_Vx = 0.05;
del_Vy = 0.05;
del_Vz = 0.05;

%% Initializing linearization point for ts
for i = 1:n
	td(i) = 0.3;
end

l = 0;
no_of_iter = 4;

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
    if(l~=0)
    	for i = 1:n
    	   %% Linearized criteria for minimizing overlap through linear approximation
    	   -crit_thresh <= linconstraintconstFunc(obs_sig, drone_sig, Ox(i), Oy(i), Oz(i), x(i), y(i), z(i), td(i)) +  ...
    	   					(Px(i) - x(i))*linconstraintPxFunc(obs_sig, drone_sig, Ox(i), Oy(i), Oz(i), x(i), y(i), z(i), td(i)) +  ...
    	   						(Py(i) - y(i))*linconstraintPyFunc(obs_sig, drone_sig, Ox(i), Oy(i), Oz(i), x(i), y(i), z(i), td(i)) +  ...
    	   							(Pz(i) - z(i))*linconstraintPzFunc(obs_sig, drone_sig, Ox(i), Oy(i), Oz(i), x(i), y(i), z(i), td(i)) +  ...
    	   								(t(i) - td(i))*linconstrainttFunc(obs_sig, drone_sig, Ox(i), Oy(i), Oz(i), x(i), y(i), z(i), td(i))  <=  crit_thresh;
    	   
    	   %% Making U1 greater than certain number
    	   u1_thresh_upper >= linU1constFunc(obs_sig, drone_sig, Ox(i), Oy(i), Oz(i), x(i), y(i), z(i), td(i)) +  ...
						    	   (Px(i) - x(i))*linU1PxFunc(obs_sig, drone_sig, Ox(i), Oy(i), Oz(i), x(i), y(i), z(i), td(i)) +  ...
    	   							 (Py(i) - y(i))*linU1PyFunc(obs_sig, drone_sig, Ox(i), Oy(i), Oz(i), x(i), y(i), z(i), td(i)) +  ...
    	   								(Pz(i) - z(i))*linU1PzFunc(obs_sig, drone_sig, Ox(i), Oy(i), Oz(i), x(i), y(i), z(i), td(i)) +  ...
    	  				 					(t(i) - td(i))*linU1tFunc(obs_sig, drone_sig, Ox(i), Oy(i), Oz(i), x(i), y(i), z(i), td(i)) >= u1_thresh_low;
    	   
    	end
    end
    %%%%%%%% Constraints %%%%%%%%
    Vx >= 0;
    Vy >= 0;
    Vz >= 0;
    
    %% Constraint on t
    0 <= t <= 1;
    
    %% The velocity during initialization shouldn't be much higher
    0 <= Vx(1) <= 0.3;
    0 <= Vy(1) <= 0.3;
    0 <= Vz(1) <= 0.3;
    
    %	% Ensuring that subsequent X vel commands do not have difference more than del_Vx
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
%         confidence_low = confidence_low + 0.05;
%         P_u1_low = Overlap(int16(confidence_low*length(Conf)))/2;
%         CP_u1_low = 1 - P_u1_low; %% Cumulative probability of U1
%         u1_thresh_low = norminv(CP_u1_low);
    end
    
    l = l + 1;
end %% End of while loop


v = VideoWriter('Drone_graph.avi', 'Motion JPEG AVI');
v.FrameRate = 10;
open(v);

hfig = figure('units','normalized','outerposition',[0 0 1 1]);
hold on;
pause(3);

for i = 1:(n-1)
	j_max = 20;	
	i
	for j = 1:j_max
		j;
		%% Defining mean and covariances
		mu1 = [Px(i) + (Px(i+1) - Px(i))*j/j_max;  Py(i) + (Py(i+1) - Py(i))*j/j_max;  Pz(i)  + (Pz(i+1) - Pz(i))*j/j_max];
		sig1 = drone_sig;
		mu2 = [Ox(i) + (Ox(i+1) - Ox(i))*j/j_max;  Oy(i) + (Oy(i+1) - Oy(i))*j/j_max;  Oz(i) + (Oz(i+1) - Oz(i))*j/j_max];
		sig2 = obs_sig;

		%% Getting the intersection area
		overlap_area = Gaussian_overlap_func(mu1, sig1, mu2, sig2);

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
		conf_interval = Conf(conf_indices(1));
		mahaDis = chi2inv(conf_interval, 2);

		[V1, D1] = eig(sig1);
		eig1 = diag(V1*D1*V1');
		[V2, D2] = eig(sig2);
		eig2 = diag(V2*D2*V2');

		%% Getting lengths of semi major axis based on confidence values

		%% Radius of ellipse along X-direction for distribution 1
		rad_X1 = sqrt(mahaDis*eig1(1)); 

		%% Radius of ellipse along Y-direction for distribution 1
		rad_Y1 = sqrt(mahaDis*eig1(2)); 

		%% Radius of ellipse along Z-direction for distribution 1
		rad_Z1 = sqrt(mahaDis*eig1(3)); 


		%% Radius of ellipse along X-direction for distribution 2
		rad_X2 = sqrt(mahaDis*eig2(1)); 
	
		%% Radius of ellipse along Y-direction for distribution 2
		rad_Y2 = sqrt(mahaDis*eig2(2)); 
	
		%% Radius of ellipse along Z-direction for distribution 2
		rad_Z2 = sqrt(mahaDis*eig2(3)); 

        %%Outermost Ellipse with confidence_e
		mahaDis_e = chi2inv(confidence_e, 2);

		[V1_e, D1_e] = eig(sig1);
		eig1_e = diag(V1_e*D1_e*V1_e');
		[V2_e, D2_e] = eig(sig2);
		eig2_e = diag(V2_e*D2_e*V2_e');

		%% Getting lengths of semi major axis based on confidence_e values

		%% Radius of outermost ellipse along X-direction for distribution 1
		rad_X1_e = sqrt(mahaDis_e*eig1_e(1)); 

		%% Radius of ellipse along Y-direction for distribution 1
		rad_Y1_e = sqrt(mahaDis_e*eig1_e(2)); 

		%% Radius of ellipse along Z-direction for distribution 1
		rad_Z1_e = sqrt(mahaDis_e*eig1_e(3)); 


		%% Radius of ellipse along X-direction for distribution 2
		rad_X2_e = sqrt(mahaDis_e*eig2_e(1)); 
	
		%% Radius of ellipse along Y-direction for distribution 2
		rad_Y2_e = sqrt(mahaDis_e*eig2_e(2)); 
	
		%% Radius of ellipse along Z-direction for distribution 2
		rad_Z2_e = sqrt(mahaDis_e*eig2_e(3)); 

        
		%% Drawing confidence ellipses around the obstacle
		h1 = subplot(1,1,1);
		set(gca,'XLim', [floor(min([Px, Ox]))-1, floor(max([Px, Ox]))+1], 'YLim', [floor(min([Py, Oy]))-1, floor(max([Py, Oy]))+1], 'ZLim',[floor(min([Pz, Oz])) - 1, floor(max([Pz, Oz]))+1]);
		hold on;
		view(81.6000, 28.4000);
		hold on;
		grid on;
		hold on;
		%% Drawing confidence ellipses around the obstacle
		title(strcat(' Total time = ', num2str(tot_time), ' Confidence of avoidance = ', num2str(conf_interval), '.', ' Min confidence: ', num2str(confidence_low), ' Max confidence: ', num2str(confidence_upper)));
		hold on;		
		
		scatter3(mu1(1),  mu1(2), mu1(3),1000, [0, 0, 0], 'filled');
		%scatter3(Ox(i), Oy(i),Oz(i), 100,  [1 - i/n, 1/2, i/n], 'filled');
		hold on;
		[el1_x, el1_y, el1_z] = ellipsoid(mu1(1),  mu1(2), mu1(3), rad_X1, rad_Y1, rad_Z1, 50);
		surf(el1_x, el1_y, el1_z, 'FaceAlpha', 0.4, 'EdgeColor','none', 'FaceColor', [0.9 0.5 0.1]);
		hold on;
        [el1_x_e, el1_y_e, el1_z_e] = ellipsoid(mu1(1),  mu1(2), mu1(3), rad_X1_e, rad_Y1_e, rad_Z1_e, 50);
		surf(el1_x_e, el1_y_e, el1_z_e, 'FaceAlpha', 0.4, 'EdgeColor','none', 'FaceColor', [0.9 0.5 0.1]);
		hold on;

		scatter3(mu2(1),  mu2(2), mu2(3), 1000,  [1, 0, 1], 'filled');
		hold on;	
		[el2_x, el2_y, el2_z] = ellipsoid(mu2(1),  mu2(2), mu2(3), rad_X2, rad_Y2, rad_Z2, 50);
		surf(el2_x, el2_y, el2_z, 'FaceAlpha', 0.4, 'EdgeColor','none', 'FaceColor', [0.1 0.5 0.9]);
		hold on;
        [el2_x_e, el2_y_e, el2_z_e] = ellipsoid(mu2(1),  mu2(2), mu2(3), rad_X2_e, rad_Y2_e, rad_Z2_e, 50);
		surf(el2_x_e, el2_y_e, el2_z_e, 'FaceAlpha', 0.4, 'EdgeColor','none', 'FaceColor', [0.1 0.5 0.9]);
		hold on;
		
		M(i*j) = getframe(hfig);
		writeVideo(v,M(i*j));
		delete(h1);
	end
end
close(v);

%movie2avi(M,'Drone_trajectory.avi');


