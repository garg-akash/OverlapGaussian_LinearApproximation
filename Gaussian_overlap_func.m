function [ overlap ] = Gaussian_overlap_func( mu1, sig1, mu2, sig2 )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

%% Defining precision
prec = 0.00000000001;

%% Incremeant factor for t
incr = 0.5;

%% Initializing the value of t
t = 0.5;

%% Initializing crit
crit = 0.8;

%% Linear approximation algo begins
while(crit > prec || crit < -prec)

	b = inv(t*sig1 + (1 - t)*sig2)*(mu2 - mu1);
	crit = transpose(b)*(t^2*sig1 - (1 - t)^2*sig2)*b;
	if crit > prec
		t = t - incr;
	end
	if crit < -prec
		t = t + incr;
	end
	incr = incr/2;
end


c1 = transpose(b)*mu1 + t*transpose(b)*sig1*b;
c2 = transpose(b)*mu2 - (1-t)*transpose(b)*sig2*b;

u1 = (c1 - transpose(b)*mu1)/sqrt(transpose(b)*sig1*b);
u2 = (transpose(b)*mu2 - c1)/sqrt(transpose(b)*sig2*b);

P_u1 = 1 - normcdf(u1);
P_u2 = 1 - normcdf(u2);
overlap = P_u1 + P_u2;

end

