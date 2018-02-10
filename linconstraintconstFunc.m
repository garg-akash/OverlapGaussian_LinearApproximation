function [ linconstraintconst] = linconstraintconstFunc(sig_obs, sig_drone, Ox, Oy, Oz, Pxd, Pyd, Pzd, td)
% Returns the constant term after linearization
% 

%% Extracting a1,b1,c1,b1,d1,e1,c1,e1,f1,a2,b2,c2,b2,d2,e2,c2,e2,f2,ox,oy,oz,pxdash,pydash,pzdash,tdash
%% Getting obstacle uncertainty variables
a1 =sig_obs(1,1);
b1 =sig_obs(1,2);
c1 =sig_obs(1,3);
d1 =sig_obs(2,2);
e1 =sig_obs(2,3);
f1 =sig_obs(3,3);

%% Getting drone uncertainty variables
a2 =sig_drone(1,1);
b2 =sig_drone(1,2);
c2 =sig_drone(1,3);
d2 =sig_drone(2,2);
e2 =sig_drone(2,3);
f2 =sig_drone(3,3);

%% Obstacle position
ox = Ox;
oy = Oy;
oz = Oz;

%% Linearization points
pxdash = Pxd;
pydash = Pyd;
pzdash = Pzd;
tdash = td;

linconstraintconst=((f2+f1.*tdash+(-1).*f2.*tdash).*((-1).*(b2+b1.*tdash+(-1).*b2.*tdash) ...
  .^2+(a2+a1.*tdash+(-1).*a2.*tdash).*(d2+d1.*tdash+(-1).*d2.*tdash))+(-1) ...
  .*(e2+e1.*tdash+(-1).*e2.*tdash).*((-1).*(b2+b1.*tdash+(-1).*b2.*tdash) ...
  .*(c2+c1.*tdash+(-1).*c2.*tdash)+(a2+a1.*tdash+(-1).*a2.*tdash).*(e2+ ...
  e1.*tdash+(-1).*e2.*tdash))+(c2+c1.*tdash+(-1).*c2.*tdash).*((-1).*(c2+ ...
  c1.*tdash+(-1).*c2.*tdash).*(d2+d1.*tdash+(-1).*d2.*tdash)+(b2+b1.* ...
  tdash+(-1).*b2.*tdash).*(e2+e1.*tdash+(-1).*e2.*tdash))).^(-1).*(((-1).* ...
  oz+pzdash).*((-1).*(c2+c1.*tdash+(-1).*c2.*tdash).*(d2+d1.*tdash+(-1).* ...
  d2.*tdash)+(b2+b1.*tdash+(-1).*b2.*tdash).*(e2+e1.*tdash+(-1).*e2.* ...
  tdash))+((-1).*oy+pydash).*((c2+c1.*tdash+(-1).*c2.*tdash).*(e2+e1.* ...
  tdash+(-1).*e2.*tdash)+(-1).*(b2+b1.*tdash+(-1).*b2.*tdash).*(f2+f1.* ...
  tdash+(-1).*f2.*tdash))+((-1).*ox+pxdash).*((-1).*(e2+e1.*tdash+(-1).* ...
  e2.*tdash).^2+(d2+d1.*tdash+(-1).*d2.*tdash).*(f2+f1.*tdash+(-1).*f2.* ...
  tdash))).*(((-1).*c2.*((-1)+tdash).^2+c1.*tdash.^2).*(((-1).*oz+pzdash) ...
  .*((-1).*(b2+b1.*tdash+(-1).*b2.*tdash).^2+(a2+a1.*tdash+(-1).*a2.* ...
  tdash).*(d2+d1.*tdash+(-1).*d2.*tdash))+((-1).*oy+pydash).*((b2+b1.* ...
  tdash+(-1).*b2.*tdash).*(c2+c1.*tdash+(-1).*c2.*tdash)+(-1).*(a2+a1.* ...
  tdash+(-1).*a2.*tdash).*(e2+e1.*tdash+(-1).*e2.*tdash))+((-1).*ox+ ...
  pxdash).*((-1).*(c2+c1.*tdash+(-1).*c2.*tdash).*(d2+d1.*tdash+(-1).*d2.* ...
  tdash)+(b2+b1.*tdash+(-1).*b2.*tdash).*(e2+e1.*tdash+(-1).*e2.*tdash))) ...
  .*((f2+f1.*tdash+(-1).*f2.*tdash).*((-1).*(b2+b1.*tdash+(-1).*b2.*tdash) ...
  .^2+(a2+a1.*tdash+(-1).*a2.*tdash).*(d2+d1.*tdash+(-1).*d2.*tdash))+(-1) ...
  .*(e2+e1.*tdash+(-1).*e2.*tdash).*((-1).*(b2+b1.*tdash+(-1).*b2.*tdash) ...
  .*(c2+c1.*tdash+(-1).*c2.*tdash)+(a2+a1.*tdash+(-1).*a2.*tdash).*(e2+ ...
  e1.*tdash+(-1).*e2.*tdash))+(c2+c1.*tdash+(-1).*c2.*tdash).*((-1).*(c2+ ...
  c1.*tdash+(-1).*c2.*tdash).*(d2+d1.*tdash+(-1).*d2.*tdash)+(b2+b1.* ...
  tdash+(-1).*b2.*tdash).*(e2+e1.*tdash+(-1).*e2.*tdash))).^(-1)+((-1).* ...
  b2.*((-1)+tdash).^2+b1.*tdash.^2).*((f2+f1.*tdash+(-1).*f2.*tdash).*(( ...
  -1).*(b2+b1.*tdash+(-1).*b2.*tdash).^2+(a2+a1.*tdash+(-1).*a2.*tdash).*( ...
  d2+d1.*tdash+(-1).*d2.*tdash))+(-1).*(e2+e1.*tdash+(-1).*e2.*tdash).*(( ...
  -1).*(b2+b1.*tdash+(-1).*b2.*tdash).*(c2+c1.*tdash+(-1).*c2.*tdash)+(a2+ ...
  a1.*tdash+(-1).*a2.*tdash).*(e2+e1.*tdash+(-1).*e2.*tdash))+(c2+c1.* ...
  tdash+(-1).*c2.*tdash).*((-1).*(c2+c1.*tdash+(-1).*c2.*tdash).*(d2+d1.* ...
  tdash+(-1).*d2.*tdash)+(b2+b1.*tdash+(-1).*b2.*tdash).*(e2+e1.*tdash+( ...
  -1).*e2.*tdash))).^(-1).*(((-1).*oz+pzdash).*((b2+b1.*tdash+(-1).*b2.* ...
  tdash).*(c2+c1.*tdash+(-1).*c2.*tdash)+(-1).*(a2+a1.*tdash+(-1).*a2.* ...
  tdash).*(e2+e1.*tdash+(-1).*e2.*tdash))+((-1).*oy+pydash).*((-1).*(c2+ ...
  c1.*tdash+(-1).*c2.*tdash).^2+(a2+a1.*tdash+(-1).*a2.*tdash).*(f2+f1.* ...
  tdash+(-1).*f2.*tdash))+((-1).*ox+pxdash).*((c2+c1.*tdash+(-1).*c2.* ...
  tdash).*(e2+e1.*tdash+(-1).*e2.*tdash)+(-1).*(b2+b1.*tdash+(-1).*b2.* ...
  tdash).*(f2+f1.*tdash+(-1).*f2.*tdash)))+((-1).*a2.*((-1)+tdash).^2+a1.* ...
  tdash.^2).*((f2+f1.*tdash+(-1).*f2.*tdash).*((-1).*(b2+b1.*tdash+(-1).* ...
  b2.*tdash).^2+(a2+a1.*tdash+(-1).*a2.*tdash).*(d2+d1.*tdash+(-1).*d2.* ...
  tdash))+(-1).*(e2+e1.*tdash+(-1).*e2.*tdash).*((-1).*(b2+b1.*tdash+(-1) ...
  .*b2.*tdash).*(c2+c1.*tdash+(-1).*c2.*tdash)+(a2+a1.*tdash+(-1).*a2.* ...
  tdash).*(e2+e1.*tdash+(-1).*e2.*tdash))+(c2+c1.*tdash+(-1).*c2.*tdash).* ...
  ((-1).*(c2+c1.*tdash+(-1).*c2.*tdash).*(d2+d1.*tdash+(-1).*d2.*tdash)+( ...
  b2+b1.*tdash+(-1).*b2.*tdash).*(e2+e1.*tdash+(-1).*e2.*tdash))).^(-1).*( ...
  ((-1).*oz+pzdash).*((-1).*(c2+c1.*tdash+(-1).*c2.*tdash).*(d2+d1.*tdash+ ...
  (-1).*d2.*tdash)+(b2+b1.*tdash+(-1).*b2.*tdash).*(e2+e1.*tdash+(-1).* ...
  e2.*tdash))+((-1).*oy+pydash).*((c2+c1.*tdash+(-1).*c2.*tdash).*(e2+e1.* ...
  tdash+(-1).*e2.*tdash)+(-1).*(b2+b1.*tdash+(-1).*b2.*tdash).*(f2+f1.* ...
  tdash+(-1).*f2.*tdash))+((-1).*ox+pxdash).*((-1).*(e2+e1.*tdash+(-1).* ...
  e2.*tdash).^2+(d2+d1.*tdash+(-1).*d2.*tdash).*(f2+f1.*tdash+(-1).*f2.* ...
  tdash))))+((f2+f1.*tdash+(-1).*f2.*tdash).*((-1).*(b2+b1.*tdash+(-1).* ...
  b2.*tdash).^2+(a2+a1.*tdash+(-1).*a2.*tdash).*(d2+d1.*tdash+(-1).*d2.* ...
  tdash))+(-1).*(e2+e1.*tdash+(-1).*e2.*tdash).*((-1).*(b2+b1.*tdash+(-1) ...
  .*b2.*tdash).*(c2+c1.*tdash+(-1).*c2.*tdash)+(a2+a1.*tdash+(-1).*a2.* ...
  tdash).*(e2+e1.*tdash+(-1).*e2.*tdash))+(c2+c1.*tdash+(-1).*c2.*tdash).* ...
  ((-1).*(c2+c1.*tdash+(-1).*c2.*tdash).*(d2+d1.*tdash+(-1).*d2.*tdash)+( ...
  b2+b1.*tdash+(-1).*b2.*tdash).*(e2+e1.*tdash+(-1).*e2.*tdash))).^(-1).*( ...
  ((-1).*oz+pzdash).*((b2+b1.*tdash+(-1).*b2.*tdash).*(c2+c1.*tdash+(-1).* ...
  c2.*tdash)+(-1).*(a2+a1.*tdash+(-1).*a2.*tdash).*(e2+e1.*tdash+(-1).* ...
  e2.*tdash))+((-1).*oy+pydash).*((-1).*(c2+c1.*tdash+(-1).*c2.*tdash).^2+ ...
  (a2+a1.*tdash+(-1).*a2.*tdash).*(f2+f1.*tdash+(-1).*f2.*tdash))+((-1).* ...
  ox+pxdash).*((c2+c1.*tdash+(-1).*c2.*tdash).*(e2+e1.*tdash+(-1).*e2.* ...
  tdash)+(-1).*(b2+b1.*tdash+(-1).*b2.*tdash).*(f2+f1.*tdash+(-1).*f2.* ...
  tdash))).*(((-1).*e2.*((-1)+tdash).^2+e1.*tdash.^2).*(((-1).*oz+pzdash) ...
  .*((-1).*(b2+b1.*tdash+(-1).*b2.*tdash).^2+(a2+a1.*tdash+(-1).*a2.* ...
  tdash).*(d2+d1.*tdash+(-1).*d2.*tdash))+((-1).*oy+pydash).*((b2+b1.* ...
  tdash+(-1).*b2.*tdash).*(c2+c1.*tdash+(-1).*c2.*tdash)+(-1).*(a2+a1.* ...
  tdash+(-1).*a2.*tdash).*(e2+e1.*tdash+(-1).*e2.*tdash))+((-1).*ox+ ...
  pxdash).*((-1).*(c2+c1.*tdash+(-1).*c2.*tdash).*(d2+d1.*tdash+(-1).*d2.* ...
  tdash)+(b2+b1.*tdash+(-1).*b2.*tdash).*(e2+e1.*tdash+(-1).*e2.*tdash))) ...
  .*((f2+f1.*tdash+(-1).*f2.*tdash).*((-1).*(b2+b1.*tdash+(-1).*b2.*tdash) ...
  .^2+(a2+a1.*tdash+(-1).*a2.*tdash).*(d2+d1.*tdash+(-1).*d2.*tdash))+(-1) ...
  .*(e2+e1.*tdash+(-1).*e2.*tdash).*((-1).*(b2+b1.*tdash+(-1).*b2.*tdash) ...
  .*(c2+c1.*tdash+(-1).*c2.*tdash)+(a2+a1.*tdash+(-1).*a2.*tdash).*(e2+ ...
  e1.*tdash+(-1).*e2.*tdash))+(c2+c1.*tdash+(-1).*c2.*tdash).*((-1).*(c2+ ...
  c1.*tdash+(-1).*c2.*tdash).*(d2+d1.*tdash+(-1).*d2.*tdash)+(b2+b1.* ...
  tdash+(-1).*b2.*tdash).*(e2+e1.*tdash+(-1).*e2.*tdash))).^(-1)+((-1).* ...
  d2.*((-1)+tdash).^2+d1.*tdash.^2).*((f2+f1.*tdash+(-1).*f2.*tdash).*(( ...
  -1).*(b2+b1.*tdash+(-1).*b2.*tdash).^2+(a2+a1.*tdash+(-1).*a2.*tdash).*( ...
  d2+d1.*tdash+(-1).*d2.*tdash))+(-1).*(e2+e1.*tdash+(-1).*e2.*tdash).*(( ...
  -1).*(b2+b1.*tdash+(-1).*b2.*tdash).*(c2+c1.*tdash+(-1).*c2.*tdash)+(a2+ ...
  a1.*tdash+(-1).*a2.*tdash).*(e2+e1.*tdash+(-1).*e2.*tdash))+(c2+c1.* ...
  tdash+(-1).*c2.*tdash).*((-1).*(c2+c1.*tdash+(-1).*c2.*tdash).*(d2+d1.* ...
  tdash+(-1).*d2.*tdash)+(b2+b1.*tdash+(-1).*b2.*tdash).*(e2+e1.*tdash+( ...
  -1).*e2.*tdash))).^(-1).*(((-1).*oz+pzdash).*((b2+b1.*tdash+(-1).*b2.* ...
  tdash).*(c2+c1.*tdash+(-1).*c2.*tdash)+(-1).*(a2+a1.*tdash+(-1).*a2.* ...
  tdash).*(e2+e1.*tdash+(-1).*e2.*tdash))+((-1).*oy+pydash).*((-1).*(c2+ ...
  c1.*tdash+(-1).*c2.*tdash).^2+(a2+a1.*tdash+(-1).*a2.*tdash).*(f2+f1.* ...
  tdash+(-1).*f2.*tdash))+((-1).*ox+pxdash).*((c2+c1.*tdash+(-1).*c2.* ...
  tdash).*(e2+e1.*tdash+(-1).*e2.*tdash)+(-1).*(b2+b1.*tdash+(-1).*b2.* ...
  tdash).*(f2+f1.*tdash+(-1).*f2.*tdash)))+((-1).*b2.*((-1)+tdash).^2+b1.* ...
  tdash.^2).*((f2+f1.*tdash+(-1).*f2.*tdash).*((-1).*(b2+b1.*tdash+(-1).* ...
  b2.*tdash).^2+(a2+a1.*tdash+(-1).*a2.*tdash).*(d2+d1.*tdash+(-1).*d2.* ...
  tdash))+(-1).*(e2+e1.*tdash+(-1).*e2.*tdash).*((-1).*(b2+b1.*tdash+(-1) ...
  .*b2.*tdash).*(c2+c1.*tdash+(-1).*c2.*tdash)+(a2+a1.*tdash+(-1).*a2.* ...
  tdash).*(e2+e1.*tdash+(-1).*e2.*tdash))+(c2+c1.*tdash+(-1).*c2.*tdash).* ...
  ((-1).*(c2+c1.*tdash+(-1).*c2.*tdash).*(d2+d1.*tdash+(-1).*d2.*tdash)+( ...
  b2+b1.*tdash+(-1).*b2.*tdash).*(e2+e1.*tdash+(-1).*e2.*tdash))).^(-1).*( ...
  ((-1).*oz+pzdash).*((-1).*(c2+c1.*tdash+(-1).*c2.*tdash).*(d2+d1.*tdash+ ...
  (-1).*d2.*tdash)+(b2+b1.*tdash+(-1).*b2.*tdash).*(e2+e1.*tdash+(-1).* ...
  e2.*tdash))+((-1).*oy+pydash).*((c2+c1.*tdash+(-1).*c2.*tdash).*(e2+e1.* ...
  tdash+(-1).*e2.*tdash)+(-1).*(b2+b1.*tdash+(-1).*b2.*tdash).*(f2+f1.* ...
  tdash+(-1).*f2.*tdash))+((-1).*ox+pxdash).*((-1).*(e2+e1.*tdash+(-1).* ...
  e2.*tdash).^2+(d2+d1.*tdash+(-1).*d2.*tdash).*(f2+f1.*tdash+(-1).*f2.* ...
  tdash))))+(((-1).*oz+pzdash).*((-1).*(b2+b1.*tdash+(-1).*b2.*tdash).^2+( ...
  a2+a1.*tdash+(-1).*a2.*tdash).*(d2+d1.*tdash+(-1).*d2.*tdash))+((-1).* ...
  oy+pydash).*((b2+b1.*tdash+(-1).*b2.*tdash).*(c2+c1.*tdash+(-1).*c2.* ...
  tdash)+(-1).*(a2+a1.*tdash+(-1).*a2.*tdash).*(e2+e1.*tdash+(-1).*e2.* ...
  tdash))+((-1).*ox+pxdash).*((-1).*(c2+c1.*tdash+(-1).*c2.*tdash).*(d2+ ...
  d1.*tdash+(-1).*d2.*tdash)+(b2+b1.*tdash+(-1).*b2.*tdash).*(e2+e1.* ...
  tdash+(-1).*e2.*tdash))).*((f2+f1.*tdash+(-1).*f2.*tdash).*((-1).*(b2+ ...
  b1.*tdash+(-1).*b2.*tdash).^2+(a2+a1.*tdash+(-1).*a2.*tdash).*(d2+d1.* ...
  tdash+(-1).*d2.*tdash))+(-1).*(e2+e1.*tdash+(-1).*e2.*tdash).*((-1).*( ...
  b2+b1.*tdash+(-1).*b2.*tdash).*(c2+c1.*tdash+(-1).*c2.*tdash)+(a2+a1.* ...
  tdash+(-1).*a2.*tdash).*(e2+e1.*tdash+(-1).*e2.*tdash))+(c2+c1.*tdash+( ...
  -1).*c2.*tdash).*((-1).*(c2+c1.*tdash+(-1).*c2.*tdash).*(d2+d1.*tdash+( ...
  -1).*d2.*tdash)+(b2+b1.*tdash+(-1).*b2.*tdash).*(e2+e1.*tdash+(-1).*e2.* ...
  tdash))).^(-1).*(((-1).*f2.*((-1)+tdash).^2+f1.*tdash.^2).*(((-1).*oz+ ...
  pzdash).*((-1).*(b2+b1.*tdash+(-1).*b2.*tdash).^2+(a2+a1.*tdash+(-1).* ...
  a2.*tdash).*(d2+d1.*tdash+(-1).*d2.*tdash))+((-1).*oy+pydash).*((b2+b1.* ...
  tdash+(-1).*b2.*tdash).*(c2+c1.*tdash+(-1).*c2.*tdash)+(-1).*(a2+a1.* ...
  tdash+(-1).*a2.*tdash).*(e2+e1.*tdash+(-1).*e2.*tdash))+((-1).*ox+ ...
  pxdash).*((-1).*(c2+c1.*tdash+(-1).*c2.*tdash).*(d2+d1.*tdash+(-1).*d2.* ...
  tdash)+(b2+b1.*tdash+(-1).*b2.*tdash).*(e2+e1.*tdash+(-1).*e2.*tdash))) ...
  .*((f2+f1.*tdash+(-1).*f2.*tdash).*((-1).*(b2+b1.*tdash+(-1).*b2.*tdash) ...
  .^2+(a2+a1.*tdash+(-1).*a2.*tdash).*(d2+d1.*tdash+(-1).*d2.*tdash))+(-1) ...
  .*(e2+e1.*tdash+(-1).*e2.*tdash).*((-1).*(b2+b1.*tdash+(-1).*b2.*tdash) ...
  .*(c2+c1.*tdash+(-1).*c2.*tdash)+(a2+a1.*tdash+(-1).*a2.*tdash).*(e2+ ...
  e1.*tdash+(-1).*e2.*tdash))+(c2+c1.*tdash+(-1).*c2.*tdash).*((-1).*(c2+ ...
  c1.*tdash+(-1).*c2.*tdash).*(d2+d1.*tdash+(-1).*d2.*tdash)+(b2+b1.* ...
  tdash+(-1).*b2.*tdash).*(e2+e1.*tdash+(-1).*e2.*tdash))).^(-1)+((-1).* ...
  e2.*((-1)+tdash).^2+e1.*tdash.^2).*((f2+f1.*tdash+(-1).*f2.*tdash).*(( ...
  -1).*(b2+b1.*tdash+(-1).*b2.*tdash).^2+(a2+a1.*tdash+(-1).*a2.*tdash).*( ...
  d2+d1.*tdash+(-1).*d2.*tdash))+(-1).*(e2+e1.*tdash+(-1).*e2.*tdash).*(( ...
  -1).*(b2+b1.*tdash+(-1).*b2.*tdash).*(c2+c1.*tdash+(-1).*c2.*tdash)+(a2+ ...
  a1.*tdash+(-1).*a2.*tdash).*(e2+e1.*tdash+(-1).*e2.*tdash))+(c2+c1.* ...
  tdash+(-1).*c2.*tdash).*((-1).*(c2+c1.*tdash+(-1).*c2.*tdash).*(d2+d1.* ...
  tdash+(-1).*d2.*tdash)+(b2+b1.*tdash+(-1).*b2.*tdash).*(e2+e1.*tdash+( ...
  -1).*e2.*tdash))).^(-1).*(((-1).*oz+pzdash).*((b2+b1.*tdash+(-1).*b2.* ...
  tdash).*(c2+c1.*tdash+(-1).*c2.*tdash)+(-1).*(a2+a1.*tdash+(-1).*a2.* ...
  tdash).*(e2+e1.*tdash+(-1).*e2.*tdash))+((-1).*oy+pydash).*((-1).*(c2+ ...
  c1.*tdash+(-1).*c2.*tdash).^2+(a2+a1.*tdash+(-1).*a2.*tdash).*(f2+f1.* ...
  tdash+(-1).*f2.*tdash))+((-1).*ox+pxdash).*((c2+c1.*tdash+(-1).*c2.* ...
  tdash).*(e2+e1.*tdash+(-1).*e2.*tdash)+(-1).*(b2+b1.*tdash+(-1).*b2.* ...
  tdash).*(f2+f1.*tdash+(-1).*f2.*tdash)))+((-1).*c2.*((-1)+tdash).^2+c1.* ...
  tdash.^2).*((f2+f1.*tdash+(-1).*f2.*tdash).*((-1).*(b2+b1.*tdash+(-1).* ...
  b2.*tdash).^2+(a2+a1.*tdash+(-1).*a2.*tdash).*(d2+d1.*tdash+(-1).*d2.* ...
  tdash))+(-1).*(e2+e1.*tdash+(-1).*e2.*tdash).*((-1).*(b2+b1.*tdash+(-1) ...
  .*b2.*tdash).*(c2+c1.*tdash+(-1).*c2.*tdash)+(a2+a1.*tdash+(-1).*a2.* ...
  tdash).*(e2+e1.*tdash+(-1).*e2.*tdash))+(c2+c1.*tdash+(-1).*c2.*tdash).* ...
  ((-1).*(c2+c1.*tdash+(-1).*c2.*tdash).*(d2+d1.*tdash+(-1).*d2.*tdash)+( ...
  b2+b1.*tdash+(-1).*b2.*tdash).*(e2+e1.*tdash+(-1).*e2.*tdash))).^(-1).*( ...
  ((-1).*oz+pzdash).*((-1).*(c2+c1.*tdash+(-1).*c2.*tdash).*(d2+d1.*tdash+ ...
  (-1).*d2.*tdash)+(b2+b1.*tdash+(-1).*b2.*tdash).*(e2+e1.*tdash+(-1).* ...
  e2.*tdash))+((-1).*oy+pydash).*((c2+c1.*tdash+(-1).*c2.*tdash).*(e2+e1.* ...
  tdash+(-1).*e2.*tdash)+(-1).*(b2+b1.*tdash+(-1).*b2.*tdash).*(f2+f1.* ...
  tdash+(-1).*f2.*tdash))+((-1).*ox+pxdash).*((-1).*(e2+e1.*tdash+(-1).* ...
  e2.*tdash).^2+(d2+d1.*tdash+(-1).*d2.*tdash).*(f2+f1.*tdash+(-1).*f2.* ...
  tdash))));

end

