clear;
clc;

pkg load symbolic

syms s t

% shape functions and gradN

N = (1-s)*(1-t);
N(2) = (1+s)*(1-t);
N(3) = (1+s)*(1+t);
N(4) = (1-s)*(1+t);
N = (1/4)*N;

dN = diff(N,s);
dN(2,:) = diff(N,t);

% Mesh info

conn = zeros(2,4);
conn(1,:) = [1 2 5 6];
conn(2,:) = [2 3 4 5];
nelem = 2;
nnode = 4;

xy = zeros(6,2);
xy(1,:) = [0.0 0.0];
xy(2,:) = [2.0 0.0];
xy(3,:) = [3.0 0.0];
xy(4,:) = [3.5 3.0];
xy(5,:) = [2.0 2.0];
xy(6,:) = [0.0 2.0];
npoin = 6;

% Energy flux data

u = zeros(6,2);
E = zeros(6,1);
pr = zeros(6,1);
u(1,:) = [1.0 1.0];
u(2,:) = [1.5 1.5];
u(3,:) = [1.0 1.0];
u(4,:) = [1.0 1.0];
u(5,:) = [1.5 1.5];
u(6,:) = [1.0 1.0];
E = 0.5;
pr = 0.5;
q = u*(E+pr);

% Gaussian integ. table

xgp = zeros(4,2);
xgp(1,:) = [-1.0/sqrt(3.0) -1.0/sqrt(3.0)];
xgp(2,:) = [1.0/sqrt(3.0) -1.0/sqrt(3.0)];
xgp(3,:) = [1.0/sqrt(3.0) 1.0/sqrt(3.0)];
xgp(4,:) = [-1.0/sqrt(3.0) 1.0/sqrt(3.0)];
ngaus = 4;

% Jacobian info

J1 = transpose(xy(conn(1,:),:))*transpose(dN);
detJ1 = det(J1);

J2 = transpose(xy(conn(2,:),:))*transpose(dN);
detJ2 = det(J2);

% Cartesian derivatives

dxdN1 = simplify(inv(J1)*dN);
dxdN2 = simplify(inv(J2)*dN);

% Integration

R = zeros(6,1);
R1 = zeros(4,1);
R2 = zeros(4,1);

for igaus = 1:ngaus
  Ngp = double(subs(N,{s t},{xgp(igaus,1) xgp(igaus,2)}));
  dxdN1gp = double(subs(dxdN1,{s t},{xgp(igaus,1) xgp(igaus,2)}));
  dxdN2gp = double(subs(dxdN2,{s t},{xgp(igaus,1) xgp(igaus,2)}));
  detJ1gp = double(subs(detJ1,{s t},{xgp(igaus,1) xgp(igaus,2)}));
  detJ2gp = double(subs(detJ2,{s t},{xgp(igaus,1) xgp(igaus,2)}));
  R1 = R1+detJ1gp*transpose(Ngp)*(dxdN1gp(1,:)*q(conn(1,:),1)+dxdN1gp(2,:)*q(conn(1,:),2));
  R2 = R2+detJ2gp*transpose(Ngp)*(dxdN2gp(1,:)*q(conn(2,:),1)+dxdN2gp(2,:)*q(conn(2,:),2));
end

R(conn(1,:),1) = R(conn(1,:),1) + R1;
R(conn(2,:),1) = R(conn(2,:),1) + R2;
