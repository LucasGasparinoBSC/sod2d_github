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
ndime = 2;

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
q = zeros(6,2);
pr = zeros(6,1);
u(1,:) = [1.0 1.0];
u(2,:) = [1.5 1.5];
u(3,:) = [1.0 1.0];
u(4,:) = [1.0 1.0];
u(5,:) = [1.5 1.5];
u(6,:) = [1.0 1.0];
q = 1.2*u;
pr(:,:) = 0.5;

% Gaussian integ. table

xgp = zeros(4,2);
xgp(1,:) = [-1.0/sqrt(3.0) -1.0/sqrt(3.0)];
xgp(2,:) = [1.0/sqrt(3.0) -1.0/sqrt(3.0)];
xgp(3,:) = [1.0/sqrt(3.0) 1.0/sqrt(3.0)];
xgp(4,:) = [-1.0/sqrt(3.0) 1.0/sqrt(3.0)];
ngaus = 4;

R = zeros(npoin,ndime);
for ielem = 1:nelem
  Re = zeros(nnode,ndime);
  ind = conn(ielem,:);           % Nodal index
  elq = q(ind,:);                % Element momentum
  elu = u(ind,:);                % Element velocity
  elpr = pr(ind,1);             % Element pressure
  J = xy(ind,:)'*transpose(dN); % Element Jacobian
  J = simplify(J);
  detJ = simplify(det(J));      % Jacobian determinant
  H = inv(J);                   % Jacobian inverse
  dNcar = simplify(H*dN);       % Cartesian derivatives
  aux = zeros(ndime,ndime);
  %
  % Compute symbolical divergence term
  %
  div = [s;s];
  div(:,:) = 0;
  for jnode = 1:nnode
    for idime = 1:ndime
      for jdime = 1:ndime
        aux = elq(jnode,idime)*elu(jnode,jdime);
        div(idime,1) = div(idime,1)+dNcar(jdime,jnode)*aux;
      endfor
    endfor
  endfor
  for igaus = 1:ngaus
    divgp = double(subs(simplify(div),{s t},{xgp(igaus,1) xgp(igaus,2)}));
    gpvol = double(subs(detJ,{s t},{xgp(igaus,1) xgp(igaus,2)}));
    Ngp = double(subs(simplify(N),{s t},{xgp(igaus,1) xgp(igaus,2)}));
    gpcar = double(subs(simplify(dNcar),{s t},{xgp(igaus,1) xgp(igaus,2)}));
    for inode = 1:nnode
      for idime = 1:ndime
        aux = 0;
        for jnode = 1:nnode
          aux = aux+gpcar(idime,jnode)*elpr(jnode,1);
        endfor
        Re(inode,idime) = Re(inode,idime) + gpvol*Ngp(inode)*(divgp(idime)+aux);
      endfor
    endfor
  endfor
  R(ind,1:ndime) = R(ind,1:ndime)+Re(:,1:ndime);
endfor