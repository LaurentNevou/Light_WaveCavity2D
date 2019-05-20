function[Ez,f0_z]=WC2D_TM_Ez_FEM_f(x,y,eps,nmodes,f0_guess,f0_min,f0_max)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

c=2.99792458e8;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Nx=length(x);
Ny=length(y);
Nxy=Nx*Ny;

dx=x(2)-x(1);
dy=y(2)-y(1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%% Building of the operators %%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%% First derivative X %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% DX1(Ny=5,Nx=4); %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                                  %
%   0   0   0   0   0 |1/2  0   0   0   0 | 0   0   0   0   0 | 0   0   0   0   0  %
%   0   0   0   0   0 | 0  1/2  0   0   0 | 0   0   0   0   0 | 0   0   0   0   0  %
%   0   0   0   0   0 | 0   0  1/2  0   0 | 0   0   0   0   0 | 0   0   0   0   0  %
%   0   0   0   0   0 | 0   0   0  1/2  0 | 0   0   0   0   0 | 0   0   0   0   0  %
%   0   0   0   0   0 | 0   0   0   0  1/2| 0   0   0   0   0 | 0   0   0   0   0  %
%   -----------------------------------------------------------------------------  %
% -1/2  0   0   0   0 | 0   0   0   0   0 |1/2  0   0   0   0 | 0   0   0   0   0  %
%   0 -1/2  0   0   0 | 0   0   0   0   0 | 0  1/2  0   0   0 | 0   0   0   0   0  %
%   0   0 -1/2  0   0 | 0   0   0   0   0 | 0   0  1/2  0   0 | 0   0   0   0   0  %
%   0   0   0 -1/2  0 | 0   0   0   0   0 | 0   0   0  1/2  0 | 0   0   0   0   0  %
%   0   0   0   0 -1/2| 0   0   0   0   0 | 0   0   0   0  1/2| 0   0   0   0   0  %
%   -----------------------------------------------------------------------------  %
%   0   0   0   0   0 -1/2  0   0   0   0 | 0   0   0   0   0 |1/2  0   0   0   0  %
%   0   0   0   0   0 | 0 -1/2  0   0   0 | 0   0   0   0   0 | 0  1/2  0   0   0  %
%   0   0   0   0   0 | 0   0 -1/2  0   0 | 0   0   0   0   0 | 0   0  1/2  0   0  %
%   0   0   0   0   0 | 0   0   0 -1/2  0 | 0   0   0   0   0 | 0   0   0  1/2  0  %
%   0   0   0   0   0 | 0   0   0   0 -1/2| 0   0   0   0   0 | 0   0   0   0  1/2 %
%   -----------------------------------------------------------------------------  %
%   0   0   0   0   0 | 0   0   0   0   0 -1/2  0   0   0   0 | 0   0   0   0   0  %
%   0   0   0   0   0 | 0   0   0   0   0 | 0 -1/2  0   0   0 | 0   0   0   0   0  %
%   0   0   0   0   0 | 0   0   0   0   0 | 0   0 -1/2  0   0 | 0   0   0   0   0  %
%   0   0   0   0   0 | 0   0   0   0   0 | 0   0   0 -1/2  0 | 0   0   0   0   0  %
%   0   0   0   0   0 | 0   0   0   0   0 | 0   0   0   0 -1/2| 0   0   0   0   0  %
%                                                                                  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%% Second derivative X %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% DX2(Ny=5,Nx=4); %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                                  %
%  -2   0   0   0   0 | 1   0   0   0   0 | 0   0   0   0   0 | 0   0   0   0   0  %
%   0  -2   0   0   0 | 0   1   0   0   0 | 0   0   0   0   0 | 0   0   0   0   0  %
%   0   0  -2   0   0 | 0   0   1   0   0 | 0   0   0   0   0 | 0   0   0   0   0  %
%   0   0   0  -2   0 | 0   0   0   1   0 | 0   0   0   0   0 | 0   0   0   0   0  %
%   0   0   0   0  -2 | 0   0   0   0   1 | 0   0   0   0   0 | 0   0   0   0   0  %
%   -----------------------------------------------------------------------------  %
%   1   0   0   0   0 |-2   0   0   0   0 | 1   0   0   0   0 | 0   0   0   0   0  %
%   0   1   0   0   0 | 0  -2   0   0   0 | 0   1   0   0   0 | 0   0   0   0   0  %
%   0   0   1   0   0 | 0   0  -2   0   0 | 0   0   1   0   0 | 0   0   0   0   0  %
%   0   0   0   1   0 | 0   0   0  -2   0 | 0   0   0   1   0 | 0   0   0   0   0  %
%   0   0   0   0   1 | 0   0   0   0  -2 | 0   0   0   0   1 | 0   0   0   0   0  %
%   -----------------------------------------------------------------------------  %
%   0   0   0   0   0 | 1   0   0   0   0 |-2   0   0   0   0 | 1   0   0   0   0  %
%   0   0   0   0   0 | 0   1   0   0   0 | 0  -2   0   0   0 | 0   1   0   0   0  %
%   0   0   0   0   0 | 0   0   1   0   0 | 0   0  -2   0   0 | 0   0   1   0   0  %
%   0   0   0   0   0 | 0   0   0   1   0 | 0   0   0  -2   0 | 0   0   0   1   0  %
%   0   0   0   0   0 | 0   0   0   0   1 | 0   0   0   0  -2 | 0   0   0   0   1  %
%   -----------------------------------------------------------------------------  %
%   0   0   0   0   0 | 0   0   0   0   0 | 1   0   0   0   0 |-2   0   0   0   0  %
%   0   0   0   0   0 | 0   0   0   0   0 | 0   1   0   0   0 | 0  -2   0   0   0  %
%   0   0   0   0   0 | 0   0   0   0   0 | 0   0   1   0   0 | 0   0  -2   0   0  %
%   0   0   0   0   0 | 0   0   0   0   0 | 0   0   0   1   0 | 0   0   0  -2   0  %
%   0   0   0   0   0 | 0   0   0   0   0 | 0   0   0   0   1 | 0   0   0   0  -2  %
%                                                                                  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Axy=ones(1,(Nx-1)*Ny);

DX1 =     (-0.5)*diag(Axy,-Ny) + (+0.5)*diag(Axy,Ny);
DX2 = (-2)*diag(ones(1,Ny*Nx)) + (1)*diag(Axy,-Ny) + (1)*diag(Axy,Ny);

%%%%%%%%%%%%%%%%%%%%%%%%%%%% First derivative Y %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% DY1(Ny=5,Nx=4); %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                                  %
%   0  1/2  0   0   0 | 0   0   0   0   0 | 0   0   0   0   0 | 0   0   0   0   0  %
% -1/2  0  1/2  0   0 | 0   0   0   0   0 | 0   0   0   0   0 | 0   0   0   0   0  %
%   0 -1/2  0  1/2  0 | 0   0   0   0   0 | 0   0   0   0   0 | 0   0   0   0   0  %
%   0   0 -1/2  0  1/2| 0   0   0   0   0 | 0   0   0   0   0 | 0   0   0   0   0  %
%   0   0   0 -1/2  0 | 0   0   0   0   0 | 0   0   0   0   0 | 0   0   0   0   0  %
%   -----------------------------------------------------------------------------  %
%   0   0   0   0   0 | 0  1/2  0   0   0 | 0   0   0   0   0 | 0   0   0   0   0  %
%   0   0   0   0   0 -1/2  0  1/2  0   0 | 0   0   0   0   0 | 0   0   0   0   0  %
%   0   0   0   0   0 | 0 -1/2  0  1/2  0 | 0   0   0   0   0 | 0   0   0   0   0  %
%   0   0   0   0   0 | 0   0 -1/2  0  1/2| 0   0   0   0   0 | 0   0   0   0   0  %
%   0   0   0   0   0 | 0   0   0 -1/2  0 | 0   0   0   0   0 | 0   0   0   0   0  %
%   -----------------------------------------------------------------------------  %
%   0   0   0   0   0 | 0   0   0   0   0 | 0  1/2  0   0   0 | 0   0   0   0   0  %
%   0   0   0   0   0 | 0   0   0   0   0 -1/2  0  1/2  0   0 | 0   0   0   0   0  %
%   0   0   0   0   0 | 0   0   0   0   0 | 0 -1/2  0  1/2  0 | 0   0   0   0   0  %
%   0   0   0   0   0 | 0   0   0   0   0 | 0   0 -1/2  0  1/2| 0   0   0   0   0  %
%   0   0   0   0   0 | 0   0   0   0   0 | 0   0   0 -1/2  0 | 0   0   0   0   0  %
%   -----------------------------------------------------------------------------  %
%   0   0   0   0   0 | 0   0   0   0   0 | 0   0   0   0   0 | 0  1/2  0   0   0  %
%   0   0   0   0   0 | 0   0   0   0   0 | 0   0   0   0   0 -1/2  0  1/2  0   0  %
%   0   0   0   0   0 | 0   0   0   0   0 | 0   0   0   0   0 | 0 -1/2  0  1/2  0  %
%   0   0   0   0   0 | 0   0   0   0   0 | 0   0   0   0   0 | 0   0 -1/2  0  1/2 %
%   0   0   0   0   0 | 0   0   0   0   0 | 0   0   0   0   0 | 0   0   0 -1/2  0  %
%                                                                                  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%% Second derivative Y %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% DY2(Ny=5,Nx=4); %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                                  %
%  -2   1   0   0   0 | 0   0   0   0   0 | 0   0   0   0   0 | 0   0   0   0   0  %
%   1  -2   1   0   0 | 0   0   0   0   0 | 0   0   0   0   0 | 0   0   0   0   0  %
%   0   1  -2   1   0 | 0   0   0   0   0 | 0   0   0   0   0 | 0   0   0   0   0  %
%   0   0   1  -2   1 | 0   0   0   0   0 | 0   0   0   0   0 | 0   0   0   0   0  %
%   0   0   0   1  -2 | 0   0   0   0   0 | 0   0   0   0   0 | 0   0   0   0   0  %
%   -----------------------------------------------------------------------------  %
%   0   0   0   0   0 |-2   1   0   0   0 | 0   0   0   0   0 | 0   0   0   0   0  %
%   0   0   0   0   0 | 1  -2   1   0   0 | 0   0   0   0   0 | 0   0   0   0   0  %
%   0   0   0   0   0 | 0   1  -2   1   0 | 0   0   0   0   0 | 0   0   0   0   0  %
%   0   0   0   0   0 | 0   0   1  -2   1 | 0   0   0   0   0 | 0   0   0   0   0  %
%   0   0   0   0   0 | 0   0   0   1  -2 | 0   0   0   0   0 | 0   0   0   0   0  %
%   -----------------------------------------------------------------------------  %
%   0   0   0   0   0 | 0   0   0   0   0 |-2   1   0   0   0 | 0   0   0   0   0  %
%   0   0   0   0   0 | 0   0   0   0   0 | 1  -2   1   0   0 | 0   0   0   0   0  %
%   0   0   0   0   0 | 0   0   0   0   0 | 0   1  -2   1   0 | 0   0   0   0   0  %
%   0   0   0   0   0 | 0   0   0   0   0 | 0   0   1  -2   1 | 0   0   0   0   0  %
%   0   0   0   0   0 | 0   0   0   0   0 | 0   0   0   1  -2 | 0   0   0   0   0  %
%   -----------------------------------------------------------------------------  %
%   0   0   0   0   0 | 0   0   0   0   0 | 0   0   0   0   0 |-2   1   0   0   0  %
%   0   0   0   0   0 | 0   0   0   0   0 | 0   0   0   0   0 | 1  -2   1   0   0  %
%   0   0   0   0   0 | 0   0   0   0   0 | 0   0   0   0   0 | 0   1  -2   1   0  %
%   0   0   0   0   0 | 0   0   0   0   0 | 0   0   0   0   0 | 0   0   1  -2   1  %
%   0   0   0   0   0 | 0   0   0   0   0 | 0   0   0   0   0 | 0   0   0   1  -2  %
%                                                                                  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

AA=ones(1,Nx*Ny);
BB=ones(1,Nx*Ny-1);
BB(Ny:Ny:end)=0;

DY1=diag(BB,-1)*(-0.5) + diag(BB,1)*(0.5);
DY2=(-2)*diag(AA) + (1)*diag(BB,-1) + (1)*diag(BB,1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%% Building of the Hamiltonien %%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

H0 = DX2/dx^2 + DY2/dy^2 ;

Hz = -diag(1./eps(:)) * H0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%% Diagonalisation of the Hamiltonien %%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Guess = (2*pi*f0_guess/c)^2;

Hz=sparse(Hz);
[psi_z,k0_z] = eigs(Hz,nmodes,'SM');   %% eigen values are ordered
%[psi_z,k0_z] = eigs(Hz,nmodes,Guess);

f0_z=sqrt(diag(k0_z)) * c /2/pi;
lambda_z= 2*pi ./ sqrt(diag(k0_z)) * 1e6;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%% Filtering and reshaping the Wavefunction %%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

idx1=real(f0_z)>f0_min;
idx2=real(f0_z)<f0_max;
%idx3=imag(f0_z)==0;
%idx=logical( idx1.*idx2.*idx3);
idx=logical( idx1.*idx2);

f0_z=f0_z(idx);
psi_z=psi_z(:,idx);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% here is a small patch due to differences between Octave and Matlab
% Matlab order the eigen values while Octave reverse it

if f0_z(2)<f0_z(1)
  psi_z=psi_z(:,end:-1:1);
  f0_z=f0_z(end:-1:1);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Ez = psi_z ;

nmodes=length(f0_z);
Ez = reshape(Ez,[Ny,Nx,nmodes]);

end
