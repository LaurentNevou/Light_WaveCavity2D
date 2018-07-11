%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%% last update 23Mai2018, lne %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all
close all
clc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Constants %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

c=2.99792458e8;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% model activation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 0 for turn off
% 1 for turn on

FE_Method=1;            % Diagonalization of the Hamiltonian (FEM)
PWE_Method=1;           % Plane Wave Expansion (PWE)

% Take care, the PWE method doesn t take into account the Poisson equation, div(D)=ro
% while the FEM does. Therefore, solutions are slightly different in TE

TE=1;
TM=0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%% Solving Parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

f0_guess= c/5e-6;     %% Guess of the frequency solutions (Hz), f=c/lambda
f0_min  = c/50e-6;    %% filter the solutions where the frequency is superior than (Hz), f=c/lambda
f0_max  = c/0.8e-6;   %% filter the solutions where the frequency is inferior than (Hz), f=c/lambda
nmodes=15;            %% number of solutions asked 

AbsorbingBoundaryCondition=0;     %% 0 or 1 (not sure it is working well...)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%% Plotting parameter %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

AAbs=0;               %% Plot abs(E)
RReal=1;              %% Plot real(E)
IImag=0;              %% Plot imag(E)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%% Optical index definition %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Two vectors (x and y) and one 2D-matrix n must be defined with homogeneous grid
% x, y [meter]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Nx=51;                  %% Meshing point in x-direction
Ny=43;                  %% Meshing point in y-direction

Dx=1E-6;                %% map X [m]
Dy=1E-6;                %% map Y [m]

x = linspace(-Dx, Dx, Nx);
y = linspace(-Dy, Dy, Ny);

dx = x(2)-x(1);
dy = y(2)-y(1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Choose between the next 2 structures or build your own on an homogeneous grid!

n1=1; n2=3;

Lx=1.3e-6; Ly=1.2e-6;
[n,eps]=epsBox_f(x,y,Lx,Ly,n1,n2,AbsorbingBoundaryCondition);

%Rx=0.5e-6; Ry=0.6e-6;
%[n,eps]=epsElipse_f(x,y,Rx,Ry,n1,n2,AbsorbingBoundaryCondition);

%break
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%% NOTHING TO CHANGE ANYMORE!!! %%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%% Selection of the model %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
f1_xy=[];f2_xy=[];f1_z=[];f2_z=[];Ex1=[];Ey1=[];Ez1=[];Exy2=[];Ez2=[];
display('=======================================')

if FE_Method==1
    
    if length(x)*length(y)*2>1e4
      Nxy=length(x)*length(y);
      display(strcat('Warning: Take care, H=',num2str(Nxy),'x',num2str(Nxy),'x2 elements'))
    end
    if TE==1
      tic
      [Ex1,Ey1,f1_xy] = WC2D_TE_Exy_FEM_f(x,y,eps,nmodes,f0_guess,f0_min,f0_max);
      display(strcat('-> Finite Elements method TE=',num2str(toc),'sec'))
    end
    if TM==1
      tic
      [Ez1,f1_z]     = WC2D_TM_Ez_FEM_f(x,y,eps,nmodes,f0_guess,f0_min,f0_max);
      display(strcat('-> Finite Elements method TM=',num2str(toc),'sec'))
    end
    
end

if PWE_Method==1
    Nx = 64 ;        % number of points on the x grid % has to be a power of 2 (32,64,128,256,512,...) (smaller => faster)
    Ny = 64 ;        % number of points on the y grid % has to be a power of 2 (32,64,128,256,512,...) (smaller => faster)
    NGx = 15;%Nx/2-1  ;    % number of harmonics % must be at least 2 times -1 smaller than Nx (smaller => faster)
    NGy = 15;%Ny/2-1  ;    % number of harmonics % must be at least 2 times -1 smaller than Ny (smaller => faster)
    
    if TE==1
      tic
      [Exy2,f2_xy] = WC2D_TE_Exy_PWE_f(x,y,eps,nmodes,f0_guess,f0_min,f0_max,Nx,Ny,NGx,NGy);
      display(strcat('-> PWE method TE=',num2str(toc),'sec'))
    end
    if TM==1
      tic
      [Ez2,f2_z]     = WC2D_TM_Ez_PWE_f(x,y,eps,nmodes,f0_guess,f0_min,f0_max,Nx,Ny,NGx,NGy);
      display(strcat('-> PWE method TM=',num2str(toc),'sec'))
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%% Display Results %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if TE==1
  fTE=nan(nmodes,2);
  fTE(1:length(f1_xy),1)=f1_xy;
  fTE(1:length(f2_xy),2)=f2_xy;

  lambdaTE=nan(nmodes,2);
  lambdaTE(1:length(f1_xy),1)=c./f1_xy;
  lambdaTE(1:length(f2_xy),2)=c./f2_xy;
  
  display('=======================================')
  display('TE: Results lambda(um):')
  lambdaTE*1e6
  %display('=======================================')
  %display('TE: Results f(THz):')
  %fTE*1e-12
  
  if AAbs==1
    EEx1=abs(Ex1);EEy1=abs(Ey1);
    EExy2=abs(Exy2);
  end
  if RReal==1
    EEx1=real(Ex1);EEy1=real(Ey1);
    EExy2=real(Exy2);
  end
  if IImag==1
    EEx1=imag(Ex1);EEy1=imag(Ey1);
    EExy2=imag(Exy2);
  end
end

if TM==1
  fTM=nan(nmodes,2);
  fTM(1:length(f1_z),1)=f1_z;
  fTM(1:length(f2_z),2)=f2_z;

  lambdaTM=nan(nmodes,2);
  lambdaTM(1:length(f1_z),1)=c./f1_z;
  lambdaTM(1:length(f2_z),2)=c./f2_z;
  display('=======================================')
  display('TM: Results lambda(um):')
  lambdaTM*1e6
  %display('=======================================')
  %display('TM: Results f(THz):')
  %fTM*1e-12
  
  if AAbs==1
    EEz1=abs(Ez1);EEz2=abs(Ez2);
  end
  if RReal==1
    EEz1=real(Ez1);EEz2=real(Ez2);  
  end
  if IImag==1
    EEz1=imag(Ez1);EEz2=imag(Ez2);
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% figures %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

X0fig = 50    ; Y0fig = 50;
Wfig  = 1500  ; Hfig  = 800;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if FE_Method==1 && TE==1
  
    figure('Name',strcat('FEM method: Electrical Field'),'position',[X0fig Y0fig Wfig Hfig])
    colormap(jet)
    
for i=1:6%length(neff)
    
    subplot(3,6,i)
    hold on
    A=EEx1(:,:,i);
    if sum(A(:)==0)
      pcolor(x*1e6,y*1e6,A)
    else
      pcolor(x*1e6,y*1e6,A/max(abs(A(:))))
    end
    contour(x*1e6,y*1e6,abs(n),1,'linewidth',2,'linecolor','w')
    
    if RReal==1 || IImag==1
        caxis([-1 1])
    elseif AAbs==1
        caxis([0 1])
    end
    
    if AAbs==1
      title(  strcat('\lambda=',  num2str( lambdaTE(i,1)*1e6 ,'%.2f'),'um,  TE:Abs(Ex)'  )    )
    elseif RReal==1
      title(  strcat('\lambda=',  num2str( lambdaTE(i,1)*1e6 ,'%.2f'),'um,  TE:Re(Ex)'  )    )
    elseif IImag==1
      title(  strcat('\lambda=',  num2str( lambdaTE(i,1)*1e6 ,'%.2f'),'um,  TE:Im(Ex)'  )    )
    end
    shading flat
    xlabel('x (um)')
    ylabel('y (um)')
    colorbar
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    subplot(3,6,i+6)
    hold on
    A=EEy1(:,:,i);
    if sum(A(:)==0)
      pcolor(x*1e6,y*1e6,A)
    else
      pcolor(x*1e6,y*1e6,A/max(abs(A(:))))
    end
    contour(x*1e6,y*1e6,abs(n),1,'linewidth',2,'linecolor','w')
    
    if RReal==1 || IImag==1
        caxis([-1 1])
    elseif AAbs==1
        caxis([0 1])
    end
    
    if AAbs==1
      title(  strcat('\lambda=',  num2str( lambdaTE(i,1)*1e6 ,'%.2f'),'um,  TE:Abs(Ey)'  )    )
    elseif RReal==1
      title(  strcat('\lambda=',  num2str( lambdaTE(i,1)*1e6 ,'%.2f'),'um,  TE:Re(Ey)'  )    )
    elseif IImag==1
      title(  strcat('\lambda=',  num2str( lambdaTE(i,1)*1e6 ,'%.2f'),'um,  TE:Im(Ey)'  )    )
    end
    shading flat
    xlabel('x (um)')
    ylabel('y (um)')
    colormap(jet)
    colorbar
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    subplot(3,6,i+12)
    hold on
    A=sqrt(EEx1(:,:,i).^2 + EEy1(:,:,i).^2 );
    if sum(A(:)==0)
      pcolor(x*1e6,y*1e6,A)
    else
      pcolor(x*1e6,y*1e6,A/max(abs(A(:))))
    end
    contour(x*1e6,y*1e6,abs(n),1,'linewidth',2,'linecolor','w')
    caxis([0 1])
    
    if AAbs==1
      title(  strcat('Abs(sqrt(Ex^2+Ey^2))'  )    )
    elseif RReal==1
      title(  strcat('Re(sqrt(Ex^2+Ey^2))'  )    )
    elseif IImag==1
      title(  strcat('Im(sqrt(Ex^2+Ey^2))'  )    )
    end
    shading flat
    xlabel('x (um)')
    ylabel('y (um)')
    colorbar
    
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if FE_Method==1 && TM==1
  
    figure('Name',strcat('FEM method: Electrical Field'),'position',[X0fig Y0fig Wfig Hfig])
    colormap(jet)
    
for i=1:length(lambdaTM)
    
    subplot(3,6,i)
    hold on
    A=EEz1(:,:,i);
    if sum(A(:)==0)
      pcolor(x*1e6,y*1e6,A)
    else
      pcolor(x*1e6,y*1e6,A/max(abs(A(:))))
    end
    contour(x*1e6,y*1e6,abs(n),1,'linewidth',2,'linecolor','w')
    
    if RReal==1 || IImag==1
        caxis([-1 1])
    elseif AAbs==1
        caxis([0 1])
    end
    
    if AAbs==1
      title(  strcat('\lambda=',  num2str( lambdaTM(i,1)*1e6 ,'%.2f'),'um,  TM:Abs(Ez)'  )    )
    elseif RReal==1
      title(  strcat('\lambda=',  num2str( lambdaTM(i,1)*1e6 ,'%.2f'),'um,  TM:Re(Ez)'  )    )
    elseif IImag==1
      title(  strcat('\lambda=',  num2str( lambdaTM(i,1)*1e6 ,'%.2f'),'um,  TM:Im(Ez)'  )    )
    end
    shading flat
    xlabel('x (um)')
    ylabel('y (um)')
    colorbar
    
end
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if PWE_Method==1 && TE==1
  
    figure('Name',strcat('PWE method: Electrical Field'),'position',[X0fig Y0fig Wfig Hfig])
    colormap(jet)
    
for i=1:6%length(neff)
    
    subplot(3,6,i)
    hold on
    A=EExy2(:,:,i);
    if sum(A(:)==0)
      pcolor(x*1e6,y*1e6,A)
    else
      pcolor(x*1e6,y*1e6,A/max(abs(A(:))))
    end
    contour(x*1e6,y*1e6,abs(n),1,'linewidth',2,'linecolor','w')
    
    if RReal==1 || IImag==1
        caxis([-1 1])
    elseif AAbs==1
        caxis([0 1])
    end
    
    if AAbs==1
      title(  strcat('\lambda=',  num2str( lambdaTE(i,2)*1e6 ,'%.2f'),'um,  TE:Abs(Exy)'  )    )
    elseif RReal==1
      title(  strcat('\lambda=',  num2str( lambdaTE(i,2)*1e6 ,'%.2f'),'um,  TE:Re(Exy)'  )    )
    elseif IImag==1
      title(  strcat('\lambda=',  num2str( lambdaTE(i,2)*1e6 ,'%.2f'),'um,  TE:Im(Exy)'  )    )
    end
    shading flat
    xlabel('x (um)')
    ylabel('y (um)')
    colorbar   
    
end
end

if PWE_Method==1 && TM==1
  
    figure('Name',strcat('PWE method: Electrical Field'),'position',[X0fig Y0fig Wfig Hfig])
    colormap(jet)
    
for i=1:6%length(neff)
    
    subplot(3,6,i)
    hold on
    A=EEz2(:,:,i);
    if sum(A(:)==0)
      pcolor(x*1e6,y*1e6,A)
    else
      pcolor(x*1e6,y*1e6,A/max(abs(A(:))))
    end
    contour(x*1e6,y*1e6,abs(n),1,'linewidth',2,'linecolor','w')
    
    if RReal==1 || IImag==1
        caxis([-1 1])
    elseif AAbs==1
        caxis([0 1])
    end
    
    if AAbs==1
      title(  strcat('\lambda=',  num2str( lambdaTM(i,2)*1e6 ,'%.2f'),'um,  TM:Abs(Ez)'  )    )
    elseif RReal==1
      title(  strcat('\lambda=',  num2str( lambdaTM(i,2)*1e6 ,'%.2f'),'um,  TM:Re(Ez)'  )    )
    elseif IImag==1
      title(  strcat('\lambda=',  num2str( lambdaTM(i,2)*1e6 ,'%.2f'),'um,  TM:Im(Ez)'  )    )
    end
    shading flat
    xlabel('x (um)')
    ylabel('y (um)')
    colorbar
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    
end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

