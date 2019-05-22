function[n,eps]=epsBox_f(x,y,Lx,Ly,n1,n2,AbsorbingBoundaryCondition)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%% Optical index definition %%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[X,Y]=meshgrid(x,y);

x0=0;y0=0;

idx  = 1 > abs((X-x0)/Lx*2);
idy  = 1 > abs((Y-y0)/Ly*2);
idXY = idx.*idy ;

n = n2*idXY + n1*(1-idXY);  %% ridge optical index


if AbsorbingBoundaryCondition==1
    LOSS=1e-4;
    n(:,1)                = n(:,1)               + LOSS*i;
    n(:,end)              = n(:,end)             + LOSS*i;
    n(1,2:end-1)          = n(1,2:end-1)         + LOSS*i;
    n(end,2:end-1)        = n(end,2:end-1)       + LOSS*i;
end

eps=n.^2;

%break
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure('Name','Optical index','position',[10 50 1600 800])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

subplot(1,1,1,'fontsize',15)
hold on;grid on;

pcolor(x*1e6,y*1e6,abs(squeeze(n)))

colormap(jet)
colorbar

%xlim([-1 1]*Dx/2*1e6)
%ylim([-1 1]*Dz/2*1e6)

xlabel('x (um)')
ylabel('y (um)')
