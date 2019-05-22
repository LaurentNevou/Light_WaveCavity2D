function[n,eps]=epsRing_f(x,y,Rx_min,Ry_min,Rx_max,Ry_max,n1,n2,AbsorbingBoundaryCondition)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%% Optical index definition %%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[X,Y]=meshgrid(x,y);

x0=0;y0=0;


idXY1 =  ((X-x0)/Rx_min).^2 + ((Y-y0)/Ry_min).^2  < 1  ;  %% internal elipse
idXY2 =  ((X-x0)/Rx_max).^2 + ((Y-y0)/Ry_max).^2  < 1  ;  %% external elipse

n = n2*idXY2 + n1*(1-idXY2) ;  
n(idXY1)=n1;


if AbsorbingBoundaryCondition==1
    LOSS=1e-2;
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

figure('Name','Optical index','position',[10 -100 1000 700])

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
