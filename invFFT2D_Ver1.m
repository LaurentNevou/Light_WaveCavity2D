function [Vxyback] = invFFT2D_Ver1(Vk2Dcut,Ny,Nx)

Nkx=length(Vk2Dcut(1,:));
Nky=length(Vk2Dcut(:,1));

Nx1=Nx/2-floor(Nkx/2);
Nx2=Nx/2+ceil(Nkx/2);
Ny1=Ny/2-floor(Nky/2);
Ny2=Ny/2+ceil(Nky/2);

Vk2Dback=zeros(Ny,Nx);
size(Ny1+1:Ny2);
size(Nx1+1:Nx2);
size(Vk2Dcut);
size(Vk2Dback);

Vk2Dback( Ny1+1:Ny2 , Nx1+1:Nx2)=Vk2Dcut;

size(Vk2Dback);

Vk2=ifftshift(Vk2Dback);

Vxyback=ifft2(Vk2);

end