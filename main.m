clc
clear
close all
%%
load('P.mat');
load('Px.mat');
load('Py.mat');
load('Pz.mat');
load('Vg.mat');
load('Vm.mat');
load('Wzg.mat');
load('Wzm.mat');
xyzg=load('grav.txt','-ascii');
xyzm=load('mag.txt','-ascii');
%%
g=xyzg(:,4);
T=xyzm(:,4);
xmin=min(xyzg(:,1));
xmax=max(xyzg(:,1));
ymin=min(xyzg(:,2));
ymax=max(xyzg(:,2));
zmin=-10000;
zmax=0;
nx=60;
ny=40;
nz=20;
N=nx*ny*nz;
Itermax=1500;
lx=(xmax-xmin)/(nx);
ly=(ymax-ymin)/(ny);
lz=(zmax-zmin)/(nz);
x=(xmin+lx/2):lx:(xmax-lx/2);
y=(ymin+ly/2):ly:(ymax-ly/2);
z=(zmax-lz/2):-lz:(zmin+lz/2);
tic
[Cg,Cm]=inversion(Vg,Vm,g,T,Itermax,Wzg,Wzm,P,Px,Py,Pz);
toc
save('Cg.mat','Cg');
save('Cm.mat','Cm');
rho=0.001*P*Cg;
r0=Vg*Cg-g;
display(sqrt(r0'*r0/length(g)));
m=0.06.*P*Cm;
r0=Vm*Cm-T;
display(sqrt(r0'*r0/length(T)));
rr=reshape(rho,nz,ny,nx);
[xxx,yyy,zzz]=meshgrid(x,y,z);
%%
rrr=permute(rr,[2 3 1]);
figure(1)
a=slice(xxx,yyy,zzz,rrr,[],y(round(ny/2)),[]);
set(a,'EdgeColor','none')
colorbar
colormap jet
axis([min(x),max(x),min(y),max(y),min(z),max(z)])
shading interp
hold on
grid off
set(gca,'PlotBoxAspectRatio',[4 2 1]);
load('mod1.mat');
ck=[0,1,1];
for i=1:length(mod)
    l=length(mod(i).K(:,1));
    for j=1:l
        xk=[mod(i).MD(mod(i).K(j,1),1),mod(i).MD(mod(i).K(j,2),1),mod(i).MD(mod(i).K(j,3),1)];
        yk=[mod(i).MD(mod(i).K(j,1),2),mod(i).MD(mod(i).K(j,2),2),mod(i).MD(mod(i).K(j,3),2)];
        zk=[mod(i).MD(mod(i).K(j,1),3),mod(i).MD(mod(i).K(j,2),3),mod(i).MD(mod(i).K(j,3),3)];
        fill3(xk,yk,zk,ck,'FaceAlpha',0,'LineWidth',1,'EdgeColor',[0 0 0])
    end
    hold on
end
mm=reshape(m,nz,ny,nx);
[xxx,yyy,zzz]=meshgrid(x,y,z);
mmm=permute(mm,[2 3 1]);
figure(2)
a=slice(xxx,yyy,zzz,mmm,[],y(round(ny/2)),[]);
set(a,'EdgeColor','none')
colorbar
colormap jet
axis([min(x),max(x),min(y),max(y),min(z),max(z)])
shading interp
hold on
grid off
set(gca,'PlotBoxAspectRatio',[4 2 1]);
load('mod1.mat');
ck=[0,1,1];
for i=1:length(mod)-1
    l=length(mod(i).K(:,1));
    for j=1:l
        xk=[mod(i).MD(mod(i).K(j,1),1),mod(i).MD(mod(i).K(j,2),1),mod(i).MD(mod(i).K(j,3),1)];
        yk=[mod(i).MD(mod(i).K(j,1),2),mod(i).MD(mod(i).K(j,2),2),mod(i).MD(mod(i).K(j,3),2)];
        zk=[mod(i).MD(mod(i).K(j,1),3),mod(i).MD(mod(i).K(j,2),3),mod(i).MD(mod(i).K(j,3),3)];
        fill3(xk,yk,zk,ck,'FaceAlpha',0,'LineWidth',1,'EdgeColor',[0 0 0])
    end
    hold on
end









