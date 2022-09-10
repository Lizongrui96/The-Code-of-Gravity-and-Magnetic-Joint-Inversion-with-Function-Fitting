clc
clear
close all
xmin=0;
xmax=30000;
ymin=0;
ymax=20000;
nx=61;
ny=41;
lx=(xmax-xmin)/(nx-1);
ly=(ymax-ymin)/(ny-1);
x=xmin:lx:xmax;
y=ymin:ly:ymax;
[xx,yy]=meshgrid(x,y);
zz=0*cos(4*pi*(xx-(xmax+xmin)/2)/(xmax-xmin)).*cos(pi*(yy-(ymax+ymin)/2)/(ymax-ymin));%地形为平面，且=0
xyzg=[xx(:),yy(:),zz(:),zz(:)*0];
nxt=61;
nyt=41;
lxt=(xmax-xmin)/(nxt-1);
lyt=(ymax-ymin)/(nyt-1);
xt=xmin:lxt:xmax;
yt=ymin:lyt:ymax;
[xxt,yyt]=meshgrid(xt,yt);
zzt=0*cos(4*pi*(xxt-(xmax+xmin)/2)/(xmax-xmin)).*cos(pi*(yyt-(ymax+ymin)/2)/(ymax-ymin));
xyzt=[xxt(:),yyt(:),zzt(:)];
outg='grav.txt';
outc='mag.txt';
outtp='topo.txt';
num=3;
II0=60;
DD0=10;
nnn=1;
mod(nnn).xyz=[7500,10000,-4000];
mod(nnn).size0=[4000,5000,4000];
t=0;
for i=-1:2:1
    for j=-1:2:1
        for k=-1:2:1
            t=t+1;
            mod(nnn).MD(t,:)=mod(nnn).xyz+[i,(j),k].*mod(nnn).size0/2;
        end
    end
end
mod(nnn).rho=0.5;
mod(nnn).M=5;
mod(nnn).II=60;
mod(nnn).DD=10;
mod(nnn).s=delaunayTriangulation(mod(nnn).MD);
[mod(nnn).K,~] = convexHull(mod(nnn).s);
nnn=2;
mod(nnn).xyz=[15000,10000,-5000];
mod(nnn).size0=[5000,5000,5000];
t=0;
for i=-1:2:1
    for j=-1:2:1
        for k=-1:2:1
            t=t+1;
            mod(nnn).MD(t,:)=mod(nnn).xyz+[i,j,k].*mod(nnn).size0/2;
        end
    end
end
mod(nnn).rho=0.5;
mod(nnn).M=5;
mod(nnn).II=60;
mod(nnn).DD=10;
mod(nnn).s=delaunayTriangulation(mod(nnn).MD);
[mod(nnn).K,~] = convexHull(mod(nnn).s);
nnn=3;
mod(nnn).xyz=[23500,10000,-4000];
mod(nnn).size0=[4000,5000,4000];
t=0;
for i=-1:2:1
    for j=-1:2:1
        for k=-1:2:1
            t=t+1;
            mod(nnn).MD(t,:)=mod(nnn).xyz+[i,j,k].*mod(nnn).size0/2;
        end
    end
end
mod(nnn).rho=0.5;
mod(nnn).M=0;
mod(nnn).II=60;
mod(nnn).DD=10;
mod(nnn).s=delaunayTriangulation(mod(nnn).MD);
[mod(nnn).K,~] = convexHull(mod(nnn).s);
save('mod1.mat','mod');
xyzm=xyzg;
xyzm(:,4)=0;
for k=1:num
    for i=1:length(xyzg(:,1))
            xyz0=xyzg(i,1:3);
            xyzg(i,4)=xyzg(i,4)+mod(k).rho*1000*grav_fun('Vz',xyz0,mod(k).xyz,mod(k).size0);
            xyzm(i,4)=xyzm(i,4)+mod(k).M*mag_fun('Vt',xyz0,mod(k).xyz,mod(k).size0,mod(k).II,-mod(k).DD-90,II0,-DD0-90);
    end
end
g=reshape(xyzg(:,4),ny,nx);
gm=reshape(xyzm(:,4),ny,nx);
figure(1)
surf(xx,yy,zz,g)
alpha(0.8)
colorbar
colormap jet
hold on
shading interp
ck=[1,0,0];
for i=1:num
    l=length(mod(i).K(:,1));
    for j=1:l
        xk=[mod(i).MD(mod(i).K(j,1),1),mod(i).MD(mod(i).K(j,2),1),mod(i).MD(mod(i).K(j,3),1)];
        yk=[mod(i).MD(mod(i).K(j,1),2),mod(i).MD(mod(i).K(j,2),2),mod(i).MD(mod(i).K(j,3),2)];
        zk=[mod(i).MD(mod(i).K(j,1),3),mod(i).MD(mod(i).K(j,2),3),mod(i).MD(mod(i).K(j,3),3)];
        fill3(xk,yk,zk,ck,'FaceAlpha',1,'LineWidth',2,'EdgeColor',[0,0,0])
        fill3(xk,yk,-10000*ones(size(zk)),ck,'FaceAlpha',0,'LineWidth',2,'EdgeColor',[0,0,0])
        fill3(xk,ymax*ones(size(yk)),zk,ck,'FaceAlpha',0,'LineWidth',2,'EdgeColor',[0,0,0])
        fill3(xmax*ones(size(xk)),yk,zk,ck,'FaceAlpha',0,'LineWidth',2,'EdgeColor',[0,0,0])
        for k=1:8
            xxk=mod(i).MD(k,1);yyk=mod(i).MD(k,2);zzk=mod(i).MD(k,3);
            plot3([xxk,xmax],[yyk,yyk],[zzk,zzk],'--','Color',[0.4,0.4,0.4],'LineWidth',1.5)
            plot3([xxk,xxk],[yyk,ymax],[zzk,zzk],'--','Color',[0.4,0.4,0.4],'LineWidth',1.5)
            plot3([xxk,xxk],[yyk,yyk],[-10000,zzk],'--','Color',[0.4,0.4,0.4],'LineWidth',1.5)
        end
    end
    hold on
end
axis([xmin xmax ymin ymax -10000 max(max(zz))])
set(gca,'PlotBoxAspectRatio',[3 2 1]);
box on
grid on
%%
figure(777)
surf(xx,yy,zz,gm)
alpha(0.8)
colorbar
colormap jet
hold on
shading interp
ck=[1,0,0];
for i=1:num-1
    l=length(mod(i).K(:,1));
    for j=1:l
        xk=[mod(i).MD(mod(i).K(j,1),1),mod(i).MD(mod(i).K(j,2),1),mod(i).MD(mod(i).K(j,3),1)];
        yk=[mod(i).MD(mod(i).K(j,1),2),mod(i).MD(mod(i).K(j,2),2),mod(i).MD(mod(i).K(j,3),2)];
        zk=[mod(i).MD(mod(i).K(j,1),3),mod(i).MD(mod(i).K(j,2),3),mod(i).MD(mod(i).K(j,3),3)];
        fill3(xk,yk,zk,ck,'FaceAlpha',1,'LineWidth',2,'EdgeColor',[0,0,0])
        fill3(xk,yk,-10000*ones(size(zk)),ck,'FaceAlpha',0,'LineWidth',2,'EdgeColor',[0,0,0])
        fill3(xk,ymax*ones(size(yk)),zk,ck,'FaceAlpha',0,'LineWidth',2,'EdgeColor',[0,0,0])
        fill3(xmax*ones(size(xk)),yk,zk,ck,'FaceAlpha',0,'LineWidth',2,'EdgeColor',[0,0,0])
        for k=1:8
            xxk=mod(i).MD(k,1);yyk=mod(i).MD(k,2);zzk=mod(i).MD(k,3);
            plot3([xxk,xmax],[yyk,yyk],[zzk,zzk],'--','Color',[0.4,0.4,0.4],'LineWidth',1.5)
            plot3([xxk,xxk],[yyk,ymax],[zzk,zzk],'--','Color',[0.4,0.4,0.4],'LineWidth',1.5)
            plot3([xxk,xxk],[yyk,yyk],[-10000,zzk],'--','Color',[0.4,0.4,0.4],'LineWidth',1.5)
        end
    end
    hold on
end
axis([xmin xmax ymin ymax -10000 max(max(zz))])
set(gca,'PlotBoxAspectRatio',[3 2 1]);
box on
grid on
%%
save(outg,'xyzg','-ascii');
save(outc,'xyzm','-ascii');
save(outtp,'xyzt','-ascii');

























