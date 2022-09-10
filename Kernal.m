clc
clear
close all
xyzg=load('grav.txt','-ascii');
outg='Gg.mat';
outm='Gm.mat';
outvg='Vg.mat';
outvm='Vm.mat';
outWg='Wzg.mat';
outWm='Wzm.mat';
nobs=length(xyzg(:,1));
xmin=min(xyzg(:,1));
xmax=max(xyzg(:,1));
ymin=min(xyzg(:,2));
ymax=max(xyzg(:,2));
zmin=-10000;
zmax=0;
nx1=12;
ny1=8;
nz1=1;
nx2=5;
ny2=5;
nz2=20;
nx=nx1*nx2;
ny=ny1*ny2;
nz=nz1*nz2;
N1=nx1*ny1*nz1;
N2=nx2*ny2*nz2;
N=nx*ny*nz;
s=3;
II0=60;
DD0=10;
II=60;
DD=10;
lx=(xmax-xmin)/(nx);
ly=(ymax-ymin)/(ny);
lz=(zmax-zmin)/(nz);
size0=[lx,ly,lz];
x=(xmin+lx/2):lx:(xmax-lx/2);
y=(ymin+ly/2):ly:(ymax-ly/2);
z=(zmax-lz/2):-lz:(zmin+lz/2);
Gg=zeros(nobs,N);
Gm=zeros(nobs,N);
Wz0=zeros(1,N);
tt=1;
xyz=zeros(N,3);
for i=1:nx
    for j=1:ny
        for k=1:nz
            xyz(tt,:)=[x(i),y(j),z(k)];
            tt=tt+1;
        end
    end
end
for t=1:nobs
    xyz0=xyzg(t,1:3);
    for tt=1:N
        Gg(t,tt)=grav_fun('Vz',xyz0,xyz(tt,:),size0);%÷ÿ¡¶
        Gm(t,tt)=mag_fun('Vt',xyz0,xyz(tt,:),size0,II,-DD-90,II0,-DD0-90);%¥≈%
    end
    disp(t);
end
for tt=1:N
    Wz0(1,tt)=(-xyz(tt,3)-zmin);
end
Wdg=sum(Gg.^2).^(0.5);
Wdm=sum(Gm.^2).^(0.5);
Wzg=sparse(1:N,1:N,Wdg);
Wzm=sparse(1:N,1:N,Wdm);
save(outg,'Gg');
save(outm,'Gm');
%%
st=0;
for tx=0:s
    for ty=0:s-tx
        for tz=0:s-tx-ty
            st=st+1;
        end
    end
end
P=zeros(N,N1*st);
Px=zeros(N,N1*st);
Py=zeros(N,N1*st);
Pz=zeros(N,N1*st);
%% P
tt=0;
for i=1:nx1
    for ii=1:nx2
        for j=1:ny1
            for jj=1:ny2
                for k=1:nz1
                    for kk=1:nz2
                        tt=tt+1;
                        t=k+(j-1)*nz1+(i-1)*nz1*ny1-1;
                        ttt=t*st+1;
                        for tx=0:s
                            for ty=0:s-tx
                                for tz=0:s-tx-ty
                                    P(tt,ttt)=(((ii-nx2/2-0.5))^tx)*(((jj-ny2/2-0.5))^ty)*(((kk-nz2))^tz);
                                    ttt=ttt+1;
                                end
                            end
                        end
                    end
                end
            end
        end
    end
end
P=sparse(P);
%% Px
tt=0;
for i=1:nx1
    for ii=1:nx2
        for j=1:ny1
            for jj=1:ny2
                for k=1:nz1
                    for kk=1:nz2
                        tt=tt+1;
                        t=k+(j-1)*nz1+(i-1)*nz1*ny1-1;
                        ttt=t*st+1;
                        for tx=0:s
                            for ty=0:s-tx
                                for tz=0:s-tx-ty
                                    if(tx==0)
                                        Px(tt,ttt)=0;
                                    else
                                        Px(tt,ttt)=tx*(((ii-nx2/2-0.5))^(tx-1))*(((jj-ny2/2-0.5))^ty)*(((kk-nz2))^tz);
                                    end
                                    ttt=ttt+1;
                                end
                            end
                        end
                    end
                end
            end
        end
    end
end
Px=sparse(Px);
%% Py
tt=0;
for i=1:nx1
    for ii=1:nx2
        for j=1:ny1
            for jj=1:ny2
                for k=1:nz1
                    for kk=1:nz2
                        tt=tt+1;
                        t=k+(j-1)*nz1+(i-1)*nz1*ny1-1;
                        ttt=t*st+1;
                        for tx=0:s
                            for ty=0:s-tx
                                for tz=0:s-tx-ty
                                    if(ty==0)
                                        Py(tt,ttt)=0;
                                    else
                                        Py(tt,ttt)=ty*(((ii-nx2/2-0.5))^tx)*(((jj-ny2/2-0.5))^(ty-1))*(((kk-nz2))^tz);
                                    end
                                    ttt=ttt+1;
                                end
                            end
                        end
                    end
                end
            end
        end
    end
end
Py=sparse(Py);
%% Pz
tt=0;
for i=1:nx1
    for ii=1:nx2
        for j=1:ny1
            for jj=1:ny2
                for k=1:nz1
                    for kk=1:nz2
                        tt=tt+1;
                        t=k+(j-1)*nz1+(i-1)*nz1*ny1-1;
                        ttt=t*st+1;
                        for tx=0:s
                            for ty=0:s-tx
                                for tz=0:s-tx-ty
                                    if(tz==0)
                                        Pz(tt,ttt)=0;
                                    else
                                        Pz(tt,ttt)=tz*(((ii-nx2/2-0.5))^tx)*(((jj-ny2/2-0.5))^ty)*(((kk-nz2))^(tz-1));
                                    end
                                    ttt=ttt+1;
                                end
                            end
                        end
                    end
                end
            end
        end
    end
end
Pz=sparse(Pz);
save('P.mat','P');
save('Px.mat','Px');
save('Py.mat','Py');
save('Pz.mat','Pz');
Vg=Gg*P;
Vm=Gm*P;
save(outvg,'Vg');
save(outvm,'Vm');
% Wzg=sparse(1:N,1:N,Wz0.^-1);
% Wzm=sparse(1:N,1:N,Wz0.^-1.5);
save(outWg,'Wzg');
save(outWm,'Wzm');


% a=diag(Wzg1);
% aa=diag(Wzg);