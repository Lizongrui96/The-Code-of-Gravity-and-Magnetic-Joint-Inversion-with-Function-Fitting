clc
clear
% close all
%% 读取数据
xyzg=load('zyg_c_7_6141.txt','-ascii');
% xyzg=load('g_lz_ys.txt','-ascii');
% cccc='213986';
cccc='49314';
cccc1=['MD' cccc '.mat'];
cccc2=['MDp' cccc '.mat'];
cccc3=['Gg' cccc '.mat'];
cccc4=['P' cccc '-443.mat'];
load(cccc1,'-mat');
load(cccc2,'-mat');
load(cccc3,'-mat');
load(cccc4,'-mat');
% load('Px2-242416.mat');
% load('Py2-242416.mat');
% load('Pz2-242416.mat');
% V=G*P;
nobs=length(xyzg(:,1));
N=length(MDp(:,1));
g=xyzg(:,4);
%% 输入参数
%迭代次数
nd=2000;
%% 反演
Wd0=sum(Gg.^2).^(0.5);
% C=inversion(V,g,nd,Wd0,P);
% [Dx1,Dy1,Dz1]=Difference1(nx,ny,nz);
% P=speye(49520);
C=inversion6(Gg,g,nd,Wd0,P);
Cg=C;
save(['Cgys1' cccc '.mat'],'Cg');
% C=inversion4(G,g,nd,Wd0,P,Dx1,Dy1,Dz1);
% C=inversion3(G,g,nd,Wd0,P);
%% 
% Gw=zeros(size(G));
% for i=1:nobs
%     Gw(i,:)=G(i,:).*Wd0;
% end
% Vw=Gw'*V;
% B=Gw'*g;
% C=inversion1(Vw,B,nd,Wd0,P);
m=P*Cg;
rr=Gg*m-g;
% nobs=length(G(:,1));
display(sqrt(rr'*rr/nobs));
%%
gg=Gg*m;
ggg=gg-g;
gggg=reshape(ggg,41,61);
figure(3)
imagesc(gggg)
colorbar
colormap jet
%%
res=zeros(N,4);
for i=1:N
    %% 坐标
    xyz=zeros(1,3);
    for ci=1:4
        xyz=xyz+MD(MDp(i,ci),:)/4;
    end
    res(i,:)=[xyz,m(i,1)];
end
save('g49520.txt','res','-ascii');
%%
% interp3
%%
% load('mod.mat');
% ck=[0,1,1];
% for i=1:length(mod)
%     l=length(mod(i).K(:,1));
%     %     trisurf(mod(i).K,mod(i).MD(:,1),mod(i).MD(:,2),mod(i).MD(:,3),'FaceAlpha',1,'FaceColor',[0,1,1],'LineWidth',0.1,'EdgeColor',[0,0.9,1])
%     for j=1:l
%         xk=[mod(i).MD(mod(i).K(j,1),1),mod(i).MD(mod(i).K(j,2),1),mod(i).MD(mod(i).K(j,3),1)];
%         yk=[mod(i).MD(mod(i).K(j,1),2),mod(i).MD(mod(i).K(j,2),2),mod(i).MD(mod(i).K(j,3),2)];
%         zk=[mod(i).MD(mod(i).K(j,1),3),mod(i).MD(mod(i).K(j,2),3),mod(i).MD(mod(i).K(j,3),3)];
%         fill3(xk,yk,zk,ck,'FaceAlpha',0,'LineWidth',2,'EdgeColor',[0 0 0])
%     end
%     hold on
% end



