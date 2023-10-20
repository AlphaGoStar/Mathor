close all
clear all
nmse_aod=0;
nmse_aoa=0;
nmse_dop=0;
for loop=1:1
    clearvars -except nmse_aod nmse_aoa nmse_dop loop
%% 参数配置
Nr=4;   %接收天线数量
Nt=Nr;  %发射天线数量
M=240;  %子载波数量
dic_size=800;   %字典大小
snr=30; %信噪比
Nu=2;   %用户数量
Lobj=4; %每个用户的多径数
Nrru=2;
Y_allpath=zeros(M,Nr);
numberofssb=4;  %SSB的数量，每个SSB有4个OFDM符号(此4是ssb数量而不是OFDM符号数量，符号数量为4*numberofssb
numberofsym=4;  %一个SSB使用的符号数量
Y=[];   %一条多径的接收信号
Y_allssb=[];    %全部ssb的接收信号,维度为子载波和一个ssb的符号长度
Y_allpath=zeros(M,Nr);   %全部多径的接收信号
Y_allpath_sym=[];       %四个ofdm符号排起来
act_aoa=(pi/2)*ones(1,Lobj*Nu)-(randperm(179,Lobj*Nu)/180)*pi;  %实际到达角
act_aod=(pi/2)*ones(1,Lobj*Nu)-(randperm(179,Lobj*Nu)/180)*pi;  %实际离去角
act_delay=randperm(dic_size-1,Lobj*Nu)/dic_size;    %实际的延迟（在网格）
act_dop=randperm(800,Lobj*Nu)/3200; %实际的多普勒
phase=[ 1 exp(1j*pi/(1)) exp(1j*pi/(2)) exp(1j*pi/3)];  %分别对准0度，±90度，60度，30度方向
for np=1:length(act_aoa)    %阵列响应向量矩阵
    Amtx(:,np)=exp(1i*pi*(0:Nt-1)*sin(act_aod(np))).'/sqrt(Nt);
    Amrx(:,np)=exp(1i*pi*(0:Nr-1)*sin(act_aoa(np))).'/sqrt(Nr);
end
qpsk=[-1 1];    %QPSK调制的字典
for row=1:M
    for column=1:numberofssb*numberofsym
        sig(row,column)=(qpsk(randi([1 2]))+qpsk(randi([1 2]))*1j)/sqrt(2); %经过QPSK调制的发射信号   M为子载波数，列数为符号数
    end
end
sig=1*sig;  %功率放大
%% 生成波束成形以及字典矩阵
for diff=1:2    %信号交替发送不同形式的信号
    for ssb=1:numberofssb
        for i=1:Nt
            w(ssb,i,diff)=phase(ssb)^i/sqrt(4); %不同ssb的波束成形矩阵
        end
        if(diff==1)
            w(ssb,Nt,diff)=0;   %奇数最后一个天线不发
        else
            w(ssb,2:Nt,diff)=w(ssb,1:Nt-1,diff);    %保持一致
            w(ssb,1,diff)=0;    %偶数第一个天线不发
        end
    end
end
for i=1:M
    dftmtx(i,:)=exp(-1i*2*pi*((0:dic_size-1)/dic_size)*i)/sqrt(M);  %延迟字典矩阵
end
%% 生成接收信号
for nssb=1:numberofssb  %SSB个数
    for nsy=1:numberofsym %一个ssb的符号数

        for k=1:M
            for nrru=1:Nrru
                for np=1:Lobj*Nu    %多径数
                    if mod(nrru,2)==0
                        Y(k,:)=sig(k,(nssb-1)*numberofsym+nsy)*conj(w(nssb,:,1+(mod(nsy,2)==0)))*Amtx(:,np)*exp(-1i*k*act_delay(np)*2*pi)*(exp(1i*((nssb-1)*numberofsym+nsy)*2*pi*act_dop(np)))*Amrx(:,np).';    %第K个子载波上的第n个多径的接收信号

                    else

                        Y(k,:)=sig(k,(nssb-1)*numberofsym+nsy)*w(nssb,:,1+(mod(nsy,2)==0))*Amtx(:,np)*exp(-1i*k*act_delay(np)*2*pi)*(exp(1i*((nssb-1)*numberofsym+nsy)*2*pi*act_dop(np)))*Amrx(:,np).';    %第K个子载波上的第n个多径的接收信号
                    end
                    Y_allpath(k,:)=Y_allpath(k,:)+Y(k,:);  %多个多径信号加起来
                end
            end
        end

        Y_allpath=1./sig(:,(nssb-1)*numberofsym+nsy).*Y_allpath./4;
        Y_allpath_sym=[Y_allpath_sym Y_allpath];    %4个ofdm符号排起来
        Y_allpath=zeros(M,Nr);   %初始化Y_allpath
        Y_allpath=zeros(M,Nr);  %初始化Y_all_rru
    end
    Y_allssb=[Y_allssb Y_allpath_sym];  %多个SSB信号排起来
    Y_allpath_sym=[];       %初始化Y_allpath_sym
end
Y_addnoise=awgn(Y_allssb,snr);  %加噪声
%% 估计
smu=[];
for i=1:numberofssb
    Y=zeros(M,numberofssb*Nr)+Y_allssb(:,(i-1)*numberofssb*Nr+1:i*numberofssb*Nr);
[smu]=[smu uamp_sbl_mmv(dftmtx,Y_addnoise)];    %使用UAMP算法
end

vapsilon=smu.*conj(smu);    %计算能量
power=mean(abs(vapsilon),2);    %计算能量
est_delay=find(abs(power)>1e-1);  %寻找符合要求的索引
for block=1:numberofssb
    block_power(:,block)=sum(smu(est_delay,(block-1)*Nr*numberofsym+1:block*Nr*numberofsym).*conj(smu(est_delay,(block-1)*Nr*numberofsym+1:block*Nr*numberofsym)),2);   %找出每一行能量最大的那一块，即波束对的准的那一个SSB
end
[~,index]=max(block_power,[],2);    %找出最大的那一块
% 到达角估计
for i=1:length(est_delay)
    for j=1:numberofsym
        s(j)=sum(smu(est_delay(i),(index(i)-1)*Nr*numberofsym+(j-1)*Nr+1:(index(i)-1)*Nr*numberofsym+((j-1)*Nr+1+Nr-2)).*conj(smu(est_delay(i),(index(i)-1)*Nr*numberofsym+(j-1)*Nr+2:(index(i)-1)*Nr*numberofsym+((j-1)*Nr+2+Nr-2))));    %块内共轭相加提高信噪比
    end
    est_aoa(i)=-angle(sum(s))/pi;
    s(:)=0; %初始化，进入下一次估计循环
end
s(:)=0; %清空，为下次使用
% 多普勒估计
for i=1:length(est_delay)
    for j=1:fix(numberofsym/2)
        s(j)=sum(smu(est_delay(i),(index(i)-1)*Nr*numberofsym+(j-1)*Nr+1:(index(i)-1)*Nr*numberofsym+j*Nr).*conj(smu(est_delay(i),(index(i)-1)*Nr*numberofsym+(j+1)*Nr+1:(index(i)-1)*Nr*numberofsym+(j+2)*Nr)));   %跨符号作共轭充分利用数据提高信噪比
    end
    est_dop(i)=-angle(sum(s))/(4*pi); 
end
% 发射角估计
for i=1:length(est_delay)
    used_sym=[smu(est_delay(i),(index(i)-1)*Nr*numberofsym+1:(index(i)-1)*Nr*numberofsym+Nr) ...
        smu(est_delay(i),(index(i)-1)*Nr*numberofsym+Nr+1:(index(i)-1)*Nr*numberofsym+Nr*2).*exp(-1j*est_dop(i)*(2*pi)) ...
        smu(est_delay(i),(index(i)-1)*Nr*numberofsym+2*Nr+1:(index(i)-1)*Nr*numberofsym+Nr*3).*exp(-1j*est_dop(i)*(4*pi)) ... %乘上多余的多普勒，使其除了发射角项不一样以外，其他都一样
        smu(est_delay(i),(index(i)-1)*Nr*numberofsym+3*Nr+1:(index(i)-1)*Nr*numberofsym+Nr*4).*exp(-1j*est_dop(i)*(6*pi))];
    for j=1:numberofsym-1
        if(mod(j,2)~=0)
            s(j)=sum(conj(used_sym((j-1)*Nr+1:j*Nr).*conj(used_sym(j*Nr+1:(j+1)*Nr))));    %共轭求到达角
        else
            s(j)=sum(used_sym((j-1)*Nr+1:j*Nr).*conj(used_sym(j*Nr+1:(j+1)*Nr)));
        end
    end
    est_aod(i)=angle(sum(s))/pi;
    s(:)=0;         %初始化，进入下一次估计循环
    used_sym=[];    %初始化，进入下一次估计循环
end
%% 画图
% 延迟-到达角
figure
scatter(sin(act_aoa),act_delay*dic_size,'r','o')   %真实延迟放大至网格
hold on;
scatter(est_aoa,est_delay-1,'r','+')
xlabel('AOA')
ylabel('Delay ')
% 到达角-多普勒
figure
scatter(sin(act_aoa),act_dop*100,'r','o')   
hold on;
scatter(est_aoa,est_dop*100,'r','+')
xlabel('AOA (Radian)')
ylabel('dop ')
%离去角-多普勒
figure
scatter(sin(act_aod),act_dop*100,'r','o')   
hold on;
scatter(est_aod,est_dop*100,'r','+')
xlabel('AOD (Radian)')
ylabel('dop ')
figure
scatter(sin(act_aod),act_delay*dic_size,'r','o')   
hold on;
scatter(est_aod,est_delay-1,'r','+')
xlabel('AOD (Radian)')
ylabel('dop ')
nmse_aod=nmse_aod+norm(sort(abs(est_aod))-sort(abs(sin(act_aod))))/norm(abs(sin(act_aod))); %发射角的NMSE
nmse_dop=nmse_dop+norm(sort(abs(est_dop))-sort(abs(sin(act_dop))))/norm(abs(sin(act_dop))); %多普勒的NMSE  
nmse_aoa=nmse_aoa+norm(sort(abs(est_aoa))-sort(abs(sin(act_aoa))))/norm(abs(sin(act_aoa))); %到达角的NMSE
loop    %显示进度 
end
nmse_dop=nmse_dop/loop;
nmse_aod=nmse_aod/loop;
nmse_aoa=nmse_aoa/loop;