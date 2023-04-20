clear all
%%%%%%%%%%%%% 前言 %%%%%%%%%%%%%%%%%
%此代码由安徽工业大学电气信息工程学院汪某人基于详解matlab通信系统建模与仿真（电子工业出版社，刘学勇）源码基础上修改而来，从而实现基于卷积码，qpsk，ofdm的noma通信系统设计

%%%%%%%%%%%%% 参数设置部分 %%%%%%%%%%%%%%%%%

Nsp=52;             %系统子载波数（不包括直流载波）
Nfft=64;            % FFT 长度
Ncp=16;             % 循环前缀长度
Ns=Nfft+Ncp;        % 1个完整OFDM符号长度
noc=53;             % 包含直流载波的总的子载波数
Nd=6;               % 每帧包含的OFDM符号数(不包括训练符号)
M1=4;               % QPSK调制
sr=250000;          % OFDM符号速率
EbNo=0:2:30;      	% 归一化信噪比
Nfrm=10000;                         % 每种信噪比下的仿真帧数
%多普勒效应 暂未考虑
ts=1/sr/Ns;                         % OFDM符号抽样时间间隔
t=0:ts:(Ns*(Nd+1)*Nfrm-1)*ts;       % 抽样时刻
fd=10;                             % 最大多普勒频移
%卷积码参数设定
L=7;                %卷积码约束长度
tblen=6*L;          %Viterbi译码器回溯深度
stage = 3;          % m序列的阶数
ptap1 = [1 3];      % m序列的寄存器连接方式
regi1 = [1 1 1];    % m序列的寄存器初始值
% 两用户功率参数
Rp_db=0;                                      %用户1与用户2的功率分配比dB值
Rp=10^(Rp_db/10);                              %用户1与用户2的功率分配比，Rp=p_1/p_2,（假设用户1分得的功率比用户2的大）
p_u1=Rp/(1+Rp);                                %p_1与p_2分别为用户1与用户2的功率分配系数，并假设 p_1+p_2=1
p_u2=1/(1+Rp);
%% 瑞利信道生成
h=rayleigh(fd,t);                   % 生成单径Rayleigh衰落信道
h1=sqrt(2/3)*h;
h2=sqrt(1/3)*rayleigh(fd,t);
h2=[zeros(1,4) h2(1:end-4)];


%%
%训练符号频域数据,采用802.11a中的长训练符号数据
Preamble=[1 1 -1 -1 1 1 -1 1 -1 1 1 1 1 1 1 -1 -1 1 1 -1 1 -1 1 1 1 1 ...
    1 -1 -1 1 1 -1 1 -1 1 -1 -1 -1 -1 -1 1 1 -1 -1 1 -1 1 -1 1 1 1 1];     %52比特数据
Preamble1=zeros(1,Nfft);
Preamble1(2:27)=Preamble(27:end);                   % 前导重排后的数据
Preamble1(39:end)=Preamble(1:26);

preamble1=ifft(Preamble1);                          % 训练符号的时域数据
preamble1=[preamble1(Nfft-Ncp+1:end) preamble1];    % 训练符号加入循环前缀

%% %%%%%%%%%%%%%%%%%%% 仿真循环 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for ii=1:length(EbNo)
    %% 用户产生随机序列
    data_tx1=randsrc(1,Nsp*Nd*Nfrm,[0:1]);       %用户3
    data_tx2=randsrc(1,Nsp*Nd*Nfrm,[0:1]);       %用户4
    data_tx1yanzhen=reshape(data_tx1,52,[]);     %统计误码率需使用的变量
    data_tx2yanzhen=reshape(data_tx2,52,[]);
    %% 信道编码（卷积码、或交织器）
%卷积码：前向纠错非线性码
%交织：使突发错误最大限度的分散化
    trellis = poly2trellis(7,[133 171]);            % (2,1,7)卷积编码，码率1/2
    code_data=convenc(data_tx1,trellis);            % 输入（1*3120000），输出（1*6240000）  
    trellis = poly2trellis(7,[133 171]);    
    code_data2=convenc(data_tx2,trellis);   
    %% 验证卷积码是否正确？ 正确的话可以删除
    trellis = poly2trellis(7,[133 171]);
    yanzheng_data = vitdec(code_data,trellis,tblen,'trunc','hard');                 
    yanzheng_shuchu=reshape(yanzheng_data,52,[]);
    %验证正确√
%% qpsk调制
   % data1=pskmod(msg1,M1,pi/4);                       % 用户1QPSK调制
    %data2=pskmod(msg2,M1,pi/4);                       % 用户2QPSK调制
%用户3数据转换与qpsk调制
    data_temp1=reshape(code_data,2,[])';              % 
    data_temp2=bi2de(data_temp1);                     % 二进制转化为十进制，方便下一步qpsk调制
    data_temp3=pskmod(data_temp2,M1,pi/M1);           % 4PSK调制，y = pskmod(x,M,ini_phase) specifies the initial phase of the PSK-modulated signal.
    data_temp4=reshape(data_temp3,52,[]) ;            % 用户3（卷积码）qpsk调制，生成卷积码过后的数据52*60000，马上进行重排
%用户4数据转换与qpsk调制
    data_temp11=reshape(code_data2,2,[])';            %
    data_temp22=bi2de(data_temp11);                   % 
    data_temp33=pskmod(data_temp22,M1,pi/M1);         % 
    data_temp44=reshape(data_temp33,52,[]) ;        
    
    
    

    %星座图验证
%    data_xingzuotu1=reshape(data1,1,[]);
  %  data_xingzuotu2=reshape(data2,1,[]);
   % data_xingzuotu3=reshape(data_temp4,1,[]);
    %% 用户叠加
    
    data_diejia=sqrt(p_u1)*data_temp4+sqrt(p_u2)*data_temp44;
    %% ifft数据重排
    %data3=zeros(Nfft,Nd*Nfrm);                  % 根据FFT要求，对数据重排
    %data4=zeros(Nfft,Nd*Nfrm);                  
    data5=zeros(Nfft,Nd*Nfrm); 
    data6=zeros(Nfft,Nd*Nfrm);
    data7=zeros(Nfft,Nd*Nfrm);
    %用户1重排数据
   % data3(2:27,:)=data1(27:end,:);              % qpsk的前导后重排的数据
   % data3(39:end,:)=data1(1:26,:);              % 进行data1的重排数据
    %用户2重排数据
   % data4(2:27,:)=data2(27:end,:);              % qpsk的前导后重排的数据
   % data4(39:end,:)=data2(1:26,:);              % 进行data1的重排数据
   
    %用户3（使用卷积码）重排数据
    data5(2:27,:)=data_temp4(27:end,:);         % qpsk的前导后重排的数据
    data5(39:end,:)=data_temp4(1:26,:);         % 进行data1的重排数据  
    %用户4（使用卷积码）重排数据
    data6(2:27,:)=data_temp44(27:end,:);         % qpsk的前导后重排的数据
    data6(39:end,:)=data_temp44(1:26,:);         % 进行data1的重排数据  
    %叠加信号（使用卷积码）重排数据
    data7(2:27,:)=data_diejia(27:end,:);         % qpsk的前导后重排的数据
    data7(39:end,:)=data_diejia(1:26,:);         % 进行data1的重排数据  
    
   % clear data1 data2 data11;                  % 清除不需要的临时变量
   
    %data3=ifft(data3);                          % 用户1IFFT变换
    %data4=ifft(data4);                          % 用户2TFFT变换
    data5=ifft(data5);                          % 用户3（卷积码）ifft变换
    data6=ifft(data6);                          % 用户4（卷积码）ifft变换
    data7=ifft(data7);
   
    %% ifft之后加入循环前缀
   % data3=[data3(Nfft-Ncp+1:end,:);data3];      % 加入循环前缀，循环前缀长度为16，插入后矩阵变为80*60000
    %data4=[data4(Nfft-Ncp+1:end,:);data4];
    data5=[data5(Nfft-Ncp+1:end,:);data5];
    data6=[data6(Nfft-Ncp+1:end,:);data6];
    data7=[data7(Nfft-Ncp+1:end,:);data7];

    %spow1=norm(data3,'fro').^2/(Nsp*Nd*Nfrm);   %norm计算qpsk数据矩阵的范数，计算符号能量
    %spow2=norm(data4,'fro').^2/(Nsp*Nd*Nfrm);   %norm计算qpsk数据矩阵的范数，计算符号能量
    spow3=norm(data5,'fro').^2/(Nsp*Nd*Nfrm);   %norm计算qpsk数据矩阵的范数，计算符号能量
    spow4=norm(data6,'fro').^2/(Nsp*Nd*Nfrm);   %norm计算qpsk数据矩阵的范数，计算符号能量
    spow5=norm(data7,'fro').^2/(Nsp*Nd*Nfrm);   %norm计算qpsk数据矩阵的范数，计算符号能量
    
    %data10=zeros(Ns,(Nd+1)*Nfrm);
    %data11=zeros(Ns,(Nd+1)*Nfrm);
    data12=zeros(Ns,(Nd+1)*Nfrm);
    data13=zeros(Ns,(Nd+1)*Nfrm);
    data14=zeros(Ns,(Nd+1)*Nfrm);
    %% 每隔7帧加入训练符号
    for indx=1:Nfrm
       
        %data10(:,(indx-1)*(Nd+1)+1)=preamble1.';                                 %用户1添加训练序列
        %data10(:,(indx-1)*(Nd+1)+2:indx*(Nd+1))=data3(:,(indx-1)*Nd+1:indx*Nd);  %用户1添加符号信息
        
        %data11(:,(indx-1)*(Nd+1)+1)=preamble1.';                                 %用户2    
        %data11(:,(indx-1)*(Nd+1)+2:indx*(Nd+1))=data4(:,(indx-1)*Nd+1:indx*Nd);  %
        
        data12(:,(indx-1)*(Nd+1)+1)=preamble1.';                                 %用户3（卷积码）    
        data12(:,(indx-1)*(Nd+1)+2:indx*(Nd+1))=data5(:,(indx-1)*Nd+1:indx*Nd);  %
        
        data13(:,(indx-1)*(Nd+1)+1)=preamble1.';                                 %用户4（卷积码）    
        data13(:,(indx-1)*(Nd+1)+2:indx*(Nd+1))=data6(:,(indx-1)*Nd+1:indx*Nd);  %
        
        data14(:,(indx-1)*(Nd+1)+1)=preamble1.';                                 %叠加信号（卷积码）    
        data14(:,(indx-1)*(Nd+1)+2:indx*(Nd+1))=data7(:,(indx-1)*Nd+1:indx*Nd);  %
    end
        
                                
    %% 并串转换
   % data10=reshape(data10,1,Ns*(Nd+1)*Nfrm);               % 用户1
    %data11=reshape(data11,1,Ns*(Nd+1)*Nfrm);               % 用户2
    data12=reshape(data12,1,Ns*(Nd+1)*Nfrm);               % 用户3
    data13=reshape(data13,1,Ns*(Nd+1)*Nfrm);               % 用户4
    data14=reshape(data14,1,Ns*(Nd+1)*Nfrm);               % 叠加用户
     
    %data51=zeros(1,length(data10));                     %用户1
    %data511=zeros(1,length(data11));                    %用户2
    data71=zeros(1,length(data12));                     %用户3（卷积码）
    data81=zeros(1,length(data13));                     %用户4（卷积码）           
    data91=zeros(1,length(data14));                     %叠加用户
    
    %data51(5:end)=data10(1:end-4);                          %第二径比第一径延迟四个采样点
    %data511(5:end)=data11(1:end-4);
    data71(5:end)=data12(1:end-4);
    data81(5:end)=data13(1:end-4);
    data91(5:end)=data14(1:end-4);
    %%  根据Eb/No计算噪声标准差  
   % sigma1=sqrt(1/2*spow1/log2(M1)*10.^(-EbNo(ii)/10)); 
    %sigma2=sqrt(1/2*spow2/log2(M1)*10.^(-EbNo(ii)/10));
    sigma3=sqrt(1/2*spow3/log2(M1)*10.^(-EbNo(ii)/10));
    sigma4=sqrt(1/2*spow4/log2(M1)*10.^(-EbNo(ii)/10));
    sigma5=sqrt(1/2*spow5/log2(M1)*10.^(-EbNo(ii)/10));
    %% 信道
    for indx=1:Nfrm
        %dd_user1=data10((indx-1)*Ns*(Nd+1)+1:indx*Ns*(Nd+1));
        %dd_user2=data11((indx-1)*Ns*(Nd+1)+1:indx*Ns*(Nd+1));
        dd_user3=data12((indx-1)*Ns*(Nd+1)+1:indx*Ns*(Nd+1));
        dd_user4=data13((indx-1)*Ns*(Nd+1)+1:indx*Ns*(Nd+1));
        dd_user5=data14((indx-1)*Ns*(Nd+1)+1:indx*Ns*(Nd+1));
        
        %dd_user11=data51((indx-1)*Ns*(Nd+1)+1:indx*Ns*(Nd+1));
        %dd_user22=data511((indx-1)*Ns*(Nd+1)+1:indx*Ns*(Nd+1));
        dd_user33=data71((indx-1)*Ns*(Nd+1)+1:indx*Ns*(Nd+1));
        dd_user44=data81((indx-1)*Ns*(Nd+1)+1:indx*Ns*(Nd+1));
        dd_user55=data91((indx-1)*Ns*(Nd+1)+1:indx*Ns*(Nd+1));
        
        % 当前帧的单径信道参数
        hh=h((indx-1)*Ns*(Nd+1)+1:indx*Ns*(Nd+1));      %暂时不用，本系统在两径信道产生
        % 当前帧的两径信道参数
        hh1=h1((indx-1)*Ns*(Nd+1)+1:indx*Ns*(Nd+1));    
        hh2=h2((indx-1)*Ns*(Nd+1)+1:indx*Ns*(Nd+1));
        
        % 信号通过2径衰落信道，并加入高斯白噪声
        %r11=hh1.*dd_user1+hh2.*dd_user11+sigma1*(randn(1,length(dd_user1))+1i*randn(1,length(dd_user1)));   
        %r31=hh1.*dd_user2+hh2.*dd_user22+sigma2*(randn(1,length(dd_user2))+1i*randn(1,length(dd_user2))); 
        r41=hh1.*dd_user3+hh2.*dd_user33+sigma3*(randn(1,length(dd_user3))+1i*randn(1,length(dd_user3))); 
        r51=hh1.*dd_user4+hh2.*dd_user44+sigma4*(randn(1,length(dd_user4))+1i*randn(1,length(dd_user4))); 
        r61=hh1.*dd_user5+hh2.*dd_user55+sigma5*(randn(1,length(dd_user5))+1i*randn(1,length(dd_user5))); 
        %% 接收端串并变换
        %俩径信道传播
        %r11=reshape(r11,Ns,Nd+1);
        %r31=reshape(r31,Ns,Nd+1);
        r41=reshape(r41,Ns,Nd+1);
        r51=reshape(r51,Ns,Nd+1);
        r61=reshape(r61,Ns,Nd+1);
     
                
        % 双径信道移除循环前缀
       % r11=r11(Ncp+1:end,:);  
        %r31=r31(Ncp+1:end,:);  
        r41=r41(Ncp+1:end,:);
        r51=r51(Ncp+1:end,:);
        r61=r61(Ncp+1:end,:);
        %% 接收端 fft运算

        % 双径信道fft运算    
        %R11=fft(r11);
       % R31=fft(r31);
        R41=fft(r41);
        R51=fft(r51);
        R61=fft(r61);
        % 双径信道数据重排
        %R11=[R11(39:end,:);R11(2:27,:)];
        %R31=[R31(39:end,:);R31(2:27,:)];
        R41=[R41(39:end,:);R41(2:27,:)];
        R51=[R51(39:end,:);R51(2:27,:)];
        R61=[R61(39:end,:);R61(2:27,:)];
       %% 接收端 信道估计
        % 信道估计 
       % HH11=(Preamble.')./R11(:,1);            
       % HH31=(Preamble.')./R31(:,1);
        HH41=(Preamble.')./R41(:,1);
        HH51=(Preamble.')./R51(:,1);
        HH61=(Preamble.')./R61(:,1);
       
        
        %HH11=HH11*ones(1,Nd);                       %用户1信道估计
        %HH31=HH31*ones(1,Nd);                       %用户2信道估计
        HH41=HH41*ones(1,Nd);                       %用户3（卷积码）
        HH51=HH51*ones(1,Nd); 
        HH61=HH61*ones(1,Nd); 
       %% 信道补偿     
        %x1=R11(:,2:end).*HH11;                      %用户1的2径信道补偿       
       % x2=R31(:,2:end).*HH31;                      %用户2的2径信道补偿
        x3=R41(:,2:end).*HH41;  
        x4=R51(:,2:end).*HH51; 
        x5=R61(:,2:end).*HH61; 
        
        %% OMA数据解调
        %x1=pskdemod(x1,M1,pi/4);  
        %x2=pskdemod(x2,M1,pi/4);
        x3=pskdemod(x3,M1,pi/4);
        x4=pskdemod(x4,M1,pi/4);
        x5_3=pskdemod(x5,M1,pi/4);
       %%  NOMA解调  用户3
        %用户3qpsk解调完成
      	x5_3= reshape(x5_3,[],1);
        x5_3= de2bi(x5_3);
        x5_3= reshape(x5_3.',1,[]);                                        %叠加信号用户3完成
        x5_decode = vitdec(x5_3,trellis,tblen,'trunc','hard');             %用户3比特流信息
        x5_3_bit=reshape(x5_decode.',Nsp,[]);
         %% NOMA解调 用户4
        x5_4=convenc(x5_decode ,trellis);
        x5_4=reshape(x5_4,2,[])';                                            % 将加密后的数据重构为一个[]*2矩阵，m=4
    	x5_4=bi2de(x5_4);                                                    % 二进制转化为十进制，方便下一步qpsk调制
        x5_4=reshape(x5_4,1,[]);
        x5_4_temp=pskmod(x5_4,M1,pi/M1);                                  
       x5_4_temp=reshape(x5_4_temp.',Nsp,[]);
        %叠加信号相减
        x5_4_temp2=x5-x5_4_temp;
        x5_4_temp2=reshape(x5_4_temp2,1,[]);
        x5_4_temp3=pskdemod(x5_4_temp2,M1,pi/M1);
        x5_4_temp4= reshape(x5_4_temp3,[],1);
        x5_4out= de2bi(x5_4_temp4);
        x5_4out2= reshape(x5_4out.',1,[]);%叠加信号用户3完成
        x5_4decode = vitdec(x5_4out2,trellis,tblen,'trunc','hard');    
        x5_4_bit=reshape(x5_4decode.',Nsp,[]);
        
        
        
        
        %% OMA 用户3解调
        x3 = reshape(x3,[],1);
        x3 = de2bi(x3);
        x3 = reshape(x3',1,[]);
        % 维特比译码
        trellis = poly2trellis(7,[133 171]);
        x3 = vitdec(x3,trellis,tblen,'trunc','hard');                  %硬判决
        x3=reshape(x3,52,[]);
        %% OMA 用户4解调
        
        x4= reshape(x4,[],1);
        x4 = de2bi(x4);
        x4= reshape(x4',1,[]);
        % 维特比译码
        trellis = poly2trellis(7,[133 171]);
        x4 = vitdec(x4,trellis,tblen,'trunc','hard');                  %硬判决
        x4=reshape(x4,52,[]);
       %% 统计一帧中的错误比特数
        %理论情况下
        [neb1(indx),temp]=biterr(x3, data_tx1yanzhen(:,(indx-1)*Nd+1:indx*Nd),log2(M1));             %用户1的两径信道误比特数
        [neb2(indx),temp]=biterr(x4,data_tx2yanzhen(:,(indx-1)*Nd+1:indx*Nd),log2(M1));             %用户2的两径信道误比特数
        [neb3(indx),temp]=biterr(x5_3_bit,data_tx1yanzhen(:,(indx-1)*Nd+1:indx*Nd),log2(M1));    %用户3的两径信道误比特数                                  
        [neb4(indx),temp]=biterr(x5_4_bit,data_tx2yanzhen(:,(indx-1)*Nd+1:indx*Nd),log2(M1));            %用户2的单信道误比特数
        
    end
    %% 误比特率
    ber1(ii)=sum(neb1)/(Nsp*log2(M1)*Nd*Nfrm);                             %用户1误码率     
    ber2(ii)=sum(neb2)/(Nsp*log2(M1)*Nd*Nfrm);                             %用户2误码率                                  
    ber3(ii)=sum(neb3)/(Nsp*log2(M1)*Nd*Nfrm);                             %用户3误码率    
    ber4(ii)=sum(neb4)/(Nsp*log2(M1)*Nd*Nfrm);                             %用户3误码率  
                               
            
end
%% 星座图显示


%figure(1);stem(msg);
%msg1 = pskmod(msg,4,pi/4); % 4psk调制  初始相位为 pi/4
%scatterplot(msg1); axis([-1.2,1.2,-1.2,1.2]);% 画星座图
%hold on;

%% 作图
figure(2)
semilogy(EbNo,ber1,'-ro',EbNo,ber2,'-rv',EbNo,ber3,'-r*',EbNo,ber4,'-rd')                   
grid on
title('NOMA与OFDM系统误比特率对比')
legend('OFDM用户3 ','OFDM用户4   ','NOMA用户3','NOMA用户4')
xlabel('Eb/No[dB]')
ylabel('BER')

    
    