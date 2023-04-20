function [h]=rayleigh(fd,t)
%�ó������øĽ���jakesģ��������������ƽ̹������˥���ŵ�
%Yahong R.Zheng and Chengshan Xiao "Improved Models for 
%the Generation of Multiple Uncorrelated Rayleigh Fading Waveforms" 
%IEEE Commu letters, Vol.6, NO.6, JUNE 2002
%�������˵����
%  fd���ŵ�����������Ƶ�� ��λHz     
%  t :�źŵĳ���ʱ�����У����������λs  
%  hΪ����������ŵ���������һ��ʱ�亯�������� 

    %��������䲨��Ŀ
    N=40; 

    wm=2*pi*fd;
    %ÿ���޵����䲨��Ŀ��������Ŀ
    N0=N/4;
    %�ŵ�������ʵ��
    Tc=zeros(1,length(t));
    %�ŵ��������鲿
    Ts=zeros(1,length(t));
    %��һ������ϵ��
    P_nor=sqrt(1/N0);
    %�������·���ľ��ȷֲ������λ
    theta=2*pi*rand(1,1)-pi;
    for ii=1:N0
          %��i�����䲨������� 
            alfa(ii)=(2*pi*ii-pi+theta)/N;
            %��ÿ�����ز�������(-pi,pi)֮����ȷֲ��������λ
            fi_tc=2*pi*rand(1,1)-pi;
            fi_ts=2*pi*rand(1,1)-pi;
            %����弤��Ӧ����
            Tc=Tc+cos(cos(alfa(ii))*wm*t+fi_tc);
            Ts=Ts+cos(sin(alfa(ii))*wm*t+fi_ts);
    end;
    %�˹�һ������ϵ���õ����亯��
   h=P_nor*(Tc+j*Ts );