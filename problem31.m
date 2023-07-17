% ������
% �� 6�����ůס������ŷֱ�Ϊ 1-6��Ϊ�������������¶�Ϊ-20�棬���ڳ�ʼ�¶����¿������ھ��ȷֲ�������ѡ��һ����ů�豸���صĳ�ʼ״̬
clc;clear;
% ��������
C_in = 1.1e6; % ���ڿ�����Ч���� J/��
C_wall = 1.86e8; % ǽ���Ч���� J/��
R_1 = 1.2e-3; % ���ڿ�����ǽ���ڲ��Ч���� ��/W
R_2 = 9.2e-3; % ǽ���������������Ч���� ��/W
P_N = 8e3; % ���ů�豸����� W

% ʱ������
dt =1; % ʱ�䲽�� 1s
T = 24*3600; % ��ʱ�� 24h
N = T/dt; % ʱ�䲽�� 
time = 0:dt:T-dt; % ʱ������ h

M = 6; % ס����Ŀ
T_out_20 = -20; % �����¶� ��
T_in_20 = zeros(M,N); % �����¶� ��
T_wall_20 = zeros(M,N); % ǽ���¶� ��
S_20 = zeros(M,N); % ����״̬ 
P_heat_20 = zeros(M,N); % ���ȹ��� W

up_time_20 = zeros(M,N); % �����ϵ��ɳ���ʱ�� s
down_time_20 = zeros(M,N); % �����µ��ɳ���ʱ�� s

T_in0_20 = linspace(18,22,M); % ���ڳ�ʼ�¶� ��
T_wall0_20 = [15.365 15.365 15.365 15.365 15.365 15.365]; % ǽ���ʼ��̬�¶� ��
S0_20 = randi([0 1],1,M); % ���س�ʼ״̬

% ѭ�����㲻ͬס���µĽ��
for i = 1:M
    
    T_in_20(i,1) = T_in0_20(i); % ��ʼ����
    T_wall_20(i,1) = T_wall0_20(i); % ǽ���ʼ��̬�¶� ��
    S_20(i,1) = S0_20(i); % ��ʼ״̬Ϊ������ر�
    P_heat_20(i,1) = S_20(i,1)*P_N; % ��ʼ���ȹ���
    
    for j = 2:N
        
        T_wall_20(i,j) = T_wall_20(i,j-1) + dt*((T_in_20(i,j-1)-T_wall_20(i,j-1))/(C_wall*R_1)-(T_wall_20(i,j-1)-T_out_20)/(C_wall*R_2)); % ŷ�������ǽ���¶�
        
        if 18 < T_in_20(i,j-1) < 22 % �¿��߼��жϿ���״̬
            S_20(i,j) = S_20(i,j-1);
        end
        
        if T_in_20(i,j-1) >= 22 
            S_20(i,j) = 0;
        end
        
        if T_in_20(i,j-1) <= 18 
            S_20(i,j) = 1;
        end
        
        P_heat_20(i,j) = S_20(i,j)*P_N; % �������ȹ���
        
        T_in_20(i,j) = T_in_20(i,j-1) + dt*(P_heat_20(i,j)/C_in-(T_in_20(i,j-1)-T_wall_20(i,j-1))/(C_in*R_1)); % ŷ������������¶�
        
        if S_20(i,j) == 0 % �ر�״̬�¿������ϵ���
            
            T_in_up = T_in_20(i,j); % ���ϵ���ʱ�������¶� ��
            T_wall_up = T_wall_20(i,j); % ���ϵ���ʱ��ǽ���¶� ��
            P_heat_up = P_N; % ���ϵ���ʱ�����ȹ��� W
            
            while T_in_up < 22
                
                up_time_20(i,j) = up_time_20(i,j) + dt; %�ۼ����ϵ��ڿɳ���ʱ��
                
                T_wall_up = T_wall_up + dt*((T_in_up-T_wall_up)/(C_wall*R_1)-(T_wall_up-T_out_20)/(C_wall*R_2));%ŷ�������ǽ���¶�
                
                T_in_up = T_in_up + dt*(P_heat_up/C_in-(T_in_up-T_wall_up)/(C_in*R_1)); % ŷ������������¶�
            end
        end
         if S_20(i,j) == 1 % ����״̬�¿������µ���
            
            T_in_down = T_in_20(i,j); % ���µ���ʱ�������¶� ��
            T_wall_down = T_wall_20(i,j); % ���µ���ʱ��ǽ���¶� ��
            P_heat_down = 0; % ���µ���ʱ�����ȹ��� W
            
            while T_in_down > 18
                
                down_time_20(i,j) = down_time_20(i,j) + dt; %�ۼ����µ��ڿɳ���ʱ��
                
                T_wall_down = T_wall_down + dt*((T_in_down-T_wall_down)/(C_wall*R_1)-(T_wall_down-T_out_20)/(C_wall*R_2)); %ŷ�������ǽ���¶�
                
                T_in_down = T_in_down + dt*(P_heat_down/C_in-(T_in_down-T_wall_down)/(C_in*R_1)); %ŷ������������¶�
                
            end
            
         end
        
    end
    
end

% ���ƽ��
figure(1)
for i=1:M
    subplot(2,3,i)
    plot(time/3600,T_in_20(i,:))
    xlabel('ʱ��/h')
    ylabel('�¶�/��')
    title(['ס��',num2str(i),'�������¶�'])
    legend('�����¶�')
end
figure(2)
for i=1:M
    subplot(2,3,i)
    plot(time/3600,T_wall_20(i,:))
    xlabel('ʱ��/h')
    ylabel('�¶�/��')
    title(['ס��',num2str(i),'��ǽ���¶�'])
    legend('ǽ���¶�')
end

figure(3)
for i=1:M
    subplot(3,2,i)
    plot(time/3600,S_20(i,:))
    xlabel('ʱ��/h')
    ylabel('����״̬')
    title(['ס��',num2str(i),'�ĵ��ů�豸����״̬'])
    set(gca,'YTick',[0 1]) % ����y��̶�Ϊ0��1
end

figure(4)
for i=1:M
    subplot(2,3,i)
    plot(time/3600,up_time_20(i,:)/3600)
    hold on
    plot(time/3600,down_time_20(i,:)/3600)
    hold off
    xlabel('ʱ��/h')
    ylabel('�ɳ���ʱ��/h')
    title(['�û�',num2str(i),'�Ĺ����ϵ����µ��ɳ���ʱ��'])
    legend('���ϵ���','���µ���')
end

figure(5) % 24Сʱ��6��ס�������õ繦������
plot(time/3600,sum(P_heat_20(:,:))/1000)
xlabel('ʱ��/h')
ylabel('���õ繦��/kW')
title('24Сʱ��6��ס�������õ繦������')
