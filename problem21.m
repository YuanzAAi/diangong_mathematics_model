% �ڶ��ʵڣ�1��С��
% �Ե���ס�����ů����Ϊ���������¶�Ϊ-15�棬���ڳ�ʼ�¶�Ϊ20��
% ���ů�豸���صĳ�ʼ״̬Ϊ�������������ס�����ů����������24h��ʱ�㣨��� 1min�������ϵ����µ��Ŀɳ���ʱ�䣬�����Ƽ�������
clc;clear;
% ��������
C_in = 1.1e6; % ���ڿ�����Ч���� J/��
C_wall = 1.86e8; % ǽ���Ч���� J/��
R_1 = 1.2e-3; % ���ڿ�����ǽ���ڲ��Ч���� ��/W
R_2 = 9.2e-3; % ǽ���������������Ч���� ��/W
P_N = 8e3; % ���ů�豸����� W
T_in0 = 20; % ���ڳ�ʼ�¶� ��

% ʱ������
dt =1; % ʱ�䲽�� 1min
T = 24*3600; % ��ʱ�� 24h
N = T/dt; % ʱ�䲽�� 
time = 0:dt:T-dt; % ʱ������ h

T_out_15 = -15; % �����¶� ��
T_in_15 = zeros(1,N); % �����¶� ��
T_wall_15 = zeros(1,N); % ǽ���¶� ��
S_15 = zeros(1,N); % ����״̬ 
P_heat_15 = zeros(1,N); % ���ȹ��� W

T_in_15(1) = T_in0; % ��ʼ����
T_wall_15(1) = 15.88; % ǽ���ʼ��̬�¶� ��
S_15(1) = 1; % ��ʼ״̬Ϊ����
P_heat_15(1) = S_15(1)*P_N; % ��ʼ���ȹ���

up_time = zeros(1,N); % �����ϵ��ɳ���ʱ�� s
down_time = zeros(1,N); % �����µ��ɳ���ʱ�� s

for j = 2:N
    
    T_wall_15(j) = T_wall_15(j-1) + dt*((T_in_15(j-1)-T_wall_15(j-1))/(C_wall*R_1)-(T_wall_15(j-1)-T_out_15)/(C_wall*R_2)); % ŷ�������ǽ���¶�
    
    if 18 < T_in_15(j-1) < 22 % �¿��߼��жϿ���״̬
        S_15(j) = S_15(j-1);
    end
    
    if T_in_15(j-1) >= 22 
        S_15(j) = 0;
    end
    
    if T_in_15(j-1) <= 18 
        S_15(j) = 1;
    end
    
    P_heat_15(j) = S_15(j)*P_N; % �������ȹ���
    
    T_in_15(j) = T_in_15(j-1) + dt*(P_heat_15(j)/C_in-(T_in_15(j-1)-T_wall_15(j-1))/(C_in*R_1)); % ŷ������������¶�
    
    if S_15(j) == 0 % �ر�״̬�¿��Բ������ϵ���
        
        T_in_up = T_in_15(j); % ���ϵ���ʱ�������¶� ��
        T_wall_up = T_wall_15(j); % ���ϵ���ʱ��ǽ���¶� ��
        P_heat_up = P_N; % ���ϵ���ʱ�����ȹ��� W
        
        while T_in_up < 22
            
            up_time(j) = up_time(j) + dt; % �ۼ����ϵ��ڿɳ���ʱ��
            
            T_wall_up = T_wall_up + dt*((T_in_up-T_wall_up)/(C_wall*R_1)-(T_wall_up-T_out_15)/(C_wall*R_2)); % ŷ�������ǽ���¶�
            
            T_in_up = T_in_up + dt*(P_heat_up/C_in-(T_in_up-T_wall_up)/(C_in*R_1)); % ŷ������������¶�
            
        end
        
    end
    
    if S_15(j) == 1 % ����״̬�¿��Բ������µ���
        
        T_in_down = T_in_15(j); % ���µ���ʱ�������¶� ��
        T_wall_down = T_wall_15(j); % ���µ���ʱ��ǽ���¶� ��
        P_heat_down = 0; % ���µ���ʱ�����ȹ��� W
        
        while T_in_down > 18
            
            down_time(j) = down_time(j) + dt; %�ۼ����µ��ڿɳ���ʱ��
            
            T_wall_down = T_wall_down + dt*((T_in_down-T_wall_down)/(C_wall*R_1)-(T_wall_down-T_out_15)/(C_wall*R_2)); %ŷ�������ǽ���¶�
            
            T_in_down = T_in_down + dt*(P_heat_down/C_in-(T_in_down-T_wall_down)/(C_in*R_1)); %ŷ������������¶�
            
        end
        
    end
    
end

% ���ƽ��
figure(1)
plot(time/60,up_time/60)
hold on
plot(time/60,down_time/60)
hold off
xlabel('ʱ��/min')
ylabel('�ɳ���ʱ��/min')
title('�����¶�Ϊ-15��ʱ�Ĺ����ϵ����µ��ɳ���ʱ��')
legend('���ϵ���','���µ���')