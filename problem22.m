% �ڶ��ʵڣ�2��С��
% �Բ�ͬ�������¶ȼ�����ů���ɹ����ϵ����µ��Ŀɳ���ʱ�䣬��������ͬ�����¶ȶԹ����ϵ����µ����Ե�Ӱ�졣
clc;clear;
% ��������
C_in = 1.1e6; % ���ڿ�����Ч���� J/��
C_wall = 1.86e8; % ǽ���Ч���� J/��
R_1 = 1.2e-3; % ���ڿ�����ǽ���ڲ��Ч���� ��/W
R_2 = 9.2e-3; % ǽ���������������Ч���� ��/W
P_N = 8e3; % ���ů�豸����� W
T_in0 = 20; % ���ڳ�ʼ�¶� ��
T_wall0 = [17.32 16.865 16.38 15.88 15.365 14.85]; % ǽ���ʼ��̬�¶� ��

% ʱ������
dt =1; % ʱ�䲽�� 1s
T = 24*3600; % ��ʱ�� 24h
N = T/dt; % ʱ�䲽�� 
time = 0:dt:T-dt; % ʱ������ h

T_out = [0 -5 -10 -15 -20 -25]; % �����¶� ��
T_in = zeros(length(T_out),N); % �����¶� ��
T_wall = zeros(length(T_out),N); % ǽ���¶� ��
S = zeros(length(T_out),N); % ����״̬ 
P_heat = zeros(length(T_out),N); % ���ȹ��� W

up_time = zeros(length(T_out),N); % �����ϵ��ɳ���ʱ�� s
down_time = zeros(length(T_out),N); % �����µ��ɳ���ʱ�� s

% ѭ�����㲻ͬ�����¶��µĽ��
for i = 1:length(T_out)
    T_in(i,1) = T_in0; % ��ʼ����
    T_wall(i,1) = T_wall0(1,i); % ǽ���ʼ��̬�¶� ��
    S(i,1) = 1; % ��ʼ״̬Ϊ����
    P_heat(i,1) = S(i,1)*P_N; % ��ʼ���ȹ���
    
    for j = 2:N
        
        T_wall(i,j) = T_wall(i,j-1) + dt*((T_in(i,j-1)-T_wall(i,j-1))/(C_wall*R_1)-(T_wall(i,j-1)-T_out(i))/(C_wall*R_2)); % ŷ�������ǽ���¶�
        
        if 18 < T_in(i,j-1) < 22 % �¿��߼��жϿ���״̬
            S(i,j) = S(i,j-1);
        end
        
        if T_in(i,j-1) >= 22 
            S(i,j) = 0;
        end
        
        if T_in(i,j-1) <= 18 
            S(i,j) = 1;
        end
        
        P_heat(i,j) = S(i,j)*P_N; % �������ȹ���
        
        T_in(i,j) = T_in(i,j-1) + dt*(P_heat(i,j)/C_in-(T_in(i,j-1)-T_wall(i,j-1))/(C_in*R_1)); % ŷ������������¶�
        
        if S(i,j) == 0 % �ر�״̬�¿������ϵ���
            
            T_in_up = T_in(i,j); % ���ϵ���ʱ�������¶� ��
            T_wall_up = T_wall(i,j); % ���ϵ���ʱ��ǽ���¶� ��
            P_heat_up = P_N; % ���ϵ���ʱ�����ȹ��� W
            
            while T_in_up < 22
                
                up_time(i,j) = up_time(i,j) + dt; % �ۼ����ϵ��ڿɳ���ʱ��
                
                T_wall_up = T_wall_up + dt*((T_in_up-T_wall_up)/(C_wall*R_1)-(T_wall_up-T_out(i))/(C_wall*R_2)); % ŷ�������ǽ���¶�
                
                T_in_up = T_in_up + dt*(P_heat_up/C_in-(T_in_up-T_wall_up)/(C_in*R_1)); % ŷ������������¶�
                
            end
            
        end
        
        if S(i,j) == 1 % ����״̬�¿������µ���
            
            T_in_down = T_in(i,j); % ���µ���ʱ�������¶� ��
            T_wall_down = T_wall(i,j); % ���µ���ʱ��ǽ���¶� ��
            P_heat_down = 0; % ���µ���ʱ�����ȹ��� W
            
            while T_in_down > 18
                
                down_time(i,j) = down_time(i,j) + dt; %�ۼ����µ��ڿɳ���ʱ��
                
                T_wall_down = T_wall_down + dt*((T_in_down-T_wall_down)/(C_wall*R_1)-(T_wall_down-T_out(i))/(C_wall*R_2)); %ŷ�������ǽ���¶�
                
                T_in_down = T_in_down + dt*(P_heat_down/C_in-(T_in_down-T_wall_down)/(C_in*R_1)); %ŷ������������¶�
                
            end
            
        end
        
    end
    
end

% ���ƽ��
figure(1)
for i=1:length(T_out)
    subplot(2,3,i)
    plot(time/60,up_time(i,:)/60)
    hold on
    plot(time/60,down_time(i,:)/60)
    hold off
    xlabel('ʱ��/min')
    ylabel('�ɳ���ʱ��/min')
    title(['�����¶�Ϊ',num2str(T_out(i)),'��ʱ�Ĺ����ϵ����µ��ɳ���ʱ��'])
    legend('���ϵ���','���µ���')
end