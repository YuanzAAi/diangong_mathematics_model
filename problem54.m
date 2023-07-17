% ������
% �Ե��ůסլ�� 600 ��ס��Ϊ�����������ס����ʼ�����¶����¿������ھ��ȷֲ����ڱ� 1 ��ʾ�ĸ�����ƽ���¶��£�����ѡ��һ����ů�豸���صĳ�ʼ״̬
clc;clear;
% ��������
C_in = 1.1e6; % ���ڿ�����Ч���� J/��
C_wall = 1.86e8; % ǽ���Ч���� J/��
R_1 = 1.2e-3; % ���ڿ�����ǽ���ڲ��Ч���� ��/W
R_2 = 9.2e-3; % ǽ���������������Ч���� ��/W
P_N = 8e3; % ���ů�豸����� W

M = 600; % ס����Ŀ
T_out = [0 -5 -10 -15 -20 ]; % �����¶� ��

% ʱ������
days = [30 40 40 40 30]; % ����������Ӧ0��-5��-10��-15��-20��
dt =1; % ʱ�䲽�� 1h
if length(T_out) == 1 || length(T_out) == 5
    T = 24*30; 
elseif length(T_out) == 2 || length(T_out) == 3 || length(T_out) == 4
    T = 24*40;
end
N = T/dt; % ʱ�䲽�� 
time = 0:dt:T-dt; % ʱ������ h

T_in = zeros(M,N,length(T_out)); % �����¶� ��
T_wall = zeros(M,N,length(T_out)); % ǽ���¶� ��
S = zeros(M,N,length(T_out)); % ����״̬ 
P_heat = zeros(M,N,length(T_out)); % ���ȹ��� W

up_index = zeros(M,N,length(T_out)); % �ɲ����ϵ��ĵ��ů�豸���
down_index = zeros(M,N,length(T_out)); % �ɲ����µ��ĵ��ů�豸���

up_power = zeros(1,N,length(T_out)); % �ܿ��ϵ����� W
down_power = zeros(1,N,length(T_out)); % �ܿ��µ����� W

T_in0 = unifrnd(18,22,1,M); % ���ڳ�ʼ�¶� ��,���������18��22֮��ľ��ȷֲ�����ֵ
T_wall0 = [17.32 16.865 16.38 15.88 15.365 14.85]; % ǽ���ʼ��̬�¶� �� 
S0 = randi([0 1],1,M); % ���س�ʼ״̬ �������0��1

% ѭ�����㲻ͬס���µĽ��
for k = 1:length(T_out) % ������ͬ�������¶�
    
    for i = 1:M
        
        T_in(i,1,k) = T_in0(i); % ��ʼ����
        T_wall(i,1,k) = T_wall0(k); % ǽ���ʼ��̬�¶� ��
        S(i,1,k) = S0(i); % ��ʼ״̬Ϊ������ر�
        P_heat(i,1,k) = S(i,1,k)*P_N; % ��ʼ���ȹ���
        
        for j = 2:N
            
            T_wall(i,j,k) = T_wall(i,j-1,k) + dt*((T_in(i,j-1,k)-T_wall(i,j-1,k))/(C_wall*R_1)-(T_wall(i,j-1,k)-T_out(k))/(C_wall*R_2)); % ŷ�������ǽ���¶�
            
            if 18 < T_in(i,j-1,k) < 22 % �¿��߼��жϿ���״̬
                S(i,j,k) = S(i,j-1,k);
            end
            
            if T_in(i,j-1,k) >= 22 
                S(i,j,k) = 0;
            end
            
            if T_in(i,j-1,k) <= 18 
                S(i,j,k) = 1;
            end
            
            P_heat(i,j,k) = S(i,j,k)*P_N; %�������ȹ���
            
            T_in(i,j,k) = T_in(i,j-1,k) + dt*(P_heat(i,j,k)/C_in-(T_in(i,j-1,k)-T_wall(i,j-1,k))/(C_in*R_1)); % ŷ������������¶�
            
            if S(i,j,k) == 0 % �ر�״̬�¿������ϵ���
                
                T_in_up = T_in(i,j,k); % ���ϵ���ʱ�������¶� ��
                T_wall_up = T_wall(i,j,k); % ���ϵ���ʱ��ǽ���¶� ��
                P_heat_up = P_N; % ���ϵ���ʱ�����ȹ��� W
                up_index(i,j,k) = i;% ����������ϵ������¼��ס�������
                
                while T_in_up < 22
                    
                    T_wall_up = T_wall_up + dt*((T_in_up-T_wall_up)/(C_wall*R_1)-(T_wall_up-T_out(k))/(C_wall*R_2));%ŷ�������ǽ���¶�
                    
                    T_in_up = T_in_up + dt*(P_heat_up/C_in-(T_in_up-T_wall_up)/(C_in*R_1)); % ŷ������������¶�
                end
                
                
            end
            
            if S(i,j,k) == 1 % ����״̬�¿������µ���
                
                T_in_down = T_in(i,j,k); % ���µ���ʱ�������¶� ��
                T_wall_down = T_wall(i,j,k); % ���µ���ʱ��ǽ���¶� ��
                P_heat_down = 0; % ���µ���ʱ�����ȹ��� W
                down_index(i,j,k) = i; % ����������µ��ڣ����¼��ס�������
                
                while T_in_down > 18
                    
                    T_wall_down = T_wall_down + dt*((T_in_down-T_wall_down)/(C_wall*R_1)-(T_wall_down-T_out(k))/(C_wall*R_2)); %ŷ�������ǽ���¶�
                    
                    T_in_down = T_in_down + dt*(P_heat_down/C_in-(T_in_down-T_wall_down)/(C_in*R_1));%ŷ������������¶�
                    
                end
                
                
            end
            
        end
        
    end
    
    up_power(:,:,k) = sum(up_index(:,:,k)~=0,1)*P_N; % �����ܿ��ϵ����� w
    down_power(:,:,k) = sum(down_index(:,:,k)~=0,1)*P_N; % �����ܿ��µ����� w
    
end

figure(1)
for k=1:length(T_out)
    subplot(3,2,k)
    plot(time,up_power(:,:,k)/1000)
    xlabel('ʱ��/h')
    ylabel('�ܿ��ϵ�����/kW')
    title(['�����¶�Ϊ',num2str(T_out(k)),'��ʱ',',����Ϊ',num2str(days(1,k))])
end

figure(2)
for k=1:length(T_out)
    subplot(3,2,k)
    plot(time,down_power(:,:,k)/1000)
    xlabel('ʱ��/h')
    ylabel('�ܿ��µ�����/kW'),
    title(['�����¶�Ϊ',num2str(T_out(k)),'��ʱ',',����Ϊ',num2str(days(1,k))])
end
