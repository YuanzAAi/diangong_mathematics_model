dbstop if error %MATLAB����������
clc;clear;
% ��������
C_in = 1.1e6; % ���ڿ�����Ч���� J/��
C_wall = 1.86e8; % ǽ���Ч���� J/��
R_1 = 1.2e-3; % ���ڿ�����ǽ���ڲ��Ч���� ��/W
R_2 = 9.2e-3; % ǽ���������������Ч���� ��/W
P_N = 8e3; % ���ů�豸����� W
T_in0 = 20; % ���ڳ�ʼ�¶� ��
T_wall0 = [17.32 16.865 16.38 15.88 15.365 14.85]; % ǽ���ʼ��̬�¶� ��
T_out = [0 -5 -10 -15 -20 -25]; % �����¶� ��
price_peak = 0.56; % ��ʱ��� Ԫ/kWh
price_valley = 0.32; % ��ʱ��� Ԫ/kWh

% ʱ������
dt =1; % ʱ�䲽�� 1s
T = 24*3600; % ��ʱ�� 24h
N = T/dt; % ʱ�䲽�� 
time = 0:dt:T-dt; % ʱ������ h

% ��ʼ������
T_in = zeros(length(T_out),N); % �����¶� ��
T_wall = zeros(length(T_out),N); % ǽ���¶� ��
S = zeros(length(T_out),N); % ����״̬ 
P_heat = zeros(length(T_out),N); % ���ȹ��� W
E_day = zeros(length(T_out),1); % ���õ��� kWh
P_avg = zeros(length(T_out),1); % ��ƽ���õ繦�� kW
cost_day = zeros(length(T_out),1); % ���õ�ɱ� Ԫ

% ѭ�����㲻ͬ�����¶��µĽ��
for i = 1:length(T_out)
    T_in(i,1) = T_in0; % ��ʼ����
    T_wall(i,1) = T_wall0(1,i);
    S(i,1) = 1; % ��ʼ״̬Ϊ����
    P_heat(i,1) = S(i,1)*P_N; % ��ʼ���ȹ���
    
    for j = 2:N
        
        T_wall(i,j) = T_wall(i,j-1) + dt*((T_in(i,j-1)-T_wall(i,j-1))/(C_wall*R_1)-(T_wall(i,j-1)-T_out(i))/(C_wall*R_2)); % ŷ�������ǽ���¶�
        
        if 18 < T_in(i,j-1) < 22
            S(i,j) = S(i,j-1);
        end
        
        if T_in(i,j-1) >= 22 % �¿��߼��жϿ���״̬
            S(i,j) = 0;
        end
        
        
        if T_in(i,j-1) <= 18 
            S(i,j) = 1;
        end
        
        P_heat(i,j) = S(i,j)*P_N; % �������ȹ���
        
        T_in(i,j) = T_in(i,j-1) + dt*(P_heat(i,j)/C_in-(T_in(i,j-1)-T_wall(i,j-1))/(C_in*R_1)); % ŷ������������¶�
        
        if j <= 8*3600/dt || j > 21*3600/dt % �жϷ��ʱ��
            price = price_valley/3600;
        else
            price = price_peak/3600;
        end
        
        E_day(i) = E_day(i) + P_heat(i,j)*dt/(1000*3600); % �ۼ����õ���
        cost_day(i) = cost_day(i) + P_heat(i,j)*dt/1000*price; % �ۼ����õ�ɱ�
        
    end
    
    P_avg(i) = E_day(i)/(T/3600); % ������ƽ���õ繦��
    
end

% ���ƽ��
figure(1)
for i = 1:length(T_out)
    subplot(2,3,i)
    plot(time,T_in(i,:))
    xlabel('ʱ��/s')
    ylabel('�����¶�/��')
    title(['�����¶�Ϊ',num2str(T_out(i)),'��ʱ�������¶�'])
    ylim([18 22]) % ����y�᷶ΧΪ18��22
    set(gca,'YTick',[18 22]) % ����y��̶�Ϊ18��22
end

figure(2)
for i = 1:length(T_out)
    subplot(2,3,i)
    plot(time,T_wall(i,:))
    xlabel('ʱ��/s')
    ylabel('ǽ���¶�/��')
    title(['�����¶�Ϊ',num2str(T_out(i)),'��ʱ��ǽ���¶�'])
end

figure(3)
for i = 1:length(T_out)
    subplot(2,3,i)
    plot(time,S(i,:))
    xlabel('ʱ��/s')
    ylabel('����״̬')
    title(['�����¶�Ϊ',num2str(T_out(i)),'��ʱ�ĵ��ů�豸����״̬'])
    ylim([0 1]) % ����y�᷶ΧΪ0��1
    set(gca,'YTick',[0 1]) % ����y��̶�Ϊ0��1
end

figure(4)
for i = 1:length(T_out)
    subplot(2,3,i)
    plot(time,P_heat(i,:)/1000)
    xlabel('ʱ��/s')
    ylabel('���ȹ���/kW')
    title(['�����¶�Ϊ',num2str(T_out(i)),'��ʱ�ĵ��ů�豸���ȹ���'])
    ylim([0 8]) % ����y�᷶ΧΪ0��8
    set(gca,'YTick',[0 8]) % ����y��̶�Ϊ0��8
end

% ������
disp('�� 1 ����ס�����ů�����õ���Ϊ������ͳ�ƽ�������ڳ�ʼ�¶�Ϊ 20�棩')
disp(['�����¶�/��',' ','ƽ������ʱ��/min',' ','ƽ������ʱ��/min',' ','����/min',' ','ƽ��ռ�ձ�/%',' ','���õ���/kWh',' ','��ƽ���õ繦��/kW',' ','���õ�ɱ�/Ԫ'])
for i = 1:length(T_out)
    
    % ����ƽ������ʱ����ƽ������ʱ�������ں�ƽ��ռ�ձ�
    rise_time = 0; % ��������ʱ��
    fall_time = 0; % ���ν���ʱ��
    rise_count = 0; % ���´���
    fall_count = 0; % ���´���
    duty_cycle = 0; % ռ�ձ�
    
    for j = 2:N
        
        if S(i,j-1) == 0 && S(i,j) == 1 % ��ʼ����
            rise_count = rise_count + 1;
            rise_time = rise_time + dt;
            duty_cycle = duty_cycle + dt;
        end
        
        if S(i,j-1) == 1 && S(i,j) == 0 % ��ʼ����
            fall_count = fall_count + 1;
            fall_time = fall_time + dt;
        end
        
        if S(i,j-1) == 0 && S(i,j) == 0 % ��������
            fall_time = fall_time + dt;
        end
        
        if S(i,j-1) == 1 && S(i,j) == 1 % ��������
            rise_time = rise_time + dt;
            duty_cycle = duty_cycle + dt;
        end
        
    end
    
    avg_rise_time = (rise_time/rise_count)/60; % ƽ������ʱ�� min
    avg_fall_time = (fall_time/fall_count)/60; % ƽ������ʱ�� min
    period = (rise_time+fall_time)/(rise_count+fall_count)/60; % ���� min
    avg_duty_cycle = duty_cycle/T*100; % ƽ��ռ�ձ� 
    
    disp([num2str(T_out(i)),' ',num2str(avg_rise_time),' ',num2str(avg_fall_time),' ',num2str(period),' ',num2str(avg_duty_cycle),' ',num2str(E_day(i)),' ',num2str(P_avg(i)),' ',num2str(cost_day(i))])
    
end
% ��һ�ʵڣ�3��С��
% ���蹩ů��Ϊ 180 �죬

T_out_avg = [0 -5 -10 -15 -20]; % ����ƽ���¶� ��
days = [30 40 40 40 30]; % ��������
E_season = zeros(length(T_out_avg),1); % ��ů���õ��� kWh
cost_season = zeros(length(T_out_avg),1); % ��ů���õ�ɱ� Ԫ

for i = 1:length(T_out_avg)
    E_season(i) = E_day(i)*days(i); % ���㹩ů���õ���
    cost_season(i) = cost_day(i)*days(i); % ���㹩ů���õ�ɱ�
end

E_total = sum(E_season); % ���㹩ů�����õ���
cost_total = sum(cost_season); % ���㹩ů���ܳɱ�

% ������
disp('�� 2 ��ů�ڵ���ס���õ������õ�ɱ�ͳ�ƽ��')
disp(['����ƽ���¶�/��',' ','��������',' ','�õ���/kWh',' ','��ů�ɱ�/Ԫ'])
for i = 1:length(T_out_avg)
    disp([num2str(T_out_avg(i)),' ',num2str(days(i)),' ',num2str(E_season(i)),' ',num2str(cost_season(i))])
end
disp(['��ů�����õ���/kWh',' ','��ů���ܳɱ�/Ԫ'])
disp([num2str(E_total),' ',num2str(cost_total)])
