% ������
% �Ե��ůסլ�� 600 ��ס��Ϊ�����������ס����ʼ�����¶����¿������ھ��ȷֲ����ڱ� 1 ��ʾ�ĸ�����ƽ���¶��£�����ѡ��һ����ů�豸���صĳ�ʼ״̬
clc;clear;
% ��������
C_in = 1.1e6; % ���ڿ�����Ч���� J/��
C_wall = 1.86e8; % ǽ���Ч���� J/��
R_1 = 1.2e-3; % ���ڿ�����ǽ���ڲ��Ч���� ��/W
R_2 = 9.2e-3; % ǽ���������������Ч���� ��/W
P_N = 8e3; % ���ů�豸����� W

% ʱ������
dt =1; % ʱ�䲽�� 1min
T = 24*60; % ��ʱ�� 24h
N = T/dt; % ʱ�䲽�� 
time = 0:dt:T-dt; % ʱ������ h

M = 600; % ס����Ŀ
T_out = [0 -5 -10 -15 -20 -25]; % �����¶� ��
T_in = zeros(M,N,length(T_out)); % �����¶� ��
T_wall = zeros(M,N,length(T_out)); % ǽ���¶� ��
S = zeros(M,N,length(T_out)); % ����״̬ 
P_heat = zeros(M,N,length(T_out)); % ���ȹ��� W

up_time = zeros(M,N,length(T_out)); % �����ϵ��ɳ���ʱ�� min
down_time = zeros(M,N,length(T_out)); % �����µ��ɳ���ʱ�� min
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
                    
                    up_time(i,j,k) = up_time(i,j,k) + dt; %�ۼ����ϵ��ڿɳ���ʱ��
                    
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
                    
                    down_time(i,j,k) = down_time(i,j,k) + dt; %�ۼ����µ��ڿɳ���ʱ��
                    
                    T_wall_down = T_wall_down + dt*((T_in_down-T_wall_down)/(C_wall*R_1)-(T_wall_down-T_out(k))/(C_wall*R_2)); %ŷ�������ǽ���¶�
                    
                    T_in_down = T_in_down + dt*(P_heat_down/C_in-(T_in_down-T_wall_down)/(C_in*R_1));%ŷ������������¶�
                    
                end
                
                
            end
            
        end
        
    end
    
    up_power(:,:,k) = sum(up_index(:,:,k)~=0,1)*P_N; % �����ܿ��ϵ����� w
    down_power(:,:,k) = sum(down_index(:,:,k)~=0,1)*P_N; % �����ܿ��µ����� w
    
end
% ���ƽ��
figure(1)
for k=1:length(T_out)
    subplot(3,2,k)
    plot(time,up_index(:,:,k),'o')
    xlabel('ʱ��/min')
    ylabel('�ɲ����ϵ��ĵ��ů�豸���')
    ylim([1 600])  % ����y��̶ȷ�ΧΪ1��600
    title(['�����¶�Ϊ',num2str(T_out(k)),'��ʱ'])
end

figure(2)
for k=1:length(T_out)
    subplot(3,2,k)
    plot(time,down_index(:,:,k),'o')
    xlabel('ʱ��/min')
    ylabel('�ɲ����µ��ĵ��ů�豸���')
    ylim([1 600])  % ����y��̶ȷ�ΧΪ1��600
    title(['�����¶�Ϊ',num2str(T_out(k)),'��ʱ'])
end

figure(3)
for k=1:length(T_out)
    subplot(3,2,k)
    plot(time,up_power(:,:,k)/1000)
    xlabel('ʱ��/min')
    ylabel('�ܿ��ϵ�����/kW')
    title(['�����¶�Ϊ',num2str(T_out(k)),'��ʱ'])
end

figure(4)
for k=1:length(T_out)
    subplot(3,2,k)
    plot(time,down_power(:,:,k)/1000)
    xlabel('ʱ��/min')
    ylabel('�ܿ��µ�����/kW')
    title(['�����¶�Ϊ',num2str(T_out(k)),'��ʱ'])
end

figure(5) % 24Сʱ��600��ס�������õ繦������
for k=1:length(T_out)
    subplot(3,2,k)
    plot(time,sum(P_heat(:,:,k))/1000)
    xlabel('ʱ��/min')
    ylabel('���õ繦��/kW')
    title('24Сʱ��600��ס�������õ繦������')
end

%������
% ��ʱ��ۼ��������񲹳�����

peak_price=0.56;%��ʱ��� Ԫ/kWh

valley_price=0.32;%��ʱ��� Ԫ/kWh

peak_compensation=1.30;%���岹���۸� Ԫ/kWh

valley_compensation=0.65;%��Ȳ����۸� Ԫ/kWh

peak_start=16*60;%���忪ʼʱ�� min

peak_end=20*60;%�������ʱ�� min

valley_start=0*60;%��ȿ�ʼʱ�� min

valley_end=4*60;%��Ƚ���ʱ�� min

% ����ͳɱ���ʼ��

benefit_per=zeros(M,N,length(T_out));% ÿ��ÿʱ���� Ԫ

cost_per=zeros(M,N,length(T_out));% ÿ��ÿʱ�ɱ� Ԫ

benefit=zeros(1,length(T_out));%������ Ԫ

cost=zeros(1,length(T_out));%�ܳɱ� Ԫ

% ���������ں�������¶ȱ仯�����ʶ�ָ�꣬����������ͳɱ�

T_in_new = zeros(M,N,length(T_out)); % ������ں�������¶� ��

T_wall_new = zeros(M,N,length(T_out)); % ������ں��ǽ���¶� ��

S_new = zeros(M,N,length(T_out)); % ������ں�Ŀ���״̬

P_heat_new = zeros(M,N,length(T_out)); % ������ں�����ȹ��� W

deviation = zeros(1,length(T_out)); % �¶�ƫ��� ��

unstability = zeros(1,length(T_out)); % �¶��ȶ��� ��/min

change_num = zeros(M, N, length(T_out)); % ��ʱ�����ڲ���������ڵ��¿�����״̬�����仯�ĵ��ů�豸����

peak_power = zeros(1,length(T_out));

valley_power = zeros(1,length(T_out));

for k = 1:length(T_out)
    
    peak_power(1,k) = max(down_power(1,16*60+1:20*60,k)); % ����ʱ�ο��ṩ�ĳ���������µ��ڹ���ֵ W
    
    valley_power(1,k) = max(up_power(1,1:4*60,k)); % ���ʱ�ο��ṩ�ĳ���������ϵ��ڹ���ֵ W


for i = 1:M
    
    T_in_new(i,1,k) = T_in0(i); % ��ʼ����
    
    T_wall_new(i,1,k) = T_wall0(k); % ǽ���ʼ��̬�¶� ��
    
    S_new(i,1,k) = S0(i); % ��ʼ״̬Ϊ������ر�
    
    P_heat_new(i,1,k) = S_new(i,1,k)*P_N; % ��ʼ���ȹ���
    
    for j = 2:N
        
        T_wall_new(i,j,k) = T_wall_new(i,j-1,k) + dt*((T_in_new(i,j-1,k)-T_wall_new(i,j-1,k))/(C_wall*R_1)-(T_wall_new(i,j-1,k)-T_out(k))/(C_wall*R_2));%ŷ�������ǽ���¶�
        
        if (peak_start <= j && j <= peak_end && down_index(i,j,k) == i) || T_in_new(i,j-1,k) >= 22 % ����ʱ���ҿ������µ��ڻ��������¶ȴ��ڵ����¿���������
            
            S_new(i,j,k) = 0; % �رյ��ů�豸
            
            
        elseif (valley_start <= j && j <= valley_end && up_index(i,j,k) == i) || T_in_new(i,j-1,k) <= 18 % ���ʱ���ҿ������ϵ��ڻ��������¶�С�ڵ����¿���������
            
            S_new(i,j,k) = 1; % �������ů�豸
            
        else
            
            S_new(i,j,k) = S(i,j,k); % �����������ԭ״̬
            
        end
        
        if peak_start <= j && j <= peak_end && S_new(i,j-1,k) ~= S_new(i,j,k) % ������ʱ�����������״̬�����仯
            
            change_num(i,j,k) = 1; % ����Ӧ��λ�ü�һ
        
        elseif valley_start <= j && j <= valley_end && S_new(i,j-1,k) ~= S_new(i,j,k)% �����ʱ�����������״̬�����仯
            
            change_num(i,j,k) = 1; % ����Ӧ��λ�ü�һ
            
        end
        
        P_heat_new(i,j,k) = S_new(i,j,k)*P_N; %�������ȹ���
        
        T_in_new(i,j,k) = T_in_new(i,j-1,k) + dt*(P_heat_new(i,j,k)/C_in-(T_in_new(i,j-1,k)-T_wall_new(i,j-1,k))/(C_in*R_1));%ŷ������������¶�
        
        if   4*60 < j && j < 8*60  %��ƽʱ���У�04:00����08:00��21:00����24:00Ϊ��ʱ��08:00����16:00��20:00����21:00Ϊ��ʱ��
        
                 cost_per(i,j,k) = cost_per(i,j,k) + P_heat_new(i,j,k)/1000*valley_price*dt/60; % ƽʱ�ɱ�
                 
        elseif   8*60 <= j && j < 16*60
        
                 cost_per(i,j,k) = cost_per(i,j,k) + P_heat_new(i,j,k)/1000*peak_price*dt/60; % ƽʱ�ɱ�
        
        elseif   20*60 < j && j <= 21*60
        
                 cost_per(i,j,k) = cost_per(i,j,k) + P_heat_new(i,j,k)/1000*peak_price*dt/60; % ƽʱ�ɱ�
        
        elseif   21*60 < j && j <= 24*60
        
                 cost_per(i,j,k) = cost_per(i,j,k) + P_heat_new(i,j,k)/1000*valley_price*dt/60; % ƽʱ�ɱ�
        end
            
        deviation(k) = deviation(k) + abs(T_in_new(i,j,k)-T_in0(1,i));%�ۼ��¶�ƫ���,���ʼ�����¶ȵĲ�ľ���ֵ
        
        unstability(k) = unstability(k) + abs(T_in_new(i,j,k)-T_in_new(i,j-1,k))/dt;%�ۼ��¶Ȳ��ȶ���
        
    end
    
end
end

%����������ȵĳɱ�������
for k = 1:length(T_out)
    for i = 1:M
        in_peak = false; % ����Ƿ�������ʱ��
        in_valley = false; % ����Ƿ������ʱ��
        j = 2; % ʱ������
        
        while peak_start <= j && j <= peak_end
            cost_per(i,j,k) = cost_per(i,j,k) + P_heat_new(i,j,k)/1000*valley_price*dt/60; %�����ڵĹ��ȳɱ������ܳɱ���
            if S_new(i, j-1, k) == 1 && S_new(i, j, k) == 0 && ~in_peak % ����Ƿ����㿪ʼ���������
                in_peak = true; % ��������ʱ��
            end
            
            if in_peak % �����������ʱ��
                benefit_per(i, j, k) = benefit_per(i, j, k) + P_N/1000 * peak_compensation * dt / 60; % ��������
                
                if S_new(i, j, k) == 1 % ����Ƿ�����������������
                    in_peak = false; % ��������ʱ��
                end
            end
            
            j = j + 1; % �ƶ�����һ������
        end
        
        while  valley_start <= j && j <= valley_end
            if S_new(i, j-1, k) == 0 && S_new(i, j, k) == 1 && ~in_valley % ����Ƿ����㿪ʼ��ȵ�����
                in_valley = true;
            end
             if in_valley % ����������ʱ��
                benefit_per(i, j, k) = benefit_per(i, j, k) + P_N/1000 * valley_compensation * dt / 60; % �������
                cost_per(i,j,k) = cost_per(i,j,k) + P_N/1000*valley_price*dt/60; % ��ȳɱ�(���Ҫ����)
                
                if S_new(i, j, k) == 0 % ����Ƿ����������ȵ�����
                    in_peak = false; % ��������ʱ��
                end
             end
             j = j + 1; % �ƶ�����һ������
        end
    end  
end
%������������ܳɱ�

for k = 1:length(T_out)
    
    benefit = sum(sum(benefit_per, 1), 2); % ���;
    benefit = squeeze(benefit); % ȥ����һά��
    
    cost = sum(sum(cost_per, 1), 2); % ���;
    cost = squeeze(cost);
    
    change_permin = sum(change_num, 1); % �Ե�һά���
end

%�ı������״��Ϊ��ͼ׼��
benefit = benefit';
cost = cost';

% ���ƽ��

figure(6)

% ������ں�������¶�����

for k=1:length(T_out)

subplot(3,2,k)

plot(time,T_in(:,:,k),'b',time,T_in_new(:,:,k),'r')

xlabel('ʱ��/min')

ylabel('�����¶�/��')

legend('���������','�������')

title(['�����¶�Ϊ',num2str(T_out(k)),'��ʱ'])
end



% ����ͳɱ���״ͼ
figure(7)
bar([benefit' cost'])
xlabel('�����¶�/��')
ylabel('���/Ԫ')
legend('����','�ɱ�')
set(gca,'xticklabel',num2str(T_out'))
title('����������ȵ�����ͳɱ�')


%����������ȶ����ʶȵ�Ӱ���ͼ
figure(8)

subplot(2,1,1)

bar(deviation', 'b')  % ʹ����ɫ��blue��������״ͼ

xlabel('�����¶�/��')
ylabel('ָ��ֵ')
legend('�¶�ƫ���')
set(gca,'xticklabel',num2str(T_out'))
title('����������ȶ����ʶȵ�Ӱ��')

subplot(2,1,2)

bar(unstability', 'g')  % ʹ����ɫ��green��������״ͼ

xlabel('�����¶�/��')
ylabel('ָ��ֵ')
legend('�¶Ȳ��ȶ���')
set(gca,'xticklabel',num2str(T_out'))
title('����������ȶ����ʶȵ�Ӱ��')

% �����ʵڣ�3��С��
% ���蹩ů��Ϊ 180 ��
%Ҫ����������¶��¸�סլ�����ů���ɲ���������ȵ��ϵ����µ����ʣ���Ҫ����ʱ�����

days = [30 40 40 40 30]; % ����������Ӧ0��-5��-10��-15��-20��
cost_per_season = zeros(1,5); %���ȳɱ�
benefit_per_season = zeros(1,5); %��������
for k = 1:length(T_out)-1
    cost_per_season(1,k) = days(k)*cost(1,k);  %���ȳɱ�����
    benefit_per_season(1,k) = days(k)*benefit(1,k); %�����������
end
cost_year = sum(cost_per_season(1,:))/(180/365); %���ܳɱ�����
benefit_year = sum(benefit_per_season(1,:)); %�����������
benefit_year_person = benefit_year / M; %ƽ��ÿ����������
cost_savingratio = (benefit_year/cost_year)*100;%��ʡ�Ĺ��ȳɱ��ٷֱ�
