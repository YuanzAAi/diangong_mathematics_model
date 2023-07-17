dbstop if error %MATLAB调试神命令
clc;clear;
% 参数设置
C_in = 1.1e6; % 室内空气等效热容 J/℃
C_wall = 1.86e8; % 墙体等效热容 J/℃
R_1 = 1.2e-3; % 室内空气和墙体内侧等效热阻 ℃/W
R_2 = 9.2e-3; % 墙体外侧和室外空气等效热阻 ℃/W
P_N = 8e3; % 电采暖设备额定功率 W
T_in0 = 20; % 室内初始温度 ℃
T_wall0 = [17.32 16.865 16.38 15.88 15.365 14.85]; % 墙体初始稳态温度 ℃
T_out = [0 -5 -10 -15 -20 -25]; % 室外温度 ℃
price_peak = 0.56; % 峰时电价 元/kWh
price_valley = 0.32; % 谷时电价 元/kWh

% 时间设置
dt =1; % 时间步长 1s
T = 24*3600; % 总时间 24h
N = T/dt; % 时间步数 
time = 0:dt:T-dt; % 时间序列 h

% 初始化变量
T_in = zeros(length(T_out),N); % 室内温度 ℃
T_wall = zeros(length(T_out),N); % 墙体温度 ℃
S = zeros(length(T_out),N); % 开关状态 
P_heat = zeros(length(T_out),N); % 制热功率 W
E_day = zeros(length(T_out),1); % 日用电量 kWh
P_avg = zeros(length(T_out),1); % 日平均用电功率 kW
cost_day = zeros(length(T_out),1); % 日用电成本 元

% 循环计算不同室外温度下的结果
for i = 1:length(T_out)
    T_in(i,1) = T_in0; % 初始条件
    T_wall(i,1) = T_wall0(1,i);
    S(i,1) = 1; % 初始状态为开启
    P_heat(i,1) = S(i,1)*P_N; % 初始制热功率
    
    for j = 2:N
        
        T_wall(i,j) = T_wall(i,j-1) + dt*((T_in(i,j-1)-T_wall(i,j-1))/(C_wall*R_1)-(T_wall(i,j-1)-T_out(i))/(C_wall*R_2)); % 欧拉法求解墙体温度
        
        if 18 < T_in(i,j-1) < 22
            S(i,j) = S(i,j-1);
        end
        
        if T_in(i,j-1) >= 22 % 温控逻辑判断开关状态
            S(i,j) = 0;
        end
        
        
        if T_in(i,j-1) <= 18 
            S(i,j) = 1;
        end
        
        P_heat(i,j) = S(i,j)*P_N; % 计算制热功率
        
        T_in(i,j) = T_in(i,j-1) + dt*(P_heat(i,j)/C_in-(T_in(i,j-1)-T_wall(i,j-1))/(C_in*R_1)); % 欧拉法求解室内温度
        
        if j <= 8*3600/dt || j > 21*3600/dt % 判断峰谷时段
            price = price_valley/3600;
        else
            price = price_peak/3600;
        end
        
        E_day(i) = E_day(i) + P_heat(i,j)*dt/(1000*3600); % 累计日用电量
        cost_day(i) = cost_day(i) + P_heat(i,j)*dt/1000*price; % 累计日用电成本
        
    end
    
    P_avg(i) = E_day(i)/(T/3600); % 计算日平均用电功率
    
end

% 绘制结果
figure(1)
for i = 1:length(T_out)
    subplot(2,3,i)
    plot(time,T_in(i,:))
    xlabel('时间/s')
    ylabel('室内温度/℃')
    title(['室外温度为',num2str(T_out(i)),'℃时的室内温度'])
    ylim([18 22]) % 设置y轴范围为18到22
    set(gca,'YTick',[18 22]) % 设置y轴刻度为18和22
end

figure(2)
for i = 1:length(T_out)
    subplot(2,3,i)
    plot(time,T_wall(i,:))
    xlabel('时间/s')
    ylabel('墙体温度/℃')
    title(['室外温度为',num2str(T_out(i)),'℃时的墙体温度'])
end

figure(3)
for i = 1:length(T_out)
    subplot(2,3,i)
    plot(time,S(i,:))
    xlabel('时间/s')
    ylabel('开关状态')
    title(['室外温度为',num2str(T_out(i)),'℃时的电采暖设备开关状态'])
    ylim([0 1]) % 设置y轴范围为0到1
    set(gca,'YTick',[0 1]) % 设置y轴刻度为0和1
end

figure(4)
for i = 1:length(T_out)
    subplot(2,3,i)
    plot(time,P_heat(i,:)/1000)
    xlabel('时间/s')
    ylabel('制热功率/kW')
    title(['室外温度为',num2str(T_out(i)),'℃时的电采暖设备制热功率'])
    ylim([0 8]) % 设置y轴范围为0到8
    set(gca,'YTick',[0 8]) % 设置y轴刻度为0和8
end

% 输出结果
disp('表 1 典型住户电采暖负荷用电行为特征量统计结果（室内初始温度为 20℃）')
disp(['室外温度/℃',' ','平均升温时长/min',' ','平均降温时长/min',' ','周期/min',' ','平均占空比/%',' ','日用电量/kWh',' ','日平均用电功率/kW',' ','日用电成本/元'])
for i = 1:length(T_out)
    
    % 计算平均升温时长、平均降温时长、周期和平均占空比
    rise_time = 0; % 单次升温时长
    fall_time = 0; % 单次降温时长
    rise_count = 0; % 升温次数
    fall_count = 0; % 降温次数
    duty_cycle = 0; % 占空比
    
    for j = 2:N
        
        if S(i,j-1) == 0 && S(i,j) == 1 % 开始升温
            rise_count = rise_count + 1;
            rise_time = rise_time + dt;
            duty_cycle = duty_cycle + dt;
        end
        
        if S(i,j-1) == 1 && S(i,j) == 0 % 开始降温
            fall_count = fall_count + 1;
            fall_time = fall_time + dt;
        end
        
        if S(i,j-1) == 0 && S(i,j) == 0 % 继续降温
            fall_time = fall_time + dt;
        end
        
        if S(i,j-1) == 1 && S(i,j) == 1 % 继续升温
            rise_time = rise_time + dt;
            duty_cycle = duty_cycle + dt;
        end
        
    end
    
    avg_rise_time = (rise_time/rise_count)/60; % 平均升温时长 min
    avg_fall_time = (fall_time/fall_count)/60; % 平均降温时长 min
    period = (rise_time+fall_time)/(rise_count+fall_count)/60; % 周期 min
    avg_duty_cycle = duty_cycle/T*100; % 平均占空比 
    
    disp([num2str(T_out(i)),' ',num2str(avg_rise_time),' ',num2str(avg_fall_time),' ',num2str(period),' ',num2str(avg_duty_cycle),' ',num2str(E_day(i)),' ',num2str(P_avg(i)),' ',num2str(cost_day(i))])
    
end
% 第一问第（3）小问
% 假设供暖期为 180 天，

T_out_avg = [0 -5 -10 -15 -20]; % 室外平均温度 ℃
days = [30 40 40 40 30]; % 持续天数
E_season = zeros(length(T_out_avg),1); % 供暖期用电量 kWh
cost_season = zeros(length(T_out_avg),1); % 供暖期用电成本 元

for i = 1:length(T_out_avg)
    E_season(i) = E_day(i)*days(i); % 计算供暖期用电量
    cost_season(i) = cost_day(i)*days(i); % 计算供暖期用电成本
end

E_total = sum(E_season); % 计算供暖期总用电量
cost_total = sum(cost_season); % 计算供暖期总成本

% 输出结果
disp('表 2 供暖期典型住户用电量和用电成本统计结果')
disp(['室外平均温度/℃',' ','持续天数',' ','用电量/kWh',' ','供暖成本/元'])
for i = 1:length(T_out_avg)
    disp([num2str(T_out_avg(i)),' ',num2str(days(i)),' ',num2str(E_season(i)),' ',num2str(cost_season(i))])
end
disp(['供暖期总用电量/kWh',' ','供暖期总成本/元'])
disp([num2str(E_total),' ',num2str(cost_total)])
