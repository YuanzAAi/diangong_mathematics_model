% 第三问
% 以 6个电采暖住户（序号分别为 1-6）为例，假设室内初始温度在温控区间内均匀分布，自行选定一组电采暖设备开关的初始状态
clc;clear;
% 参数设置
C_in = 1.1e6; % 室内空气等效热容 J/℃
C_wall = 1.86e8; % 墙体等效热容 J/℃
R_1 = 1.2e-3; % 室内空气和墙体内侧等效热阻 ℃/W
R_2 = 9.2e-3; % 墙体外侧和室外空气等效热阻 ℃/W
P_N = 8e3; % 电采暖设备额定功率 W

% 时间设置
dt =1; % 时间步长 1min
T = 24*60; % 总时间 24h
N = T/dt; % 时间步数 
time = 0:dt:T-dt; % 时间序列 h

M = 6; % 住户数目，算600户的可以改为600，记得把绘图那里的坐标轴也改一下
T_out = [0 -5 -10 -15 -20 -25]; % 室外温度 ℃
T_in = zeros(M,N,length(T_out)); % 室内温度 ℃
T_wall = zeros(M,N,length(T_out)); % 墙体温度 ℃
S = zeros(M,N,length(T_out)); % 开关状态 
P_heat = zeros(M,N,length(T_out)); % 制热功率 W

up_time = zeros(M,N,length(T_out)); % 功率上调可持续时间 min
down_time = zeros(M,N,length(T_out)); % 功率下调可持续时间 min

up_index = zeros(M,N,length(T_out)); % 可参与上调的电采暖设备序号
down_index = zeros(M,N,length(T_out)); % 可参与下调的电采暖设备序号

up_power = zeros(1,N,length(T_out)); % 总可上调功率 W
down_power = zeros(1,N,length(T_out)); % 总可下调功率 W

T_in0 = linspace(18,22,M); % 室内初始温度 ℃
T_wall0 = [17.32 16.865 16.38 15.88 15.365 14.85]; % 墙体初始稳态温度 ℃
S0 = randi([0 1],1,M); % 开关初始状态 

% 循环计算不同住户下的结果
for k = 1:length(T_out) % 遍历不同的室外温度
    
    for i = 1:M
        
        T_in(i,1,k) = T_in0(i); % 初始条件
        T_wall(i,1,k) = T_wall0(k); % 墙体初始稳态温度 ℃
        S(i,1,k) = S0(i); % 初始状态为开启或关闭
        P_heat(i,1,k) = S(i,1,k)*P_N; % 初始制热功率
        
        for j = 2:N
            
            T_wall(i,j,k) = T_wall(i,j-1,k) + dt*((T_in(i,j-1,k)-T_wall(i,j-1,k))/(C_wall*R_1)-(T_wall(i,j-1,k)-T_out(k))/(C_wall*R_2)); % 欧拉法求解墙体温度
            
            if 18 < T_in(i,j-1,k) < 22 % 温控逻辑判断开关状态
                S(i,j,k) = S(i,j-1,k);
            end
            
            if T_in(i,j-1,k) >= 22 
                S(i,j,k) = 0;
            end
            
            if T_in(i,j-1,k) <= 18 
                S(i,j,k) = 1;
            end
            
            P_heat(i,j,k) = S(i,j,k)*P_N; %计算制热功率
            
            T_in(i,j,k) = T_in(i,j-1,k) + dt*(P_heat(i,j,k)/C_in-(T_in(i,j-1,k)-T_wall(i,j-1,k))/(C_in*R_1)); % 欧拉法求解室内温度
            
            if S(i,j,k) == 0 % 关闭状态下可以向上调节
                
                T_in_up = T_in(i,j,k); % 向上调节时的室内温度 ℃
                T_wall_up = T_wall(i,j,k); % 向上调节时的墙体温度 ℃
                P_heat_up = P_N; % 向上调节时的制热功率 W
                up_index(i,j,k) = i;% 如果可以向上调节则记录该住户的序号
                
                while T_in_up < 22
                    
                    up_time(i,j,k) = up_time(i,j,k) + dt; %累计向上调节可持续时间
                    
                    T_wall_up = T_wall_up + dt*((T_in_up-T_wall_up)/(C_wall*R_1)-(T_wall_up-T_out(k))/(C_wall*R_2));%欧拉法求解墙体温度
                    
                    T_in_up = T_in_up + dt*(P_heat_up/C_in-(T_in_up-T_wall_up)/(C_in*R_1)); % 欧拉法求解室内温度
                end
                
                
            end
            
            if S(i,j,k) == 1 % 开启状态下可以向下调节
                
                T_in_down = T_in(i,j,k); % 向下调节时的室内温度 ℃
                T_wall_down = T_wall(i,j,k); % 向下调节时的墙体温度 ℃
                P_heat_down = 0; % 向下调节时的制热功率 W
                down_index(i,j,k) = i; % 如果可以向下调节，则记录该住户的序号
                
                while T_in_down > 18
                    
                    down_time(i,j,k) = down_time(i,j,k) + dt; %累计向下调节可持续时间
                    
                    T_wall_down = T_wall_down + dt*((T_in_down-T_wall_down)/(C_wall*R_1)-(T_wall_down-T_out(k))/(C_wall*R_2)); %欧拉法求解墙体温度
                    
                    T_in_down = T_in_down + dt*(P_heat_down/C_in-(T_in_down-T_wall_down)/(C_in*R_1)); %欧拉法求解室内温度
                    
                end
                
                
            end
            
        end
        
    end
    
up_power(:,:,k) = sum(up_index(:,:,k)~=0,1)*P_N; % 计算总可上调功率 w
down_power(:,:,k) = sum(down_index(:,:,k)~=0,1)*P_N; % 计算总可下调功率 w
    
end

% 绘制结果
figure(1)
for k=1:length(T_out)
    subplot(3,2,k)
    plot(time,up_index(:,:,k),'o')
    xlabel('时间/min')
    ylabel('可参与上调的电采暖设备序号')
    ylim([1 6]) % 设置y轴范围为1到6
    set(gca,'YTick',[ 1 2 3 4 5 6]) % 设置y轴刻度为1到6
    title(['室外温度为',num2str(T_out(k)),'℃时'])
end

figure(2)
for k=1:length(T_out)
    subplot(3,2,k)
    plot(time,down_index(:,:,k),'o')
    xlabel('时间/min')
    ylabel('可参与下调的电采暖设备序号')
    ylim([1 6]) % 设置y轴范围为1到6
    set(gca,'YTick',[ 1 2 3 4 5 6]) % 设置y轴刻度为1到6
    title(['室外温度为',num2str(T_out(k)),'℃时'])
end

figure(3)
for k=1:length(T_out)
    subplot(3,2,k)
    plot(time,up_power(:,:,k)/1000)
    xlabel('时间/min')
    ylabel('总可上调功率/kW')
    title(['室外温度为',num2str(T_out(k)),'℃时'])
end

figure(4)
for k=1:length(T_out)
    subplot(3,2,k)
    plot(time,down_power(:,:,k)/1000)
    xlabel('时间/min')
    ylabel('总可下调功率/kW')
    title(['室外温度为',num2str(T_out(k)),'℃时'])
end