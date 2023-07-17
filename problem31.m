% 第三问
% 以 6个电采暖住户（序号分别为 1-6）为例，假设室外温度为-20℃，室内初始温度在温控区间内均匀分布，自行选定一组电采暖设备开关的初始状态
clc;clear;
% 参数设置
C_in = 1.1e6; % 室内空气等效热容 J/℃
C_wall = 1.86e8; % 墙体等效热容 J/℃
R_1 = 1.2e-3; % 室内空气和墙体内侧等效热阻 ℃/W
R_2 = 9.2e-3; % 墙体外侧和室外空气等效热阻 ℃/W
P_N = 8e3; % 电采暖设备额定功率 W

% 时间设置
dt =1; % 时间步长 1s
T = 24*3600; % 总时间 24h
N = T/dt; % 时间步数 
time = 0:dt:T-dt; % 时间序列 h

M = 6; % 住户数目
T_out_20 = -20; % 室外温度 ℃
T_in_20 = zeros(M,N); % 室内温度 ℃
T_wall_20 = zeros(M,N); % 墙体温度 ℃
S_20 = zeros(M,N); % 开关状态 
P_heat_20 = zeros(M,N); % 制热功率 W

up_time_20 = zeros(M,N); % 功率上调可持续时间 s
down_time_20 = zeros(M,N); % 功率下调可持续时间 s

T_in0_20 = linspace(18,22,M); % 室内初始温度 ℃
T_wall0_20 = [15.365 15.365 15.365 15.365 15.365 15.365]; % 墙体初始稳态温度 ℃
S0_20 = randi([0 1],1,M); % 开关初始状态

% 循环计算不同住户下的结果
for i = 1:M
    
    T_in_20(i,1) = T_in0_20(i); % 初始条件
    T_wall_20(i,1) = T_wall0_20(i); % 墙体初始稳态温度 ℃
    S_20(i,1) = S0_20(i); % 初始状态为开启或关闭
    P_heat_20(i,1) = S_20(i,1)*P_N; % 初始制热功率
    
    for j = 2:N
        
        T_wall_20(i,j) = T_wall_20(i,j-1) + dt*((T_in_20(i,j-1)-T_wall_20(i,j-1))/(C_wall*R_1)-(T_wall_20(i,j-1)-T_out_20)/(C_wall*R_2)); % 欧拉法求解墙体温度
        
        if 18 < T_in_20(i,j-1) < 22 % 温控逻辑判断开关状态
            S_20(i,j) = S_20(i,j-1);
        end
        
        if T_in_20(i,j-1) >= 22 
            S_20(i,j) = 0;
        end
        
        if T_in_20(i,j-1) <= 18 
            S_20(i,j) = 1;
        end
        
        P_heat_20(i,j) = S_20(i,j)*P_N; % 计算制热功率
        
        T_in_20(i,j) = T_in_20(i,j-1) + dt*(P_heat_20(i,j)/C_in-(T_in_20(i,j-1)-T_wall_20(i,j-1))/(C_in*R_1)); % 欧拉法求解室内温度
        
        if S_20(i,j) == 0 % 关闭状态下可以向上调节
            
            T_in_up = T_in_20(i,j); % 向上调节时的室内温度 ℃
            T_wall_up = T_wall_20(i,j); % 向上调节时的墙体温度 ℃
            P_heat_up = P_N; % 向上调节时的制热功率 W
            
            while T_in_up < 22
                
                up_time_20(i,j) = up_time_20(i,j) + dt; %累计向上调节可持续时间
                
                T_wall_up = T_wall_up + dt*((T_in_up-T_wall_up)/(C_wall*R_1)-(T_wall_up-T_out_20)/(C_wall*R_2));%欧拉法求解墙体温度
                
                T_in_up = T_in_up + dt*(P_heat_up/C_in-(T_in_up-T_wall_up)/(C_in*R_1)); % 欧拉法求解室内温度
            end
        end
         if S_20(i,j) == 1 % 开启状态下可以向下调节
            
            T_in_down = T_in_20(i,j); % 向下调节时的室内温度 ℃
            T_wall_down = T_wall_20(i,j); % 向下调节时的墙体温度 ℃
            P_heat_down = 0; % 向下调节时的制热功率 W
            
            while T_in_down > 18
                
                down_time_20(i,j) = down_time_20(i,j) + dt; %累计向下调节可持续时间
                
                T_wall_down = T_wall_down + dt*((T_in_down-T_wall_down)/(C_wall*R_1)-(T_wall_down-T_out_20)/(C_wall*R_2)); %欧拉法求解墙体温度
                
                T_in_down = T_in_down + dt*(P_heat_down/C_in-(T_in_down-T_wall_down)/(C_in*R_1)); %欧拉法求解室内温度
                
            end
            
         end
        
    end
    
end

% 绘制结果
figure(1)
for i=1:M
    subplot(2,3,i)
    plot(time/3600,T_in_20(i,:))
    xlabel('时间/h')
    ylabel('温度/℃')
    title(['住户',num2str(i),'的室内温度'])
    legend('室内温度')
end
figure(2)
for i=1:M
    subplot(2,3,i)
    plot(time/3600,T_wall_20(i,:))
    xlabel('时间/h')
    ylabel('温度/℃')
    title(['住户',num2str(i),'的墙体温度'])
    legend('墙体温度')
end

figure(3)
for i=1:M
    subplot(3,2,i)
    plot(time/3600,S_20(i,:))
    xlabel('时间/h')
    ylabel('开关状态')
    title(['住户',num2str(i),'的电采暖设备开关状态'])
    set(gca,'YTick',[0 1]) % 设置y轴刻度为0和1
end

figure(4)
for i=1:M
    subplot(2,3,i)
    plot(time/3600,up_time_20(i,:)/3600)
    hold on
    plot(time/3600,down_time_20(i,:)/3600)
    hold off
    xlabel('时间/h')
    ylabel('可持续时间/h')
    title(['用户',num2str(i),'的功率上调、下调可持续时间'])
    legend('向上调节','向下调节')
end

figure(5) % 24小时内6个住户的总用电功率曲线
plot(time/3600,sum(P_heat_20(:,:))/1000)
xlabel('时间/h')
ylabel('总用电功率/kW')
title('24小时内6个住户的总用电功率曲线')
