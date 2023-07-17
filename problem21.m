% 第二问第（1）小问
% 以单个住户电采暖负荷为对象，室外温度为-15℃，室内初始温度为20℃
% 电采暖设备开关的初始状态为开启，计算典型住户电采暖负荷在日内24h各时点（间隔 1min）功率上调、下调的可持续时间，并绘制计算结果。
clc;clear;
% 参数设置
C_in = 1.1e6; % 室内空气等效热容 J/℃
C_wall = 1.86e8; % 墙体等效热容 J/℃
R_1 = 1.2e-3; % 室内空气和墙体内侧等效热阻 ℃/W
R_2 = 9.2e-3; % 墙体外侧和室外空气等效热阻 ℃/W
P_N = 8e3; % 电采暖设备额定功率 W
T_in0 = 20; % 室内初始温度 ℃

% 时间设置
dt =1; % 时间步长 1min
T = 24*3600; % 总时间 24h
N = T/dt; % 时间步数 
time = 0:dt:T-dt; % 时间序列 h

T_out_15 = -15; % 室外温度 ℃
T_in_15 = zeros(1,N); % 室内温度 ℃
T_wall_15 = zeros(1,N); % 墙体温度 ℃
S_15 = zeros(1,N); % 开关状态 
P_heat_15 = zeros(1,N); % 制热功率 W

T_in_15(1) = T_in0; % 初始条件
T_wall_15(1) = 15.88; % 墙体初始稳态温度 ℃
S_15(1) = 1; % 初始状态为开启
P_heat_15(1) = S_15(1)*P_N; % 初始制热功率

up_time = zeros(1,N); % 功率上调可持续时间 s
down_time = zeros(1,N); % 功率下调可持续时间 s

for j = 2:N
    
    T_wall_15(j) = T_wall_15(j-1) + dt*((T_in_15(j-1)-T_wall_15(j-1))/(C_wall*R_1)-(T_wall_15(j-1)-T_out_15)/(C_wall*R_2)); % 欧拉法求解墙体温度
    
    if 18 < T_in_15(j-1) < 22 % 温控逻辑判断开关状态
        S_15(j) = S_15(j-1);
    end
    
    if T_in_15(j-1) >= 22 
        S_15(j) = 0;
    end
    
    if T_in_15(j-1) <= 18 
        S_15(j) = 1;
    end
    
    P_heat_15(j) = S_15(j)*P_N; % 计算制热功率
    
    T_in_15(j) = T_in_15(j-1) + dt*(P_heat_15(j)/C_in-(T_in_15(j-1)-T_wall_15(j-1))/(C_in*R_1)); % 欧拉法求解室内温度
    
    if S_15(j) == 0 % 关闭状态下可以参与向上调节
        
        T_in_up = T_in_15(j); % 向上调节时的室内温度 ℃
        T_wall_up = T_wall_15(j); % 向上调节时的墙体温度 ℃
        P_heat_up = P_N; % 向上调节时的制热功率 W
        
        while T_in_up < 22
            
            up_time(j) = up_time(j) + dt; % 累计向上调节可持续时间
            
            T_wall_up = T_wall_up + dt*((T_in_up-T_wall_up)/(C_wall*R_1)-(T_wall_up-T_out_15)/(C_wall*R_2)); % 欧拉法求解墙体温度
            
            T_in_up = T_in_up + dt*(P_heat_up/C_in-(T_in_up-T_wall_up)/(C_in*R_1)); % 欧拉法求解室内温度
            
        end
        
    end
    
    if S_15(j) == 1 % 开启状态下可以参与向下调节
        
        T_in_down = T_in_15(j); % 向下调节时的室内温度 ℃
        T_wall_down = T_wall_15(j); % 向下调节时的墙体温度 ℃
        P_heat_down = 0; % 向下调节时的制热功率 W
        
        while T_in_down > 18
            
            down_time(j) = down_time(j) + dt; %累计向下调节可持续时间
            
            T_wall_down = T_wall_down + dt*((T_in_down-T_wall_down)/(C_wall*R_1)-(T_wall_down-T_out_15)/(C_wall*R_2)); %欧拉法求解墙体温度
            
            T_in_down = T_in_down + dt*(P_heat_down/C_in-(T_in_down-T_wall_down)/(C_in*R_1)); %欧拉法求解室内温度
            
        end
        
    end
    
end

% 绘制结果
figure(1)
plot(time/60,up_time/60)
hold on
plot(time/60,down_time/60)
hold off
xlabel('时间/min')
ylabel('可持续时间/min')
title('室外温度为-15℃时的功率上调、下调可持续时间')
legend('向上调节','向下调节')