% 第四问
% 以电采暖住宅区 600 个住户为分析对象，设各住户初始室内温度在温控区间内均匀分布，在表 1 所示的各室外平均温度下，自行选定一组电采暖设备开关的初始状态
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

M = 600; % 住户数目
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

T_in0 = unifrnd(18,22,1,M); % 室内初始温度 ℃,随机生成在18到22之间的均匀分布的数值
T_wall0 = [17.32 16.865 16.38 15.88 15.365 14.85]; % 墙体初始稳态温度 ℃ 
S0 = randi([0 1],1,M); % 开关初始状态 随机生成0或1

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
                    
                    T_in_down = T_in_down + dt*(P_heat_down/C_in-(T_in_down-T_wall_down)/(C_in*R_1));%欧拉法求解室内温度
                    
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
    ylim([1 600])  % 设置y轴刻度范围为1到600
    title(['室外温度为',num2str(T_out(k)),'℃时'])
end

figure(2)
for k=1:length(T_out)
    subplot(3,2,k)
    plot(time,down_index(:,:,k),'o')
    xlabel('时间/min')
    ylabel('可参与下调的电采暖设备序号')
    ylim([1 600])  % 设置y轴刻度范围为1到600
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

figure(5) % 24小时内600个住户的总用电功率曲线
for k=1:length(T_out)
    subplot(3,2,k)
    plot(time,sum(P_heat(:,:,k))/1000)
    xlabel('时间/min')
    ylabel('总用电功率/kW')
    title('24小时内600个住户的总用电功率曲线')
end

%第五问
% 分时电价及辅助服务补偿机制

peak_price=0.56;%峰时电价 元/kWh

valley_price=0.32;%谷时电价 元/kWh

peak_compensation=1.30;%削峰补偿价格 元/kWh

valley_compensation=0.65;%填谷补偿价格 元/kWh

peak_start=16*60;%削峰开始时间 min

peak_end=20*60;%削峰结束时间 min

valley_start=0*60;%填谷开始时间 min

valley_end=4*60;%填谷结束时间 min

% 收益和成本初始化

benefit_per=zeros(M,N,length(T_out));% 每人每时收益 元

cost_per=zeros(M,N,length(T_out));% 每人每时成本 元

benefit=zeros(1,length(T_out));%总收益 元

cost=zeros(1,length(T_out));%总成本 元

% 计算参与调节后的室内温度变化和舒适度指标，并计算收益和成本

T_in_new = zeros(M,N,length(T_out)); % 参与调节后的室内温度 ℃

T_wall_new = zeros(M,N,length(T_out)); % 参与调节后的墙体温度 ℃

S_new = zeros(M,N,length(T_out)); % 参与调节后的开关状态

P_heat_new = zeros(M,N,length(T_out)); % 参与调节后的制热功率 W

deviation = zeros(1,length(T_out)); % 温度偏离度 ℃

unstability = zeros(1,length(T_out)); % 温度稳定性 ℃/min

change_num = zeros(M, N, length(T_out)); % 各时点由于参与电网调节导致开、关状态发生变化的电采暖设备数量

peak_power = zeros(1,length(T_out));

valley_power = zeros(1,length(T_out));

for k = 1:length(T_out)
    
    peak_power(1,k) = max(down_power(1,16*60+1:20*60,k)); % 削峰时段可提供的持续最大向下调节功率值 W
    
    valley_power(1,k) = max(up_power(1,1:4*60,k)); % 填谷时段可提供的持续最大向上调节功率值 W


for i = 1:M
    
    T_in_new(i,1,k) = T_in0(i); % 初始条件
    
    T_wall_new(i,1,k) = T_wall0(k); % 墙体初始稳态温度 ℃
    
    S_new(i,1,k) = S0(i); % 初始状态为开启或关闭
    
    P_heat_new(i,1,k) = S_new(i,1,k)*P_N; % 初始制热功率
    
    for j = 2:N
        
        T_wall_new(i,j,k) = T_wall_new(i,j-1,k) + dt*((T_in_new(i,j-1,k)-T_wall_new(i,j-1,k))/(C_wall*R_1)-(T_wall_new(i,j-1,k)-T_out(k))/(C_wall*R_2));%欧拉法求解墙体温度
        
        if (peak_start <= j && j <= peak_end && down_index(i,j,k) == i) || T_in_new(i,j-1,k) >= 22 % 削峰时段且可以向下调节或者室内温度大于等于温控区间上限
            
            S_new(i,j,k) = 0; % 关闭电采暖设备
            
            
        elseif (valley_start <= j && j <= valley_end && up_index(i,j,k) == i) || T_in_new(i,j-1,k) <= 18 % 填谷时段且可以向上调节或者室内温度小于等于温控区间下限
            
            S_new(i,j,k) = 1; % 开启电采暖设备
            
        else
            
            S_new(i,j,k) = S(i,j,k); % 其他情况保持原状态
            
        end
        
        if peak_start <= j && j <= peak_end && S_new(i,j-1,k) ~= S_new(i,j,k) % 在削峰时段如果开、关状态发生变化
            
            change_num(i,j,k) = 1; % 在相应的位置加一
        
        elseif valley_start <= j && j <= valley_end && S_new(i,j-1,k) ~= S_new(i,j,k)% 在填谷时段如果开、关状态发生变化
            
            change_num(i,j,k) = 1; % 在相应的位置加一
            
        end
        
        P_heat_new(i,j,k) = S_new(i,j,k)*P_N; %计算制热功率
        
        T_in_new(i,j,k) = T_in_new(i,j-1,k) + dt*(P_heat_new(i,j,k)/C_in-(T_in_new(i,j-1,k)-T_wall_new(i,j-1,k))/(C_in*R_1));%欧拉法求解室内温度
        
        if   4*60 < j && j < 8*60  %（平时段中，04:00——08:00；21:00——24:00为谷时，08:00——16:00；20:00——21:00为峰时）
        
                 cost_per(i,j,k) = cost_per(i,j,k) + P_heat_new(i,j,k)/1000*valley_price*dt/60; % 平时成本
                 
        elseif   8*60 <= j && j < 16*60
        
                 cost_per(i,j,k) = cost_per(i,j,k) + P_heat_new(i,j,k)/1000*peak_price*dt/60; % 平时成本
        
        elseif   20*60 < j && j <= 21*60
        
                 cost_per(i,j,k) = cost_per(i,j,k) + P_heat_new(i,j,k)/1000*peak_price*dt/60; % 平时成本
        
        elseif   21*60 < j && j <= 24*60
        
                 cost_per(i,j,k) = cost_per(i,j,k) + P_heat_new(i,j,k)/1000*valley_price*dt/60; % 平时成本
        end
            
        deviation(k) = deviation(k) + abs(T_in_new(i,j,k)-T_in0(1,i));%累计温度偏离度,与初始室内温度的差的绝对值
        
        unstability(k) = unstability(k) + abs(T_in_new(i,j,k)-T_in_new(i,j-1,k))/dt;%累计温度不稳定性
        
    end
    
end
end

%计算削峰填谷的成本和收益
for k = 1:length(T_out)
    for i = 1:M
        in_peak = false; % 标记是否处于削峰时段
        in_valley = false; % 标记是否处于填谷时段
        j = 2; % 时间索引
        
        while peak_start <= j && j <= peak_end
            cost_per(i,j,k) = cost_per(i,j,k) + P_heat_new(i,j,k)/1000*valley_price*dt/60; %削峰期的供热成本，算总成本用
            if S_new(i, j-1, k) == 1 && S_new(i, j, k) == 0 && ~in_peak % 检查是否满足开始削峰的条件
                in_peak = true; % 进入削峰时段
            end
            
            if in_peak % 如果处于削峰时段
                benefit_per(i, j, k) = benefit_per(i, j, k) + P_N/1000 * peak_compensation * dt / 60; % 削峰收益
                
                if S_new(i, j, k) == 1 % 检查是否满足结束削峰的条件
                    in_peak = false; % 结束削峰时段
                end
            end
            
            j = j + 1; % 移动到下一个分钟
        end
        
        while  valley_start <= j && j <= valley_end
            if S_new(i, j-1, k) == 0 && S_new(i, j, k) == 1 && ~in_valley % 检查是否满足开始填谷的条件
                in_valley = true;
            end
             if in_valley % 如果处于填谷时段
                benefit_per(i, j, k) = benefit_per(i, j, k) + P_N/1000 * valley_compensation * dt / 60; % 填谷收益
                cost_per(i,j,k) = cost_per(i,j,k) + P_N/1000*valley_price*dt/60; % 填谷成本(填谷要开机)
                
                if S_new(i, j, k) == 0 % 检查是否满足结束填谷的条件
                    in_peak = false; % 结束削峰时段
                end
             end
             j = j + 1; % 移动到下一个分钟
        end
    end  
end
%计算总收益和总成本

for k = 1:length(T_out)
    
    benefit = sum(sum(benefit_per, 1), 2); % 求和;
    benefit = squeeze(benefit); % 去掉单一维度
    
    cost = sum(sum(cost_per, 1), 2); % 求和;
    cost = squeeze(cost);
    
    change_permin = sum(change_num, 1); % 对第一维求和
end

%改变矩阵形状，为作图准备
benefit = benefit';
cost = cost';

% 绘制结果

figure(6)

% 参与调节后的室内温度曲线

for k=1:length(T_out)

subplot(3,2,k)

plot(time,T_in(:,:,k),'b',time,T_in_new(:,:,k),'r')

xlabel('时间/min')

ylabel('室内温度/℃')

legend('不参与调节','参与调节')

title(['室外温度为',num2str(T_out(k)),'℃时'])
end



% 收益和成本柱状图
figure(7)
bar([benefit' cost'])
xlabel('室外温度/℃')
ylabel('金额/元')
legend('收益','成本')
set(gca,'xticklabel',num2str(T_out'))
title('参与削峰填谷的收益和成本')


%参与削峰填谷对舒适度的影响绘图
figure(8)

subplot(2,1,1)

bar(deviation', 'b')  % 使用蓝色（blue）绘制柱状图

xlabel('室外温度/℃')
ylabel('指标值')
legend('温度偏离度')
set(gca,'xticklabel',num2str(T_out'))
title('参与削峰填谷对舒适度的影响')

subplot(2,1,2)

bar(unstability', 'g')  % 使用绿色（green）绘制柱状图

xlabel('室外温度/℃')
ylabel('指标值')
legend('温度不稳定性')
set(gca,'xticklabel',num2str(T_out'))
title('参与削峰填谷对舒适度的影响')

% 第五问第（3）小问
% 假设供暖期为 180 天
%要估算各室外温度下该住宅区电采暖负荷参与削峰填谷的上调、下调功率，需要把总时间更换

days = [30 40 40 40 30]; % 持续天数对应0，-5，-10，-15，-20℃
cost_per_season = zeros(1,5); %季度成本
benefit_per_season = zeros(1,5); %季度收益
for k = 1:length(T_out)-1
    cost_per_season(1,k) = days(k)*cost(1,k);  %季度成本计算
    benefit_per_season(1,k) = days(k)*benefit(1,k); %季度收益计算
end
cost_year = sum(cost_per_season(1,:))/(180/365); %年总成本计算
benefit_year = sum(benefit_per_season(1,:)); %年总收益计算
benefit_year_person = benefit_year / M; %平均每户的年收益
cost_savingratio = (benefit_year/cost_year)*100;%节省的供热成本百分比
