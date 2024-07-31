clc;clear;
[TrainData,TestData] = Re_data_gene11(0,0,0,0,0);

Nn = 3000;  %训练样本个数
my = 5;     %观测样本y的维数
n_ns = 2;     %非平稳因子的个数
n_s = 3;     %平稳因子的个数
Nt = 1600;      %测试样本个数
time = 800;    %产生故障时间

[Ps, Pn, As, An, ssa_results] = ssa(TrainData', n_s);

%2.slice
N = 50;
sample_N = floor((sample-1)/N) + 1 ;                %每个epoch中sample数
sample_N_last = sample - sample_N * (N-1);          %最后epoch中sample数
epoch = cell(N,1);                                %存储每个epoch所含的smple    sample_N * dim
for j = 1:(N-1)
    epoch{j,1} = TrainData((j-1)*sample_N+1:j*sample_N,:);
end
epoch{N,1} = TrainData((N-1)*sample_N+1:(N-1)*sample_N+sample_N_last,:);


miu = cell(N,1);                             
covi = cell(N,1);                    
sum_miu = zeros(dim,1);
sum_cov = zeros(dim,dim);
for i = 1:N
    miu{i,1} = mean(epoch{i,1})';               %    dim*1
    sum_miu = sum_miu + miu{i,1};
    covi{i,1} = cov(epoch{i,1});
    sum_cov = sum_cov + (covi{i,1})^(1/2);
end

miu_average = sum_miu./N;
cov_average_r = (sum_cov./N)^2;

%我加的
T1 = zeros(Nn+Nt,1);
Bs_w_r = Ps;
covx_average = cov(TrainData);


for i = 1:Nn

    ss1 = Ps*TrainData(i,:)';        %3*1
    T1(i) = ((ss1-Bs_w_r*miu_average)'*(Bs_w_r*covx_average*Bs_w')^(-1)*(ss1-Bs_w_r*miu_average));%^(1/2); 
    %T1(k) = ((ss1-Bs_w*mean(TrainData)')'*(Bs_w*cov(TrainData)*Bs_w')^(-1)*(ss1-Bs_w*mean(TrainData)'))%^(1/2); 
    
end
%whos T1;
range = KDE_fcn(T1(1:Nn),0.99);
CTR1 = range(2);

for i = 1:Nt

    ss1 = Ps*TestData(i,:)';        %3*1
   % T1(i) = ((ss1-Bs_w*miu_average)'*(Bs_w*cov(TrainData)*Bs_w')^(-1)*(ss1-Bs_w*miu_average))^(1/2); 
    T1(i+Nn) = ((ss1-Ps*mean(TrainData).')'*(Ps*cov(TrainData)*Ps')^(-1)*(ss1-Ps*mean(TrainData).'));
end
error1 = 0;
for i = 1:time
    if T1(i+Nn) > CTR1%ceshi 
        error1 = error1+1;
    end
   
end
FAR1 = error1/time*100;


error1 = 0;

for i = time+1:Nt
    if T1(i+Nn) > CTR1
        error1 = error1+1;
    end
    
end
FDR1 = error1/(Nt - time)*100;





%8.求平稳非平稳源
sst =  TestData * Ps';             %估计的源
ssn =  TestData * Pn';


%8.绘制平稳源数据图
figure('Name','平稳信号');
subplot(3,1,1);
plot(1:Nt,sst(:,1),'k');
hold on;
title('平稳信号变化图');
ylabel('s1');

subplot(3,1,2);
plot(1:Nt,sst(:,2),'k');
hold on;

ylabel('s2');

subplot(3,1,3);
plot(1:Nt,sst(:,3),'k');
hold on;

xlabel('采样数');
ylabel('x3');

%9.绘制非平稳源数据图
figure('Name','非平稳信号');
subplot(2,1,1);
plot(1:Nt,ssn(:,1),'k');
hold on;

title('非平稳信号变化图');
ylabel('u1');

subplot(2,1,2);
plot(1:Nt,ssn(:,2),'k');
hold on;

xlabel('采样数');
ylabel('u2');



%10.绘制训练数据图
figure('Name','训练数据');
subplot(5,1,1);
plot(1:Nn,TrainData(:,1),'k');
hold on;
title('统计量变化图');
ylabel('x1');

subplot(5,1,2);
plot(1:Nn,TrainData(:,2),'k');
hold on;
ylabel('x2');

subplot(5,1,3);
plot(1:Nn,TrainData(:,3),'k');
hold on;
ylabel('x3');

subplot(5,1,4);
plot(1:Nn,TrainData(:,4),'k');
hold on;
ylabel('x4');

subplot(5,1,5);
plot(1:Nn,TrainData(:,5),'k');
hold on;
xlabel('采样数');
ylabel('x5');

range = KDE_fcn(T1(1:Nn),0.99);
CTR1 = range(2);









%我加的
T1 = zeros(Nn+Nt,1);

covx_average = cov(TrainData);

for i = 1:Nn

    ss1 = Ps*TrainData(i,:)';        %3*1
    T1(i) = ((ss1-Ps*mean(TrainData).')'*(Ps*covx_average*Ps')^(-1)*(ss1-Ps*mean(TrainData).'));%^(1/2); 
    
end

range = KDE_fcn(T1(1:Nn),0.99);
CTR1 = range(2);

for i = 1:Nt

    ss1 = Ps*TestData(i,:)';        %3*1 
    T1(i+Nn) = ((ss1-Ps*mean(TrainData).')'*(Ps*cov(TrainData)*Ps')^(-1)*(ss1-Ps*mean(TrainData).'));
end
error1 = 0;
for i = 1:time
    if T1(i+Nn) > CTR1%ceshi 
        error1 = error1+1;
    end
   
end
FAR1 = error1/time*100;


error1 = 0;

for i = time+1:Nt
    if T1(i+Nn) > CTR1
        error1 = error1+1;
    end
    
end
FDR1 = error1/(Nt - time)*100;


%11检测量绘图
%SSA
figure('Name','检测量变化图');
subplot(2,1,1);
plot(1:Nt,T1(Nn+1:Nn+Nt),'k');
hold on;
title(sprintf('误报率:%.1f%%    检测率:%.1f%%',FAR1,FDR1));
xlabel('采样时间');
ylabel('T^2');
hold on;
line([0,Nt],[CTR1,CTR1],'LineStyle','--','Color','r');
line([time,time],[0,20],'LineStyle','--','Color','b');
legend('SSA的T^2','控制限','故障发生');


