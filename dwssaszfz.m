clc;clear;
[TrainData,TestData] = Re_data_gene11(0,2.7,0,0,0);
Nn = 3000;  %训练样本个数
my = 5;     %观测样本y的维数
n_ns = 2;     %非平稳因子的个数
n_s = 3;     %平稳因子的个数
Nt = 1600;      %测试样本个数
time = 800;    %产生故障时间


% 创建矩阵
old_matrix1 = TrainData;
% 确定新矩阵的大小
new_rows1 = size(old_matrix1, 1) - 2;
new_cols1 = size(old_matrix1, 2) * 3 ;
% 初始化新的矩阵
new_matrix1 = zeros(new_rows1, new_cols1);
% 进行迭代
for i = 1:new_rows1
    for j = 1:new_cols1
        if j <= size(old_matrix1, 2)%如果j小于旧矩阵的列数
            new_matrix1(i, j) = old_matrix1(i, j);
        elseif j <= size(old_matrix1, 2) * 2 
            new_matrix1(i, j) = old_matrix1(i+1, j-size(old_matrix1, 2));
        else
            new_matrix1(i, j) = old_matrix1(i+2, j-2*size(old_matrix1, 2));
        end
    end
end


%对测试集
% 创建矩阵
old_matrix2 = TestData;
% 确定新矩阵的大小
new_rows2 = size(old_matrix2, 1) - 2;
new_cols2 = size(old_matrix2, 2) * 3 ;
% 初始化新的矩阵
new_matrix2 = zeros(new_rows2, new_cols2);
% 进行迭代
for i = 1:new_rows2
    for j = 1:new_cols2
        if j <= size(old_matrix2, 2)%如果j小于旧矩阵的列数
            new_matrix2(i, j) = old_matrix2(i, j);
        elseif j <= size(old_matrix2, 2) * 2 
            new_matrix2(i, j) = old_matrix2(i+1, j-size(old_matrix2, 2));
        else
            new_matrix2(i, j) = old_matrix2(i+2, j-2*size(old_matrix2, 2));
        end
    end
end

Nn = 2998;  %训练样本个数
my = 15;     %观测样本y的维数
n_ns = 6;     %非平稳因子的个数
n_s = 9;     %平稳因子的个数
Nt = 1598;      %测试样本个数
time = 800;    %产生故障时间

[sample,dim] = size(new_matrix1);
  
x_mean = mean(new_matrix1);
new_matrix1 = new_matrix1 - repmat(x_mean,sample,1);


%2.slice
N = 50;
sample_N = floor((sample-1)/N) + 1 ;                %每个epoch中sample数
sample_N_last = sample - sample_N * (N-1);          %最后epoch中sample数
epoch = cell(N,1);                                %存储每个epoch所含的smple    sample_N * dim
for j = 1:(N-1)
    epoch{j,1} = new_matrix1((j-1)*sample_N+1:j*sample_N,:);
end
epoch{N,1} = new_matrix1((N-1)*sample_N+1:(N-1)*sample_N+sample_N_last,:);


%3.epoch 
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

%3.1. frechet mean
K = 100;
covf = cov_average_r;
theta = 1;
k = 1;
while theta >= 1e-4 && k < K
    T = zeros(dim,dim);
    for i = 1:N
        T = T + covf^(-1/2)*((covf^(1/2)*covi{i,1}*covf^(1/2))^(1/2))*covf^(-1/2);
    end
    T = T/N;
    cov_next = T * covf * T;
    theta = (trace((covf-cov_next)'*(covf-cov_next)))^(1/2);
    k = k+1;
    covf = cov_next;
end

cov_frechet = cov_next;
P1 = zeros(dim,dim);
P2 = zeros(dim,dim);
P3 = zeros(dim,dim);
P4 = zeros(dim,dim);
P5 = zeros(dim,dim);

for i = 1:N                                              
    P1 = P1+miu{i,1}*miu{i,1}';
    P2 = P2+covi{i,1};
    P3 = P3+(covi{i,1}^0.5*cov_average_r*covi{i,1}^0.5)^(1/2);
end
P4 = miu_average*miu_average';
P5 = cov_average_r;
S = (P1+P2-2*P3)./N -P4+P5;


%5.特征值分解
[E_vector,E_value] = eig(S);                         
for i = 1:dim
    E_vector(:,i) = E_vector(:,i)./(E_vector(:,i)'*E_vector(:,i))^(1/2);
end                                                                     %使特征向量模为1
[E_value_sort,index] = sort(diag(E_value),'ascend');                    %用排列函数将特征值和特征向量升序排列
E_value_sort
E_vector_sort = E_vector(:,index)
% E_vector_sort(:,1)' * E_vector_sort(:,1)


%6.求解Bn，Bn
Bs_w = E_vector_sort(:,1:n_s)';
E_vector_sort_rev = fliplr(E_vector_sort);
Bn_w = E_vector_sort_rev(:,1:n_ns)';
%改进
T1 = zeros(Nn+Nt,1);
Bs_w_r = Bs_w;
covx_average = cov(new_matrix1);


for i = 1:Nn

    ss1 = Bs_w*new_matrix1(i,:)';        %3*1
    T1(i) = ((ss1-Bs_w_r*miu_average)'*(Bs_w_r*covx_average*Bs_w')^(-1)*(ss1-Bs_w_r*miu_average));%^(1/2); 
    %T1(k) = ((ss1-Bs_w*mean(TrainData)')'*(Bs_w*cov(TrainData)*Bs_w')^(-1)*(ss1-Bs_w*mean(TrainData)'))%^(1/2); 
    
end
%whos T1;
range = KDE_fcn(T1(1:Nn),0.99);
CTR1 = range(2);

for i = 1:Nt

    ss1 = Bs_w*new_matrix2(i,:)';        %3*1
   % T1(i) = ((ss1-Bs_w*miu_average)'*(Bs_w*cov(TrainData)*Bs_w')^(-1)*(ss1-Bs_w*miu_average))^(1/2); 
    T1(i+Nn) = ((ss1-Bs_w*mean(new_matrix1).')'*(Bs_w*cov(new_matrix1)*Bs_w')^(-1)*(ss1-Bs_w*mean(new_matrix1).'));
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

%7.求真实子空间和投影
A = [0.2,-0.7;0.3,0.1;-0.7,0.4;0.1,-0.6;-0.5,0.5];
B = [-0.6,0.2,0.5;0.1,0.4,-0.7;0.7,-0.3,0.8;-0.5,0.2,0.1;0.3,-0.4,0.6];
A_A = [B,A];                             %真实的A
B_r = A_A^(-1);                          %真实的B

Bs_r = B_r(1:3,:);                        %真实的Bs
Bn_r = (B_r(4:5,:));                        %真实的Bn


%8.求平稳非平稳源
sst =  new_matrix2 * Bs_w';             %估计的源
ssn =  new_matrix2 * Bn_w';
%sst_r =new_matrix2 * Bs_r';            %真实的源
%ssn_r = new_matrix2 * Bn_r';


%8.绘制平稳源数据图
figure('Name','平稳信号');
subplot(3,3,1);
plot(1:Nt,sst(:,1),'k');
hold on;
%plot(1:Nt,sst_r(:,1),'r');
%hold on;
title('平稳信号变化图');
ylabel('s1');

subplot(3,3,2);
plot(1:Nt,sst(:,2),'k');
hold on;
%plot(1:Nt,sst_r(:,2),'r');
%hold on;
ylabel('s2');

subplot(3,3,3);
plot(1:Nt,sst(:,3),'k');
hold on;
%plot(1:Nt,sst_r(:,3),'r');
%hold on;
xlabel('采样数');
ylabel('x3');

subplot(3,3,4);
plot(1:Nt,sst(:,4),'k');
hold on;
%plot(1:Nt,sst_r(:,3),'r');
%hold on;
xlabel('采样数');
ylabel('x3');

subplot(3,3,5);
plot(1:Nt,sst(:,5),'k');
hold on;
%plot(1:Nt,sst_r(:,3),'r');
%hold on;
xlabel('采样数');
ylabel('x5');

subplot(3,3,6);
plot(1:Nt,sst(:,6),'k');
hold on;
%plot(1:Nt,sst_r(:,3),'r');
%hold on;
xlabel('采样数');
ylabel('x6');

subplot(3,3,7);
plot(1:Nt,sst(:,7),'k');
hold on;
%plot(1:Nt,sst_r(:,3),'r');
%hold on;
xlabel('采样数');
ylabel('x7');

subplot(3,3,8);
plot(1:Nt,sst(:,8),'k');
hold on;
%plot(1:Nt,sst_r(:,3),'r');
%hold on;
xlabel('采样数');
ylabel('x8');

subplot(3,3,9);
plot(1:Nt,sst(:,9),'k');
hold on;
%plot(1:Nt,sst_r(:,3),'r');
%hold on;
xlabel('采样数');
ylabel('x9');
ylabel('s2');


%9.绘制非平稳源数据图
figure('Name','非平稳信号');
subplot(2,3,1);
plot(1:Nt,ssn(:,1),'k');
hold on;
%plot(1:Nt,ssn_r(:,1),'r');
%hold on;
title('非平稳信号变化图');
ylabel('u1');

subplot(2,3,2);
plot(1:Nt,ssn(:,2),'k');
hold on;
%plot(1:Nt,ssn_r(:,2),'r');
%hold on;
xlabel('采样数');
ylabel('u2');

subplot(2,3,3);
plot(1:Nt,ssn(:,3),'k');
hold on;
%plot(1:Nt,ssn_r(:,2),'r');
%hold on;
xlabel('采样数');
ylabel('u3');

subplot(2,3,4);
plot(1:Nt,ssn(:,4),'k');
hold on;
%plot(1:Nt,ssn_r(:,2),'r');
%hold on;
xlabel('采样数');
ylabel('u4');

subplot(2,3,5);
plot(1:Nt,ssn(:,5),'k');
hold on;
%plot(1:Nt,ssn_r(:,2),'r');
%hold on;
xlabel('采样数');
ylabel('u5');

subplot(2,3,6);
plot(1:Nt,ssn(:,6),'k');
hold on;
%plot(1:Nt,ssn_r(:,2),'r');
%hold on;
xlabel('采样数');
ylabel('u6');


%10.绘制训练数据图
figure('Name','训练数据');
subplot(5,1,1);
plot(1:Nn,new_matrix1(:,1),'k');
hold on;
title('统计量变化图');
ylabel('x1');

subplot(5,1,2);
plot(1:Nn,new_matrix1(:,2),'k');
hold on;
ylabel('x2');

subplot(5,1,3);
plot(1:Nn,new_matrix1(:,3),'k');
hold on;
ylabel('x3');

subplot(5,1,4);
plot(1:Nn,new_matrix1(:,4),'k');
hold on;
ylabel('x4');

subplot(5,1,5);
plot(1:Nn,new_matrix1(:,5),'k');
hold on;
xlabel('采样数');
ylabel('x5');

range = KDE_fcn(T1(1:Nn),0.99);
CTR1 = range(2);

%11检测量绘图
%WSSA
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
legend('DWSSA的T^2','控制限','故障发生');
