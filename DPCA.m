clc
close all
X_train = load('D:\Desktop\TE_process\d00.dat'); %载入训练数据
X_test =load('D:\Desktop\TE_process\d02_te.dat'); %载入测试数据

X_train = double(X_train');%注意此处有个转置，是因为TE数据集的d00数据和其他数据行列颠倒了。只要保证和样本数据维度一致即可。
X_test = double(X_test);

%对训练集
% 创建矩阵
old_matrix1 = X_train;
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
old_matrix2 = X_test;
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


X_mean = mean(new_matrix1);  %求均值                           
X_std = std(new_matrix1);    %求标准差                      

[X_row,X_col] = size(new_matrix1);                                                    
new_matrix1=(new_matrix1-repmat(X_mean,X_row,1))./repmat(X_std,X_row,1);%标准化处理

SXtrain = cov(new_matrix1);%求协方差矩阵
[T,lm] = eig(SXtrain);%求特征值及特征向量                                                               
D = flipud(diag(lm));%将特征值从大到小排列，                            

num = 1;                                         
while sum(D(1:num))/sum(D) < 0.9   %若累计贡献率小于90%则增加主元个数
num = num +1;
end                                                 
P = T(:,X_col-num+1:X_col); %取对应的特征向量                           
                     
TU=num*(X_row-1)*(X_row+1)*finv(0.99,num,X_row - num)/(X_row*(X_row - num));%求置信度为99%时的T2统计控制限  

for i = 1:3
    th(i) = sum((D(num+1:X_col)).^i);
end
h0 = 1 - 2*th(1)*th(3)/(3*th(2)^2);
ca = norminv(0.99,0,1);
QU = th(1)*(h0*ca*sqrt(2*th(2))/th(1) + 1 + th(2)*h0*(h0 - 1)/th(1)^2)^(1/h0); %置信度为99%的Q统计控制限                          
%模型测试
n = size(new_matrix2,1);
new_matrix2=(new_matrix2-repmat(X_mean,n,1))./repmat(X_std,n,1);%标准化处理
%求T2统计量，Q统计量
[r,y] = size(P*P');
I = eye(r,y); 
T2 = zeros(n,1);
Q = zeros(n,1);
for i = 1:n
T2(i)=new_matrix2(i,:)*P*inv(lm(156-num+1:156,156-num+1:156))*P'*new_matrix2(i,:)';                                           
    Q(i) = new_matrix2(i,:)*(I - P*P')*new_matrix2(i,:)';                                                                                    
end
%绘图
for i=1:n
    if T2(i,1)>TU
        break
    end
end

figure('Name','DPCA');
subplot(3,1,1);
    plot(1:i,T2(1:i),'k');
    hold on;
    plot(i:n,T2(i:n),'k');
   title('统计量变化图');
    xlabel('采样数');
    ylabel('T2');
    hold on;
    line([0,n],[TU,TU],'LineStyle','--','Color','g');
   for i=1:n
    if Q(i,1)>QU
        break
    end
   end 

subplot(3,1,2);
    plot(1:i,Q(1:i),'k');
    hold on;
     plot(i:n,Q(i:n),'k');
    title('统计量变化图');
    xlabel('采样数');
    ylabel('SPE');
    hold on;
    line([0,n],[QU,QU],'LineStyle','--','Color','g');
    
%贡献图
%4.计算每个变量对Q的贡献
e = new_matrix2(600,:)*(I - P*P');%选取第600个样本来检测哪个变量出现问题。
contq = e.^2;
[mxq_bor,mxq]=max(contq);
contqq=zeros(1,156);
contqq(1,mxq)=contq(1,mxq);
%5. 绘制贡献图
subplot(3,1,3);
    bar(contq,'g');
     hold on;
   bar(contqq,'r');
    xlabel('变量号');
    ylabel('SPE贡献率 %');
%检测结果分析
%T2故障检测率
  k1=0;
for i=160:958
   if T2(i)>TU
       k1=k1+1;
   end
end
FDRT=k1/(958-160);
%T2误报率
k1=0;
for i=1:159
   if T2(i)>TU
       k1=k1+1;
   end
end
FART=k1/159;
%Q故障检测率
k1=0;
for i=160:958
   if Q(i)>QU
       k1=k1+1;
   end
end
FDRQ=k1/(958-160);
%Q故障误报率
k1=0;
for i=1:159
   if Q(i)>QU
       k1=k1+1;
   end
end
FARQ=k1/159;

