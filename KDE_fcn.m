function range=KDE_fcn(x,a)
pd = fitdist(x,'Kernel');
t=fmincon(@(x)icdf(pd,x)-icdf(pd,x-a),a,[],[],[],[],[a],[1]);   %a是你设定的置信概率
range=[icdf(pd,t-a), icdf(pd,t)];  %最短置信区间