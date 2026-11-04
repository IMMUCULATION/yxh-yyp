function [freqout,specout,lcmax,pcmax]=CZT(sig,L1,L2,m,L)
%% 输入：
% sig 重采样后的测量路信号
% L1 L2 范围
% m 细分数
% L 参考光纤长度
%% 输出：
% freqout： 横轴，距离
% specout：强度谱 

%%
w = exp(-1j*2*pi*(L2-L1)/(m*L)); 
  a = exp(1j*2*pi*L1/L); 
specout=czt(sig,m,w,a);
[pc,lc]=max(abs(specout));
lcmax=(lc-1)/m*(L2-L1)+L1;
pcmax=20*log10(pc);

freqout=(L1+(L2-L1)/m*(0:m-1))';

