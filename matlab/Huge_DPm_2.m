close all
clear all
Laux=24.964926;%fuzhu光纤光程

NUM=25;%每个位置点数

for l=1:3
for j=1:10 
mar_zero=[];
martmp=[];nor_smo=[];nortmp=[];tmp=[];tmp_mul=[];mar_peak=[];
fid=fopen(['H:\tsinghua\laser lidar\USS8802 SDK光盘 21.3\VC x64\Release\',num2str(l),'\huge_',num2str(j-1),'.bin']);%aux辅助干涉
[huge0,huge]=fread(fid,'int16');
fclose(fid);
clear CH1 CH2 CH3 CH4
n = length(huge0)/4;
idx = 1:4:length(huge0);
CH1 = huge0(idx);
CH2 = huge0(idx + 1);
CH3 = huge0(idx + 2);
CH4 = huge0(idx + 3);

MAXTAG=length(CH1)/64/65536;%位置个数
SEG1=zeros(NUM,MAXTAG);
SEG2=zeros(NUM,MAXTAG);
SEG=zeros(NUM,MAXTAG);

SEG11 = zeros(NUM, MAXTAG);
SEG12 = zeros(NUM, MAXTAG);
SEG21 = zeros(NUM, MAXTAG);
SEG22 = zeros(NUM, MAXTAG);
Signal_up   = cell(NUM, MAXTAG);
Signal_down = cell(NUM, MAXTAG);
AuxIndex_up   = cell(NUM, MAXTAG);
AuxIndex_down = cell(NUM, MAXTAG);
Spectrum     = cell(NUM, MAXTAG, 2, 2);
SpectrumAxis = cell(NUM, MAXTAG, 2, 2);
Window = zeros(NUM, MAXTAG, 4);
CyclePeak = zeros(1, MAXTAG);

% figure,plot(CH1(1:length(CH1)/MAXTAG))







for tag=1:MAXTAG


xnfs0=CH1(length(CH1)*(tag-1)/MAXTAG+1:length(CH1)*(tag)/MAXTAG);
xns0=CH3(length(CH3)*(tag-1)/MAXTAG+1:length(CH3)*(tag)/MAXTAG);
xnfs0=xnfs0-mean(xnfs0);
xns0=xns0-mean(xns0);


stp=150252;

[~,fp]=max(xnfs0(1:stp));
CyclePeak(tag) = fp;
str=fp+0.05*stp; sto=fp+0.45*stp; 
str1=fp+0.55*stp; sto1=fp+0.95*stp;
for i=1:2
if(i==1)
L1=7.39;L2=7.41; %CZT范围
else
L1=6.645;L2=6.665; %CZT范围
end
for k=1:25
xnfs=xnfs0(str+stp*(k-1):sto+stp*(k-1));
xnfs1=xnfs0(str1+stp*(k-1):sto1+stp*(k-1));
xns=xns0(str+stp*(k-1):sto+stp*(k-1));
xns1=xns0(str1+stp*(k-1):sto1+stp*(k-1));


[xns_cal_peak,mar_peak]=ReSample20260908(xnfs,xns,3);%！
[xns_cal_peak1,mar_peak1]=ReSample20260908(xnfs1,xns1,3);%！


% 两个目标共用这些信号，只记录一次。
% 只在保存副本时转换类型，不改变后续计算变量。
if i == 1
    Signal_up{k,tag}   = single(xns_cal_peak(:));
    Signal_down{k,tag} = single(xns_cal_peak1(:));

    AuxIndex_up{k,tag}   = int32(mar_peak(:));
    AuxIndex_down{k,tag} = int32(mar_peak1(:));

    Window(k,tag,1) = str  + stp*(k-1);  % 上扫起点表达式
    Window(k,tag,2) = sto  + stp*(k-1);  % 上扫终点表达式
    Window(k,tag,3) = str1 + stp*(k-1);  % 下扫起点表达式
    Window(k,tag,4) = sto1 + stp*(k-1);  % 下扫终点表达式
end


% [freqout,specout0,lcmax,pcmax]=CZT(xns_cal_peak,L1,L2,80000,Laux);
% pktmp=peakget(freqout,L1,L2,20*log10(2*abs(specout0/length(xns_cal_peak))));
% RS1(k,j)=pktmp; 
% [~, peakIndex] = min(abs(freqout - pktmp));
% localIndex = max(1,peakIndex-512): ...
%              min(length(freqout),peakIndex+512);
% 
% Spectrum{k,tag,i,1} = single(specout0(localIndex));
% SpectrumAxis{k,tag,i,1} = freqout(localIndex);
% 
% 
% [freqout1,specout1,lcmax,pcmax]=CZT(xns_cal_peak1,L1,L2,80000,Laux);
% pktmp=peakget(freqout1,L1,L2,20*log10(2*abs(specout1/length(xns_cal_peak1))));
% RS2(k,j)=pktmp; %原始信号 ,用作对比
% [~, peakIndex] = min(abs(freqout1 - pktmp));
% localIndex = max(1,peakIndex-512): ...
%              min(length(freqout1),peakIndex+512);
% 
% Spectrum{k,tag,i,2} = single(specout1(localIndex));
% SpectrumAxis{k,tag,i,2} = freqout1(localIndex);
% 
% 
% RS(k,j)=RS1(k,j)/2+RS2(k,j)/2;
% span(j)=(length(xns_cal_peak)+length(xns_cal_peak1))/2*3e8/2/Laux;


end
% if(i==1)
% SEG11(:,tag)=RS1(:,j);
% SEG12(:,tag)=RS2(:,j);
% SEG1(:,tag)=RS(:,j);
% else
% SEG21(:,tag)=RS1(:,j);
% SEG22(:,tag)=RS2(:,j);
% SEG2(:,tag)=RS(:,j);
% 
% end
end
% figure, hold on
% plot(xnfs)
% 
% ttmp=1:length(xnfs);
% plot(ttmp(mar_peak),xnfs(mar_peak),'ro')
% figure, hold on
% plot(xnfs1)
% 
% ttmp=1:length(xnfs1);
% plot(ttmp(mar_peak1),xnfs1(mar_peak1),'ro')
end
PositionIndex = 0:MAXTAG-1;
NominalDisplacement_m = -PositionIndex * 1e-6;

SampleRate_Hz = 75e6;
NominalModulation_Hz = 500;


% 距离、参数和索引
% save(['group',num2str(l),'_scan',num2str(j-1),'_info.mat'], ...
%     'l', 'j', ...
%     'PositionIndex', 'NominalDisplacement_m', ...
%     'SEG1', 'SEG2', 'SEG11', 'SEG12', 'SEG21', 'SEG22', ...
%     'Window', 'CyclePeak', ...
%     'Laux', 'NUM', 'stp', 'MAXTAG', ...
%     'SampleRate_Hz', 'NominalModulation_Hz', '-v7');

% 重采样信号和辅助极值索引
save(['group',num2str(l),'_scan',num2str(j-1),'_signals.mat'], ...
    'Signal_up', 'Signal_down', ...
     '-v7');
save(['group',num2str(l),'_scan',num2str(j-1),'_auxindexs.mat'], ...
    'AuxIndex_up', 'AuxIndex_down', '-v7');
% 局部复数谱及其距离坐标
% save(['group',num2str(l),'_scan',num2str(j-1),'_spectra.mat'], ...
%     'Spectrum', 'SpectrumAxis', '-v7');

end
end
% sspe=20*log10(2*abs(fft(xns_cal_peak/length(xns_cal_peak))));
% sspe=sspe(1:length(sspe)/2);
% ffre=(0:length(sspe)-1)/length(sspe)*Laux/2;
% figure,plot(ffre,sspe)
% set(0,'defaultfigurecolor','w')
% xlim([0 Laux/2])
% set(gca,'FontSize',25,'FontName', 'Times New Roman')
% xlabel("Distance(m)","FontSize",25,'FontName', 'Times New Roman');
% ylabel('Power(dB)','FontSize',25,'FontName', 'Times New Roman')
% 
% figure,plot(freqout,20*log10(2*abs(specout0/length(xns_cal_peak))))
% set(gca,'FontSize',25,'FontName', 'Times New Roman')
% xlabel("Distance(m)","FontSize",25,'FontName', 'Times New Roman');
% ylabel('Power(dB)','FontSize',25,'FontName', 'Times New Roman')

% res1=mean(SEG1);
% res2=mean(SEG2);
% figure,plot(res1-res2)
% figure,plot(res1)
% figure,plot(res2)
