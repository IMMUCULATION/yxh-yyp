function [Rsresults,mar_peak]=ReSample20260908(xnfs,xns,Nsmo)
if Nsmo==1
%     nortmp=xnfs-mean(xnfs);
nortmp=xnfs;
else
%    nortmp=smooth(xnfs-mean(xnfs),Nsmo);
nortmp=smooth(xnfs,Nsmo);
end

% nor_smo=nortmp./abs(hilbert(nortmp));
nor_smo=nortmp;
tmp1=[0;nor_smo(1:(length(nor_smo)-1))];
tmp_mul=nor_smo.*tmp1;


negtiveN = find((tmp_mul<0));
%%直接极值
% for i=1:(length(negtiveN))
%    
%       
%         martmp(i)=negtiveN(i);
%    
% end


%%过零后极值
for i=1:(length(negtiveN)-1)
    if (nor_smo(negtiveN(i)))<0
        [aa,bb]=min(nor_smo(negtiveN(i):negtiveN(i+1)));
        martmp(i)=bb+negtiveN(i)-1;
    else
        [aa,bb]=max(nor_smo(negtiveN(i):negtiveN(i+1)));
        martmp(i)=bb+negtiveN(i)-1;
    end
end

%%极值+过零
% for i=1:(length(negtiveN)-1)
%     if (nor_smo(negtiveN(i)))<0
%         [aa,bb]=min(nor_smo(negtiveN(i):negtiveN(i+1)));
%         martmp(2*i-1)=bb+negtiveN(i)-1;
%     else
%         [aa,bb]=max(nor_smo(negtiveN(i):negtiveN(i+1)));
%         martmp(2*i-1)=bb+negtiveN(i)-1;
%     end
% end
% for i=1:(length(negtiveN)-1)-1
%     martmp(2*i)= round(martmp(2*i-1)/2+martmp(2*i+1)/2);
% end

%%极值后过零
% for i=1:(length(negtiveN)-1)
%     if (nor_smo(negtiveN(i)))<0
%         [aa,bb]=min(nor_smo(negtiveN(i):negtiveN(i+1)));
%         martmp0(i)=bb+negtiveN(i)-1;
%     else
%         [aa,bb]=max(nor_smo(negtiveN(i):negtiveN(i+1)));
%         martmp0(i)=bb+negtiveN(i)-1;
%     end
% end
% 
% for i=1:(length(negtiveN)-1)-1
%     martmp(i)= round(martmp0(i)/2+martmp0(i+1)/2);
% end

mar_peak=martmp';
Rsresults=xns(mar_peak);