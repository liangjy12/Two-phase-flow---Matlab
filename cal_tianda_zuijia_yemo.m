%批量绘制递归图，需要调用函数RP3.m
function cal_tianda_zuijia_yemo()
%clear all;
load('name_str');
load('data_5000_thickness.mat');
data_num = numel(name_str);

for k = 1:data_num
    tic
    data0 = name_str(k);
    data0 = char(data0)
    [non1,t_new,non2,m_new] = CC_new(eval(data0));
    fid1=fopen('m_t.txt','a');
    fprintf(fid1,'%s: m=%d, t=%d, td1=%d, tw=%d \r\n\r\n',data0,m_new,t_new,non1,non2);
    fclose(fid1);
    
    RP3(eval(data0),m_new,t_new,data0);
    %attractor(t_new,eval(data0),data0);
    %RQA(data0);
    t_runtime = toc
end
end

function [td1,td2,tw,m_new] = CC_new(x)
td1=0;
td2=0;
tw=0;
N = length(x);
sigma = std(x);
tmax = 100;
m = 2:5;
r = 1:4;
St = 0;
for tt = 1:tmax
    L = floor(N/tt);
    for mm = m
        for rr = r
            rj = sigma*rr/2;
            for ss = 1:tt;
                xsub=x(ss:tt:N);
                xsub=xsub(1:L);
                Cmrt(ss) = getC(xsub,mm,rj);
                C1rt(ss) = getC(xsub,1,rj);
            end
            St(tt,mm,rr) = (sum(Cmrt)-sum(C1rt.^mm))/tt;
        end
    end
S_t = sum(sum(St,3),2)/16;
%[tt,mm,rr]
end
for i = 2:tmax
    if S_t(i-1)*S_t(i)<=0
        td1 = i;
        break;
    end
end
dSt = max(St,[],3)-min(St,[],3);
dS_t = sum(dSt,2)/4;
Scor = dS_t + abs(S_t);
IndMin=find(diff(sign(diff(dS_t')))>0)+1;
if ~isempty(IndMin)
    td2 = IndMin(1);
end
[t1,tw] = min(Scor);
m_new = ceil(tw/td2)+1;
end

function Cmrt = getC(x,m,r)
    p = length(x) - (m-1);
    phase_new = zeros(m,p);
    for i=1:m
        phase_new(i,:)=x([((i-1)+1):1:((i-1)+p)]);
    end
    phase_new = phase_new';
    dd = pdist(phase_new,'chebychev');
    theta = sum(dd<=r);
    Cmrt = 2*theta/(p*(p-1));
end

%批量绘制递归图时使用的函数
function RP3(x,m,t,data0)
%tic
clf;
p = length(x) - (m-1)*t;   %嵌入相空间向量个数为p， p=N0 - (m-1)t, N0为时间序列x数据点数

%取p行m列的矩阵y来表示重构的相空间
y = zeros(m,p);
for i=1:m
    y(i,:)=x([((i-1)*t+1):1:((i-1)*t+p)]);
end
y = y';

%计算各行之间的距离，并生成一个p*p的距离矩阵
D = pdist2(y,y,'chebychev');
clear y;

%计算阈值
r = 0.2*std(x);
D2 = D<r;          %阈值判断，小于阈值r则置1，大于r则置0
clear D;
D2 = double(D2);

%生成每一行均为1-p的矩阵
H = 1:1:p;
DD = ones(p,1)*H;
% parfor k = 1:1:p
%     DD(k,:) = H; 
% end
yaxis = DD.*(D2);   %将矩阵置为对应位置的列数
clear DD;
clear D2;

%h = inputname(1);
h = data0;
filename = sprintf('%s%s',h,'_y');   %将纵坐标矩阵yaxis保存在filename中
save(filename,'yaxis');
xaxis = 1:1:p;
plot(xaxis,yaxis,'k.');
set(gcf,'Name',h);
title(h);
saveas(gcf,h,'bmp');
%t = toc
end

function RQA(name_data1)
    fid2=fopen('RQA.txt','a');
    name_y = sprintf('%s%s',name_data1,'_y');    
    fprintf(fid2,'Recurrence Quqntification Analysis of %s:\r\n\r\n',name_y);
    load(name_y);
    N = length(yaxis);
    
    %非零点数，即递归点数
    nzero_num = numel(yaxis(yaxis~=0));
    
    %递归率（Recurrence Rate，RR)计算
    RR = nzero_num/(N^2);
    fprintf(fid2,'RR = %d\r\n',RR);
    
    %确定性（Determinism，DET）计算,l_min = 3
    pl = 0;
    m = 0;
    for d = -(N-3):1:(N-3)
        l = 0;
        if d == 0
            continue;
        end
        dia = diag(yaxis,d);
        for k = 1:(N-abs(d))
            if dia(k) ~=0
                l = l+1;
            elseif dia(k) ==0 & l<3
                l = 0;
            elseif dia(k) ==0 & l>=3
                m = m+1;
                pl(m) = l;
                l = 0;
            end
        end
        if l>=3
            m = m+1;
            pl(m) = l;
        end
    end
    DET = sum(pl)/nzero_num;
    fprintf(fid2,'DET = %d\r\n',DET);
    
    %比率（Ratio）计算
    Ratio = DET/RR;
    fprintf(fid2,'Ratio = %d\r\n',Ratio);
    
    %平均对角线长度（Average diagonal line length，L）
    L = sum(pl)/m;
    fprintf(fid2,'L = %d\r\n',L);
    
    %最大对角线长度（Maximum length，Lmax）
    Lmax = max(pl);
    fprintf(fid2,'Lmax = %d\r\n',Lmax);
    
    %分叉性（Divergence，DIV）
    DIV = 1/Lmax;
    fprintf(fid2,'DIV = %d\r\n',DIV);
    
    %熵（Entropy，ENTR）计算
    if pl==0
        fprintf(fid2,'ENTR: pl = 0, have no diag\r\n\r\n');
    elseif pl~=0
    pl_without_re = unique(pl);
    ENTR = 0;
    for q = pl_without_re
        ppl(q) = numel(pl(pl==q))/m;
        ENTR = ENTR + ppl(q)*log(ppl(q));
    end
    ENTR = -ENTR;
    fprintf(fid2,'ENTR = %d\r\n\r\n',ENTR);
    end
    
    clearvars -except name_str_y fid2;
fclose(fid2);
end

function attractor(t_new,data0,name_data0)
%clearvars -except tt;
clf;
%tt = 8;
    xy = data0;
    for i = 1:(5000-t_new)
        y(i) = xy(i+t_new);
    end
    x = xy(1:(5000-t_new))';
    plot(x,y);
    
    name_attractor = sprintf('%s%s',name_data0,'_attractor');
    set(gcf,'Name',name_attractor);
    title(name_attractor);
    xlabel('x(t)');
    label0 = sprintf('x(t+%d)',t_new);
    ylabel(label0);
    saveas(gcf,name_attractor,'bmp');
end