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