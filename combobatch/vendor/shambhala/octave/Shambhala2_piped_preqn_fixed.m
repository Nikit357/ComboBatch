% Shambhala2_piped_preqn_fixed.m
% Combined variant: no quantilenorm (input was pre-QN'd Python-side) AND
% CuBlock_fixed with precomputed gene clusters.
%
% Added by ComboBatch. Upstream had the two speed-ups only as separate scripts,
% and its script selection let the fixed-clusters branch overwrite the pre-QN
% branch, so the approved A+E combination quantile-normalized twice (once in
% Python, once here). This script is the missing fourth combination.
%
% Required --eval variables: NH, NP, k, FIXED_CLUSTERS (G×1 int32, labels 1..k).

inData=readExpressionData('/dev/stdin','log2');
Exp = inData.Samples;
SYMBOL = inData.GeneList;
SN = inData.SamplesName;
NS = length(SN);

for ( i = 1:NH )
    message = sprintf('Harmonizing sample %d out of %d',i,NH);
    disp(message);

    i0 = i;
    for ( jjj = 1:NP )
        i0 = [i0 (jjj+NH)];
    end

    EXP = Exp(:,i0);
    % quantilenorm deliberately omitted: input was pre-QN'd by Python

    dataN = CuBlock_fixed(real(EXP),[],k,FIXED_CLUSTERS);

    log2e = log2(exp(1));
    DataN = dataN/log2e;

    EXPN = exp(DataN)-1;
    vecN = EXPN(:,1);

    if ( i == 1 )
        OUT = vecN;
    else
        OUT = [OUT vecN];
    end

end

nG=size(OUT,1);
nS=size(OUT,2);
for i=1:nG
    fprintf(1,'%s',SYMBOL{i,1});
    for j=1:nS
        fprintf(1,' %f',OUT(i,j));
    end
    fprintf(1, '\n');
end
