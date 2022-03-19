function Aq = decodefilter(codeA,cb1,cb2)
    for i = 1:size(codeA, 1)
        
        idx1 = codeA(i, 1);
        idx2 = codeA(i, 2);
        
        lsf = cb1(idx1, :);
        residual = cb2(idx2, :);
        
        output = sort(lsf + residual);
        
        Aq(i, :) = lsf2poly(output);
    end
end