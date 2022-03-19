function codeA = encodefilter(A,cb1,cb2)
    for i = 1:length(A)  
        
        lsf = poly2lsf(A(i,:));
       
        distance = sum((lsf' - cb1).^2, 2);
        distance = distance/length(A);
        idx1 = find(distance == min(distance));
        
        residual = lsf' - cb1(idx1,:);
        
        distance = sum((residual - cb2).^2, 2);
        distance = distance/length(A);
        idx2 = find(distance == min(distance));
        
        codeA(i, :) = [idx1, idx2] ;
    end
end