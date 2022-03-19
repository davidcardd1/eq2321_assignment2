function outq = sq_dec(idx, n_bits, xmax, m)
    L = 2^n_bits - 1;
    del = 2*xmax/L;
    outq = del*idx + del/2 - m;
end