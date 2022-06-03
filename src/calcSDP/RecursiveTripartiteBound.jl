using Combinatorics

function Z(m1,m2)
    return floor(m1/2)*floor((m1-1)/2)*floor(m2/2)*floor((m2-1)/2)
end

function A(m1, m2, m3)
    return floor(m1/2)*floor((m1-1)/2)*floor(m2/2)*floor((m2-1)/2)+floor(m3/2)*floor((m3-1)/2)*floor(m1*m2/2) + floor(m1/2)*floor((m1-1)/2)*floor(m3/2)*floor((m3-1)/2)+floor(m2/2)*floor((m2-1)/2)*floor(m1*m3/2) + floor(m3/2)*floor((m3-1)/2)*floor(m2/2)*floor((m2-1)/2)+floor(m1/2)*floor((m1-1)/2)*floor(m3*m2/2)
end

function recurA(m1, m2, m3)
    return Z(m1, m2) + Z(m3, m1+m2)
end

## 
using Plots
plot([A(m,m,m) for m = 1:30])
plot!([recurA(m,m,m)*3/2 for m = 1:30])