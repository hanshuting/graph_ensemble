function C = padcat_hor(A,B)

s1 = size(A,1);
s2 = size(B,1);

if s2 > s1
    A = [A; zeros(s2-s1,size(A,2))];
else
    B = [B; zeros(s1-s2,size(B,2))];
end

C = [A B];
        