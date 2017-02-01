T = 100;
D = 1000;
l2 = log(2);

He = @(x) x.*log(x) + (1 - x).*log(1 - x);
H2 = @(x) l2*(x.*log2(x) + (1 - x).*log2(1 - x));

for t = 1:T
    p = rand(D,D);
    assertElementsAlmostEqual(He(p), H2(p));
end
