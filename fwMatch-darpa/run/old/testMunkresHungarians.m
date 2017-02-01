function testMunkresHungarians(N, T)

    for t = 1:T
        A = 2*rand(N, N) - 1;
        [ah, vh] = hungarian(A);
        [am, vm] = munkres(A');
        
        assertElementsAlmostEqual(vh, vm);
        assertElementsAlmostEqual(ah, am);        
    end


end

