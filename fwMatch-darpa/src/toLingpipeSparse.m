function v = toLingpipeSparse(x)
    D = length(x);
    [~, k, v] = find(x);
    v = com.aliasi.matrix.SparseFloatVector(int32(k - 1), v, D);    
end
