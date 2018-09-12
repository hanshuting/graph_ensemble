cd Differentiation
%mex trw_bprop_scheduled_helper1.cpp -I../eigen3/ CXXFLAGS="\$CFLAGS -O3"
%mex trw_bprop_scheduled_helper2.cpp -I../eigen3/ CXXFLAGS="\$CFLAGS -O3"
mex trw_bprop_scheduled_helper1.cpp -I../eigen3/ CXXFLAGS="\$CFLAGS -fopenmp" LDFLAGS="\$LDFLAGS -fopenmp"
mex trw_bprop_scheduled_helper2.cpp -I../eigen3/ CXXFLAGS="\$CFLAGS -fopenmp" LDFLAGS="\$LDFLAGS -fopenmp"
cd ..

% cd Differentiation
% mex trw_bprop_helper1.cpp -I../eigen3/
% mex trw_bprop_helper2.cpp -I../eigen3/
% mex meanfield_bprop_helper1.cpp -I../eigen3/
% mex meanfield_bprop_helper2.cpp -I../eigen3/
% cd ..
% 
% cd Inference
% mex meanfield_helper.cpp -I../eigen3/
% mex trw_helper.cpp      -I../eigen3/
% cd ..
% 
% cd Losses
% mex pseudo_helper.cpp -I../eigen3
% cd ..
% 
% % public version doesn't include sampling code
% % cd Sampling
% % mex gibbs_helper.cpp -I../eigen
% % cd ..
% 
% cd Vision
% mex reflectim_helper.cpp -I../eigen3
% cd ..
% 
