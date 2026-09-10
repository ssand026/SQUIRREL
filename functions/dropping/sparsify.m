function [out] = sparsify(M)
% Checks if the memory/computation costs improve before converting the input to sparse. 
arguments
	M (:,:) double {mustBeNonsparse}
end
flM = M;
spM = sparse(flM);

full_mem = whos('flM').bytes;
sparse_mem = whos('spM').bytes;
%-----------------------------------------------------------
test_time = 0.25;
is_sparse = 0;
if (sparse_mem < full_mem)
    fl_time = 0; 
    sp_time = 0;
    fl_vec = zeros(size(flM,2),1);
    sp_vec = zeros(size(flM,2),1);
    
    iter = 0;
    tic;
    while (toc<test_time) || ((iter<=2)&&(toc<2*test_time))
        test_vect = randn(size(flM,2),1);
        time_a = toc;
        fl_vec = flM*test_vect;
        time_b = toc;
        sp_vec = spM*test_vect;
        time_c = toc;
        fl_time = fl_time + (time_b-time_a);
        sp_time = sp_time + (time_c-time_b);
        iter = iter+1;
    end
    if sp_time <= fl_time; is_sparse = 1; end
end
%-----------------------------------------------------------
% output sparse matrix if tests were successful
if is_sparse; out=spM; else; out=flM; end
end