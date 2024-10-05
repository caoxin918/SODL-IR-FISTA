function [D, X, Xout] = ODL(Y, G, lambda, opts, method)
% * Solving the following problem:
%  (D, X) = \arg\min_{D,X} 0.5||Y - DX||_F^2 + lambda||X||_1
% * Syntax: `[D, X] = ODL(Y, k, lambda, opts, sc_method)`
%   - INPUT: 
%     + `Y`: collection of samples.4/7/2016 7:35:39 PM
%     + `k`: number of atoms in the desired dictionary.
%     + `lambda`: norm 1 regularization parameter.
%     + `opts`: option.
%     + `sc_method`: sparse coding method used in the sparse coefficient update. Possible values:
%       * `'fista'`: using FISTA algorithm. See also [`fista`](#fista).
%       * `'spams'`: using SPAMS toolbox [[12]](#fn_spams). 
%   - OUTPUT:
%     + `D, X`: as in the problem.
% -----------------------------------------------
% Author: Tiep Vu, thv102@psu.edu, 4/7/2016
%         (http://www.personal.psu.edu/thv102/)
% -----------------------------------------------
	% if nargin == 0
		addpath(fullfile('DICTOL-master', 'utils'));
        addpath(fullfile('DICTOL-master', 'build_spams'));
        % d      = 50; % data dimension
        % N      = 100; % number of samples 
        % k      = 50; % dictionary size 
        % lambda = 0.1;
        % Y      = normc(rand(d, N));        
        %
        % 500次迭代一般足够
		opts.max_iter      = 500;
		opts.show_progress = 0;
		opts.check_grad    = false;  
		opts.tol           = 1e-8;  
		opts.verbose     = true;
		
	% end 
	%%
	opts = initOpts(opts);
	%%
	%% ========= initial D ==============================
	%D = PickDfromY(Y, [0, size(Y,2)], k);
    D = G;
    for ii = 1:size(G,2)
        D(:,ii)=D(:,ii)/norm(D(:,ii));
    end
    X = zeros(size(D,2), size(Y,2));
    Xout = X;
    if opts.verbose 
        fprintf('cost: %f', ODL_cost(Y, D, X, lambda));
    end 
    optsX = opts;
	optsX.max_iter = 200;
	optsX.tol      = 1e-8;
    optsD = opts;
	optsD.max_iter = 100;
	optsD.tol      = 1e-5;
	iter = 0;
	while iter < opts.max_iter
		iter = iter + 1;
		%% ========= sparse coding step ==============================
		% X = lasso_fista(Y, D, X, lambda, optsX);
       	X = IR_FISTA(D, Y, lambda);
        % 将负数变为0
        X(X<0)=0;
        Xout = Xout + X;
        if opts.verbose 
			costX = ODL_cost(Y, D, X, lambda);
			fprintf('iter: %3d, costX = %5f\n', iter, costX)
		end 
		%% ========= dictionary update step ==============================
		F = X*X'; E = Y*X';
		D = ODL_updateD(D, E, F, optsD);
		%效果没有上面好
        %D = KSVD(Y, D, X);
        if opts.verbose 
			costD = ODL_cost(Y, D, X, lambda);
			fprintf('iter: %3d, costD = %5f\n', iter, costD)
		end 
    end
    % 优化输出，只保留大于某个数的元素
    Xout = adjustOutput(Xout, 10);
    sourceNum = 1;
    sparsity = 4;
    threshold1 = 1;
    threshold2 = 200;
    % 再稀疏阶段，获得稀疏后的向量和离群点坐标
    % Xout =  re_sparse(Xout, sourceNum, sparsity,threshold1, threshold2);

	%%
	if nargin == 0
		pause;
	end
end

function cost = ODL_cost(Y, D, X, lambda)
% -----------------------------------------------
% Author: Tiep Vu, thv102@psu.edu, 04/07/2016
%         (http://www.personal.psu.edu/thv102/)
% -----------------------------------------------
	cost = 0.5*normF2(Y - D*X) + lambda*norm1(X);
end 

function X = adjustOutput(X, num)
    X(X<num)=0;
end
