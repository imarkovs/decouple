function decoupledPolynomial = DecouplePolynomial(coupledPolynomial)
% This function decouples a coupled polynomial function and gives back a
% decoupled representation of this polynomial, using tensor
% decomposition techniques.
% Written by Gabriel Hollander, Vrije Universiteit Brussel, Dept. ELEC, 10/06/2016

[coupledPolynomial, modelInfo] = checkInputStructure(coupledPolynomial);

% Construct the Jacobian tensor and decompose it using CPD
[WeWeights, VeWeights, HeWeights, ystarWeights, relerr, it, ~, ~, fullWeightMatrix_var, ~, Jstar, equations_H, exit_crit, costfunction_relerr] = ...
    CPdecoupleInOperatingPoints_weighted(coupledPolynomial.uOperatingPoints, coupledPolynomial.approximatedCoeffs, coupledPolynomial.monomials, coupledPolynomial.varSym, modelInfo, coupledPolynomial.covarianceMatrix, coupledPolynomial.CPDtype, coupledPolynomial.oldLoad, coupledPolynomial.lambda, [coupledPolynomial.stepSizeTol_wTx, coupledPolynomial.maxNwTx, coupledPolynomial.wTxCommandWindowLog]);

% reconstructing ge from tensor decomposition with the Nway toolbox and {u*,y*} data
[geWeights, x] = ComputeInternalMappings(modelInfo, WeWeights, VeWeights, HeWeights, coupledPolynomial.uOperatingPoints, ystarWeights);
GeWeights = cell(coupledPolynomial.r, 1); for k = 1 : coupledPolynomial.r, GeWeights{k} = fliplr(CoefficientList(geWeights(k), x(k), coupledPolynomial.degreeMax).'); end

% take the appropriate relative error to use in the output structure
if isempty(costfunction_relerr), relative_error = relerr(end); else relative_error = costfunction_relerr(end); end

decoupledPolynomial = struct(...
    'type', ['decoupling with ' coupledPolynomial.CPDtype ' weight'], ...
    'We', WeWeights, ...
    'Ve', VeWeights, ... %     'He', HeWeights, ...
    'Ge', PadLeft(GeWeights), ...
    'iteration_count', it, ...
    'relerr', relative_error);
end

function [coupledPolynomial, modelInfo] = checkInputStructure(coupledPolynomial)
n = 2; m = size(coupledPolynomial.approximatedCoeffs,2);
if ~isfield(coupledPolynomial, 'covarianceMatrix'), coupledPolynomial.covarianceMatrix = []; coupledPolynomial.CPDtype = 'no'; 
elseif isempty(coupledPolynomial.covarianceMatrix)
    coupledPolynomial.CPDtype = 'no';
end
if ~isfield(coupledPolynomial, 'NOperatingPoints'), coupledPolynomial.NOperatingPoints = 50; end
if ~isfield(coupledPolynomial, 'uOperatingPoints'), coupledPolynomial.uOperatingPoints = rand(n, coupledPolynomial.NOperatingPoints); end
switch size(coupledPolynomial.approximatedCoeffs,1)
    case 3,  coupledPolynomial.degreeMax = 1;
    case 6,  coupledPolynomial.degreeMax = 2;
    case 10, coupledPolynomial.degreeMax = 3;
    case 15, coupledPolynomial.degreeMax = 4;
    case 21, coupledPolynomial.degreeMax = 5;
    case 28, coupledPolynomial.degreeMax = 6;
    case 36, coupledPolynomial.degreeMax = 7;
    case 45, coupledPolynomial.degreeMax = 8;
end
modelInfo.degreeMax = coupledPolynomial.degreeMax; modelInfo.branchesNumber = coupledPolynomial.r; 
[coupledPolynomial.monomials, coupledPolynomial.varSym]= Monomials('u', n, coupledPolynomial.degreeMax);
if ~isfield(coupledPolynomial, 'oldLoad'), coupledPolynomial.oldLoad = randOldLoad(m, n, coupledPolynomial.r, coupledPolynomial.NOperatingPoints, 2); end
if ~isfield(coupledPolynomial, 'lambda'), coupledPolynomial.lambda = 1; end
if ~isfield(coupledPolynomial, 'stepSizeTol_wTx'), coupledPolynomial.stepSizeTol_wTx = 1e-06; end
if ~isfield(coupledPolynomial, 'maxNwTx'), coupledPolynomial.maxNwTx = 20; end
if ~isfield(coupledPolynomial, 'wTxCommandWindowLog'), coupledPolynomial.wTxCommandWindowLog = 50; end

end

function [ge, x] = ComputeInternalMappings(modelInfo, We, Ve, He, ustar, ystar)
% Last version: 22/10/2014, updated by Gabriel Hollander

% [ge,gVTu,feu] = compute_internal_mappings( model, data )
% INPUTS
% model.V       =    CPD factor V
% model.W       =    CPD factor W
% model.H       =    CPD factor H 
% model.d       =    max degree
% model.x       =    symbolic vars internal
% model.u       =    symbolic vars inputs
% data.ustar    =    experimental data inputs    
% data.ystar    =    experimental data outputs
% 
% requires vdmvec.m
%
% Philippe Dreesen, 2014
% Vrije Universiteit Brussel, Dept. ELEC
%

d = modelInfo.degreeMax;
r = modelInfo.branchesNumber; 
x = sym('x', [r 1]);
re = size( Ve, 2 );
n = size( We, 1 );

% xstare
xstare = Ve.'*ustar;

% 1) fitting polynomials through he_i(x_i) = ge_i'(x_i)
he = sym( zeros( re, 1 ));
ge0 = sym( zeros( re, 1 ));

for i = 1:re,
    % fit degree d-1 polynomial
    hetmp = polyfit( xstare(i,:)', He(:,i), d-1 );
    vx = vdmvec(x(i), d-1);
    he(i) = flipud( hetmp(:) ).' * vx;
    
    % reconstructing g_i(x_i) (by integrating)
    ge0(i) = int( he(i));
    ge0(i) = ge0(i) - subs( ge0(i), x(i), 0);       % remove constant term
    
end

% 2) determining integration intercepts
K = 1:1; 
ystarlong = zeros( n * length(K), 1 );
Wlong = repmat( We, [length(K), 1 ]);
for ki = 1:length(K),
    k = K(ki);
    g0Veustar = zeros( re, 1 );
    for j=1:re,
        g0Veustar(j) = subs( ge0(j), x(j), xstare(j,k));
    end
    ystarlong((ki-1)*n+1:ki*n,1) = ystar(:,k) - We*g0Veustar;
end

% determine int constants from LS problem
gec = Wlong\ystarlong;
ge = ge0 + gec;

end

function [We, Ve, He, ystar, relerr, it, weightTensor_var, weightTensor, fullWeightMatrix_var, fullWeightMatrix, Jstar, equations_H, exit_crit, costfunction_relerr] = CPdecoupleInOperatingPoints_weighted(operatingPoints, coefficients, monomials, varSym, modelInfo, covarianceMatrix, diag_or_blockdiag, oldLoad, lambda_weightForSecondSetOfEqs, optionsForCPD)
% This function is based on the working code of Philippe Dreessen and
% decouples the tensor made from the jacobian matrices in the operating
% points. In this version of the function, we use the Nway toolbox instead
% of tensorlab.
% Usage:
%         CPdecoupleInOperatingPoints(operatngPoints, coefficients, monomials)
% Written by Gabriel Hollander, Vrije Universiteit Brussel, Dept. ELEC
% 15/10/2014

% Is there weighting info or not?
if isempty(covarianceMatrix)
    withWeighting = 'noWeighting';
else
    switch diag_or_blockdiag
        case 'diag'
            withWeighting = 'withWeighting_diag';
        case 'blockdiag'
            withWeighting = 'withWeighting_blockdiag';
        case 'full'
            withWeighting = 'withWeighting_fullmatrix';
    end
end

% Get the general variables about input, output, branches
    m = size(operatingPoints, 1);       % the number of inputs
    n = size(coefficients, 2);          % the number of outputs
    r = modelInfo.branchesNumber;       % the number of branches
    N = size(operatingPoints, 2);       % the number of operating points
    lengthMonomials = length(monomials);
    
% The coupled function for each output
    f = sym(zeros(n, 1));               
    for k = 1 : n
        f(k) = sum(coefficients(:,k) .* monomials);
    end
    f_mf = matlabFunction( f, 'vars', varSym);
    
% The symbolic Jacobian matrix
    Ju = jacobian(f, varSym);
    Ju_mf = matlabFunction( Ju, 'vars', varSym );

% Initialization of the outputs and of the jacobian tensor
    ystar = zeros(n,N);
    Jstar = zeros(n,m,N);

% Evaluating f and Jf in ustar(k)    
    for k = 1 : N,
        uTemporary = num2cell(operatingPoints(:,k));
        Jstar(:,:,k) = Ju_mf(uTemporary{:});
        ystar(:,k) = f_mf(uTemporary{:});
    end

% Check if the rank estimation is correct
    rkest = rankest(Jstar);
%     if (rkest ~= r), 
%         display(['Rank estitation = ' num2str(rkest) ' is not correct: overriding by correct one = ' num2str(r)]); 
%     end
    re = r;
    
% Perform the Canonical Polyadic Decomposition on the Jacobian matrix
% special case of 2x2xN tensor: no need to do iterations, so maxIterations = 1
% Use my own implentation of the CPD.
    switch withWeighting
        case 'noWeighting'      % when no weights are used
            [WVH, it, relerr] = weightedcpd(Jstar, oldLoad, [], optionsForCPD);
            weightTensor_var = []; weightTensor = []; fullWeightMatrix_var = []; fullWeightMatrix = []; equations_H = []; exit_crit = []; costfunction_relerr = [];
        case 'withWeighting_diag'    % when weights are used, generate the diagonal weighting matrix and use it in the parafac decomposition
            monomialsDerived = sym(zeros(m, lengthMonomials));
            monomialsDerived_mf = cell(m, 1);
            for j_iter = 1 : m
                monomialsDerived(j_iter, :) = diff(monomials, varSym(j_iter)).';
                monomialsDerived_mf{j_iter} = matlabFunction(monomialsDerived(j_iter, :), 'vars', varSym);
            end
            
            weightTensor_var = zeros(size(Jstar));
            for k_iter = 1 : N
                for i_iter = 1 : n
                    for j_iter = 1 : m
                        uTemp = num2cell(operatingPoints(:, k_iter));
                        weightTensor_var(i_iter, j_iter, k_iter) = monomialsDerived_mf{j_iter}(uTemp{:}) *  covarianceMatrix((i_iter-1)*lengthMonomials+1 : i_iter*lengthMonomials, (i_iter-1)*lengthMonomials+1 : i_iter*lengthMonomials) * monomialsDerived_mf{j_iter}(uTemp{:}).';
                    end
                end
            end
            inverted_weightTensor_var = 1 ./ weightTensor_var;
            weightTensor = zeros(m*n, m*n, N);
            for opPoint_iter = 1 : N
               weightTensor(:,:,opPoint_iter) = diag(reshape(inverted_weightTensor_var(:,:,opPoint_iter), [m*n, 1]));
            end
            [WVH, it, relerr] = weightedcpd(Jstar, oldLoad, weightTensor, optionsForCPD);
            fullWeightMatrix_var = []; fullWeightMatrix = []; equations_H = []; exit_crit = []; costfunction_relerr = [];
        case 'withWeighting_blockdiag'
            % We define the derivatives of the monomials, to be used in the matrix A (row j contains the j-th derivative)
                monomialsDerived = sym(zeros(m, lengthMonomials));
                for j_iter = 1 : m
                    monomialsDerived(j_iter, :) = diff(monomials, varSym(j_iter)).';
                end
            % We define the matrix A, used in the propagation of uncertainty (wiki). In our case, this matrix has a lot of zeros, except where there is information about a certain derivative.
                A = sym(zeros(m*n, lengthMonomials*n));
                for row_iter = 1 : m*n
                    function_no = mod(row_iter, n);
                    if function_no == 0, function_no = n; end
                    derivative_to = ceil(row_iter/n); % this gives the column number of the tensor of Jacobians (the j-th derivative)
                    A(row_iter, (function_no-1)*lengthMonomials + 1: function_no*lengthMonomials) = monomialsDerived(derivative_to, :);
                end
            % A_mf is a functionized version of the matrix A, so we can use it with the different operating points using cells
              A_mf = matlabFunction(A, 'vars', varSym);
            % Construct the tensor of covariances: (m*n) x (m*n) x N
                weightTensor_var = zeros(m*n, m*n, N);
                weightTensor = zeros(m*n, m*n, N);
                for opPoint_iter = 1 : N
                    uTemp = num2cell(operatingPoints(:, opPoint_iter));
                    weightTensor_var(:,:, opPoint_iter) =  A_mf(uTemp{:})*covarianceMatrix*(A_mf(uTemp{:})).';
                    weightTensor(:,:, opPoint_iter) = pinv(weightTensor_var(:,:, opPoint_iter));
                end
                [WVH, it, relerr] = weightedcpd(Jstar, oldLoad, weightTensor, optionsForCPD);
                fullWeightMatrix_var = []; fullWeightMatrix = []; equations_H = []; exit_crit = []; costfunction_relerr = [];
        case 'withWeighting_fullmatrix'
            % We define the derivatives of the monomials, to be used in the matrix A (row j contains the j-th derivative)
                monomialsDerived = sym(zeros(m, lengthMonomials));
                for j_iter = 1 : m
                    monomialsDerived(j_iter, :) = diff(monomials, varSym(j_iter)).';
                end
            % We define the matrix A, used in the propagation of uncertainty (wiki). In our case, this matrix has a lot of zeros, except where there is information about a certain derivative.
                A = sym(zeros(m*n, lengthMonomials*n));
                for row_iter = 1 : m*n
                    function_no = mod(row_iter, n);
                    if function_no == 0, function_no = n; end
                    derivative_to = ceil(row_iter/n); % this gives the column number of the tensor of Jacobians (the j-th derivative)
                    A(row_iter, (function_no-1)*lengthMonomials + 1: function_no*lengthMonomials) = monomialsDerived(derivative_to, :);
                end
            % A_mf is a functionized version of the matrix A, so we can use it with the different operating points using cells
                A_mf = matlabFunction(A, 'vars', varSym);
            % Construct the full covariance matrix at once for all the operating points. For this, gather all the N different (m*n) x (l*n) A-matrices.
                A_cell = cell(N,1);
                for opPoint_iter = 1 : N
                    uTemp = num2cell(operatingPoints(:, opPoint_iter));
                    A_cell{opPoint_iter} = A_mf(uTemp{:});
                end
                A_combined = cell2mat(A_cell);
                fullWeightMatrix_var = A_combined*covarianceMatrix*(A_combined.');
             % Make the fullWeightMatrix_var symmetric (it's due to rounding-off errors in the order of e-20)
                fullWeightMatrix_var = (fullWeightMatrix_var + fullWeightMatrix_var.')/2;
                rank_fullWMatrix_var = rank(fullWeightMatrix_var);
%                 if n*lengthMonomials ~= rank_fullWMatrix_var, cprintf('red', '   The rank of the full weighting matrix is %d instead of %d\n', rank_fullWMatrix_var, n*lengthMonomials); end
%                 [U_A, ~, ~] = svd(A_combined);
%                 Aperp = (U_A(:,rank_fullWMatrix_var + 1:end)).';
                [WVH, it, relerr, costfunction_relerr, exit_crit] = wTxCPD(Jstar, oldLoad, fullWeightMatrix_var, lambda_weightForSecondSetOfEqs, optionsForCPD);
                weightTensor_var = []; weightTensor = []; fullWeightMatrix = []; equations_H = [];
    end

% Give the output to the function
    We = WVH{1};                    
    Ve = WVH{2};                    
    He = WVH{3};
end

function out = MakeBlockDiag(matrix, times)
% This function makes a block diagonal matrix by repeating the given matrix
% a given number of times.
% Written by Gabriel Hollander, Vrije Universiteit Brussel, Dept. ELEC, 30/3/2015
matrixCell = cell(1, times);
for iter = 1 : times
    matrixCell{iter} = matrix;
end
out = blkdiag(matrixCell{:});

end

function out = randOldLoad(m, n, r, NOperatingPoints, maxOldLoadElement)
% This function randomizes the oldLoad which is given to the Parafac
% algorithm as an initizal value.
% Written by Gabriel Hollander, Vrije Universiteit Brussel, Dept. ELEC, 18/3/2015
out = {RandomReal([-maxOldLoadElement, maxOldLoadElement], [n, r]), ...
       RandomReal([-maxOldLoadElement, maxOldLoadElement], [m, r]), ...
       RandomReal([-maxOldLoadElement, maxOldLoadElement], [NOperatingPoints, r])}; 
end

function out = reshapeWeightForMode(weightMatrix, mode, Ninputs, Noutputs, NOperatingPonits)
% This function takes the weight tensor used for the decoupling of
% multivariate polynomials in the decoupling procedure. It reshapes the
% weight tensor to be used in the Altenating Least Squares iterations for
% the weighted CP decomposition.
% It takes two inputs: the weight tensor and the mode used in the
% Alternating Least Squares.
% Written by Gabriel Hollander, dept. ELEC, VUB, 26/03/2015

m = Ninputs;
n = Noutputs;
N = NOperatingPonits;

% define a tensor of the appropriate dimensions in order to find the good permutation for each mode
T = reshape(1:m*n*N, [n,m,N]);
I = eye(m*n*N);

switch mode
    case 1
        reorder = vec(reshape(T, [n, m*N]).');
        P = I(reorder, :);      % due to the code of Nway toolbox (and Kolda)
        out = P*weightMatrix*P.';
    case 2
        reorder = vec(reshape(permute(T, [2, 1, 3]), [m, n*N]).');
        P = I(reorder, :);    % due to Kolda
        out = P*weightMatrix*P.';
    case 3
        reorder = vec(reshape(permute(T, [3, 1, 2]), [N, m*n]).');
        P = I(reorder, :);    % due to the code of Nway toolbox (and Kolda)
        out = P*weightMatrix*P.'; % or transposed, doens't matter: symmetric weight matrix
end
end

function vecske=vdmvec(xi,di,nosym)
% Build a Vandermonde structured (power product) vector.  
%
% SIGNATURE
% vecout = vdmvec(xi,di,nosym)
%
% DESCRIPTION
% Build a power product vector of a given symbolic variable xi containing
% the powers zero up to di. 
% 
% INPUTS
%    xi         =    symbolic input variable
%    di         =    desired max degree
%    nosym      =    flag for NOT symbolic result (1 means not symbolic)
%
% OUTPUTS
%    vecout     =    Vandermonde vector of size di+1 x 1 containing the
%                    powers xi^d, with d=0,1,...,di
%
% EXAMPLE
% >> x=sym('x');
% >> vdmvec(x,4)
% 
% ans =
%  
%    1
%    x
%  x^2
%  x^3
%  x^4
%  
% CALLS
% sym 
%
% AUTHOR
%     Philippe Dreesen (philippe.dreesen@gmail.com)
%     April 2014
%

if (nargin==2 || nosym==0),
    vecske = sym(ones(di+1,1));
else,
    vecske=ones(di+1,1);
end

for dd = 1:di,
    vecske(dd+1:end) = vecske(dd+1:end)*xi;
end
end

function v = vec( x )

% VEC   Vectorize.
%    VEC(X), where X is a vector, matrix, or N-D array, returns a column vector
%    containing all of the elements of X; i.e., VEC(X)=X(:). 
%
% Copyright 2005-2013 CVX Research, Inc.
% See the file COPYING.txt for full copyright information.
% The command 'cvx_where' will show where this file is located.

v = reshape( x, numel( x ), 1 );


end

function [weightedCPD, iter, relerr] = weightedcpd(tensor, oldLoad, weightTensorOrMatrix, options)
% My personal implentation of the CPD using alternating least squares.
% Written by Gabriel Hollander, Vrije Universiteit Brussel, Dept. ELEC, 21/04/2015

% general variables: threshold error and maximum number of iterations
if isempty(options), stepSizeTol_wTx = 1e-2; maxIter = 5000; wTxCommandWindowLog = 0; else stepSizeTol_wTx = options(1); maxIter = options(2); wTxCommandWindowLog = options(3); end

m = size(tensor, 2); n = size(tensor, 1); N = size(tensor, 3); iter = 1;
tensor_1 = reshape(tensor,  [n, m*N]);
tensor_2 = reshape(permute(tensor, [2, 1, 3]), [m, n*N]);
tensor_3 = reshape(permute(tensor, [3, 1, 2]), [N, m*n]);

% initialization for the CPD algorithm
W_cpd = oldLoad{1}; V_cpd = oldLoad{2}; H_cpd = oldLoad{3};
weightedCPD = {W_cpd, V_cpd, H_cpd};
relerr = zeros(maxIter, 1); relerr(iter) = norm(vec(cpdgen(weightedCPD) - tensor))/norm(vec(tensor));     % code from Philippe Dreesen
go_on = 1;

if isempty(weightTensorOrMatrix) % no weighting tensor is given, do unweighted CPD
    while and(iter < maxIter, go_on)
        if mod(iter, wTxCommandWindowLog) == 0, disp(['   wTx: no weight: iteration ' num2str(iter)]); end
        iter = iter + 1; W_cpd_old = W_cpd; H_cpd_old = H_cpd; V_cpd_old = V_cpd;
        % update of factor W
        Z = kr(H_cpd, V_cpd);                  % Khatri-Rao product
        ZtZ=Z'*Z; ZtX=Z'*(tensor_1)';          % code of Nway toolbox
        W_cpd = (pinv(ZtZ)*ZtX)';              % written on 21/04/2015
        % update of factor H
        Z = kr(V_cpd, W_cpd);                  % Khatri-Rao product
        ZtZ=Z'*Z; ZtX=Z'*(tensor_3)';          % code of Nway toolbox
        H_cpd = (pinv(ZtZ)*ZtX)';              % written on 21/04/2015
        % update of factor V
        Z = kr(H_cpd, W_cpd);                  % Khatri-Rao product
        ZtZ=Z'*Z; ZtX=Z'*(tensor_2)';          % code of Nway toolbox
        V_cpd = (pinv(ZtZ)*ZtX)';              % written on 21/04/2015
        weightedCPD = {W_cpd, V_cpd, H_cpd};
        relerr(iter) = norm(vec(cpdgen(weightedCPD) - tensor))/norm(vec(tensor));
%         step_size = sqrt(norm(W_cpd - W_cpd_old, 'fro') + norm(H_cpd - H_cpd_old, 'fro') + norm(V_cpd - V_cpd_old, 'fro')) / sqrt(norm(W_cpd, 'fro') + norm(H_cpd, 'fro') + norm(V_cpd, 'fro'));          % old step_size: problem with possible permutations of the factors, that's why we use the new step_size here under
        step_size = sum( cpderr( {W_cpd_old,V_cpd_old,H_cpd_old} , weightedCPD ));          % update of the step_size 2016/02/16
        if step_size < stepSizeTol_wTx, go_on = 0; end
    end
    relerr = relerr(1: iter);       % cut off the error vector after the loop has stopped
else    % use weighted CPD, using diagonal, block diagonal or full weighting matrix
    if size(weightTensorOrMatrix,3) > 1 % when a weighting tensor (3 dimensions) is given, use weighted CPD with diagonal or block diagonal weighting matrix
        % define the block diagonal weight matrix using cells of every slice
        Wcell = cell(1,N); for cell_iter = 1:N, Wcell{cell_iter} = weightTensorOrMatrix(:,:,cell_iter); end
        weightMatrix = blkdiag(Wcell{:});
%         weightMatrix
    else    % when a weighting matrix is given (so the third dimension = 1), then use a weighted CPD with full weighting matrix
        weightMatrix = weightTensorOrMatrix;
    end
%     disp('   wTx: sqrtm 1 starting...');
    weightMatrix_1 = reshapeWeightForMode(weightMatrix, 1, m, n, N); weightMatrix_1_sqrtm = sqrtm(weightMatrix_1);
%     disp('   wTx: sqrtm 2 starting...');
    weightMatrix_2 = reshapeWeightForMode(weightMatrix, 2, m, n, N); weightMatrix_2_sqrtm = sqrtm(weightMatrix_2);
%     disp('   wTx: sqrtm 3 starting...');
    weightMatrix_3 = reshapeWeightForMode(weightMatrix, 3, m, n, N); weightMatrix_3_sqrtm = sqrtm(weightMatrix_3);
    while and(iter < maxIter, go_on)
        if mod(iter, wTxCommandWindowLog) == 0, disp(['   wTx: weight: iteration ' num2str(iter)]); end
        iter = iter + 1; W_cpd_old = W_cpd; H_cpd_old = H_cpd; V_cpd_old = V_cpd;
        
        HV_blockmatrix = MakeBlockDiag(kr(H_cpd, V_cpd), n);
        W_cpd = reshape((weightMatrix_1_sqrtm*HV_blockmatrix) \ (weightMatrix_1_sqrtm*vec(tensor_1.')), size(W_cpd.')).';

        VW_blockmatrix = MakeBlockDiag(kr(V_cpd, W_cpd), N); VW_blockmatrix_sparse = sparse(VW_blockmatrix); weightMatrix_3_sqrtm_sparse = sparse(weightMatrix_3_sqrtm);
        H_cpd = reshape( (weightMatrix_3_sqrtm_sparse*VW_blockmatrix_sparse) \ (weightMatrix_3_sqrtm_sparse*vec(tensor_3.')), size(H_cpd.')).';

        HW_blockmatrix = MakeBlockDiag(kr(H_cpd, W_cpd), m); 
        V_cpd = reshape((weightMatrix_2_sqrtm*HW_blockmatrix) \ (weightMatrix_2_sqrtm*vec(tensor_2.')), size(V_cpd.')).';

        weightedCPD = {W_cpd, V_cpd, H_cpd};
        relerr(iter) = norm(vec(cpdgen(weightedCPD) - tensor))/norm(vec(tensor));
%         step_size = sqrt(norm(W_cpd - W_cpd_old, 'fro') + norm(H_cpd - H_cpd_old, 'fro') + norm(V_cpd - V_cpd_old, 'fro')) / sqrt(norm(W_cpd, 'fro') + norm(H_cpd, 'fro') + norm(V_cpd, 'fro'));         % old step_size: problem with possible permutations of the factors, that's why we use the new step_size here under
        step_size = sum( cpderr( {W_cpd_old,V_cpd_old,H_cpd_old} , weightedCPD ));          % update of the step_size 2016/02/16
        if step_size < stepSizeTol_wTx, go_on = 0; end
    end
    relerr = relerr(1:iter);
end
end

function [weightedCPD, iter, relerr, costfunction_relerr, exit_crit] = wTxCPD(tensor, oldLoad, covarianceMatrix, lambda, options)
% My second personal implentation of the CPD using alternating least squares. This implementaiton doesn't invert the full covariance matrix, but uses its SVD decomposition in order to achieve this goal.
% Written by Gabriel Hollander, Vrije Universiteit Brussel, Dept. ELEC, 28/05/2015

% general variables: threshold error and maximum number of iterations
if isempty(options), stepSizeTol_wTx = 1e-2; maxIter = 500; wTxCommandWindowLog = 0; else stepSizeTol_wTx = options(1); maxIter = options(2); wTxCommandWindowLog = options(3); end
maxCondit = 1e14;

% Defining the maintenance: dimensions, permutations, SVD 
n = size(tensor, 1); m = size(tensor, 2); N = size(tensor, 3); rankCovMatrix = rank(covarianceMatrix); iter = 1;
tensor_1 = reshape(tensor,  [n, m*N]);
tensor_2 = reshape(permute(tensor, [2, 1, 3]), [m, n*N]);
tensor_3 = reshape(permute(tensor, [3, 1, 2]), [N, m*n]);
[V,D,~] = svd(covarianceMatrix);                                          % because the covariance is symmetric, U = V, so we only take U and D
D1sqrtinv = diag(sqrt(1./diag(D(1:rankCovMatrix, 1:rankCovMatrix))));     % take only the non-zero singular values, invert them and take the sqrt
V1 = V(:, 1:rankCovMatrix);                                               % take only the first rankCovMatrix columns of V
V2 = V(:, (rankCovMatrix+1):end);                                         % take all the columns of V except the first rankCovMatrix 

T = reshape(1:m*n*N, [n,m,N]); I = eye(m*n*N);
reorder_1 = vec(reshape(T, [n, m*N]).');
P_1 = I(reorder_1, :);      % due to the code of Nway toolbox (and Kolda)
reorder_2 = vec(reshape(permute(T, [2, 1, 3]), [m, n*N]).');
P_2 = I(reorder_2, :);    % due to Kolda, so using a MATLAB-friendly matricizing order (following the three dimensions), and not a "mathematical"-friendly way (with cyclic permutations)
reorder_3 = vec(reshape(permute(T, [3, 1, 2]), [N, m*n]).');
P_3 = I(reorder_3, :);    % due to the code of Nway toolbox (and Kolda)

Q_1 = D1sqrtinv * V1.' * P_1.'; Q_2 = D1sqrtinv * V1.' * P_2.'; Q_3 = D1sqrtinv * V1.' * P_3.';

% initialization for the CPD algorithm
W_cpd = oldLoad{1}; V_cpd = oldLoad{2}; H_cpd = oldLoad{3}; norm_scalings = ones(size(V_cpd, 2), 1);
weightedCPD = {W_cpd, V_cpd, H_cpd};
relerr = NaN(maxIter, 1); relerr(iter) = norm(vec(cpdgen(weightedCPD) - tensor))/norm(vec(tensor));     % code from Philippe Dreesen
costfunction_relerr = NaN(maxIter, 1); costfunction_relerr(iter) = norm([D1sqrtinv * V1.' * vec(cpdgen(weightedCPD) - tensor); V2.' * vec(cpdgen(weightedCPD) - tensor)])/norm(vec(tensor));
go_on = 1;

while and(iter < maxIter, go_on)
    if mod(iter, wTxCommandWindowLog) == 0, disp(['   wTx: iteration ' num2str(iter)]); end
    iter = iter + 1;
    W_cpd_old = W_cpd; H_cpd_old = H_cpd; V_cpd_old = V_cpd;

%   Update of the factor W
    HV_blockmatrix = MakeBlockDiag(kr(H_cpd, V_cpd), n);
    W_cpd = reshape([(Q_1*HV_blockmatrix); lambda * V2.' * P_1.' * HV_blockmatrix] \ [Q_1*vec(tensor_1.'); lambda * V2.' * P_1.' * vec(tensor_1.')], size(W_cpd.')).';            % last update, with the weight factor lamda, 25/03/2016
    for witer = 1 : size(W_cpd, 2), 
        W_cpd(:,witer) = W_cpd(:,witer) / norm(W_cpd(:,witer)); 
    end

%   Update of the factor H
    VW_blockmatrix = MakeBlockDiag(kr(V_cpd, W_cpd), N); VW_blockmatrix_sparse = sparse(VW_blockmatrix);
    H_cpd = reshape([Q_3*VW_blockmatrix_sparse; lambda * V2.' * P_3.' * VW_blockmatrix_sparse] \ [Q_3*vec(tensor_3.'); lambda * V2.' * P_3.' * vec(tensor_3.')], size(H_cpd.')).';            % last update, with the weight factor lamda, 25/03/2016
    for hiter = 1 : size(H_cpd, 2), 
        H_cpd(:,hiter) = H_cpd(:,hiter) / norm(H_cpd(:,hiter)); 
    end    
    
%     equations_H = [Q_3*VW_blockmatrix_sparse; V2.' * VW_blockmatrix_sparse];

%   Update of the factor V
    HW_blockmatrix = MakeBlockDiag(kr(H_cpd, W_cpd), m); 
    V_cpd = reshape([Q_2 * HW_blockmatrix; lambda * V2.' * P_2.' * HW_blockmatrix] \ [Q_2 * vec(tensor_2.'); lambda * V2.' * P_2.' * vec(tensor_2.')], size(V_cpd.')).';            % last update, with the weight factor lamda, 25/03/2016

    weightedCPD = {W_cpd, V_cpd, H_cpd};
    relerr(iter) = norm(vec(cpdgen(weightedCPD) - tensor))/norm(vec(tensor));
    
    costfunction_relerr(iter) = norm([D1sqrtinv * V1.' * vec(cpdgen(weightedCPD) - tensor); V2.' * vec(cpdgen(weightedCPD) - tensor)])/norm(vec(tensor));
%     step_size = sqrt(norm(W_cpd - W_cpd_old, 'fro') + norm(H_cpd - H_cpd_old, 'fro') + norm(V_cpd - V_cpd_old*diag(norm_scalings), 'fro')) / sqrt(norm(W_cpd, 'fro') + norm(H_cpd, 'fro') + norm(V_cpd,'fro'));                   % old step_size: problem with possible permutations of the factors, that's why we use the new step_size here under
    step_size = sum( cpderr( {W_cpd_old,V_cpd_old*diag(norm_scalings),H_cpd_old} , weightedCPD ));          % update of the step_size 2016/02/16
    if step_size < stepSizeTol_wTx, go_on = 0; exit_crit = 1; end
    if any(vec(isnan([W_cpd; H_cpd; V_cpd])))   % if there is a serious ill-conditioned matrix W_cpd, H_cpd and/or V_cpd, then return one step in the iteration
        go_on = 0; exit_crit = 2; iter = iter - 1; W_cpd = W_cpd_old; H_cpd = H_cpd_old; V_cpd = V_cpd_old*diag(norm_scalings); weightedCPD = {W_cpd, V_cpd, H_cpd};
    end
    
%   If we don't stop the algorithm, then rescale the columns of the factor V (as we did in any case for the factors W and H)
    if and( iter < (maxIter - 1), go_on)
        norm_scalings = NaN(size(V_cpd, 2), 1);
        for viter = 1 : size(V_cpd, 2), 
            norm_scalings(viter) = norm(V_cpd(:,viter));
            V_cpd(:,viter) = V_cpd(:,viter) / norm_scalings(viter); 
        end
    end

end
if go_on == 1, exit_crit = 3; end;
relerr = relerr(1:(iter));
costfunction_relerr = costfunction_relerr(1:(iter));
% plot(costfunction_relerr(50:end))             % plot the cost function at the end of the iterations

% Make a normalization on the V factor, and accomodate with the two other factors W and H
for viter = 1 : size(V_cpd,2)
    scaledFactor = norm(V_cpd(:,viter));
    V_cpd(:,viter) = V_cpd(:,viter) / scaledFactor;
    W_cpd(:,viter) = W_cpd(:,viter) * sqrt(scaledFactor);
    H_cpd(:,viter) = H_cpd(:,viter) * sqrt(scaledFactor);
end
weightedCPD = {W_cpd, V_cpd, H_cpd};

end
