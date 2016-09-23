cwd; clear; clc;
% Initialize:
rng('shuffle');
m = 2; n = 2; r = 2; degreeMax = 3;
Nsample = 1000; NOperatingPoints = 50; NoptimizedValidation = 1000; Ntransients = 2;
V = randn(m, r); for columns = 1 : r, V(:, columns) = V(:, columns)/norm(V(:, columns)); end
signs = [1 -1]; 
G = [RandomReal([0.05 0.2], [r 1]) RandomReal([0.1 0.4], [r 1]) RandomReal([0.5 1], [r 1]) zeros(r,1)].*signs(RandomInteger([1 length(signs)], [r, 4]));
W = randn(n, r); for rows = 1 : n, W(rows, :) = W(rows, :)/norm(W(rows, :)); end
[bInput, aInput] = CreateRandomChebyFilters([2 4],[0 20],[0.05 0.2 ],m); [bOutput, aOutput] = CreateRandomButterFilters([1 4],[0.1 0.3],n);
validationScales = [1 0.8 0.5];
operatingPointsSample = sort(randsample(Nsample, NOperatingPoints).');
lambda_weightForSecondSetOfEqs = 1;
uMultisine = RandomPhaseMultisine(Nsample, 1/8); 
uMultisineOpt = RandomPhaseMultisine(Nsample, 1/8);
for scale_iter = 1 : length(validationScales), uMultisineVal(scale_iter,:) = RandomPhaseMultisine(NoptimizedValidation, 1/8); end
added_noise_sigma_y = 0.25;
lambda_weightForSecondSetOfEqs = 1;
options{1} = 0; %  CPD initialization: 0  = DTLD/GRAM for initialization,  2  = random orthogonalized values for initialization
options{2} = NaN; % log during CPD in command window: NaN = almost none
options{3} = randOldLoad(m, n, r, NOperatingPoints, 2);             % randomize the initizal value for the CPD
options{4} = 10; % max number of functions calls in optimization
options{5} = 'off'; % choose 'off' or 'iter' for logging during optimization
options{6} = 1e-6;      % max relative error for the wTx CPD
options{7} = 2000;       % max number of iterations for the wTx CPD
options{8} = 1000;        % show log in command window at every N-th iteration of the wTX CPD
options{9} = 8000;      % max number of functions calls in optimization with wTx

added_noise_on_output = randn(1, Nsample); 
added_noise_on_outputOptim = randn(1, Nsample);
for scale_iter = 1 : length(validationScales), added_noise_on_outputVal(scale_iter,:) = randn(1, Nsample); end

generalVariables         = makeGeneralVariables(Nsample, NOperatingPoints, NoptimizedValidation, Ntransients, operatingPointsSample, uMultisine, uMultisineOpt, uMultisineVal, V, G, W, bInput, aInput, bOutput, aOutput, added_noise_sigma_y, added_noise_on_output, added_noise_on_outputOptim, added_noise_on_outputVal, validationScales, lambda_weightForSecondSetOfEqs, options);
approximatedPoly         = approximatePolynomial(generalVariables);
truePoly                 = truePolynomial(generalVariables);

% with weighted toolbox (wTx)
% cpdTrue_wTx              = CPD_weighted_True(generalVariables, approximatedPoly, truePoly);

uOperatingPoints = approximatedPoly.uOperatingPoints;
monomials = approximatedPoly.monomials; varSym = approximatedPoly.varSym;
modelInfo = generalVariables.modelInfo;
approximatedCoeffs = approximatedPoly.approximatedCoeffs;
% covarianceMatrix = approximatedPoly.covarianceMatrix;
covarianceMatrix = CreatePositiveDefiniteMatrix(20);
r = generalVariables.r; 
degreeMax = generalVariables.degreeMax;
oldLoad = generalVariables.oldLoad;
lambda_weightForSecondSetOfEqs = generalVariables.lambda_weightForSecondSetOfEqs;
stepSizeTol_wTx = generalVariables.stepSizeTol_wTx;
maxNwTx = generalVariables.maxNwTx;
wTxCommandWindowLog = generalVariables.wTxCommandWindowLog;




coupledPolynomial.approximatedCoeffs = [0 0; approximatedCoeffs];
coupledPolynomial.covarianceMatrix = covarianceMatrix; 
coupledPolynomial.CPDtype = 'blockdiag';
coupledPolynomial.r = 2;

