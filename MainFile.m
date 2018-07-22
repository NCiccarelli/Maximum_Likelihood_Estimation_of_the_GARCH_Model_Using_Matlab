 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% We report Matlab code for Maximum Likelihood estimation of the GARCH model; moreover, we report a Monte Carlo simulation which shows that the Maximum Likelihood estimator converges to the true parameters. We use the t5-student innovation for the GARCH process.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% number of repetitions for the Monte Carlo simulation:
repetitions = 2000; 

% initialize matrices for coefficients (ML estimator):
MLE_alpha= NaN(1,repetitions);
MLE_beta= NaN(1,repetitions);



for w = 1: repetitions

    
%% sample size for the GARCH process is set to 1000:    
N=1000;


%% generate the t5-student innovation for the GARCH process 
t5_distr = trnd(5,[N,1])';
% generate t5-student innovation with mean zero and variance 1: 
epsilon =  (t5_distr  - mean(t5_distr) )./ (sqrt(var(t5_distr)  ));


%% parameter values (GARCH volatility equation): 
beta0= 0.8;
alfa0=0.1;


%% generate the data: 
h=zeros(1,N);
y=zeros(1,N);
h(1)= 1./ (1-alfa0-beta0);
y(1)= ((h(1)).^0.5) * (epsilon(1)); 

for  i=1: (N-1);
    h(i+1)= 1+ beta0 .* h(i) + alfa0 .* ((y(i)).^2);
    y(i+1)=  ( (h(i+1)).^0.5) * (epsilon(i+1));
end


%% starting values for the estimation of parameters (maximum likelihood estimation): 
startingvalues = [0.0001;0.0001; 0.0001];
lowerbound = [  0 ;  0  ; 0];
upperbound = [  1;  1; 10];
  
%% estimate the parameters using MLE (maximum likelihood estimation)
% see the MLE_t5_NEW function for more details:
likelihood = @(x) MLE_t5_NEW(x(1),x(2),x(3), y,h );

options    = optimset('fmincon');
options    = optimset(options, 'TolFun', 1e-006);
options    = optimset(options, 'LargeScale', 'off');
options    = optimset(options, 'MaxFunEvals', 1000);
options    = optimset(options, 'MaxIter', 400);
options = optimset('maxfunevals',20000);

[PARMLE, fval] = fmincon(likelihood,startingvalues,[],[],[],[],lowerbound ,upperbound , @mycon ,options);

MLE_alpha(1,w)= PARMLE(1,1);
MLE_beta(1,w)= PARMLE(2,1);
MLE_sigma(1,w)= PARMLE(3,1);





w
end

% means and standard deviation of the GARCH parameters using the ML
% estimator: 
Mean_alpha_MLE= mean(MLE_alpha,2)
Mean_beta_MLE= mean(MLE_beta,2)
Std_dev_alpha_MLE =std(MLE_alpha)
Std_dev_beta_MLE= std(MLE_beta) 
Mean_sigma_MLE= mean(MLE_sigma,2)




