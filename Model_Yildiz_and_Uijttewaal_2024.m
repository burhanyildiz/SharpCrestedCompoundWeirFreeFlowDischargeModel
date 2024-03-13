%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Sharp-Crested Rectangular Compound Weir Flow Rate Calculation 
% Yildiz and Uijttewaal (2024) DOI:
% All lenghts are in meters, flow rates in m^3/s
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Clear weir plate heights at each notch (Modify it according to your weir configuration)
config = [0.2 0.2 0.15 0.1 0.05 0 0 0.2 0.2 0.2 0.2 0.2];
% Threshold amounts at the weir footings, if any (Modify it according to your setup)
thresholds = [0.00358 0.00207 0.00288 0.0029 0.00228 0.0031 0.00272 0.00284 0.0023 0.00186 0.00176 0.0021];
% Upstream flow depth measurements at 6 cases (Modify it according to your data)
hu = [0.1662625; 0.152875; 0.141455; 0.1259625; 0.108725; 0.08695];
% Weir width (Modify it according to your setup)
b=0.23;
% Full width including weir and buttress (bful=b when no buttresses are
% used) (Modify it according to your setup)
bful=0.25;

% Calculation lines are given below. Only modify lines 30 and 31 if the stated
% condition given at line 28 is valid
% Number of weir openings, N
N=size(config,2);
% Number of cases, M
M=size(hu,1);
% Total weir height
P = config + thresholds;

% Eqn.3 RHS coefficient calculations. If (bful-b)/bful)> 0.1 arrange these
% two lines according Kindsvater&Carter method
S=0.602-(10*(bful-b)/bful)*(0.602-0.599);
T=0.075-(10*(bful-b)/bful)*(0.075-0.064);

% Kindsvater&Carter method correction parameters against viscosity and
% surface tension
Kb=0.015;
Kh=0.001;
% Initial zero matrixes for related parameters
Cd=zeros(M,N);
Q1=zeros(M,N);
stDev=zeros(M,1);
sumQ=zeros(M,1);
Qd=zeros(M,1);
CL=zeros(M,10);
Q1_update=zeros(M,N);
% Main calculation loops for discharge coefficients, correction amounts and flow
% rates
for i=1:M
    for j=1:N
        if P(j)<0.01 % assigns empty weir coefficient
            Cd(i,j)=0.66;
        else
            Cd(i,j)=S+T*(hu(i)-P(j))/P(j); % Eqn. (3)
        end
        if hu(i)<P(j) % non-flowing notch condition
            Q1(i,j)=0;
        else
            Q1(i,j)= (2/3)*(b+Kb)*sqrt(19.62)*Cd(i,j)*((hu(i)+Kh-P(j))^1.5); % Eqn. (1)
        end
    end
    Q_model_uncorrected=sum(Q1,2); % Flow rate over the weir when uncorrected model was used (For comparison)
    % Before entering the CL iteration, first values of CL are estimated
    % using the uncorrected flow rates at each notch
    stDev(i)=std(Q1(i,:),1); % Standard deviation of the flow rate values at each notch
    sumQ(i)=sum(Q1(i,:)); % Total flow rate over compound weir (Uncorrected yet)
    Qd(i)= stDev(i)/(N*b*(9.81^0.5)*(hu(i)^1.5)); % Qd* calculation Eqn. (11)
    k=2; % CL iteration matrix column starting number
    CL(i,k)=1.01-16.22*Qd(i); % Eqn 13 for the first element in CL convergence matrix
    tol=0.0001; % Allowed tolerance in CL convergence
    % CL iteration is applied in the while loop given below
    while abs(CL(i,k)-CL(i,k-1))>tol
        Q1_update(i,:)=CL(i,k).*Q1(i,:); % Intermediate parameter for flow rates at each i
        stDev(i)=std(Q1_update(i,:),1); % standard deviation of the flow rates at each notch
        sumQ(i)=sum(Q1_update(i,:)); % Total flow rate over compound weir (intermediate step)
        Qd(i)= stDev(i)/(N*b*(9.81^0.5)*(hu(i)^1.5)); % Qd* calculation Eqn. (11)
        CL(i,k+1)=1.01-16.22*Qd(i); % Eqn. (13)
        k=k+1;
    end
end
[row_indices, col_indices] = find(CL ~= 0); % Estimating the last non-zero element of CL matrix
last_non_zero_indices = accumarray(row_indices, col_indices, [], @max);
CL_end = zeros(size(CL, 1), 1); % Converged CL matrix empty form
for row = 1:size(CL, 1)
    linear_index = sub2ind(size(CL), row, last_non_zero_indices(row));
    CL_end(row) = CL(linear_index); % Converged CL values at each case
    if CL_end(row)<0.9 % Applying correction to the model
        Q1_update(row,:)=CL_end(row)*Q1(row,:);
    else
        Q1_update(row,:)=Q1(row,:);
    end
end
Q_model=sum(Q1_update,2); % Flow rate over the weir (Result)
% The measured flow rates at each case for comparison
Q_measured_C8 = [0.07036; 0.0592; 0.05044; 0.03973; 0.02956; 0.01983];