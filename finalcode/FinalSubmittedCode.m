clear;clc;close all;clear all;

num_assets =2;num_periods = 50;CL = 0.9;alph=1-CL;covScaleFactor = 3;

mu_set = randi([1 350], [num_assets,1])
Sigma = generateSPDmatrix(num_assets,covScaleFactor);
%Sigma = diag(diag(Sigma))  %un-comment to imply independent exp. returns
S = Sigma*(num_periods-1);

mu_i = mvnrnd(mu_set, Sigma, num_periods)';

%% calculation of sample values and ordering of the eigen decomposition

mu_bar = mean(mu_i')';              %p x 1 column vector of sample means
estS =(mu_i-mu_bar)*(mu_i-mu_bar)'; % p x p  scatter matrix of mu_i
estSigma = estS/(num_periods-1);    %p x p  sample covariance matrix

[T, Lambda] = eig(estSigma);        %T is p x p eigenvector matrix
                                    %Lambda is p x p diagonal matrix of 
                                    %eigen val.
 
%sort eigenvectors according to descending eigenvalues                                    
[misc, perm] = sort(diag(Lambda), 'descend');
Lambda = Lambda(perm, perm);T=T(:,perm); 

%% Computation of critical values and confidence regions
if(num_assets==1)
    muSEM = sqrt(estSigma/num_periods);
    ts = tinv([alph/2 1-alph/2] , num_periods-1)
    muDeviation = (ts)*muSEM;
    muCI = mu_bar + muDeviation
    
    chisq = chi2inv([alph/2 1-alph/2], num_periods-1)
    varSigmaCoeff = [num_periods-1 num_periods-1]./chisq;
    varSigmaCI = flip(varSigmaCoeff)*estSigma;
else
    critF = num_assets*(num_periods-1)/(num_periods-num_assets)*finv(CL,num_periods-1,num_assets) 
end

%% Plotting the confidence regions 



if(num_assets == 3) 
figure('Position',[400 400 450 450])            
hold on
title(['3D Confidence Area for \mu, CL = ' num2str(CL)]  )
q1=quiver3(mu_bar(1),mu_bar(2),mu_bar(3),T(1,1)*sqrt(Lambda(1,1)),T(1,2)*sqrt(Lambda(1,1)),T(1,3)*sqrt(Lambda(1,1)),'LineWidth',3);
q2=quiver3(mu_bar(1),mu_bar(2),mu_bar(3),T(2,1)*sqrt(Lambda(2,2)),T(2,2)*sqrt(Lambda(2,2)),T(2,3)*sqrt(Lambda(2,2)),'LineWidth',3);
q3=quiver3(mu_bar(1),mu_bar(2),mu_bar(3),T(3,1)*sqrt(Lambda(3,3)),T(3,2)*sqrt(Lambda(3,3)),T(3,3)*sqrt(Lambda(3,3)),'LineWidth',3);
plot3(mu_i(1,:),mu_i(2,:),mu_i(3,:), '.') ;
xlabel('\mu_1');
ylabel('\mu_2');
zlabel('\mu_3');
H3 = plot_gaussian_ellipsoid(mu_bar,estSigma,critF/num_periods);
set(H3,'FaceAlpha',0.075, 'edgecolor', 'k', 'facecolor', [.5 .9 .1]);
grid on;
ax = gca;
set(q1, 'ShowArrowHead','off')
set(q2, 'ShowArrowHead','off')
set(q3, 'ShowArrowHead','off')
ax.Clipping = 'off' 
hold off            


elseif(num_assets==2)
figure('Position',[400 400 400 350]) 
hold on
title(['2D Confidence Area for \mu, CL = ' num2str(CL)]  )
critF = num_assets*(num_periods-1)/(num_periods-num_assets)*finv(CL,num_periods-1,num_assets)
q1 = quiver(mu_bar(1),mu_bar(2),T(1,1)*sqrt(Lambda(1,1)),T(1,2)*sqrt(Lambda(1,1)),'LineWidth',3)
q2 = quiver(mu_bar(1),mu_bar(2),T(2,1)*sqrt(Lambda(2,2)),T(2,2)*sqrt(Lambda(2,2)),'LineWidth',3)
data1 = plot(mu_i(1,:),mu_i(2,:),'.');
xlim([min(min(mu_i))-10 max(max((mu_i)))+10])                  
ylim([min(min(mu_i))-10 max(max((mu_i)))+10])    
xlabel('\mu_1');
ylabel('\mu_2');
H2 = plot_gaussian_ellipsoid(mu_bar,estSigma,critF/num_periods);  
alpha(.5)
set(q1, 'ShowArrowHead','off')
set(q2, 'ShowArrowHead','off')
hold off
grid on



%%
elseif(num_assets==1)
figure('Position',[400 400 400 350]) 
hold on
title(['Confidence Interval for \mu, CL = ' num2str(CL)]  )

p1 = [num_periods/2 mu_bar ];
p2 = [num_periods/2 mu_bar+muDeviation(1)];
p3 = [num_periods/2 mu_bar+muDeviation(2)];
dp = p2-p1;
dp2= p3-p1;

q1 = quiver(p1(1),p1(2),dp(1),dp(2),0,'LineWidth',3);
q2 = quiver(p1(1),p1(2),dp2(1),dp2(2),0,'LineWidth',3);
%data1 = plot(mu_i(1,:),'.');
xlim([0 num_periods]);
ylim([p2(2) p3(2)]);
xlabel('Period');
ylabel('\mu_i values');
set(q1, 'ShowArrowHead','off')
set(q2, 'ShowArrowHead','off')
axis 'auto x'
axis 'auto y'
hold off

figure('Position',[400 400 400 350]) 
hold on
title(['Confidence Interval for \sigma^2, CL = ' num2str(CL)])
estvarSigma = 2*(estSigma)^(2)/(num_periods-1)

ps1 = [num_periods/2 estSigma];
ps2 = [num_periods/2 varSigmaCI(1)];
ps3 = [num_periods/2 varSigmaCI(2)];
dps = ps2-ps1;
dps2= ps3-ps1;
q3 = quiver(ps1(1),ps1(2),dps(1),dps(2),0,'LineWidth',3);
q4 = quiver(ps1(1),ps1(2),dps2(1),dps2(2),0,'LineWidth',3);
xlim([0 num_periods]);
ylim([ps2(2) ps3(2)]);
ylabel('\sigma^2 values');
ax2Dsig = gca;
set(ax2Dsig, 'XTickLabel',[])
axis 'auto y'
hold off
grid on
set(q3, 'ShowArrowHead','off')
set(q4, 'ShowArrowHead','off')
ax = gca;
ax.Clipping = 'off';



end



return    
%% ignore below this line

%% Generating a random sample of asset returns given the sampling dist. of
%  the covariance
%  this section of the code is to demonstrate the use of the sampling dist.
%  of the sample covariance (Wishart dist.) to generate similar expected
%  returns vectors as was done by the mvnrnd function at the beginning of
%  the script


%  first the scatter matrices are generated

%genS1 = wishrnd(estS/(num_periods-1),num_periods-1,chol(estS/(num_periods-1)));
%genS2 = wishrnd(estS/(num_periods-1),num_periods-1,chol(estS/(num_periods-1)));
%genS3 = wishrnd(estS/(num_periods-1),num_periods-1,chol(estS/(num_periods-1)));

%  then the covariance matrices are calculated via S/n-1

%genSigma1=genS1/(num_periods-1);
%genSigma2=genS2/(num_periods-1);
%genSigma3=genS3/(num_periods-1);
%Z = randn(num_assets,num_periods);
%X1 = chol(genSigma1);
%X2 = chol(genSigma2);        %cholesky decomp. was used to confirm PD 
%X3 = chol(genSigma3);        %nature of generated covariance matrices

%gen_mu1 = X1*Z + mu_bar;
%gen_mu2 = X2*Z + mu_bar;     %defining new datasets with covariance mat.
%gen_mu3 = X3*Z + mu_bar;     %generated from the Wishart dist.
%norm(Sigma-estSigma);
%norm(S-estS)/(num_periods-1);
%norm(S-genS1)/(num_periods-1);
%norm(S-genS2)/(num_periods-1);     %verifying the magnitude of difference
%norm(S-genS3)/(num_periods-1);     %between the cov. matrices

%norm(gen_mu1-mu_i);
%norm(gen_mu2-mu_i);      %verifying magnitude of difference between the
%norm(gen_mu3-mu_i);      %generated expected returns

