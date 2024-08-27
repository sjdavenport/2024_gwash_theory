X = Xgen( 1000, 2000, 0.95, '012');
mean(std(X,0,1))

corrval = 0;
for I = 1:size(X,2)
    corrval = corrval + corr(X(:,1), X(:,2));
end
corrval = corrval/size(X,2)

%%
rho = 0.99;
X = Xgen( 1000, 2000, rho, 'equi');
Y = Xgen( 1000, 2000, rho, 'equi');

Z = (X > 0) + (Y>0);

corrval = 0;
lag = 1;
for I = 1:(size(Z,2)-lag)
    corrval = corrval + corr(Z(:,1), Z(:,1+lag));
end
corrval = corrval/(size(Z,2) -lag)

std(Z,0,1)/(1/sqrt(2))
mean(Z,1)