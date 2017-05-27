%---------------------------------------------------------------------%
%Q1
data = csvread('Q1.csv');
data = sortrows(data, 1);
%%a.
monthlysimplenetreturn = data(2:end,2) ./ data(1:end-1,2) - 1;
figure()
plot(monthlysimplenetreturn)
xlabel('time');
ylabel('return');
title('monthly simple net return')

%%b.
monthlygrossreturn = data(2:end,2) ./ data(1:end-1,2);
plot(monthlygrossreturn)
xlabel('time');
ylabel('return');
title('monthly gross return')
%%c.
continuouslycompoundedreturn = log(data(2:end,2)) - log(data(1:end-1,2)) - 1;
plot(continuouslycompoundedreturn)
xlabel('time');
ylabel('return');
title('continuously compounded return')
%%d.
meanret = mean(monthlysimplenetreturn);
medianret = median(monthlysimplenetreturn);
stdret = std(monthlysimplenetreturn);
fivepercentileret = prctile(monthlysimplenetreturn, 5);
ninetyfivepercentileret = prctile(monthlysimplenetreturn, 95);
skewnessret = skewness(monthlysimplenetreturn);
kurtosisret = kurtosis(monthlysimplenetreturn);
statistictable = [meanret;medianret;stdret;fivepercentileret;ninetyfivepercentileret; skewnessret; kurtosisret];
csvwrite('Q1table.csv', statistictable);
%---------------------------------------------------------------------%

%---------------------------------------------------------------------%
%Q2
data = csvread('Q2.csv');
retdata = data(2:end, 2:end) ./ data(1:end-1, 2:end) - 1;
meanret = mean(retdata);
stdret = std(retdata);
statistictable = [meanret', stdret'];
csvwrite('Q2table.csv', statistictable);
%---------------------------------------------------------------------%

%---------------------------------------------------------------------%
%Q3
Rm = csvread('Q3.csv');
Rtdata = csvread('Q2.csv');
Rtdata = Rtdata(2:end, 2:end) ./ Rtdata(1:end-1, 2:end) - 1;
Rm = Rm(2:end, 2) ./ Rm(1:end-1, 2) - 1;
Rf = 0.035 / 12 * ones(size(Rtdata, 1), 1);
Rme = Rm - Rf;
for i = 1 : size(Rtdata, 2)
    y = Rtdata(:,i) - Rf;
    x = [ones(size(Rte,1), 1), Rme];
    [B,BINT,R,RINT,STATS] = regress(y, x);
    Intercept(i) = B(1);
    Beta(i) = B(2);
    Residual(i,:) = R;
end
statistictable = [Intercept', Beta'];
csvwrite('Q3table.csv', statistictable);
%---------------------------------------------------------------------%

%---------------------------------------------------------------------%
%Q4
a = Intercept';
sigma = Residual * Residual' / T;
T = size(Rtdata, 1);
J0 = T * (1 + mean(Rm)^2 / var(Rm))^(-1) * a' * sigma^(-1) * a;
%---------------------------------------------------------------------%

%---------------------------------------------------------------------%
%Q5
S = (mean(Rm) - Rf(1)) / std(Rm);
N = size(Rtdata, 2);
J1 = (T - N - 1) / N * a' * sigma^(-1) * a / (1 + S^2);
%---------------------------------------------------------------------%
%---------------------------------------------------------------------%
%Q6. Momemtum
data = csvread('Q6_2.csv');
data(:,3) = data(:,3) + 1;
timelist = data(:,4) * 100 + data(:,5);
timelist = unique(timelist);
timerank = [];
for i = 1 : length(data)
    timerank(i) = find(timelist == data(i,4) * 100 + data(i,5));
end
data(:,7) = timerank;
stocklist = unique(data(:,1));
totaldata = [];
recpos = 0;
data(:,8) = zeros(length(data),1);

for i = 1 : length(stocklist);
    i
    slicedata = data(data(:,1) == stocklist(i), :);
    if size(slicedata,1) < 12
        recpos = recpos + size(slicedata,1);
        continue
    end
    for j = 1 : size(slicedata,1) - 12
        ret = cumprod(slicedata(j:j+10, 3));
        ret = ret(end);
        data(j + 12 + recpos, 8) = ret;
    end
    recpos = recpos + size(slicedata,1);
end

retdata = data(data(:,8) ~= 0, :);
datelist = unique(retdata(:, 2));

timerank = unique(timerank);
groupret = [];
for i = 1 : length(timerank)
    slicedata = retdata(retdata(:, 7) == timerank(i), :);
    slicedata = sortrows(slicedata, 8);
    idx = 1;
    retlist = [];
    for j = 1 : 9
        monthlyret = slicedata(idx : idx + floor(size(slicedata,1) / 10) - 1, 3);
        retlist = [retlist, mean(monthlyret)];
        idx = idx + floor(size(slicedata,1) / 10);
    end
    monthlyret = slicedata(idx : end, 3);
    retlist = [retlist, mean(monthlyret)];
    groupret = [groupret; retlist];
end
groupret = groupret(~isnan(groupret(:,1)), :);
cumret = cumprod(groupret);
figure()
plot(cumret);
legend('1','2','3','4','5','6','7','8','9','10', 'Location','NorthWest');
xlabel('time')
ylabel('return')
title('Returns For Momentum Groups')

%Fama-French Test
threefactor = csvread('threefactor.csv', 1, 0);
factoridx = ismember(threefactor(:,1), datelist);
threefactor = threefactor(factoridx, :);
rm = threefactor(:,2);
smb = threefactor(:,3);
hml = threefactor(:,4);
x = [ones(size(groupret(:, 1), 1),1), rm, smb,hml];
fftable = [];
for i = 1 : 10
    y = groupret(:, i) - 1;
    [B,BINT,R,RINT,STATS] = regress(y, x);
    SSE = R' * R;
    MSE = SSE / (size(x,1) - size(x,2));
    se = sqrt(MSE*diag(inv(x'*x)));
    t = B ./ se;
    fftable(:,i) = [B; t; STATS(1)]
end
csvwrite('momemtumfftest.csv', fftable);
%---------------------------------------------------------------------%
%---------------------------------------------------------------------%
%Q7.Reversal
data = csvread('Q6_2.csv');
data(:,3) = data(:,3) + 1;
timelist = data(:,4) * 100 + data(:,5);
timelist = unique(timelist);
timerank = [];
for i = 1 : length(data)
    timerank(i) = find(timelist == data(i,4) * 100 + data(i,5));
end
data(:,7) = timerank;
stocklist = unique(data(:,1));
totaldata = [];
recpos = 0;
data(:,8) = zeros(length(data),1);

for i = 1 : length(stocklist);
    i
    slicedata = data(data(:,1) == stocklist(i), :);
    if size(slicedata,1) < 48
        recpos = recpos + size(slicedata,1);
        continue
    end
    for j = 1 : size(slicedata,1) - 48
        ret = cumprod(slicedata(j:j+35, 3));
        ret = ret(end);
        data(j + 48 + recpos, 8) = ret;
    end
    recpos = recpos + size(slicedata,1);
end

retdata = data(data(:,8) ~= 0, :);
datelist = unique(retdata(:, 2));

timerank = unique(timerank);
groupret = [];
for i = 1 : length(timerank)
    slicedata = retdata(retdata(:, 7) == timerank(i), :);
    slicedata = sortrows(slicedata, 8);
    idx = 1;
    retlist = [];
    for j = 1 : 9
        monthlyret = slicedata(idx : idx + floor(size(slicedata,1) / 10) - 1, 3);
        retlist = [retlist, mean(monthlyret)];
        idx = idx + floor(size(slicedata,1) / 10);
    end
    monthlyret = slicedata(idx : end, 3);
    retlist = [retlist, mean(monthlyret)];
    groupret = [groupret; retlist];
end
groupret = groupret(~isnan(groupret(:,1)), :);
cumret = cumprod(groupret);
figure()
plot(cumret);
legend('1','2','3','4','5','6','7','8','9','10', 'Location','NorthWest');
xlabel('time')
ylabel('return')
title('Returns For Reversal Groups')

%Fama-French Test
threefactor = csvread('threefactor.csv', 1, 0);
factoridx = ismember(threefactor(:,1), datelist);
threefactor = threefactor(factoridx, :);
rm = threefactor(:,2);
smb = threefactor(:,3);
hml = threefactor(:,4);
x = [ones(size(groupret(:, 1), 1),1), rm, smb,hml];
fftable = [];
for i = 1 : 10
    y = groupret(:, i) - 1;
    [B,BINT,R,RINT,STATS] = regress(y, x);
    SSE = R' * R;
    MSE = SSE / (size(x,1) - size(x,2));
    se = sqrt(MSE*diag(inv(x'*x)));
    t = B ./ se;
    fftable(:,i) = [B; t; STATS(1)]
end
csvwrite('reversalfftest.csv', fftable);
%---------------------------------------------------------------------%
%---------------------------------------------------------------------%
%Q8
r104=zeros(length(EAPS4(:,1)),length(EAPS4(1,:)));
for i=1:length(EAPS4(:,1))
    for j=2:length(EAPS4(:,1))
        r104(j,i)=EAPS4(j,i)/EAPS4(j-1,i)-1;
    end
end
 omg = cov(r104);
 coeff = pcacov(omg);
omg = cov(r104);
 coeff = pcacov(omg);
N=10;
compvar=zeros(1,N);
for i = 1:N
    compvar(i)=coeff(:,i)'*omg*coeff(:,i);
end
totvar = sum(compvar);
prctvar=compvar/totvar;
plot(prctvar)
