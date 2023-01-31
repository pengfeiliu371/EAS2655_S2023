% Single line comment

% Mutiple line comment

%{
We plot the time series of a synthetic data generated and
the histogram to see the distribution of monthly temperatures 
for the month of January for the Disneyland
%}

%% Section 1
% importing data (avg. monthly temperatures from 1970 to 2020 for all months)
yrs= 1970:2019;
dat= 66+ 2 .*rand(50,12);

%% Section 2
% Plotting avg. monthly temperatures across years (1970- 2020)

figure(1);
plot(yrs, dat(:,1));
ylabel("Temperature (deg F)");
xlabel("years")
title("Jan monthly avg. temp (1970-2020)")

%% Section 3
%{
Plotting histogram for the month of January to view the distribution of temp.
%}

figure(2);
hist(dat(:,1));
ylabel("count")
xlabel("Temperature (deg F)");
title("Histogram for Jan monthly temp (1970-2020)");
%%
%{
Publish: 
https://www.mathworks.com/help/matlab/matlab_prog/publishing-matlab-code.html

Unusal blank spaces: 
https://www.mathworks.com/matlabcentral/answers/880303-how-can-i-remove-blank-spaces-in-generated-pdf-using-publish-command
%}