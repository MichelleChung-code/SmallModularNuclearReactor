% Example of how to do a cashflow diagram

% Cashflows and Dates should be read in from a CSV in the future when we
% get valid numvers

clc

cashflow_dates = [20200110, 20210502, 20231001];
cashflow_dates = datenum(num2str(cashflow_dates'),'yyyymmdd')';

cashflow_amts = [-10, 100, 200]; 

figure;
show_pts = 0:1:length(cashflow_amts);
[h, axes_handle] = cfplot(cashflow_dates, cashflow_amts, 'Groups', 'off', 'ShowAmnt', 1, 'DateSpacing', 1, 'DateFormat', 1)