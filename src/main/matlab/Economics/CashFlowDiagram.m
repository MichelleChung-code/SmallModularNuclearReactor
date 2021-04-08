% Example of how to do a cashflow diagram

% Cashflows and Dates should be read in from a CSV in the future when we
% get valid numvers

clc
close all

from_csv = readtable('Cashflows.csv');

cashflow_dates1 = table2array(from_csv(1:21,1));
cashflow_dates1 = datenum(num2str(cashflow_dates1),'yyyymmdd')';

cashflow_amts1 = table2array(from_csv(1:21,2)); 

figure(1);
[h, axes_handle] = cfplot(cashflow_dates1, cashflow_amts1', 'Groups', 'off', 'ShowAmnt', 1, 'DateSpacing', 1, 'DateFormat', 11)
set(gca,'YTickLabel',[]);
title('Cash Flow Diagram for Year 0 to 20')
xlabel('Year')
ylabel('Cashflow Amounts ($MM)')
savefig('Cashflow_Diagram1.fig')

cashflow_dates2 = table2array(from_csv(22:41,1));
cashflow_dates2 = datenum(num2str(cashflow_dates2),'yyyymmdd')';

cashflow_amts2 = table2array(from_csv(22:41,2)); 

figure(2)
[h, axes_handle] = cfplot(cashflow_dates2, cashflow_amts2', 'Groups', 'off', 'ShowAmnt', 1, 'DateSpacing', 1, 'DateFormat', 11)
set(gca,'YTickLabel',[]);
title('Cash Flow Diagram for Year 21 to 40')
xlabel('Year')
ylabel('Cashflow Amounts ($MM)')
savefig('Cashflow_Diagram2.fig')
