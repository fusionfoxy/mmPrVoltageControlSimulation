% For plotting the consumption of single households
clc; clear all; close all;

load('consData/house1to120days7.mat');

xdate = datenum(Data.time,'yy-mm-dd HH:MM:SS');

days = 1;
fromH = 115;
toH = 117;
N = days*24*4;

figure
plot(xdate(1:N),Data.HouseP(1:N,fromH:toH)*1000)
grid
datetick('x','HH:MM')
title('Consumption of 3 Houses. 15 min. sampling')
xlabel('Time [hrs:min]')
ylabel('Power [W]')

figure
plot(xdate(1:N),sum(Data.HouseP(1:N,fromH:toH),2))
datetick('x','HH:MM')
xlabel('Time [hrs:min]')
ylabel('Power [kW]')
