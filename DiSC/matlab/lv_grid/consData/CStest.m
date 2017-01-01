clc
clear
close all


load house1to120days2


stairs(0:15:15*(size(Data.time,1)-1),Data.HouseP(:,1:3))

ylabel('Consumption [kW]')
xlabel('Time [min]')