% A = importdata('/Users/bill/Desktop/data_mining/data_visualization/Programming\ Assignment\ 1\ Data\ New/ExcelFormattedGISTEMPData2TXT.txt',' ',1);
A = importdata('/Users/bill/Desktop/data_mining/data_visualization/Programming Assignment 1 Data New/ExcelFormattedGISTEMPData2TXT.txt',' ',1);
Glob = importdata('/Users/bill/Desktop/data_mining/data_visualization/Programming Assignment 1 Data New/Glob.txt',' ',1);
Year = importdata('/Users/bill/Desktop/data_mining/data_visualization/Programming Assignment 1 Data New/Year.txt',' ',1);
NHem = importdata('/Users/bill/Desktop/data_mining/data_visualization/Programming Assignment 1 Data New/NHem.txt',' ',1);
SHem = importdata('/Users/bill/Desktop/data_mining/data_visualization/Programming Assignment 1 Data New/SHem.txt',' ',1);
N2490 = importdata('/Users/bill/Desktop/data_mining/data_visualization/Programming Assignment 1 Data New/24N90N.txt',' ',1);
N24S24 = importdata('/Users/bill/Desktop/data_mining/data_visualization/Programming Assignment 1 Data New/24N24S.txt',' ',1);
S2490 = importdata('/Users/bill/Desktop/data_mining/data_visualization/Programming Assignment 1 Data New/90S24S.txt',' ',1);
% disp(A.colheaders{1,2})
% disp(A.data(:,2))
plot(Year.data,Glob.data,'color',[0,0.2,0.8],'LineWidth',2)
hold on
% plot(Year.data,NHem.data,'color',[1,0.4,0.1])
% plot(Year.data,SHem.data, 'color',[1,0.1,0.9])
plot(Year.data,N2490.data, 'color',[1,0.2,0.1])
plot(Year.data,N24S24.data,'color',[0.4,0.1,0.3])
% plot(Year.data,N24S24.data,'color',[0.8,0.3,0.7])
plot(Year.data,S2490.data,'color',[0.2,0.9,0.4])
% legend('global','north hemisphere','south hemisphere','24N-90N','24N-24S','90S-24S','location','northwest')
legend('global','24N-90N','24N-24S','90S-24S','location','northwest')
xlabel('year')
ylabel('deviation')
ylh = get(gca,'ylabel');
gyl = get(ylh);                                                         % Object Information 'Position',ylp, 
ylp = get(ylh, 'Position');
set(ylh, 'Rotation',0, 'VerticalAlignment','middle', 'HorizontalAlignment','right')
ylim([-80,120])
title('Deviation of global and regional annual temperature')
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 1, 0.96]);
saveas(gcf,'/Users/bill/Desktop/data_mining/data_visualization/regionalTemperature.png')