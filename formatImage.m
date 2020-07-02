clear all;
fname = 'axicon_10.fig';
folder_path = 'D:/GoogleDrive/Thesis/MATLAB/';
fig = [folder_path,fname];
inf = openfig(fig);

C = gca;

%truesize(figure(1),[500,1500]);
title('');
xlabel('Distance (\mum)');
xticks([0,100,200,300,400,500]);
ylabel('x (\mum)');
yticks([-50,-25,0,25,50]);
set(C,'FontSize',20);
set(figure(1),'position',[342,438,1400,400]);
text(C,400,40,['\beta = 10',char(176)],'FontSize',36,'FontWeight','bold','color','w');
savefig([folder_path,fname]);
saveas(gcf,[folder_path,'axicon_10.png']);