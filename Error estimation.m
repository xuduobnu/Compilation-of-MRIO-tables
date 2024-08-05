%% This code is used to calculate the estimation error of the technical matrix A
% Created by: Xu Duo, Beijing Normal University
% Created date: 27/07/2024
clear
clc
% Read official technical matrix
load A1997.mat
load A2002.mat
load A2007.mat
load A2012.mat
load A2017.mat

cd('C:\Users\ASUS\Desktop\new MRIO')
%% Matrix SSIM for 2003-2006
filename = 'Solved matrix from 03 to 06.xlsx'
Data03to06 = importdata(filename);
% Based on 2002
A2003_02 = Data03to06.data.x2003(2:1210,2:1171)*inv(diag(Data03to06.data.x2003(1213,2:1171)));
A2003_02(isnan(A2003_02)) = 0;
A2004_02 = Data03to06.data.x2004(2:1210,2:1171)*inv(diag(Data03to06.data.x2004(1213,2:1171)));
A2004_02(isnan(A2004_02)) = 0;
A2005_02 = Data03to06.data.x2005(2:1210,2:1171)*inv(diag(Data03to06.data.x2005(1213,2:1171)));
A2005_02(isnan(A2005_02)) = 0;
A2006_02 = Data03to06.data.x2006(2:1210,2:1171)*inv(diag(Data03to06.data.x2006(1213,2:1171)));
A2006_02(isnan(A2006_02)) = 0;
% Based on 2007
A2003_07 = Data03to06.data.x2003(1218:2426,2:1171)*inv(diag(Data03to06.data.x2003(2429,2:1171)));
A2003_07(isnan(A2003_07)) = 0;
A2004_07 = Data03to06.data.x2004(1218:2426,2:1171)*inv(diag(Data03to06.data.x2004(2429,2:1171)));
A2004_07(isnan(A2004_07)) = 0;
A2005_07 = Data03to06.data.x2005(1218:2426,2:1171)*inv(diag(Data03to06.data.x2005(2429,2:1171)));
A2005_07(isnan(A2005_07)) = 0;
A2006_07 = Data03to06.data.x2006(1218:2426,2:1171)*inv(diag(Data03to06.data.x2006(2429,2:1171)));
A2006_07(isnan(A2006_07)) = 0;

% SSIM
[ssimval, ssimmap] = ssim(A2007,A2006_07);
fprintf('SSIM: %f\n', ssimval);
figure;
imshow(ssimmap, []);  
title(['SSIM Index Map with Mean SSIM Value: ', num2str(ssimval)], 'FontSize', 12); 
axis on;
imageSize = size(ssimmap);
xticks(0:39:imageSize(2)-1);
yticks(0:39:imageSize(1)-1);
xticklabels(arrayfun(@num2str, 0:39:imageSize(2)-1, 'UniformOutput', false));
yticklabels(arrayfun(@num2str, 0:39:imageSize(1)-1, 'UniformOutput', false));
set(gca, 'TickDir', 'out'); 
colorbar;
c = colorbar;
c.Label.String = 'SSIM Value'; 
c.Label.FontSize = 12; 
c.Ticks = [0.1, 0.3, 0.5, 0.7, 0.9, 1]; 
set(gca, 'LooseInset', get(gca, 'TightInset')); 

%% Matrix SSIM for 2008-2011
filename = 'Solved matrix from 2008 to 2011.xlsx'
Data08to11 = importdata(filename);
% Based on 2007
A2008_07 = Data08to11.data.x2008(2:1210,2:1171)*inv(diag(Data08to11.data.x2008(1213,2:1171)));
A2008_07(isnan(A2008_07)) = 0;
A2009_07 = Data08to11.data.x2009(2:1210,2:1171)*inv(diag(Data08to11.data.x2009(1213,2:1171)));
A2009_07(isnan(A2009_07)) = 0;
A2010_07 = Data08to11.data.x2010(2:1210,2:1171)*inv(diag(Data08to11.data.x2010(1213,2:1171)));
A2010_07(isnan(A2010_07)) = 0;
A2011_07 = Data08to11.data.x2011(2:1210,2:1171)*inv(diag(Data08to11.data.x2011(1213,2:1171)));
A2011_07(isnan(A2011_07)) = 0;
% Based on 2012
A2008_12 = Data08to11.data.x2008(1218:2426,2:1171)*inv(diag(Data08to11.data.x2008(2429,2:1171)));
A2008_12(isnan(A2008_12)) = 0;
A2009_12 = Data08to11.data.x2009(1218:2426,2:1171)*inv(diag(Data08to11.data.x2009(2429,2:1171)));
A2009_12(isnan(A2009_12)) = 0;
A2010_12 = Data08to11.data.x2010(1218:2426,2:1171)*inv(diag(Data08to11.data.x2010(2429,2:1171)));
A2010_12(isnan(A2010_12)) = 0;
A2011_12 = Data08to11.data.x2011(1218:2426,2:1171)*inv(diag(Data08to11.data.x2011(2429,2:1171)));
A2011_12(isnan(A2011_12)) = 0;

[ssimval, ssimmap] = ssim(A2011_12,A2012);
fprintf('SSIM: %f\n', ssimval);
figure;
imshow(ssimmap, []);  
title(['SSIM Index Map with Mean SSIM Value: ', num2str(ssimval)], 'FontSize', 12); 
axis on;
imageSize = size(ssimmap);
xticks(0:39:imageSize(2)-1);
yticks(0:39:imageSize(1)-1);
xticklabels(arrayfun(@num2str, 0:39:imageSize(2)-1, 'UniformOutput', false));
yticklabels(arrayfun(@num2str, 0:39:imageSize(1)-1, 'UniformOutput', false));
set(gca, 'TickDir', 'out'); 
colorbar;
c = colorbar;
c.Label.String = 'SSIM Value';
c.Label.FontSize = 12; 
c.Ticks = [0.1, 0.3, 0.5, 0.7, 0.9, 1]; 
set(gca, 'LooseInset', get(gca, 'TightInset')); 


%% Matrix SSIM for 2013-2016
filename = 'Solved matrix from 2013 to 2016.xlsx'
Data13to16 = importdata(filename);
% Based on 2012
A2013_12 = Data13to16.data.x2013(2:1210,2:1171)*inv(diag(Data13to16.data.x2013(1213,2:1171)));
A2013_12(isnan(A2013_12)) = 0;
A2014_12 = Data13to16.data.x2014(2:1210,2:1171)*inv(diag(Data13to16.data.x2014(1213,2:1171)));
A2014_12(isnan(A2014_12)) = 0;
A2015_12 = Data13to16.data.x2015(2:1210,2:1171)*inv(diag(Data13to16.data.x2015(1213,2:1171)));
A2015_12(isnan(A2015_12)) = 0;
A2016_12 = Data13to16.data.x2016(2:1210,2:1171)*inv(diag(Data13to16.data.x2016(1213,2:1171)));
A2016_12(isnan(A2016_12)) = 0;
% Based on 2017
A2013_17 = Data13to16.data.x2013(1218:2426,2:1171)*inv(diag(Data13to16.data.x2013(2429,2:1171)));
A2013_17(isnan(A2013_17)) = 0;
A2014_17 = Data13to16.data.x2014(1218:2426,2:1171)*inv(diag(Data13to16.data.x2014(2429,2:1171)));
A2014_17(isnan(A2014_17)) = 0;
A2015_17 = Data13to16.data.x2015(1218:2426,2:1171)*inv(diag(Data13to16.data.x2015(2429,2:1171)));
A2015_17(isnan(A2015_17)) = 0;
A2016_17 = Data13to16.data.x2016(1218:2426,2:1171)*inv(diag(Data13to16.data.x2016(2429,2:1171)));
A2016_17(isnan(A2016_17)) = 0;

[ssimval, ssimmap] = ssim(A2016_17,A2017);
fprintf('SSIM: %f\n', ssimval);
figure;
imshow(ssimmap, []);  
title(['SSIM Index Map with Mean SSIM Value: ', num2str(ssimval)], 'FontSize', 12); 
axis on;
imageSize = size(ssimmap);
xticks(0:39:imageSize(2)-1);
yticks(0:39:imageSize(1)-1);
xticklabels(arrayfun(@num2str, 0:39:imageSize(2)-1, 'UniformOutput', false));
yticklabels(arrayfun(@num2str, 0:39:imageSize(1)-1, 'UniformOutput', false));
set(gca, 'TickDir', 'out'); 
colorbar;
c = colorbar;
c.Label.String = 'SSIM Value'; 
c.Label.FontSize = 12; 
c.Ticks = [0.1, 0.3, 0.5, 0.7, 0.9, 1]; 
set(gca, 'LooseInset', get(gca, 'TightInset')); 


%% Matrix SSIM for 1995,1996,2018-2020
filename = 'Solved matrix for other years.xlsx'
Dataothers = importdata(filename);

A1995_97 = Dataothers.data.x1995(2:1210,2:1171)*inv(diag(Dataothers.data.x1995(1213,2:1171)));
A1995_97(isnan(A1995_97)) = 0;

A1996_97 = Dataothers.data.x1996(2:1210,2:1171)*inv(diag(Dataothers.data.x1996(1213,2:1171)));
A1996_97(isnan(A1996_97)) = 0;

A2018_17 = Dataothers.data.x2018(2:1210,2:1171)*inv(diag(Dataothers.data.x2018(1213,2:1171)));
A2018_17(isnan(A2018_17)) = 0;

A2019_17 = Dataothers.data.x2019(2:1210,2:1171)*inv(diag(Dataothers.data.x2019(1213,2:1171)));
A2019_17(isnan(A2019_17)) = 0;

A2020_17 = Dataothers.data.x2020(2:1210,2:1171)*inv(diag(Dataothers.data.x2020(1213,2:1171)));
A2020_17(isnan(A2020_17)) = 0;

[ssimval, ssimmap] = ssim(A2020_17,A2017);
fprintf('SSIM: %f\n', ssimval);
figure;
imshow(ssimmap, []);  
title(['SSIM Index Map with Mean SSIM Value: ', num2str(ssimval)], 'FontSize', 12); 
axis on;
imageSize = size(ssimmap);
xticks(0:39:imageSize(2)-1);
yticks(0:39:imageSize(1)-1);
xticklabels(arrayfun(@num2str, 0:39:imageSize(2)-1, 'UniformOutput', false));
yticklabels(arrayfun(@num2str, 0:39:imageSize(1)-1, 'UniformOutput', false));
set(gca, 'TickDir', 'out'); 
colorbar;
c = colorbar;
c.Label.String = 'SSIM Value';
c.Label.FontSize = 12;
c.Ticks = [0.1, 0.3, 0.5, 0.7, 0.9, 1];
set(gca, 'LooseInset', get(gca, 'TightInset')); 


%% Matrix SSIM for 1998-2001
filename = 'Solved matrix from 1998 to 2001.xlsx'
Data98to01 = importdata(filename);
% Based on 1997
A1998_97 = Data98to01.data.x1998(2:1210,2:1171)*inv(diag(Data98to01.data.x1998(1213,2:1171)));
A1998_97(isnan(A1998_97)) = 0;

A1999_97 = Data98to01.data.x1999(2:1210,2:1171)*inv(diag(Data98to01.data.x1999(1213,2:1171)));
A1999_97(isnan(A1999_97)) = 0;

A2000_97 = Data98to01.data.x2000(2:1210,2:1171)*inv(diag(Data98to01.data.x2000(1213,2:1171)));
A2000_97(isnan(A2000_97)) = 0;

A2001_97 = Data98to01.data.x2001(2:1210,2:1171)*inv(diag(Data98to01.data.x2001(1213,2:1171)));
A2001_97(isnan(A2001_97)) = 0;
% Based on 2002
A1998_02 = Data98to01.data.x1998(1218:2426,2:1171)*inv(diag(Data98to01.data.x1998(2429,2:1171)));
A1998_02(isnan(A1998_02)) = 0;

A1999_02 = Data98to01.data.x1999(1218:2426,2:1171)*inv(diag(Data98to01.data.x1999(2429,2:1171)));
A1999_02(isnan(A1999_02)) = 0;

A2000_02 = Data98to01.data.x2000(1218:2426,2:1171)*inv(diag(Data98to01.data.x2000(2429,2:1171)));
A2000_02(isnan(A2000_02)) = 0;

A2001_02 = Data98to01.data.x2001(1218:2426,2:1171)*inv(diag(Data98to01.data.x2001(2429,2:1171)));
A2001_02(isnan(A2001_02)) = 0;

[ssimval, ssimmap] = ssim(A1998_97,A1997);
fprintf('SSIM: %f\n', ssimval);
figure;
imshow(ssimmap, []);  
title(['SSIM Index Map with Mean SSIM Value: ', num2str(ssimval)], 'FontSize', 12); 
axis on;
imageSize = size(ssimmap);
xticks(0:39:imageSize(2)-1);
yticks(0:39:imageSize(1)-1);
xticklabels(arrayfun(@num2str, 0:39:imageSize(2)-1, 'UniformOutput', false));
yticklabels(arrayfun(@num2str, 0:39:imageSize(1)-1, 'UniformOutput', false));
set(gca, 'TickDir', 'out'); 
colorbar;
c = colorbar;
c.Label.String = 'SSIM Value'; 
c.Label.FontSize = 12; 
c.Ticks = [0.1, 0.3, 0.5, 0.7, 0.9, 1]; 
set(gca, 'LooseInset', get(gca, 'TightInset')); 



% MSE
MSE = mean((A2020_17(:) - A2017(:)).^2);
% Correlation Coefficient
CorrCoeff = corr2(A2020_17, A2017);
% NMSE
NMSE = MSE / mean(A2017(:).^2);
% Cosine Similarity
CosineSimilarity = dot(A2020_17(:), A2017(:)) / (norm(A2020_17(:)) * norm(A2017(:)));
fprintf('MSE: %f\n', MSE);
fprintf('Correlation Coefficient: %f\n', CorrCoeff);
fprintf('NMSE: %f\n', NMSE);
fprintf('Cosine Similarity: %f\n', CosineSimilarity);









