% Capital_carbon
% Created by: Xu Duo, Beijing Normal University
% Created date: 27/07/2024
%% data preparation
clear
clc
r = 30;%province
s = 39;%sector
ts = 2020-1994;

% Read MRIO data, integrate it to IOT
% Initialize structure array IOT
% Pre allocated structure array, containing 26 years
IOT(26).Y = []; % Storage of 1170 element columns for final use except for exports inventory changes
IOT(26).Z = []; % intermediate matrix
IOT(26).x = []; % total output
fileY = 'final use and export.xlsx'; 
fileZ = 'RAS solved matrix.xlsx'; 
filex = 'RAS solved matrix.xlsx'; 
% 指定数据范围
dataRangeY = 'B2434:DR3603'; % RC,UC,GC,GFCF+export
dataRangeZ = 'B2:ASA1171'; 
dataRangex = 'B1213:ASA1213';
years = 1995:2020; 
for idx = 1:length(years)
    year = years(idx);
    sheetName = sprintf('%d', year); % sheet name "x1995"
    % 从 Excel 文件中读取数据
    IOT(idx).Y = readmatrix(fileY, 'Sheet', sheetName, 'Range', dataRangeY);
    IOT(idx).Z = readmatrix(fileZ, 'Sheet', sheetName, 'Range', dataRangeZ);
    IOT(idx).x = readmatrix(filex, 'Sheet', sheetName, 'Range', dataRangex);
end

% Eliminate price fluctuations
filename = 'consistency matrix.xlsx';
CoormatrixData = importdata(filename);
cpi = CoormatrixData.data.cpi(1:26,end);
for idx = 1:length(IOT)
    factor = cpi(idx); 
    IOT(idx).Y = IOT(idx).Y / factor; 
    IOT(idx).Z = IOT(idx).Z / factor;
    IOT(idx).x = IOT(idx).x / factor; 
end

save('IOT.mat',"IOT")

%% Carbon intensity
filename = 'Solved carbon emission.xlsx'
CarbonData = importdata(filename);
years = 1997:2020;
fields = arrayfun(@(year) sprintf('x%d', year), years, 'UniformOutput', false);
CarbonDataT = cell(1, length(fields));% carbon emission from 1997 to 2020, 28 sectors
for i = 1:length(fields)
    field = fields{i};
    data = CarbonData.data.(field);
    data1 = data(1:25, 3:32);
    data2 = data(27:29, 3:32);
    matrix = [data1; data2];
    CarbonDataT{i} = matrix;
end

% Integrate total output, calculate intensity
transformedMatrices = cell(1, 24);
for i = 3:26
    x_matrix = IOT(i).x; 
    transformedMatrices{i-2} = reshape(x_matrix, [39, 30]);
end

B = cell(1, 24);
for i = 1:24
    transformedMatrix = transformedMatrices{i};
    resultMatrix = zeros(28, 30);
    for j = 1:25
        resultMatrix(j, :) = transformedMatrix(j, :);
    end
    resultMatrix(26, :) = transformedMatrix(27, :);
    resultMatrix(27, :) = transformedMatrix(26, :)+transformedMatrix(28, :);
    resultMatrix(28, :) = sum(transformedMatrix(29:39, :), 1);
    B{i} = resultMatrix;
end

%  Original co2 intensity C
C = cell(1, 24);
for i = 1:24
    matrix_A = CarbonDataT{i};
    matrix_B = B{i};
    result_matrix = (matrix_A ./ matrix_B)*1000000;%
    result_matrix(isnan(result_matrix) | isinf(result_matrix)) = 0;
    C{i} = result_matrix; 
end

CarbonIntensity = cell(1, 26);
CarbonIntensity{1} = C{1};
CarbonIntensity{2} = C{1};
for i = 3:26
    CarbonIntensity{i} = C{i-2};
end

for i = 1:26
    expanded_matrix = zeros(39, 30);
    expanded_matrix(1:25, :) = CarbonIntensity{i}(1:25, :);
    expanded_matrix(26, :) = CarbonIntensity{i}(27, :);
    expanded_matrix(27, :) = CarbonIntensity{i}(26, :);
    expanded_matrix(28, :) = CarbonIntensity{i}(27, :);
    expanded_matrix(29:39, :) = repmat(CarbonIntensity{i}(end, :), [11, 1]);
    CarbonIntensity{i} = expanded_matrix;
end
save('CarbonIntensity.mat','CarbonIntensity')

%% Calculate annual capital depreciation
load ki_chn.mat % Rural investment
load ki_chnU.mat % Urban investment
load depr.mat % Asset depreciation rate
load c_axp.mat % Matrix of correspondence between assets and sectors

ki_chnT = ki_chn + ki_chnU;
cd_ss = zeros(r,ts,ts,s,s);% first ts is the depreciation year, second ts is the capital formation year
for m  = 1:r
    for n = 1:ts % each capital consuming year
        for t = 1:n % each capital investment year
            Yk = IOT(t).Y(:,90+m);
            Ki = squeeze(ki_chnT(t,m,:,:))';
            gfcf_a = Ki;%39*2
            gfcf_a(isnan(gfcf_a))=0;
            gfcf_p = sum(reshape(Yk,[s,r]),2);
            gfcf_p(gfcf_p<0) = 0;
            cfc_v = IOT(n).cfc(:,(m-1)*s+1:m*s);

            % converting newly capital investment data based on IOT: from assets to the capital goods producing sectors
            c_axp_v1 = c_axp; % products producing assets: adjusted concordance - considering the distribution of sectors producing the same assets based on the production structure described in Yk 
            for a = 1:size(Ki,2)
                temp = find(c_axp(a,:)==1);
                c_axp_v1(a,temp) = c_axp(a,temp).*gfcf_p(temp)'./sum(gfcf_p(temp));
            end
            clear temp;
            c_axp_v1(isnan(c_axp_v1)) = 0; % 2*39
            gfcf_ap = gfcf_a*c_axp_v1; % 39*39
            gfcf_ap = gfcf_ap.*gfcf_p'./sum(gfcf_ap,'omitnan'); % scale to Yk
            gfcf_ap(isnan(gfcf_ap))=0; 
            clear c_axp_v1
            
            % adjust when Ki = 0 whereas gfcf of sector is not 0
            for j=1:s
                if gfcf_p(j)>0 && sum(gfcf_ap(:,j))==0
                    gfcf_ap(:,j) = gfcf_p(j).*(sum(gfcf_a,2)./sum(sum(gfcf_a)))';
                end
            end
            
            cd=zeros(size(Ki,1),s); %39*39
            for j = 1:s
                if sum(gfcf_ap(:,j))>0
                    temp = find(c_axp(:,j)>0);
                    temp1 = zeros(size(Ki,2),1); % 2*1, address: one sector producing more than one assets
                    if any(temp)
                        temp1(temp)=sum(gfcf_a(:,temp))./sum(sum(gfcf_a(:,temp)));
                        if sum(temp1,'omitnan')==0
                            temp1=repmat(1/size(Ki,2),[size(Ki,2),1]);
                        end
                    else
                        temp1= sum(gfcf_a,1)./sum(sum(gfcf_a,1));
                    end
                    cdr = ((1-depr).^(n-t)).*depr;
                    cd(:,j) = sum(gfcf_ap(:,j)*c_axp(:,j)'*diag(temp1).*cdr,2);
                    
                    clear cdr temp temp1
                end
            end               
            cd_ss(m,n,t,:,:) = cd'; 
            clear Yk impt_y Ki gfcf_a gfcf_p cfc_v gfcf_ap cd cd_n_t
        end
    end
end

save('cd_ss.mat','cd_ss')


%% supply chain-wide CO2 emissions embodied in capital depreciation
f_co2 = zeros(ts,s,r);
for i = 1:26
    f_co2(i, :, :) = CarbonIntensity{i};
end

Fk = zeros(ts,ts,r,s*r); % supply chain-wide carbon emissions of capital consumption, 1st for capital depreciation year; 2nd for capital investment year
Fka = zeros(ts,ts,r,s*r); % supply chain-wide carbon emissions of capital consumption, estimated by current technology assumption 
PBEk_p_s_cd  = zeros(ts,ts,r*s,r*s); % production based carbon emissions of capital depreciation
%3rd dimension: basic production sectors
%4rd dimension:regarding the capital depreciation by capital using sectors
for t = 1:ts % production year of the depreciated capital goods
    A = IOT(t).Z./IOT(t).x;
    A(isnan(A)) = 0;
    A(isinf(A)) = 0;
    L = inv(eye(r*s)-A);
    
    f0 = reshape(squeeze(f_co2(t,:,:)),r*s,1)';% 1*1170
    for n = t:ts % year of final consumption
        f0n = reshape(squeeze(f_co2(n,:,:)),r*s,1)'; % 1*1170
        for i = 1:r
            gfcf_p = IOT(t).Y(:,90+i); %i: investing province
%             gfcf_impt = Imp_Yk(t).Imp_Yk(:,i); % imported capital 
            if sum(sum(cd_ss(i,n,t,:,:),'omitnan'))>0 % capital goods produced in year t and depreciated in year n, invested by province i
                temp = squeeze(cd_ss(i,n,t,:,:)); % rows (dimension 4): capital goods producing sectors; columns: capital goods consuming sectors
                temp(isnan(temp))=0;
                temp0 = reshape(gfcf_p,[s,r]); % capital goods purchased in year t by province i; columns: province where the capital goods were purchased from
                temp0 = temp0./sum(temp0,2);
                temp0(isnan(temp0)) = 0;
                temp0(isinf(temp0)) = 0;
                temp1 = zeros(r*s,s);
                for m = 1:r % 'finished' capital goods producer
                    for j = 1:s % sectors that produced the 'finished' capital goods
                        temp1((m-1)*s+j,:) = temp0(j,m).*temp(j,:); % for the 'same' capital goods produced in the same year, assuming the cross-province distributions of depreciation and investment are idential
                    end
                end
                clear temp temp0
                %                     test_k1 =  test_k1+sum(sum(temp1,'omitnan'));
                
                temp2 = f0'.*L*temp1;
                for m = 1:r % m: provinces producing for the capital goods; i: provinces using the capital goods; columns: capital using sectors in province i
                    Fk(n,t,m,(i-1)*s+1:i*s) = sum(temp2((m-1)*s+1:m*s,:));
                end
                PBEk_p_s_cd(n,t,:,(i-1)*s+(1:s)) = temp2;
                
                temp2a = f0n'.*L*temp1;
                for m = 1:r % m: countries producing for the capital goods; i: countries using the capital goods; columns: capital using sectors in country i
                    Fka(n,t,m,(i-1)*s+1:i*s) = sum(temp2a((m-1)*s+1:m*s,:));
                end
                clear temp1 temp2 temp2a
            end
            clear gfcf_p
        end
        clear f0n
    end
    
    clear L A f0
    t
end

save('FkCO2.mat','Fk')
save('FkaCO2.mat','Fka')
save('PBEk_p_s_cdCO2.mat','PBEk_p_s_cd','-v7.3')

%% conventional consumption-based carbon emissions
for n = 1:ts
    for m = 1:r
        Yc(n).Yc(:,m) = IOT(n).Y(:,m)+IOT(n).Y(:,30+m)+IOT(n).Y(:,60+m);
        Yk(n).Yk(:,m) = IOT(n).Y(:,90+m); % GFCF列
    end
end

EF_c = zeros(ts,r,s,r,s); % conventional consumption-based carbon emissions of final consumption; 2nd-3rd: province of final consumption and sectors; 4th-5th: where impacts occurred
EF_gfcf = zeros(ts,r,s,r,s); % conventional consumption-based carbon emissions of gross fixed capital formation; 2nd-3rd: province of GFCF and sectors; 4th-5th: where impacts occurred
EF_expt = zeros(ts,r,s,s);
for n = 1:ts % year of capital investment
    tic
    A = IOT(n).Z./IOT(n).x;
    A(isnan(A)) = 0;
    A(isinf(A)) = 0;
    L = inv(eye(r*s)-A);
    
    yc = zeros(r*s,r*s);
    yk = zeros(r*s,r*s);
    c_m = repmat(eye(s),r,1);
    for m = 1:r
        yc(:,(m-1)*s+(1:s)) = Yc(n).Yc(:,m).*c_m;
        yk(:,(m-1)*s+(1:s)) = Yk(n).Yk(:,m).*c_m;
    end
    expt = IOT(n).Y(:,end).*c_m;
    
    f0 = reshape(squeeze(f_co2(n,:,:)),s*r,1)';
    temp_ef_c = f0'.*L*yc;
    temp_ef_gfcf = f0'.*L*yk;
    temp_ef_expt = f0'.*L*expt;
    for m = 1:r % region where EP occur
        for i = 1:r % region where final goods consumed
            EF_c(n,i,:,m,:) = temp_ef_c((m-1)*s+(1:s),(i-1)*s+(1:s))';
            EF_gfcf(n,i,:,m,:) = temp_ef_gfcf((m-1)*s+(1:s),(i-1)*s+(1:s))';
        end
        EF_expt(n,m,:,:) = temp_ef_expt((m-1)*s+(1:s),:);
    end
    clear temp_ef_c temp_ef_gfcf f0 temp_ef_expt
        n
    toc
    clear A L y
end

save('EF_c_CO2.mat','EF_c')
save('EF_gfcf_CO2.mat','EF_gfcf')
save('EF_expt_CO2.mat','EF_expt')

%% Capital-associate carbon emissions
PBEk_p_s_cdNew = squeeze(sum(PBEk_p_s_cd,1));%26*1170*1170
xnew = struct('x', cell(1, 26)); 
for i = 1:length(IOT)
    xnew(i).x = 1 ./ IOT(i).x;
end
for i = 1:length(xnew)
    invalidIndices = isnan(xnew(i).x) | isinf(xnew(i).x);
    xnew(i).x(invalidIndices) = 0;
end

fknnew = zeros(26,1170,1170);
five = ones(1, 1170);
for i = 1:26
    fknnew(i,:,:) = squeeze(PBEk_p_s_cdNew(i,:,:))*diag(xnew(i).x);
end
fknnew(isnan(fknnew))=0;
fknnew(isinf(fknnew))=0;
fknnew2 = zeros(26,1170);
for i = 1:26
    fknnew2(i,:) = five*squeeze(fknnew(i,:,:));
end

for n = 1:ts
    for m = 1:r
        Yc(n).Yc(:,m) = IOT(n).Y(:,m)+IOT(n).Y(:,30+m)+IOT(n).Y(:,60+m);
        Yk(n).Yk(:,m) = IOT(n).Y(:,90+m); % GFCF列
    end
end

EFk_c = zeros(ts,r,s,r,s); 
EFk_gfcf = zeros(ts,r,s,r,s);
EFk_expt = zeros(ts,r,s,s);

for n = 1:ts % year of capital investment
    tic
    A = IOT(n).Z./IOT(n).x;
    A(isnan(A)) = 0;
    A(isinf(A)) = 0;
    L = inv(eye(r*s)-A);
    
    yc = zeros(r*s,r*s);
    yk = zeros(r*s,r*s);
    c_m = repmat(eye(s),r,1);
    for m = 1:r
        yc(:,(m-1)*s+(1:s)) = Yc(n).Yc(:,m).*c_m;
        yk(:,(m-1)*s+(1:s)) = Yk(n).Yk(:,m).*c_m;
    end
    expt = IOT(n).Y(:,end).*c_m;
    
    f0 = fknnew2(n,:);
    temp_ef_c = f0'.*L*yc;
    temp_ef_gfcf = f0'.*L*yk;
    temp_ef_expt = f0'.*L*expt;
    for m = 1:r % region where EP occur
        for i = 1:r % region where final goods consumed
            EFk_c(n,i,:,m,:) = temp_ef_c((m-1)*s+(1:s),(i-1)*s+(1:s))';
            EFk_gfcf(n,i,:,m,:) = temp_ef_gfcf((m-1)*s+(1:s),(i-1)*s+(1:s))';
        end
        EFk_expt(n,m,:,:) = temp_ef_expt((m-1)*s+(1:s),:);
    end
    clear temp_ef_c temp_ef_gfcf f0 temp_ef_expt
        n
    toc
    clear A L y
end

save('EFk_c_CO2.mat','EFk_c')
save('EFk_gfcf_CO2.mat','EFk_gfcf')
save('EFk_expt_CO2.mat','EFk_expt')
