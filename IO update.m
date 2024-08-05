clear all
clc
%% This program is used for updating the intermediate matrix in MRIO tables from 1995 to 2020
% Created by: Xu Duo, Beijing Normal University
% Created date: 27/07/2024
% Read the processed data from 1998 to 2016 first
filename = 'Solved total output and value added from 98 to 16.xlsx';
OtVadata = importdata(filename);

% Read official MRIO tables data
filename2 = '3 1997IOtable.xlsx';
filename3 = '4 2002IOtable.xlsx';
filename4 = '5 2007IOtable.xlsx';
filename5 = '6 2012IOtable.xlsx';
filename6 = '7 2017IOtable.xlsx';
sheet2 = 'MRIO1997';
sheet3 = 'MRIO2002';
sheet4 = 'MRIO2007';
sheet5 = 'MRIO2012';
sheet6 = 'MRIO2017';
IO1997data = importdata(filename2,sheet2);
IO2002data = importdata(filename3,sheet3);
IO2007data = importdata(filename4,sheet4);
IO2012data = importdata(filename5,sheet5);
IO2017data = importdata(filename6,sheet6);

%% Obtain intermediate use vectors based on different base tables
% 1997 is the base table
interuse1998 = IO1997data.data.MRIO19970x65B0(1:1209,1171).*(1333003307/sum(IO1997data.data.MRIO19970x65B0(1:1209,1171)));
interuse1999 = IO1997data.data.MRIO19970x65B0(1:1209,1171).*(1425126921/sum(IO1997data.data.MRIO19970x65B0(1:1209,1171)));
interuse2000 = IO1997data.data.MRIO19970x65B0(1:1209,1171).*(1583933029/sum(IO1997data.data.MRIO19970x65B0(1:1209,1171)));
interuse2001 = IO1997data.data.MRIO19970x65B0(1:1209,1171).*(1747856542/sum(IO1997data.data.MRIO19970x65B0(1:1209,1171)));

% 2002 is the base table
interuse2003 = IO2002data.data.MRIO20020x65B0(1:1209,1171).*(2392466167/sum(IO2002data.data.MRIO20020x65B0(1:1209,1171)));
interuse2004 = IO2002data.data.MRIO20020x65B0(1:1209,1171).*(2863173434/sum(IO2002data.data.MRIO20020x65B0(1:1209,1171)));
interuse2005 = IO2002data.data.MRIO20020x65B0(1:1209,1171).*(3306536852/sum(IO2002data.data.MRIO20020x65B0(1:1209,1171)));
interuse2006 = IO2002data.data.MRIO20020x65B0(1:1209,1171).*(3961049458/sum(IO2002data.data.MRIO20020x65B0(1:1209,1171)));

% 2007 is the base table
interuse2008 = IO2007data.data.MRIO20070x65B0(1:1209,1171).*(6020998361/sum(IO2007data.data.MRIO20070x65B0(1:1209,1171)));
interuse2009 = IO2007data.data.MRIO20070x65B0(1:1209,1171).*(6551132852/sum(IO2007data.data.MRIO20070x65B0(1:1209,1171)));
interuse2010 = IO2007data.data.MRIO20070x65B0(1:1209,1171).*(7911013648/sum(IO2007data.data.MRIO20070x65B0(1:1209,1171)));
interuse2011 = IO2007data.data.MRIO20070x65B0(1:1209,1171).*(9372792332/sum(IO2007data.data.MRIO20070x65B0(1:1209,1171)));
 
% 2012 is the base table
interuse2013 = IO2012data.data.MRIO20120x65B0(1:1209,1171).*(12089273813/sum(IO2012data.data.MRIO20120x65B0(1:1209,1171)));
interuse2014 = IO2012data.data.MRIO20120x65B0(1:1209,1171).*(12404535364/sum(IO2012data.data.MRIO20120x65B0(1:1209,1171)));
interuse2015 = IO2012data.data.MRIO20120x65B0(1:1209,1171).*(13248468553/sum(IO2012data.data.MRIO20120x65B0(1:1209,1171)));
interuse2016 = IO2012data.data.MRIO20120x65B0(1:1209,1171).*(13856555170/sum(IO2012data.data.MRIO20120x65B0(1:1209,1171)));

% 2002 is the base table
interuse1998 = IO2002data.data.MRIO20020x65B0(1:1209,1171).*(1421103719/sum(IO2002data.data.MRIO20020x65B0(1:1209,1171)));
interuse1999 = IO2002data.data.MRIO20020x65B0(1:1209,1171).*(1507590206/sum(IO2002data.data.MRIO20020x65B0(1:1209,1171)));
interuse2000 = IO2002data.data.MRIO20020x65B0(1:1209,1171).*(1667710755/sum(IO2002data.data.MRIO20020x65B0(1:1209,1171)));
interuse2001 = IO2002data.data.MRIO20020x65B0(1:1209,1171).*(1822515741/sum(IO2002data.data.MRIO20020x65B0(1:1209,1171)));

% 2007 is the base table
interuse2003 = IO2007data.data.MRIO20070x65B0(1:1209,1171).*(2440776199/sum(IO2007data.data.MRIO20070x65B0(1:1209,1171)));
interuse2004 = IO2007data.data.MRIO20070x65B0(1:1209,1171).*(2946603948/sum(IO2007data.data.MRIO20070x65B0(1:1209,1171)));
interuse2005 = IO2007data.data.MRIO20070x65B0(1:1209,1171).*(3472078209/sum(IO2007data.data.MRIO20070x65B0(1:1209,1171)));
interuse2006 = IO2007data.data.MRIO20070x65B0(1:1209,1171).*(4122558119/sum(IO2007data.data.MRIO20070x65B0(1:1209,1171)));

% 2012 is the base table
interuse2008 = IO2012data.data.MRIO20120x65B0(1:1209,1171).*(6501488186/sum(IO2012data.data.MRIO20120x65B0(1:1209,1171)));
interuse2009 = IO2012data.data.MRIO20120x65B0(1:1209,1171).*(7094166929/sum(IO2012data.data.MRIO20120x65B0(1:1209,1171)));
interuse2010 = IO2012data.data.MRIO20120x65B0(1:1209,1171).*(8573060766/sum(IO2012data.data.MRIO20120x65B0(1:1209,1171)));
interuse2011 = IO2012data.data.MRIO20120x65B0(1:1209,1171).*(10140169979/sum(IO2012data.data.MRIO20120x65B0(1:1209,1171)));

% 2017 is the base table
interuse2013 = IO2017data.data.MRIO20170x65B0(1:1209,1171).*(13005172764/sum(IO2017data.data.MRIO20170x65B0(1:1209,1171)));
interuse2014 = IO2017data.data.MRIO20170x65B0(1:1209,1171).*(14617079179/sum(IO2017data.data.MRIO20170x65B0(1:1209,1171)));
interuse2015 = IO2017data.data.MRIO20170x65B0(1:1209,1171).*(14385502913/sum(IO2017data.data.MRIO20170x65B0(1:1209,1171)));
interuse2016 = IO2017data.data.MRIO20170x65B0(1:1209,1171).*(15096992231/sum(IO2017data.data.MRIO20170x65B0(1:1209,1171)));

%% Read solved intermediate input vector
% 1997 is the base table
interinput1998 = reshape(OtVadata.data.x0x4E240x7C7B0x68C00x9A8C(1:39,1:30),1,1170);
interinput1999 = reshape(OtVadata.data.x0x4E240x7C7B0x68C00x9A8C(43:81,1:30),1,1170);
interinput2000 = reshape(OtVadata.data.x0x4E240x7C7B0x68C00x9A8C(85:123,1:30),1,1170);
interinput2001 = reshape(OtVadata.data.x0x4E240x7C7B0x68C00x9A8C(127:165,1:30),1,1170);
 
% 2002 is the base table
interinput2003 = reshape(OtVadata.data.x0x4E240x7C7B0x68C00x9A8C(169:207,34:63),1,1170);
interinput2004 = reshape(OtVadata.data.x0x4E240x7C7B0x68C00x9A8C(211:249,34:63),1,1170);
interinput2005 = reshape(OtVadata.data.x0x4E240x7C7B0x68C00x9A8C(253:291,34:63),1,1170);
interinput2006 = reshape(OtVadata.data.x0x4E240x7C7B0x68C00x9A8C(295:333,34:63),1,1170);
 
% 2007 is the base table
interinput2008 = reshape(OtVadata.data.x0x4E240x7C7B0x68C00x9A8C(169:207,67:96),1,1170);
interinput2009 = reshape(OtVadata.data.x0x4E240x7C7B0x68C00x9A8C(211:249,67:96),1,1170);
interinput2010 = reshape(OtVadata.data.x0x4E240x7C7B0x68C00x9A8C(253:291,67:96),1,1170);
interinput2011 = reshape(OtVadata.data.x0x4E240x7C7B0x68C00x9A8C(295:333,67:96),1,1170);
 
% 2012 is the base table
interinput2013 = reshape(OtVadata.data.x0x4E240x7C7B0x68C00x9A8C(169:207,100:129),1,1170);
interinput2014 = reshape(OtVadata.data.x0x4E240x7C7B0x68C00x9A8C(211:249,100:129),1,1170);
interinput2015 = reshape(OtVadata.data.x0x4E240x7C7B0x68C00x9A8C(253:291,100:129),1,1170);
interinput2016 = reshape(OtVadata.data.x0x4E240x7C7B0x68C00x9A8C(295:333,100:129),1,1170);

% 2002 is the base table
interinput1998 = reshape(OtVadata.data.x0x4E240x7C7B0x68C00x9A8C(1:39,34:63),1,1170);
interinput1999 = reshape(OtVadata.data.x0x4E240x7C7B0x68C00x9A8C(43:81,34:63),1,1170);
interinput2000 = reshape(OtVadata.data.x0x4E240x7C7B0x68C00x9A8C(85:123,34:63),1,1170);
interinput2001 = reshape(OtVadata.data.x0x4E240x7C7B0x68C00x9A8C(127:165,34:63),1,1170);

% 2007 is the base table
interinput2003 = reshape(OtVadata.data.x0x4E240x7C7B0x68C00x9A8C(1:39,67:96),1,1170);
interinput2004 = reshape(OtVadata.data.x0x4E240x7C7B0x68C00x9A8C(43:81,67:96),1,1170);
interinput2005 = reshape(OtVadata.data.x0x4E240x7C7B0x68C00x9A8C(85:123,67:96),1,1170);
interinput2006 = reshape(OtVadata.data.x0x4E240x7C7B0x68C00x9A8C(127:165,67:96),1,1170);

% 2012 is the base table
interinput2008 = reshape(OtVadata.data.x0x4E240x7C7B0x68C00x9A8C(1:39,100:129),1,1170);
interinput2009 = reshape(OtVadata.data.x0x4E240x7C7B0x68C00x9A8C(43:81,100:129),1,1170);
interinput2010 = reshape(OtVadata.data.x0x4E240x7C7B0x68C00x9A8C(85:123,100:129),1,1170);
interinput2011 = reshape(OtVadata.data.x0x4E240x7C7B0x68C00x9A8C(127:165,100:129),1,1170);

% 2017 is the base table
interinput2013 = reshape(OtVadata.data.x0x4E240x7C7B0x68C00x9A8C(1:39,133:162),1,1170);
interinput2014 = reshape(OtVadata.data.x0x4E240x7C7B0x68C00x9A8C(43:81,133:162),1,1170);
interinput2015 = reshape(OtVadata.data.x0x4E240x7C7B0x68C00x9A8C(85:123,133:162),1,1170);
interinput2016 = reshape(OtVadata.data.x0x4E240x7C7B0x68C00x9A8C(127:165,133:162),1,1170);

%% Read solved total input vector
% 1997 is the base table
Ot1998 = reshape(OtVadata.data.x0x4E240x7C7B0x603B0x4EA70x503C(1:39,1:30),1,1170);
Ot1999 = reshape(OtVadata.data.x0x4E240x7C7B0x603B0x4EA70x503C(43:81,1:30),1,1170);
Ot2000 = reshape(OtVadata.data.x0x4E240x7C7B0x603B0x4EA70x503C(85:123,1:30),1,1170);
Ot2001 = reshape(OtVadata.data.x0x4E240x7C7B0x603B0x4EA70x503C(127:165,1:30),1,1170);

% 2002 is the base table
Ot2003 = reshape(OtVadata.data.x0x4E240x7C7B0x603B0x4EA70x503C(169:207,34:63),1,1170);
Ot2004 = reshape(OtVadata.data.x0x4E240x7C7B0x603B0x4EA70x503C(211:249,34:63),1,1170);
Ot2005 = reshape(OtVadata.data.x0x4E240x7C7B0x603B0x4EA70x503C(253:291,34:63),1,1170);
Ot2006 = reshape(OtVadata.data.x0x4E240x7C7B0x603B0x4EA70x503C(295:333,34:63),1,1170);
 
% 2007 is the base table
Ot2008 = reshape(OtVadata.data.x0x4E240x7C7B0x603B0x4EA70x503C(169:207,67:96),1,1170);
Ot2009 = reshape(OtVadata.data.x0x4E240x7C7B0x603B0x4EA70x503C(211:249,67:96),1,1170);
Ot2010 = reshape(OtVadata.data.x0x4E240x7C7B0x603B0x4EA70x503C(253:291,67:96),1,1170);
Ot2011 = reshape(OtVadata.data.x0x4E240x7C7B0x603B0x4EA70x503C(295:333,67:96),1,1170);
 
% 2012 is the base table
Ot2013 = reshape(OtVadata.data.x0x4E240x7C7B0x603B0x4EA70x503C(169:207,100:129),1,1170);
Ot2014 = reshape(OtVadata.data.x0x4E240x7C7B0x603B0x4EA70x503C(211:249,100:129),1,1170);
Ot2015 = reshape(OtVadata.data.x0x4E240x7C7B0x603B0x4EA70x503C(253:291,100:129),1,1170);
Ot2016 = reshape(OtVadata.data.x0x4E240x7C7B0x603B0x4EA70x503C(295:333,100:129),1,1170);

% 2002 is the base table
Ot1998 = reshape(OtVadata.data.x0x4E240x7C7B0x603B0x4EA70x503C(1:39,34:63),1,1170);
Ot1999 = reshape(OtVadata.data.x0x4E240x7C7B0x603B0x4EA70x503C(43:81,34:63),1,1170);
Ot2000 = reshape(OtVadata.data.x0x4E240x7C7B0x603B0x4EA70x503C(85:123,34:63),1,1170);
Ot2001 = reshape(OtVadata.data.x0x4E240x7C7B0x603B0x4EA70x503C(127:165,34:63),1,1170);

% 2007 is the base table
Ot2003 = reshape(OtVadata.data.x0x4E240x7C7B0x603B0x4EA70x503C(1:39,67:96),1,1170);
Ot2004 = reshape(OtVadata.data.x0x4E240x7C7B0x603B0x4EA70x503C(43:81,67:96),1,1170);
Ot2005 = reshape(OtVadata.data.x0x4E240x7C7B0x603B0x4EA70x503C(85:123,67:96),1,1170);
Ot2006 = reshape(OtVadata.data.x0x4E240x7C7B0x603B0x4EA70x503C(127:165,67:96),1,1170);

% 2012 is the base table
Ot2008 = reshape(OtVadata.data.x0x4E240x7C7B0x603B0x4EA70x503C(1:39,100:129),1,1170);
Ot2009 = reshape(OtVadata.data.x0x4E240x7C7B0x603B0x4EA70x503C(43:81,100:129),1,1170);
Ot2010 = reshape(OtVadata.data.x0x4E240x7C7B0x603B0x4EA70x503C(85:123,100:129),1,1170);
Ot2011 = reshape(OtVadata.data.x0x4E240x7C7B0x603B0x4EA70x503C(127:165,100:129),1,1170);

% 2017 is the base table
Ot2013 = reshape(OtVadata.data.x0x4E240x7C7B0x603B0x4EA70x503C(1:39,133:162),1,1170);
Ot2014 = reshape(OtVadata.data.x0x4E240x7C7B0x603B0x4EA70x503C(43:81,133:162),1,1170);
Ot2015 = reshape(OtVadata.data.x0x4E240x7C7B0x603B0x4EA70x503C(85:123,133:162),1,1170);
Ot2016 = reshape(OtVadata.data.x0x4E240x7C7B0x603B0x4EA70x503C(127:165,133:162),1,1170);

% Read official technical matrix
load A1997.mat
load A2002.mat
load A2007.mat
load A2012.mat
load A2017.mat

%% RAS
[m,n]=size(A2017);
U1=[];
R1=[];
CsIx=A2017.*(Ot2016); 
for i=1:m
    U1(i)=sum(CsIx(i,:));% row total value firstly
    R1(i)=interuse2016(i)/U1(i);% row adjustment coefficient firstly
end
R1(isnan(R1)) = 0;
R1(isinf(R1)) = 0;
Z_R1=diag(R1)*CsIx;%first intermediate matrix after row adjust
for j=1:n
    V1(j)=sum(Z_R1(:,j));% column total value firstly
    S1(j)=interinput2016(j)/V1(j);% column adjustment coefficient firstly
end
S1(isnan(S1)) = 0;
S1(isinf(S1)) = 0;
Z_S1=Z_R1*diag(S1);%first intermediate matrix after column adjust


%% Set the number of iterations, store the errors for each iteration, and output the adjusted intermediate matrix
count=1;
countin=1;
Erro=[];
Uro=[];
Cro=[];
Ux=[];
Rx=[];
Vx=[];
Sx=[];
Z_Sx=Z_S1;
while count<100
    for i=1:m
        Ux(i)=sum(Z_Sx(i,:));
        Rx(i)=interuse2016(i)/Ux(i);
        if isnan(Rx(i)) || isinf(Rx(i))
            Rx(i) = 0;
        end
        Uro(i)=(Rx(i)-1)^2;%row error
    end
    Zrx=diag(Rx)*Z_Sx;
    for j=1:n
        Vx(j)=sum(Zrx(:,j));
        Sx(j)=interinput2016(j)/Vx(j);
        if isnan(Sx(j)) || isinf(Sx(j))
            Sx(j) = 0;
        end
        Cro(j)=(Sx(j)-1)^2;%column error
    end
    Erro(count)=sum(Uro)+sum(Cro);
    Z_Sx=Zrx*diag(Sx);
    count=count+1;
end
[ExPv,Exp]=min(Erro);%minimum error value and its location
Z_Sx=Z_S1;
t=1;
while t<=Exp
    for i=1:m
        Ux(i)=sum(Z_Sx(i,:));
        Rx(i)=interuse2016(i)/Ux(i);
        if isnan(Rx(i)) || isinf(Rx(i))
            Rx(i) = 0;
        end
    end
    Zrx=diag(Rx)*Z_Sx;
    for j=1:n
        Vx(j)=sum(Zrx(:,j));
        Sx(j)=interinput2016(j)/Vx(j);
        if isnan(Sx(j)) || isinf(Sx(j))
            Sx(j) = 0;
        end
    end
    Z_Sx=Zrx*diag(Sx);% intermediate matrix with the minimum error
    t=t+1;
end

%% Read weighted value added 
Valaddnew98 = reshape(OtVadata.data.x0x4E240x7C7B0x589E0x52A00x503C(169:207,1:30),1,1170);
Valaddnew99 = reshape(OtVadata.data.x0x4E240x7C7B0x589E0x52A00x503C(211:249,1:30),1,1170);
Valaddnew00 = reshape(OtVadata.data.x0x4E240x7C7B0x589E0x52A00x503C(253:291,1:30),1,1170);
Valaddnew01 = reshape(OtVadata.data.x0x4E240x7C7B0x589E0x52A00x503C(295:333,1:30),1,1170);

Valaddnew03 = reshape(OtVadata.data.x0x4E240x7C7B0x589E0x52A00x503C(169:207,67:96),1,1170);
Valaddnew04 = reshape(OtVadata.data.x0x4E240x7C7B0x589E0x52A00x503C(211:249,67:96),1,1170);
Valaddnew05 = reshape(OtVadata.data.x0x4E240x7C7B0x589E0x52A00x503C(253:291,67:96),1,1170);
Valaddnew06 = reshape(OtVadata.data.x0x4E240x7C7B0x589E0x52A00x503C(295:333,67:96),1,1170);

Valaddnew08 = reshape(OtVadata.data.x0x4E240x7C7B0x589E0x52A00x503C(169:207,133:162),1,1170);
Valaddnew09 = reshape(OtVadata.data.x0x4E240x7C7B0x589E0x52A00x503C(211:249,133:162),1,1170);
Valaddnew10 = reshape(OtVadata.data.x0x4E240x7C7B0x589E0x52A00x503C(253:291,133:162),1,1170);
Valaddnew11 = reshape(OtVadata.data.x0x4E240x7C7B0x589E0x52A00x503C(295:333,133:162),1,1170);

Valaddnew13 = reshape(OtVadata.data.x0x4E240x7C7B0x589E0x52A00x503C(169:207,199:228),1,1170);
Valaddnew14 = reshape(OtVadata.data.x0x4E240x7C7B0x589E0x52A00x503C(211:249,199:228),1,1170);
Valaddnew15 = reshape(OtVadata.data.x0x4E240x7C7B0x589E0x52A00x503C(253:291,199:228),1,1170);
Valaddnew16 = reshape(OtVadata.data.x0x4E240x7C7B0x589E0x52A00x503C(295:333,199:228),1,1170);

% Read unweighted value added
Valaddnew98_97 = reshape(OtVadata.data.x0x4E240x7C7B0x589E0x52A00x503C(1:39,1:30),1,1170);
Valaddnew99_97 = reshape(OtVadata.data.x0x4E240x7C7B0x589E0x52A00x503C(43:81,1:30),1,1170);
Valaddnew00_97 = reshape(OtVadata.data.x0x4E240x7C7B0x589E0x52A00x503C(85:123,1:30),1,1170);
Valaddnew01_97 = reshape(OtVadata.data.x0x4E240x7C7B0x589E0x52A00x503C(127:165,1:30),1,1170);

Valaddnew98_02 = reshape(OtVadata.data.x0x4E240x7C7B0x589E0x52A00x503C(1:39,34:63),1,1170);
Valaddnew99_02 = reshape(OtVadata.data.x0x4E240x7C7B0x589E0x52A00x503C(43:81,34:63),1,1170);
Valaddnew00_02 = reshape(OtVadata.data.x0x4E240x7C7B0x589E0x52A00x503C(85:123,34:63),1,1170);
Valaddnew01_02 = reshape(OtVadata.data.x0x4E240x7C7B0x589E0x52A00x503C(127:165,34:63),1,1170);

Valaddnew03_02 = reshape(OtVadata.data.x0x4E240x7C7B0x589E0x52A00x503C(1:39,67:96),1,1170);
Valaddnew04_02 = reshape(OtVadata.data.x0x4E240x7C7B0x589E0x52A00x503C(43:81,67:96),1,1170);
Valaddnew05_02 = reshape(OtVadata.data.x0x4E240x7C7B0x589E0x52A00x503C(85:123,67:96),1,1170);
Valaddnew06_02 = reshape(OtVadata.data.x0x4E240x7C7B0x589E0x52A00x503C(127:165,67:96),1,1170);

Valaddnew03_07 = reshape(OtVadata.data.x0x4E240x7C7B0x589E0x52A00x503C(1:39,100:129),1,1170);
Valaddnew04_07 = reshape(OtVadata.data.x0x4E240x7C7B0x589E0x52A00x503C(43:81,100:129),1,1170);
Valaddnew05_07 = reshape(OtVadata.data.x0x4E240x7C7B0x589E0x52A00x503C(85:123,100:129),1,1170);
Valaddnew06_07 = reshape(OtVadata.data.x0x4E240x7C7B0x589E0x52A00x503C(127:165,100:129),1,1170);

Valaddnew08_07 = reshape(OtVadata.data.x0x4E240x7C7B0x589E0x52A00x503C(1:39,133:162),1,1170);
Valaddnew09_07 = reshape(OtVadata.data.x0x4E240x7C7B0x589E0x52A00x503C(43:81,133:162),1,1170);
Valaddnew10_07 = reshape(OtVadata.data.x0x4E240x7C7B0x589E0x52A00x503C(85:123,133:162),1,1170);
Valaddnew11_07 = reshape(OtVadata.data.x0x4E240x7C7B0x589E0x52A00x503C(127:165,133:162),1,1170);

Valaddnew08_12 = reshape(OtVadata.data.x0x4E240x7C7B0x589E0x52A00x503C(1:39,166:195),1,1170);
Valaddnew09_12 = reshape(OtVadata.data.x0x4E240x7C7B0x589E0x52A00x503C(43:81,166:195),1,1170);
Valaddnew10_12 = reshape(OtVadata.data.x0x4E240x7C7B0x589E0x52A00x503C(85:123,166:195),1,1170);
Valaddnew11_12 = reshape(OtVadata.data.x0x4E240x7C7B0x589E0x52A00x503C(127:165,166:195),1,1170);

Valaddnew13_12 = reshape(OtVadata.data.x0x4E240x7C7B0x589E0x52A00x503C(1:39,199:228),1,1170);
Valaddnew14_12 = reshape(OtVadata.data.x0x4E240x7C7B0x589E0x52A00x503C(43:81,199:228),1,1170);
Valaddnew15_12 = reshape(OtVadata.data.x0x4E240x7C7B0x589E0x52A00x503C(85:123,199:228),1,1170);
Valaddnew16_12 = reshape(OtVadata.data.x0x4E240x7C7B0x589E0x52A00x503C(127:165,199:228),1,1170);

Valaddnew13_17 = reshape(OtVadata.data.x0x4E240x7C7B0x589E0x52A00x503C(1:39,232:261),1,1170);
Valaddnew14_17 = reshape(OtVadata.data.x0x4E240x7C7B0x589E0x52A00x503C(43:81,232:261),1,1170);
Valaddnew15_17 = reshape(OtVadata.data.x0x4E240x7C7B0x589E0x52A00x503C(85:123,232:261),1,1170);
Valaddnew16_17 = reshape(OtVadata.data.x0x4E240x7C7B0x589E0x52A00x503C(127:165,232:261),1,1170);

%% 1995、1996、2018、2019、2020 RAS
filename = 'Solved total output and value added before 97 after 17.xlsx';
OtVadata = importdata(filename);

interinput1995 = reshape(OtVadata.data.x0x68C00x9A8C(1:39,1:30),1,1170);
interinput1996 = reshape(OtVadata.data.x0x68C00x9A8C(43:81,1:30),1,1170);

Ot1995 = reshape(OtVadata.data.x0x65B00x603B0x4EA70x503C(1:39,1:30),1,1170);
Ot1996 = reshape(OtVadata.data.x0x65B00x603B0x4EA70x503C(43:81,1:30),1,1170);

interuse1995 = IO1997data.data.MRIO19970x65B0(1:1209,1171).*(976488835/sum(IO1997data.data.MRIO19970x65B0(1:1209,1171)));
interuse1996 = IO1997data.data.MRIO19970x65B0(1:1209,1171).*(1116591381/sum(IO1997data.data.MRIO19970x65B0(1:1209,1171)));

Valaddnew95 = reshape(OtVadata.data.x0x65B00x589E0x52A00x503C(1:39,1:30),1,1170);
Valaddnew96 = reshape(OtVadata.data.x0x65B00x589E0x52A00x503C(43:81,1:30),1,1170);


interinput2018 = reshape(OtVadata.data.x0x68C00x9A8C(1:39,34:63),1,1170);
interinput2019 = reshape(OtVadata.data.x0x68C00x9A8C(43:81,34:63),1,1170);
interinput2020 = reshape(OtVadata.data.x0x68C00x9A8C(85:123,34:63),1,1170);

Ot2018 = reshape(OtVadata.data.x0x65B00x603B0x4EA70x503C(1:39,34:63),1,1170);
Ot2019 = reshape(OtVadata.data.x0x65B00x603B0x4EA70x503C(43:81,34:63),1,1170);
Ot2020 = reshape(OtVadata.data.x0x65B00x603B0x4EA70x503C(85:123,34:63),1,1170);

interuse2018 = IO2017data.data.MRIO20170x65B0(1:1209,1171).*(17006540538/sum(IO2017data.data.MRIO20170x65B0(1:1209,1171)));
interuse2019 = IO2017data.data.MRIO20170x65B0(1:1209,1171).*(17759629512/sum(IO2017data.data.MRIO20170x65B0(1:1209,1171)));
interuse2020 = IO2017data.data.MRIO20170x65B0(1:1209,1171).*(17833373656/sum(IO2017data.data.MRIO20170x65B0(1:1209,1171)));

Valaddnew18 = reshape(OtVadata.data.x0x65B00x589E0x52A00x503C(1:39,34:63),1,1170);
Valaddnew19 = reshape(OtVadata.data.x0x65B00x589E0x52A00x503C(43:81,34:63),1,1170);
Valaddnew20 = reshape(OtVadata.data.x0x65B00x589E0x52A00x503C(85:123,34:63),1,1170);

[m,n]=size(A2017);
U1=[];
R1=[];
CsIx=A2017.*(Ot2020); 
for i=1:m
    U1(i)=sum(CsIx(i,:));
    R1(i)=interuse2020(i)/U1(i);
end
R1(isnan(R1)) = 0;
R1(isinf(R1)) = 0;
Z_R1=diag(R1)*CsIx;
for j=1:n
    V1(j)=sum(Z_R1(:,j));
    S1(j)=interinput2020(j)/V1(j);
end
S1(isnan(S1)) = 0;
S1(isinf(S1)) = 0;
Z_S1=Z_R1*diag(S1);

count=1;
countin=1;
Erro=[];
Uro=[];
Cro=[];
Ux=[];
Rx=[];
Vx=[];
Sx=[];
Z_Sx=Z_S1;
while count<100
    for i=1:m
        Ux(i)=sum(Z_Sx(i,:));
        Rx(i)=interuse2020(i)/Ux(i);
        if isnan(Rx(i)) || isinf(Rx(i))
            Rx(i) = 0;
        end
        Uro(i)=(Rx(i)-1)^2;
    end
    Zrx=diag(Rx)*Z_Sx;
    for j=1:n
        Vx(j)=sum(Zrx(:,j));
        Sx(j)=interinput2020(j)/Vx(j);
        if isnan(Sx(j)) || isinf(Sx(j))
            Sx(j) = 0;
        end
        Cro(j)=(Sx(j)-1)^2;
    end
    Erro(count)=sum(Uro)+sum(Cro);
    Z_Sx=Zrx*diag(Sx);
    count=count+1;
end
[ExPv,Exp]=min(Erro);
Z_Sx=Z_S1;
t=1;
while t<=Exp
    for i=1:m
        Ux(i)=sum(Z_Sx(i,:));
        Rx(i)=interuse2020(i)/Ux(i);
        if isnan(Rx(i)) || isinf(Rx(i))
            Rx(i) = 0;
        end
    end
    Zrx=diag(Rx)*Z_Sx;
    for j=1:n
        Vx(j)=sum(Zrx(:,j));
        Sx(j)=interinput2020(j)/Vx(j);
        if isnan(Sx(j)) || isinf(Sx(j))
            Sx(j) = 0;
        end
    end
    Z_Sx=Zrx*diag(Sx);
    t=t+1;
end

