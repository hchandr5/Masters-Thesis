File=xlsread('ACTIVISg2000_load_time_series_MVAR.csv');   % One day load data
Input=xlsread('System_Mapping1.xlsx'); 
num=1;
RDM=[1:366];

% TD= 'TrainingDays.xlsx'
% xlswrite(TD,RDM);


hour = 15;
for hour=15:15
    
    hour
    
P118=[];
Pval118=[];
Ptest118=[];

for num=1:99
    num
ent=Input(num,5);
Time = File(:,1);
[m1,n1]= size(Time);

Tindex=0;
Temp=[];
PData=[];
cnt=0;
for i=1:m1
    k= rem(i,24);
    if  (k==0)
        k=24;
    end
    if (Time(i,1)==0)
        cnt=cnt+1;
    end
        Pdata(cnt,k)= File(i,ent+4);
end
Temp(k)=Tindex;
Tindex=Tindex+(1/24);


D=Pdata;
D=[];    
ind=hour;
%while(i<25)
    D(:,1)=Pdata(:,ind);
% x=17.5:0.01:22; 
DR=[];
drind=0;
[d1,d2]=size(D);

for i=1:length(RDM)
        drind=drind+1;
        DR(drind,1)=D(RDM(i));
end

D=DR;

[a,a2]=min(D);
[b,b2]=max(D);
% 
x= (0.9*a):0.01:(b*1.1);

[d1,d2] = size(D);
D1=sort(D);


alfa = 0.05;
eps = sqrt((1/(2*d1))*log(2/0.05));
pdf_each_hour = (1/d1)*ones(1,d1);
LB=[];
UB=[];
for i = 1:d1
    cdf_each_hour(1,i) = sum(pdf_each_hour(1:i));
    LB(1,i) = max (cdf_each_hour(1,i)-eps, 0);
    UB(1,i) = min (cdf_each_hour(1,i)+eps, 1);
end

pd5=fitdist(D, 'Kernel');
y5=pdf(pd5,x);
pd10=fitdist(D, 'Kernel');
y10=pdf(pd10,x);
bd= pd5.BandWidth;
[h,p5] = adtest(D,'Distribution',pd5);

viol=0;

if(p5 > 0.05)
while ((p5 > 0.05) && (viol <= (0.05*d1)))
    pd5old=pd5;
    pd10old=pd10;
    y5old=y5;
    y10old=y10;
    bdold=bd;
    violold=viol;
    p5old=p5;

    bd=bd+0.05;
    pd5=fitdist(D, 'Kernel', 'Bandwidth',bd);
    [h,p5] = adtest(D,'Distribution',pd5);

    
y5=pdf(pd5,x);

pd10=fitdist(D, 'Kernel','Bandwidth',bd);
y10=cdf(pd10,x);




id=1;
viol=0;
i=1;
l=length(x);
x(l+1) = 10*b;
l1=length(D1);
D1(l1+1)=8*b;
while(i<= l)
    if(abs(x(i)-D1(id)) <0.005)
         diff1= y10(i)-LB(id);
         diff2 = UB(id)-y10(i);
         if ((diff1 <=0) || (diff2 <=0))
             viol=viol+1;
         end
         id=id+1;
         i=i-1;
    end
    i=i+1;
end
 
x(l+1)=[];
D1(l1+1)=[];

end

pd5=pd5old;
pd10=pd10old;
y5=y5old;
y10=y10old;
bd=bdold;

else
    
   while ((p5 < 0.05) && (bd >0) && (viol <= (0.05*d1)))
    pd5old=pd5;
    pd10old=pd10;
    y5old=y5;
    y10old=y10;
    bdold=bd;
    violold=viol;
    p5old=p5;

    bd=bd-0.05;
    if(bd <0)
        bd = bd + 0.05 - 0.005;
    end
    pd5=fitdist(D, 'Kernel', 'Bandwidth',bd);
    [h,p5] = adtest(D,'Distribution',pd5);

    
y5=pdf(pd5,x);

pd10=fitdist(D, 'Kernel','Bandwidth',bd);
y10=cdf(pd10,x);




id=1;
viol=0;
i=1;
l=length(x);
x(l+1) = 10*b;
l1=length(D1);
D1(l1+1)=8*b;
while(i<= l)
    if(abs(x(i)-D1(id)) <0.005)
         diff1= y10(i)-LB(id);
         diff2 = UB(id)-y10(i);
         if ((diff1 <=0) || (diff2 <=0))
             viol=viol+1;
         end
         id=id+1;
         i=i-1;
    end
    i=i+1;
end
 
x(l+1)=[];
D1(l1+1)=[];

end
    
end
%viol=violold;
%p5=p5old;



% figure(3)
% plot(x, y10,'b', D1,LB,'--', D1, UB, '--')
% 
% figure(4)
% plot(x, y5)
N=100;
%randomNumbers = ksdensity(D, rand(N,1), pd1, 'icdf');
mpc = case118;
busno = Input(num,2);
if(mpc.bus(busno,3)~=0)
    if(mpc.bus(busno,4)==0)
        P118(:,num) = zeros(85000,1);
        Pval118(:,num) = zeros(27000,1);
        Ptest118(:,num) = zeros(1000,1);
    else
       

R=[];
for i=1:175000
    R(i)= random(pd5);
end



avg = mean(R);



R1=[];
R1 = transpose(R) - ((avg - Input(num,4)*ones(175000,1)));

for i=1:175000
    if(R1(i,1) <= 0)
        while(R1(i,1) <= 0)
            R(i) = random(pd5);
            avg = mean(R);
            R1 = [];
            R1 = transpose(R) - ((avg - Input(num,4)*ones(175000,1)));
        end
    end
end

Rf = [];
cnt=1;
m=1;

while(cnt<=85000) && (m <= 175000)
    if(R1(m,1) >=0)
        Rf(cnt,1) = R1(m,1);
        cnt = cnt+1;
    end
    m = m+1;
end

P118(:,num) = Rf(:,1);

RVal=[];
for i=1:48000
    Rval(i)= random(pd5);
end

avg = mean(Rval);



R1val=[];
R1val = transpose(Rval) - ((avg - Input(num,4)*ones(48000,1)));

for i=1:48000
    if(R1val(i,1) <= 0)
        while(R1val(i,1) <= 0)
            Rval(i) = random(pd5);
            avg = mean(Rval);
            R1val = [];
            R1val = transpose(Rval) - ((avg - Input(num,4)*ones(48000,1)));
        end
    end
end

Rfval = [];
cnt=1;
m=1;
while(cnt<=27000) && (m <= 48000)
    if(R1val(m,1) >=0)
        Rfval(cnt,1) = R1val(m,1);
        cnt = cnt +1;
    end
    m = m+1;
end

Pval118(:,num) = Rfval(:,1);

RTest=[];
for i=1:2000
    Rtest(i)= random(pd5);
end

avg = mean(Rtest);



R1test=[];
R1test = transpose(Rtest) - ((avg - Input(num,4)*ones(2000,1)));

for i=1:2000
    if(R1test(i,1) <= 0)
        while(R1test(i,1) <= 0)
            Rtest(i) = random(pd5);
            avg = mean(Rtest);
            R1test = [];
            R1test = transpose(Rtest) - ((avg - Input(num,4)*ones(2000,1)));
        end
    end
end

Rftest = [];
cnt=1;
m=1;

while(cnt<=1000) && (m <= 2000)
    if(R1test(m,1) >=0)
        Rftest(cnt,1) = R1test(m,1);
        cnt = cnt + 1;
    end
    m = m+1;
end

Ptest118(:,num) = Rftest(:,1);

    end
end
end

filename= sprintf('LQTrainDat%d.csv', hour)
csvwrite(filename,P118);

filename1= sprintf('LQValDat%d.csv', hour)
csvwrite(filename1,Pval118);

filename2 = sprintf('LQTestDat%d.csv', hour)
csvwrite(filename2,Ptest118);

end
% R1=[];
% R1 = D*2;
% %R1=transpose(R1);
% 
% 
% pd11= fitdist(R1, 'Kernel', 'Bandwidth', bd);
% pd11= fitdist(R1, 'Lognormal');
% 
% pd12= fitdist(R1, 'Kernel', 'Bandwidth', bd);
% y12=cdf(pd11,x2);
% 
% x2 = 2*a: 0.01*2: b*2;
% 
% y11=pdf(pd11,x2);
% figure(4)
% plot(x,y5,'b',x2, y11,'r')
% 
% [h,p11] = adtest(R1,'Distribution',pd11);
% display(p11)
