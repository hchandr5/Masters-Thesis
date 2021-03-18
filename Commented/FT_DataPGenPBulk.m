tic
File=xlsread('ACTIVISg2000_load_time_series_MW.csv');   % One day real load data
Input=xlsread('System_Mapping1.xlsx');  % File with mapping information between the IEEE 118 bus system and Texas system
num=1;

% Code to calculate the number of loads begins

mpc = case118;
[mbus,nbus] = size(mpc.bus)
cnt_load = 0;
for i=1:mbus
if(mpc.bus(i,3)~=0)
    cnt_load = cnt_load+1; %cnt_load represents the number of loads
end
end

% Ends

prompt = 'Enter the hour number (From):';
hour_num_from = input(prompt);

prompt = 'Enter the hour number (To):';
hour_num_to = input(prompt);

for hour = hour_num_from:hour_num_to
    hour
    
P118=[];
Pval118=[];
Ptest118=[];

for num=1:cnt_load
    num
ent=Input(num,5); % 'ent' represents the Texas system load number for which the distribution will be obtained
Time = File(:,1);
[m1,n1]= size(Time);

% Code to obtain the required load data from the input file and store it
% variable Pdata - Begins

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
        Pdata(cnt,k)= File(i,ent+4); % The array Pdata is used to store the data in the 366*24 format for load 'ent'
end
Temp(k)=Tindex;
Tindex=Tindex+(1/24);

% Code to obtain the required load data from the input file and store it
% variable Pdata - Ends

D=Pdata;
D=[];    
ind=hour;
    D(:,1)=Pdata(:,ind); % Vector D stores the load data for a particular hour

[d1,d2]=size(D);

[a,a2]=min(D);
[b,b2]=max(D);
 
x= (0.9*a):0.01:(b*1.1); % Range of values for the distribution

[d1,d2] = size(D);
D1=sort(D); % D1 rdenotes the sorted load array

% Code for obtaining lower and upper bounds of DKW Inequality test - Begins

alfa = 0.05;
eps = sqrt((1/(2*d1))*log(2/0.05));
pdf_each_hour = (1/d1)*ones(1,d1);
LB=[];
UB=[];
for i = 1:d1
    cdf_each_hour(1,i) = sum(pdf_each_hour(1:i));
    LB(1,i) = max (cdf_each_hour(1,i)-eps, 0); % LB denotes the array to store lower bound values for DKW inequality test
    UB(1,i) = min (cdf_each_hour(1,i)+eps, 1); % UB denotes the array to store upper bound values for DKW inequality test
end

% Code for obtaining lower and upper bounds of DKW Inequality test - Ends

pd5=fitdist(D, 'Kernel');
y5=pdf(pd5,x); % Obtaining the pdf using a Kernel Density Estimator
bd= pd5.BandWidth; % 'bd' initially stores the default Matlab bandwidth value
[h,p5] = adtest(D,'Distribution',pd5); % 'p5' denotes the p value of the AD test
pd10=fitdist(D, 'Kernel');
y10=cdf(pd10,x);
viol=0;

% Code for obtaining the distribution by increasing bandwidth  - Begins
if(p5 > 0.05)
while ((p5 > 0.05) && (viol <= (0.05*d1)))
    % Storing the values of variables in current iteration before varying the bandwidth
    pd5old=pd5;
    y5old=y5;
    bdold=bd;
    violold=viol;
    p5old=p5;
    pd10old=pd10;
    y10old=y10;

    bd=bd+0.05; % Modifying badnwidth
    pd5=fitdist(D, 'Kernel', 'Bandwidth',bd);
    [h,p5] = adtest(D,'Distribution',pd5);

    
y5=pdf(pd5,x); % Obtaining new pdf with modified bandwidth

pd10=fitdist(D, 'Kernel','Bandwidth',bd);

y10=cdf(pd10,x); % Obtaining new cdf with modified bandwidth


% Code for checking DKW inequality- Begins

id=1;
viol=0;
i=1;
l=length(x);
x(l+1) = 10*b;
l1=length(D1);
D1(l1+1)=8*b; % Dummy value stored for the code to run
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

% Code for checking DKW inequality- Ends
 
x(l+1)=[];
D1(l1+1)=[];

end

% The final values of variables will be the previous iteration values
% before it quits the while loop
pd5=pd5old;
pd10=pd10old;
y5=y5old;
y10=y10old;
bd=bdold;

% Code for obtaining the distribution by increasing bandwidth  - Ends

% Code for obtaining the distribution by decreasing bandwidth - Begins
else
    while ((p5 < 0.05) && (viol <= (0.05*d1)))
    pd5old=pd5;
    pd10old=pd10;
    y5old=y5;
    y10old=y10;
    bdold=bd;
    violold=viol;
    p5old=p5;

    bd=bd-0.05;
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
   
% Code for obtaining the distribution by decreasing bandwidth - Ends
end


% Code for sampling discrete load data from the distribution - Begins

R=[];
for i=1:175000
    R(i)= random(pd5); % Load Points for training data set
end




avg = mean(R);



R1=[];
R1 = transpose(R) - ((avg - Input(num,3)*ones(175000,1))); % Load Shifting (The shifted load value could be negative)

for i=1:175000  % A loop to partially remove negative values after shifting
    if(R1(i,1) <= 0)
        while(R1(i,1) <= 0)
            R(i) = random(pd5);
            avg = mean(R);
            R1 = [];
            R1 = transpose(R) - ((avg - Input(num,3)*ones(175000,1))); % Load Shifting after recalculation of average
        end
    end
end % There could still be negative values of load even after this loop

Rf = [];
cnt_train=1; % 'cnt_train' represents the variable for number of training samples
m=1;

while((cnt_train<=85000) && (m <= 175000)) 
    if(R1(m,1) >=0)
        Rf(cnt_train,1) = R1(m,1); % 'Rf' represents the array to store training load samples
        cnt_train = cnt_train+1;
    end
    m = m+1;
end

P118(:,num) = Rf(:,1); % 'num' represents the load number for IEEE 118 bus system

RVal=[];
for i=1:48000
    Rval(i)= random(pd5); % Load Points for validation data set
end

avg = mean(Rval);



R1val=[];
R1val = transpose(Rval) - ((avg - Input(num,3)*ones(48000,1))); % Load Shifting (The shifted load value could be negative)
 
for i=1:48000 % A loop to partially remove negative values after shifting
    if(R1val(i,1) <= 0)
        while(R1val(i,1) >= 0)
            Rval(i) = random(pd5);
            avg = mean(Rval);
            R1val = [];
            R1val = transpose(Rval) - ((avg - Input(num,3)*ones(48000,1))); % Load Shifting  after recalculation of average
        end
    end
end % There could still be negative values of load even after this loop

Rfval = [];
cnt_val=1;  % 'cnt_val' represents the variable for number of validation samples
m=1; 
while((cnt_val<=27000) && (m <= 48000))
    if(R1val(m,1) >=0)
        Rfval(cnt,1) = R1val(m,1); % 'Rfval' represents the array to store validation load samples
        cnt_val = cnt_val +1;
    end
    m = m+1;
end 

Pval118(:,num) = Rfval(:,1); % 'num' represents the load number for IEEE 118 bus system

RTest=[];
for i=1:2100
    Rtest(i)= random(pd5); % Load Points for testing data set
end

avg = mean(Rtest);



R1test=[];
R1test = transpose(Rtest) - ((avg - Input(num,3)*ones(2100,1))); % Load Shifting (The shifted load value could be negative)

for i=1:2100 % A loop to partially remove negative values after shifting
    if(R1test(i,1) <= 0)
        while(R1test(i,1) >= 0)
            Rtest(i) = random(pd5);
            avg = mean(Rtest);
            R1test = [];
            R1test = transpose(Rtest) - ((avg - Input(num,3)*ones(2100,1))); % Load Shifting  after recalculation of average
        end
    end
end  % There could still be negative values of load even after this loop

Rftest = [];
cnt_test=1; % 'cnt_test' represents the variable for number of testing samples
m=1;

while((cnt_test<=1000) && (m <= 2100))
    if(R1test(m,1) >=0)
        Rftest(cnt_test,1) = R1test(m,1); % 'Rftest' represents the array to store testing load samples
        cnt_test = cnt_test + 1;
    end
    m = m+1;
end

Ptest118(:,num) = Rftest(:,1); % 'num' represents the load number for IEEE 118 bus system

end

% Code for sampling discrete load data from the distribution - Ends

% Code for storing the load data - Begins

filename= sprintf('LPTrainDat%d.csv', hour)
xlswrite(filename,P118);

filename1= sprintf('LPValDat%d.csv', hour)
xlswrite(filename1,Pval118);

filename2 = sprintf('LPTestDat%d.csv', hour)
xlswrite(filename2,Ptest118);

% Code for storing the load data - Ends
end

time =toc;