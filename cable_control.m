% initialization section

%path(path,'../../hhfiles')

nodes=200;
% arrays for currents from upstream and downstream neighbors
left = zeros(1,nodes);
right = zeros(1,nodes);

% timing parameters

cm=1.0; dt=0.002; tfinal = 8.0; 
stim=200.0;
% stimulate cell with 1msec pulse at 200 uA/cm^2 strength
tstimon = [1];
tstimoff=tstimon+1;

% add monitor for time and step for time monitor
tmonitor = 0.0;
tmonstep = 0.5;

% spacing parameters

% space step 100 micron
dx = 0.01;
% cable radius 20 micron
a= 0.002;
% specific intracellular resistivity(Ri) in kohm-cm
Ri = 0.1;
% mesh ratio
gamma = (a*dt)/(2*Ri*cm*dx*dx);
% Assume total current (I) can be neglected
% Core-conductor resistances
ri = 8; % intracellular core-conductor resistance of 8 Mohm/cm 
re = 2; % extracellular core-conductor resistance of 2 Mohm/cm. 
% ri = Ri/(pi()*a*a);
% re = Re/(pi()*(b-a)*(b-a));

% spacewrite variable for debugging
spacewrite =tfinal; 
% initialize vectors
times = 0:dt:tfinal;
allvm=zeros(length(times),nodes);

% voltage parameters

vrest=-67; tc=27;
vm=vrest*ones(1,nodes); vlast=vrest*ones(1,nodes);
dvdt = zeros(1,nodes);

% sodium parameters

ina=zeros(1,nodes); 
nai=61; nae=485; gnabar=120.; ena=nernst(nai,nae,tc,1);
alfam=zeros(1,nodes); betam=zeros(1,nodes);
for index=1:nodes
    [alfam(index) betam(index)]=mgate(vlast(index),vrest);
end
m0=zeros(1,nodes); 
m0=alfam./(alfam+betam);
    
alfah=zeros(1,nodes); betah=zeros(1,nodes);
for index=1:nodes
    [alfah(index) betah(index)]=hgate(vlast(index),vrest);
end
h0=zeros(1,nodes); 
h0=alfah./(alfah+betah);
    
% potassium parameters

ik=zeros(1,nodes); 

ki=280; ke=10; gkbar=36.; ek=nernst(ki,ke,tc,1);
alfan=zeros(1,nodes); betan=zeros(1,nodes); 
for index=1:nodes
    [alfan(index) betan(index)]=ngate(vlast(index),vrest);
%    [alfan(index) betan(index)]=ngate_bh(vlast(index),vrest);
end
n0=zeros(1,nodes); 
n0=alfan./(alfan+betan);

% chloride parameters

icl=zeros(1,nodes); 
cli=50; cle=400; gclbar=0.3; ecl=nernst(cli,cle,tc,-1);

iion=zeros(1,nodes); ij=zeros(1,nodes);

k=1;
for time=0.0:dt:tfinal
    
    for index=1:nodes
        [alfam(index) betam(index)]=mgate(vm(index),vrest);
        m(index)=gateupdate(alfam(index),betam(index),m0(index),dt); 
        [alfah(index) betah(index)]=hgate(vm(index),vrest);
        h(index)=gateupdate(alfah(index),betah(index),h0(index),dt); 
        [alfan(index) betan(index)]=ngate(vm(index),vrest);
%        [alfan(index) betan(index)]=ngate_bh(vm(index),vrest);
        n(index)=gateupdate(alfan(index),betan(index),n0(index),dt); 

        ina(index)=gnabar*h(index)*m(index)*m(index)*m(index)*(vm(index)-ena);
        ik(index)=gkbar*n(index)*n(index)*n(index)*n(index)*(vm(index)-ek);
        icl(index)=gclbar*(vm(index)-ecl);
    end
    
    m0=m; n0=n; h0=h;
    iion=ina+ik+icl;

    % stimulation at one end
    
    for tindex=1:length(tstimon)
        if((time>tstimon(tindex)) && (time<tstimoff(tindex)))
            iion(1)=iion(1)-stim;
            %iion(200)=iion(200)-stim;
        end
    end
    
    % implement boundary conditions at cable ends
    
    left(1) = 0.0;
    right(1) = gamma*(vlast(2)-vlast(1));
    left(nodes) = gamma*(vlast(nodes-1)-vlast(nodes));
    right(nodes) = 0.0;
    for index=2:nodes-1
        left(index) = gamma*(vlast(index-1)-vlast(index));
        right(index) = gamma*(vlast(index+1)-vlast(index));
    end
    
    % vm update

    vm = vlast + left + right -dt/cm*(iion);
    
    allvm(k,:)=vm;
    k=k+1;
    
    vlast=vm;    
    
    % Please store the value for Vm at every point on the cable on the last
    % time step of the simulation in a file named vm_01.dat.
    if time ==tfinal
        dlmwrite('vm_01.dat',vm);
    end
    
end


% Plot Vm as a function of position before you calculate the
% second derivative of Vm with respect to the x direction (d2Vm/dx2). 
vm = dlmread('vm_01.dat');
position = [0:nodes-1]/100;
%subplot(3,1,1);
%plot(position,vm);title('Vm at each node: Control');
%ylabel('Vm(V)');xlabel('Position (cm)');

% find mesh ratio (tolerance=0.02) 
gamma;

cells = nodes;
% Initialize arrays for 1st and 2nd derivatives of transmembrane potiential over space
dvmdx(cells)=0;
d2vmdx2(cells) = 0;
% Calculate the derivative of transmembrane potiential over space
for i = 1:cells-1
    dvmdx(i) =(vm(i+1)-vm(i))/dx;
end
% plot dVm/dx as a function of position before you answer the following questions.
%subplot(3,1,2);
%figure;plot(position,dvmdx);title('dVm/dx');

% Calculate the second derivative of transmembrane potiential over space
for i = 1:cells-2
    d2vmdx2(i) = (dvmdx(i+1)-dvmdx(i))/dx;
end
%subplot(3,1,3);
%figure;plot(position,d2vmdx2);title('d2Vm/dx2');
% Find the cell peak to end depolarization and start repolarization
[dep_max_vmdx i_dep_max_vm]= max(vm);
[rep_min_vmdx i_rep_min_vm]= min(vm(60:120));
% The peak maximum value for (d2Vm/dx2) in the part of the model that is
% undergoing depolarization. (tolerance=200) 
% depolarization:position 1.1-1.3
[dep_max_d2vmdx2 i_dep_max_d2vmdx2]= max(d2vmdx2(i_dep_max_vm:cells));
dep_max_d2vmdx2;
% The peak minimum value for (d2Vm/dx2) in the part of the model that is
% undergoing depolarization. (tolerance=500) 
[dep_min_d2vmdx2 i_dep_min_d2vmdx2]= min(d2vmdx2(i_dep_max_vm:cells));
dep_min_d2vmdx2;
% The peak maximum value for (d2Vm/dx2) in the part of the model that is
% undergoing repolarization. (tolerance=200) 
% repolarization: .6-1.2
[rep_max_d2vmdx2 i_rep_max_d2vmdx2]= max(d2vmdx2(i_rep_min_vm:i_dep_max_vm));
rep_max_d2vmdx2;
% The peak minimum value for (d2Vm/dx2) in the part of the model that is
% undergoing repolarization. (tolerance=200) 
[rep_min_d2vmdx2 i_rep_min_d2vmdx2]= min(d2vmdx2(i_rep_min_vm:i_dep_max_vm-6));
rep_min_d2vmdx2;

% use the average for the gradient calculated between the cell and its 
% neighbor to the left, and between the cell and its neighbor to the right.

gradR(cells) = 0;
for i = 1:cells-2
    gradR(i) = (vm(i+1)-vm(i))/dx;
end

gradL(cells)=0;
for i = 2:cells-1
    gradL(i-1) = (vm(i)-vm(i-1))/dx;
end

gradient(cells)=0;
for i = 2:cells-1
    gradient(i) = (gradL(i-1)+gradR(i-1))/2;
end
dlmwrite('gradvm_02.dat',gradient);
%figure;plot(position,gradient, position, dvmdx);title('Gradient (dVm/dx)');
%  find the x value at which gradient is at its most negative (peak minimum)
[gradient_min, i_gradient_min]= min(gradient);
x_gradient_min = position(i_gradient_min);

% The ratio between the gradient at a position that is 0.2 mm away from the
% peak minimum site and the gradient at the peak minimum site in the +x
% direction. (tolerance=0.10) 0.2mm = 0.02cm
x = x_gradient_min+0.02; i_x = find(x==position);
gradient_x = gradient(i_x);
ratio = gradient_x/gradient_min;
% The ratio between the gradient at a position that is 0.5 mm away from the
% peak minimum site and the gradient at the peak minimum site in the -x
% direction. (tolerance=0.10)
x = x_gradient_min-0.05; i_x = find(x==position);
gradient_x = gradient(i_x);
ratio = gradient_x/gradient_min;
% The ratio between the gradient at a position that is 2.0 mm away from the
% peak minimum site and the gradient at the peak minimum site in the -x
% direction. (tolerance=0.05)
x = x_gradient_min-0.2; i_x = find(x==position);
gradient_x = gradient(i_x);
ratio = gradient_x/gradient_min;
% The ratio between the gradient at a position that is 5.0 mm away from the
% peak minimum site and the gradient at the peak minimum site in the -x
% direction. (tolerance=0.05)
x = x_gradient_min-0.5; i_x = find(x==position);
gradient_x = gradient(i_x);
ratio = gradient_x/gradient_min;

% 1102 prob3: assume total current (I) can be neglected, an intracellular 
% core-conductor resistance of 8 Mohm/cm and an extracellular 
% core-conductor resistance of 2 Mohm/cm.
ri = 8;re = 2;
% Calculate the intracellular potential gradient 
gradphii = (ri/(ri+re))*dvmdx;
% Save the intracellular potential gradient in a  file (gradphii_03.dat)
dlmwrite('gradphii_03.dat',gradphii);
% plot(position,gradphii,position,dvmdx);title('Intracellular potential gradient');

% What is the ratio for the gradient of intracellular potential to the
% gradient of transmembrane potential at a position that is 0.2 mm away
% from the peak minimum site for the gradient in transmembrane potential in
% the +x direction. (tolerance=0.02)
% gradphii = ri/(ri+re)*dvmdx and same location so gradphii/dvmdx = ri/(ri+re)
ratio = gradphii/dvmdx;
ratio = ri/(ri+re);
% What is the ratio for the gradient of intracellular potential to the
% gradient of transmembrane potential at a position that is 0.5 mm away
% from the peak minimum site for the gradient in transmembrane potential in
% the -x direction. (tolerance=0.02)
ratio = ri/(ri+re);
% What is the ratio for the gradient of intracellular potential to the
% gradient of transmembrane potential at a position that is 2.0 mm away
% from the peak minimum site for the gradient in transmembrane potential in
% the -x direction. (tolerance=0.02)
ri = 8;re = 2;
ratio = ri/(ri+re);
% What is the ratio for the gradient of intracellular potential to the
% gradient of transmembrane potential at a position that is 5.0 mm away
% from the peak minimum site for the gradient in transmembrane potential in
% the -x direction. (tolerance=0.02)
ratio = ri/(ri+re);

% 1102 prob4: assume total current (I) can be neglected, an intracellular 
% core-conductor resistance of 8 Mohm/cm and an extracellular 
% core-conductor resistance of 2 Mohm/cm.
% Calculate the extracellular potential gradient 
ri = 8;re = 2;
gradphie = (-re/(ri+re))*dvmdx;
% Save the intracellular potential gradient in a  file(gradphie_04.dat)
dlmwrite('gradphie_04.dat',gradphie);
plot(position,gradphie,position,gradphii);title('Extracellular potential gradient');
legend('Extracellular potential gradient','Intracellular potential gradient');

% What is the ratio for the gradient of extracellular potential to the
% gradient of transmembrane potential at a position that is 0.2 mm away
% from the peak minimum site for the gradient in transmembrane potential in
% the +x direction.(tolerance=0.02)
% phie = -re/(ri+re)*dvmdx and same location so phie/dvmdx = -re/(ri+re)
ratio = gradphie/dvmdx;
ratio = -re/(ri+re);
% What is the ratio for the gradient of extracellular potential to the
% gradient of transmembrane potential at a position that is 0.5 mm away
% from the peak minimum site for the gradient in transmembrane potential in
% the -x direction. (tolerance=0.02)
ratio = -re/(ri+re);
% What is the ratio for the gradient of extracellular potential to the
% gradient of transmembrane potential at a position that is 2.0 mm away
% from the peak minimum site for the gradient in transmembrane potential in
% the -x direction. (tolerance=0.02)
ratio = -re/(ri+re);
% What is the ratio for the gradient of extracellular potential to the
% gradient of transmembrane potential at a position that is 5.0 mm away
% from the peak minimum site for the gradient in transmembrane potential in
% the -x direction. (tolerance=0.02)
ratio = -re/(ri+re);

% 1102 prob5: use intracellular potential gradient (gradphii03.dat) to find
% the core-conductor transmembrane current (im) with units microA/cm

% Calculate im as a function of position
% Calculate the derivative of gradphii wrt position(dgradphii/dx)
dgradphiidx(cells) = 0;
for i = 1:cells-1
    dgradphiidx(i) = (gradphii(i+1)-gradphii(i))/dx;
end
im = dgradphiidx/(ri*1000); % divide by 1000 to get units microA/cm

% store im in im_05.dat file
dlmwrite('im_05.dat', im);
%plot im over position
plot(position,im);title('transmembrane current (im)');
%  find the x value at which im is at its most negative (peak minimum)
[im_min, i_im_min]= min(im);
x_im_min = position(i_im_min);
% What is the value for im at a position that is 0.2 mm away (+x direction)
% from the peak minimum site for im? (absolute tolerance=0.6)
x = x_im_min+0.02; i_x = find(x==position);
im_x = im(i_x);sprintf('%.4f',im_x);
% Is the current at this site outward (O) ir inward (I)?  I
% What is the value for im at a position that is 0.5 mm away (-x direction)
% from the peak minimum site for im? (absolute tolerance=0.2)
x = x_im_min-0.05; i_x = int8(100*x);
im_x = im(i_x);sprintf('%.4f',im_x);
% Is the current at this site outward (O) ir inward (I)?  I
% What is the value for im at a position that is 2.0 mm away (-x direction)
% from the peak minimum site for im?(absolute tolerance=0.02)
x = x_im_min-0.2; i_x = find(x==position);
im_x = im(i_x);sprintf('%.4f',im_x);
% What is the value for im at a position that is 5.0 mm away (-x direction)
% from the peak minimum site for im?  (absolute tolerance=0.22)
x = x_im_min-0.5; i_x = int8(100*x);
im_x = im(i_x);sprintf('%.4f',im_x);
% Is the current at this site outward (O) ir inward (I)?  O

% 1102 prob 6 use the extracellular potential gradient (gradphie_04.dat) to
% find the core-conductor transmembrane current (im), which has units of
% microA/cm. Please compare the im measured from the extracellular space
% with the im measured from the intracellular space (im_05.dat). Assume an
% extracellular core-conductor resistance of 2 Mohm/cm.
% Calculate im as a function of position
% Calculate the derivative of gradphii wrt position(dgradphii/dx)
dgradphiedx(cells) = 0;
for i = 1:cells-1
    dgradphiedx(i) = (gradphie(i+1)-gradphie(i))/dx;
end
ime = dgradphiedx/(re*1000); % divide by 1000 to get units microA/cm

% store im in the file im_06.dat
dlmwrite('im_06.dat', ime);
%  find the x value at which ime is at its most negative (peak minimum)
[ime_min, i_ime_min]= min(ime);
x_ime_min = position(i_ime_min);
% compare to imi
imi = dlmread('im_05.dat');
[imi_min, i_imi_min]= min(imi);
x_imi_min = position(i_imi_min);
%plot im over position
plot(position,ime,position,imi);title('transmembrane current (im)');
% What is the ratio between im calculated from the extracellular space and
% im calculated from the intracellular space at a position that is 0.2 mm
% away (+x direction) from the peak minimum site for im? (tolerance=0.02)
x = x_ime_min+0.02; i_x = find(x==position);
Ie_x = ime(i_x);%sprintf('%.4f',im_xe);
Ii_x = imi(i_x);%sprintf('%.4f',im_xi);
ratio = Ie_x/Ii_x;
% What is the ratio between im calculated from the extracellular space and
% im calculated from the intracellular space at a position that is 0.5 mm
% away (-x direction) from the peak minimum site for im? (tolerance=0.02)
x = x_ime_min-0.05; i_x = find(x==position);
Ie_x = ime(i_x);%sprintf('%.4f',im_xe);
Ii_x = imi(i_x);%sprintf('%.4f',im_xi);
ratio = Ie_x/Ii_x;
% What is the ratio between im calculated from the extracellular space and
% im calculated from the intracellular space at a position that is 5.0 mm
% away (-x direction) from the peak minimum site for im?  (tolerance=0.02)
xe = x_ime_min+0.5; i_xe = find(xe==position);
Ie_x = ime(i_xe);%sprintf('%.4f',im_xe);
xi = x_imi_min+0.5; i_xi = find(xi==position);
Ii_x = ime(i_xi);%sprintf('%.4f',im_xi);
ratio = Ie_x/Ii_x;

% 1102 prob 7: use the intracellular (gradphii_03.dat) and
% extracellular (gradphie_04.dat) potential gradients to measure the
% intracellular (iintra_07.dat) and extracellular (iextra_07.dat) currents
gradphii = dlmread('gradphii_03.dat');
gradphie = dlmread('gradphie_04.dat');
% Assume an intracellular core conductor resistance of 8 Mohm/cm
% extracellular core-conductor resistance of 2 Mohm/cm.
ri = 8;re = 2; 
% Measure the intracellular (iintra_07.dat) and extracellular (iextra_07.dat) currents
Ii = -gradphii/ri; dlmwrite('iintra_07.dat',Ii);
Ie = -gradphie/re; dlmwrite('iextra_07.dat',Ie);
%  find the x value at which ime is at its most negative (peak minimum)
[Ie_min, i_Ie_min]= min(Ie);
x_Ie_min = position(i_Ie_min);
% What is the ratio between the intracellular and extracellular currents at
% a position that is 0.2 mm away (+x direction) from the peak minimum site
% for the extracellular current?(tolerance=0.02)
x = x_Ie_min+0.02; i_x = find(x==position);
Ie_x = Ie(i_x);%sprintf('%.4f',im_xe);
Ii_x = Ii(i_x);%sprintf('%.4f',im_xi);
ratio = Ie_x/Ii_x
% What is the ratio between the intracellular and extracellular currents at
% a position that is 0.5 mm away (-x direction) from the peak minimum site
% for the extracellular current?
%  (tolerance=0.02)
% 
% What is the ratio between the intracellular and extracellular currents at
% a position that is 2.0 mm away (-x direction) from the peak minimum site
% for the extracellular current?
%  (tolerance=0.02)
% 
% What is the ratio between the intracellular and extracellular currents at
% a position that is 5.0 mm away (-x direction) from the peak minimum site
% for the extracellular current?
%  (tolerance=0.02)
