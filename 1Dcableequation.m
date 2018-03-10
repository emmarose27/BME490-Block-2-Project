%path(path,'../../hhfiles')
% filePtr=fopen('hh.dat','w');
% 
nodes=200;
% filePtr1=fopen('node1_03.dat','w');
% filePtr2=fopen('node2_03.dat','w');

cm=1.0; dt=0.002; tfinal=20.;
stim=200.0; 
tstimon=[1];
tstimoff=tstimon+1;
dx=0.01;
a=0.002;
Ri=0.1;
p=0:dt:tfinal;
allvm=zeros(length(p),nodes);
% voltage parameters

vrest=-67; tc=27;
vm=vrest*ones(1,nodes); vlast=vrest*ones(1,nodes);
ina=zeros(1,nodes); 
nai=61; nae=485; gnabar=120.; ena=nernsthw0831(nai,nae,tc,1);
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
ki=280; ke=10; gkbar=36.; ek=nernsthw0831(ki,ke,tc,1);
alfan=zeros(1,nodes); betan=zeros(1,nodes); 
for index=1:nodes
    [alfan(index) betan(index)]=ngate(vlast(index),vrest);
end
n0=zeros(1,nodes); 
n0=alfan./(alfan+betan);

% chloride parameters

icl=zeros(1,nodes); 
cli=50; cle=400; gclbar=0.3; ecl=nernsthw0831(cli,cle,tc,-1);

iion=zeros(1,nodes); % ij=zeros(1,nodes);

% tmonitor=0.0
% tmonstep=dt
gamma=(a*dt)/(2*Ri*cm*dx^2);
k=1;
for time=0.0:dt:tfinal
%     if (time==tmonitor)
%         time
%         tmonitor=tmonitor+tmonstep;
%     end
%     disp(time);
    %q=find(p==time);
    
    % ionic currents

    for index=1:nodes
        [alfam(index) betam(index)]=mgate(vm(index),vrest);
        m(index)=gateupdate(alfam(index),betam(index),m0(index),dt); 
        [alfah(index) betah(index)]=hgate(vm(index),vrest);
        h(index)=gateupdate(alfah(index),betah(index),h0(index),dt); 
        [alfan(index) betan(index)]=ngate(vm(index),vrest);
        n(index)=gateupdate(alfan(index),betan(index),n0(index),dt); 

        ina(index)=gnabar*h(index)*m(index)*m(index)*m(index)*(vm(index)-ena);
        ik(index)=gkbar*n(index)*n(index)*n(index)*n(index)*(vm(index)-ek);
        icl(index)=gclbar*(vm(index)-ecl);
    end
    
    m0=m; n0=n; h0=h;
    iion=ina+ik+icl;
    
    for tindex=1:length(tstimon)
        if((time>tstimon(tindex)) && (time<tstimoff(tindex)))
           iion(1)=iion(1)-stim;
%            iion(200)=iion(200)-stim;
        end
    end
    
    % junctional current
    left(1)=0;
    right(1)=gamma*(vlast(2)-vlast(1));
    left(nodes)=gamma*(vlast(nodes-1)-vlast(nodes));
    right(nodes)=0;
    for ind=2:nodes-1;
        left(ind)=gamma*(vlast(ind-1)-vlast(ind));
        right(ind)=gamma*(vlast(ind+1)-vlast(ind));
    end
    
    % vm update
    vm=vlast+left+right-dt/cm*(iion);
    
    
    allvm(k,:)=vm;
    k=k+1;
%       fprintf(filePtr1,'%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\n',time,vm(1),vlast(1),iion(1),ina(1),ik(1),icl(1),h(1),m(1),n(1));
%       fprintf(filePtr2,'%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\n',time,vm(2),vlast(2),iion(2),ina(2),ik(2),icl(2),h(2),m(2),n(2));
    vlast=vm;
    
end

allvm=allvm';
time=p;
alldvmdt=zeros(nodes,length(p));
for i=1:nodes;
    for j=2:length(time);
        alldvmdt(i,j-1)=(allvm(i,j)-allvm(i,j-1))/(time(j)-time(j-1));
    end
end

alldvmdt=alldvmdt';
a=max(alldvmdt);
% aoriginal=max(alldvmdt);
% plot(aoriginal);
% hold on;
dvmdtmax=max(max(alldvmdt))
dvmdtmin=min(max(alldvmdt))

% b=linspace(aoriginal(14),aoriginal(96),96-14);
% a(15:96)=b(1:82)
meandvmdt=mean(a)
plot(a);
% plot a is the plot of all of the maxdvmdts

% lastdvdm14=(alldvmdt(:,14));
% lastdvdm14=lastdvdm14';
% at1=time(find(lastdvdm14==a(1,14)));
% lastdvdm100=(alldvmdt(:,100));
% lastdvdm100=lastdvdm100';
% at2=time(find(lastdvdm100==a(1,100)));
% velocity=(1.8-0.14)*0.01/((at2-at1)*0.001) %%% how to find this?


%% 
for i=1:nodes

    for j=1:length(time)-1

        alldvmdt(i,j)=(allvm(i,j+1)-allvm(i,j))/dt;

    end

    [maxDVm(i),maxDVtime] = max(alldvmdt(i,:));

    DVtime(i) = time(maxDVtime);

end

 

position = 0:dx:(dx*nodes)-dx;

 

[dvmdtmax,indd] = max(a);

dvmdtmin = min(a);

meandvmdt = mean(a);

 

P = position(120)-position(80);

T = DVtime(120)-DVtime(80);

 

velocity = (P*10)/T;
