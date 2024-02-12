clear all;
close all;
clc;
tic;
L=10000;
D=-1000;
Aa=.001;
Ss=.0001;
Kx=5; % m/day
Kz=1;  % m/day
Ani=Kx/Kz;
mloop=10;
nloop=10;
ploop=mloop;
qloop=nloop;
HL=50;
HR=50;
mm=7;
ww=pi/365;
t=0.25*365; % days
alp=sin(ww*t)^2;
[x,z]=meshgrid(0:200:L,0:-100:D);

for m=1:mloop
    for n=1:nloop

        Nm=(m-1)*pi/L;
        f1=@(x)(cos(Nm*x).*cos(Nm*x));
        Bm=(integral(f1,0,L))^.5;
        Nn=(1-2*n)*pi/(2*D);
        f2=@(z)(sin(Nn*z).*sin(Nn*z));
        Bn=(integral(f2,0,D))^.5;

        for p=1:ploop
            for q=1:qloop

                Np=(p-1)*pi/L;
                f3=@(x)(cos(Np*x).*cos(Np*x));
                Bp=(integral(f3,0,L))^.5;
                Nq=(1-2*q)*pi/(2*D);
                f4=@(z)(sin(Nq*z).*sin(Nq*z));
                Bq=(integral(f4,0,D))^.5;

                % Matrix for time-derivative term
                f5=@(x)(cos(Nm*x).*cos(Np*x));
                w1=integral(f5,0,L);
                f6=@(z)(sin(Nn*z).*sin(Nq*z));
                w2=integral(f6,0,D);
                aa=w1/(Bm*Bp);
                ab=w2/(Bn*Bq);
                A((m-1)*nloop+n,(p-1)*qloop+q)=aa*ab;

                % Matrix for second order x-derivative term
                f5=@(x)(cos(Nm*x).*cos(Np*x));
                w1=integral(f5,0,L);
                f6=@(z)(sin(Nn*z).*sin(Nq*z).*exp(Aa*z));
                w2=integral(f6,0,D);
                aa=w1/(Bm*Bp);
                ab=w2/(Bn*Bq);
                B((m-1)*nloop+n,(p-1)*qloop+q)=Kx*aa*ab*Np^2;

                % Matrix for second order z-derivative term
                f5=@(x)(cos(Nm*x).*cos(Np*x));
                w1=integral(f5,0,L);
                f6=@(z)(sin(Nn*z).*sin(Nq*z).*exp(Aa*z));
                w2=integral(f6,0,D);
                aa=w1/(Bm*Bp);
                ab=w2/(Bn*Bq);
                C((m-1)*nloop+n,(p-1)*qloop+q)=Kz*aa*ab*Nq^2;

                % Matrix for first order z-derivative term
                f5=@(x)(cos(Nm*x).*cos(Np*x));
                w1=integral(f5,0,L);
                f6=@(z)(sin(Nn*z).*cos(Nq*z).*exp(Aa*z));
                w2=integral(f6,0,D);
                aa=w1/(Bm*Bp);
                ab=w2/(Bn*Bq);
                E((m-1)*nloop+n,(p-1)*qloop+q)=Kz*Aa*aa*ab*Nq;
            end
        end
    end
end

F=(B+C-E)/Ss;

aR=((Aa/2)^2+(Ani*(pi/L)^2))^.5;
aL=((Aa/2)^2+(Ani*(mm*pi/L)^2))^.5;

aa=HR*cos(pi*x/L).*exp(-Aa*z/2);
ab=(Aa*sinh(aR*(z-D)))+(2*aR*cosh(aR*(z-D)));
ac=(Aa*sinh(aR*D))-(2*aR*cosh(aR*D));
ad=aa.*(ab/ac);
ae=HL*alp*cos(mm*pi*x/L).*exp(-Aa*z/2);
af=(Aa*sinh(aL*(z-D)))+(2*aL*cosh(aL*(z-D)));
ag=(Aa*sinh(aL*D))-(2*aL*cosh(aL*D));
ah=ae.*(af/ag);
hst=(HR+(HL*alp))+ad+ah;

%equating to g(m,n)
for m=1:mloop
    for n=1:nloop

        Nm=(m-1)*pi/L;
        f1=@(x)(cos(Nm*x).*cos(Nm*x));
        Bm=(integral(f1,0,L))^.5;
        Nn=(1-2*n)*pi/(2*D);
        f2=@(z)(sin(Nn*z).*sin(Nn*z));
        Bn=(integral(f2,0,D))^.5;

        f3=@(x)(cos(Nm*x));
        w3=integral(f3,0,L);
        ba=w3/Bm;
        f4=@(z)(sin(Nn*z));
        w4=integral(f4,0,D);
        bb=w4/Bn;
        bc=ba*bb;

        f7=@(x)(cos(Nm*x).*cos(mm*pi*x/L));
        w7=integral(f7,0,L);
        bg=w7/Bm;
        f8=@(z)(sin(Nn*z).*exp(-Aa*z/2).*((Aa*sinh(aL*(z-D)))+(2*aL*cosh(aL*(z-D)))));
        w8=integral(f8,0,D);
        bh=w8/(Bn*ag);
        bi=bg*bh;
        gmn((m-1)*nloop+n,1)=-ww*(bc+bi)*HL;
    end
end

da=((2*ww)*eye(nloop*mloop))^2+F^2;
db=2*ww*eye(nloop*mloop)*expm(-F*t);
dc=F*sin(2*ww*t);
dd=2*ww*eye(nloop*mloop)*cos(2*ww*t);
de=db+dc-dd;
df=(da^-1)*de;
Tmn=df*gmn;


% transient solution in terms of U
st=0;
for m=1:mloop
    for n=1:nloop

        Nm=(m-1)*pi/L;
        f1=@(x)(cos(Nm*x).*cos(Nm*x));
        Bm=(integral(f1,0,L))^.5;
        Nn=(1-2*n)*pi/(2*D);
        f2=@(z)(sin(Nn*z).*sin(Nn*z));
        Bn=(integral(f2,0,D))^.5;

        ca=cos(Nm*x).*sin(Nn*z)/(Bm*Bn);
        cb=Tmn((m-1)*nloop+n,1)*ca;
        st=st+cb;
    end
end

hh=hst+st;

disp(hh)

toc;        