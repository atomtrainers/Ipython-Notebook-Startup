import matplotlib as plt;
from scipy.interpolate import interp1d

def GriffAnalysis(fname, Fsamp = 1000.,f11=10.,f22=50.):
	
	fpath= 'D:\\Users\\thadlabs\\Documents\\Chambo LabVIEW Data\\DATA\\GriffithsData\\';
	file=fpath+fname;
	data = transpose(loadtxt(file));
	GriffNum=fname.strip('ShotNoise.dat')
	
	diffTransImp = 20.E3;
	sumTransImp = 1.E3;
	
	diffData = data[0]/diffTransImp;
	sumData = data[1]/sumTransImp;
	sz = size(sumData);
	
	IavgSum = mean(sumData);
	IavgSumuA=round(IavgSum*1.E6)
	
	IavgDiff = mean(diffData);
	
	angle=round(1E6*IavgDiff/(2*IavgSum)); #angle in micro-radians
	
	sP,freq = mlab.psd(sumData, NFFT = sz,Fs = Fsamp);  #A**2/Hz
	sumPSD = sqrt(sP);
	
	dP,freq = mlab.psd(diffData, NFFT = sz,Fs = Fsamp);  #A/rHz
	diffPSD = sqrt(dP);
	
	floor=mean(diffPSD[2*f11:2*f22]);
	
	#CMR = -20*log10(diffPSD[2*cmrf]/sumPSD[2*cmrf]);
	#return[CMR]
	
	npoints = size(freq)
	#I0 = <Isum>
	PSN = sqrt(2*1.6E-19*IavgSum) + zeros(npoints);
	
	figure(1,figsize=[10,8]);
	tick_params(labelsize=18)
	loglog(freq,sumPSD,'k-',label = 'Sum');
	loglog(freq,diffPSD,'r-',label = 'Difference');
	loglog(freq,PSN,'g-', linewidth = 3, label = 'PSN');
	legend(fontsize = 18);
	grid('on');
	xlabel('freq,Hz',fontsize=18);
	ylabel('A/rHz',fontsize=18);
	suptitle('Griffiths Circuit, Run %s, $\phi$ = %.0f $\mu$rad, $I_{sum}$ = %.0f $\mu$A' %(GriffNum,angle,IavgSumuA),fontsize=18)
	show();
	
	return[floor]
	
def SimpSubt(fnameSum='ShotNoise0000.dat',fnameDiff='ShotNoise0000.dat',senseSum=500E-6,senseDiff=5E-6,Fsamp = 1000.,f11 = 10., f22 = 50.):
	fpath= 'D:\\Users\\thadlabs\\Documents\\Chambo LabVIEW Data\\DATA\\GriffithsData\\';
	fileSum=fpath+fnameSum;
	fileDiff=fpath+fnameDiff;
	
	dataSum = transpose(loadtxt(fileSum));
	dataDiff = transpose(loadtxt(fileDiff));
	
	SumNum=fnameSum.strip('ShotNoise.dat')
	DiffNum=fnameDiff.strip('ShotNoise.dat')
	
	diffTransImp = 1/senseDiff;
	sumTransImp = 1/senseSum;
	
	sumData = dataSum[0]/sumTransImp;
	sz = size(sumData);
	
	IavgSum = mean(sumData);
	IavgSumuA=round(IavgSum*1.E6)
	
	diffData = dataDiff[0]/diffTransImp;
	IavgDiff = mean(diffData[0]);
	
	angle=round(1E6*IavgDiff/(2*IavgSum)); #Angle in micro-radians

	sP,freq = mlab.psd(sumData, NFFT = sz,Fs = Fsamp);  #A**2/Hz
	sumPSD = sqrt(sP);                                  #A/rHz
	
	dP,freq = mlab.psd(diffData, NFFT = sz,Fs = Fsamp);  
	diffPSD = sqrt(dP);
	
	floor=mean(diffPSD[2*f11:2*f22]);
	
	npoints = size(freq)
	#I0 = <Isum>
	PSN = sqrt(2*1.6E-19*IavgSum) + zeros(npoints);
	
	
	figure(1,figsize=[10,8]);
	tick_params(labelsize=18)
	
	loglog(freq,sumPSD,'k-',label = 'Sum, %s'%SumNum);
	loglog(freq,diffPSD,'r-',label = 'Difference, %s'%DiffNum);
	loglog(freq,PSN,'g-', linewidth = 3, label = 'PSN');
	legend(fontsize = 18);
	grid('on');
	xlabel('freq,Hz',fontsize=18);
	ylabel('A/rHz',fontsize=18);
	suptitle('Passive Subtraction, $\phi$ = %.0f $\mu$rad, $I_{sum}$ = %.0f $\mu$A' %(angle,IavgSumuA),fontsize=18)
	show();
	
def SimpSubt2016(loglogplot=1,fnameSum='ShotNoise0000.dat',fnameDiff='ShotNoise0000.dat',senseSum=500E-6,senseDiff=5E-6,Fsamp = 1000.,fmax=200):
	fpath= 'D:\\Users\\thadlabs\\Documents\\Chambo LabVIEW Data\\DATA\\GriffithsData\\';
	fileSum=fpath+fnameSum;
	fileDiff=fpath+fnameDiff;
	
	dataSum = transpose(loadtxt(fileSum));
	dataDiff = transpose(loadtxt(fileDiff));
	
	SumNum=fnameSum.strip('ShotNoise.dat')
	DiffNum=fnameDiff.strip('ShotNoise.dat')
	
	diffTransImp = 1/senseDiff;
	sumTransImp = 1/senseSum;
	
	sumData = dataSum[0]/sumTransImp;
	sz = size(sumData);
	
	IavgSum = mean(sumData);
	IavgSumuA=round(IavgSum*1.E6)
	
	diffData = dataDiff[0]/diffTransImp;
	IavgDiff = mean(diffData[0]);
	
	angle=round(1E6*IavgDiff/(2*IavgSum)); #Angle in micro-radians

	sP,freq = mlab.psd(sumData, NFFT = sz,Fs = Fsamp);  #A**2/Hz
	sumPSD = sqrt(sP);                                  #A/rHz
	
	dP,freq = mlab.psd(diffData, NFFT = sz,Fs = Fsamp); #A**2/Hz
	diffPSD = sqrt(dP);									#A/rHz
	
	npoints = size(freq)
	PSN = sqrt(2*1.6E-19*IavgSum) + zeros(npoints);
	
	figure(1,figsize=[10,8]);
	tick_params(labelsize=18)
	if loglogplot==1:
		loglog(freq,sumPSD,'k-',label = 'Sum, %s'%SumNum);
		loglog(freq,diffPSD,'r-',label = 'Difference, %s'%DiffNum);
		loglog(freq,PSN,'g-', linewidth = 3, label = 'PSN');
	else:
		semilogy(freq,sumPSD,'k-',label = 'Sum, %s'%SumNum);
		semilogy(freq,diffPSD,'r-',label = 'Difference, %s'%DiffNum);
		semilogy(freq,PSN,'g-', linewidth = 3, label = 'PSN');
		legend(fontsize = 18);
	xlim(0,fmax)
	grid('on',which='both');
	xlabel('freq,Hz',fontsize=18);
	ylabel('A/rHz',fontsize=18);
	suptitle('Passive Subtraction, $\phi$ = %.0f $\mu$rad, $I_{sum}$ = %.0f $\mu$A' %(angle,IavgSumuA),fontsize=18);
	legend();
	show();