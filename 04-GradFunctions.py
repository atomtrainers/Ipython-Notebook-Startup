from nptdms import TdmsFile

import numpy as np
from scipy.fftpack import rfft, irfft, fftfreq, fft
from scipy import fft
from scipy.optimize import curve_fit
from scipy import signal;

import warnings
warnings.filterwarnings("ignore")

def GetCompEnvGrad():
	'''Sets directory of gradiometer TDMS files'''
	cname = os.getenv('computername');

	if cname =='RAMSESE-II':
		fpath = 'C:\\Users\\sulai\\Documents\\LabVIEW Data\\DATA\\GradiometerData\\';
	elif cname =='THADLABS':
		fpath = 'D:\\Sulai\\Analysis\\CommonCalib\\';
	return fpath
	
def LoadGradData(filenum='0000'):
	'''loads data from TMDS file'''
	fpath=GetCompEnvGrad();
	
	fname=fpath+'CommonCalib'+str(filenum)+'.tdms'
	
	tdms_file=TdmsFile(fname)
	Ch1 = tdms_file.object('MagCalib','Ch1'); Data1 = Ch1.data;
	Ch2 = tdms_file.object('MagCalib','Ch2'); Data2 = Ch2.data;
	Ch3 = tdms_file.object('MagCalib','Ch3'); Data3 = Ch3.data;
	Ch4 = tdms_file.object('MagCalib','Ch4'); Data4 = Ch4.data;
	
	time=Ch1.time_track()
	
	return [Data1,Data2,Data3,Data4,time]
	
def LoadGradData4(date='2016.01.01',filenum='0000'):
	'''loads 4-channel data from TMDS file'''
	fpath=GetCompEnvGrad();
	
	fname=fpath+date+'\\CommonCalib'+str(filenum)+'.tdms'
	
	tdms_file=TdmsFile(fname)
	Ch1 = tdms_file.object('MagCalib','Ch1'); Data1 = Ch1.data;
	Ch2 = tdms_file.object('MagCalib','Ch2'); Data2 = Ch2.data;
	Ch3 = tdms_file.object('MagCalib','Ch3'); Data3 = Ch3.data;
	Ch4 = tdms_file.object('MagCalib','Ch4'); Data4 = Ch4.data;
	ChRef = tdms_file.object('MagCalib','ChRef'); DataRef= ChRef.data;
	
	time=Ch1.time_track()
	
	return [DataRef,Data1,Data2,Data3,Data4,time] #such that the reference is index 0, Channel 1 is index 1, Channel 2 is index 2, etc
	
def plotTDMS(filenumber='0000',ch='yyyy',xbot=0,xtop=55000,ybot=-10,ytop=10):
	'''plots all four channels from TDMS file. The "ch" variable allows you go select which channels to plot.'''
	data=LoadGradData(filenum=filenumber);
	figure(fignum,figsize=(14,8))
	if ch[0]=='y':
		plot(data[0],'y-',label=filenumber+" AI0");
	if ch[1]=='y':
		plot(data[1],'b-',label=filenumber+" AI1");
	if ch[2]=='y':
		plot(data[2],'m-',label=filenumber+" AI2");
	if ch[3]=='y':
		plot(data[3],'g-',label=filenumber+" AI3");
	xlim(xmin=xbot,xmax=xtop);
	ylim(ymin=ybot,ymax=ytop);
	grid()
	legend();
	
def plotTDMS4(date='2016.01.01',filenumber='0000',ch='yyyy',xbot=0,xtop=55000,ybot=-10,ytop=10,fignum=1):
	'''plots reference and four channels from TDMS file. The "ch" variable allows you go select which channels to plot.'''
	data=LoadGradData4(date,filenum=filenumber);
	figure(fignum,figsize=(14,8))
	plot(data[0],'k-',label=filenumber+" Ref");
	if ch[0]=='y':
		plot(data[1],'k-',label=filenumber+" AI0");
	if ch[1]=='y':
		plot(data[2],'b-',label=filenumber+" AI1");
	if ch[2]=='y':
		plot(data[3],'m-',label=filenumber+" AI2");
	if ch[3]=='y':
		plot(data[4],'g-',label=filenumber+" AI3");
	xlim(xmin=xbot,xmax=xtop);
	ylim(ymin=ybot,ymax=ytop);
	grid()
	legend();
	
def PlotRawGradData(filenum='0000',):
	'''plots data from TDMS file'''
	[D1,D2,D3,D4,time]=LoadGradData(filenum);
	
	figure(1, figsize = [16,10]);
	clf();
	
	subplot(211);
	plot(D1, 'y-', label = 'Ch1');
	plot(D2, 'b-', label = 'Ch2');
	grid('on');
	legend();
	
	subplot(212);
	plot(D4, 'k-', label = 'Applied Signal');
	grid('on');
	legend();
	
def PlotRawGradData234(date='2016.01.01',filenum='0000',ch='yyyy',fignum=1):
	'''plots data from TDMS file'''
	data=LoadGradData4(date,filenum);
	
	figure(fignum, figsize = [16,10]);
	clf();
	
	subplot(211);
	if ch[0]=='y':
		plot(data[1], 'y-', label = 'Ch1');
	if ch[1]=='y':
		plot(data[2], 'b-', label = 'Ch2');
	if ch[2]=='y':
		plot(data[3], 'm-', label = 'Ch3');
	if ch[3]=='y':
		plot(data[4], 'g-', label = 'Ch4');
	grid('on');
	legend();
	
	subplot(212);
	plot(data[0], 'k-', label = 'Applied Signal');
	grid('on');
	legend();
	
def CreateCalibFilters(filenum='0000',fsamp=2500,Rout=1000,AdderFac=1,CoilFac=5.9E-6,DivFac=100,RefStart=200,RefStop=24800,fcutoff=150):
	'''creates calibration filters from raw data. For use with data taken before 2016.'''
	#import data
	[D1,D2,D3,D4,time]=LoadGradData(filenum);
	
	#converts voltage recorded on reference signal to fT, CoilFac is in [T/A]
	sfactor=AdderFac*CoilFac*1E15/(Rout*DivFac);
	
	#selects data
	Sig1 = D1[RefStart:RefStop];
	Sig2 = D2[RefStart:RefStop];
	Ref = D4[RefStart:RefStop]*sfactor;
	
	#AC-couples signals
	Ref = Ref-mean(Ref[-100:-1]);
	Sig1 = Sig1-mean(Sig1[-100:-1]);
	Sig2 = Sig2-mean(Sig2[-100:-1]);
	
	#convert to frequency space
	Dref = fft(Ref);
	Dsig1 = fft(Sig1);
	Dsig2 = fft(Sig2);
	
	#construct calibration between grad signals and reference
	filter1 = Dref/Dsig1;
	filter2 = Dref/Dsig2;
	
	#Window data using fcutoff
	freq = fftfreq(Dref.size, 1./fsamp);
	W = (abs(freq)<fcutoff);
	filter1w = filter1*W;
	filter2w = filter2*W;
	
	#reconstruct calibrated signals
	Dsig1r = Dsig1*filter1;
	Sig1r = ifft(Dsig1r);
	Dsig2r = Dsig2*filter2;
	Sig2r = ifft(Dsig2r);
	
	return [freq,Ref,Sig1,Sig2,filter1,filter2,filter1w,filter2w,Sig1r,Sig2r]

def CreateCalibFilters1(filenum='0000',fsamp=2500,Rout=1000,AdderFac=1,Dir='Y',DivFac=100,RefStart=200,RefStop=24800,fcutoff=150):
	'''creates calibration filters from raw data - modified to use the freq-dependent calibration of the room Y coils. For use with data taken before 2016.'''
	#import data
	[D1,D2,D3,D4,time]=LoadGradData(filenum);
	
	#selects data
	Sig1 = D1[RefStart:RefStop];
	Sig2 = D2[RefStart:RefStop];
	Ref = D4[RefStart:RefStop]; #different from original function - does not apply scale factor yet
	
	#AC-couples signals
	Ref = Ref-mean(Ref[-100:-1]);
	Sig1 = Sig1-mean(Sig1[-100:-1]);
	Sig2 = Sig2-mean(Sig2[-100:-1]);
	
	#convert to frequency space
	Dref = fft(Ref);
	Dsig1 = fft(Sig1);
	Dsig2 = fft(Sig2);
	
	#create frequency spectrum from reference
	freq = fftfreq(Dref.size, 1./fsamp);
	
	if Dir=='Y':
		sfactor1=8.15E-6*1E15/(Rout*DivFac); #DC coil calibration. Converts V to fT
		sfactor2=.0752*abs(freq)**.668+1; #frequency-dependent part of the calibration.
	if Dir=='X':
		sfactor1=5.16E-6; #DC coil calibration
		sfactor2=1; #We haven't measured the frequency-dependent part of the calibration for the x coils. Assume for now that there isn't one.
		
	Dref = Dref*sfactor1/sfactor2 #Converts V to fT, and applies frequency-dependent calibration part
	Ref=ifft(Dref) #converts back to time-space, with the frequency-dependent calibration of the room coils included
	
	#construct calibration between grad signals and reference
	filter1 = Dref/Dsig1;
	filter2 = Dref/Dsig2;
	
	#Window data using fcutoff
	W = (abs(freq)<fcutoff);
	filter1w = filter1*W;
	filter2w = filter2*W;
	
	#reconstruct calibrated signals
	Dsig1r = Dsig1*filter1;
	Sig1r = ifft(Dsig1r);
	Dsig2r = Dsig2*filter2;
	Sig2r = ifft(Dsig2r);
	
	return [freq,Ref,Sig1,Sig2,filter1,filter2,filter1w,filter2w,Sig1r,Sig2r]
	
def CreateCalibFilters2(date='2016.01.01',filenum='0000',Chan=[1,2],fsamp=2500,Rout=1000,Dir='Y',DivFac=100,RefStart=200,RefStop=24800,fcutoff=150):
	'''creates calibration filters from two channels of raw data. For use with data taken 2016.01.01 onward.'''
	
	data=LoadGradData4(date,filenum); #import data
	
	#Reference
	
	Ref = data[0][RefStart:RefStop]; #Selects reference data used for filter, in V
	Ref = Ref-mean(Ref[-100:-1]); #AC-couples reference data
	Dref = fft(Ref); #Converts to frequency space, in V
	freq = fftfreq(Dref.size, 1./fsamp); #Creates frequency spectrum from reference
	
	if Dir=='Y':
		sfactor1=8.15E-6*1E15/(Rout*DivFac); #DC coil calibration. The CoilFac is the current->field conversion at DC
		sfactor2=.0752*abs(freq)**.668+1; #frequency-dependent part of the calibration.
	if Dir=='X':
		sfactor1=5.16E-6; #DC coil calibration
		sfactor2=1; #We haven't measured the frequency-dependent part of the calibration for the x coils. Assume for now that there isn't one.

	Dref = Dref*sfactor1/sfactor2 #Converts V to fT, and applies frequency-dependent calibration part
	Ref=ifft(Dref) #converts back to time-space, with the frequency-dependent calibration of the room coils included
	
	W = (abs(freq)<fcutoff); #Creates "window" that can limit the frequency size of filter (*)
	
	#First Data Channel
	
	index1=Chan[0]; #Channel that recorded first set of sensor data
	Data1=data[index1]; #Selects data
	Sig1 = Data1[RefStart:RefStop]; #Selects data used for filter
	Sig1 = Sig1-mean(Sig1[-100:-1]); #AC-couples data
	Dsig1 = fft(Sig1); #Transforms data to frequency space
	filter1 = Dref/Dsig1; #Creates calibration filter
	filter1w = filter1*W; #Creates "windowed" calibration filter given the window in (*) above
	Dsig1r = Dsig1*filter1; #Applies calibration filter to original data, this is the reconstructed signal in frequency space, in fT
	Sig1r = ifft(Dsig1r); #Transforms back to time space; this is the reconstructed signal in time space, in fT
	
	#Second Data Channel
	
	index2=Chan[1]; #Channel that recorded second set of sensor data
	Data2=data[index2]; #Selects data
	Sig2 = Data2[RefStart:RefStop]; #Selects data used for filter
	Sig2 = Sig2-mean(Sig2[-100:-1]); #AC-couples data
	Dsig2 = fft(Sig2); #Transforms data to frequency space
	filter2 = Dref/Dsig2; #Creates calibration filter
	filter2w = filter2*W; #Creates "windowed" calibration filter given the window in (*) above
	Dsig2r = Dsig2*filter2; #Applies calibration filter to original data, this is the reconstructed signal in frequency space,  in fT
	Sig2r = ifft(Dsig2r); #Transforms back to time space; this is the reconstructed signal in time space, in fT
	
	return [freq,Ref,index1,Sig1,filter1,filter1w,Sig1r,index2,Sig2,filter2,filter2w,Sig2r]
	
def PlotReconstructed(RunInfo,fignum=1):
	'''Plots reconstructed calibration signals'''
	#plots raw data
	figure(fignum,figsize = [16,10]);
	rawplot=subplot(211)
	plot(RunInfo[1], 'k-', label = 'Ref',);
	#plot(RunInfo[2], 'y-', label = 'Ch1');
	#plot(RunInfo[3], 'b-', label = 'Ch2');
	ylabel('Applied Reference Signal [fT]');
	legend();
	grid();
	
	#plots reconstructed data
	subplot(212)
	plot(RunInfo[1], 'k-', label = 'reference');
	plot(RunInfo[8].real, 'y-', label = 'reconstructed - Ch1');
	plot(RunInfo[9], 'b-', label = 'reconstructed - Ch2');
	grid();
	legend();

def PlotReconstructed2(RunInfo,fignum=1):
	'''Plots reconstructed calibration signals'''
	#plots reference signal
	figure(fignum,figsize = [16,10]);
	suptitle('Calibration Reconstruction',fontsize=18)
	rawplot=subplot(211)
	plot(RunInfo[1], 'k-', label = 'Reference',);
	ylabel('Applied Reference Signal [fT]');
	legend();
	grid();
	
	ColorMap=['k-','y-','b-','m-','-g']
	colorRef=ColorMap[0]
	
	label1='Reconstructed - Ch. '+str(RunInfo[2]);
	color1=ColorMap[RunInfo[2]];
	
	label2='Reconstructed - Ch. '+str(RunInfo[7]);
	color2=ColorMap[RunInfo[7]];
	
	#plots reconstructed data
	subplot(212);
	plot(RunInfo[6].real, color1, label = label1);
	plot(RunInfo[11].real, color2, label = label2);
	grid();
	legend();

def SelectNoiseData(filenum='0000',NoiseStart=40000,NoiseLength=12000,filter1='a',filter2='b'):
	'''Selects noise portion of data, takes FFT, and applies filters, use with data stored in old format (2015)'''
	#import data
	[D1,D2,D3,D4,time]=LoadGradData(filenum);
	
	#select data of interest
	Sig1 = D1[NoiseStart:NoiseStart+NoiseLength]; 
	Sig2 = D2[NoiseStart:NoiseStart+NoiseLength];
	
	#AC-couple
	Sig1 = Sig1-mean(Sig1)
	Sig2 = Sig2-mean(Sig2)
	
	#fourier transform
	Dsig1 = fft(Sig1)
	Dsig2 = fft(Sig2)
	
	#apply filters and transform
	Dsig1r = Dsig1*filter1
	Sig1r = ifft(Dsig1r)
	#Sig1r = Sig1r[5000:20000]
	Dsig2r = Dsig2*filter2
	Sig2r = ifft(Dsig2r)
	#Sig2r = Sig2r[5000:20000]
	
	figure(1, figsize = [12,10]);
	clf();
	
	subplot(211);
	plot(Sig1r, 'y-', label = 'Ch1recon');
	plot(10**5*Sig1, 'm-', label = 'Ch1raw');
	grid('on');
	legend();
	
	subplot(212);
	plot(Sig2r, 'b-', label = 'Ch2recon');
	plot(10**5*Sig1, 'g-', label = 'Ch2raw');
	grid('on');
	legend();
	
	return [Sig1r.real,Sig2r.real]

def SelectNoiseData2(date='2016.01.01',filenum='0000',Chan=[1,2],NoiseStart=40000,NoiseLength=12000,filter1='a',filter2='b',fignum=1):
	'''Selects noise portion of data, takes FFT, and applies filters, use with data stored in old format (2016)'''
	
	data=LoadGradData4(date,filenum); #import data
	
	#First set of noise data
	
	index1=Chan[0]; #Channel that recorded first set of sensor data
	Data1=data[index1]; #Selects data entire dataset
	Sig1 = Data1[NoiseStart:NoiseStart+NoiseLength]; #Selects data that will be analyzed 
	Sig1 = Sig1-mean(Sig1); #AC-couples data
	Dsig1 = fft(Sig1); #Transforms to frequency space
	Dsig1r = Dsig1*filter1; #Applies filter to data, in doing so converts to fT
	Sig1r = ifft(Dsig1r); #Transforms data back to time space, in fT
	
	#Second set of noise data
	
	index2=Chan[1]; #Channel that recorded first set of sensor data
	Data2=data[index2]; #Selects data entire dataset
	Sig2 = Data2[NoiseStart:NoiseStart+NoiseLength]; #Selects data that will be analyzed 
	Sig2 = Sig2-mean(Sig2); #AC-couples data
	Dsig2 = fft(Sig2); #Transforms to frequency space
	Dsig2r = Dsig2*filter2; #Applies filter to data, in doing so converts to fT
	Sig2r = ifft(Dsig2r); #Transforms data back to time space, in fT
	
	ColorMap=['k-','y-','b-','m-','-g']
	colorRaw=ColorMap[0]
	
	label1A='Raw - Ch. '+str(index1)
	label1B='Reconstructed - Ch. '+str(index1);
	color1=ColorMap[index1];
	
	label2A='Raw - Ch. '+str(index2)
	label2B='Reconstructed - Ch. '+str(index2);
	color2=ColorMap[index2];
	
	figure(fignum, figsize = [16,10]);
	suptitle('Noise Reconstruction',fontsize=18)
	chA=subplot(211);
	plot(Sig1, colorRaw, label = label1A);
	ylabel('Signal [V]');
	xlabel('point');
	legend(loc=2);
	grid('on');
	chA2=chA.twinx();
	plot(Sig1r, color1, label = label1B);
	ylabel('Reconstructed Signal [fT]');
	legend();
	
	chB=subplot(212);
	plot(Sig2, colorRaw, label = label2A);
	ylabel('Signal [V]');
	xlabel('point');
	legend(loc=2);
	grid('on');
	chB2=chB.twinx();
	plot(Sig2r, color2, label = label2B);
	ylabel('Reconstructed Signal [fT]');
	legend();
	
	return [Sig1,Sig1r.real,Sig2,Sig2r.real]

def	FilterCompare(filter1,filter2,freq,freqmax):
	f1shift=fftshift(filter1);
	f2shift=fftshift(filter2);
	freqshift=fftshift(freq);
	
	angle1=angle(f1shift);
	angle2=angle(f2shift);
	
	figure(1, figsize=[14,10])
	subplots_adjust(left=None, bottom=None, right=None, top=None, wspace=0.32, hspace=0.32)
	
	subplot(221)
	semilogy(freqshift,abs(f1shift),'b-');
	semilogy(freqshift,abs(f2shift),'g-');
	grid('on');
	xlim(-freqmax,freqmax);
	ylabel('Gain');
	
	subplot(222)
	plot(freqshift,angle1,'bo');
	plot(freqshift,angle2,'go');
	grid('on');
	xlim(-freqmax,freqmax);
	ylabel('Phase');
	
	subplot(223)
	plot(freqshift,abs(f1shift)-abs(f2shift), 'r-');
	xlim(-freqmax,freqmax);
	ylim(-2e5,2e5);
	grid('on');
	ylabel('Gain Difference');
	xlabel('Frequency [Hz]');
	
	subplot(224)
	plot(freqshift,angle1- angle2, 'r-');
	grid('on');
	xlim(-100,100);
	ylim(-1e-0,1e-0);
	ylabel('phase diff');
	xlabel('Frequency [Hz]');
	
def getCh1Ch2PCAv2(ch1,ch2,LPF = 100., HPF = 1.,bs =1,rScale = 1, note1 = '', note2 = '', tmsize = 5000,Fsamp = 2500,minY=10,cmMax = 100,fmax=200):
	'''Plots raw time-series from two channels, as well as PCA-transformed data (ie, sum and difference). Also plots frequency-space FFTs of raw and PCA data, as well as CMRR of the PCA data'''
	
	#set Nyquist frequency
	Nyq = Fsamp/2.0;

	#apply lowpass and highpass filters
	bLP,aLP = butter(2,[HPF/Nyq, LPF/Nyq], btype = 'bandpass');
	C1 = filtfilt(bLP,aLP,ch1);
	C2 = filtfilt(bLP,aLP,ch2);

	#chop off edge effects
	trimFactor = Fsamp/2.
	sz = size(C1); sz2 = sz-trimFactor
	C1,C2 = [C1[trimFactor:sz2],C2[trimFactor:sz2]]

	#create PCA eigenvectors
	data = rScale*array([C1,C2])
	w,v=eig(cov(data))	
	a=dot(inv(v),data)

	sz = size(C1)

	#create PSD of raw data
	P1a,freq = mlab.psd(data[0], NFFT =sz, Fs = Fsamp); P1a = sqrt(P1a)
	P1b,freq = mlab.psd(data[1], NFFT =sz, Fs = Fsamp); P1b = sqrt(P1b)

	#create PSD of PCA-transformed data
	P2a,freq = mlab.psd(a[0], NFFT =sz, Fs = Fsamp); P2a = sqrt(P2a)
	P2b,freq = mlab.psd(a[1], NFFT =sz, Fs = Fsamp); P2b = sqrt(P2b)

	#bin data
	P1a = binit(P1a,bs);
	P1b = binit(P1b,bs);
	freq = binit(freq,bs);

	P2a = binit(P2a,bs);
	P2b = binit(P2b,bs);

	timepoints=size(data[0]);
	time=linspace(1,timepoints,timepoints)/Fsamp; #time in seconds

	#Plots data
	figure(1, figsize=[14,10])
	subplots_adjust(left=None, bottom=None, right=None, top=None, wspace=None, hspace=0.32)

	#subplots 621 and 623 are the time series of the raw data
	subplot(621);
	plot(time,1e-3*data[0],'b-');
	grid('on');
	xlim(xmax = ln1)
	ylabel('field, pT');
	title('Raw Data')

	subplot(623)
	plot(1e-3*data[1],'g-')
	grid('on');
	xlabel('time, ms');
	ylabel('field, pT');
	xlim(xmax = ln1)

	#subplot 323 is the power spectral density of the raw data
	subplot(323)
	semilogy(freq,P1a);
	semilogy(freq,P1b);
	grid('on',which = 'both');
	xlim(xmax = fmax)
	grid('on');
	ylim(ymin = minY)
	xlabel('freq,Hz');
	ylabel('PSD, fT/rHz');
	#ylim(ymin = .1);

	#subplots 622 and 624 are the time series of the PCA-transformed data
	subplot(622);
	plot(1e-3*a[0],'b-');
	xlim(xmax = ln1)
	grid('on')
	ylabel('field, pT');
	title('PCA')

	subplot(624)
	plot(1e-3*a[1],'g-')
	xlim(xmax = ln1)
	grid('on');
	xlabel('time, ms');
	ylabel('field, pT');

	#subplot 324 is the power spectral density of the PCA-transformed data
	subplot(324)
	semilogy(freq,P2a);
	semilogy(freq,P2b);
	grid('on',which = 'both');
	xlim(xmax = fmax)
	grid('on');
	xlabel('freq,Hz');
	ylabel('PSD, fT/rHz');
	ylim(ymin = minY);

	#calculates the ratio of the two PCA-transformed PSDs (ideally, sum and difference).
	PratioA = P2b/P2a;
	PratioB = P2a/P2b;

	#Depending on how the the PCA assigned its sum and difference, 
	#either PratioA or PratioB will be the CMRR
	if mean(PratioA[0:20])>mean(PratioB[0:20]):
		Pratio = PratioA;
	else:
		Pratio = PratioB;

	#subplot 326 is the CMRR vs frequency
	subplot(326)
	plot(freq,Pratio, 'r-');
	grid('on',which = 'both');
	xlim(xmax = fmax)
	ylim(ymax = cmMax)
	xlabel('freq,Hz');
	ylabel('PSD ratio')

	#Adds a notation to the figure
	fignotes = subplot(325)
	fignotes.axis('off');
	fs = 16
	frame1 = gca()
	frame1.axes.get_xaxis().set_visible(False)
	frame1.axes.get_yaxis().set_visible(False)

	bptext = 'Band Pass: ' + str(HPF) +' - '+str(LPF) +  ' Hz';
	evtext = 'EigenVector = '+str(v[0]);


	text(0.1,0.8,bptext, fontsize = fs)
	text(0.1,0.5,note1, fontsize = fs)
	text(0.1,0.2,note2, fontsize = fs)
	text(0.1,0.0,evtext, fontsize = 12)

def getCh1Ch2GRAv2(ch1,ch2,LPF = 100., HPF = 1.,bs =1, p2 = 10, note1 = '', note2 = '',Fsamp=2500,rScale=1,tmsize = 10,gmode = [0.707,0.707],minY=1,cmMax=200,fmax=200,fignum=1):
	'''Basically the same function as getCh1Ch2PCAv2, except forces both eigenvalues to be 1/sqrt(2); essentially straight subtraction'''
	
	#set Nyquist frequency
	Nyq = Fsamp/2.0;

	#apply lowpass and highpass filters
	bLP,aLP = butter(2,[HPF/Nyq, LPF/Nyq], btype = 'bandpass');
	C1 = filtfilt(bLP,aLP,ch1);
	C2 = filtfilt(bLP,aLP,ch2);

	#chop off edge effects
	trimFactor = Fsamp/2.
	sz = size(C1); sz2 = sz-trimFactor
	C1,C2 = [C1[trimFactor:sz2],C2[trimFactor:sz2]]

	data = rScale*array([C1,C2])
	w,v=eig(cov(data))
	#print 'The eigenvalues are %.2f and %.2f' % (w[0],w[1]) #eigenvalues
	#print 'The eigenvectors are' 
	#print v #eigenvectors 
	
	#sets the two eigenvalues to the user-defined ones, by default [1/sqrt(2), 1/sqrt(2)]
	g1 = gmode[0]; g2 = gmode[1]
	gMatrix = array([[g1,g2],[g2,-g1]])
	v = gMatrix;
	a=dot(inv(gMatrix),data)


	sz = size(C1);

	#calculates PSD of raw data
	P1a,freq = mlab.psd(data[0], NFFT =sz, Fs = Fsamp); P1a = sqrt(P1a)
	P1b,freq = mlab.psd(data[1], NFFT =sz, Fs = Fsamp); P1b = sqrt(P1b)

	#calculates PSD of transformed data
	P2a,freq = mlab.psd(a[0], NFFT =sz, Fs = Fsamp); P2a = sqrt(P2a)
	P2b,freq = mlab.psd(a[1], NFFT =sz, Fs = Fsamp); P2b = sqrt(P2b)


	#fmax = 2*LPF

	#bins data
	P1a = binit(P1a,bs);
	P1b = binit(P1b,bs);
	freq = binit(freq,bs);


	P2a = binit(P2a,bs);
	P2b = binit(P2b,bs);

	timepoints=size(data[0]);
	time=linspace(1,timepoints,timepoints)/Fsamp; #time in seconds

	#plots data
	figure(1, figsize=[14,10])
	subplots_adjust(left=None, bottom=None, right=None, top=None, wspace=None, hspace=0.32)

	#subplots 621 and 623 are the time series of the raw data
	subplot(621);
	plot(time,1e-3*data[0],'b-');
	grid('on');
	xlim(xmax = tmsize)
	ylabel('field, pT');
	title('Raw Data')

	subplot(623)
	plot(time,1e-3*data[1],'g-')
	grid('on');
	xlabel('time, s');
	ylabel('field, pT');
	xlim(xmax = tmsize)

	#subplot 323 is the PSD of the raw data
	subplot(323)
	semilogy(freq,P1a);
	semilogy(freq,P1b);
	grid('on',which = 'both');
	xlim(xmax = fmax)
	grid('on');
	xlabel('freq,Hz');
	ylabel('PSD, fT/rHz');
	ylim(ymin = .1);

	#subplots 622 and 624 are the time series of the transformed data
	subplot(622);
	plot(time,1e-3*a[0],'b-');
	xlim(xmax = tmsize)
	grid('on')
	ylabel('field, pT');
	title('PCA')

	subplot(624)
	plot(time,1e-3*a[1],'g-')
	xlim(xmax = tmsize)
	grid('on');
	xlabel('time, s');
	ylabel('field, pT');

	#subplot 324 is the PSDs of the transformed data
	subplot(324)
	semilogy(freq,P2a);
	semilogy(freq,P2b);
	grid('on',which = 'both');
	xlim(xmax = fmax)
	grid('on');
	xlabel('freq,Hz');
	ylabel('PSD, fT/rHz');
	ylim(ymin = .1);

	#calculate ratios of transformed data
	PratioA = P2b/P2a;
	PratioB = P2a/P2b;

	#depending on how the data transformed, one of the 
	#[PratioA,PratioB] is the [sum,difference]. The ratio
	#of the two is the CMRR
	if mean(PratioA)>mean(PratioB):
		Pratio = PratioA;
	else:
		Pratio = PratioB;

	#subplot 326 is the CMRR vs frequency
	subplot(326)
	plot(freq,Pratio, 'r-');
	grid('on',which = 'both');
	xlim(xmax = fmax);
	ylim(ymax = cmMax);
	xlabel('freq,Hz');
	ylabel('PSD ratio')

	#adds an annotation
	fignotes = subplot(325)
	fignotes.axis('off');
	fs = 16
	frame1 = gca()
	frame1.axes.get_xaxis().set_visible(False)
	frame1.axes.get_yaxis().set_visible(False)

	bptext = 'Band Pass: ' + str(HPF) +' - '+str(LPF) +  ' Hz';
	evtext = 'EigenVector = '+str(v[0]);
	
	
	text(0.1,0.8,bptext, fontsize = fs)
	text(0.1,0.5,note1, fontsize = fs)
	text(0.1,0.2,note2, fontsize = fs)
	text(0.1,0.0,evtext, fontsize = 12)
	
	#figtitle = day+'Run '+run+' Noise '+noise;
	#suptitle(figtitle);

	show();
	
	#returns time, raw time series (data[0] and data[1]), PCA time series (a[0] and a[1]),
	#fft frequency, raw PSDs (P1a and P1b), PCA PSDs (P2a and P2b) 
	return time,data[0],data[1],a[0],a[1],freq,P1a,P1b,P2a,P2b,Pratio

def GradSubtract(ch1,ch2,LPF = 100., HPF = 1.,bs =1,Fsamp=2500,tmsize = 10,PCAGRA='GRA',cmMax=200,fmax=200, note1 = '', note2 = '',fignum=1):
	'''Combines the getCh1Ch2GRAv2 and getCh1Ch2PCAv2 functions. Select using the PCAGRA variable'''
	
	#set Nyquist frequency
	Nyq = Fsamp/2.0;

	#apply lowpass and highpass filters
	bLP,aLP = butter(2,[HPF/Nyq, LPF/Nyq], btype = 'bandpass');
	C1 = filtfilt(bLP,aLP,ch1);
	C2 = filtfilt(bLP,aLP,ch2);

	#chop off edge effects
	trimFactor = Fsamp/2.
	sz = size(C1); sz2 = sz-trimFactor
	C1,C2 = [C1[trimFactor:sz2],C2[trimFactor:sz2]]

	data = array([C1,C2])
	
	if PCAGRA=='GRA':
		g1 = 0.707; g2 = 0.707
		v = array([[g1,g2],[g2,-g1]])
		a=dot(inv(v),data)
	elif PCAGRA=='PCA':
		w,v=eig(cov(data))	
		a=dot(inv(v),data)

	sz = size(C1);

	#calculates PSD of raw data
	P1a,freq = mlab.psd(data[0], NFFT =sz, Fs = Fsamp); P1a = sqrt(P1a)
	P1b,freq = mlab.psd(data[1], NFFT =sz, Fs = Fsamp); P1b = sqrt(P1b)

	#calculates PSD of transformed data
	P2a,freq = mlab.psd(a[0], NFFT =sz, Fs = Fsamp); P2a = sqrt(P2a)
	P2b,freq = mlab.psd(a[1], NFFT =sz, Fs = Fsamp); P2b = sqrt(P2b)


	#fmax = 2*LPF

	#bins data
	P1a = binit(P1a,bs);
	P1b = binit(P1b,bs);
	freq = binit(freq,bs);


	P2a = binit(P2a,bs);
	P2b = binit(P2b,bs);

	timepoints=size(data[0]);
	time=linspace(1,timepoints,timepoints)/Fsamp; #time in seconds

	#plots data
	figure(fignum, figsize=[16,10])
	suptitle('Noise Subtraction ('+PCAGRA+')',fontsize=18)
	subplots_adjust(left=None, bottom=None, right=None, top=None, wspace=None, hspace=0.32)

	#subplots 621 and 623 are the time series of the raw data
	subplot(621);
	plot(time,1e-3*data[0],'b-');
	grid('on');
	xlim(xmax = tmsize)
	ylabel('field, pT');
	title('Raw Data')

	subplot(623)
	plot(time,1e-3*data[1],'g-')
	grid('on');
	xlabel('time, s');
	ylabel('field, pT');
	xlim(xmax = tmsize)

	#subplot 323 is the PSD of the raw data
	subplot(323)
	semilogy(freq,P1a);
	semilogy(freq,P1b);
	grid('on',which = 'both');
	xlim(xmax = fmax)
	grid('on');
	xlabel('freq,Hz');
	ylabel('PSD, fT/rHz');
	ylim(ymin = .1);

	#subplots 622 and 624 are the time series of the transformed data
	subplot(622);
	plot(time,1e-3*a[0],'b-');
	xlim(xmax = tmsize)
	grid('on')
	ylabel('field, pT');
	title('PCA')

	subplot(624)
	plot(time,1e-3*a[1],'g-')
	xlim(xmax = tmsize)
	grid('on');
	xlabel('time, s');
	ylabel('field, pT');

	#subplot 324 is the PSDs of the transformed data
	subplot(324)
	semilogy(freq,P2a);
	semilogy(freq,P2b);
	grid('on',which = 'both');
	xlim(xmax = fmax)
	grid('on');
	xlabel('freq,Hz');
	ylabel('PSD, fT/rHz');
	ylim(ymin = .1);

	#calculate ratios of transformed data
	PratioA = P2b/P2a;
	PratioB = P2a/P2b;

	#depending on how the data transformed, one of the 
	#[PratioA,PratioB] is the [sum,difference]. The ratio
	#of the two is the CMRR
	if mean(PratioA)>mean(PratioB):
		Pratio = PratioA;
	else:
		Pratio = PratioB;

	#subplot 326 is the CMRR vs frequency
	subplot(326)
	plot(freq,Pratio, 'r-');
	grid('on',which = 'both');
	xlim(xmax = fmax);
	ylim(ymax = cmMax);
	xlabel('freq,Hz');
	ylabel('PSD ratio')

	#adds an annotation
	fignotes = subplot(325)
	fignotes.axis('off');
	fs = 16
	frame1 = gca()
	frame1.axes.get_xaxis().set_visible(False)
	frame1.axes.get_yaxis().set_visible(False)

	bptext = 'Band Pass: ' + str(HPF) +' - '+str(LPF) +  ' Hz';
	evtext = 'EigenVector = '+str(v[0]);
	
	
	text(0.1,0.8,bptext, fontsize = fs)
	text(0.1,0.5,note1, fontsize = fs)
	text(0.1,0.2,note2, fontsize = fs)
	text(0.1,0.0,evtext, fontsize = 12)
	
	#figtitle = day+'Run '+run+' Noise '+noise;
	#suptitle(figtitle);

	show();
	
	#returns time, raw time series (data[0] and data[1]), PCA time series (a[0] and a[1]),
	#fft frequency, raw PSDs (P1a and P1b), PCA PSDs (P2a and P2b), and CMRR (Pratio)
	return time,data[0],data[1],a[0],a[1],freq,P1a,P1b,P2a,P2b,Pratio	

def plotMagRun4Grad(day='2015.01.01',runNum='00',noiseNum = '00',noiseNum2 = '01', noiseNum3 = '02',bs = 4,fmax = 200,fignum = 1,ver = 'v16',leg=1,rScale=5,phaselim=-180,direc = 'Y', el = 1,Isum1=100E-6,calib1=2E-6,Isum2=100E-6,calib2=2E-6):
	[basepath, analysispath] = getCompEnv();
	basepath = basepath + ver + '\\';
	
	dayDir = analysispath+day
	a = os.path.isdir(dayDir);
	if a==0:
		os.mkdir(dayDir);
	
	runDir = analysispath+day+'\\Run'+runNum
	a = os.path.isdir(runDir);
	if a==0:
		os.mkdir(runDir);
	
	fs = 14;
	
	Chan = 1;
	path = basepath +str(day) + '\\' + 'run_' + str(runNum);
	
	fPSD,fPSD2,fPSD3,fRESP,figTitle = getRunSpecs(day,path,noiseNum,noiseNum2,noiseNum3,Chan,direc,runNum)
	
	
	if ver == 'v12':
		fPSD = path + '\\'+'noise_'+str(noiseNum)+'\\'+'magnetic-field-PSD_'+str(Chan)+'.txt';
		fPSD2 = path + '\\'+'noise_'+str(noiseNum2)+'\\'+'magnetic-field-PSD_'+str(Chan)+'.txt';
		fPSD3 = path + '\\'+'noise_'+str(noiseNum3)+'\\'+'magnetic-field-PSD_'+str(Chan)+'.txt';
		fRESP =  path+'\\'+'response'+'\\'+'amp-n-phase_response_'+str(Chan)+'.txt';
	
	data = loadtxt(fPSD);
	data = transpose(data);
	freq = binit(data[0],bs);
	PS = binit(data[1],bs);
	dataResp = loadtxt(fRESP);
	dataResp1 = transpose(dataResp)
	
	data = loadtxt(fPSD2);
	data = transpose(data);
	freq = binit(data[0],bs);
	PS2 = binit(data[1],bs);
	
	data = loadtxt(fPSD3);
	data = transpose(data);
	freq = binit(data[0],bs);
	PS3 = binit(data[1],bs);
	
	shotnoise1=sqrt(4*Isum1*1.6E-19)/calib1/dataResp1[1];
	
	figs=figure(figsize = (16,8))
	clf();
	
	ticklblsize=18
	axlblsize=20
	titsize=20
	legsize=18
	
	matplotlib.rc('xtick',labelsize=ticklblsize)
	matplotlib.rc('ytick',labelsize=ticklblsize)
	
	subplot(131)
	semilogy(freq,PS,'k-',label='Total Noise');
	semilogy(freq,PS2,'r-',label='Probe Noise');
	if el == 1:
		semilogy(freq,PS3,'b-',label='Electronic Noise');
	semilogy(dataResp1[0],shotnoise1,'g-', linewidth = 3,label='photon shot noise');
	grid(which='both');
	axis([0,fmax,0.1,10000]);
	ylabel('Power Spectral Density [fT/rHz]',fontsize=axlblsize);
	xlabel('Frequency [Hz]',fontsize=axlblsize);
	title('Channel 1',fontsize=titsize)
	#~ legend();
	#~ legend();
	#text(25,2000,'Ch1',fontsize=fs)
	
	ns1 = PS;
	
	Chan = 2;
	path = basepath +str(day) + '\\' + 'run_' + str(runNum);
	fPSD,fPSD2,fPSD3,fRESP,figTitle = getRunSpecs(day,path,noiseNum,noiseNum2,noiseNum3,Chan,direc,runNum)
	
	if ver == 'v12':
		fPSD = path + '\\'+'noise_'+str(noiseNum)+'\\'+'magnetic-field-PSD_'+str(Chan)+'.txt';
		fPSD2 = path + '\\'+'noise_'+str(noiseNum2)+'\\'+'magnetic-field-PSD_'+str(Chan)+'.txt';
		fPSD3 = path + '\\'+'noise_'+str(noiseNum3)+'\\'+'magnetic-field-PSD_'+str(Chan)+'.txt';
		fRESP =  path+'\\'+'response'+'\\'+'amp-n-phase_response_'+str(Chan)+'.txt';
	
	
	data = loadtxt(fPSD);
	data = transpose(data);
	freq = binit(data[0],bs);
	PS = binit(data[1],bs);
	dataResp = loadtxt(fRESP);
	dataResp2 = transpose(dataResp)
	
	data = loadtxt(fPSD2);
	data = transpose(data);
	freq = binit(data[0],bs);
	PS2 = binit(data[1],bs);
	
	data = loadtxt(fPSD3);
	data = transpose(data);
	freq = binit(data[0],bs);
	PS3 = binit(data[1],bs);
	
	shotnoise2=sqrt(4*Isum1*1.6E-19)/calib2/dataResp2[1];
	
	subplot(132)
	semilogy(freq,PS,'k-',label='Total Noise');
	semilogy(freq,PS2,'r-',label='Probe Noise');
	if el == 1:
		semilogy(freq,PS3,'b-',label='Electronic Noise');
	semilogy(dataResp2[0],shotnoise2,'g-', linewidth = 3,label='photon shot noise');
	grid(which='both');
	axis([0,fmax,0.1,10000])
	xlabel('Frequency [Hz]',fontsize=axlblsize);
	# ylabel('PSD, fT/rHz');
	title('Channel 2',fontsize=titsize);
	
	if leg==1:
		legend(fontsize=legsize);
	#text(100,100,'Ch2',fontsize=fs)
		
	ns2 = PS;
		
	shotnoise2=sqrt(4*Isum2*1.6E-19)/(calib2*dataResp2[1]);
	
	xR = subplot(233);
	plot(dataResp1[0],1E6*dataResp1[1],'y-', linewidth = 3, label = 'Ch1');
	plot(dataResp2[0],1E6*dataResp2[1],'b-', linewidth = 3, label = 'Ch2');
	axis([0,fmax,0,rScale]);
	xR.yaxis.set_label_position("right");
	ylabel('Response [$\mu$V/fT]',fontsize=axlblsize);
	grid(which='both');
	legend(fontsize=legsize);
	
	xR2 = subplot(236);
	plot(dataResp1[0],dataResp1[2] - dataResp1[2][0],'y-', linewidth = 3, label = 'Ch1');
	plot(dataResp2[0],dataResp2[2] - dataResp2[2][0],'b-', linewidth = 3, label = 'Ch2');
	axis([0,fmax,phaselim,10])
	xR2.yaxis.set_label_position("right");
	xlabel('Frequency [Hz]',fontsize=axlblsize);
	ylabel('Phase [deg]',fontsize=axlblsize);
	grid(which='both');
	#legend()
		
	#figTitle = 'Individual Channel Performance'
	fname1 = dayDir+'\\Mag4'+str(direc)+str(day)+'_'+'Run'+str(runNum)
	fname2 = runDir+'\\Mag4'+str(direc)+str(day)+'_'+'Run'+str(runNum)
	#suptitle(figTitle, fontsize=18);
	#show()
	figs.savefig(fname2+'.pdf', format = 'pdf')
	figs.savefig(fname1+'.png')
	
	return figs
	
def GradiometerComplete(date='2016.01.01',filenum='0000',Chan=[1,2],fsamp=2500,Rout=1000,Dir='Y',DivFac=100,RefStart=200,RefStop=24800,fcutoff=150,NoiseStart=40000,NoiseLength=12000,PCALPF=200,PCAHPF=1,bs=1,tmsize=10,PCAGRA='GRA',cmMax=200,fmax=120, note1 = '', note2 = ''):
	'''Runs a complete batch of gradiometer functions. For use with data acquired January 2016 and on.'''
	calibfilter=CreateCalibFilters2(date=date,filenum=filenum,Chan=Chan,fsamp=fsamp,Rout=Rout,Dir=Dir,DivFac=DivFac,RefStart=RefStart,RefStop=RefStop,fcutoff=fcutoff);
	
	PlotReconstructed2(calibfilter,fignum=1);
	
	magnoise=SelectNoiseData2(date=date,filenum=filenum,Chan=Chan,NoiseStart=NoiseStart,NoiseLength=NoiseLength,filter1=calibfilter[5],filter2=calibfilter[10],fignum=2);
	
	GradSubtract(ch1=magnoise[1],ch2=magnoise[3],LPF=PCALPF,HPF=PCAHPF,bs=bs,Fsamp=fsamp,tmsize=tmsize,PCAGRA=PCAGRA,cmMax=cmMax,fmax=fmax, note1 = note1, note2 = note2,fignum=3)
	