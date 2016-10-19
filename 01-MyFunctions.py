import os;
from scipy.signal import *
from scipy.optimize import curve_fit
from numpy import *
from pylab import *
import mdp
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider, Button, RadioButtons

def AddIndexExt(BaseName,extn):
    """ This function looks for the last file with 'BaseName', and increments the number by 1. It then returns the new file name as a string. If the file does not exist, it is created with an index of '0'"""
    
    extn = '.'+extn;
    cmd = 'dir '+BaseName+'*/B';
        
    lst = [lst[:-1] for lst in os.popen(cmd)]
    if lst==[]:   #File with such a base doesnt exist
        newBase = BaseName+str(0)
        newFname = BaseName + str(0)+ extn;
           
    
    else:
        fn = lst[-1];
        p1 = len(BaseName);
        p2 = fn.find(".");
            
        lastnum = int(fn[p1:p2]);
        if lastnum == '':
            newBase = BaseName+str(0)
            newFname = BaseName+str(0) + extn;
        else:
            newBase = BaseName+str(lastnum+1)    
            newFname = BaseName + str(lastnum+1)+ extn;
        
    return [newBase,newFname];
    

def getCurrentNoise(day, run, noiseSum,noisediff,calibSum,calibDiff,Chan = 3):
    rcParams["font.size"] = "12"
    dataSum = plotMagRun3(day,run,noiseSum,Chan = Chan, fignum = 1,showfig=0)
    dataProbe = plotMagRun3(day,run,noisediff,Chan = Chan, fignum = 1,showfig=0)
    close(1)

    dSum = calibSum*dataSum[2];

    dAv = 0.5*mean(dSum);  #Mean current on single photodiode
    dAv = abs(dAv)
    PSN = sqrt(4*1.6E-19*dAv)

    dProbe = calibDiff*dataProbe[2];

    Pss,freq = mlab.psd(dSum, NFFT = size(dSum), Fs = 1000.)
    Ppr,freq = mlab.psd(dProbe, NFFT = size(dProbe), Fs = 1000.)

    Pss = sqrt(Pss); Ppr = sqrt(Ppr); 
    Psn = freq*0+PSN;

    bf = binit(freq,5);
    cmr = Pss/Ppr;
    cmrr = 20*log10(binit(cmr,5));

    figure()
    subplot(211)
    loglog(freq,abs(Pss),'b-',label = 'Sum');
    loglog(freq,abs(Ppr),'r-',label = 'Difference');
    loglog(freq,Psn,'g-',linewidth = 3,label = 'Photon Shot Noise');
    legend(loc=3,fontsize='small')
    ylabel('A/rHz')
    xlabel('freq,Hz');
    grid('on')
    xlim(xmin =0.1,xmax = 500)
    
  
    subplot(212)
    loglog(bf,abs(cmrr), 'k-',label='CMRR')
    xlim(xmin =0.1,xmax = 500)
    grid('on', which='both')
    xlabel('freq, Hz');
    ylabel('CMRR, dB');
    legend(loc=3)
    
def get1ChanData(day,runNum,noiseNum,Chan = 1,LPF = 80, HPF = 2,f120=120,f60=60, ver = 'v16',calib = 0, lineNotch = 1,fan=0):
	#The calib flag (0 or 1) determines whether it is the raw or calibrated data that is returned.
	
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
	
	noiseNum = int(noiseNum);
		
	if  noiseNum<10:
		noiseNum = '0'+str(noiseNum);

	path = basepath +str(day) + '\\' + 'run_' + str(runNum);
	
	if calib == 0:
		ftser1 = path + '\\'+'noise_'+str(noiseNum)+'\\'+'Y_noise_'+str(Chan)+'.bin';

	if calib == 1:
		ftser1 = path + '\\'+'noise_'+str(noiseNum)+'\\'+'Y_calibrated_time_series_fT_'+str(Chan)+'.bin';


	dataTS1 = fromfile(ftser1,dtype = '<d')

	Nyq = 500.
	#100 Hz low pass
	bLP,aLP = butter(4,[HPF/Nyq, LPF/Nyq], btype = 'bandpass');
	b120,a120 = butter(2,[118./Nyq,122./Nyq],btype='bandstop');
	b60,a60 = butter(4,[59./Nyq,61./Nyq],btype='bandstop');
	bfan,afan = butter(2,[4./Nyq,5./Nyq],btype='bandstop');

	dataTSf1 = filtfilt(bLP,aLP,dataTS1);
	
	if lineNotch==1:
		dataTSf1 = filtfilt(b120,a120,dataTSf1);
		dataTSf1 = filtfilt(b60,a60,dataTSf1);

	if fan==1:
		dataTSf1 = filtfilt(bfan,afan,dataTSf1)


	return [dataTS1,dataTSf1]
		 
    
def plotVPSD(day='2016.10.01',run='00',noise='00',Chan = 1, bs = 1,logplot = 1, fmax = 100, tracelabel = '',scale=1):
	'''plots the time series and PSD of a noise run from the main magnetometer program. Uses the raw Y_noise data, so no magnetometer response calibration is performed. The scale variable allows for comparisons of datasets with differing I-V settings, etc.'''
	a=get1ChanData(day,run,noise,Chan = Chan,ver = 'v16',calib = 0);
	v=a[0]
	sz= size(v);
	dcv=mean(v)
	time=arange(0,sz/1000,.001); #assumes 1000 Hz acquisition rate.
	
	v=v-dcv #AC couples data
	
	Pxx,freq = mlab.psd(v,NFFT = sz,Fs = 1000.); #Calculates power spectrum of unfiltered, AC-coupled data
	
	Px = sqrt(Pxx)*scale; #calculates PSD, scales output so it's possible to compare to other data sets.
	
	freq = binit(freq,bs); Px = binit(Px,bs); #bins data
	
	matplotlib.rcParams['figure.figsize'] = (12.0,9.0)
	
	TSplot=subplot(211);
	plot(time,v,label=tracelabel);
	grid('on');
	xlabel('time [s]');
	ylabel('Signal [V]');
	legend();
	
	
	PSDplot=subplot(212);
	semilogy(freq,Px,label = tracelabel);
	grid('on',which='both');
	xlabel('freq,Hz');ylabel('PSD [1/rHz]');
	xlim(xmax = fmax);
	#ylim(ymax = ymax,ymin = 1E-5);
	if tracelabel !='':
		legend();
	
	return [time,v,freq,Px]

def plotCurrPSD(day='2016.10.01',run='00',noise='00',Chan = 1, bs = 1,calib=2E-6,fmax = 100,tracelabel = ''):
	'''plots the time series and PSD of a noise run from the main magnetometer program. Uses the raw Y_noise data, so no magnetometer response calibration is performed.'''
	a=get1ChanData(day,run,noise,Chan = Chan,ver = 'v16',calib = 0);
	v=a[0]
	I=v*calib
	sz=size(I);
	dcI=mean(I)
	time=arange(0,sz/1000,.001); #assumes 1000 Hz acquisition rate.
	
	I=I-dcI #AC couples data
	
	Pxx,freq = mlab.psd(I,NFFT = sz,Fs = 1000.); #Calculates power spectrum of unfiltered, AC-coupled data
	
	Px = sqrt(Pxx); #calculates PSD, scales output so it's possible to compare to other data sets.
	
	freq = binit(freq,bs); Px = binit(Px,bs); #bins data
	
	matplotlib.rcParams['figure.figsize'] = (12.0,9.0)
	
	TSplot=subplot(211);
	plot(time,1E6*I,label=tracelabel);
	grid("on");
	xlabel('time [s]');
	ylabel('Signal [uA]');
	#legend();
	
	
	PSDplot=subplot(212);
	semilogy(freq,Px,label = tracelabel);
	grid("on",which='both');
	xlabel('freq,Hz');ylabel('PSD [A/rHz]');
	xlim(xmax = fmax);
	if tracelabel !='':
		legend();
	
	return [time,I,freq,Px]

def Fluxgate(day='2016.10.01',run='00',noise='00',Chan = 1, bs = 1,direc='x',logplot = 1, fmax = 100, tracelabel = '',scale=1):
	'''Analysis for data collected from the fluxgate by the main program. Uses the raw voltage data and converts to field using the calibratio that Bob Wyllie worked out.'''
	a=get1ChanData(day,run,noise,Chan = Chan,ver = 'v16',calib = 0);
	v=a[0]
	
	if direc=='x':
		calib=.723909
	elif direc=='y':
		calib=.718403
	elif direc=='z':
		calib=.743970
	
	B=v*1000/calib/143*1000. #converts volts to millivolts, uses calibration calculated by Bob
							 #The factor of 1000 converts uT to nT
	
	sz= size(B);
	dcB=mean(B)
	time=arange(0,sz/1000,.001); #assumes 1000 Hz acquisition rate.
	
	B=B-dcB #AC couples data
	
	Pxx,freq = mlab.psd(B,NFFT = sz,Fs = 1000.); #Calculates power spectrum of unfiltered, AC-coupled data
	
	Px = sqrt(Pxx)*scale; #calculates PSD, scales output so it's possible to compare to other data sets.
	
	freq = binit(freq,bs); Px = binit(Px,bs); #bins data
	
	matplotlib.rcParams['figure.figsize'] = (12.0,9.0)
	
	TSplot=subplot(211);
	plot(time,B,label=tracelabel);
	grid('on');
	xlabel('time [s]');
	ylabel('Signal [nT]');
	legend();
	
	
	PSDplot=subplot(212);
	semilogy(freq,Px,label = tracelabel);
	grid('on',which='both');
	xlabel('freq,Hz');ylabel('PSD [nT/rHz]');
	xlim(xmax = fmax);
	#ylim(ymax = ymax,ymin = 1E-5);
	if tracelabel !='':
		legend();
	
	return [time,B,freq,Px]
	
def CoilNoise(day='2016.10.01',runNum='00',noiseNum='00',bs=1,Chan=1,calib=0,SRScalib=5E-6,dBdI=43E9,fmax=100,label=''):
		
	Vt = get1ChanData(day,runNum,noiseNum,Chan=Chan,calib=0); #Imports raw voltage time-series V(t) from specified file
	sz = size(Vt[0]);
	Vxx,freq = mlab.psd(Vt[0],NFFT = sz,Fs = 1000.); #Creates the power spectrum in V^2/Hz
	Vx = sqrt(Vxx); #computes psd in V/rHz

	freq = binit(freq,bs); Vx = binit(Vx,bs); #bins frequency and psd
	
	Ix=Vx*SRScalib; #Converts the voltage psd to a current psd
	Bx=Ix*dBdI; #converts to magnetic field psd [fT/rHz] given the coil calibration dB/dI
	
	matplotlib.rcParams['figure.figsize'] = (12.0,9.0)
	
	#sets label on trace, gives default label name if none provided
	if label=='':
		tracelabel='Run '+runNum+', Noise '+noiseNum+', Chan '+str(Chan)
	else:
		tracelabel=label
		
	#plots current psd
	Iplot=subplot(211);
	semilogy(freq,Ix,label=tracelabel);
	grid('on',which='both')
	ylabel('Current PSD [A/rHz]');
	xlim(xmax=fmax);
	ylim(ymax=10*max(Ix));
	legend();
	
	#plots magnetic field psd
	Bplot=subplot(212);
	semilogy(freq,Bx)
	grid('on',which='both');
	xlabel('Frequency [Hz]'); ylabel('Mag. Field PSD [fT/rHz]');
	xlim(xmax=fmax);
	ylim(ymax=10*max(Bx));
	legend();

def NullFields(x=0,y=0,z=0,Rx=1000,Ry=1000,Rz=1000,Supply='Sulai',Coils='Shell'):
	if Coils=='Shell':
		dBdIx=21000; dBdIy=43000; dBdIz=21000; #Shell coil calibrations in nT/A
	elif Coils=='Room':
		dBdIx=5160; dBdIy=8150; dBdIz=5220; #Room coil calibrations (@DC) in nT/A

	if Supply=='Sulai':
		Vx=5*x/2**15; #normalizes 15 bit output value, computes output voltage
		Vy=5*y/2**15;
		Vz=5*z/2**15;
	elif Supply=='Wyllie':
		Vx=x/2; #When read from DMM, mon voltage is actually twice the voltage output to the coil.
		Vy=y/2;
		Vz=z/2;
		
	Bx=round(Vx/Rx*dBdIx,1);
	By=round(Vy/Ry*dBdIy,1);
	Bz=round(Vz/Rz*dBdIz,1);
	Btot=round(sqrt(Bx**2+By**2+Bz**2),1)
	
	print(Coils,'coil fields:')
	print('[Bx, By, Bz] =',[Bx,By,Bz],'nT')
	print('Btotal =',Btot,'nT')
	
	
def getRoomBNoise(fdata,calib, bcalib = 1):
    v2 = fdata;
    Curr2 = calib*v2  #Amps

    P2,freq = mlab.psd(Curr2, NFFT = size(Curr2),Fs = 1000); P2a = sqrt(P2);

    bs = 10
    f2 = binit(freq,bs);
    P2Amp = binit(P2a,bs);
    subplot(211)
    semilogy(f2,P2Amp);
    grid('on', which = 'both');
    xlabel('freq,Hz');
    ylabel('A/rHz');
    #ylim(ymax = 1E-6);
    xlim(xmin = 0.1, xmax = 100)
    

    Bnoise = bcalib*P2Amp
    
    subplot(212)
    semilogy(f2,Bnoise, 'r-');
    grid('on', which = 'both');
    xlabel('freq,Hz');
    ylabel('fT/rHz');
    #ylim(ymax = 100);
    xlim(xmin = 0.1, xmax = 100)

def getCh1Ch2XCorr(day,run,noise,calib = 0,LPF = 100., HPF = 1.,bs =1,rScale = 1, p2 = 10):
    

    d1 = get1ChanData(day,run,noise,calib = calib, Chan = 1);
    d2 = get1ChanData(day,run,noise,calib = calib, Chan = 2);

    ch1,ch2 = [d1[0],d2[0]]

    Nyq = 500.
    bLP,aLP = butter(4,[HPF/Nyq, LPF/Nyq], btype = 'bandpass');
    C1 = filtfilt(bLP,aLP,ch1);
    C2 = filtfilt(bLP,aLP,ch2);
    
    sz = size(C1); sz2 = sz-500
    C1,C2 = [C1[500:sz2],C2[500:sz2]]
    
    data = rScale*array([C1,C2])
    w,v=eig(cov(data))
    #print 'The eigenvalues are %.2f and %.2f' % (w[0],w[1]) #eigenvalues
    #print 'The eigenvectors are' 
    #print v #eigenvectors 

    a=dot(inv(v),data)



    sz = 2**p2;

    P1a,freq = mlab.psd(data[0], NFFT =sz, Fs = 1000); P1a = sqrt(P1a)
    P1b,freq = mlab.psd(data[1], NFFT =sz, Fs = 1000); P1b = sqrt(P1b)

    P2a,freq = mlab.psd(a[0], NFFT =sz, Fs = 1000); P2a = sqrt(P2a)
    P2b,freq = mlab.psd(a[1], NFFT =sz, Fs = 1000); P2b = sqrt(P2b)


    fmax = 2*LPF

    P1a = binit(P1a,bs);
    P1b = binit(P1b,bs);
    freq = binit(freq,bs);


    P2a = binit(P2a,bs);
    P2b = binit(P2b,bs);



    figure(1)
    subplot(211);
    plot(ch1,'b-');
    oset = 2*max(d1[0]);
    plot(ch2+oset,'g-')
    grid('on');
    xlabel('time, ms');
    ylabel('field, fT');


    subplot(212)
    semilogy(freq,P1a);
    semilogy(freq,P1b);
    grid('on',which = 'both');
    xlim(xmax = fmax)
    grid('on');
    suptitle('Raw Data')
    xlabel('freq,Hz');
    ylabel('PSD, fT/rHz');
    ylim(ymin = .1);
    

    #CREATING SPACE BETWEEN FIGURES
    figure()
    subplot(811)
    axis('off')
    show();
    #########################

    figure(2)
    subplot(211);
    plot(a[0],'b-');

    plot(oset + a[1],'g-')
    grid('on');
    xlabel('time, ms');
    ylabel('field, fT');


    subplot(212)
    semilogy(freq,P2a);
    semilogy(freq,P2b);
    grid('on',which = 'both');
    xlim(xmax = fmax)
    suptitle('X-Correlated')
    grid('on');
    xlabel('freq,Hz');
    ylabel('PSD, fT/rHz');
    ylim(ymin = .1);

    #CREATING SPACE BETWEEN FIGURES
    figure()
    subplot(811)
    axis('off')
    show();
    #########################


    PratioA = P2b/P2a;
    PratioB = P2a/P2b;
    
    if mean(PratioA)>mean(PratioB):
	    Pratio = PratioA;
    else:
	    Pratio = PratioB;

    figure(3);
    subplot(313)
    plot(freq,Pratio, 'r-');
    grid('on',which = 'both');
    xlim(xmax = fmax)
    #ylim(ymax = 10)
    xlabel('freq,Hz');
    ylabel('PSD ratio')
    show();


    
def getCh1Ch2PCA(day,run,noise,LPF = 100., HPF = 1.,bs =1,rScale = 1, p2 = 10, calib = 1, note1 = '', note2 = '', tmsize = 5000,Ch = [1,2], Fsamp = 1000):
    
	Nyq = Fsamp/2.0;
	
	
	d1 = get1ChanData(day,run,noise,calib = 1, Chan = Ch[0]);
	d2 = get1ChanData(day,run,noise,calib = 1, Chan = Ch[1]);

	ch1,ch2 = [d1[0],d2[0]]


	bLP,aLP = butter(2,[HPF/Nyq, LPF/Nyq], btype = 'bandpass');
	C1 = filtfilt(bLP,aLP,ch1);
	C2 = filtfilt(bLP,aLP,ch2);

	
	trimFactor = Fsamp/2.
	sz = size(C1); sz2 = sz-trimFactor
	C1,C2 = [C1[trimFactor:sz2],C2[trimFactor:sz2]]

	data = rScale*array([C1,C2])
	w,v=eig(cov(data))
	#print 'The eigenvalues are %.2f and %.2f' % (w[0],w[1]) #eigenvalues
	#print 'The eigenvectors are' 
	#print v #eigenvectors 

	a=dot(inv(v),data)



	sz = 2**p2;

	P1a,freq = mlab.psd(data[0], NFFT =sz, Fs = Fsamp); P1a = sqrt(P1a)
	P1b,freq = mlab.psd(data[1], NFFT =sz, Fs = Fsamp); P1b = sqrt(P1b)

	P2a,freq = mlab.psd(a[0], NFFT =sz, Fs = Fsamp); P2a = sqrt(P2a)
	P2b,freq = mlab.psd(a[1], NFFT =sz, Fs = Fsamp); P2b = sqrt(P2b)


	fmax = 2*LPF

	P1a = binit(P1a,bs);
	P1b = binit(P1b,bs);
	freq = binit(freq,bs);


	P2a = binit(P2a,bs);
	P2b = binit(P2b,bs);

	ln1 = tmsize;

	figure(1, figsize=[14,10])

	subplots_adjust(left=None, bottom=None, right=None, top=None, wspace=None, hspace=0.32)


	subplot(621);
	plot(1e-3*data[0],'b-');
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

	subplot(323)
	semilogy(freq,P1a);
	semilogy(freq,P1b);
	grid('on',which = 'both');
	xlim(xmax = fmax)
	grid('on');

	xlabel('freq,Hz');
	ylabel('PSD, fT/rHz');
	ylim(ymin = .1);

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





	subplot(324)
	semilogy(freq,P2a);
	semilogy(freq,P2b);
	grid('on',which = 'both');
	xlim(xmax = fmax)
	grid('on');
	xlabel('freq,Hz');
	ylabel('PSD, fT/rHz');
	ylim(ymin = .1);

	PratioA = P2b/P2a;
	PratioB = P2a/P2b;

	if mean(PratioA)>mean(PratioB):
		Pratio = PratioA;
	else:
		Pratio = PratioB;


	subplot(326)
	plot(freq,Pratio, 'r-');
	grid('on',which = 'both');
	xlim(xmax = fmax)
	ylim(ymax = 1000)
	xlabel('freq,Hz');
	ylabel('PSD ratio')

	fignotes = subplot(325)
	fignotes.axis('off');
	fs = 16
	frame1 = gca()
	frame1.axes.get_xaxis().set_visible(False)
	frame1.axes.get_yaxis().set_visible(False)

	bptext = 'Band Pass: ' + str(HPF) +' - '+str(LPF) +  ' Hz';1.0
	evtext = 'EigenVector = '+str(v[0]);
	
	
	text(0.1,0.8,bptext, fontsize = fs)
	text(0.1,0.5,note1, fontsize = fs)
	text(0.1,0.2,note2, fontsize = fs)
	text(0.1,0.0,evtext, fontsize = 12)
	
	figtitle = day+'Run '+run+' Noise '+noise;
	suptitle(figtitle);

	show();
	return [data[0], data[1]]


#Functions for displaying the magnetic fields applied by Compensating Fields on Shells
	
def ChanVal4(V):
    print(round(V[0],2),round(V[1],2),round(V[2],2));
    print(round(V[3],2),round(V[4],2),round(V[5],2));
    print(round(V[6],2),round(V[7],2),round(V[8],2));
    print(round(V[9],2),round(V[10],2),round(V[11],2));


    
    
def getBout(EncVals,OutResVals):

    Bcal = 1E-6*array([17,25.5,17,17,25.5,17,17,25.5,17,17,25.5,17])  #Tesla/amp

    Vout = 5*array(EncVals)/(2**15);
    Iout = Vout/OutResVals;  #Amps
    Bout = 1E9*Bcal*Iout

    ChanVal4(Bout)
	

def getCh1Ch2GRA(day,run,noise,calib = 1,LPF = 100., HPF = 1.,bs =1,rScale = 1, p2 = 10, note1 = '', note2 = '', tmsize = 5000,gmode = [0.707,0.707],Ch = [1,2]):
    
	
	d1 = get1ChanData(day,run,noise,calib = calib, Chan = Ch[0]);
	d2 = get1ChanData(day,run,noise,calib = calib, Chan = Ch[1]);

	ch1,ch2 = [d1[0],d2[0]]

	Nyq = 500.
	bLP,aLP = butter(4,[HPF/Nyq, LPF/Nyq], btype = 'bandpass');
	C1 = filtfilt(bLP,aLP,ch1);
	C2 = filtfilt(bLP,aLP,ch2);

	sz = size(C1); sz2 = sz-500
	C1,C2 = [C1[500:sz2],C2[500:sz2]]

	data = rScale*array([C1,C2])
	w,v=eig(cov(data))
	#print 'The eigenvalues are %.2f and %.2f' % (w[0],w[1]) #eigenvalues
	#print 'The eigenvectors are' 
	#print v #eigenvectors 
	g1 = gmode[0]; g2 = gmode[1]
	gMatrix = array([[g1,g2],[g2,-g1]])
	v = gMatrix;
	a=dot(inv(gMatrix),data)



	sz = 2**p2;

	P1a,freq = mlab.psd(data[0], NFFT =sz, Fs = 1000); P1a = sqrt(P1a)
	P1b,freq = mlab.psd(data[1], NFFT =sz, Fs = 1000); P1b = sqrt(P1b)

	P2a,freq = mlab.psd(a[0], NFFT =sz, Fs = 1000); P2a = sqrt(P2a)
	P2b,freq = mlab.psd(a[1], NFFT =sz, Fs = 1000); P2b = sqrt(P2b)


	fmax = 2*LPF

	P1a = binit(P1a,bs);
	P1b = binit(P1b,bs);
	freq = binit(freq,bs);


	P2a = binit(P2a,bs);
	P2b = binit(P2b,bs);

	ln1 = tmsize;

	figure(1, figsize=[14,10])

	subplots_adjust(left=None, bottom=None, right=None, top=None, wspace=None, hspace=0.32)


	subplot(621);
	plot(1e-3*data[0],'b-');
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

	subplot(323)
	semilogy(freq,P1a);
	semilogy(freq,P1b);
	grid('on',which = 'both');
	xlim(xmax = fmax)
	grid('on');

	xlabel('freq,Hz');
	ylabel('PSD, fT/rHz');
	ylim(ymin = .1);

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





	subplot(324)
	semilogy(freq,P2a);
	semilogy(freq,P2b);
	grid('on',which = 'both');
	xlim(xmax = fmax)
	grid('on');
	xlabel('freq,Hz');
	ylabel('PSD, fT/rHz');
	ylim(ymin = .1);

	PratioA = P2b/P2a;
	PratioB = P2a/P2b;

	if mean(PratioA)>mean(PratioB):
		Pratio = PratioA;
	else:
		Pratio = PratioB;


	subplot(326)
	plot(freq,Pratio, 'r-');
	grid('on',which = 'both');
	xlim(xmax = fmax)
	#ylim(ymax = 10)
	xlabel('freq,Hz');
	ylabel('PSD ratio')

	fignotes = subplot(325)
	fignotes.axis('off');
	fs = 16
	frame1 = gca()
	frame1.axes.get_xaxis().set_visible(False)
	frame1.axes.get_yaxis().set_visible(False)

	bptext = 'Band Pass: ' + str(HPF) +' - '+str(LPF) +  ' Hz';1.0
	evtext = 'EigenVector = '+str(v[0]);
	
	
	text(0.1,0.8,bptext, fontsize = fs)
	text(0.1,0.5,note1, fontsize = fs)
	text(0.1,0.2,note2, fontsize = fs)
	text(0.1,0.0,evtext, fontsize = 12)
	
	figtitle = day+'Run '+run+' Noise '+noise;
	suptitle(figtitle);

	show();
	
def MagRelRes(x,a1,b1,c1,a2,b2):
	return a1/sqrt(b1**2+(x-c1)**2)+a2/sqrt(b2**2+(x-c2)**2);
	#return a1+b1*x+c1*x**2


def MagRelResEven(x, a1, b1,c1):
    if x>=0:
        return a1+b1*x+c1*x**2
    else:
        return a1-b1*x+c1*x**2
    
def MagPhRel(x, a1, b1):
    return -a1*arctan(b1*x)
	
def MagPhRelOdd(x, a1, b1, c1):
    if x>=0:
        return a1+b1*x+c1*x**2
    else:
        return a1+b1*x-c1*x**2
    
    
def MagPhRelEven(x, a1, b1, c1):
    if x>=0:
        return a1+b1*x+c1*x**2
    else:
        return a1-b1*x+c1*x**2


def getCh1Ch2PCAv2(ch1,ch2,LPF = 100., HPF = 1.,bs =1,rScale = 1, note1 = '', note2 = '', tmsize = 1000,Fsamp = 1000,minY=10,fmax=100,cmMax = 1000):
	Nyq = Fsamp/2.0;

	bLP,aLP = butter(2,[HPF/Nyq, LPF/Nyq], btype = 'bandpass');
	C1 = filtfilt(bLP,aLP,ch1);
	C2 = filtfilt(bLP,aLP,ch2);


	trimFactor = Fsamp/2.
	sz = size(C1); sz2 = sz-trimFactor
	C1,C2 = [C1[trimFactor:sz2],C2[trimFactor:sz2]]

	data = rScale*array([C1,C2])
	w,v=eig(cov(data))
	#print 'The eigenvalues are %.2f and %.2f' % (w[0],w[1]) #eigenvalues
	#print 'The eigenvectors are' 
	#print v #eigenvectors 

	a=dot(inv(v),data)



	sz = size(C1)

	P1a,freq = mlab.psd(data[0], NFFT =sz, Fs = Fsamp); P1a = sqrt(P1a)
	P1b,freq = mlab.psd(data[1], NFFT =sz, Fs = Fsamp); P1b = sqrt(P1b)

	P2a,freq = mlab.psd(a[0], NFFT =sz, Fs = Fsamp); P2a = sqrt(P2a)
	P2b,freq = mlab.psd(a[1], NFFT =sz, Fs = Fsamp); P2b = sqrt(P2b)




	P1a = binit(P1a,bs);
	P1b = binit(P1b,bs);
	freq = binit(freq,bs);


	P2a = binit(P2a,bs);
	P2b = binit(P2b,bs);

	ln1 = tmsize;

	figure(1, figsize=[14,10])

	subplots_adjust(left=None, bottom=None, right=None, top=None, wspace=None, hspace=0.32)


	subplot(621);
	plot(1e-3*data[0],'b-');
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





	subplot(324)
	semilogy(freq,P2a);
	semilogy(freq,P2b);
	grid('on',which = 'both');
	xlim(xmax = fmax)
	grid('on');
	xlabel('freq,Hz');
	ylabel('PSD, fT/rHz');
	ylim(ymin = minY);

	PratioA = P2b/P2a;
	PratioB = P2a/P2b;

	if mean(PratioA[0:20])>mean(PratioB[0:20]):
		Pratio = PratioA;
	else:
		Pratio = PratioB;


	subplot(326)
	plot(freq,Pratio, 'r-');
	grid('on',which = 'both');
	xlim(xmax = fmax)
	ylim(ymax = cmMax)
	xlabel('freq,Hz');
	ylabel('PSD ratio')

	fignotes = subplot(325)
	fignotes.axis('off');
	fs = 16
	frame1 = gca()
	frame1.axes.get_xaxis().set_visible(False)
	frame1.axes.get_yaxis().set_visible(False)

	bptext = 'Band Pass: ' + str(HPF) +' - '+str(LPF) +  ' Hz';1.0
	evtext = 'EigenVector = '+str(v[0]);


	text(0.1,0.8,bptext, fontsize = fs)
	text(0.1,0.5,note1, fontsize = fs)
	text(0.1,0.2,note2, fontsize = fs)
	text(0.1,0.0,evtext, fontsize = 12)
		
