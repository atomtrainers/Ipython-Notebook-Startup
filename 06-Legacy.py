import os;
from scipy.signal import *
from scipy.optimize import curve_fit
from numpy import *
from pylab import *
import mdp
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider, Button, RadioButtons
import time

def plotMagRun(day,runNum,noiseNum,noiseNum2,exNoise='00',Chan=1,bs = 1,fignum = 1, ver='v16',rScale=1,rgain = 1,fmax = 200, minpsd = 0.05, figLabel=["Total Noise","Probe Noise","Electronic Noise"]):
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

    
    if noiseNum<10:
        noiseNum = '0'+str(noiseNum);
        
    fname1 = runDir+'\\MagRun2'+str(day)+'_'+'Run'+str(runNum)+'Chan'+str(Chan)
    fname2 = dayDir+'\\MagRun2'+str(day)+'_'+'Run'+str(runNum)+'Chan'+str(Chan)
    

    path = basepath +str(day) + '\\' + 'run_' + str(runNum);
    fPSD = path + '\\'+'noise_'+str(noiseNum)+'\\'+'Y_magnetic-field-PSD_'+str(Chan)+'.txt';
    fPSD2 = path + '\\'+'noise_'+str(noiseNum2)+'\\'+'Y_magnetic-field-PSD_'+str(Chan)+'.txt';
    fPSD3 = path + '\\'+'noise_'+str(exNoise)+'\\'+'Y_magnetic-field-PSD_'+str(Chan)+'.txt';
    fRESP =  path+'\\'+'response'+'\\'+'Y_amp-n-phase_response_'+str(Chan)+'.txt';
    
    if ver == 'v12':
        fPSD = path + '\\'+'noise_'+str(noiseNum)+'\\'+'magnetic-field-PSD_'+str(Chan)+'.txt';
        fPSD2 = path + '\\'+'noise_'+str(noiseNum2)+'\\'+'magnetic-field-PSD_'+str(Chan)+'.txt';
        fPSD3 = path + '\\'+'noise_'+str(exNoise)+'\\'+'magnetic-field-PSD_'+str(Chan)+'.txt';
        fRESP =  path+'\\'+'response'+'\\'+'amp-n-phase_response_'+str(Chan)+'.txt';
    
    
    figTitle = str(day) + '\\' + 'run_' + str(runNum)+'\\'+'noise_'+str(noiseNum) + ': Channel ' + str(Chan);

    data = loadtxt(fPSD);
    data = transpose(data);
    freq = binit(data[0],bs);
    PS = binit(data[1],bs);
    dataResp = loadtxt(fRESP);
    dataResp = transpose(dataResp)

    data = loadtxt(fPSD2);
    data = transpose(data);
    freq = binit(data[0],bs);
    PS2 = binit(data[1],bs);

    data = loadtxt(fPSD3);
    data = transpose(data);
    freq = binit(data[0],bs);
    PS3 = binit(data[1],bs);

    figure(fignum)
    clf();
    subplot(121)
    semilogy(freq,PS,'k-',label=figLabel[0]);
    semilogy(freq,PS2,'r-',label=figLabel[1]);
    if exNoise!='00':
	    semilogy(freq,PS3,'b-',label=figLabel[2]);
    grid(which='both');
    axis([0,fmax,minpsd,10000])
    xlabel('frequency, Hz');
    ylabel('PSD, fT/rHz');
    legend();


    axR = subplot(122);
    plot(dataResp[0],1E5*rgain*dataResp[1],'k-');
    axis([0,fmax,0,rScale]);
    xlabel('frequency, Hz');
    axR.yaxis.set_label_position("right")
    ylabel('Resp, 1E-5*V/fT');
    grid(which='both');
    suptitle(figTitle)
    #show();
    savefig(fname1+'.pdf', format = 'pdf')
    savefig(fname2+'.png')

def plot4LaserScan(date = '', laser = 'Probe', fnum=0, fignum = 1):
	if fnum<10:
		fnum = '000'+str(fnum);
	elif fnum<100:
		fnum = '00'+str(fnum);
	elif fnum<1000:
		fnum = '0'+str(fnum);
	
	[basepath, analysispath] = getCompEnv();
	date2=date[0:4]+'.'+date[4:6]+'.'+date[6:8]
	dayDir = analysispath+date2
	a = os.path.isdir(dayDir);
	if a==0:
		os.mkdir(dayDir);

	ftit = date+laser+'Scan'+fnum;
	
	datafile = 'C:\\Users\\sulai\\Documents\\LabVIEW Data\\DATA\\LaserScanLog\\'+ftit+'.dat'
	
	data = loadtxt(datafile);
	
	data0p=np.delete(data[0],0)
	data1p=np.delete(data[1],0)
	data2p=np.delete(data[2],0)
	data3p=np.delete(data[3],0)
	data4p=np.delete(data[4],0)
	
	xbot=round(min(data0p),1);
	xtop=round(max(data0p),1);
	
	pngname=dayDir+'\\'+date2+'_'+laser+'Scan'+str(fnum)+'.png'
	
	figure(fignum, figsize = [16,8])
	clf();
	plot(data0p,data1p, 'y-', linewidth = '3',label='Channel 1');
	plot(data0p,data2p, 'b-', linewidth = '3',label='Channel 2');
	plot(data0p,data3p, 'm-', linewidth = '3',label='Channel 3');
	plot(data0p,data4p, 'g-', linewidth = '3',label='Channel 4');
	title(ftit,fontsize=16)
	xlim(xbot,xtop)
	xticks(np.arange(xbot, xtop, 0.2))
	xlabel('TEC Temp ($ \propto \Delta$) [kOhm]',fontsize=14);
	ylabel('Signal [V]',fontsize=14);
	ylim(ymin=-.1)
	grid(which='both');
	legend(loc=0);
	savefig(pngname)
	
	return data0p,data1p,data2p,data3p,data4p;
	
def LaserScanCompare(laser='Probe',day1 = '20160101', fnum1=0,Chan1='1',s1=1,lbl1='-',day2 = '20160101',fnum2=0,Chan2='1',s2=1,lbl2='-',day3='-',fnum3=0,Chan3='1',s3=1,lbl3='-',day4='-',fnum4=0,Chan4='1',s4=1,lbl4='-',day5='-',fnum5=0,Chan5='1',s5=1,lbl5='-',legendloc=0):
	'''Plots multiple runs taken with the LaserScan.vi program. The "s" inputs allow you to scale the data in case you change the applied field or other parameters. By default, only two traces are plotted, though more can be added manually'''
	def AddZeroes(fnum):
		if fnum<10:
			fnum = '000'+str(fnum);
		elif fnum<100:
			fnum = '00'+str(fnum);
		elif fnum<1000:
			fnum = '0'+str(fnum);
		return fnum
	
	fnum1=AddZeroes(fnum1);
	fnum2=AddZeroes(fnum2);
	fnum3=AddZeroes(fnum3);
	fnum4=AddZeroes(fnum4);
	fnum5=AddZeroes(fnum5);
	
	def LabelMaker(day,laser,fnum,lbl):
		if lbl=='-':
			label=day+laser+'Scan'+str(fnum)
		else:
			label=lbl
		return label

	label1=LabelMaker(day1,laser,fnum1,lbl1);
	label2=LabelMaker(day2,laser,fnum2,lbl2);
	label3=LabelMaker(day3,laser,fnum3,lbl3);
	label4=LabelMaker(day4,laser,fnum4,lbl4);
	label5=LabelMaker(day5,laser,fnum5,lbl5);
	
	def LoadLaserScan(day,fnum,Chan):
		file = day+laser+'Scan'+str(fnum);
		path = 'C:\\Users\\sulai\\Documents\\LabVIEW Data\\DATA\\LaserScanLog\\'+file+'.dat';
		data=loadtxt(path);
		detun=np.delete(data[0],0);
		sig=np.delete(data[Chan],0);
		
		return detun,sig
	
	figure(1, figsize = [16,10])
	clf();
	
	[detun1,sig1]=LoadLaserScan(day1,fnum1,Chan1);
	plot(detun1,s1*sig1,'b-',linewidth=3,label=label1);
	
	[detun2,sig2]=LoadLaserScan(day2,fnum2,Chan2);
	plot(detun2,s2*sig2,'g-',linewidth=3,label=label2);
	
	if day3 != '-':
		[detun3,sig3]=LoadLaserScan(day3,fnum3,Chan3);
		plot(detun3,s3*sig3,'r-',linewidth=3,label=label3);
		
	if day4 != '-':
		[detun4,sig4]=LoadLaserScan(day4,fnum4,Chan4);
		plot(detun4,s4*sig4,'c-',linewidth=3,label=label4);
		
	if day5 != '-':
		[detun5,sig5]=LoadLaserScan(day5,fnum5,Chan5);
		plot(detun5,s5*sig5,'y-',linewidth=3,label=label5);
	
	xlim(xmin=min(detun1)-.1,xmax=max(detun1)+.1)
	xlabel('TEC Temp ($ \propto \Delta$) [kOhm]',fontsize=14);
	ylabel('Signal [V]',fontsize=14);
	ylim(ymin = -.1)
	grid(which='both');
	legend(loc=legendloc);
	show()
	
def plotMagRun4(day='2016.01.01',runNum='00',noiseNum = '00',noiseNum2 = '01', noiseNum3 = '02',bs = 4,fmax = 200,fignum = 1,ver = 'v16',leg=1,rScale=50,direc = 'Y'):
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


	
	figs=figure(figsize = (16,10))
	clf();
	subplot(231)
	semilogy(freq,PS,'k-',label='Total Noise');
	if noiseNum2 != noiseNum:
		semilogy(freq,PS2,'r-',label='Probe Noise');
	if noiseNum3 != noiseNum:
		semilogy(freq,PS3,'b-',label='Electronic Noise');
	grid(which='both');
	axis([0,fmax,0.1,10000])
	# xlabel('frequency, Hz');
	ylabel('PSD, fT/rHz');
	#~ legend();
	#~ legend();
	text(25,2000,'Ch1',fontsize=fs)

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



	subplot(232)
	semilogy(freq,PS,'k-',label='Total Noise');
	if noiseNum2 != noiseNum:
		semilogy(freq,PS2,'r-',label='Probe Noise');
	if noiseNum3 != noiseNum:
		semilogy(freq,PS3,'b-',label='Electronic Noise');
	grid(which='both');
	axis([0,fmax,0.1,10000])
	# xlabel('frequency, Hz');
	# ylabel('PSD, fT/rHz');

	if leg==1:
		legend();
	text(25,2000,'Ch2',fontsize=fs)


	ns2 = PS;


	Chan = 3;
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
	dataResp3 = transpose(dataResp)

	data = loadtxt(fPSD2);
	data = transpose(data);
	freq = binit(data[0],bs);
	PS2 = binit(data[1],bs);

	data = loadtxt(fPSD3);
	data = transpose(data);
	freq = binit(data[0],bs);
	PS3 = binit(data[1],bs);


	subplot(234)
	semilogy(freq,PS,'k-',label='Total Noise');
	if noiseNum2 != noiseNum:
		semilogy(freq,PS2,'r-',label='Probe Noise');
	if noiseNum3 != noiseNum:
		semilogy(freq,PS3,'b-',label='Electronic Noise');
	grid(which='both');
	axis([0,fmax,0.1,10000])
	xlabel('frequency, Hz');
	ylabel('PSD, fT/rHz');
	text(25,2000,'Ch3',fontsize=fs)



	ns3 = PS;

	Chan = 4;
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
	dataResp4 = transpose(dataResp)

	data = loadtxt(fPSD2);
	data = transpose(data);
	freq = binit(data[0],bs);
	PS2 = binit(data[1],bs);

	data = loadtxt(fPSD3);
	data = transpose(data);
	freq = binit(data[0],bs);
	PS3 = binit(data[1],bs);



	subplot(235)
	semilogy(freq,PS,'k-',label='Total Noise');
	if noiseNum2 != noiseNum:
		semilogy(freq,PS2,'r-',label='Probe Noise');
	if noiseNum3 != noiseNum:
		semilogy(freq,PS3,'b-',label='Electronic Noise');
	grid(which='both');
	axis([0,fmax,0.1,10000])
	xlabel('frequency, Hz');
	# ylabel('PSD, fT/rHz');
	#~ legend();
	text(25,2000,'Ch4',fontsize=fs)

	ns4 = PS;





	xR = subplot(233);
	plot(dataResp1[0],1E6*dataResp1[1],'y-', linewidth = 3, label = 'Ch1');
	plot(dataResp2[0],1E6*dataResp2[1],'b-', linewidth = 3, label = 'Ch2');
	plot(dataResp3[0],1E6*dataResp3[1],'m-', linewidth = 3, label = 'Ch3');
	plot(dataResp4[0],1E6*dataResp4[1],'g-', linewidth = 3, label = 'Ch4');
	axis([0,fmax,0,rScale]);
	xR.yaxis.set_label_position("right")
	ylabel('Resp, $\mu$V/fT');
	grid(which='both');
	legend()

	xR2 = subplot(236);
	plot(dataResp1[0],dataResp1[2] - dataResp1[2][0],'y-', linewidth = 3, label = 'Ch1');
	plot(dataResp2[0],dataResp2[2]-dataResp2[2][0],'b-', linewidth = 3, label = 'Ch2');
	plot(dataResp3[0],dataResp3[2]- dataResp3[2][0],'m-', linewidth = 3, label = 'Ch3');
	plot(dataResp4[0],dataResp4[2]-dataResp4[2][0],'g-', linewidth = 3, label = 'Ch4');
	axis([0,fmax,-180,180]);
	xR2.yaxis.set_label_position("right")
	ylabel('Phase (deg)');
	grid(which='both');
	#~ legend()



	figTitle = str(day) + '\\' + 'run_' + str(runNum)+ ' B-'+str(direc)+':Channels 1 - 4';
	fname1 = dayDir+'\\Mag4'+str(direc)+str(day)+'_'+'Run'+str(runNum);
	fname2 = runDir+'\\Mag4'+str(direc)+str(day)+'_'+'Run'+str(runNum);
	suptitle(figTitle);
	print(fname1);
	show();
	figs.savefig(fname2+'.pdf', format = 'pdf');
	figs.savefig(fname1+'.png');
	
def plotMag1Chan(day,runNum,noiseNum,Chan=1,bs = 1,fignum = 1, ver='v16',rScale=1,rgain = 1,fmax = 200, psn = 0,cfit = 0,calib = 500E-9,Ipr = 100E-6,fitfunc=0):
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


    if noiseNum<10:
        noiseNum = '0'+str(noiseNum);

    fname1 = runDir+'\\MagRun2PS'+str(day)+'_'+'Run'+str(runNum)+'Chan'+str(Chan)
    fname2 = dayDir+'\\MagRun2PS'+str(day)+'_'+'Run'+str(runNum)+'Chan'+str(Chan)


    path = basepath +str(day) + '\\' + 'run_' + str(runNum);
    fPSD = path + '\\'+'noise_'+str(noiseNum)+'\\'+'Y_magnetic-field-PSD_'+str(Chan)+'.txt';
    fRESP =  path+'\\'+'response'+'\\'+'Y_amp-n-phase_response_'+str(Chan)+'.txt';
    
    if psn==1:
        figTitle = str(day) + '\\' + 'run_' + str(runNum)+'\\'+'noise_'+str(noiseNum) + ': Channel ' + str(Chan) + '  (' + str(1E6*Ipr) + ' uAmp per PD)';
    if psn==0:
        figTitle = str(day) + '\\' + 'run_' + str(runNum)+'\\'+'noise_'+str(noiseNum) + ': Channel ' + str(Chan) ;

    data = loadtxt(fPSD);
    data = transpose(data);
    freq = binit(data[0],bs);
    PS = binit(data[1],bs);
    dataResp = loadtxt(fRESP);
    dataResp = transpose(dataResp)

    shotnoiseB=0;



    fig = figure(fignum)
    clf();

    subplot(121)
    semilogy(freq,PS,'k-',label='Total Noise');
    grid(which='both');
    axis([0,fmax,0.05,10000])
    xlabel('frequency (Hz)');
    ylabel('PSD (fT/rHz)');
    legend();




    ax1 = subplot(122);
    plot(dataResp[0],1E5*rgain*dataResp[1],'k-', label = 'meas.');
    xlim(xmax = fmax);
    xlabel('frequency (Hz)');
    grid(which='both');
    suptitle(figTitle)

    phSc = 1E-5*calib/(4*Ipr);

    fig.subplots_adjust(wspace = 0.32)
    
    pfit = [0,0,0,0];
    if psn ==1:
        
        dx,dy = dataResp[0],1E5*dataResp[1];

        
        shotnoiseA = sqrt(4*Ipr*1.6E-19);   #measured in A/rHz
        shotnoiseV = shotnoiseA/calib    #measured in V/rHz
        
        
        shotnoiseB = shotnoiseV/dataResp[1];  #measured in fT/rHz
        
        if cfit == 1:
            if fitfunc==0:
                pfit,perr = curve_fit(func,dx,dy, p0 = [rScale,20.,rScale/3,30.]);
                fa,fb,fc,fd = pfit;

                yfit = func(dx,fa,fb,fc,fd)
                shotnoiseB = 1E5*shotnoiseV/yfit;  #measured in fT/rHz
            
            if fitfunc==1:
                pfit,perr = curve_fit(func1,dx,dy);
                fa,fb = pfit;

                yfit = func1(dx,fa,fb)
                shotnoiseB = 1E5*shotnoiseV/yfit;  #measured in fT/rHz
            
            
            figure(fignum);
            subplot(122)
            plot(dx,yfit,'r-', label = 'fit');
            legend()
            
        a = ax1.axis();

        ax2 = ax1.twinx()
        ax2.axis([a[0],a[1],a[2],1E7*phSc*a[3]])
        ax1.set_ylabel('Mag. Resp. (1E-5*V/fT)');
        ax2.set_ylabel('Mag. Resp. (1E-7*rad/fT)');

             
        
        
        subplot(121);
        semilogy(dataResp[0],shotnoiseB,'g-', linewidth = 3,label='photon shot noise');
        legend();
        
        
        
        show();

    fig.savefig(fname1+'.pdf', format = 'pdf');
    fig.savefig(fname2+'.png');
    return [freq,PS,dataResp[0],1E5*rgain*dataResp[1]];

def plotMagRun4a(day,runNum,noiseNum = '00',bs = 2,fmax = 200,fignum = 1,ver = 'v16',leg=1,rScale=5,direc = 'Y', el = 1):
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

	dataResp1 = [array([0,0]),array([0,0]),array([0,0])];
	dataResp2 =[array([0,0]),array([0,0]),array([0,0])];
	dataResp3 =[array([0,0]),array([0,0]),array([0,0])];
	dataResp4 = [array([0,0]),array([0,0]),array([0,0])];

	Chan = 1;
	path = basepath +str(day) + '\\' + 'run_' + str(runNum);
	
	fPSD,fRESP,figTitle = getSingleRunSpecs(day,path,noiseNum,Chan,direc,runNum)
	
	if os.path.exists(fPSD)==True:
		data = loadtxt(fPSD);
		data = transpose(data);
		freq = binit(data[0],bs);
		ns1 = binit(data[1],bs);
		dataResp = loadtxt(fRESP);
		dataResp1 = transpose(dataResp)

			
		figs=figure(figsize = (11,8.5))
		clf();
		subplot(231)
		semilogy(freq,ns1,'k-',label='Total Noise');

		grid(which='both');
		axis([0,fmax,0.1,10000])
		# xlabel('frequency, Hz');
		ylabel('PSD, fT/rHz');
		#~ legend();
		#~ legend();
		text(25,2000,'Ch1',fontsize=fs)
		
	elif os.path.exists(fPSD)==False:
		figs=figure(figsize = (11,8.5))
		clf();
		subplot(231)

	Chan = 2;
	path = basepath +str(day) + '\\' + 'run_' + str(runNum);
	
	fPSD,fRESP,figTitle = getSingleRunSpecs(day,path,noiseNum,Chan,direc,runNum)

	if os.path.exists(fPSD):
		data = loadtxt(fPSD);
		data = transpose(data);
		freq = binit(data[0],bs);
		ns2 = binit(data[1],bs);
		dataResp = loadtxt(fRESP);
		dataResp2 = transpose(dataResp)

		
		subplot(232)
		semilogy(freq,ns2,'k-',label='Total Noise');

		grid(which='both');
		axis([0,fmax,0.1,10000])
		text(25,2000,'Ch2',fontsize=fs)
	elif os.path.exists(fPSD)==False:
		subplot(232)


	Chan = 3;
	path = basepath +str(day) + '\\' + 'run_' + str(runNum);
	fPSD,fRESP,figTitle = getSingleRunSpecs(day,path,noiseNum,Chan,direc,runNum)

	if os.path.exists(fPSD):
	
		data = loadtxt(fPSD);
		data = transpose(data);
		freq = binit(data[0],bs);
		ns3 = binit(data[1],bs);
		dataResp = loadtxt(fRESP);
		dataResp3 = transpose(dataResp)

		
		subplot(234)
		semilogy(freq,ns3,'k-',label='Total Noise');

		grid(which='both');
		axis([0,fmax,0.1,10000])
		xlabel('frequency, Hz');
		ylabel('PSD, fT/rHz');
		text(25,2000,'Ch3',fontsize=fs)
	elif os.path.exists(fPSD)==False:
		subplot(234)


	Chan = 4;
	path = basepath +str(day) + '\\' + 'run_' + str(runNum);
	fPSD,fRESP,figTitle = getSingleRunSpecs(day,path,noiseNum,Chan,direc,runNum)
	
	if os.path.exists(fPSD):
			
		data = loadtxt(fPSD);
		data = transpose(data);
		freq = binit(data[0],bs);
		ns4 = binit(data[1],bs);
		dataResp = loadtxt(fRESP);
		dataResp4 = transpose(dataResp)


		subplot(235)
		semilogy(freq,ns4,'k-',label='Total Noise');
		
		grid(which='both');
		axis([0,fmax,0.1,10000])
		xlabel('frequency, Hz');
		# ylabel('PSD, fT/rHz');
		#~ legend();
		text(25,2000,'Ch4',fontsize=fs)
	elif os.path.exists(fPSD)==False:
		subplot(235)

	xR = subplot(233);
	plot(dataResp1[0],1E5*dataResp1[1],'y-', linewidth = 3, label = 'Ch1');
	plot(dataResp2[0],1E5*dataResp2[1],'b-', linewidth = 3, label = 'Ch2');
	plot(dataResp3[0],1E5*dataResp3[1],'m-', linewidth = 3, label = 'Ch3');
	plot(dataResp4[0],1E5*dataResp4[1],'g-', linewidth = 3, label = 'Ch4');
	axis([0,fmax,0,rScale]);
	xR.yaxis.set_label_position("right")
	ylabel('Resp, 1E-5*V/fT');
	grid(which='both');
	legend()

	xR2 = subplot(236);
	plot(dataResp1[0],dataResp1[2] - dataResp1[2][0],'y-', linewidth = 3, label = 'Ch1');
	plot(dataResp2[0],dataResp2[2]-dataResp2[2][0],'b-', linewidth = 3, label = 'Ch2');
	plot(dataResp3[0],dataResp3[2]- dataResp3[2][0],'m-', linewidth = 3, label = 'Ch3');
	plot(dataResp4[0],dataResp4[2]-dataResp4[2][0],'g-', linewidth = 3, label = 'Ch4');
	axis([0,fmax,-180,180]);
	xR2.yaxis.set_label_position("right")
	ylabel('Phase (deg)');
	grid(which='both');
	#~ legend()



	figTitle = str(day) + '\\' + 'run_' + str(runNum)+ ' B-'+str(direc)+':Channels 1 - 4;'
	fname1 = dayDir+'\\Mag4'+str(direc)+str(day)+'_'+'Run'+str(runNum)
	fname2 = runDir+'\\Mag4'+str(direc)+str(day)+'_'+'Run'+str(runNum)
	suptitle(figTitle);
	print(fname1)
	show()
	figs.savefig(fname2+'.pdf', format = 'pdf')
	figs.savefig(fname1+'.png')

def plotMagRun2a(day,runNum,noiseNum = '00',bs = 2,fmax = 200,fignum = 1,ver = 'v16',leg=1,rScale=5,direc = 'Y', el = 1,relScale = 1,Ch = [1,2]):
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
	col = ['black','gold','blue','magenta','green']
	
	dataResp1 = [array([0,0]),array([0,0]),array([0,0])];
	dataResp2 =[array([0,0]),array([0,0]),array([0,0])];
	
	Chan = Ch[0];
	path = basepath +str(day) + '\\' + 'run_' + str(runNum);
	
	fPSD,fRESP,figTitle = getSingleRunSpecs(day,path,noiseNum,Chan,direc,runNum)
	
	if os.path.exists(fPSD)==True:
		data = loadtxt(fPSD);
		data = transpose(data);
		freq = binit(data[0],bs);
		ns1 = relScale*binit(data[1],bs);
		dataResp = loadtxt(fRESP);
		dataResp1 = transpose(dataResp)

			
		figs=figure(figsize = (11,8.5))
		clf();
		subplot(121)
		semilogy(freq,ns1,color = col[Ch[0]],label='Ch'+str(Chan));
		legend()
		grid(which='both');
		axis([0,fmax,1.0,10000])
		# xlabel('frequency, Hz');
		ylabel('PSD, fT/rHz');
		#~ legend();
		#~ legend();
		grid('on', which = 'both')
		
	elif os.path.exists(fPSD)==False:
		figs=figure(figsize = (11,8.5))
		clf();
		subplot(121)

	Chan = Ch[1];
	path = basepath +str(day) + '\\' + 'run_' + str(runNum);
	
	fPSD,fRESP,figTitle = getSingleRunSpecs(day,path,noiseNum,Chan,direc,runNum)

	if os.path.exists(fPSD):
		data = loadtxt(fPSD);
		data = transpose(data);
		freq = binit(data[0],bs);
		ns2 = relScale*binit(data[1],bs);
		dataResp = loadtxt(fRESP);
		dataResp2 = transpose(dataResp)

		
		subplot(121)
		semilogy(freq,ns2,color = col[Ch[1]],label='Ch '+str(Chan));
		legend();
		grid(which='both');
		axis([0,fmax,1.0,10000])
		grid('on', which = 'both')
	elif os.path.exists(fPSD)==False:
		subplot(121)





	xR = subplot(222);
	plot(dataResp1[0],1E5*dataResp1[1],color = col[Ch[0]], linewidth = 3, label = 'Ch'+str(Ch[0]));
	plot(dataResp2[0],1E5*dataResp2[1],color = col[Ch[1]], linewidth = 3, label = 'Ch'+str(Ch[1]));
	axis([0,fmax,0,rScale]);
	xR.yaxis.set_label_position("right")
	ylabel('Resp, 1E-5*V/fT');
	grid(which='both');
	legend()

	xR2 = subplot(224);
	plot(dataResp1[0],dataResp1[2] - dataResp1[2][0],color = col[Ch[0]], linewidth = 3, label = 'Ch'+str(Ch[0]));
	plot(dataResp2[0],dataResp2[2]-dataResp2[2][0],color = col[Ch[1]], linewidth = 3, label = 'Ch'+str(Ch[1]));
	#axis([0,fmax,-180,180]);
	xlim(xmax = fmax)
	xR2.yaxis.set_label_position("right")
	ylabel('Phase (deg)');
	grid(which='both');
	#~ legend()



	figTitle = str(day) + '\\' + 'run_' + str(runNum)+ ' B-'+str(direc)+':Channels' + str(Ch[0])+ ' and '+ str(Ch[1]);
	fname1 = dayDir+'\\Mag2'+str(direc)+str(day)+'_'+'Run'+str(runNum)
	fname2 = runDir+'\\Mag2'+str(direc)+str(day)+'_'+'Run'+str(runNum)
	suptitle(figTitle);
	print(fname1)
	show()
	figs.savefig(fname2+'.pdf', format = 'pdf')
	figs.savefig(fname1+'.png');

