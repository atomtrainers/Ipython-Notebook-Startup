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


def binit(Vec, bs):
    NVec = []
    nstep = np.int(np.size(Vec)/bs);

    for j in range (0,nstep-1):
        ss = np.sum(Vec[bs*j:bs*j+bs])
        NVec.append(ss/bs);
    return np.array(NVec)

def MagRunHelp():
	
	print("List of Frequently used functions")
	print("")
	f1 = str("1.  plotMagRun(day,runNum,noiseNum,noiseNum2,noiseNum3='00',Chan=1,bs = 1,fignum = 1, ver='v16',rScale=1,rgain = 1,fmax = 200, minpsd = 0.05, figLabel=['Total Noise','Probe Noise','Electronic Noise'])");

	f2 = str("2. (day,runNum,noiseNum,noiseNum2,noiseNum3='00',Chan=1,bs = 1,fignum = 1, ver='v16',rScale=1,rgain = 1,fmax = 200, psn = 0,cfit = 0,calib = 500E-9,Ipr = 100E-6,fitfunc=0)");

	f3 = str("3. plotMagRun4(day,runNum,noiseNum = '00',noiseNum2 = '01', noiseNum3 = '02',bs = 2,fmax = 200,fignum = 1,ver = 'v16',leg=1,rScale=5,direc = 'Y', el = 1)");

	f4 = str("4. plotMagRun3(day,runNum,noiseNum,Chan=1,bs = 1,fignum = 1, ver='v16')");

	f5 = str("5. plot4ProbeScan(fnum, fignum = 1)");

	f6 = str("6. plotMag1Chan(day,runNum,noiseNum,Chan=1,bs = 1,fignum = 1, ver='v16',rScale=1,rgain = 1,fmax = 200, psn = 0,cfit = 0,calib = 500E-9,Ipr = 100E-6,fitfunc=0)");

	f7 = str("7. plot4LaserScan(fnum, day = '', laser = 'Probe', fignum = 1)");

	f8 = str("8. plotMag1Chan(day,runNum,noiseNum,Chan=1,bs = 1,fignum = 1, ver='v16',rScale=1,rgain = 1,fmax = 200, psn = 0,cfit = 0,calib = 500E-9,Ipr = 100E-6,fitfunc=0)")
	 



	print(f1);
	print("")
	print(f2);
	print("")
	print(f3);
	print("")
	print(f4);
	print("")
	print(f5);
	print("")
	print(f6);
	print("")
	print(f7);
	print("")
	print(f8);

def getCompEnv(datasource='WIMR'):
	"""set the directories of the computer accordingly. Returns [basepath, and analysispath]"""
	cname = os.getenv('computername');

	if cname =='RAMSESE-II':
		b_path = 'C:\\Users\\sulai\\Documents\\LabVIEW Data\\DATA\\';
		a_path = 'C:\\analysis\\';
	elif cname =='THADLABS':
		if datasource=='WIMR':
			b_path = 'D:\\Users\\thadlabs\\Documents\\LabVIEW Data\\DATA\\';
			a_path = 'D:\\Analysis\\';
		elif datasource=='Chambo':
			b_path = 'D:\\Users\\thadlabs\\Documents\\Chambo LabVIEW Data\\DATA\\';
			a_path = 'D:\\Analysis\\';

	return [b_path,a_path];
    	
def getSingleRunSpecs(day,path,noiseNum,Chan,direc,runNum):
	if direc=='Y':
		f1 = path + '\\'+'noise_'+str(noiseNum)+'\\'+'Y_magnetic-field-PSD_'+str(Chan)+'.txt';
		f2 =  path+'\\'+'response'+'\\'+'Y_amp-n-phase_response_'+str(Chan)+'.txt';
		f3 = str(day) + '\\' + 'run_' + str(runNum)+'\\'+'Y - noise_'+str(noiseNum) + ': Channel ' + str(Chan);
	if direc=='X':
		f1 = path + '\\'+'noise_'+str(noiseNum)+'\\'+'X_magnetic-field-PSD_'+str(Chan)+'.txt';
		f2 =  path+'\\'+'response'+'\\'+'X_amp-n-phase_response_'+str(Chan)+'.txt';
		f3 = str(day) + '\\' + 'run_' + str(runNum)+'\\'+'X - noise_'+str(noiseNum) + ': Channel ' 
	return [f1,f2,f3]

def getRunSpecs(day,path,noiseNum,noiseNum2,noiseNum3,Chan,direc,runNum):
	if direc=='Y':
		f1 = path + '\\'+'noise_'+str(noiseNum)+'\\'+'Y_magnetic-field-PSD_'+str(Chan)+'.txt';
		f2 = path + '\\'+'noise_'+str(noiseNum2)+'\\'+'Y_magnetic-field-PSD_'+str(Chan)+'.txt';
		f3 = path + '\\'+'noise_'+str(noiseNum3)+'\\'+'Y_magnetic-field-PSD_'+str(Chan)+'.txt';
		f4 =  path+'\\'+'response'+'\\'+'Y_amp-n-phase_response_'+str(Chan)+'.txt';
		f5 = str(day) + '\\' + 'run_' + str(runNum)+'\\'+'Y - noise_'+str(noiseNum) + ': Channel ' + str(Chan);
	if direc=='X':
		f1 = path + '\\'+'noise_'+str(noiseNum)+'\\'+'X_magnetic-field-PSD_'+str(Chan)+'.txt';
		f2 = path + '\\'+'noise_'+str(noiseNum2)+'\\'+'X_magnetic-field-PSD_'+str(Chan)+'.txt';
		f3 = path + '\\'+'noise_'+str(noiseNum3)+'\\'+'X_magnetic-field-PSD_'+str(Chan)+'.txt';
		f4 =  path+'\\'+'response'+'\\'+'X_amp-n-phase_response_'+str(Chan)+'.txt';
		f5 = str(day) + '\\' + 'run_' + str(runNum)+'\\'+'X - noise_'+str(noiseNum) + ': Channel ' 
	return [f1,f2,f3,f4,f5]
	
def OneLorentz(x, a,b):
	return a*b**2/(x**2+b**2)
	
def TwoLorentz(x, a,b,c,d):
	return a*b**2/(x**2+b**2)+c*d**2/(x**2+d**2)

def OneSqrtLorentz(x, a,b):
	return a*b/(x**2+b**2)**.5

def plotMagRunPSN(loc='WIMR',day='2016.10.01',runNum='00',noiseNum='00',noiseNum2='01',noiseNum3='02',Chan=1,direc='Y',ver='v16',fignum = 1,bs = 1,nLim=[0.,200.,.05,10000.],psn=1,Ipr = 100E-6,calib = 2E-6,cfit = 0,fitfunc=1,dispFig=1,ZM=0,leg=1,lbl1=0,lbl2=0,lbl3=0):
	[basepath, analysispath] = getCompEnv(loc);
	basepath = basepath + ver + '\\';

	dayDir = analysispath+day
	

	a = os.path.isdir(dayDir);
	if a==0:
		os.mkdir(dayDir);

	runDir = analysispath+day+'\\Run'+runNum
	a = os.path.isdir(runDir);
	if a==0:
		os.mkdir(runDir);


	#    if int(noiseNum)<10:
	#       noiseNum = '0'+str(noiseNum);

	fname1 = runDir+'\\MagRun2PSN_'+str(day)+'_'+'R'+str(runNum)+'N'+str(noiseNum)+str(noiseNum2)+str(noiseNum3)+'_Chan'+str(Chan)+direc;
	fname2 = dayDir+'\\MagRun2PSN_'+str(day)+'_'+'R'+str(runNum)+'N'+str(noiseNum)+str(noiseNum2)+str(noiseNum3)+'_Chan'+str(Chan)+direc;


	path = basepath +str(day) + '\\' + 'run_' + str(runNum);
	fPSD = path + '\\'+'noise_'+str(noiseNum)+'\\'+'Y_magnetic-field-PSD_'+str(Chan)+'.txt';
	fPSD2 = path + '\\'+'noise_'+str(noiseNum2)+'\\'+'Y_magnetic-field-PSD_'+str(Chan)+'.txt';
	fPSD3 = path + '\\'+'noise_'+str(noiseNum3)+'\\'+'Y_magnetic-field-PSD_'+str(Chan)+'.txt';
	fRESP =  path+'\\'+'response'+'\\'+'Y_amp-n-phase_response_'+str(Chan)+'.txt';

	if direc=='X':
		fPSD = path + '\\'+'noise_'+str(noiseNum)+'\\'+'X_magnetic-field-PSD_'+str(Chan)+'.txt';
		fPSD2 = path + '\\'+'noise_'+str(noiseNum2)+'\\'+'X_magnetic-field-PSD_'+str(Chan)+'.txt';
		fPSD3 = path + '\\'+'noise_'+str(noiseNum3)+'\\'+'X_magnetic-field-PSD_'+str(Chan)+'.txt';
		fRESP =  path+'\\'+'response'+'\\'+'X_amp-n-phase_response_'+str(Chan)+'.txt';
		
	if psn==1 and lbl1==0: #calculates the "average" unbalance during the total noise run
		totalTS=fromfile(path+'\\'+'noise_'+str(noiseNum)+'\\'+'Y_noise_'+str(Chan)+'.bin',dtype='<d'); #imports total noise time series (in V)
		tmeanV=mean(abs(totalTS)) #calculates avg absolute value
		tanglerad=tmeanV*calib/(4*Ipr) #calculates angle in radians, I_diff/(2*I_sum)
		tangleurad=int(tanglerad*1E6) #converts to micro-radians and rounds

	if psn==1 and lbl2==0: #calculates the "average" unbalance during the probe noise run
		probeTS=fromfile(path+'\\'+'noise_'+str(noiseNum2)+'\\'+'Y_noise_'+str(Chan)+'.bin',dtype='<d'); #imports probe noise time series (in V)
		pmeanV=mean(abs(probeTS)) #calculates avg absolute value
		panglerad=pmeanV*calib/(4*Ipr) #calculates angle in radians, I_diff/(2*I_sum)
		pangleurad=int(panglerad*1E6) #converts to micro-radians and rounds

	if ZM==0:
		if psn==1:
			figTitle = str(day) + ' - ' + 'Run' + str(runNum)+' - '+'Noise'+str(noiseNum)+','+str(noiseNum2)+','+str(noiseNum3)+': Channel '+str(Chan)+direc+'  (' + str(round(1E6*Ipr)) + ' $\mu$A per PD)';
		if psn==0:
			figTitle = str(day) + ' - ' + 'Run' + str(runNum)+' - '+'Noise'+str(noiseNum)+','+str(noiseNum2)+','+str(noiseNum3)+': Channel '+str(Chan)+direc;
	if ZM==1:
		if psn==1:
			figTitle = str(day) + ' - ' + 'Run' + str(runNum)+' - '+'Noise'+str(noiseNum)+','+str(noiseNum2)+','+str(noiseNum3)+': Channel '+str(Chan)+direc+' [Z-mode] (' + str(round(1E6*Ipr)) + ' $\mu$A per PD)';
		if psn==0:
			figTitle = str(day) + ' - ' + 'Run' + str(runNum)+' - '+'Noise'+str(noiseNum)+','+str(noiseNum2)+','+str(noiseNum3)+': Channel '+str(Chan)+direc+' [Z-mode]';	
	
	data = loadtxt(fPSD);
	data = transpose(data);
	freq = binit(data[0],bs);
	PS = binit(data[1],bs);
	dataResp = loadtxt(fRESP);
	dataResp = transpose(dataResp)

	if noiseNum2 != noiseNum:
		data = loadtxt(fPSD2);
		data = transpose(data);
		freq2 = binit(data[0],bs);
		PS2 = binit(data[1],bs);


	if noiseNum3 != noiseNum:
		data = loadtxt(fPSD3);
		data = transpose(data);
		freq3 = binit(data[0],bs);
		PS3 = binit(data[1],bs);

	shotnoiseB=0;

	if Chan==1:
		respcolor='y-';
	elif Chan==2:
		respcolor='b-';
	elif Chan==3:
		respcolor='m-';
	elif Chan==4:
		respcolor='g-';
	else:
		respcolor='k-';
	
	fig = figure(fignum,figsize=(16.0,6.0))

	#Assigns standard labels of 'total noise', 'probe noise', and 'electronic noise' unless user specifies otherwise
	if lbl1==0:
		if psn==1:
			label1='Total Noise, $\phi$ = '+str(tangleurad)+' $\mu$rad'
		else:
			label1='Total Noise'
	else:
		label1=lbl1	
	if lbl2==0:
		if psn==1:
			label2='Probe Noise, $\phi$ = '+str(pangleurad)+' $\mu$rad'
		else:
			label2='Probe Noise'
	else:
		label2=lbl2
	if lbl3==0:
		label3='Electronic Noise'
	else:
		label3=lbl3
		
	subplot(121)
	semilogy(freq,PS,'k-',label=label1);
	if noiseNum2!=noiseNum:
		semilogy(freq2,PS2,'r-',label=label2);
	if noiseNum3!=noiseNum and noiseNum3!=noiseNum2:
		semilogy(freq3,PS3,'b-',label=label3);
		
	grid(which='both');
	axis(nLim)
	xlabel('Frequency [Hz]',fontsize=14); xticks(fontsize=14)
	ylabel('PSD [fT/rHz]',fontsize=14); yticks(fontsize=14)
	if leg==1:
		legend();

	ax1 = subplot(122);
	plot(dataResp[0],1E6*dataResp[1],respcolor,linewidth=2, label = 'meas.');
	xlim(xmax = nLim[1]);
	xlabel('Frequency [Hz]',fontsize=14); xticks(fontsize=14)
	ylabel('Mag. Resp. [$\mu$V/fT]',fontsize=14); yticks(fontsize=14)
	grid('on');
	suptitle(figTitle,fontsize=16)

	phSc = 1E-6*calib/(4*Ipr);

	fig.subplots_adjust(wspace = 0.32)

	pfit = [0,0,0,0];
	if psn ==1:
		
		shotnoiseA = sqrt(4*Ipr*1.6E-19);   #measured in A/rHz
		shotnoiseV = shotnoiseA/calib    #measured in V/rHz
		shotnoiseB = shotnoiseV/dataResp[1];  #measured in fT/rHz
		
		if cfit == 1:
			dx,dy = dataResp[0],1E6*dataResp[1];
		
			textcoordx=nLim[1]/2.5;
			textcoordy=max(dy)/1.5;
			
			if fitfunc==1:
				pfit,perr = curve_fit(OneLorentz,dx,dy,p0 = [max(dy),30.]);
				fa,fb = pfit;
				width=str(round(fb));
				amp=str(round(fa));
				yfit = OneLorentz(dx,fa,fb);
				shotnoiseB = 1E6*shotnoiseV/yfit;  #measured in fT/rHz
				
			if fitfunc==2:
				pfit,perr = curve_fit(TwoLorentz,dx,dy, p0 = [max(dy),30.,max(dy)/3,30.]);
				fa,fb,fc,fd = pfit;
				width=str(round(fb))+', '+str(round(fd));
				amp=str(round(fa))+', '+str(round(fc));
				yfit = TwoLorentz(dx,fa,fb,fc,fd)
				shotnoiseB = 1E6*shotnoiseV/yfit;  #measured in fT/rHz
				
			if fitfunc==3:
				pfit,perr = curve_fit(OneSqrtLorentz,dx,dy,p0 = [max(dy),30.]);
				fa,fb = pfit;
				width=str(round(fb));
				amp=str(round(fa));
				yfit = OneSqrtLorentz(dx,fa,fb);
				shotnoiseB = 1E6*shotnoiseV/yfit;  #measured in fT/rHz
			
			figure(fignum);
			subplot(122)
			plot(dx,yfit,'r-', label = 'fit');
			legend()
			#text(textcoordx,textcoordy,'Linewidth = '+width+' Hz',fontsize=14);
			#text(textcoordx,textcoordy*.9,'Amplitudes(s) = '+amp+' $\mu$V/fT',fontsize=14);
			
		a = ax1.axis();

		ax2 = ax1.twinx()
		ax2.axis([a[0],a[1],a[2],1E7*phSc*a[3]])
		ax2.set_ylabel('Mag. Resp. (1E-7*rad/fT)',fontsize=14); yticks(fontsize=14)

			 
		
		
		subplot(121);
		semilogy(dataResp[0],shotnoiseB,linestyle='-',color=(0,1,0), linewidth = 3,label='photon shot noise');
		legend();
		
		
	if dispFig!=1:
		close(fignum);
	show();

	fig.savefig(fname1+'.pdf', format = 'pdf');
	fig.savefig(fname2+'.png');
	return [freq,PS,dataResp[0],1E5*dataResp[1]];

def plotMagRunZM(loc='WIMR',day='2016.10.01',runNum='00',noiseNum1='00',noiseNum2='01',noiseNum3='02',Chan=1,bs = 1,fignum=1,ver='v16',fmax=200,psn=1,calib=2E-6,Ipr=100E-6,leg=1,lbl1=0,lbl2=0,lbl3=0):
	'''Plots the X and Y noise from a single channel Z-mode run, as well as the X and Y response'''
	[basepath, analysispath] = getCompEnv(loc);
	basepath = basepath + ver + '\\';

	dayDir = analysispath+day
	
	a = os.path.isdir(dayDir);
	if a==0:
		os.mkdir(dayDir);

	runDir = analysispath+day+'\\Run'+runNum
	a = os.path.isdir(runDir);
	if a==0:
		os.mkdir(runDir);

	fname1 = runDir+'\\MagRun2PSN_'+str(day)+'_'+'R'+str(runNum)+'N'+str(noiseNum1)+str(noiseNum2)+str(noiseNum3)+'_Chan'+str(Chan)+'XY';
	fname2 = dayDir+'\\MagRun2PSN_'+str(day)+'_'+'R'+str(runNum)+'N'+str(noiseNum1)+str(noiseNum2)+str(noiseNum3)+'_Chan'+str(Chan)+'XY';


	path = basepath +str(day) + '\\' + 'run_' + str(runNum);
	frespY = path + '\\'+'response'+'\\'+'Y_amp-n-phase_response_'+str(Chan)+'.txt';
	fpsd1Y = path + '\\'+'noise_'+str(noiseNum1)+'\\'+'Y_magnetic-field-PSD_'+str(Chan)+'.txt';
	fpsd2Y = path + '\\'+'noise_'+str(noiseNum2)+'\\'+'Y_magnetic-field-PSD_'+str(Chan)+'.txt';
	fpsd3Y = path + '\\'+'noise_'+str(noiseNum3)+'\\'+'Y_magnetic-field-PSD_'+str(Chan)+'.txt';

	frespX = path + '\\'+'response'+'\\'+'X_amp-n-phase_response_'+str(Chan)+'.txt';
	fpsd1X = path + '\\'+'noise_'+str(noiseNum1)+'\\'+'X_magnetic-field-PSD_'+str(Chan)+'.txt';
	fpsd2X = path + '\\'+'noise_'+str(noiseNum2)+'\\'+'X_magnetic-field-PSD_'+str(Chan)+'.txt';
	fpsd3X = path + '\\'+'noise_'+str(noiseNum3)+'\\'+'X_magnetic-field-PSD_'+str(Chan)+'.txt';

	dataRespY = transpose(loadtxt(frespY));
	dataRespX = transpose(loadtxt(frespX));
	
	data = transpose(loadtxt(fpsd1Y));
	freq = binit(data[0],bs);
	PS1Y = binit(data[1],bs);
	
	data = transpose(loadtxt(fpsd1X));
	PS1X = binit(data[1],bs);

	if noiseNum2 != noiseNum1:
		data = transpose(loadtxt(fpsd2Y));
		PS2Y = binit(data[1],bs);
		data = transpose(loadtxt(fpsd2X));
		PS2X = binit(data[1],bs);

	if noiseNum3 != noiseNum1 and noiseNum3 != noiseNum2:
		data = transpose(loadtxt(fpsd3Y));
		PS3Y = binit(data[1],bs);
		data = transpose(loadtxt(fpsd3X));
		PS3X = binit(data[1],bs);
		
	if psn==1 and lbl1==0: #calculates the "average" unbalance during the total noise run
		totalTS=fromfile(path+'\\'+'noise_'+str(noiseNum1)+'\\'+'Y_noise_'+str(Chan)+'.bin',dtype='<d'); #imports Y total noise time series (in V)
		tmeanV=mean(abs(totalTS)) #calculates avg absolute value
		tanglerad=tmeanV*calib/(4*Ipr) #calculates angle in radians, I_diff/(2*I_sum)
		tangleurad=int(tanglerad*1E6) #converts to micro-radians and rounds
    
	if psn==1 and lbl2==0: #calculates the "average" unbalance during the probe noise run
		probeTS=fromfile(path+'\\'+'noise_'+str(noiseNum2)+'\\'+'Y_noise_'+str(Chan)+'.bin',dtype='<d'); #imports Y probe noise time series (in V)
		pmeanV=mean(abs(probeTS)) #calculates avg absolute value
		panglerad=pmeanV*calib/(4*Ipr) #calculates angle in radians, I_diff/(2*I_sum)
		pangleurad=int(panglerad*1E6) #converts to micro-radians and rounds

		
	respcolor=['k-','y-','b-','m-','g-'];
	respcolor2=['k--','y--','b--','m--','g--'];
		
	#Assigns standard labels of 'total noise', 'probe noise', and 'electronic noise' unless user specifies otherwise
	if lbl1==0:
		label1='Total Noise, $\phi$ = '+str(tangleurad)+' $\mu$rad'
	else:
		label1=lbl1	
	if lbl2==0:
		label2='Probe Noise, $\phi$ = '+str(pangleurad)+' $\mu$rad'
	else:
		label2=lbl2
	if lbl3==0:
		label3='Electronic Noise'
	else:
		label3=lbl3
	
	fig = figure(fignum,figsize=(16.0,6.0))
	
	noisey=subplot(131);
	semilogy(freq,PS1Y,'k-',label=label1);
	if noiseNum2!=noiseNum1:
		semilogy(freq,PS2Y,'r-',label=label2);
	if noiseNum3!=noiseNum1 and noiseNum3!=noiseNum2:
		semilogy(freq,PS3Y,'b-',label=label3);
	grid('on',which='both');
	axis([0,fmax,0.05,10000]);
	xlabel('Frequency [Hz]',fontsize=14); xticks(fontsize=14);
	ylabel('PSD [fT/rHz]',fontsize=14); yticks(fontsize=14);
	title('Y Noise',fontsize=16)
	if leg==1:
		legend();
	
	noisex=subplot(132);
	semilogy(freq,PS1X,'k-',label=label1);
	if noiseNum2!=noiseNum1:
		semilogy(freq,PS2X,'r-',label=label2);
	if noiseNum3!=noiseNum1 and noiseNum3!=noiseNum2:
		semilogy(freq,PS3X,'b-',label=label3);
	grid('on',which='both');
	axis([0,fmax,0.05,10000]);
	xlabel('Frequency [Hz]',fontsize=14); xticks(fontsize=14);
	yticks(fontsize=14);
	title('X Noise',fontsize=16)
	
		
	resp=subplot(133);
	plot(dataRespY[0],1E6*dataRespY[1],respcolor[Chan],linewidth=2,label='Y response');
	plot(dataRespX[0],1E6*dataRespX[1],respcolor2[Chan],linewidth=2,label='X response');
	xlim(xmax = fmax);
	xlabel('Frequency [Hz]',fontsize=14); xticks(fontsize=14);
	ylabel('Mag. Resp. [$\mu$V/fT]',fontsize=14); yticks(fontsize=14);
	if psn==0:
		resp.yaxis.set_label_position('right');
	grid('on');
	title('Response',fontsize=16);
	legend();
	
	if psn==1:
		shotnoiseA = sqrt(4*Ipr*1.6E-19);   #measured in A/rHz
		shotnoiseV = shotnoiseA/calib    #measured in V/rHz
		shotnoiseBy = shotnoiseV/dataRespY[1];  #measured in fT/rHz
		shotnoiseBx = shotnoiseV/dataRespX[1];  #measured in ft/rHz
		
		subplot(131);
		semilogy(dataRespY[0],shotnoiseBy,linestyle='-',color=(0,1,0), linewidth = 3,label='photon shot noise');
		subplot(132);
		semilogy(dataRespX[0],shotnoiseBx,linestyle='-',color=(0,1,0), linewidth = 3,label='photon shot noise');
		
		subplot(133);
		a = resp.axis();

		phSc = 1E-6*calib/(4*Ipr);
		ax2 = resp.twinx();
		ax2.axis([a[0],a[1],a[2],1E7*phSc*a[3]]);
		ax2.set_ylabel('Mag. Resp. (1E-7*rad/fT)',fontsize=14); yticks(fontsize=14);

	if psn==0:
		figtitle=day+' - Run '+runNum+' - Noise '+str(noiseNum1)+','+str(noiseNum2)+','+str(noiseNum3)+': Channel '+str(Chan)+'XY';
	elif psn==1:
		figtitle=day+' - Run '+runNum+' - Noise '+str(noiseNum1)+','+str(noiseNum2)+','+str(noiseNum3)+': Channel '+str(Chan)+'XY'+'  (' + str(round(1E6*Ipr)) + ' $\mu$A per PD)';

	suptitle(figtitle,fontsize=16);
	
	fig.savefig(fname1+'.pdf', format = 'pdf');
	fig.savefig(fname2+'.png');

def plotMagRun5(loc='WIMR',day='2016.10.01',runNum='00',noiseNum1 = '00',noiseNum2 = '01', noiseNum3 = '02',Chan=1,bs = 2,fmax = 100,fignum = 1,ver = 'v16',leg=1,rScale=50,direc = 'Y'):
	'''Plots noise, amplitude and phase response for a single run, single channel'''
	[basepath, analysispath] = getCompEnv(loc);
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

	path = basepath +str(day) + '\\' + 'run_' + str(runNum);
	
	fPSD,fPSD2,fPSD3,fRESP,figTitle = getRunSpecs(day,path,noiseNum1,noiseNum2,noiseNum3,Chan,direc,runNum)
	
	if ver == 'v12':
		fPSD = path + '\\'+'noise_'+str(noiseNum1)+'\\'+'magnetic-field-PSD_'+str(Chan)+'.txt';
		fPSD2 = path + '\\'+'noise_'+str(noiseNum2)+'\\'+'magnetic-field-PSD_'+str(Chan)+'.txt';
		fPSD3 = path + '\\'+'noise_'+str(noiseNum3)+'\\'+'magnetic-field-PSD_'+str(Chan)+'.txt';
		fRESP =  path+'\\'+'response'+'\\'+'amp-n-phase_response_'+str(Chan)+'.txt';

	data = loadtxt(fPSD);
	data = transpose(data);
	freq = binit(data[0],bs);
	PS = binit(data[1],bs);

	data = loadtxt(fPSD2);
	data = transpose(data);
	freq = binit(data[0],bs);
	PS2 = binit(data[1],bs);

	data = loadtxt(fPSD3);
	data = transpose(data);
	freq = binit(data[0],bs);
	PS3 = binit(data[1],bs);
	
	dataResp = loadtxt(fRESP);
	dataResp1 = transpose(dataResp)
	
	color=['k-','y-','b-','m-','g-']
	
	figs=figure(figsize = (16,8))
	clf();
	subplot(121)
	semilogy(freq,PS,'k-',label='Total Noise');
	if noiseNum2 != noiseNum1:
		semilogy(freq,PS2,'r-',label='Probe Noise');
	if noiseNum3 != noiseNum1:
		semilogy(freq,PS3,'b-',label='Electronic Noise');
	grid(which='both');
	axis([0,fmax,0.1,10000])
	xlabel('Frequency [Hz]');
	ylabel('PSD [fT/rHz]');
	legend();

	xR = subplot(222);
	plot(dataResp1[0],1E6*dataResp1[1],color[Chan], linewidth = 3, label = 'Chan '+str(Chan));
	xlim(0,fmax)
	#axis([0,fmax,0,rScale]);
	xR.yaxis.set_label_position("right")
	ylabel('Resp, $\mu$V/fT');
	grid(which='both');
	legend()

	xR2 = subplot(224);
	plot(dataResp1[0],dataResp1[2]-dataResp1[2][0],color[Chan], linewidth = 3);
	axis([0,fmax,-180,180]);
	xR2.yaxis.set_label_position("right")
	ylabel('Phase (deg)');
	xlabel('Frequency [Hz]')
	grid(which='both');
	#~ legend()

	figTitle = str(day) + '\\' + 'run_' + str(runNum)+ ' B-'+str(direc)+':Channel '+str(Chan);
	fname1 = dayDir+'\\MagRun5_'+str(day)+'_R'+str(runNum)+'_N'+str(noiseNum1)+str(noiseNum2)+str(noiseNum3)+'_Chan'+str(Chan)+'_'+str(direc);
	fname2 = runDir+'\\MagRun5_'+str(day)+'_R'+str(runNum)+'_N'+str(noiseNum1)+str(noiseNum2)+str(noiseNum3)+'_Chan'+str(Chan)+'_'+str(direc);
	suptitle(figTitle);
	show();
	figs.savefig(fname2+'.pdf', format = 'pdf');
	figs.savefig(fname1+'.png');

def plotMagRun4PSN(loc='WIMR',day='2016.10.01',runNum='00',noiseNum='00',noiseNum2='01',noiseNum3='02',bs=2,fmax=200,fignum=1,ver='v16',leg=1,direc = 'Y',psn=0,Ipr=[1E-4,1E-4,1E-4,1E-4],calib=[2E-6,1E-6,2E-6,2E-6],ZM=0,lbl=[0,0,0]):
	[basepath, analysispath] = getCompEnv(loc);
	basepath = basepath + ver + '\\';

	dayDir = analysispath+day
	a = os.path.isdir(dayDir);
	if a==0:
		os.mkdir(dayDir);

	runDir = analysispath+day+'\\Run'+runNum
	a = os.path.isdir(runDir);
	if a==0:
		os.mkdir(runDir);

	fs = 18;
	labelfs=14
	
	path = basepath +str(day) + '\\' + 'run_' + str(runNum);
	
	if lbl[0]==0:
		lbl0='Total Noise'
	else:
		lbl0=lbl[0]
		
	if lbl[1]==0:
		lbl1='Probe Noise'
	else:
		lbl1=lbl[1]
		
	if lbl[2]==0:
		lbl2='Electronic Noise'
	else:
		lbl2=lbl[2]
	
	Chan = 1;
	#File paths for noise measurements and response
	fPSD,fPSD2,fPSD3,fRESP,figTitle = getRunSpecs(day,path,noiseNum,noiseNum2,noiseNum3,Chan,direc,runNum)
	
	#Response
	dataResp1 = transpose(loadtxt(fRESP));

	#Total magnetic noise
	data = transpose(loadtxt(fPSD));
	freq = binit(data[0],bs);
	PSa = binit(data[1],bs);
	
	#Probe (technical) noise
	data = transpose(loadtxt(fPSD2));
	freq = binit(data[0],bs);
	PSb = binit(data[1],bs);

	#Electronic noise
	data = transpose(loadtxt(fPSD3));
	freq = binit(data[0],bs);
	PSc = binit(data[1],bs);
	
	figs=figure(figsize = (16,10));
	clf();
	
	#Plots power spectra for Channel 1
	n1=subplot(231);
	semilogy(freq,PSa,'k-',label=lbl0);
	if noiseNum2 != noiseNum:
		semilogy(freq,PSb,'r-',label=lbl1);
	if noiseNum3 != noiseNum and noiseNum3 != noiseNum2:
		semilogy(freq,PSc,'b-',label=lbl2);
	if psn==1:
		shotnoiseB=sqrt(4*Ipr[0]*1.6E-19)/calib[0]/dataResp1[1];
		semilogy(dataResp1[0],shotnoiseB,linestyle='-',color=(0,1,0),linewidth=3,label='Photon Shot Noise');
	grid(which='both');
	axis([0,fmax,0.1,10000]);
	xticks(fontsize=labelfs)
	ylabel('PSD [fT/rHz]',fontsize=labelfs); yticks(fontsize=labelfs)
	text(2,4000,'Ch1',fontsize=fs);

	Chan = 2;
	#File paths for noise measurements and response
	fPSD,fPSD2,fPSD3,fRESP,figTitle = getRunSpecs(day,path,noiseNum,noiseNum2,noiseNum3,Chan,direc,runNum)

	#Response
	dataResp2 = transpose(loadtxt(fRESP));
	
	#Total magnetic noise
	data = transpose(loadtxt(fPSD));
	freq = binit(data[0],bs);
	PSa = binit(data[1],bs);

	#Probe (technical) noise
	data = transpose(loadtxt(fPSD2));
	freq = binit(data[0],bs);
	PSb = binit(data[1],bs);

	#Electronic noise
	data = transpose(loadtxt(fPSD3));
	freq = binit(data[0],bs);
	PSc = binit(data[1],bs);
	
	#Plots power spectra for Channel 2
	n2=subplot(232);
	semilogy(freq,PSa,'k-',label=lbl0);
	if noiseNum2 != noiseNum:
		semilogy(freq,PSb,'r-',label=lbl1);
	if noiseNum3 != noiseNum and noiseNum3 != noiseNum2:
		semilogy(freq,PSc,'b-',label=lbl2);
	if psn==1:
		shotnoiseB=sqrt(4*Ipr[1]*1.6E-19)/calib[1]/dataResp2[1];
		semilogy(dataResp2[0],shotnoiseB,linestyle='-',color=(0,1,0),linewidth=3,label='Photon Shot Noise');
	grid(which='both');
	axis([0,fmax,0.1,10000]);
	xticks(fontsize=labelfs);
	yticks(fontsize=labelfs);
	text(2,4000,'Ch2',fontsize=fs);
	if leg==1:
		legend();
	
	Chan = 3;
	#File paths for noise measurements and response
	fPSD,fPSD2,fPSD3,fRESP,figTitle = getRunSpecs(day,path,noiseNum,noiseNum2,noiseNum3,Chan,direc,runNum)
	
	#Response
	dataResp3 = transpose(loadtxt(fRESP));
	
	#Total magnetic noise
	data = transpose(loadtxt(fPSD));
	freq = binit(data[0],bs);
	PSa = binit(data[1],bs);

	#Probe (technical) noise
	data = transpose(loadtxt(fPSD2));
	freq = binit(data[0],bs);
	PSb = binit(data[1],bs);

	#Electronic noise
	data = transpose(loadtxt(fPSD3));
	freq = binit(data[0],bs);
	PSc = binit(data[1],bs);
	
	#Plots power spectra for Channel 3
	n3=subplot(234);
	semilogy(freq,PSa,'k-',label=lbl0);
	if noiseNum2 != noiseNum:
		semilogy(freq,PSb,'r-',label=lbl1);
	if noiseNum3 != noiseNum and noiseNum3 != noiseNum2:
		semilogy(freq,PSc,'b-',label=lbl2);
	if psn==1:
		shotnoiseB=sqrt(4*Ipr[2]*1.6E-19)/calib[2]/dataResp3[1];
		semilogy(dataResp3[0],shotnoiseB,linestyle='-',color=(0,1,0),linewidth=3,label='Photon Shot Noise');
	grid(which='both');
	axis([0,fmax,0.1,10000]);
	xlabel('Frequency [Hz]',fontsize=labelfs); xticks(fontsize=labelfs);
	ylabel('PSD [fT/rHz]',fontsize=labelfs); yticks(fontsize=labelfs);
	text(2,4000,'Ch3',fontsize=fs);
	
	Chan = 4;
	#File paths for noise measurements and response
	fPSD,fPSD2,fPSD3,fRESP,figTitle = getRunSpecs(day,path,noiseNum,noiseNum2,noiseNum3,Chan,direc,runNum)
	
	#Response
	dataResp4 = transpose(loadtxt(fRESP));

	#Total magnetic noise
	data = transpose(loadtxt(fPSD));
	freq = binit(data[0],bs);
	PSa = binit(data[1],bs);

	#Probe (technical) noise
	data = transpose(loadtxt(fPSD2));
	freq = binit(data[0],bs);
	PSb = binit(data[1],bs);

	#Electronic noise
	data = transpose(loadtxt(fPSD3));
	freq = binit(data[0],bs);
	PSc = binit(data[1],bs);
	
	#Plots power spectra for Channel 4
	n3=subplot(235);
	semilogy(freq,PSa,'k-',label=lbl0);
	if noiseNum2 != noiseNum:
		semilogy(freq,PSb,'r-',label=lbl1);
	if noiseNum3 != noiseNum and noiseNum3 != noiseNum2:
		semilogy(freq,PSc,'b-',label=lbl2);
	if psn==1:
		shotnoiseB=sqrt(4*Ipr[3]*1.6E-19)/calib[3]/dataResp4[1];
		semilogy(dataResp4[0],shotnoiseB,linestyle='-',color=(0,1,0),linewidth=3,label='Photon Shot Noise');
	grid(which='both');
	axis([0,fmax,0.1,10000]);
	xlabel('Frequency [Hz]',fontsize=labelfs); xticks(fontsize=labelfs);
	yticks(fontsize=labelfs);
	text(2,4000,'Ch4',fontsize=fs);
	
	ampresp=subplot(233);
	if psn==0:
		plot(dataResp1[0],1E6*dataResp1[1],'y-', linewidth = 3, label = 'Ch1'); #response in V/fT
		plot(dataResp2[0],1E6*dataResp2[1],'b-', linewidth = 3, label = 'Ch2');
		plot(dataResp3[0],1E6*dataResp3[1],'m-', linewidth = 3, label = 'Ch3');
		plot(dataResp4[0],1E6*dataResp4[1],'g-', linewidth = 3, label = 'Ch4');
		ampresp.yaxis.set_label_position("right");
		xticks(fontsize=labelfs);
		ylabel('Response [$\mu$V/fT]',fontsize=labelfs); yticks(fontsize=labelfs);
	elif psn==1:
		Aresp1=dataResp1[1]*calib[0]/(4*Ipr[0]); #response in rad/fT
		plot(dataResp1[0],1E7*Aresp1,'y-', linewidth = 3, label = 'Ch1');
		Aresp2=dataResp2[1]*calib[1]/(4*Ipr[1]);
		plot(dataResp2[0],1E7*Aresp2,'b-', linewidth = 3, label = 'Ch2');
		Aresp3=dataResp3[1]*calib[2]/(4*Ipr[2]);
		plot(dataResp3[0],1E7*Aresp3,'m-', linewidth = 3, label = 'Ch3');
		Aresp4=dataResp4[1]*calib[3]/(4*Ipr[3]);
		plot(dataResp4[0],1E7*Aresp4,'g-', linewidth = 3, label = 'Ch4');
		ampresp.yaxis.set_label_position("right");
		xticks(fontsize=labelfs);
		ylabel('Response [$10^{-7}$ rad/fT]',fontsize=labelfs); yticks(fontsize=labelfs);
	xlim(0,fmax);
	ylim(ymin=0);
	grid('on');
	legend(loc=0);
	

	phaseresp = subplot(236);
	plot(dataResp1[0],dataResp1[2]-dataResp1[2][0],'y-', linewidth = 3, label = 'Ch1');
	plot(dataResp2[0],dataResp2[2]-dataResp2[2][0],'b-', linewidth = 3, label = 'Ch2');
	plot(dataResp3[0],dataResp3[2]-dataResp3[2][0],'m-', linewidth = 3, label = 'Ch3');
	plot(dataResp4[0],dataResp4[2]-dataResp4[2][0],'g-', linewidth = 3, label = 'Ch4');
	axis([0,fmax,-180,180]);
	phaseresp.yaxis.set_label_position("right");
	xlabel('Frequency [Hz]',fontsize=labelfs); xticks(fontsize=labelfs);
	ylabel('Phase [deg]',fontsize=labelfs); yticks(fontsize=labelfs);
	grid(which='both');

	#show();
	
	direclower=direc.lower()
	
	if ZM==0:
		figTitle = str(day) + ' - Run ' + str(runNum)+ ' - Noise '+str(noiseNum)+','+str(noiseNum2)+','+str(noiseNum3)+' - B'+direclower+' [DC-SERF]: Ch. 1 - 4';
	if ZM==1 and direc=='X':	
		figTitle = str(day) + ' - Run ' + str(runNum)+ ' - Noise '+str(noiseNum)+','+str(noiseNum2)+','+str(noiseNum3)+' - B'+direclower+' [Z-mode, 1f demodulated]: Ch. 1 - 4';
	if ZM==1 and direc=='Y':	
		figTitle = str(day) + ' - Run ' + str(runNum)+ ' - Noise '+str(noiseNum)+','+str(noiseNum2)+','+str(noiseNum3)+' - B'+direclower+' [Z-mode, DC lowpassed]: Ch. 1 - 4';
		
	fname1 = dayDir+'\\Mag4PSN_'+str(day)+'_'+'R'+str(runNum)+'N'+str(noiseNum)+str(noiseNum2)+str(noiseNum3)+'_'+str(direc);
	fname2 = runDir+'\\Mag4PSN_'+str(day)+'_'+'R'+str(runNum)+'N'+str(noiseNum)+str(noiseNum2)+str(noiseNum3)+'_'+str(direc);
	suptitle(figTitle,fontsize=18);
	figs.savefig(fname2+'.pdf', format = 'pdf');
	figs.savefig(fname1+'.png');
	
def plotMagRun2PSN(loc='WIMR',day='2016.10.01',runNum='00',noiseNum='00',noiseNum2='01',noiseNum3='02',Channels=[1,2],bs=4,fmax=200,fignum=1,ver='v16',leg=1,direc = 'Y',psn=0,Ipr=[1E-4,1E-4],calib=[1E-6,1E-6]):
	[basepath, analysispath] = getCompEnv(loc);
	basepath = basepath + ver + '\\';

	dayDir = analysispath+day
	a = os.path.isdir(dayDir);
	if a==0:
		os.mkdir(dayDir);

	runDir = analysispath+day+'\\Run'+runNum
	a = os.path.isdir(runDir);
	if a==0:
		os.mkdir(runDir);

	fs = 18;
	labelfs=14
	
	path = basepath +str(day) + '\\' + 'run_' + str(runNum);
	
	Chan = Channels[0];
	#File paths for noise measurements and response
	fPSD,fPSD2,fPSD3,fRESP,figTitle = getRunSpecs(day,path,noiseNum,noiseNum2,noiseNum3,Chan,direc,runNum)
	
	#Response
	dataResp1 = transpose(loadtxt(fRESP));

	#Total magnetic noise
	data = transpose(loadtxt(fPSD));
	freq = binit(data[0],bs);
	PSa = binit(data[1],bs);
	
	#Probe (technical) noise
	data = transpose(loadtxt(fPSD2));
	freq = binit(data[0],bs);
	PSb = binit(data[1],bs);

	#Electronic noise
	data = transpose(loadtxt(fPSD3));
	freq = binit(data[0],bs);
	PSc = binit(data[1],bs);
	
	figs=figure(figsize = (16,6));
	clf();
	
	#Plots power spectra for first channel
	n1=subplot(131);
	semilogy(freq,PSa,'k-',label='Total Noise');
	if noiseNum2 != noiseNum:
		semilogy(freq,PSb,'r-',label='Probe Noise');
	if noiseNum3 != noiseNum and noiseNum3 != noiseNum2:
		semilogy(freq,PSc,'b-',label='Electronic Noise');
	if psn==1:
		shotnoiseB=sqrt(4*Ipr[0]*1.6E-19)/calib[0]/dataResp1[1];
		semilogy(dataResp1[0],shotnoiseB,linestyle='-',color=(0,1,0),linewidth=3,label='Photon Shot Noise');
	grid(which='both');
	axis([0,fmax,0.1,10000]);
	xticks(fontsize=labelfs)
	ylabel('PSD [fT/rHz]',fontsize=labelfs); yticks(fontsize=labelfs)
	text(2,4000,'Ch1',fontsize=fs);

	Chan = Channels[1];
	#File paths for noise measurements and response
	fPSD,fPSD2,fPSD3,fRESP,figTitle = getRunSpecs(day,path,noiseNum,noiseNum2,noiseNum3,Chan,direc,runNum)

	#Response
	dataResp2 = transpose(loadtxt(fRESP));
	
	#Total magnetic noise
	data = transpose(loadtxt(fPSD));
	freq = binit(data[0],bs);
	PSa = binit(data[1],bs);

	#Probe (technical) noise
	data = transpose(loadtxt(fPSD2));
	freq = binit(data[0],bs);
	PSb = binit(data[1],bs);

	#Electronic noise
	data = transpose(loadtxt(fPSD3));
	freq = binit(data[0],bs);
	PSc = binit(data[1],bs);
	
	#Plots power spectra for second channel
	n2=subplot(132);
	semilogy(freq,PSa,'k-',label='Total Noise');
	if noiseNum2 != noiseNum:
		semilogy(freq,PSb,'r-',label='Probe Noise');
	if noiseNum3 != noiseNum and noiseNum3 != noiseNum2:
		semilogy(freq,PSc,'b-',label='Electronic Noise');
	if psn==1:
		shotnoiseB=sqrt(4*Ipr[1]*1.6E-19)/calib[1]/dataResp2[1];
		semilogy(dataResp2[0],shotnoiseB,linestyle='-',color=(0,1,0),linewidth=3,label='Photon Shot Noise');
	grid(which='both');
	axis([0,fmax,0.1,10000]);
	xticks(fontsize=labelfs);
	yticks(fontsize=labelfs);
	text(2,4000,'Ch2',fontsize=fs);
	if leg==1:
		legend();
	
	respcolor=['k-','y-','b-','m-','g-'];
	
	ampresp=subplot(233);
	if psn==0:
		plot(dataResp1[0],1E6*dataResp1[1],respcolor[Channels[0]], linewidth = 3, label = 'Ch'+str(Channels[0])); #response in V/fT
		plot(dataResp2[0],1E6*dataResp2[1],respcolor[Channels[1]], linewidth = 3, label = 'Ch'+str(Channels[1]));
		ampresp.yaxis.set_label_position("right");
		xticks(fontsize=labelfs);
		ylabel('Response [$\mu$V/fT]',fontsize=labelfs); yticks(fontsize=labelfs);
	elif psn==1:
		Aresp1=dataResp1[1]*calib[0]/(4*Ipr[0]); #response in rad/fT
		plot(dataResp1[0],1E7*Aresp1,respcolor[Channels[0]], linewidth = 3, label = 'Ch'+str(Channels[0]));
		Aresp2=dataResp2[1]*calib[1]/(4*Ipr[1]);
		plot(dataResp2[0],1E7*Aresp2,respcolor[Channels[1]], linewidth = 3, label = 'Ch'+str(Channels[1]));
		ampresp.yaxis.set_label_position("right");
		xticks(fontsize=labelfs);
		ylabel('Response [$10^{-7}$ rad/fT]',fontsize=labelfs); yticks(fontsize=labelfs);
	xlim(0,fmax);
	ylim(ymin=0);
	grid('on');
	legend(loc=0);
	

	phaseresp = subplot(236);
	plot(dataResp1[0],dataResp1[2]-dataResp1[2][0],respcolor[Channels[0]], linewidth = 3, label = 'Ch'+str(Channels[0]));
	plot(dataResp2[0],dataResp2[2]-dataResp2[2][0],respcolor[Channels[1]], linewidth = 3, label = 'Ch'+str(Channels[1]));
	axis([0,fmax,-180,180]);
	phaseresp.yaxis.set_label_position("right");
	xlabel('Frequency [Hz]',fontsize=labelfs); xticks(fontsize=labelfs);
	ylabel('Phase [deg]',fontsize=labelfs); yticks(fontsize=labelfs);
	grid('on');

	#show();

	figTitle = str(day) + ' Run ' + str(runNum)+ ' Noise '+str(noiseNum)+','+str(noiseNum2)+','+str(noiseNum3)+' B-'+str(direc)+': Channels '+str(Channels[0])+', '+str(Channels[1]);
	fname1 = dayDir+'\\Mag2PSN_'+str(direc)+'_'+str(day)+'_'+'R'+str(runNum)+'N'+str(noiseNum)+str(noiseNum2)+str(noiseNum3)+'_'+'Ch'+str(Channels[0])+str(Channels[1]);
	fname2 = runDir+'\\Mag2PSN_'+str(direc)+'_'+str(day)+'_'+'R'+str(runNum)+'N'+str(noiseNum)+str(noiseNum2)+str(noiseNum3)+'_'+'Ch'+str(Channels[0])+str(Channels[1]);
	suptitle(figTitle,fontsize=18);
	figs.savefig(fname2+'.pdf', format = 'pdf');
	figs.savefig(fname1+'.png');

def plotMagRun3(loc='WIMR',day='2016.10.01',runNum='00',noiseNum='00',Chan=1,direc='Y',bs = 1,fignum = 1, ver='v16', fsamp=1000., ds=0, fsampd=1000, trim=.01, HPF = 1, LPF = 80, fan=0, raw = 0, t0 = 0, dt = 10, showfig = 1,ZM=0):
	'''Single Channel. Plots PSD, response curve, entire time series, and time series starting at x0 and ending at x0+dx. Applies bandpass filter from HPF to LPF, and notch filters at 60 Hz and 120 Hz. trim chops off a fraction of the data at each end of the original time series'''
	[basepath, analysispath] = getCompEnv(loc);
	basepath = basepath + ver + '\\';
	dayDir = analysispath+day

	a = os.path.isdir(dayDir);
	if a==0:
		os.mkdir(dayDir);

	runDir = analysispath+day+'\\Run'+runNum
	a = os.path.isdir(runDir);
	if a==0:
		os.mkdir(runDir);

	fname1 = runDir+'\\MagRun3_'+str(day)+'_'+'R'+str(runNum)+'N'+str(noiseNum)+'_Chan'+str(Chan)+direc
	fname2 = dayDir+'\\MagRun3_'+str(day)+'_'+'R'+str(runNum)+'N'+str(noiseNum)+'_Chan'+str(Chan)+direc

	noiseNum = int(noiseNum)
	if noiseNum<10:
		noiseNum = '0'+str(noiseNum);
	
	path = basepath +str(day) + '\\' + 'run_' + str(runNum);
	fPSD = path + '\\'+'noise_'+str(noiseNum)+'\\'+direc+'_magnetic-field-PSD_'+str(Chan)+'.txt';
	ftser = path + '\\'+'noise_'+str(noiseNum)+'\\'+direc+'_calibrated_time_series_fT_'+str(Chan)+'.bin';
	ftserRaw = path + '\\'+'noise_'+str(noiseNum)+'\\'+direc+'_noise_'+str(Chan)+'.bin';
	fRESP =  path+'\\'+'response'+'\\'+direc+'_amp-n-phase_response_'+str(Chan)+'.txt';

	if ZM==0:
		figTitle = str(day)+' Run '+str(runNum)+' Noise '+str(noiseNum)+': Channel '+str(Chan)+' B-'+str(direc)+' [DC-SERF]';
	if ZM==1:
		if direc=='Y':
			figTitle = str(day)+' Run '+str(runNum)+' Noise '+str(noiseNum)+': Channel '+str(Chan)+' B-'+str(direc)+' [Z-mode (DC lowpassed)]';
		if direc=='X':
			figTitle = str(day)+' Run '+str(runNum)+' Noise '+str(noiseNum)+': Channel '+str(Chan)+' B-'+str(direc)+' [Z-mode (1f demodulated)]';

	#loads and bins PSD
	data = transpose(loadtxt(fPSD));
	freq = binit(data[0],bs);
	PS = binit(data[1],bs);
	
	#loads response function
	dataResp = transpose(loadtxt(fRESP));
	
	#loads calibrated time series, chops off ends, and downsamples if required
	dataTS = fromfile(ftser,dtype = '<d');
	sz=size(dataTS);
	chop=sz*trim
	dataTS = dataTS[chop:sz-chop]
	if ds==1:
		r=fsamp/fsampd;
		dataTS=binit(dataTS,r);
		fsamp=fsampd #sets fsamp (which will be used again) equal to the downsampled rate
		
	#loads uncalibrated raw time series and chops off ends
	dataTSRaw = fromfile(ftserRaw,dtype = '<d');
	szRaw=size(dataTSRaw);
	chopRaw=szRaw*trim
	dataTSRaw = dataTSRaw[chopRaw:szRaw-chopRaw]

	#creates array of times based on size of time series and sampling rate;
	sz=size(dataTS);
	if ds==0:
		tx=linspace(0,(sz-1)/fsamp,sz);
	elif ds==1:
		tx=linspace(0,(sz-1)/fsampd,sz);
		
	Nyq=fsamp/2;

	b60,a60 = butter(4,[59./Nyq,61./Nyq],btype='bandstop'); #60 Hz notch filter
	b120,a120 = butter(2,[118./Nyq,122./Nyq],btype='bandstop'); #120 Hz notch filter 
	bLP,aLP = butter(2,[HPF/Nyq, LPF/Nyq], btype = 'bandpass'); #bandpass filter for region of interest
	
	#applies filters to data
	dataTSf = filtfilt(bLP,aLP,dataTS);
	dataTSf = filtfilt(b120,a120,dataTSf);
	dataTSf = filtfilt(b60,a60,dataTSf);
	
	if fan==1:
		bfan,afan = butter(2,[5.5/Nyq,6.5/Nyq],btype='bandstop'); #fan notch filter - no longer needed after fan was removed
		dataTSf = filtfilt(bfan,afan,dataTSf);

	figure(fignum,figsize=(16,10))
	subplot(221);
	semilogy(freq,PS,'k-');
	grid('on',which='both');
	axis([0,200,0.1,10000])
	xlabel('Frequency [Hz]');
	ylabel('PSD [fT/rHz]');
	
	axR = subplot(222);
	ColorMap=['k-','y-','b-','m-','g-'];
	plot(dataResp[0],1E6*dataResp[1],ColorMap[Chan],linewidth=3);
	axis([0,200,0,1.1E6*max(dataResp[1])]);
	xlabel('Frequency [Hz]');
	axR.yaxis.set_label_position("right")
	ylabel('Resp [$\mu$V/fT]');
	grid('on');
	
	#plots entire calibrated time series, unfiltered
	subplot(413)
	plot(tx,1E-3*dataTS, 'k-');
	ylabel('Field [pT]');
	grid('on');
	
	#gets indicies for range of data plotted in 414
	x0=t0*fsamp
	xmax=(t0+dt)*fsamp
	
	subplot(414)
	#plots calibrated time series of interest, filtered
	if(raw==0):
		plot(tx[x0:xmax],1E-3*dataTSf[x0:xmax], ColorMap[Chan]);
		xlim(t0,t0+dt);
		ylabel('Field [pT]')
		grid('on');
	
	elif (raw==1):
	#plots raw, uncalibrated time series of interest, unfiltered
		plot(tx[x0:xmax],dataTSRaw[x0:xmax], ColorMap[Chan])
		grid('on')

	elif(raw==2):
	#plots raw, uncalibrated time series of interest, filtered
		dataTSR = filtfilt(b120,a120,dataTSRaw);
		dataTSR = filtfilt(b60,a60,dataTSR);
		dataTSR = filtfilt(bLP,aLP,dataTSR);
		plot(tx[x0:xmax],dataTSR[x0:xmax], ColorMap[Chan])
		grid('on')

	suptitle(figTitle,fontsize=18)
	savefig(fname1+'.pdf')
	savefig(fname2+'.png')

	if showfig ==1:
		show();	
	
	return dataTS, dataTSf, dataTSRaw;
	
def get4ChanData(loc='WIMR',day='2016.10.01',runNum='00',noiseNum='00',LPF = 80, HPF = 2,f120=120,f60=60, ver = 'v16'):
	[basepath, analysispath] = getCompEnv(loc);
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
	
	path = basepath +str(day) + '\\' + 'run_' + str(runNum);
	ftser1 = path + '\\'+'noise_'+str(noiseNum)+'\\'+'Y_calibrated_time_series_fT_'+str(1)+'.bin';
	ftser2 = path + '\\'+'noise_'+str(noiseNum)+'\\'+'Y_calibrated_time_series_fT_'+str(2)+'.bin';
	ftser3 = path + '\\'+'noise_'+str(noiseNum)+'\\'+'Y_calibrated_time_series_fT_'+str(3)+'.bin';
	ftser4 = path + '\\'+'noise_'+str(noiseNum)+'\\'+'Y_calibrated_time_series_fT_'+str(4)+'.bin';
	
	
	dataTS1 = fromfile(ftser1,dtype = '<d')
	dataTS2 = fromfile(ftser2,dtype = '<d')
	dataTS3 = fromfile(ftser3,dtype = '<d')
	dataTS4 = fromfile(ftser4,dtype = '<d')
	
	Nyq = 500.
	#100 Hz low pass
	b120,a120 = butter(2,[118./Nyq,122./Nyq],btype='bandstop');

	bLP,aLP = butter(4,[HPF/Nyq, LPF/Nyq], btype = 'bandpass');
	bLP,aLP = butter(4,[HPF/Nyq, LPF/Nyq], btype = 'bandpass');
	bf,af = butter(2,[4./Nyq,5./Nyq],btype='bandstop');
	b60,a60 = butter(4,[59./Nyq,61./Nyq],btype='bandstop');

	dataTSf1 = filtfilt(bLP,aLP,dataTS1);
	dataTSf1 = filtfilt(b120,a120,dataTSf1);
	dataTSf1 = filtfilt(b60,a60,dataTSf1);
	

	dataTSf2 = filtfilt(bLP,aLP,dataTS2);
	dataTSf2 = filtfilt(b120,a120,dataTSf2);
	dataTSf2 = filtfilt(b60,a60,dataTSf2);
	
	
	
	dataTSf3 = filtfilt(bLP,aLP,dataTS3);
	dataTSf3 = filtfilt(b120,a120,dataTSf3);
	dataTSf3 = filtfilt(b60,a60,dataTSf3);
	
	
	
	dataTSf4 = filtfilt(bLP,aLP,dataTS4);
	dataTSf4 = filtfilt(b120,a120,dataTSf4);
	dataTSf4 = filtfilt(b60,a60,dataTSf4);
	
	
	return [dataTSf1,dataTSf2,dataTSf3,dataTSf4]

def PlotTimeSeries(loc='WIMR',day='2016.10.01',runNum='00',noiseNum='00',Chan=1,inv=0,ver='v16',leg=1,direc='Y',filt=0,HPF=1,LPF=80,fan=0,t0=0,tlen=1,ymax=0,fsamp=1000):
	[basepath, analysispath] = getCompEnv(loc);
	basepath = basepath + ver + '\\';
	dayDir = analysispath+day
	

	a = os.path.isdir(dayDir);
	if a==0:
		os.mkdir(dayDir);

	runDir = analysispath+day+'\\Run'+runNum
	a = os.path.isdir(runDir);
	if a==0:
		os.mkdir(runDir);

	path = basepath +str(day) + '\\' + 'run_' + runNum + '\\' + 'noise_' + noiseNum + '\\';
	
	#sampling rate = 1000; fNyq = 500 Hz;
	#To create a notch filter at 60 Hz; .116 to .124 (58 Hz to 62 Hz)
	Nyq = fsamp/2

	b120,a120 = butter(2,[118./Nyq,122./Nyq],btype='bandstop');
	bLP,aLP = butter(4,[HPF/Nyq, LPF/Nyq], btype = 'bandpass');
	bf,af = butter(2,[4./Nyq,7./Nyq],btype='bandstop');
	b60,a60 = butter(4,[59./Nyq,61./Nyq],btype='bandstop');
	
	channel=Chan[0]
	filepath=path+direc+'_Calibrated_time_series_fT_'+str(channel)+'.bin'
	dataTS1 = fromfile(filepath,dtype = '<d')
	
	dataTSf1 = filtfilt(bLP,aLP,dataTS1);
	dataTSf1 = filtfilt(b120,a120,dataTSf1);
	dataTSf1 = filtfilt(b60,a60,dataTSf1);
	
	if fan==1:
			dataTSf1 = filtfilt(bf,af,dataTSf1)
	
	#allows inversion of data to align all maternal QRS peaks
	if inv[0]==1:
		dataTSf1=-1*dataTSf1
	
	lw=1
	
	color=['k-','y-','b-','m-','g-']
	
	time=linspace(0,len(dataTSf1)-1,len(dataTSf1))
	time=time/fsamp
	
	figure(fignum,figsize=(14,10));
	plot(time,dataTSf1,color[channel],linewidth=lw);
	xlim(t0,t0+tlen);
	if ymax==0:
		ylim(-1.1*max(dataTSf1),1.1*max(dataTSf1))
	else:
		ylim(-1*ymax,ymax)
	xlabel('Time [s]')
	ylabel('Field [fT]')
	
	grid();

def PlotNTimeSeries(loc='WIMR',day='2016.10.01',runNum='00',noiseNum='00',Chan=[1,2,3,4],inv=[0,0,0,0],bs=4,fignum=1,ver='v16',leg=1,direc='Y',LPF=80,HPF=1,fan=0,t0=0,tlen=1,ymax=0,fsamp=1000):
	[basepath, analysispath] = getCompEnv(loc);
	basepath = basepath + ver + '\\';
	dayDir = analysispath+day
	

	a = os.path.isdir(dayDir);
	if a==0:
		os.mkdir(dayDir);

	runDir = analysispath+day+'\\Run'+runNum
	a = os.path.isdir(runDir);
	if a==0:
		os.mkdir(runDir);

	path = basepath +str(day) + '\\' + 'run_' + runNum + '\\' + 'noise_' + noiseNum + '\\';
	
	#sampling rate = 1000; fNyq = 500 Hz;
	#To create a notch filter at 60 Hz; .116 to .124 (58 Hz to 62 Hz)
	Nyq = fsamp/2

	b120,a120 = butter(2,[118./Nyq,122./Nyq],btype='bandstop');
	bLP,aLP = butter(4,[HPF/Nyq, LPF/Nyq], btype = 'bandpass');
	bf,af = butter(2,[4./Nyq,7./Nyq],btype='bandstop');
	b60,a60 = butter(4,[59./Nyq,61./Nyq],btype='bandstop');
	
	channel=Chan[0]
	filepath=path+direc+'_Calibrated_time_series_fT_'+str(channel)+'.bin'
	dataTS1 = fromfile(filepath,dtype = '<d')
	
	dataTSf1 = filtfilt(bLP,aLP,dataTS1);
	dataTSf1 = filtfilt(b120,a120,dataTSf1);
	dataTSf1 = filtfilt(b60,a60,dataTSf1);
	
	if fan==1:
			dataTSf1 = filtfilt(bf,af,dataTSf1)
	
	#allows inversion of data to align all maternal QRS peaks
	if inv[0]==1:
		dataTSf1=-1*dataTSf1
	
	lw=1
	
	color=['k-','y-','b-','m-','g-']
	
	time=linspace(0,len(dataTSf1)-1,len(dataTSf1))
	time=time/fsamp
	
	figure(fignum,figsize=(14,10));
	plot(time,dataTSf1,color[channel],linewidth=lw);
	xlim(t0,t0+tlen);
	if ymax==0:
		ylim(-1.1*max(dataTSf1),1.1*max(dataTSf1))
	else:
		ylim(-1*ymax,ymax)
	xlabel('Time [s]')
	ylabel('Field [fT]')
	
	grid();

	
	if len(Chan) > 1:
		channel=Chan[1]
		filepath=path+direc+'_Calibrated_time_series_fT_'+str(channel)+'.bin'
		dataTS2 = fromfile(filepath,dtype = '<d')
		
		##sampling rate = 1000; fNyq = 500 Hz;
		##To create a notch filter at 60 Hz; .116 to .124 (58 Hz to 62 Hz)
		#Nyq = fsamp/2
	    #
		#b120,a120 = butter(2,[118./Nyq,122./Nyq],btype='bandstop');
		#bLP,aLP = butter(4,[HPF/Nyq, LPF/Nyq], btype = 'bandpass');
		#bf,af = butter(2,[4./Nyq,5./Nyq],btype='bandstop');
		#b60,a60 = butter(4,[59./Nyq,61./Nyq],btype='bandstop');
	
		dataTSf2 = filtfilt(bLP,aLP,dataTS2);
		dataTSf2 = filtfilt(b120,a120,dataTSf2);
		dataTSf2 = filtfilt(b60,a60,dataTSf2);
		
		if fan==1:
			dataTSf2 = filtfilt(bf,af,dataTSf2)
		
		if inv[1]==1:
			dataTSf2=-1*dataTSf2
			
		plot(time,dataTSf2,color[channel-1],linewidth=lw);
		
	if len(Chan) > 2:
		channel=Chan[2]
		filepath=path+direc+'_Calibrated_time_series_fT_'+str(channel)+'.bin'
		dataTS3 = fromfile(filepath,dtype = '<d')
		
		##sampling rate = 1000; fNyq = 500 Hz;
		##To create a notch filter at 60 Hz; .116 to .124 (58 Hz to 62 Hz)
		#Nyq = fsamp/2
	    #
		#b120,a120 = butter(2,[118./Nyq,122./Nyq],btype='bandstop');
		#bLP,aLP = butter(4,[HPF/Nyq, LPF/Nyq], btype = 'bandpass');
		#bf,af = butter(2,[4./Nyq,5./Nyq],btype='bandstop');
		#b60,a60 = butter(4,[59./Nyq,61./Nyq],btype='bandstop');
	
		dataTSf3 = filtfilt(bLP,aLP,dataTS3);
		dataTSf3 = filtfilt(b120,a120,dataTSf3);
		dataTSf3 = filtfilt(b60,a60,dataTSf3);
		
		if fan==1:
			dataTSf3 = filtfilt(bf,af,dataTSf3)
		
		if inv[2]==1:
			dataTSf3=-1*dataTSf3
			
		plot(time,dataTSf3,color[channel-1],linewidth=lw);

	if len(Chan) > 3:
		channel=Chan[3]
		filepath=path+direc+'_Calibrated_time_series_fT_'+str(channel)+'.bin'
		dataTS4 = fromfile(filepath,dtype = '<d')
		
		##sampling rate = 1000; fNyq = 500 Hz;
		##To create a notch filter at 60 Hz; .116 to .124 (58 Hz to 62 Hz)
		#Nyq = fsamp/2
	    #
		#b120,a120 = butter(2,[118./Nyq,122./Nyq],btype='bandstop');
		#bLP,aLP = butter(4,[HPF/Nyq, LPF/Nyq], btype = 'bandpass');
		#bf,af = butter(2,[4./Nyq,5./Nyq],btype='bandstop');
		#b60,a60 = butter(4,[59./Nyq,61./Nyq],btype='bandstop');
	
		dataTSf4 = filtfilt(bLP,aLP,dataTS4);
		dataTSf4 = filtfilt(b120,a120,dataTSf4);
		dataTSf4 = filtfilt(b60,a60,dataTSf4);
		
		if fan==1:
			dataTSf4 = filtfilt(bf,af,dataTSf4)

		if inv[3]==1:
			dataTSf4=-1*dataTSf4
			
		plot(time,dataTSf4,color[channel-1],linewidth=lw);
		
	return dataTS1,dataTSf1

def PlotNPSDs(loc='WIMR',day='2016.10.01',runNum='00',noiseNum='00',Chan=[1,2,3,4],bs=4,fignum=1,ver='v16',leg=1,direc='Y',fmax=200):
	[basepath, analysispath] = getCompEnv(loc);
	basepath = basepath + ver + '\\';
	dayDir = analysispath+day
	
	a = os.path.isdir(dayDir);
	if a==0:
		os.mkdir(dayDir);

	runDir = analysispath+day+'\\Run'+runNum
	a = os.path.isdir(runDir);
	if a==0:
		os.mkdir(runDir);
	
	channel=Chan[0]
	path = basepath +str(day) + '\\' + 'run_' + str(runNum) + '\\' + 'noise_' + noiseNum + '\\';
	fPSD1 = path+direc+'_magnetic-field-PSD_'+str(channel)+'.txt';
	data1=transpose(loadtxt(fPSD1))
	
	freq1 = binit(data1[0],bs);
	PSD1 = binit(data1[1],bs);
	
	fig = figure(fignum,figsize=(16.0,10.0))
	suptitle('Power Spectral Densities, '+day+' Run '+runNum+' Noise '+noiseNum,fontsize=18)
	
	subplot(221)
	semilogy(freq1,PSD1,'k-',label='Channel'+str(channel));
	xlim(0,fmax);
	ylim(1,1E5)
	grid(which='both');
	title('Channel '+str(channel));
	
	if len(Chan) > 1:
		channel=Chan[1]
		path = basepath +str(day) + '\\' + 'run_' + str(runNum) + '\\' + 'noise_' + noiseNum + '\\';
		fPSD2 = path+direc+'_magnetic-field-PSD_'+str(channel)+'.txt';
		data2=transpose(loadtxt(fPSD2))
		
		freq2 = binit(data2[0],bs);
		PSD2 = binit(data2[1],bs);
		
		subplot(222);
		semilogy(freq2,PSD2,'k-',label='Channel'+str(channel));
		xlim(0,fmax);
		ylim(1,1E5)
		grid(which='both');
		title('Channel '+str(channel));
	
	if len(Chan) > 2:
		channel=Chan[2]
		path = basepath +str(day) + '\\' + 'run_' + str(runNum) + '\\' + 'noise_' + noiseNum + '\\';
		fPSD3 = path+direc+'_magnetic-field-PSD_'+str(channel)+'.txt';
		data3=transpose(loadtxt(fPSD3))
		
		freq3 = binit(data3[0],bs);
		PSD3 = binit(data3[1],bs);
		
		subplot(223);
		semilogy(freq3,PSD3,'k-',label='Channel'+str(channel));
		xlim(0,fmax);
		ylim(1,1E5)
		grid(which='both');
		title('Channel '+str(channel));
		
	if len(Chan) > 3:
		channel=Chan[3]
		path = basepath +str(day) + '\\' + 'run_' + str(runNum) + '\\' + 'noise_' + noiseNum + '\\';
		fPSD4 = path+direc+'_magnetic-field-PSD_'+str(channel)+'.txt';
		data4=transpose(loadtxt(fPSD4))
		
		freq4 = binit(data4[0],bs);
		PSD4 = binit(data4[1],bs);
		
		subplot(224);
		semilogy(freq4,PSD4,'k-',label='Channel'+str(channel));
		xlim(0,fmax);
		ylim(1,1E5)
		grid(which='both');
		title('Channel '+str(channel));

def magFilt(d1,fsamp=1000,fanfreq = 6.0, fanfilt = 0):
	Nyq = fsamp/2

	b0,a0 = butter(4,[1./Nyq,80./Nyq],btype='bandpass');
	b60,a60 = butter(4,[59./Nyq,61/Nyq],btype='bandstop');
	bf,af = butter(2,[(fanfreq-0.75)/Nyq,(fanfreq+0.75)/Nyq],btype='bandstop');

	d1 = filtfilt(b0,a0,d1);
	d1 = filtfilt(b0,a0,d1);
	d1 = filtfilt(b0,a0,d1);
	d1 = filtfilt(b60,a60,d1);
	#~ d1 = filtfilt(b60,a60,d1);
	if fanfilt == 1:
		d1 = filtfilt(bf,af,d1);
		#~ d1 = filtfilt(bf,af,d1);
	
	return d1
	
def plotMagRunSlider(loc='WIMR',day='2016.10.01',runNum='00',noiseNum='00',Chan=1,direc='Y',fsamp=1000.,fignum=1,ver='v16',Yscale=2000,dX=1000,lnw=1.5):
	[basepath, analysispath] = getCompEnv(loc);
	basepath = basepath + ver + '\\';
	path = basepath +str(day) + '\\' + 'run_' + str(runNum);
	
	figTitle = str(day) + '\\' + 'run_' + str(runNum)+'\\'+'noise_'+str(noiseNum);
	
	#time series file
	ftser1 = path + '\\'+'noise_'+str(noiseNum)+'\\'+direc+'_calibrated_time_series_fT_'+str(Chan)+'.bin';
	
	#import time series
	d1 = fromfile(ftser1,dtype = '<d')
	
	#filter time series
	Nyq=fsamp/2
	b0,a0 = butter(2,[1./Nyq,80./Nyq],btype='bandpass');
	b60,a60 = butter(4,[59./Nyq,61/Nyq],btype='bandstop');
	
	d1 = filtfilt(b0,a0,d1);
	d1 = filtfilt(b60,a60,d1);
	
	#plot data
	figure(fignum,figsize=[16,8]);
	clf();
	
	colormap=['k-','y-','b-','m-','g-'];
	
	ax1 = subplot(111);
	l1, = plot(d1, colormap[Chan],lw=lnw);
	axis([0, dX,-Yscale,Yscale]);
	grid('on');
	
	suptitle(figTitle);
	
	xMin = axes([0.15, 0.01, 0.65, 0.03])
	sX = Slider(xMin,'Xmin', 0, size(d1), valinit=0)
	
	def update(val):
		x0 = sX.val
		ax1.set_xlim(x0,x0+dX);        
		ax1.set_ylim(auto=True);
		draw()
	
	sX.on_changed(update)
	show()

def plotMagRunSlider2(loc='WIMR',day='2016.10.01',runNum='00',noiseNum = '00',Chan=[1,2,3,4],direc='Y',Fsamp=1000., fignum = 1,ver = 'v16',Yscale=2000,dX = 1000,lw=1.5):
	[basepath, analysispath] = getCompEnv(loc);
	basepath = basepath + ver + '\\';
	path = basepath +str(day) + '\\' + 'run_' + str(runNum);
	
	ftserdata=[]
	
	for channel in Chan:
		ftser=path+'\\noise_'+str(noiseNum)+'\\'+direc+'_calibrated_time_series_fT_'+str(channel)+'.bin';
		d=fromfile(ftser,dtype = '<d');
		d=magFilt(d,fsamp=Fsamp,fanfreq=6,fanfilt=0)
		ftserdata.append(d);
	
	colormap=['k-','y-','b-','m-','g-']
	
	timeseries=figure(fignum,figsize=[16,8]);
	
	xMin  = axes([0.15, 0.01, 0.65, 0.03])
	sX = Slider(xMin,'Xmin', 0, size(ftserdata[0])-dX, valinit=0);
	
	if len(Chan)==1:
		ax1=subplot(111);
		l1, = plot(ftserdata[0],colormap[Chan[0]],linewidth=lw);
		axis([0, dX,-Yscale,Yscale]);
		grid('on');
		
		def update(val):
			x0 = sX.val
			ax1.set_xlim(x0,x0+dX); 			
			ax1.set_ylim(auto=True);
			
			draw()

		sX.on_changed(update)
		show()
		
	elif len(Chan)==2:
		ax1=subplot(211);
		plot(ftserdata[0],colormap[Chan[0]],linewidth=lw);
		axis([0, dX,-Yscale,Yscale]);
		grid('on');
		
		ax2=subplot(212);
		plot(ftserdata[1],colormap[Chan[1]],linewidth=lw);
		axis([0, dX,-Yscale,Yscale]);
		grid('on');
		
		def update(val):
			x0 = sX.val
			ax1.set_xlim(x0,x0+dX); 			
			ax1.set_ylim(auto=True);
			
			ax2.set_xlim(x0,x0+dX); 			
			ax2.set_ylim(auto=True);
			
			draw()

		sX.on_changed(update)
		show()
		
	elif len(Chan)==3:
		ax1=subplot(311);
		plot(ftserdata[0],colormap[Chan[0]],linewidth=lw);
		axis([0, dX,-Yscale,Yscale]);
		grid('on');
		
		ax2=subplot(312);
		plot(ftserdata[1],colormap[Chan[1]],linewidth=lw);
		axis([0, dX,-Yscale,Yscale]);
		grid('on');
		
		ax3=subplot(313);
		plot(ftserdata[2],colormap[Chan[2]],linewidth=lw);
		axis([0, dX,-Yscale,Yscale]);
		grid('on');
		
		def update(val):
			x0 = sX.val
			ax1.set_xlim(x0,x0+dX); 			
			ax1.set_ylim(auto=True);
			
			ax2.set_xlim(x0,x0+dX); 			
			ax2.set_ylim(auto=True);
			
			ax3.set_xlim(x0,x0+dX); 			
			ax3.set_ylim(auto=True);
			
			draw()

		sX.on_changed(update);
		show();
		
	elif len(Chan)==4:
		ax1=subplot(411);
		plot(ftserdata[0],colormap[Chan[0]],linewidth=lw);
		axis([0, dX,-Yscale,Yscale]);
		grid('on');
		
		ax2=subplot(412);
		plot(ftserdata[1],colormap[Chan[1]],linewidth=lw);
		axis([0, dX,-Yscale,Yscale]);
		grid('on');
		
		ax3=subplot(413);
		plot(ftserdata[2],colormap[Chan[2]],linewidth=lw);
		axis([0, dX,-Yscale,Yscale]);
		grid('on');
		
		ax4=subplot(414);
		plot(ftserdata[3],colormap[Chan[3]],linewidth=lw);
		axis([0, dX,-Yscale,Yscale]);
		grid('on');
		
		def update(val):
			x0 = sX.val
			ax1.set_xlim(x0,x0+dX); 			
			ax1.set_ylim(auto=True);
			
			ax2.set_xlim(x0,x0+dX); 			
			ax2.set_ylim(auto=True);
			
			ax3.set_xlim(x0,x0+dX); 			
			ax3.set_ylim(auto=True);

			ax4.set_xlim(x0,x0+dX); 			
			ax4.set_ylim(auto=True);
			
			draw()

		sX.on_changed(update)
		show()
		
def getMagRunSlider4(loc='WIMR',day='2016.10.01',runNum='00',noiseNum='00',bs = 2,notch = 0, ffreq = 6, fignum = 1,ver = 'v15',Yscale=2000,Grad=0,dX = 1000,lnw=1.5):
    [basepath, analysispath] = getCompEnv(loc);
    basepath = basepath + ver + '\\';
    path = basepath +str(day) + '\\' + 'run_' + str(runNum);
    
    figTitle = str(day) + '\\' + 'run_' + str(runNum)+'\\'+'noise_'+str(noiseNum);

    ftser1 = path + '\\'+'noise_'+str(noiseNum)+'\\'+'Y_calibrated_time_series_fT_'+str(1)+'.bin';
    ftser2 = path + '\\'+'noise_'+str(noiseNum)+'\\'+'Y_calibrated_time_series_fT_'+str(2)+'.bin';
    ftser3 = path + '\\'+'noise_'+str(noiseNum)+'\\'+'Y_calibrated_time_series_fT_'+str(3)+'.bin';
    ftser4 = path + '\\'+'noise_'+str(noiseNum)+'\\'+'Y_calibrated_time_series_fT_'+str(4)+'.bin';

    d1 = fromfile(ftser1,dtype = '<d')
    d2 = fromfile(ftser2,dtype = '<d')
    d3 = fromfile(ftser3,dtype = '<d')
    d4 = fromfile(ftser4,dtype = '<d')


    d1 = magFilt(d1,fsamp=1000,fanfreq = ffreq, fanfilt = notch)
    d2 = magFilt(d2,fsamp=1000,fanfreq = ffreq, fanfilt = notch)
    d3 = magFilt(d3,fsamp=1000,fanfreq = ffreq, fanfilt = notch)
    d4 = magFilt(d4,fsamp=1000,fanfreq = ffreq, fanfilt = notch)

    DS = np.array([d1,d2,d3,d4])

    a = shape(DS);
    if a[0]!=4:
        DS = transpose(DS);


    figure(fignum,figsize=[16,8]);
    clf()

    ax1 = subplot(411)
    l1, = plot(DS[0], 'y-',lw=lnw)
    axis([0, 10000,-Yscale,Yscale])
    grid('on')

    ax2 = subplot(412)
    l2, = plot(DS[1], 'b-',lw=lnw)
    axis([0, 10000,-Yscale,Yscale])
    grid('on')

    ax3 = subplot(413)
    l3, = plot(DS[2], 'm-',lw=lnw)
    axis([0, 10000,-Yscale,Yscale])
    grid('on')

    ax4 = subplot(414)
    l4, = plot(DS[3],'g-', lw=lnw)
    axis([0, 10000,-Yscale,Yscale])
    grid('on')
    
    suptitle(figTitle)


    
    xMin  = axes([0.15, 0.01, 0.65, 0.03])
    sX = Slider(xMin,'Xmin', 0, size(DS[0]), valinit=0)

    def update(val):
        x0 = sX.val
        
        ax1.set_xlim(x0,x0+dX); 
        ax2.set_xlim(x0,x0+dX);
        ax3.set_xlim(x0,x0+dX);
        ax4.set_xlim(x0,x0+dX);
        
        
        ax1.set_ylim(auto=True);
        ax2.set_ylim(auto=True);
        ax3.set_ylim(auto=True);
        ax4.set_ylim(auto=True);
        
        draw()

    sX.on_changed(update)
    show()
	
def AmpPhaseCompare(loc='WIMR',day='2016.10.01',run='00',Chan=1,direc='Y',label='',ver='v16',fmax=200,radresp=0,calib=2E-6,Ipr=100E-6):
	'''Plots amplitude and phase response for one run. Desgined to be "stackable": Run multiple fucntions sequentially and all will be plotted on same plot.'''
	import time
	today=time.strftime('%Y.%m.%d')
	
	[basepath, analysispath] = getCompEnv(loc);
	basepath = basepath + ver + '\\';
		
	daydir=analysispath+today+'\\'
	
	AmpPhase=transpose(loadtxt(basepath+day+'\\'+'run_'+run+'\\response\\'+direc+'_amp-n-phase_response_'+str(Chan)+'.txt'))
	
	matplotlib.rcParams['figure.figsize'] = (16,8.0)
	
	if radresp==1:
		AmpPhase[1]=AmpPhase[1]*calib/(4*Ipr)
		subplot(121);
		plot(AmpPhase[0],1E7*AmpPhase[1],linewidth=2,label=label);
		xlabel('Frequency [Hz]');
		ylabel('Response x 10$^{-7}$ [rad/fT]')
		xlim(0,fmax);
		grid('on')
		if label!='':
			legend()
	elif radresp==0:
		subplot(121);
		plot(AmpPhase[0],1E6*AmpPhase[1],linewidth=2,label=label);
		xlabel('Frequency [Hz]');
		ylabel('Response [$\mu$V/fT]')
		xlim(0,fmax);
		grid('on')
		if label!='':
			legend()
	
	subplot(122);
	plot(AmpPhase[0],AmpPhase[2],linewidth=2,label=label);
	xlabel('Frequency [Hz]');
	ylabel('Phase [deg]')
	xlim(0,fmax);
	#ylim(-180,180);
	grid('on')
	if label!='':
		legend()

def NoiseCompare(loc='WIMR',day='2016.10.01',run='00',noise='00',bs=1,Chan=1,direc='Y',label='',ver='v16',fmax=200,radresp=0,calib=2E-6,Ipr=100E-6,psn=0):
	'''Plots response and noise (in fT/rHz) for one run. Desgined to be "stackable": Run multiple fucntions in one input box and all will be plotted on same plot.'''
	import time
	today=time.strftime('%Y.%m.%d')
	
	[basepath, analysispath] = getCompEnv(loc);
	basepath = basepath + ver + '\\';
		
	daydir=analysispath+today+'\\'
	
	AmpPhase=transpose(loadtxt(basepath+day+'\\'+'run_'+run+'\\response\\'+direc+'_amp-n-phase_response_'+str(Chan)+'.txt'))
	Noise=transpose(loadtxt(basepath+day+'\\'+'run_'+run+'\\noise_'+noise+'\\'+direc+'_magnetic-field-PSD_'+str(Chan)+'.txt'))
	
	freq = binit(Noise[0],bs);
	PS = binit(Noise[1],bs);
	
	if psn==1:
		shotnoiseA = sqrt(4*Ipr*1.6E-19);   #measured in A/rHz
		shotnoiseV = shotnoiseA/calib    #measured in V/rHz
		shotnoiseB = shotnoiseV/AmpPhase[1];  #measured in fT/rHz
	
	fs=14
	matplotlib.rcParams['figure.figsize'] = (16,8.0);
	
	if radresp==1:
		AmpPhase[1]=AmpPhase[1]*calib/(4*Ipr)
		subplot(121);
		plot(AmpPhase[0],1E7*AmpPhase[1],linewidth=2,label=label);
		xlabel('Frequency [Hz]',fontsize=fs); xticks(fontsize=fs);
		ylabel('Response x 10$^{-7}$ [rad/fT]',fontsize=fs); yticks(fontsize=fs);
		xlim(0,fmax);
		grid('on')
		if label!='':
			legend()
	elif radresp==0:
		subplot(121);
		plot(AmpPhase[0],1E6*AmpPhase[1],linewidth=2,label=label);
		xlabel('Frequency [Hz]',fontsize=fs); xticks(fontsize=fs);
		ylabel('Response [$\mu$V/fT]',fontsize=fs); yticks(fontsize=fs);
		xlim(0,fmax);
		grid('on')
		if label!='':
			legend()
	
	subplot(122);
	semilogy(freq,PS,linewidth=1,label=label);
	if psn==1:
		semilogy(AmpPhase[0],shotnoiseB,linestyle='-',color=(0,1,0), linewidth = 3,label='photon shot noise');
	xlabel('Frequency [Hz]',fontsize=fs); xticks(fontsize=fs);
	ylabel('PSD [fT/rHz]',fontsize=fs); yticks(fontsize=fs);
	xlim(0,fmax);
	grid('on',which='both')
	if label!='':
		legend()