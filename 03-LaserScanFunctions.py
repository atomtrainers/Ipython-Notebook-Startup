import os;
	
def ResponseCompare(dayA='2015.10.01',runA='00',ChanA=1,dirA='Y',calibA=5E-6,IprA=100E-6,dayB='2015.01.01',runB='00',ChanB=1,dirB='Y',calibB=5E-6,IprB=100E-6,ver='v16',fmax=200):
	[basepath, analysispath] = getCompEnv();
	basepath = basepath + ver + '\\';

	pathA = basepath +str(dayA) + '\\' + 'run_' + str(runA);
	fRESPA =  pathA+'\\'+'response'+'\\'+'Y_amp-n-phase_response_'+str(ChanA)+'.txt';
	dataA = transpose(loadtxt(fRESPA));
	
	pathB = basepath +str(dayB) + '\\' + 'run_' + str(runB);
	fRESPB =  pathB+'\\'+'response'+'\\'+'Y_amp-n-phase_response_'+str(ChanB)+'.txt';
	dataB = transpose(loadtxt(fRESPB));
	
	labelA=dayA+' - Run '+runA+' - Chan '+str(ChanA);
	labelB=dayB+' - Run '+runB+' - Chan '+str(ChanB);
	
	respAvolts=dataA[1]; #response in V/fT
	respAamps=respAvolts*calibA; #response in A/fT
	respAradians=respAamps/(4*IprA); #response in rad/fT
	
	respBvolts=dataB[1]; #response in V/fT
	respBamps=respBvolts*calibB; #response in A/fT
	respBradians=respBamps/(4*IprB); #response in rad/fT
	
	figure(1,figsize=(8,6))
	plot(dataA[0],1E7*respAradians,'r-',linewidth=2, label = labelA);
	plot(dataB[0],1E7*respBradians,'b-',linewidth=2, label = labelB);
	xlim(xmax = fmax);
	xlabel('frequency (Hz)');
	ylim(ymax = 1.1E7*max(max(respAradians),max(respBradians)));
	ylabel('Response [1E-7 rad/fT]');
	grid(which='both');
	legend();
		
def LaserScan2(date = '2016.10.01', laser = 'Probe', fnum=0, fignum = 1, gr='res'):
	'''Plots four channels of data taken with the Laser Scan v1.1 program. the gr variable will set the x grid with either the resistance value ('res') or the approximate detuning ('detun').'''
	if fnum<10:
		fnum = '000'+str(fnum);
	elif fnum<100:
		fnum = '00'+str(fnum);
	elif fnum<1000:
		fnum = '0'+str(fnum);
	
	[basepath, analysispath] = getCompEnv();
	dayDir = analysispath+date;
	a = os.path.isdir(dayDir);
	if a==0:
		os.mkdir(dayDir);

	ftit = date+' - '+laser+'Scan'+fnum;
	
	datafile = 'C:\\Users\\sulai\\Documents\\LabVIEW Data\\DATA\\LaserScanData\\'+date+'\\'+laser+'Scan'+fnum+'.dat';
	
	data = loadtxt(datafile);
	
	data0p=np.delete(data[0],0);
	data1p=np.delete(data[1],0);
	data2p=np.delete(data[2],0);
	data3p=np.delete(data[3],0);
	data4p=np.delete(data[4],0);
	
	xbot=round(min(data0p),1);
	xtop=round(max(data0p),1);
	
	pngname=dayDir+'\\'+date+'_'+laser+'Scan'+str(fnum)+'.png';
	
	figure(fignum, figsize = [16,10]);
	ax1=subplot(111);
	plot(data0p,data1p, 'y-', linewidth = '3',label='Channel 1');
	plot(data0p,data2p, 'b-', linewidth = '3',label='Channel 2');
	plot(data0p,data3p, 'm-', linewidth = '3',label='Channel 3');
	plot(data0p,data4p, 'g-', linewidth = '3',label='Channel 4');
	title(ftit, fontsize=20, y=1.08);
	xlim(xbot,xtop);
	xticks(np.arange(xbot, xtop, 0.2),fontsize=16);
	xlabel('TEC Temp[kOhm]',fontsize=16);
	ylabel('Signal [V]',fontsize=16); yticks(fontsize=16);
	ylim(ymin=-.1);
	if gr=='res':
		grid('on');
	legend(loc=0);
	
	a = ax1.axis();
	
	ax2=ax1.twiny();
	if laser=='Probe':
		xmin=(a[0]-7.02)*80; xmax=(a[1]-7.05)*80;
	elif laser=='Pump':
		xmin=(a[0]-8.32)*80; xmax=(a[1]-8.42)*80;
	ax2.axis([xmin,xmax,a[2],a[3]]);
	ax2.set_xlabel('Detuning [GHz]',fontsize=16); yticks(fontsize=14);
	xticks(fontsize=16);
	if gr=='detun':
		grid('on');
	
	savefig(pngname)
	
	
	return data0p,data1p,data2p,data3p,data4p;

def LaserScan1Chan(date='2016.10.01',laser='Probe',fnum=0,Chan=1,curr=0,srscalib=100E-6,ymax=14,step=.2,sc=1,lbl='-'):
	'''Plots 1 channel data from the LaserScan VI. Designed to be "stackable". Call many instances in one notebook command line and the data will be shown on same plot. The curr flag converts voltage readings to currents using the gain on the SRS. The sc parameter allows for scaling.'''
	if fnum<10:
			fnum = '000'+str(fnum);
	elif fnum<100:
		fnum = '00'+str(fnum);
	elif fnum<1000:
		fnum = '0'+str(fnum);
	
	if lbl=='-':
			label=date+'_'+laser+'Scan'+str(fnum)
	else:
			label=lbl
	
	datafile = 'C:\\Users\\sulai\\Documents\\LabVIEW Data\\DATA\\LaserScanData\\'+date+'\\'+laser+'Scan'+fnum+'.dat'
	data = loadtxt(datafile);
	
	detun=np.delete(data[0],0); #deletes sometimes-problematic first point of datasets
	sig=np.delete(data[Chan],0);
	
	if curr==1:
		sig=1E6*srscalib*sig #converts to micro-amps
	
	xbot=round(min(detun),1);
	xtop=round(max(detun),1);
	
	matplotlib.rcParams['figure.figsize'] = (16.0,10.0)
	
	plot(detun,sc*sig, linewidth = '3',label=lbl);
	xlim(xbot,xtop)
	xticks(np.arange(xbot, xtop, step),fontsize=16)
	xlabel('TEC Temp ($ \propto \Delta$) [kOhm]',fontsize=16); 
	if curr==1:
		ylabel('Signal [$\mu$A]',fontsize=16); yticks(fontsize=16);
	else:
		ylabel('Signal [V]',fontsize=16); yticks(fontsize=16);

	ylim(ymin=-.1,ymax=ymax)
	grid('on');
	legend(loc=0);
	
def AbsScanWithRef(day,refnum,absnum):
	def AddZeroes(fnum):
		if fnum<10:
			fnum = '000'+str(fnum);
		elif fnum<100:
			fnum = '00'+str(fnum);
		elif fnum<1000:
			fnum = '0'+str(fnum);
		return fnum
	
	refnum=AddZeroes(refnum);
	absnum=AddZeroes(absnum);
	
	def LoadLaserScan2(day,fnum):
		file = day+'ProbeScan'+str(fnum);
		path = 'C:\\Users\\sulai\\Documents\\LabVIEW Data\\DATA\\LaserScanLog\\'+file+'.dat';
		data=loadtxt(path);
		detun=np.delete(data[0],0);
		sig1=np.delete(data[1],0);
		sig2=np.delete(data[2],0);
		sig3=np.delete(data[3],0);
		sig4=np.delete(data[4],0);
		
		return detun,sig1,sig2,sig3,sig4
		
	refdetun,ref1,ref2,ref3,ref4=LoadLaserScan2(day,refnum)
	absdetun,abs1,abs2,abs3,abs4=LoadLaserScan2(day,absnum)
	
	norm1=abs1/ref1
	norm2=abs2/ref2
	norm3=abs3/ref3
	norm4=abs4/ref4
	
	figure(1, figsize = [16,10])
	clf();
	plot(refdetun,norm1, 'y-', linewidth = '3');
	plot(refdetun,norm2, 'b-', linewidth = '3');
	plot(refdetun,norm3, 'm-', linewidth = '3');
	plot(refdetun,norm4, 'g-', linewidth = '3');
	xlim(xmin=min(refdetun)-.1,xmax=max(refdetun)+.1)
	xlabel('TEC Temp ($ \propto \Delta$) [kOhm]',fontsize=14);
	ylabel('Transmission',fontsize=14);
	ylim(ymin = -.1)
	grid(which='both');
	show()
	
	