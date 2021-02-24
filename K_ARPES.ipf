#pragma rtGlobals=3		// Use modern global access method and strict wave access.

//	K_ARPES for ARPES simulation in k-space with rotation
//			by Takashi Kono

Function K_a2k(inputwave,kdim)
	string inputwave
	variable kdim  // 0 : return kx, 1 : return ky, 2 : return kz
	
	// make/N=22 inputwave=NaN
	
	// inputwave[0] = alpha (deg)
	// inputwave[1] = beta (deg)
	//
	// inputwave[2] = analyser rotation (deg)
	// inputwave[3] = polar (deg)
	// inputwave[4] = tilt (deg)
	// inputwave[5] = azimuth (deg)
	//
	// inputwave[6] = polar normal pos (deg)
	// inputwave[7] = tilt normal pos (deg)
	// inputwave[8] = azimuth normal pos (deg)
	//
	// inputwave[9] = polar ofs (deg)
	// inputwave[10] = tilt ofs (deg)
	// inputwave[11] = azimuth ofs (deg)
	//
	// inputwave[12] = alpha reverse ; ON : -1 , OFF = 1
	// inputwave[13] = beta reverse ; ON : -1 , OFF = 1
	// inputwave[14] = polar reverse ; ON : -1 , OFF = 1
	// inputwave[15] = tilt reverse ; ON : -1 , OFF = 1
	// inputwave[16] = azimuth reverse ; ON : -1 , OFF = 1
	//
	// inputwave[17] = Anlular notation ; Polar-ang : 0 , Tilt-ang : 1
	//
	// inputwave[18] = Photon energy (eV)
	// inputwave[19] = Inner potential (eV)
	// inputwave[20] = Work function (eV)
	// inputwave[21] = Binding energy (eV)
	
	
	if(Exists(inputwave))
		wave iw=$inputwave
		
		string nf=GetDataFolder(1)
		
		SetDataFolder root:K_ARPES:misc
		variable ng=0
		
		if(numtype(iw[0])==2)
			print "Alpha parameter is necessary."
			ng=1
		else
			variable ang=iw[0]
		endif
		
		if(numtype(iw[1])==0)
			variable/g v_bet=iw[1]
		endif
		
		if(numtype(iw[2])==0)
			variable/g v_anal_rot=iw[2]
		endif
		
		if(numtype(iw[3])==0)
			variable/g v_po=iw[3]
		endif
		
		if(numtype(iw[4])==0)
			variable/g v_ph=iw[4]
		endif
		
		if(numtype(iw[5])==0)
			variable/g v_az=iw[5]
		endif
		
		if(numtype(iw[6])==0)
			variable/g v_poo=iw[6]
		endif
		
		if(numtype(iw[7])==0)
			variable/g v_pho=iw[7]
		endif
		
		if(numtype(iw[8])==0)
			variable/g v_azo=iw[8]
		endif
		
		if(numtype(iw[9])==0)
			variable/g v_po_ofs=iw[9]
		endif
		
		if(numtype(iw[10])==0)
			variable/g v_ph_ofs=iw[10]
		endif
		
		if(numtype(iw[11])==0)
			variable/g v_az_ofs=iw[11]
		endif
		
		if(numtype(iw[12])==0)
			if(iw[12]==1 || iw[12]==-1)
				variable/g v_revth=iw[12]
			else
				print "Reverse parameter shuld be 1 or -1."
				ng=1
			endif
		endif
		
		if(numtype(iw[13])==0)
			if(iw[13]==1 || iw[13]==-1)
				variable/g v_revbet=iw[13]
			else
				print "Reverse parameter shuld be 1 or -1."
				ng=1
			endif
		endif
		
		if(numtype(iw[14])==0)
			if(iw[14]==1 || iw[14]==-1)
				variable/g v_revpo=iw[14]
			else
				print "Reverse parameter shuld be 1 or -1."
				ng=1
			endif
		endif
		
		if(numtype(iw[15])==0)
			if(iw[15]==1 || iw[15]==-1)
				variable/g v_revph=iw[15]
			else
				print "Reverse parameter shuld be 1 or -1."
				ng=1
			endif
		endif
		
		if(numtype(iw[16])==0)
			if(iw[16]==1 || iw[16]==-1)
				variable/g v_revaz=iw[16]
			else
				print "Reverse parameter shuld be 1 or -1."
				ng=1
			endif
		endif
		
		if(numtype(iw[17])==0)
			if(iw[17]==1 || iw[17]==0)
				variable/g v_analtype=iw[17]
			else
				print "Angular notation shuld be 0 or 1."
				ng=1
			endif
		endif
		
		if(numtype(iw[18])==0)
			if(iw[18]>0)
				variable/g v_hn=iw[18]
			else
				print "Photon energy shuld be lager than 0."
				ng=1
			endif
		endif
		
		if(numtype(iw[19])==0)
			if(iw[19]>0)
				variable/g v_V0=iw[19]
			else
				print "Inner potential shuld be lager than 0."
				ng=1
			endif
		endif
		
		if(numtype(iw[20])==0)
			if(iw[20]>0)
				variable/g v_W=iw[20]
			else
				print "Work function shuld be lager than 0."
				ng=1
			endif
		endif
		
		if(numtype(iw[21])==0)
			variable/g v_EB=iw[21]
		endif
		
		if(strlen(WinList("K_ARPES_p", ";", "WIN:64")))
			K_panel_update()
		endif
		
		if(ng==0)
			
			K_calc_normal()
			K_calc_emis(ang)
			
			nvar smh=v_smh
			variable/g v_hn
			nvar V0=v_V0
			variable/g v_W
			variable/g v_EB
			variable EK=v_hn-v_W-v_EB
			wave vw=$("vec")
		//	print vw[0],vw[1],vw[2]
			variable kx,ky,kz
			
			if(EK<0)
				print "WARNING : EK < 0"
			endif
			
			kx=smh*sqrt(EK)*vw[0]
			ky=smh*sqrt(EK)*vw[1]
			kz=smh*sqrt(EK*vw[2]^2+V0)
			
			SetDataFolder $nf
			if(kdim==0)
				return kx
			elseif(kdim==1)
				return ky
			elseif(kdim==2)
				return kz
			endif
		endif
		SetDataFolder $nf
	else
		print "I can't find your ["+inputwave+"] wave."
	endif
End

Function K_a2k_example(p,t,pn,tn,an,po,to,ao) // k conv 1D wave
	variable p,t,pn,tn,an,po,to,ao
	string iws="input"
	wave iw=$iws
	
	string pws="para"
	Make/O/N=22 $pws
	wave pw=$pws
	
	string ows="output"
	Make/O/N=(DimSize(iw,0),3) $ows
	wave ow=$ows
	
	pw=NaN
	
	pw[3]=p
	pw[4]=t
	
	pw[6]=pn
	pw[7]=tn
	pw[8]=an
	
	pw[9]=po
	pw[10]=to
	pw[11]=ao
	
	
	variable i
	for(i=0;i<DimSize(iw,0);i+=1)
		if(NumType(iw[i][0])==0)
			pw[0]=iw[i][0]
			pw[5]=iw[i][1]
			
			ow[i][0]=K_a2k(pws,0)
			ow[i][1]=K_a2k(pws,1)
			ow[i][2]=K_a2k(pws,2)
		endif
	endfor
End


Macro K_ARPES_start()
	if(strlen(WinList("K_ARPES_p", ";", "WIN:64")))
		DoWindow/F K_ARPES_p
	else
		String nf=GetDataFolder(1)
		NewDataFolder/O root:K_ARPES
		NewDataFolder/O root:K_ARPES:misc
		NewDataFolder/O root:K_ARPES:misc:automap
		NewDataFolder/O root:K_ARPES:misc:rot_matrix
		NewDataFolder/O root:K_ARPES:Curves
		NewDataFolder/S/O root:K_ARPES:global
		make/T/O/N=0 $("win_list")
		SetDataFolder root:K_ARPES:misc
//		variable/g v_mstar=1
		variable/g v_hn
		if(v_hn==0)
			v_hn=50
		endif
		variable/g v_V0
		if(v_V0==0)
			v_V0=12
		endif
		variable/g v_W
		if(v_W==0)
			v_W=4.5
		endif
		variable/g v_EB
		variable/g v_th_s
		variable/g v_th_e
		if(v_th_s==0 && v_th_e==0)
			variable/g v_th_s=-15
			variable/g v_th_e=15
		endif
		variable/g v_po=0
		variable/g v_ph=0
		variable/g v_az=0
		variable/g v_poo=0
		variable/g v_pho=0
		variable/g v_azo=0
		
		variable/g v_po_ofs=0
		variable/g v_ph_ofs=0
		variable/g v_az_ofs=0
		
		variable/g v_bet=0
		variable/g v_alp_emis_set=0
		variable/g v_bet_emis_set=0
		variable/g v_app=0
		variable/g v_map=0
		
		variable/g v_kxmin_w
		variable/g v_kxmax_w
		variable/g v_kymin_w
		variable/g v_kymax_w
		variable/g v_kzmin_w
		variable/g v_kzmax_w
		
		variable/g v_kxwid_w
		variable/g v_kywid_w
		variable/g v_kzwid_w
		
		variable/g v_cross_z
		variable/g v_cross_y
		variable/g v_cross_x
		
		variable/g v_pkf
		variable/g v_zeta
		variable/g v_eta
		
		if(v_cross_x==0 && v_cross_y==0 && v_cross_z==0)
			variable/g v_cross_x=1
		endif
		
		variable/g v_smh=0.512316 // sqrt(2m)/hbar
		variable/g v_th_num=100
		variable/g v_kzfixed=0
		string/g s_winname="window"
		string wnwn=s_winname
		make/O/N=3 $("vec")
		make/O/N=2 $("winposi")
		
		winposi[0]=675.75
		winposi[1]=82.25
		
		variable/g v_revth=1
		variable/g v_revpo=1
		variable/g v_revph=1
		variable/g v_revaz=1
		variable/g v_revbet=1
		
		variable/g v_anal_rot
		
		SetDatafolder root:K_ARPES:misc:automap
		if(!Exists("amap_parameter"))
			K_make_apws()
		endif
		variable/g v_amapval=1
		string/g s_amapval="none"
		
		SetDataFolder root:K_ARPES:global
		variable/g v_kconv=0
		variable/g v_EF=0
		variable/g v_kx_n=100
		variable/g v_ky_n=100
		variable/g v_kconv_vol=0
		variable/g v_kconv_en=0
		variable/g v_outr=-inf
		string/g s_fintime=""
		
		K_make_window()
		K_make_color_list()
		SetDataFolder $nf
		
		K_ARPES_p()
		K_panel_update_cross()
		K_set_status(0,"")
	endif
End

Function K_mani_all_zero()
	string nf=GetDataFolder(1)
	SetDataFolder root:K_ARPES:misc
	variable/g v_po=0
	variable/g v_ph=0
	variable/g v_az=0
	variable/g v_poo=0
	variable/g v_pho=0
	variable/g v_azo=0
	variable/g v_po_ofs=0
	variable/g v_ph_ofs=0
	variable/g v_az_ofs=0
	SetDataFolder $nf
End

Function K_mani_nml2val()
	string nf=GetDataFolder(1)
	SetDataFolder root:K_ARPES:misc
	
	SetDataFolder $nf
End

Function K_make_apws()
	make/O/N=4 $("amap_parameter")
	make/O/N=4 $("amap_parameter_polar")
	wave apw=$("amap_parameter_polar")
	apw[0]=-15
	apw[1]=1
	apw[2]=15
	apw[3]=31
	make/O/N=4 $("amap_parameter_tilt")
	wave apw=$("amap_parameter_tilt")
	apw[0]=-15
	apw[1]=1
	apw[2]=15
	apw[3]=31
	make/O/N=4 $("amap_parameter_azimuth")
	wave apw=$("amap_parameter_azimuth")
	apw[0]=0
	apw[1]=5
	apw[2]=360
	apw[3]=73
	make/O/N=4 $("amap_parameter_beta")
	wave apw=$("amap_parameter_beta")
	apw[0]=-15
	apw[1]=1
	apw[2]=15
	apw[3]=31
	make/O/N=4 $("amap_parameter_hn")
	wave apw=$("amap_parameter_hn")
	apw[0]=400
	apw[1]=10
	apw[2]=600
	apw[3]=21
	make/O/N=4 $("amap_parameter_EB")
	wave apw=$("amap_parameter_EB")
	apw[0]=0
	apw[1]=0.1
	apw[2]=3
	apw[3]=31
End

Window K_ARPES_p() : Panel
	PauseUpdate; Silent 1		// building window...
	NewPanel /K=1 /W=(524,61,823,706) as "K_ARPES panel"
	ModifyPanel cbRGB=(65535,54611,49151), fixedSize=1
	SetDrawLayer UserBack
	DrawText 146,197,"Alpha (deg)"
	DrawText 169,293,"\\f01Energy\\f00 (eV)"
	DrawText 123,62,"Val.@NML"
	DrawText 73,63,"Value"
	DrawText 240,62,"Rev."
	DrawText 20,174,"Emis. Ang."
	DrawText 32,196,"Slit ll"
	DrawText 32,219,"Slit L"
	DrawText 251,198,"Rev."
	DrawText 23,42,"\\f01Manipulator"
	DrawText 143,161,"\\f01Analyzer"
	DrawText 187,62,"Offsets"
	SetVariable setvar_status,pos={224.00,0.00},size={66.00,14.00},bodyWidth=30,title="Status :"
	SetVariable setvar_status,frame=0,fStyle=0,valueBackColor=(65535,54611,49151)
	SetVariable setvar_status,limits={-inf,inf,0},value= _STR:"Idle",noedit= 1
	SetVariable setvar_step,pos={182.00,27.00},size={95.00,14.00},bodyWidth=46,proc=K_SetVarProc_panelstep,title="Panel step"
	SetVariable setvar_step,limits={0,inf,1},value= _NUM:1
	GroupBox Curve,pos={8.00,7.00},size={282.00,389.00},title="Experimental setting"
	SetVariable setvar_po,pos={30.00,67.00},size={86.00,14.00},bodyWidth=60,proc=K_SetVarProc_calc,title="Polar"
	SetVariable setvar_po,valueBackColor=(65535,65534,49151)
	SetVariable setvar_po,limits={-90,90,1},value= root:K_ARPES:misc:v_po
	SetVariable setvar_ph,pos={38.00,88.00},size={78.00,14.00},bodyWidth=60,proc=K_SetVarProc_calc,title="Tilt"
	SetVariable setvar_ph,valueBackColor=(65535,65534,49151)
	SetVariable setvar_ph,limits={-90,90,1},value= root:K_ARPES:misc:v_ph
	SetVariable setvar_az,pos={35.00,109.00},size={81.00,14.00},bodyWidth=60,proc=K_SetVarProc_calc,title="Azi."
	SetVariable setvar_az,valueBackColor=(65535,65534,49151)
	SetVariable setvar_az,value= root:K_ARPES:misc:v_az
	SetVariable setvar_po_nml,pos={128.00,66.00},size={50.00,14.00},bodyWidth=50,proc=K_SetVarProc_calc,title=" "
	SetVariable setvar_po_nml,limits={-90,90,1},value= root:K_ARPES:misc:v_poo
	SetVariable setvar_ph_nml,pos={128.00,88.00},size={50.00,14.00},bodyWidth=50,proc=K_SetVarProc_calc,title=" "
	SetVariable setvar_ph_nml,limits={-90,90,1},value= root:K_ARPES:misc:v_pho
	SetVariable setvar_az_nml,pos={128.00,109.00},size={50.00,14.00},bodyWidth=50,proc=K_SetVarProc_calc,title=" "
	SetVariable setvar_az_nml,value= root:K_ARPES:misc:v_azo
	SetVariable setvar_po_ofs,pos={184.00,66.00},size={50.00,14.00},bodyWidth=50,proc=K_SetVarProc_calc,title=" "
	SetVariable setvar_po_ofs,limits={-90,90,1},value= root:K_ARPES:misc:v_po_ofs
	SetVariable setvar_ph_ofs,pos={184.00,88.00},size={50.00,14.00},bodyWidth=50,proc=K_SetVarProc_calc,title=" "
	SetVariable setvar_ph_ofs,limits={-90,90,1},value= root:K_ARPES:misc:v_ph_ofs
	SetVariable setvar_az_ofs,pos={184.00,109.00},size={50.00,14.00},bodyWidth=50,proc=K_SetVarProc_calc,title=" "
	SetVariable setvar_az_ofs,value= root:K_ARPES:misc:v_az_ofs
	CheckBox check_revpo,pos={244.00,69.00},size={15.00,16.00},proc=K_CheckProc_revpo,title=""
	CheckBox check_revpo,value= 0
	CheckBox check_revph,pos={244.00,91.00},size={15.00,16.00},proc=K_CheckProc_revph,title=""
	CheckBox check_revph,value= 0
	CheckBox check_revaz,pos={244.00,112.00},size={15.00,16.00},proc=K_CheckProc_revaz,title=""
	CheckBox check_revaz,value= 0
	SetVariable setvar_alp,pos={60.00,179.00},size={50.00,18.00},bodyWidth=50,title=" "
	SetVariable setvar_alp,limits={-inf,inf,0},value= root:K_ARPES:misc:v_alp_emis_set
	SetVariable setvar_bet,pos={60.00,201.00},size={50.00,18.00},bodyWidth=50,title=" "
	SetVariable setvar_bet,limits={-inf,inf,0},value= root:K_ARPES:misc:v_bet_emis_set
	Button button_find,pos={79.00,156.00},size={50.00,20.00},proc=K_ButtonProc_sample,title="Set"
	SetVariable setvar_analyzer,pos={143.00,162.00},size={104.00,18.00},bodyWidth=49,proc=K_SetVarProc_ar,title="Rotation :"
	SetVariable setvar_analyzer,limits={-inf,inf,90},value= root:K_ARPES:misc:v_anal_rot
	SetVariable setvar_th_s,pos={161.00,197.00},size={90.00,14.00},bodyWidth=60,proc=K_SetVarProc_thsnum,title="start :"
	SetVariable setvar_th_s,limits={-90,90,1},value= root:K_ARPES:misc:v_th_s
	SetVariable setvar_th_e,pos={166.00,221.00},size={86.00,14.00},bodyWidth=60,proc=K_SetVarProc_thenum,title="end :"
	SetVariable setvar_th_e,limits={-90,90,1},value= root:K_ARPES:misc:v_th_e
	CheckBox check_revth,pos={256.00,210.00},size={15.00,16.00},proc=K_CheckProc_revth,title=""
	CheckBox check_revth,value= 0
	SetVariable setvar_bet_offs,pos={156.00,245.00},size={95.00,14.00},bodyWidth=45,proc=K_SetVarProc_calc,title="Beta (deg):"
	SetVariable setvar_bet_offs,valueBackColor=(65535,65534,49151)
	SetVariable setvar_bet_offs,limits={-90,90,1},value= root:K_ARPES:misc:v_bet
	CheckBox check_revbet,pos={256.00,248.00},size={15.00,16.00},proc=K_CheckProc_revbet,title=""
	CheckBox check_revbet,value= 0
	SetVariable setvar_hn,pos={177.00,295.00},size={81.00,14.00},bodyWidth=60,proc=K_SetVarProc_calc,title="hn :"
	SetVariable setvar_hn,valueBackColor=(65535,65534,49151)
	SetVariable setvar_hn,limits={0,inf,1},value= root:K_ARPES:misc:v_hn
	SetVariable setvar_V0,pos={176.00,315.00},size={81.00,14.00},bodyWidth=60,proc=K_SetVarProc_calc,title="V0 :"
	SetVariable setvar_V0,value= root:K_ARPES:misc:v_V0
	SetVariable setvar_W,pos={180.00,335.00},size={77.00,14.00},bodyWidth=60,proc=K_SetVarProc_calc,title="W :"
	SetVariable setvar_W,value= root:K_ARPES:misc:v_W
	SetVariable setvar_EB,pos={177.00,355.00},size={80.00,14.00},bodyWidth=60,proc=K_SetVarProc_calc,title="EB :"
	SetVariable setvar_EB,value= root:K_ARPES:misc:v_EB
	GroupBox Mapping,pos={17.00,272.00},size={135.00,118.00},title="Mapping"
	PopupMenu popup_amapval,pos={28.00,290.00},size={103.00,23.00},proc=K_PopMenuProc_amap_variable,title="Variable :"
	PopupMenu popup_amapval,mode=1,popvalue="none",value= #"\"none;Polar;Tilt;Azimuth;Beta;hn;EB\""
	SetVariable setvar_map_st,pos={48.00,311.00},size={70.00,14.00},bodyWidth=40,disable=2,proc=K_SetVarProc_amap_parameter,title="Start :"
	SetVariable setvar_map_st,value= root:K_ARPES:misc:automap:amap_parameter[0]
	SetVariable setvar_map_step,pos={48.00,331.00},size={69.00,14.00},bodyWidth=40,disable=2,proc=K_SetVarProc_amap_parameter,title="Step :"
	SetVariable setvar_map_step,value= root:K_ARPES:misc:automap:amap_parameter[1]
	SetVariable setvar_map_en,pos={51.00,351.00},size={66.00,14.00},bodyWidth=40,disable=2,proc=K_SetVarProc_amap_parameter,title="End :"
	SetVariable setvar_map_en,value= root:K_ARPES:misc:automap:amap_parameter[2]
	SetVariable setvar_map_num,pos={60.00,367.00},size={60.00,14.00},bodyWidth=60,disable=2,title=" "
	SetVariable setvar_map_num,format="%d points",frame=0
	SetVariable setvar_map_num,limits={-inf,inf,0},value= root:K_ARPES:misc:automap:amap_parameter[3],noedit= 1
	TabControl tab_mode,pos={29.00,399.00},size={240.00,20.00},proc=K_TabProc_panelmode
	TabControl tab_mode,tabLabel(0)="Simulation mode"
	TabControl tab_mode,tabLabel(1)="Data analyze mode",value= 0
	GroupBox run_setting,pos={9.00,422.00},size={279.00,182.00},title="Plot setting"
	SetVariable setvar_winname,pos={24.00,443.00},size={135.00,14.00},bodyWidth=100,proc=K_SetVarProc_restore,title="Name :"
	SetVariable setvar_winname,value= root:K_ARPES:misc:s_winname
	PopupMenu popup_winlist,pos={159.00,443.00},size={147.00,23.00},proc=K_PopMenuProc_winlist,title="< select from exists"
	PopupMenu popup_winlist,mode=0,value= #"K_winlist_str()"
	SetVariable setvar_show_pr,pos={25.00,465.00},size={100.00,14.00},frame=0
	SetVariable setvar_show_pr,value= _STR:"Plot range (/A)",noedit= 1
	SetVariable setvar_show_min,pos={38.00,479.00},size={40.00,14.00},frame=0
	SetVariable setvar_show_min,value= _STR:"min.",noedit= 1
	SetVariable setvar_kxmin,pos={21.00,496.00},size={58.00,14.00},bodyWidth=40,proc=K_SetVarProc_SetAxix_w,title="kx "
	SetVariable setvar_kxmin,limits={-inf,inf,0.5},value= root:K_ARPES:misc:v_kxmin_w
	SetVariable setvar_kymin,pos={22.00,518.00},size={57.00,14.00},bodyWidth=40,proc=K_SetVarProc_SetAxix_w,title="ky "
	SetVariable setvar_kymin,limits={-inf,inf,0.5},value= root:K_ARPES:misc:v_kymin_w
	SetVariable setvar_kzmin,pos={22.00,539.00},size={57.00,14.00},bodyWidth=40,proc=K_SetVarProc_SetAxix_w,title="kz "
	SetVariable setvar_kzmin,limits={-inf,inf,0.5},value= root:K_ARPES:misc:v_kzmin_w
	SetVariable setvar_show_max,pos={80.00,479.00},size={40.00,14.00},frame=0
	SetVariable setvar_show_max,value= _STR:"max.",noedit= 1
	SetVariable setvar_kxmax,pos={80.00,496.00},size={40.00,14.00},bodyWidth=40,proc=K_SetVarProc_SetAxix_w,title=" "
	SetVariable setvar_kxmax,limits={-inf,inf,0.5},value= root:K_ARPES:misc:v_kxmax_w
	SetVariable setvar_kymax,pos={80.00,518.00},size={40.00,14.00},bodyWidth=40,proc=K_SetVarProc_SetAxix_w,title=" "
	SetVariable setvar_kymax,limits={-inf,inf,0.5},value= root:K_ARPES:misc:v_kymax_w
	SetVariable setvar_kzmax,pos={80.00,539.00},size={40.00,14.00},bodyWidth=40,proc=K_SetVarProc_SetAxix_w,title=" "
	SetVariable setvar_kzmax,limits={-inf,inf,0.5},value= root:K_ARPES:misc:v_kzmax_w
	SetVariable setvar_show_wid,pos={120.00,480.00},size={40.00,14.00},frame=0
	SetVariable setvar_show_wid,value= _STR:"wid.",noedit= 1
	SetVariable setvar_kxwid,pos={122.00,497.00},size={25.00,14.00},bodyWidth=25,proc=K_SetVarProc_SetAxix_w,title=" "
	SetVariable setvar_kxwid,limits={-inf,inf,0},value= root:K_ARPES:misc:v_kxwid_w
	SetVariable setvar_kywid,pos={122.00,519.00},size={25.00,14.00},bodyWidth=25,proc=K_SetVarProc_SetAxix_w,title=" "
	SetVariable setvar_kywid,limits={-inf,inf,0},value= root:K_ARPES:misc:v_kywid_w
	SetVariable setvar_kzwid,pos={122.00,540.00},size={25.00,14.00},bodyWidth=25,proc=K_SetVarProc_SetAxix_w,title=" "
	SetVariable setvar_kzwid,limits={-inf,inf,0},value= root:K_ARPES:misc:v_kzwid_w
	SetVariable setvar_show_fix,pos={149.00,480.00},size={40.00,14.00},frame=0
	SetVariable setvar_show_fix,value= _STR:"fix.",noedit= 1
	CheckBox check_fixkx,pos={152.00,499.00},size={15.00,16.00},proc=K_CheckProc_fixkz,title=""
	CheckBox check_fixkx,variable= root:K_ARPES:misc:v_kxfixed
	CheckBox check_fixky,pos={152.00,521.00},size={15.00,16.00},proc=K_CheckProc_fixkz,title=""
	CheckBox check_fixky,variable= root:K_ARPES:misc:v_kyfixed
	CheckBox check_fixkz,pos={152.00,542.00},size={15.00,16.00},proc=K_CheckProc_fixkz,title=""
	CheckBox check_fixkz,variable= root:K_ARPES:misc:v_kzfixed
	Button button_select_wave,pos={110.00,444.00},size={153.00,20.00},disable=3,proc=K_ButtonProc_select_wave,title=""
	CheckBox check_new,pos={185.00,519.00},size={36.00,16.00},proc=K_CheckProc_app,title="\\f01New"
	CheckBox check_new,variable= root:K_ARPES:misc:v_app
	CheckBox check_map,pos={185.00,538.00},size={44.00,16.00},proc=K_CheckProc_map,title="\\f01Sweep"
	CheckBox check_map,variable= root:K_ARPES:misc:v_map
	Button button_init_window,pos={24.00,571.00},size={100.00,20.00},proc=K_ButtonProc_make_window,title="Init. window"
	SetVariable setvar_show_ARPES_data,pos={35.00,445.00},size={70.00,14.00},disable=3
	SetVariable setvar_show_ARPES_data,frame=0,value= _STR:"ARPES data :",noedit= 1
	SetVariable setvar_kconv_ef,pos={44.00,467.00},size={92.00,14.00},bodyWidth=40,disable=3,proc=K_SetVarProc_kconv_ef,title="Ekin of EF :"
	SetVariable setvar_kconv_ef,limits={-inf,inf,0},value= root:K_ARPES:global:v_EF
	SetVariable setvar_kconv_ek,pos={34.00,487.00},size={101.00,14.00},bodyWidth=40,disable=3,proc=K_SetVarProc_kconv_write_ktw,title="Ekin of data :"
	SetVariable setvar_kconv_ek,limits={-inf,inf,0},value= root:K_ARPES:global:v_ek
	SetVariable setvar_kconv_ef_data,pos={42.00,491.00},size={93.00,14.00},bodyWidth=40,disable=3,proc=K_SetVarProc_kconv_ef_data,title="EF of data :"
	SetVariable setvar_kconv_ef_data,limits={-inf,inf,0},value= root:K_ARPES:global:v_EF_data
	SetVariable setvar_kconv_eofs,pos={27.00,514.00},size={107.00,14.00},bodyWidth=40,disable=3,proc=K_SetVarProc_kconv_eofs,title="Energy offset :"
	SetVariable setvar_kconv_eofs,limits={-inf,inf,0},value= root:K_ARPES:global:v_eofs
	PopupMenu popup_ydim,pos={174.00,469.00},size={89.00,23.00},disable=3,proc=K_PopMenuProc_img,title="Y dim: "
	PopupMenu popup_ydim,mode=1,popvalue="Eng.",value= #"\"Eng.;Map.\""
	CheckBox check_kconv_vol,pos={48.00,541.00},size={62.00,16.00},disable=3,proc=K_CheckProc_kconv_vol,title="Map to 3D"
	CheckBox check_kconv_vol,variable= root:K_ARPES:global:v_kconv_vol
	SetVariable setvar_kconv_en_str,pos={184.00,466.00},size={56.00,14.00},disable=3,title="Energy :"
	SetVariable setvar_kconv_en_str,frame=0
	SetVariable setvar_kconv_en_str,limits={-inf,inf,0},value= _STR:"",noedit= 1
	SetVariable setvar_kconv_en,pos={231.00,467.00},size={40.00,14.00},bodyWidth=40,disable=3,proc=K_SetVarProc_kconv_en_r,title=" "
	SetVariable setvar_kconv_en,valueBackColor=(65535,32768,32768)
	SetVariable setvar_kconv_en,limits={-inf,inf,0},value= root:K_ARPES:global:v_kconv_en
	SetVariable setvar_kconv_en_s,pos={200.00,484.00},size={70.00,14.00},bodyWidth=40,disable=3,proc=K_SetVarProc_kconv_en_r,title="Start :"
	SetVariable setvar_kconv_en_s,limits={-inf,inf,0},value= root:K_ARPES:global:v_e_s
	SetVariable setvar_kconv_en_e,pos={204.00,504.00},size={66.00,14.00},bodyWidth=40,disable=3,proc=K_SetVarProc_kconv_en_r,title="End :"
	SetVariable setvar_kconv_en_e,limits={-inf,inf,0},value= root:K_ARPES:global:v_e_e
	SetVariable setvar_kconv_kxn,pos={179.00,523.00},size={91.00,14.00},bodyWidth=40,disable=3,proc=K_SetVarProc_kconv_write_ktw,title="kx points :"
	SetVariable setvar_kconv_kxn,limits={2,inf,0},value= root:K_ARPES:global:v_kx_n
	SetVariable setvar_kconv_kyn,pos={180.00,541.00},size={90.00,14.00},bodyWidth=40,disable=3,proc=K_SetVarProc_kconv_write_ktw,title="ky points :"
	SetVariable setvar_kconv_kyn,limits={2,inf,0},value= root:K_ARPES:global:v_ky_n
	SetVariable setvar_outr,pos={34.00,561.00},size={137.00,14.00},bodyWidth=40,disable=3,proc=K_SetVarProc_kconv_write_ktw,title="Fill out of range with "
	SetVariable setvar_outr,limits={-inf,inf,0},value= root:K_ARPES:global:v_outr
	Button button_go,pos={176.00,561.00},size={100.00,28.00},proc=K_ButtonProc_go,title="\\Zr200\\f01PLOT"
	Button button_init,pos={76.00,612.00},size={150.00,20.00},proc=K_ButtonProc_init,title="Initialize K_ARPES"
	Button button_init,fColor=(65535,16385,16385)
	Button button_yz,pos={218.00,477.00},size={20.00,20.00},proc=K_ButtonProc_cross,title=" "
	Button button_xy,pos={235.00,477.00},size={20.00,20.00},proc=K_ButtonProc_cross,title=" "
	Button button_xy,fColor=(0,65535,65535)
	Button button_xz,pos={235.00,494.00},size={20.00,20.00},proc=K_ButtonProc_cross,title=" "
	Button button_m_zero,pos={25.00,132.00},size={50.00,20.00},proc=K_ButtonProc_m_zero,title="Zero"
	SetVariable setvar_zeta,pos={29.00,244.00},size={39.00,18.00},bodyWidth=30,proc=K_SetVarProc_calc,title="ζ"
	SetVariable setvar_zeta,limits={-90,90,0},value= root:K_ARPES:misc:v_zeta
	SetVariable setvar_eta,pos={72.00,244.00},size={41.00,18.00},bodyWidth=30,proc=K_SetVarProc_calc,title="η"
	SetVariable setvar_eta,limits={-90,90,0},value= root:K_ARPES:misc:v_eta
	CheckBox check_pkf,pos={23.00,225.00},size={120.00,15.00},proc=K_CheckProc_pkf,title="Photon Momentum"
	CheckBox check_pkf,variable= root:K_ARPES:misc:v_pkf
EndMacro

Function K_panel_update()
	string nf=GetDataFolder(1)
	SetDataFolder root:K_ARPES:misc
	nvar revth=v_revth
	nvar revbet=v_revbet
	nvar revpo=v_revpo
	nvar revph=v_revph
	nvar revaz=v_revaz
	nvar at=v_analtype
	
	if(revth==1)
		CheckBox check_revth win=K_ARPES_p, value=0
	elseif(revth==-1)
		CheckBox check_revth win=K_ARPES_p, value=1
	endif
	if(revbet==1)
		CheckBox check_revbet win=K_ARPES_p, value=0
	elseif(revbet==-1)
		CheckBox check_revbet win=K_ARPES_p, value=1
	endif
	if(revpo==1)
		CheckBox check_revpo win=K_ARPES_p, value=0
	elseif(revpo==-1)
		CheckBox check_revpo win=K_ARPES_p, value=1
	endif
	if(revph==1)
		CheckBox check_revph win=K_ARPES_p, value=0
	elseif(revph==-1)
		CheckBox check_revph win=K_ARPES_p, value=1
	endif
	if(revaz==1)
		CheckBox check_revaz win=K_ARPES_p, value=0
	elseif(revaz==-1)
		CheckBox check_revaz win=K_ARPES_p, value=1
	endif
	if(at==0)
		Slider slider_notation win=K_ARPES_p, value=1
	elseif(at==1)
		Slider slider_notation win=K_ARPES_p, value=0
	endif
	SetDataFolder $nf
End

Function K_refresh_panel()
	DoWindow/F K_ARPES_p
	SetVariable setvar_status,pos={224.00,0.00},size={66.00,14.00},bodyWidth=30,title="Status :"
	SetVariable setvar_status,frame=0,fStyle=0,valueBackColor=(65535,54611,49151)
	SetVariable setvar_status,limits={-inf,inf,0},value= _STR:"Idle",noedit= 1
	SetVariable setvar_step,pos={182.00,27.00},size={95.00,14.00},bodyWidth=46,proc=K_SetVarProc_panelstep,title="Panel step"
	SetVariable setvar_step,limits={0,inf,1},value= _NUM:1
	GroupBox Curve,pos={8.00,7.00},size={282.00,389.00},title="Experimental setting"
	SetVariable setvar_po,pos={30.00,67.00},size={86.00,14.00},bodyWidth=60,proc=K_SetVarProc_calc,title="Polar"
	SetVariable setvar_po,valueBackColor=(65535,65534,49151)
	SetVariable setvar_po,limits={-90,90,1},value= root:K_ARPES:misc:v_po
	SetVariable setvar_ph,pos={38.00,88.00},size={78.00,14.00},bodyWidth=60,proc=K_SetVarProc_calc,title="Tilt"
	SetVariable setvar_ph,valueBackColor=(65535,65534,49151)
	SetVariable setvar_ph,limits={-90,90,1},value= root:K_ARPES:misc:v_ph
	SetVariable setvar_az,pos={35.00,109.00},size={81.00,14.00},bodyWidth=60,proc=K_SetVarProc_calc,title="Azi."
	SetVariable setvar_az,valueBackColor=(65535,65534,49151)
	SetVariable setvar_az,value= root:K_ARPES:misc:v_az
	SetVariable setvar_po_nml,pos={128.00,66.00},size={50.00,14.00},bodyWidth=50,proc=K_SetVarProc_calc,title=" "
	SetVariable setvar_po_nml,limits={-90,90,1},value= root:K_ARPES:misc:v_poo
	SetVariable setvar_ph_nml,pos={128.00,88.00},size={50.00,14.00},bodyWidth=50,proc=K_SetVarProc_calc,title=" "
	SetVariable setvar_ph_nml,limits={-90,90,1},value= root:K_ARPES:misc:v_pho
	SetVariable setvar_az_nml,pos={128.00,109.00},size={50.00,14.00},bodyWidth=50,proc=K_SetVarProc_calc,title=" "
	SetVariable setvar_az_nml,value= root:K_ARPES:misc:v_azo
	SetVariable setvar_po_ofs,pos={184.00,66.00},size={50.00,14.00},bodyWidth=50,proc=K_SetVarProc_calc,title=" "
	SetVariable setvar_po_ofs,limits={-90,90,1},value= root:K_ARPES:misc:v_po_ofs
	SetVariable setvar_ph_ofs,pos={184.00,88.00},size={50.00,14.00},bodyWidth=50,proc=K_SetVarProc_calc,title=" "
	SetVariable setvar_ph_ofs,limits={-90,90,1},value= root:K_ARPES:misc:v_ph_ofs
	SetVariable setvar_az_ofs,pos={184.00,109.00},size={50.00,14.00},bodyWidth=50,proc=K_SetVarProc_calc,title=" "
	SetVariable setvar_az_ofs,value= root:K_ARPES:misc:v_az_ofs
	CheckBox check_revpo,pos={244.00,69.00},size={15.00,16.00},proc=K_CheckProc_revpo,title=""
	CheckBox check_revpo,value= 0
	CheckBox check_revph,pos={244.00,91.00},size={15.00,16.00},proc=K_CheckProc_revph,title=""
	CheckBox check_revph,value= 0
	CheckBox check_revaz,pos={244.00,112.00},size={15.00,16.00},proc=K_CheckProc_revaz,title=""
	CheckBox check_revaz,value= 0
	SetVariable setvar_alp,pos={60.00,179.00},size={50.00,18.00},bodyWidth=50,title=" "
	SetVariable setvar_alp,limits={-inf,inf,0},value= root:K_ARPES:misc:v_alp_emis_set
	SetVariable setvar_bet,pos={60.00,201.00},size={50.00,18.00},bodyWidth=50,title=" "
	SetVariable setvar_bet,limits={-inf,inf,0},value= root:K_ARPES:misc:v_bet_emis_set
	Button button_find,pos={79.00,156.00},size={50.00,20.00},proc=K_ButtonProc_sample,title="Set"
	SetVariable setvar_analyzer,pos={143.00,162.00},size={104.00,18.00},bodyWidth=49,proc=K_SetVarProc_ar,title="Rotation :"
	SetVariable setvar_analyzer,limits={-inf,inf,90},value= root:K_ARPES:misc:v_anal_rot
	SetVariable setvar_th_s,pos={161.00,197.00},size={90.00,14.00},bodyWidth=60,proc=K_SetVarProc_thsnum,title="start :"
	SetVariable setvar_th_s,limits={-90,90,1},value= root:K_ARPES:misc:v_th_s
	SetVariable setvar_th_e,pos={166.00,221.00},size={86.00,14.00},bodyWidth=60,proc=K_SetVarProc_thenum,title="end :"
	SetVariable setvar_th_e,limits={-90,90,1},value= root:K_ARPES:misc:v_th_e
	CheckBox check_revth,pos={256.00,210.00},size={15.00,16.00},proc=K_CheckProc_revth,title=""
	CheckBox check_revth,value= 0
	SetVariable setvar_bet_offs,pos={156.00,245.00},size={95.00,14.00},bodyWidth=45,proc=K_SetVarProc_calc,title="Beta (deg):"
	SetVariable setvar_bet_offs,valueBackColor=(65535,65534,49151)
	SetVariable setvar_bet_offs,limits={-90,90,1},value= root:K_ARPES:misc:v_bet
	CheckBox check_revbet,pos={256.00,248.00},size={15.00,16.00},proc=K_CheckProc_revbet,title=""
	CheckBox check_revbet,value= 0
	SetVariable setvar_hn,pos={177.00,295.00},size={81.00,14.00},bodyWidth=60,proc=K_SetVarProc_calc,title="hn :"
	SetVariable setvar_hn,valueBackColor=(65535,65534,49151)
	SetVariable setvar_hn,limits={0,inf,1},value= root:K_ARPES:misc:v_hn
	SetVariable setvar_V0,pos={176.00,315.00},size={81.00,14.00},bodyWidth=60,proc=K_SetVarProc_calc,title="V0 :"
	SetVariable setvar_V0,value= root:K_ARPES:misc:v_V0
	SetVariable setvar_W,pos={180.00,335.00},size={77.00,14.00},bodyWidth=60,proc=K_SetVarProc_calc,title="W :"
	SetVariable setvar_W,value= root:K_ARPES:misc:v_W
	SetVariable setvar_EB,pos={177.00,355.00},size={80.00,14.00},bodyWidth=60,proc=K_SetVarProc_calc,title="EB :"
	SetVariable setvar_EB,value= root:K_ARPES:misc:v_EB
	GroupBox Mapping,pos={17.00,272.00},size={135.00,118.00},title="Mapping"
	PopupMenu popup_amapval,pos={28.00,290.00},size={103.00,23.00},proc=K_PopMenuProc_amap_variable,title="Variable :"
	PopupMenu popup_amapval,mode=1,popvalue="none",value= #"\"none;Polar;Tilt;Azimuth;Beta;hn;EB\""
	SetVariable setvar_map_st,pos={48.00,311.00},size={70.00,14.00},bodyWidth=40,disable=2,proc=K_SetVarProc_amap_parameter,title="Start :"
	SetVariable setvar_map_st,value= root:K_ARPES:misc:automap:amap_parameter[0]
	SetVariable setvar_map_step,pos={48.00,331.00},size={69.00,14.00},bodyWidth=40,disable=2,proc=K_SetVarProc_amap_parameter,title="Step :"
	SetVariable setvar_map_step,value= root:K_ARPES:misc:automap:amap_parameter[1]
	SetVariable setvar_map_en,pos={51.00,351.00},size={66.00,14.00},bodyWidth=40,disable=2,proc=K_SetVarProc_amap_parameter,title="End :"
	SetVariable setvar_map_en,value= root:K_ARPES:misc:automap:amap_parameter[2]
	SetVariable setvar_map_num,pos={60.00,367.00},size={60.00,14.00},bodyWidth=60,disable=2,title=" "
	SetVariable setvar_map_num,format="%d points",frame=0
	SetVariable setvar_map_num,limits={-inf,inf,0},value= root:K_ARPES:misc:automap:amap_parameter[3],noedit= 1
	TabControl tab_mode,pos={29.00,399.00},size={240.00,20.00},proc=K_TabProc_panelmode
	TabControl tab_mode,tabLabel(0)="Simulation mode"
	TabControl tab_mode,tabLabel(1)="Data analyze mode",value= 0
	GroupBox run_setting,pos={9.00,422.00},size={279.00,182.00},title="Plot setting"
	SetVariable setvar_winname,pos={24.00,443.00},size={135.00,14.00},bodyWidth=100,proc=K_SetVarProc_restore,title="Name :"
	SetVariable setvar_winname,value= root:K_ARPES:misc:s_winname
	PopupMenu popup_winlist,pos={159.00,443.00},size={147.00,23.00},proc=K_PopMenuProc_winlist,title="< select from exists"
	PopupMenu popup_winlist,mode=0,value= #"K_winlist_str()"
	SetVariable setvar_show_pr,pos={25.00,465.00},size={100.00,14.00},frame=0
	SetVariable setvar_show_pr,value= _STR:"Plot range (/A)",noedit= 1
	SetVariable setvar_show_min,pos={38.00,479.00},size={40.00,14.00},frame=0
	SetVariable setvar_show_min,value= _STR:"min.",noedit= 1
	SetVariable setvar_kxmin,pos={21.00,496.00},size={58.00,14.00},bodyWidth=40,proc=K_SetVarProc_SetAxix_w,title="kx "
	SetVariable setvar_kxmin,limits={-inf,inf,0.5},value= root:K_ARPES:misc:v_kxmin_w
	SetVariable setvar_kymin,pos={22.00,518.00},size={57.00,14.00},bodyWidth=40,proc=K_SetVarProc_SetAxix_w,title="ky "
	SetVariable setvar_kymin,limits={-inf,inf,0.5},value= root:K_ARPES:misc:v_kymin_w
	SetVariable setvar_kzmin,pos={22.00,539.00},size={57.00,14.00},bodyWidth=40,proc=K_SetVarProc_SetAxix_w,title="kz "
	SetVariable setvar_kzmin,limits={-inf,inf,0.5},value= root:K_ARPES:misc:v_kzmin_w
	SetVariable setvar_show_max,pos={80.00,479.00},size={40.00,14.00},frame=0
	SetVariable setvar_show_max,value= _STR:"max.",noedit= 1
	SetVariable setvar_kxmax,pos={80.00,496.00},size={40.00,14.00},bodyWidth=40,proc=K_SetVarProc_SetAxix_w,title=" "
	SetVariable setvar_kxmax,limits={-inf,inf,0.5},value= root:K_ARPES:misc:v_kxmax_w
	SetVariable setvar_kymax,pos={80.00,518.00},size={40.00,14.00},bodyWidth=40,proc=K_SetVarProc_SetAxix_w,title=" "
	SetVariable setvar_kymax,limits={-inf,inf,0.5},value= root:K_ARPES:misc:v_kymax_w
	SetVariable setvar_kzmax,pos={80.00,539.00},size={40.00,14.00},bodyWidth=40,proc=K_SetVarProc_SetAxix_w,title=" "
	SetVariable setvar_kzmax,limits={-inf,inf,0.5},value= root:K_ARPES:misc:v_kzmax_w
	SetVariable setvar_show_wid,pos={120.00,480.00},size={40.00,14.00},frame=0
	SetVariable setvar_show_wid,value= _STR:"wid.",noedit= 1
	SetVariable setvar_kxwid,pos={122.00,497.00},size={25.00,14.00},bodyWidth=25,proc=K_SetVarProc_SetAxix_w,title=" "
	SetVariable setvar_kxwid,limits={-inf,inf,0},value= root:K_ARPES:misc:v_kxwid_w
	SetVariable setvar_kywid,pos={122.00,519.00},size={25.00,14.00},bodyWidth=25,proc=K_SetVarProc_SetAxix_w,title=" "
	SetVariable setvar_kywid,limits={-inf,inf,0},value= root:K_ARPES:misc:v_kywid_w
	SetVariable setvar_kzwid,pos={122.00,540.00},size={25.00,14.00},bodyWidth=25,proc=K_SetVarProc_SetAxix_w,title=" "
	SetVariable setvar_kzwid,limits={-inf,inf,0},value= root:K_ARPES:misc:v_kzwid_w
	SetVariable setvar_show_fix,pos={149.00,480.00},size={40.00,14.00},frame=0
	SetVariable setvar_show_fix,value= _STR:"fix.",noedit= 1
	CheckBox check_fixkx,pos={152.00,499.00},size={15.00,16.00},proc=K_CheckProc_fixkz,title=""
	CheckBox check_fixkx,variable= root:K_ARPES:misc:v_kxfixed
	CheckBox check_fixky,pos={152.00,521.00},size={15.00,16.00},proc=K_CheckProc_fixkz,title=""
	CheckBox check_fixky,variable= root:K_ARPES:misc:v_kyfixed
	CheckBox check_fixkz,pos={152.00,542.00},size={15.00,16.00},proc=K_CheckProc_fixkz,title=""
	CheckBox check_fixkz,variable= root:K_ARPES:misc:v_kzfixed
	Button button_select_wave,pos={110.00,444.00},size={153.00,20.00},disable=3,proc=K_ButtonProc_select_wave,title=""
	CheckBox check_new,pos={185.00,519.00},size={36.00,16.00},proc=K_CheckProc_app,title="\\f01New"
	CheckBox check_new,variable= root:K_ARPES:misc:v_app
	CheckBox check_map,pos={185.00,538.00},size={44.00,16.00},proc=K_CheckProc_map,title="\\f01Sweep"
	CheckBox check_map,variable= root:K_ARPES:misc:v_map
	Button button_init_window,pos={24.00,571.00},size={100.00,20.00},proc=K_ButtonProc_make_window,title="Init. window"
	SetVariable setvar_show_ARPES_data,pos={35.00,445.00},size={70.00,14.00},disable=3
	SetVariable setvar_show_ARPES_data,frame=0,value= _STR:"ARPES data :",noedit= 1
	SetVariable setvar_kconv_ef,pos={44.00,467.00},size={92.00,14.00},bodyWidth=40,disable=3,proc=K_SetVarProc_kconv_ef,title="Ekin of EF :"
	SetVariable setvar_kconv_ef,limits={-inf,inf,0},value= root:K_ARPES:global:v_EF
	SetVariable setvar_kconv_ek,pos={34.00,487.00},size={101.00,14.00},bodyWidth=40,disable=3,proc=K_SetVarProc_kconv_write_ktw,title="Ekin of data :"
	SetVariable setvar_kconv_ek,limits={-inf,inf,0},value= root:K_ARPES:global:v_ek
	SetVariable setvar_kconv_ef_data,pos={42.00,491.00},size={93.00,14.00},bodyWidth=40,disable=3,proc=K_SetVarProc_kconv_ef_data,title="EF of data :"
	SetVariable setvar_kconv_ef_data,limits={-inf,inf,0},value= root:K_ARPES:global:v_EF_data
	SetVariable setvar_kconv_eofs,pos={27.00,514.00},size={107.00,14.00},bodyWidth=40,disable=3,proc=K_SetVarProc_kconv_eofs,title="Energy offset :"
	SetVariable setvar_kconv_eofs,limits={-inf,inf,0},value= root:K_ARPES:global:v_eofs
	PopupMenu popup_ydim,pos={174.00,469.00},size={89.00,23.00},disable=3,proc=K_PopMenuProc_img,title="Y dim: "
	PopupMenu popup_ydim,mode=1,popvalue="Eng.",value= #"\"Eng.;Map.\""
	CheckBox check_kconv_vol,pos={48.00,541.00},size={62.00,16.00},disable=3,proc=K_CheckProc_kconv_vol,title="Map to 3D"
	CheckBox check_kconv_vol,variable= root:K_ARPES:global:v_kconv_vol
	SetVariable setvar_kconv_en_str,pos={184.00,466.00},size={56.00,14.00},disable=3,title="Energy :"
	SetVariable setvar_kconv_en_str,frame=0
	SetVariable setvar_kconv_en_str,limits={-inf,inf,0},value= _STR:"",noedit= 1
	SetVariable setvar_kconv_en,pos={231.00,467.00},size={40.00,14.00},bodyWidth=40,disable=3,proc=K_SetVarProc_kconv_en_r,title=" "
	SetVariable setvar_kconv_en,valueBackColor=(65535,32768,32768)
	SetVariable setvar_kconv_en,limits={-inf,inf,0},value= root:K_ARPES:global:v_kconv_en
	SetVariable setvar_kconv_en_s,pos={200.00,484.00},size={70.00,14.00},bodyWidth=40,disable=3,proc=K_SetVarProc_kconv_en_r,title="Start :"
	SetVariable setvar_kconv_en_s,limits={-inf,inf,0},value= root:K_ARPES:global:v_e_s
	SetVariable setvar_kconv_en_e,pos={204.00,504.00},size={66.00,14.00},bodyWidth=40,disable=3,proc=K_SetVarProc_kconv_en_r,title="End :"
	SetVariable setvar_kconv_en_e,limits={-inf,inf,0},value= root:K_ARPES:global:v_e_e
	SetVariable setvar_kconv_kxn,pos={179.00,523.00},size={91.00,14.00},bodyWidth=40,disable=3,proc=K_SetVarProc_kconv_write_ktw,title="kx points :"
	SetVariable setvar_kconv_kxn,limits={2,inf,0},value= root:K_ARPES:global:v_kx_n
	SetVariable setvar_kconv_kyn,pos={180.00,541.00},size={90.00,14.00},bodyWidth=40,disable=3,proc=K_SetVarProc_kconv_write_ktw,title="ky points :"
	SetVariable setvar_kconv_kyn,limits={2,inf,0},value= root:K_ARPES:global:v_ky_n
	SetVariable setvar_outr,pos={34.00,561.00},size={137.00,14.00},bodyWidth=40,disable=3,proc=K_SetVarProc_kconv_write_ktw,title="Fill out of range with "
	SetVariable setvar_outr,limits={-inf,inf,0},value= root:K_ARPES:global:v_outr
	Button button_go,pos={176.00,561.00},size={100.00,28.00},proc=K_ButtonProc_go,title="\\Zr200\\f01PLOT"
	Button button_init,pos={76.00,612.00},size={150.00,20.00},proc=K_ButtonProc_init,title="Initialize K_ARPES"
	Button button_init,fColor=(65535,16385,16385)
	Button button_yz,pos={218.00,477.00},size={20.00,20.00},proc=K_ButtonProc_cross,title=" "
	Button button_xy,pos={235.00,477.00},size={20.00,20.00},proc=K_ButtonProc_cross,title=" "
	Button button_xy,fColor=(0,65535,65535)
	Button button_xz,pos={235.00,494.00},size={20.00,20.00},proc=K_ButtonProc_cross,title=" "
	Button button_m_zero,pos={25.00,132.00},size={50.00,20.00},proc=K_ButtonProc_m_zero,title="Zero"
	SetVariable setvar_zeta,pos={29.00,244.00},size={39.00,18.00},bodyWidth=30,proc=K_SetVarProc_calc,title="ζ"
	SetVariable setvar_zeta,limits={-90,90,0},value= root:K_ARPES:misc:v_zeta
	SetVariable setvar_eta,pos={72.00,244.00},size={41.00,18.00},bodyWidth=30,proc=K_SetVarProc_calc,title="η"
	SetVariable setvar_eta,limits={-90,90,0},value= root:K_ARPES:misc:v_eta
	CheckBox check_pkf,pos={23.00,225.00},size={120.00,15.00},proc=K_CheckProc_pkf,title="Photon Momentum"
	CheckBox check_pkf,variable= root:K_ARPES:misc:v_pkf
	
	string nf=GetDataFolder(1)
	SetDataFolder root:K_ARPES:misc
	variable/g v_analtype=0
	K_panel_update_cross()
	SetDataFolder root:K_ARPES:global
	variable/g v_kconv
	variable kc=v_kconv
	SetDataFolder root:K_ARPES:misc:automap
	variable/g v_amapval
	K_set_amapval(v_amapval)
	if(kc==1)
		K_show_targwave()
	endif
	K_set_status(0,"")
	SetDataFolder $nf
End

Function K_ARPES_calc_kline()
	K_find_curve()
	K_calc_normal()
	K_ARPES_calc_kline_calc()
	K_calcAxisRange_w()
	K_save_misc()
End

Function K_cord_k2a(kx,ky,kz,adim)
	variable kx,ky,kz,adim
	// ang = 0:alpha, 1:beta
	// type = 0:Polar-angular notation 1:Tilt-angular notation
	
	variable alp_r,bet_r
	
	string nf=GetDataFolder(1)
	SetDataFolder root:K_ARPES:misc
	variable/g v_analtype
	variable type=v_analtype
	
	SetDataFolder root:K_ARPES:misc:rot_matrix
	make/O/N=3 $("temp")
	wave tw=$("temp")
	tw[0]=kx
	tw[1]=ky
	tw[2]=kz
	
	wave arw=$("azofs_rm")
	
	MatrixMultiply arw,tw
	
	wave nw=$("M_product")
	
	kx=nw[0]
	ky=nw[1]
	kz=nw[2]
	
	KillWaves tw,nw
	
	SetDataFolder $nf
	
	
	if(type==0)
		if(kx==0 && ky==0)
			alp_r=0
			bet_r=0
		else
			alp_r=(kx/sqrt(kx^2+ky^2))*acos(kz/sqrt(kx^2+ky^2+kz^2))
			bet_r=(ky/sqrt(kx^2+ky^2))*acos(kz/sqrt(kx^2+ky^2+kz^2))
		endif
	elseif(type==1)
		if(kx==0 && kz==0)
			alp_r=0
			bet_r=0
		else
			alp_r=(acos(kx/sqrt(kx^2+ky^2+kz^2))-pi/2)*-1
			bet_r=sign(kz)*acos(-ky/sqrt(ky^2+kz^2))-pi/2
		endif
	endif
	
	if(adim==0)
		return alp_r
	elseif(adim==1)
		return bet_r
	endif
	
End

Function K_calc_normal()
	variable kx,ky,kz
	
	String nf=GetDataFolder(1)
	SetDataFolder root:K_ARPES:misc
	
	variable alp_r,bet_r
	
	K_make_manip_rot_matrix()
	SetDataFolder root:K_ARPES:misc:rot_matrix
	wave aw=$("NML")
	
	kx=aw[0]
	ky=aw[1]
	kz=aw[2]
	
	alp_r=K_cord_k2a(kx,ky,kz,0)
	bet_r=K_cord_k2a(kx,ky,kz,1)
	
	SetDataFolder root:K_ARPES:misc
	
	variable/g v_alp_N=alp_r*180/pi
	variable/g v_bet_N=bet_r*180/pi
	
	SetDataFolder $nf
End

Function K_calc_emis(alp)
	variable alp
	variable vx,vy,vz
	variable vx1,vy1,vz1
	string nf=GetDataFolder(1)
	SetDataFolder root:K_ARPES:misc
	string vws="vec"
	string pkws="photon_k"
	make/O/N=3 $vws
	make/O/N=3 $pkws
	wave pkw=$pkws
	variable alp_r=alp*pi/180
	variable/g v_bet
	variable/g v_revbet
	variable bet_r=v_revbet*v_bet*pi/180
	variable/g v_az
	variable/g v_azo
	variable rot_r=(v_az-v_azo)*pi/180
	nvar v_kconv=root:K_ARPES:global:v_kconv
	variable/g v_analtype
	variable type=v_analtype
	
	if(type==0)
		if(alp_r==0 && bet_r==0)
			vx1=0
			vy1=0
			vz1=1
		else
			vx1=alp_r/sqrt(alp_r^2+bet_r^2)*sin(sqrt(alp_r^2+bet_r^2))
			vy1=bet_r/sqrt(alp_r^2+bet_r^2)*sin(sqrt(alp_r^2+bet_r^2))
			vz1=cos(sqrt(alp_r^2+bet_r^2))
		endif
	elseif(type==1)
		vx1=sin(alp_r)
		vy1=cos(alp_r)*sin(bet_r)
		vz1=cos(alp_r)*cos(bet_r)
	endif
	
	variable/g v_anal_rot
	variable ar_r=v_anal_rot*pi/180
	
	nvar hn=v_hn
	variable pkx,pky,pkz
	nvar W=v_W
	nvar smh=v_smh
	nvar pkf=v_pkf
	nvar vzeta=v_zeta
	nvar eta=v_eta
	
	variable pk=K_ev2ai(hn)
	
	if(pkf)
		pkw[0]=pk*sin(vzeta/180*pi)
		pkw[1]=pk*cos(vzeta/180*pi)*sin(eta/180*pi)
		pkw[2]=pk*cos(vzeta/180*pi)*cos(eta/180*pi)
	else
		pkw[0]=0
		pkw[1]=0
		pkw[2]=0
	endif
	
	SetDataFolder root:K_ARPES:misc:rot_matrix
	wave aw=$("NML")
	SetDataFolder root:K_ARPES:misc
	wave vw=$vws
	vw[0]=vx1-pkx
	vw[1]=vy1-pky
	vw[2]=vz1-pkz
	wave vw=$vws
	SetDataFolder root:K_ARPES:misc:rot_matrix
	MatrixMultiply $("man_rot_inv"),vw
	wave mpw=$("M_product")
	vw=mpw
	MatrixMultiply $("man_rot_inv"),pkw
	wave mpw=$("M_product")
	pkw=mpw
	
	Killwaves mpw
	
	vx=vw[0]
	vy=vw[1]
	vz=vw[2]
	
//	vx=vx1*cos(ar_r)-vy1*sin(ar_r)
//	vy=vx1*sin(ar_r)+vy1*cos(ar_r)
//	vz=vz1
	
//	vw[0]=vx
//	vw[1]=vy
//	vw[2]=vz
	
	
	variable alp_e_r,bet_e_r
	
	if(v_kconv==0)
		alp_e_r=K_cord_k2a(vx,vy,vz,0)
		bet_e_r=K_cord_k2a(vx,vy,vz,1)
		
		SetDataFolder root:K_ARPES:misc
		
		variable/g v_alp_emis=alp_e_r*180/pi
		variable/g v_bet_emis=bet_e_r*180/pi
	endif
	SetDataFolder $nf
End

Function K_make_manip_rot_matrix()
	string nf=GetDataFolder(1)
	
	K_make_Rar()
	K_make_Rpoo()
	K_make_Rtio()
	K_make_Raz()
	K_make_Rti()
	K_make_Rpo()
	K_make_Razofs()
	
	SetDataFolder root:K_ARPES:misc:rot_matrix
	wave pom=$("po_rm")
	wave tim=$("ti_rm")
	wave azm=$("az_rm")
	wave tiom=$("tio_rm")
	wave poom=$("poo_rm")
	wave azom=$("azofs_rm")
	string nws="NML"
	Make/O/N=3 $nws
	wave nw=$nws
	
	MatrixMultiply pom,tim,azm,tiom,poom,azom
	Duplicate/O $("M_product"),man_rot
	Killwaves $("M_product")
	
	MatrixInverse $("man_rot")
	Duplicate/O $("M_Inverse"),man_rot_inv
	Killwaves $("M_Inverse")
	
	nw=0
	nw[2]=1
	
	wave nw=$nws
	MatrixMultiply $("man_rot_inv"),nw
	Duplicate/O $("M_product"),NML
	Killwaves $("M_product")
	
	K_sep_Rmat()
	
	SetDataFolder $nf
End

Function K_sep_Rmat()
	string nf=GetDataFolder(1)
	SetDataFolder root:K_ARPES:misc:automap
	variable/g v_amapval
	variable amap=v_amapval
	SetDataFolder root:K_ARPES:misc:rot_matrix
	wave arm=$("ar_rm")
	wave arm_i=$("ar_rm_inv")
	wave pom=$("po_rm")
	wave tim=$("ti_rm")
	wave azm=$("az_rm")
	wave tiom=$("tio_rm")
	wave poom=$("poo_rm")
	
	if(amap==2)
	
		MatrixMultiply arm_i,tim,azm,tiom,poom
		Duplicate/O $("M_product"),R1
		KillWaves $("M_product")
		MatrixInverse $("R1")
		Duplicate/O $("M_inverse"),R1_inv
		KillWaves $("M_inverse")
		
		Duplicate/O arm,R2
		MatrixInverse $("R2")
		Duplicate/O $("M_inverse"),R2_inv
		KillWaves $("M_inverse")
	
	elseif(amap==3)
	
		MatrixMultiply arm_i,azm,tiom,poom
		Duplicate/O $("M_product"),R1
		KillWaves $("M_product")
		MatrixInverse $("R1")
		Duplicate/O $("M_inverse"),R1_inv
		KillWaves $("M_inverse")
		
		MatrixMultiply pom,arm
		Duplicate/O $("M_product"),R2
		KillWaves $("M_product")
		MatrixInverse $("R2")
		Duplicate/O $("M_inverse"),R2_inv
		KillWaves $("M_inverse")
		
	elseif(amap==4)
	
		MatrixMultiply arm_i,tiom,poom
		Duplicate/O $("M_product"),R1
		KillWaves $("M_product")
		MatrixInverse $("R1")
		Duplicate/O $("M_inverse"),R1_inv
		KillWaves $("M_inverse")
		
		MatrixMultiply pom,tim,arm
		Duplicate/O $("M_product"),R2
		KillWaves $("M_product")
		MatrixInverse $("R2")
		Duplicate/O $("M_inverse"),R2_inv
		KillWaves $("M_inverse")
	
	endif
	
	
	SetDataFolder $nf
End

Function K_ARPES_calc_kline_calc()
	variable EK,th,th_step,kx,ky,kz,kx1,ky1,kz1
	
	String nf=GetDataFolder(1)
	SetDataFolder root:K_ARPES:global
	variable/g v_kconv
	variable kconv=v_kconv
	SetDataFolder root:K_ARPES:misc
	Variable/g v_th_num
	variable th_num=v_th_num
	Variable/g v_smh
	variable smh=v_smh
	variable/g v_hn
	variable hn=v_hn
	variable/g v_V0
	variable V0=v_V0
	variable/g v_W
	variable W=v_W
	variable/g v_EB
	variable EB=v_EB
	variable/g v_th_s
	variable th_s=v_th_s
	variable/g v_th_e
	variable th_e=v_th_e
	string/g s_thsnum
	string thsnum=s_thsnum
	string/g s_thenum
	string thenum=s_thenum
	variable/g v_alp_N
	variable alp_N=v_alp_N
	variable/g v_bet_N
	variable bet_N=v_bet_N
	variable/g v_revth
	variable revth=v_revth
	
	string/g s_winname
	string winna=s_winname
	
	variable/g v_curve_num
	variable cn=v_curve_num
	
	if(th_e==th_s)
		th_step=1
		th_num=0
	else
		th_step=(th_e-th_s)/th_num
	endif
	
	EK=hn-W-EB
	
	variable pk=smh*hn
	
	
	string cxws="C"+num2str(cn)+"kx"
	string cyws="C"+num2str(cn)+"ky"
	string czws="C"+num2str(cn)+"kz"
	
	wave winw=$winna
	
	SetDataFolder root:K_ARPES:Curves
	SetDataFolder $winna
	
	Make/O/N=(th_num+1) $cxws
	Make/O/N=(th_num+1) $cyws
	Make/O/N=(th_num+1) $czws
	
	wave cxw=$cxws
	wave cyw=$cyws
	wave czw=$czws
	
	note/K cxw,K_make_note()
	note/K cyw,K_make_note()
	note/K czw,K_make_note()
	
	variable counter=0
	
	SetDataFolder root:K_ARPES:misc
	
	wave vw=$("vec")
	
	for(th=th_s;th<th_e+th_step/2;th+=th_step)
		variable thr=revth*th
		
		K_calc_emis(thr)
		wave pkw=$("photon_k")
		
		kx=smh*sqrt(EK)*vw[0]-pkw[0]
		ky=smh*sqrt(EK)*vw[1]-pkw[1]
		kz=smh*sqrt(EK*vw[2]^2+V0)-pkw[2]
		
		cxw[counter]=kx
		cyw[counter]=ky
		czw[counter]=kz
		
		counter=counter+1
	endfor
	
	K_calc_emis((th_s+th_e)/2)
	
	variable/g v_alp_emis
	variable/g v_bet_emis
	variable/g v_alp_emis_set=round(v_alp_emis*1000)/1000
	variable/g v_bet_emis_set=round(v_bet_emis*1000)/1000
	
	SetDataFolder $nf
End

Function K_set_emis_ang_window()
	string nf=GetDataFolder(1)
	SetDataFolder root:K_ARPES:misc
	variable/g v_th_s
	variable/g v_th_e
	K_calc_normal()
	K_calc_emis((v_th_s+v_th_e)/2)
	variable/g v_alp_emis
	variable/g v_bet_emis
	variable/g v_alp_emis_set=round(v_alp_emis*1000)/1000
	variable/g v_bet_emis_set=round(v_bet_emis*1000)/1000
	SetDataFolder $nf
End


Function/S K_make_note()
	String nf=GetDataFolder(1)
	SetDataFolder root:K_ARPES:misc
	Variable/g v_th_num
	variable th_num=v_th_num
	variable/g v_hn
	variable hn=v_hn
	variable/g v_V0
	variable V0=v_V0
	variable/g v_W
	variable W=v_W
	variable/g v_EB
	variable EB=v_EB
	variable/g v_th_s
	variable th_s=v_th_s
	variable/g v_th_e
	variable th_e=v_th_e
	variable/g v_bet
	variable bet=v_bet
	variable/g v_po
	variable po=v_po
	variable/g v_ph
	variable ph=v_ph
	variable/g v_az
	variable az=v_az
	variable/g v_poo
	variable poo=v_poo
	variable/g v_pho
	variable pho=v_pho
	variable/g v_azo
	variable azo=v_azo
	variable/g v_po_ofs
	variable poof=v_po_ofs
	variable/g v_ph_ofs
	variable phof=v_ph_ofs
	variable/g v_az_ofs
	variable azof=v_az_ofs
	variable/g v_revpo
	variable revpo=v_revpo
	variable/g v_revph
	variable revph=v_revph
	variable/g v_revaz
	variable revaz=v_revaz
	variable/g v_revth
	variable revth=v_revth
	variable/g v_analtype
	variable type=v_analtype
	variable/g v_anal_rot
	variable ar=v_anal_rot
	
	
	string notes="\nMANIPULATOR\n\n"
	notes+="Value\n"
	notes+="Polor = "+num2str(po)+" (deg)\n"
	notes+="Tilt = "+num2str(ph)+" (deg)\n"
	notes+="Azimuth = "+num2str(az)+" (deg)\n"
	notes+="\n"
	notes+="Val.@NML\n"
	notes+="Polor = "+num2str(poo)+" (deg)\n"
	notes+="Tilt = "+num2str(pho)+" (deg)\n"
	notes+="Azimuth = "+num2str(azo)+" (deg)\n"
	notes+="\n"
	notes+="Val. offsets\n"
	notes+="Polor = "+num2str(poof)+" (deg)\n"
	notes+="Tilt = "+num2str(phof)+" (deg)\n"
	notes+="Azimuth = "+num2str(azof)+" (deg)\n"
	notes+="\n"
	notes+="Reverse\n"
	notes+="Polor : "+K_revonoff_str(revpo)+"\n"
	notes+="Tilt : "+K_revonoff_str(revph)+"\n"
	notes+="Azimuth : "+K_revonoff_str(revaz)+"\n"
	notes+="____________________\n\n"
	
	notes+="ANALYZER\n\n"
	notes+="Analyzer rot. = "+num2str(ar)+" (deg)\n"
	notes+="Anlyzer type : "+K_analtype_str(type)+"\n"
	notes+="Alpha range = ["+num2str(th_s)+":"+num2str(th_e)+"] (deg)\n"
	notes+="Beta = "+num2str(Bet)+" (deg)\n"
	notes+="____________________\n\n"
	
	notes+="ENERGY\n\n"
	notes+="hn = "+num2str(hn)+" (eV)\n"
	notes+="V0 = "+num2str(V0)+" (eV)\n"
	notes+="W = "+num2str(W)+" (eV)\n"
	notes+="EB = "+num2str(EB)+" (eV)\n"
	
	SetDataFolder $nf
	
	return notes
End

Function/S K_make_note_map()
	String nf=GetDataFolder(1)
	SetDataFolder root:K_ARPES:misc
	Variable/g v_th_num
	variable th_num=v_th_num
	variable/g v_hn
	variable hn=v_hn
	variable/g v_V0
	variable V0=v_V0
	variable/g v_W
	variable W=v_W
	variable/g v_EB
	variable EB=v_EB
	variable/g v_th_s
	variable th_s=v_th_s
	variable/g v_th_e
	variable th_e=v_th_e
	variable/g v_bet
	variable bet=v_bet
	variable/g v_po
	variable po=v_po
	variable/g v_ph
	variable ph=v_ph
	variable/g v_az
	variable az=v_az
	variable/g v_poo
	variable poo=v_poo
	variable/g v_pho
	variable pho=v_pho
	variable/g v_azo
	variable azo=v_azo
	variable/g v_po_ofs
	variable poof=v_po_ofs
	variable/g v_ph_ofs
	variable phof=v_ph_ofs
	variable/g v_az_ofs
	variable azof=v_az_ofs
	variable/g v_revpo
	variable revpo=v_revpo
	variable/g v_revph
	variable revph=v_revph
	variable/g v_revaz
	variable revaz=v_revaz
	variable/g v_revth
	variable revth=v_revth
	variable/g v_analtype
	variable type=v_analtype
	variable/g v_anal_rot
	variable ar=v_anal_rot
	
	SetDataFolder root:K_ARPES:misc:automap
	variable/g v_amapval
	variable amap=v_amapval
	string/g s_amapval
	string amaps=s_amapval
	wave apw=$("amap_parameter")
	variable a_s=apw[0]
	variable a_st=apw[1]
	variable a_e=apw[2]
	string uni
	if(amap<6)
		uni=" (deg)"
	else
		uni=" (eV)"
	endif
	
	string notes="\nMAPPING\n\n"
	notes+="Variable : "+amaps+"\n"
	notes+="Start = "+num2str(a_s)+uni+"\n"
	notes+="Step = "+num2str(a_st)+uni+"\n"
	notes+="End = "+num2str(a_e)+uni+"\n"
	notes+="\nMANIPULATOR\n\n"
	notes+="Value\n"
	if(amap!=2)
		notes+="Polor = "+num2str(po)+" (deg)\n"
	endif
	if(amap!=3)
		notes+="Tilt = "+num2str(ph)+" (deg)\n"
	endif
	if(amap!=4)
		notes+="Azimuth = "+num2str(az)+" (deg)\n"
	endif
	notes+="\n"
	notes+="Val.@NML\n"
	notes+="Polor = "+num2str(poo)+" (deg)\n"
	notes+="Tilt = "+num2str(pho)+" (deg)\n"
	notes+="Azimuth = "+num2str(azo)+" (deg)\n"
	notes+="\n"
	notes+="Val. offsets\n"
	notes+="Polor = "+num2str(poof)+" (deg)\n"
	notes+="Tilt = "+num2str(phof)+" (deg)\n"
	notes+="Azimuth = "+num2str(azof)+" (deg)\n"
	notes+="\n"
	notes+="Reverse\n"
	notes+="Polor : "+K_revonoff_str(revpo)+"\n"
	notes+="Tilt : "+K_revonoff_str(revph)+"\n"
	notes+="Azimuth : "+K_revonoff_str(revaz)+"\n"
	notes+="____________________\n\n"
	
	notes+="ANALYZER\n\n"
	notes+="Analyzer rot. = "+num2str(ar)+" (deg)\n"
	notes+="Anlyzer type : "+K_analtype_str(type)+"\n"
	notes+="Alpha range = ["+num2str(th_s)+":"+num2str(th_e)+"] (deg)\n"
	if(amap!=5)
		notes+="Beta = "+num2str(Bet)+" (deg)\n"
	endif
	notes+="____________________\n\n"
	
	notes+="ENERGY\n\n"
	if(amap!=6)
		notes+="hn = "+num2str(hn)+" (eV)\n"
	endif
	notes+="V0 = "+num2str(V0)+" (eV)\n"
	notes+="W = "+num2str(W)+" (eV)\n"
	notes+="EB = "+num2str(EB)+" (eV)\n"
	
	SetDataFolder $nf
	
	return notes
End

Function/S K_make_note_kconv()
	String nf=GetDataFolder(1)
	SetDataFolder root:K_ARPES:global
	nvar v=v_kconv_vol
	string/g s_kconv_target
	string ws=s_kconv_target
	nvar en=v_kconv_en
	nvar v=v_kconv_vol
	nvar ef=v_EF
	nvar efd=v_EF_data
	nvar es=v_e_s
	nvar ee=v_e_e
	
	string notes="\nk-conv.\n\n"
	
	notes+="data : "+ws+"\n"
	notes+="\n"
	notes+="Ekin at EF = "+num2str(ef)+"\n"
	notes+="EF of data = "+num2str(efd)+"\n"
	if(v==0)
		notes+="Energy = "+num2str(en)+"\n"
	elseif(v==1)
		notes+="Energy range = ["+num2str(es)+":"+num2str(ee)+"] (eV)\n"
	endif
	
	notes+="____________________\n\n"
	
	
	return notes
End

Function/S K_revonoff_str(val)
	variable val
	string out
	
	if(val==1)
		out="OFF"
	elseif(val==-1)
		out="ON"
	endif
	
	return out
End

Function/S K_analtype_str(val)
	variable val
	string out
	
	if(val==1)
		out="Tilt-anlular notation"
	elseif(val==0)
		out="Polar-anlular notation"
	endif
	
	return out
End

Function K_find_normal_pos()
	string nf=GetDataFolder(1)
	SetDataFolder root:K_ARPES:misc
	variable/g v_alp_emis_set
	variable/g v_bet_emis_set
	
	variable alp_r,bet_r,alp,bet
	variable/g v_th_s
	variable/g v_th_e
	variable th=(v_th_e+v_th_s)/2
	variable th_r=th*pi/180
	
	variable po,ph
	
	variable pom,phm //min
	variable alp_emis,bet_emis
	variable alp_emis_m,bet_emis_m
	variable dmin=10000
	variable kx_s,ky_s,kz
	
	SetDataFolder root:K_ARPES:misc:rot_matrix
	wave aw=$("NML")
	SetDataFolder root:K_ARPES:misc
	
	for(po=-80;po<=80;po+=10)
		for(ph=-80;ph<=80;ph+=10)
			variable/g v_po=po
			variable/g v_ph=ph
			
			K_calc_normal()
			
			K_calc_emis(th)
			
			variable/g v_alp_emis
			variable/g v_bet_emis
			
			if(((v_alp_emis-v_alp_emis_set)^2+(v_bet_emis-v_bet_emis_set)^2)<dmin)
				dmin=(v_alp_emis-v_alp_emis_set)^2+(v_bet_emis-v_bet_emis_set)^2
				pom=po
				phm=ph
			endif
		endfor
	endfor
	
	for(po=pom-10;po<=pom+10;po+=1)
		for(ph=phm-10;ph<=phm+10;ph+=1)
			variable/g v_po=po
			variable/g v_ph=ph
			
			K_calc_normal()
			
			K_calc_emis(th)
			
			variable/g v_alp_emis
			variable/g v_bet_emis
			
			if(((v_alp_emis-v_alp_emis_set)^2+(v_bet_emis-v_bet_emis_set)^2)<dmin)
				dmin=(v_alp_emis-v_alp_emis_set)^2+(v_bet_emis-v_bet_emis_set)^2
				pom=po
				phm=ph
			endif
		endfor
	endfor
	
	for(po=pom-1;po<=pom+1;po+=0.1)
		for(ph=phm-1;ph<=phm+1;ph+=0.1)
			variable/g v_po=po
			variable/g v_ph=ph
			
			K_calc_normal()
			
			K_calc_emis(th)
			
			variable/g v_alp_emis
			variable/g v_bet_emis
			
			if(((v_alp_emis-v_alp_emis_set)^2+(v_bet_emis-v_bet_emis_set)^2)<dmin)
				dmin=(v_alp_emis-v_alp_emis_set)^2+(v_bet_emis-v_bet_emis_set)^2
				pom=po
				phm=ph
			endif
		endfor
	endfor
	
	for(po=pom-0.1;po<=pom+0.1;po+=0.01)
		for(ph=phm-0.1;ph<=phm+0.1;ph+=0.1)
			variable/g v_po=po
			variable/g v_ph=ph
			
			K_calc_normal()
			
			K_calc_emis(th)
			
			variable/g v_alp_emis
			variable/g v_bet_emis
			
			if(((v_alp_emis-v_alp_emis_set)^2+(v_bet_emis-v_bet_emis_set)^2)<dmin)
				dmin=(v_alp_emis-v_alp_emis_set)^2+(v_bet_emis-v_bet_emis_set)^2
				pom=po
				phm=ph
			endif
		endfor
	endfor

	for(po=pom-0.01;po<=pom+0.01;po+=0.001)
		for(ph=phm-0.01;ph<=phm+0.01;ph+=0.001)
			variable/g v_po=po
			variable/g v_ph=ph
			
			K_calc_normal()
			
			K_calc_emis(th)
			
			variable/g v_alp_emis
			variable/g v_bet_emis
			
			if(((v_alp_emis-v_alp_emis_set)^2+(v_bet_emis-v_bet_emis_set)^2)<dmin)
				dmin=(v_alp_emis-v_alp_emis_set)^2+(v_bet_emis-v_bet_emis_set)^2
				pom=po
				phm=ph
			endif
		endfor
	endfor
	
	for(po=pom-0.001;po<=pom+0.001;po+=0.0001)
		for(ph=phm-0.001;ph<=phm+0.001;ph+=0.0001)
			variable/g v_po=po
			variable/g v_ph=ph
			
			K_calc_normal()
			
			K_calc_emis(th)
			
			variable/g v_alp_emis
			variable/g v_bet_emis
			
			if(((v_alp_emis-v_alp_emis_set)^2+(v_bet_emis-v_bet_emis_set)^2)<dmin)
				dmin=(v_alp_emis-v_alp_emis_set)^2+(v_bet_emis-v_bet_emis_set)^2
				pom=po
				phm=ph
				
				alp_emis_m=v_alp_emis
				bet_emis_m=v_bet_emis
			endif
		endfor
	endfor
	
	variable/g v_po=round(pom*1000)/1000
	variable/g v_ph=round(phm*1000)/1000
	
	variable/g v_alp_emis_set=alp_emis
	variable/g v_bet_emis_set=bet_emis
	
	SetDataFolder $nf
End

Function K_Rodrigues_matrix(x,y,z,rot,r)
	variable x,y,z,rot,r
	string nf=GetDataFolder(1)
	SetDataFolder root:K_ARPES:misc:rot_matrix
	string ms="rotation_matrix"
	make/O/N=(3,3) $ms
	wave mw=$ms
	variable x1,y1,z1,rot_r
	if(r==0)
		rot_r=rot/180*pi
	elseif(r==1)
		rot_r=rot
	endif
	
	x1=x/sqrt(x^2+y^2+z^2)
	y1=y/sqrt(x^2+y^2+z^2)
	z1=z/sqrt(x^2+y^2+z^2)
	
	x=x1
	y=y1
	z=z1
	
	mw[0][0]=x^2*(1-cos(rot_r))+cos(rot_r)
	mw[0][1]=x*y*(1-cos(rot_r))-z*sin(rot_r)
	mw[0][2]=z*x*(1-cos(rot_r))+y*sin(rot_r)
	
	mw[1][0]=x*y*(1-cos(rot_r))+z*sin(rot_r)
	mw[1][1]=y^2*(1-cos(rot_r))+cos(rot_r)
	mw[1][2]=y*z*(1-cos(rot_r))-x*sin(rot_r)
	
	mw[2][0]=z*x*(1-cos(rot_r))-y*sin(rot_r)
	mw[2][1]=y*z*(1-cos(rot_r))+x*sin(rot_r)
	mw[2][2]=z^2*(1-cos(rot_r))+cos(rot_r)
	
	
	//					 / mw[0][0] mw[0][1] mw[0][2] \ 
	// R(x,y,z,rot)=	|  mw[1][0] mw[1][1] mw[1][2]  |
	//					 \ mw[2][0] mw[2][1] mw[2][2] /
	
	SetDataFolder $nf
End


//	Function K_rev_Rmat(ax,ay,az,a1,a2,a3,b1,b2,b3)
//		variable ax,ay,az,a1,a2,a3,b1,b2,b3
//		
//		variable ax1,ay1,az1
//		
//		ax1=ax/sqrt(ax^2+ay^2+az^2)
//		ay1=ay/sqrt(ax^2+ay^2+az^2)
//		az1=az/sqrt(ax^2+ay^2+az^2)
//		
//		ax=ax1
//		ay=ay1
//		az=az1
//		
//		variable xx,yy,zz
//		
//		xx=(ax*a1+ay*a2+az*a3)*ax-b1
//		yy=a1-ax*(ax*a1+ay*a2+az*a3)
//		zz=-az*a2+ay*a3
//		
//		variable cx1=(-xx*yy+sqrt((-xx^2+yy^2+zz^2)*zz^2))/(yy^2+zz^2)
//		variable cx2=(-xx*yy-sqrt((-xx^2+yy^2+zz^2)*zz^2))/(yy^2+zz^2)
//		
//		
//		xx=(ax*a1+ay*a2+az*a3)*ay-b2
//		yy=a2-ay*(ax*a1+ay*a2+az*a3)
//		zz=-ax*a3+az*a1
//		
//		variable cy1=(-xx*yy+sqrt((-xx^2+yy^2+zz^2)*zz^2))/(yy^2+zz^2)
//		variable cy2=(-xx*yy-sqrt((-xx^2+yy^2+zz^2)*zz^2))/(yy^2+zz^2)
//		
//		xx=(ax*a1+ay*a2+az*a3)*az-b3
//		yy=a3-az*(ax*a1+ay*a2+az*a3)
//		zz=-ay*a1+ax*a2
//		
//		variable cz1=(-xx*yy+sqrt((-xx^2+yy^2+zz^2)*zz^2))/(yy^2+zz^2)
//		variable cz2=(-xx*yy-sqrt((-xx^2+yy^2+zz^2)*zz^2))/(yy^2+zz^2)
//		
//		print cx1,cx2,cy1,cy2,cz1,cz2
//		print acos(cx1)*180/pi,acos(cx2)*180/pi,acos(cy1)*180/pi,acos(cy2)*180/pi,acos(cz1)*180/pi,acos(cz2)*180/pi
//	End

Function K_make_Rar()
	String nf=GetDataFolder(1)
	
	SetDataFolder root:K_ARPES:misc
	variable/g v_anal_rot
	variable anal_rot=v_anal_rot
	variable ar_r=anal_rot/180*pi
	
	SetDataFolder root:K_ARPES:misc:rot_matrix
	
	K_Rodrigues_matrix(0,0,-1,ar_r,1)
	
	wave rm1=$("rotation_matrix")
	
	Duplicate/O rm1,$("ar_rm")
	
	MatrixInverse rm1
	Duplicate/O $("M_Inverse"),$("ar_rm_inv")
	KillWaves $("M_Inverse")
	
	SetDataFolder $nf
End

Function K_make_Razofs()
	String nf=GetDataFolder(1)
	
	SetDataFolder root:K_ARPES:misc
	variable/g v_az_ofs
	variable azo=v_az_ofs*pi/180
	
	SetDataFolder root:K_ARPES:misc:rot_matrix
	
	K_Rodrigues_matrix(0,0,-1,azo,1)
	
	wave rm1=$("rotation_matrix")
	
	Duplicate/O rm1,$("azofs_rm")
	
	SetDataFolder $nf
End

Function K_make_Rpo()
	String nf=GetDataFolder(1)
	
	SetDataFolder root:K_ARPES:misc
	variable/g v_po
	variable/g v_po_ofs
	variable/g v_revpo
	variable po_r=(v_po-v_po_ofs)/180*pi*v_revpo
	
	SetDataFolder root:K_ARPES:misc:rot_matrix
	wave ar=$("ar_rm")
	
	variable axx,axy,axz,axx1,axy1,axz1
	
	axx=0
	axy=-1
	axz=0
	
	axx1=axx*ar[0][0]+axy*ar[0][1]+axz*ar[0][2]
	axy1=axx*ar[1][0]+axy*ar[1][1]+axz*ar[1][2]
	axz1=axx*ar[2][0]+axy*ar[2][1]+axz*ar[2][2]
	
	K_Rodrigues_matrix(axx1,axy1,axz1,po_r,1)
	
	wave rm1=$("rotation_matrix")
	
	Duplicate/O rm1,$("po_rm")
	
	SetDataFolder $nf
End

Function K_make_Rti()
	String nf=GetDataFolder(1)
	
	SetDataFolder root:K_ARPES:misc
	variable/g v_ph
	variable/g v_ph_ofs
	variable/g v_revph
	variable ph_r=(v_ph-v_ph_ofs)/180*pi*v_revph
	
	SetDataFolder root:K_ARPES:misc:rot_matrix
	wave ar=$("ar_rm")
	
	variable axx,axy,axz,axx1,axy1,axz1
	
	axx=1
	axy=0
	axz=0
	
	axx1=axx*ar[0][0]+axy*ar[0][1]+axz*ar[0][2]
	axy1=axx*ar[1][0]+axy*ar[1][1]+axz*ar[1][2]
	axz1=axx*ar[2][0]+axy*ar[2][1]+axz*ar[2][2]
	
	K_Rodrigues_matrix(axx1,axy1,axz1,ph_r,1)
	
	wave rm1=$("rotation_matrix")
	
	Duplicate/O rm1,$("ti_rm")
	
	SetDataFolder $nf
End

Function K_make_Raz()
	String nf=GetDataFolder(1)
	
	SetDataFolder root:K_ARPES:misc
	variable/g v_az
	variable/g v_azo
	variable/g v_revaz
	variable az=v_az-v_azo
	variable az_r=az/180*pi*v_revaz
	
	SetDataFolder root:K_ARPES:misc:rot_matrix
	wave ar=$("ar_rm")
	
	variable axx,axy,axz,axx1,axy1,axz1
	
	axx=0
	axy=0
	axz=-1
	
	axx1=axx*ar[0][0]+axy*ar[0][1]+axz*ar[0][2]
	axy1=axx*ar[1][0]+axy*ar[1][1]+axz*ar[1][2]
	axz1=axx*ar[2][0]+axy*ar[2][1]+axz*ar[2][2]
	
	K_Rodrigues_matrix(axx1,axy1,axz1,az_r,1)
	
	wave rm1=$("rotation_matrix")
	
	Duplicate/O rm1,$("az_rm")
	
	SetDataFolder $nf
End

Function K_make_Rpoo()
	String nf=GetDataFolder(1)
	
	SetDataFolder root:K_ARPES:misc
	variable/g v_poo
	variable/g v_po_ofs
	variable/g v_revpo
	variable poo_r=(v_poo-v_po_ofs)/180*pi*v_revpo
	
	SetDataFolder root:K_ARPES:misc:rot_matrix
	wave ar=$("ar_rm")
	
	variable axx,axy,axz,axx1,axy1,axz1
	
	axx=0
	axy=1
	axz=0
	
	axx1=axx*ar[0][0]+axy*ar[0][1]+axz*ar[0][2]
	axy1=axx*ar[1][0]+axy*ar[1][1]+axz*ar[1][2]
	axz1=axx*ar[2][0]+axy*ar[2][1]+axz*ar[2][2]
	
	K_Rodrigues_matrix(axx1,axy1,axz1,poo_r,1)
	
	wave rm1=$("rotation_matrix")
	
	Duplicate/O rm1,$("poo_rm")
	
	SetDataFolder $nf
End

Function K_make_Rtio()
	String nf=GetDataFolder(1)
	
	SetDataFolder root:K_ARPES:misc
	variable/g v_pho
	variable/g v_ph_ofs
	variable/g v_revph
	variable pho_r=(v_pho-v_ph_ofs)/180*pi*v_revph
	
	SetDataFolder root:K_ARPES:misc:rot_matrix
	wave ar=$("ar_rm")
	
	variable axx,axy,axz,axx1,axy1,axz1
	
	axx=-1
	axy=0
	axz=0
	
	axx1=axx*ar[0][0]+axy*ar[0][1]+axz*ar[0][2]
	axy1=axx*ar[1][0]+axy*ar[1][1]+axz*ar[1][2]
	axz1=axx*ar[2][0]+axy*ar[2][1]+axz*ar[2][2]
	
	K_Rodrigues_matrix(axx1,axy1,axz1,pho_r,1)
	
	wave rm1=$("rotation_matrix")
	
	Duplicate/O rm1,$("tio_rm")
	
	SetDataFolder $nf
End

Function K_make_window()
	string nf=GetDataFolder(1)
	SetDataFolder root:K_ARPES:misc
	string/g s_winname
	string winn=s_winname
	make/O/N=(0,3) $("curves")
	
	SetDataFolder root:K_ARPES:Curves
	if(DataFolderExists(winn))
		KillDataFolder $winn
	endif
	NewDataFolder/O $winn
	SetDataFolder $nf
	K_winlist()
End

Function K_winlist()
	string nf=GetDataFolder(1)
	SetDataFolder root:K_ARPES:global
	string wlws="win_list"
	wave/T wlw=$wlws
	SetDataFolder root:K_ARPES:misc
	string/g s_winname
	wave wp=$("winposi")
	
	variable nown=Dimsize(wlw,0)
	variable i
	string wnwn
	variable exis=0
	SetDataFolder root:K_ARPES:global
	for(i=0;i<nown;i+=1)
		wnwn=wlw[i]
		if(stringmatch(wnwn, s_winname))
			exis=1
		endif
	endfor
	
	if(exis==0)
		SetDataFolder root:K_ARPES:global
		InsertPoints nown,1,$wlws
		wlw[nown]=s_winname
		SetDataFolder root:K_ARPES:misc
		wp=wp+10
	endif
	
	SetDataFolder $nf
End

Function K_calcAxisRange_w()
	K_calcAxisRange()
	K_CheckAxis_w()
	K_select_crosssection()
	K_SetAxis_w()
End

Function K_calcAxisRange()
	string nf=GetDataFolder(1) 
	
	SetDataFolder root:K_ARPES:misc
	variable/g v_kxfixed
	variable kxfixed=v_kxfixed
	variable/g v_kyfixed
	variable kyfixed=v_kyfixed
	variable/g v_kzfixed
	variable kzfixed=v_kzfixed
	
	string/g s_winname
	string hlws=s_winname
	
	string BZzs,nds
	variable cnmax,cnmin
	
	SetDataFolder root:K_ARPES:global
	wave hlw=$hlws
	
	variable i
	variable kxmaxw,kymaxw,kzmaxw
	variable kxmax=-100,kymax=-100,kzmax=-100
	variable kxmaxs,kymaxs,kzmaxs
	
	variable kxminw,kyminw,kzminw
	variable kxmin=100,kymin=100,kzmin=100
	variable kxmins,kymins,kzmins
	
	SetDataFolder root:K_ARPES:Curves:$hlws
	i=0
	do
		string kxcws="C"+num2str(i)+"kx"
		string kycws="C"+num2str(i)+"ky"
		string kzcws="C"+num2str(i)+"kz"
		
		wave kxcw=$kxcws
		wave kycw=$kycws
		wave kzcw=$kzcws
		
		if(Exists(kxcws) && Exists(kycws) && Exists(kzcws))
			
			kxminw=WaveMin(kxcw)
			kxmaxw=WaveMax(kxcw)
			kyminw=WaveMin(kycw)
			kymaxw=WaveMax(kycw)
			kzminw=WaveMin(kzcw)
			kzmaxw=WaveMax(kzcw)
			
			if(kxmaxw>kxmax)
				kxmax=kxmaxw
			endif
			if(kymaxw>kymax)
				kymax=kymaxw
			endif
			if(kzmaxw>kzmax)
				kzmax=kzmaxw
			endif
			
			if(kxminw<kxmin)
				kxmin=kxminw
			endif
			if(kyminw<kymin)
				kymin=kyminw
			endif
			if(kzminw<kzmin)
				kzmin=kzminw
			endif
		
		elseif(!Exists(kxcws) && !Exists(kycws) && !Exists(kzcws))
			break
		else
			print "ERRPR: please initialize window"
			break
		endif
		
		i+=1
	while(1<2)
	
	if(kxmax<round(kxmax))
		kxmaxs=round(kxmax)
	else
		kxmaxs=round(kxmax)+1
	endif
	if(kxmin<=round(kxmin))
		kxmins=round(kxmin)-1
	else
		kxmins=round(kxmin)
	endif
	
	if(kymax<round(kymax))
		kymaxs=round(kymax)
	else
		kymaxs=round(kymax)+1
	endif
	if(kymin<=round(kymin))
		kymins=round(kymin)-1
	else
		kymins=round(kymin)
	endif
	
	if(kzmax<round(kzmax))
		kzmaxs=round(kzmax)
	else
		kzmaxs=round(kzmax)+1
	endif
	if(kzmin<=round(kzmin))
		kzmins=round(kzmin)-1
	else
		kzmins=round(kzmin)
	endif
	
	SetDataFolder root:K_ARPES:misc

	if(kxfixed==0)
		variable/g v_kxmax_w=kxmaxs
		variable/g v_kxmin_w=kxmins
	endif
	if(kyfixed==0)
		variable/g v_kymax_w=kymaxs
		variable/g v_kymin_w=kymins
	endif
	if(kzfixed==0)
		variable/g v_kzmax_w=kzmaxs
		variable/g v_kzmin_w=kzmins
	endif
	
	make/O/N=(2,3) $("krange")
	wave krw=$("krange")
	krw[0][0]=kxmin
	krw[1][0]=kxmax
	krw[0][1]=kymin
	krw[1][1]=kymax
	krw[0][2]=kzmin
	krw[1][2]=kzmax
	
	SetDataFolder $nf
End

Function K_CheckAxis_w()
	string nf=GetDataFolder(1)
	SetDataFolder root:K_ARPES:misc
	
	variable/g v_kxmax_w
	variable kxa=v_kxmax_w
	variable/g v_kymax_w
	variable kya=v_kymax_w
	variable/g v_kzmax_w
	variable kza=v_kzmax_w
	
	variable/g v_kxmin_w
	variable kxi=v_kxmin_w
	variable/g v_kymin_w
	variable kyi=v_kymin_w
	variable/g v_kzmin_w
	variable kzi=v_kzmin_w
	
	if(kxa<kxi)
		variable/g v_kxmax_w=kxi
		variable/g v_kxmin_w=kxa
	endif
	
	if(kya<kyi)
		variable/g v_kymax_w=kyi
		variable/g v_kymin_w=kya
	endif
	
	if(kza<kzi)
		variable/g v_kzmax_w=kzi
		variable/g v_kzmin_w=kza
	endif
	
	variable/g v_kxwid_w=kxa-kxi
	variable/g v_kywid_w=kya-kyi
	variable/g v_kzwid_w=kza-kzi
	
	SetDataFolder $nf
End

Function K_set_winrange_wid()
	string nf=GetDataFolder(1)
	SetDataFolder root:K_ARPES:misc
	variable/g v_kxwid_w
	variable/g v_kywid_w
	variable/g v_kzwid_w
	
	if(v_kxwid_w==0)
		v_kxwid_w=2
	endif
	if(v_kywid_w==0)
		v_kywid_w=2
	endif
	if(v_kzwid_w==0)
		v_kzwid_w=2
	endif
	
	variable/g v_kxmax_w
	variable/g v_kymax_w
	variable/g v_kzmax_w
	
	variable/g v_kxmin_w
	variable/g v_kymin_w
	variable/g v_kzmin_w
	
	variable kxc=(v_kxmax_w+v_kxmin_w)/2
	variable kyc=(v_kymax_w+v_kymin_w)/2
	variable kzc=(v_kzmax_w+v_kzmin_w)/2
	
	variable/g v_kxmax_w=kxc+abs(v_kxwid_w)/2
	variable/g v_kymax_w=kyc+abs(v_kywid_w)/2
	variable/g v_kzmax_w=kzc+abs(v_kzwid_w)/2
	
	variable/g v_kxmin_w=kxc-abs(v_kxwid_w)/2
	variable/g v_kymin_w=kyc-abs(v_kywid_w)/2
	variable/g v_kzmin_w=kzc-abs(v_kzwid_w)/2
	
	SetDataFolder $nf
End


Function K_SetAxis_w()
	string win
	String nf=GetDataFolder(1)

	SetDataFolder root:K_ARPES:misc
	variable/g v_kxmax_w
	variable kxa=v_kxmax_w
	variable/g v_kxmin_w
	variable kxi=v_kxmin_w
	variable/g v_kymax_w
	variable kya=v_kymax_w
	variable/g v_kymin_w
	variable kyi=v_kymin_w
	variable/g v_kzmax_w
	variable kza=v_kzmax_w
	variable/g v_kzmin_w
	variable kzi=v_kzmin_w
	string/g s_winname
	
	nvar vz=v_cross_z
	nvar vy=v_cross_y
	nvar vx=v_cross_x
	
	variable kxw=kxa-kxi
	variable kyw=kya-kyi
	variable kzw=kza-kzi
	
	if(kxw>0 && kyw>0 && kzw>0)
		string wnn=s_winname
		
		if(WinType(wnn)==1)
			
			ModifyGraph /W=$wnn standoff=0,mirror=2;DelayUpdate
			ModifyGraph /W=$wnn lblPosMode=2;DelayUpdate
			if(vx==1)
				SetAxis /W=$wnn xy_x kxi,kxa;DelayUpdate
				SetAxis /W=$wnn xy_y kyi,kya;DelayUpdate
				ModifyGraph /W=$wnn freePos(xy_x)={vz*kzw/(kyw+kzw),kwFraction};DelayUpdate
				ModifyGraph /W=$wnn freePos(xy_y)={vy*kzw/(kxw+kzw),kwFraction};DelayUpdate
				ModifyGraph /W=$wnn axisEnab(xy_x)={vy*kzw/(kxw+kzw),1};DelayUpdate
				ModifyGraph /W=$wnn axisEnab(xy_y)={vz*kzw/(kyw+kzw),1};DelayUpdate
				ModifyGraph /W=$wnn gridEnab(xy_x)={vz*kzw/(kyw+kzw),1};DelayUpdate
				ModifyGraph /W=$wnn gridEnab(xy_y)={vy*kzw/(kxw+kzw),1};DelayUpdate
				ModifyGraph /W=$wnn lblPosMode(xy_y)=1,lblPosMode(xy_x)=1;DelayUpdate
				if(vz==1)
					Label /W=$wnn xy_x "";DelayUpdate
					ModifyGraph /W=$wnn noLabel(xy_x)=1;DelayUpdate
				else
					Label /W=$wnn xy_x "\\f02k\\Bx\SS\\M\\f00 (/A)";DelayUpdate
				endif
				if(vy==1)
					Label /W=$wnn xy_y "";DelayUpdate
					ModifyGraph /W=$wnn noLabel(xy_y)=1;DelayUpdate
				else
					Label /W=$wnn xy_y "\\f02k\\By\SS\\M\\f00 (/A)";DelayUpdate
				endif
			endif
			
			if(vy==1)
				SetAxis /W=$wnn yz_y kyi,kya;DelayUpdate
				SetAxis /W=$wnn yz_z kzi,kza;DelayUpdate
				ModifyGraph /W=$wnn freePos(yz_y)={0,kwFraction};DelayUpdate
				ModifyGraph /W=$wnn freePos(yz_z)={vz*kzw/(kyw+kzw),kwFraction};DelayUpdate
				ModifyGraph /W=$wnn axisEnab(yz_y)={vz*kzw/(kyw+kzw),1};DelayUpdate
				ModifyGraph /W=$wnn axisEnab(yz_z)={0,kzw/(kxw+kzw)};DelayUpdate
				ModifyGraph /W=$wnn gridEnab(yz_y)={0,kzw/(kxw+kzw)};DelayUpdate
				ModifyGraph /W=$wnn gridEnab(yz_z)={vz*kzw/(kyw+kzw),1};DelayUpdate
				Label /W=$wnn yz_y "\\f02k\\By\SS\\M\\f00 (/A)";DelayUpdate
				Label /W=$wnn yz_z "\\f02k\\Bz\SS\\M\\f00 (/A)";DelayUpdate
				if(vz==1)
					ModifyGraph /W=$wnn lblPosMode(yz_z)=4,lblPos(yz_z)=40;DelayUpdate
				endif
			endif
			
			if(vz==1)
				SetAxis /W=$wnn xz_x kxi,kxa;DelayUpdate
				SetAxis /W=$wnn xz_z kzi,kza;DelayUpdate
				ModifyGraph /W=$wnn freePos(xz_x)={0,kwFraction};DelayUpdate
				ModifyGraph /W=$wnn freePos(xz_z)={0,kwFraction};DelayUpdate
				ModifyGraph /W=$wnn axisEnab(xz_x)={vy*kzw/(kxw+kzw),1};DelayUpdate
				ModifyGraph /W=$wnn axisEnab(xz_z)={0,1-vx*kyw/(kyw+kzw)};DelayUpdate
				ModifyGraph /W=$wnn gridEnab(xz_x)={0,1-vx*kyw/(kyw+kzw)};DelayUpdate
				ModifyGraph /W=$wnn gridEnab(xz_z)={0,1-vy*kzw/(kxw+kzw)};DelayUpdate
				if(vy==0)
					ModifyGraph /W=$wnn mirror(xz_z)=2;DelayUpdate
				else
					ModifyGraph /W=$wnn mirror(xz_z)=0;DelayUpdate
				endif
				Label /W=$wnn xz_z "\\f02k\\Bz\SS\\M\\f00 (/A)";DelayUpdate
				Label /W=$wnn xz_x "\\f02k\\Bx\SS\\M\\f00 (/A)";DelayUpdate
			endif
			
			if(vx==1)
				ModifyGraph /W=$wnn width={Plan,1,xy_x,xy_y};DelayUpdate
			elseif(vz==1)
				ModifyGraph /W=$wnn width={Plan,1,xz_x,xz_z};DelayUpdate
			else
				ModifyGraph /W=$wnn width={Plan,1,yz_z,yz_y};DelayUpdate
			endif
			ModifyGraph /W=$wnn tick=3,grid=1;DelayUpdate
			ModifyGraph /W=$wnn margin(bottom)=34,margin(right)=42,margin(left)=42
		endif
	else
		SetDataFolder root:K_ARPES:misc
		if(kxw<=0)
			variable/g v_kxmax_w=kxi+1
		endif
		if(kyw<=0)
			variable/g v_kymax_w=kyi+1
		endif
		if(kzw<=0)
			variable/g v_kzmax_w=kzi+1
		endif
		
		K_SetAxis_w()
	endif
	
	SetDataFolder $nf
End

Function K_append_curve()
	String fldrSav0= GetDataFolder(1)
	SetDataFolder root:K_ARPES:misc
	variable/g v_hn
	variable hn=v_hn
	variable/g v_V0
	variable V0=v_V0
	variable/g v_W
	variable W=v_W
	variable/g v_EB
	variable EB=v_EB
	variable/g v_th_s
	variable th_s=v_th_s
	variable/g v_th_e
	variable th_e=v_th_e
	string/g s_thsnum
	string thsnum=s_thsnum
	string/g s_thenum
	string thenum=s_thenum
	
	variable/g v_po
	variable po=v_po
	variable/g v_ph
	variable ph=v_ph
	variable/g v_az
	variable az=v_az
	variable/g v_poo
	variable poo=v_poo
	variable/g v_pho
	variable pho=v_pho
	variable/g v_azo
	variable azo=v_azo
	
	string/g s_winname
	string hlws=s_winname
	variable/g v_map
	variable map=v_map
	
	variable/g v_labhn
	variable labhn=v_labhn
	variable/g v_labV0
	variable labV0=v_labV0
	variable/g v_labW
	variable labW=v_labW
	variable/g v_labEB
	variable labEB=v_labEB
	variable/g v_labth
	variable labth=v_labth
	variable/g v_labph
	variable labph=v_labph
	
	SetDataFolder root:K_ARPES:global
	
	string clws="color_list"
	wave clw=$clws
	
	SetDataFolder root:K_ARPES:misc
	variable/g v_flag_app
	variable flag_app=v_flag_app
	variable/g v_curve_num
	variable cn=v_curve_num
	variable/g v_curve_num_map
	variable cnmap=v_curve_num_map
	wave winp=$("winposi")
	variable wpx=winp[0]
	variable wpy=winp[1]
	variable cr,cb,cg
	wave cvsw=$("curves")
	
	if(map==0)
		cr=clw[mod(cn,10)][0]
		cg=clw[mod(cn,10)][1]
		cb=clw[mod(cn,10)][2]
	elseif(map==1)
		cr=clw[mod(cnmap,10)][0]
		cg=clw[mod(cnmap,10)][1]
		cb=clw[mod(cnmap,10)][2]
	endif
		
	SetDataFolder root:K_ARPES:Curves
	SetDataFolder $hlws
	
	string cwxs="C"+num2str(cn)+"kx"
	string cwys="C"+num2str(cn)+"ky"
	string cwzs="C"+num2str(cn)+"kz"
	
	wave cwx=$cwxs
	wave cwy=$cwys
	wave cwz=$cwzs
	variable gr
		
	if(flag_app==1 || WinType(hlws)==0)
		variable px,py,pz
		px=cwx[0]
		py=cwy[0]
		pz=cwz[0]
		
		if(WinType(hlws)==0)
		
			Display/W=(wpx,wpy,wpx+400,wpy+400)/K=1/N=$hlws/R=xz_z/B=xz_x $cwzs vs $cwxs as hlws
			AppendToGraph/L=yz_y/B=yz_z $cwys vs $cwzs
			AppendToGraph/L=xy_y/B=xy_x/VERT $cwxs vs $cwys
			ModifyGraph tick(xy_x)=1
			ModifyGraph tick(xy_y)=1
			ModifyGraph standoff=0
			ModifyGraph grid=1
			ModifyGraph mirror(xz_z)=0,mirror(yz_z)=2,mirror(xy_y)=2,mirror(xy_x)=2
			ModifyGraph noLabel(xy_y)=1,noLabel(xy_x)=1
			ModifyGraph lblPosMode(xz_z)=3,lblPosMode(xz_x)=1,lblPosMode(yz_y)=1,lblPosMode(yz_z)=3
			ModifyGraph lblPos(xz_z)=50,lblPos(yz_z)=40
			ModifyGraph freePos(xz_x)={0,kwFraction}
			ModifyGraph freePos(yz_y)={0,kwFraction}
			if(th_s==th_e)
				ModifyGraph mode($cwxs)=3,marker($cwxs)=42
				ModifyGraph mode($cwys)=3,marker($cwys)=42
				ModifyGraph mode($cwzs)=3,marker($cwzs)=42
			else
				ModifyGraph mode($cwxs)=0,lsize($cwxs)=1.2
				ModifyGraph mode($cwys)=0,lsize($cwxs)=1.2
				ModifyGraph mode($cwzs)=0,lsize($cwxs)=1.2
			endif
		else
			DoWindow/F $hlws
			
			AppendToGraph/R=xz_z/B=xz_x $cwzs vs $cwxs
			AppendToGraph/L=yz_y/B=yz_z $cwys vs $cwzs
			AppendToGraph/L=xy_y/B=xy_x/VERT $cwxs vs $cwys
			if(th_s==th_e)
				ModifyGraph mode($cwxs)=3,marker($cwxs)=42
				ModifyGraph mode($cwys)=3,marker($cwys)=42
				ModifyGraph mode($cwzs)=3,marker($cwzs)=42
			else
				ModifyGraph mode($cwxs)=0,lsize($cwxs)=1.2
				ModifyGraph mode($cwys)=0,lsize($cwys)=1.2
				ModifyGraph mode($cwzs)=0,lsize($cwzs)=1.2
			endif
			
			ModifyGraph/W=$hlws rgb($cwxs)=(cr,cg,cb)
			ModifyGraph/W=$hlws rgb($cwys)=(cr,cg,cb)
			ModifyGraph/W=$hlws rgb($cwzs)=(cr,cg,cb)
		endif
		
		variable s=DimSize(cvsw,0)
		InsertPoints s,1,cvsw
		cvsw[s][0]=cr
		cvsw[s][1]=cg
		cvsw[s][2]=cb
		
	endif
	
	
	
	SetDataFolder $fldrSav0
EndMacro

Function K_select_crosssection()
	string nf=GetDataFolder(1)
	SetDataFolder root:K_ARPES:global
	wave clw=$("color_list")
	SetDataFolder root:K_ARPES:misc
	string/g s_winname
	string winn=s_winname
	wave cvsw=$("curves")
	
	variable/g v_cross_z
	variable/g v_cross_y
	variable/g v_cross_x
	variable vz=v_cross_z
	variable vy=v_cross_y
	variable vx=v_cross_x
	
	variable i,cr,cg,cb
	string trl,trs,trlx="",trly="",trlz=""
	if(WinType(winn)!=0)
		SetDataFolder root:K_ARPES:Curves:$winn
		trl=TraceNameList(winn,";",1)
		i=0
		do
			trs=StringFromList(i,trl)
			if(!StrLen(trs)>0)
				break
			endif
			if(StringMatch(trs,"*kx"))
				RemoveFromGraph /W=$winn $trs
				trs=ReplaceString("kx",trs,"")
				trlx=AddListItem(trs,trlx)
			endif
			if(StringMatch(trs,"*ky"))
				RemoveFromGraph /W=$winn $trs
				trs=ReplaceString("ky",trs,"")
				trly=AddListItem(trs,trly)
			endif
			if(StringMatch(trs,"*kz"))
				RemoveFromGraph /W=$winn $trs
				trs=ReplaceString("kz",trs,"")
				trlz=AddListItem(trs,trlz)
			endif
			i+=1
		while(1)
		
		trlx=SortList(trlx,";",16)
		trly=SortList(trly,";",16)
		trlz=SortList(trlz,";",16)
		
		if(StrLen(trlx)>0)
			trl=trlx
		elseif(StrLen(trly)>0)
			trl=trly
		else
			trl=trlz
		endif
		
		i=0
		do
			trs=StringFromList(i,trl)
			if(!StrLen(trs)>0)
				break
			endif
			
			cr=cvsw[i][0]
			cg=cvsw[i][1]
			cb=cvsw[i][2]
			
			trs=ReplaceString("kx",trs,"")
			if(vx==1)
				AppendToGraph/W=$winn/L=xy_y/B=xy_x/VERT $(trs+"kx") vs $(trs+"ky")
				if(DimSize($(trs+"kx"),0)>1)
					ModifyGraph /W=$winn mode($(trs+"kx"))=0,lsize($(trs+"kx"))=1.2
				else
					ModifyGraph /W=$winn mode($(trs+"kx"))=3,marker($(trs+"kx"))=42
				endif
				ModifyGraph /W=$winn rgb($(trs+"kx"))=(cr,cg,cb)
			endif
			if(vy==1)
				AppendToGraph/W=$winn/L=yz_y/B=yz_z $(trs+"ky") vs $(trs+"kz")
				if(DimSize($(trs+"ky"),0)>1)
					ModifyGraph /W=$winn mode($(trs+"ky"))=0,lsize($(trs+"ky"))=1.2
				else
					ModifyGraph /W=$winn mode($(trs+"ky"))=3,marker($(trs+"ky"))=42
				endif
				ModifyGraph /W=$winn rgb($(trs+"ky"))=(cr,cg,cb)
			endif
			if(vz==1)
				AppendToGraph/W=$winn/R=xz_z/B=xz_x $(trs+"kz") vs $(trs+"kx")
				if(DimSize($(trs+"kz"),0)>1)
					ModifyGraph /W=$winn mode($(trs+"kz"))=0,lsize($(trs+"kz"))=1.2
				else
					ModifyGraph /W=$winn mode($(trs+"kz"))=3,marker($(trs+"kz"))=42
				endif
				ModifyGraph /W=$winn rgb($(trs+"kz"))=(cr,cg,cb)
			endif
			i+=1
		while(1)
		K_save_misc()
	endif
	
	SetDataFolder $nf
End

Function K_change_mode()
	string nf=GetDataFolder(1)
	SetDataFolder root:K_ARPES:misc
	
	string/g s_winname
	string hlws=s_winname
	variable/g v_th_s
	variable th_s=v_th_s
	variable/g v_th_e
	variable th_e=v_th_e
	variable/g v_curve_num
	variable cn=v_curve_num
	
	nvar vx=v_cross_x
	nvar vy=v_cross_y
	nvar vz=v_cross_z
	
	string cwxs="C"+num2str(cn)+"kx"
	string cwys="C"+num2str(cn)+"ky"
	string cwzs="C"+num2str(cn)+"kz"
//	
//	if(th_s==th_e)
//		if(vx==1)
//			ModifyGraph /W=$hlws mode($cwxs)=3,marker($cwxs)=42
//		endif
//		if(vy==1)
//			ModifyGraph /W=$hlws mode($cwys)=3,marker($cwys)=42
//		endif
//		if(vz==1)
//			ModifyGraph /W=$hlws mode($cwzs)=3,marker($cwzs)=42
//		endif
//	else
//		if(vx==1)
//			ModifyGraph /W=$hlws mode($cwxs)=0
//		endif
//		if(vy==1)
//			ModifyGraph /W=$hlws mode($cwys)=0
//		endif
//		if(vz==1)
//			ModifyGraph /W=$hlws mode($cwzs)=0
//		endif
//	endif
	
	SetDataFolder $nf
End

Function K_make_color_list()
	string nf=GetDataFolder(1)
	SetDataFolder root:K_ARPES:global
	string ws="color_list"
	make/O/N=(10,3) $ws
	wave w=$ws
	
	// red
	w[0][0]=65535 // R
	w[0][1]=0 // G
	w[0][2]=0 // B
	
	// blue
	w[1][0]=0 // R
	w[1][1]=0 // G
	w[1][2]=65535 // B
	
	// dark green
	w[2][0]=26205 // R
	w[2][1]=52428 // G
	w[2][2]=1 // B
	
	// orange
	w[3][0]=65535 // R
	w[3][1]=43690 // G
	w[3][2]=1 // B
	
	// purple
	w[4][0]=29524 // R
	w[4][1]=1 // G
	w[4][2]=58982 // B
	
	// cyan
	w[5][0]=0 // R
	w[5][1]=65535 // G
	w[5][2]=65535 // B
	
	// yellow
	w[6][0]=52428 // R
	w[6][1]=52425 // G
	w[6][2]=1 // B
	
	// magenta
	w[7][0]=65535 // R
	w[7][1]=0 // G
	w[7][2]=65535 // B
	
	// black
	w[8][0]=0 // R
	w[8][1]=0 // G
	w[8][2]=0 // B
	
	// green
	w[9][0]=0 // R
	w[9][1]=65535 // G
	w[9][2]=0 // B
	
	SetDataFolder $nf
End

Function K_find_curve()
	string nf=GetDataFolder(1)
	SetDataFolder root:K_ARPES:misc
	variable/g v_hn
	variable hn=v_hn
	variable/g v_V0
	variable V0=v_V0
	variable/g v_W
	variable W=v_W
	variable/g v_EB
	variable EB=v_EB
	variable/g v_th_s
	variable th_s=v_th_s
	variable/g v_th_e
	variable th_e=v_th_e
	
	variable/g v_po
	variable po=v_po
	variable/g v_ph
	variable ph=v_ph
	variable/g v_az
	variable az=v_az
	variable/g v_poo
	variable poo=v_poo
	variable/g v_pho
	variable pho=v_pho
	variable/g v_azo
	variable azo=v_azo
	
	variable/g v_app
	variable app=v_app
	variable/g v_map
	variable map=v_map
	
	string/g s_winname
	string hlws=s_winname
	
	SetDataFolder root:K_ARPES:Curves:$hlws
	
	variable i,flag_app,cn
	string cxws
	
	flag_app=1
	
	i=0
	do
		cxws="C"+num2str(i)+"kx"
		if(Exists(cxws))
			if(StringMatch(note($cxws),K_make_note()))
//				flag_app=0
				cn=i
			endif
		else
			break
		endif
		i+=1
	while(1)
	
	if(flag_app==1)
		cn=i
	endif
	
	if(app==0)
		if(i==0)
			flag_app=1
			cn=0
		else
			flag_app=0
			cn=i-1
		endif
	endif
	
	SetDataFolder root:K_ARPES:misc
	variable/g v_flag_app=flag_app
	variable/g v_curve_num=cn
	if(map==0)
		variable/g v_curve_num_map=cn
	endif
	SetDataFolder $nf
End

Function K_SetVarProc_calc_kline(sva) : SetVariableControl
	STRUCT WMSetVariableAction &sva

	switch( sva.eventCode )
		case 1: // mouse up
		case 2: // Enter key
		case 3: // Live update
			Variable dval = sva.dval
			String sval = sva.sval
			K_ARPES_calc_kline()
			break
		case -1: // control being killed
			break
	endswitch

	return 0
End

Function K_SetVarProc_SetAxix_w(sva) : SetVariableControl
	STRUCT WMSetVariableAction &sva

	switch( sva.eventCode )
		case 1: // mouse up
		case 2: // Enter key
		case 3: // Live update
			Variable dval = sva.dval
			String sval = sva.sval
			string nf=GetDataFolder(1)
			if(StringMatch(sva.ctrlName,"setvar_kxwid")||StringMatch(sva.ctrlName,"setvar_kywid")||StringMatch(sva.ctrlName,"setvar_kzwid"))
				K_set_winrange_wid()
			endif
			K_CheckAxis_w()
			SetDataFolder root:K_ARPES:misc
			string/g s_winname
			if(WinType(s_winname)==1)
				K_SetAxis_w()
			endif
			SetDataFolder $nf
			break
		case -1: // control being killed
			break
	endswitch

	return 0
End

Function K_set_thxnum(th_s,th_e)
	variable th_s,th_e
	
	string nf=GetDataFolder(1)
	SetDataFolder root:K_ARPES:misc
	
	if(th_s<0)
		string/g s_thsnum="m"+num2str(th_s*-1)
	elseif(th_s>0)
		string/g s_thsnum="p"+num2str(th_s)
	elseif(th_s==0)
		string/g s_thsnum=num2str(th_s)
	endif
	
	if(th_e<0)
		string/g s_thenum="m"+num2str(th_e*-1)
	elseif(th_e>0)
		string/g s_thenum="p"+num2str(th_e)
	elseif(th_e==0)
		string/g s_thenum=num2str(th_e)
	endif
	
	SetDataFolder $nf
End
Function K_CheckProc_fixkz(cba) : CheckBoxControl
	STRUCT WMCheckboxAction &cba

	switch( cba.eventCode )
		case 2: // mouse up
			Variable checked = cba.checked
			string nf=GetDataFolder(1)
			SetDataFolder root:K_ARPES:misc
			string/g s_winname
			if(checked==0 && WinType(s_winname)==1)
				K_calcAxisRange_w()
			endif
			SetDataFolder $nf
			break
		case -1: // control being killed
			break
	endswitch

	return 0
End








Function K_ButtonProc_ksp(ba) : ButtonControl
	STRUCT WMButtonAction &ba

	switch( ba.eventCode )
		case 2: // mouse up
			// click code here
			string nf=GetDataFolder(1)
			SetDataFolder root:K_ARPES:misc
			variable/g v_lc_a=pi
			variable/g v_lc_b=pi
			variable/g v_lc_c=pi
			SetDataFolder $nf
			break
		case -1: // control being killed
			break
	endswitch

	return 0
End


Function K_CheckProc_app(cba) : CheckBoxControl
	STRUCT WMCheckboxAction &cba

	switch( cba.eventCode )
		case 2: // mouse up
			Variable checked = cba.checked
			string nf=GetDataFolder(1)
			SetDataFolder root:K_ARPES:misc
			if(checked==0)
				variable/g v_map=0
			endif
			SetDataFolder $nf
			break
		case -1: // control being killed
			break
	endswitch

	return 0
End


Function K_CheckProc_map(cba) : CheckBoxControl
	STRUCT WMCheckboxAction &cba

	switch( cba.eventCode )
		case 2: // mouse up
			Variable checked = cba.checked
			string nf=GetDataFolder(1)
			SetDataFolder root:K_ARPES:misc
			if(checked==1)
				variable/g v_app=1
			endif
			SetDataFolder $nf
			break
		case -1: // control being killed
			break
	endswitch

	return 0
End


Function K_CheckProc_revpo(cba) : CheckBoxControl
	STRUCT WMCheckboxAction &cba

	switch( cba.eventCode )
		case 2: // mouse up
			Variable checked = cba.checked
			string nf=GetDatafolder(1)
			SetDataFolder root:K_ARPES:misc:automap
			variable/g v_amapval
			variable amap=v_amapval
			SetDataFolder root:K_ARPES:misc
			string/g s_winname
			variable/g v_app
			variable/g v_map
			variable/g v_smh
			nvar v_kconv=root:K_ARPES:global:v_kconv
			
			if(checked==1)
				variable/g v_revpo=-1
			elseif(checked==0)
				variable/g v_revpo=1
			endif
			
			K_set_emis_ang_window()
			
			if(!WinType(s_winname)==0 && v_kconv==0)
				if(amap==1)
					if(v_app==0)
						K_ARPES_calc_kline()
						K_append_curve()
					elseif(v_map==1)
						K_ARPES_calc_kline()
						K_append_curve()
					endif
				else
					if(v_app==0)
						K_auto_mapping(1)
					endif
				endif
			endif
			
			SetDataFolder $nf
			break
		case -1: // control being killed
			break
	endswitch

	return 0
End

Function K_CheckProc_revaz(cba) : CheckBoxControl
	STRUCT WMCheckboxAction &cba

	switch( cba.eventCode )
		case 2: // mouse up
			Variable checked = cba.checked
			string nf=GetDatafolder(1)
			SetDataFolder root:K_ARPES:misc:automap
			variable/g v_amapval
			variable amap=v_amapval
			SetDataFolder root:K_ARPES:misc
			string/g s_winname
			variable/g v_app
			variable/g v_map
			variable/g v_smh
			nvar v_kconv=root:K_ARPES:global:v_kconv
			
			if(checked==1)
				variable/g v_revaz=-1
			elseif(checked==0)
				variable/g v_revaz=1
			endif
			
			K_set_emis_ang_window()
			
			if(!WinType(s_winname)==0 && v_kconv==0)
				if(amap==1)
					if(v_app==0)
						K_ARPES_calc_kline()
						K_append_curve()
					elseif(v_map==1)
						K_ARPES_calc_kline()
						K_append_curve()
					endif
				else
					if(v_app==0)
						K_auto_mapping(1)
					endif
				endif
			endif
			
			SetDataFolder $nf
			break
		case -1: // control being killed
			break
	endswitch

	return 0
End

Function K_CheckProc_revph(cba) : CheckBoxControl
	STRUCT WMCheckboxAction &cba

	switch( cba.eventCode )
		case 2: // mouse up
			Variable checked = cba.checked
			string nf=GetDatafolder(1)
			SetDataFolder root:K_ARPES:misc:automap
			variable/g v_amapval
			variable amap=v_amapval
			SetDataFolder root:K_ARPES:misc
			string/g s_winname
			variable/g v_app
			variable/g v_map
			variable/g v_smh
			nvar v_kconv=root:K_ARPES:global:v_kconv
			
			if(checked==1)
				variable/g v_revph=-1
			elseif(checked==0)
				variable/g v_revph=1
			endif
			
			K_set_emis_ang_window()
			
			if(!WinType(s_winname)==0 && v_kconv==0)
				if(amap==1)
					if(v_app==0)
						K_ARPES_calc_kline()
						K_append_curve()
					elseif(v_map==1)
						K_ARPES_calc_kline()
						K_append_curve()
					endif
				else
					if(v_app==0)
						K_auto_mapping(1)
					endif
				endif
			endif
			
			SetDataFolder $nf
			break
		case -1: // control being killed
			break
	endswitch

	return 0
End


Function K_CheckProc_revth(cba) : CheckBoxControl
	STRUCT WMCheckboxAction &cba

	switch( cba.eventCode )
		case 2: // mouse up
			Variable checked = cba.checked
			string nf=GetDatafolder(1)
			SetDataFolder root:K_ARPES:misc:automap
			variable/g v_amapval
			variable amap=v_amapval
			SetDataFolder root:K_ARPES:misc
			string/g s_winname
			variable/g v_app
			variable/g v_map
			variable/g v_smh
			nvar v_kconv=root:K_ARPES:global:v_kconv
			
			if(checked==1)
				variable/g v_revth=-1
			elseif(checked==0)
				variable/g v_revth=1
			endif
			
			K_set_emis_ang_window()
			
			if(!WinType(s_winname)==0 && v_kconv==0)
				if(amap==1)
					if(v_app==0)
						K_ARPES_calc_kline()
						K_append_curve()
					elseif(v_map==1)
						K_ARPES_calc_kline()
						K_append_curve()
					endif
				else
					if(v_app==0)
						K_auto_mapping(1)
					endif
				endif
			endif
			
			SetDataFolder $nf
			break
		case -1: // control being killed
			break
	endswitch

	return 0
End

Function K_CheckProc_revbet(cba) : CheckBoxControl
	STRUCT WMCheckboxAction &cba

	switch( cba.eventCode )
		case 2: // mouse up
			Variable checked = cba.checked
			string nf=GetDatafolder(1)
			SetDataFolder root:K_ARPES:misc:automap
			variable/g v_amapval
			variable amap=v_amapval
			SetDataFolder root:K_ARPES:misc
			string/g s_winname
			variable/g v_app
			variable/g v_map
			variable/g v_smh
			nvar v_kconv=root:K_ARPES:global:v_kconv
			
			if(checked==1)
				variable/g v_revbet=-1
			elseif(checked==0)
				variable/g v_revbet=1
			endif
			
			K_set_emis_ang_window()
			
			if(!WinType(s_winname)==0 && v_kconv==0)
				if(amap==1)
					if(v_app==0)
						K_ARPES_calc_kline()
						K_append_curve()
					elseif(v_map==1)
						K_ARPES_calc_kline()
						K_append_curve()
					endif
				else
					if(v_app==0)
						K_auto_mapping(1)
					endif
				endif
			endif
			
			SetDataFolder $nf
			break
		case -1: // control being killed
			break
	endswitch

	return 0
End

Function K_ButtonProc_sample(ba) : ButtonControl
	STRUCT WMButtonAction &ba

	switch( ba.eventCode )
		case 2: // mouse up
			// click code here
			string  nf=GetDataFolder(1)
			SetDataFolder root:K_ARPES:misc
			string/g s_winname
			
			K_find_normal_pos()
			if(!WinType(s_winname)==0)
				K_ARPES_calc_kline()
				K_append_curve()
			endif
			
			SetDataFolder $nf
			break
		case -1: // control being killed
			break
	endswitch

	return 0
End


Function/S K_winlist_str()
	string nf=GetDataFolder(1)
	SetDataFolder root:K_ARPES:global
	string wlws="win_list"
	wave/T wlw=$wlws
	variable n
	string str=""
	
	if(DimSize(wlw,0)>0)
		for(n=0;n<DimSize(wlw,0)-1;n+=1)
			str+=wlw[n]+";"
		endfor
		
		str+=wlw[DimSize(wlw,0)-1]
	else
		str="none"
	endif
	
	SetDataFolder $nf
	
	return str
End




Function K_ButtonProc_init(ba) : ButtonControl
	STRUCT WMButtonAction &ba

	switch( ba.eventCode )
		case 2: // mouse up
			// click code here
			SetDataFolder root:
			string nf=GetDataFolder(1)
			
			SetDataFolder root:K_ARPES:global
			wave/T wl=$("win_list")
			variable iend=DimSize(wl,0)
			variable i
			for(i=0;i<iend;i+=1)
				string wins=wl[i]
				if(WinType(wins)==1)
					KillWindow $wins
				endif
			endfor
			
			if(WinType("K_ARPES_p")==7)
				Killwindow K_ARPES_p
			endif
			if(DataFolderExists("root:K_ARPES:Curves"))
				KillDataFolder root:K_ARPES:Curves
			endif
			if(DataFolderExists("root:K_ARPES:global"))
				KillDataFolder root:K_ARPES:global
			endif
			if(DataFolderExists("root:K_ARPES:misc_saved"))
				KillDataFolder root:K_ARPES:misc_saved
			endif
			
			execute "K_ARPES_start()"
			SetDataFolder $nf
			
			Print " "
			Print "   --- K_ARPES has been initialized ---"
			Print " "
			
			break
		case -1: // control being killed
			break
	endswitch

	return 0
End

Function K_kconv_prep(tws)
	string tws
	wave tw1=$tws
	nvar img=root:K_ARPES:global:v_kconv_img
	nvar outr=root:K_ARPES:global:v_outr
	
	Duplicate/O tw1,$(tws+"_dupl")
	tws+="_dupl"
	
	wave tw=$tws
	
	InsertPoints/M=1 0,1,tw
	InsertPoints/M=1 DimSize(tw,1),1,tw
	SetScale/P y,DimOffset(tw1,1)-DimDelta(tw1,1),DimDelta(tw1,1),tw
	
	if(img!=1)
		InsertPoints/M=2 0,1,tw
		InsertPoints/M=2 DimSize(tw,2),1,tw
		SetScale/P z,DimOffset(tw1,2)-DimDelta(tw1,2),DimDelta(tw1,2),tw
	endif
	
	variable be,bx,by
	if(img==1)
		for(be=0;be<Dimsize(tw,0);be+=1)
			tw[be][0]=outr
			tw[be][DimSize(tw,1)-1]=outr
		endfor
	else
		for(be=0;be<Dimsize(tw,0);be+=1)
			for(bx=0;bx<DimSize(tw,1);bx+=1)
				tw[be][bx][0]=outr
				tw[be][bx][DimSize(tw,2)-1]=outr
			endfor
			for(by=0;by<DimSize(tw,2);by+=1)
				tw[be][0][by]=outr
				tw[be][DimSize(tw,1)-1][by]=outr
			endfor
		endfor
	endif
	
	variable e_s=DimOffset(tw,0)
	variable e_st=DimDelta(tw,0)
	variable e_n=DimSize(tw,0)
	variable t1_s=DimOffset(tw,1)
	variable t1_st=DimDelta(tw,1)
	variable t1_n=DimSize(tw,1)
	variable t1_e=t1_s+(t1_n-1)*t1_st
	if(img!=1)
		variable t2_s=DimOffset(tw,2)
		variable t2_st=DimDelta(tw,2)
		variable t2_n=DimSize(tw,2)
		variable t2_e=t2_s+(t2_n-1)*t2_st
	endif
	
	string nf=GetDataFolder(1)
	SetDataFolder root:K_ARPES:misc
	string/g s_winname
	string logs_winn=s_winname
	string/g s_winname="k_saved"
	K_save_misc()
	SetDataFolder root:K_ARPES:misc_saved:k_saved
	string/g s_winname=logs_winn
	SetDataFolder root:K_ARPES:misc
	string/g s_winname="k"
	variable/g v_th_s=t1_s
	variable/g v_th_e=t1_e
	
	if(img!=1)
		SetDataFolder root:K_ARPES:misc:automap
		wave amp=$("amap_parameter")
		variable s=amp[0]
		variable st=amp[1]
		variable e=amp[2]
		variable n=amp[3]
		if(img!=1)
			amp[0]=s-st
			amp[1]=st
			amp[2]=e+st
			amp[3]=n+2
		endif
	endif
	
	K_make_window()
	K_save_misc()
	SetDataFolder $nf
End

Function K_kconv_finish(tws)
	string tws
	string nf=GetDataFolder(1)
	nvar img=root:K_ARPES:global:v_kconv_img
	
	wave tw=$(tws+"_dupl")
	KillWaves tw
	
	SetDataFolder root:K_ARPES:misc:automap
	wave amp=$("amap_parameter")
	variable s=amp[0]
	variable st=amp[1]
	variable e=amp[2]
	variable n=amp[3]
	if(img!=1)
		amp[0]=s+st
		amp[1]=st
		amp[2]=e-st
		amp[3]=n-2
	endif
	string notes
	if(img==1)
		notes=K_make_note_kconv()+K_make_note()
	else
		notes=K_make_note_kconv()+K_make_note_map()
	endif
	note $(tws+"_k"),notes
	K_save_misc()
	K_restore_misc("k_saved")
	DoWindow/F K_ARPES_p
	K_refresh_panel()
	K_change_panelmode(1)
	K_show_targwave()
	SetDataFolder $nf
End

Function K_kconv_plot(tws)
	string tws
	
	wave tw=$(tws+"_dupl")
	
	string nf=GetDataFolder(1)
	nvar amap=root:K_ARPES:misc:automap:v_amapval
	nvar ken=root:K_ARPES:global:v_kconv_en
	SetDataFolder root:K_ARPES:global
	variable/g v_EF
	variable/g v_kconv_en
	variable/g v_eofs
	variable eofs=v_eofs
	variable ek=v_kconv_en+eofs
	nvar img=v_kconv_img
	
	SetDataFolder root:K_ARPES:misc
	variable/g v_th_num
	variable savevth=v_th_num
	variable/g v_th_num=DimSize(tw,1)-1
	variable/g v_EB=0
	variable es=DimOffset(tw,0)
	variable ee=DimOffset(tw,0)+DimDelta(tw,0)*(DimSize(tw,0)-1)
	
	nvar kxn=root:K_ARPES:global:v_kx_n
	
	if(kxn*2<DimSize(tw,1))
		variable/g v_th_num=kxn*2
	endif
	
	if(amap!=1)
		K_auto_mapping(0)
	endif
	
	SetDataFolder root:K_ARPES:Curves:k
	wave ckx=$("Ckx")
	wave cky=$("Cky")
	wave ckz=$("Ckz")
	
	SetDataFolder root:K_ARPES:misc
//	if(v_kconv_vol==1)
		if(!Exists("krange"))
			make/O/N=(2,3) $("krange")
		endif
		wave krw=$("krange")
		
		variable kxmax1,kxmin1,kymax1,kymin1,kzmax1,kzmin1
		variable kxmax2,kxmin2,kymax2,kymin2,kzmax2,kzmin2
		
		
		if(img==2)
			kxmax1=WaveMax(ckx)*sqrt(ek)
			kxmin1=WaveMin(ckx)*sqrt(ek)
			kymax1=WaveMax(cky)*sqrt(ek)
			kymin1=WaveMin(cky)*sqrt(ek)
			kzmax1=WaveMax(ckz)*sqrt(ek)
			kzmin1=WaveMin(ckz)*sqrt(ek)
			
			krw[0][0]=kxmin1
			krw[1][0]=kxmax1
			
			krw[0][1]=kymin1
			krw[1][1]=kymax1
			
			krw[0][2]=kzmin1
			krw[1][2]=kzmax1
			
		else
			kxmax1=WaveMax(ckx)*sqrt(eofs+ee)
			kxmin1=WaveMin(ckx)*sqrt(eofs+ee)
			kymax1=WaveMax(cky)*sqrt(eofs+ee)
			kymin1=WaveMin(cky)*sqrt(eofs+ee)
			kzmax1=WaveMax(ckz)*sqrt(eofs+ee)
			kzmin1=WaveMin(ckz)*sqrt(eofs+ee)
			
			kxmax2=WaveMax(ckx)*sqrt(eofs+es)
			kxmin2=WaveMin(ckx)*sqrt(eofs+es)
			kymax2=WaveMax(cky)*sqrt(eofs+es)
			kymin2=WaveMin(cky)*sqrt(eofs+es)
			kzmax2=WaveMax(ckz)*sqrt(eofs+es)
			kzmin2=WaveMin(ckz)*sqrt(eofs+es)
			
			variable f
			
			f=0
			
			krw[0][1]=Min(kymin1,kymin2)
			krw[1][1]=Max(kymax1,kymax2)
			
			krw[0][2]=Min(kzmin1,kzmin2)
			krw[1][2]=Max(kzmax1,kzmax2)
			
			if(NumType(kxmax1)!=0)
				krw[1][0]=kxmax2
				f+=1
			endif
			
			if(NumType(kxmax2)!=0)
				krw[1][0]=kxmax1
				f+=1
			endif
			
			if(NumType(kxmin1)!=0)
				krw[0][0]=kxmin2
				f+=2
			endif
			
			if(NumType(kxmin2)!=0)
				krw[0][0]=kxmin1
				f+=2
			endif
			
			if(f<2)
				krw[0][0]=Min(kxmin1,kxmin2)
			endif
			
			if(mod(f,2)==0)
				krw[1][0]=Max(kxmax1,kxmax2)
			endif
			
			f=0
			
			if(NumType(kymax1)!=0)
				krw[1][1]=kymax2
				f+=1
			endif
			
			if(NumType(kymax2)!=0)
				krw[1][1]=kymax1
				f+=1
			endif
			
			if(NumType(kymin1)!=0)
				krw[0][1]=kymin2
				f+=2
			endif
			
			if(NumType(kymin2)!=0)
				krw[0][1]=kymin1
				f+=2
			endif
			
			if(f<2)
				krw[0][1]=Min(kymin1,kymin2)
			endif
			
			if(mod(f,2)==0)
				krw[1][1]=Max(kymax1,kymax2)
			endif
			
			f=0
			
			if(NumType(kzmax1)!=0)
				krw[1][2]=kzmax2
				f+=1
			endif
			
			if(NumType(kzmax2)!=0)
				krw[1][2]=kzmax1
				f+=1
			endif
			
			if(NumType(kzmin1)!=0)
				krw[0][2]=kzmin2
				f+=2
			endif
			
			if(NumType(kzmin2)!=0)
				krw[0][2]=kzmin1
				f+=2
			endif
			
			if(f<2)
				krw[0][2]=Min(kzmin1,kzmin2)
			endif
			
			if(mod(f,2)==0)
				krw[1][2]=Max(kzmax1,kzmax2)
			endif
			
		endif
//	endif
	
	variable/g v_th_num=savevth
	SetDataFolder $nf
End

Function K_kconv_calc_img(tws)
	string tws
	wave tw=$(tws+"_dupl")
	string nf=GetDataFolder(1)
	SetDataFolder root:K_ARPES:global
	variable/g v_prog1=0
	NVAR prog1=v_prog1
	variable prog
	nvar k_n=v_kx_n
	variable/g v_EF
	variable ef=v_EF
	nvar efd=v_EF_data
	variable/g v_eofs
	variable eofs=v_eofs
	
	nvar outr=v_outr
	
	variable e_s=DimOffset(tw,0)
	variable e_st=DimDelta(tw,0)
	variable e_n=DimSize(tw,0)
	variable e_e=e_s+e_st*(e_n-1)
	
	SetDataFolder root:K_ARPES:Curves:k
	wave ckx=$("C0kx")
	Duplicate ckx,$("C0kx_e")
	wave ckxe=$("C0kx_e")
	
	wave cky=$("C0ky")
	Duplicate cky,$("C0ky_e")
	wave ckye=$("C0ky_e")
	
	variable kmax1,kmin1,kmax2,kmin2
	
	kmax1=WaveMax(ckx)*sqrt(e_e+eofs)
	kmax2=WaveMax(ckx)*sqrt(e_s+eofs)
	kmin1=WaveMin(ckx)*sqrt(e_e+eofs)
	kmin2=WaveMin(ckx)*sqrt(e_s+eofs)
	
	
	variable k_s
	variable k_e
	variable f=0
	
	if(NumType(kmax1)!=0)
		k_e=kmax2
		f+=1
	endif
	
	if(NumType(kmax2)!=0)
		k_e=kmax1
		f+=1
	endif
	
	if(NumType(kmin1)!=0)
		k_s=kmin2
		f+=2
	endif
	
	if(NumType(kmin2)!=0)
		k_s=kmin1
		f+=2
	endif
	
	if(f<2)
		k_s=Min(kmin1,kmin2)
	endif
	
	if(mod(f,2)==0)
		k_e=Max(kmax1,kmax2)
	endif
	
	variable k_st=(k_e-k_s)/(k_n-1)
	
	string kws=tws+"_k"
	
	make/O/N=(e_n,k_n) $kws
	wave kw=$kws
	kw=outr
	
	SetDataFolder root:K_ARPES:global
	string vwys="viewer_y"
	make/O/N=(1) $vwys
	wave vwy=$vwys
	
	string vwxs="viewer_x"
	make/O/N=(1) $vwxs
	wave vwx=$vwxs
	variable vwc=0
	
	K_Progress_panel_v()
	
//	ValDisplay prog win=K_Prog_p, highColor= (65535,0,0)
	DoUpdate /W=K_Prog_p
	
	variable stt=datetime
	variable dtt=datetime
	string nt=Secs2Time(stt,3)
	variable e,k,ek,e_i,k_i,i,kximin,dk,dkmin
	variable kx,ky,a,a_i,c,ch=1
	for(e=e_s;e<e_e+e_st/2;e+=e_st)
		ek=e+eofs
		if(ek>0)
			ckxe=ckx*sqrt(ek)
			ckye=cky*sqrt(ek)
			e_i=round((e-e_s)/e_st)
			for(k=k_s;k<k_e+k_st/2;k+=k_st)
				k_i=round((k-k_s)/k_st)
				dkmin=inf
				for(i=0;i<DimSize(ckxe,0);i+=1)
					dk=(ckxe[i]-k)^2
					if(dk<dkmin)
						dkmin=dk
						kximin=i
					endif
				endfor
				
				kx=ckxe[kximin]
				ky=ckye[kximin]
				a=K_kconv_k2ang(ek,kx,ky,0,1)
				a_i=round((a-DimOffset(tw,1))/DimDelta(tw,1))
				
				
				if(a_i>0 && a_i<DimSize(tw,1))
					kw[e_i][k_i]=tw[e_i][a_i]
				endif
				
				c+=1
				if(!StringMatch(secs2time(datetime,3),nt))
					if(ch==1)
						prog=c/(e_n*k_n)*100
					else
						prog=round(c/(e_n*k_n)*200)/2
					endif
					if(prog>prog1)
						SetDataFolder root:K_ARPES:global
						string/g s_fintime=Secs2Time(stt+(datetime-stt)*100/prog,3)
						if(ch==1)
							SetVariable setvar_fin disable=0,win=K_Prog_p
							string/g s_fintime=s_fintime+" (rough estimate)"
						endif
						if(!prog<0.5)
							if(ch==1)
								vwy[vwc]=(prog-prog1)/(datetime-dtt)
								vwx[vwc]=0
								vwc+=1
							endif
							InsertPoints vwc,1,vwx
							InsertPoints vwc,1,vwy
							vwy[vwc]=(prog-prog1)/(datetime-dtt)
							vwx[vwc]=prog
							vwc+=1
							dtt=datetime
							DoUpdate /W=K_Prog_p
						endif
						prog1=prog
						ch=0
					endif
					nt=secs2time(datetime,3)
				endif
			endfor
		endif
	endfor
	
//	ValDisplay prog win=K_Prog_p, highColor= (2,39321,1)
	
	KillWindow K_Prog_p
	
	SetScale/P x,e_s-efd,e_st,"E - EF (eV)",kw
	SetScale/P y,k_s,k_st,"kx (/A)",kw
	
	SetDataFolder $nf
End

Function K_kconv_calc(tws,vol,e,k1_n,k2_n,krange_f)
	string tws
	variable vol,e,k1_n,k2_n,krange_f // krange_f 0: use panel value // 1: use krange wave
	
	wave tw=$(tws+"_dupl")
	
	string nws=tws+"_k"
	
	print "k conv start "+secs2date(datetime,0)+" "+secs2time(datetime,3)
	K_set_status(1,"")
	
	string nf=GetDataFolder(1)
	SetDataFolder root:K_ARPES:global
	
	nvar img=v_kconv_img
	variable/g v_prog1=0
	NVAR prog1=v_prog1
	
	variable/g v_EF
	variable ef=v_EF
	nvar efd=v_EF_data
	variable/g v_eofs
	variable eofs=v_eofs
	variable/g v_kdim2=0
	variable kdim2=v_kdim2
	
	variable/g v_e_s
	variable/g v_e_st
	variable/g v_e_e
	variable/g v_e_n
	
	nvar outr=v_outr
	
	variable e_s=v_e_s
	variable e_st=v_e_st
	variable e_e=v_e_e
	variable e_n=round((v_e_e-v_e_s)/v_e_st)+1
	variable t1_s=DimOffset(tw,1)
	variable t1_st=DimDelta(tw,1)
	variable t1_n=DimSize(tw,1)
	
	variable t1_e=t1_s+(t1_n-1)*t1_st
	
	variable/g v_th_s=t1_s
	variable/g v_th_e=t1_e
	
	SetDataFolder root:K_ARPES:misc:automap
	variable/g v_amapval
	variable amap=v_amapval
	
	wave apw=$("amap_parameter")
	variable t2_st=apw[1]
	variable t2_s=apw[0]//-t2_st
	variable t2_e=apw[2]//+t2_st
	variable t2_n=apw[3]//+2
	
	variable k1_st,k2_st,e_i,k1_i,k2_i,t1_i,t2_i
	variable k1_s,k1_e,k2_s,k2_e
	
	variable k1,k2,ek,ax,ay,meb,eb
	if(amap!=1)
		SetDataFolder root:K_ARPES:global
		string vwys="viewer_y"
		make/O/N=(1) $vwys
		wave vwy=$vwys
		
		string vwxs="viewer_x"
		make/O/N=(1) $vwxs
		wave vwx=$vwxs
		variable vwc=0
		
		K_Progress_panel_v()
		
//		ValDisplay prog win=K_Prog_p, highColor= (65535,0,0)
		DoUpdate /W=K_Prog_p
		
		SetDataFolder root:K_ARPES:misc:
		wave krw=$("krange")
		if(kdim2==0)
			k1_s=krw[0][0]
			k1_e=krw[1][0]
			
			k2_s=krw[0][1]
			k2_e=krw[1][1]
		endif
		if(amap==6)
			k1_s=krw[0][0]
			k1_e=krw[1][0]
			
			k2_s=krw[0][2]
			k2_e=krw[1][2]
		endif
		
		k1_st=(k1_e-k1_s)/(k1_n-1)
		k2_st=(k2_e-k2_s)/(k2_n-1)
		
		if(vol==0)
			e_s=e
			e_e=e
			e_st=1
			e_n=1
			
			SetDataFolder $nf
			make/O/N=(k1_n,k2_n) $nws
			wave nw=$nws
		else
			SetDataFolder $nf
			make/O/N=(e_n,k1_n,k2_n) $nws
			wave nw=$nws
		endif
		
		variable dtt=datetime
		variable stt=datetime
		nw=outr
		
//		make/O/N=(k1_n,k2_n,2) $("ang_table")
//		wave at=$("ang_table")
//		at=0
//		make/O/N=(k1_n,k2_n,2) $("ang_table_e")
//		wave ate=$("ang_table_e")
//		ate=0
		
		variable in,c,prog,ch=1,e_ni
		string nt=Secs2Time(datetime,3)
		string ck1s
		string ck2s
		if(kdim2==0)
			ck1s="Ckx"
			ck2s="Cky"
		elseif(kdim2==0)
			ck1s="Ckx"
			ck2s="Ckz"
		endif
		for(e=e_s;e<e_e+e_st/2;e+=e_st)
			ek=e+eofs
			if(amap==6)
				eb=e*-1
				ek=eb
			endif
			if(img==2)
				e_i=0
			else
				e_i=round((e-DimOffset(tw,0))/DimDelta(tw,0))
			endif
			SetDataFolder root:K_ARPES:Curves:k
			
			
			wave ck1=$ck1s
			wave ck2=$ck2s
			Duplicate/O ck1,$(ck1s+"_e")
			Duplicate/O ck2,$(ck2s+"_e")
			wave ck1e=$(ck1s+"_e")
			wave ck2e=$(ck2s+"_e")
			ck1e=ck1*sqrt(ek)
			ck2e=ck2*sqrt(ek)
			
			k2_i=0
			for(k2=k2_s;k2<k2_e+k2_st/2;k2+=k2_st)
				k1_i=0
				for(k1=k1_s;k1<k1_e+k1_st/2;k1+=k1_st)
					ax=K_kconv_k2ang(ek,k1,k2,0,amap)
					ay=K_kconv_k2ang(ek,k1,k2,1,amap)
					//print k1,k2,ax,ay
					t1_i=round((ax-t1_s)/t1_st)
					t2_i=round((ay-t2_s)/t2_st)
					if(t1_i >= 0 && t1_i < t1_n && t2_i >= 0 && t2_i < t2_n)
						if(vol==0)
							nw[k1_i][k2_i]=tw[e_i][t1_i][t2_i]
						else
							nw[e_ni][k1_i][k2_i]=tw[e_i][t1_i][t2_i]
						endif
					else
						if(vol==0)
							nw[k1_i][k2_i]=outr
						else
							nw[e_ni][k1_i][k2_i]=outr
						endif
					endif
					
					c+=1
					if(!StringMatch(secs2time(datetime,3),nt))
						if(ch==1)
							prog=c/(e_n*k1_n*k2_n)*100
						else
							prog=round(c/(e_n*k1_n*k2_n)*200)/2
						endif
						if(prog>prog1)
							SetDataFolder root:K_ARPES:global
							string/g s_fintime=Secs2Date(stt+(datetime-stt)*100/prog,0)+" "+Secs2Time(stt+(datetime-stt)*100/prog,3)
							if(ch==1)
								SetVariable setvar_fin disable=0,win=K_prog_p
								string/g s_fintime=s_fintime+" (rough estimate)"
							endif
							if(!prog<0.5)
								if(ch==1)
									vwy[vwc]=(prog-prog1)/(datetime-dtt)
									vwx[vwc]=0
									vwc+=1
								endif
								InsertPoints vwc,1,vwx
								InsertPoints vwc,1,vwy
								vwy[vwc]=(prog-prog1)/(datetime-dtt)
								vwx[vwc]=prog
								vwc+=1
								dtt=datetime
								DoUpdate /W=K_Prog_p
							endif
							prog1=prog
							ch=0
						endif
						nt=secs2time(datetime,3)
					endif
					
					k1_i+=1
				endfor
				k2_i+=1
			endfor
			e_ni+=1
		endfor
		
//		KillWaves at,ate
		
//		ValDisplay prog1 win=K_Prog_p, highColor= (2,39321,1)
		
		KillWindow K_Prog_p
		
		string k1str,k2str
		if(kdim2==0)
			k1str="kx"
			k2str="ky"
		elseif(kdim2==1)
			k1str="kx"
			k2str="kz"
		endif
		
		if(vol==0)
			SetScale/P x,k1_s,k1_st,k1str+" (/A)",nw
			SetScale/P y,k2_s,k2_st,k2str+" (/A)",nw
		else
			SetScale/P x,e_s-efd,e_st,"E - EF (eV)",nw
			SetScale/P y,k1_s,k1_st,k1str+" (/A)",nw
			SetScale/P z,k2_s,k2_st,k2str+" (/A)",nw
		endif
		print "Finish "+secs2date(datetime,0)+" "+secs2time(datetime,3)
		K_set_status(0,"")
	else
		print "Please set mapping parameter."
		K_set_status(3,"")
	endif
	
	SetDataFolder $nf
End

Function K_sim_ARPES()
	string kws=GetBrowserSelection(0)
	string nf=GetDataFolder(1)
	
	wave kw=$kws
	
	variable e_s=DimOffset(kw,0)
	variable e_st=DimDelta(kw,0)
	variable e_n=DimSize(kw,0)
	variable kx_s=DimOffset(kw,1)
	variable kx_st=DimDelta(kw,1)
	variable kx_n=DimSize(kw,1)
	variable ky_s=DimOffset(kw,2)
	variable ky_st=DimDelta(kw,2)
	variable ky_n=DimSize(kw,2)
	
	variable e_e=e_s+(e_n-1)*e_st
	
	SetDataFolder root:K_ARPES:misc
	string/g s_winname
	string wn=s_winname
	variable/g v_th_s
	variable/g v_th_e
	variable/g v_th_num
	
	variable ax_s=v_th_s
	variable ax_e=v_th_e
	variable ax_n=v_th_num+1
	variable ax_st=(ax_e-ax_s)/(ax_n-1)
	
	variable/g v_hn
	variable/g v_W
	variable ef=v_hn-v_W
	
	SetDataFolder root:K_ARPES:misc:automap
	wave apw=$("amap_parameter")
	variable ay_s=apw[0]
	variable ay_e=apw[2]
	variable ay_n=apw[3]
	variable ay_st=apw[1]
	
	SetDataFolder $nf
	make/O/N=(e_n,ax_n,ay_n) $(kws+"_ang")
	wave nw=$(kws+"_ang")
	nw=0

	SetDataFolder root:K_ARPES:Curves:$wn
	wave ckx=$("Ckx")
	wave cky=$("Cky")
	Duplicate/O ckx,$("Ckx_e")
	Duplicate/O cky,$("Cky_e")
	wave ckxe=$("Ckx_e")
	wave ckye=$("Cky_e")
	variable e,ax,ay,kx,ky
	for(e=e_s;e<e_e+e_st/2;e+=e_st)
		ckxe=ckx*sqrt(e+ef)
		ckye=cky*sqrt(e+ef)
		for(ay=0;ay<ay_n;ay+=1)
			for(ax=0;ax<ax_n;ax+=1)
				kx=round((ckxe[ax][ay]-kx_s)/kx_st)
				ky=round((ckye[ax][ay]-ky_s)/ky_st)
				if(kx>=0 && kx<kx_n && ky>=0 && ky<ky_n)
					nw[round((e-e_s)/e_st)][ax][ay]=kw[round((e-e_s)/e_st)][kx][ky]
				else
					nw[round((e-e_s)/e_st)][ax][ay]=NaN
				endif
			endfor
		endfor
	endfor
	
	SetScale/P x,e_s,e_st,nw
	SetScale/P y,ax_s,ax_st,"ax",nw
	SetScale/P z,ay_s,ay_st,"ay",nw
	
	SetDataFolder $nf
End

Function K_kconv_k2ang(ek,k1,k2,adim,amap)
	variable ek,k1,k2,adim,amap // adim 0: ax, 1: ay
	
	string nf=GetDataFolder(1)
	variable ax,ay
	
	SetDataFolder root:K_ARPES:global
	variable/g v_kdim2
	variable kdim2=v_kdim2
	
	SetDataFolder root:K_ARPES:misc
	variable/g v_smh
	variable smh=v_smh
	variable/g v_th_s
	variable/g v_th_e
	variable th_s=v_th_s
	variable th_e=v_th_e
	variable/g v_bet
	variable bet=v_bet
	variable/g v_V0
	variable V0=v_V0
	variable/g v_W
	variable W=v_W
	
	v0=0 // 2020.07.01
	
	variable kx,ky,kz,krength2,eb,hn
	if(kdim2==0)
		kx=k1
		ky=k2
		kz=sqrt((ek+V0)*smh^2-k1^2-k2^2)
	elseif(kdim2==1)
//		kx=k1
//		ky=sqrt(ek*smh^2-k1^2-k2^2)
//		kz=k2
	endif
	
	if(amap==6)
		kx=k1
		ky=0
		kz=k2
		krength2=(kx^2+kz^2)/smh^2
		eb=ek
		hn=krength2-V0+W+eb
	endif
	
	
	SetDataFolder root:K_ARPES:misc:rot_matrix
	string mmxs="man_rot"
	wave mmx=$mmxs
	wave mmx_1=$("R1")
	wave mmx_1i=$("R1_inv")
	wave mmx_2=$("R2")
	wave mmx_2i=$("R2_inv")
	wave rm=$("rotation_matrix")
	string kvws="kvec"
	make/O/N=3 $kvws
	wave kvw=$kvws
	
	kvw[0]=kx
	kvw[1]=ky
	kvw[2]=kz
	
	SetDataFolder root:K_ARPES:Curves:k
	if(amap!=1)
		wave cxw=$("Ckx_e")
	endif
	
	SetDataFolder root:K_ARPES:misc:rot_matrix
	
	variable tt=0
	
	variable i,d,dmin,alp,alpmin
	if(amap==1) // no mapping
		MatrixMultiply mmx,kvw
		wave nkvw=$("M_product")
		ax=K_cord_k2a(nkvw[0],nkvw[1],nkvw[2],0)*180/pi
	elseif(amap==5) // beta map DA30
		MatrixMultiply mmx,kvw
		wave nkvw=$("M_product")
		if(adim==0)
			ax=K_cord_k2a(nkvw[0],nkvw[1],nkvw[2],0)*180/pi
		else
			ay=K_cord_k2a(nkvw[0],nkvw[1],nkvw[2],1)*180/pi
		endif
	elseif(amap==6)
		MatrixMultiply mmx,kvw
		wave nkvw=$("M_product")
		if(abs(nkvw[1])!=0)
			print "warning"
		endif
		if(adim==0)
			ax=K_cord_k2a(nkvw[0],nkvw[1],nkvw[2],0)*180/pi
		else
			ay=hn
		endif
	elseif(amap==3 && tt==1)
//		MatrixMultiply mmx_1,kvw
//		wave nkvw=$("M_product")
//		
//		dmin=inf
//		for(i=0;i<DimSize(cxw,0);i+=1)
//			alp=(th_s+i*(th_e-th_s)/(DimSize(cxw,0)-1))*pi/180
//			bet=bet*pi/180
//			d=((mmx_2i[0][0]*alp+mmx_2i[0][1]*bet)/sqrt(alp^2+bet^2)*sin(sqrt(alp^2+bet^2))+mmx_2i[0][2]*cos(sqrt(alp^2+bet^2))-nkvw[0])^2
//			if(d<dmin)
//				dmin=d
//				alpmin=alp
//			endif
//		endfor
//		
//		ax=alpmin*180/pi
//		if(adim==0)
//			print dmin
//		endif
//		if(adim==1)
//			K_q_run(nkvw[0],nkvw[1],nkvw[2],alpmin/sqrt(alpmin^2+bet^2)*sin(sqrt(alpmin^2+bet^2)),bet/sqrt(alpmin^2+bet^2)*sin(sqrt(alpmin^2+bet^2)),cos(sqrt(alpmin^2+bet^2)))
//			wave qw=$("q_mw")
//			
//			MatrixMultiply mmx_2i,qw,mmx_1i
//			Duplicate/O $("M_product"),rm
//			KillWaves $("M_product")
//			
//			ay=acos(rm[1][1])*180/pi
//		endif
		
	else
		SetDataFolder root:K_ARPES:misc:automap
		wave awp=$("amap_parameter")
		SetDataFolder root:K_ARPES:Curves:k
		if(kdim2==0)
			wave c1w=$("Ckx_e")
			wave c2w=$("Cky_e")
		elseif(kdim2==1)
			wave c1w=$("Ckx_e")
			wave c2w=$("Ckz_e")
		endif
		wave cdw=$("Cdim")
		wave alw=$("amap_log")
		
		MultiThread cdw=(c1w-k1)^2+(c2w-k2)^2  //  GEKIOSO
		FindValue/V=(WaveMin(cdw)) cdw // OSOI
		
		variable minrow,mincol
		
		//iv8
//		minrow=V_row
//		mincol=V_col
		
		//iv7
		mincol=Floor(V_value/DimSize(cdw,0))
		minrow=V_value-mincol*DimSize(cdw,0)
		
		ax=th_s+minrow*(th_e-th_s)/(DimSize(cdw,0)-1)
		ay=awp[0]+mincol*awp[1]
		
	endif
	
	SetDataFolder $nf
	
	if(adim==0)
		return ax
	elseif(adim==1)
		return ay
	endif
End

Function K_SetVarProc_thsnum(sva) : SetVariableControl
	STRUCT WMSetVariableAction &sva

	switch( sva.eventCode )
		case 1: // mouse up
		case 2: // Enter key
		case 3: // Live update
			Variable dval = sva.dval
			String sval = sva.sval
			string nf=GetDataFolder(1)
			SetDataFolder root:K_ARPES:misc:automap
			variable/g v_amapval
			variable amap=v_amapval
			SetDataFolder root:K_ARPES:misc
			
			variable/g v_th_s,v_th_e
			nvar v_kconv=root:K_ARPES:global:v_kconv
			
			if(dval>v_th_e)
				variable/g v_th_e=v_th_s
			endif
			
			execute "K_set_thxnum("+num2str(v_th_s)+","+num2str(v_th_e)+")"
			
			K_set_emis_ang_window()
			
			string/g s_winname
			variable/g v_app
			
			if(!WinType(s_winname)==0 && v_kconv==0)
				if(amap==1)
					if(v_app==0)
						K_ARPES_calc_kline()
						K_append_curve()
						K_change_mode()
					endif
				else
					K_auto_mapping(1)
				endif
			endif
			
			SetDataFolder $nf
			break
		case -1: // control being killed
			break
	endswitch

	return 0
End

Function K_SetVarProc_thenum(sva) : SetVariableControl
	STRUCT WMSetVariableAction &sva

	switch( sva.eventCode )
		case 1: // mouse up
		case 2: // Enter key
		case 3: // Live update
			Variable dval = sva.dval
			String sval = sva.sval
			string nf=GetDataFolder(1)
			SetDataFolder root:K_ARPES:misc:automap
			variable/g v_amapval
			variable amap=v_amapval
			SetDataFolder root:K_ARPES:misc
			
			variable/g v_th_s,v_th_e
			nvar v_kconv=root:K_ARPES:global:v_kconv
			
			if(dval<v_th_s)
				variable/g v_th_s=v_th_e
			endif
			
			execute "K_set_thxnum("+num2str(v_th_s)+","+num2str(v_th_e)+")"
			
			K_set_emis_ang_window()
			
			string/g s_winname
			variable/g v_app
			
			if(!WinType(s_winname)==0 && v_kconv==0)
				if(amap==1)
					if(v_app==0)
						K_ARPES_calc_kline()
						K_append_curve()
						K_change_mode()
					endif
				else
					K_auto_mapping(1)
				endif
			endif
			
			SetDataFolder $nf
			break
		case -1: // control being killed
			break
	endswitch

	return 0
End



Function K_Progress_panel()
	PauseUpdate; Silent 1		// building window...
	NewPanel /N=K_Prog_p /K=1 /W=(41,78,317,133) as "Calculation Progress"
	ValDisplay prog,pos={33.00,17.00},size={200.00,17.00},title="Progress:"
	ValDisplay prog,limits={0,100,0},barmisc={0,0},mode= 3
	ValDisplay prog,value= #"root:K_ARPES:global:v_prog1"
	SetVariable setvar_prog,pos={142.00,15.00},size={50.00,18.00},title=" "
	SetVariable setvar_prog,format="%.1f %",frame=0
	SetVariable setvar_prog,limits={-inf,inf,0},value= root:K_ARPES:global:v_prog1,noedit= 1
	SetVariable setvar_fin,pos={19.00,32.00},size={220.00,18.00},bodyWidth=150,title="Finish time : "
	SetVariable setvar_fin,frame=0
	SetVariable setvar_fin,limits={-inf,inf,0},value= root:K_ARPES:global:s_fintime,noedit= 1
End

Function K_Progress_panel_v()
	PauseUpdate; Silent 1		// building window...
	NewPanel /N=K_Prog_p /K=1 /W=(110,77,386,184) as "Calculation Progress"
	SetVariable setvar_prog,pos={4.00,81.00},size={50.00,18.00},title=" "
	SetVariable setvar_prog,format="%.1f %",frame=0
	SetVariable setvar_prog,limits={-inf,inf,0},value= root:K_ARPES:global:v_prog1,noedit= 1
	SetVariable setvar_fin,pos={52.00,81.00},size={220.00,18.00},bodyWidth=150,title="Finish time : "
	SetVariable setvar_fin,frame=0
	SetVariable setvar_fin,limits={-inf,inf,0},value= root:K_ARPES:global:s_fintime,noedit= 1
	SetVariable setvar_fin disable=1
	String fldrSav0= GetDataFolder(1)
	SetDataFolder root:K_ARPES:global
	Display/W=(13,10,268,81)/HOST=#  $("viewer_y") vs $("viewer_x")
	SetDataFolder fldrSav0
	ModifyGraph margin(left)=14,margin(bottom)=28,margin(top)=5,margin(right)=11,wbRGB=(61166,61166,61166)
	ModifyGraph gbRGB=(61166,61166,61166)
	ModifyGraph mode=7
	ModifyGraph rgb=(1,52428,26586)
	ModifyGraph hbFill=2
	ModifyGraph tick(bottom)=3
	ModifyGraph mirror=2
	ModifyGraph nticks(left)=0
	ModifyGraph noLabel(left)=0
	ModifyGraph standoff=1
	SetAxis left 0,*
	SetAxis bottom 0,100
	Label left "Speed"
	RenameWindow #,G0
	SetActiveSubwindow ##
End



Function K_SetVarProc_ar(sva) : SetVariableControl
	STRUCT WMSetVariableAction &sva

	switch( sva.eventCode )
		case 1: // mouse up
		case 2: // Enter key
		case 3: // Live update
			Variable dval = sva.dval
			String sval = sva.sval
			string nf=GetDataFolder(1)
			SetDataFolder root:K_ARPES:misc:automap
			variable/g v_amapval
			variable amap=v_amapval
			SetDataFolder root:K_ARPES:misc
			string/g s_winname
			variable/g v_app
			variable/g v_map
			variable/g v_smh
			nvar v_kconv=root:K_ARPES:global:v_kconv
			
			if(dval>=0)
				variable/g v_anal_rot=mod(dval,360)
			elseif(dval<0)
				variable/g v_anal_rot=359+mod(dval+1,360)
			endif
			
			K_set_emis_ang_window()
			
			if(!WinType(s_winname)==0 && v_kconv==0)
				if(amap==1)
					if(v_app==0)
						K_ARPES_calc_kline()
						K_append_curve()
					elseif(v_map==1)
						K_ARPES_calc_kline()
						K_append_curve()
					endif
				else
					if(amap==2 && mod(dval,180)==0)
						SetDataFolder root:K_ARPES:misc:automap
						variable/g v_amapval=3
						string/g s_amapval="Tilt"
						DoWindow/F K_ARPES_p
						PopupMenu popup_amapval mode=3
						SetVariable setvar_po disable=0
						SetVariable setvar_ph disable=2
					elseif(amap==3 && mod(dval,180)==90)
						SetDataFolder root:K_ARPES:misc:automap
						variable/g v_amapval=2
						string/g s_amapval="Polar"
						DoWindow/F K_ARPES_p
						PopupMenu popup_amapval mode=2
						SetVariable setvar_po disable=2
						SetVariable setvar_ph disable=0
					endif
					if(v_app==0)
						K_auto_mapping(1)
					endif
				endif
			endif
			
			SetDataFolder $nf
			break
		case -1: // control being killed
			break
	endswitch

	return 0
End





Function K_save_misc()
	string nf=GetDataFolder(1)
	NewDataFolder/O root:K_ARPES:misc_saved
	SetDataFolder root:K_ARPES:misc
	string/g s_winname
	string wn=s_winname
	
	SetDataFolder root:K_ARPES:misc_saved
	if(DataFolderExists(wn))
		KillDataFolder root:K_ARPES:misc_saved:$wn
	endif
	DuplicateDataFolder root:K_ARPES:misc root:K_ARPES:misc_saved:$wn
	SetDataFolder $nf
End

Function K_restore_misc(wn)
	string wn
	string nf=GetDataFolder(1)
	NewDataFolder/O/S root:K_ARPES:misc_saved
	if(DataFolderExists(wn))
		SetDataFolder root:K_ARPES
		if(DataFolderExists("misc"))
			KillDataFolder root:K_ARPES:misc
		endif
		DuplicateDataFolder root:K_ARPES:misc_saved:$wn root:K_ARPES:misc
	endif
	
	SetDataFolder $nf
End

Function K_SetVarProc_restore(sva) : SetVariableControl
	STRUCT WMSetVariableAction &sva

	switch( sva.eventCode )
		case 1: // mouse up
		case 2: // Enter key
		case 3: // Live update
			Variable dval = sva.dval
			String sval = sva.sval
			
			string nf=GetDataFolder(1)
			
			SetDataFolder root:K_ARPES:global
			wave/T wlw=$("win_list")
			
			variable i=0,flag=0
			if(DimSize(wlw,0)>0)
				do
					if(StringMatch(sval,wlw[i]))
						flag=1
						break
					endif
					i+=1
				while(i<DimSize(wlw,0))
			else
				flag=1
			endif
			
			if(flag==1)
				K_restore_misc(sval)
				K_refresh_panel()
			endif
			
			SetDataFolder $nf
			
			break
		case -1: // control being killed
			break
	endswitch

	return 0
End





Function K_auto_mapping(plot)
	variable plot
	string nf=GetDataFolder(1)
	K_save_misc()
	K_set_status(1,"")
	SetDataFolder root:K_ARPES:misc:automap
	variable/g v_amapval
	variable val=v_amapval
	
	if(val>1)
		
		wave ampw=$("amap_parameter")
		variable st=ampw[1]
		variable s=ampw[0]
		variable e=ampw[2]
		
		SetDataFolder root:K_ARPES:misc
		variable/g v_app
		variable save_v_app=v_app
		variable/g v_map
		variable save_v_map=v_map
		string/g s_winname
		string winn=s_winname
		variable/g v_hn
		variable/g v_W
		variable/g v_EB
		variable ek=v_hn-v_W-v_EB
		nvar vx=v_cross_x
		nvar vy=v_cross_y
		nvar vz=v_cross_z
		SetDataFolder root:K_ARPES:global
		variable/g v_kconv
		variable kconv=v_kconv
		
		variable cns,plsiz
		string logws="amap_log"
		SetDataFolder root:K_ARPES:Curves:$winn
		if(!Exists(logws) || kconv==1)
			make/O/N=0 $logws
			wave logw=$logws
			cns=0
			plsiz=0
		else
			wave logw=$logws
			plsiz=DimSize(logw,0)
			if(save_v_app==0)
				cns=logw[0]
				make/O/N=0 $logws
			elseif(save_v_app==1)
				cns=logw[DimSize(logw,0)-1]+1
				make/O/N=0 $logws
			endif
		endif
		wave logw=$logws
		
		SetDataFolder root:K_ARPES:misc
		variable save_para
		if(save_v_app==1)
			variable/g v_curve_num_map
			variable/g v_curve_num_map=v_curve_num_map+1
		endif
		
		switch(val)
			case 2:
				variable/g v_po
				save_para=v_po
				break
			case 3:
				variable/g v_ph
				save_para=v_ph
				break
			case 4:
				variable/g v_az
				save_para=v_az
				break
			case 5:
				variable/g v_bet
				save_para=v_bet
				break
			case 6:
				variable/g v_hn
				save_para=v_hn
				break
			case 7:
				variable/g v_EB
				save_para=v_EB
				break
		endswitch
		
		variable i,c,f_app,tr
		string trl,trs
		
		c=cns
		for(i=s;i<e+st/2;i+=st)
			switch(val)
				case 2:
					variable/g v_po=i
					break
				case 3:
					variable/g v_ph=i
					break
				case 4:
					variable/g v_az=i
					break
				case 5:
					variable/g v_bet=i
					break
				case 6:
					variable/g v_hn=i
					break
				case 7:
					variable/g v_EB=i
					break
			endswitch
			
			variable/g v_map=1
			variable/g v_app=1
			
			variable/g v_curve_num=c
			
			if(WinType(winn)==1)
				if(save_v_app==1)
					f_app=1
				elseif(save_v_app==0)
					if(plsiz>c-cns)
						f_app=0
					else
						f_app=1
					endif
				else
				endif
			else
				f_app=1
			endif
			
			SetDataFolder root:K_ARPES:Curves:$winn
			InsertPoints c,1,logw
			wave logw=$logws
			SetDataFolder root:K_ARPES:misc
			variable/g v_curve_num
			logw[c-cns]=v_curve_num
			
			if(plot==0)
				f_app=0
			endif
			
			K_calc_normal()
			K_ARPES_calc_kline_calc()
			trl=TraceNameList(winn,";",1)
			tr=0
			do
				trs=StringFromList(tr,trl)
				if(!StrLen(trs)>0)
					break
				endif
				if(StringMatch(trs,"C"+num2str(c)+"*"))
					f_app=0
					break
				endif
				tr+=1
			while(1)
			if(f_app==0)
			elseif(f_app==1)
				variable/g v_flag_app=1
				K_append_curve()
			endif
			if(plot==1)
				K_change_mode()
			endif
			
			c+=1
		endfor
		
		variable kc
		if(save_v_app==0)
			for(i=c-cns;i<plsiz;i+=1)
				DoWindow/F $winn
				
				kc=i+cns
				if(vx==1)
					RemoveFromGraph $("C"+num2str(kc)+"kx")
				endif
				if(vy==1)
					RemoveFromGraph $("C"+num2str(kc)+"ky")
				endif
				if(vz==1)
					RemoveFromGraph $("C"+num2str(kc)+"kz")
				endif
				SetDataFolder root:K_ARPES:Curves:$winn
				KillWaves $("C"+num2str(kc)+"kx"),$("C"+num2str(kc)+"ky"),$("C"+num2str(kc)+"kz")
			endfor
		endif
		
		K_calcAxisRange_w()
		
		if(c>0)
			SetDataFolder root:K_ARPES:Curves:$winn
			string cxwlist="",cywlist="",czwlist=""
			for(i=0;i<DimSize(logw,0);i+=1)
				cxwlist+="C"+num2str(logw[i])+"kx;"
				cywlist+="C"+num2str(logw[i])+"ky;"
				czwlist+="C"+num2str(logw[i])+"kz;"
			endfor
			Concatenate/O cxwlist,$("Ckx")
			Concatenate/O cywlist,$("Cky")
			Concatenate/O czwlist,$("Ckz")
			
			wave ckx=$("Ckx")
			wave cky=$("Cky")
			wave ckz=$("Ckz")
			
			ckx=ckx/sqrt(ek)
			cky=cky/sqrt(ek)
			ckz=ckz/sqrt(ek)
			
			Duplicate/O ckx,Cdim
			wave dw=$("Cdim")
			dw=0
			K_set_status(0,"")
		else
			K_set_status(3,"Please set mapping parameters.")
		endif
		
		SetDataFolder root:K_ARPES:misc
		
		switch(val)
			case 2:
				variable/g v_po=save_para
				break
			case 3:
				variable/g v_ph=save_para
				break
			case 4:
				variable/g v_az=save_para
				break
			case 5:
				variable/g v_bet=save_para
				break
			case 6:
				variable/g v_hn=save_para
				break
			case 7:
				variable/g v_EB=save_para
				break
		endswitch
		
		variable/g v_map=save_v_map
		variable/g v_app=save_v_app
	endif
	
	SetDataFolder $nf
End

Function K_ButtonProc_amap(ba) : ButtonControl
	STRUCT WMButtonAction &ba

	switch( ba.eventCode )
		case 2: // mouse up
			// click code here
			K_auto_mapping(1)
			break
		case -1: // control being killed
			break
	endswitch

	return 0
End

Function K_set_amapval(val)
	variable val
	string nf=GetDataFolder(1)
	SetDataFolder root:K_ARPES:global
	variable/g v_kconv
	variable kconv=v_kconv
	SetDataFolder root:K_ARPES:misc:automap
	
	string amapList=";none;Polar;Tilt;Azimuth;Beta;hn;EB"
	
	string/g s_amapval=StringFromList(val,amapList)
	variable/g v_amapval=val
	
	if(val==1)
		DoWindow/F K_ARPES_p
		SetVariable setvar_map_st disable=2
		SetVariable setvar_map_step disable=2
		SetVariable setvar_map_en disable=2
		SetVariable setvar_map_num disable=2
		
		SetVariable setvar_alp disable=0
		SetVariable setvar_bet disable=0
		if(kconv==0)
			CheckBox check_map disable=0
		endif
		Button button_find disable=0
	else
		DoWindow/F K_ARPES_p
		SetVariable setvar_map_st disable=0
		SetVariable setvar_map_step disable=0
		SetVariable setvar_map_en disable=0
		SetVariable setvar_map_num disable=0
		
		CheckBox check_map disable=3
		SetVariable setvar_alp disable=2
		SetVariable setvar_bet disable=2
		Button button_find disable=2
		
		if(kconv==0)
			Duplicate/O $("amap_parameter_"+StringFromList(val,amapList)),$("amap_parameter")
		endif
	endif

	DoWindow/F K_ARPES_p
	SetVariable setvar_po disable=0
	SetVariable setvar_ph disable=0
	SetVariable setvar_az disable=0
	SetVariable setvar_bet_offs disable=0
	SetVariable setvar_hn disable=0
	SetVariable setvar_EB disable=0
	PopupMenu popup_amapval mode=val
	
	switch(val)
		case 2:
			SetVariable setvar_po disable=2
			break
		case 3:
			SetVariable setvar_ph disable=2
			break
		case 4:
			SetVariable setvar_az disable=2
			break
		case 5:
			SetVariable setvar_bet_offs disable=2
			break
		case 6:
			SetVariable setvar_hn disable=2
			break
		case 7:
			SetVariable setvar_EB disable=2
			break
	endswitch
	SetDataFolder $nf
End

Function K_kconv_img_CE()
	string nf=GetDataFolder(1)
	SetDataFolder root:K_ARPES:global
	string/g s_kconv_target
	string tws=s_kconv_target
	variable/g v_kx_n
	variable/g v_ky_n
	variable/g v_kconv_ek
	variable en=v_kconv_ek
	variable/g v_kconv_eofs=0
	string winn="k"
	
	if(StringMatch(tws,""))
		print "Please set ARPES data."
	else
		if(Exists(tws)==1)
			wave tw=$tws
			
			string vtws=tws+"_vol"
			Duplicate/O tw,$vtws
			Redimension/N=(-1,-1,1) $vtws
			ImageTransform/G=2 transposeVol $vtws
			Duplicate/O $("M_VolumeTranspose") $vtws
			KillWaves $("M_VolumeTranspose")
			
			K_kconv_prep(vtws)
			
			SetDataFolder root:K_ARPES:global
			if(Exists(winn))
				KillWaves $winn
			endif
			SetDataFolder $nf
			if(WinType(winn)==1)
				KillWindow $winn
			endif
			K_make_window()
			
			K_kconv_plot(vtws)
			K_kconv_calc(vtws,0,en,v_kx_n,v_ky_n,1)
			K_kconv_finish(vtws)
			
			Duplicate/O $(vtws+"_k"),$(tws+"_k")
			KillWaves $(vtws+"_k"),$vtws
		else
			print "Please check ARPES data."
		endif
	endif
	SetDataFolder $nf
End

Function K_kconv_vol()
	string nf=GetDataFolder(1)
	SetDataFolder root:K_ARPES:global
	string/g s_kconv_target
	string tws=s_kconv_target
	variable/g v_kx_n
	variable/g v_ky_n
	variable/g v_kconv_vol
	variable/g v_kconv_en
	variable en=v_kconv_en
	variable/g v_EF
	variable ef=v_EF
	string winn="k"
	
	if(StringMatch(tws,""))
		print "Please set ARPES data."
	else
		if(Exists(tws)==1)
			wave tw=$tws
			variable e_i=round((en-DimOffset(tw,0))/DimDelta(tw,0))
			if((e_i >=0 && e_i < DimSize(tw,0)) || v_kconv_vol==1)
				K_kconv_prep(tws)
				
				SetDataFolder root:K_ARPES:global
				if(Exists(winn))
					KillWaves $winn
				endif
				SetDataFolder $nf
				if(WinType(winn)==1)
					KillWindow $winn
				endif
				K_make_window()
				
				K_kconv_plot(tws)
				K_kconv_calc(tws,v_kconv_vol,en,v_kx_n,v_ky_n,1)
				K_kconv_finish(tws)
			else
			print "E - EB is out of range."
			endif
		else
			print "Please check ARPES data."
		endif
	endif
	SetDataFolder $nf
End

Function K_kconv_img_Ek()
	string nf=GetDataFolder(1)
	SetDataFolder root:K_ARPES:global
	string/g s_kconv_target
	string tws=s_kconv_target
	variable/g v_kx_n
	variable k_n=v_kx_n
	variable/g v_ky_n
	variable/g v_kconv_img
	variable img=v_kconv_img
	variable/g v_EF
	variable ef=v_EF
	string winn="k"
	string vtws=tws+"_vol"
	
	if(StringMatch(tws,""))
		print "Please set ARPES data."
	else
		if(Exists(tws)==1)
			K_kconv_prep(tws)
			SetDataFolder root:K_ARPES:misc
			variable/g v_th_num
			variable stn=v_th_num
			wave tw=$(tws+"_dupl")
			if(DimSize(tw,1)>k_n*2)
				variable/g v_th_num=k_n*2
			else
				variable/g v_th_num=DimSize(tw,1)-1
			endif
			variable/g v_curve_num=0
			K_calc_normal()
			K_ARPES_calc_kline_calc()
			variable/g v_th_num=stn
			nvar hn=root:K_ARPES:misc:v_hn
			nvar W=root:K_ARPES:misc:v_W
			nvar EB=root:K_ARPES:misc:v_EB
			variable ek=hn-W-EB
			SetDataFolder root:K_ARPES:Curves:k
			wave ckx=$("C0kx")
			wave cky=$("C0ky")
			ckx=ckx/sqrt(ek)
			cky=cky/sqrt(ek)
			K_kconv_calc_img(tws)
			K_kconv_finish(tws)
			
		else
			print "Please check ARPES data."
		endif
	endif
	SetDataFolder $nf
End


Function K_show_targwave()
	string nf=GetDataFolder(1)
	SetDataFolder root:K_ARPES:global
	string/g s_kconv_target
	string targ=s_kconv_target
	K_Set_Status(0,"")
	variable ng=1
	string targws="kconv_target_waves"
	if(Exists(targ)==1)
		wave tw=$targ
		SetDataFolder root:K_ARPES:global
		if(!WaveExists($targws))
			make/T/N=(0,11) $targws
		endif
		wave/T targw=$targws
		
		variable tex=K_ktw_search(targ)
		
		DoWindow/F K_ARPES_p
		Button button_select_wave title="\\JL"+targ
		SetVariable setvar_kconv_eofs disable=0
		SetVariable setvar_kconv_ef_data disable=0
		
		if(WaveDims(tw)==3)
			SetDataFolder root:K_ARPES:misc:automap
			wave apw=$("amap_parameter")
			apw[0]=DimOffset(tw,2)
			apw[1]=DimDelta(tw,2)
			apw[2]=DimOffset(tw,2)+DimDelta(tw,2)*(DimSize(tw,2)-1)
			apw[3]=DimSize(tw,2)
			
			K_kconv_panel_targw_eofs(targ)
			K_kconv_check_enrs_panel()
			
			SetDataFolder root:K_ARPES:global
			variable/g v_kconv_vol=0
			K_write_ktw()
			tex=K_ktw_search(targ)
			if(tex<0)
				variable/g v_e_s=DimOffset(tw,0)
				variable/g v_e_st=DimDelta(tw,0)
				variable/g v_e_n=DimSize(tw,0)
				variable/g v_e_e=v_e_s+v_e_st*(v_e_n-1)
			else
				variable/g v_e_s=str2num(targw[tex][5])
				variable/g v_e_st=str2num(targw[tex][6])
				variable/g v_e_n=DimSize(tw,0)
				variable/g v_e_e=v_e_s+v_e_st*(v_e_n-1)
			endif
			
			variable/g v_kconv_img=0
			
			DoWindow/F K_ARPES_p
			CheckBox check_kconv_vol disable=0
			SetVariable setvar_kconv_en disable=0
			SetVariable setvar_kconv_en_s disable=3
			SetVariable setvar_kconv_en_e disable=3
			SetVariable setvar_kconv_en_str,disable=0
			SetVariable setvar_kconv_kxn disable=0,title="kx points :"
			SetVariable setvar_kconv_kyn disable=0
			SetVariable setvar_kconv_ef disable=0
			SetVariable setvar_kconv_ek disable=3
			PopupMenu popup_ydim disable=3
			ng=0
		elseif(WaveDims(tw)==2)
			SetDataFolder root:K_ARPES:misc:automap
			wave apw=$("amap_parameter")
			apw[0]=DimOffset(tw,1)
			apw[1]=DimDelta(tw,1)
			apw[2]=DimOffset(tw,1)+DimDelta(tw,1)*(DimSize(tw,1)-1)
			apw[3]=DimSize(tw,1)
			SetDataFolder root:K_ARPES:global
			variable/g v_kdim2=0
			DoWindow/F K_ARPES_p
			CheckBox check_kconv_vol disable=3
			SetVariable setvar_kconv_en disable=3
			SetVariable setvar_kconv_en_str,disable=3
			SetVariable setvar_kconv_kxn disable=0,title="k points :"
			SetVariable setvar_kconv_kyn disable=3
			SetVariable setvar_kconv_ef disable=0
			SetVariable setvar_kconv_ek disable=3
			PopupMenu popup_ydim disable=0,mode=1
			Button button_go disable=0
			K_set_amapval(1)
			SetVariable setvar_EB disable=2
			CheckBox check_map disable=3
			K_kconv_panel_targw_eofs(targ)
			SetDataFolder root:K_ARPES:global
			K_write_ktw()
			tex=K_ktw_search(targ)
			wave/T ktw=$("kconv_target_waves")
			K_kconv_2D_img(str2num(ktw[tex][10]))
			ng=0
		endif
		SetDataFolder root:K_ARPES:global
		K_write_ktw()
		tex=K_ktw_search(targ)
		variable/g v_kx_n=str2num(targw[tex][8])
		variable/g v_ky_n=str2num(targw[tex][9])
	endif
	if(ng==1)
		DoWindow/F K_ARPES_p
		Button button_select_wave title=""
		string/g s_kconv_target=""
		CheckBox check_kconv_vol disable=3
		SetVariable setvar_kconv_en disable=3
		SetVariable setvar_kconv_en_str,disable=3
		SetVariable setvar_kconv_kxn disable=3
		SetVariable setvar_kconv_kyn disable=3
		SetVariable setvar_kconv_ef disable=3
		PopupMenu popup_ydim disable=3
		SetVariable setvar_kconv_eofs disable=3
		SetVariable setvar_kconv_ef_data disable=3
	endif
	SetDataFolder $nf
End

Function K_kconv_panel_targw_eofs(tws)
	string tws // abs
	wave tw=$tws
	string nf=GetDataFolder(1)
	SetDataFolder root:K_ARPES:global
	wave/t ktw=$("kconv_target_waves")
	variable tex=K_ktw_search(tws)
	
	if(tex<0)
		if(DimOffset(tw,0)<0)
			variable/g v_f_ef=0
			nvar hn=root:K_ARPES:misc:v_hn
			nvar W=root:K_ARPES:misc:v_W
			variable/g v_EF=hn-W
			variable/g v_eofs=v_EF
			variable/g v_EF_data=0
			variable/g v_kconv_en=0
		else
			variable/g v_f_ef=1
			nvar hn=root:K_ARPES:misc:v_hn
			nvar W=root:K_ARPES:misc:v_w
			variable/g v_EF=hn-W
			variable/g v_eofs=0
			variable/g v_EF_data=v_EF
			variable/g v_kconv_en=v_EF_data
		endif
	else
			variable/g v_EF=str2num(ktw[tex][1])
			variable/g v_eofs=str2num(ktw[tex][4])
			variable/g v_EF_data=str2num(ktw[tex][2])
			variable/g v_kconv_en=str2num(ktw[tex][7])
	endif
	SetDataFolder $nf
End

Function K_kconv_check_enrs_panel()
	string nf=GetDataFolder(1)
	SetDataFolder root:K_ARPES:global
	variable/g v_kconv_en
	variable/g v_e_e
	variable/g v_e_s
	variable ng=0
	
	DoWindow/F K_ARPES_p
	if(K_kconv_check_enr(v_kconv_en)==2)
		SetVariable setvar_kconv_en valueBackColor=(65535,32768,32768)
		ng=1
	else
		SetVariable setvar_kconv_en valueBackColor=0
	endif
	
	if(K_kconv_check_enr(v_e_e)==2)
		SetVariable setvar_kconv_en_e valueBackColor=(65535,32768,32768)
		ng=1
	else
		SetVariable setvar_kconv_en_e valueBackColor=0
	endif
	
	if(K_kconv_check_enr(v_e_s)==2)
		SetVariable setvar_kconv_en_s valueBackColor=(65535,32768,32768)
		ng=1
	else
		SetVariable setvar_kconv_en_s valueBackColor=0
	endif
	
	if(ng==1)
		Button button_go disable=2
	else
		Button button_go disable=0
	endif
	
	SetDataFolder $nf
End

Function K_kconv_check_enr(val)
	variable val
	string nf=GetDataFolder(1)
	SetDataFolder root:K_ARPES:global
	string/g s_kconv_target
	string targ=s_kconv_target
	
	variable out=0
	if(Exists(targ))
		wave tw=$targ
		
		variable e_s=DimOffset(tw,0)
		variable e_e=DimOffset(tw,0)+DimDelta(tw,0)*(DimSize(tw,0)-1)
		
		if(val>=e_s && val<=e_e)
			out=1
		else
			out=2
		endif
	endif
	
	SetDataFolder $nf
	return out
End

Function K_CheckProc_kconv_vol(cba) : CheckBoxControl
	STRUCT WMCheckboxAction &cba

	switch( cba.eventCode )
		case 2: // mouse up
			Variable checked = cba.checked
			if(checked==0)
				DoWindow/F K_ARPES_p
				SetVariable setvar_kconv_en,disable=0
				SetVariable setvar_kconv_en_s,disable=3
				SetVariable setvar_kconv_en_e,disable=3
			elseif(checked==1)
				DoWindow/F K_ARPES_p
				SetVariable setvar_kconv_en,disable=3
				SetVariable setvar_kconv_en_s,disable=0
				SetVariable setvar_kconv_en_e,disable=0
			endif
			
			break
		case -1: // control being killed
			break
	endswitch

	return 0
End

Function K_PopMenuProc_notation(pa) : PopupMenuControl
	STRUCT WMPopupAction &pa

	switch( pa.eventCode )
		case 2: // mouse up
			Variable popNum = pa.popNum
			String popStr = pa.popStr

			break
		case -1: // control being killed
			break
	endswitch

	return 0
End


Function K_SetVarProc_panelstep(sva) : SetVariableControl
	STRUCT WMSetVariableAction &sva

	switch( sva.eventCode )
		case 1: // mouse up
		case 2: // Enter key
		case 3: // Live update
			Variable dval = sva.dval
			String sval = sva.sval
			
			string nf=GetDataFolder(1)
			
			variable step
			
			if(dval>0)
				step=dval
			else
				step=1
				SetVariable setvar_step value= _NUM:1
			endif
			DoWindow/F K_ARPES_p

			SetVariable setvar_po,limits={-90,90,step}
			SetVariable setvar_ph,limits={-90,90,step}
			SetVariable setvar_az,limits={-inf,inf,step}
			SetVariable setvar_po_ofs,limits={-90,90,step}
			SetVariable setvar_ph_ofs,limits={-90,90,step}
			SetVariable setvar_az_ofs,limits={-inf,inf,step}
			
			SetVariable setvar_hn,limits={0,inf,step}
			SetVariable setvar_V0,limits={0,inf,step}
			SetVariable setvar_W,limits={0,inf,step}
			SetVariable setvar_EB,limits={0,inf,step}
			SetVariable setvar_th_s,limits={-90,90,step}
			SetVariable setvar_th_e,limits={-90,90,step}
			SetVariable setvar_bet_offs limits={-90,90,step}
			
			SetDataFolder $nf
			break
		case -1: // control being killed
			break
	endswitch

	return 0
End

Function K_TabProc_panelmode(tca) : TabControl
	STRUCT WMTabControlAction &tca

	switch( tca.eventCode )
		case 2: // mouse up
			Variable tab = tca.tab
			
			K_change_panelmode(tab)
			
			break
		case -1: // control being killed
			break
	endswitch

	return 0
End

Function K_change_panelmode(val)
	variable val // 0: sim, 1: kconv
	string nf=GetDataFolder(1)
	SetDataFolder root:K_ARPES:misc:automap
	variable/g v_amapval
	variable amap=v_amapval
	SetDataFolder root:K_ARPES:global
	variable/g v_kconv=val
	K_Set_Status(0,"")
	DoWindow/F K_ARPES_p
	TabControl tab_mode value=val
	if(val==0)
		
		variable/g v_kconv_vol=0
		string/g s_kconv_target=""
		K_show_targwave()
	
		Button button_go title="\\Zr200\\f01PLOT",disable=0
		GroupBox run_setting title="Plot setting"
		
		///
		SetVariable setvar_winname disable=0
		PopupMenu popup_winlist disable=0
		SetVariable setvar_show_min disable=0
		SetVariable setvar_show_max disable=0
		SetVariable setvar_show_wid disable=0
		SetVariable setvar_show_fix disable=0
		SetVariable setvar_kxmin disable=0
		SetVariable setvar_kxmax disable=0
		SetVariable setvar_kymin disable=0
		SetVariable setvar_kymax disable=0
		SetVariable setvar_kzmin disable=0
		SetVariable setvar_kzmax disable=0
		SetVariable setvar_kxwid disable=0
		SetVariable setvar_kywid disable=0
		SetVariable setvar_kzwid disable=0
		CheckBox check_fixkx disable=0
		CheckBox check_fixky disable=0
		CheckBox check_fixkz disable=0
		Button button_init_window disable=0
		CheckBox check_new disable=0
		CheckBox check_map disable=0
		SetVariable setvar_show_pr disable=0
		SetVariable setvar_EB disable=0
		Button button_xy disable=0
		Button button_yz disable=0
		Button button_xz disable=0
		
		///
		SetVariable setvar_show_ARPES_data disable=3
		Button button_select_wave disable=3
		CheckBox check_kconv_vol disable=3
		SetVariable setvar_kconv_en disable=3
		SetVariable setvar_kconv_kxn disable=3
		SetVariable setvar_kconv_kyn disable=3
		SetVariable setvar_kconv_ef disable=3
		SetVariable setvar_kconv_ek disable=3
		SetVariable setvar_kconv_eofs disable=3
		SetVariable setvar_kconv_ef_data disable=3
		PopupMenu popup_ydim disable=3
		SetVariable setvar_kconv_en_s,disable=3
		SetVariable setvar_kconv_en_e,disable=3
		SetVariable setvar_kconv_en_str,disable=3
		SetVariable setvar_outr disable=3
		
		///
		K_set_amapval(amap)
		
	elseif(val==1)
		Button button_go title="\\Zr200\\f01RUN"
		GroupBox run_setting title="k conv setting"
		
		///
		SetVariable setvar_winname disable=3
		PopupMenu popup_winlist disable=3
		SetVariable setvar_show_min disable=3
		SetVariable setvar_show_max disable=3
		SetVariable setvar_show_wid disable=3
		SetVariable setvar_show_fix disable=3
		SetVariable setvar_kxmin disable=3
		SetVariable setvar_kxmax disable=3
		SetVariable setvar_kymin disable=3
		SetVariable setvar_kymax disable=3
		SetVariable setvar_kzmin disable=3
		SetVariable setvar_kzmax disable=3
		SetVariable setvar_kxwid disable=3
		SetVariable setvar_kywid disable=3
		SetVariable setvar_kzwid disable=3
		CheckBox check_fixkx disable=3
		CheckBox check_fixky disable=3
		CheckBox check_fixkz disable=3
		Button button_init_window disable=3
		CheckBox check_new disable=3
		CheckBox check_map disable=3
		SetVariable setvar_show_pr disable=3
		SetVariable setvar_EB disable=2
		Button button_xy disable=3
		Button button_yz disable=3
		Button button_xz disable=3
		
		///
		SetVariable setvar_show_ARPES_data disable=0
		Button button_select_wave disable=0
		PopupMenu popup_ydim mode=1
		SetVariable setvar_outr disable=0
				
		///
		
	endif
	SetDataFolder $nf
End


Function K_ButtonProc_kconv(ba) : ButtonControl
	STRUCT WMButtonAction &ba

	switch( ba.eventCode )
		case 2: // mouse up
			// click code here
			string nf=GetDataFolder(1)

			SetDataFolder $nf
			break
		case -1: // control being killed
			break
	endswitch

	return 0
End


Function K_set_status(val,message)
	variable val
	string message
	
	DoWindow/F K_ARPES_p
	
	switch(val)
	case 0:
		SetVariable setvar_status value= _STR:"Idle",valueBackColor=(65535,54611,49151),help={message}
		break
	case 1:
		SetVariable setvar_status value= _STR:"Busy",valueBackColor=(65535,65534,49151),help={message}
		break
	case 2:
		SetVariable setvar_status value= _STR:"Warning",valueBackColor=(65535,54607,32768),help={message}
		break
	case 3:
		SetVariable setvar_status value= _STR:"Error",valueBackColor=(65535,49151,55704),help={message}
		break
	endswitch
End


Function K_SetVarProc_calc(sva) : SetVariableControl
	STRUCT WMSetVariableAction &sva

	switch( sva.eventCode )
		case 1: // mouse up
		case 2: // Enter key
		case 3: // Live update
			Variable dval = sva.dval
			String sval = sva.sval
			string  nf=GetDataFolder(1)
			SetDataFolder root:K_ARPES:misc:automap
			variable/g v_amapval
			variable amap=v_amapval
			nvar v_kconv=root:K_ARPES:global:v_kconv
			SetDataFolder root:K_ARPES:misc
			string/g s_winname
			variable/g v_app
			variable/g v_map
			variable/g v_smh
			
			K_set_emis_ang_window()
			
			if(!WinType(s_winname)==0 && v_kconv==0)
				if(amap==1)
					if(v_app==0)
						K_ARPES_calc_kline()
						K_append_curve()
					elseif(v_map==1)
						K_ARPES_calc_kline()
						K_append_curve()
					endif
				else
					if(v_app==0)
						K_auto_mapping(1)
					endif
				endif
				K_calcAxisRange_w()
			endif
			
			SetDataFolder $nf
			break
		case -1: // control being killed
			break
	endswitch

	return 0
End

Function K_PopMenuProc_amap_variable(pa) : PopupMenuControl
	STRUCT WMPopupAction &pa

	switch( pa.eventCode )
		case 2: // mouse up
			Variable popNum = pa.popNum
			String popStr = pa.popStr
			
			K_set_status(0,"")
			
			K_set_amapval(popNum)
			
			K_save_misc()
			break
		case -1: // control being killed
			break
	endswitch

	return 0
End


Function K_SliderProc_notation(sa) : SliderControl
	STRUCT WMSliderAction &sa

	switch( sa.eventCode )
		case -3: // Control received keyboard focus
		case -2: // Control lost keyboard focus
		case -1: // Control being killed
			break
		default:
			if( sa.eventCode & 1 ) // value set
				Variable curval = sa.curval
				string nf=GetDataFolder(1)
				SetDataFolder root:K_ARPES:misc:automap
				variable/g v_amapval
				variable amap=v_amapval
				SetDataFolder root:K_ARPES:global
				variable/g v_kconv
				variable kconv=v_kconv
				SetDataFolder root:K_ARPES:misc
				variable/g v_analtype=(curval-1)*-1
				string/g s_winname
				variable/g v_app
				variable/g v_map
				
				if(!WinType(s_winname)==0 && kconv==0)
					if(amap==1)
						if(v_app==0)
							K_ARPES_calc_kline()
							K_append_curve()
						elseif(v_map==1)
							K_ARPES_calc_kline()
							K_append_curve()
						endif
					else
						K_auto_mapping(1)
					endif
				endif
				
				SetDataFolder $nf
			endif
			break
	endswitch

	return 0
End



Function K_ButtonProc_make_window(ba) : ButtonControl
	STRUCT WMButtonAction &ba

	switch( ba.eventCode )
		case 2: // mouse up
			// click code here
			
			string nf=GetDataFolder(1)
			SetDataFolder root:K_ARPES:misc
			string/g s_winname
			string winn=s_winname
			variable/g v_map=0
			variable/g v_app=0
			variable/g v_curve_num=0
			variable/g v_curve_num_map=0
			SetDataFolder $nf
			if(WinType(winn)==1)
				KillWindow $winn
			endif
			K_make_window()
			SetDataFolder root:K_ARPES:Curves
			if(DataFolderExists(winn))
				KillDataFolder $winn
			endif
			SetDataFolder root:K_ARPES:global
			wave/T wl=$("win_list")
			variable i,di=nan
			for(i=0;i<DimSize(wl,0);i+=1)
				if(StringMatch(winn,wl[i]))
					di=i
				endif
				if(NumType(di)==0)
					DeletePoints di,1,wl
				endif
			endfor
			SetDataFolder $nf
			break
		case -1: // control being killed
			break
	endswitch

	return 0
End



Function K_SetVarProc_amap_parameter(sva) : SetVariableControl
	STRUCT WMSetVariableAction &sva
	
	switch( sva.eventCode )
		case 1: // mouse up
		case 2: // Enter key
		case 3: // Live update
			Variable dval = sva.dval
			String sval = sva.sval
			string nf=GetDataFolder(1)
			SetDataFolder root:K_ARPES:global
			variable/g v_kconv
			variable kconv=v_kconv
			SetDataFolder root:K_ARPES:misc:automap
			
			string/g s_amapval
			wave pw=$("amap_parameter")
			wave spw=$("amap_parameter_"+s_amapval)
			
			if(kconv==0)
				if(pw[1]==0)
					pw[1]=1
				endif
				
				pw[3]=round((pw[2]-pw[0])/pw[1]*10000)/10000+1
				
				if(round(pw[3])!=pw[3])
					pw[2]=pw[0]+round(pw[3])*pw[1]
					pw[3]=round(pw[3])
				endif
				
				spw=pw
				
			elseif(kconv==1)
				
				if(StringMatch(sva.ctrlName,"setvar_map_st"))
					pw[1]=(pw[2]-pw[0])/(pw[3]-1)
				elseif(StringMatch(sva.ctrlName,"setvar_map_step"))
					pw[2]=pw[0]+pw[1]*(pw[3]-1)
				elseif(StringMatch(sva.ctrlName,"setvar_map_en"))
					pw[1]=(pw[2]-pw[0])/(pw[3]-1)
				endif
				
			endif
			
			K_save_misc()
			
			SetDataFolder $nf
			break
		case -1: // control being killed
			break
	endswitch

	return 0
End

Function K_PopMenuProc_winlist(pa) : PopupMenuControl
	STRUCT WMPopupAction &pa

	switch( pa.eventCode )
		case 2: // mouse up
			Variable popNum = pa.popNum
			String popStr = pa.popStr
			string nf=GetDataFolder(1)
			
			K_restore_misc(popstr)
			K_refresh_panel()
			
			K_change_panelmode(0)
			K_refresh_panel()
			DoWindow/F $popstr
			
			SetDataFolder $nf
			break
		case -1: // control being killed
			break
	endswitch

	return 0
End


//Function K_q_run(ax1,ay1,az1,bx1,by1,bz1)
//	variable ax1,ay1,az1,bx1,by1,bz1
//	string nf=GetDataFolder(1)
//	SetDataFolder root:K_ARPES:misc:rot_matrix
//	K_q_make_q(ax1,ay1,az1,bx1,by1,bz1)
//	K_q_make_m()
//	SetDataFolder $nf
//End
//
//
//Function K_q_make_q(ax1,ay1,az1,bx1,by1,bz1)
//	variable ax1,ay1,az1,bx1,by1,bz1
//	
//	string vw1s="qv1_tempo"
//	string vw2s="qv2_tempo"
//	string vw3s="qv3_tempo"
//	string vw4s="qv4_tempo"
//	string vw5s="qv5_tempo"
//	string vw6s="qv"
//	make/O/N=3 $vw1s
//	make/O/N=3 $vw2s
//	wave vw1=$vw1s
//	wave vw2=$vw2s
//	make/O/N=4 $vw6s
//	wave vw6=$vw6s
//	
//	vw1[0]=ax1/sqrt(ax1^2+ay1^2+az1^2)
//	vw1[1]=ay1/sqrt(ax1^2+ay1^2+az1^2)
//	vw1[2]=az1/sqrt(ax1^2+ay1^2+az1^2)
//	
//	vw2[0]=bx1/sqrt(bx1^2+by1^2+bz1^2)
//	vw2[1]=by1/sqrt(bx1^2+by1^2+bz1^2)
//	vw2[2]=bz1/sqrt(bx1^2+by1^2+bz1^2)
//	
//	Cross/DEST=$vw3s vw1,vw2
//	wave vw3=$vw3s
//	variable vw3l=sqrt(vw3[0]^2+vw3[1]^2+vw3[2]^2)
//	vw3=vw3/vw3l
//	vw3l=-vw3l
//	
//	variable de=0.0002
//	variable dp=MatrixDot(vw1,vw2)
//	if(-de<vw3l || dp>1)
//		if(dp<de-1)
//			Make/O/N=3 $vw4s
//			wave vw4=$vw4s
//			vw4[0]=-vw1[1]
//			vw4[1]=vw1[2]
//			vw4[2]=vw1[0]
//			
//			Cross/DEST=$vw5s vw4,vw1
//			wave vw5=$vw5s
//			ax1=vw5[0]/sqrt(vw5[0]^2+vw5[1]^2+vw5[2]^2)
//			ay1=vw5[1]/sqrt(vw5[0]^2+vw5[1]^2+vw5[2]^2)
//			az1=vw5[2]/sqrt(vw5[0]^2+vw5[1]^2+vw5[2]^2)
//			
//			vw5[0]=ax1
//			vw5[1]=ay1
//			vw5[2]=az1
//			
//			vw6[0]=0
//			vw6[1]=vw5[0]
//			vw6[2]=vw5[1]
//			vw6[3]=vw5[2]
//			
//		else
//			vw6[0]=1
//			vw6[1]=0
//			vw6[2]=0
//			vw6[3]=0
//		endif
//	else
//		vw6[0]=sqrt((1+dp)/2)
//		vw6[1]=vw3[0]*sqrt((1-dp)/2)
//		vw6[2]=vw3[1]*sqrt((1-dp)/2)
//		vw6[3]=vw3[2]*sqrt((1-dp)/2)
//	endif
//	
//	KillWaves vw1,vw2,vw3,vw4,vw5
//End
//
//Function K_q_make_m()
//	string qws="qv"
//	wave qw=$qws
//	variable q0,q1,q2,q3
//	q0=qw[0]
//	q1=qw[1]
//	q2=qw[2]
//	q3=qw[3]
//	
//	qw[0]=q1
//	qw[1]=q2
//	qw[2]=q3
//	qw[3]=q0
//	
//	string nws="q_mw"
//	make/O/N=(3,3) $nws
//	wave nw=$nws
//	
//	nw[0][0]=qw[3]*qw[3]+qw[0]*qw[0]-qw[1]*qw[1]-qw[2]*qw[2]
//	nw[0][1]=2*qw[0]*qw[1]-2*qw[3]*qw[2]
//	nw[0][2]=2*qw[0]*qw[2]+2*qw[3]*qw[1]
////	nw[0][3]=0
//	
//	nw[1][0]=2*qw[0]*qw[1]+2*qw[3]*qw[2]
//	nw[1][1]=qw[3]*qw[3]-qw[0]*qw[0]+qw[1]*qw[1]-qw[2]*qw[2]
//	nw[1][2]=2*qw[1]*qw[2]-2*qw[3]*qw[0]
////	nw[1][3]=0
//
//	nw[2][0]=2*qw[0]*qw[2]-2*qw[3]*qw[1]
//	nw[2][1]=2*qw[1]*qw[2]+2*qw[3]*qw[0]
//	nw[2][2]=qw[3]*qw[3]-qw[0]*qw[0]-qw[1]*qw[1]+qw[2]*qw[2]
////	nw[2][3]=0
//
////	nw[3][0]=0
////	nw[3][1]=0
////	nw[3][2]=0
//	
//	nw=nw/(qw[3]*qw[3]+qw[0]*qw[0]+qw[1]*qw[1]+qw[2]*qw[2])
//	KillWaves qw
//End



///// Make Sample Wave
Function K_msw_Volumesample1()
	string sws="Volume_sample1"
	variable en,ax,ay
	variable en_s=-2
	variable en_n=101
	variable en_e=0.2
	
	variable ax_s=-15
	variable ax_n=201
	variable ax_e=15
	
	variable ay_s=-15
	variable ay_n=101
	variable ay_e=15
	
	make/O/N=(en_n,ax_n,ay_n) $sws
	wave sw=$sws
	sw=0
	
	variable en_st=(en_e-en_s)/(en_n-1)
	variable ax_st=(ax_e-ax_s)/(ax_n-1)
	variable ay_st=(ay_e-ay_s)/(ay_n-1)
	
	variable ew,ewid=2
	
	variable en_i,ax_i,ay_i
	for(ay_i=0;ay_i<ay_n;ay_i+=1)
		ay=ay_s+ay_i*ay_st
		for(ax_i=0;ax_i<ax_n;ax_i+=1)
			ax=ax_s+ax_i*ax_st
			en=(ax^2+ay^2)/100*1^2-2
			en_i=round((en-en_s)/en_st)
			
			for(ew=-ewid;ew<=ewid;ew+=1)
				if(en_i+ew>0 && en_i+ew<en_n)
					sw[en_i+ew][ax_i][ay_i]+=1
				endif
			endfor
		endfor
	endfor
	
	SetScale/P x,en_s,en_st,sw
	SetScale/P y,ax_s,ax_st,sw
	SetScale/P z,ay_s,ay_st,sw
End


Function K_ButtonProc_select_wave(ba) : ButtonControl
	STRUCT WMButtonAction &ba

	switch( ba.eventCode )
		case 2: // mouse up
			// click code here
			string nf=GetDataFolder(1)
			SetDataFolder root:K_ARPES:global
			if(Exists(GetBrowserSelection(0))==1)
				string/g s_kconv_target=GetBrowserSelection(0)
			endif
			K_show_targwave()
			SetDataFolder $nf
			break
		case -1: // control being killed
			break
	endswitch

	return 0
End

Function K_PopMenuProc_img(pa) : PopupMenuControl
	STRUCT WMPopupAction &pa

	switch( pa.eventCode )
		case 2: // mouse up
			Variable popNum = pa.popNum
			String popStr = pa.popStr
			K_kconv_2D_img(popNum)
			break
		case -1: // control being killed
			break
	endswitch

	return 0
End

Function K_kconv_2D_img(popNum)
	variable popNum
	string nf=GetDataFolder(1)
	nvar hn=root:K_ARPES:misc:v_hn
	nvar W=root:K_ARPES:misc:v_W
	SetDataFolder root:K_ARPES:global
	variable/g v_kconv_img=popNum
	string/g s_kconv_target
	string targ=s_kconv_target
	
	DoWindow/F K_ARPES_p
	if(popNum==1)
		SetVariable setvar_kconv_kxn disable=0,title="k points :"
		SetVariable setvar_kconv_kyn disable=3
		SetVariable setvar_kconv_ef_data disable=0
		SetVariable setvar_kconv_ek disable=3
		SetVariable setvar_kconv_eofs disable=0
		PopupMenu popup_ydim mode=1
		
		K_kconv_panel_targw_eofs(targ)
		
	elseif(popNum==2)
		SetVariable setvar_kconv_kxn disable=0,title="kx points :"
		SetVariable setvar_kconv_kyn disable=0
		SetVariable setvar_kconv_ef disable=0
		SetVariable setvar_kconv_ef_data disable=3
		SetVariable setvar_kconv_eofs disable=3
		SetVariable setvar_kconv_ek disable=0
		PopupMenu popup_ydim mode=2
		variable tex=K_ktw_search(targ)
		wave/T ktw=$("kconv_target_waves")
		if(tex<0)
			variable/g v_EF=hn-W
			variable/g v_ek=hn-W
		else
			variable/g v_EF=str2num(ktw[tex][1])
			variable/g v_ek=str2num(ktw[tex][3])
		endif
	endif
	K_write_ktw()
	SetDataFolder $nf
End

Function K_ButtonProc_go(ba) : ButtonControl
	STRUCT WMButtonAction &ba

	switch( ba.eventCode )
		case 2: // mouse up
			// click code here
			K_set_status(1,"")
			string  nf=GetDataFolder(1)
			SetDataFolder root:K_ARPES:misc:automap
			variable/g v_amapval
			variable amap=v_amapval
			SetDataFolder root:K_ARPES:global
			variable/g v_kconv
			variable kconv=v_kconv
			SetDataFolder root:K_ARPES:misc
			
			if(kconv==0)
				
				string/g s_winname
				string winn=s_winname
				
				if(!WinType(winn))
					SetDataFolder root:K_ARPES:misc
					variable/g v_map=0
					variable/g v_app=0
					
					if(WinType(winn)==1)
						KillWindow $winn
					endif
					K_make_window()
				endif
				
				if(amap==1)
					K_ARPES_calc_kline()
					K_append_curve()
					K_calcAxisRange_w()
				else
					K_auto_mapping(1)
				endif
				
			else
				
				SetDataFolder root:K_ARPES:global
				string/g s_kconv_target
				string targ=s_kconv_target
				
				if(WaveDims($targ)==2)
					variable/g v_kconv_img
					if(v_kconv_img==1)
						K_kconv_img_Ek()
					elseif(v_kconv_img==2)
						if(amap>1 && amap<6)
							K_kconv_img_CE()
						elseif(amap==6)
							print "hn map no kaiseki ha kansei site imasenn."
							print "hn map wo kaiseki sitai baai ha, kouno ni renraku site kudasai."
						else
							print "Please check mapping parameters."
						endif
					endif
				elseif(WaveDims($targ)==3)
					variable/g v_kconv_img=0
					if(amap>1 && amap<7)
						K_kconv_vol()
					elseif(amap==6)
						print "hn map no kaiseki ha kansei site imasenn."
						print "hn map wo kaiseki sitai baai ha, kouno ni renraku site kudasai."
					else
						print "Please check mapping parameters."
					endif
				else
					K_set_status(3,"Please check ARPES data.")
					print "Please check ARPES data."
				endif
			
			endif
			
			SetDataFolder $nf
			K_set_status(0,"")
			break
		case -1: // control being killed
			break
	endswitch

	return 0
End






Function K_SetVarProc_kconv_ef(sva) : SetVariableControl
	STRUCT WMSetVariableAction &sva

	switch( sva.eventCode )
		case 1: // mouse up
		case 2: // Enter key
		case 3: // Live update
			Variable dval = sva.dval
			String sval = sva.sval
			string nf=GetDataFolder(1)
			SetDataFolder root:K_ARPES:misc
			variable/g v_hn
			variable/g v_W=v_hn-dval
			SetDataFolder root:K_ARPES:global
			variable/g v_eofs
			variable/g v_EF_data=dval-v_eofs
			SetDataFolder $nf
			K_write_ktw()
			break
		case -1: // control being killed
			break
	endswitch

	return 0
End


Function K_pm_init()
	NewDataFolder/O/S root:K_ARPES:pm_panel
	make/O/N=(3,3,3) $("panel_value")
	variable/g v_p=0
	variable/g v_t=0
	variable/g v_a=0
	variable/g v_pn=0
	variable/g v_tn=0
	variable/g v_an=0
	variable/g v_po=0
	variable/g v_to=0
	variable/g v_ao=0
End

Window K_pm() : Panel
	K_pm_init()
	PauseUpdate; Silent 1		// building window...
	NewPanel /K=1 /W=(484,145,1282,360) as "Parameter set"
	ModifyPanel fixedSize=1
	SetDrawLayer UserBack
	DrawText 21,37,"Polar"
	DrawText 105,37,"Tilt"
	DrawText 179,37,"Azimuth"
	DrawText 278,37,"Polar"
	DrawText 362,37,"Tilt"
	DrawText 436,37,"Azimuth"
	DrawText 546,37,"Polar"
	DrawText 630,37,"Tilt"
	DrawText 704,37,"Azimuth"
	SetVariable setvar_p,pos={31.00,39.00},size={50.00,18.00},proc=K_pm_SetVarProc_update,title=" "
	SetVariable setvar_p,limits={-inf,inf,0},value= root:K_ARPES:pm_panel:v_p
	SetVariable setvar_pmax,pos={49.00,60.00},size={50.00,18.00},proc=K_pm_SetVarProc_maxmin,title=" "
	SetVariable setvar_pmax,limits={-inf,inf,0},value= root:K_ARPES:pm_panel:panel_value[0]
	SetVariable setvar_pmin,pos={50.00,182.00},size={50.00,18.00},proc=K_pm_SetVarProc_maxmin,title=" "
	SetVariable setvar_pmin,limits={-inf,inf,0},value= root:K_ARPES:pm_panel:panel_value[1]
	Slider slider_p,pos={20.00,64.00},size={22.00,136.00},proc=K_pm_SliderProc_update
	Slider slider_p,limits={0,1,0},variable= root:K_ARPES:pm_panel:v_p,ticks= 0
	SetVariable setvar_t,pos={110.00,39.00},size={50.00,18.00},proc=K_pm_SetVarProc_update,title=" "
	SetVariable setvar_t,limits={-inf,inf,0},value= root:K_ARPES:pm_panel:v_t
	SetVariable setvar_tmax,pos={128.00,60.00},size={50.00,18.00},proc=K_pm_SetVarProc_maxmin,title=" "
	SetVariable setvar_tmax,limits={-inf,inf,0},value= root:K_ARPES:pm_panel:panel_value[0][0][1]
	SetVariable setvar_tmin,pos={129.00,182.00},size={50.00,18.00},proc=K_pm_SetVarProc_maxmin,title=" "
	SetVariable setvar_tmin,limits={-inf,inf,0},value= root:K_ARPES:pm_panel:panel_value[1][0][1]
	Slider slider_t,pos={99.00,64.00},size={22.00,136.00},proc=K_pm_SliderProc_update
	Slider slider_t,limits={0,1,0},variable= root:K_ARPES:pm_panel:v_t,ticks= 0
	SetVariable setvar_a,pos={190.00,39.00},size={50.00,18.00},proc=K_pm_SetVarProc_update,title=" "
	SetVariable setvar_a,limits={-inf,inf,0},value= root:K_ARPES:pm_panel:v_a
	SetVariable setvar_pmax2,pos={208.00,60.00},size={50.00,18.00},proc=K_pm_SetVarProc_maxmin,title=" "
	SetVariable setvar_pmax2,limits={-inf,inf,0},value= root:K_ARPES:pm_panel:panel_value[0][0][2]
	SetVariable setvar_amin,pos={209.00,182.00},size={50.00,18.00},proc=K_pm_SetVarProc_maxmin,title=" "
	SetVariable setvar_amin,limits={-inf,inf,0},value= root:K_ARPES:pm_panel:panel_value[1][0][2]
	Slider slider_a,pos={179.00,64.00},size={22.00,136.00},proc=K_pm_SliderProc_update
	Slider slider_a,limits={0,1,0},variable= root:K_ARPES:pm_panel:v_a,ticks= 0
	SetVariable setvar_pn,pos={287.00,39.00},size={50.00,18.00},proc=K_pm_SetVarProc_update,title=" "
	SetVariable setvar_pn,limits={-inf,inf,0},value= root:K_ARPES:pm_panel:v_pn
	SetVariable setvar_pnmax,pos={305.00,60.00},size={50.00,18.00},proc=K_pm_SetVarProc_maxmin,title=" "
	SetVariable setvar_pnmax,limits={-inf,inf,0},value= root:K_ARPES:pm_panel:panel_value[0][1]
	SetVariable setvar_pnmin,pos={306.00,182.00},size={50.00,18.00},proc=K_pm_SetVarProc_maxmin,title=" "
	SetVariable setvar_pnmin,limits={-inf,inf,0},value= root:K_ARPES:pm_panel:panel_value[1][1]
	Slider slider_pn,pos={276.00,64.00},size={22.00,136.00},proc=K_pm_SliderProc_update
	Slider slider_pn,limits={0,1,0},variable= root:K_ARPES:pm_panel:v_pn,ticks= 0
	SetVariable setvar_tn,pos={366.00,39.00},size={50.00,18.00},proc=K_pm_SetVarProc_update,title=" "
	SetVariable setvar_tn,limits={-inf,inf,0},value= root:K_ARPES:pm_panel:v_tn
	SetVariable setvar_tnmax,pos={384.00,60.00},size={50.00,18.00},proc=K_pm_SetVarProc_maxmin,title=" "
	SetVariable setvar_tnmax,limits={-inf,inf,0},value= root:K_ARPES:pm_panel:panel_value[0][1][1]
	SetVariable setvar_tnmin,pos={385.00,182.00},size={50.00,18.00},proc=K_pm_SetVarProc_maxmin,title=" "
	SetVariable setvar_tnmin,limits={-inf,inf,0},value= root:K_ARPES:pm_panel:panel_value[1][1][1]
	Slider slider_tn,pos={355.00,64.00},size={22.00,136.00},proc=K_pm_SliderProc_update
	Slider slider_tn,limits={0,2,0},variable= root:K_ARPES:pm_panel:v_tn,ticks= 0
	SetVariable setvar_an,pos={446.00,39.00},size={50.00,18.00},proc=K_pm_SetVarProc_update,title=" "
	SetVariable setvar_an,limits={-inf,inf,0},value= root:K_ARPES:pm_panel:v_an
	SetVariable setvar_anmax,pos={464.00,60.00},size={50.00,18.00},proc=K_pm_SetVarProc_maxmin,title=" "
	SetVariable setvar_anmax,limits={-inf,inf,0},value= root:K_ARPES:pm_panel:panel_value[0][1][2]
	SetVariable setvar_anmin,pos={465.00,182.00},size={50.00,18.00},proc=K_pm_SetVarProc_maxmin,title=" "
	SetVariable setvar_anmin,limits={-inf,inf,0},value= root:K_ARPES:pm_panel:panel_value[1][1][2]
	Slider slider_an,pos={435.00,64.00},size={22.00,136.00},proc=K_pm_SliderProc_update
	Slider slider_an,limits={0,1,0},variable= root:K_ARPES:pm_panel:v_an,ticks= 0
	SetVariable setvar_po,pos={554.00,39.00},size={50.00,18.00},proc=K_pm_SetVarProc_update,title=" "
	SetVariable setvar_po,limits={-inf,inf,0},value= root:K_ARPES:pm_panel:v_po
	SetVariable setvar_pomax,pos={572.00,60.00},size={50.00,18.00},proc=K_pm_SetVarProc_maxmin,title=" "
	SetVariable setvar_pomax,limits={-inf,inf,0},value= root:K_ARPES:pm_panel:panel_value[0][2]
	SetVariable setvar_pomin,pos={573.00,182.00},size={50.00,18.00},proc=K_pm_SetVarProc_maxmin,title=" "
	SetVariable setvar_pomin,limits={-inf,inf,0},value= root:K_ARPES:pm_panel:panel_value[1][2]
	Slider slider_po,pos={543.00,64.00},size={22.00,136.00},proc=K_pm_SliderProc_update
	Slider slider_po,limits={0,1,0},variable= root:K_ARPES:pm_panel:v_po,ticks= 0
	SetVariable setvar_to,pos={633.00,39.00},size={50.00,18.00},proc=K_pm_SetVarProc_update,title=" "
	SetVariable setvar_to,limits={-inf,inf,0},value= root:K_ARPES:pm_panel:v_to
	SetVariable setvar_tomax,pos={651.00,60.00},size={50.00,18.00},proc=K_pm_SetVarProc_maxmin,title=" "
	SetVariable setvar_tomax,limits={-inf,inf,0},value= root:K_ARPES:pm_panel:panel_value[0][2][1]
	SetVariable setvar_tomin,pos={652.00,182.00},size={50.00,18.00},proc=K_pm_SetVarProc_maxmin,title=" "
	SetVariable setvar_tomin,limits={-inf,inf,0},value= root:K_ARPES:pm_panel:panel_value[1][2][1]
	Slider slider_to,pos={622.00,64.00},size={22.00,136.00},proc=K_pm_SliderProc_update
	Slider slider_to,limits={0,1,0},variable= root:K_ARPES:pm_panel:v_to,ticks= 0
	SetVariable setvar_ao,pos={713.00,39.00},size={50.00,18.00},proc=K_pm_SetVarProc_update,title=" "
	SetVariable setvar_ao,limits={-inf,inf,0},value= root:K_ARPES:pm_panel:v_ao
	SetVariable setvar_aomax,pos={731.00,60.00},size={50.00,18.00},proc=K_pm_SetVarProc_maxmin,title=" "
	SetVariable setvar_aomax,limits={-inf,inf,0},value= root:K_ARPES:pm_panel:panel_value[0][2][2]
	SetVariable setvar_aomin,pos={732.00,182.00},size={50.00,18.00},proc=K_pm_SetVarProc_maxmin,title=" "
	SetVariable setvar_aomin,limits={-inf,inf,0},value= root:K_ARPES:pm_panel:panel_value[1][2][2]
	Slider slider_ao,pos={702.00,64.00},size={22.00,136.00},proc=K_pm_SliderProc_update
	Slider slider_ao,limits={0,1,0},variable= root:K_ARPES:pm_panel:v_ao,ticks= 0
	GroupBox group_v,pos={8.00,3.00},size={256.00,205.00},title="Value"
	GroupBox group_n,pos={269.00,3.00},size={256.00,205.00},title="Normal"
	GroupBox group_o,pos={535.00,3.00},size={256.00,205.00},title="Offset"
	K_pm_maxmin()
EndMacro

Function K_pm_maxmin()
	string nf=GetDataFolder(1)
	wave pw=$("root:K_ARPES:pm_panel:panel_value")
	
	
	Slider slider_p,limits={pw[1][0][0],pw[0][0][0],pw[2][0][0]},win=K_pm
	Slider slider_t,limits={pw[1][0][1],pw[0][0][1],pw[2][0][1]},win=K_pm
	Slider slider_a,limits={pw[1][0][2],pw[0][0][2],pw[2][0][2]},win=K_pm
	
	Slider slider_pn,limits={pw[1][1][0],pw[0][1][0],pw[2][1][0]},win=K_pm
	Slider slider_tn,limits={pw[1][1][1],pw[0][1][1],pw[2][1][1]},win=K_pm
	Slider slider_an,limits={pw[1][1][2],pw[0][1][2],pw[2][1][2]},win=K_pm
	
	Slider slider_po,limits={pw[1][2][0],pw[0][2][0],pw[2][2][0]},win=K_pm
	Slider slider_to,limits={pw[1][2][1],pw[0][2][1],pw[2][2][1]},win=K_pm
	Slider slider_ao,limits={pw[1][2][2],pw[0][2][2],pw[2][2][2]},win=K_pm
	
	SetDataFolder $nf
End

Function K_pm_SetVarProc_maxmin(sva) : SetVariableControl
	STRUCT WMSetVariableAction &sva

	switch( sva.eventCode )
		case 1: // mouse up
		case 2: // Enter key
		case 3: // Live update
			Variable dval = sva.dval
			String sval = sva.sval
			K_pm_maxmin()
			break
		case -1: // control being killed
			break
	endswitch

	return 0
End

Function K_pm_SetVarProc_update(sva) : SetVariableControl
	STRUCT WMSetVariableAction &sva

	switch( sva.eventCode )
		case 1: // mouse up
		case 2: // Enter key
		case 3: // Live update
			Variable dval = sva.dval
			String sval = sva.sval
				K_pm_function()
			break
		case -1: // control being killed
			break
	endswitch

	return 0
End

Function K_pm_SliderProc_update(sa) : SliderControl
	STRUCT WMSliderAction &sa

	switch( sa.eventCode )
		case -3: // Control received keyboard focus
		case -2: // Control lost keyboard focus
		case -1: // Control being killed
			break
		default:
			if( sa.eventCode & 1 ) // value set
				Variable curval = sa.curval
				K_pm_function()
			endif
			break
	endswitch

	return 0
End

Function K_pm_function()
	string nf=GetDataFolder(1)
	SetDataFolder root:K_ARPES:pm_panel
	nvar p=v_p
	nvar t=v_t
	nvar pn=v_pn
	nvar tn=v_tn
	nvar an=v_an
	nvar po=v_po
	nvar to=v_to
	nvar ao=v_ao
	SetDataFolder $nf
	
	K_a2k_example(p,t,pn,tn,an,po,to,ao)
End

Function K_ButtonProc_cross(ba) : ButtonControl
	STRUCT WMButtonAction &ba
	
	switch( ba.eventCode )
		case 2: // mouse up
			// click code here
			string nf=GetDataFolder(1)
			SetDataFolder root:K_ARPES:misc
			variable/g v_cross_x
			variable/g v_cross_y
			variable/g v_cross_z
			
			variable vx=v_cross_x
			variable vy=v_cross_y
			variable vz=v_cross_z
			if(StringMatch(ba.ctrlName,"button_xy"))
				if(vx==0)
					variable/g v_cross_x=1
				elseif(vx==1)
					if(vz==1)
						variable/g v_cross_x=0
						variable/g v_cross_y=0
					elseif(vz==0)
						variable/g v_cross_x=0
						variable/g v_cross_z=1
					endif
				endif
			endif
			if(StringMatch(ba.ctrlName,"button_yz"))
				if(vy==0)
					variable/g v_cross_y=1
					variable/g v_cross_x=1
					variable/g v_cross_z=1
				elseif(vy==1)
					variable/g v_cross_y=0
				endif
			endif
			if(StringMatch(ba.ctrlName,"button_xz"))
				if(vz==0)
					variable/g v_cross_z=1
				elseif(vz==1)
					if(vx==1)
						variable/g v_cross_z=0
						variable/g v_cross_y=0
					elseif(vx==0)
						variable/g v_cross_z=0
						variable/g v_cross_x=1
					endif
				endif
			endif
			K_select_crosssection()
			K_SetAxis_w()
			K_panel_update_cross()
			SetDataFolder $nf
			break
		case -1: // control being killed
			break
	endswitch

	return 0
End

Function K_panel_update_cross()
	string nf=GetDataFolder(1)
	SetDataFolder root:K_ARPES:misc
	nvar vx=v_cross_x
	nvar vy=v_cross_y
	nvar vz=v_cross_z
	
	if(vx==0)
		Button button_xy fColor=(0,0,0),win=K_ARPES_p
	elseif(vx==1)
		Button button_xy fColor=(0,65535,65535),win=K_ARPES_p
	endif
	
	if(vy==0)
		Button button_yz fColor=(0,0,0),win=K_ARPES_p
	elseif(vy==1)
		Button button_yz fColor=(0,65535,65535),win=K_ARPES_p
	endif
	
	if(vz==0)
		Button button_xz fColor=(0,0,0),win=K_ARPES_p
	elseif(vz==1)
		Button button_xz fColor=(0,65535,65535),win=K_ARPES_p
	endif
	
	SetDataFolder $nf
End

Function K_ButtonProc_m_zero(ba) : ButtonControl
	STRUCT WMButtonAction &ba

	switch( ba.eventCode )
		case 2: // mouse up
			// click code here
			K_mani_all_zero()
			break
		case -1: // control being killed
			break
	endswitch

	return 0
End

Function K_write_ktw()
	string nf=GetDataFolder(1)
	
	SetDataFolder root:K_ARPES:misc
	nvar hn=v_hn
	nvar W=v_W
	SetDataFolder root:K_ARPES:global
	wave/T ktw=$("kconv_target_waves")
	svar targ=s_kconv_target
	nvar ef=v_EF
	nvar ef_data=v_EF_data
	nvar eofs=v_eofs
	nvar e_e=v_e_e
	nvar e_st=v_e_st
	nvar e_s=v_e_s
	nvar kx_n=v_kx_n
	nvar ky_n=v_ky_n
	nvar outr=v_outr
	nvar ek=v_ek
	nvar kconv_en=v_kconv_en
	nvar img=v_kconv_img
	
	variable tex=K_ktw_search(targ)
	
	if(tex<0)
		tex=DimSize(ktw,0)
		InsertPoints tex,1,ktw
		ktw[tex][0]=targ
	endif
	
	if(NumType(ek)==2)
		ek=hn-W
	endif
	
	ktw[tex][1]=num2str(ef)
	ktw[tex][2]=num2str(ef_data)
	ktw[tex][3]=num2str(ek)
	ktw[tex][4]=num2str(eofs)
	ktw[tex][5]=num2str(e_s)
	ktw[tex][6]=num2str(e_st)
	ktw[tex][7]=num2str(kconv_en)
	ktw[tex][8]=num2str(kx_n)
	ktw[tex][9]=num2str(ky_n)
	ktw[tex][10]=num2str(img)
  
	SetDataFolder $nf
End

Function K_SetVarProc_kconv_eofs(sva) : SetVariableControl
	STRUCT WMSetVariableAction &sva

	switch( sva.eventCode )
		case 1: // mouse up
		case 2: // Enter key
		case 3: // Live update
			Variable dval = sva.dval
			String sval = sva.sval
			string nf=GetDataFolder(1)
			SetDataFolder root:K_ARPES:global
			variable/g v_EF
			variable/g v_EF_data=v_EF-dval
			SetDataFolder $nf
			K_write_ktw()
			break
		case -1: // control being killed
			break
	endswitch

	return 0
End

Function K_SetVarProc_kconv_ef_data(sva) : SetVariableControl
	STRUCT WMSetVariableAction &sva

	switch( sva.eventCode )
		case 1: // mouse up
		case 2: // Enter key
		case 3: // Live update
			Variable dval = sva.dval
			String sval = sva.sval
			string nf=GetDataFolder(1)
			SetDataFolder root:K_ARPES:global
			variable/g v_EF
			variable/g v_eofs=v_EF-dval
			variable/g v_kconv_en=dval
			SetDataFolder $nf
			K_write_ktw()
			break
		case -1: // control being killed
			break
	endswitch

	return 0
End

Function K_CheckProc_pkf(cba) : CheckBoxControl
	STRUCT WMCheckboxAction &cba

	switch( cba.eventCode )
		case 2: // mouse up
			Variable checked = cba.checked
			string  nf=GetDataFolder(1)
			SetDataFolder root:K_ARPES:misc:automap
			variable/g v_amapval
			variable amap=v_amapval
			nvar v_kconv=root:K_ARPES:global:v_kconv
			SetDataFolder root:K_ARPES:misc
			string/g s_winname
			variable/g v_app
			variable/g v_map
			variable/g v_smh
			
			K_set_emis_ang_window()
			
			if(!WinType(s_winname)==0 && v_kconv==0)
				if(amap==1)
					if(v_app==0)
						K_ARPES_calc_kline()
						K_append_curve()
					elseif(v_map==1)
						K_ARPES_calc_kline()
						K_append_curve()
					endif
				else
					if(v_app==0)
						K_auto_mapping(1)
					endif
				endif
				K_calcAxisRange_w()
			endif
			
			SetDataFolder $nf

			break
		case -1: // control being killed
			break
	endswitch

	return 0
End

Function K_SetVarProc_kconv_en_r(sva) : SetVariableControl
	STRUCT WMSetVariableAction &sva

	switch( sva.eventCode )
		case 1: // mouse up
		case 2: // Enter key
		case 3: // Live update
			Variable dval = sva.dval
			String sval = sva.sval
			K_kconv_check_enrs_panel()
			K_write_ktw()
			break
		case -1: // control being killed
			break
	endswitch

	return 0
End

Function K_SetVarProc_kconv_write_ktw(sva) : SetVariableControl
	STRUCT WMSetVariableAction &sva

	switch( sva.eventCode )
		case 1: // mouse up
		case 2: // Enter key
		case 3: // Live update
			Variable dval = sva.dval
			String sval = sva.sval
			K_write_ktw()
			break
		case -1: // control being killed
			break
	endswitch

	return 0
End


Static Function K_ktw_search(targ)
	string targ
	string nf=GetDataFolder(1)
	SetDataFolder root:K_ARPES:global
	variable i,tex=-1
	if(WaveExists($("kconv_target_waves")))
		wave/T ktw=$("kconv_target_waves")
		for(i=0;i<DimSize(ktw,0);i+=1)
			if(StringMatch(ktw[i][0], targ))
				tex=i
			endif
		endfor
		if(tex<0)
			tex=-1
		endif
	else
		tex=-1
	endif
	SetDataFolder $nf
	return tex
End

Function K_ev2ai(ev)
	variable ev
	variable hc=1.23984193e4 // [ev*AA]
	variable ai=ev/hc
	return ai
End
