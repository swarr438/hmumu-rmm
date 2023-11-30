/***************************************
This is a ROOT macro used to generate RMMs
Run it with: root -q filename.C
****************************************/

#include<vector>
#define MMU 105.7 

// Reproduction of S.Chekanov\'s codes using ROOT
// return angle product
double getAngle(const TLorentzVector p1, const TLorentzVector p2){
	double HL;
	if (p1.Et()>0 && p2.Et()>0) {
		double y1=p1.Rapidity();
		double y2=p2.Rapidity();
		HL=0.1*( TMath::CosH( 0.5*(y2-y1) )-1);
	} else HL=0;
	return HL;
}

// get masses
double getMass(const TLorentzVector p1, const TLorentzVector p2){
	TLorentzVector pp=p1+p2;
	double xmass=pp.M();
	if (xmass<0) xmass=0;
	//if (xmass>1) cout << "xmass=" << xmass << " CMS=" << CMS << endl;
	return xmass;
}

double getHL(const TLorentzVector p1){
	double HL;
	if (p1.Et()>0) {
		double y=p1.Rapidity();
		HL=0.1*(TMath::CosH(y)-1);
	} else HL=0;
	return HL;
}

// return transverse mass using experimental method 
double getMT(const TLorentzVector met, const TLorentzVector jet) {
	double ss= (jet.Et()+met.Et()) * (jet.Et() +met.Et()) -
                ( (jet.Px()+met.Px())*(jet.Px() +met.Px()) ) - ( (jet.Py()+met.Py())*(jet.Py() +met.Py()) );
     	double Mt_exact=0;
     	if (ss>0) Mt_exact= TMath::Sqrt(ss);
	return Mt_exact;
}

// map2rmm
double** map2rmm(const double CMS, const int maxN, const int maxTypes,
		const TLorentzVector missing,
		const vector<TLorentzVector> jets,
		const vector<TLorentzVector> bjets,
		const vector<TLorentzVector> muons) {

	const int mSize = maxN*maxTypes+1;
	double** outMatrix = 0;
	outMatrix = new double*[mSize];
	for (int h=0; h<mSize; h++) outMatrix[h] = new double[mSize];

	unsigned int INCR = 0;
	for(int i=0; i<mSize; i++)
		for(int j=0; j<mSize; j++)  outMatrix[i][j] = 0;
	outMatrix[0][0] = missing.Et()/CMS;
	
	// jets	
	for (unsigned int k1=0; k1<maxN; k1++) {
		TLorentzVector jet1 = jets.at(k1);
		if (missing.Et()>0) outMatrix[0][k1+INCR*maxN+1] = getMT(missing, jet1)/CMS;
		outMatrix[k1+INCR*maxN+1][0] = getHL(jet1);
				
		// diagonal
		if (k1 == 0) {
			outMatrix[k1+INCR*maxN+1][k1+INCR*maxN+1] = jet1.Et()/CMS;
		} else if (k1>0)  {
			TLorentzVector jet3 = jets.at(k1-1);
			double imbalance = 0;
			if (jet3.Et()!=0) imbalance = (jet3.Et() - jet1.Et())/(jet3.Et() + jet1.Et());
			outMatrix[k1+INCR*maxN+1][k1+INCR*maxN+1] = imbalance;
		}
		
		// nondiagonal
		for (unsigned int k2 = 0; k2<maxN; k2++) {
			TLorentzVector jet2 = jets.at(k2);
			if (k1<k2)    outMatrix[k1+INCR*maxN+1][k2+INCR*maxN+1] = getMass(jet1,jet2)/CMS;
			if (k1>k2)    outMatrix[k1+INCR*maxN+1][k2+INCR*maxN+1] = getAngle(jet1,jet2);
		}
	}
	
	// bjets
	INCR = 1;	
	for (unsigned int k1=0; k1<maxN; k1++) {
		TLorentzVector bet1 = bjets.at(k1);
		if (missing.Et()>0) outMatrix[0][k1+INCR*maxN+1] = getMT(missing, bet1)/CMS;
		outMatrix[k1+INCR*maxN+1][0] = getHL(bet1);
		//cout << bet1.Et()<<" "<<bet1.Py()<<" "<<getHL(bet1)<<endl;
				
		// diagonal
		if (k1 == 0) {
			outMatrix[k1+INCR*maxN+1][k1+INCR*maxN+1] = bet1.Et()/CMS;
		} else if (k1>0)  {
			TLorentzVector bet3 = bjets.at(k1-1);
			double imbalance = 0;
			if (bet3.Et()!=0) imbalance = (bet3.Et() - bet1.Et())/(bet3.Et() + bet1.Et());
			outMatrix[k1+INCR*maxN+1][k1+INCR*maxN+1] = imbalance;
		}
		
		// nondiagonal
		for (unsigned int k2 = 0; k2<maxN; k2++) {
			TLorentzVector bet2 = bjets.at(k2);
			if (k1<k2)    outMatrix[k1+INCR*maxN+1][k2+INCR*maxN+1] = getMass(bet1,bet2)/CMS;
			if (k1>k2)    outMatrix[k1+INCR*maxN+1][k2+INCR*maxN+1] = getAngle(bet1,bet2);
		}
	}

	// muons	
	INCR = 2;
	for (unsigned int k1=0; k1<maxN; k1++) {
		TLorentzVector mu1 = muons.at(k1);
		if (missing.Et()>0) outMatrix[0][k1+INCR*maxN+1] = getMT(missing, mu1)/CMS;
		outMatrix[k1+INCR*maxN+1][0] = getHL(mu1);
				
		// diagonal
		if (k1 == 0) {
			outMatrix[k1+INCR*maxN+1][k1+INCR*maxN+1] = mu1.Et()/CMS;
		} else if (k1>0)  {
			TLorentzVector mu3 = muons.at(k1-1);
			double imbalance = (mu3.Et() - mu1.Et())/(mu3.Et() + mu1.Et());
			outMatrix[k1+INCR*maxN+1][k1+INCR*maxN+1] = imbalance;
		}
		
		// nondiagonal
		for (unsigned int k2 = 0; k2<maxN; k2++) {
			TLorentzVector mu2 = muons.at(k2);
			if (k1<k2)    outMatrix[k1+INCR*maxN+1][k2+INCR*maxN+1] = getMass(mu1,mu2)/CMS;
			if (k1>k2)    outMatrix[k1+INCR*maxN+1][k2+INCR*maxN+1] = getAngle(mu1,mu2);
		}
	}
	
	// cross elements
	unsigned int INCR_shift = 0; // for jets=1
	for (unsigned int k1=0;k1<maxN;k1++) {
		TLorentzVector p1 = jets.at(k1);
		
		INCR = 1;
		for (unsigned int k2=0; k2<maxN; k2++) {
			TLorentzVector p2 = bjets.at(k2);
			outMatrix[k1+INCR_shift*maxN+1][k2+INCR*maxN+1] = getMass(p1,p2)/CMS;
			outMatrix[k2+INCR*maxN+1][k1+INCR_shift*maxN+1] = getAngle(p1,p2);
		}
	
		INCR = 2;
		for (unsigned int k2=0; k2<maxN; k2++) {
			TLorentzVector p2 = muons.at(k2);
			outMatrix[k1+INCR_shift*maxN+1][k2+INCR*maxN+1] = getMass(p1,p2)/CMS;
			outMatrix[k2+INCR*maxN+1][k1+INCR_shift*maxN+1] = getAngle(p1,p2);
		}
	}

	INCR_shift = 1; // for bjets=2
	for (unsigned int k1=0;k1<maxN;k1++) {
		TLorentzVector p1 = bjets.at(k1);
	
		INCR = 2;
		for (unsigned int k2=0; k2<maxN; k2++) {
			TLorentzVector p2 = muons.at(k2);
			outMatrix[k1+INCR_shift*maxN+1][k2+INCR*maxN+1] = getMass(p1,p2)/CMS;
			outMatrix[k2+INCR*maxN+1][k1+INCR_shift*maxN+1] = getAngle(p1,p2);
		}
	}

	return outMatrix;
}

//sort by Et
bool compare(const TLorentzVector & l, const TLorentzVector & r)
{
    return l.Et() > r.Et();
}


//read event indices from skimmed data
int readindex(vector<long>& eventnumber, const string x) {
	ifstream file(x.c_str());
	if (file.is_open()) {
	    	string line;
		int i=0;
		while (getline(file, line)) {
			if (i==0) i++;
			else {
				eventnumber.push_back(stol(line.c_str()));
				i++;
			}
		}
	}

	file.close();
	
	return 0;
}

// read mc id
vector<string> readdsid(const string channel) {
	string x = channel+"/index/dsid.txt"; 
	ifstream file(x.c_str());
	vector<string> dsid;
	if (file.is_open()) {
	    	string line;
		while (getline(file, line)) {
			dsid.push_back(line);
		}
	}

	file.close();
	
	return dsid;
}

// main
void rmmf() {
	string channel; // "ttbar";
	cout << "Channel Name.:\t";
	cin >> channel;
	vector<string> snum = {"a","d","e"};
	vector<string> dsid = readdsid(channel);
	// auto tree = new TChain("DiMuonNtuple");
	
	ofstream oset(channel+".txt");
	int nsig=0;
	int filenum = dsid.size();
	int setnum = 100000;
	int loopnum = setnum / filenum + 1;
	float wt = 0;	

	for (int j=0;j<filenum;j++) {
		for (int i=0;i<snum.size();i++) {				

			string datapath = channel+"/mc16"+snum[i]+"."+dsid[j];
			string oripath = "/eos/atlas/atlascerngroupdisk/phys-higgs/HSG2/Hmumu/common_ntuples/v23/mc16"+snum[i]+"."+dsid[j]+".root";
			cout << "----------------------------------" << endl;
			cout << "Processing " << oripath << endl;
			
			auto file = TFile::Open(oripath.c_str());
			if (!file || file->IsZombie()) {
				std::cerr << "Error opening file" << endl;
				exit(-1);
			}

			auto tree = file->Get<TTree>("DiMuonNtuple");

			// read event number from index
			size_t pos = 0;
			string delimiter = "/";
			string token;
			while ((pos = datapath.find(delimiter)) != std::string::npos) {
				token = datapath.substr(0, pos);
				datapath.erase(0, pos + delimiter.length());
			}
			string name = token+"/index/"+datapath+".root.index";
			vector<long> eventnumber {};
			readindex(eventnumber,name);
	
			// randomly rearrange elements
			// unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
			// std::sort(eventnumber.begin(),eventnumber.end());
			// std::shuffle(eventnumber.begin(),eventnumber.end(),std::default_random_engine(seed));

			int inen = eventnumber.size();
			// cout << "----------------------------------" << endl;
					
			// Disable everything but the branches we need
			tree->SetBranchStatus("*", false);
			tree->SetBranchStatus("Jets_E_Lead",true);
			tree->SetBranchStatus("Jets_PT_Lead",true);
			tree->SetBranchStatus("Jets_Eta_Lead",true);
			tree->SetBranchStatus("Jets_Phi_Lead",true);
			tree->SetBranchStatus("Jets_E_Sub",true);
			tree->SetBranchStatus("Jets_PT_Sub",true);
			tree->SetBranchStatus("Jets_Eta_Sub",true);
			tree->SetBranchStatus("Jets_Phi_Sub",true);
			tree->SetBranchStatus("Jets_PartonID_Lead",true);
			tree->SetBranchStatus("Jets_PartonID_Sub",true);
			tree->SetBranchStatus("Muons_Minv_MuMu",true);
			tree->SetBranchStatus("Muons_PT_Lead",true);
			tree->SetBranchStatus("Muons_Eta_Lead",true);
			tree->SetBranchStatus("Muons_Phi_Lead",true);
			tree->SetBranchStatus("Muons_PT_Sub",true);
			tree->SetBranchStatus("Muons_Eta_Sub",true);
			tree->SetBranchStatus("Muons_Phi_Sub",true);
			tree->SetBranchStatus("Event_MET",true);
			tree->SetBranchStatus("Event_MET_Phi",true);
			tree->SetBranchStatus("GlobalWeight",true);
			tree->SetBranchStatus("EventInfo_EventNumber",true);

			// tree->Print();
			int maxN = 2;
			int maxT = 3;
			int maxS = maxN*maxT+1;
			double CMS = 13000;

			// data for jets
			float_t jl_pt, jl_eta, jl_phi, jl_ee;
			float_t js_pt, js_eta, js_phi, js_ee;
			Float_t jl_id, js_id;
			// double jlpx=0, jlpy=0, jlpy=0, jlee=0;
			// double jspx=0, jspy=0, jspy=0; jsee=0;
			tree->SetBranchAddress("Jets_E_Lead",&jl_ee);
			tree->SetBranchAddress("Jets_PT_Lead",&jl_pt);
			tree->SetBranchAddress("Jets_Eta_Lead",&jl_eta);
			tree->SetBranchAddress("Jets_Phi_Lead",&jl_phi);
			tree->SetBranchAddress("Jets_E_Sub",&js_ee);
			tree->SetBranchAddress("Jets_PT_Sub",&js_pt);
			tree->SetBranchAddress("Jets_Eta_Sub",&js_eta);
			tree->SetBranchAddress("Jets_Phi_Sub",&js_phi);

			tree->SetBranchAddress("Jets_PartonID_Lead",&jl_id);
			tree->SetBranchAddress("Jets_PartonID_Sub",&js_id);

			// data for muons
			float_t ml_pt, ml_eta, ml_phi, m_mm;
			float_t ms_pt, ms_eta, ms_phi;
			// double mlpx=0, mlpy=0, mlpy=0, mlee=0;
			// double mspx=0, mspy=0, mspy=0; msee=0;
			tree->SetBranchAddress("Muons_Minv_MuMu",&m_mm);
			tree->SetBranchAddress("Muons_PT_Lead",&ml_pt);
			tree->SetBranchAddress("Muons_Eta_Lead",&ml_eta);
			tree->SetBranchAddress("Muons_Phi_Lead",&ml_phi);
			tree->SetBranchAddress("Muons_PT_Sub",&ms_pt);
			tree->SetBranchAddress("Muons_Eta_Sub",&ms_eta);
			tree->SetBranchAddress("Muons_Phi_Sub",&ms_phi);

			// data for missing
			float_t met, met_phi, mc_wt;
			tree->SetBranchAddress("Event_MET",&met);
			tree->SetBranchAddress("Event_MET_Phi",&met_phi);
			tree->SetBranchAddress("GlobalWeight",&mc_wt);

			// data for events
			ULong64_t evenum;
			// UInt_t ttH, VH, DiMuon;
			// tree->SetBranchAddress("PassesttHSelection",&ttH);
			// tree->SetBranchAddress("PassesVHSelection",&VH);
			// tree->SetBranchAddress("PassesDiMuonSelection",&DiMuon);
			tree->SetBranchAddress("EventInfo_EventNumber",&evenum);
			
			tree->BuildIndex("EventInfo_EventNumber", "0");

			Int_t nen = (Int_t)tree->GetEntries();
			double projsig[7][7]{};

			cout << "Total Evenets:\t" << nen << endl; 
			cout << "Indexed Events:\t" << inen << endl;
			// cout << "Sample Events:\t" << loopnum <<endl;
			
			int n = 0;
			for (int i=0;i<inen;i++) {
				// create empty variables				
				TLorentzVector missing;
				vector<TLorentzVector> jets(maxT), bjets(maxT), muons(maxT);
				int bnum=0;				

				// int randnum = i+nen%(i+1+rand()%100);
				tree->GetEntryWithIndex(eventnumber[i]);
				if (m_mm>120 && m_mm<130 && ml_pt>27 && ms_pt>15) {
					if (jl_pt!=-9999 && jl_eta!=-9999 && jl_phi!=-9999 
							&& jl_ee!=-9999 && js_pt!= -9999 && js_eta!=-9999 
							&& js_phi!=-9999 && js_ee!=-9999) {
						
						/* Adjust the conditions for different channels. 
						Here is an example of ggF channel, where >=1 b-jet is required.*/
						if (jl_id==5 || js_id==5) {

							if (jl_id==5) bjets[0].SetPtEtaPhiE(jl_pt,jl_eta,jl_phi,jl_ee);
							else jets[0].SetPtEtaPhiE(jl_pt,jl_eta,jl_phi,jl_ee);

							if (js_id==5) bjets[1].SetPtEtaPhiE(js_pt,js_eta,js_phi,js_ee);
							else jets[1].SetPtEtaPhiE(js_pt,js_eta,js_phi,js_ee);
							
							if (jl_id==5) bnum+=1;
							if (js_id==5) bnum+=1;

							muons[0].SetPtEtaPhiM(ml_pt,ml_eta,ml_phi,MMU);
							muons[1].SetPtEtaPhiM(ms_pt,ms_eta,ms_phi,MMU);
							
							missing.SetPtEtaPhiM(met,0,met_phi,0);
							
							sort(jets.begin(),jets.end(), compare);
							sort(bjets.begin(),bjets.end(), compare);
							sort(muons.begin(),muons.end(), compare);
						
							double** projsep = map2rmm(CMS,maxN,maxT,missing,jets,bjets,muons);

							/* Set Mjj cells of jets to 0
							projsep[1][1] = 0;
							projsep[1][1+maxN] = 0;
							projsep[1+maxN][1+maxN] = 0;
							projsep[1][2] = 0;
							projsep[1+maxN][2+maxN] = 0;
							*/						
							
							// oset << bnum << " ";
							for (int x=0;x<maxS;x++) {
								for (int y=0;y<maxS;y++) {
									float dia = (bnum+0.5)/3;
									if (x==y) oset << dia << " ";
									else oset << projsep[x][y] << " ";
								}
							}

							oset << endl;
							// oset << mc_wt << endl;

							n++;
							nsig++;
							wt += mc_wt;
							if (nsig%10000==0) cout << "Finished " << nsig <<" events..."<<endl;						
							if (nsig>=setnum) goto endpart; // disable it if you want to loop over all files
						}						
					}
			
				}
	
			}

		cout << "Effective Events: " << n << endl;			
		}
	}

	endpart:
		cout << "Total Weight: " << wt << endl;
		oset.close();
}
