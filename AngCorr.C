// Angular correlations and mixing ratios follow the definitions in
// K.S. Krane and R.M. Steffen, Phys. Rev. C 2, 724 (1970)
// DOI: https://doi.org/10.1103/PhysRevC.2.724 

// L1f + L2f -> multipolarities of first transition
// L1s + L2s -> multipolarities of second transition
// J1, J2, J3 -> spins of first, intermediate, last levels
// delta1, delta2 -> mixing ratios of first and second transitions
// Q2, Q4 -> Q coefficients accounting for detector geometry

double J1, J2, J3;
double L1f, L2f, L1s, L2s;
double delta1, delta2;

bool gamma_cascade_defined = false;

double Q2 = 0.948; 
double Q4 = 0.862;
double delta_min=-50;
double delta_max=50;
double delta_step=0.01;



cout << "\nUsage of this root macro:\n" << endl;
cout << "SetGammaCascade(J1, J2, J3, L1f, L2f, L1s, L2s, delta1, delta2)\n";
cout << "\tJ1, J2, J3       ---> level spins in the order\n";
cout << "\tL1f, L2f, delta1 ---> multipolarities and mixing ratio of the first gamma\n";
cout << "\tL1s, L2s, delta2 ---> multipolarities and mixing ratio of the second gamma\n\n";
cout << "GetDelta(datafile.txt, whichdelta [1|2], arctan [0|1], logy [0|1])\n";
cout << "\tdatafile.txt     ---> txt file with cos(theta), W, Werr\n";
cout << "\twhichdelta       ---> [1|2] minimize delta1 or delta2\n";
cout << "\tarctan           ---> [0|1] plot chi2 vs arctan(delta)\n";
cout << "\tlogy             ---> [0|1] set log scale on chi2\n\n";



//===============================================================================================
// AVAILABLE METHODS
//===============================================================================================

// Print SetGammaCascade(...) usage
void PrintSetGammaCascadeUsage() {

	cout << "SetGammaCascade(J1, J2, J3, L1f, L2f, L1s, L2s, delta1, delta2)\n";
	cout << "\tJ1, J2, J3       ---> level spins in the order\n";
	cout << "\tL1f, L2f, delta1 ---> multipolarities and mixing ratio of the first gamma\n";
	cout << "\tL1s, L2s, delta2 ---> multipolarities and mixing ratio of the second gamma\n\n";	
}



// Define the gamma cascade
void SetGammaCascade(double j1, double j2, double j3,
                     double l1f, double l2f, double l1s, double l2s,
                     double d1, double d2) {
    
    J1 = j1; J2 = j2; J3 = j3;
    L1f = l1f; L2f = l2f; L1s = l1s; L2s = l2s;
    delta1 = d1; delta2 = d2;
	gamma_cascade_defined = true;

	cout << endl;
	cout << "-------------------- J = " << J1 << endl;
	cout << "          |          " << endl;
	cout << "  Gamma 1 | L = " << L1f << "+" << L2f << endl;
	cout << "          |          " << endl;
	cout << "         \\|/         " << endl;
	cout << "-------------------- J = " << J2 << endl;
	cout << "          |          " << endl;
	cout << "  Gamma 2 | L = " << L1s << "+" << L2s << endl;
	cout << "          |          " << endl;
	cout << "         \\|/         " << endl;
    cout << "-------------------- J = " << J3 << endl << endl;
}



// Change Q values accounting for detector geometry
void SetQvalues(double q2, double q4) { 
   
    Q2 = q2; Q4 =q4;
	cout << endl;
	cout << "Set Q2 = " << Q2 << endl;
	cout << "Set Q4 = " << Q4 << endl << endl;
}



// Change delta search range and step
void SetDeltaRange(double dmin, double dmax, double dstep) {

	delta_min = dmin;
	delta_max = dmax;
	delta_step = dstep;
	cout << endl;	
	cout << "Scanning delta between " << dmin << " and " << dmax << ", step = " << dstep;	
	cout << endl << endl;
}



// Evaluate Wigner 3-j symbol
double Symbol3j(double j1, double j2, double j3, double m1, double m2, double m3) {

    int L1 = int(round(2*j1));
    int L2 = int(round(2*j2));
    int L3 = int(round(2*j3));
    int M1 = int(round(2*m1));
    int M2 = int(round(2*m2));
    int M3 = int(round(2*m3));
    double symbol3j = ROOT::Math::wigner_3j(L1, L2, L3, M1, M2, M3);
    return symbol3j;
}



// Evaluate Wigner 6-j symbol
double Symbol6j(double j1, double j2, double j3, double m1, double m2, double m3) {

    int L1 = int(round(2*j1));
    int L2 = int(round(2*j2));
    int L3 = int(round(2*j3));
    int M1 = int(round(2*m1));
    int M2 = int(round(2*m2));
    int M3 = int(round(2*m3));
    double symbol6j = ROOT::Math::wigner_6j(L1, L2, L3, M1, M2, M3);
    return symbol6j;
}



// Evaluate Fk coefficient
double Fcoeff(double k, double l1, double l2, double j1, double j2) {

    double symbol3j = Symbol3j(l1, l2, k, 1, -1, 0);
    double symbol6j = Symbol6j(l1, l2, k, j2, j2, j1);
    int phase_exp = int(j1 + j2 - 1);
    double phase_sign = (phase_exp % 2 == 0) ? 1.0 : -1.0;
    double F = phase_sign * sqrt((2.*k+1.) * (2.*l1+1.) * (2.*l2+1.) * (2.*j2+1.)) * symbol3j * symbol6j;
    return F;
}



// Evaluate orientation coefficient Bk
double Bcoeff(double k, double l1, double l2, double j1, double j2, double delta) {

    int phase_exp = int(l1 + l2);
    double phase_sign = (phase_exp % 2 == 0) ? 1.0 : -1.0;
    double B = ( Fcoeff(k, l1, l1, j1, j2) +
                 phase_sign * 2. * delta * Fcoeff(k, l1, l2, j1, j2) +
                 delta*delta * Fcoeff(k, l2, l2, j1, j2) ) / (1. + delta*delta);
    return B;
}



// Evaluate distribution coefficient Ak
double Acoeff(double k, double l1, double l2, double j2, double j3, double delta) {

    double A = ( Fcoeff(k, l1, l1, j3, j2) +
                 2. * delta * Fcoeff(k, l1, l2, j3, j2) +
                 delta*delta * Fcoeff(k, l2, l2, j3, j2) ) / (1. + delta*delta);

    return A;
}



// Evaluate correlation coefficient Akk
double Akkcoeff(double k, double l1f, double l2f, double l1s, double l2s,
                double j1, double j2, double j3, double delta1, double delta2) {

    double Akk = Bcoeff(k, l1f, l2f, j1, j2, delta1) * Acoeff(k, l1s, l2s, j2, j3, delta2);
    return Akk;
}



// Main function which runs the chi2 minimization to get the mixing ratio
void GetDelta(const string& datafile, int whichdelta, bool arctan=1, bool logy=1) {
	
	if (!gamma_cascade_defined) {
		cerr << "\nERROR: gamma cascade not defined";
		cerr << "\nSet it through SetGammaCascade(J1, J2, J3, L1f, L2f, L1s, L2s, delta1, delta2)\n\n";
		PrintSetGammaCascadeUsage();
		return;
	}

	// Read cos(theta), W, Werr from file
    TCanvas* can = new TCanvas("can", "AngCorr", 1200, 500);
	can->Divide(2,1);
    TGraphErrors* graph = new TGraphErrors();
    graph->SetName("graph");
    
	ifstream in;
    in.open(datafile);
    if (!in.is_open()) {
        cerr << "Error while opening data file " << datafile << endl;
        return;
    }

    int i = 0;
    double theta, w, werr;
    while (in >> theta >> w >> werr) {
        graph->SetPoint(i, theta, w);
        graph->SetPointError(i, 0., werr);
        ++i;
    }
    in.close();



    // Fit with A0(1 + A2*Q2*P2(cos(theta)) + A4*Q4*P4(cos(theta)))
	TGraph* chi2Scan = new TGraph();
    chi2Scan->SetName("chi2Scan");
    
	TF1* func = new TF1("func", "[0]*(1. + [1]*[2]*ROOT::Math::legendre(2, x) + [3]*[4]*ROOT::Math::legendre(4, x))", -1., 1.);
    func->FixParameter(2, Q2);
    func->FixParameter(4, Q4);

    double* X = graph->GetX();
    double* Y = graph->GetY();
    double* Xerr = graph->GetEX();
    double* Yerr = graph->GetEY();

    int k = 0;
    double chi2_bestglob = 1e+9;
    double delta_bestglob;
    double A0_best, A22_best, A44_best;
    int ndof = graph->GetN() -1;

    for (double delta = delta_min; delta <= delta_max; delta+=delta_step) {
        double chi2 = 0.;
        if (whichdelta==1) delta1=delta;		// Minimize delta of first transition
        else if (whichdelta==2) delta2=delta;	// Minimize delta of second transition
        double A22 = Akkcoeff(2., L1f, L2f, L1s, L2s, J1, J2, J3, delta1, delta2);
        double A44 = Akkcoeff(4., L1f, L2f, L1s, L2s, J1, J2, J3, delta1, delta2);
        func->SetParameter(0, 1.); 
        func->FixParameter(1, A22);
        func->FixParameter(3, A44);
        graph->Fit(func, "NRQ");				// Optimize the A0 parameter only
        for (int i=0; i<graph->GetN(); ++i) 
            chi2 += pow( (Y[i] - func->Eval(X[i])) / Yerr[i], 2);
        chi2Scan->SetPoint(k, delta, chi2);
        if (chi2 < chi2_bestglob) {
            chi2_bestglob = chi2;
            delta_bestglob = delta;
            A0_best = func->GetParameter(0);
            A22_best = A22;
            A44_best = A44;
        }
        ++k;
    }

	cout << endl;
	cout << "########################################" << endl;
	cout << "Fit parameters (chi2 global minimum):" << endl;
	cout << "A0  = " << A0_best << endl;
	cout << "A22 = " << A22_best << endl;
    cout << "A44 = " << A44_best << endl;
	cout << "########################################" << endl;



	// Interpolate chi2Scan to get minima position
	vector<double> minima;
	TSpline3* spline = new TSpline3("spline", chi2Scan);
	for (double x=delta_min; x<delta_max; x+=delta_step) {
		double x1 = x;
		double x2 = x+delta_step;
		double d1 = spline->Derivative(x1);
		double d2 = spline->Derivative(x2);
		if (d1*d2<0) {
			double x0 = (x1+x2)/2;
			double dx = 1e-3;
			double d2 = (spline->Derivative(x0+dx) - spline->Derivative(x0-dx)) / (2*dx);
			if (d2>0) 
				minima.push_back(x0);
		}
	}



	// Fit near chi2Scan minima to get delta with errors
	int num_minima = minima.size();
	for (int i=0; i<num_minima; ++i) {
		double delta = minima[i];
		double chi2 = spline->Eval(minima[i]);
		double target = chi2+1.;
		double x_left = delta_min;
		double x_right = delta_max;
		for (int j=0; j<chi2Scan->GetN(); ++j) {
			double x, y;
			chi2Scan->GetPoint(j, x, y);
			if (x<minima[i] && y>=target && x>x_left) x_left = x;
			if (x>minima[i] && y>=target && x<x_right) x_right = x;
		}
		double err_left = delta - x_left;
		double err_right = x_right - delta;
		cout << "DELTA = " << delta << " (+" << err_right << ", -" << err_left << ")" << endl;
		cout << "CHI2 = " << chi2 << ", NDOF = " << ndof << endl;
		cout << "P-VALUE = " << TMath::Prob(chi2, ndof) << endl;
		cout << "########################################" << endl;
	}
	cout << endl;



	// Plot the experimental point with angular correlation fit
    can->SetCrosshair();
    can->cd(1);
    graph->SetTitle("Angular correlation fit");
	graph->GetYaxis()->SetTitle("W(#theta_{#gamma#gamma})");
	graph->GetXaxis()->SetTitle("cos(#theta_{#gamma#gamma})");
    graph->SetMarkerStyle(20);
    graph->SetMarkerColor(kBlue);
    func->SetParameter(0, A0_best);
    func->SetParameter(1, A22_best);
    func->SetParameter(3, A44_best);
    func->SetLineWidth(2);
    func->SetLineColor(kRed);
    gPad->SetGridx();
	gPad->SetGridy();
	graph->Draw("AP");
    func->Draw("same");



	// Plot the chi2 scan as a function of delta (or ATan(delta))
	if (arctan) {
		double x, y;
		for (int i=0; i<chi2Scan->GetN(); ++i) {
			chi2Scan->GetPoint(i, x, y);
			chi2Scan->SetPoint(i, TMath::ATan(x), y);
		}
	} 
	can->cd(2);
    chi2Scan->SetTitle("#chi^{2} scan");
	chi2Scan->GetYaxis()->SetTitle("#chi^{2}");
	if (!arctan) chi2Scan->GetXaxis()->SetTitle("#delta");
	if (arctan) chi2Scan->GetXaxis()->SetTitle("tan^{-1}(#delta)");
    if (logy) gPad->SetLogy();
	chi2Scan->SetLineWidth(2);
    chi2Scan->SetLineColor(kBlack);
	gPad->SetGridx();
	gPad->SetGridy();
    chi2Scan->Draw("APL");
	
}
