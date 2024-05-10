struct RBF_fn {
	virtual double rbf(double r) = 0;
};

struct RBF_interp {
	int dim, n;
	const Matdouble &pts;
	const Vecdouble &vals;
//	Vecdouble w;
	double* w = new double[n];
	RBF_fn &fn;
	Bool norm;

	RBF_interp(Matdouble_I &ptss, Vecdouble_I &valss, RBF_fn &func, Bool nrbf=false)
	: dim(ptss.ncols()), n(ptss.nrows()) , pts(ptss), vals(valss),
	w(n), fn(func), norm(nrbf) {
		int i,j;
		double sum;
		Matdouble rbf(n,n);
		Vecdouble rhs(n);
		for (i=0;i<n;i++) {
			sum = 0.;
			for (j=0;j<n;j++) {
				sum += (rbf[i][j] = fn.rbf(rad(&pts[i][0],&pts[j][0])));
			}
			if (norm) rhs[i] = sum*vals[i];
			else rhs[i] = vals[i];
		}
		LUdcmp lu(rbf);
		lu.solve(rhs,w);
	}

	double interp(Vecdouble_I &pt) {
		double fval, sum=0., sumw=0.;
		if (pt.size() != dim) throw("RBF_interp bad pt size");
		for (int i=0;i<n;i++) {
			fval = fn.rbf(rad(&pt[0],&pts[i][0]));
			sumw += w[i]*fval;
			sum += fval;
		}
		return norm ? sumw/sum : sumw;
	}

	double rad(const double *p1, const double *p2) {
		double sum = 0.;
		for (int i=0;i<dim;i++) sum += SQR(p1[i]-p2[i]);
		return sqrt(sum);
	}
};
struct RBF_multiquadric : RBF_fn {
	double r02;
	RBF_multiquadric(double scale=1.) : r02(SQR(scale)) {}
	double rbf(double r) { return sqrt(SQR(r)+r02); }
};

struct RBF_thinplate : RBF_fn {
	double r0;
	RBF_thinplate(double scale=1.) : r0(scale) {}
	double rbf(double r) { return r <= 0. ? 0. : SQR(r)*log(r/r0); }
};

struct RBF_gauss : RBF_fn {
	double r0;
	RBF_gauss(double scale=1.) : r0(scale) {}
	double rbf(double r) { return exp(-0.5*SQR(r/r0)); }
};

struct RBF_inversemultiquadric : RBF_fn {
	double r02;
	RBF_inversemultiquadric(double scale=1.) : r02(SQR(scale)) {}
	double rbf(double r) { return 1./sqrt(SQR(r)+r02); }
};
struct Shep_interp {
	int dim, n;
	const Matdouble &pts;
	const Vecdouble &vals;
	double pneg;

	Shep_interp(Matdouble_I &ptss, Vecdouble_I &valss, double p=2.)
	: dim(ptss.ncols()), n(ptss.nrows()) , pts(ptss),
	vals(valss), pneg(-p) {}

	double interp(Vecdouble_I &pt) {
		double r, w, sum=0., sumw=0.;
		if (pt.size() != dim) throw("RBF_interp bad pt size");
		for (int i=0;i<n;i++) {
			if ((r=rad(&pt[0],&pts[i][0])) == 0.) return vals[i];
			sum += (w = pow(r,pneg));
			sumw += w*vals[i];
		}
		return sumw/sum;
	}

	double rad(const double *p1, const double *p2) {
		double sum = 0.;
		for (int i=0;i<dim;i++) sum += SQR(p1[i]-p2[i]);
		return sqrt(sum);
	}
};
